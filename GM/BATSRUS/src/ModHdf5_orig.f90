module ModHdf5

  use ModProcMH, ONLY: iProc, nProc, iComm 
  use hdf5
  use ModHdf5Utils

  implicit none
  SAVE

  ! These parameters are needed if there is no HDF5 module
  !integer,parameter:: Hsize_t=4,hid_t=4,H5T_NATIVE_INTEGER=4,H5T_NATIVE_CHARACTER=1
  !integer,parameter:: H5P_DEFAULT_F=8,H5S_ALL_F=8,H5FD_MPIO_COLLECTIVE_F=8,H5S_SELECT_SET_F=8
  !integer,parameter:: H5P_DATASET_XFER_F=8,H5T_NATIVE_DOUBLE=8,H5F_ACC_TRUNC_F=8,H5P_FILE_ACCESS_F=1

  private ! except
  public write_var_hdf5
  public write_cut_var_hdf5
  public write_plot_hdf5
  public hdf5_sph_plot

  real,    allocatable :: PlotVarIdx(:,:,:,:,:)
  real,    allocatable :: plotXYZNodes(:,:,:,:,:)

  integer, allocatable :: UsedBlocks(:)
  integer, allocatable :: usednodes(:)
  integer, allocatable :: UsedBlocksTemp(:)

  integer, parameter:: lNameVar = 10
  integer, parameter:: lNameH5 = lNameVar + 1

  character (len=lNameH5), allocatable :: UnknownNameArray(:)

  integer :: nBlocksUsed, nBlkUsedGlobal, iBlk ,arrSize(3), nodeArrSize(3)
  integer :: iSizeGlobal, jSizeGlobal, kSizeGlobal, nPlotDim, plotDim(3) 
  logical :: collectiveWrite
contains
  subroutine  hdf5_sph_plot(i,j,x,y,z,theta_out,phi_out,PointVar)
    real, intent(in) :: x,y,z, theta_out, phi_out,PointVar(:)
    integer, intent(in) :: i,j 

  end subroutine hdf5_sph_plot
  subroutine write_cut_var_hdf5(iFile, plotType, iBlock, H5Index,nPlotVar,PlotVar, &
       xMin, xMax, yMin, yMax, zMin, zMax, DxBlock, DyBlock, DzBlock,&
       isNonCartesian,nCell,H5Advance)

    ! Save all cells within plotting range, for each processor

    use ModProcMH
    use ModMain, ONLY: nI, nJ, nK, &
         x_, y_, z_, Phi_, nBlockMax, unusedBLK
    use ModGeometry, ONLY: Dx_BLK, Dy_BLK, Dz_BLK,&
         XyzStart_BLK
    use ModIO
    use ModNumConst
    use ModCovariant, ONLY: is_axial_geometry          
    use ModMpi
    use BATL_lib, ONLY : Xyz_DNB

    implicit none

    ! Arguments

    integer, intent(in)   :: iFile, iBlock, H5Index
    integer, intent(in)   :: nPlotVar
    real,    intent(in)   :: PlotVar(-1:nI+2,-1:nJ+2,-1:nK+2,nPlotVar)
    real,    intent(in)   :: xMin,xMax,yMin,yMax,zMin,zMax
    logical, intent(in) :: isNonCartesian
    character (len=*), intent(in) :: plotType
    real,    intent(inout):: DxBlock,DyBlock,DzBlock
    integer, intent(out)  :: nCell 
    logical, intent(out) :: H5advance

    ! Local variables
    ! Indices and coordinates
    integer :: i, j, k, iMin, iMax, jMin, jMax, kMin, kMax
    integer :: iMin1, iMax1, jMin1, jMax1, kMin1, kMax1
    integer ::  iMaxN,jMaxN, kMaxN,iMaxN1,jMaxN1, kMaxN1
    integer :: ii, jj, kk, iSize, jSize, kSize
    integer :: error
    real ::Dx
    real :: xMin1,xMax1,yMin1,yMax1,zMin1,zMax1
    real :: DxBlockOut

    real :: ySqueezed

    real, parameter:: cHalfMinusTiny=cHalf*(cOne-cTiny)

    character(len=*), parameter :: NameSub = 'write_plot_idl'


    if (.not. allocated(PlotVarIdx)) then
       iMin = huge(iMin)
       iMax = -huge(iMax)
       jMin = huge(jMin)
       jMax = -huge(jMax)
       kMin = huge(kMin)
       kMax = -huge(kMax)
       iSize = 0
       jSize = 0
       kSize = 0
       do iBlk = 1, nBlockMax
          xMin1 = xMin - cHalfMinusTiny*Dx_BLK(iBlk)
          xMax1 = xMax + cHalfMinusTiny*Dx_BLK(iBlk)
          yMin1 = yMin - cHalfMinusTiny*Dy_BLK(iBlk)
          yMax1 = yMax + cHalfMinusTiny*Dy_BLK(iBlk)
          zMin1 = zMin - cHalfMinusTiny*Dz_BLK(iBlk)
          zMax1 = zMax + cHalfMinusTiny*Dz_BLK(iBlk)

          if(is_axial_geometry())then                 
             ! Make sure that angles around 3Pi/2 are moved to Pi/2 for x=0 cut
             ySqueezed = mod(xyzStart_BLK(Phi_,iBlk),cPi)
             ! Make sure that small angles are moved to Pi degrees for y=0 cut
             if(ySqueezed < 0.25*cPi .and. &
                  abs(yMin+yMax-cTwoPi) < cTiny .and. yMax-yMin < 0.01) &
                  ySqueezed = ySqueezed + cPi
          else                                          
             ySqueezed = xyzStart_BLK(y_,iBlk)
          end if

          if( unusedBLK(iBLK) .or. xyzStart_BLK(x_,iBlk) > xMax1.or.&
               xyzStart_BLK(x_,iBlk)+(nI-1)*Dx_BLK(iBlk) < xMin1.or.&
               ySqueezed > yMax1.or.&
               ySqueezed+(nJ-1)*Dy_BLK(iBlk) < yMin1.or.&  
               xyzStart_BLK(z_,iBlk) > zMax1.or.&
               xyzStart_BLK(z_,iBlk)+(nK-1)*Dz_BLK(iBlk) < zMin1) then

             CYCLE

          else

             Dx = plot_Dx(1,iFile)
             DxBlock=Dx_BLK(iBlk); DyBlock=Dy_BLK(iBlk); DzBlock=Dz_BLK(iBlk)

             ! Calculate index limits of cells inside cut
             i = max(1 ,floor((xMin1-xyzStart_BLK(x_,iBlk))/DxBlock)+2)
             ii = min(nI,floor((xMax1-xyzStart_BLK(x_,iBlk))/DxBlock)+1)

             j = max(1 ,floor((yMin1-ySqueezed)/DyBlock)+2)
             jj = min(nJ,floor((yMax1-ySqueezed)/DyBlock)+1)

             k = max(1 ,floor((zMin1-xyzStart_BLK(z_,iBlk))/DzBlock)+2)
             kk = min(nK,floor((zMax1-xyzStart_BLK(z_,iBlk))/DzBlock)+1)

             if (i < iMin)&
                  iMin = i
             if (ii > iMax)&
                  iMax = ii
             if ((ii-i + 1) > iSize)&
                  iSize = ii-i+1

             if (j < jMin)&
                  jMin = j
             if (jj > jMax)&
                  jMax = jj
             if ((jj-j+1) > jSize)&
                  jSize = jj-j+1

             if (k < kMin)&
                  kMin = k
             if (kk > kMax)&
                  kMax = kk
             if ((kk-k+1) > kSize)&
                  kSize = kk-k+1

          end if
       end do
       call MPI_AllReduce(iSize, iSizeGlobal, 1,&
            MPI_INTEGER, MPI_MAX, iComm, error)
       call MPI_AllReduce(jSize, jSizeGlobal, 1,&
            MPI_INTEGER, MPI_MAX, iComm, error)
       call MPI_AllReduce(kSize, kSizeGlobal, 1,&
            MPI_INTEGER, MPI_MAX, iComm, error)
       nBlocksUsed = 0
       allocate(UsedBlocksTemp(nBlockMax))
       if (isNonCartesian) then
          if (plotType=='x=0') then
             !              write (*,*) "allocated x=0"
             nPlotDim = 2
             plotDim = (/1,3,0/)
             arrSize = (/iSizeGlobal, kSizeGlobal, 1/)
             nodeArrSize = (/iSizeGlobal+1, kSizeGlobal+1, 1/)
          else if(plotType=="y=0)") then
             !              write (*,*) "allocated y=0"
             nPlotDim = 2
             plotDim = (/1,3,0/)
             arrSize = (/iSizeGlobal, kSizeGlobal, 1/)
             nodeArrSize = (/iSizeGlobal+1, kSizeGlobal+1, 1/)
          else if (plotType=='z=0' )then
             !              write (*,*) "allocated z=0"
             nPlotDim = 2
             plotDim = (/1,2,0/)
             arrSize = (/iSizeGlobal, jSizeGlobal, 1/)
             nodeArrSize = (/iSizeGlobal+1, jSizeGlobal+1, 1/)
          end if
       else
          if (iSizeGlobal == 1 .and. jSizeGlobal == 1) then
             nPlotDim = 1
             plotDim = (/3,0,0/)
             arrSize = (/kSizeGlobal,1,1/)
             nodeArrSize = (/kSizeGlobal+1,1,1/)
          else if (iSizeGlobal == 1 .and. kSizeGlobal == 1) then
             nPlotDim = 1
             plotDim = (/2,0,0/)
             arrSize = (/jSizeGlobal,1,1/)
             nodeArrSize = (/jSizeGlobal+1,1,1/)
          else if (jSizeGlobal == 1 .and. kSizeGlobal == 1) then
             nPlotDim = 1
             plotDim = (/1,0,0/)
             arrSize = (/iSizeGlobal,1,1/)
             nodeArrSize = (/iSizeGlobal+1,1,1/)
          else if (iSizeGlobal == 1) then
             nPlotDim = 2
             plotDim = (/2,3,0/)
             arrSize = (/jSizeGlobal,kSizeGlobal,1/)
             nodeArrSize = (/jSizeGlobal+1,kSizeGlobal+1,1/)
          else if (jSizeGlobal == 1) then
             nPlotDim = 2
             plotDim = (/1,3,0/)
             arrSize = (/iSizeGlobal,kSizeGlobal,1/)
             nodeArrSize = (/iSizeGlobal+1,kSizeGlobal+1,1/)
          else if (kSizeGlobal == 1) then
             nPlotDim = 2
             plotDim = (/1,2,0/)
             arrSize = (/iSizeGlobal,jSizeGlobal,1/)
             nodeArrSize = (/iSizeGlobal+1,jSizeGlobal+1,1/)
          else
             nPlotDim = 3
             plotDim = (/1,2,3/)
             arrSize = (/iSizeGlobal,jSizeGlobal,kSizeGlobal/)
             nodeArrSize = (/iSizeGlobal+1,jSizeGlobal+1,kSizeGlobal+1/)
          end if
       end if
       allocate(PlotVarIdx(arrSize(1),arrSize(2), arrSize(3), nBlockMax, nPlotVar))
       PlotVarIdx=-50
       if (isNonCartesian) allocate(plotXYZNodes(nPlotDim, &
            nodeArrSize(1),nodeArrSize(2),nodeArrSize(3),nBlockMax))
    end if

    xMin1 = xMin - cHalfMinusTiny*Dx_BLK(iBlock)
    xMax1 = xMax + cHalfMinusTiny*Dx_BLK(iBlock)
    yMin1 = yMin - cHalfMinusTiny*Dy_BLK(iBlock)
    yMax1 = yMax + cHalfMinusTiny*Dy_BLK(iBlock)
    zMin1 = zMin - cHalfMinusTiny*Dz_BLK(iBlock)
    zMax1 = zMax + cHalfMinusTiny*Dz_BLK(iBlock)

    nCell = 0
    if(isNonCartesian)then                 
       ! Make sure that angles around 3Pi/2 are moved to Pi/2 for x=0 cut
       ySqueezed = mod(xyzStart_BLK(Phi_,iBlock),cPi)
       ! Make sure that small angles are moved to Pi degrees for y=0 cut
       if(ySqueezed < 0.25*cPi .and. &
            abs(yMin+yMax-cTwoPi) < cTiny .and. yMax-yMin < 0.01) &
            ySqueezed = ySqueezed + cPi
    else                                          
       ySqueezed = xyzStart_BLK(y_,iBlock)
    end if
    ! 
    !   if(DoTestMe)then
    !      write(*,*) NameSub, 'xMin1,xMax1,yMin1,yMax1,zMin1,zMax1=',&
    !           xMin1,xMax1,yMin1,yMax1,zMin1,zMax1
    !      write(*,*) NameSub, 'xyzStart_BLK=',iBlock,xyzStart_BLK(:,iBlock)
    !      write(*,*) NameSub, 'ySqueezed =',ySqueezed
    !      write(*,*) NameSub, 'xyzEnd=',xyzStart_BLK(x_,iBlock)+(nI-1)*Dx_BLK(iBlock),&
    !           ySqueezed+(nJ-1)*Dy_BLK(iBlock),&
    !           xyzStart_BLK(z_,iBlock)+(nK-1)*Dz_BLK(iBlock)
    !   end if
    ! 
    ! If block is fully outside of cut then cycle
    if(  xyzStart_BLK(x_,iBlock) > xMax1.or.&
         xyzStart_BLK(x_,iBlock)+(nI-1)*Dx_BLK(iBlock) < xMin1.or.&
         ySqueezed > yMax1.or.&
         ySqueezed+(nJ-1)*Dy_BLK(iBlock) < yMin1.or.&  
         xyzStart_BLK(z_,iBlock) > zMax1.or.&
         xyzStart_BLK(z_,iBlock)+(nK-1)*Dz_BLK(iBlock) < zMin1) then
       H5Advance = .false.
       RETURN
    end if

    Dx = plot_Dx(1,iFile)
    DxBlock=Dx_BLK(iBlock); DyBlock=Dy_BLK(iBlock); DzBlock=Dz_BLK(iBlock)

    ! Calculate index limits of cells inside cut
    iMin = max(1 ,floor((xMin1-xyzStart_BLK(x_,iBlock))/DxBlock)+2)
    iMax = min(nI,floor((xMax1-xyzStart_BLK(x_,iBlock))/DxBlock)+1)

    jMin = max(1 ,floor((yMin1-ySqueezed)/DyBlock)+2)
    jMax = min(nJ,floor((yMax1-ySqueezed)/DyBlock)+1)

    kMin = max(1,floor((zMin1-xyzStart_BLK(z_,iBlock))/DzBlock)+2)
    kMax = min(nK,floor((zMax1-xyzStart_BLK(z_,iBlock))/DzBlock)+1)

    iMin1 = max(0 ,floor((xMin1-xyzStart_BLK(x_,iBlock))/DxBlock)+2)
    iMax1 = min(nI+1,floor((xMax1-xyzStart_BLK(x_,iBlock))/DxBlock)+1)

    kMin1 = max(0,floor((zMin1-xyzStart_BLK(z_,iBlock))/DzBlock)+2)
    kMax1 = min(nK+1,floor((zMax1-xyzStart_BLK(z_,iBlock))/DzBlock)+1)

    jMin1 = max(0 ,floor((yMin1-ySqueezed)/DyBlock)+2)
    jMax1 = min(nJ+1,floor((yMax1-ySqueezed)/DyBlock)+1)

    iMaxN = iMax
    jMaxN = jMax
    kMaxN = kMax


    if (isNonCartesian) then
       if ( plotType(1:3)=='x=0' .or. plotType(1:3)=='y=0') then
          iMaxN = iMax + 1
          kMaxN = kMax + 1

          iMaxN1 = iMax1 + 1
          kMaxN1 = kMax1 + 1

       else if ( plotType(1:3)=='z=0') then
          iMaxN = iMax + 1
          jMaxN = jmax + 1
          iMaxN1 = iMax1 + 1
          jMaxN1 = jmax1 + 1


       end if
    else
       do ii = 1, nPlotDim
          if (plotDim(ii) == 1) then
             iMaxN = iMax + 1
          else if (plotDim(ii) == 2) then
             jMaxN = jMax + 1
          else if (plotDim(ii) == 3) then
             kMaxN = kMax + 1
          end if
       end do
    end if


    !   if(DoTestMe)then
    !      write(*,*) NameSub, 'iMin,iMax,jMin,jMax,kMin,kMax=',&
    !           iMin,iMax,jMin,jMax,kMin,kMax
    !      write(*,*) NameSub, 'DxBlock,x1,y1,z1',DxBlock,xyzStart_BLK(:,iBlock)
    !      write(*,*) NameSub, 'ySqueezed  =',ySqueezed
    !      write(*,*) NameSub, 'xMin1,xMax1=',xMin1,xMax1
    !      write(*,*) NameSub, 'yMin1,yMax1=',yMin1,yMax1
    !      write(*,*) NameSub, 'zMin1,zMax1=',zMin1,zMax1
    !   end if
    ! 
    DxBlockOut = DxBlock
    ! if (plot_dimensional(iFile))DxBlockOut = DxBlockOut*No2Io_V(UnitX_)
    !  if (nPlotDim == 2) then
    !     if (plotDim(1) == 2) then
    !         iMin = jMin
    !         iMax = jMax
    !     end if
    !     if (plotDim(2) == 3) then
    !         jMin = kMin
    !         jMax = kMax
    !     end if
    !     kMin = 1
    !     kMax = 1
    !   end if
    ! 

    ! this seems to do y=0 and z=0 plots okay but x=0 plots appear to be broken
    !I have only tried one test case so the problem be bigger than just x=0 cuts
    if (isNonCartesian) then
       do kk=kMin,kMaxN; do jj=jMin,jMaxN; do ii=iMin,iMaxN
          if (plotType=='x=0') then
             i = ii - iMin + 1
             j = kk - kMin + 1
             k = jj - jMin + 1

             plotXYZNodes(1,i, j, k, H5Index) = &
                  Xyz_DNB(2,ii, jj, kk, iBlock)

             plotXYZNodes(2,i, j, k, H5Index) = &
                  Xyz_DNB(3,ii, jj, kk, iBlock)
          else if (plotType=='y=0') then
             i = ii - iMin + 1
             j = kk - kMin + 1
             k = jj - jMin + 1

             plotXYZNodes(1,i, j, k, H5Index) = &
                  xyz_dnb(1,ii, jj, kk, iblock)
             plotXYZNodes(2,i, j, k, H5Index) = &
                  Xyz_DNB(3,ii, jj, kk, iBlock)
             ! 
          else if (plotType=='z=0') then
             i = ii - iMin + 1
             j = jj - jMin + 1
             k = kk - kMin + 1

             plotXYZNodes(1,i, j, k, H5Index) = &
                  Xyz_DNB(1,ii, jj, kk, iBlock)
             plotXYZNodes(2,i, j, k, H5Index) = &
                  Xyz_DNB(2,ii, jj, kk, iBlock)

          else if (nPlotDim == 3) then
             write (*,*) "NPLOTDIM==3"
          end if
          where(abs(PlotXYZNodes(:,:,:,:,H5Index))<1.e-10) &
               PlotXYZNodes(:,:,:,:,H5Index) = 0.
       end do; end do; end do
       do kk=kMin,kMax; do jj=jMin,jMax; do ii=iMin,iMax
          if ( plotType=='x=0' .or. plotType=='y=0') then
             !              write (*,*) "x=0 y=0"
             i = ii - iMin + 1
             j = kk - kMin + 1
             k = jj - jMin + 1
             PlotVarIdx(i, j, k, H5Index, 1:nPlotVar) = &
                  .5*(PlotVar(ii, jMin1, kk, 1:nPlotVar)+PlotVar(ii, jMax1, kk, 1:nPlotVar))
          else if ( plotType=='z=0') then
             !              write (*,*) 'z=0'
             i = ii - iMin + 1
             j = jj - jMin + 1
             k = kk - kMin + 1
             PlotVarIdx(i, j, k, H5Index, 1:nPlotVar) = &
                  .5*(PlotVar(ii, jj, kMin1, 1:nPlotVar)+PlotVar(ii, jj, kMax1, 1:nPlotVar))
          else if (nPlotDim == 3) then
             write (*,*) "NPLOTDIM==3"
          end if
          nCell = nCell+1
       end do; end do; end do
    else
       do kk=kMin,kMax; do jj=jMin,jMax; do ii=iMin,iMax

          if (nPlotDim == 2) then
             if (kSizeGlobal == 1) then
                i = ii - iMin + 1
                j = jj - jMin + 1
                k = kk - kMin + 1
                PlotVarIdx(i, j, k, H5Index, 1:nPlotVar) = &
                     .5*(plotvar(ii, jj, kmin1, 1:nplotvar)+plotvar(ii, jj, kmax1, 1:nplotvar))
             else if (jSizeGlobal == 1) then
                i = ii - iMin + 1
                j = kk - kMin + 1
                k = jj - jMin + 1
                PlotVarIdx(i, j, k, H5Index, 1:nPlotVar) = &
                     .5*(plotvar(ii, jMin1 ,kk, 1:nplotvar)+plotvar(ii, jMax1, kk, 1:nplotvar))
             else if (iSizeGlobal == 1) then
                i = jj - jMin + 1
                j = kk - kMin + 1
                k = ii - iMin + 1
                PlotVarIdx(i, j, k, H5Index, 1:nPlotVar) = &
                     .5*(plotvar(iMin1, jj ,kk, 1:nplotvar)+plotvar(iMax1, jj, kk, 1:nplotvar))
             end if
          else if (nPlotDim == 1) then
             if (plotDim(1) == 1) then
                i = ii - iMin + 1
                j = jj - jMin + 1
                k = kk - kMin + 1
                PlotVarIdx(i, j, k, H5Index, 1:nPlotVar) = &
                     plotvar(ii, jj ,kk, 1:nplotvar) 
             else if (plotDim(1) == 2) then
                i = jj - jMin + 1
                j = ii - iMin + 1
                k = kk - kMin + 1
                PlotVarIdx(i, j, k, H5Index, 1:nPlotVar) = &
                     plotvar(ii, jj ,kk, 1:nplotvar) 
             else if (plotDim(1) == 3) then
                i = kk - kMin + 1
                j = ii - iMin + 1
                k = kk - kMin + 1
                PlotVarIdx(i, j, k, H5Index, 1:nPlotVar) = &
                     plotvar(ii, jj ,kk, 1:nplotvar) 

             end if
          else
             i = ii - iMin + 1
             j = jj - jMin + 1
             k = kk - kMin + 1
             PlotVarIdx(i, j, k, H5Index, 1:nPlotVar) = &
                  plotvar(ii, jj ,kk, 1:nplotvar) 
          end if

          nCell = nCell+1
       end do; end do; end do
    end if
    H5Advance = .true.
    UsedBlocksTemp(H5Index) = iBlock
    nBlocksUsed = H5Index

  end subroutine write_cut_var_hdf5


  !=================================================================
  !================================================================


  subroutine write_var_hdf5(PlotVar, nPlotVar, H5Index, iBlock,&
       isNonCartesian)
    use ModMain, only : nBlockMax, nI,nJ,nK
    use BATL_lib, ONLY : nDim, Xyz_DNB


    !  Handles adding one block of variable data that for plotting
    !  in the next plotfile.

    integer, intent(in) :: nPlotVar, H5Index, iBlock
    real, intent(in) :: PlotVar(-1:nI+2,-1:nJ+2,-1:nK+2,nPlotVar)
    logical, intent(in) :: isNonCartesian
    integer:: ii, jj, kk

    !----------------------------------------------------------------------
    ! Allocate the storage array if it has not been allocated for
    ! this plotfile yet

!!! Use nBlock if possible
    !!NODE ORDER TEMPORARILY SWITCHED
    if (.not. allocated(PlotVarIdx)) then
       allocate(PlotVarIdx(nI, nJ, nK, nBlockMax, nPlotVar))
       allocate(UsedBlocksTemp(nBlockMax))
       iSizeGlobal = nI
       jSizeGlobal = nJ
       kSizeGlobal = nK
       arrSize = (/nI,nJ,nK/)
       nPlotDim = nDim
       plotDim = (/1,2,3/)
       nBlocksUsed = 0

       nodeArrSize(1) = nI + 1
       if (nDim .ge. 2) then
          nodeArrSize(2) = nJ + 1
       else
          nodeArrSize(2) = nJ
       end if
       if (nDim == 3) then
          nodeArrSize(3) = nK + 1
       else
          nodeArrSize(3) = nK
       end if
       if (isNonCartesian) &
            allocate(plotXYZNodes(nPlotDim,1:1+nI,1:1+nJ,1:1+nK,nBlockMax))
    end if
    if (isNonCartesian) then  
       do kk = 1, nodeArrSize(3); do jj = 1, nodeArrSize(2); do ii = 1, nodeArrSize(1)
          plotXYZNodes(1:nPlotDim,ii, jj, kk, H5Index) = &
               Xyz_DNB(1:nPlotDim,ii, jj, kk, iBlock)

       end do; end do; end do

    end if
    !translates iBlk to index for plotfile where unused blocks are skipped  
!!! Efficient loop ordering: do k = ...; do j = ...; do i = ...

    do kk = 1, nK; do jj = 1, nJ; do ii = 1, nI
       PlotVarIdx(ii, jj, kk, H5Index, 1:nPlotVar) = &
            PlotVar(ii,jj, kk,1:nPlotVar)


    end do; end do; end do
    UsedBlocksTemp(H5Index) = iBlock
    nBlocksUsed = H5Index


  end subroutine write_var_hdf5
  !===============================================================           
  !================================================================

  subroutine write_plot_hdf5(filename, plotType, plotVarNames, plotVarUnits,&
       nPlotVar,isCutFile, nonCartesian,plot_dimensional, xmin, xmax, &
       ymin, ymax, zmin, zmax)!, nBLKcells)

    use BATL_tree, only: iNode_B, iTree_IA, Coord0_,&
         Level_,iMortonNode_A
    use BATL_lib, only : CoordMin_DB, CoordMax_DB
    !  use ModIO
    use ModNumConst
    use ModMpi
    use BATL_lib, ONLY : CoordMin_DB, MaxDim
    use BATL_mpi, only : barrier_mpi


    integer,                 intent(in):: nPlotVar
    character(len=80),       intent(in):: filename
    character(len=lNameVar), intent(in):: plotVarNames(nPlotVar)
    character(len=lNameVar),  intent(in):: plotVarUnits(nPlotVar)
    character (len=*), intent(in) :: plotType
    real, intent(in)  ::  xMin, xMax, yMin, yMax, zMin, zMax
    !    integer,                 intent(in):: nBLKcells
    logical, intent(in) :: isCutFile, plot_dimensional, nonCartesian
    integer :: offsetPerProc(0:nProc-1), blocksPerProc(0:nProc-1)
    integer :: error, procIdx, iLen, iVar, nBlocksUsedMax 
    integer(HID_T) :: fileID
    integer :: labelLeng, i, nBlocksUsedMin
    integer :: hdfMajor, hdfMinor, hdfRelease
    real, allocatable :: coordinates(:,:), bbox(:,:,:),blockDataExtents(:,:)
    integer, allocatable :: procNum(:)
    real :: varMin, varMax
    real :: globalVarMin(nPlotVar), globalVarMax(nPlotVar)
    real :: minCoords(3), maxCoords(3)
    character (len=lnamevar+4) :: dsname
    real, parameter:: cHalfMinusTiny=cHalf*(cOne-cTiny)




    !--------------------------------------------------------------------------
    !     !If you make changes to the hdf part of this code, uncomment this and 
    !     !the writes at the end and run the code before commiting it to
    !     !make sure that it doesn't cause the writes to to be really slow.
    !     if (iProc == 0) write (*,*) '===================================='
    !     if (iProc == 0) write (*,*) 'Writing BATL hdf5 parallel plotfile'
    !     if (iProc == 0) write (*,*) '------------------------------------'
    !     if (iProc == 0) write (*,*) '  opening hdf5 file'

    call barrier_mpi()
    do iVar = 1, nPlotVar
       if (nBlocksUsed > 0) then
          varMin = minval(PlotVarIdx(:,:,:,1:nBlocksUsed,iVar))
          varMax = maxval(PlotVarIdx(:,:,:,1:nBlocksUsed,iVar))
       else
          varMin = huge(varMin)
          varMax = -huge(varMax)
       end if
       call MPI_allreduce(varMin, globalVarMin(iVar), 1, MPI_REAL,&
            MPI_MIN, iComm, error)

       call MPI_allreduce(varMax, globalVarMax(iVar), 1, MPI_REAL,&
            MPI_MAX, iComm, error)

    end do


    !initialize the Fortran Hdf5 library
    call h5open_f(error)

    !determine if we want to do collective writes. Collective write
    !is faster but all processors must call h5dwrite, even in cuts where 
    !some processors write no data.  In newer versions of hdf5 a null
    !write can be called but attempting to do so on older versions will
    !cause the code to crash. We should be able to do collective for any
    !hdf5 1.8.x where all procs write data but that didn't work for me 
    !for some reason so for now collective write mode is restricted to
    !hdf5 version 1.8.8+
    !     write (*,*) "write plot error 3"v
    call h5get_libversion_f(hdfMajor, hdfMinor, hdfRelease, error)   
    if (hdfmajor > 1) then
       collectiveWrite = .true.
    else if (hdfmajor == 1 .and. hdfMinor > 8) then
       collectiveWrite = .true.
    else if (hdfmajor == 1 .and. hdfMinor == 8) then
       if (hdfRelease .ge. 8 ) then!.or. allProcsWrite) then
          collectiveWrite = .true.
       else
          collectiveWrite = .false.
       end if
    else 
       collectiveWrite = .false.
    end if

    ! Open the hdf5 file and write data to it
    !---------------------------------------------------------------------    
    fileID = -1
    call hdf5_init_file(fileID, filename, iComm)
    if(fileID == -1) then
       if (iProc == 0) write (*,*)  "Error: unable to initialize file"
       call stop_mpi("unable to initialize hdf5 file")
    end if
    !get processor block offsets for writes
    call MPI_Allgather(nBlocksUsed, 1, MPI_INTEGER, blocksPerProc, &
         1, MPI_INTEGER, iComm, error)
    nBlocksUsedMax = maxval(blocksPerProc)
    nBlocksUsedMin = minval(blocksPerProc)
    offsetPerProc(0) = 0
    nBlkUsedGlobal = blocksPerProc(0) 
    if (nProc > 1) then
       do procIdx = 1, nProc-1
          offsetPerProc(procIdx) = offsetPerProc(procIdx-1) &
               + blocksPerProc(procIdx-1)
          nBlkUsedGlobal = nBlkUsedGlobal + blocksPerProc(procIdx)
       end do
    end if
    !The UsedBlocksTemp array will not serve our purpose because it will
    !have trailing negative numbers.  Move it to an array of the correct
    !size and use it to create an array of used Nodes.
    allocate(UsedBlocks(nBlocksUsed))
    UsedBlocks(1:nBlocksUsed) = UsedBlocksTemp(1:nBlocksUsed)
    deallocate(UsedBlocksTemp)


    if (nBlocksUsed > 0) then
       allocate(UsedNodes(nBlocksUsed))
       UsedNodes(1:nBlocksUsed) = iNode_B(UsedBlocks(1:nBlocksUsed))
    else 
       allocate(UsedNodes(1))
       UsedNodes(1) = 0
    end if

    allocate(unknownNameArray(nPlotVar))
    do iVar = 1, nPlotVar
       UnknownNameArray(iVar) = plotVarNames(iVar)
       !The VisIt plugin needs null padded names.
       labelLeng = len_trim(UnknownNameArray(iVar))
       do iLen = labelLeng + 1,lNameH5 
          UnknownNameArray(iVar)(iLen:iLen) = CHAR(0)
       end do
    end do
    call write_plot_string(nPlotVar,lNameH5, UnknownNameArray,"plotVarNames",&
         fileID)

    ! The plugn doesn't need these. However if does use this info to speed some 
    ! things up if it is there.  Perhaps this could be an option in PARAM.in
    allocate (blockDataExtents(2,nBlocksUsed))
    do iVar = 1, nPlotVar
       do iBlk = 1, nBlocksUsed
          blockDataExtents(1, iBlk) = minval(PlotVarIdx(:,:,:,iBlk,iVar))
          blockDataExtents(2, iBlk) = maxval(PlotVarIdx(:,:,:,iBlk,iVar))
       end do
       dsname = trim(plotVarNames(iVar))//"_Ext"
       labelLeng = len_trim(UnknownNameArray(iVar)) + 4
       do iLen = labelLeng + 1,lNameH5+3
          dsname(iLen:iLen) = CHAR(0)
       end do
       !        if (nBlocksUsed == 0 ) then
       !           if (collectiveWrite)&
       !                call writeHdf5Rank2Real(fileID, reshape((/0.0,0.0,0.0,0.0/),(/2,2/)), nBlocksUsed,&
       !                nBlkUsedGlobal, offsetPerProc(iProc), dsname, 2, collectiveWrite)
       !        else
       call writeHdf5Rank2Real(fileID, blockDataExtents, nBlocksUsed,&
            nBlkUsedGlobal, offsetPerProc(iProc), dsname, 2, collectiveWrite)   
       !       end if
    end do

    deallocate (blockDataExtents)

    do iVar = 1, nPlotVar
       call writeHdf5Rank4Real(fileID, PlotVarIdx(:,:,:,1:nBlocksUsed,iVar), nBlocksUsed,&
            nBlkUsedGlobal, offsetPerProc(iProc), UnknownNameArray(iVar), &
            arrSize(1), arrSize(2), arrSize(3), .true., globalVarMin(iVar), globalVarMax(iVar),&
            collectiveWrite)
    end do
    allocate(bbox(2,nPlotDim,nBlocksUsed))    
    if (nonCartesian) then

       call writeHdf5Rank4Real(fileID, plotXYZNodes(1,:,:,:,1:nBlocksUsed), nBlocksUsed,&
            nBlkUsedGlobal, offsetPerProc(iProc), 'NodesX', &
            nodeArrSize(1), nodeArrSize(2), nodeArrSize(3), .false., 0.0, 0.0, collectiveWrite)
       if (nPlotDim .ge. 2) &
            call writeHdf5Rank4Real(fileID, plotXYZNodes(2,:,:,:,1:nBlocksUsed), nBlocksUsed,&
            nBlkUsedGlobal, offsetPerProc(iProc), 'NodesY', &
            nodeArrSize(1), nodeArrSize(2), nodeArrSize(3), .false., 0.0, 0.0,collectiveWrite)
       if (nPlotDim == 3) &
            call writeHdf5Rank4Real(fileID, plotXYZNodes(3,:,:,:,1:nBlocksUsed), nBlocksUsed,&
            nBlkUsedGlobal, offsetPerProc(iProc), 'NodesZ', &
            nodeArrSize(1), nodeArrSize(2), nodeArrSize(3), .false., 0.0, 0.0,collectiveWrite)


       do iBlk = 1, nBlocksUsed
          bbox(1,1,iBlk) = minval(plotXYZNodes(1,:,:,:,iBlk))
          bbox(2,1,iBlk) = maxval(plotXYZNodes(1,:,:,:,iBlk))
          if (nPlotDim .ge. 2) then
             bbox(1,2,iBlk) = minval(plotXYZNodes(2,:,:,:,iBlk))
             bbox(2,2,iBlk) = maxval(plotXYZNodes(2,:,:,:,iBlk))
          end if
          if (nPlotDim == 3) then
             bbox(1,3,iBlk) = minval(plotXYZNodes(3,:,:,:,iBlk))
             bbox(2,3,iBlk) = maxval(plotXYZNodes(3,:,:,:,iBlk))
          end if
       end do
       deallocate(plotXYZNodes)


    else
       bbox(1,:,:) = CoordMin_DB(plotDim(1:nPlotDim), UsedBlocks)
       bbox(2,:,:) = CoordMax_DB(plotDim(1:nPlotDim), UsedBlocks) 
    end if
    call writeHdf5Rank3Real(fileID, bbox, nBlocksUsed,&
         nBlkUsedGlobal, offsetPerProc(iProc), "bounding box", 2, nPlotDim,collectiveWrite)
    ! 
    !     call MPI_allreduce(minval(bbox(1,1,:)), minCoords(1), 1, MPI_REAL,&
    !          MPI_MIN, iComm, error)
    !     call MPI_allreduce(maxval(bbox(2,1,:)), maxCoords(1), 1, MPI_REAL,&
    !          MPI_MAX, iComm, error)
    !     if (nPlotDim .ge. 2) then
    !        call MPI_allreduce(minval(bbox(1,2,:)), minCoords(2), 1, MPI_REAL,&
    !             MPI_MIN, iComm, error)
    !        call MPI_allreduce(maxval(bbox(2,2,:)), maxCoords(2), 1, MPI_REAL,&
    !             MPI_MAX, iComm, error)
    !     end if
    !     if (nPlotDim == 3) then
    !        call MPI_allreduce(minval(bbox(1,3,:)), minCoords(3), 1, MPI_REAL,&
    !             MPI_MIN, iComm, error)
    !        call MPI_allreduce(maxval(bbox(2,3,:)), maxCoords(3), 1, MPI_REAL,&
    !             MPI_MAX, iComm, error)
    !     end if
    allocate(coordinates(nPlotDim, nBlocksUsed))
    do iBlk = 1, nBlocksUsed
       do i = 1, nPlotDim
          coordinates(i, iBlk) = .5*(bbox(1,i,iBlk) + bbox(2,i,iBlk))
       end do
    end do
    deallocate(bbox)

    call writeHdf5Rank2Real(fileID, coordinates, nBlocksUsed,&
         nBlkUsedGlobal, offsetPerProc(iProc), "coordinates", nPlotDim, collectiveWrite)
    deallocate(coordinates)

    if (nBlocksUsed > 0) then
       call writeHdf5Rank2Integer(fileID, iTree_IA(plotDim(1:nPlotDim)+Coord0_,UsedNodes), nBlocksUsed,&
            nBlkUsedGlobal, offsetPerProc(iProc), "MinLogicalExtents", nPlotDim, collectiveWrite)
    else if (collectiveWrite .and. nBlocksUsed == 0) then
       call writeHdf5Rank2Integer(fileID, reshape((/0,0,0,0/),(/2,2/)), nBlocksUsed,&
            nBlkUsedGlobal, offsetPerProc(iProc), "MinLogicalExtents", nPlotDim, collectiveWrite)
    end if
    !
    do iVar = 1, nPlotVar
       UnknownNameArray(iVar) = plotVarUnits(iVar)
       !The VisIt plugin needs null padded names.
       labelLeng = len_trim(UnknownNameArray(iVar))
       do iLen = labelLeng + 1,lNameH5 
          UnknownNameArray(iVar)(iLen:iLen) = CHAR(0)
       end do
    end do
    call write_plot_string(nPlotVar, lNameH5, UnknownNameArray, "plotVarUnits",&
         fileID)

    deallocate(PlotVarIdx)
    deallocate(unknownNameArray)

    allocate(bbox(2,nPlotDim, nBlocksUsed))

    if (.not. isCutFile .and. nBlocksUsed > 0) then
       call writeHdf5Rank1Integer(fileID, iMortonNode_A(UsedNodes), &
            nBlocksUsed, nBlkUsedGlobal, offsetPerProc(iProc),"iMortonNode_A",collectiveWrite)
    else if (.not. isCutFile .and. collectiveWrite) then
       call writeHdf5Rank1Integer(fileID, (/0/), &
            nBlocksUsed, nBlkUsedGlobal, offsetPerProc(iProc),"iMortonNode_A", collectiveWrite)
    end if

    if (nBlocksUsed > 0 ) then
       call writeHdf5Rank1Integer(fileID, iTree_IA(Level_,UsedNodes), &
            nBlocksUsed, nBlkUsedGlobal, offsetPerProc(iProc),"refine level", collectiveWrite)
    else if (collectiveWrite) then
       call writeHdf5Rank1Integer(fileID, (/0/), &
            nBlocksUsed, nBlkUsedGlobal, offsetPerProc(iProc),"refine level", collectiveWrite)
    end if

    deallocate(usedNodes)
    allocate(procNum(nBlocksUsed))
    procNum = iProc
    call writeHdf5Rank1Integer(fileID, procNum,nBlocksUsed, &
         nBlkUsedGlobal, offsetPerProc(iProc),"Processor Number", collectiveWrite)
    deallocate(procNum)

    deallocate(usedBlocks)
    allocate(UnknownNameArray(nPlotDim))

    if (plotType=='x=0') then
       UnknownNameArray(1) = "Y-Axis"
       UnknownNameArray(2) = "Z-Axis"
    else
       do iVar = 1, nPlotDim

          if (plotDim(iVar) == 1) then
             UnknownNameArray(iVar) = "X-Axis"
          else if (plotDim(iVar) == 2) then
             UnknownNameArray(iVar) = "Y-Axis"
          else if (plotDim(iVar) == 3) then
             UnknownNameArray(iVar) = "Z-Axis"
          end if
          !The VisIt plugin needs null padded names.
       end do
    end if
    do iVar = 1, nPlotDim
       labelLeng = len_trim(UnknownNameArray(iVar))
       do iLen = labelLeng + 1,lNameH5 
          UnknownNameArray(iVar)(iLen:iLen) = CHAR(0)
       end do
    end do

    call write_plot_string(nPlotDim, lNameH5, UnknownNameArray, "Axis Labels",&
         fileID)
    deallocate(UnknownNameArray)
    call write_integer_sim_metadata(fileID, nPlotVar)
    call write_real_sim_metadata(FileID,plot_dimensional)
    call write_integer_plot_metadata(fileID, nPlotVar, isCutFile)
    if (plotType=='x=0') then
       call write_real_plot_metadata(fileid,plot_dimensional, .true.)
    else
       call write_real_plot_metadata(FileID,plot_dimensional, .false.)
    end if

    !     call write_real_plot_metadata(FileID,minCoords(1),maxCoords(1),minCoords(2),maxCoords(2),&
    !          minCoords(3),maxCoords(3))
    !     call write_real_plot_metadata(FileID,xMin,xMax,yMin,yMax,&
    !          zMin,zMax)


    !   !if one of the hdf routines leaves a dataspace, dataset, proterty list, etc... open
    !   !the code will hang on h5garbage_collect_f
    !    if (iProc == 0) write (*,*) '  closing file'
    !Close the file
    call h5garbage_collect_f(error)
    call h5fclose_f(fileID,error)
    !closing the hdf5 interface
    call h5close_f(error)
    if (error == -1) write (*,*) 'h5fclose_f failed'

    !     if (iProc == 0) write (*,*) '===================================='
  end subroutine write_plot_hdf5

  !======================================================================
  !==========================================================================

  subroutine write_real_plot_metadata(FileID,plot_dimensional, isXZero)

    use ModMain, only : Time_Simulation, CodeVersion
    use ModProcMH, only : iProc
    use BATL_lib, only : nLevel, nDim
    use ModGeometry, only : x1,x2,y1,y2,z1,z2
    use ModPhysics, ONLY : No2Io_V, UnitX_


    integer (HID_T), intent(in) :: fileID
    integer(HID_T) :: dataset, dataspace
    integer(HSIZE_T) :: dimens1D(1)
    integer :: error, rank, iData,i
    logical, intent (in) :: plot_dimensional, isXZero

    real :: RealMetaData(7)
    !    allocate(attName(nAtts))
    iData = 1
    !    attName(1) = 'Simulation Time'
    RealMetaData(iData) = Time_Simulation
    iData = iData + 1
    if (isXZero) then
       if (plot_dimensional) then
          RealMetaData(iData) = y1*No2Io_V(UnitX_)
          iData = iData + 1
          RealMetaData(iData) = y2*No2Io_V(UnitX_)
          iData = iData + 1
          RealMetaData(iData) = z1*No2Io_V(UnitX_)
          iData = iData + 1
          RealMetaData(iData) = z2*No2Io_V(UnitX_)
          iData = iData + 1
          RealMetaData(iData) = 0.0
          iData = iData + 1
          RealMetaData(iData) = 0.0
          iData = iData + 1
       else
          RealMetaData(iData) = y1!*No2Io_V(UnitX_)
          iData = iData + 1
          RealMetaData(iData) = y2!*No2Io_V(UnitX_)
          iData = iData + 1
          RealMetaData(iData) = z1!*No2Io_V(UnitX_)
          iData = iData + 1
          RealMetaData(iData) = z2!*No2Io_V(UnitX_)
          iData = iData + 1
          RealMetaData(iData) = 0.0
          iData = iData + 1
          RealMetaData(iData) = 0.0
          iData = iData + 1
       end if
    else

       if (plot_dimensional) then
          do i=1,3
             if (plotDim(i) == 1) then
                RealMetaData(iData) = x1*No2Io_V(UnitX_)
                iData = iData + 1
                RealMetaData(iData) = x2*No2Io_V(UnitX_)
                iData = iData + 1
             else if(plotDim(i) == 2) then
                RealMetaData(iData) = y1*No2Io_V(UnitX_)
                iData = iData + 1
                RealMetaData(iData) = y2*No2Io_V(UnitX_)
                iData = iData + 1
             else if(plotDim(i) == 3) then
                RealMetaData(iData) = z1*No2Io_V(UnitX_)
                iData = iData + 1
                RealMetaData(iData) = z2*No2Io_V(UnitX_)
                iData = iData + 1
             else if (plotDim(i) == 0) then
                RealMetaData(iData) = 0.0
                iData = iData + 1
                RealMetaData(iData) = 0.0
                iData = iData + 1
             end if
          end do
       else
          do i=1,3
             if (plotDim(i) == 1) then
                RealMetaData(iData) = x1
                iData = iData + 1
                RealMetaData(iData) = x2
                iData = iData + 1
             else if(plotDim(i) == 2) then
                RealMetaData(iData) = y1
                iData = iData + 1
                RealMetaData(iData) = y2
                iData = iData + 1
             else if(plotDim(i) == 3) then
                RealMetaData(iData) = z1
                iData = iData + 1
                RealMetaData(iData) = z2
                iData = iData + 1
             else if (plotDim(i) == 0) then
                RealMetaData(iData) = 0.0
                iData = iData + 1
                RealMetaData(iData) = 0.0
                iData = iData + 1
             end if
          end do
       end if
    end if
    iData = iData -1

    !-------------------------------------------------------------------
    !write the real Metadata
    rank = 1
    dimens1D = iData

    call h5screate_simple_f(rank, dimens1D, dataspace, error) 
    call h5dcreate_f(fileID, "Real Plot Metadata",&
         H5T_NATIVE_DOUBLE,dataspace, dataset, error)


    if (iProc == 0)&
         call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, &
         RealMetaData,dimens1D, error, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
    if (error == -1)&
         call stop_mpi(&
         "error in subroutine init_grid. Error marker 3")
    call write_hdf5_attribute("Code Version", dataset, H5T_NATIVE_DOUBLE,&
         realAtt=CodeVersion)

    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)
  end subroutine write_real_plot_metadata
  !===========================================================================

  subroutine write_real_sim_metadata(FileID,plot_dimensional)

    use ModMain, only : Time_Simulation, CodeVersion
    use BATL_lib, only : nLevel, nDim
    use ModGeometry, only : x1,x2,y1,y2,z1,z2
    use ModPhysics, ONLY : No2Io_V, UnitX_

    integer (HID_T), intent(in) :: fileID
    logical, intent (in) :: plot_dimensional
    integer(HID_T) :: dataset, dataspace
    integer(HSIZE_T) :: dimens1D(1)
    integer :: error, rank, iData


    real :: RealMetaData(7)
    !    allocate(attName(nAtts))
    iData = 1
    !    attName(1) = 'Simulation Time'
    RealMetaData(iData) = Time_Simulation
    iData = iData + 1   
    !    attName(2) = 'xmin'

    if (plot_dimensional) then
       RealMetaData(iData) = x1*No2Io_V(UnitX_)
       iData = iData + 1
       !    attName(3) = 'xmax'
       RealMetaData(iData) = x2*No2Io_V(UnitX_)
       iData = iData + 1
       !    attName(4) = 'ymin'
       RealMetaData(iData) = y1*No2Io_V(UnitX_)
       iData = iData + 1
       !    attName(3) = 'ymax'
       RealMetaData(iData) = y2*No2Io_V(UnitX_)
       iData = iData + 1
       !    attName(6) = 'zmin'
       RealMetaData(iData) = z1*No2Io_V(UnitX_)
       iData = iData + 1
       !    attName(3) = 'zmax'
       RealMetaData(iData) = z2*No2Io_V(UnitX_)
    else
       RealMetaData(iData) = x1
       iData = iData + 1
       !    attName(3) = 'xmax'
       RealMetaData(iData) = x2
       iData = iData + 1
       !    attName(4) = 'ymin'
       RealMetaData(iData) = y1
       iData = iData + 1
       !    attName(3) = 'ymax'
       RealMetaData(iData) = y2
       iData = iData + 1
       !    attName(6) = 'zmin'
       RealMetaData(iData) = z1
       iData = iData + 1
       !    attName(3) = 'zmax'
       RealMetaData(iData) = z2
    end if
    !-------------------------------------------------------------------
    !write the real Metadata
    rank = 1
    dimens1D = iData

    call h5screate_simple_f(rank, dimens1D, dataspace, error) 
    call h5dcreate_f(fileID, "Real Simulation Metadata",&
         H5T_NATIVE_DOUBLE,dataspace, dataset, error)


    !   if (iProc == 0)&
    call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, &
         RealMetaData,dimens1D, error, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
    if (error == -1)&
         call stop_mpi(&
         "error in subroutine init_grid. Error marker 3")
    call write_hdf5_attribute("Code Version", dataset, H5T_NATIVE_DOUBLE,&
         realAtt=CodeVersion)

    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)
  end subroutine write_real_sim_metadata

  !==========================================================================

  subroutine write_integer_plot_metadata(fileID,nPlotVar,isCutFile)

    use ModProcMH, only : nProc
    use BATL_lib, only : nDimAmr, nLevel, IsPeriodic_D 
    use ModMain, only : n_step
    use ModGeometry, ONLY : TypeGeometry
    integer (HID_T), intent(in) :: fileID, nPlotVar
    logical, intent(in) :: isCutFile
    integer(HID_T) :: dataset, dataspace
    integer(HSIZE_T) :: dimens1D(1)

    integer :: error, rank, i, iData
    integer :: IntegerMetaData(16)
    integer, parameter  :: FFV = 1 ! file format version
    iData = 1

    IntegerMetaData(iData) = FFV
    !    attName(1) = 'File Format Version'
    iData = iData + 1
    IntegerMetaData(iData) = n_step 
    !    attName(2) = "Time Step"
    iData = iData + 1
    IntegerMetaData(iData) = nPlotDim
    !    attName(3) = 'nDim'
    iData = iData + 1
    IntegerMetaData(iData) = nDimAmr
    !    attName(3) = 'nDim'
    iData = iData + 1
    IntegerMetaData(iData) = nBlkUsedGlobal
    !    attName(4) = 'globalNumBlocks'
    iData = iData + 1    
    IntegerMetaData(iData) = nProc
    !    attName(5) = 'numProcessors'
    iData = iData + 1

    IntegerMetaData(iData) = nLevel
    !    attName(7) = 'nLevel'
    iData = iData + 1
    do i = 1, 3
       if (plotDim(i) == 1) then
          IntegerMetaData(iData) = iSizeGlobal
          iData = iData + 1
       else if (plotDim(i) == 2) then
          IntegerMetaData(iData) = jSizeGlobal
          iData = iData + 1
       else if (plotDim(i) == 3) then
          IntegerMetaData(iData) = kSizeGlobal
          iData = iData + 1
       else if (plotDim(i) == 0) then
          IntegerMetaData(iData) = 0
          iData = iData + 1

       end if
    end do

    !     IntegerMetaData(iData) = iSizeGlobal
    !     !    attName(8) = 'nxb' 
    !     iData = iData + 1
    !     IntegerMetaData(iData) = jSizeGlobal
    !     !    attName(9) = 'nyb'
    !     iData = iData + 1
    !     IntegerMetaData(iData) = kSizeGlobal
    !     !    attName(10) = 'nzb'
    ! 
    !HDF5 has no boolean datatype so any logical parameters should go here.

    if(TypeGeometry=='cartesian') then
       IntegerMetaData(iData) = 0
    elseif(TypeGeometry=='rz') then
       IntegerMetaData(iData) = 1
    elseif(TypeGeometry=='roundcube') then
       IntegerMetaData(iData) = 2
    elseif(TypeGeometry=='cylindrical') then
       IntegerMetaData(iData) = 3
    elseif(TypeGeometry=='cylindrical_lnr') then
       IntegerMetaData(iData) = 4
    elseif(TypeGeometry=='cylindrical_genr') then
       IntegerMetaData(iData) = 5
    elseif(TypeGeometry=='spherical') then
       IntegerMetaData(iData) = 6
    elseif(TypeGeometry=='spherical_lnr') then
       IntegerMetaData(iData) = 7
    elseif(TypeGeometry=='spherical_genr') then
       IntegerMetaData(iData) = 8
    elseif(TypeGeometry=='rlonlat') then
       IntegerMetaData(iData) = 9
    elseif(TypeGeometry=='rlonlat_lnr') then
       IntegerMetaData(iData) = 10
    elseif(TypeGeometry=='rlonlat_genr') then
       IntegerMetaData(iData) = 11
    endif
    iData = iData + 1
    !as of 2/3/2012 this is not implimented in the plugin but it probably
    !should be in the future
    do i = 1, 3
       if(IsPeriodic_D(i)) then
          IntegerMetaData(iData) = 1
       else
          IntegerMetaData(iData) = 0
       endif
       iData = iData + 1
    end do

    if (isCutFile) then
       integerMetaData(iData) = 1
    else
       integerMetaData(iData) = 0
    end if
    iData = iData + 1
    integerMetaData(iData) = nPlotVar

    !-------------------------------------------------------------------
    !write the integer Metadata
    rank = 1
    dimens1D = iData

    call h5screate_simple_f(rank, dimens1D, dataspace, error) 
    call h5dcreate_f(fileID, "Integer Plot Metadata",&
         H5T_NATIVE_INTEGER,dataspace, dataset, error)

    !    if (iProc == 0)&   
    call h5dwrite_f(dataset,H5T_NATIVE_INTEGER,IntegerMetaData,dimens1D,&
         error, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
    if (error == -1)&
         call stop_mpi(&
         "error in subroutine init_grid. Error marker 3")

    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)
  end subroutine write_integer_plot_metadata

  !=========================================================================

  subroutine write_integer_sim_metadata(fileID, nPlotVar)
    !Not read by plugin at this time  Only exists so that one looking at
    !the file may know something about the simulation that created it

    use ModProcMH, only : nProc
    use BATL_lib, only : nDim, nDimAmr, nLevel
    use ModMain, only : n_step, nI, nJ, nK
    integer (HID_T), intent(in) :: fileID, nPlotVar
    integer(HID_T) :: dataset, dataspace
    integer(HSIZE_T) :: dimens1D(1)

    integer :: error, rank, iData
    integer :: IntegerMetaData(8)
    iData = 1

    IntegerMetaData(iData) = nI
    !    attName(8) = 'nxb' 
    iData = iData + 1
    IntegerMetaData(iData) = nJ
    !    attName(9) = 'nyb'
    iData = iData + 1
    IntegerMetaData(iData) = nK
    !    attName(10) = 'nzb'
    iData = iData + 1
    IntegerMetaData(iData) = nDim
    !    attName(3) = 'nDim'
    iData = iData + 1
    IntegerMetaData(iData) = nDimAmr
    !    attName(3) = 'nDim'
    iData = iData + 1
    IntegerMetaData(iData) = nProc
    !    attName(5) = 'numProcessors'
    iData = iData + 1

    IntegerMetaData(iData) = nLevel
    !    attName(7) = 'nLevel'
    iData = iData + 1

    IntegerMetaData(iData) = n_step 
    !    attName(2) = "Time Step"
    !    iData = iData + 1

    !HDF5 has no boolean datatype so any logical parameters should go here.
    !  write character attributes for key.
    !     iData = iData + 1
    !     if(TypeGeometry=='cartesian') then
    !        IntegerMetaData(iData) = 0
    !     elseif(TypeGeometry=='rz') then
    !        IntegerMetaData(iData) = 1
    !     elseif(TypeGeometry=='roundcube') then
    !        IntegerMetaData(iData) = 2
    !     elseif(TypeGeometry=='cylindrical') then
    !        IntegerMetaData(iData) = 3
    !     elseif(TypeGeometry=='cylindrical_lnr') then
    !        IntegerMetaData(iData) = 4
    !     elseif(TypeGeometry=='cylindrical_genr') then
    !        IntegerMetaData(iData) = 5
    !     elseif(TypeGeometry=='spherical') then
    !        IntegerMetaData(iData) = 6
    !     elseif(TypeGeometry=='spherical_lnr') then
    !        IntegerMetaData(iData) = 7
    !     elseif(TypeGeometry=='spherical_genr') then
    !        IntegerMetaData(iData) = 8
    !     elseif(TypeGeometry=='rlonlat') then
    !        IntegerMetaData(iData) = 9
    !     elseif(TypeGeometry=='rlonlat_lnr') then
    !        IntegerMetaData(iData) = 10
    !     elseif(TypeGeometry=='rlonlat_genr') then
    !        IntegerMetaData(iData) = 11
    !     endif
    ! 
    !     !as of 2/3/2012 this is not implimented in the plugin but it probably
    !     !should be in the future
    !     do i = 1, 3
    !        if(IsPeriodic_D(i)) then
    !           IntegerMetaData(iData+i) = 1
    !        else
    !           IntegerMetaData(iData+i) = 0
    !        endif
    !     end do
    !     iData = iData + 3



    !-------------------------------------------------------------------
    !write the integer Metadata
    rank = 1
    dimens1D = iData

    call h5screate_simple_f(rank, dimens1D, dataspace, error) 
    call h5dcreate_f(fileID, "Integer Simulation Metadata",&
         H5T_NATIVE_INTEGER,dataspace, dataset, error)

    !    if (iProc == 0)&   
    call h5dwrite_f(dataset,H5T_NATIVE_INTEGER,IntegerMetaData,dimens1D,&
         error, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
    if (error == -1)&
         call stop_mpi(&
         "error in subroutine init_grid. Error marker 3")

    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)
  end subroutine write_integer_sim_metadata

  !=====================================================================
end module ModHdf5
