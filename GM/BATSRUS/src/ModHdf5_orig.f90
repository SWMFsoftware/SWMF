module ModHdf5

  use ModProcMH, ONLY: iProc, nProc, iComm 
  use hdf5

  implicit none
  SAVE


  private ! except
  public write_var_hdf5
  public write_cut_var_hdf5
  public write_plot_hdf5

  real,  allocatable :: PlotVarIdx(:,:,:,:,:)
  real,    allocatable :: plotXYZNodes(:,:,:,:,:)

  integer, allocatable :: UsedBlocks(:)
  integer, allocatable :: usednodes(:)
  integer, allocatable :: UsedBlocksTemp(:)

  integer, parameter:: lNameVar = 10
  integer, parameter:: lNameH5 = lNameVar + 1

  character (len=lNameH5), allocatable :: UnknownNameArray(:)

  integer :: nBlocksUsed, nBlkUsedGlobal, iBlk ,arrSize(3), nodeArrSize(3)
  integer :: iSizeGlobal, jSizeGlobal, kSizeGlobal, nPlotDim, plotDim(3) 
  real, allocatable :: Position_DA(:,:)
  logical :: collectiveWrite
contains

  subroutine write_cut_var_hdf5(iFile, iBlock, H5Index,nPlotVar,PlotVar, &
       xMin, xMax, yMin, yMax, zMin, zMax, DxBlock, DyBlock, DzBlock,&
       isNonCartesian,nCell,H5Advance)

    ! Save all cells within plotting range, for each processor

    use ModProcMH
    use ModMain, ONLY: nI, nJ, nK, PROCtest, BLKtest, test_string, &
         x_, y_, z_, Phi_, nBlockMax, unusedBLK
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK, Dx_BLK, Dy_BLK, Dz_BLK,&
         x1, x2, y1, y2, z1, z2, XyzStart_BLK, XyzMin_D, XyzMax_D
    use ModCovariant, ONLY: is_axial_geometry          
    use ModPhysics, ONLY : No2Io_V, UnitX_
    use ModIO
    use ModNumConst
    use ModMpi
    use ModNodes, ONLY: NodeX_NB,NodeY_NB,NodeZ_NB
    use BATL_lib, ONLY : nDim, Xyz_DNB

    implicit none

    ! Arguments

    integer, intent(in)   :: iFile, iBlock, H5Index
    integer, intent(in)   :: nPlotVar
    real,    intent(in)   :: PlotVar(-1:nI+2,-1:nJ+2,-1:nK+2,nPlotVar)
    real,    intent(in)   :: xMin,xMax,yMin,yMax,zMin,zMax
    logical, intent(in) :: isNonCartesian
    real,    intent(inout):: DxBlock,DyBlock,DzBlock
    integer, intent(out)  :: nCell 
    logical, intent(out) :: H5advance

    ! Local variables
    ! Indices and coordinates
    integer :: iVar, i, j, k, i2, j2, k2, iMin, iMax, jMin, jMax, kMin, kMax
    integer ::  iMaxN,jMaxN, kMaxN
    integer :: ii, jj, kk, iSize, jSize, kSize, d
    integer :: nRestrict, nRestrictX, nRestrictY, nRestrictZ
    integer :: sizePlotBlkI, sizePlotBlkJ, sizePlotBlkK, error
    real :: x,y,z,Dx,Restrict
    real :: xMin1,xMax1,yMin1,yMax1,zMin1,zMax1
    real :: DxBlockOut
    real :: Plot_V(nPlotVarMax)

    real :: ySqueezed

    real, parameter:: cHalfMinusTiny=cHalf*(cOne-cTiny)

    character(len=*), parameter :: NameSub = 'write_plot_hdf'
    logical :: DoTest, DoTestMe


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

       if (iSizeGlobal == 1 .and. jSizeGlobal == 1) then
          nPlotDim = 1
          plotDim = (/3,0,0/)
          arrSize = (/kSizeGlobal,1,1/)
       else if (iSizeGlobal == 1 .and. kSizeGlobal == 1) then
          nPlotDim = 1
          plotDim = (/2,0,0/)
          arrSize = (/jSizeGlobal,1,1/)
       else if (jSizeGlobal == 1 .and. kSizeGlobal == 1) then
          nPlotDim = 1
          plotDim = (/1,0,0/)
          arrSize = (/iSizeGlobal,1,1/)
       else if (iSizeGlobal == 1) then
          nPlotDim = 2
          plotDim = (/2,3,0/)
          arrSize = (/jSizeGlobal,kSizeGlobal,1/)
       else if (jSizeGlobal == 1) then
          nPlotDim = 2
          plotDim = (/1,3,0/)
          arrSize = (/iSizeGlobal,kSizeGlobal,1/)
       else if (kSizeGlobal == 1) then
          nPlotDim = 2
          plotDim = (/1,2,0/)
          arrSize = (/iSizeGlobal,jSizeGlobal,1/)
       else
          nPlotDim = 3
          plotDim = (/1,2,3/)
          arrSize = (/iSizeGlobal,jSizeGlobal,kSizeGlobal/)
       end if
       allocate(PlotVarIdx(arrSize(1), arrSize(2), arrSize(3),&
            nBlockMax, nPlotVar))
       !       allocate(PlotVarIdx(nI, nJ, nK, nBlockMax, nPlotVar))
       PlotVarIdx=-50
       nodeArrSize(1:3) = arrSize(1:3)
       do ii = 1, nPlotDim
          nodeArrSize(ii) = arrSize(ii) + 1
       end do
       if (isNonCartesian) allocate(plotXYZNodes(nPlotDim, &
            nodeArrSize(1),nodeArrSize(2),nodeArrSize(3),nBlockMax))
       !             allocate(plotXYZNodes(nodeArrSize(1),nodeArrSize(2),&
       !                 nodeArrSize(3),nPlotDim,nBlockMax))


       !       allocate(cutIndexLimits(6,nBlockMax))

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

    kMin = max(1 ,floor((zMin1-xyzStart_BLK(z_,iBlock))/DzBlock)+2)
    kMax = min(nK,floor((zMax1-xyzStart_BLK(z_,iBlock))/DzBlock)+1)

    iMaxN = iMax
    jMaxN = jMax
    kMaxN = kMax

    do ii = 1, nPlotDim
       if (plotDim(ii) == 1) then
          iMaxN = iMax + 1
       else if (plotDim(ii) == 2) then
          jMaxN = jMax + 1
       else if (plotDim(ii) == 3) then
          kMaxN = kMax + 1
       end if
    end do


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
    do kk=kMin,kMax; do jj=jMin,jMax; do ii=iMin,iMax

       if (nPlotDim == 2) then
          if (kSizeGlobal == 1) then
             i = ii - iMin + 1
             j = jj - jMin + 1
             k = kk - kMin + 1
          else if (jSizeGlobal == 1) then
             i = ii - iMin + 1
             j = kk - kMin + 1
             k = jj - jMin + 1
          else if (iSizeGlobal == 1) then
             i = jj - jMin + 1
             j = kk - kMin + 1
             k = ii - iMin + 1
          end if
       else if (nPlotDim == 1) then
          if (plotDim(1) == 1) then
             i = ii - iMin + 1
             j = jj - jMin + 1
             k = kk - kMin + 1
          else if (plotDim(1) == 2) then
             i = jj - jMin + 1
             j = ii - iMin + 1
             k = kk - kMin + 1
          else if (plotDim(1) == 3) then
             i = kk - kMin + 1
             j = ii - iMin + 1
             k = kk - kMin + 1
          end if
       else
          i = ii - iMin + 1
          j = jj - jMin + 1
          k = kk - kMin + 1
       end if
       PlotVarIdx(i, j, k, H5Index, 1:nPlotVar) = &
            PlotVar(ii, jj, kk, 1:nPlotVar)
       nCell = nCell+1
    end do; end do; end do


    ! this seems to do y=0 and z=0 plots okay but x=0 plots appear to be broken
    !I have only tried one test case so the problem be bigger than just x=0 cuts
    if (isNonCartesian) then
       do kk=kMin,kMaxN; do jj=jMin,jMaxN; do ii=iMin,iMaxN

          if (nPlotDim == 2) then
             if (kSizeGlobal == 1) then
                i = ii - iMin + 1
                j = jj - jMin + 1
                k = kk - kMin + 1
             else if (jSizeGlobal == 1) then
                i = ii - iMin + 1
                j = kk - kMin + 1
                k = jj - jMin + 1
             else if (iSizeGlobal == 1) then
                i = jj - jMin + 1
                j = kk - kMin + 1
                k = ii - iMin + 1
             end if
          else if (nPlotDim == 1) then
             if (plotDim(1) == 1) then
                i = ii - iMin + 1
                j = jj - jMin + 1
                k = kk - kMin + 1
             else if (plotDim(1) == 2) then
                i = jj - jMin + 1
                j = ii - iMin + 1
                k = kk - kMin + 1
             else if (plotDim(1) == 3) then
                i = kk - kMin + 1
                j = ii - iMin + 1
                k = kk - kMin + 1
             end if
          else
             i = ii - iMin + 1
             j = jj - jMin + 1
             k = kk - kMin + 1
          end if
          do d=1,nPlotVar
             plotXYZNodes(d,i, j, k, H5Index) = &
                  Xyz_DNB(plotDim(d),ii, jj, kk, iBlock)
          end do

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
    use BATL_lib, ONLY : nDim, Xyz_DNB, Xyz_DGB
    use ModNodes, ONLY: NodeX_NB,NodeY_NB,NodeZ_NB
    use ModCoordTransform, only : xyz_to_sph, sph_to_xyz
    use ModSort, only : sort_quick_func


    !  Handles adding one block of variable data that for plotting
    !  in the next plotfile.

    integer, intent(in) :: nPlotVar, H5Index, iBlock
    real, intent(in) :: PlotVar(-1:nI+2,-1:nJ+2,-1:nK+2,nPlotVar)
    logical, intent(in) :: isNonCartesian
    integer:: i,j,k, ii, jj, kk, iPlot, error, ijkIndex(3,nI*nJ*nK), &
         iCell, cellIndex(nI*nJ*nK)

    !----------------------------------------------------------------------
    ! Allocate the storage array if it has not been allocated for
    ! this plotfile yet

!!! Use nBlock if possible

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

    do iPlot = 1, nPlotVar
       do kk = 1, nK; do jj = 1, nJ; do ii = 1, nI
          PlotVarIdx(ii, jj, kk, H5Index, iPlot) = &
               PlotVar(ii, jj, kk, iPlot)


       end do; end do; end do
    end do
    UsedBlocksTemp(H5Index) = iBlock
    nBlocksUsed = H5Index


  end subroutine write_var_hdf5
  !===============================================================           

  logical function sort_test(i, j)
    use BATL_lib, ONLY : nDim 

    integer, intent(in) :: i,j
    if(ndim == 3) then
       sort_test = 0
       if(Position_DA(3,i) < Position_DA(3,j)) then
          sort_test = .false.
       elseif(Position_DA(3,j) < Position_DA(3,i)) then
          sort_test = .true.
       else
          if(Position_DA(2,i) < Position_DA(2,j)) then
             sort_test = .false.
          elseif(Position_DA(2,j) < Position_DA(2,i)) then
             sort_test = .true.
          else
             if(Position_DA(1,i) < Position_DA(1,j)) then
                sort_test = .false.
             elseif(Position_DA(1,j) < Position_DA(1,i)) then
                sort_test = .true.
             end if
          end if
       end if
    elseif(nDim == 2) then
       if(Position_DA(2,i) < Position_DA(2,j)) then
          sort_test = .false.
       elseif(Position_DA(2,j) < Position_DA(2,i)) then
          sort_test = .true.
       else
          if(Position_DA(1,i) < Position_DA(1,j)) then
             sort_test = .false.
          elseif(Position_DA(1,j) < Position_DA(1,i)) then
             sort_test = .true.
          end if
       end if
    elseif(nDim == 1) then
       if(Position_DA(1,i) < Position_DA(1,j)) then
          sort_test = .false.
       elseif(Position_DA(1,j) < Position_DA(1,i)) then
          sort_test = .true.
       end if

    else
       sort_test = .false.
    end if
  end function sort_test
  !=================================================================
  !================================================================

  subroutine write_plot_hdf5(filename, plotVarNames, plotVarUnits,&
       nPlotVar,isCutFile, nonCartesian,plot_dimensional, xmin, xmax, &
       ymin, ymax, zmin, zmax)!, nBLKcells)

    use BATL_tree, only: iNode_B, iTree_IA, Coord1_, Coord2_, Coord3_,&
         Level_,iMortonNode_A, MaxNode, get_tree_position, nNodeUsed,&
         Block_, Proc_,iNodeMorton_I, Unused_BP
    use BATL_lib, only : CoordMin_DB, CoordMax_DB
    use ModMpi
    use ModMain, only : nI,nJ,nK, nBlockMax
    use ModGeometry, ONLY : x_BLK, y_BLK, z_BLK
    use BATL_lib, ONLY : nDim, MaxDim
    use ModSort, only : sort_quick_func
    use BATL_lib, ONLY : CoordMin_DB

    integer,                 intent(in):: nPlotVar
    character(len=80),       intent(in):: filename
    character(len=lNameVar), intent(in):: plotVarNames(nPlotVar)
    character(len=lNameVar),  intent(in):: plotVarUnits(nPlotVar)
    character(len=lNameVar) :: AxisLabels(3)!, AxisLabelsNull(3)
    real, intent(in)  ::  xMin, xMax, yMin, yMax, zMin, zMax
    !    integer,                 intent(in):: nBLKcells
    logical, intent(in) :: isCutFile,plot_dimensional, nonCartesian
    integer :: offsetPerProc(0:nProc-1), blocksPerProc(0:nProc-1)
    integer :: fileID, error, procIdx, lastCoord, iNode, iLen, iVar 
    integer :: labelLeng, i, ii, PosMinType, iMorton
    integer :: hdfMajor, hdfMinor, hdfRelease, minBlocksUsed
    real, allocatable :: coordinates(:,:), bbox(:,:,:),blockDataExtents(:,:)
    integer, allocatable :: minLogicalExtents(:,:)
    integer, allocatable :: procNum(:)
    integer, allocatable :: VisItIndx(:)
    integer, allocatable :: NumUnsedBlkTo(:,:)
    integer(HSIZE_T) :: openObjs
    real :: PositionMin_D(MaxDim), PositionMax_D(MaxDim), varMin, varMax
    real :: globalVarMin(nPlotVar), globalVarMax(nPlotVar)
    real :: minCoords(3), maxCoords(3)
    character (len=lnamevar+4) :: dsname
    character (len=4) :: tempstr


    !--------------------------------------------------------------------------

    if (iProc == 0) write (*,*) '===================================='
    if (iProc == 0) write (*,*) 'Writing BATL hdf5 parallel plotfile'
    if (iProc == 0) write (*,*) '------------------------------------'
    if (iProc == 0) write (*,*) '  opening hdf5 file'

    do iVar = 1, nPlotVar
       if (nBlocksUsed > 0) then
          varMin = minval(PlotVarIdx(:,:,:,:,iVar))
          varMax = maxval(PlotVarIdx(:,:,:,:,iVar))
       else
          varMin = huge(varMin)
          varMax = -huge(varMax)
       end if
       call MPI_allreduce(varMin, globalVarMin(iVar), 1, MPI_REAL,&
            MPI_MIN, iComm, error)

       call MPI_allreduce(varMax, globalVarMax(iVar), 1, MPI_REAL,&
            MPI_MAX, iComm, error)

       call MPI_Allgather(nBlocksUsed, 1, MPI_INTEGER, blocksPerProc, &
            1, MPI_INTEGER, iComm, error)
    end do


    !initialize the Fortran Hdf5 library
    call h5open_f(error)

    !determine if we want to do collective writes. Collective write
    !is faster but all processors must call h5dwrite, even in cuts where 
    !some processors write no data.  In newer versions of hdf5 a null
    !write can be called but attempting to do so on older versions will
    !cause the code to crash.  

    call h5get_libversion_f(hdfMajor, hdfMinor, hdfRelease, error)   
    if (hdfMajor > 1) then
       collectiveWrite = .true.
    else if (hdfMajor == 1) then
       if (hdfMinor .ge. 8) then
          collectiveWrite = .true.
       else if (hdfMinor == 6 .or. hdfMinor == 7) then
          call MPI_allreduce(nBlocksUsed, minBlocksUsed, 1, MPI_INTEGER, &
               MPI_MIN, iComm, error)
          if (minBlocksUsed == 0) then
             collectiveWrite = .false.
          else 
             collectiveWrite = .true.
          end if
       else
          collectiveWrite = .false.
       end if
    else
       collectiveWrite = .false.
    end if



    ! Open the hdf5 file and write data to it
    !---------------------------------------------------------------------    
    fileID = -1
    call hdf5_init_file(fileID, filename)
    if(fileID == -1) then
       if (iProc == 0) write (*,*)  "Error: unable to initialize file"
       call stop_mpi("unable to initialize hdf5 file")
    end if
    !get processor block offsets for writes
    call MPI_Allgather(nBlocksUsed, 1, MPI_INTEGER, blocksPerProc, &
         1, MPI_INTEGER, iComm, error)

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

    allocate(UsedNodes(nBlocksUsed))
    UsedNodes(1:nBlocksUsed) = iNode_B(UsedBlocks(1:nBlocksUsed))
    allocate(NumUnsedBlkTo(nBlockMax,0:nProc-1))
    NumUnsedBlkTo = 0
    do procIdx = 0, nProc-1
       if(Unused_BP(1, procIdx))&
            NumUnsedBlkTo(1, procIdx) = 1
       do iBlk = 2, nBlockMax
          if (Unused_BP(iBlk,procIdx)) then
             NumUnsedBlkTo(iBlk,procIdx)=NumUnsedBlkTo(iBlk-1,ProcIdx)&
                  +1
          else
             NumUnsedBlkTo(iBlk,procIdx)=NumUnsedBlkTo(iBlk-1,procIdx)
          end if
       end do
    end do
    allocate(Position_DA(nPlotDim,nNodeUsed))
    do iMorton = 1, nNodeUsed
       iNode = iNodeMorton_I(iMorton)
       iBlk = iTree_IA(Block_,iNode)
       procIdx = iTree_IA(Proc_,iNode) 
       i = offsetPerProc(procIdx) + iBlk - NumUnsedBlkTo(iBlk,procIdx)
       call get_tree_position(iNode, PositionMin_D, PositionMax_D)
       Position_DA(:,i) = PositionMin_D(1:nPlotDim)
    end do
    deallocate(NumUnsedBlkTo)

    allocate(VisItIndx(nBlkUsedGlobal))        
    call sort_quick_func(nBlkUsedGlobal, sort_test, VisItIndx)
    deallocate(Position_DA)
    call writeHdf5Rank1Integer(fileID,VisItIndx, nBlkUsedGlobal,&
         nBlkUsedGlobal, 0,"VisIt Index")

    deallocate(VisItIndx)

    allocate(unknownNameArray(nPlotVar))
    do iVar = 1, nPlotVar
       UnknownNameArray(iVar) = plotVarNames(iVar)
       !The VisIt plugin needs null padded names.
       labelLeng = len_trim(UnknownNameArray(iVar))
       do iLen = labelLeng + 1,lNameH5 
          UnknownNameArray(iVar)(iLen:iLen) = CHAR(0)
       end do
    end do
    call write_plot_string(nPlotVar,UnknownNameArray,"plotVarNames",&
         fileID)


    ! The plugn doesn't need these. However if does use this info to speed some 
    ! things up if it is there.  Perhaps this could be an option in PARAM.in
    allocate (blockDataExtents(2,nBlocksUsed))

    do iVar = 1, nPlotVar
       do iBlk = 1, nBlocksUsed
          blockDataExtents(1, iBlk) = minval(PlotVarIdx(:,:,:,iBlk,iVar))
          blockDataExtents(2, iBlk) = maxval(PlotVarIdx(:,:,:,iBlk,iVar))
       end do

       dsname = UnknownNameArray(iVar)
       do iLen = 1,lNameH5
          if (UnknownNameArray(iVar)(iLen:iLen) == CHAR(0)) then
             labelLeng = iLen
             exit
          end if
       end do

       tempstr = "_Ext"
       dsname(labelLeng:3+labelLeng) = tempstr(1:4)

       call writeHdf5Rank2Real(fileID, blockDataExtents, nBlocksUsed,&
            nBlkUsedGlobal, offsetPerProc(iProc), trim(dsname), 2)

    end do

    deallocate (blockDataExtents)

    do iVar = 1, nPlotVar
       call writeHdf5Rank4Real(fileID, PlotVarIdx(:,:,:,1:nBlocksUsed,iVar), nBlocksUsed,&
            nBlkUsedGlobal, offsetPerProc(iProc), UnknownNameArray(iVar), &
            arrSize(1), arrSize(2), arrSize(3), .true., globalVarMin(iVar), globalVarMax(iVar))
    end do

    allocate(bbox(2,nPlotDim,nBlocksUsed))    
    if (nonCartesian) then

       call writeHdf5Rank4Real(fileID, plotXYZNodes(1,:,:,:,1:nBlocksUsed), nBlocksUsed,&
            nBlkUsedGlobal, offsetPerProc(iProc), 'NodesX', &
            nodeArrSize(1), nodeArrSize(2), nodeArrSize(3), .false., 0.0, 0.0)
       if (nPlotDim .ge. 2) &
            call writeHdf5Rank4Real(fileID, plotXYZNodes(2,:,:,:,1:nBlocksUsed), nBlocksUsed,&
            nBlkUsedGlobal, offsetPerProc(iProc), 'NodesY', &
            nodeArrSize(1), nodeArrSize(2), nodeArrSize(3), .false., 0.0, 0.0)
       if (nPlotDim == 3) &
            call writeHdf5Rank4Real(fileID, plotXYZNodes(3,:,:,:,1:nBlocksUsed), nBlocksUsed,&
            nBlkUsedGlobal, offsetPerProc(iProc), 'NodesZ', &
            nodeArrSize(1), nodeArrSize(2), nodeArrSize(3), .false., 0.0, 0.0)


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
         nBlkUsedGlobal, offsetPerProc(iProc), "bounding box", 2, nPlotDim)

    call MPI_allreduce(minval(bbox(1,1,:)), minCoords(1), 1, MPI_REAL,&
         MPI_MIN, iComm, error)
    call MPI_allreduce(maxval(bbox(2,1,:)), maxCoords(1), 1, MPI_REAL,&
         MPI_MAX, iComm, error)
    if (nPlotDim .ge. 2) then
       call MPI_allreduce(minval(bbox(1,2,:)), minCoords(2), 1, MPI_REAL,&
            MPI_MIN, iComm, error)
       call MPI_allreduce(maxval(bbox(2,2,:)), maxCoords(2), 1, MPI_REAL,&
            MPI_MAX, iComm, error)
    end if
    if (nPlotDim == 3) then
       call MPI_allreduce(minval(bbox(1,3,:)), minCoords(3), 1, MPI_REAL,&
            MPI_MIN, iComm, error)
       call MPI_allreduce(maxval(bbox(2,3,:)), maxCoords(3), 1, MPI_REAL,&
            MPI_MAX, iComm, error)
    end if
    allocate(coordinates(nPlotDim, nBlocksUsed))
    do iBlk = 1, nBlocksUsed
       do i = 1, nPlotDim
          coordinates(i, iBlk) = .5*(bbox(1,i,iBlk) + bbox(2,i,iBlk))
       end do
    end do
    deallocate(bbox)
    call writeHdf5Rank2Real(fileID, coordinates, nBlocksUsed,&
         nBlkUsedGlobal, offsetPerProc(iProc), "coordinates", nPlotDim)
    deallocate(coordinates)

    do iVar = 1, nPlotVar
       UnknownNameArray(iVar) = plotVarUnits(iVar)
       !The VisIt plugin needs null padded names.
       labelLeng = len_trim(UnknownNameArray(iVar))
       do iLen = labelLeng + 1,lNameH5 
          UnknownNameArray(iVar)(iLen:iLen) = CHAR(0)
       end do
    end do
    call write_plot_string(nPlotVar, UnknownNameArray, "plotVarUnits",&
         fileID)

    deallocate(PlotVarIdx)
    deallocate(unknownNameArray)

    allocate(bbox(2,nPlotDim, nBlocksUsed))
    allocate(minLogicalExtents(nPlotDim, nBlocksUsed))

    do iBlk = 1, nBlocksUsed
       do i = 1,nPlotDim
          ii = Coord1_+i-1
          minLogicalExtents(i,iBlk) = iTree_IA(ii,iNode)&
               -1.0
       end do
    end do

    call writeHdf5Rank2Integer(fileID, minLogicalExtents, nBlocksUsed,&
         nBlkUsedGlobal, offsetPerProc(iProc), "MinLogicalExtents", nPlotDim)

    deallocate(minLogicalExtents)

    if (.not. isCutFile)& 
         call writeHdf5Rank1Integer(fileID, iMortonNode_A(UsedNodes), &
         nBlocksUsed, nBlkUsedGlobal, offsetPerProc(iProc),"iMortonNode_A")

    call writeHdf5Rank1Integer(fileID, iTree_IA(Level_,UsedNodes), &
         nBlocksUsed, nBlkUsedGlobal, offsetPerProc(iProc),"refine level")
    deallocate(usedNodes)

    allocate(procNum(nBlocksUsed))
    procNum = iProc
    call writeHdf5Rank1Integer(fileID, procNum,nBlocksUsed, &
         nBlkUsedGlobal, offsetPerProc(iProc),"Processor Number")
    deallocate(procNum)

    deallocate(usedBlocks)

    allocate(UnknownNameArray(nPlotDim))
    do iVar = 1, nPlotDim
       if (plotDim(iVar) == 1) then
          UnknownNameArray(iVar) = "X-Axis"
       else if (plotDim(iVar) == 2) then
          UnknownNameArray(iVar) = "Y-Axis"
       else if (plotDim(iVar) == 3) then
          UnknownNameArray(iVar) = "Z-Axis"
       end if
       !The VisIt plugin needs null padded names.
       labelLeng = len_trim(UnknownNameArray(iVar))
       do iLen = labelLeng + 1,lNameH5 
          UnknownNameArray(iVar)(iLen:iLen) = CHAR(0)
       end do
    end do
    call write_plot_string(nPlotDim, UnknownNameArray, "Axis Labels",&
         fileID)
    deallocate(UnknownNameArray)

    call write_integer_sim_metadata(fileID, nPlotVar)
    call write_real_sim_metadata(FileID,plot_dimensional)
    call write_integer_plot_metadata(fileID, nPlotVar, isCutFile)
    call write_real_plot_metadata(FileID,minCoords(1),maxCoords(1),minCoords(2),maxCoords(2),&
         minCoords(3),maxCoords(3))
    if (iProc == 0) write (*,*) '  closing file'
    !Close the file



    call h5garbage_collect_f(error)
    call h5fclose_f(fileID,error)
    !closing the hdf5 interface
    call h5close_f(error)
    if (error == -1) write (*,*) 'h5fclose_f failed'

    if (iProc == 0) write (*,*) '===================================='

    ! Deallocate the variable array
  end subroutine write_plot_hdf5
  !======================================================================
  !=====================================================================

  subroutine writeHdf5Rank3Real(fileID, dataBuff, localNumBlocks,&
       globalNumBlocks, localOffset, description, nIplot, nJplot)

    implicit none

    integer :: rank
    integer, intent(in) :: nIplot,nJplot
    integer :: error
    integer(HID_T) :: dataset, plist_id
    integer(HSIZE_T) :: dimens3D(3), dimens1D(1)
    integer(HID_T) :: dataspace, attributeSpace, attribute
    integer(HID_T) :: memspace
    integer(HSIZE_T) :: start3D(3)
    integer(HSIZE_T) :: stride3D(3)
    integer(HSIZE_T) :: count3D(3)
    integer(HID_T), intent(in) :: fileID
    integer, intent(in) :: localNumBlocks, globalNumBlocks
    real, intent(in) :: dataBuff(:,:,:)
    integer, intent(in) :: localOffset
    character (len=*), intent(in) :: description
    integer, pointer :: nullPointer => NULL( )
    integer(HSIZE_T) :: nullSize(1) = 1
    integer(HSIZE_T) :: openObjs, one

    !Set the dimensions of the dataset

    rank = 3
    dimens3D(1) = nIplot
    dimens3D(2) = nJplot
    dimens3D(3) = globalNumBlocks
    call h5open_f(error)
    if (error == -1) &
         call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 1")


    call h5screate_simple_f(rank, dimens3D, dataspace, error) 
    !create the dataset
    call h5dcreate_f(fileID, description, H5T_NATIVE_DOUBLE, dataspace, dataset, error)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    if (collectiveWrite) then
       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    else
       !          call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
    end if

    if (error == -1) &
         call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 2")

    if (localNumBlocks == 0) then
       one = 0
       call h5screate_simple_f(1, (/one/), memspace, error)
       if (collectiveWrite) then
          call h5sselect_none_f(memspace,error)
          call h5sselect_none_f(dataspace,error)
          nullSize(1:rank) = 0
          call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, nullPointer, nullSize, error,&
               mem_space_id = memspace,file_space_id = dataspace, &
               xfer_prp = plist_id)
       end if
       call h5sclose_f(memspace,error)
       call h5pclose_f(plist_id, error)
       call h5sclose_f(dataspace, error)
       call h5dclose_f(dataset, error)
       return
    end if


    start3D(1:2) = 0
    start3D(3) = localOffset 
    stride3D(:) = 1
    count3D(1:2) = dimens3D(1:2)
    count3D(3) = localNumBlocks


    if (localNumBlocks > 0 ) then    
       !create the hyperslab.  This will differ on the different processors
       call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start3D, count3D, error, &
            stride = stride3D)
       !         !create the memory space
       if (error == -1)& 
            call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 3")

       dimens3D(3) = localNumBlocks
       call h5screate_simple_f(rank, dimens3D, memspace, error)
       !Write the data
       if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 3")

       call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, dataBuff, dimens3D, error, &
            mem_space_id = memspace, file_space_id = dataspace, xfer_prp = plist_id)
    end if
    if (error == -1) &
         call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 5")

    call h5sclose_f(memspace,error)
    call h5pclose_f(plist_id, error)
    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)

    if (error == -1)& 
         call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 6")

  end subroutine writeHdf5Rank3Real



  !======================================================================
  !=====================================================================

  subroutine writeHdf5Rank1Integer(fileID, dataBuff, localNumBlocks,&
       globalNumBlocks, localOffset, description)

    implicit none

    integer :: rank, dimens
    integer :: error
    integer(HID_T) :: dataset, plist_id
    integer(HSIZE_T) :: dimens1D(2)
    integer(HID_T) :: dataspace
    integer(HID_T) :: memspace
    integer(HSIZE_T) :: start1D(2)
    integer(HSIZE_T) :: stride1D(2)
    integer(HSIZE_T) :: count1D(2)
    integer(HID_T), intent(in) :: fileID
    integer, intent(in) :: localNumBlocks, globalNumBlocks
    integer, intent(in) :: dataBuff(:)
    integer, intent(in) :: localOffset
    character (len=*), intent(in) :: description
    integer, pointer :: nullPointer => NULL( )
    integer(HSIZE_T) :: nullSize(1) = 1
    !Set the dimensions of the dataset
    rank = 1
    dimens1D(1) = globalNumBlocks
    call h5open_f(error)
    if (error == -1) &
         
         call stop_mpi("error in subroutine writeHdf5Rank1Real. Error marker 1")


    call h5screate_simple_f(rank, dimens1D, dataspace, error) 
    !create the dataset
    call h5dcreate_f(fileID, description, H5T_NATIVE_INTEGER, dataspace, dataset, error)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    if (collectiveWrite) then
       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    else
       !nothing = INDEPENDENT  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
    end if

    if (localNumBlocks == 0) then
       if (collectiveWrite) then
          call h5screate_f(H5S_NULL_F, memspace, error)
          call h5sselect_none_f(dataspace,error)
          call h5dwrite_f(dataset, H5T_NATIVE_INTEGER, nullPointer, nullSize, error,&
               file_space_id = dataspace, xfer_prp = plist_id)
          call h5sclose_f(memspace,error)
       end if
       call h5pclose_f(plist_id, error)
       call h5sclose_f(dataspace, error)
       call h5dclose_f(dataset, error)
       return
    end if

    if (error == -1) &
         call stop_mpi("error in subroutine writeHdf5Rank1Real. Error marker 2")

    start1D(1) = localOffset 
    stride1D(:) = 1
    count1D(1) = localNumBlocks

    if (localNumBlocks > 0 ) then    
       !create the hyperslab.  This will differ on the different processors
       call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start1D, count1D, error, &
            stride = stride1D)
       !         !create the memory space
       if (error == -1)& 
            call stop_mpi("error in subroutine writeHdf5Rank1Real. Error marker 3")

       dimens1D(1) = localNumBlocks
       call h5screate_simple_f(rank, dimens1D, memspace, error)
       !Write the data
       if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Rank1Real. Error marker 4")

       call h5dwrite_f(dataset, H5T_NATIVE_INTEGER, dataBuff, dimens1D, error, &
            mem_space_id = memspace, file_space_id = dataspace, xfer_prp = plist_id)


    end if
    if (error == -1) &
         call stop_mpi("error in subroutine writeHdf5Rank1Real. Error marker 5")

    call h5sclose_f(memspace,error)
    call h5pclose_f(plist_id, error)
    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)

    if (error == -1)& 
         call stop_mpi("error in subroutine writeHdf5Rank1Real. Error marker 6")

  end subroutine writeHdf5Rank1Integer

  !=====================================================================
  !=====================================================================

  subroutine writeHdf5Rank2Integer(fileID, dataBuff, localNumBlocks,&
       globalNumBlocks, localOffset, description, dimens)

    implicit none

    integer :: rank, dimens
    integer :: error
    integer(HID_T) :: dataset, plist_id
    integer(HSIZE_T) :: dimens2D(2)
    integer(HID_T) :: dataspace
    integer(HID_T) :: memspace
    integer(HSIZE_T) :: start2D(2)
    integer(HSIZE_T) :: stride2D(2)
    integer(HSIZE_T) :: count2D(2)
    integer(HID_T), intent(in) :: fileID
    integer, intent(in) :: localNumBlocks, globalNumBlocks
    integer, intent(in) :: dataBuff(:,:)
    integer, intent(in) :: localOffset
    character (len=*), intent(in) :: description
    integer, pointer :: nullPointer => NULL( )
    integer(HSIZE_T) :: one, nullSize(1) = 1

    !Set the dimensions of the dataset
    rank = 2
    dimens2D(1) = dimens
    dimens2D(2) = globalNumBlocks
    call h5open_f(error)
    if (error == -1) &
         
         call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 1")


    call h5screate_simple_f(rank, dimens2D, dataspace, error) 
    !create the dataset
    call h5dcreate_f(fileID, description, H5T_NATIVE_INTEGER, dataspace, dataset, error)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    if (collectiveWrite) then
       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    else
       !         call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
    end if

    if (localNumBlocks == 0) then
       if (collectiveWrite) then
          call h5screate_f(H5S_NULL_F, memspace, error)
          call h5sselect_none_f(dataspace,error)
          call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, nullPointer, nullSize, error,&
               file_space_id = dataspace, xfer_prp = plist_id)
          call h5sclose_f(memspace,error)
       end if
       call h5pclose_f(plist_id, error)
       call h5sclose_f(dataspace, error)
       call h5dclose_f(dataset, error)
       return
    end if

    if (error == -1) &
         call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 2")

    start2d(1) = 0
    start2D(2) = localOffset 
    stride2D(:) = 1
    count2D(1) = dimens2D(1)
    count2D(2) = localNumBlocks

    if (localNumBlocks > 0 ) then    
       !create the hyperslab.  This will differ on the different processors
       call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start2D, count2D, error, &
            stride = stride2D)
       !         !create the memory space
       if (error == -1)& 
            call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 3")

       dimens2D(2) = localNumBlocks
       call h5screate_simple_f(rank, dimens2D, memspace, error)
       !Write the data
       if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 4")

       call h5dwrite_f(dataset, H5T_NATIVE_INTEGER, dataBuff, dimens2D, error, &
            mem_space_id = memspace, file_space_id = dataspace, xfer_prp = plist_id)
    end if
    if (error == -1) &
         call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 5")

    call h5sclose_f(memspace,error)
    call h5pclose_f(plist_id, error)
    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)

    if (error == -1)& 
         call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 6")

  end subroutine writeHdf5Rank2Integer

  !======================================================================
  !=====================================================================

  subroutine writeHdf5Rank2Real(fileID, dataBuff, localNumBlocks,&
       globalNumBlocks, localOffset, description, dimens)

    implicit none

    integer :: rank, dimens
    integer :: error
    integer(HID_T) :: dataset, plist_id
    integer(HSIZE_T) :: dimens2D(2)
    integer(HID_T) :: dataspace
    integer(HID_T) :: memspace
    integer(HSIZE_T) :: start2D(2)
    integer(HSIZE_T) :: stride2D(2)
    integer(HSIZE_T) :: count2D(2)
    integer(HID_T), intent(in) :: fileID
    integer, intent(in) :: localNumBlocks, globalNumBlocks
    real, intent(in) :: dataBuff(:,:)
    integer, intent(in) :: localOffset
    character (len=*), intent(in) :: description
    integer, pointer :: nullPointer => NULL( )
    integer(HSIZE_T) :: one, nullSize(1) = 1

    !Set the dimensions of the dataset
    rank = 2
    dimens2D(1) = dimens
    dimens2D(2) = globalNumBlocks
    call h5open_f(error)
    if (error == -1) &
         
         call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 1")


    call h5screate_simple_f(rank, dimens2D, dataspace, error) 
    !create the dataset
    call h5dcreate_f(fileID, description, H5T_NATIVE_DOUBLE, dataspace, dataset, error)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    if (collectiveWrite) then
       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    else
       !         call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
    end if


    if (localNumBlocks == 0) then
       one = 0
       call h5screate_simple_f(1, (/one/), memspace, error)
       if (collectiveWrite) then
          call h5sselect_none_f(memspace,error)
          call h5sselect_none_f(dataspace,error)
          nullSize(1:rank) = 0
          call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, nullPointer, nullSize, error,&
               mem_space_id = memspace,file_space_id = dataspace, &
               xfer_prp = plist_id)
       end if
       call h5sclose_f(memspace,error)
       call h5pclose_f(plist_id, error)
       call h5sclose_f(dataspace, error)
       call h5dclose_f(dataset, error)
       return
    end if
    if (error == -1) &
         call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 2")

    start2d(1) = 0
    start2D(2) = localOffset 
    stride2D(:) = 1
    count2D(1) = dimens2D(1)
    count2D(2) = localNumBlocks

    if (localNumBlocks > 0 ) then    
       !create the hyperslab.  This will differ on the different processors
       call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start2D, count2D, error, &
            stride = stride2D)
       !         !create the memory space
       if (error == -1)& 
            call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 3")

       dimens2D(2) = localNumBlocks
       call h5screate_simple_f(rank, dimens2D, memspace, error)
       !Write the data
       if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 4")

       call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, dataBuff, dimens2D, error, &
            mem_space_id = memspace, file_space_id = dataspace, xfer_prp = plist_id)
    end if
    if (error == -1) &
         call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 5")

    call h5sclose_f(memspace,error)
    call h5pclose_f(plist_id, error)
    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)

    if (error == -1)& 
         call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 6")

  end subroutine writeHdf5Rank2Real

  !======================================================================
  !=====================================================================
  subroutine hdf5_init_file(fileID, filename)

    use ModMpi

    character (len=80), intent(in) :: filename
    integer :: error, accTemplate
    integer, intent(inout) :: fileID
    INTEGER :: mpierror       ! MPI error flag
    INTEGER :: dxpList
!!!    INTEGER :: mpi_size, mpi_rank, comm

    !Initalize currentBlock andplotArrayBegin for this file

    !    call h5open_f(error)                    

    !create MPI info Object

    !Create file access propertty list
    call h5pcreate_f(H5P_FILE_ACCESS_F, accTemplate, error)
    CALL h5pset_fapl_mpio_f(accTemplate, iComm, MPI_INFO_NULL, error)
    ! Create the file collectively.

    CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, fileID, error,&
         access_prp = accTemplate)
    CALL h5pclose_f(accTemplate, error)
    if (error == -1) &
         call stop_mpi(&
         "error in subroutine hdf5_init_file. Error marker 1")

    if (error == -1) fileID = -1
  end subroutine hdf5_init_file

  !=====================================================================
  !=====================================================================

  subroutine write_plot_string(nPlotVar,dataBuff, description, fileID)

    integer, intent(in) :: nPlotVar
    character (len=11), intent(in) :: dataBuff(nPlotVar)
    integer (HID_T), intent(in) :: fileID
    character (len=*), intent(in) :: description
    integer (HID_T) :: dataType, dataset, dataspace
    integer (HSIZE_T) :: sizeString, dimens1D(1) 

    integer :: error, rank, iLen, iVar, labelleng

    rank = 1  
    dimens1d(1) = nPlotVar

    sizeString = lNameH5


    call h5tcopy_f(H5T_NATIVE_CHARACTER, datatype, error)
    call h5tset_size_f(datatype, sizeString, error)
    if (error == -1) &
         call stop_mpi(&
         "error in subroutine writeHdf5Header. Error marker 10")
    call h5screate_simple_f(rank, dimens1D, dataspace, error)
    call h5dcreate_f(&
         fileID, description, datatype, dataspace, dataset, error)
    if (error == -1) &
         call stop_mpi("error in subroutine writeHdf5Header. Error marker 11")

    if (iProc == 0) &
         call h5dwrite_f(dataset, datatype, UnknownNameArray, dimens1D, error, H5S_ALL_F,&
         H5S_ALL_F, H5P_DEFAULT_F)
    if (error == -1) &
         call stop_mpi("error in subroutine writeHdf5Header. Error marker 11")



    call h5tclose_f(datatype, error)
    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)

    if (error == -1)& 
         call stop_mpi(&
         "error in subroutine writeHdf5Header. Error marker 13")
  end subroutine write_plot_string

  !=====================================================================
  !=====================================================================

  subroutine writeHdf5Rank4Real(fileID, dataBuff, localNumBlocks,&
       globalNumBlocks, localOffset, description, nIplot, nJplot, nKplot,&
       minMax ,vMin, vMax)

    implicit none

    integer :: rank
    integer, intent(in) :: nIplot,nJplot,nKplot
    integer :: error
    integer(HID_T) :: dataset, plist_id
    integer(HSIZE_T) :: dimens4D(4), dimens1D(1)
    integer(HID_T) :: dataspace, attributeSpace, attribute
    integer(HID_T) :: memspace
    integer(HSIZE_T) :: start4D(4)
    integer(HSIZE_T) :: stride4D(4)
    integer(HSIZE_T) :: count4D(4)
    integer(HID_T), intent(in) :: fileID
    integer, intent(in) :: localNumBlocks, globalNumBlocks
    real, intent(in) :: dataBuff(:,:,:,:)
    integer, intent(in) :: localOffset
    character (len=*), intent(in) :: description
    logical, intent(in) :: minMax
    real, intent(in) :: vMin, vMax
    integer, pointer :: nullPointer => NULL( )
    integer(HSIZE_T) :: nullSize(1) = 1
    integer(HSIZE_T) :: openObjs, one

    !Set the dimensions of the dataset

    rank = 4
    dimens4D(1) = nIplot
    dimens4D(2) = nJplot
    dimens4D(3) = nKplot
    dimens4D(4) = globalNumBlocks
    call h5open_f(error)
    if (error == -1) &
         call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 1")


    call h5screate_simple_f(rank, dimens4D, dataspace, error) 
    !create the dataset
    call h5dcreate_f(fileID, description, H5T_NATIVE_DOUBLE, dataspace, dataset, error)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    if (collectiveWrite) then
       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    else
       !          call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
    end if

    if (error == -1) &
         call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 2")

    if (minMax) then
       !Minimum
       dimens1D=(1)
       call h5screate_simple_f(1, dimens1D, attributeSpace, error)

       call h5acreate_f(dataset, "minimum", H5T_NATIVE_DOUBLE, attributeSpace, attribute,&
            error, H5P_DEFAULT_F)
       call h5awrite_f(attribute, H5T_NATIVE_DOUBLE, vMin, dimens1D, error)
       call h5sclose_f(attributeSpace, error)
       call h5aclose_f(attribute, error)
       !Maximum
       dimens1D=(1)
       call h5screate_simple_f(1, dimens1D, attributeSpace, error)

       call h5acreate_f(dataset, "maximum", H5T_NATIVE_DOUBLE, attributeSpace, attribute,&
            error, H5P_DEFAULT_F)
       call h5awrite_f(attribute, H5T_NATIVE_DOUBLE, vMax, dimens1D, error)
       call h5sclose_f(attributeSpace, error)
       call h5aclose_f(attribute, error)

    end if

    if (localNumBlocks == 0) then
       one = 0
       call h5screate_simple_f(1, (/one/), memspace, error)
       if (collectiveWrite) then
          call h5sselect_none_f(memspace,error)
          call h5sselect_none_f(dataspace,error)
          nullSize(1:rank) = 0
          call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, nullPointer, nullSize, error,&
               mem_space_id = memspace,file_space_id = dataspace, &
               xfer_prp = plist_id)
       end if
       call h5sclose_f(memspace,error)
       call h5pclose_f(plist_id, error)
       call h5sclose_f(dataspace, error)
       call h5dclose_f(dataset, error)
       return
    end if


    start4D(1:3) = 0
    start4D(4) = localOffset 
    stride4D(:) = 1
    count4D(1:3) = dimens4D(1:3)
    count4D(4) = localNumBlocks


    if (localNumBlocks > 0 ) then    
       !create the hyperslab.  This will differ on the different processors
       call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start4D, count4D, error, &
            stride = stride4D)
       !         !create the memory space
       if (error == -1)& 
            call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 3")

       dimens4D(4) = localNumBlocks
       call h5screate_simple_f(rank, dimens4D, memspace, error)
       !Write the data
       if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 4")

       call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, dataBuff, dimens4D, error, &
            mem_space_id = memspace, file_space_id = dataspace, xfer_prp = plist_id)
    end if
    if (error == -1) &
         call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 5")

    call h5sclose_f(memspace,error)
    call h5pclose_f(plist_id, error)
    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)

    if (error == -1)& 
         call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 6")

  end subroutine writeHdf5Rank4Real



  !======================================================================
  !=====================================================================

  subroutine write_real_plot_metadata(FileID,xmin,xmax,ymin,ymax, zmin,zmax)

    use ModMain, only : Time_Simulation, CodeVersion
    use ModProcMH, only : iProc
    use BATL_lib, only : nLevel, nDim
    use ModPhysics, ONLY : No2Io_V, UnitX_
    use ModIO, only : plot_range
    use ModGeometry, only : XyzMin_D, XyzMax_D, x1,x2,y1,y2,z1,z2

    integer (HID_T), intent(in) :: fileID
    integer(HID_T) :: dataset, dataspace, attribute, attributeSpace
    integer :: levelFactorArray(0:nLevel, nBlocksUsed)
    integer :: minLogicalExtents(ndim,nBlocksUsed)
    integer(HSIZE_T) :: dimens1D(1)
    integer :: error, rank, iData, i

    real, intent(in) :: xmin,xmax,ymin,ymax, zmin,zmax

    real :: RealMetaData(7)
    !    allocate(attName(nAtts))
    iData = 1
    !    attName(1) = 'Simulation Time'
    RealMetaData(iData) = Time_Simulation
    iData = iData + 1
    RealMetaData(iData) = xMin
    iData = iData + 1
    RealMetaData(iData) = xMax
    iData = iData + 1
    RealMetaData(iData) = yMin
    iData = iData + 1
    RealMetaData(iData) = yMax
    iData = iData + 1
    RealMetaData(iData) = zMin
    iData = iData + 1
    RealMetaData(iData) = zMax

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

    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)
  end subroutine write_real_plot_metadata
  !===========================================================================

  subroutine write_real_sim_metadata(FileID,plot_dimensional)

    use ModMain, only : Time_Simulation, CodeVersion
    use ModProcMH, only : iProc
    use BATL_lib, only : nLevel, nDim
    use ModGeometry, only : XyzMin_D, XyzMax_D, x1,x2,y1,y2,z1,z2
    use ModPhysics, ONLY : No2Io_V, UnitX_

    integer (HID_T), intent(in) :: fileID
    logical, intent (in) :: plot_dimensional
    integer(HID_T) :: dataset, dataspace, attribute, attributeSpace
    integer :: levelFactorArray(0:nLevel, nBlocksUsed)
    integer :: minLogicalExtents(ndim,nBlocksUsed)
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
    dimens1D=(1)
    call h5screate_simple_f(1, dimens1D, attributeSpace, error)

    call h5acreate_f(dataset, "Code Version", H5T_NATIVE_DOUBLE,&
         attributeSpace, attribute,&
         error, H5P_DEFAULT_F)
    call h5awrite_f(attribute, H5T_NATIVE_DOUBLE, CodeVersion, dimens1D, error)
    call h5sclose_f(attributeSpace, error)
    call h5aclose_f(attribute, error)


    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)
  end subroutine write_real_sim_metadata

  !==========================================================================

  subroutine write_integer_plot_metadata(fileID,nPlotVar,isCutFile)

    use ModProcMH, only : nProc, iProc
    use BATL_lib, only : nDim, nDimAmr, nLevel, IsPeriodic_D 
    use ModMain, only : n_step
    use ModGeometry, ONLY : TypeGeometry
    integer (HID_T), intent(in) :: fileID, nPlotVar
    logical, intent(in) :: isCutFile
    integer(HID_T) :: dataset, dataspace
    integer(HSIZE_T) :: dimens1D(1)

    integer :: error, rank, i, iData
    integer :: IntegerMetaData(16)
    integer :: FFV = 1 ! file format version
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

    use ModProcMH, only : nProc, iProc
    use BATL_lib, only : nDim, nDimAmr, nLevel, IsPeriodic_D 
    use ModMain, only : n_step, nI, nJ, nK
    use ModGeometry, ONLY : TypeGeometry
    integer (HID_T), intent(in) :: fileID, nPlotVar
    integer(HID_T) :: dataset, dataspace
    integer(HSIZE_T) :: dimens1D(1)

    integer :: error, rank, i, iData
    integer :: IntegerMetaData(8)
    integer :: FFV = 1 ! file format version
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
