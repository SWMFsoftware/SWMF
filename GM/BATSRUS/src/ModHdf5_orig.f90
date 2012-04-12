module ModHdf5

  use ModProcMH, ONLY: iProc, nProc, iComm 
  use hdf5

  implicit none
  SAVE


  private ! except
  public write_var_hdf5
  public write_plot_hdf5
  public write_cut_var_hdf5
  real,  allocatable :: PlotVarIdx(:,:,:,:,:)
  real, allocatable :: plotXYZNodes(:,:,:,:,:)
  real, allocatable :: plotXYZNodes_DNB(:,:,:,:,:)
  integer, allocatable :: UsedBlocks(:)
  integer, allocatable :: usednodes(:)
  integer, allocatable :: UsedBlocksTemp(:)
  integer, allocatable :: cutIndexLimits(:,:)

  integer, parameter:: lNameVar = 10
  integer, parameter:: lNameH5 = lNameVar + 1

  character (len=lNameH5), allocatable :: UnknownNameArray(:)

  integer :: nBlocksUsed, nBlkUsedGlobal, iBlk ,arrSize(3), nodeArrSize(3)
  integer :: iSizeGlobal, jSizeGlobal, kSizeGlobal, nPlotDim, plotDim(3) 
  integer(HSIZE_T) :: nBlocksUsedHSIZE_T, nBlkUsedGlobalHSIZE_T

  real, allocatable :: Position_DA(:,:)

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

  character(len=*), parameter :: NameSub = 'write_plot_idl'
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
            cycle
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
        if (isNonCartesian) &
            allocate(plotXYZNodes(nPlotDim,nodeArrSize(1),nodeArrSize(2),nodeArrSize(3),nBlockMax))
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
  !===============================================================
  !=============================================================== 
 
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
    integer:: i,j,k, ii, jj, kk, iPlot, error, ijkIndex(3,nI*nJ*nK), iCell, cellIndex(nI*nJ*nK)

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
  !=================================================================
  !================================================================
!   subroutine fill_nodeXYZ
!     ! Fill array with position (optionally dimensioned)
!     if (plot_dimensional) then
!        NodeXYZ_N(:,:,:,1:3)=PlotXYZNodes_NBI(:,:,:,iBLK,1:3)*No2Io_V(UnitX_)
!     else
!        NodeXYZ_N(:,:,:,1:3)=PlotXYZNodes_NBI(:,:,:,iBLK,1:3)
!     end if
!   end subroutine fill_nodeXYZ
  !=================================================================
  !================================================================

  subroutine write_plot_hdf5(filename, plotVarNames, plotVarUnits,&
    nPlotVar, isCutFile, isnonCartesian,plot_dimensional)

    !REORDER TO REDUCE MEM CONSUMPTION AND CONSIDER MAKING MORE THINGS
    !ALLOCATABLE

    use BATL_tree, only: iNode_B, iTree_IA, Coord0_, Coord1_, Coord2_,&
        Coord3_,Level_,iMortonNode_A, MaxNode, get_tree_position,&
        nNodeUsed,Block_, Proc_,iNodeMorton_I, Unused_BP
    use ModMpi
    use ModMain, only : nI,nJ,nK, nBlockMax, UseBatl, nBLK, nBlock, UnusedBlk
    use ModGeometry, ONLY : x_BLK, y_BLK, z_BLK, TypeGeometry
    use BATL_lib, ONLY : nDim, MaxDim, CoordMin_DB, CoordMax_DB, message_pass_node, calc_error_amr_criteria, Xyz_DNB,&
       message_pass_node

    use ModSort, only : sort_quick_func
    use ModProcMH, only : nProc, iProc, iComm
    use ModNodes, ONLY: NodeX_NB,NodeY_NB,NodeZ_NB
    use ModPhysics, ONLY : No2Io_V, UnitX_
    use BATL_geometry, only : coord_to_xyz
    use ModCoordTransform, only : xyz_to_sph, sph_to_xyz

    integer,                 intent(in):: nPlotVar
    character(len=80),       intent(in):: filename
    character(len=lNameVar), intent(in):: plotVarNames(nPlotVar)
    character(len=lNameVar),  intent(in):: plotVarUnits(nPlotVar)
    real :: xMin, xMax, yMin, yMax, zMin, zMax

    integer :: offsetPerProc(0:nProc-1), blocksPerProc(0:nProc-1)
    integer :: fileID, error, procIdx, lastCoord, iNode, iLen, iVar 
    integer :: labelLeng, PosMinType, iMorton, ii, jj, kk,i, j,k,b
    integer :: sizeplotblkI, sizeplotblkK,sizeplotblkJ,maxBlocksPerProc
    integer :: minPlotI, maxPlotI, minPlotJ, maxPlotJ, minPlotK, maxPlotK
    real, allocatable :: coordinates(:,:)
    integer, allocatable :: minLogicalExtents(:,:),UsedNodesProc(:,:)
    integer, allocatable :: procNum(:)
    integer, allocatable :: VisItIndx(:), VisItIdxLocal(:)
    integer, allocatable :: NumUnusedBlkTo(:,:)
    real :: PositionMin_D(MaxDim), PositionMax_D(MaxDim),PositionMin(MaxDim)
    real, allocatable :: globalVarMin(:), globalVarMax(:)
    real, allocatable :: blockDataExtents(:,:)
    real :: varMin, varMax, coordinate(nDim), tempReal
    logical, intent(in) :: isCutFile, isnonCartesian,plot_dimensional
    integer(HID_T) ::accTemplate, dataset, dataspace, memspace, plist_id
    integer(HSIZE_T) :: rank, localSize(5),globalSize(5),offset(5),TempHSTInt
    character (len=lnamevar+4) :: dsname
    character (len=4) :: tempstr
    
    logical :: hasData
    
    
    !------------------------------------------------------------------
    if (iProc == 0) write (*,*) '===================================='
    if (iProc == 0) write (*,*) 'Writing BATL hdf5 parallel plotfile'
    if (iProc == 0) write (*,*) '------------------------------------'
    if (iProc == 0) write (*,*) '  opening hdf5 file'


    !initialize the Fortran Hdf5 library
    call h5open_f(error)
    ! Open the hdf5 file and write data to it
    offset = 0
    !---------------------------------------------------------------------    

    fileID = -1
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

    if(fileID == -1) then
       if (iProc == 0) write (*,*)  "Error: unable to initialize file"
       call stop_mpi("unable to initialize hdf5 file")
    end if

    !collect min and max variable values before writing because otherwise
    !MPI_allreduce will hang on cuts with uneven bloc-proc distribution
    allocate(globalVarMin(nPlotVar)); allocate(globalVarMax(nPlotVar))
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
    
    if (nBlocksUsed > 0) then
        hasData = .true.
    else
        hasData = .false.
    end if

    !get processor block offsets for writes
!    write (*,*) "write_plot_hdf5 error marker 2"
    offsetPerProc(0) = 0
    nBlkUsedGlobal = blocksPerProc(0) 
    if (nProc > 1) then
       do procIdx = 1, nProc-1
          offsetPerProc(procIdx) = offsetPerProc(procIdx-1) &
               + blocksPerProc(procIdx-1)
          nBlkUsedGlobal = nBlkUsedGlobal + blocksPerProc(procIdx)
       end do
    end if
    nBlocksUsedHSIZE_T = nBlocksUsed
    nBlkUsedGlobalHSIZE_T = nBlkUsedGlobal
    !The UsedBlocksTemp array will not serve our purpose because it will
    !have trailing negative numbers.  Move it to an array of the correct
    !size and use it to create an array of used Nodes.
!    if (.not. isCutFile) then
        allocate(UsedBlocks(nBlocksUsed))
        UsedBlocks(1:nBlocksUsed) = UsedBlocksTemp(1:nBlocksUsed)
        deallocate(UsedBlocksTemp)
!    end if


    maxBlocksPerProc = maxval(blocksPerProc)
    allocate(UsedNodesProc(maxBlocksPerProc,0:nProc-1))
  
   if (nBlocksUsed > 0) then
        allocate(UsedNodes(nBlocksUsed))
        UsedNodes(1:nBlocksUsed) = iNode_B(usedBlocks(1:nBlocksUsed))
        call MPI_Allgather(UsedNodes, nBlocksUsed, MPI_INTEGER, &
            UsedNodesProc, maxBlocksPerProc, MPI_INTEGER, iComm, error)
   else
    !because if this block has to send something or MPI_Allgather fails
       allocate(UsedNodes(1))
       call MPI_Allgather(UsedNodes, 1, MPI_INTEGER, &
           UsedNodesProc, maxBlocksPerProc, MPI_INTEGER, iComm, error)

    end if
     allocate(Position_DA(nPlotDim,nNodeUsed))
        
        do procIdx = 0, nProc-1
                do iBlk = 1, blocksPerProc(procIdx) 
                    iNode = UsedNodesProc(iBlk,procIdx)
                    i = offsetPerProc(procIdx) + iBlk
                    if (blocksPerProc(procIdx) > 0) then               
                        call get_tree_position(iNode, PositionMin, PositionMax_D)
                        call coord_to_xyz(PositionMin,PositIonMax_D) 
                        call xyz_to_sph(PositionMax_D,PositionMin_D)
                        
                        do ii = 1, nPlotDim
                            Position_DA(ii,i) = PositionMin_D(plotDim(ii))
                        end do
                    end if
                end do
        end do
    
    allocate(VisItIndx(nBlkUsedGlobal))        

    call sort_quick_func(nBlkUsedGlobal, sort_test_block, VisItIndx)

    deallocate(Position_DA)

    allocate(VisItIdxLocal(nBlocksUsed))
    !We just want the results for local proc
    VisItIdxLocal(1:nBlocksUsed) = VisItIndx(offsetPerProc(iProc)+1:&
    nBlocksUsed+offsetPerProc(iProc))

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
    localSize(4) = 2
    localSize(5) = nBlocksUsed
    globalSize(4) = 2
    globalSize(5) = nBlkUsedGlobal
    offset(5) = offsetPerProc(iProc)
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


       if (hasData) then
            call writeHdf5DataSetHS(fileID,H5T_NATIVE_DOUBLE,2,&
            trim(dsname),localSize(4:5),globalSize(4:5),&
            offset(4:5), .true., .true.,hasData,dataspace,memspace,&
            dataset,plist_id,R2Real = blockDataExtents)
       else 
            call writeHdf5DataSetHS(fileID,H5T_NATIVE_DOUBLE,2,&
            trim(dsname),localSize(4:5),globalSize(4:5),&
            offset(4:5), .true., .true.,hasData,dataspace,&
            memspace,dataset,plist_id)
       end if

    end do

    deallocate (blockDataExtents)

    localSize(2:4) = arrSize(1:3)
    localSize(5) = nBlocksUsed
    globalSize(2:4) = arrSize(1:3)
    globalSize(5) = nBlkUsedGlobal

    do iVar = 1, nPlotVar
       if (nBlocksUsed > 0) then
            call writeHdf5DataSetHS(fileID,H5T_NATIVE_DOUBLE,4,&
            UnknownNameArray(iVar),localSize(2:5),globalSize(2:5),&
            offset(2:5), .false., .true.,hasData,dataspace,memspace,&
            dataset,plist_id,R4Real=PlotVarIdx(:,:,:,1:nBlocksUsed,iVar))
       else 
            call writeHdf5DataSetHS(fileID,H5T_NATIVE_DOUBLE,4,&
            UnknownNameArray(iVar),localSize(2:5),globalSize(2:5),&
            offset(2:5), .false., .true.,hasData,dataspace,&
            memspace,dataset,plist_id)
       end if


        call write_hdf5_attribute('minimum', dataset, H5T_NATIVE_DOUBLE,&
            realAtt=globalVarMin(iVar))
 
        call write_hdf5_attribute('maximum', dataset, H5T_NATIVE_DOUBLE,&
            realAtt=globalVarMax(iVar))
        call hdf5_close(memspace, plist_id, dataspace, dataset)

            
       
    end do
    deallocate(globalVarMin); deallocate(globalVarMax)
    deallocate(PlotVarIdx)

! 
  if (isnonCartesian) then
    
    localSize(5) = nBlocksUsed
    globalSize(5) = nBlkUsedGlobal
    localSize(2:4) = nodeArrSize(1:3)
    globalSize(2:4) = nodeArrSize(1:3)
    localSize(1) = nPlotDim
    globalSize(1) = nPlotDim
    offset(1:4) = 0
    offset(5) = offsetPerProc(iProc)

    call writeHdf5DataSetHS(fileID, H5T_NATIVE_DOUBLE,4,'NodesX',&
    localSize(2:5),globalSize(2:5),offset(2:5),.true.,.true.,&
    hasData, dataspace, memspace, dataset, plist_id,&
    R4Real = plotXYZNodes(1,:,:,:,1:nBlocksUsed))
    if (nPlotDim .GE. 2) &
     call writeHdf5DataSetHS(fileID, H5T_NATIVE_DOUBLE,4,'NodesY',&
    localSize(2:5),globalSize(2:5),offset(2:5),.true.,.true.,&
    hasData, dataspace, memspace, dataset, plist_id,&
    R4Real = plotXYZNodes(2,:,:,:,1:nBlocksUsed))
    if (nPlotDim == 3) &
    call writeHdf5DataSetHS(fileID, H5T_NATIVE_DOUBLE,4,'NodesZ',&
    localSize(2:5),globalSize(2:5),offset(2:5),.true.,.true.,&
    hasData, dataspace, memspace, dataset, plist_id,&
    R4Real = plotXYZNodes(3,:,:,:,1:nBlocksUsed))
    

   end if
  !else


    allocate(minLogicalExtents(nPlotDim, nBlocksUsed))

    lastCoord = Coord1_ + nPlotDim-1
    if (nBlocksUsed > 0) then
    do iBLK = 1, nBlocksUsed
       iNode = UsedNodes(iBLK)
       do ii = 1, nPlotDim
       minLogicalExtents(ii,iBlk) = iTree_IA(Coord0_+plotDim(ii),iNode)&
            -1.0
        end do
    end do
    else
        minLogicalExtents = 0
    end if
    if (isCutFile) then
        do b = 1, nPlotDim
            ii = minVal(minLogicalExtents(b,:))
            minLogicalExtents(b,1:nBlocksUsed) = &
                minLogicalExtents(b,1:nBlocksUsed) - ii
        end do
    end if
    !
    localSize(1) = 1
    localSize(2) = nPlotDim
    localSize(3) = nBlocksUsed
    globalSize(1) = 2
    globalSize(2) = nPlotDim
    globalSize(3) = nBlkUsedGlobal
    offset(3) = offsetPerProc(iProc)
    call writeHdf5DataSetHS(fileID,H5T_NATIVE_INTEGER,2,&
    "MinLogicalExtents",localSize(2:3), globalSize(2:3),offset(2:3),&
     .true., .true., hasData, dataspace,memspace, dataset,&
    plist_id,R2Int=minLogicalExtents)

    deallocate(minLogicalExtents)

    !end if



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

    deallocate(unknownNameArray)


    xmin = huge(xmin)
    ymin = huge(ymin)
    zmin = huge(zmin)
    xmax = -huge(xmax)
    ymax = -huge(ymax)
    zmax = -huge(zmax)
    localSize(1) = 1
    localSize(2) = nPlotDim
    localSize(3) = nBlocksUsed
    globalSize(1) = 2
    globalSize(2) = nPlotDim
    globalSize(3) = nBlkUsedGlobal
    offset(3) = offsetPerProc(iProc)
    if (isNonCartesian) then

        if (nBlocksUsed > 0) then

           allocate(coordinates(nPlotDim, nBlocksUsed))
            do iBlk = 1, nBlocksUsed
                    !!!FIX!!
                    !to speed things up only test outer nodes
                    do i = 1, nDim
                        coordinate(i) = huge(coordinate(i))
                    end do
                    do kk=1, nodeArrSize(3); do jj=1, nodeArrSize(2); do ii=1, nodeArrSize(1)
                        if (plotXYZNodes(1,ii,jj,kk,iBlk) < coordinate(1))&
                            coordinate(1) = plotXYZNodes(1,ii,jj,kk,iBlk)
                        if (nPlotDim .ge. 2) then
                            if (plotXYZNodes(2,ii,jj,kk,iBlk) < coordinate(2))&
                                coordinate(2) = plotXYZNodes(2,ii,jj,kk,iBlk)
                        end if
                        if (nPlotDim == 3) then
                            if (plotXYZNodes(3,ii,jj,kk,iBlk) < coordinate(3))&
                                coordinate(3) = plotXYZNodes(3,ii,jj,kk,iBlk)
                        end if
                    end do; end do; end do
                    coordinates(1:nPlotDim,iBlk) = coordinate(1:nPlotDim)
                    do i = 1,nDim
                        if (plotDim(i) == 1) then
                            if (coordinate(i) < xmin)&
                                xmin = coordinate(i)
                        else if (plotDim(i) == 2) then
                            if (coordinate(i) < ymin)&
                                ymin = coordinate(i)
                        else if (plotDim(i) == 3) then
                            if (coordinate(i) < zmin)&
                                zmin = coordinate(i)
                        else
                            if (nPlotDim == 2) then
                                if (plotDim(1) == 1 .and. plotDim(2) == 2) then
                                     if (coordinate(i) < zmin)&
                                        zmin = coordinate(i)
                                else if (plotDim(1) == 1 .and. plotDim(2) == 3) then
                                     if (coordinate(i) < ymin)&
                                        ymin = coordinate(i)
                                else if (plotDim(1) == 2 .and. plotDim(2) == 3) then
                                     if (coordinate(i) < xmin)&
                                        xmin = coordinate(i)
                                end if
                            end if
                        end if
                    end do
                end do
!                             else !1D
!                               if (plotDim(1) == 1) then
!                                 if (i == 2) then
!                                     if (coordinate(i) < ymin)&
!                                         ymin = coordinate(i)
!                                 else
!                                     if (coordinate(i) < zmin)&
!                                     zmin = coordinate(i)
!                                 end if
!                               else if (plotDim(1) == 2) then
!                                 if (i == 2) then
!                                     if (coordinate(i) < xmin)&
!                                         xmin = coordinate(i)
!                                 else
!                                     if (coordinate(i) < zmin)&
!                                         zmin = coordinate(i)
!                                 end if
!                              if (plotDim(1) == 3) then
!                                 if (i == 2) then
!                                     if (coordinate(i) < xmin) &
!                                         xmin = coordinate(i)
!                                 else
!                                     if (coordinate(i) < ymin)&
!                                         ymin = coordinate(i)
!                                 end if
!                             end if
!                         end if
!                     end if
!                 end do
!             end do
! 
           call writeHdf5DataSetHS(fileID,H5T_NATIVE_DOUBLE,3,"bounding box",&
            localSize(1:3), globalSize(1:3), offset(1:3),.false.,&
            .true. , hasData,dataspace,memspace, dataset,plist_id,R2Real=&
            coordinates)
            do iBlk = 1, nBlocksUsed
                    !!!FIX!!
                    !to speed things up only test outer nodes
                    do i = 1, nDim
                        coordinate(i) = -huge(coordinate(i))
                    end do
                    do kk=1, nodeArrSize(3); do jj=1, nodeArrSize(2); do ii=1, nodeArrSize(1)
                        if (plotXYZNodes(1,ii,jj,kk,iBlk) > coordinate(1))&
                            coordinate(1) = plotXYZNodes(1,ii,jj,kk,iBlk)
                        if (nPlotDim .ge. 2) then
                            if (plotXYZNodes(2,ii,jj,kk,iBlk) > coordinate(2))&
                                coordinate(2) = plotXYZNodes(2,ii,jj,kk,iBlk)
                        end if
                        if (nPlotDim == 3) then
                            if (plotXYZNodes(3,ii,jj,kk,iBlk) > coordinate(3))&
                                coordinate(3) = plotXYZNodes(3,ii,jj,kk,iBlk)
                        end if
                    end do; end do; end do
                    coordinates(1:nPlotDim,iBlk) = coordinate(1:nPlotDim)
                    do i = 1,nDim
                        if (plotDim(i) == 1) then
                            if (coordinate(i) > xmax)&
                                xmax = coordinate(i)
                        else if (plotDim(i) == 2) then
                            if (coordinate(i) > ymax)&
                                ymax = coordinate(i)
                        else if (plotDim(i) == 3) then
                            if (coordinate(i) > zmax)&
                                zmax = coordinate(i)
                        else
                            if (nPlotDim == 2) then
                                if (plotDim(1) == 1 .and. plotDim(2) == 2) then
                                     if (coordinate(i) > zmax)&
                                        zmax = coordinate(i)
                                else if (plotDim(1) == 1 .and. plotDim(2) == 3) then
                                     if (coordinate(i) > ymax)&
                                        ymax = coordinate(i)
                                else if (plotDim(1) == 2 .and. plotDim(2) == 3) then
                                     if (coordinate(i) > xmax)&
                                        xmax = coordinate(i)
                                end if
                            end if
                        end if
                    end do
                end do
!                             else !1D
!                               if (plotDim(1) == 1) then
!                                 if (i == 2) then
!                                     if (coordinate(i) > ymax)&
!                                         ymax = coordinate(i)
!                                 else
!                                     if (coordinate(i) > zmax)&
!                                     zmax = coordinate(i)
!                                 end if
!                               else if (plotDim(1) == 2) then
!                                 if (i == 2) then
!                                     if (coordinate(i) > xmax)&
!                                         xmax = coordinate(i)
!                                 else
!                                     if (coordinate(i) > zmax)&
!                                         zmax = coordinate(i)
!                                 end if
!                              if (plotDim(1) == 3) then
!                                 if (i == 2) then
!                                     if (coordinate(i) > xmax) &
!                                         xmax = coordinate(i)
!                                 else
!                                     if (coordinate(i) > ymax)&
!                                         ymax = coordinate(i)
!                                 end if
!                             end if
!                         end if
!                     end if
!                 end do
!             end do
! 
            offset(1) = 1
           call writeHdf5DataSetHS(fileID,H5T_NATIVE_DOUBLE,3,"bounding box",&
            localSize(1:3), globalSize(1:3),offset(1:3),.true.,&
            .false. ,hasData, dataspace,memspace, dataset,plist_id,R2Real=&
            coordinates)
            deallocate(coordinates)

        else
            call writeHdf5DataSetHS(fileID,H5T_NATIVE_DOUBLE,3,"bounding box",&
            localSize(1:3), globalSize(1:3),offset(1:3),.false.,&
            .true. ,hasData, dataspace,memspace, dataset,plist_id)
            offset(1) = 1
             call writeHdf5DataSetHS(fileID,H5T_NATIVE_DOUBLE,3,"bounding box",&
            localSize(1:3), globalSize(1:3),offset(1:3),.true.,&
            .false. ,hasData, dataspace,memspace, dataset,plist_id)
        end if
        deallocate(plotXYZNodes)
    else
        if (nBlocksUsed > 0) then
            if (plot_dimensional) then
                call writeHdf5DataSetHS(fileID,H5T_NATIVE_DOUBLE,3,"bounding box",&
                localSize(1:3), globalSize(1:3),offset(1:3),.false.,&
                .true. , hasData,dataspace,memspace, dataset,plist_id,R2Real=&
                CoordMin_DB(plotDim(1:nPlotDim),UsedBlocks)*No2Io_V(UnitX_))
                offset(1) = 1
                call writeHdf5DataSetHS(fileID,H5T_NATIVE_DOUBLE,3,"bonding box",&
                localSize(1:3), globalSize(1:3),offset(1:3),.true.,&
                .false. ,hasData, dataspace,memspace, dataset,plist_id,R2Real=&
                CoordMax_DB(plotDim(1:nPlotDim),UsedBlocks)*No2Io_V(UnitX_))
                xmax = maxval(CoordMax_DB(1,UsedBlocks))*No2Io_V(UnitX_)
                xmin = minval(CoordMin_DB(1,UsedBlocks))*No2Io_V(UnitX_)
                ymax = maxval(CoordMax_DB(2,UsedBlocks))*No2Io_V(UnitX_)
                ymin = minval(CoordMin_DB(2,UsedBlocks))*No2Io_V(UnitX_)
                zmax = maxval(CoordMax_DB(3,UsedBlocks))*No2Io_V(UnitX_)
                zmin = minval(CoordMin_DB(3,UsedBlocks))*No2Io_V(UnitX_)
 
            !deallocate(coordinates)) 
            else
                call writeHdf5DataSetHS(fileID,H5T_NATIVE_DOUBLE,3,"bounding box",&
                localSize(1:3), globalSize(1:3), offset(1:3),.false.,&
                .true. , hasData,dataspace,memspace, dataset,plist_id,R2Real=&
                CoordMin_DB(plotDim(1:nPlotDim),UsedBlocks))
                offset(1) = 1
                call writeHdf5DataSetHS(fileID,H5T_NATIVE_DOUBLE,3,"bonding box",&
                localSize(1:3), globalSize(1:3),offset(1:3),.true.,&
                .false. ,hasData, dataspace,memspace, dataset,plist_id,R2Real=&
                CoordMax_DB(plotDim(1:nPlotDim),UsedBlocks))
            end if
                xmax = maxval(CoordMax_DB(1,UsedBlocks))
                xmin = minval(CoordMin_DB(1,UsedBlocks))
                ymax = maxval(CoordMax_DB(2,UsedBlocks))
                ymin = minval(CoordMin_DB(2,UsedBlocks))
                zmax = maxval(CoordMax_DB(3,UsedBlocks))
                zmin = minval(CoordMin_DB(3,UsedBlocks))
 
        else
            call writeHdf5DataSetHS(fileID,H5T_NATIVE_DOUBLE,3,"bounding box",&
            localSize(1:3), globalSize(1:3),offset(1:3),.false.,&
            .true. ,hasData,dataspace,memspace, dataset,plist_id)

            offset(1) = 1
            call writeHdf5DataSetHS(fileID,H5T_NATIVE_DOUBLE,3,"bonding box",&
            localSize(1:3), globalSize(1:3),offset(1:3),.true.,&
            .false. ,hasData, dataspace,memspace, dataset,plist_id)
        end if
 
    end if
    tempReal = xmin
    call MPI_allreduce(tempReal, xmin, 1, MPI_REAL,&
        MPI_MIN, iComm, error)
    tempReal = ymin
    call MPI_allreduce(tempReal, ymin, 1, MPI_REAL,&
        MPI_MIN, iComm, error)
    tempReal = zmin
    call MPI_allreduce(tempReal, zmin, 1, MPI_REAL,&
        MPI_MIN, iComm, error)
     tempReal = xmax
    call MPI_allreduce(tempReal, xmax, 1, MPI_REAL,&
        MPI_MAX, iComm, error)
    tempReal = ymax
    call MPI_allreduce(tempReal, ymax, 1, MPI_REAL,&
        MPI_MAX, iComm, error)
    tempReal = zmax
    call MPI_allreduce(tempReal, zmax, 1, MPI_REAL,&
        MPI_MAX, iComm, error)

! 

    allocate(coordinates(nPlotDim,nBlocksUsed))
    if (nBlocksUsed > 0) then
    do iBlk = 1,nBlocksUsed
       ii = 1
       if (iSizeGlobal .ne. 1) then
            if (plot_dimensional) then
                coordinates(ii,iBlk) = 0.5*(x_BLK(nI,nJ,nK,UsedBlocks(iBLK))+&
                x_BLK(1,1,1,UsedBlocks(iBLK)))*No2Io_V(UnitX_)
            else
                coordinates(ii,iBlk) = 0.5*(x_BLK(nI,nJ,nK,UsedBlocks(iBLK))+&
                x_BLK(1,1,1,UsedBlocks(iBLK)))
            end if
            ii = ii+1
       end if
       if (jSizeGlobal .ne. 1) then
            if (plot_dimensional) then
                coordinates(ii,iBlk) = 0.5*(y_BLK(nI,nJ,nK,&
                UsedBlocks(iBLK))+ y_BLK(1,1,1,UsedBlocks(iBLK)))*No2Io_V(UnitX_)
            else
                coordinates(ii,iBlk) = 0.5*(y_BLK(nI,nJ,nK,&
                UsedBlocks(iBLK))+ y_BLK(1,1,1,UsedBlocks(iBLK)))
            end if

            ii = ii + 1
        end if
       if (kSizeGlobal .ne. 1) then
            if (plot_dimensional) then
                coordinates(ii,iBlk) = 0.5*&
                (z_BLK(nI,nJ,nK,UsedBlocks(iBLK))+&
                z_BLK(1,1,1,UsedBlocks(iBLK)))*No2Io_V(UnitX_)
            else
                coordinates(ii,iBlk) = 0.5*&
                (z_BLK(nI,nJ,nK,UsedBlocks(iBLK))+&
                z_BLK(1,1,1,UsedBlocks(iBLK)))
            end if
        end if
    end do
    else
        coordinates = 0
    end if
    !



        call writeHdf5DataSetHS(fileID,H5T_NATIVE_DOUBLE,2,"coordinates",&
        localSize(2:3), globalSize(2:3),offset(2:3),.true.,&
        .true. , hasData, dataspace,memspace, dataset,plist_id,R2Real=&
        coordinates)


    deallocate(coordinates)

    
    if (nBlocksUsed > 0) then
    call writeHdf5DataSetHS(fileID,H5T_NATIVE_INTEGER,1,"iMortonNode_A",&
    localSize(3), globalSize(3), offset(3),.true.,&
    .true.,hasData, dataspace,memspace,dataset,plist_id,&
    R1Int=iMortonNode_A(usedNodes))
    else
    call writeHdf5DataSetHS(fileID,H5T_NATIVE_INTEGER,1,"iMortonNode_A",&
    localSize(3), globalSize(3), offset(3),.true.,&
    .true.,hasData,dataspace,memspace,dataset,plist_id)
    end if

    deallocate(usedBlocks)
    if (nBlocksUsed > 0) then
    call writeHdf5DataSetHS(fileID,H5T_NATIVE_INTEGER,1,"refine level",&
    localSize(3), globalSize(3), offset(3),.true.,&
    .true.,hasData,dataspace,memspace,dataset,plist_id,R1Int=&
    iTree_IA(Level_,UsedNodes))
     else
     call writeHdf5DataSetHS(fileID,H5T_NATIVE_INTEGER,1,"refine level",&
    localSize(3), globalSize(3), offset(3),.true.,&
    .true.,hasData,dataspace,memspace,dataset,plist_id)
     end if
     deallocate(usedNodes)
    
    allocate(procNum(nBlocksUsed))
    procNum = iProc+1
    call writeHdf5DataSetHS(fileID,H5T_NATIVE_INTEGER,1,"Processor Number",&
    localSize(3), globalSize(3), offset(3), .true.,&
    .true.,hasData,dataspace,memspace,dataset,plist_id,R1Int=procNum)

   deallocate(procNum)
    call writeHdf5DataSetHS(fileID,H5T_NATIVE_INTEGER,1,"VisIt Index",&
    localSize(3), globalSize(3), offset(3), .true.,&
    .true.,hasData,dataspace,memspace,dataset,plist_id,R1Int=VisItIdxLocal)


    deallocate(VisItIdxLocal)
    call write_integer_sim_metadata(fileID, nPlotVar)
    call write_real_sim_metadata(FileID,plot_dimensional)
    call write_integer_plot_metadata(fileID, nPlotVar,isCutFile)

    call write_real_plot_metadata(FileID,xmin,xmax,ymin,ymax, zmin,zmax)
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

  !=================================================================
  !===============================================================``=


  subroutine writeHdf5DataSetHS(fileID,datatype,rank,description,localSize,&
    globalSize, offset, closeDataset,openDataset,hasData,dataspace,&
    memspace,dataset,plist_id,R1Int, R1Real,R2Int, R2Real, R3Int, &
    R3Real,R4Int,R4Real,R5Real)
    use ModProcMH, only : nProc, iProc
    implicit none

    integer, intent(in) :: rank 
    integer :: error
!     integer(HSIZE_T), intent(in) :: localDimens(4), globalDimens(4)
!     integer(HSIZE_T), intent(in) :: startArr(4), strideArr(4), countArr(4)
    integer(HID_T), intent(inout) :: dataspace
    integer(HID_T), intent(inout) :: memspace
    integer(HID_T), intent(inout) :: dataset, plist_id
    integer(HID_T), intent(in) :: fileID, datatype
    integer(HSIZE_T), intent(in) :: localSize(rank), globalSize(rank), offset(rank)
    integer(HSIZE_T) :: nullSize(rank), one, strideArr(rank)
    character (len=*), intent(in) :: description
    logical, intent(in) :: closeDataset
    logical, intent(in) :: openDataset, hasData
    real, intent(in), optional :: R1Real(:)
    integer, intent(in), optional :: R1Int(:)
    real, intent(in), optional :: R2Real(:,:)
    integer, intent(in), optional :: R2Int(:,:)
    real, intent(in), optional :: R3Real(:,:,:)
    integer, intent(in), optional :: R3Int(:,:,:)
    real, intent(in), optional :: R4Real(:,:,:,:)
    integer, intent(in), optional :: R4Int(:,:,:,:)
    real, intent(in), optional :: R5Real(:,:,:,:,:)
    integer, pointer :: nullPointer => NULL( )
    
    strideArr(1:rank) = 1
    !Set the dimensions of the dataset

    call h5open_f(error)
    if (error == -1) &
         
         call stop_mpi("error in subroutine writeHdf5Rank4Real. Error marker 1")
    
    if (openDataset) then  
        call h5screate_simple_f(rank, globalSize, dataspace, error) 
        !create the dataset

        call h5dcreate_f(fileID, description, datatype, dataspace,&
        dataset, error)

       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)

       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    end if

        if (error == -1) &
             call stop_mpi("error in subroutine writeHdf5Rank4Real. Error marker 4")




    if (hasData) then
    call h5screate_simple_f(rank, localSize, memspace, error)

!      call h5sselect_elements_f(dataspace, H5S_SELECT_SET_F, rank, &
!          numPoints, pSelect, error)
    call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, localSize, error, &
            stride = strideArr)
    

     else
!    call h5screate_f(H5S_NULL_F, memspace, error)
!    if (nBlocksUsed == 0) then
! !    call h5sset_extent_none_f(dataspace, error)
!     call h5sset_extent_none_f(memspace, error)
    one = 0
    call h5screate_simple_f(1, (/one/), memspace, error)

    call h5sselect_none_f(memspace,error)
    call h5sselect_none_f(dataspace,error)
    nullSize(1:rank) = 0
    call h5dwrite_f(dataset, datatype, nullPointer, nullSize, error,&
    mem_space_id = memspace,file_space_id = dataspace, &
    xfer_prp = plist_id)
    if (closeDataset)&
        call hdf5_close(memspace, plist_id, dataspace, dataset) 
    return
    end if

!    if (nBlocksUsed > 0) then
    if (Present(R1Real)) then
        call h5dwrite_f(dataset, datatype, R1Real, localSize, error,&
        mem_space_id = memspace,file_space_id = dataspace, &
        xfer_prp = plist_id)    
    elseif(Present(R1Int)) then
        call h5dwrite_f(dataset, datatype, R1Int, localSize, error,&
         mem_space_id = memspace,file_space_id = dataspace,&
         xfer_prp = plist_id)
    elseif(Present(R2Real)) then
        CALL H5dwrite_f(dataset, datatype, R2Real, localSize, error, &
            mem_space_id=memspace, file_space_id=dataspace,&
         xfer_prp = plist_id)
    elseif(Present(R2Int)) then
        call h5dwrite_f(dataset, datatype, R2Int, localSize, error,&
        mem_space_id = memspace,file_space_id = dataspace, &
        xfer_prp = plist_id)
    elseif(Present(R3Real)) then
        call h5dwrite_f(dataset, datatype, R3Real, localSize, error,&
        mem_space_id = memspace,file_space_id = dataspace, &
        xfer_prp = plist_id)
    elseif(Present(R3Int)) then
        call h5dwrite_f(dataset, datatype, R3Int, localSize, error,&
        mem_space_id = memspace,file_space_id = dataspace, &
        xfer_prp = plist_id)
    elseif(Present(R4Real)) then
        call h5dwrite_f(dataset, datatype, R4Real, localSize, error,&
        mem_space_id = memspace,file_space_id = dataspace, &
        xfer_prp = plist_id)
    elseif(Present(R4Int)) then
        call h5dwrite_f(dataset, datatype, R4Int, localSize, error,&
        mem_space_id = memspace,file_space_id = dataspace, &
        xfer_prp = plist_id)
    elseif(Present(R5Real)) then
        call h5dwrite_f(dataset, datatype, R5Real, localSize, error,&
        mem_space_id = memspace,file_space_id = dataspace, &
        xfer_prp = plist_id)
    end if

!     else 
!          call h5dwrite_f(dataset, datatype, char(0), localSize, error,&
!         mem_space_id = memspace,file_space_id = dataspace, &
!         xfer_prp = plist_id)
!    
!     end if
    if (closeDataset)&
        call hdf5_close(memspace, plist_id, dataspace, dataset) 
       
  end subroutine writeHdf5DataSetHS

  !=================================================================
  !===============================================================``=

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
    do i=1,3
        if (plotDim(i) == 1) then 
            RealMetaData(iData) = xMin
            iData = iData + 1
            RealMetaData(iData) = xMax
            iData = iData + 1
         else if (plotDim(i) == 2) then
            RealMetaData(iData) = yMin
            iData = iData + 1
            RealMetaData(iData) = yMax
            iData = iData + 1
        else if (plotDim(i) == 3) then
            RealMetaData(iData) = zMin
            iData = iData + 1
            RealMetaData(iData) = zMax
            iData = iData + 1
        else if (plotDim(i) == 0) then
             RealMetaData(iData) = 0
            iData = iData + 1
            RealMetaData(iData) = 0
            iData = iData + 1
        end if
    end do
    iData = iData - 1
    
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
  !=================================================================
  !================================================================

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
  end subroutine write_real_sim_metadata
  !=================================================================
  !===============================================================``=

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

    if (iProc == 0)&   
         call h5dwrite_f(dataset,H5T_NATIVE_INTEGER,IntegerMetaData,dimens1D,&
         error, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
    if (error == -1)&
         call stop_mpi(&
         "error in subroutine init_grid. Error marker 3")

    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)
  end subroutine write_integer_plot_metadata

  !=================================================================
  !===============================================================``=
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

    if (iProc == 0)&   
         call h5dwrite_f(dataset,H5T_NATIVE_INTEGER,IntegerMetaData,dimens1D,&
         error, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
    if (error == -1)&
         call stop_mpi(&
         "error in subroutine init_grid. Error marker 3")

    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)
  end subroutine write_integer_sim_metadata



  !=====================================================================
  !=====================================================================

  subroutine write_hdf5_attribute(attName, dataset, datatype, realAtt,&
    intAtt)
    integer(HID_T), intent(in) :: dataset, datatype
    character (len=*), intent(in) :: attName
    real, optional, intent(in) :: realAtt
    real, optional, intent(in) :: intAtt
    integer(HID_T) :: attributeSpace, attribute
    integer(HSIZE_T) :: dimens1D(1)
    integer :: error

       dimens1D=(1)
       call h5screate_simple_f(1, dimens1D, attributeSpace, error)

       call h5acreate_f(dataset, attname, datatype,&
            attributeSpace, attribute,&
            error, H5P_DEFAULT_F)
       if (Present(realAtt)) then
           call h5awrite_f(attribute, datatype, realAtt, dimens1D, error)
       elseif (present(intAtt)) then
           call h5awrite_f(attribute, datatype, realAtt, dimens1D, error)
       endif
       call h5sclose_f(attributeSpace, error)
       call h5aclose_f(attribute, error)

  end subroutine write_hdf5_attribute
  !=====================================================================
  !=====================================================================

  !======================================================================
  !=====================================================================
  subroutine hdf5_close(memspace, plist_id, dataspace, dataset) 
    integer(HID_T), intent(in) ::memspace, plist_id, dataspace, dataset
    integer :: error
    call h5sclose_f(memspace,error)
    call h5pclose_f(plist_id, error)
    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)
  end subroutine hdf5_close
  !======================================================================
   !=====================================================================

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
 ! 
  !=================================================================
  !================================================================

!   subroutine get_var_extrema_hdf5( iVar, globalVarMin, globalVarMax)
! 
!     use ModMpi
! 
!     ! Returns the global variable extrema for each unkown.
! 
!     integer, intent(in) :: iVar
!     real, intent(out) :: globalVarMin, globalVarMax
! 
!     real :: varMin, varMax
!     integer :: error
!     !------------------------------------------------------------------
! 
! 
!     varMin = huge(varMin)
!     varMax = -huge(varMax)
!     if (nBlocksUsed > 0) then
!     varMin = minval(PlotVarIdx(:,:,:,:,iVar))
!     varMax = maxval(PlotVarIdx(:,:,:,:,iVar))
!     else
! 
!     end if
! 
! 
!   end subroutine get_var_extrema_hdf5
!   !===============================================================
!   !=============================================================== 
!        
!   logical function sort_test_cells(i, j)
! 
!     integer, intent(in) :: i,j
!     if(ndim == 3) then
!         sort_test_cells = 0
!         if(Position_DA(1,i) < Position_DA(1,j)) then
!            sort_test_cells = .false.
!         elseif(Position_DA(1,j) < Position_DA(1,i)) then
!            sort_test_cells = .true.
!         else
!            if(Position_DA(2,i) < Position_DA(2,j)) then
!               sort_test_cells = .false.
!            elseif(Position_DA(2,j) < Position_DA(2,i)) then
!               sort_test_cells = .true.
!            else
!               if(Position_DA(3,i) < Position_DA(3,j)) then
!                  sort_test_cells = .false.
!               elseif(Position_DA(3,j) < Position_DA(3,i)) then
!                  sort_test_cells = .true.
!               end if
!            end if
!         end if
!     elseif(nDim == 2) then
!        if(Position_DA(1,i) < Position_DA(1,j)) then
!           sort_test_cells = .false.
!        elseif(Position_DA(1,j) < Position_DA(1,i)) then
!           sort_test_cells = .true.
!        else
!           if(Position_DA(2,i) < Position_DA(2,j)) then
!              sort_test_cells = .false.
!           elseif(Position_DA(2,j) < Position_DA(2,i)) then
!              sort_test_cells = .true.
!           end if
!        end if
!     elseif(nDim == 1) then
!         if(Position_DA(1,i) < Position_DA(1,j)) then
!             sort_test_cells = .false.
!         elseif(Position_DA(1,j) < Position_DA(1,i)) then
!              sort_test_cells = .true.
!         end if
! 
!     else
!         sort_test_cells = .false.
!     end if
!     end function sort_test_cells
!   !===============================================================
!   !=============================================================== 
!        
  logical function sort_test_block(i, j)
!    use BATL_lib, only : nDim
    integer, intent(in) :: i,j
!    if(nDim == 3) then
     if(nPlotDim == 3) then  
        sort_test_block = 0
        if(Position_DA(1,i) < Position_DA(1,j)) then
           sort_test_block = .false.
        elseif(Position_DA(1,j) < Position_DA(1,i)) then
           sort_test_block = .true.
        else
           if(Position_DA(2,i) < Position_DA(2,j)) then
              sort_test_block = .false.
           elseif(Position_DA(2,j) < Position_DA(2,i)) then
              sort_test_block = .true.
           else
              if(Position_DA(3,i) < Position_DA(3,j)) then
                 sort_test_block = .false.
              elseif(Position_DA(3,j) < Position_DA(3,i)) then
                 sort_test_block = .true.
              end if
           end if
        end if
!    elseif(nDim == 2) then
    elseif(nPlotDim == 2) then
       if(Position_DA(2,i) < Position_DA(2,j)) then
          sort_test_block = .false.
       elseif(Position_DA(2,j) < Position_DA(2,i)) then
          sort_test_block = .true.
       else
          if(Position_DA(1,i) < Position_DA(1,j)) then
             sort_test_block = .false.
          elseif(Position_DA(1,j) < Position_DA(1,i)) then
             sort_test_block = .true.
          end if
       end if
!    elseif(nDim == 1) then
    elseif(nPlotDim == 1) then
        if(Position_DA(1,i) < Position_DA(1,j)) then
            sort_test_block = .false.
        elseif(Position_DA(1,j) < Position_DA(1,i)) then
             sort_test_block = .true.
        end if

    else
        sort_test_block = .false.
    end if
    end function sort_test_block

end module ModHdf5
