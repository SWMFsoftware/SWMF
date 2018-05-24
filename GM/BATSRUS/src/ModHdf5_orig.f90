!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModHdf5

  use BATL_lib, ONLY: &
       test_start, test_stop

  use ModSize, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK
  use hdf5
  use ModHdf5Utils
  use ModNumConst
  use BATL_size, ONLY: nDim
  use ModMpi
  implicit none
  SAVE

  private ! except

  public:: write_var_hdf5
  public:: write_plot_hdf5
  public:: init_hdf5_plot

  real, parameter:: cHalfMinusTiny = 0.5*(1.0 - 1e-6)

  real, allocatable :: PlotVarIdx(:,:,:,:,:)
  real, allocatable :: PlotXYZNodes(:,:,:,:,:)
  real, allocatable :: CellData(:,:)

  integer, allocatable :: UsedBlocks(:)
  integer, allocatable :: usednodes(:)
  integer, allocatable :: BlockProcNum(:)
  integer, allocatable :: BlocksPerProc(:)
  integer, allocatable :: MinLogicalExtents(:,:)
  integer, allocatable :: CellProc(:,:,:)
  integer(HSIZE_T), allocatable :: CoordArray(:,:)

  integer(HID_T), allocatable  :: VarDataset(:), VarPlistId(:), VarDataspace(:)

  integer(HID_T) :: fileID
  character (len=lNameH5), allocatable :: UnknownNameArray(:)

  integer ::  iBlock
  integer :: iSizeGlobal, jSizeGlobal, kSizeGlobal, iPlotDim(3)
  integer :: nMessagesNeeded, nMessagesTotal, nBlocksThetaPhi(2), iProcOffset,nBlocksUsed
  integer(HSIZE_T) :: nCellsLocal,nCellsPerBlock(3), nBlkUsedGlobal, nPlotDim
  integer(HSIZE_T) :: nNodeCellsPerBlock(3), nLogicalDim
  integer(HSIZE_T) :: iInteger8, iInteger8a
  integer(HID_T) :: iInteger4

  integer :: iplotComm,nPlotProc,iPlotProc ! one comunicator for each plot

contains
  !============================================================================
  subroutine stop_hdf5(String)
    character(len=*), intent(in):: String
    integer:: iError, nError

    !--------------------------------------------------------------------------
    write(*,*)'ERROR in ModHdf5: ', String
    call MPI_abort(MPI_COMM_WORLD, nError, iError)
    stop

  end subroutine stop_hdf5
  !============================================================================
  subroutine init_hdf5_plot(iFile, plotType, nPlotVar, &
       xMin, xMax, yMin, yMax, zMin, zMax, DxBlock, DyBlock, DzBlock,&
       isNonCartesian, NotACut)

    ! Save all cells within plotting range, for each processor

    use ModProcMH
    use ModMain, ONLY: nI, nJ, nK, &
         x_, y_, z_, Phi_, nBlockMax, Unused_B
    use ModGeometry, ONLY: CellSize_DB,&
         XyzStart_BLK
    use ModIO
    use BATL_lib, ONLY: Xyz_DNB, IsRLonLat, IsCylindrical, CoordMin_DB,&
         CoordMax_DB, CoordMin_D, CoordMax_D, DiLevelNei_IIIB

    ! Arguments

    integer, intent(in) :: iFile
    integer, intent(in) :: nPlotVar
    real,    intent(in) :: xMin,xMax,yMin,yMax,zMin,zMax
    logical, intent(in) :: isNonCartesian, NotACut
    character(len=*), intent(in) :: plotType
    real,    intent(inout):: DxBlock,DyBlock,DzBlock

    ! Local variables
    ! Indices and Coordinates
    integer :: i, j, k, iMin, iMax, jMin, jMax, kMin, kMax
    integer :: iMin1, iMax1, jMin1, jMax1, kMin1, kMax1
    integer :: iMaxN, jMaxN, kMaxN, iMaxN1, jMaxN1, kMaxN1
    integer :: ii, jj, kk, iSize, jSize, kSize
    integer :: Error,i1, j1, k1,ii1, jj1, kk1
    real    :: Dx
    real    :: xMin1, xMax1, yMin1, yMax1, zMin1, zMax1
    real    :: ySqueezed

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'init_hdf5_plot'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    iMin = huge(iMin)
    iMax = -huge(iMax)
    jMin = huge(jMin)
    jMax = -huge(jMax)
    kMin = huge(kMin)
    kMax = -huge(kMax)
    iSize = 0
    jSize = 0
    kSize = 0
    nBlocksUsed = 0
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE

       xMin1 = xMin - cHalfMinusTiny*CellSize_DB(x_,iBlock)
       xMax1 = xMax + cHalfMinusTiny*CellSize_DB(x_,iBlock)
       yMin1 = yMin - cHalfMinusTiny*CellSize_DB(y_,iBlock)
       yMax1 = yMax + cHalfMinusTiny*CellSize_DB(y_,iBlock)
       zMin1 = zMin - cHalfMinusTiny*CellSize_DB(z_,iBlock)
       zMax1 = zMax + cHalfMinusTiny*CellSize_DB(z_,iBlock)

       if(IsRLonLat .or. IsCylindrical)then
          ! Make sure that angles around 3Pi/2 are moved to Pi/2 for x=0 cut
          ySqueezed = mod(xyzStart_BLK(Phi_,iBlock), cPi)
          ! Make sure that small angles are moved to Pi degrees for y=0 cut
          if(ySqueezed < 0.25*cPi .and. &
               abs(yMin+yMax-cTwoPi) < cTiny .and. yMax-yMin < 0.01) &
               ySqueezed = ySqueezed + cPi
       else
          ySqueezed = xyzStart_BLK(y_,iBlock)
       end if

       if( xyzStart_BLK(x_,iBlock) > xMax1.or.&
            xyzStart_BLK(x_,iBlock)+(nI-1)*CellSize_DB(x_,iBlock) < xMin1.or.&
            ySqueezed > yMax1.or.&
            ySqueezed+(nJ-1)*CellSize_DB(y_,iBlock) < yMin1.or.&
            xyzStart_BLK(z_,iBlock) > zMax1.or.&
            xyzStart_BLK(z_,iBlock)+(nK-1)*CellSize_DB(z_,iBlock) < zMin1) CYCLE
       if (NotACut) then
          iMin = 1
          iMax = nI
          jMin = 1
          jMax = nJ
          kMin = 1
          kMax = nK
          nBlocksUsed = nBlocksUsed + 1
       else

          Dx = plot_Dx(1,iFile)
          DxBlock=CellSize_DB(x_,iBlock); DyBlock=CellSize_DB(y_,iBlock); DzBlock=CellSize_DB(z_,iBlock)

          ! Calculate index limits of cells inside cut
          i = max(1 ,floor((xMin1-xyzStart_BLK(x_,iBlock))/DxBlock)+2)
          ii = min(nI,floor((xMax1-xyzStart_BLK(x_,iBlock))/DxBlock)+1)

          j = max(1 ,floor((yMin1-ySqueezed)/DyBlock)+2)
          jj = min(nJ,floor((yMax1-ySqueezed)/DyBlock)+1)

          k = max(1 ,floor((zMin1-xyzStart_BLK(z_,iBlock))/DzBlock)+2)
          kk = min(nK,floor((zMax1-xyzStart_BLK(z_,iBlock))/DzBlock)+1)

          i1 = max(0 ,floor((xMin1-xyzStart_BLK(x_,iBlock))/DxBlock)+2)
          ii1 = min(nI+1,floor((xMax1-xyzStart_BLK(x_,iBlock))/DxBlock)+1)

          j1 = max(0 ,floor((yMin1-ySqueezed)/DyBlock)+2)
          jj1 = min(nJ+1,floor((yMax1-ySqueezed)/DyBlock)+1)

          k1 = max(0,floor((zMin1-xyzStart_BLK(z_,iBlock))/DzBlock)+2)
          kk1 = min(nK+1,floor((zMax1-xyzStart_BLK(z_,iBlock))/DzBlock)+1)

          !         if (i1==0 .and.  i==1 .and. ii==1) then ! Cut falls between iMinCell and iMinCell - 1
          !             if (DiLevelNei_IIIB(-1, 0, 0,iBlock) == -1) cycle
          !             ! Neighbor is more refined so it writes
          !         end if
          !
          !         if (j1==0 .and. j==1 .and. jj == 1) then
          !             ! Neighbor is more refined so it writes
          !             if (DiLevelNei_IIIB(0, -1, 0,iBlock) == -1) cycle
          !         end if
          !
          !         if (k1==0 .and. k==1 .and. kk == 1) then
          !             ! Neighbor is more refined so it writes
          !             if (DiLevelNei_IIIB(0, 0, -1,iBlock) == -1) cycle
          !         end if
          !
          !         if (ii1==nI+1 .and.  ii==nI .and. i == nI) then ! Cut falls between iMaxCell and iMaxCell 1
          !             if (DiLevelNei_IIIB(1, 0, 0,iBlock) == -1) then
          !                 cycle
          !             else if (CoordMax_DB(1, iBlock) /= CoordMax_D(1)) then
          !                 cycle
          !
          !             end if
          !         end if
          !
          !         if (jj1==nJ+1 .and. jj==nJ .and. j == nJ) then
          !             if (DiLevelNei_IIIB(0, 1, 0,iBlock) == -1) then
          !                 cycle
          !             else if (CoordMax_DB(2, iBlock) /= CoordMax_D(2)) then
          !                 cycle
          !
          !             end if
          !         end if
          !
          !         if (kk1==nK+1 .and. kk==nK .and. k==nK) then
          !             if (DiLevelNei_IIIB(0, 0, 1,iBlock) == -1) then
          !                 cycle
          !             else if (CoordMax_DB(3, iBlock) /= CoordMax_D(3)) then
          !                 cycle
          !
          !             end if
          !         end if
          !
          ! Neighbor is more refined so it writes

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

          nBlocksUsed = nBlocksUsed + 1
       end if
    end do

    ! Make a comincator for only procs which has data to save
    if(nBlock == 0) then
       call MPI_Comm_split(iComm,MPI_UNDEFINED,iProc,iplotComm,Error)
    else
       call MPI_Comm_split(iComm,1,iProc,iplotComm,Error)
    endif

    if(MPI_COMM_NULL == iplotComm ) RETURN
    call MPI_COMM_RANK (iplotComm, iPlotProc, Error)
    call MPI_COMM_SIZE (iplotComm, nPlotProc, Error)

    allocate(blocksPerProc(0:nPlotProc-1))
    call MPI_Allgather(nBlocksUsed, 1, MPI_INTEGER, blocksPerProc, &
         1, MPI_INTEGER, iplotComm, Error)
    if (iPlotProc == 0) then
       iProcOffset = 0
    else
       iProcOffset = sum(blocksPerProc(0:iPlotProc-1))
    end if
    nBlkUsedGlobal = sum(blocksPerProc(0:nPlotProc-1))
    deallocate(blocksPerProc)

    allocate(UsedBlocks(nBlocksUsed))
    if(NotACut) then
       iSizeGlobal = nI
       jSizeGlobal = nJ
       kSizeGlobal = nK
    else
       call MPI_AllReduce(iSize, iSizeGlobal, 1,&
            MPI_INTEGER, MPI_MAX, iplotComm, Error)
       call MPI_AllReduce(jSize, jSizeGlobal, 1,&
            MPI_INTEGER, MPI_MAX, iplotComm, Error)
       call MPI_AllReduce(kSize, kSizeGlobal, 1,&
            MPI_INTEGER, MPI_MAX, iplotComm, Error)
    end if

    if (isNonCartesian) then
       if (NotACut) then
          nPlotDim = nDim
          iPlotDim = 0
          nCellsPerBlock = 1
          nNodeCellsPerBlock = 1
          do i=1,nPlotDim
             iPlotDim(i) = i
             if (i == 1) then
                nCellsPerBlock(i) = iSizeGlobal
                nNodeCellsPerBlock(i) = iSizeGlobal+1
             else if (i == 2) then
                nCellsPerBlock(i) = jSizeGlobal
                nNodeCellsPerBlock(i) = jSizeGlobal+1
             else if (i == 3) then
                nCellsPerBlock(i) = kSizeGlobal
                nNodeCellsPerBlock(i) = kSizeGlobal+1
             end if
          end do
       else if (plotType=='x=0') then
          !              write (*,*) "allocated x=0"
          nPlotDim = 2
          iPlotDim = (/1,3,0/)
          nCellsPerBlock = (/iSizeGlobal, kSizeGlobal, 1/)
          nNodeCellsPerBlock = (/iSizeGlobal+1, kSizeGlobal+1, 1/)
       else if(plotType=="y=0)") then
          !              write (*,*) "allocated y=0"
          nPlotDim = 2
          iPlotDim = (/1,3,0/)
          nCellsPerBlock = (/iSizeGlobal, kSizeGlobal, 1/)
          nNodeCellsPerBlock = (/iSizeGlobal+1, kSizeGlobal+1, 1/)
       else if (plotType=='z=0' )then
          !              write (*,*) "allocated z=0"
          nPlotDim = 2
          iPlotDim = (/1,2,0/)
          nCellsPerBlock = (/iSizeGlobal, jSizeGlobal, 1/)
          nNodeCellsPerBlock = (/iSizeGlobal+1, jSizeGlobal+1, 1/)
       end if
    else
       if (NotACut) then
          nPlotDim = nDim
          iPlotDim = 0
          nCellsPerBlock = 1
          nNodeCellsPerBlock = 1
          do i=1,nPlotDim
             iPlotDim(i) = i
             if (i == 1) then
                nCellsPerBlock(i) = iSizeGlobal
                nNodeCellsPerBlock(i) = iSizeGlobal+1
             else if (i == 2) then
                nCellsPerBlock(i) = jSizeGlobal
                nNodeCellsPerBlock(i) = jSizeGlobal+1
             else if (i == 3) then
                nCellsPerBlock(i) = kSizeGlobal
                nNodeCellsPerBlock(i) = kSizeGlobal+1
             end if
          end do

       else if (iSizeGlobal == 1 .and. jSizeGlobal == 1) then
          nPlotDim = 1
          iPlotDim = (/3,0,0/)
          nCellsPerBlock = (/kSizeGlobal,1,1/)
          nNodeCellsPerBlock = (/kSizeGlobal+1,1,1/)
       else if (iSizeGlobal == 1 .and. kSizeGlobal == 1) then
          nPlotDim = 1
          iPlotDim = (/2,0,0/)
          nCellsPerBlock = (/jSizeGlobal,1,1/)
          nNodeCellsPerBlock = (/jSizeGlobal+1,1,1/)
       else if (jSizeGlobal == 1 .and. kSizeGlobal == 1) then
          nPlotDim = 1
          iPlotDim = (/1,0,0/)
          nCellsPerBlock = (/iSizeGlobal,1,1/)
          nNodeCellsPerBlock = (/iSizeGlobal+1,1,1/)
       else if (iSizeGlobal == 1) then
          nPlotDim = 2
          iPlotDim = (/2,3,0/)
          nCellsPerBlock = (/jSizeGlobal,kSizeGlobal,1/)
          nNodeCellsPerBlock = (/jSizeGlobal+1,kSizeGlobal+1,1/)
       else if (jSizeGlobal == 1) then
          nPlotDim = 2
          iPlotDim = (/1,3,0/)
          nCellsPerBlock = (/iSizeGlobal,kSizeGlobal,1/)
          nNodeCellsPerBlock = (/iSizeGlobal+1,kSizeGlobal+1,1/)
       else if (kSizeGlobal == 1) then
          nPlotDim = 2
          iPlotDim = (/1,2,0/)
          nCellsPerBlock = (/iSizeGlobal,jSizeGlobal,1/)
          nNodeCellsPerBlock = (/iSizeGlobal+1,jSizeGlobal+1,1/)
       end if
    end if

    allocate(PlotVarIdx(nCellsPerBlock(1),nCellsPerBlock(2), nCellsPerBlock(3), nBlocksUsed, nPlotVar))
    PlotVarIdx=-50
    if (isNonCartesian) allocate(PlotXYZNodes(nPlotDim, &
         nNodeCellsPerBlock(1),nNodeCellsPerBlock(2),nNodeCellsPerBlock(3),nBlocksUsed))

    call test_stop(NameSub, DoTest)

  end subroutine init_hdf5_plot
  !============================================================================

  subroutine write_var_hdf5(iFile, plotType, iBlock, H5Index,nPlotVar,PlotVar, &
       xMin, xMax, yMin, yMax, zMin, zMax, DxBlock, DyBlock, DzBlock,&
       isNonCartesian, NotACut, nCell,H5Advance)

    ! Save all cells within plotting range, for each processor

    use ModMain, ONLY: nI, nJ, nK, &
         x_, y_, z_, Phi_, nBlockMax, Unused_B
    use ModGeometry, ONLY: CellSize_DB,&
         XyzStart_BLK
    use ModIO
    use ModMpi
    use BATL_lib, ONLY : Xyz_DNB, IsRLonLat, IsCylindrical, CoordMin_DB,&
         CoordMax_DB, CoordMin_D,CoordMax_D,DiLevelNei_IIIB

    ! Arguments

    integer, intent(in)   :: iFile, iBlock, H5Index
    integer, intent(in)   :: nPlotVar
    real,    intent(in)   :: PlotVar(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nPlotVar)
    real,    intent(in)   :: xMin,xMax,yMin,yMax,zMin,zMax
    logical, intent(in) :: isNonCartesian, NotACut
    character (len=*), intent(in) :: plotType
    real,    intent(inout):: DxBlock,DyBlock,DzBlock
    integer, intent(out)  :: nCell
    logical, intent(out) :: H5advance

    ! Local variables
    ! Indices and Coordinates
    integer :: i, j, k, iMin, iMax, jMin, jMax, kMin, kMax
    integer :: iMin1, iMax1, jMin1, jMax1, kMin1, kMax1
    integer ::  iMaxN,jMaxN, kMaxN,iMaxN1,jMaxN1, kMaxN1
    integer :: ii, jj, kk, iSize, jSize, kSize
    integer :: Error,i1, j1, k1,ii1, jj1, kk1
    real ::Dx
    real :: xMin1,xMax1,yMin1,yMax1,zMin1,zMax1

    real :: ySqueezed

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'write_var_hdf5'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    if(MPI_COMM_NULL == iplotComm ) RETURN

    xMin1 = xMin - cHalfMinusTiny*CellSize_DB(x_,iBlock)
    xMax1 = xMax + cHalfMinusTiny*CellSize_DB(x_,iBlock)
    yMin1 = yMin - cHalfMinusTiny*CellSize_DB(y_,iBlock)
    yMax1 = yMax + cHalfMinusTiny*CellSize_DB(y_,iBlock)
    zMin1 = zMin - cHalfMinusTiny*CellSize_DB(z_,iBlock)
    zMax1 = zMax + cHalfMinusTiny*CellSize_DB(z_,iBlock)

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
    !   if(DoTest)then
    !      write(*,*) NameSub, 'xMin1,xMax1,yMin1,yMax1,zMin1,zMax1=',&
    !           xMin1,xMax1,yMin1,yMax1,zMin1,zMax1
    !      write(*,*) NameSub, 'xyzStart_BLK=',iBlock,xyzStart_BLK(:,iBlock)
    !      write(*,*) NameSub, 'ySqueezed =',ySqueezed
    !      write(*,*) NameSub, 'xyzEnd=',xyzStart_BLK(x_,iBlock)+(nI-1)*CellSize_DB(x_,iBlock),&
    !           ySqueezed+(nJ-1)*CellSize_DB(y_,iBlock),&
    !           xyzStart_BLK(z_,iBlock)+(nK-1)*CellSize_DB(z_,iBlock)
    !   end if
    !
    ! If block is fully outside of cut then cycle
    if(  xyzStart_BLK(x_,iBlock) > xMax1.or.&
         xyzStart_BLK(x_,iBlock)+(nI-1)*CellSize_DB(x_,iBlock) < xMin1.or.&
         ySqueezed > yMax1.or.&
         ySqueezed+(nJ-1)*CellSize_DB(y_,iBlock) < yMin1.or.&
         xyzStart_BLK(z_,iBlock) > zMax1.or.&
         xyzStart_BLK(z_,iBlock)+(nK-1)*CellSize_DB(z_,iBlock) < zMin1) then
       H5Advance = .false.
       RETURN
    end if

    Dx = plot_Dx(1,iFile)
    DxBlock=CellSize_DB(x_,iBlock); DyBlock=CellSize_DB(y_,iBlock); DzBlock=CellSize_DB(z_,iBlock)

    ! Calculate index limits of cells inside cut

    iMin = max(1 ,floor((xMin1-xyzStart_BLK(x_,iBlock))/DxBlock)+2)
    iMax = min(nI,floor((xMax1-xyzStart_BLK(x_,iBlock))/DxBlock)+1)

    jMin = max(1 ,floor((yMin1-ySqueezed)/DyBlock)+2)
    jMax = min(nJ,floor((yMax1-ySqueezed)/DyBlock)+1)

    kMin = max(1,floor((zMin1-xyzStart_BLK(z_,iBlock))/DzBlock)+2)
    kMax = min(nK,floor((zMax1-xyzStart_BLK(z_,iBlock))/DzBlock)+1)

    iMin1 = max(-1 ,floor((xMin1-xyzStart_BLK(x_,iBlock))/DxBlock)+2)
    iMax1 = min(nI+1,floor((xMax1-xyzStart_BLK(x_,iBlock))/DxBlock)+1)

    jMin1 = max(-1 ,floor((yMin1-ySqueezed)/DyBlock)+2)
    jMax1 = min(nJ+1,floor((yMax1-ySqueezed)/DyBlock)+1)

    kMin1 = max(-1,floor((zMin1-xyzStart_BLK(z_,iBlock))/DzBlock)+2)
    kMax1 = min(nK+1,floor((zMax1-xyzStart_BLK(z_,iBlock))/DzBlock)+1)

    !    if (iMin1==0 .and.  iMin==1 .and. iMax==1) then ! Cut falls between iMinCell and iMinCell - 1
    !         if (DiLevelNei_IIIB(-1, 0, 0,iBlock) == -1) then
    !             h5Advance = .false.
    !             return
    !         end if
    !         ! Neighbor is more refined so it writes
    !     end if
    !
    !     if (jMin1==0 .and. jMin==1 .and. jMax == 1) then
    !         ! Neighbor is more refined so it writes
    !         if (DiLevelNei_IIIB(0, -1, 0,iBlock) == -1) then
    !             h5Advance = .false.
    !             return
    !         end if
    !     end if
    !     if (kMin1==0 .and. kMin==1 .and. kMax == 1) then
    !         ! Neighbor is more refined so it writes
    !         if (DiLevelNei_IIIB(0, 0, -1,iBlock) == -1) then
    !             h5Advance = .false.
    !             return
    !         end if
    !     end if
    !     if (iMax1==nI+1 .and.  iMax==nI .and. iMin == nI) then ! Cut falls between iMaxCell and iMaxCell 1
    !         if (DiLevelNei_IIIB(1, 0, 0,iBlock) == -1) then
    !             h5Advance = .false.
    !             return
    !         else if (CoordMax_DB(1, iBlock) /= CoordMax_D(1)) then
    !             h5Advance = .false.
    !             return
    !         end if
    !     end if
    !
    !     if (jMax1==nJ+1 .and. jMax==nJ .and. jMin == nJ) then
    !         if (DiLevelNei_IIIB(0, 1, 0,iBlock) == -1) then
    !             h5Advance = .false.
    !             return
    !         else if (CoordMax_DB(2, iBlock) /= CoordMax_D(2)) then
    !             h5Advance = .false.
    !             return
    !         end if
    !     end if
    !
    !     if (kMax1==nK+1 .and. kMax==nK .and. kMin==nK) then
    !         if (DiLevelNei_IIIB(0, 0, 1,iBlock) == -1) then
    !             h5Advance = .false.
    !             return
    !         else if (CoordMax_DB(3, iBlock) /= CoordMax_D(3)) then
    !             h5Advance = .false.
    !             return
    !         end if
    !     end if

    iMaxN = iMax
    jMaxN = jMax
    kMaxN = kMax

    !   if(DoTest)then
    !      write(*,*) NameSub, 'iMin,iMax,jMin,jMax,kMin,kMax=',&
    !           iMin,iMax,jMin,jMax,kMin,kMax
    !      write(*,*) NameSub, 'DxBlock,x1,y1,z1',DxBlock,xyzStart_BLK(:,iBlock)
    !      write(*,*) NameSub, 'ySqueezed  =',ySqueezed
    !      write(*,*) NameSub, 'xMin1,xMax1=',xMin1,xMax1
    !      write(*,*) NameSub, 'yMin1,yMax1=',yMin1,yMax1
    !      write(*,*) NameSub, 'zMin1,zMax1=',zMin1,zMax1
    !   end if
    !

    if (isNonCartesian) then
       if ( plotType=='x=0' .or. plotType=='y=0') then
          iMaxN = iMax + 1
          kMaxN = kMax + 1

          iMaxN1 = iMax1 + 1
          kMaxN1 = kMax1 + 1

       else if ( plotType=='z=0') then
          iMaxN = iMax + 1
          jMaxN = jmax + 1
          iMaxN1 = iMax1 + 1
          jMaxN1 = jmax1 + 1
       else if (NotACut) then
          do ii = 1, nPlotDim
             if (iPlotDim(ii) == 1) then
                iMaxN = iMax + 1
             else if (iPlotDim(ii) == 2) then
                jMaxN = jMax + 1
             else if (iPlotDim(ii) == 3) then
                kMaxN = kMax + 1
             end if
          end do
       end if
    else
       do ii = 1, nPlotDim
          if (iPlotDim(ii) == 1) then
             iMaxN = iMax + 1
          else if (iPlotDim(ii) == 2) then
             jMaxN = jMax + 1
          else if (iPlotDim(ii) == 3) then
             kMaxN = kMax + 1
          end if
       end do
    end if

    if (isNonCartesian) then
       do kk=kMin,kMaxN; do jj=jMin,jMaxN; do ii=iMin,iMaxN
          if (nPlotDim == 3) then
             i = ii - iMin + 1
             j = jj - jMin + 1
             k = kk - kMin + 1

             PlotXYZNodes(1,i, j, k, H5Index) = &
                  Xyz_DNB(1,ii, jj, kk, iBlock)
             PlotXYZNodes(2,i, j, k, H5Index) = &
                  Xyz_DNB(2,ii, jj, kk, iBlock)
             PlotXYZNodes(3,i, j, k, H5Index) = &
                  Xyz_DNB(3,ii, jj, kk, iBlock)

          else if (plotType=='x=0') then
             i = ii - iMin + 1
             j = kk - kMin + 1
             k = jj - jMin + 1

             PlotXYZNodes(1,i, j, k, H5Index) = &
                  Xyz_DNB(2,ii, jj, kk, iBlock)

             PlotXYZNodes(2,i, j, k, H5Index) = &
                  Xyz_DNB(3,ii, jj, kk, iBlock)
          else if (plotType=='y=0') then
             i = ii - iMin + 1
             j = kk - kMin + 1
             k = jj - jMin + 1

             PlotXYZNodes(1,i, j, k, H5Index) = &
                  xyz_dnb(1,ii, jj, kk, iblock)
             PlotXYZNodes(2,i, j, k, H5Index) = &
                  Xyz_DNB(3,ii, jj, kk, iBlock)
             !
          else if (plotType=='z=0') then
             i = ii - iMin + 1
             j = jj - jMin + 1
             k = kk - kMin + 1

             PlotXYZNodes(1,i, j, k, H5Index) = &
                  Xyz_DNB(1,ii, jj, kk, iBlock)
             PlotXYZNodes(2,i, j, k, H5Index) = &
                  Xyz_DNB(2,ii, jj, kk, iBlock)

          else
             i = ii - iMin + 1
             j = jj - jMin + 1
             k = kk - kMin + 1

             PlotXYZNodes(1,i, j, k, H5Index) = &
                  Xyz_DNB(1,ii, jj, kk, iBlock)
             PlotXYZNodes(2,i, j, k, H5Index) = &
                  Xyz_DNB(2,ii, jj, kk, iBlock)
             PlotXYZNodes(3,i, j, k, H5Index) = &
                  Xyz_DNB(3,ii, jj, kk, iBlock)

          end if
          where(abs(PlotXYZNodes(:,:,:,:,H5Index))<1.e-10) &
               PlotXYZNodes(:,:,:,:,H5Index) = 0.
       end do; end do; end do
       do kk=kMin,kMax; do jj=jMin,jMax; do ii=iMin,iMax
          if (nPlotDim == 3 .or. NotACut) then
             i = ii - iMin + 1
             j = jj - jMin + 1
             k = kk - kMin + 1
             PlotVarIdx(i, j, k, H5Index, 1:nPlotVar) = &
                  PlotVar(ii, jj, kk, 1:nPlotVar)
          else if ( plotType=='x=0' .or. plotType=='y=0') then
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
          end if
          nCell = nCell+1
       end do; end do; end do
    else
       do kk=kMin,kMax; do jj=jMin,jMax; do ii=iMin,iMax

          if (nPlotDim < 3) then
             if (PlotType == 'z=0') then
                i = ii - iMin + 1
                j = jj - jMin + 1
                k = kk - kMin + 1
                PlotVarIdx(i, j, k, H5Index, 1:nPlotVar) = &
                     .5*(plotvar(ii, jj, kmin1, 1:nplotvar)+plotvar(ii, jj, kmax1, 1:nplotvar))
             else if (PlotType == 'y=0') then
                i = ii - iMin + 1
                j = kk - kMin + 1
                k = jj - jMin + 1
                PlotVarIdx(i, j, k, H5Index, 1:nPlotVar) = &
                     .5*(plotvar(ii, jMin1 ,kk, 1:nplotvar)+plotvar(ii, jMax1, kk, 1:nplotvar))
             else if (PlotType == 'x=0') then
                i = jj - jMin + 1
                j = kk - kMin + 1
                k = ii - iMin + 1
                PlotVarIdx(i, j, k, H5Index, 1:nPlotVar) = &
                     .5*(plotvar(iMin1, jj ,kk, 1:nplotvar)+plotvar(iMax1, jj, kk, 1:nplotvar))
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
    UsedBlocks(h5Index) = iBlock
    call test_stop(NameSub, DoTest, iBlock)

  end subroutine write_var_hdf5
  !============================================================================

  subroutine write_plot_hdf5(filename, plotType, plotVarNames, plotVarUnits,&
       nPlotVar, NotACut, nonCartesian, IsSphPlot, plot_dimensional, &
       xmin, xmax, ymin, ymax, zmin, zmax)

    use ModNumConst
    use ModMpi
    use BATL_lib,  ONLY: CoordMin_DB, CoordMax_DB, MaxDim, &
         iNode_B, iTree_IA, Coord0_, Level_, iMortonNode_A
    use BATL_mpi,  ONLY: barrier_mpi

    integer,                 intent(in):: nPlotVar
    character(len=80),       intent(in):: filename
    character(len=lnamevar), intent(in):: plotvarnames(nplotvar)
    character(len=lNameVar),  intent(in):: plotVarUnits(nPlotVar)
    character (len=*), intent(in) :: plotType
    real, intent(in)  ::  xMin, xMax, yMin, yMax, zMin, zMax
    logical, intent(in) :: NotACut, plot_dimensional, nonCartesian, IsSphPlot
    integer :: offsetPerProc(0:nPlotProc-1)
    integer :: Error, procIdx, iLen, iVar, nBlocksUsedMax

    integer :: labelLeng, i, ii, nBlocksUsedMin
    integer(HID_T) :: iIntHIDT

    real, allocatable :: Coordinates(:,:), BoundingBox(:,:,:),BlockDataExtents(:,:)
    integer, allocatable :: procNum(:)
    real :: VarMin, VarMax
    real :: GlobalVarMin(nPlotVar), GlobalVarMax(nPlotVar)
    real :: minCoords(3), maxCoords(3)
    character (len=lnamevar+4) :: DatasetNameTemp

    logical :: isCutFile

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'write_plot_hdf5'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    if(MPI_COMM_NULL == iplotComm ) RETURN

    if (NotACut) then
       isCutFile = .false.
    else
       isCutFile = .true.
    end if

    if (nBlocksUsed > 0) then
       allocate(UsedNodes(nBlocksUsed))
       UsedNodes(1:nBlocksUsed) = iNode_B(UsedBlocks(1:nBlocksUsed))
    else
       allocate(UsedNodes(1))
       UsedNodes(1) = 0
    end if
    do iVar = 1, nPlotVar
       if (nBlocksUsed > 0) then
          VarMin = minval(PlotVarIdx(:,:,:,:,iVar))
          VarMax = maxval(PlotVarIdx(:,:,:,:,iVar))
       else
          VarMin = huge(VarMin)
          VarMax = -huge(VarMax)
       end if
       call MPI_allreduce(VarMin, GlobalVarMin(iVar), 1, MPI_REAL,&
            MPI_MIN, iplotComm, Error)

       call MPI_allreduce(VarMax, GlobalVarMax(iVar), 1, MPI_REAL,&
            MPI_MAX, iplotComm, Error)

    end do

    ! Open the hdf5 file
    call open_hdf5_file(fileid, filename, iplotComm)
    if(fileid == -1) then
       if (iPlotProc == 0) write (*,*)  "Error: unable to initialize file"
       call stop_hdf5("unable to initialize hdf5 file")
    end if

    allocate(unknownNameArray(nPlotVar))

    call pad_string_with_null(nPlotVar, int(lNameH5), PlotVarNames, UnknownNameArray)
    iInteger8=nPlotVar
    iInteger8a=lNameH5
    call write_hdf5_data(FileID, "plotVarNames", 1,  (/iInteger8/),&
         CharacterData=UnknownNameArray, nStringChars=iInteger8a)

    ! The plugn doesn't need these. However it does use this info to speed some
    ! things up if it is there.  Perhaps this could be an option in PARAM.in
    allocate (BlockDataExtents(2,nBlocksUsed))
    do iVar = 1, nPlotVar
       do iBlock = 1, nBlocksUsed
          BlockDataExtents(1, iBlock) = minval(PlotVarIdx(:,:,:,iBlock,iVar))
          BlockDataExtents(2, iBlock) = maxval(PlotVarIdx(:,:,:,iBlock,iVar))
       end do
       DatasetNameTemp = trim(plotVarNames(iVar))//"_Ext"
       labelLeng = len_trim(UnknownNameArray(iVar)) + 4
       do iLen = labelLeng + 1,lNameH5+3
          DatasetNameTemp(iLen:iLen) = CHAR(0)
       end do
       if (nBlocksUsed == 0 ) then
          iInteger8 = 2
          call write_hdf5_data(FileID, DatasetNameTemp, 2, (/iInteger8,nBlkUsedGlobal/),&
               nOffsetLocal=iProcOffset, nBlocksLocalIn=nBlocksUsed,&
               Rank2RealData=reshape((/0.0,0.0,0.0,0.0/),(/2,2/)))
       else
          iInteger8 = 2
          call write_hdf5_data(FileID, DatasetNameTemp, 2, (/iInteger8,nBlkUsedGlobal/),&
               nOffsetLocal=iProcOffset, nBlocksLocalIn=nBlocksUsed,&
               Rank2RealData=BlockDataExtents)
       end if

    end do

    deallocate (BlockDataExtents)

    do iVar = 1, nPlotVar
       call write_hdf5_data(FileID, UnknownNameArray(iVar), 4, (/nCellsPerBlock(1), nCellsPerBlock(2),&
            nCellsPerBlock(3),nBlkUsedGlobal/),nOffsetLocal=iProcOffset, nBlocksLocalIn=nBlocksUsed,&
            Rank4RealData=PlotVarIdx(:,:,:,:,iVar),&
            RealAttribute1=GlobalVarMin(iVar), RealAttribute2=GlobalVarMax(iVar), &
            NameRealAttribute1="minimum", NameRealAttribute2="maximum")
    end do

    allocate(BoundingBox(2,nPlotDim,nBlocksUsed))
    if (nonCartesian) then
       call write_hdf5_data(FileID, "NodesX", 4, (/nNodeCellsPerBlock(1), nNodeCellsPerBlock(2),&
            nNodeCellsPerBlock(3),nBlkUsedGlobal/),nOffsetLocal=iProcOffset, nBlocksLocalIn=nBlocksUsed,&
            Rank4RealData=PlotXYZNodes(1,:,:,:,:))
       if (nPlotDim >= 2) then
          call write_hdf5_data(FileID, "NodesY", 4, (/nNodeCellsPerBlock(1), nNodeCellsPerBlock(2),&
               nNodeCellsPerBlock(3),nBlkUsedGlobal/),nOffsetLocal=iProcOffset, nBlocksLocalIn=nBlocksUsed,&
               Rank4RealData=PlotXYZNodes(2,:,:,:,:))
          if (nPlotDim == 3) &
               call write_hdf5_data(FileID, "NodesZ", 4, (/nNodeCellsPerBlock(1), nNodeCellsPerBlock(2),&
               nNodeCellsPerBlock(3),nBlkUsedGlobal/),nOffsetLocal=iProcOffset, nBlocksLocalIn=nBlocksUsed,&
               Rank4RealData=PlotXYZNodes(3,:,:,:,:))
       end if
       do iBlock = 1, nBlocksUsed
          BoundingBox(1,1,iBlock) = minval(PlotXYZNodes(1,:,:,:,iBlock))
          BoundingBox(2,1,iBlock) = maxval(PlotXYZNodes(1,:,:,:,iBlock))
          if (nPlotDim >= 2) then
             BoundingBox(1,2,iBlock) = minval(PlotXYZNodes(2,:,:,:,iBlock))
             BoundingBox(2,2,iBlock) = maxval(PlotXYZNodes(2,:,:,:,iBlock))
          end if
          if (nPlotDim == 3) then
             BoundingBox(1,3,iBlock) = minval(PlotXYZNodes(3,:,:,:,iBlock))
             BoundingBox(2,3,iBlock) = maxval(PlotXYZNodes(3,:,:,:,iBlock))
          end if
       end do
       deallocate(PlotXYZNodes)

    else
       BoundingBox(1,:,:) = CoordMin_DB(iPlotDim(1:nPlotDim), UsedBlocks)
       BoundingBox(2,:,:) = CoordMax_DB(iPlotDim(1:nPlotDim), UsedBlocks)
    end if
    iInteger8=2
    call write_hdf5_data(FileID, "bounding box", 3, (/iInteger8,nPlotDim,nBlkUsedGlobal/),&
         nOffsetLocal=iProcOffset, nBlocksLocalIn=nBlocksUsed,&
         Rank3RealData=BoundingBox)
    allocate(Coordinates(nPlotDim, nBlocksUsed))
    do iBlock = 1, nBlocksUsed
       do i = 1, nPlotDim
          Coordinates(i, iBlock) = .5*(BoundingBox(1,i,iBlock) + BoundingBox(2,i,iBlock))
       end do
    end do
    deallocate(BoundingBox)
    call MPI_Barrier(iPlotComm,Error)
    call write_hdf5_data(FileID, "coordinates", 2, (/nPlotDim,nBlkUsedGlobal/),&
         nOffsetLocal=iProcOffset, nBlocksLocalIn=nBlocksUsed,&
         Rank2RealData=Coordinates)
    deallocate(Coordinates)

    call MPI_Barrier(iPlotComm,Error)
    if (nBlocksUsed > 0) then
       allocate(MinLogicalExtents(nPlotDim, nBlocksUsed))
       do iVar = 1,nPlotDim
          MinLogicalExtents(:,:) = iTree_IA(iPlotDim(1:nPlotDim)+Coord0_,UsedNodes)*nCellsPerBlock(iVar)
       end do
    else
       allocate(MinLogicalExtents(nPlotDim, 1))
       MinLogicalExtents = huge(MinLogicalExtents(1,1))
    end if

    do iVar = 1,nPlotDim
       i = minval(MinLogicalExtents(iVar,:))
       call MPI_allreduce(i, ii, 1, MPI_INTEGER,&
            MPI_MIN, iplotComm, Error)
       do iBlock = 1, nBlocksUsed
          MinLogicalExtents(iVar,iBlock) = MinLogicalExtents(iVar,iBlock) - ii
       end do
    end do

    call  write_hdf5_data(FileID, "MinLogicalExtents", 2, (/nPlotDim,nBlkUsedGlobal/),&
         nOffsetLocal=iProcOffset, nBlocksLocalIn=nBlocksUsed,&
         Rank2IntegerData=MinLogicalExtents)
    deallocate(MinLogicalExtents)

    call pad_string_with_null(nPlotVar, int(lNameH5), PlotVarUnits, UnknownNameArray)

    iInteger8=nPlotVar
    iInteger8a=lNameH5
    call write_hdf5_data(FileID, "plotVarUnits", 1,  (/iInteger8/),&
         CharacterData=UnknownNameArray, nStringChars=iInteger8a)

    call MPI_Barrier(iPlotComm,Error)
    deallocate(PlotVarIdx)
    deallocate(unknownNameArray)

    call MPI_Barrier(iPlotComm,Error)
    if (NotACut .and. nBlocksUsed > 0) then
       call  write_hdf5_data(FileID, "iMortonNode_A", 1, (/nBlkUsedGlobal/),&
            nOffsetLocal=iProcOffset, nBlocksLocalIn=nBlocksUsed,&
            Rank1IntegerData=iMortonNode_A(UsedNodes))
    else if (NotACut) then
       call  write_hdf5_data(FileID, "iMortonNode_A", 1, (/nBlkUsedGlobal/),&
            nOffsetLocal=iProcOffset, nBlocksLocalIn=nBlocksUsed,&
            Rank1IntegerData=(/0/))
    end if
    !
    call MPI_Barrier(iPlotComm,Error)
    if (nBlocksUsed > 0 ) then
       call  write_hdf5_data(FileID, "refine level", 1, (/nBlkUsedGlobal/),&
            nOffsetLocal=iProcOffset, nBlocksLocalIn=nBlocksUsed,&
            Rank1IntegerData=iTree_IA(Level_,UsedNodes))
       deallocate(usedNodes)
    else
       call  write_hdf5_data(FileID, "refine level", 1, (/nBlkUsedGlobal/),&
            nOffsetLocal=iProcOffset, nBlocksLocalIn=nBlocksUsed,&
            Rank1IntegerData=(/0/))
       deallocate(usedNodes)
    end if

    call MPI_Barrier(iPlotComm,Error)
    allocate(procNum(nBlocksUsed))
    procNum = iPlotProc
    call  write_hdf5_data(FileID, "Processor Number", 1, (/nBlkUsedGlobal/),&
         nOffsetLocal=iProcOffset, nBlocksLocalIn=nBlocksUsed,&
         Rank1IntegerData=procNum)
    deallocate(procNum)
    deallocate(usedBlocks)
    allocate(UnknownNameArray(nPlotDim))

    call MPI_Barrier(iPlotComm,Error)
    if (plotType=='x=0') then
       UnknownNameArray(1) = "Y-Axis"
       if (nPlotDim == 2) &
            UnknownNameArray(2) = "Z-Axis"
    else
       do iVar = 1, nPlotDim

          if (iPlotDim(iVar) == 1) then
             UnknownNameArray(iVar) = "X-Axis"
          else if (iPlotDim(iVar) == 2) then
             UnknownNameArray(iVar) = "Y-Axis"
          else if (iPlotDim(iVar) == 3) then
             UnknownNameArray(iVar) = "Z-Axis"
          end if
       end do
    end if

    call pad_string_with_null(int(nPlotDim), int(lNameH5), UnknownNameArray, UnknownNameArray)
    iInteger8 = lNameH5
    call write_hdf5_data(FileID, "Axis Labels", 1,  (/nPlotDim/),&
         CharacterData=UnknownNameArray, nStringChars=iInteger8)
    
    deallocate(UnknownNameArray)
    iIntHIDT = nPlotVar ! Convert to correct data type.
    call write_integer_sim_metadata(fileID, iIntHIDT)
    call write_real_sim_metadata(FileID,plot_dimensional)
    call write_integer_plot_metadata(fileID, iIntHIDT, isCutFile)
    if (plotType=='x=0') then
       call write_real_plot_metadata(fileid,plot_dimensional, .true.)
    else
       call write_real_plot_metadata(FileID,plot_dimensional, .false.)
    end if

    call MPI_Barrier(iPlotComm,Error)
    call close_hdf5_file(FileID)

    call MPI_Comm_Free(iPlotComm,Error)

    call test_stop(NameSub, DoTest)

  end subroutine write_plot_hdf5
  !============================================================================

  subroutine write_real_plot_metadata(FileID,plot_dimensional, isXZero)

    use ModMain, only : Time_Simulation, CodeVersion
    ! use ModProcMH, only : iProc
    use BATL_lib, only : nLevelMax, nDim
    use ModGeometry, only : x1,x2,y1,y2,z1,z2
    use ModPhysics, ONLY : No2Io_V, UnitX_
    integer(HID_T), intent(in) :: FileID
    integer(HSIZE_T) :: iData
    integer :: i
    logical, intent (in) :: plot_dimensional, isXZero

    real :: RealMetaData(7)

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'write_real_plot_metadata'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    if(MPI_COMM_NULL == iplotComm ) RETURN

    iData = 1
    !    attName(1) = 'Simulation Time'
    RealMetaData(iData) = Time_Simulation
    iData = iData + 1
    if (isXZero) then
       RealMetaData(iData) = y1
       iData = iData + 1
       RealMetaData(iData) = y2
       iData = iData + 1
       if (nPlotDim == 1) then
          RealMetaData(iData) = 0.0
          iData = iData + 1
          RealMetaData(iData) = 1.0
       else
          RealMetaData(iData) = z1
          iData = iData + 1
          RealMetaData(iData) = z2
       end if
       iData = iData + 1
       RealMetaData(iData) = 0.0
       iData = iData + 1
       RealMetaData(iData) = 0.0
       iData = iData + 1
    else
       do i=1,3
          if (iPlotDim(i) == 1) then
             RealMetaData(iData) = x1
             iData = iData + 1
             RealMetaData(iData) = x2
          else if(iPlotDim(i) == 2) then
             RealMetaData(iData) = y1
             iData = iData + 1
             RealMetaData(iData) = y2
          else if(iPlotDim(i) == 3) then
             RealMetaData(iData) = z1
             iData = iData + 1
             RealMetaData(iData) = z2
          else if (i == 2 .and. nPlotDim == 1) then
             RealMetaData(iData) = 0.0
             iData = iData + 1
             RealMetaData(iData) = 1.0
          else if (iPlotDim(i) == 0) then
             RealMetaData(iData) = 0.0
             iData = iData + 1
             RealMetaData(iData) = 0.0
          end if
          iData = iData + 1
       end do
    end if
    iData = iData -1

    !-------------------------------------------------------------------
    ! write the real Metadata
    call  write_hdf5_data(FileID, "Real Plot Metadata", 1, (/iData/),&
         Rank1RealData=RealMetaData)

    call test_stop(NameSub, DoTest)
  end subroutine write_real_plot_metadata
  !============================================================================

  subroutine write_real_sim_metadata(FileID,plot_dimensional)

    use ModMain, only : Time_Simulation, CodeVersion
    use BATL_lib, only : nLevelMax, nDim
    use ModGeometry, only : x1,x2,y1,y2,z1,z2
    use ModPhysics, ONLY : No2Io_V, UnitX_

    integer (HID_T), intent(in) :: fileID
    logical, intent (in) :: plot_dimensional
    integer :: rank
    integer(HSIZE_T) :: iData

    real :: RealMetaData(7)
    !    allocate(attName(nAtts))

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'write_real_sim_metadata'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    if(MPI_COMM_NULL == iplotComm ) RETURN

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
       !         write (*,*) "x1,x2,y1,y2,z1,z2",x1,x2,y1,y2,z1,z2
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
    ! write the real Metadata
    call  write_hdf5_data(FileID, "Real Simulation Metadata", 1, (/iData/),&
         Rank1RealData=RealMetaData, RealAttribute1=CodeVersion,&
         NameRealAttribute1="CodeVersion")

    call test_stop(NameSub, DoTest)
  end subroutine write_real_sim_metadata
  !============================================================================

  subroutine write_integer_plot_metadata(fileID,nPlotVar,isCutFile)

    use BATL_lib, only : nDimAmr, nLevelMax, IsPeriodic_D
    use ModMain, only : n_step
    use ModGeometry, ONLY : TypeGeometry
    integer (HID_T), intent(in) :: fileID, nPlotVar
    logical, intent(in) :: isCutFile

    integer(HSIZE_T) ::  iData
    integer :: i,IntegerMetaData(16)
    integer, parameter  :: FFV = 1 ! file format version

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'write_integer_plot_metadata'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    if(MPI_COMM_NULL == iplotComm ) RETURN

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
    IntegerMetaData(iData) = nPlotProc
    !    attName(5) = 'numProcessors'
    iData = iData + 1

    IntegerMetaData(iData) = nLevelMax
    !    attName(7) = 'nLevelMax'
    iData = iData + 1
    do i = 1, 3
       if (iPlotDim(i) == 1) then
          IntegerMetaData(iData) = iSizeGlobal
       else if (iPlotDim(i) == 2) then
          IntegerMetaData(iData) = jSizeGlobal
       else if (iPlotDim(i) == 3) then
          IntegerMetaData(iData) = kSizeGlobal
       else if (iPlotDim(i) == 0) then
          IntegerMetaData(iData) = 0
       end if
       iData = iData + 1
    end do

    if(TypeGeometry=='cartesian') then
       IntegerMetaData(iData) = CartesianPlot_
    elseif(TypeGeometry=='rz') then
       IntegerMetaData(iData) = RzPlot_
    elseif(TypeGeometry=='roundcube') then
       IntegerMetaData(iData) = RoundCubePlot_
    elseif(TypeGeometry=='cylindrical') then
       IntegerMetaData(iData) = CylindricalPlot_
    elseif(TypeGeometry=='cylindrical_lnr') then
       IntegerMetaData(iData) = CylindricalLnrPlot_
    elseif(TypeGeometry=='cylindrical_genr') then
       IntegerMetaData(iData) = CylindricalGenrPlot_
    elseif(TypeGeometry=='spherical') then
       IntegerMetaData(iData) = SphericalPlot_
    elseif(TypeGeometry=='spherical_lnr') then
       IntegerMetaData(iData) = SphericalLnrPlot_
    elseif(TypeGeometry=='spherical_genr') then
       IntegerMetaData(iData) = SphericalGenrPlot_
    elseif(TypeGeometry=='rlonlat') then
       IntegerMetaData(iData) = rLonLatPlot_
    elseif(TypeGeometry=='rlonlat_lnr') then
       IntegerMetaData(iData) = rLonLatLnrPlot_
    elseif(TypeGeometry=='rlonlat_genr') then
       IntegerMetaData(iData) = rLonLatGenrPlot_
    endif
    iData = iData + 1
    ! as of 2/3/2012 this is not implimented in the plugin but it probably
    ! should be in the future
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
    ! write the integer Metadata
    call  write_hdf5_data(FileID, "Integer Plot Metadata", 1, (/iData/),&
         Rank1IntegerData=IntegerMetaData)

    call test_stop(NameSub, DoTest)

  end subroutine write_integer_plot_metadata
  !============================================================================

  subroutine write_integer_sim_metadata(fileID, nPlotVar)
    ! Not read by plugin at this time  Only exists so that one looking at
    ! the file may know something about the simulation that created it

    use BATL_lib, only : nDim, nDimAmr, nLevelMax
    use ModMain, only : n_step, nI, nJ, nK
    integer (HID_T), intent(in) :: fileID, nPlotVar

    integer(HSIZE_T) :: iData
    integer :: IntegerMetaData(8)
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'write_integer_sim_metadata'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    if(MPI_COMM_NULL == iplotComm ) RETURN

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
    IntegerMetaData(iData) = nPlotProc
    !    attName(5) = 'numProcessors'
    iData = iData + 1

    IntegerMetaData(iData) = nLevelMax
    !    attName(7) = 'nLevelMax'
    iData = iData + 1

    IntegerMetaData(iData) = n_step
    !    attName(2) = "Time Step"
    !    iData = iData + 1

    !-------------------------------------------------------------------
    ! write the integer Metadata
    call  write_hdf5_data(FileID, "Integer Sim Metadata", 1, (/iData/),&
         Rank1IntegerData=IntegerMetaData)
    call test_stop(NameSub, DoTest)

  end subroutine write_integer_sim_metadata
  !============================================================================

  real function minmod(x,y)
    real, intent(in) :: x,y
    !--------------------------------------------------------------------------
    minmod = max(0.0, min(abs(x), sign(1.0,x)*y))
  end function minmod
  !============================================================================

end module ModHdf5
!==============================================================================
