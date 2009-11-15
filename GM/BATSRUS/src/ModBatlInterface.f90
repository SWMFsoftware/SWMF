module ModBatlInterface

  implicit none

contains
  !===========================================================================
  subroutine set_batsrus_grid

    use BATL_lib, ONLY: nNodeUsed, nBlock, MaxBlock, Unused_B, Unused_BP, &
         iProc, iComm

    use ModAmr, ONLY: UnusedBlock_BP

    use ModMain, ONLY: nBlockAll, nBlockBats => nBlock, nBlockMax, UnusedBlk
    
    use ModGeometry, ONLY: dx_BLK, MinDxValue, MaxDxValue

    use ModAdvance, ONLY: iTypeAdvance_B, iTypeAdvance_BP, &
         SkippedBlock_, ExplBlock_
    use ModMpi

    integer:: iBlock, iError
    real   :: DxMin, DxMax
    !-------------------------------------------------------------------------

    nBlockAll  = nNodeUsed
    nBlockBats = nBlock
    call MPI_ALLREDUCE(nBlock, nBlockMax, 1, MPI_INTEGER, MPI_MAX, &
         iComm, iError)

    UnusedBlk      = Unused_B
    UnusedBlock_BP = Unused_BP

    where(Unused_BP)
       iTypeAdvance_BP = SkippedBlock_
    elsewhere
       iTypeAdvance_BP = ExplBlock_
    end where
    iTypeAdvance_B = iTypeAdvance_BP(:,iProc)

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       call set_batsrus_block(iBlock)
    end do

    DxMin = minval(dx_BLK, MASK=(.not.Unused_B))
    DxMax = maxval(dx_BLK, MASK=(.not.Unused_B))
    call MPI_allreduce(DxMin, minDXvalue,1,MPI_REAL,MPI_MIN,iComm,iError)
    call MPI_allreduce(DxMax, maxDXvalue,1,MPI_REAL,MPI_MAX,iComm,iError)

  end subroutine set_batsrus_grid
  !===========================================================================
  subroutine set_batsrus_block(iBlock)

    use BATL_lib, ONLY: CellSize_DB, CoordMin_DB, DiLevelNei_IIIB, iProc
    use ModGeometry, ONLY: dx_BLK, dy_BLK, dz_BLK, XyzStart_BLK
    use ModParallel, ONLY: BLKneighborLEV,  neiLEV, &
         neiLeast, neiLwest, neiLsouth, neiLnorth, neiLbot, neiLtop, &
         neiBeast, neiBwest, neiBsouth, neiBnorth, neiBbot, neiBtop, &
         neiPeast, neiPwest, neiPsouth, neiPnorth, neiPbot, neiPtop

    integer, intent(in):: iBlock
    !-------------------------------------------------------------------------
    BLKneighborLEV(:,:,:,iBlock) = DiLevelNei_IIIB(:,:,:,iBlock)

    neiLeast(iBlock)  = BLKneighborLEV(-1,0,0,iBlock)
    neiLwest(iBlock)  = BLKneighborLEV(+1,0,0,iBlock)
    neiLsouth(iBlock) = BLKneighborLEV(0,-1,0,iBlock)
    neiLnorth(iBlock) = BLKneighborLEV(0,+1,0,iBlock)
    neiLbot(iBlock)   = BLKneighborLEV(0,0,-1,iBlock)
    neiLtop(iBlock)   = BLKneighborLEV(0,0,+1,iBlock)

    neiLEV(1,iBlock)  = neiLeast(iBlock)
    neiLEV(2,iBlock)  = neiLwest(iBlock)
    neiLEV(3,iBlock)  = neiLsouth(iBlock)
    neiLEV(4,iBlock)  = neiLnorth(iBlock)
    neiLEV(5,iBlock)  = neiLbot(iBlock)
    neiLEV(6,iBlock)  = neiLtop(iBlock)

    ! To avoid problems in ModFaceValue::correct_monotone_restrict
    ! we pretend that the neighbor is the block itself, which is used.
    neiBeast(:,iBlock)  = iBlock
    neiBwest(:,iBlock)  = iBlock
    neiBsouth(:,iBlock) = iBlock
    neiBnorth(:,iBlock) = iBlock
    neiBtop(:,iBlock)   = iBlock
    neiBbot(:,iBlock)   = iBlock

    neiPeast(:,iBlock)  = iProc
    neiPwest(:,iBlock)  = iProc
    neiPsouth(:,iBlock) = iProc
    neiPnorth(:,iBlock) = iProc
    neiPtop(:,iBlock)   = iProc
    neiPbot(:,iBlock)   = iProc

    XyzStart_BLK(:,iBlock) = CoordMin_DB(:,iBlock) + 0.5*CellSize_DB(:,iBlock)

    dx_BLK(iBlock) = CellSize_DB(1,iBlock)
    dy_BLK(iBlock) = CellSize_DB(2,iBlock)
    dz_BLK(iBlock) = CellSize_DB(3,iBlock)

    call fix_block_geometry(iBlock)

  end subroutine set_batsrus_block

end module ModBatlInterface
