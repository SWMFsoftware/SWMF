!^CFG COPYRIGHT UM
!------------------------------------------------------------------------
! Front end routine to collect refined blocks from each processor into
!   global list before calling parallel_refine
!------------------------------------------------------------------------
subroutine refine_grid(refBLK)
  use ModProcMH
  use ModMain, ONLY : nBLK
  use ModAMR, ONLY : refine_list
  use ModMpi
  implicit none

  logical, intent(in) :: refBLK(nBLK)

  integer :: iError

  call MPI_ALLGATHER(refBLK, nBLK, MPI_LOGICAL, &
       refine_list, nBLK, MPI_LOGICAL, iComm, iError)

  call parallel_refine

end subroutine refine_grid

!------------------------------------------------------------------------
! parallel_refine assumes that "refine_list" is filled
!   Execute refinement in parallel
!------------------------------------------------------------------------
subroutine parallel_refine
  use ModProcMH
  use ModMain, ONLY : iNewGrid,nBlock,nBlockMax,unusedBLK,nBLK,lVerbose, &
       BlkTest
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK,R_BLK,dx_BLK,dy_BLK,dz_BLK,&
       XyzStart_BLK
  use ModAMR, ONLY : local_cube,local_cubeBLK,availableBLKs,refine_list
  use ModAdvance, ONLY : iTypeAdvance_B, SkippedBlock_
  use ModIO, ONLY: write_prefix, iUnitOut
  use ModOctree
  implicit none

  integer :: iPE, iBLK, icube, nPEsRefBlk, PEsRefBlk(9), minPE, loc(1)
  integer :: ii, nlevel

  logical :: UseOldBlock
  logical :: oktest, oktest_me

  type (adaptive_block_ptr) :: refine_block_ptr

  character(len=*), parameter:: NameSub = 'parallel_refine'
  !---------------------------------------------------------------------------
  UseOldBlock = .false.

  call set_oktest(NameSub,oktest,oktest_me)

  if(oktest_me) write(*,*)NameSub,' refine_list(Test), dx=',&
       refine_list(BlkTest, iProc+1), dx_Blk(BlkTest)

  if (.not. any(refine_list)) RETURN
  if(iProc==0 .and. lVerbose>0)then
     call write_prefix; write (iUnitOut,*) 'parallel_refine', &
          ' starting to refine ',count(refine_list),' block(s) on PE 0'
  end if
  ! Change grid identification number
  iNewGrid = mod( iNewGrid+1, 10000)

  nullify(refine_block_ptr % ptr)

  do iPE=0,nProc-1; do iBLK=1,nBLK
     ! check each solution block for refinement
     if (.not. refine_list(iBLK,iPE+1)) CYCLE

     refine_block_ptr % ptr => global_block_ptrs(iBLK, iPE+1) % ptr

     ! ensure flagged block is allocated
     if (.not. associated(refine_block_ptr % ptr)) then  
        ! ERROR, block not associated, turn refine flag false
        refine_list(iBLK,iPE+1) = .false.
        CYCLE
     end if

     ! refine only solution blocks in use
     if (.not. refine_block_ptr % ptr % used) then    
        ! ERROR, block not used, turn refine flag false
        refine_list(iBLK,iPE+1) = .false.
        CYCLE
     end if

     if(availableBLKs(0,iPE)+7<=nBLK)then
        ! if staying on iPE is possible, stay.
        minPE=iPE
     else
        ! find PE with minimum blocks
        loc=minloc(availableBLKs(0,0:nProc-1)); minPE=loc(1)-1
     end if

     ! fill local_cube and local_cubeBLK
     if(availableBLKs(0,minPE)+7<=nBLK) then
        do icube=1,8
           if(icube==1 .and. (UseOldBlock .or. iPE==minPE)) then
              ! use current block as first of new 8
              local_cube   (1) = iPE
              local_cubeBLK(1) = iBLK
              CYCLE
           end if
           local_cube   (icube) = minPE
           local_cubeBLK(icube) = availableBLKs(availableBLKs(0,minPE),minPE)
           availableBLKs(0,minPE) = availableBLKs(0,minPE)+1
        end do
        if(.not.UseOldBlock .and. iPE/=minPE) then
           ! put old original block back on available list
           availableBLKs(0,iPE) = availableBLKs(0,iPE)-1
           availableBLKs(availableBLKs(0,iPE),iPE) = iBLK
        end if
     else
        ! split new blocks on processors with the least blocks
        local_cube(1) = iPE
        local_cubeBLK(1) = iBLK
        do icube=2,8
           loc=minloc(availableBLKs(0,0:nProc-1)); minPE=loc(1)-1

           if(availableBLKs(0,minPE)>nBLK) then
              write(*,*) "refine_grid: PE = ",iProc, &
                   " available nBLK=",nBLK
              call stop_mpi("ERROR in refine_grid: "// &
                   "insufficient number of blocks !")
           end if
           local_cube   (icube) = minPE
           local_cubeBLK(icube) = &
                availableBLKs(availableBLKs(0,minPE),minPE)
           availableBLKs(0,minPE) = availableBLKs(0,minPE)+1
        end do
     end if

     nPEsRefBlk = 1
     PEsRefBlk = -1
     PEsRefBlk(1) = iPE
     do icube = 1, 8
        ! update nBlock and nBlockMax

        nBlockMax = max(nBlockMax, local_cubeBLK(icube))
        nBlock    = nBlockMax

        ! find out if new block is on same PE as any other previous new blocks
        do ii = 1,nPEsRefBlk
           if(local_cube(icube) == PEsRefBlk(ii)) EXIT
           if(ii == nPEsRefBlk) then
              nPEsRefBlk = nPEsRefBlk + 1
              PEsRefBlk(nPEsRefBlk) = local_cube(icube)
           end if
        end do

        ! set block ready to use
        if (iProc == local_cube(icube)) &
             unusedBLK(local_cubeBLK(icube)) = .false.
     end do

!     if(oktest .and. iProc == iPE) then
!        write(*,'(a,i4,/,a,i1,/,a,9i4,/,2(a,8i4,/))') &
!             "R--> refine_grid: PE = ",iProc, &
!             "R-     nPEsRefBlk    = ",nPEsRefBlk, &
!             "R-     PEsRefBlk     = ",PEsRefBlk, &
!             "R-     local_cube    = ",local_cube, &
!             "R-     local_cubeBLK = ",local_cubeBLK
!     end if

     call create_refined_soln_blocks(nPEsRefBlk, PEsRefBlk, iPE,iBLK)
     call refine_octree_block(refine_block_ptr, local_cube, &
          local_cubeBLK, iPE,iBLK)

     if (iProc==iPE) then
        do icube=1,8
           if (local_cube(icube)==iPE .and. &
                local_cubeBLK(icube)==iBLK) EXIT
           if (icube==8) then
              dx_BLK(iBLK) = -777777.
              dy_BLK(iBLK) = -777777.
              dz_BLK(iBLK) = -777777.

              xyzStart_BLK(:,iBLK) = -777777.

              x_BLK(:,:,:,iBLK) = -777777.
              y_BLK(:,:,:,iBLK) = -777777.
              z_BLK(:,:,:,iBLK) = -777777.
              R_BLK(:,:,:,iBLK) = -777777.
              unusedBLK(iBLK) = .true.
              iTypeAdvance_B(iBLK) = SkippedBlock_
           end if
        end do
     end if
  end do; end do

  if(oktest_me)write(*,*)NameSub,' final dx=',dx_BLK(BlkTest)

end subroutine parallel_refine
