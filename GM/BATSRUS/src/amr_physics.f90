!^CFG COPYRIGHT UM
subroutine amr_physics
  use ModProcMH
  use ModMain, ONLY : nBLK,nBlockALL,lVerbose
  use ModAMR
  use ModOctree
  use ModIO, ONLY: write_myname, write_prefix, iUnitOut
  use ModMpi
  implicit none

  integer :: nDesiredMax
  integer :: idir, i1,i2, nDesired, nRefined, nCoarsened, currentBlocks
  integer :: i,j,k,n, iBLK, nSort, nRefineMax, nCoarsenMax, Itmp
  integer :: iError, SortIndex(4)
  real :: percentCoarsenMax=100.0
  real :: refine_criteria(4, nBLK), Rtmp
  logical :: unique, noNeighbor, stopRefinement, done
  logical :: oktest=.false., oktest_me=.false.
  type (adaptive_block_ptr) :: BlockPtr, tmpBlockPtr
  !---------------------------------------------------------------------------
  call set_oktest('amr',oktest,oktest_me)
  if(oktest .and. iProc == 0) then
     call write_myname; write(*,*) '| amr_physics beginning ...'
     call write_myname; write(*,*) '| '
  end if

  !--------------------------------------------------------------------
  ! Initial error checking and initializations
  !
  if(maxTotalBlocks > nBLK*nProc) then
     if(iProc == 0) then
        call write_myname; write(*,*) &
             '| maxTotalBlocks lowered to',nBLK*nProc,' from',maxTotalBlocks
     end if
     maxTotalBlocks = nBLK*nProc
  end if
  nCoarsened=0
  nRefined=0
  call set_body_flag

  !--------------------------------------------------------------------
  ! Compute and gather criteria for each block
  !
  call amr_criteria(refine_criteria)
  call MPI_ALLGATHER(refine_criteria(1,1), 4*nBLK, MPI_REAL, &
       refine_criteria_list(1,1,1), 4*nBLK, MPI_REAL, iComm, iError)

  !--------------------------------------------------------------------
  ! Sort gathered criteria
  !

  !Create lists to sort
  SortB=-1
  SortP=-1
  SortC=-1.0
  k=0
  do i=1,nBLK; do j=1,nProc
     BlockPtr%ptr => global_block_ptrs(i,j)%ptr
     if(associated(BlockPtr%ptr)) then ! ensure block is allocated
        k=k+1
        do n=1,nRefineCrit
           SortB(n,k)=i
           SortP(n,k)=j-1
           SortC(n,k)=refine_criteria_list(n,i,j)
        end do
     end if
  end do; end do
  nSort=k
  if(oktest .and. iProc == 0) then
     call write_myname; write(*,*) '| Blocks to sort=',nSort
  end if

  !Sort criteria lists
  do n=1,nRefineCrit
     do j=1,nSort-1
        do i=1,nSort-j
           if(SortC(n,i) > SortC(n,i+1)) then
              Rtmp         = SortC(n,i+1)
              SortC(n,i+1) = SortC(n,i)
              SortC(n,i)   = Rtmp

              Itmp         = SortB(n,i+1)
              SortB(n,i+1) = SortB(n,i)
              SortB(n,i)   = Itmp

              Itmp         = SortP(n,i+1)
              SortP(n,i+1) = SortP(n,i)
              SortP(n,i)   = Itmp
           end if
        end do
     end do
  end do

  !Put max value into next SortC entry for creating list for coarsening later
  SortC(:,nSort+1) = SortC(:,nSort)

  !--------------------------------------------------------------------
  ! Create list of possible refined blocks from sorted lists
  !
  SortIndex=nSort+1
  ListToRefine=-1
  k=0
  do i=nSort,1,-1
     do n=1,nRefineCrit
        if(SortIndex(n)<=i) CYCLE
        SortIndex(n)=min(i,SortIndex(n))

        !Check to see if block locked for refining
        BlockPtr%ptr => global_block_ptrs(SortB(n,i),SortP(n,i)+1)%ptr
        if(BlockPtr%ptr%LEV >= BlockPtr%ptr%LEVmax) CYCLE

        ! Do not refine blocks with criterion below the minimum
        if(SortC(n,i) < RefineCritMin_I(n)) CYCLE 
       
        unique=.true.
        do j=1,k
           if ( ListToRefine(1,j)==SortB(n,i) .and. &
                ListToRefine(2,j)==SortP(n,i)) unique=.false.
        end do
        if(unique) then
           k=k+1
           ListToRefine(1,k) = SortB(n,i)
           ListToRefine(2,k) = SortP(n,i)
        else
        end if

        !check to see if next is symmetric point and add
        do i1=1,i-1
           ! Exit loop if criteria is too different
           if(abs(SortC(n,i)-SortC(n,i-i1)) > (SortC(n,i)*1.E-2)) EXIT
           ! Cycle loop if not symmetric point
           if(abs(refine_criteria_list(4,SortB(n,i   ),SortP(n,i   )+1)- &
                  refine_criteria_list(4,SortB(n,i-i1),SortP(n,i-i1)+1)) > &
                 (refine_criteria_list(4,SortB(n,i   ),SortP(n,i   )+1)*1.E-6)) CYCLE

           !Check to see if block locked for refining
           BlockPtr%ptr => global_block_ptrs(SortB(n,i-i1),SortP(n,i-i1)+1)%ptr
           if(BlockPtr%ptr%LEV >= BlockPtr%ptr%LEVmax) CYCLE


           unique=.true.
           do j=1,k
              if ( ListToRefine(1,j)==SortB(n,i-i1) .and. &
                   ListToRefine(2,j)==SortP(n,i-i1) ) unique=.false.
           end do
           if(unique) then
              k=k+1
              ListToRefine(1,k) = SortB(n,i-i1)
              ListToRefine(2,k) = SortP(n,i-i1)
              if(SortIndex(n)==i-i1+1) SortIndex(n)=i-i1
           else
           end if
        end do
     end do
  end do
  nRefineMax=k
  if(oktest .and. iProc == 0) then
     call write_myname; write(*,*) '| Blocks possible to refine=',nRefineMax
  end if
  !--------------------------------------------------------------------
  ! Create new ordering for coarsening (coarsening needs AND of criteria)
  !
  !Create list to sort
  SortB=-1
  SortP=-1
  SortC(:,1:nSort) = -1.0
  k=0
  do i=1,nBLK; do j=1,nProc
     BlockPtr%ptr => global_block_ptrs(i,j)%ptr
     if(associated(BlockPtr%ptr)) then ! ensure block is allocated
        k=k+1
        SortB(1,k)=i
        SortP(1,k)=j-1
        SortC(1,k)=0.
        do n=1,nRefineCrit
           SortC(1,k)=SortC(1,k)+refine_criteria_list(n,i,j)/SortC(n,nSort+1)
        end do
     end if
  end do; end do

  !Sort criteria list
  n=1
  do j=1,nSort-1
     do i=1,nSort-j
        if(SortC(n,i) > SortC(n,i+1)) then
           Rtmp         = SortC(n,i+1)
           SortC(n,i+1) = SortC(n,i)
           SortC(n,i)   = Rtmp

           Itmp         = SortB(n,i+1)
           SortB(n,i+1) = SortB(n,i)
           SortB(n,i)   = Itmp

           Itmp         = SortP(n,i+1)
           SortP(n,i+1) = SortP(n,i)
           SortP(n,i)   = Itmp
        end if
     end do
  end do

  !--------------------------------------------------------------------
  ! Create list of possible coarsened blocks
  !
  nDesiredMax = int((percentCoarsenMax/100.)*nBlockALL)
  SortIndex=0
  ListToCoarsen=-1
  k=0
  n=1
  do i=1,nSort
     if(SortIndex(1)>=i) CYCLE
     SortIndex(1)=max(i,SortIndex(1))

     ! Do not coarsen blocks with criterion above maximum
     if(SortC(n,i) > CoarsenCritMax) CYCLE

     if(i>nDesiredMax)then
        EXIT    !Stop if passing nDesiredMax
     end if

     !Check to see if block locked for coarsening
     BlockPtr%ptr => global_block_ptrs(SortB(n,i),SortP(n,i)+1)%ptr
     if(BlockPtr%ptr%LEV <= BlockPtr%ptr%LEVmin) CYCLE

     unique=.true.
     do j=1,k
        if ( ListToCoarsen(1,j)==SortB(n,i) .and. &
             ListToCoarsen(2,j)==SortP(n,i)) unique=.false.
     end do
     if(unique) then
        !Remove if part of multilevel block
        if(associated(BlockPtr%ptr%parent)) then ! ensure parent exists
           if  (BlockPtr%ptr%parent%child1%used .and. &
                BlockPtr%ptr%parent%child2%used .and. &
                BlockPtr%ptr%parent%child3%used .and. &
                BlockPtr%ptr%parent%child4%used .and. &
                BlockPtr%ptr%parent%child5%used .and. &
                BlockPtr%ptr%parent%child6%used .and. &
                BlockPtr%ptr%parent%child7%used .and. &
                BlockPtr%ptr%parent%child8%used) then
              !All siblings have no children
           else
              !A sibling has children, can't coarsen
              unique=.false.
           end if
        else
           !Root cell, no parent, can't coarsen
           unique=.false.
        end if
     end if
     if(unique) then
        !Check for more refined neighbors
        do idir=1,27
           call findTreeNeighbor(BlockPtr,tmpBlockPtr,idir,noNeighbor)
           if(.not. noNeighbor) then
              if(associated(tmpBlockPtr%ptr)) then
                 if(associated(tmpBlockPtr%ptr%child1)) unique=.false.
              end if
           end if
        end do
     end if
     if(unique) then
        k=k+1
        ListToCoarsen(1,k) = SortB(n,i)
        ListToCoarsen(2,k) = SortP(n,i)
        !check to see if next is symmetric point and add
        do i1=1,nSort-i
           ! Exit loop if criteria is too different
           if(abs(SortC(n,i)-SortC(n,i+i1)) > (SortC(n,i)*1.E-2)) EXIT
           ! Cycle loop if not symmetric point
           if(abs(refine_criteria_list(4,SortB(n,i   ),SortP(n,i   )+1)- &
                  refine_criteria_list(4,SortB(n,i+i1),SortP(n,i+i1)+1)) > &
                 (refine_criteria_list(4,SortB(n,i   ),SortP(n,i   )+1)*1.E-6)) CYCLE

           !Check to see if block locked for coarsening
           BlockPtr%ptr => global_block_ptrs(SortB(n,i+i1),SortP(n,i+i1)+1)%ptr
           if(BlockPtr%ptr%LEV <= BlockPtr%ptr%LEVmin) CYCLE


           unique=.true.
           do j=1,k
              if ( ListToCoarsen(1,j)==SortB(n,i+i1) .and. &
                   ListToCoarsen(2,j)==SortP(n,i+i1) ) unique=.false.
           end do
           if(unique) then
              !Remove if part of multilevel block
              if(associated(BlockPtr%ptr%parent)) then ! ensure parent exists
                 if  (BlockPtr%ptr%parent%child1%used .and. &
                      BlockPtr%ptr%parent%child2%used .and. &
                      BlockPtr%ptr%parent%child3%used .and. &
                      BlockPtr%ptr%parent%child4%used .and. &
                      BlockPtr%ptr%parent%child5%used .and. &
                      BlockPtr%ptr%parent%child6%used .and. &
                      BlockPtr%ptr%parent%child7%used .and. &
                      BlockPtr%ptr%parent%child8%used) then
                    !All siblings have no children
                 else
                    !A sibling has children, can't coarsen
                    unique=.false.
                 end if
              else
                 !Root cell, no parent, can't coarsen
                 unique=.false.
              end if
           end if
           if(unique) then
              !Check for more refined neighbors
              do idir=1,27
                 call findTreeNeighbor(BlockPtr,tmpBlockPtr,idir,noNeighbor)
                 if(.not. noNeighbor) then
                    if(associated(tmpBlockPtr%ptr)) then
                       if(associated(tmpBlockPtr%ptr%child1)) unique=.false.
                    end if
                 end if
              end do
           end if
           if(unique) then
              k=k+1
              ListToCoarsen(1,k) = SortB(n,i+i1)
              ListToCoarsen(2,k) = SortP(n,i+i1)
              if(SortIndex(1)==i+i1-1) SortIndex(1)=i+i1
           end if
        end do
     end if
  end do
  nCoarsenMax=k
  if(oktest .and. iProc == 0) then
     call write_myname; write(*,*) '| Blocks possible to coarsen=',nCoarsenMax
  end if
  !--------------------------------------------------------------------
  ! Set all refinement and coarsening flags false
  !
  call clear_amr_flags

  !--------------------------------------------------------------------
  ! Flag blocks for coarsening
  !
  nDesired = int((percentCoarsen/100.)*nBlockALL)
  if(oktest .and. iProc == 0) then
     call write_myname; write(*,*) &
          '| Coarsening:',percentCoarsen,nDesired,nBlockALL,'(',percentCoarsenMax,')'
  end if
  if(nDesired>0)then
     do i=1,nCoarsenMax
        BlockPtr%ptr => global_block_ptrs(ListToCoarsen(1,i),ListToCoarsen(2,i)+1)%ptr
        call FlagBlockCoarsen
        if(nCoarsened+7 > nDesired) then
           ! check to see if a few more blocks should be coarsened to preserve symmetry
           do i1=1,nCoarsenMax-i
              if(abs(refine_criteria_list(4,ListToCoarsen(1,i   ),ListToCoarsen(2,i   )+1)- &
                   refine_criteria_list(4,ListToCoarsen(1,i+i1),ListToCoarsen(2,i+i1)+1) ) > &
                   (refine_criteria_list(4,ListToCoarsen(1,i   ),ListToCoarsen(2,i   )+1)*1.E-6)) EXIT
              BlockPtr%ptr => global_block_ptrs(ListToCoarsen(1,i+i1),ListToCoarsen(2,i+i1)+1)%ptr
              call FlagBlockCoarsen
           end do
           EXIT
        end if
     end do

     !Clean up extra unused coarsening flags
     do i=1,nCoarsenMax
        BlockPtr%ptr => global_block_ptrs(ListToCoarsen(1,i),ListToCoarsen(2,i)+1)%ptr
        if(associated(BlockPtr%ptr%parent)) then ! ensure parent exists
           if  (BlockPtr%ptr%parent%child1%coarsen .and. &
                BlockPtr%ptr%parent%child2%coarsen .and. &
                BlockPtr%ptr%parent%child3%coarsen .and. &
                BlockPtr%ptr%parent%child4%coarsen .and. &
                BlockPtr%ptr%parent%child5%coarsen .and. &
                BlockPtr%ptr%parent%child6%coarsen .and. &
                BlockPtr%ptr%parent%child7%coarsen .and. &
                BlockPtr%ptr%parent%child8%coarsen) then
              !All siblings coarsened, leave there
           else
              BlockPtr%ptr%coarsen = .false.
           end if
        end if
     end do

     !Check to see if flagging would result in body blocks at different levels
     do
        done = .true.
        do i=1,nBLK; do j=1,nProc
           BlockPtr%ptr => global_block_ptrs(i,j)%ptr
           if(associated(BlockPtr%ptr)) then ! ensure block is allocated
              if(BlockPtr%ptr%body .and. .not.BlockPtr%ptr%coarsen) then
                 do idir=1,27
                    call findTreeNeighbor(BlockPtr,tmpBlockPtr,idir,noNeighbor)
                    if(.not. noNeighbor)then
                       if(associated(tmpBlockPtr%ptr)) then
                          if(tmpBlockPtr%ptr%body .and. tmpBlockPtr%ptr%coarsen) then
                             tmpBlockPtr%ptr%parent%child1%coarsen = .false.
                             tmpBlockPtr%ptr%parent%child2%coarsen = .false.
                             tmpBlockPtr%ptr%parent%child3%coarsen = .false.
                             tmpBlockPtr%ptr%parent%child4%coarsen = .false.
                             tmpBlockPtr%ptr%parent%child5%coarsen = .false.
                             tmpBlockPtr%ptr%parent%child6%coarsen = .false.
                             tmpBlockPtr%ptr%parent%child7%coarsen = .false.
                             tmpBlockPtr%ptr%parent%child8%coarsen = .false.
                             nCoarsened=nCoarsened-7
                             currentBlocks = nBlockALL - nCoarsened + nRefined
                             done = .false.
                          end if
                       end if
                    end if
                 end do
              end if
           end if
        end do; end do
        if(done) EXIT
     end do
  end if

  if(oktest .and. iProc == 0) then
     i1=0; i2=0; call count_amr_flags(i1,i2)
     call write_myname; write(*,*) '| amr_flags after coarsening: ',i1,i2
  end if

  !--------------------------------------------------------------------
  ! Flag blocks for refining
  !
  stopRefinement = .false.
  nDesired = int((percentRefine/100.)*nBlockALL)
  if(oktest .and. iProc == 0) then
     call write_myname; write(*,*) &
          '| Refining:',percentRefine,nDesired,nBlockALL
  end if
  currentBlocks = nBlockALL - nCoarsened + nRefined
  if(nDesired>0)then
     do i=1,nRefineMax
        BlockPtr%ptr => global_block_ptrs(ListToRefine(1,i),ListToRefine(2,i)+1)%ptr
        call FlagBlockRefine
        if(stopRefinement) EXIT
        if(currentBlocks+7 > maxTotalBlocks) EXIT
        if(nRefined+7 > nDesired) then
           ! check to see if a few more blocks should be refined to preserve symmetry
           do i1=1,nRefineMax-i
              if(abs(refine_criteria_list(4,ListToRefine(1,i   ),ListToRefine(2,i   )+1)- &
                   refine_criteria_list(4,ListToRefine(1,i+i1),ListToRefine(2,i+i1)+1) ) > &
                   (refine_criteria_list(4,ListToRefine(1,i   ),ListToRefine(2,i   )+1)*1.E-6)) EXIT
              BlockPtr%ptr => global_block_ptrs(ListToRefine(1,i+i1),ListToRefine(2,i+i1)+1)%ptr
              call FlagBlockRefine
              if(stopRefinement) EXIT
           end do
           EXIT
        end if
     end do
  end if
  if(oktest .and. iProc == 0) then
     i1=0; i2=0; call count_amr_flags(i1,i2)
     call write_myname; write(*,*) '| amr_flags after refining:   ',i1,i2
  end if

  !--------------------------------------------------------------------
  ! Write final summary
  !
  if(iProc == 0.and.(lVerbose>0.or.oktest)) then
     call write_prefix; write(iUnitOut,*) '| amr_physics summary:'
     call write_prefix; write(iUnitOut,*) '| '
     call write_prefix; write(iUnitOut,'(a,i6)') ' |     blocks:',nBlockALL
     call write_prefix; write(iUnitOut,'(a,i6,a,f7.2,a,a,f7.2,a)') &
          ' |    removed:', &
          nCoarsened,'  (',100.*real(nCoarsened)/real(nBlockALL),'% change', &
          ', goal=',percentCoarsen,' %)'
     call write_prefix; write(iUnitOut,'(a,i6,a,f7.2,a,a,f7.2,a)') &
          ' |      added:', &
          nRefined  ,'  (',100.*real(nRefined  )/real(nBlockALL),'% change', &
          ', goal=',percentRefine,' %)'
     currentBlocks = nBlockALL - nCoarsened + nRefined
     call write_prefix; write(iUnitOut,'(a,i6,a,f7.2,a)') &
          ' |  new total:', &
          currentBlocks,'  (', &
          100.*real(currentBlocks-nBlockALL)/real(nBlockALL), &
          '% change)'
     call write_prefix; write(iUnitOut,*) '| '
     call write_prefix; write(iUnitOut,'(a,i6)') &
          ' | max blocks:',maxTotalBlocks
  end if
  if(oktest .and. iProc == 0 .and. lVerbose>0) then
     call write_prefix; write(iUnitOut,*) '| '
     call write_prefix; write(iUnitOut,*) '|     ... amr_physics ending.'
  end if

  !--------------------------------------------------------------------
  ! Copy coarsening and refine flags to lists and do coarsening and refinement
  !
  do i=1,nBLK; do j=1,nProc
     BlockPtr%ptr => global_block_ptrs(i,j)%ptr
     if(associated(BlockPtr%ptr)) then ! ensure block is allocated
        refine_list(i,j)  = BlockPtr%ptr%refine
        coarsen_list(i,j) = BlockPtr%ptr%coarsen
     end if
  end do; end do
  call parallel_coarsen
  call parallel_refine

Contains

  subroutine FlagBlockRefine
    if(BlockPtr%ptr%refine) RETURN  !already flagged
    BlockPtr%ptr%refine = .true.
    nRefined=nRefined+7
    currentBlocks = nBlockALL - nCoarsened + nRefined
    if(BlockPtr%ptr%coarsen) then
       BlockPtr%ptr%parent%child1%coarsen = .false.
       BlockPtr%ptr%parent%child2%coarsen = .false.
       BlockPtr%ptr%parent%child3%coarsen = .false.
       BlockPtr%ptr%parent%child4%coarsen = .false.
       BlockPtr%ptr%parent%child5%coarsen = .false.
       BlockPtr%ptr%parent%child6%coarsen = .false.
       BlockPtr%ptr%parent%child7%coarsen = .false.
       BlockPtr%ptr%parent%child8%coarsen = .false.
       nCoarsened=nCoarsened-7
       currentBlocks = nBlockALL - nCoarsened + nRefined
    end if
    call fix_octree_refine_flags(BlockPtr, nCoarsened, nRefined, &
         currentBlocks, stopRefinement)
    if(stopRefinement) then
       BlockPtr%ptr%refine = .false.
       nRefined=nRefined-7
       currentBlocks = nBlockALL - nCoarsened + nRefined
    end if
  end subroutine FlagBlockRefine

  subroutine FlagBlockCoarsen
    BlockPtr%ptr%coarsen = .true.
    if(associated(BlockPtr%ptr%parent)) then ! ensure parent exists
       if  (BlockPtr%ptr%parent%child1%coarsen .and. &
            BlockPtr%ptr%parent%child2%coarsen .and. &
            BlockPtr%ptr%parent%child3%coarsen .and. &
            BlockPtr%ptr%parent%child4%coarsen .and. &
            BlockPtr%ptr%parent%child5%coarsen .and. &
            BlockPtr%ptr%parent%child6%coarsen .and. &
            BlockPtr%ptr%parent%child7%coarsen .and. &
            BlockPtr%ptr%parent%child8%coarsen) then
          !All siblings coarsened, add to list
          nCoarsened=nCoarsened+7
          currentBlocks = nBlockALL - nCoarsened + nRefined
       end if
    end if
  end subroutine FlagBlockCoarsen

end subroutine amr_physics

!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
subroutine clear_amr_flags
  use ModParallel, ONLY : proc_dims
  use ModOctree
  implicit none

  integer :: i,j,k
  type (adaptive_block_ptr) :: octree

  do k=1,proc_dims(3); do j=1,proc_dims(2); do i=1,proc_dims(1)
     octree % ptr => octree_roots(i, j, k) % ptr
     call clear_octree_amr_flags(octree)
  end do; end do; end do

end subroutine clear_amr_flags

recursive subroutine clear_octree_amr_flags(octree)
  use ModOctree
  implicit none

  integer :: icube
  type (adaptive_block_ptr) :: octree, child

  if(associated(octree % ptr)) then
     octree % ptr % refine  = .false.
     octree % ptr % coarsen = .false.
     if(.not. (octree % ptr % used) ) then
        do icube = 1,8
           select case (icube)
           case (1)
              child % ptr => octree % ptr % child1
           case (2)
              child % ptr => octree % ptr % child2
           case (3)
              child % ptr => octree % ptr % child3
           case (4)
              child % ptr => octree % ptr % child4
           case (5)
              child % ptr => octree % ptr % child5
           case (6)
              child % ptr => octree % ptr % child6
           case (7)
              child % ptr => octree % ptr % child7
           case (8)
              cHild % ptr => octree % ptr % child8
           end select
           call clear_octree_amr_flags(child)
        end do
     end if
  end if
end subroutine clear_octree_amr_flags

!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
subroutine count_amr_flags(count1,count2)
  use ModParallel, ONLY : proc_dims
  use ModOctree
  implicit none

  integer, intent(inout) :: count1, count2
  integer :: i,j,k
  type (adaptive_block_ptr) :: octree

  do k=1,proc_dims(3); do j=1,proc_dims(2); do i=1,proc_dims(1)
     octree % ptr => octree_roots(i, j, k) % ptr
     call count_octree_amr_flags(octree, count1, count2)
  end do; end do; end do

end subroutine count_amr_flags

recursive subroutine count_octree_amr_flags(octree, count1, count2)
  use ModOctree
  implicit none

  integer, intent(inout) :: count1, count2
  integer :: icube
  type (adaptive_block_ptr) :: octree, child

  if(associated(octree % ptr)) then
     if(octree % ptr % coarsen) count1=count1+1
     if(octree % ptr % refine) count2=count2+1
     if(.not. (octree % ptr % used) ) then
        do icube = 1,8
           select case (icube)
           case (1)
              child % ptr => octree % ptr % child1
           case (2)
              child % ptr => octree % ptr % child2
           case (3)
              child % ptr => octree % ptr % child3
           case (4)
              child % ptr => octree % ptr % child4
           case (5)
              child % ptr => octree % ptr % child5
           case (6)
              child % ptr => octree % ptr % child6
           case (7)
              child % ptr => octree % ptr % child7
           case (8)
              cHild % ptr => octree % ptr % child8
           end select
           call count_octree_amr_flags(child, count1, count2)
        end do
     end if
  end if
end subroutine count_octree_amr_flags

!------------------------------------------------------------------------
! This routine takes a pointer to a block (assuming it has just been
!   flagged for refinement) and checks to see if any of it's neighbors
!   also needs to be refined (or not coarsened)
!------------------------------------------------------------------------
recursive subroutine fix_octree_refine_flags(inBlockPtr, nCoarsened, nRefined, &
     currentBlocks, stopRefinement)
  use ModMain
  use ModAMR, ONLY : maxTotalBlocks
  use ModOctree
  implicit none

  integer, intent(inout) :: nCoarsened, nRefined, currentBlocks
  logical, intent(inout) :: stopRefinement
  integer :: idir, iLEV1, iLEV2
  logical :: noNeighbor
  type (adaptive_block_ptr) :: inBlockPtr, outBlockPtr
  !---------------------------------------------------------------------------

  if(associated(inBlockPtr % ptr)) then
     !Check for valid neighbor level changes
     iLEV1 = inBlockPtr%ptr%LEV + 1
     do idir=1,27
        call findTreeNeighbor(inBlockPtr,outBlockPtr,idir,noNeighbor)
        if(.not. noNeighbor) then
           if(associated(outBlockPtr%ptr)) then
              if(.not. associated(outBlockPtr%ptr%child1)) then
                 iLEV2 = outBlockPtr%ptr%LEV
                 if(outBlockPtr%ptr%refine)  iLEV2=iLEV2+1
                 if(outBlockPtr%ptr%coarsen) iLEV2=iLEV2-1
                 if(iLEV1-iLEV2 > 1 .and. outBlockPtr%ptr%coarsen) then
                    !Turn off coarsening
                    nCoarsened=nCoarsened-7
                    currentBlocks = nBlockALL - nCoarsened + nRefined
                    outBlockPtr%ptr%parent%child1%coarsen = .false.
                    outBlockPtr%ptr%parent%child2%coarsen = .false.
                    outBlockPtr%ptr%parent%child3%coarsen = .false.
                    outBlockPtr%ptr%parent%child4%coarsen = .false.
                    outBlockPtr%ptr%parent%child5%coarsen = .false.
                    outBlockPtr%ptr%parent%child6%coarsen = .false.
                    outBlockPtr%ptr%parent%child7%coarsen = .false.
                    outBlockPtr%ptr%parent%child8%coarsen = .false.
                    iLEV2=iLEV2+1
                 end if
                 if(  (iLEV1-iLEV2 > 1) .or. &
                      (iLEV1-iLEV2 > 0 .and. &
                      inBlockPtr%ptr%body .and. outBlockPtr%ptr%body) ) then
                    !Flag block to refine and check it's neighbors
                    outBlockPtr%ptr%refine = .true.
                    nRefined=nRefined+7
                    currentBlocks = nBlockALL - nCoarsened + nRefined
                    if(currentBlocks > maxTotalBlocks) then
                       stopRefinement = .true.
                    end if


                    if(.not. stopRefinement) then
                       call fix_octree_refine_flags(outBlockPtr, nCoarsened, nRefined, &
                            currentBlocks, stopRefinement)
                    end if
                    if(stopRefinement) then
                       outBlockPtr%ptr%refine = .false.
                       nRefined=nRefined-7
                       currentBlocks = nBlockALL - nCoarsened + nRefined
                       return
                    end if
                 end if
              end if
           end if
        end if
     end do
  end if
end subroutine fix_octree_refine_flags

subroutine fixCheck
  use ModProcMH
  use ModMain, ONLY : nBLK
  use ModOctree
  use ModIO, ONLY: write_myname
  implicit none

  integer :: idir, inPE,inBLK, iLEV1,iLEV2, iCount, maxLev, curLev
  logical :: noNeighbor
  type (adaptive_block_ptr) :: inBlockPtr,outBlockPtr,tmpBlockPtr

  iCount = 0
  !find max level
  maxLev=0
  do inPE = 1,nProc
     do inBLK = 1,nBLK
        inBlockPtr%ptr => global_block_ptrs(inBLK,inPE)%ptr
        if(associated(inBlockPtr%ptr)) &
             maxLev = max(maxLev, inBlockPtr%ptr%LEV)
     end do
  end do
  !loop from max level down to zero
  do curLev = maxLev,1,-1
     do inPE = 1,nProc; do inBLK = 1,nBLK
        inBlockPtr%ptr => global_block_ptrs(inBLK,inPE)%ptr
        if(associated(inBlockPtr%ptr)) then
           if(inBlockPtr%ptr%LEV == curLev) then
              if(inBlockPtr%ptr%refine .and. inBlockPtr%ptr%coarsen) then
                 iCount = iCount+1
                 if(iProc==0) then
                    call write_myname; write(*,*) &
                         'fixCheck: found block coarsened and refined'
                 end if
              end if
              iLEV1 = inBlockPtr%ptr%LEV
              if(inBlockPtr%ptr%refine)  iLEV1 = iLEV1 + 1
              if(inBlockPtr%ptr%coarsen) iLEV1 = iLEV1 - 1
              do idir=1,27
                 call findTreeNeighbor(inBlockPtr,outBlockPtr,idir,noNeighbor)
                 if(.not. noNeighbor)then
                    if(associated(outBlockPtr%ptr)) then
                       if(.not. associated(outBlockPtr%ptr%child1)) then
                          iLEV2 = outBlockPtr%ptr%LEV
                          if(outBlockPtr%ptr%refine)  iLEV2 = iLEV2 + 1
                          if(outBlockPtr%ptr%coarsen) iLEV2 = iLEV2 - 1
                          if(  (iLEV1-iLEV2 > 1) .or. &
                               (iLEV1-iLEV2 > 0 .and. &
                               inBlockPtr%ptr%body .and. outBlockPtr%ptr%body) ) then
                             iCount = iCount +1
                             if(iProc==0) then
                                call write_myname; write(*,*) &
                                     'fixCheck: found large block level jump'
                             end if
                          end if
                       end if
                    end if
                 end if
              end do
           end if
        end if
     end do; end do
  end do
  if(iCount > 0) then
     if(iProc == 0) then
        write (*,*) '    FixCheck found blocks to fix', iCount
     end if
  end if

end subroutine fixCheck
