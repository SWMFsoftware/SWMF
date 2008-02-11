!^CFG COPYRIGHT UM
subroutine build_octree_roots
  use ModProcMH
  use ModMain, ONLY : nBLK
  use ModParallel, ONLY : proc_dims
  use ModOctree
  use ModMpi
  implicit none

  type (adaptive_block_ptr) :: octree
  integer :: i, j, k, iPE, iBLK, iLEV, iLEVmin, iLEVmax, iError

  !\  
  ! Allocate and initialize octree block data structure.
  !/
  allocate( octree_roots(proc_dims(1),proc_dims(2),proc_dims(3)), &
       stat=iError )
  if (iError > 0) write(*,*) "build_octree_roots: PE = ",iProc, &
       " allocation error for octree_roots"

  allocate( global_block_ptrs(nBLK,nProc), stat=iError )
  if (iError > 0) write(*,*) "build_octree_roots: PE = ",iProc, &
       " allocation error for global_block_ptrs"

  allocate( blocknumber_ptrs(nBLK*nProc), stat=iError )
  if (iError > 0) write(*,*) "build_octree_roots: PE = ",iProc, &
       " allocation error for blocknumber_ptrs"

  do iPE = 1, nProc; do iBLK = 1, nBLK
     nullify (global_block_ptrs(iBLK,iPE) % ptr)
  end do; end do

  do i = 1, nProc*nBLK
     nullify (blocknumber_ptrs(i) % ptr)
  end do

  ! All octree roots go to processor 0
  iBLK = 0
  do k = 1, proc_dims(3)
     do j = 1, proc_dims(2)
        do i = 1, proc_dims(1)
           iBLK = iBLK+1
           iPE  = 0
           iLEV = 0
           iLEVmin = 0
           iLEVmax = 999
           nullify(octree % ptr)
           call initialize_octree_block(octree, iPE, iBLK, iLEV, iLEVmin, iLEVmax)
           octree_roots(i, j, k) % ptr => octree  % ptr
           global_block_ptrs(iBLK, iPE+1) % ptr => octree % ptr
        end do
     end do
  end do

end subroutine build_octree_roots

subroutine initialize_octree_block(octree, iPE, iBLK, iLEV, iLEVmin, iLEVmax)
  use ModOctree
  implicit none

  type (adaptive_block_ptr) :: octree
  integer, intent(in) :: iPE, iBLK, iLEV, iLEVmin, iLEVmax

  integer :: ierror

  if (associated(octree % ptr)) deallocate (octree % ptr)
  allocate ( octree % ptr, stat=ierror )
  if (ierror > 0) write(*,*) "initialize_octree_block: ", &
       & " allocation error for octree"
  nullify (octree % ptr % parent)
  nullify (octree % ptr % child1)
  nullify (octree % ptr % child2)
  nullify (octree % ptr % child3)
  nullify (octree % ptr % child4)
  nullify (octree % ptr % child5)
  nullify (octree % ptr % child6)
  nullify (octree % ptr % child7)
  nullify (octree % ptr % child8)

  octree % ptr % number  = 0
  octree % ptr % child_number = 0
  octree % ptr % PE      = iPE
  octree % ptr % BLK     = iBLK
  octree % ptr % LEV     = iLEV
  octree % ptr % LEVmin  = iLEVmin
  octree % ptr % LEVmax  = iLEVmax
  octree % ptr % used    = .true.
  octree % ptr % refine  = .false.
  octree % ptr % coarsen = .false.
  octree % ptr % body    = .false.
  octree % ptr % IsExtraBoundaryOrPole    = .false.   
  octree % ptr % IsOuterBoundary    = .false.  

end subroutine initialize_octree_block

subroutine refine_octree_block(octree, iPEs, iBLKs, fromPE, fromBLK)
  use ModOctree
  implicit none

  type (adaptive_block_ptr) :: octree
  integer, intent(in) :: fromPE, fromBLK
  integer, intent(in), dimension(8) :: iPEs, iBLKs

  integer :: iLEV, iLEVmin, iLEVmax
  type (adaptive_block_ptr) :: child

  if (associated(octree % ptr)) then
     iLEV = octree % ptr % LEV + 1
     iLEVmin = octree % ptr % LEVmin
     iLEVmax = octree % ptr % LEVmax

     nullify(child % ptr)
     call initialize_octree_block(child, iPEs(1), iBLKs(1), iLEV, iLEVmin, iLEVmax)
     octree % ptr % child1 =>  child % ptr
     child % ptr % child_number = 1

     nullify(child % ptr)
     call initialize_octree_block(child, iPEs(2), iBLKs(2), iLEV, iLEVmin, iLEVmax)
     octree % ptr % child2 =>  child % ptr
     child % ptr % child_number = 2

     nullify(child % ptr)
     call initialize_octree_block(child, iPEs(3), iBLKs(3), iLEV, iLEVmin, iLEVmax)
     octree % ptr % child3 =>  child % ptr
     child % ptr % child_number = 3

     nullify(child % ptr)
     call initialize_octree_block(child, iPEs(4), iBLKs(4), iLEV, iLEVmin, iLEVmax)
     octree % ptr % child4 =>  child % ptr
     child % ptr % child_number = 4

     nullify(child % ptr)
     call initialize_octree_block(child, iPEs(5), iBLKs(5), iLEV, iLEVmin, iLEVmax)
     octree % ptr % child5 =>  child % ptr
     child % ptr % child_number = 5

     nullify(child % ptr)
     call initialize_octree_block(child, iPEs(6), iBLKs(6), iLEV, iLEVmin, iLEVmax)
     octree % ptr % child6 =>  child % ptr
     child % ptr % child_number = 6

     nullify(child % ptr)
     call initialize_octree_block(child, iPEs(7), iBLKs(7), iLEV, iLEVmin, iLEVmax)
     octree % ptr % child7 =>  child % ptr
     child % ptr % child_number = 7

     nullify(child % ptr)
     call initialize_octree_block(child, iPEs(8), iBLKs(8), iLEV, iLEVmin, iLEVmax)
     octree % ptr % child8 =>  child % ptr
     child % ptr % child_number = 8

     octree % ptr % child1 % parent => octree % ptr
     octree % ptr % child2 % parent => octree % ptr
     octree % ptr % child3 % parent => octree % ptr
     octree % ptr % child4 % parent => octree % ptr
     octree % ptr % child5 % parent => octree % ptr
     octree % ptr % child6 % parent => octree % ptr
     octree % ptr % child7 % parent => octree % ptr
     octree % ptr % child8 % parent => octree % ptr

     octree % ptr % used    = .false.
     octree % ptr % refine  = .false.
     octree % ptr % coarsen = .false.
     octree % ptr % body    = .false.
     octree % ptr % IsExtraBoundaryOrPole = .false. 
     octree % ptr % IsOuterBoundary = .false. 
     octree % ptr % number  = -1
     octree % ptr % PE      = -1
     octree % ptr % BLK     = -1

     if(fromBLK>=0 .and. fromPE>=0) &
          nullify(global_block_ptrs(fromBLK, fromPE+1) % ptr)

     global_block_ptrs(iBLKs(1), iPEs(1)+1) % ptr => &
          octree % ptr % child1
     global_block_ptrs(iBLKs(2), iPEs(2)+1) % ptr => &
          octree % ptr % child2
     global_block_ptrs(iBLKs(3), iPEs(3)+1) % ptr => &
          octree % ptr % child3
     global_block_ptrs(iBLKs(4), iPEs(4)+1) % ptr => &
          octree % ptr % child4
     global_block_ptrs(iBLKs(5), iPEs(5)+1) % ptr => &
          octree % ptr % child5
     global_block_ptrs(iBLKs(6), iPEs(6)+1) % ptr => &
          octree % ptr % child6
     global_block_ptrs(iBLKs(7), iPEs(7)+1) % ptr => &
          octree % ptr % child7
     global_block_ptrs(iBLKs(8), iPEs(8)+1) % ptr => &
          octree % ptr % child8
  end if

end subroutine refine_octree_block

subroutine move_octree_block(iBlockFrom, iProcFrom, iBlockTo, iProcTo)
  use ModProcMH
  use ModOctree
  implicit none

  integer, intent(in) :: iBlockFrom, iProcFrom, iBlockTo, iProcTo
  type (adaptive_block_ptr) :: octree

  octree % ptr => global_block_ptrs(iBlockFrom, iProcFrom+1) % ptr
  if (associated(octree % ptr)) then
     octree % ptr % PE  = iProcTo
     octree % ptr % BLK = iBlockTo
  else
     write(*,*)'iProc,iBlockFrom,iProcFrom,iBlockTo,iProcTo', &
          iProc,iBlockFrom,iProcFrom,iBlockTo,iProcTo
     call stop_mpi("ERROR in move_octree_block")
  end if
  global_block_ptrs(iBlockTo, iProcTo+1) % ptr => octree % ptr
  nullify (global_block_ptrs(iBlockFrom, iProcFrom+1) % ptr)

end subroutine move_octree_block

subroutine coarsen_octree_block(octree, iPEs, iBLKs)
  use ModOctree
  implicit none

  type (adaptive_block_ptr) :: octree
  integer, intent(in), dimension(8) :: iPEs, iBLKs

  type (adaptive_block_ptr) :: child

  if (associated(octree % ptr)) then
     deallocate (octree % ptr % child1)
     deallocate (octree % ptr % child2)
     deallocate (octree % ptr % child3)
     deallocate (octree % ptr % child4)
     deallocate (octree % ptr % child5)
     deallocate (octree % ptr % child6)
     deallocate (octree % ptr % child7)
     deallocate (octree % ptr % child8)
     !     nullify (octree % ptr % child1)
     !     nullify (octree % ptr % child2)
     !     nullify (octree % ptr % child3)
     !     nullify (octree % ptr % child4)
     !     nullify (octree % ptr % child5)
     !     nullify (octree % ptr % child6)
     !     nullify (octree % ptr % child7)
     !     nullify (octree % ptr % child8)

     octree % ptr % number  = -2
     octree % ptr % PE      = iPEs(1)
     octree % ptr % BLK     = iBLKs(1)
     !     octree % ptr % LEV     = octree % ptr % LEV
     octree % ptr % used    = .true.
     octree % ptr % refine  = .false.
     octree % ptr % coarsen = .false.
     octree % ptr % body    = .false.
     octree % ptr % IsExtraBoundaryOrPole = .false.   
     octree % ptr % IsOuterBoundary = .false.  

     global_block_ptrs(iBLKs(1), iPEs(1)+1) % ptr => octree % ptr

     nullify (global_block_ptrs(iBLKs(2), iPEs(2)+1) % ptr)
     nullify (global_block_ptrs(iBLKs(3), iPEs(3)+1) % ptr)
     nullify (global_block_ptrs(iBLKs(4), iPEs(4)+1) % ptr)
     nullify (global_block_ptrs(iBLKs(5), iPEs(5)+1) % ptr)
     nullify (global_block_ptrs(iBLKs(6), iPEs(6)+1) % ptr)
     nullify (global_block_ptrs(iBLKs(7), iPEs(7)+1) % ptr)
     nullify (global_block_ptrs(iBLKs(8), iPEs(8)+1) % ptr)

  end if

end subroutine coarsen_octree_block

!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
subroutine set_levels
  use ModAMR

  call set_min_level(min_block_level)
  call set_max_level(max_block_level)
  if (fix_body_level) call set_body_level

end subroutine set_levels

!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
subroutine set_min_level(iLEVmin)
  use ModParallel, ONLY : proc_dims
  use ModOctree
  implicit none

  integer, intent(in) :: iLEVmin
  integer :: i,j,k
  type (adaptive_block_ptr) :: octree

  do k=1,proc_dims(3); do j=1,proc_dims(2); do i=1,proc_dims(1)
     octree % ptr => octree_roots(i, j, k) % ptr
     call set_octree_min_level(octree, iLEVmin)
  end do; end do; end do

end subroutine set_min_level

recursive subroutine set_octree_min_level(octree, iLEVmin)
  use ModOctree
  implicit none

  integer, intent(in) :: iLEVmin
  integer :: icube
  type (adaptive_block_ptr) :: octree, child

  if (associated(octree % ptr)) then
     if (iLEVmin < 0) then
        octree % ptr % LEVmin = octree % ptr % LEV
     else
        octree % ptr % LEVmin = iLEVmin
     end if
     if (.not. (octree % ptr % used) ) then
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
              child % ptr => octree % ptr % child8
           end select
           call set_octree_min_level(child, iLEVmin)
        end do
     end if
  end if
end subroutine set_octree_min_level

!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
subroutine set_max_level(iLEVmax)
  use ModParallel, ONLY : proc_dims
  use ModOctree
  implicit none

  integer, intent(in) :: iLEVmax
  integer :: i,j,k
  type (adaptive_block_ptr) :: octree

  do k=1,proc_dims(3); do j=1,proc_dims(2); do i=1,proc_dims(1)
     octree % ptr => octree_roots(i, j, k) % ptr
     call set_octree_max_level(octree, iLEVmax)
  end do; end do; end do

end subroutine set_max_level

recursive subroutine set_octree_max_level(octree, iLEVmax)
  use ModOctree
  implicit none

  integer, intent(in) :: iLEVmax
  integer :: icube
  type (adaptive_block_ptr) :: octree, child

  if (associated(octree % ptr)) then
     if (iLEVmax < 0) then
        octree % ptr % LEVmax = octree % ptr % LEV
     else
        octree % ptr % LEVmax = iLEVmax
     end if
     if (.not. (octree % ptr % used) ) then
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
              child % ptr => octree % ptr % child8
           end select
           call set_octree_max_level(child, iLEVmax)
        end do
     end if
  end if
end subroutine set_octree_max_level

!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
subroutine set_body_level
  use ModParallel, ONLY : proc_dims
  use ModOctree
  implicit none

  integer :: i,j,k
  type (adaptive_block_ptr) :: octree

  call set_body_flag
  do k=1,proc_dims(3); do j=1,proc_dims(2); do i=1,proc_dims(1)
     octree % ptr => octree_roots(i, j, k) % ptr
     call set_octree_body_level(octree)
  end do; end do; end do

end subroutine set_body_level

recursive subroutine set_octree_body_level(octree)
  use ModOctree
  implicit none

  integer :: icube
  type (adaptive_block_ptr) :: octree, child

  if (associated(octree % ptr)) then
     if (octree % ptr % body) then
        octree % ptr % LEVmin = octree % ptr % LEV
        octree % ptr % LEVmax = octree % ptr % LEV
     end if
     if (.not. (octree % ptr % used) ) then
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
              child % ptr => octree % ptr % child8
           end select
           call set_octree_body_level(child)
        end do
     end if
  end if
end subroutine set_octree_body_level

!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
subroutine set_body_flag
  use ModProcMH
  use ModSize,ONLY:nBLK
  use ModMain, ONLY : DoFixOuterBoundary,East_,Top_,&       
       ExtraBc_,DoFixExtraBoundaryOrPole,TypeBC_I     
  use ModGeometry, ONLY :IsBoundaryBlock_IB,DoFixExtraBoundary_B
  use ModGeometry, ONLY : BodyFlg_B
  use ModOctree
  use ModMpi
  implicit none

  integer :: inPE,inBLK, iError
  logical,dimension(nBLK) :: DoFixBoundary_B
  logical, dimension(:,:), allocatable :: tmp_logical_list
  type (adaptive_block_ptr) :: inBlockPtr

  allocate( tmp_logical_list(nBLK,nProc), stat=ierror )
  if (ierror > 0) write(*,*) "set_body_flag: PE = ",iProc, &
       " allocation error: tmp_logical_list"
  tmp_logical_list = .false.

  !collect BodyFlags to tmp_logical_list
  call MPI_ALLGATHER(BodyFlg_B(1), nBLK, MPI_LOGICAL, &
       tmp_logical_list(1,1), nBLK, MPI_LOGICAL, iComm, iError)

  !set body flag
  do inPE = 1,nProc
     do inBLK = 1,nBLK
        inBlockPtr%ptr => global_block_ptrs(inBLK,inPE)%ptr
        if (associated(inBlockPtr%ptr)) &
             inBlockPtr%ptr%body = tmp_logical_list(inBLK,inPE)
     end do
  end do

  if(DoFixExtraBoundaryOrPole)then

     tmp_logical_list = .false.

     call MPI_ALLGATHER(DoFixExtraBoundary_B(1), nBLK, MPI_LOGICAL, &
          tmp_logical_list(1,1), nBLK, MPI_LOGICAL, iComm, iError)
     
     !set body flag
     do inPE = 1,nProc
        do inBLK = 1,nBLK
           inBlockPtr%ptr => global_block_ptrs(inBLK,inPE)%ptr
           if (associated(inBlockPtr%ptr)) &
                inBlockPtr%ptr%IsExtraBoundaryOrPole = &
                tmp_logical_list(inBLK,inPE)
        end do
     end do
  end if
  if(DoFixOuterBoundary)then
     do inBLK=1,nBLK
        DoFixBoundary_B(inBLK)=any(IsBoundaryBlock_IB(East_:Top_,inBLK).and.&
        TypeBC_I(East_:Top_)/='float')
     end do
     tmp_logical_list = .false.

     !collect BodyFlags to tmp_logical_list
     call MPI_ALLGATHER(DoFixBoundary_B(1), nBLK, MPI_LOGICAL, &
          tmp_logical_list(1,1), nBLK, MPI_LOGICAL, iComm, iError)
     
     !set body flag
     do inPE = 1,nProc
        do inBLK = 1,nBLK
           inBlockPtr%ptr => global_block_ptrs(inBLK,inPE)%ptr
           if (associated(inBlockPtr%ptr)) &
                inBlockPtr%ptr%IsOuterBoundary = tmp_logical_list(inBLK,inPE)
        end do
     end do
  end if

  deallocate(tmp_logical_list)

end subroutine set_body_flag
!^CFG IF DEBUGGING BEGIN
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
subroutine analyze_neighbors
  use ModProcMH
  use ModMain, ONLY : nBlockMax, unusedBLK
  use ModParallel, ONLY : NOBLK
  use ModOctree
  use ModIO, ONLY : iUnitOut, write_prefix
  use ModMpi
  implicit none

  integer, dimension(-1:2) :: NMsame,NMtotal
  integer :: i,j,k, n, inBLK,inPE, returnedLEV, NMblk,NMholes, iError
  integer, dimension(4) :: returnedPE,returnedBLK,returnedCHILD
  logical :: topBlockFound
  type (adaptive_block_ptr) :: octree
  integer, save, dimension(:), allocatable :: PEblk,PEholes,PEsame,PEtotal
  logical :: oktest, oktest_me

  !---------------------------------------------------------------------------
  call set_oktest('analyze_neighbors',oktest,oktest_me)
  if(.not.oktest) return

  if(.not.allocated(PEblk))then
     allocate( PEblk(nProc), stat=iError )
     if(iError>0)write(*,*) "analyze_neighbors: PE = ",iProc," allocation error: PEblk"
  end if
  PEblk = 0

  if(.not.allocated(PEholes))then
     allocate( PEholes(nProc), stat=iError )
     if(iError>0)write(*,*) "analyze_neighbors: PE = ",iProc," allocation error: PEholes"
  end if
  PEholes = 0

  if(.not.allocated(PEsame))then
     allocate( PEsame(nProc), stat=iError )
     if(iError>0)write(*,*) "analyze_neighbors: PE = ",iProc," allocation error: PEsame"
  end if
  PEsame = 0

  if(.not.allocated(PEtotal))then
     allocate( PEtotal(nProc), stat=iError )
     if(iError>0)write(*,*) "analyze_neighbors: PE = ",iProc," allocation error: PEtotal"
  end if
  PEtotal = 0

  !loop over blocks adding up information
  NMblk=0
  NMsame=0
  NMtotal=0
  inPE = iProc
  topBlockFound = .false.
  NMholes=0
  if(iProc==0)then
     call write_prefix; write(iUnitOut,*) &
          'Starting analyze_neighbors: nBlockMax=',nBlockMax
  end if
  do inBLK = nBlockMax,1,-1
     if(.not.unusedBLK(inBLK)) topBlockFound = .true.
     if(topBlockFound .and. unusedBLK(inBLK)) NMholes=NMholes+1
     if (unusedBLK(inBLK)) CYCLE
     NMblk=NMblk+1
     do i=-1,1; do j=-1,1; do k=-1,1
        if((abs(i)+abs(j)+abs(k)) == 1) then
           call treeNeighbor(inPE,inBLK,i,j,k, &
                returnedPE,returnedBLK,returnedCHILD,returnedLEV)
           do n=1,4
              if (returnedBLK(n) == NOBLK) CYCLE
              NMtotal(returnedLEV)=NMtotal(returnedLEV)+1
              if (returnedPE(n) == iProc) NMsame(returnedLEV)=NMsame(returnedLEV)+1
           end do
        end if
     end do; end do; end do
  end do

  NMtotal(2)=NMtotal(-1)+NMtotal(0)+NMtotal(1)
  NMsame (2)=NMsame (-1)+NMsame (0)+NMsame (1)

!!$  !each processor write its values
!!$  do i=-1,1
!!$     write(*,'(i2,a,i3,a,i3,a,i3,a,i6,a,i6,a,f6.2,a)') i, &
!!$          ' PE=',iProc,' Blocks:',NMblk,' Holes:',NMholes, &
!!$          '  Total Messages:',NMtotal(i),'  Copies:',NMsame(i), &
!!$          ' (',100.*real(NMsame(i))/(real(NMtotal(i))+1.E-8),'%)'
!!$  end do
!!$  write(*,'(2x,a,i3,a,i3,a,i3,a,i6,a,i6,a,f6.2,a)') ' PE=',iProc,' Blocks:',NMblk,&
!!$       ' Holes:',NMholes,'  Total Messages:',NMtotal(2),'  Copies:',NMsame(2), &
!!$       ' (',100.*real(NMsame(2))/real(NMtotal(2)),'%)'

  !processor 0 write summary
  call MPI_allgather(NMblk,     1,MPI_INTEGER,PEblk(1),  1,MPI_INTEGER,iComm,iError)
  call MPI_allgather(NMholes,   1,MPI_INTEGER,PEholes(1),1,MPI_INTEGER,iComm,iError)
  call MPI_allgather(NMsame(2), 1,MPI_INTEGER,PEsame(1), 1,MPI_INTEGER,iComm,iError)
  call MPI_allgather(NMtotal(2),1,MPI_INTEGER,PEtotal(1),1,MPI_INTEGER,iComm,iError)
  if(iProc==0) then
     do i=1,nProc
        write(*,'(2x,a,i3,a,i3,a,i3,a,i6,a,i6,a,f6.2,a)') ' PE=',i-1,' Blocks:',PEblk(i),&
             ' Holes:',PEholes(i),'  Total Messages:',PEtotal(i),'  Copies:',PEsame(i), &
             ' (',100.*real(PEsame(i))/real(PEtotal(i)),'%)'
     end do
     write(*,*)' '
  end if

end subroutine analyze_neighbors
!^CFG END DEBUGGING
