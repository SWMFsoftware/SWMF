!^CFG COPYRIGHT UM
subroutine build_octree_roots
  use ModProcMH
  use ModMain, ONLY : nBLK,nBlockAll
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
  nBlockAll=iBLK
end subroutine build_octree_roots

subroutine initialize_octree_block(octree, iPE, iBLK, iLEV, iLEVmin, iLEVmax)
  use ModOctree
  implicit none

  type (adaptive_block_ptr) :: octree
  integer, intent(in) :: iPE, iBLK, iLEV, iLEVmin, iLEVmax

  integer :: ierror,iChild

  if (associated(octree % ptr)) deallocate (octree % ptr)
  allocate ( octree % ptr, stat=ierror )
  if (ierror > 0) write(*,*) "initialize_octree_block: ", &
       & " allocation error for octree"
  nullify (octree % ptr % parent%ptr)
  do iChild=1,8
     nullify (octree % ptr % child(iChild)%ptr)
  end do

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
  use ModMain,ONLY:nBlockAll
  implicit none

  type (adaptive_block_ptr) :: octree
  integer, intent(in) :: fromPE, fromBLK
  integer, intent(in), dimension(8) :: iPEs, iBLKs

  integer :: iLEV, iLEVmin, iLEVmax,iChild
  type (adaptive_block_ptr) :: child

  if (associated(octree % ptr)) then
     iLEV = octree % ptr % LEV + 1
     iLEVmin = octree % ptr % LEVmin
     iLEVmax = octree % ptr % LEVmax
     do iChild=1,8
        nullify(child % ptr)
        call initialize_octree_block(child, &
             iPEs(iChild), iBLKs(iChild), iLEV, iLEVmin, iLEVmax)
        octree % ptr % child(iChild)%ptr =>  child % ptr
        child % ptr % child_number = iChild
        
        octree % ptr % child(iChild)%ptr % parent%ptr => octree % ptr
     end do
     
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

     do iChild=1,8
        global_block_ptrs(iBLKs(iChild), iPEs(iChild)+1) % ptr => &
             octree % ptr % child(iChild)%ptr
     end do
     nBlockAll=nBlockAll+7
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
  use ModMain,ONLY:nBlockAll
  implicit none

  type (adaptive_block_ptr) :: octree
  integer, intent(in), dimension(8) :: iPEs, iBLKs
  integer::iChild

  type (adaptive_block_ptr) :: child

  if (associated(octree % ptr)) then
     do iChild=1,8
        deallocate (octree % ptr % child(iChild)%ptr)
     end do
     !     nullify (octree % ptr % child(1)%ptr)
     !     nullify (octree % ptr % child(2)%ptr)
     !     nullify (octree % ptr % child(3)%ptr)
     !     nullify (octree % ptr % child(4)%ptr)
     !     nullify (octree % ptr % child(5)%ptr)
     !     nullify (octree % ptr % child(6)%ptr)
     !     nullify (octree % ptr % child(7)%ptr)
     !     nullify (octree % ptr % child(8)%ptr)

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
     do iChild=2,8
        nullify (global_block_ptrs(iBLKs(iChild), iPEs(iChild)+1) % ptr)
     end do
     nBlockAll=nBlockAll-7
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
           child % ptr => octree % ptr % child(iCube)%ptr
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
           child % ptr => octree % ptr % child(iCube)%ptr
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
           child % ptr => octree % ptr % child(iCube)%ptr
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
  call MPI_ALLGATHER(BodyFlg_B, nBLK, MPI_LOGICAL, &
       tmp_logical_list, nBLK, MPI_LOGICAL, iComm, iError)

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

     call MPI_ALLGATHER(DoFixExtraBoundary_B, nBLK, MPI_LOGICAL, &
          tmp_logical_list, nBLK, MPI_LOGICAL, iComm, iError)
     
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
     call MPI_ALLGATHER(DoFixBoundary_B, nBLK, MPI_LOGICAL, &
          tmp_logical_list, nBLK, MPI_LOGICAL, iComm, iError)
     
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
  call MPI_allgather(NMblk,     1,MPI_INTEGER,PEblk,  1,MPI_INTEGER,&
       iComm,iError)
  call MPI_allgather(NMholes,   1,MPI_INTEGER,PEholes,1,MPI_INTEGER,&
       iComm,iError)
  call MPI_allgather(NMsame(2), 1,MPI_INTEGER,PEsame, 1,MPI_INTEGER,&
       iComm,iError)
  call MPI_allgather(NMtotal(2),1,MPI_INTEGER,PEtotal,1,MPI_INTEGER,&
       iComm,iError)
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
