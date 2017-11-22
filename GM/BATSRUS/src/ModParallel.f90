!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModParallel

  use BATL_lib, ONLY: &
       test_start, test_stop
  use ModSize, ONLY: MaxBlock
  use BATL_tree, ONLY: Unset_
  implicit none
  save

  !\
  ! Neighbor solution block refinement levels
  ! ( 0=neighbors at same level,
  !  -1=neighbors at lower level,
  !  +1=neighbors at higher level,
  !  NOBLK=no neighbors).
  !/
  integer, parameter :: NOBLK = Unset_

  integer, allocatable :: neiLtop(:)
  integer, allocatable :: neiLbot(:)
  integer, allocatable :: neiLeast(:)
  integer, allocatable :: neiLwest(:)
  integer, allocatable :: neiLnorth(:)
  integer, allocatable :: neiLsouth(:)

  integer, allocatable :: neiLEV(:,:)

  !\
  ! Neighbor processor and block numbers (a value of NOBLK
  ! means not used).  As only one level change is permitted
  ! between neighboring solution blocks, there are either 1 or 4
  ! neighboring blocks in each of the six directions.
  !/
  integer, allocatable :: neiPtop(:,:)
  integer, allocatable :: neiPbot(:,:)
  integer, allocatable :: neiPeast(:,:)
  integer, allocatable :: neiPwest(:,:)
  integer, allocatable :: neiPnorth(:,:)
  integer, allocatable :: neiPsouth(:,:)
  integer, allocatable :: neiBtop(:,:)
  integer, allocatable :: neiBbot(:,:)
  integer, allocatable :: neiBeast(:,:)
  integer, allocatable :: neiBwest(:,:)
  integer, allocatable :: neiBnorth(:,:)
  integer, allocatable :: neiBsouth(:,:)

  integer, allocatable :: neiPE(:,:,:)
  integer, allocatable :: neiBLK(:,:,:)

  ! used by mpi_allgatherv for a more efficient replacment of mpi_allgather
  integer, allocatable :: nBlockMax_P(:), MaxBlockDisp_P(:)

contains
  !============================================================================

  subroutine init_mod_parallel

    ! allocate and initialize the displacement and the maximum number
    ! of received data for the special using mpi_allgatherv instead of
    ! MPI_allgather. The displacement is always MaxBlock, the maximum
    ! received data varies together with the value of nBlockMax.

    use ModProcMH, ONLY: nProc
    use ModSize, ONLY: MaxBlock

    integer :: jProc
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'init_mod_parallel'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    if(allocated(nBlockMax_P)) RETURN

    allocate(neiLtop(MaxBlock))
    allocate(neiLbot(MaxBlock))
    allocate(neiLeast(MaxBlock))
    allocate(neiLwest(MaxBlock))
    allocate(neiLnorth(MaxBlock))
    allocate(neiLsouth(MaxBlock))
    allocate(neiLEV(1:6,MaxBlock))
    allocate(neiPtop(4,MaxBlock))
    allocate(neiPbot(4,MaxBlock))
    allocate(neiPeast(4,MaxBlock))
    allocate(neiPwest(4,MaxBlock))
    allocate(neiPnorth(4,MaxBlock))
    allocate(neiPsouth(4,MaxBlock))
    allocate(neiBtop(4,MaxBlock))
    allocate(neiBbot(4,MaxBlock))
    allocate(neiBeast(4,MaxBlock))
    allocate(neiBwest(4,MaxBlock))
    allocate(neiBnorth(4,MaxBlock))
    allocate(neiBsouth(4,MaxBlock))
    allocate(neiPE(4,1:6,MaxBlock))
    allocate(neiBLK(4,1:6,MaxBlock))

    allocate(MaxBlockDisp_P(0:nProc-1))
    do jProc = 0, nProc-1
       MaxBlockDisp_P(jProc) = jProc*MaxBlock
    end do

    allocate(nBlockMax_P(0:nProc-1))
    nBlockMax_P = 0

    call test_stop(NameSub, DoTest)
  end subroutine init_mod_parallel
  !============================================================================
  subroutine clean_mod_parallel

    if(.not.allocated(nBlockMax_P)) RETURN

    deallocate(neiLtop)
    deallocate(neiLbot)
    deallocate(neiLeast)
    deallocate(neiLwest)
    deallocate(neiLnorth)
    deallocate(neiLsouth)
    deallocate(neiLEV)
    deallocate(neiPtop)
    deallocate(neiPbot)
    deallocate(neiPeast)
    deallocate(neiPwest)
    deallocate(neiPnorth)
    deallocate(neiPsouth)
    deallocate(neiBtop)
    deallocate(neiBbot)
    deallocate(neiBeast)
    deallocate(neiBwest)
    deallocate(neiBnorth)
    deallocate(neiBsouth)
    deallocate(neiPE)
    deallocate(neiBLK)
    deallocate(MaxBlockDisp_P)
    deallocate(nBlockMax_P)

  end subroutine clean_mod_parallel

end module ModParallel
