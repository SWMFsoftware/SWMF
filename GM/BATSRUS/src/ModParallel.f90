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

  integer, dimension(MaxBlock) :: &
       neiLtop, neiLbot, neiLeast, neiLwest, neiLnorth, neiLsouth

  integer, dimension(1:6,MaxBlock):: neiLEV

  !\
  ! Neighbor processor and block numbers (a value of NOBLK
  ! means not used).  As only one level change is permitted
  ! between neighboring solution blocks, there are either 1 or 4
  ! neighboring blocks in each of the six directions.
  !/
  integer, dimension(4,MaxBlock) :: &
       neiPtop, neiPbot, neiPeast, neiPwest, neiPnorth, neiPsouth, &
       neiBtop, neiBbot, neiBeast, neiBwest, neiBnorth, neiBsouth

  integer, dimension(4,1:6,MaxBlock) :: neiPE, neiBLK

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

    allocate(nBlockMax_P(0:nProc-1), MaxBlockDisp_P(0:nProc-1))
    do jProc = 0, nProc-1
       MaxBlockDisp_P(jProc) = jProc*MaxBlock
    end do
    nBlockMax_P = 0

    call test_stop(NameSub, DoTest)
  end subroutine init_mod_parallel
  !============================================================================

end module ModParallel
!==============================================================================
