!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!This code is a copyright protected software (c) 2002- University of Michigan
module ModParallel
  use ModSize, ONLY: MaxBlock
  implicit none
  save

  !\
  ! Neighbor solution block refinement levels
  ! ( 0=neighbors at same level, 
  !  -1=neighbors at lower level,
  !  +1=neighbors at higher level,
  !  NOBLK=no neighbors).
  !/
  integer, parameter :: NOBLK=-100

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

  subroutine init_mod_parallel

    ! allocate and initialize the displacement and the maximum number
    ! of received data for the special using mpi_allgatherv instead of
    ! MPI_allgather. The displacement is always MaxBlock, the maximum
    ! received data varies together with the value of nBlockMax.

    use ModProcMH, ONLY: nProc
    use ModSize, ONLY: MaxBlock

    integer :: jProc
    !---------------------------------------------------------------


    if(allocated(nBlockMax_P)) RETURN
    
    allocate(nBlockMax_P(0:nProc-1), MaxBlockDisp_P(0:nProc-1))
    do jProc = 0, nProc-1
       MaxBlockDisp_P(jProc) = jProc*MaxBlock
    end do
    nBlockMax_P = 0

  end subroutine init_mod_parallel

end module ModParallel
