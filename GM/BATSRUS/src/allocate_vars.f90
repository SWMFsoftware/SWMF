!^CFG COPYRIGHT UM
subroutine allocate_vars
  use ModProcMH, ONLY: nProc
  use ModSize, ONLY: MaxBlock
  use ModParallel, ONLY: iBlock_A, iProc_A, iBlockRestartALL_A, &
       nBlockMax_P, MaxBlockDisp_P
  implicit none

  integer :: jProc
  integer :: MaxBlockALL

  !\
  ! allocate index arrays
  !/

  MaxBlockALL = MaxBlock * nProc
  allocate(iBlock_A(MaxBlockALL), iProc_A(MaxBlockALL), &
       iBlockRestartALL_A(MaxBlockALL))

  ! allocate and initialize the displacement and the maximum number
  ! of received data for the special using mpi_allgatherv instead of
  ! MPI_allgather. The displacement is always MaxBlock, the maximum
  ! received data varies together with the value of nBlockMax.

  allocate(nBlockMax_P(0:nProc-1), MaxBlockDisp_P(0:nProc-1))
  do jProc = 0, nProc-1
     MaxBlockDisp_P(jProc) = jProc*MaxBlock
  end do
  nBlockMax_P = 0

end subroutine allocate_vars
