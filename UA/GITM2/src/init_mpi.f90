
subroutine init_mpi

  use ModGITM
  use ModTime
  use ModMpi

  implicit none

  integer :: iError

  ! setup mpi
  call MPI_INIT(iError)
  iCommGITM = MPI_COMM_WORLD
  call MPI_COMM_RANK(iCommGITM, iProc, iError)
  call MPI_COMM_SIZE(iCommGITM, nProcs, iError)

end subroutine init_mpi
