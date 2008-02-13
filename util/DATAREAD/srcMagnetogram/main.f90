program Harmonics

  ! Transform raw magnetogram file into spherical harmonics file

  use ModMagHarmonics
  implicit none

  integer :: iError
  !----------------------------------------------------------------------------
  call MPI_INIT(iError)
  iComm = MPI_COMM_WORLD
  call MPI_COMM_RANK (iComm, iProc, iError)
  call MPI_COMM_SIZE (iComm, nProc, iError)

  call read_raw_magnetogram 
  call calc_harmonics

  call MPI_finalize(iError)
  
end program Harmonics
