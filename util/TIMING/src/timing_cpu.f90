!^CFG COPYRIGHT UM
function timing_cpu()
  use ModKind
  implicit none
  real(Real8_)           :: timing_cpu
  real(Real8_), external :: MPI_WTIME

  timing_cpu=MPI_WTIME()

end function timing_cpu

