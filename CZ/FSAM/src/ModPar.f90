module ModPar

  ! Primary parameters:
  
  integer, parameter :: inmax = 101, jnmax = 517, knmax = 773, ijknmax = 773
  real(8), parameter    :: gg = 6.67231D-8, rgas = 8.3170000D+07, kboltz = 1.381D-16, &
       mproton = 1.673D-24
  real(8), parameter    :: pi = 3.1415926535898D0, tiny = 1.0D-99, huge = 1.0D+99

  ! NOTE that inmax-5 must be divisible by nproc1
  ! jnmax-5 must be divisible by nproc2
  ! knmax-5 must be divisible by nproc
  integer, parameter :: nproc1 = 12, nproc2 = 64
  integer, parameter :: nproc = nproc1*nproc2
  integer, parameter :: in = (inmax-5)/nproc1+5, jn = (jnmax-5)/nproc2+5, &
       kns = (knmax-5)/nproc, kn = knmax, ijkn = 773
  integer :: myid, myid1, myid2

end module ModPar
