module ModPar

  ! Primary parameters:
  integer, parameter :: inmax = 101, jnmax = 101, knmax = 101
  real, parameter    :: gg = 6.67231D-8, rgas = 8.3170000D+07, kboltz = 1.381D-16, &
       mproton = 1.673D-24
  real, parameter    :: pi = 3.1415926535898D0, tiny = 1.0D-99, huge = 1.0D+99

  ! NOTE that inmax-5 must be divisible by nproc1
  ! jnmax-5 must be divisible by nproc2
  ! knmax-5 must be divisible by nproc
  integer, parameter :: nproc1 = 3, nproc2 = 2
  integer, parameter :: nproc = nproc1*nproc2
  integer, parameter :: in = (inmax-5)/nproc1+5, jn = (jnmax-5)/nproc2+5, &
       kns = (knmax-5)/nproc, kn = knmax
  integer, parameter :: ijkn = max(in,jn,kn)
  integer :: myid, myid1, myid2

end module ModPar
