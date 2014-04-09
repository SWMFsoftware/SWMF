module ModBoundary
  use ModPar,   ONLY: jn, kn
  implicit none

  integer :: niib, noib, nijb, nojb
  real    :: prbt, prtp, srbt, srtp
  
  real, dimension(1:jn,1:kn) :: siib1
  
end module ModBoundary
