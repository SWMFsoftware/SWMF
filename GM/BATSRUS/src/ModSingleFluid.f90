module ModSingleFluid

  ! Define single fluid values for multifluid parameters 

  implicit none

  integer,           parameter :: nFluid = 1
  character (len=1), parameter :: NameFluid_I(nFluid) = (/' '/)
  character (len=3), parameter :: TypeFluid_I(nFluid) = (/'ion'/)
  integer,           parameter :: iVarFluid_I(nFluid) = (/ 0 /)

end module ModSingleFluid

