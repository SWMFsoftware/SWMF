module ModSingleFluid

  ! Define single fluid values for multifluid parameters 

  implicit none

  integer,           parameter :: nFluid = 1
  integer,           parameter :: nIonFluid = 1
  logical,           parameter :: UseMultiIon = .false.
  character (len=1), parameter :: NameFluid_I(nFluid) = (/' '/)
  character (len=3), parameter :: TypeFluid_I(nFluid) = (/'ion'/)

  ! Default is proton mass, but it can be changed
  real :: MassFluid_I(nFluid) = 1.0

end module ModSingleFluid

