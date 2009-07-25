module ModSingleFluid

  ! Define single fluid values for multifluid parameters 

  implicit none

  integer,          parameter :: nFluid    = 1
  integer,          parameter :: IonFirst_ = 1
  integer,          parameter :: IonLast_  = 1
  logical,          parameter :: IsMhd     = .true.
  character(len=1), parameter :: NameFluid_I(nFluid) = (/' '/)
  
  ! Default is proton mass, but it can be changed
  real :: MassFluid_I(nFluid) = 1.0

  ! Set impossible value for indexes that are not used in most equation modules
  integer, parameter :: Pe_ = 1, Ppar_ = 1, Hyp_ = 1, Erad_ = 1
  integer, parameter :: Palfp_ = 1, Palfm_ = 1

end module ModSingleFluid

