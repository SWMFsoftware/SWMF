module ModFluid

  implicit none

  integer, parameter :: &
       Rho_   = 1,          &
       RhoUx_ = 2, Ux_ = 2, &
       RhoUy_ = 3, Uy_ = 3, &
       RhoUz_ = 4, Uz_ = 4, &
       Bx_    = 5, &
       By_    = 6, &
       Bz_    = 7, &
       p_     = 8

  ! This allows to calculate RhoUx_ as RhoU_+x_ and so on.
  integer, parameter :: U_ = Ux_ - 1, RhoU_ = RhoUx_-1, B_ = Bx_-1

end module ModFluid

!==============================================================
module ModSingleFluid

  use ModFluid
  implicit none

  integer,           parameter :: nFluid = 1
  character (len=1), parameter :: NameFluid_I(nFluid) = (/' '/)
  character (len=3), parameter :: TypeFluid_I(nFluid) = (/'ion'/)
  integer,           parameter :: iVarFluid_I(nFluid) = (/ 0 /)

end module ModSingleFluid

