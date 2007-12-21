
module ModRIM
  use ModSizeRIM

  implicit none

  real, dimension(0:nLons+1,nLats) :: &
       Latitude, Longitude, Potential, AveE, Eflux, SigmaH, SigmaP, &
       dLatitude, dLongitude

  real, allocatable :: AllLons(:)

  logical :: IsTimeAccurate
  real :: ThetaTilt, DipoleStrength

  integer, dimension(7) :: TimeArray
  integer :: nSolve=0

end module ModRIM
