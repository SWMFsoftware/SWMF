!\
! MODULE ------------------------------------------------
!/

module ModVertical

  use ModSize, only: nAlts
  use ModPlanet, only: nSpecies, nIonsAdvect

  implicit none

  real, dimension(-1:nAlts+2)   :: LogRho, Temp
  real, dimension(-1:nAlts+2,3) :: Vel_GD
  real, dimension(-1:nAlts+2,3) :: IVel

  real, dimension(-1:nAlts+2)             :: NewLogRho, NewTemp
  real, dimension(-1:nAlts+2,3)           :: NewVel_GD
  real, dimension(-1:nAlts+2,nSpecies)    :: NewLogNS, NewVertVel
  real, dimension(-1:nAlts+2,nIonsAdvect) :: NewLogINS
  real, dimension(-1:nAlts+2,nIonsAdvect) :: LogINS

  real, dimension(-1:nAlts+2, nSpecies) :: LogNS, NDensityS, VertVel

  real, dimension(0:nAlts+1) :: cMax

  real :: Heating(nAlts)
  real :: KappaTemp(1:nAlts+1)
!  real :: Kappa0(nAlts), KappaNS(nSpecies,nSpecies)
  real :: KappaNS(nAlts,nSpecies,nSpecies)

  real :: Centrifugal, Coriolis, Lat, Lon

end module ModVertical

