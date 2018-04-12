!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!\
! MODULE ------------------------------------------------
!/

module ModVertical

  use ModSizeGitm, only: nAlts
  use ModPlanet, only: nSpecies, nIonsAdvect, nSpeciesTotal, nIons

  implicit none

  real, dimension(-1:nAlts+2)   :: dAlt_F = 0.0, InvDAlt_F = 0.0, Altitude_G = 0.0, Gravity_G = 0.0
  real, dimension(-1:nAlts+2)   :: InvRadialDistance_C = 0.0, dAlt_C = 0.0, Cv_1D = 0.0
  real, dimension(-1:nAlts+2)   :: LogRho = 0.0, Temp = 0.0, MeanMajorMass_1d = 0.0, Gamma_1d = 0.0
  real, dimension(-1:nAlts+2,3) :: Vel_GD = 0.0
  real, dimension(-1:nAlts+2,3) :: IVel = 0.0
  
  real, dimension(-1:nAlts+2)             :: NewLogRho = 0.0, NewTemp = 0.0
  real, dimension(-1:nAlts+2,3)           :: NewVel_GD = 0.0
  real, dimension(-1:nAlts+2,nSpecies) :: NewLogNS = 0.0, NewVertVel = 0.0
!  real, dimension(-1:nAlts+2,nIonsAdvect) :: NewLogINS
!  real, dimension(-1:nAlts+2,nIonsAdvect) :: LogINS

  real, dimension(-1:nAlts+2,nIons) :: NewLogINS = 0.0
  real, dimension(-1:nAlts+2,nIons) :: LogINS = 0.0
!  real, dimension(1:2,nIons) :: OldLogINSc


  real, dimension(-1:nAlts+2, nSpecies) :: LogNS = 0.0, VertVel = 0.0
  real, dimension(nAlts, nSpecies) :: NDensityS_1D = 0.0

  real, dimension(0:nAlts+1) :: cMax = 0.0
  
  real :: SZAVertical = 0.0
  real :: MLatVertical = 0.0
 
  integer :: iLon1D, iLat1D, iBlock1D

  real :: Heating(nAlts)
  real, dimension(-1:nAlts+2) :: EddyCoef_1d = 0.0
  real, dimension( 0:nAlts+1) :: ViscCoef_1d = 0.0
  real, dimension( 0:nAlts+1,1:nSpecies) :: ViscCoefS_1d = 0.0
  real :: KappaTemp_1d( 0:nAlts+1) = 0.0

  real :: Centrifugal, Coriolis, Lat, Lon, dAltdlon_1D, dAltdLat_1D

end module ModVertical

