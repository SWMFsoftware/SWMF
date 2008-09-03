
module ModGITM

  use ModSizeGitm
  use ModPlanet

  implicit none

  real :: dt = 0.0

  integer :: iCommGITM, iProc, nProcs

  real, dimension(-1:nLons+2,-1:nLats+2,-1:nAlts+2,nBlocksMax) :: &
       dLonDist_GB, InvDLonDist_GB, &
       dLonDist_FB, InvDLonDist_FB, &
       dLatDist_GB, InvDLatDist_GB, &
       dLatDist_FB, InvDLatDist_FB, &
       Altitude_GB, dAlt_GB, RadialDistance_GB, InvRadialDistance_GB, &
       Gravity_GB

  ! Topography
  real, dimension(nLons,nLats,nAlts,nBlocksMax) :: dAltDLon_CB, dAltDLat_CB

  real, dimension(-1:nLons+2, nBlocksMax) :: Longitude
  real, dimension(-1:nLats+2, nBlocksMax) :: Latitude, TanLatitude, CosLatitude

  real, dimension(nLons, nBlocksMax) :: GradLonM_CB, GradLon0_CB, GradLonP_CB
  real, dimension(nLats, nBlocksMax) :: GradLatM_CB, GradLat0_CB, GradLatP_CB

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2, nBlocksMax) :: &
       Rho, Temperature, InvScaleHeight, Pressure, &
       NDensity, eTemperature, ITemperature, &
       IPressure, ePressure

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2, &
       nSpecies, nBlocksMax) :: &
       LogRhoS, LogNS, VerticalVelocity

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2, &
       nSpeciesTotal, nBlocksMax) :: &
       NDensityS

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2, &
       nIons, nBlocksMax) :: &
       IDensityS, IRIDensity

  real :: Gamma(-1:nLons+2, -1:nLats+2, -1:nAlts+2,nBlocksMax)

  real, dimension(nLons, nLats, 0:nAlts+1, nBlocksMax) :: &
       KappaTemp, Ke,dKe

  real, dimension(nLons, nLats, nBlocksMax) :: &
       SurfaceAlbedo=0.0, SurfaceTemp=0.0

  real :: cp(nLons, nLats, 0:nAlts+1,nBlocksMax)
  real :: ViscCoef(nLons, nLats, 0:nAlts+1)

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2) :: &
       MeanIonMass, MeanMajorMass

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2) :: &
       e_gyro, i_gyro

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2, 3) :: Collisions

  real :: B0(-1:nLons+2, -1:nLats+2, -1:nAlts+2, 4, nBlocksMax)
  real :: MLatitude(-1:nLons+2, -1:nLats+2, -1:nAlts+2, nBlocksMax)
  real :: MLT(-1:nLons+2, -1:nLats+2, -1:nAlts+2)
  real :: MLongitude(-1:nLons+2, -1:nLats+2, -1:nAlts+2, nBlocksMax)
  real :: DipAngle(-1:nLons+2, -1:nLats+2, -1:nAlts+2, nBlocksMax)
  real :: DecAngle(-1:nLons+2, -1:nLats+2, -1:nAlts+2, nBlocksMax)

  real :: b0_d1(-1:nLons+2, -1:nLats+2, -1:nAlts+2, 3, nBlocksMax)
  real :: b0_d2(-1:nLons+2, -1:nLats+2, -1:nAlts+2, 3, nBlocksMax)
  real :: b0_d3(-1:nLons+2, -1:nLats+2, -1:nAlts+2, 3, nBlocksMax)
  real :: b0_e1(-1:nLons+2, -1:nLats+2, -1:nAlts+2, 3, nBlocksMax)
  real :: b0_e2(-1:nLons+2, -1:nLats+2, -1:nAlts+2, 3, nBlocksMax)
  real :: b0_e3(-1:nLons+2, -1:nLats+2, -1:nAlts+2, 3, nBlocksMax)
  real :: b0_cD(-1:nLons+2, -1:nLats+2, -1:nAlts+2, nBlocksMax)
  real :: b0_Be3(-1:nLons+2, -1:nLats+2, -1:nAlts+2, nBlocksMax)

  real :: cMax_GDB(0:nLons+1, 0:nLats+1, 0:nAlts+1, 3, nBlocksMax)

  real, dimension(1:nLons, 1:nLats, 1:nAlts, 3) :: &
       IonDrag, Viscosity

  real, dimension(1:nLons, 1:nLats, 1:nAlts, nSpecies) :: &
       VerticalIonDrag

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2, nBlocksMax) :: &
       Potential

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2, 3) :: &
       ExB, EField

  real, dimension(-1:nLons+2, -1:nLats+2) :: &
       ElectronEnergyFlux, ElectronAverageEnergy

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2, 3, nBlocksMax) :: &
       Velocity, IVelocity

  real :: Emissions(nLons, nLats, nAlts, nEmissions, nBlocksMax)

  real, dimension(nLons,nLats,nAlts,nEmissionWavelengths,nBlocksMax) :: &
       vEmissionRate
  
  real, dimension(nLons,nLats,nAlts,nPhotoBins,nBlocksMax) :: &
       PhotoElectronDensity,PhotoElectronRate,PhotoEFluxU,PhotoEFluxD
  
  real, dimension(nLons,nLats,nAlts,nBlocksMax,2) :: PhotoEFluxTotal

  real, dimension(nPhotoBins)                 :: PhotoEBins
  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2) :: TempUnit

  real :: LocalTime(-1:nLons+2)

  real :: SubsolarLatitude, SubsolarLongitude 
  real :: MagneticPoleColat, MagneticPoleLon


  integer, parameter :: iEast_ = 1, iNorth_ = 2, iUp_ = 3, iMag_ = 4
  integer, parameter :: iVIN_ = 1, iVEN_ = 2, iVEI_ = 3

contains
  !=========================================================================
  subroutine init_mod_gitm
  end subroutine init_mod_gitm
  !=========================================================================
  subroutine clean_mod_gitm
  end subroutine clean_mod_gitm
  !=========================================================================
end module ModGITM
