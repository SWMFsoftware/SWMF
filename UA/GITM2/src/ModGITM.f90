
module ModGITM

  use ModSize
  use ModPlanet

  implicit none

  real :: dt

  integer :: iCommGITM, iProc, nProcs

  real, dimension(-1:nAlts+2) :: Altitude, dAlt, RadialDistance, Gravity

  real, dimension(-1:nLats+2,-1:nAlts+2,nBlocksMax) :: dLonDist_GB, dLatDist_GB

  real, dimension(-1:nLons+2, nBlocksMax) :: Longitude
  real, dimension(-1:nLats+2, nBlocksMax) :: Latitude, TanLatitude

  real, allocatable :: Rho(:,:,:,:)
  real, allocatable :: Temperature(:,:,:,:)
  real, allocatable :: InvScaleHeight(:,:,:,:)
  real, allocatable :: Pressure(:,:,:,:)
  real, allocatable :: NDensity(:,:,:,:)
  real, allocatable :: eTemperature(:,:,:,:)
  real, allocatable :: ITemperature(:,:,:,:)
  real, allocatable :: IPressure(:,:,:,:)
  real, allocatable :: ePressure(:,:,:,:)

  real, allocatable :: LogRhoS(:,:,:,:,:)
  real, allocatable :: LogNS(:,:,:,:,:)
  real, allocatable :: VerticalVelocity(:,:,:,:,:)

  real, allocatable :: NDensityS(:,:,:,:,:)

  real, allocatable :: IDensityS(:,:,:,:,:)
  real, allocatable :: IRIDensity(:,:,:,:,:)

  real, allocatable :: KappaTemp(:,:,:,:)
  real, allocatable :: Ke(:,:,:,:)
  real, allocatable :: dKe(:,:,:,:)

  real, allocatable :: cp(:,:,:,:)
  real :: ViscCoef(nLons, nLats, 0:nAlts+1)

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2) :: &
       MeanIonMass, MeanMajorMass

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2) :: &
       e_gyro, i_gyro

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2, 3) :: Collisions

  real, allocatable :: B0(:,:,:,:,:)
  real, allocatable :: MLatitude(:,:,:,:)
  real :: MLT(-1:nLons+2, -1:nLats+2, -1:nAlts+2)
  real, allocatable :: MLongitude(:,:,:,:)
  real, allocatable :: DipAngle(:,:,:,:)
  real, allocatable :: DecAngle(:,:,:,:)
  real, allocatable :: cMax_GDB(:,:,:,:,:)

  real, dimension(1:nLons, 1:nLats, 1:nAlts, 3) :: &
       IonDrag, Viscosity

  real, allocatable :: Potential(:,:,:,:)

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2, 3) :: &
       ExB, EField

  real, dimension(-1:nLons+2, -1:nLats+2) :: &
       ElectronEnergyFlux, ElectronAverageEnergy

  real, allocatable :: Velocity(:,:,:,:,:)
  real, allocatable :: IVelocity(:,:,:,:,:)

  real, allocatable :: Emissions(:,:,:,:,:)

  real :: TempUnit

  real :: LocalTime(-1:nLons+2)

  real :: SubsolarLatitude, SubsolarLongitude 
  real :: MagneticPoleColat, MagneticPoleLon


  integer, parameter :: iEast_ = 1, iNorth_ = 2, iUp_ = 3, iMag_ = 4

  integer, parameter :: iVIN_ = 1, iVEN_ = 2, iVEI_ = 3

contains
  !=========================================================================
  subroutine init_mod_gitm

    use ModUtilities, ONLY: check_allocate
    integer :: iError
    !----------------------------------------------------------------------
    if(allocated(Rho)) return
    allocate(Rho(-1:nLons+2, -1:nLats+2, -1:nAlts+2, nBlocks),stat=iError)
    call check_allocate(iError,'Rho')
    allocate(Temperature(-1:nLons+2, -1:nLats+2, -1:nAlts+2, nBlocks),stat=iError)
    call check_allocate(iError,'Temperature')
    allocate(InvScaleHeight(-1:nLons+2, -1:nLats+2, -1:nAlts+2, nBlocks),stat=iError)
    call check_allocate(iError,'InvScaleHeight')
    allocate(Pressure(-1:nLons+2, -1:nLats+2, -1:nAlts+2, nBlocks),stat=iError)
    call check_allocate(iError,'Pressure')
    allocate(NDensity(-1:nLons+2, -1:nLats+2, -1:nAlts+2, nBlocks),stat=iError)
    call check_allocate(iError,'NDensity')
    allocate(eTemperature(-1:nLons+2, -1:nLats+2, -1:nAlts+2, nBlocks),stat=iError)
    call check_allocate(iError,'eTemperature')
    allocate(ITemperature(-1:nLons+2, -1:nLats+2, -1:nAlts+2, nBlocks),stat=iError)
    call check_allocate(iError,'ITemperature')
    allocate(IPressure(-1:nLons+2, -1:nLats+2, -1:nAlts+2, nBlocks),stat=iError)
    call check_allocate(iError,'IPressure')
    allocate(ePressure(-1:nLons+2, -1:nLats+2, -1:nAlts+2, nBlocks),stat=iError)
    call check_allocate(iError,'ePressure')
    allocate(LogRhoS(-1:nLons+2, -1:nLats+2, -1:nAlts+2, &
       nSpecies, nBlocksMax),stat=iError)
    call check_allocate(iError,'LogRhoS')
    allocate(LogNS(-1:nLons+2, -1:nLats+2, -1:nAlts+2, &
       nSpecies, nBlocksMax),stat=iError)
    call check_allocate(iError,'LogNS')
    allocate(VerticalVelocity(-1:nLons+2, -1:nLats+2, -1:nAlts+2, &
       nSpecies, nBlocksMax),stat=iError)
    call check_allocate(iError,'VerticalVelocity')
    allocate(NDensityS(-1:nLons+2, -1:nLats+2, -1:nAlts+2, &
       nSpeciesTotal, nBlocksMax),stat=iError)
    call check_allocate(iError,'NDensityS')
    allocate(IDensityS(-1:nLons+2, -1:nLats+2, -1:nAlts+2, &
       nIons, nBlocksMax),stat=iError)
    call check_allocate(iError,'IDensityS')
    allocate(IRIDensity(-1:nLons+2, -1:nLats+2, -1:nAlts+2, &
       nIons, nBlocksMax),stat=iError)
    call check_allocate(iError,'IRIDensity')
    allocate(KappaTemp(nLons, nLats, 0:nAlts+1, nBlocks),stat=iError)
    call check_allocate(iError,'KappaTemp')
    allocate(Ke(nLons, nLats, 0:nAlts+1, nBlocks),stat=iError)
    call check_allocate(iError,'Ke')
    allocate(dKe(nLons, nLats, 0:nAlts+1, nBlocks),stat=iError)
    call check_allocate(iError,'dKe')
    allocate(cp(nLons,nLats,0:nAlts+1,nBlocks),stat=iError)
    call check_allocate(iError,'cp')
    allocate(B0(-1:nLons+2,-1:nLats+2,-1:nAlts+2,4,nBlocks),stat=iError)
    call check_allocate(iError,'B0')
    allocate(MLatitude(-1:nLons+2,-1:nLats+2,-1:nAlts+2,nBlocks),stat=iError)
    call check_allocate(iError,'MLatitude')
    allocate(MLongitude(-1:nLons+2,-1:nLats+2,-1:nAlts+2,nBlocks),stat=iError)
    call check_allocate(iError,'MLongitude')
    allocate(DipAngle(-1:nLons+2,-1:nLats+2,-1:nAlts+2,nBlocks),stat=iError)
    call check_allocate(iError,'DipAngle')
    allocate(DecAngle(-1:nLons+2,-1:nLats+2,-1:nAlts+2,nBlocks),stat=iError)
    call check_allocate(iError,'DecAngle')
    allocate(cMax_GDB(0:nLons+1,0:nLats+1,0:nAlts+1,3,nBlocks),stat=iError)
    call check_allocate(iError,'cMax_GDB')
    allocate(Potential(-1:nLons+2, -1:nLats+2, -1:nAlts+2, nBlocks),stat=iError)
    call check_allocate(iError,'Potential')
    allocate(Velocity(-1:nLons+2, -1:nLats+2, -1:nAlts+2, 3, nBlocks),stat=iError)
    call check_allocate(iError,'Velocity')
    allocate(IVelocity(-1:nLons+2, -1:nLats+2, -1:nAlts+2, 3, nBlocks),stat=iError)
    call check_allocate(iError,'IVelocity')
    allocate(Emissions(nLons,nLats,nAlts,nEmissions,nBlocks),stat=iError)
    call check_allocate(iError,'Emissions')
  end subroutine init_mod_gitm
  !=========================================================================
  subroutine clean_mod_gitm

    if(.not.allocated(Rho)) return
    deallocate(Rho)
    deallocate(Temperature)
    deallocate(InvScaleHeight)
    deallocate(Pressure)
    deallocate(NDensity)
    deallocate(eTemperature)
    deallocate(ITemperature)
    deallocate(IPressure)
    deallocate(ePressure)
    deallocate(LogRhoS)
    deallocate(LogNS)
    deallocate(VerticalVelocity)
    deallocate(NDensityS)
    deallocate(IDensityS)
    deallocate(IRIDensity)
    deallocate(KappaTemp)
    deallocate(Ke)
    deallocate(dKe)
    deallocate(cp)
    deallocate(B0)
    deallocate(MLatitude)
    deallocate(MLongitude)
    deallocate(DipAngle)
    deallocate(DecAngle)
    deallocate(cMax_GDB)
    deallocate(Potential)
    deallocate(Velocity)
    deallocate(IVelocity)
    deallocate(Emissions)
  end subroutine clean_mod_gitm
  !=========================================================================
end module ModGITM
