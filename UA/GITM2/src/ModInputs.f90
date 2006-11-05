
module ModInputs

  use ModConstants
  use ModPlanet, only : nSpecies
  use ModIoUnit, only : UnitTmp_

  implicit none

  integer, parameter        :: iCharLen_     = 100

  integer                   :: iOutputUnit_  = UnitTmp_
  integer                   :: iInputUnit_   = UnitTmp_
  integer                   :: iRestartUnit_ = UnitTmp_

  integer                   :: iLogFileUnit_ = 92
  logical                   :: IsOpenLogFile = .false.
  real                      :: DtLogFile = 60.0

  integer, parameter        :: nInputMaxLines = 10000
  integer                   :: nInputLines
  character (len=iCharLen_) :: cInputText(nInputMaxLines)

  character (len=iCharLen_) :: cInputFile = "UAM.in"

  character (len=iCharLen_) :: cAMIEFileSouth
  character (len=iCharLen_) :: cAMIEFileNorth

  character (len=iCharLen_) :: TypeLimiter = "minmod"

  integer, dimension(7) :: iStartTime

  real :: CPUTimeMax = 1.0e32

  logical :: UseVariableInputs = .false.
  logical :: IsFramework = .false.

  logical :: DoRestart = .false.

  integer :: iAltTest = -1
  integer :: iDebugLevel = 0
  logical :: UseBarriers = .false.
  integer :: nSteps = 10

  real :: CFL = 0.25

  integer :: nOutputTypes = 0
  integer, parameter :: nMaxOutputTypes = 50
  character (len=iCharLen_), dimension(nMaxOutputTypes) :: OutputType

  real :: DtPlot(nMaxOutputTypes)

  real :: DtRestart   = 60.0*60.0
  real :: DtReport    =  1.0*60.0
  real :: DtAurora    = 60.0*1.0
  real :: DtPotential = 60.0*1.0

  real :: f107  = 150.0
  real :: f107a = 150.0
  integer :: iModelSolar = 0

  real :: AltMin =  95.0 * 1000.0
  real :: AltMax = 500.0 * 1000.0

  logical :: UseMinMod = .false.
  logical :: UseMC = .true.

  real :: ConcentrationLatitude = 45.0
  real :: StretchingPercentage = 0.0
  real :: StretchingFactor = 1.0

  real :: TempMax = 1000.0
  real :: TempMin =  200.0
  real :: TempWidth    =  25.0*1e3
  real :: TempHeight   = 150.0*1e3

  real :: LogRho0

  integer :: nBlocksLon = 2
  integer :: nBlocksLat = 2

  logical :: Is1D = .false.
  logical :: IsFullSphere = .false.
  real    :: LonStart = 0.0
  real    :: LatStart = -pi/4.0
  real    :: LatEnd   =  pi/4.0

  logical :: UseStretchedAltitude = .true.

  logical :: UseApex = .true.
  logical :: UseMSIS = .true.
  logical :: UseIRI  = .true.
  logical :: UseStatisticalModelsOnly = .false.
  real    :: DtStatisticalModels = 3600.0

  !\
  ! These are things for the ion precipitation for the April 2002 storm:
  !/

  logical :: UseIonPrecipitation = .false.
  character (len=iCharLen_) :: IonIonizationFilename
  character (len=iCharLen_) :: IonHeatingRateFilename

  !\
  ! These are logicals to turn off and on forcing terms:
  !/

  logical :: UsePressureGradient = .true.
  logical :: UseGravity          = .true.
  logical :: UseIonDrag          = .true.
  logical :: UseViscosity        = .true.
  logical :: UseCoriolis         = .true.
  logical :: UseHorAdvection     = .true.
  logical :: UseVerAdvection     = .true.
  logical :: UseNeutralFriction  = .true.

  logical :: UseIonPressureGradient = .true.
  logical :: UseIonGravity          = .true.
  logical :: UseNeutralDrag         = .true.
  logical :: UseExB                 = .true.

  logical :: UseSolarHeating   = .true.
  logical :: UseJouleHeating   = .true.
  logical :: UseAuroralHeating = .true.
  logical :: UseNOCooling      = .true.
  logical :: UseOCooling       = .true.
  logical :: UseConduction     = .true.
  logical :: UseDiffusion      = .false.
  logical :: UseVerAdvectionT  = .true.

  real :: KappaTemp0 = 5.6e-4
  real :: KappaEddy  = 2.5
  real :: EddyDiffusionCoef = 0.0
  real :: EddyDiffusionPressure0 = 0.0
  real :: EddyDiffusionPressure1 = 0.0
  real :: Kappa1DCorrectionFactor = 45.0
  real :: Kappa1DCorrectionPower  = 1.75
  logical :: UseKappa1DCorrection = .true.

  logical :: UseIonChemistry     = .true.
  logical :: UseNeutralChemistry = .true.
  logical :: UseIonAdvection     = .true.

  !\
  ! Methods for completing chemistry
  !/

  integer, parameter :: cSteadyStateChemType_ = 1
  integer, parameter :: cImplicitChemType_    = 2
  integer, parameter :: cSubCycleChemType_    = 3
  integer, parameter :: nChemTypes_       = 3

  character (len=100), dimension(nChemTypes_) :: sChemType 

  character (len=100)                         :: sInputIonChemType
  character (len=100)                         :: sInputNeutralChemType
  integer :: iInputIonChemType, iInputNeutralChemType

  real :: LogNS0(nSpecies)

contains

  ! -------------------------------------------------------------
  ! set_strings
  !     Here we initialize all of the strings which need to be
  !     initialized in the code.
  ! -------------------------------------------------------------

  subroutine set_strings

    sChemType(cSteadyStateChemType_) = "steady"
    sChemType(cImplicitChemType_)    = "implicit"
    sChemType(cSubCycleChemType_)    = "subcycle"

  end subroutine set_strings

  ! -------------------------------------------------------------
  ! set_defaults
  ! -------------------------------------------------------------

  subroutine set_defaults

    use ModTime

    call set_strings

    sInputIonChemType     = sChemType(cSubCycleChemType_)
    sInputNeutralChemType = sChemType(cSubCycleChemType_)
    iInputIonChemType     = cSubCycleChemType_
    iInputNeutralChemType = cSubCycleChemType_

    dTAurora = 120.0

    tSimulation = 0.0

    iStartTime(1) = 1999
    iStartTime(2) =    3
    iStartTime(3) =   21
    iStartTime(4) =    0
    iStartTime(5) =    0
    iStartTime(6) =    0
    iStartTime(7) =    0

    iTimeArray = iStartTime

    call time_int_to_real(iStartTime, CurrentTime)

    nOutputTypes = 0
    DtPlot       = -1.0
    DtRestart    = -1.0

    ! Initiate Neutral Variables From MSIS90

    f107a = 150.0
    f107 = 150.0

    iModelSolar = 0

    UseApex   = .false.
    DoRestart = .false.

  end subroutine set_defaults

end module ModInputs
