
module ModPlanet

  use ModConstants

  implicit none

  integer, parameter :: nSpecies = 4
  integer, parameter :: iCO2_    = 1
  integer, parameter :: iCO_     = 2
  integer, parameter :: iO_      = 3
  integer, parameter :: iN2_     = 4

  integer, parameter :: nSpeciesTotal = 8
  integer, parameter :: iAr_ =  5
  integer, parameter :: iO2_ =  6
  integer, parameter :: iHe_ =  7
  integer, parameter :: iH_  =  8

  integer, parameter  :: iO2P_  = 1
  integer, parameter  :: iCO2P_ = 2
  integer, parameter  :: iNOP_  = 3
  integer, parameter  :: iOP_   = 4
  integer, parameter  :: ie_    = 5
  integer, parameter  :: nIons  = ie_
  integer, parameter  :: nIonsAdvect = 4

  real :: Mass(nSpeciesTotal), MassI(nIons)

  real, parameter :: GC_Mars                = 3.73                    ! m/s^2
  real, parameter :: RP_Mars                = 88800.0                 ! seconds
  real, parameter :: R_Mars                 = 3388.25*1000.0          ! meters
  real, parameter :: DP_Mars                = 0.0

  real, parameter :: Gravitational_Constant = GC_Mars
  real, parameter :: Rotation_Period        = RP_Mars
  real, parameter :: RBody                  = R_Mars
  real, parameter :: DipoleStrength         = DP_Mars

  real, parameter :: OMEGABody              = 2.00*pi/Rotation_Period  ! rad/s

  real, parameter :: HoursPerDay = Rotation_Period / 3600.0
  real, parameter :: Tilt = 25.0

  ! This is the Vernal Equanox at Midnight (hopefully!!!)
  integer, parameter :: iVernalYear   = 1998
  integer, parameter :: iVernalMonth  =    7
  integer, parameter :: iVernalDay    =   14
  integer, parameter :: iVernalHour   =   16
  integer, parameter :: iVernalMinute =    0
  integer, parameter :: iVernalSecond =    0

  real, parameter :: SunOrbit_A = 1.52
  real, parameter :: SunOrbit_B = 0.04
  real, parameter :: SunOrbit_C = 0.15
  real, parameter :: SunOrbit_D = 0.00
  real, parameter :: SunOrbit_E = 0.00

  real, parameter :: DaysPerYear = 670.0
  real, parameter :: SecondsPerYear = DaysPerYear * Rotation_Period

  logical :: IsEarth = .false.
  character (len=10) :: cPlanet = "Mars"

  real :: KappaTemp0 = 2.22e-4

contains

  subroutine init_planet

    use ModTime

    integer :: iTime(7)

    Mass(iO_)    = 15.9994 * AMU
    Mass(iCO_)   = 12.011 * AMU + Mass(iO_)
    Mass(iCO2_)  = Mass(iCO_) + Mass(iO_)
    Mass(iN2_)   = 14.00674 * AMU * 2

    Mass(iO2_)   = 2 * Mass(iO_)
    Mass(iAr_)   = 39.948 * AMU * 2
    Mass(iHe_)   = 4.0026 * AMU * 2
    Mass(iH_)    = 1.0079 * AMU * 2

    MassI(iOP_)   = Mass(iO_)
    MassI(iNOP_)  = Mass(iO_) + Mass(iN2_)/2.0
    MassI(iCO2P_) = Mass(iCO2_)
    MassI(iO2P_)  = Mass(iO2_)
    MassI(ie_) = Mass_Electron

    itime = 0
    itime(1) = iVernalYear
    itime(2) = iVernalMonth
    itime(3) = iVernalDay
    itime(4) = iVernalHour
    itime(5) = iVernalMinute
    itime(6) = iVernalSecond
    call time_int_to_real(itime, VernalTime)

  end subroutine init_planet

end module ModPlanet
