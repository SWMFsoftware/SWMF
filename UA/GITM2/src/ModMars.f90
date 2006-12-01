
module ModPlanet

  use ModConstants

  implicit none

  integer, parameter :: nSpecies = 4
  integer, parameter :: iCO2_    = 1
  integer, parameter :: iCO_     = 2
  integer, parameter :: iO_      = 3
  integer, parameter :: iN2_     = 4

  integer, parameter :: iN_  =  5
  integer, parameter :: iNO_ =  6
  integer, parameter :: iAr_ =  7
  integer, parameter :: iO2_ =  8
  integer, parameter :: iHe_ =  9
  integer, parameter :: iH_  =  10
  integer, parameter :: nSpeciesTotal = iH_

  integer, parameter  :: iO2P_  = 1
  integer, parameter  :: iCO2P_ = 2
  integer, parameter  :: iNOP_  = 3
  integer, parameter  :: iOP_   = 4
  integer, parameter  :: ie_    = 5
  integer, parameter  :: nIons  = ie_
  integer, parameter  :: nIonsAdvect = 4

  real :: Mass(nSpeciesTotal), MassI(nIons)

  real :: Vibration(nSpeciesTotal)

  ! When you want to program in emissions, you can use this...
  integer, parameter :: nEmissions = 10

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

!  real :: KappaTemp0 = 2.22e-4

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! These are totally wrong !!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! These are the numerical coefficients in Table 1 in m^2 instead of cm^2
  real, parameter, dimension(4, 4) :: Diff0 = 1.0e4 * reshape( (/ &
       ! 0      02     N2      N     NO
       !---------------------------------+
       0.00,  0.260, 0.260, 0.300, &            ! O
       0.26,  0.000, 0.181, 0.220, &            ! O2
       0.26,  0.181, 0.000, 0.220, &            ! N2
       0.30,  0.220, 0.220, 0.000 /), (/4,4/) )  ! N

  ! These are the exponents
  real, parameter, dimension(4, 4) :: DiffExp = reshape( (/ &
       ! 0      02     N2
       !---------------------------------+
       0.00,  0.75,  0.75, 0.75, &             ! O
       0.75,  0.00,  0.75, 0.75, &             ! O2
       0.75,  0.75,  0.00, 0.75, &             ! N2
       0.75,  0.75,  0.75, 0.00 /), (/4,4/) )  ! N

contains

  subroutine init_planet

    use ModTime

    integer :: iTime(7)

    Mass(iO_)    = 15.9994 * AMU
    Mass(iCO_)   = 12.011 * AMU + Mass(iO_)
    Mass(iCO2_)  = Mass(iCO_) + Mass(iO_)
    Mass(iN_)    = 14.00674 * AMU
    Mass(iN2_)   = Mass(iN_) * 2
    Mass(iNO_)   = Mass(iN_) + Mass(iO_)

    Mass(iO2_)   = 2 * Mass(iO_)
    Mass(iAr_)   = 39.948 * AMU * 2
    Mass(iHe_)   = 4.0026 * AMU * 2
    Mass(iH_)    = 1.0079 * AMU * 2

    Vibration(iCO2_)  = 9.0  ! Is this right???
    Vibration(iCO_)   = 7.0
    Vibration(iO_)    = 5.0
    Vibration(iN2_)   = 7.0

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
