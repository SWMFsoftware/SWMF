
module ModPlanet

  use ModConstants
  use ModSizeGITM, only: nAlts

  implicit none

! Modified (01/18/07) : SWB :   Aij, s-exponents for mutual diffusion
! Majors (4):  COntrol the Pressures Gradients and winds
  integer, parameter :: nSpecies = 4
  integer, parameter :: iCO2_    = 1
  integer, parameter :: iCO_     = 2
  integer, parameter :: iO_      = 3
  integer, parameter :: iN2_     = 4

! Minors (6) : Ride on background of mean winds derived from majors
  integer, parameter :: iN_  =  5
  integer, parameter :: iNO_ =  6
  integer, parameter :: iAr_ =  7
  integer, parameter :: iO2_ =  8
  integer, parameter :: iHe_ =  9
  integer, parameter :: iH_  =  10
  integer, parameter :: nSpeciesTotal = iH_

! Major Ions (4):  Most Important to MWACM code
  integer, parameter  :: iO2P_  = 1
  integer, parameter  :: iCO2P_ = 2
  integer, parameter  :: iNOP_  = 3
  integer, parameter  :: iOP_   = 4
  integer, parameter  :: ie_    = 5
  integer, parameter  :: nIons  = ie_
  integer, parameter  :: nIonsAdvect = 4

  character (len=20) :: cSpecies(nSpeciesTotal)
  character (len=20) :: cIons(nIons)

  real :: Mass(nSpeciesTotal), MassI(nIons)

  real :: Vibration(nSpeciesTotal)

  ! When you want to program in emissions, you can use this...
  integer, parameter :: nEmissions = 10

  real, parameter :: GC_Mars                = 3.73                    ! m/s^2
  real, parameter :: RP_Mars                = 88775.0                 ! seconds
  real, parameter :: R_Mars                 = 3388.25*1000.0          ! meters
  real, parameter :: DP_Mars                = 0.0

  real, parameter :: Gravitational_Constant = GC_Mars
  real, parameter :: Rotation_Period        = RP_Mars
  real, parameter :: RBody                  = R_Mars
  real, parameter :: DipoleStrength         = DP_Mars

  real, parameter :: OMEGABody              = 2.00*pi/Rotation_Period  ! rad/s

  real, parameter :: HoursPerDay = Rotation_Period / 3600.0
  real, parameter :: Tilt = 25.19

  ! This is the Vernal Equinox at Midnight (Ls = 0!!!)
  ! Earth-Mars clocks are set from this epoch at vernal equinox
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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! These are Modified for Mars by SWB: 1/18/07
  ! -- Most source state Dij = Dji (check)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! These are Aij coefficients from B&K (1973) formulation: Aij*1.0E+17
  ! Use Jared Bell's Titan GITM formulation for Mars GITM

  real, parameter, dimension(4, 4) :: Diff0 = 1.0e17 * reshape( (/ &
! integer, parameter :: iCO2_    = 1
! integer, parameter :: iCO_     = 2
! integer, parameter :: iO_      = 3
! integer, parameter :: iN2_     = 4
     !---------------------------------+
     ! i=C02      CO      O     N2  
     !---------------------------------+
!       0.0000, 0.7762, 0.2119,  0.6580, &            ! CO2
!       0.4940, 0.0000, 0.9466,  0.9280, &            ! CO
!       0.7703, 0.9481, 0.0000,  0.9690, &            ! O
!       0.6580, 0.9280, 0.9690,  0.0000 /), (/4,4/) ) ! N2
!
       0.0000, 0.7762, 0.2119,  0.6580, &            ! CO2
       0.7762, 0.0000, 0.9466,  0.9280, &            ! CO
       0.2219, 0.9466, 0.0000,  0.9690, &            ! O
       0.6580, 0.9280, 0.9690,  0.0000 /), (/4,4/) ) ! N2

  ! These are s-exponents from B&K (1973) formulation: T**s
  real, parameter, dimension(4, 4) :: DiffExp = reshape( (/ &
       !---------------------------------+
       !i= CO2    CO      O      N2
       !---------------------------------+
       0.000,  0.750,  0.750, 0.752, &             ! CO2
       0.750,  0.000,  0.750, 0.710, &             ! CO
       0.750,  0.750,  0.000, 0.774, &             ! O
       0.752,  0.710,  0.774, 0.000 /), (/4,4/) )  ! N2

!     Arrays filled in init_radcool in data statements (np = 68)
  integer, parameter :: np=68
  real,dimension(np) :: pnbr,ef1,ef2,co2vmr,o3pvmr,n2covmr

  !! Stuff for initial conditions

  real , Dimension(-1:nAlts + 2) :: newalt
  real , Dimension(-1:nAlts + 2) :: InTemp
  real , Dimension(-1:nAlts + 2) :: IneTemp
  real , Dimension(-1:nAlts + 2) :: InITemp
  real , Dimension(-1:nAlts + 2,nSpeciesTotal) :: InNDensityS 
  real , Dimension(-1:nAlts + 2,nIons) :: InIDensityS

  real, parameter:: AltMinIono=100.0 ! in km

contains

  subroutine init_planet

    use ModTime
    use ModIoUnit, only : UnitTmp_

    integer :: iTime(7), iiAlt

!   Mass = AMU * mean molecular weight  (unlike TGCM codes)

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

    cSpecies(iO_)    = "O"
    cSpecies(iO2_)   = "O!D2!N"
    cSpecies(iN_)    = "N"
    cSpecies(iN2_)   = "N!D2!N"
    cSpecies(iCO_)   = "CO"
    cSpecies(iCO2_)  = "CO!D2!N"
    cSpecies(iNO_)   = "NO"
    cSpecies(iAr_)   = "Ar"
    cSpecies(iH_)    = "H"
    cSpecies(iHe_)   = "He"

    cIons(iO2P_)   = "O!D2!U+!N"
    cIons(iCO2P_)   = "CO!D2!U+!N"
    cIons(iNOP_)   = "NO!U+!N"
    cIons(iOP_)    = "O!U+!N"
    cIons(ie_)     = "e-"

    Vibration(iCO2_)  = 7.0  ! Corrected by Bougher (01/18/07)!!!!
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

    write(*,*) 'Reading in the Mars_input.txt'
    open(UNIT = UnitTmp_, FILE = 'UA/DataIn/Mars_input.txt', &
         STATUS='OLD', ACTION = 'READ')

111 FORMAT(F6.1,1X, F8.2,1X, F8.2,1X, F8.2,1X,   &  
         ES10.3,1X, ES10.3, 1X,  ES10.3, 1X, ES10.3, 1X, &
         ES10.3, 1X, ES10.3, 1X,  ES10.3)

    InNDensityS(:,:) = 1.0e+3
    InIDensityS(:,:) = 1.0e+3

    do iiAlt = -1,nAlts+2
       read(UnitTmp_,111) &
            newalt(iiAlt), &
            InTemp(iiAlt), &
            InITemp(iiAlt), &
            IneTemp(iiAlt), &
            !
            InNDensityS(iiAlt,iCO2_), &
            InNDensityS(iiAlt,iO2_), &
            InNDensityS(iiAlt,iCO_), &
            InNDensityS(iiAlt,iN2_), &
            !
            InNDensityS(iiAlt,iO_), &
            InNDensityS(iiAlt,iAr_), &
            
            InIDensityS(iiAlt,ie_)

    end do

    close(Unit = UnitTmp_)

!open(UNIT = 56, FILE = 'UA/DataIn/Densities.out', STATUS='NEW',ACTION = 'WRITE')
!
!  do iiAlt = -1,nAlts+2
!        write(56,111) &
!      	newalt(iiAlt), &
!       	InTemp(iiAlt), &
!       	InITemp(iiAlt), &
!       	IneTemp(iiAlt), &
!
!       	InNDensityS(iiAlt,iCO2_), &
!       	InNDensityS(iiAlt,iO2_), &
!       	InNDensityS(iiAlt,iCO_), &
!       	InNDensityS(iiAlt,iN2_), &
!!
!       	InNDensityS(iiAlt,iO_), &
!       	InNDensityS(iiAlt,iAr_), &
!
!       	InIDensityS(iiAlt,ie_)
!  enddo
!  close(Unit = 56)
!
!write(*,*) 'End Writing the Densities.output'



  end subroutine init_planet

end module ModPlanet
