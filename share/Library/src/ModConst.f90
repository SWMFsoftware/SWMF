!^CFG COPYRIGHT UM
Module ModConst
  use ModNumConst
  use ModKind

  implicit none

  save

  !\  
  ! Physical and astronomical constants
  !/
  !\
  ! Physical constants
  !/

  real,parameter::  cBoltzmann  = 13.807/cE24
  ! Boltzmann constant  1.3807E-23 [J/K]

  real, parameter :: cProtonMass = 1.6726/cE27
  ! Proton mass  1.6726E-27[kg]

  real, parameter :: cElectronMass = 9.1094/cE30/cE1
  ! Electron mass  9.1094E-31[kg]

  real, parameter :: cLightSpeed =  0.29979*cE9
  ! Speed of light (true) 2.9979E+08 [m/s]

  real, parameter :: cMu = cPi*0.4/cE6
  ! Vacuum permeability Pi*4.0E-7[H/m]

  real, parameter :: cEps = cOne/cMu/(cLightSpeed**2)
  ! Vacuum permittivity 8.8542E-12[F/m]
  
  real,parameter :: cVacuumResistivity=cLightSpeed*cMu
  ! Vacuum resistance 120*cPi [Ohm]

  real, parameter :: cElectronCharge  = 0.16022/cE18
  ! Fundamental charge 1.6022E-19 [C (coulomb)]

  real, parameter :: cAvogadro =0.6022045*cE24
  ! Number of particles per mole 6.022045E+23[mole^-1]

  real, parameter :: cGravitation = 66.726/cE12
  ! Gravitation constant 6.6726E-11  (NRL 1994)

  real(Real8_), parameter :: cSecondPerYear   = 31536000.0
  real(Real8_), parameter :: cSecondPerDay    =    86400.0
  real(Real8_), parameter :: cSecondPerHour   =     3600.0
  real(Real8_), parameter :: cSecondPerMinute =       60.0

  !\
  ! Astronomical constants : r = radius, m = mass
  !/

  !The Earth:
  real,parameter :: rEarth = 6378.00*cThousand   ! = 6378.00E03 [ m]
  real,parameter :: mEarth =  5.976*cE24         ! = 5.976E+24 [kg]
  real,parameter :: RotationPeriodEarth = 24.0 * 3600.0

  real,parameter :: IonoHeightEarth = 110000.0   ! = 110,000 [m]

  ! The orbital period belongs to the TROPICAL YEAR, which is
  ! relative to the vernal equinox which is slowly moving 
  ! due to the precession of the Earth's rotation axis.
  real,parameter :: OrbitalPeriodEarth  = 365.24218967 * 24.0 * 3600.0

  ! The rotational angular velocity is relative to an inertial frame
  real,parameter :: OmegaEarth = &
       cTwoPi/RotationPeriodEarth + cTwoPi/OrbitalPeriodEarth

  real,parameter :: TiltEarth = 23.5 * cDegToRad

  ! Reference equinox time taken from 
  ! http://aa.usno.navy.mil/data/docs/EarthSeasons.html
  integer,parameter:: &
       iYearEquinoxEarth  = 2000, &
       iMonthEquinoxEarth =    3, &
       iDayEquinoxEarth   =   20, &
       iHourEquinoxEarth  =    7, &
       iMinuteEquinoxEarth=   35, &
       iSecondEquinoxEarth=    0
  real(Real8_),parameter :: FracSecondEquinoxEarth = 0.0

  ! The angle between the zero meridian and the eqinox direction at 
  ! equinox time. For Earth this can be calculated from the time of day.
  ! For other planets there is no analogous method to calculate this angle.

  real, parameter :: AngleEquinoxEarth = cTwoPi * &
       ( iHourEquinoxEarth * 3600 + iMinuteEquinoxEarth * 60 &
       + iSecondEquinoxEarth + FracSecondEquinoxEarth) / (24 * 3600)

  real,parameter :: DipoleStrengthEarth = -31100.0 * 1.0e-9  ! in Tesla !!!
  real,parameter :: bAxisThetaEarth    = 11.0    * cDegToRad
  real,parameter :: bAxisPhiEarth      = 289.1   * cDegToRad

  ! Sun

  real,parameter :: Rsun = 0.696*cE9                ! = 6.96E+08 [m]   
  real,parameter :: mSun = 1.99*cE30                ! = 1.99E+30 [kg]

  real,parameter :: tSunRot = 26.0    !Rotation period in days
  real,parameter :: RotationPeriodSun = Tsunrot * RotationPeriodEarth

  ! Astronomical Unit (1AU)

  real,parameter :: cAU = 1.4959787000000000000000*cE9*cE2

  ! Saturn

  real,parameter :: rSaturn = 60268.00*cE3          ! = 60268.00E+03 [m]
  real,parameter :: mSaturn = 0.5685*cE27           ! = 5.685E26 [kg]
  real,parameter :: RotationPeriodSaturn = 10.5 * 3600.0

  real,parameter :: rTitan = 2575.00*cE3            ! = 2575.00E3
  real,parameter :: rTitan_Orbit = 1.222*cE9        ! = 1.224E9 
  real,parameter :: OmegaTitan_Orbit = 4.56/cE6     ! = 4.56E-6

  ! Jupiter

  real,parameter :: rJupiter = 71492.00*cE3         !=71492.00E03 [m]
  real,parameter :: mJupiter = 1.8990*cE27          !=1.899E27[kg]
  real,parameter :: RotationPeriodJupiter = 9.925 * 3600.0

  real,parameter :: rIo =1821.00*cE3                !=1821.00E3
  real,parameter :: OmegaIo = 41.105/cE6            !=4.1105E-5

  ! Venus

  real,parameter :: rVenus = 6052.0*cE3
  ! Radius of the Venus 6052.00E03 [ m]   

  real,parameter :: mVenus =  4.865*cE24
  ! Mass of the Venus 4.865E+24  [kg]

  real,parameter :: RotationPeriodVenus = 5834.4

  ! Mars

  real,parameter :: rMars =  3396.00*cE3
  ! Radius of the Mars 3396.00E03 [ m]   

  real,parameter :: mMars =  0.6436*cE24
  ! Mass of the Mars 6.436*cE23  [kg]

  real,parameter :: RotationPeriodMars = 1.026*24*3600.0
  ! rotation period in hours

end module ModConst
