!^CFG COPYRIGHT UM
module ModPhysics

  use ModNumConst, ONLY: cDegToRad
  use ModConst
  use ModMain, ONLY:body2_,Top_
  use ModVarIndexes, ONLY: nVar, nFluid, IonFirst_
  implicit none
  save

  ! adiabatic index (gamma) and derived values
  real, parameter:: Gamma0 = 5./3.               ! default value
  real:: g_half = Gamma0/2.0
  real:: g = Gamma0, inv_g = 1.0/Gamma0, inv_gm1 = 1.0/(Gamma0 - 1.0)
  real:: gm1 = Gamma0 - 1.0, gm2 = Gamma0 - 2.0, gp1 = Gamma0 + 1.0

  ! electron charge in normalized units
  real:: ElectronCharge

  ! plasma parameters
  real:: AverageIonCharge         = 1.0
  real:: ElectronTemperatureRatio = 0.0
  real:: IonMassPerCharge         = 1.0

  ! thermal/total energy ratio limits for correctP
  real    :: Pratio_lo=0.01, Pratio_hi=0.1

  ! speed of light, inverse, square, inverse of square, boris correction
  real :: Clight, InvClight, C2light, Inv_C2light
  real :: Boris_Clight_Factor = 1.0

  ! normalized radiation constant (Erad = cRadiationNo*Trad**4)
  real :: cRadiationNo

  ! the dipole moment for body2
  real :: BdpBody2_D(3)=0.0, BdpDimBody2_D(3)=0.0

  !\
  ! Dipole and multipole expansion terms NOW ONLY IH SHOULD USE THESE
  !/
  real :: Bdp,Bdpx,Bdpy,Bdpz,DipoleStrengthSi =0.0 ! the dipole moment of B0
  real, dimension(1:3,1:3) :: Qqp  =0.0   ! the quadrupole moment of B0
  real, dimension(1:3,1:3,1:3) :: Oop=0   ! the octupole moment of B0

  real :: THETAtilt=0.0, &                ! tilt angle of magnetic axis
       SinThetaTilt=0.0, CosThetaTilt=1.0 ! NOW ONLY IH SHOULD USE THIS !!!

  !\
  ! The following are some notes on how to pick the Q's.  I have used the
  ! cartesian version of the quadrupole magnetic potential because it was
  ! the easiest to differentiate in our coordinate system.  The two ways 
  ! to write the potential are as follows (in SI - see Jackson pp.136-8):
  !
  !    spherical:       phi = (1/5)*mu * sum(m=-2..2) qlm * Ylm/r^3
  !    cartesian:       phi = 1/(8*pi)*mu * sum(i,j=1..3) Qqpij * xi*xj/r^5
  !
  !  the coefficients are related to each other as follows (note that I am 
  !  using the relations in Jackson.  He uses Gaussian units.  I would guess
  !  that these relations are the same for both systems but I have not checked
  !  them.  They are most usefull in getting the theta and phi dependance that 
  !  you want and are not really used to do any type of converting):
  !
  !               q22 = (1/12)*sqrt(15/2/pi)*(Qqp11-i*Qqp12 - Qqp22)
  !               q21 = -(1/3)*sqrt(15/8/pi)*(Qqp13-i*Qqp23)
  !               q20 = (1/2)*sqrt(5/4/pi)*Qqp33
  !
  !  Note that Qqp is TRACELESS.  Also note that it is symmetric.  The q's have
  !  the following property:
  !
  !                  ql(-m) = (-1)^m  *  complex_conjugate(qlm)   
  !
  !/


  !\
  ! Far field solar wind solution variables.
  !/
  real :: SW_T_dim=0.0, &
       SW_rho=0.0, SW_rho_dim=0.0, &
       SW_n=0.0,   SW_n_dim=0.0  , &
       SW_p=0.0  , SW_p_dim=0.0  , &
       SW_Ux=0.0 , SW_Ux_dim=0.0 , &
       SW_Uy=0.0 , SW_Uy_dim=0.0 , &
       SW_Uz=0.0 , SW_Uz_dim=0.0 , &
       SW_Bx=0.0 , SW_Bx_dim=0.0 , &
       SW_By=0.0 , SW_By_dim=0.0 , &
       SW_Bz=0.0 , SW_Bz_dim=0.0

  !\
  ! General Body parameters
  !/
  character (len=2) :: NamePlanetRadius = 'R ' !can be 'km' if there is no body
  real :: rPlanetSi=0.0, rBody=0.0, rCurrents=0.0
  real :: gBody=0.0
  real :: RotPeriodSi=0.0, OmegaBody=0.0

  ! The dimensional quantities are given for individual ion and neutral fluids
  real, dimension(IonFirst_:nFluid) :: &
       BodyNDim_I = 1.0, BodyTDim_I = 1.0, &
       PolarNDim_I= 1.0, PolarTDim_I= 1.0, PolarUDim_I = 0.0

  ! The normalized quantities include the total ion fluid (if present)
  real, dimension(nFluid) :: &
       BodyRho_I = 1.0, BodyP_I = 1.0, &
       PolarRho_I= 1.0, PolarP_I= 1.0, PolarU_I=0.0

  ! Polar boundary conditions are applied above this latitude only
  real :: PolarLatitude = 0.0, PolarTheta = 90.0*cDegToRad

  ! Density ratio of major and minor ions/neutrals (e.g. in the solar wind)
  real :: LowDensityRatio = 0.0001

  ! Minimum threshold for electron pressure 
  real :: PeMinSi = -1.1e5, PeMin

  !\
  ! General variables for the second body
  !/
  real :: rPlanetDimBody2=0.0, rBody2=0.0, rCurrentsBody2=0.0
  real :: xBody2=0.0, yBody2=0.0, zBody2=0.0
  real :: xBody2init=0.0, yBody2init=0.0
  real :: RhoDimBody2=0.0, tDimBody2=0.0, RhoBody2=0.0, pBody2=0.0
  real :: gBody2=0.0
  real :: OrbitPeriod=0.0

  ! Variables for two-state shock tube problems
  logical :: UseShockTube = .false.
  real :: ShockLeftState_V(nVar)=0.0, ShockRightState_V(nVar)=0.0
  real :: ShockPosition = 0.0, ShockSlope = 0.0

  ! State for the boundary conditions
  real,dimension(nVar,body2_:Top_):: FaceState_VI, CellState_VI

  !\
  ! Position of Earth for Dynamic AMR purposes, including
  ! a parameter that specifies the opening angle of the ray
  !/
  logical :: UseSunEarth=.false.
  real    :: xEarth=0.0, yEarth=0.0, zEarth=0.0, InvD2Ray

  !\
  ! Heat conduction parameters
  !/
  real :: Kappa0Heat=0.0, ExponentHeat=0.0

  !\
  ! Units for normalization of variables
  !/
  character (len=20) :: TypeNormalization="PLANETARY"

  !\
  ! Type of units used for I/O (input params, log files, plot files, etc.)
  !/
  character (len=20) :: TypeIoUnit = "PLANETARY"

  ! Named indexes for I/O variable units
  integer, parameter :: nIoUnit = 15

  integer, parameter :: UnitX_           = 1
  integer, parameter :: UnitU_           = 2
  integer, parameter :: UnitRho_         = 3
  integer, parameter :: UnitT_           = 4
  integer, parameter :: UnitN_           = 5
  integer, parameter :: UnitP_           = 6
  integer, parameter :: UnitB_           = 7
  integer, parameter :: UnitRhoU_        = 8
  integer, parameter :: UnitEnergyDens_  = 9
  integer, parameter :: UnitPoynting_    = 10
  integer, parameter :: UnitJ_           = 11
  integer, parameter :: UnitElectric_    = 12
  integer, parameter :: UnitTemperature_ = 13
  integer, parameter :: UnitDivB_        = 14
  integer, parameter :: UnitAngle_       = 15

  ! Conversion between units: e.g. VarSi = VarNo*No2Si_V(UnitVar_)
  ! The following should always be true: No2Si_V*Si2Io_V = No2Io_V
  real, dimension(nIoUnit) :: &
       Io2Si_V, Si2Io_V, Io2No_V, No2Io_V, Si2No_V, No2Si_V

  character (len=20), dimension(nIoUnit) :: &
       NameIdlUnit_V, NameTecUnit_V, NameSiUnit_V

  ! Some strange logical used in calc_heat_flux
  logical :: UseDefaultUnits = .false.

end module ModPhysics
