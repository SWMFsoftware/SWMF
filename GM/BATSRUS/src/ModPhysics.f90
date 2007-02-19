!^CFG COPYRIGHT UM
module ModPhysics
  use ModConst
  use ModMain, ONLY:body2_,Top_
  use ModVarIndexes,ONLY:nVar
  implicit none
  save

  real :: g = cOne + cTwo/cThree,&
       gm1 = cTwo/cThree,&
       gm2 = -cOne/cThree,&
       gp1 = cTwo + cTwo/cThree,&
       inv_g = cThree/(cTwo + cThree),&
       inv_gm1 = cThree/cTwo,&
       g_half = cHalf + cOne/cThree       ! gamma and derived values

  ! plasma parameters
  real:: AverageIonMass           = 1.0   ! only used if not UseMultiSpecies 
  real:: AverageIonCharge         = 1.0
  real:: ElectronTemperatureRatio = 0.0

  ! thermal/total energy ratio limits for correctP  !^CFG IF PROJECTION
  real    :: Pratio_lo=0.01, Pratio_hi=0.1          !^CFG IF PROJECTION

  ! speed of light, inverse, square, inverse of square, boris correction
  real :: Clight, InvClight, C2light, Inv_C2light, Boris_Clight_Factor   

  ! the dipole moment for body2           !^CFG IF SECONDBODY
  real :: BdpBody2_D(3),BdpDimBody2_D(3)  !^CFG IF SECONDBODY

  !\
  ! Dipole and multipole expansion terms NOW ONLY IH SHOULD USE THESE
  !/
  real :: Bdp,Bdpx,Bdpy,Bdpz,Bdp_dim=0.0  ! the dipole moment of B0
  real, dimension(1:3,1:3) :: Qqp  =0.0   ! the quadrupole moment of B0
  real, dimension(1:3,1:3,1:3) :: Oop=0   ! the octupole moment of B0

  real :: THETAtilt=0.0, &                ! tilt angle of magnetic axis
       sinTHETAtilt,cosTHETAtilt          ! NOW ONLY IH SHOULD USE THIS !!!

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
  real ::      SW_T_dim=0.0  , &
       SW_a_dim=0.0  , &
       SW_rho=0.0, SW_rho_dim=0.0, &
       SW_p=0.0  , SW_p_dim=0.0  , &
       SW_Ux=0.0 , SW_Ux_dim=0.0 , &
       SW_Uy=0.0 , SW_Uy_dim=0.0 , &
       SW_Uz=0.0 , SW_Uz_dim=0.0 , &
       SW_Bx=0.0 , SW_Bx_dim=0.0 , &
       SW_By=0.0 , SW_By_dim=0.0 , &
       SW_Bz=0.0 , SW_Bz_dim=0.0 , &
       SW_B_factor=0.0

  !\
  ! General Body parameters
  !/
  real :: rPlanet_dim, rBody, rCurrents
  real :: Body_rho_dim, Body_T_dim, Body_rho, Body_p  
  real :: gBody
  real :: rot_period_dim, OMEGAbody 

  !^CFG IF SECONDBODY BEGIN
  !\
  ! General variables for the second body
  !/
  real :: rPlanetDimBody2, rBody2, rCurrentsBody2, xBody2, yBody2, zBody2
  real :: RhoDimBody2, tDimBody2, RhoBody2, pBody2  
  real :: gBody2
  !^CFG END SECONDBODY

  ! Variables for two-state shock tube problems
  logical :: UseShockTube = .false.
  real :: ShockLeftState_V(nVar), ShockRightState_V(nVar), ShockSlope

  ! State for the boundary conditions
  real,dimension(nVar,body2_:Top_):: FaceState_VI, CellState_VI

  !\
  ! Position of Earth for Dynamic AMR purposes, including
  ! a parameter that specifies the opening angle of the ray
  !/
  logical :: UseSunEarth=.false.
  real :: xEarth,yEarth,zEarth, InvD2Ray

  !\
  ! Heat conduction parameters
  !/
  real :: Kappa0Heat,ExponentHeat

  !\
  ! Resistivity parameters
  !/
  character (len=20) :: TypeResist
  real :: TimeInitRise,TimeConstLev
  real :: Eta0Resist,Alpha0Resist,yShiftResist
  real :: Eta0AnomResist,EtaAnomMaxResist,jCritResist

  !\
  ! Ilia's reconnection & test problems
  !/
  character (len=20) :: TypeProblemDiss
  real :: EpsilonDiss,DeltaDiss,ThetaDiss
  real :: RhoDifDiss,yShiftDiss
  real :: scaleHeightDiss,scaleFactorDiss,BZ0Diss
  !\
  ! Normalization Units
  !/
  logical :: UseDefaultUnits = .false.
  real :: Grav0Diss,Beta0Diss,Length0Diss
  real :: Time0Diss,Rho0Diss,Tem0Diss

  !\
  ! Type of Units the user is using - for doing I/O
  !/
  character (len=20) :: IoUnits

  !\
  ! Units for normalization
  !/
  ! used as follows:     rho_with_dimensions = State_VGB(rho_,...)*unit_rho
  !                      State_VGB(rho_,... = rho_with_dimensions/unit_rho
  ! unitUSER_x = units are assigned by user, this is used for input and ouput
  ! unitSI_x   = units are standard SI units - m, kg, s, T, K, V, A, ...

  character (len=20) :: TypeNormalization='SOLARWIND'

  real ::  &
       unitUSER_x, unitUSER_t, unitUSER_angle,        & ! time and space
       unitUSER_rho, unitUSER_n, unitUSER_U,          & ! primitive MHD quantities 
       unitUSER_p, unitUSER_B,                        & ! primitive MHD quantities
       unitUSER_rhoU,  unitUSER_energydens,           & ! conservative MHD quantities
       unitUSER_J, unitUSER_electric, unitUSER_DivB,  & ! derived quantities
       unitUSER_temperature,                          & ! derived quantities
       unitUSER_Poynting                                ! Poynting vector
  real ::  &
       unitSI_x, unitSI_t, unitSI_angle,              & ! time and space
       unitSI_rho, unitSI_n, unitSI_U,                & ! primitive MHD quantities 
       unitSI_p, unitSI_B,                            & ! primitive MHD quantities
       unitSI_rhoU,  unitSI_energydens,               & ! conservative MHD quantities
       unitSI_J, unitSI_electric, unitSI_DivB,        & ! derived quantities
       unitSI_temperature,                            & ! derived quantities
       unitSI_Poynting                                  ! Poynting vector

  ! String representation of the units above - used for output - IDL
  character (len=20) ::  &
       unitstr_IDL_x, unitstr_IDL_t, unitstr_IDL_angle, & ! time and space
       unitstr_IDL_rho, unitstr_IDL_n, unitstr_IDL_U,   & ! primitive MHD quantities
       unitstr_IDL_p, unitstr_IDL_B,                    & ! primitive MHD quantities
       unitstr_IDL_rhoU,  unitstr_IDL_energydens,       & ! conservative MHD quantities
       unitstr_IDL_J, unitstr_IDL_electric,             & ! derived quantities
       unitstr_IDL_DivB, unitstr_IDL_temperature,       & ! derived quantities
       unitstr_IDL_Poynting                               ! Poynting vector

  ! String representation of the units above - used for output - TEC
  character (len=20) ::  &
       unitstr_TEC_x, unitstr_TEC_t, unitstr_TEC_angle, & ! time and space
       unitstr_TEC_rho, unitstr_TEC_n, unitstr_TEC_U,   & ! primitive MHD quantities
       unitstr_TEC_p, unitstr_TEC_B,                    & ! primitive MHD quantities
       unitstr_TEC_rhoU,  unitstr_TEC_energydens,       & ! conservative MHD quantities
       unitstr_TEC_J, unitstr_TEC_electric,             & ! derived quantities
       unitstr_TEC_DivB, unitstr_TEC_temperature,       & ! derived quantities
       unitstr_TEC_Poynting                               ! Poynting vector

  ! IO unit type and string for standard variables and plots types.  These are
  ! used by the user to set the values.
  integer, parameter :: nIoUnit = 13
  real, dimension(nIoUnit) :: Si2User
  character(len=20), dimension(nIoUnit) :: IoUnitStr
  integer, parameter :: UnitX_           = 1
  integer, parameter :: UnitRho_         = 2
  integer, parameter :: UnitN_           = 3
  integer, parameter :: UnitU_           = 4
  integer, parameter :: UnitT_           = 5
  integer, parameter :: UnitP_           = 6
  integer, parameter :: UnitB_           = 7
  integer, parameter :: UnitRhoU_        = 8
  integer, parameter :: UnitEnergyDens_  = 9
  integer, parameter :: UnitPoynting_    = 10
  integer, parameter :: UnitJ_           = 11
  integer, parameter :: UnitElectric_    = 12
  integer, parameter :: UnitTemperature_ = 13


end module ModPhysics
