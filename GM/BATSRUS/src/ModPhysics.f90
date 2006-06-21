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
  real ::      SW_T_dim  , &
       SW_a_dim  , &
       SW_rho, SW_rho_dim, &
       SW_p  , SW_p_dim  , &
       SW_Ux , SW_Ux_dim , &
       SW_Uy , SW_Uy_dim , &
       SW_Uz , SW_Uz_dim , &
       SW_Bx , SW_Bx_dim , &
       SW_By , SW_By_dim , &
       SW_Bz , SW_Bz_dim , &
       SW_B_factor

  real, dimension(0:1) :: &
       SW_rho_t,  &
       SW_p_t  ,  &
       SW_Ux_t ,  &
       SW_Uy_t ,  &
       SW_Uz_t ,  &
       SW_Bx_t ,  &
       SW_By_t ,  &
       SW_Bz_t ,  &
       SW_time_t
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

  ! Shocktube initial state values
  real, dimension(1:nVar) :: shock_Lstate, shock_Rstate
  real :: ShockSlope

  ! State for the boundary conditions
  real,dimension(nVar,body2_:Top_):: FaceState_VI, CellState_VI

  !\
  ! Heliosphere terms.
  !/
  real ::  PreSun, RhoSun, SspSun
  real ::  VelSun
  real ::  qSun, tHeat, rHeat, SigmaHeat

  !\
  ! CME and Arcade parameters
  !/
  logical:: UseFluxRope=.false.
  character(len=10) :: cme_type
  real :: cme_a, cme_r1, cme_r0, cme_a1, cme_alpha 
  real :: cme_rho1, cme_rho2, cme_B1_dim, cme_v_erupt
  real :: Rscl, RHOscl, rho1scl, rho2scl, SSPscl, Vscl 
  real :: B1scl, a1scl
  real :: ModulationRho,ModulationP
  real :: OrientationGL98, LatitudeGL98, LongitudeGL98
  real :: widthArc, phi0Arc, muArc
  real :: RhoArcDim, TArcDim, UzArcDim, BArcDim, ByArcDim
  real :: B0_scl, B0y_scl, Phtscl
  integer :: ExpArc

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
  ! Comet terms.
  !/
  real ::  kin,kin_in,mbar,Unr,Unr_in,ionization_rate, Qprod

  !\
  ! Units for normalization
  !/
  ! used as follows:     rho_with_dimensions = State_VGB(rho_,...)*unit_rho
  !                      State_VGB(rho_,... = rho_with_dimensions/unit_rho
  ! unitUSER_x = units are assigned by user, this is used for input and ouput mostly
  ! unitSI_x   = units are standard SI units - m, kg, s, T, K, V, A, ...
  real ::  &
       unitUSER_x, unitUSER_t, unitUSER_angle,        & ! time and space
       unitUSER_rho, unitUSER_n, unitUSER_U,          & ! primative MHD quantities 
       unitUSER_p, unitUSER_B,                        & ! primative MHD quantities
       unitUSER_rhoU,  unitUSER_energydens,           & ! conservative MHD quantities
       unitUSER_J, unitUSER_electric, unitUSER_DivB,  & ! derived quantities
       unitUSER_temperature,                          & ! derived quantities
       unitUSER_Poynting                                ! Poynting vector
  real ::  &
       unitSI_x, unitSI_t, unitSI_angle,              & ! time and space
       unitSI_rho, unitSI_n, unitSI_U,                & ! primative MHD quantities 
       unitSI_p, unitSI_B,                            & ! primative MHD quantities
       unitSI_rhoU,  unitSI_energydens,               & ! conservative MHD quantities
       unitSI_J, unitSI_electric, unitSI_DivB,        & ! derived quantities
       unitSI_temperature,                            & ! derived quantities
       unitSI_Poynting                                  ! Poynting vector

  ! String representation of the units above - used for output - IDL
  character (len=20) ::  &
       unitstr_IDL_x, unitstr_IDL_t, unitstr_IDL_angle, & ! time and space
       unitstr_IDL_rho, unitstr_IDL_n, unitstr_IDL_U,   & ! primative MHD quantities
       unitstr_IDL_p, unitstr_IDL_B,                    & ! primative MHD quantities
       unitstr_IDL_rhoU,  unitstr_IDL_energydens,       & ! conservative MHD quantities
       unitstr_IDL_J, unitstr_IDL_electric,             & ! derived quantities
       unitstr_IDL_DivB, unitstr_IDL_temperature,       & ! derived quantities
       unitstr_IDL_Poynting                               ! Poynting vector

  ! String representation of the units above - used for output - TEC
  character (len=20) ::  &
       unitstr_TEC_x, unitstr_TEC_t, unitstr_TEC_angle, & ! time and space
       unitstr_TEC_rho, unitstr_TEC_n, unitstr_TEC_U,   & ! primative MHD quantities
       unitstr_TEC_p, unitstr_TEC_B,                    & ! primative MHD quantities
       unitstr_TEC_rhoU,  unitstr_TEC_energydens,       & ! conservative MHD quantities
       unitstr_TEC_J, unitstr_TEC_electric,             & ! derived quantities
       unitstr_TEC_DivB, unitstr_TEC_temperature,       & ! derived quantities
       unitstr_TEC_Poynting                               ! Poynting vector

end module ModPhysics
