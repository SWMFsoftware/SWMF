!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModPhysics

  use ModNumConst, ONLY: cDegToRad
  use ModConst
  use ModMain, ONLY: body2_
  use ModVarIndexes, ONLY: nVar, nFluid, IonFirst_, SpeciesFirst_, SpeciesLast_
  implicit none
  save

  ! default gamma value
  real, parameter:: Gamma0 = 5./3.

  ! adiabatic index (gamma) and derived values for fluids
  real:: Gamma_I(nFluid)          = 5./3.
  real:: GammaMinus1_I(nFluid)    = 5./3. - 1.0
  real:: InvGammaMinus1_I(nFluid) = 1.0/(5./3. - 1.0)

  ! adiabatic index (gamma) and derived values for the first/total fluid
  real :: Gamma          = 5./3.
  real :: GammaMinus1    = 5./3. - 1.0
  real :: InvGammaMinus1 = 1.0/(5./3. - 1.0)

  ! adiabatic index (gamma) and derived values for electrons
  real :: GammaElectron          = 5./3.
  real :: GammaElectronMinus1    = 5./3. - 1.0
  real :: InvGammaElectronMinus1 = 1.0/(5./3. - 1.0)

  ! gamma of the waves
  real:: GammaWave = 1.5

  ! electron charge in normalized units (actually proton charge/mass)
  real:: ElectronCharge

  ! plasma parameters
  real:: AverageIonCharge         = 1.0
  real:: ElectronTemperatureRatio = 0.0
  real:: ElectronPressureRatio    = 0.0
  real:: PePerPtotal              = 0.0
  real:: IonMassPerCharge         = 1.0

  ! Ion charge for multi-species.
  real:: ChargeSpecies_I(SpeciesFirst_:SpeciesLast_) = 1.0

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
  real :: MonopoleStrength = 0.0, MonopoleStrengthSi = 0.0 ! the monopole B0
  real :: Bdp, DipoleStrengthSi=0.0            ! the dipole moment of B0
  real :: Qqp(3,3)  =0.0                       ! the quadrupole moment of B0
  real :: Oop(3,3,3)=0.0                       ! the octupole moment of B0

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

  ! Density ratio of major and minor ions/neutrals (e.g. in the solar wind)
  real :: LowDensityRatio = 0.0001

  ! Minimum threshold for electron pressure and extra internal energy
  real :: PeMinSi = -1.1e5, PeMinDim = -1.0, PeMin
  real :: TeMinDim = -1.0, TeMin
  real :: ExtraEintMinSi = 0.0, ExtraEintMin

  ! Minimum threshold for MHD density and pressure
  real :: RhoMinDim_I(nFluid) = -1.0, RhoMin_I(nFluid)
  real :: pMinDim_I(nFluid)   = -1.0, pMin_I(nFluid)
  real :: TMinDim_I(nFluid)   = -1.0, TMin_I(nFluid)

  ! Boundary pressure for subsonic outflow
  logical:: UseOutflowPressure = .false.
  real :: pOutflowSi = -1.0, pOutflow = -1.0


  ! Relaxation time for anisotropic pressure
  logical :: UseConstantTau = .false.
  real    :: TauInstabilitySi = -1.0, TauInstability
  real    :: TauGlobalSi = -1.0, TauGlobal
  
  !\
  ! General variables for the second body
  !/
  real :: rPlanetDimBody2=0.0, rBody2=0.0, rCurrentsBody2=0.0
  real :: xBody2=0.0, yBody2=0.0, zBody2=0.0
  real :: PhaseBody2=0.0, DistanceBody2=0.0
  real :: RhoDimBody2=0.0, tDimBody2=0.0, RhoBody2=0.0, pBody2=0.0
  real :: gBody2=0.0
  logical:: UseBody2Orbit = .false.
  real :: OrbitPeriod=0.0

  ! Variables for two-state shock tube problems
  logical :: UseShockTube = .false.
  real :: ShockLeftState_V(nVar)=0.0, ShockRightState_V(nVar)=0.0
  real :: ShockPosition = 0.0, ShockSlope = 0.0

  ! State for the boundary conditions
  real,dimension(nVar,body2_:6):: FaceState_VI, CellState_VI

  !\
  ! Position of Earth for Dynamic AMR purposes, including
  ! a parameter that specifies the opening angle of the ray
  !/
  logical :: UseSunEarth=.false.
  real    :: xEarth=0.0, yEarth=0.0, zEarth=0.0, InvD2Ray

  !\
  ! Units for normalization of variables
  !/
  character (len=20) :: TypeNormalization="PLANETARY"

  !\
  ! Type of units used for I/O (input params, log files, plot files, etc.)
  !/
  character (len=20) :: TypeIoUnit = "PLANETARY"

  ! Named indexes for I/O variable units
  integer, parameter :: nIoUnit = 18

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
  integer, parameter :: UnitMass_        = 16
  integer, parameter :: UnitCharge_      = 17
  integer, parameter :: UnitUnity_       = 18

  ! Conversion between units: e.g. VarSi = VarNo*No2Si_V(UnitVar_)
  ! The following should always be true: No2Si_V*Si2Io_V = No2Io_V
  real, dimension(nIoUnit) :: &
       Io2Si_V, Si2Io_V, Io2No_V, No2Io_V, Si2No_V, No2Si_V

  ! Mapping between state array indices and unit conversion array indices
  integer, dimension(nVar) :: iUnitCons_V
  integer, dimension(nVar) :: iUnitPrim_V

  character (len=20), dimension(nIoUnit) :: &
       NameIdlUnit_V, NameTecUnit_V, NameSiUnit_V

  ! Some strange logical used in calc_heat_flux
  logical :: UseDefaultUnits = .false.

  ! Use Stellar parameters
  logical :: UseStar = .false.
  real :: RadiusStar=1.0,MassStar=1.0,RotationPeriodStar=25.38

  ! Number and indexes of vector variables in State_VGB
  integer:: nVectorVar = 0
  integer, allocatable:: iVectorVar_I(:) ! Index of first components

end module ModPhysics
