!^CFG COPYRIGHT UM
!\
! set_physics_constants set normalizations, physical constants, 
! module (GM/IH/SC) dependent physical parameters
!/
subroutine set_physics_constants
  use ModProcMH
  use ModMain
  use ModPhysics
  use CON_axes,   ONLY: get_axes
  use CON_planet, ONLY: get_planet, RadiusPlanet, MassPlanet, OmegaPlanet
  use ModVarIndexes

  implicit none

  real :: MassBodySi
  real :: MassBody2Si                 !^CFG IF SECONDBODY
  real :: cosTheta, sinTheta, cosPhi, sinPhi, & 
       xx, yy, zz

  integer :: i, iBoundary

  logical :: oktest, oktest_me

  call set_oktest('set_physics_constants',oktest, oktest_me)

  !\
  ! Load body rotation rates, masses, and radii by module (GM/IH/SC)
  !/
  select case(NameThisComp)
  case('GM')
     call get_planet( &
          RadiusPlanetOut   = rPlanetSi, &
          MassPlanetOut     = MassBodySi, &
          RotationPeriodOut = RotPeriodSi)
  case("SC","IH")
     rPlanetSi = rSun
     MassBodySi = mSun
     RotPeriodSi = RotationPeriodSun
  end select
 
  ! Note for GM  !!! BATSRUS's OmegaBody is siderial (relative to the Sun)
  ! and it is DIFFERENT from SWMF's inertial OmegaPlanet defined in CON_planet !!!
  OmegaBody = cTwoPi/RotPeriodSi

  ! Second body mass is set to zero by default   !^CFG IF SECONDBODY
  MassBody2Si = 0.0                                !^CFG IF SECONDBODY

  !\
  ! Call set_units, which set the quantities for converting from
  ! normalized to  dimensional quantities and vice versa.  Also
  ! sets input and output strings for the conversion quantities
  !/
  call set_units

  if(oktest .and. iProc==0) then
     write(*,'(4a15)')'No2Io_V','NameIdlUnit_V','No2Si_V','NameSiUnit_V'
     do i=1, nIoUnit
        write(*,'(es15.6," ",a6,"        ",es15.6," ",a6)') &
             No2Io_V(i), NameIdlUnit_V(i), &
             No2Si_V(i), NameSiUnit_V(i)
     end do
  end if

  !\
  ! set the (corrected) speed of light and get normalization
  !/
  Clight      = Boris_Clight_Factor * cLightSpeed * Si2No_V(UnitU_)
  C2light     = cLIGHT**2
  InvClight   = cOne/cLight
  Inv_C2light = cOne/c2LIGHT

  !\
  ! Convert rotation, and gravity to non-dimensional values
  !/
  ! if the rotation period is less than 1 second then you made
  ! a mistake - the period is to fast
  if (abs(RotPeriodSi) > 1.) then
     OmegaBody = OmegaBody * (1.0/Si2No_V(UnitT_))
  else
     if(UseRotatingFrame)then
        write(*,*) "--------------------------------------------------"
        write(*,*) "WARNING in set_physics:                           "
        write(*,*) "You have set UseRotatingFrame = ",UseRotatingFrame
        write(*,*) "but RotPeriodSi in seconds= ",RotPeriodSi
        write(*,*) "This is too fast! Setting OmegaBody=0.0           "
        write(*,*) "--------------------------------------------------"
     end if
     OmegaBody = 0.0
  end if

  ! Note: The mass of the body is in SI units
  Gbody  = -cGravitation*MassBodySi*(Si2No_V(UnitU_)**2 * Si2No_V(UnitX_))
  !^CFG IF SECONDBODY BEGIN
  GBody2 = -cGravitation*MassBody2Si*(Si2No_V(UnitU_)**2 * Si2No_V(UnitX_))
  !^CFG END SECONDBODY

  !\
  ! Normalize solar wind values. Note: the solarwind is in I/O units
  !/
  SW_n   = SW_n_dim*Io2No_V(UnitN_)
  SW_rho = SW_n*No2Si_V(UnitN_)*AverageIonMass*cProtonMass*Si2No_V(UnitRho_)
  SW_p   = SW_rho * SW_T_dim*Io2No_V(UnitTemperature_)
  SW_Ux  = SW_Ux_dim*Io2No_V(UnitU_)
  SW_Uy  = SW_Uy_dim*Io2No_V(UnitU_)
  SW_Uz  = SW_Uz_dim*Io2No_V(UnitU_)
  SW_Bx  = SW_Bx_dim*Io2No_V(UnitB_)
  SW_By  = SW_By_dim*Io2No_V(UnitB_)
  SW_Bz  = SW_Bz_dim*Io2No_V(UnitB_)

  ! These are useful for printing out values
  SW_rho_dim = SW_rho*No2Io_V(UnitRho_)
  SW_p_dim   = SW_p*No2Io_V(UnitP_)

  Body_rho= Body_rho_dim * Io2No_V(UnitRho_)
  Body_p  = Body_rho * Body_T_dim*Io2No_V(UnitTemperature_)
  !^CFG IF SECONDBODY BEGIN
  RhoBody2= RhoDimBody2 * Io2No_V(UnitRho_)
  pBody2  = RhoBody2 * TDimBody2*Io2No_V(UnitTemperature_)
  !^CFG END SECONDBODY

  !Here the arrays of the FACE VALUE are formed
  !Initialization
  do iBoundary=body2_,Top_
     FaceState_VI(:,iBoundary)=DefaultState_V(1:nVar)
  end do

  !For bodies:

  FaceState_VI(rho_,body1_)=Body_rho
  FaceState_VI(P_,body1_)=Body_p
  
  !The following part of the code is sensitive to a particular physical
  !model. It should be modified in adding/deleting the physical effects 
  !or features

  FaceState_VI(rho_,body2_)=RhoBody2                !^CFG IF SECONDBODY
  FaceState_VI(P_,body2_)=pBody2                    !^CFG IF SECONDBODY

  
  !For Outer Boundaries
  do iBoundary=East_,Top_
     FaceState_VI(rho_, iBoundary) = SW_rho
     FaceState_VI(Ux_,  iBoundary) = SW_Ux
     FaceState_VI(Uy_,  iBoundary) = SW_Uy
     FaceState_VI(Uz_,  iBoundary) = SW_Uz
     FaceState_VI(Bx_,  iBoundary) = SW_Bx
     FaceState_VI(By_,  iBoundary) = SW_By
     FaceState_VI(Bz_,  iBoundary) = SW_Bz
     FaceState_VI(P_,   iBoundary) = SW_p
  end do

  !Cell State is used for filling the ghostcells
  CellState_VI=FaceState_VI

  do iBoundary=body2_,Top_  
     CellState_VI(rhoUx_:rhoUz_,iBoundary) = &
          FaceState_VI(Ux_:Uz_,iBoundary)*FaceState_VI(rho_,iBoundary)
  end do

  !\
  ! Now do the magnetic field stuff
  !/
  ! The *ThetaTilt and Bdp variables are not needed for GM except
  ! for reporting them in write_progress.

  ! Nondimensionalize dipole strength.
  if(NameThisComp == 'GM') then
     call get_axes(Time_Simulation, MagAxisTiltGsmOut = ThetaTilt)
     call get_planet(DipoleStrengthOut = DipoleStrengthSi)
  end if
  Bdp  = DipoleStrengthSi*Si2No_V(UnitB_)

  BdpBody2_D = BdpDimBody2_D*Io2No_V(UnitB_)              !^CFG IF SECONDBODY

  ! Compute dipole tilt variables
  if(NameThisComp=='IH')then
     ! For IH ThetaTilt should be set with the #HELIODIPOLETILT command
     ! But how is it going to rotate ?
     CosThetaTilt = cos(ThetaTilt)
     SinThetaTilt = sin(ThetaTilt)
  end if

  ! by default quadrupole and octupole terms are zero
  Qqp(1:3,1:3) = 0.00
  Oop(1:3,1:3,1:3) = 0.00

end subroutine set_physics_constants

!===========================================================================

subroutine set_units

  use ModProcMH, ONLY:iProc
  use ModMain
  use ModPhysics
  use ModVarIndexes
  use ModUser, ONLY: user_io_units, user_normalization
  implicit none

  character (len=*), parameter :: NameSub="set_units"

  logical :: oktest, oktest_me
  !-----------------------------------------------------------------------
  call set_oktest(NameSub,oktest, oktest_me)

  !\
  ! set variables for converting from normalized to SI units and back:
  !
  ! Normalized*UnitSi = Si, Si/UnitSi = Normalized
  !
  ! There are three independent variables: distance(x), velocity(u)
  ! and density(rho).  All others are built from these three.
  !/
  select case(TypeNormalization)
  case("PLANETARY")
     ! rPlanet, rPlanet/sec, amu/cm^3
     No2Si_V(UnitX_)   = rPlanetSi
     No2Si_V(UnitU_)   = rPlanetSi
     No2Si_V(UnitRho_) = 1000000*cProtonMass
  case("SOLARWIND")
     ! rPlanet, SW sound speed, SW density in amu/cm^3
     No2Si_V(UnitX_)   = rPlanetSi                             
     No2Si_V(UnitU_)   = sqrt(g*cBoltzmann*SW_T_dim/cProtonMass)
     No2Si_V(UnitRho_) = 1000000*cProtonMass*AverageIonMass*SW_n_dim
  case("NONE", "READ")
     ! Already set in MH_set_parameters
  case("USER")
     call user_normalization
  case default
     call stop_mpi(NameSub//' ERROR: unknown TypeNormalization='// &
          trim(TypeNormalization))
  end select

  !\
  ! set other normalizing SI variables from the independent ones
  !/
  No2Si_V(UnitT_)          = No2Si_V(UnitX_)/No2Si_V(UnitU_)         ! s
  No2Si_V(UnitN_)          = No2Si_V(UnitRho_)/cProtonMass           ! #/m^3
  No2Si_V(UnitP_)          = No2Si_V(UnitRho_)*No2Si_V(UnitU_)**2    ! Pa
  No2Si_V(UnitB_)          = No2Si_V(UnitU_) &
       *sqrt(cMu*No2Si_V(UnitRho_))                                  ! T
  No2Si_V(UnitRhoU_)       = No2Si_V(UnitRho_)*No2Si_V(UnitU_)       ! kg/m^2/s
  No2Si_V(UnitEnergyDens_) = No2Si_V(UnitP_)                         ! J/m^3
  No2Si_V(UnitPoynting_)   = No2Si_V(UnitEnergyDens_)*No2Si_V(UnitU_)! J/m^2/s
  No2Si_V(UnitJ_)          = No2Si_V(UnitB_)/(No2Si_V(UnitX_)*cMu)   ! A/m^2
  No2Si_V(UnitElectric_)   = No2Si_V(UnitU_)*No2Si_V(UnitB_)         ! V/m
  No2Si_V(UnitTemperature_)= (No2Si_V(UnitP_)/No2Si_V(UnitRho_)) &
       *(cProtonMass/cBoltzmann)                                     ! K 
  No2Si_V(UnitDivB_)       = No2Si_V(UnitB_)/No2Si_V(UnitX_)         ! T/m
  No2Si_V(UnitAngle_)      = 1.0                                     ! radian

  !\
  ! Set inverse conversion SI -> normalized
  !/
  Si2No_V = 1.0/No2Si_V
  
  !\
  ! set variables to go from Input/Output units to SI units:
  !
  ! Io*Io2Si_V = Si
  !
  ! Note that the input/output units are not necessarily consistent, e.g.
  ! units of distance divided by units of time does not necessarily
  ! coincide with the units of velocity.
  ! 
  ! Also load the unit name strings for IDL and TEC output
  !/

  ! As a default use SI units, so below only the differences need to be set
  Io2Si_V = 1.0
  No2Io_V = No2Si_V

  !\
  ! set string variables used for writing Tecplot output
  !/
  NameTecUnit_V(UnitX_)           = '[m]'            
  NameTecUnit_V(UnitU_)           = '[m/s]'          
  NameTecUnit_V(UnitRho_)         = '[kg/m^3]'
  NameTecUnit_V(UnitT_)           = '[s]'             
  NameTecUnit_V(UnitN_)           = '[m^-^3]'        
  NameTecUnit_V(UnitP_)           = '[Pa]'           
  NameTecUnit_V(UnitB_)           = '[T]'            
  NameTecUnit_V(UnitRhoU_)        = '[kg m^-^2 s^-^2]'
  NameTecUnit_V(UnitEnergydens_)  = '[J/m^3]'             
  NameTecUnit_V(UnitPoynting_)    = '[J m^-^2 s^-^1]'
  NameTecUnit_V(UnitJ_)           = '[A/m^2]'       
  NameTecUnit_V(UnitElectric_)    = '[V/m]'          
  NameTecUnit_V(UnitTemperature_) = '[K]'             
  NameTecUnit_V(UnitDivB_)        = '[T/m]'           
  NameTecUnit_V(UnitAngle_)       = '[rad]'
  !\
  ! set string variables used for writing IDL output
  !/
  NameIdlUnit_V(UnitX_)           = 'm'            
  NameIdlUnit_V(UnitRho_)         = 'kg/m3'
  NameIdlUnit_V(UnitU_)           = 'm/s'          
  NameIdlUnit_V(UnitT_)           = 's'             
  NameIdlUnit_V(UnitN_)           = '/m3'        
  NameIdlUnit_V(UnitP_)           = 'Pa'           
  NameIdlUnit_V(UnitB_)           = 'T'            
  NameIdlUnit_V(UnitRhoU_)        = 'kg/m2s2'
  NameIdlUnit_V(UnitEnergyDens_)  = 'J/m3'           
  NameIdlUnit_V(UnitPoynting_)    = 'J/m2s'
  NameIdlUnit_V(UnitJ_)           = 'A/m2'       
  NameIdlUnit_V(UnitElectric_)    = 'V/m'          
  NameIdlUnit_V(Unittemperature_) = 'K'             
  NameIdlUnit_V(UnitDivB_)        = 'T/m'           
  NameIdlUnit_V(UnitAngle_)       = 'rad'

  ! Store SI unit names for writing out variables in SI units
  NameSiUnit_V = NameIdlUnit_V

  select case(TypeIoUnit)
  case("SI")
     ! Already set above
  case("PLANETARY")
     Io2Si_V(UnitX_)        = rPlanetSi                       ! planetary radii
     Io2Si_V(UnitRho_)      = 1.0E6*cProtonMass               ! Mp/cm^3
     Io2Si_V(UnitN_)        = 1.0E6                           ! #/cm^3
     Io2Si_V(UnitU_)        = 1.0E3                           ! km/s
     Io2Si_V(UnitP_)        = 1.0E-9                          ! nPa
     Io2Si_V(UnitB_)        = 1.0E-9                          ! nT
     Io2Si_V(UnitJ_)        = 1.0E-6                          ! microA/m^2
     Io2Si_V(UnitElectric_) = 1.0E-3                          ! mV/m
     Io2Si_V(UnitDivB_)     = 1.0E-9/rPlanetSi                ! nT/R_planet
     Io2Si_V(UnitAngle_)    = cRadToDeg                       ! degrees
     !\
     ! set string variables used for writing output - TECPLOT
     !/
     NameTecUnit_V(UnitX_)           = '[R]'            
     NameTecUnit_V(UnitRho_)         = '[amu/cm^3]'
     NameTecUnit_V(UnitU_)           = '[km/s]'          
     NameTecUnit_V(UnitN_)           = '[cm^-^3]'        
     NameTecUnit_V(UnitP_)           = '[nPa]'           
     NameTecUnit_V(UnitB_)           = '[nT]'            
     NameTecUnit_V(UnitJ_)           = '[`mA/m^2]'       
     NameTecUnit_V(UnitElectric_)    = '[mV/m]'          
     NameTecUnit_V(UnitDivB_)        = '[nT/R]'
     NameTecUnit_V(UnitAngle_)       = '[deg]'

     !\
     ! set string variables used for writing output - IDL
     !/
     NameIdlUnit_V(UnitX_)           = 'R'
     NameIdlUnit_V(UnitRho_)         = 'Mp/cc'
     NameIdlUnit_V(UnitU_)           = 'km/s'
     NameIdlUnit_V(UnitN_)           = '/cc'
     NameIdlUnit_V(UnitP_)           = 'nPa'           
     NameIdlUnit_V(UnitB_)           = 'nT'
     NameIdlUnit_V(UnitJ_)           = 'uA/m2'
     NameIdlUnit_V(UnitElectric_)    = 'mV/m'
     NameIdlUnit_V(UnitDivB_)        = 'nT/R'
     NameIdlUnit_V(UnitAngle_)       = 'deg'

  case("HELIOSPHERIC")
     Io2Si_V(UnitX_)           = rPlanetSi                 ! R
     Io2Si_V(UnitRho_)         = 1.0E+3                    ! g/cm^3
     Io2Si_V(UnitN_)           = 1.0E+6                    ! #/cm^3
     Io2Si_V(UnitU_)           = 1.0E+3                    ! km/s
     Io2Si_V(UnitP_)           = 1.0E-1                    ! dyne/cm^2
     Io2Si_V(UnitB_)           = 1.0E-4                    ! Gauss
     Io2Si_V(UnitRhoU_)        = 1.0E+1                    ! g/cm^2/s
     Io2Si_V(UnitEnergydens_)  = 1.0E-1                    ! erg/cm^3
     Io2Si_V(UnitJ_)           = 1.0E-6                    ! uA/m^2
     Io2Si_V(UnitDivB_)        = 1.0E-2                    ! Gauss/cm
     Io2Si_V(UnitAngle_)       = cRadToDeg                 ! degrees
     !\
     ! set string variables used for writing output - TECPLOT
     !/
     NameTecUnit_V(UnitX_)           = '[R]'
     NameTecUnit_V(UnitRho_)         = '[g/cm^3]'
     NameTecUnit_V(UnitU_)           = '[km/s]'
     NameTecUnit_V(UnitN_)           = '[amu/cm^3]'
     NameTecUnit_V(UnitP_)           = '[dyne/cm^2]'
     NameTecUnit_V(UnitB_)           = '[Gauss]'
     NameTecUnit_V(UnitRhoU_)        = '[g/cm^2/s]'
     NameTecUnit_V(UnitEnergyDens_)  = '[erg/cm^3]'
     NameTecUnit_V(UnitJ_)           = '[`mA/m^2]'
     NameTecUnit_V(UnitDivB_)        = '[Gauss/cm]'
     NameTecUnit_V(UnitAngle_)       = '[deg]'
     !\
     ! set string variables used for writing output - IDL
     !/
     NameIdlUnit_V(UnitX_)           = 'R'
     NameIdlUnit_V(UnitRho_)         = 'g/cm3'
     NameIdlUnit_V(UnitU_)           = 'km/s'
     NameIdlUnit_V(UnitN_)           = 'mp/cc'
     NameIdlUnit_V(UnitP_)           = 'dyne/cm^2'
     NameIdlUnit_V(UnitB_)           = 'G'
     NameIdlUnit_V(UnitRhoU_)        = 'g/cm^2/s'
     NameIdlUnit_V(UnitEnergyDens_)  = 'erg/cm3'
     NameIdlUnit_V(UnitJ_)           = 'uA/m2'
     NameIdlUnit_V(UnitDivB_)        = 'G/cm'
     NameIdlUnit_V(UnitTemperature_) = 'K'

  case("NONE")
     ! I/O and normalized units are the same, so
     Io2Si_V = No2Si_V
     NameTecUnit_V = ''            
     NameIdlUnit_V = ''   

  case("USER")
     ! User method provides the conversion from I/O to SI units
     ! and Tecplot and IDL strings for all units differing from SI units.
     call user_io_units
  case default
     call stop_mpi(NameThisComp//': Unknown TypeIoUnit='//TypeIoUnit)
  end select

  ! Calculate the remaining unit conversions
  Si2Io_V = 1/Io2Si_V
  No2Io_V = No2Si_V*Si2Io_V
  Io2No_V = 1/No2Io_V

end subroutine set_units

!==============================================================================

subroutine init_mhd_variables

  ! Set default I/O units and unit names for the state variables 
  ! in MHD type equations

  use ModVarIndexes
  use ModPhysics

  implicit none

  integer :: iVar
  !--------------------------------------------------------------------------

  UnitUser_V(Rho_)          = No2Io_V(UnitRho_)
  UnitUser_V(RhoUx_:RhoUz_) = No2Io_V(UnitRhoU_)
  UnitUser_V(Bx_:Bz_)       = No2Io_V(UnitB_)
  UnitUser_V(p_)            = No2Io_V(UnitP_)
  UnitUser_V(Energy_)       = No2Io_V(UnitEnergyDens_)

  NameUnitUserIdl_V(rho_)          = NameTecUnit_V(UnitRho_)
  NameUnitUserIdl_V(RhoUx_:rhoUz_) = NameTecUnit_V(UnitRhoU_)
  NameUnitUserIdl_V(Bx_:Bz_)       = NameTecUnit_V(UnitB_)
  NameUnitUserIdl_V(P_)            = NameTecUnit_V(UnitP_)
  NameUnitUserIdl_V(Energy_)       = NameTecUnit_V(UnitEnergyDens_)

  NameUnitUserIdl_V(rho_)          = NameIdlUnit_V(UnitRho_)
  NameUnitUserIdl_V(RhoUx_:rhoUz_) = NameIdlUnit_V(UnitRhoU_)
  NameUnitUserIdl_V(Bx_:Bz_)       = NameIdlUnit_V(UnitB_)
  NameUnitUserIdl_V(P_)            = NameIdlUnit_V(UnitP_)
  NameUnitUserIdl_V(Energy_)       = NameIdlUnit_V(UnitEnergyDens_)

  ! By default the scalar advected variables are assumed to behave like density
  do iVar = ScalarFirst_, ScalarLast_
     UnitUser_V(iVar)        = No2Io_V(UnitRho_)
     NameUnitUserTec_V(iVar) = NameTecUnit_V(UnitRho_)
     NameUnitUserIdl_V(iVar) = NameIdlUnit_V(UnitRho_)
  end do

end subroutine init_mhd_variables
