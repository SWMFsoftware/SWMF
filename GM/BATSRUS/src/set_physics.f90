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

  real :: Mbody_dim
  real :: MBody2Dim                 !^CFG IF SECONDBODY
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
          RadiusPlanetOut   = rPlanet_Dim, &
          MassPlanetOut     = mBody_Dim, &
          RotationPeriodOut = rot_period_dim)
  case("SC","IH")
     rPlanet_Dim = rSun
     mBody_Dim = mSun
     rot_period_dim = RotationPeriodSun
  end select
 
  ! Note for GM  !!! BATSRUS's OmegaBody is siderial (relative to the Sun)
  ! and it is DIFFERENT from SWMF's inertial OmegaPlanet defined in CON_planet !!!
  OmegaBody = cTwoPi/rot_period_dim

  ! Second body mass is set to zero by default   !^CFG IF SECONDBODY
  MBody2Dim = 0.0                                !^CFG IF SECONDBODY

  !\
  ! Call set_units, which set the quantities for converting from
  ! normalized to  dimensional quantities and vice versa.  Also
  ! sets input and output strings for the conversion quantities
  !/
  call set_units

  if(oktest .and. iProc==0) then
     write(*,'(E15.6,11X,E15.6)') unitUSER_x, unitSI_x
     write(*,'(E15.6,11X,E15.6)') unitUSER_t, unitSI_t
     write(*,'(E15.6,11X,E15.6)') unitUSER_angle, unitSI_angle       
     write(*,'(E15.6,11X,E15.6)') unitUSER_rho,  unitSI_rho
     write(*,'(E15.6,11X,E15.6)') unitUSER_n, unitSI_n 
     write(*,'(E15.6,11X,E15.6)') unitUSER_U, unitSI_U          
     write(*,'(E15.6,11X,E15.6)') unitUSER_p, unitSI_p 
     write(*,'(E15.6,11X,E15.6)') unitUSER_B, unitSI_B                        
     write(*,'(E15.6,11X,E15.6)') unitUSER_rhoU, unitSI_rhoU
     write(*,'(E15.6,11X,E15.6)') unitUSER_energydens, unitSI_energydens
     write(*,'(E15.6,11X,E15.6)') unitUSER_J, unitSI_J
     write(*,'(E15.6,11X,E15.6)') unitUSER_electric, unitSI_electric
     write(*,'(E15.6,11X,E15.6)') unitUSER_DivB, unitSI_DivB
  end if

  !\
  ! set the (corrected) speed of light and get normalization
  !/
  Clight      = Boris_Clight_Factor * cLightSpeed/unitSI_U
  C2light     = cLIGHT**2
  InvClight   = cOne/cLight
  Inv_C2light = cOne/c2LIGHT

  !\
  ! Convert rotation, and gravity to non-dimensional values
  !/
  ! if the rotation period is less than 1 second then you made
  ! a mistake - the period is to fast
  if (abs(rot_period_dim) > 1.) then
     OmegaBody = OmegaBody / (1.0/unitSI_t)
  else
     if(UseRotatingFrame)then
        write(*,*) "--------------------------------------------------"
        write(*,*) "WARNING in set_physics:                           "
        write(*,*) "You have set UseRotatingFrame = ",UseRotatingFrame
        write(*,*) "but rot_period_dim in hours= ",rot_period_dim
        write(*,*) "This is too fast! Setting OmegaBody=0.0           "
        write(*,*) "--------------------------------------------------"
     end if
     OmegaBody = 0.0
  end if

  Gbody  = -cGravitation*Mbody_dim*(1/unitSI_U**2/unitSI_x)
  GBody2 = -cGravitation*MBody2Dim*(1/unitSI_U**2/unitSI_x) !^CFG IF SECONDBODY

  !\
  ! Nondimensionalize dimensional SW values - 
  !/
  SW_a_dim = unitUSER_U
  SW_p_dim = unitUSER_p*inv_g
  SW_B_factor = unitUSER_B

!  SW_rho = SW_rho_dim/unitUSER_rho
  SW_rho = 1.0    
  SW_p   = inv_g
  SW_Ux  = SW_Ux_dim/unitUSER_U
  SW_Uy  = SW_Uy_dim/unitUSER_U
  SW_Uz  = SW_Uz_dim/unitUSER_U
  SW_Bx  = SW_Bx_dim/unitUSER_B
  SW_By  = SW_By_dim/unitUSER_B
  SW_Bz  = SW_Bz_dim/unitUSER_B

  Body_rho= Body_rho_dim/unitUSER_n
  Body_p  = cBoltzmann*Body_rho_dim*&
            1.0E6*Body_T_dim/unitSI_p
  RhoBody2= RhoDimBody2/unitUSER_n               !^CFG IF SECONDBODY
  pBody2  = cBoltzmann*RhoDimBody2*&             !^CFG IF SECONDBODY
            1.0E6*TDimBody2/unitSI_p             !^CFG IF SECONDBODY

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
     call get_planet(DipoleStrengthOut = Bdp_dim)
     Bdp      = Bdp_dim/unitSI_B 
  else
     Bdp      = Bdp_dim/unitUSER_B 
  end if
  BdpBody2_D = BdpDimBody2_D/unitUSER_B                 !^CFG IF SECONDBODY

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
  use ModUser, ONLY: user_io_units
  implicit none

  character (len=*), parameter :: NameSub="set_units"

  logical :: oktest, oktest_me
  !-----------------------------------------------------------------------
  call set_oktest(NameSub,oktest, oktest_me)

  select case(TypeNormalization)
  case("NONE")
     UnitSi_x   = 1.0
     UnitSi_u   = 1.0
     UnitSi_rho = 1.0
  case("PLANETARY")
     UnitSi_x   = rPlanet_Dim                         ! rPlanet
     UnitSi_u   = rPlanet_Dim                         ! rPlanet/sec
     UnitSi_rho = cProtonMass*cMillion                ! amu/cm^3
  case("SOLARWIND")
     UnitSi_x   = rPlanet_Dim                         ! rPlanet
     UnitSI_u   = sqrt(g*cBoltzmann*SW_T_dim/cProtonMass) ! SW sound speed
     UnitSI_rho = cProtonMass*(cMillion*SW_rho_dim)   ! SW density in amu/cm^3
  case("USER")
     ! Set in set_parameters
  case default
     call stop_mpi(NameSub//' ERROR: unknown TypeNormalization='// &
          trim(TypeNormalization))
  end select

  !\
  ! Load variables used for converting from dimensional to non-dimensional 
  ! units and back.  Also load the name variables for each of the units for
  ! use in writing output.
  !/

  unitSI_angle = 180/cPi   
  unitUSER_angle = unitSI_angle
  unitstr_TEC_angle = '[degree]'            
  unitstr_IDL_angle = 'deg'

  !\
  ! set variables for converting from unitless to SI units and back.
  ! unitless*unitSI = SI, SI/UnitSI = unitless
  !
  ! Note there are three independent variables, position(x), velocity(u)
  ! and density(rho).  All others are built from these three.
  !/
  unitSI_rho  = cProtonMass*(cMillion*SW_rho_dim)         ! kg/m^3
  unitSI_U    = sqrt(g*cBoltzmann*SW_T_dim/cProtonMass)   ! m/s

  !\
  ! set other normalizing SI variables from the independent ones
  !/
  unitSI_t           = unitSI_x/unitSI_U                     ! s
  unitSI_n           = unitSI_rho/cProtonMass                ! #/m^3
  unitSI_p           = unitSI_rho*unitSI_U**2                ! Pa
  unitSI_B           = unitSI_U*sqrt(cMu*unitSI_rho)         ! T
  unitSI_rhoU        = unitSI_rho*unitSI_U                   ! kg/m^2/s
  unitSI_energydens  = unitSI_p                              ! J/m^3
  unitSI_Poynting    = unitSI_energydens*unitSI_U            ! J/m^2/s
  unitSI_J           = unitSI_B/(unitSI_x*cMu)               ! A/m^2
  unitSI_electric    = unitSI_U*unitSI_B                     ! V/m
  unitSI_DivB        = unitSI_B/unitSI_x                     ! T/m
  ! set temperature - note that the below is only strictly true for a
  ! limited number of cases.  If the only ion is proton then 
  ! Tcode = Te+Ti.
  !
  ! If the ions are "heavy" (m = A * mp) then n above is not really a 
  ! number density but is a nucleon density (#nucleons/m^3).  In this 
  ! case the temperature in the code using unitSI_temperature is really
  ! Tcode = (Ti+Te)/A.
  !
  ! For the special case where A=2 we have Tcode = 1/2*(Te+Ti).  If
  ! we assume Ti=Te then Tcode=Ti=Te.
  !
  ! Also note that if the are heavy ions of mass #*mp then you could
  ! be off in temperature by as much as a factor of #.  
  ! There is no way around this in MHD.
  unitSI_temperature = (unitSI_p/unitSI_rho)*(cProtonMass/cBoltzmann) ! Kelvin 

  !\
  ! set variables to go from a limited set of units to normalized units
  ! and back
  ! unitless*unitUSER = User units,    User Units/unitUSER = unitless
  !
  ! Note that the user units are not consistant in a unit sense.
  ! 
  ! Also load the string variables associated with the USER variables - 
  ! note that they are loaded differently for IDL and TEC output
  !/
  select case(IoUnits)
  case("SI")
     unitUSER_x           = unitSI_X              ! m
     unitUSER_rho         = unitSI_rho            ! amu/m^3          
     unitUSER_U           = unitSI_U              ! #/m^3            
     unitUSER_t           = unitSI_t              ! m/s              
     unitUSER_n           = unitSI_n              ! s                           
     unitUSER_p           = unitSI_p              ! Pa                          
     unitUSER_B           = unitSI_B              ! T                          
     unitUSER_rhoU        = unitSI_rhoU           ! kg/m^2/s                    
     unitUSER_energydens  = unitSI_energydens     ! J/m^3                       
     unitUSER_Poynting    = unitSI_Poynting       ! J/m^2/s           
     unitUSER_J           = unitSI_J              ! A/m^2                      
     unitUSER_electric    = unitSI_electric       ! V/m                        
     unitUSER_DivB        = unitSI_DivB           ! T/m
     unitUSER_temperature = unitSI_temperature    ! dimensionless
     !\
     ! set string variables used for writing output - TECPLOT
     !/
     unitstr_TEC_x           = '[m]'            
     unitstr_TEC_rho         = '[amu/m^3]'       
     unitstr_TEC_U           = '[m/s]'          
     unitstr_TEC_t           = '[s]'             
     unitstr_TEC_n           = '[m^-^3]'        
     unitstr_TEC_p           = '[Pa]'           
     unitstr_TEC_B           = '[T]'            
     unitstr_TEC_rhoU        = '[kg m^-^2 s^-^2]'
     unitstr_TEC_energydens  = '[J/m^3]'             
     unitstr_TEC_Poynting    = '[J m^-^2 s^-^1]'
     unitstr_TEC_J           = '[A/m^2]'       
     unitstr_TEC_electric    = '[V/m]'          
     unitstr_TEC_DivB        = '[T/m]'           
     unitstr_TEC_temperature = '[K]'             
     !\
     ! set string variables used for writing output - IDL
     !/
     unitstr_IDL_x           = 'm'            
     unitstr_IDL_rho         = 'amu/m3'       
     unitstr_IDL_U           = 'm/s'          
     unitstr_IDL_t           = 's'             
     unitstr_IDL_n           = '/m3'        
     unitstr_IDL_p           = 'Pa'           
     unitstr_IDL_B           = 'T'            
     unitstr_IDL_rhoU        = 'kg/m2s2'
     unitstr_IDL_energydens  = 'J/m3'           
     unitstr_IDL_Poynting    = 'J/m^2s'
     unitstr_IDL_J           = 'A/m2'       
     unitstr_IDL_electric    = 'V/m'          
     unitstr_IDL_DivB        = 'T/m'           
     unitstr_IDL_temperature = 'K'             

  case("PLANETARY")
     unitUSER_x           = unitSI_x/unitSI_x           ! planetary radii
     unitUSER_rho         = 1.0E-6*unitSI_n             ! amu/cm^3
     unitUSER_n           = 1.0E-6*unitSI_n             ! #/cm^3
     unitUSER_U           = 1.0E-3*unitSI_U             ! km/s
     unitUSER_t           = unitSI_t                    ! s
     unitUSER_p           = 1.0E+9*unitSI_p             ! nPa
     unitUSER_B           = 1.0E+9*unitSI_B             ! nT
     unitUSER_rhoU        = unitSI_rhoU                 ! kg/m^2/s
     unitUSER_energydens  = unitSI_energydens           ! J/m^3 
     unitUSER_Poynting    = unitSI_Poynting             ! J/m^2/s
     unitUSER_J           = 1.0E+6*unitSI_J             ! uA/m^2
     unitUSER_electric    = 1.0E+3*unitSI_electric      ! mV/m
     unitUSER_DivB        = 1.0E+9*unitSI_DivB*unitSI_x ! nT/planetary radii
     unitUSER_temperature = unitSI_temperature          ! Kelvin 
     !\
     ! set string variables used for writing output - TECPLOT
     !/
     unitstr_TEC_x           = '[R]'            
     unitstr_TEC_rho         = '[amu/cm^3]'       
     unitstr_TEC_U           = '[km/s]'          
     unitstr_TEC_t           = '[s]'             
     unitstr_TEC_n           = '[cm^-^3]'        
     unitstr_TEC_p           = '[nPa]'           
     unitstr_TEC_B           = '[nT]'            
     unitstr_TEC_rhoU        = '[kg m^-^2 s^-^2]'
     unitstr_TEC_energydens  = '[J/m^3]'             
     unitstr_TEC_Poynting    = '[J m^-^2 s^-^1]'
     unitstr_TEC_J           = '[`mA/m^2]'       
     unitstr_TEC_electric    = '[mV/m]'          
     unitstr_TEC_DivB        = '[nT/R]'           
     unitstr_TEC_temperature = '[K]'             
     !\
     ! set string variables used for writing output - IDL
     !/
     unitstr_IDL_x           = 'R'            
     unitstr_IDL_rho         = 'amu/cm3'       
     unitstr_IDL_U           = 'km/s'          
     unitstr_IDL_t           = 's'             
     unitstr_IDL_n           = '/cc'        
     unitstr_IDL_p           = 'nPa'           
     unitstr_IDL_B           = 'nT'            
     unitstr_IDL_rhoU        = 'kg/m2s2'
     unitstr_IDL_energydens  = 'J/m3'           
     unitstr_IDL_Poynting    = 'J/m^2s'
     unitstr_IDL_J           = 'uA/m2'       
     unitstr_IDL_electric    = 'mV/m'          
     unitstr_IDL_DivB        = 'nT/R'           
     unitstr_IDL_temperature = 'K'             

  case("HELIOSPHERIC")
     unitUSER_x           = unitSI_x/unitSI_x                ! R
     unitUSER_rho         = 1.0e-3*unitSI_rho                ! g/cm^3
     unitUSER_n           = 1.0E-6*unitSI_n                  ! #/cm^3
     unitUSER_U           = 1.0E-3*unitSI_U                  ! km/s
     unitUSER_t           = unitSI_t                         ! s
     unitUSER_p           = 1.0E+1*unitSI_p                  ! dyne/cm^2
     unitUSER_B           = 1.0E+4*unitSI_B                  ! Gauss
     unitUSER_rhoU        = 1.0E-1*unitSI_rhoU               ! g/cm^2/s
     unitUSER_Poynting    = unitSI_Poynting                  ! J/m^2/s
     unitUSER_energydens  = 1.0E+1*unitSI_energydens         ! erg/cm^3
     unitUSER_J           = 1.0E+6*unitSI_J                  ! uA/m^2
     unitUSER_electric    = unitSI_electric                  ! V/m
     unitUSER_DivB        = 1.0E+2*unitSI_DivB*unitSI_x      ! Gauss/cm
     unitUSER_temperature = unitSI_temperature               ! Kelvin
     !\
     ! set string variables used for writing output - TECPLOT
     !/
     unitstr_TEC_x           = 'R'
     unitstr_TEC_rho         = 'g/cm3'
     unitstr_TEC_U           = 'km/s'
     unitstr_TEC_t           = 's'
     unitstr_TEC_n           = 'amu/cm3'
     unitstr_TEC_p           = 'dyne/cm^2'
     unitstr_TEC_B           = 'Gauss'
     unitstr_TEC_rhoU        = 'g/cm^2/s'
     unitstr_TEC_Poynting    = '[J m^-^2 s^-^1]'
     unitstr_TEC_energydens  = 'erg/cm3'
     unitstr_TEC_J           = 'uA/m2'
     unitstr_TEC_electric    = 'V/m'
     unitstr_TEC_DivB        = 'Gauss/cm'
     unitstr_TEC_temperature = 'K'
     !\
     ! set string variables used for writing output - IDL
     !/
     unitstr_IDL_x           = 'R'
     unitstr_IDL_rho         = 'g/cm3'
     unitstr_IDL_U           = 'km/s'
     unitstr_IDL_t           = 's'
     unitstr_IDL_n           = 'amu/cm3'
     unitstr_IDL_p           = 'dyne/cm^2'
     unitstr_IDL_B           = 'Gauss'
     unitstr_IDL_rhoU        = 'g/cm^2/s'
     unitstr_IDL_Poynting    = 'J/m^2s'
     unitstr_IDL_energydens  = 'erg/cm3'
     unitstr_IDL_J           = 'uA/m2'
     unitstr_IDL_electric    = 'V/m'
     unitstr_IDL_DivB        = 'Gauss/cm'
     unitstr_IDL_temperature = 'K'

  case("NONE")
     unitUSER_x           = 1.0
     unitUSER_rho         = 1.0
     unitUSER_U           = 1.0
     unitUSER_t           = 1.0
     unitUSER_n           = 1.0
     unitUSER_p           = 1.0
     unitUSER_B           = 1.0
     unitUSER_rhoU        = 1.0
     unitUSER_energydens  = 1.0
     unitUSER_Poynting    = 1.0
     unitUSER_J           = 1.0
     unitUSER_electric    = 1.0
     unitUSER_DivB        = 1.0
     unitUSER_temperature = 1.0
     !\
     ! set string variables used for writing output - TECPLOT
     !/
     unitstr_TEC_x           = ''            
     unitstr_TEC_rho         = ''
     unitstr_TEC_U           = ''
     unitstr_TEC_t           = ''
     unitstr_TEC_n           = ''
     unitstr_TEC_p           = ''
     unitstr_TEC_B           = ''
     unitstr_TEC_rhoU        = ''
     unitstr_TEC_energydens  = ''
     unitstr_TEC_Poynting    = ''
     unitstr_TEC_J           = ''
     unitstr_TEC_electric    = ''
     unitstr_TEC_DivB        = ''
     unitstr_TEC_temperature = ''            
     !\
     ! set string variables used for writing output - IDL
     !/
     unitstr_IDL_x           = ''   
     unitstr_IDL_rho         = ''   
     unitstr_IDL_U           = ''   
     unitstr_IDL_t           = ''    
     unitstr_IDL_n           = ''  
     unitstr_IDL_p           = ''   
     unitstr_IDL_B           = ''   
     unitstr_IDL_rhoU        = ''
     unitstr_IDL_energydens  = ''    
     unitstr_IDL_Poynting    = ''
     unitstr_IDL_J           = ''  
     unitstr_IDL_electric    = ''   
     unitstr_IDL_DivB        = ''    
     unitstr_IDL_temperature = ''           

  case("USER")
     call user_io_units
     ! Users provide the conversion from SI units to their IO 
     ! units.  Since we don't use the units this way, we have to 
     ! multiply by unitSI_ to set the unitUSER variables.
     unitUSER_x           = Si2User(UnitX_)           *unitSI_x  
     unitUSER_rho         = Si2User(UnitRho_)         *unitSI_rho          
     unitUSER_n           = Si2User(UnitN_)           *unitSI_U            
     unitUSER_U           = Si2User(UnitU_)           *unitSI_t            
     unitUSER_t           = Si2User(UnitT_)           *unitSI_n            
     unitUSER_p           = Si2User(UnitP_)           *unitSI_p            
     unitUSER_B           = Si2User(UnitB_)           *unitSI_B            
     unitUSER_rhoU        = Si2User(UnitRhoU_)        *unitSI_rhoU         
     unitUSER_energydens  = Si2User(UnitEnergyDens_)  *unitSI_energydens   
     unitUSER_Poynting    = Si2User(UnitPoynting_)    *unitSI_Poynting     
     unitUSER_J           = Si2User(UnitJ_)           *unitSI_J               
     unitUSER_electric    = Si2User(UnitElectric_)    *unitSI_electric
     unitUSER_temperature = Si2User(UnitTemperature_) *unitSI_temperature  
     unitUSER_DivB        = unitSI_DivB  ! Users don't need to know this
     
     unitstr_TEC_x           = IoUnitStr(UnitX_)           
     unitstr_TEC_rho         = IoUnitStr(UnitRho_)         
     unitstr_TEC_n           = IoUnitStr(UnitN_)           
     unitstr_TEC_U           = IoUnitStr(UnitU_)           
     unitstr_TEC_t           = IoUnitStr(UnitT_)           
     unitstr_TEC_p           = IoUnitStr(UnitP_)           
     unitstr_TEC_B           = IoUnitStr(UnitB_)           
     unitstr_TEC_rhoU        = IoUnitStr(UnitRhoU_)        
     unitstr_TEC_energydens  = IoUnitStr(UnitEnergyDens_)  
     unitstr_TEC_Poynting    = IoUnitStr(UnitPoynting_)    
     unitstr_TEC_J           = IoUnitStr(UnitJ_)           
     unitstr_TEC_electric    = IoUnitStr(UnitElectric_)    
     unitstr_TEC_temperature = IoUnitStr(UnitTemperature_) 
     unitstr_TEC_DivB        = '[T/m]'

     unitstr_IDL_x           = IoUnitStr(UnitX_)           
     unitstr_IDL_rho         = IoUnitStr(UnitRho_)         
     unitstr_IDL_n           = IoUnitStr(UnitN_)           
     unitstr_IDL_U           = IoUnitStr(UnitU_)           
     unitstr_IDL_t           = IoUnitStr(UnitT_)           
     unitstr_IDL_p           = IoUnitStr(UnitP_)           
     unitstr_IDL_B           = IoUnitStr(UnitB_)           
     unitstr_IDL_rhoU        = IoUnitStr(UnitRhoU_)        
     unitstr_IDL_energydens  = IoUnitStr(UnitEnergyDens_)  
     unitstr_IDL_Poynting    = IoUnitStr(UnitPoynting_)    
     unitstr_IDL_J           = IoUnitStr(UnitJ_)           
     unitstr_IDL_electric    = IoUnitStr(UnitElectric_)    
     unitstr_IDL_temperature = IoUnitStr(UnitTemperature_) 
     unitstr_IDL_DivB        = 'T/m'

  case default
     call stop_mpi(NameThisComp//': IoUnitType='//IoUnits// ' is invalid.')

  end select

end subroutine set_units

!==============================================================================

subroutine init_mhd_variables

  use ModVarIndexes
  use ModPhysics

  implicit none

  integer :: iVar
  !--------------------------------------------------------------------------

  UnitUser_V(Rho_)     = UnitUser_Rho
  UnitUser_V(RhoUx_)   = UnitUser_RhoU
  UnitUser_V(RhoUy_)   = UnitUser_RhoU
  UnitUser_V(RhoUz_)   = UnitUser_RhoU
  UnitUser_V(Bx_)      = UnitUser_B
  UnitUser_V(By_)      = UnitUser_B
  UnitUser_V(Bz_)      = UnitUser_B
  UnitUser_V(p_)       = UnitUser_p
  UnitUser_V(Energy_)  = UnitUser_EnergyDens

  NameUnitUserTec_V(rho_)    = UnitStr_Tec_Rho
  NameUnitUserTec_V(rhoUx_)  = UnitStr_Tec_RhoU
  NameUnitUserTec_V(rhoUy_)  = UnitStr_Tec_RhoU
  NameUnitUserTec_V(rhoUz_)  = UnitStr_Tec_RhoU
  NameUnitUserTec_V(Bx_)     = UnitStr_Tec_B
  NameUnitUserTec_V(By_)     = UnitStr_Tec_B
  NameUnitUserTec_V(Bz_)     = UnitStr_Tec_B
  NameUnitUserTec_V(P_)      = UnitStr_Tec_p
  NameUnitUserTec_V(Energy_) = UnitStr_Tec_EnergyDens

  NameUnitUserIdl_V(rho_)    = UnitStr_Idl_Rho
  NameUnitUserIdl_V(rhoUx_)  = UnitStr_Idl_RhoU
  NameUnitUserIdl_V(rhoUy_)  = UnitStr_Idl_RhoU
  NameUnitUserIdl_V(rhoUz_)  = UnitStr_Idl_RhoU
  NameUnitUserIdl_V(Bx_)     = UnitStr_Idl_B
  NameUnitUserIdl_V(By_)     = UnitStr_Idl_B
  NameUnitUserIdl_V(Bz_)     = UnitStr_Idl_B
  NameUnitUserIdl_V(P_)      = UnitStr_Idl_p
  NameUnitUserIdl_V(Energy_) = UnitStr_Idl_EnergyDens

  ! By default the scalar advected variables are assumed to behave like density
  do iVar = ScalarFirst_, ScalarLast_
     UnitUser_V(iVar)        = UnitUser_Rho
     NameUnitUserTec_V(iVar) = UnitStr_Tec_Rho
     NameUnitUserIdl_V(iVar) = UnitStr_Idl_Rho
  end do

end subroutine init_mhd_variables
