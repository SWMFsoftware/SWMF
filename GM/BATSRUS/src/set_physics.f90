!^CFG COPYRIGHT UM
!\
! set_physics_constants set normalizations, physical constants, 
! problem_type dependent physical parameters
!/
subroutine set_physics_constants
  use ModProcMH
  use ModMain
  use ModPhysics
  use CON_axes,   ONLY: get_axes
  use CON_planet, ONLY: get_planet
  use ModVarIndexes
  use ModUser
  implicit none

  real :: Qqpmag, Oopmag, Gsun
  real :: Mbody_dim
  real :: MBody2Dim                 !^CFG IF SECONDBODY
  real :: cosTheta, sinTheta, cosPhi, sinPhi, & 
       xx, yy, zz

  integer :: i, iBoundary

  logical :: oktest, oktest_me

  call set_oktest('set_physics_constants',oktest, oktest_me)

  !\
  ! Load body rotation rates, masses, and radii by problem type
  !/
  select case(problem_type)
  case (problem_rotation,problem_saturn)
     unitSI_x       = Rsaturn       ! Radius - NOTE - MUST BE IN meters
     Mbody_dim      = Msaturn       ! Mass of body in kg
     rot_period_dim = RotationPeriodSaturn/3600.0
  case (problem_heliosphere, problem_cme)
     unitSI_x       = Rsun                      ! Radius - MUST BE IN meters
     Mbody_dim      = Msun                      ! Mass of body in kg
     rot_period_dim = RotationPeriodSun/3600.0  ! rotation period in hours
     !^CFG IF GLOBALHELIOSPHERE BEGIN
  case (problem_globalhelio)
     unitSI_x       = 215.0*Rsun                ! Radius - MUST BE IN meters
     Mbody_dim      = Msun                      ! Mass of body in kg
     rot_period_dim = RotationPeriodSun/3600.0  ! rotation period in hours
     !^CFG END GLOBALHELIOSPHERE
  case (problem_earth)
     unitSI_x       = Rearth        ! Radius - MUST BE IN meters
     Mbody_dim      = Mearth        ! Mass of body in kg
     rot_period_dim = RotationPeriodEarth/3600.0
  case (problem_jupiter)
     unitSI_x       = Rjupiter      ! Radius - NOTE - MUST BE IN meters
     Mbody_dim      = Mjupiter      ! Mass of body in kg
     rot_period_dim = RotationPeriodJupiter/3600.0
  case (problem_venus)
     unitSI_x       = Rvenus        ! Radius - NOTE - MUST BE IN meters
     Mbody_dim      = Mvenus        ! Mass of body in kg
     rot_period_dim = RotationPeriodVenus/3600.0
 case (problem_mars)
     unitSI_x       = Rmars        ! Radius - NOTE - MUST BE IN meters
     Mbody_dim      = Mmars        ! Mass of body in kg
     rot_period_dim = RotationPeriodMars/3600.0
  case (problem_arcade)
     Gsun           = -cGravitation*Msun/(Rsun*Rsun)
     unitSI_x       = SSPsun**2/abs(Gsun)
     SSPsun         = sqrt( cBoltzmann*TArcDim/(muArc*cProtonMass) )
  case (problem_comet)
     unitSI_x       = cE9           ! 1 million km=1E+9 m
     Mbody_dim      = 1.0           ! Mass of body in kg
     rot_period_dim = 1.0           ! rotation period in hours
  case default
     unitSI_x       = 1.0           ! Radius - NOTE - MUST BE IN meters
     Mbody_dim      = 0.00          ! Mass of body in kg
     rot_period_dim = 0.00          ! rotation period in hours
  end select

  ! Second body mass is set to zero by default   !^CFG IF SECONDBODY
  MBody2Dim = 0.0                                !^CFG IF SECONDBODY

  !\
  ! Call set_dimensions, which set the quantities for converting from
  ! normalized to  dimensional quantities and vice versa.  Also
  ! sets input and output strings for the conversion quantities
  !/
  call set_dimensions

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
  ! set the speed of light and get normalization
  !/
  cLIGHT = boris_cLIGHT_factor * cLightSpeed/unitSI_U
  c2LIGHT = cLIGHT**2
  inv_c2LIGHT = cOne/c2LIGHT

  !\
  ! Convert rotation, and gravity to non-dimensional values
  !/
  ! if the rotation period is less than 1 second then you made
  ! a mistake - the period is to fast
  if (abs(rot_period_dim) > 1./3600.) then
     OmegaBody = (cTwoPi/(rot_period_dim*3600.00)) / (1.0/unitSI_t)
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
  !^CFG IF ALWAVES BEGIN
  !Gbody  = Gbody/Tnot  !^CFG UNCOMMENT IF ALWAVES
  !^CFG END ALWAVES

  !\
  ! Nondimensionalize dimensional SW values - 
  !/
  SW_a_dim = unitUSER_U
  SW_p_dim = unitUSER_p*inv_g
  SW_B_factor = unitUSER_B

  if(problem_type==problem_globalhelio)then !^CFG IF GLOBALHELIOSPHERE BEGIN
     ! unitUSER_rho #/cm3
     SW_rho = SW_rho_dim/unitUSER_rho
  else                                      !^CFG END GLOBALHELIOSPHERE
     SW_rho = 1.0    
  end if                                    !^CFG IF GLOBALHELIOSPHERE
  SW_p   = inv_g
  SW_Ux  = SW_Ux_dim/unitUSER_U
  SW_Uy  = SW_Uy_dim/unitUSER_U
  SW_Uz  = SW_Uz_dim/unitUSER_U
  SW_Bx  = SW_Bx_dim/unitUSER_B
  SW_By  = SW_By_dim/unitUSER_B
  SW_Bz  = SW_Bz_dim/unitUSER_B

  !\
  ! Convert comet input parameters to SI
  !/
  kin=kin_in*1E-6
  Unr=Unr_in*1E3

  !\
  ! Nondimensionalize body pressure and density
  ! Note that the Body_rho_dim is not really a rho, but rather an n
  ! so the correct normalizing factor is unitUSER_n
  !/
  select case(problem_type)
  case(problem_shocktube, problem_uniform)
     Body_rho= Body_rho_dim
     Body_p  = Body_rho*Body_T_dim
     RhoBody2= RhoDimBody2                          !^CFG IF SECONDBODY
     pBody2  = RhoDimBody2*TDimBody2                !^CFG IF SECONDBODY

     !^CFG IF GLOBALHELIOSPHERE BEGIN
     !Correcting the problem of the pressure and density near the body
  case(problem_globalhelio)
     !
     ! setting Body_rho_dim to be like SW_rho_dim, look at write_plot_common
     !
     Body_rho_dim = SW_rho_dim
     Body_T_dim = SW_T_dim
     Body_rho = Body_rho_dim/unitUSER_rho
     Body_p = inv_g
     RhoBody2= Body_rho                            !^CFG IF SECONDBODY
     pBody2  = Body_p                              !^CFG IF SECONDBODY
     !^CFG END GLOBALHELIOSPHERE

  case default
     Body_rho= Body_rho_dim/unitUSER_n
     Body_p  = cBoltzmann*Body_rho_dim*&
               1.0E6*Body_T_dim/unitSI_p
     RhoBody2= RhoDimBody2/unitUSER_n               !^CFG IF SECONDBODY
     pBody2  = cBoltzmann*RhoDimBody2*&             !^CFG IF SECONDBODY
               1.0E6*TDimBody2/unitSI_p             !^CFG IF SECONDBODY
  end select

  !Here the arrays of the FACE VALUE are formed
  !Initialization
  do iBoundary=body2_,Top_
     FaceState_VI(:,iBoundary)=DefaultState_V
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
     FaceState_VI(Bx_,  iBoundary) =  SW_Bx
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

  !\
  ! Now do extra varibles that are problem dependent.
  ! As well as case by case exceptions to the above.
  !/
  select case (problem_type)
  case (problem_earth)
     !     !\
     !     ! Temporary for the paleo earth case
     !     !/
     !
     !     WARNING: CHECK IF SWMF CAN HANDLE NON-DIPOLE PLANET FIELD ALREADY
     !
     !     ! set dipole to zero
     !     Bdpx = 0.0
     !     Bdpy = 0.0
     !     Bdpz = .1*Bdp
     !
     !     ! Quadrupole strength     
     !     Qqpmag = .1*Bdp
     !     
     !     Qqp(1,1) = 1.00*Qqpmag         ! 1
     !     Qqp(1,2) = 0.00                ! 2
     !     Qqp(1,3) = 0.00                ! 3
     !     Qqp(2,1) = 0.0
     !     Qqp(2,2) = 1.00*Qqpmag         ! 4 
     !     Qqp(2,3) = 0.00                ! 5
     !     Qqp(3,1) = 0.00
     !     Qqp(3,2) = 0.00
     !     Qqp(3,3) = -(Qqp(1,1)+Qqp(2,2))

  case (problem_saturn)

  case (problem_jupiter)

  case (problem_heliosphere)
     !\
     ! Heliospheric flow problem.
     !/
     RHOsun = cProtonMass*Body_rho_dim
     Presun = cBoltzmann*RHOsun*Body_T_dim/cProtonMass
     SSPsun = sqrt(g*cBoltzmann*Body_T_dim/cProtonMass)
     Velsun = 0.75

     !\
     ! Chip's Flux Rope scaling
     !/
     B1scl   = 1.0E-4*cme_B1_dim/sqrt(cMu*RHOsun*SSPsun*SSPsun)
     a1scl   = 1.0E-4*cme_a1/sqrt(cMu*RHOsun*SSPsun*SSPsun)
     rho1scl = cme_rho1/RHOsun
     rho2scl = cme_rho2/RHOsun

     ! Dipole strength
     Bdp = 1.0E-4*Bdp_dim/sqrt(cMu*RHOsun*SSPsun*SSPsun)

     Bdpx = -sinTHETAtilt*Bdp       ! 1
     Bdpy = 0.00                    ! 2
     Bdpz = cosTHETAtilt*Bdp        ! 3 

     ! Quadrapole strength     
     Qqpmag = 0.00

     Qqp(1,1) = 1.00*Qqpmag         ! 1
     Qqp(1,2) = 0.00                ! 2
     Qqp(1,3) = 0.00                ! 3
     Qqp(2,1) = Qqp(1,2)
     Qqp(2,2) = 0.00*Qqpmag         ! 4 
     Qqp(2,3) = 0.00                ! 5
     Qqp(3,1) = Qqp(1,3)
     Qqp(3,2) = Qqp(2,3)
     Qqp(3,3) = -(Qqp(1,1)+Qqp(2,2))

     ! Octupole strength
     !====== Milestone 2 Octupole strength ==========
     !Oopmag = 0.00
     !===============================================

     !====== Milestone 3 Octupole strength ==========
     Oopmag = Bdp/9.00
     !===============================================

     Oop(1,1,1) = 0.00*Oopmag              ! 1
     Oop(1,1,2) = 0.00*Oopmag              ! 2
     Oop(1,1,3) = -1.00*Oopmag             ! 3
     Oop(1,2,1) = Oop(1,1,2)
     Oop(1,2,2) = 0.00*Oopmag              ! 4
     Oop(1,2,3) = 0.00                     ! 5
     Oop(1,3,1) = Oop(1,1,3) 
     Oop(1,3,2) = Oop(1,2,3)
     Oop(1,3,3) = -(Oop(1,1,1)+Oop(1,2,2))

     Oop(2,1,1) = Oop(1,1,2)
     Oop(2,1,2) = Oop(1,2,2)
     Oop(2,1,3) = Oop(1,2,3)
     Oop(2,2,1) = Oop(1,2,2)
     Oop(2,2,2) = 0.00*Oopmag              ! 6
     Oop(2,2,3) = -1.00*Oopmag             ! 7
     Oop(2,3,1) = Oop(1,2,3)
     Oop(2,3,2) = Oop(2,2,3)
     Oop(2,3,3) = -(Oop(2,1,1)+Oop(2,2,2))

     Oop(3,1,1) = Oop(1,1,3)
     Oop(3,1,2) = Oop(1,2,3)
     Oop(3,1,3) = Oop(1,3,3)
     Oop(3,2,1) = Oop(1,2,3)
     Oop(3,2,2) = Oop(2,2,3)
     Oop(3,2,3) = Oop(2,3,3)
     Oop(3,3,1) = Oop(1,3,3)
     Oop(3,3,2) = Oop(2,3,3)
     Oop(3,3,3) = -(Oop(3,1,1)+Oop(3,2,2))

     ! Speed of light used in Boris correction
     cLIGHT = boris_cLIGHT_factor * cLightSpeed/SSPsun
     c2LIGHT = cLIGHT**2
     inv_c2LIGHT = cOne/c2LIGHT

  case (problem_cme)
     !\
     ! CME flux rope problem.
     !/
     RHOsun = Body_rho_dim*1.0E+6*cProtonMass
     Presun = cBoltzmann/cProtonMass*RHOsun*Body_T_dim
     SSPsun = sqrt(g*cBoltzmann*Body_T_dim/cProtonMass)
     Velsun = 0.75
     !\
     ! put solution variables in nondimensional form
     !/
     Rscl    = 1.0
     RHOscl  = 1.0
     SSPscl  = 1.0
     Vscl    = cme_v_erupt/SSPsun
     Bdp     = 1.0E-4*Bdp_dim/sqrt(cMu*RHOsun*SSPsun*SSPsun)
     B1scl   = 1.0E-4*cme_B1_dim/sqrt(cMu*RHOsun*SSPsun*SSPsun)
     a1scl   = 1.0E-4*cme_a1/sqrt(cMu*RHOsun*SSPsun*SSPsun)
     rho1scl = cme_rho1/RHOsun
     rho2scl = cme_rho2/RHOsun
  case (problem_arcade)
     Phtscl  = 1.0
     RHOscl  = 1.0
     SSPscl  = 1.0
     Vscl    = UzArcDim/SSPsun
     Gbody   = Gsun*(UnitSI_x/(SSPsun*SSPsun))
     B0_scl  = 1.0E-4*BArcDim /sqrt(cMu*RhoArcDim*SSPsun*SSPsun)
     B0y_scl = 1.0E-4*ByArcDim/sqrt(cMu*RhoArcDim*SSPsun*SSPsun)
  end select
  if( UseUserSetPhysConst)call user_set_physics

end subroutine set_physics_constants

!===========================================================================

subroutine set_dimensions
  use ModProcMH, ONLY:iProc
  use ModMain
  use ModPhysics
  use ModVarIndexes
  implicit none

  logical :: oktest, oktest_me

  call set_oktest('set_dimensions',oktest, oktest_me)

  !\
  ! Load variables used for converting from dimensional to non-dimensional 
  ! units and back.  Also load the name variables for each of the units for
  ! use in writing output.
  !/

  unitSI_angle = 180/cPi   
  unitUSER_angle = unitSI_angle
  unitstr_TEC_angle = '[degree]'            
  unitstr_IDL_angle = '--'            

  !\
  ! set independent normalizing SI variables first - on a case by case basis
  !/
  select case(problem_type)
  case(problem_earth,problem_saturn,problem_jupiter,problem_rotation, &
       problem_venus,problem_comet, problem_mars)
     unitSI_rho  = cProtonMass*(cMillion*SW_rho_dim)         ! kg/m^3
     unitSI_U    = sqrt(g*cBoltzmann*SW_T_dim/cProtonMass)   ! m/s
  case(problem_arcade)
     unitSI_rho  = RHOsun                                    ! kg/m^3
     unitSI_U    = SSPsun                                    ! m/s
  case(problem_heliosphere,problem_cme)
     unitSI_x    = Rsun                                      ! m
     unitSI_rho  = cProtonMass*Body_rho_dim                  ! kg/m^3
     unitSI_U    = sqrt(g*cBoltzmann*Body_T_dim/cProtonMass) ! m/s  (SSPsun)
     !^CFG IF GLOBALHELIOSPHERE BEGIN
     ! SW_rho_dim is in units of n/cc
    case(problem_globalhelio)
       unitSI_x    = 215.0*Rsun                              ! m
       unitSI_rho  = cProtonMass*SW_rho_dim*1.0E+6           ! kg/m^3
       unitSI_U    = sqrt(g*cBoltzmann*SW_T_dim/cProtonMass) ! m/s  (SSPsun)
       !^CFG END GLOBALHELIOSPHERE
  case(problem_dissipation)
     if (.not.UseDefaultUnits) then
        unitSI_x    = Length0Diss                            ! m
        unitSI_t    = Time0Diss                              ! s
        unitSI_U    = unitSI_x/unitSI_t                      ! m/s
        unitSI_rho  = Rho0Diss                               ! kg/m3
     else
        unitSI_x    = 1.0                                    ! m
        unitSI_t    = 1.0                                    ! s
        unitSI_U    = 1.0                                    ! m/s
        unitSI_rho  = 1.0                                    ! kg/m3
     end if
  case default
     unitSI_x    = 1.0                                       ! dimensionless
     unitSI_rho  = 1.0                                       ! dimensionless
     unitSI_U    = 1.0                                       ! dimensionless
  end select
  !\
  ! set other normalizing SI variables from the independent ones
  !/
  unitSI_t           = unitSI_x/unitSI_U                     ! s
  unitSI_n           = unitSI_rho/cProtonMass                     ! (#/m^3) (amu/m^3)
  unitSI_p           = unitSI_rho*unitSI_U**2                ! Pa
  unitSI_B           = unitSI_U*sqrt(cMu*unitSI_rho)         ! T
  unitSI_rhoU        = unitSI_rho*unitSI_U                   ! kg/m^2/s
  unitSI_energydens  = unitSI_p                              ! J/m^3 - energy density
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
  ! case the temperature code using unitSI_temperature is really
  ! Tcode = (Ti+Te)/A.
  !
  ! For the special case where A=2 we have Tcode = 1/2*(Te+Ti).  If
  ! we assume Ti=Te then Tcode=Ti=Te.
  unitSI_temperature = (unitSI_p/unitSI_rho)*(cProtonMass/cBoltzmann) ! Kelvin 

  !\
  ! set USER variables used for normalization, input, and output
  ! They all have user defined coordinates and are not consistant in a unit sense.
  ! Every variable is case dependent.
  ! 
  ! Also load the string variables associated with the USER variables - 
  ! note that they are loaded differently for IDL and TEC output
  !/
  select case(problem_type)
  case(problem_earth,problem_saturn,problem_jupiter,problem_rotation, &
       problem_venus,problem_comet,problem_mars)
     ! load units
     unitUSER_x           = unitSI_x/unitSI_x                ! planetary radii
     unitUSER_rho         = 1.0E-6*unitSI_n                  ! (#/cm^3) (amu/cm^3)
     unitUSER_U           = 1.0E-3*unitSI_U                  ! km/s
     unitUSER_t           = unitSI_t                         ! s
     unitUSER_n           = 1.0E-6*unitSI_n                  ! (#/cm^3) (amu/cm^3)
     unitUSER_p           = 1.0E+9*unitSI_p                  ! nPa
     unitUSER_B           = 1.0E+9*unitSI_B                  ! nT
     unitUSER_rhoU        = unitSI_rhoU                      ! kg/m^2/s
     unitUSER_energydens  = unitSI_energydens                ! J/m^3 - energy density
     unitUSER_Poynting    = unitSI_Poynting                  ! J/m^2/s
     unitUSER_J           = 1.0E+6*unitSI_J                  ! uA/m^2
     unitUSER_electric    = 1.0E+3*unitSI_electric           ! mV/m
     unitUSER_DivB        = 1.0E+9*unitSI_DivB*unitSI_x      ! nT/planetary radii
     ! set temperature - note that the below is only strictly true for a
     ! limited number of cases.  If the only ion is proton then 
     ! Tcode = Te+Ti.
     !
     ! If the ions are "heavy" (m = A * mp) then n above is not really a 
     ! number density but is a nucleon density (#nucleons/m^3).  In this 
     ! case the temperature code using unitSI_temperature is really
     ! Tcode = (Ti+Te)/A.
     !
     ! For the special case where A=2 we have Tcode = 1/2*(Te+Ti).  If
     ! we assume Ti=Te then Tcode=Ti=Te..
     unitUSER_temperature = unitSI_temperature               ! Kelvin 

     !\
     ! set string variables used for writing output - TECPLOT
     !/
     unitstr_TEC_x           = '[R]'            
     unitstr_TEC_rho         = '[amu/cm^3]'       
     unitstr_TEC_U           = '[km/s]'          
     unitstr_TEC_t           = '[s]'             
     unitstr_TEC_n           = '[amu/cm^3]'        
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
     unitstr_IDL_n           = 'amu/cm3'        
     unitstr_IDL_p           = 'nPa'           
     unitstr_IDL_B           = 'nT'            
     unitstr_IDL_rhoU        = 'kg/m2s2'
     unitstr_IDL_energydens  = 'J/m3'           
     unitstr_IDL_Poynting    = 'J/m^2s'
     unitstr_IDL_J           = 'uA/m2'       
     unitstr_IDL_electric    = 'mV/m'          
     unitstr_IDL_DivB        = 'nT/R'           
     unitstr_IDL_temperature = 'K'             

  case(problem_heliosphere,problem_cme)
     unitUSER_x           = unitSI_x/unitSI_x                ! R
     unitUSER_rho         = 1.0e-3*unitSI_rho                ! g/cm^3
     unitUSER_U           = 1.0E-3*unitSI_U                  ! km/s
     unitUSER_t           = unitSI_t                         ! s
     unitUSER_n           = 1.0E-6*unitSI_n                  ! (#/cm^3) (amu/cm^3)
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
     
  case(problem_dissipation)
     ! define the correct temperature::
     unitSI_temperature = (unitSI_p/unitSI_rho)   ! Kelvin 
     ! defalut case - same as the SI default above
     unitUSER_x           = unitSI_X              ! dimensionless
     unitUSER_rho         = unitSI_rho            ! dimensionless
     unitUSER_U           = unitSI_U              ! dimensionless
     unitUSER_t           = unitSI_t              ! dimensionless   
     unitUSER_n           = unitSI_n              ! dimensionless                
     unitUSER_p           = unitSI_p              ! dimensionless                 
     unitUSER_B           = unitSI_B              ! dimensionless                
     unitUSER_rhoU        = unitSI_rhoU           ! dimensionless                
     unitUSER_energydens  = unitSI_energydens     ! dimensionless                
     unitUSER_J           = unitSI_J              ! dimensionless                
     unitUSER_electric    = unitSI_electric       ! dimensionless                
     unitUSER_DivB        = unitSI_DivB           ! dimensionless                
     ! set temperature - note that the below is  only strictly true for a
     ! pure proton plasma.  If the are heavy ions of mass #*mp then you could
     ! be off in temperature by as much as a factor of #.  There is no way around
     ! this in MHD.
     unitUSER_temperature = unitSI_temperature    ! dimensionless
     !\
     ! set string variables used for writing output - TECPLOT
     !/
     unitstr_TEC_x           = 'm'
     unitstr_TEC_rho         = 'kg/m^3'
     unitstr_TEC_U           = 'm/s'
     unitstr_TEC_t           = 's'
     unitstr_TEC_n           = '1/m^3'
     unitstr_TEC_p           = 'Pa'
     unitstr_TEC_B           = 'T'
     unitstr_TEC_rhoU        = 'kg/m^2/s'
     unitstr_TEC_Poynting    = 'J/m^2/s'
     unitstr_TEC_energydens  = 'J/m^3'
     unitstr_TEC_J           = 'A/m^2'
     unitstr_TEC_electric    = 'V/m'
     unitstr_TEC_DivB        = 'T/m'
     unitstr_TEC_temperature = 'K'
     !\
     ! set string variables used for writing output - IDL
     !/
     unitstr_IDL_x           = 'm'
     unitstr_IDL_rho         = 'kg/m^3'
     unitstr_IDL_U           = 'm/s'
     unitstr_IDL_t           = 's'
     unitstr_IDL_n           = '1/m^3'
     unitstr_IDL_p           = 'Pa'
     unitstr_IDL_B           = 'T'
     unitstr_IDL_rhoU        = 'kg/m^2/s'
     unitstr_IDL_Poynting    = 'J/m^2/s'
     unitstr_IDL_energydens  = 'J/m^3'
     unitstr_IDL_J           = 'A/m^2'
     unitstr_IDL_electric    = 'V/m'
     unitstr_IDL_DivB        = 'T/m'
     unitstr_IDL_temperature = 'K'
     if(oktest.and.iProc==0) then
        write(*,*) 'unitSI_t           = ',unitSI_t
        write(*,*) 'unitSI_x           = ',unitSI_x
        write(*,*) 'unitSI_U           = ',unitSI_U
        write(*,*) 'unitSI_rho         = ',unitSI_rho
        write(*,*) 'unitSI_p           = ',unitSI_p
        write(*,*) 'unitSI_energydens  = ',unitSI_energydens
        write(*,*) 'unitSI_temperature = ',unitSI_temperature
     end if
     
  case default
     ! defalut case - same as the SI default above
     unitUSER_x           = unitSI_X              
     unitUSER_rho         = unitSI_rho            
     unitUSER_U           = unitSI_U              
     unitUSER_t           = unitSI_t              
     unitUSER_n           = unitSI_n                        
     unitUSER_p           = unitSI_p                         
     unitUSER_B           = unitSI_B                        
     unitUSER_rhoU        = unitSI_rhoU                     
     unitUSER_energydens  = unitSI_energydens               
     unitUSER_Poynting    = unitSI_Poynting       
     unitUSER_J           = unitSI_J                        
     unitUSER_electric    = unitSI_electric                 
     unitUSER_DivB        = unitSI_DivB                     
     ! set temperature - note that the below is  only strictly true for a
     ! pure proton plasma.  If the are heavy ions of mass #*mp then you could
     ! be off in temperature by as much as a factor of #.  There is no way around
     ! this in MHD.
     unitUSER_temperature = unitSI_temperature  ! dimensionless
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

  end select
  unitUSERVars_V(rho_)     = unitUSER_rho
  unitUSERVars_V(rhoUx_)   = unitUSER_rhoU
  unitUSERVars_V(rhoUy_)   = unitUSER_rhoU
  unitUSERVars_V(rhoUz_)   = unitUSER_rhoU
  unitUSERVars_V(Bx_)      = unitUSER_B
  unitUSERVars_V(By_)      = unitUSER_B
  unitUSERVars_V(Bz_)      = unitUSER_B
  unitUSERVars_V(P_)   = unitUSER_p

  TypeUnitVarsTec_V(rho_)    = unitstr_TEC_rho
  TypeUnitVarsTec_V(rhoUx_)  = unitstr_TEC_rhoU
  TypeUnitVarsTec_V(rhoUy_)  = unitstr_TEC_rhoU
  TypeUnitVarsTec_V(rhoUz_)  = unitstr_TEC_rhoU
  TypeUnitVarsTec_V(Bx_)     = unitstr_TEC_B
  TypeUnitVarsTec_V(By_)     = unitstr_TEC_B
  TypeUnitVarsTec_V(Bz_)     = unitstr_TEC_B
  TypeUnitVarsTec_V(P_)  = unitstr_TEC_p

  TypeUnitVarsIdl_V(rho_)    = unitstr_IDL_rho
  TypeUnitVarsIdl_V(rhoUx_)  = unitstr_IDL_rhoU
  TypeUnitVarsIdl_V(rhoUy_)  = unitstr_IDL_rhoU
  TypeUnitVarsIdl_V(rhoUz_)  = unitstr_IDL_rhoU
  TypeUnitVarsIdl_V(Bx_)     = unitstr_IDL_B
  TypeUnitVarsIdl_V(By_)     = unitstr_IDL_B
  TypeUnitVarsIdl_V(Bz_)     = unitstr_IDL_B
  TypeUnitVarsIdl_V(P_)  = unitstr_IDL_p
end subroutine set_dimensions
