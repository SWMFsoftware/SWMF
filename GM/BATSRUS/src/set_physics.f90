!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!
subroutine set_physics_constants 

  ! set normalizations, physical constants, etc.

  use ModProcMH
  use ModAMR
  use ModMain
  use ModPhysics
  use CON_planet,  ONLY: get_planet, NamePlanet
  use ModVarIndexes
  use ModMultiFluid
  use ModAdvance, ONLY: UseElectronPressure, Pe_, UseAnisoPressure, Ppar_

  implicit none

  real :: MassBodySi
  real :: MassBody2Si
  real :: pCoef

  integer :: i, iVar, iBoundary

  character (len=*), parameter :: NameSub = "set_physics_constants"

  logical :: DoTest, DoTestMe
  !----------------------------------------------------------------------------
  call set_oktest(NameSub, DoTest, DoTestMe)

  !\
  ! Load body rotation rates, masses, and radii
  !/
  NamePlanetRadius = 'R'
  select case(NameThisComp)
  case('GM')
     call get_planet( &
          RadiusPlanetOut   = rPlanetSi, &
          MassPlanetOut     = MassBodySi, &
          RotationPeriodOut = RotPeriodSi)
     if(NamePlanet == 'NONE') then
        rPlanetSi = 1000.0  ! 1 km
        NamePlanetRadius = 'km'
     end if
  case('SC', 'IH', 'OH')
     if(UseStar) then
        rPlanetSi   = RadiusStar*rSun
        MassBodySi  = MassStar*mSun
        RotPeriodSi = RotationPeriodStar*cSecondPerDay 
     else
        rPlanetSi   = rSun
        MassBodySi  = mSun
        RotPeriodSi = RotationPeriodSun
     endif
     SW_n_dim    = BodyNDim_I(IonFirst_)  ! Needed for SOLARWIND normalization only
     SW_T_dim    = BodyTDim_I(IonFirst_)  ! Needed for SOLARWIND normalization only
  case('EE')
     MassBodySi = mSun
     RotPeriodSi = RotationPeriodSun
     rPlanetSi = 1000.0  ! 1 km
     NamePlanetRadius = 'km'
  end select
 
  ! Note for GM  !!! BATSRUS's OmegaBody is siderial (relative to the Sun)
  ! and it is DIFFERENT from SWMF's inertial OmegaPlanet 
  ! defined in CON_planet !!!
  if (RotPeriodSi == cZero) then
     OmegaBody = 0.0
  else
     OmegaBody = cTwoPi/RotPeriodSi
  end if
  ! Second body mass is set to zero by default
  MassBody2Si = 0.0


  ! Make sure that MassIon_I is consistent with MassFluid_I
  MassIon_I = MassFluid_I(IonFirst_:IonLast_)

  !\
  ! Call set_units, which set the quantities for converting from
  ! normalized to  dimensional quantities and vice versa.  Also
  ! sets input and output strings for the conversion quantities
  !/
  call set_units

  if(DoTest .and. iProc==0) then
     write(*,'(4a15)')'No2Io_V','NameIdlUnit_V','No2Si_V','NameSiUnit_V'
     do i=1, nIoUnit
        write(*,'(es15.6," ",a6,"        ",es15.6," ",a6)') &
             No2Io_V(i), NameIdlUnit_V(i), &
             No2Si_V(i), NameSiUnit_V(i)
     end do
  end if

  ! Average ion mass per charge used by Hall MHD and
  ! by ion-electron thermal exchange term
  if(TypeNormalization == 'NONE')then
     IonMassPerCharge = 1.0
  else
     IonMassPerCharge = cProtonMass/cElectronCharge &
          *Si2No_V(UnitMass_)/Si2No_V(UnitCharge_)
  end if

  ! Electron charge in normalized units (actually it is proton charge/mass)
  ! This is useful in formulas like n_s q_s (u_s - u_+) x B 
  ! in the multi-ion momentum equation. It is used in some user modules,
  ! so we keep the name.
  ElectronCharge = 1.0/IonMassPerCharge

  ! For single ion fluid the average ion mass per charge is constant
  if(.not. (UseMultiSpecies .or. UseMultiIon)) &
       IonMassPerCharge = IonMassPerCharge * MassIon_I(1)/ChargeIon_I(1)

  !\
  ! set the (corrected) speed of light and get normalization
  !/
  Clight      = Boris_Clight_Factor * cLightSpeed * Si2No_V(UnitU_)
  C2light     = cLIGHT**2
  InvClight   = cOne/cLight
  Inv_C2light = cOne/c2LIGHT

  ! normalize the radiation constant
  cRadiationNo = cRadiation &
       * Si2No_V(UnitEnergyDens_) / Si2No_V(UnitTemperature_)**4

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
        write(*,*) "WARNING in ",NameSub,":"
        write(*,*) "You have set UseRotatingFrame = ",UseRotatingFrame
        write(*,*) "but RotPeriodSi in seconds= ",RotPeriodSi
        write(*,*) "This is too fast! Setting OmegaBody=0.0           "
        write(*,*) "--------------------------------------------------"
     end if
     OmegaBody = 0.0
  end if

  if(GravityDir == 0)then
     ! Note: The mass of the body is in SI units
     Gbody  = -cGravitation*MassBodySi*(Si2No_V(UnitU_)**2 * Si2No_V(UnitX_))
  else
     ! Normalize gravitational acceleration
     Gbody  = GravitySi*(Si2No_V(UnitU_)**2 / Si2No_V(UnitX_))
  end if


  GBody2 = -cGravitation*MassBody2Si*(Si2No_V(UnitU_)**2 * Si2No_V(UnitX_))


  !\
  ! Normalize solar wind values. Note: the solarwind is in I/O units
  !/
  SW_n   = SW_n_dim*Io2No_V(UnitN_)
  SW_rho = SW_n * MassIon_I(1)
  SW_p   = SW_n * SW_T_dim*Io2No_V(UnitTemperature_)
  SW_Ux  = SW_Ux_dim*Io2No_V(UnitU_)
  SW_Uy  = SW_Uy_dim*Io2No_V(UnitU_)
  SW_Uz  = SW_Uz_dim*Io2No_V(UnitU_)
  SW_Bx  = SW_Bx_dim*Io2No_V(UnitB_)
  SW_By  = SW_By_dim*Io2No_V(UnitB_)
  SW_Bz  = SW_Bz_dim*Io2No_V(UnitB_)

  if(.not.UseElectronPressure .and. IsMhd) &
       SW_p = SW_p*(1 + ElectronPressureRatio)

  ! These are useful for printing out values
  SW_rho_dim = SW_rho*No2Io_V(UnitRho_)
  SW_p_dim   = SW_p*No2Io_V(UnitP_)

  ! The normalized quantities extend to the first MHD fluid too
  BodyRho_I(IonFirst_:) = BodyNDim_I*Io2No_V(UnitN_)*MassFluid_I
  BodyP_I(IonFirst_:)   = BodyNDim_I*Io2No_V(UnitN_)*BodyTDim_I &
       * Io2No_V(UnitTemperature_)

  PolarRho_I(IonFirst_:) = PolarNDim_I*Io2No_V(UnitN_)*MassFluid_I
  PolarP_I(IonFirst_:)   = PolarNDim_I*Io2No_V(UnitN_)*PolarTDim_I &
       * Io2No_V(UnitTemperature_)
  PolarU_I(IonFirst_:)   = PolarUDim_I*Io2No_V(UnitU_)

  if(UseMultiIon .and. IsMhd)then
     ! Add up ion fluids for total ion fluid
     BodyRho_I(1)  = sum(BodyRho_I(IonFirst_:IonLast_))
     BodyP_I(1)    = sum(BodyP_I(IonFirst_:IonLast_))
     PolarRho_I(1) = sum(PolarRho_I(IonFirst_:IonLast_))
     PolarP_I(1)   = sum(PolarP_I(IonFirst_:IonLast_))
     PolarU_I(1)   = &
          sum(PolarRho_I(IonFirst_:IonLast_)*PolarU_I(IonFirst_:IonLast_)) &
          /sum(PolarRho_I(IonFirst_:IonLast_))
  end if

  if(.not.UseElectronPressure .and. IsMhd) then
     BodyP_I(1)  = BodyP_I(1)*(1 + ElectronPressureRatio)
     PolarP_I(1) = PolarP_I(1)*(1 + ElectronPressureRatio)
  end if


  RhoBody2= RhoDimBody2 *Io2No_V(UnitN_)*MassIon_I(1)
  pBody2  = RhoBody2 * TDimBody2*Io2No_V(UnitTemperature_)



  ! Here the arrays of the FACE VALUE are formed
  ! Initialization
  do iBoundary=body2_,6
     FaceState_VI(:,iBoundary)=DefaultState_V(1:nVar)
  end do

  !For bodies:
  FaceState_VI(iRho_I, Body1_) = BodyRho_I
  FaceState_VI(iP_I,   Body1_) = BodyP_I

  if (UseMultiSpecies) then
     FaceState_VI(SpeciesFirst_, Body1_) = &
          BodyRho_I(1)*(1.0 - LowDensityRatio * (SpeciesLast_-SpeciesFirst_))
     FaceState_VI(SpeciesFirst_+1:SpeciesLast_, Body1_) = &
          LowDensityRatio*BodyRho_I(1)
  endif

  if(UseElectronPressure)then
     FaceState_VI(Pe_,Body1_) = sum(BodyP_I(IonFirst_:IonLast_))
     FaceState_VI(Pe_,Body2_) = pBody2
  end if

  if(UseAnisoPressure)then
     FaceState_VI(Ppar_,Body1_) = BodyP_I(1)
     FaceState_VI(Ppar_,Body2_) = pBody2
  end if

  !The following part of the code is sensitive to a particular physical
  !model. It should be modified in adding/deleting the physical effects 
  !or features

  FaceState_VI(rho_,body2_) = RhoBody2
  FaceState_VI(P_,body2_)   = pBody2

  !For Outer Boundaries (if SW_* are set)
  if(SW_rho > 0.0)then

     FaceState_VI(Rho_, 1:6) = SW_rho
     FaceState_VI(Ux_,  1:6) = SW_Ux
     FaceState_VI(Uy_,  1:6) = SW_Uy
     FaceState_VI(Uz_,  1:6) = SW_Uz
     FaceState_VI(Bx_,  1:6) = SW_Bx
     FaceState_VI(By_,  1:6) = SW_By
     FaceState_VI(Bz_,  1:6) = SW_Bz
     FaceState_VI(P_,   1:6) = SW_p

     if(UseElectronPressure) FaceState_VI(Pe_, 1:6) = SW_p

     if(UseAnisoPressure) FaceState_VI(Ppar_, 1:6) = SW_p

     if (UseMultiSpecies) then
        FaceState_VI(SpeciesFirst_, 1:6) = &
             SW_rho*(1 - LowDensityRatio * (SpeciesLast_-SpeciesFirst_))
        FaceState_VI(SpeciesFirst_+1:SpeciesLast_, 1:6) = &
             LowDensityRatio*Sw_rho
     endif

     if(nFluid > 1)then
        ! Ratio of total pressure and sum of ion pressures depends on UseElectronPressure
        pCoef = 1 + ElectronPressureRatio
        if(UseElectronPressure) pCoef = 1.0

        iFluid=IonFirst_
        call select_fluid
        FaceState_VI(iRho, 1:6) = &
             SW_Rho*(1-LowDensityRatio*(nFluid-IonFirst_))
        FaceState_VI(iUx,  1:6) = SW_Ux
        FaceState_VI(iUy,  1:6) = SW_Uy
        FaceState_VI(iUz,  1:6) = SW_Uz
        ! Use solar wind temperature and reduced density to get pressure
        FaceState_VI(iP,   1:6) = SW_p/pCoef*(1.0-LowDensityRatio*(nFluid-IonFirst_))

        do iFluid = IonFirst_+1, nFluid
           call select_fluid
           FaceState_VI(iRho, 1:6) = SW_Rho*LowDensityRatio
           FaceState_VI(iUx,  1:6) = SW_Ux
           FaceState_VI(iUy,  1:6) = SW_Uy
           FaceState_VI(iUz,  1:6) = SW_Uz
           ! Use solar wind temperature and reduced density to get pressure 
           FaceState_VI(iP,   1:6) = SW_p/pCoef &
                *LowDensityRatio*MassIon_I(1)/MassFluid_I(iFluid)
        end do
        ! Fix total pressure if necessary (density and temperature are kept)
        ! the min(2,nFluid) instead of 2 is needed so it compiles with gfortran
        if(UseMultiIon .and. IsMhd) FaceState_VI(P_,1:6) = &
             pCoef*sum(FaceState_VI(iP_I(min(2,nFluid):nFluid),1))
     end if
  end if

  ! Cell State is used for filling the ghostcells
  CellState_VI = FaceState_VI

  ! Convert velocity to momentum for all fluids and boundaries
  do iFluid = 1, nFluid
     call select_fluid
     CellState_VI(iRhoUx,:) = FaceState_VI(iUx,:)*FaceState_VI(iRho,:)
     CellState_VI(iRhoUy,:) = FaceState_VI(iUy,:)*FaceState_VI(iRho,:)
     CellState_VI(iRhoUz,:) = FaceState_VI(iUz,:)*FaceState_VI(iRho,:)
  end do

  if(UseOutflowPressure) pOutflow = pOutflowSi*Si2No_V(UnitP_)

  ! Minimum value for densities, pressures, energy
  RhoMin_I     = RhoMinDim_I * Io2No_V(UnitRho_)
  pMin_I       = pMinDim_I * Io2No_V(UnitP_)
  PeMin        = PeMinSi*Si2No_V(UnitP_)
  ExtraEintMin = ExtraEintMinSi*Si2No_V(UnitEnergyDens_)

  ! Relaxation parameters for anisotropic ion pressure
  TauInstability  = TauInstabilitySi*Si2No_V(UnitT_)
  TauGlobal  = TauGlobalSi*Si2No_V(UnitT_)

  !\
  ! Now do the magnetic field stuff
  !/
  ! The *ThetaTilt and Bdp variables are not needed for GM except
  ! for reporting them in write_progress.

  ! Nondimensionalize dipole strength.
  Bdp  = DipoleStrengthSi*Si2No_V(UnitB_)


  BdpBody2_D = BdpDimBody2_D*Io2No_V(UnitB_)              

  ! Saving initial coordinates of second body:
  if(UseOrbit)then
     PhaseBody2    = atan2(ybody2, xbody2)
     DistanceBody2 = sqrt(xbody2**2+ybody2**2)
  end if


  ! Compute dipole tilt variables
  ! For IH ThetaTilt should be set with the #HELIODIPOLETILT command
  ! For static dipole only, e.g. in 2D
  CosThetaTilt = cos(ThetaTilt)
  SinThetaTilt = sin(ThetaTilt)

  ! by default quadrupole and octupole terms are zero
  Qqp(1:3,1:3) = 0.0
  Oop(1:3,1:3,1:3) = 0.0

  ! Hyperbolic cleaning uses SpeedHyp velocity
  if(UseHyperbolicDivb .and. SpeedHypDim > 0)then
     SpeedHyp  = SpeedHypDim*Io2No_V(UnitU_)
     SpeedHyp2 = SpeedHyp**2
  end if

end subroutine set_physics_constants

!===========================================================================

subroutine set_units

  use ModMain
  use ModPhysics
  use ModVarIndexes
  use ModMultiFluid, ONLY: MassIon_I
  use ModUserInterface ! user_io_units, user_normalization
  implicit none

  character (len=*), parameter :: NameSub="set_units"

  logical :: DoTest, DoTestMe
  !-----------------------------------------------------------------------
  call set_oktest(NameSub,DoTest, DoTestMe)

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
     No2Si_V(UnitRho_) = 1000000*cProtonMass !AtomicMass
  case("SOLARWIND")
     ! rPlanet, SW sound speed, SW density in amu/cm^3
     No2Si_V(UnitX_)   = rPlanetSi
     if(NameThisComp=='OH')No2Si_V(UnitX_)=cAU
     No2Si_V(UnitU_)   = sqrt(g*cBoltzmann*SW_T_dim/cProtonMass/MassIon_I(1))
     No2Si_V(UnitRho_) = 1000000*cProtonMass*MassIon_I(1)*SW_n_dim
  case("NONE", "READ")
     ! Already set in set_parameters
  case("USER")
     call user_normalization
  case default
     call stop_mpi(NameSub//' ERROR: unknown TypeNormalization='// &
          trim(TypeNormalization))
  end select

  !\
  ! Set other normalizing SI variables from the independent ones.
  !
  ! For sake of convenience
  !  units of B are chosen to satisfy v_A = B/sqrt(rho)       (mu = 1)
  !  units of n are chosen to satisfy  n  = rho/(ionmass/amu) (mp = 1)
  !  units of T are chosen to satisfy  T  = p/n               (kBoltzmann = 1)
  !
  ! Note that No2Si_V(UnitN_) is NOT EQUAL TO 1/No2Si_V(UnitX_)^3 !!!
  !/
  No2Si_V(UnitT_)          = No2Si_V(UnitX_)/No2Si_V(UnitU_)         ! s
  No2Si_V(UnitN_)          = No2Si_V(UnitRho_)/cProtonMass           ! #/m^3
  No2Si_V(UnitP_)          = No2Si_V(UnitRho_)*No2Si_V(UnitU_)**2    ! Pa
  No2Si_V(UnitB_)          = No2Si_V(UnitU_) &
       *sqrt(cMu*No2Si_V(UnitRho_))                                  ! T
  No2Si_V(UnitRhoU_)       = No2Si_V(UnitRho_)*No2Si_V(UnitU_)       ! kg/m^2/s
  No2Si_V(UnitEnergyDens_) = No2Si_V(UnitP_)                         ! J/m^3
  No2Si_V(UnitPoynting_)   = No2Si_V(UnitEnergyDens_)*No2Si_V(UnitU_)! J/m^2/s
  No2Si_V(UnitJ_)          = No2Si_V(UnitB_)/( No2Si_V(UnitX_)*cMu ) ! A/m^2
  No2Si_V(UnitElectric_)   = No2Si_V(UnitU_)*No2Si_V(UnitB_)         ! V/m
  No2Si_V(UnitTemperature_)= No2Si_V(UnitP_) &
       /( No2Si_V(UnitN_)*cBoltzmann )                               ! K 
  No2Si_V(UnitDivB_)       = No2Si_V(UnitB_)/No2Si_V(UnitX_)         ! T/m
  No2Si_V(UnitAngle_)      = 1.0                                     ! radian
  No2Si_V(UnitMass_)       = No2Si_V(UnitRho_)*No2Si_V(UnitX_)**3    ! kg
  No2Si_V(UnitCharge_)     = No2Si_V(UnitJ_)/No2Si_V(UnitU_) &       ! C
       *No2Si_V(UnitX_)**3
  No2Si_V(UnitUnity_)      = 1.0                                     ! Fallback conversion for undefined units

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
     Io2Si_V(UnitDivB_)     = 1.0E-9/rPlanetSi                ! nT/R_planet
     Io2Si_V(UnitRho_)      = 1.0E6*cProtonMass               ! Mp/cm^3
     Io2Si_V(UnitN_)        = 1.0E6                           ! #/cm^3
     Io2Si_V(UnitU_)        = 1.0E3                           ! km/s
     Io2Si_V(UnitP_)        = 1.0E-9                          ! nPa
     Io2Si_V(UnitB_)        = 1.0E-9                          ! nT
     Io2Si_V(UnitJ_)        = 1.0E-6                          ! microA/m^2
     Io2Si_V(UnitElectric_) = 1.0E-3                          ! mV/m
     Io2Si_V(UnitAngle_)    = cRadToDeg                       ! degrees
     !\
     ! set string variables used for writing output - TECPLOT
     !/
     NameTecUnit_V(UnitX_)           = '['//trim(NamePlanetRadius)//']'
     NameTecUnit_V(UnitRho_)         = '[amu/cm^3]'
     NameTecUnit_V(UnitU_)           = '[km/s]'          
     NameTecUnit_V(UnitN_)           = '[cm^-^3]'        
     NameTecUnit_V(UnitP_)           = '[nPa]'           
     NameTecUnit_V(UnitB_)           = '[nT]'            
     NameTecUnit_V(UnitJ_)           = '[`mA/m^2]'       
     NameTecUnit_V(UnitElectric_)    = '[mV/m]'          
     NameTecUnit_V(UnitDivB_)        = '[nT/'//trim(NamePlanetRadius)//']'
     NameTecUnit_V(UnitAngle_)       = '[deg]'

     !\
     ! set string variables used for writing output - IDL
     !/
     NameIdlUnit_V(UnitX_)           = NamePlanetRadius
     NameIdlUnit_V(UnitRho_)         = 'Mp/cc'
     NameIdlUnit_V(UnitU_)           = 'km/s'
     NameIdlUnit_V(UnitN_)           = '/cc'
     NameIdlUnit_V(UnitP_)           = 'nPa'           
     NameIdlUnit_V(UnitB_)           = 'nT'
     NameIdlUnit_V(UnitJ_)           = 'uA/m2'
     NameIdlUnit_V(UnitElectric_)    = 'mV/m'
     NameIdlUnit_V(UnitDivB_)        = 'nT/'//NamePlanetRadius
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

  use ModProcMH,  ONLY: iProc
  use ModVarIndexes
  use ModPhysics
  use ModMultiFluid
  use ModAdvance, ONLY: UseElectronPressure, UseAnisoPressure, UseIdealEos
  use ModMain,    ONLY: UseB
  use ModUtilities, ONLY: lower_case
  implicit none

  integer :: iVar
  character (len=len(NameVar_V)) :: NameVar
  character (len=*), parameter :: NameSub="init_mhd_variables"
  !--------------------------------------------------------------------------

  ! Set mapping from state variable indices to unit conversion array indices
  do iVar = 1, nVar
     NameVar = NameVar_V(iVar)
     call lower_case(NameVar)
     call extract_fluid_name(NameVar)
     select case(NameVar)
     case('rho', &
          'h', 'hp', 'hpsw', 'h2p', 'sw', 'iono', &
          'o', 'op', 'o2p',  'he', 'ohp', 'n',  &
          'cop', 'co2p',  'h2O', 'h2Op', 'h3Op',  &
          'mp', 'lp', 'mhcp', 'hhcp', 'hnip', &
          'neu1', 'neu2', 'neu3', 'neu4', &
          'pui1', 'pui2', 'pui3', 'pui4')
        iUnitCons_V(iVar) = UnitRho_
        iUnitPrim_V(iVar) = UnitRho_
     case('mx','my','mz','ux','uy','uz')
        iUnitCons_V(iVar) = UnitRhoU_
        iUnitPrim_V(iVar) = UnitU_
     case('bx','by','bz')
        iUnitCons_V(iVar) = UnitB_
        iUnitPrim_V(iVar) = UnitB_
     case('p','pe', 'ppar')
        iUnitCons_V(iVar) = UnitP_
        iUnitPrim_V(iVar) = UnitP_
     case('e', 'ew', 'ehot', 'eint')
        iUnitCons_V(iVar) = UnitEnergyDens_
        iUnitPrim_V(iVar) = UnitEnergyDens_
     case default
        if(WaveFirst_ >= iVar .and. iVar <= WaveLast_)then
           iUnitCons_V(iVar) = UnitEnergyDens_
           iUnitPrim_V(iVar) = UnitEnergyDens_
        else
           if(iProc == 0) write(*,*) NameSub,': ',NameVar, &
                ' was not found to set iUnitState_V,', &
                ' using UnitUnity_ instead.'
           iUnitCons_V(iVar) = UnitUnity_
           iUnitPrim_V(iVar) = UnitUnity_
        end if
     end select
  end do

  if(UseB)then
     UnitUser_V(Bx_:Bz_)        = No2Io_V(UnitB_)
     NameUnitUserTec_V(Bx_:Bz_) = NameTecUnit_V(UnitB_)
     NameUnitUserIdl_V(Bx_:Bz_) = NameIdlUnit_V(UnitB_)
  end if

  do iFluid = 1, nFluid
     call select_fluid
     UnitUser_V(iRho)          = No2Io_V(UnitRho_)
     UnitUser_V(iRhoUx:iRhoUz) = No2Io_V(UnitRhoU_)
     UnitUser_V(iP)            = No2Io_V(UnitP_)
     UnitUser_V(iEnergy)       = No2Io_V(UnitEnergyDens_)

     NameUnitUserTec_V(iRho)          = NameTecUnit_V(UnitRho_)
     NameUnitUserTec_V(iRhoUx:iRhoUz) = NameTecUnit_V(UnitRhoU_)
     NameUnitUserTec_V(iP)            = NameTecUnit_V(UnitP_)
     NameUnitUserTec_V(iEnergy)       = NameTecUnit_V(UnitEnergyDens_)

     NameUnitUserIdl_V(iRho)          = NameIdlUnit_V(UnitRho_)
     NameUnitUserIdl_V(iRhoUx:iRhoUz) = NameIdlUnit_V(UnitRhoU_)
     NameUnitUserIdl_V(iP)            = NameIdlUnit_V(UnitP_)
     NameUnitUserIdl_V(iEnergy)       = NameIdlUnit_V(UnitEnergyDens_)
  end do

  ! By default the scalar advected variables are assumed to behave like density
  do iVar = ScalarFirst_, ScalarLast_
     UnitUser_V(iVar)        = No2Io_V(UnitRho_)
     NameUnitUserTec_V(iVar) = NameTecUnit_V(UnitRho_)
     NameUnitUserIdl_V(iVar) = NameIdlUnit_V(UnitRho_)
  end do

  if(UseElectronPressure)then
     UnitUser_V(Pe_)        = No2Io_V(UnitP_)
     NameUnitUserTec_V(Pe_) = NameTecUnit_V(UnitP_)
     NameUnitUserIdl_V(Pe_) = NameIdlUnit_V(UnitP_)
  end if

  if(UseAnisoPressure)then
     UnitUser_V(Ppar_)        = No2Io_V(UnitP_)
     NameUnitUserTec_V(Ppar_) = NameTecUnit_V(UnitP_)
     NameUnitUserIdl_V(Ppar_) = NameIdlUnit_V(UnitP_)
  end if

  if(WaveLast_ > 1)then
     ! Set the unit and unit name for the wave energy variable
     UnitUser_V(WaveFirst_:WaveLast_)        = No2Io_V(UnitEnergyDens_)
     NameUnitUserTec_V(WaveFirst_:WaveLast_) = NameTecUnit_V(UnitEnergyDens_)
     NameUnitUserIdl_V(WaveFirst_:WaveLast_) = NameIdlUnit_V(UnitEnergyDens_)
  end if

  if(Ew_ > 1)then
     ! Set the unit and unit name for the total wave energy variable
     UnitUser_V(Ew_)        = No2Io_V(UnitEnergyDens_)
     NameUnitUserTec_V(Ew_) = NameTecUnit_V(UnitEnergyDens_)
     NameUnitUserIdl_V(Ew_) = NameIdlUnit_V(UnitEnergyDens_)
  end if

  if(.not.UseIdealEos)then
     UnitUser_V(ExtraEint_)        = No2Io_V(UnitEnergyDens_)
     NameUnitUserTec_V(ExtraEint_) = NameTecUnit_V(UnitEnergyDens_)
     NameUnitUserIdl_V(ExtraEint_) = NameIdlUnit_V(UnitEnergyDens_)
  end if

  if(Hyp_ > 1)then
     ! Set the scalar field Phi used in hyperbolic cleaning
     UnitUser_V(Hyp_) = No2Io_V(UnitB_)*No2Io_V(UnitU_)
     NameUnitUserTec_V(Hyp_) = &
          trim(NameTecUnit_V(UnitB_)) // trim(NameTecUnit_V(UnitU_))
    
     NameUnitUserIdl_V(Hyp_) = &
          trim(NameIdlUnit_V(UnitB_)) // trim(NameIdlUnit_V(UnitU_))
  end if

  if(SignB_ > 1)then
     UnitUser_V(SignB_)        = 1.0
     NameUnitUserTec_V(SignB_) = ''
     NameUnitUserIdl_V(SignB_) = '1'
  end if

  if(Ehot_ > 1)then
     ! Set the unit and unit name for the extra internal energy due
     ! to suprathermal electrons
     UnitUser_V(Ehot_)        = No2Io_V(UnitEnergyDens_)
     NameUnitUserTec_V(Ehot_) = NameTecUnit_V(UnitEnergyDens_)
     NameUnitUserIdl_V(Ehot_) = NameIdlUnit_V(UnitEnergyDens_)
  end if

end subroutine init_mhd_variables
