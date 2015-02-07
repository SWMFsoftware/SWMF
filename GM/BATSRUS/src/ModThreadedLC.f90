!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModThreadedLC
  use ModFieldLineThread, ONLY: &
       BoundaryThreads, BoundaryThreads_B, AWHeating_, APlusBC_, &
       LengthPAvrSi_, UHeat_, HeatFluxLength_, DHeatFluxXOverU_, &
       RadCool2Si, Use1DModel, DoInit_, Done_, Enthalpy_, Heat_

  use ModCoronalHeating, ONLY:PoyntingFluxPerBSi, PoyntingFluxPerB, &
                              LPerpTimesSqrtBSi, LPerpTimesSqrtB

  use ModAdvance,    ONLY: UseElectronPressure, UseIdealEos
  use ModPhysics,    ONLY: AverageIonCharge
  use ModConst,      ONLY: rSun, mSun, cBoltzmann, cAtomicMass, cGravitation
  implicit none

  real :: TeFraction, TiFraction, PeFraction
  real,allocatable :: Te_G(:,:,:)
  real,allocatable,dimension(:)::ReflCoef_I, APlus_I, AMinus_I, &
       TeSi_I, Cons_I, PeSi_I, Res_I, U_I, L_I, M_I, Xi_I, &
       VaLog_I,  DCons_I, DXi_I, ResHeating_I, ResCooling_I, ResEnthalpy_I,&
       ResHeatCond_I, SpecHeat_I, IntEnergy_I
  

  real, parameter:: TeSiMin = 3.5e5    ![K]
  real           :: TeMin, ConsMin

  real :: cCoolingPerPe2, RhoNoDimCoef
  !\
  ! For transforming conservative to TeSi and back 
  !/
  real, parameter:: cTwoSevenths = 2.0/7.0
  !\
  ! Gravitation potential
  !/
  real, parameter:: cGravPot = cGravitation*mSun*cAtomicMass/&
       (cBoltzmann*rSun)
  !\
  ! The use of this constant:
  !/
  !\
  !   Hydrostatic equilibrium in an isothermal corona: 
  !    d(N_i*k_B*(Z*T_e +T_i) )/dr=G*M_sun*N_I*M_i*d(1/r)/dr
  ! => N_i*Te\propto exp(cGravPot/TeSi*(M_i[amu]/(1+Z))*\Delta(R_sun/r)) 
  !/
  !\
  ! The plasma properties dependent coefficient needed to evaluate the 
  ! eefect of gravity on the hydrostatic equilibrium
  !/
  real    :: GravHydroStat != cGravPot*MassIon_I(1)/(AverageIonCharge + 1)
  !\
  !energy flux needed to raise the mass flux rho*u to the heliocentric 
  !distance r equals: rho*u*G*Msun*(1/R_sun -1/r)=
  !=k_B*N_i*M_i(amu)u*cGravPot*(1-R_sun/r)=
  !=P_e/T_e*cGravPot*u(M_i[amu]/Z)*(1/R_sun -1/r)

  integer:: iP
  integer, parameter:: Impl_=3
contains
  !=========================================================================
  subroutine init_threaded_lc
    use BATL_lib, ONLY:  MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    use ModMultiFluid,   ONLY: MassIon_I
    use ModFieldLineThread, ONLY: check_tr_table, & 
         nPointThreadMax, HeatCondParSi
    use ModPhysics,            ONLY: UnitTemperature_, Si2No_V, UnitEnergyDens_
    use ModVarIndexes,         ONLY: Pe_, p_
    !-------------------
    allocate(Te_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)); Te_G = 0.0

    allocate(ReflCoef_I(0:nPointThreadMax)); ReflCoef_I = 0.0
    allocate(   APlus_I(0:nPointThreadMax));    APlus_I = 0.0
    allocate(  AMinus_I(0:nPointThreadMax));   AMinus_I = 0.0

    allocate(   TeSi_I(nPointThreadMax));     TeSi_I = 0.0
    allocate(   Cons_I(nPointThreadMax));     Cons_I = 0.0
    allocate( PeSi_I(nPointThreadMax));   PeSi_I = 0.0
    allocate(    Res_I(nPointThreadMax));      Res_I = 0.0
    allocate(   SpecHeat_I(nPointThreadMax));   SpecHeat_I = 0.0
    allocate(  IntEnergy_I(nPointThreadMax));  IntEnergy_I = 0.0

    allocate( ResHeating_I(nPointThreadMax)); ResHeating_I = 0.0
    allocate( ResCooling_I(nPointThreadMax)); ResCooling_I = 0.0
    allocate(ResEnthalpy_I(nPointThreadMax));ResEnthalpy_I = 0.0
    allocate(ResHeatCond_I(nPointThreadMax));ResHeatCond_I = 0.0

    allocate(      U_I(nPointThreadMax));        U_I = 0.0
    allocate(      L_I(nPointThreadMax));        L_I = 0.0
    allocate(      M_I(nPointThreadMax));        M_I = 0.0
    allocate(     Xi_I(0:nPointThreadMax));     Xi_I = 0.0
    allocate(  VaLog_I(nPointThreadMax));    VaLog_I = 0.0
    allocate(  DCons_I(nPointThreadMax));    DCons_I = 0.0
    allocate(    DXi_I(nPointThreadMax));      DXi_I = 0.0

    !\
    ! TeFraction is used for ideal EOS:
    !/
    if(UseElectronPressure)then
       ! Pe = ne*Te (dimensionless) and n=rho/ionmass
       ! so that Pe = ne/n *n*Te = (ne/n)*(rho/ionmass)*Te
       ! TeFraction is defined such that Te = Pe/rho * TeFraction
       TiFraction = MassIon_I(1)
       TeFraction = MassIon_I(1)/AverageIonCharge
       iP = Pe_
       PeFraction = 1.0
    else
       ! p = n*T + ne*Te (dimensionless) and n=rho/ionmass
       ! so that p=rho/massion *T*(1+ne/n Te/T)
       ! TeFraction is defined such that Te = p/rho * TeFraction
       TiFraction = MassIon_I(1) &
            /(1 + AverageIonCharge)
       TeFraction = TiFraction
       iP = p_
       PeFraction = AverageIonCharge/(1.0 + AverageIonCharge) 
    end if

    !\
    ! Therefore Te = TeFraction*State_V(iP)/State_V(Rho_)
    ! Pe = PeFraction*State_V(iP)
    !/

    TeMin = TeSiMin*Si2No_V(UnitTemperature_)
    call check_tr_table
    !\
    !   Hydrostatic equilibrium in an isothermal corona: 
    !    d(N_i*k_B*(Z*T_e +T_i) )/dr=G*M_sun*N_I*M_i*d(1/r)/dr
    ! => N_i*Te\propto exp(cGravPot/TeSi*(M_i[amu]/(1+Z))*\Delta(R_sun/r)) 
    !/
    GravHydroStat = cGravPot*MassIon_I(1)/(AverageIonCharge + 1)

    ConsMin = cTwoSevenths*HeatCondParSi*TeSiMin**3.50
    cCoolingPerPe2 = RadCool2Si/(cBoltzmann*cBoltzmann*AverageIonCharge)
    RhoNoDimCoef = Si2No_V(UnitEnergyDens_)/PeFraction*&
         TeFraction/Si2No_V(UnitTemperature_)
  end subroutine init_threaded_lc
  !================================
  !\
  ! Main routine:
  ! solves MHD equations along thread, to link the state above the
  ! inner boundary of the global solar corona to the photosphere
  !/
  subroutine solve_boundary_thread(j, k, iBlock, &
       iAction, TeSiIn, PeSiIn, USiIn, AMinorIn, &
       DTeOverDsSiOut, PeSiOut, RhoNoDimOut, AMajorOut)
    !\
    ! USE:
    !/
    use ModFieldLineThread, ONLY: HeatCondParSi, solve_a_plus_minus
    use ModPhysics,      ONLY: inv_gm1, No2Si_V, UnitX_,Si2No_V, &
                               UnitB_
    use ModMultiFluid,   ONLY: MassIon_I
    use ModLookupTable,  ONLY: i_lookup_table, interpolate_lookup_table
    use ModMain,         ONLY: BlkTest, ProcTest, jTest, kTest
    use ModProcMH,       ONLY: iProc
    !\
    !The initial version: pressure is constant along the thread,
    !reflection and dissipation of the major wave is ignored
    !/
    !INPUT:
    !\
    !Cell and block indexes for the boundary point
    !/
    integer,intent(in):: j, k, iBlock, iAction 
    !\
    ! Parameters of the state in the true cell near the boundary:
    ! TeSiIn: Temperature in K 
    ! USiIn:  Velocity progection on the magnetic field direction.
    ! It is positive if the wind blows outward the Sun.
    ! AMinorIn: for the wave propagating toward the Sun 
    !            EnergyDensity = (\Pi/B)\sqrt{\rho} AMinor**2 
    !/
    real,   intent(in):: TeSiIn, PeSiIn, USiIn, AMinorIn
    !\
    !OUTPUT:
    !DTeOverDsSiOut: Temperature derivative along the thread, at the end point
    !                Used to find the electron temperature in the ghostcell
    !PeSiOut: The electron pressure 
    !AMajorOut: For the wave propagating outward the Sun
    !            EnergyDensity = (\Pi/B)\sqrt{\rho} AMajor**2  
    !/
    real,  intent(out):: DTeOverDsSiOut, PeSiOut, RhoNoDimOut, AMajorOut

    !\
    ! Arrays and table numbers needed to use lookup table
    !/ 
    real    :: Value_V(4), AWValue_V(2), ValCooling(1)
    integer :: iTable, iTableAW, iTableRadcool
    !\
    !---------Used in 0D analytical model-----------------------
    !/
    !\
    !Dimmensionless length (related to the wave dissipation length)
    !/
    real    :: AWLength
    !\
    ! Total contribution from the AW heating to the TR energetics
    !/
    real    :: AWHeating
    !\
    !---------Used in 1D numerical model------------------------
    !/
    !\
    ! Number of TEMPERATURE nodes (first one is on the top of TR
    ! the last one is in the physical cell of the SC model
    !/
    integer :: nPoint 
    integer        :: nIter=5
    !\
    ! Enthalpy correction coefficient
    !/
    real    :: FluxConst, EnthalpyCorrection, DEnthalpyCorrOverDU
    !\
    ! Limited USiIn
    !/
    real :: USiLtd

    logical :: DoTest, DoTestMe

    character(len=*), parameter :: NameSub = 'solve_boundary_thread'
    !-------------------------------------------------------------------------
    if(iBlock==BLKtest.and.iProc==PROCtest.and.j==jTest.and.k==kTest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
       if(DoTestMe)then
          write(*,*)'TeSiIn=       ',TeSiIn,' K '
          write(*,*)'PeSiIn = ', PeSiIn
          write(*,*)'NeSiIn = ', PeSiIn/(TeSiIn*cBoltzmann)
          write(*,*)'AMinorIn=     ', AMinorIn
          write(*,*)'USiIn=        ',USiIn,' m/s'
          write(*,*)'Thread Length=', &
               BoundaryThreads_B(iBlock)% LengthSi_III(0,j,k) &
               ,' m = ',  Si2No_V(UnitX_)*&
               BoundaryThreads_B(iBlock)% LengthSi_III(0,j,k),' Rs'
          write(*,*)'iAction=',iAction
       end if
    else
       DoTest=.false.; DoTestMe=.false.
    endif
    iTable = i_lookup_table('TR')
    if(iTable<=0)call CON_stop('TR table is not set')

    call interpolate_lookup_table(iTable, TeSiIn, 1.0e8, Value_V, &
         DoExtrapolate=.false.)
    USiLtd = sign(min(abs(USiIn),0.05*Value_V(UHeat_)),USiIn)
    !\
    ! First value is now the product of the thread length in meters times
    ! a geometric mean pressure, so that
    !/
    PeSiOut = Value_V(LengthPAvrSi_)*sqrt(AverageIonCharge)/&
         BoundaryThreads_B(iBlock)% LengthSi_III(0,j,k)
    if(DoTestMe)write(*,*)'Pressure 0D (SI) = ',PeSiOut
    if(.not.Use1DModel)then
       RhoNoDimOut = RhoNoDimCoef*PeSiIn/TeSiIn
       
       !Dimmensionless length (related to the wave dissipation length)
       AWLength = BoundaryThreads_B(iBlock)% DXi_III(0,j,k)*&
            sqrt(sqrt(RhoNoDimOut))
       !\
       !Calculate Alfven waves for the given thread length and BC 
       !for ingoing wave
       !/ 
       iTableAW = i_lookup_table('AW_TR')
       if(iTableAW <=0 )call CON_stop('AW_TR table is not set')
       call interpolate_lookup_table(iTableAW, AWLength, AMinorIn, AWValue_V, &
            DoExtrapolate=.false.)

       AWHeating = AWValue_V(AWHeating_)*PoyntingFluxPerBSi*&
            BoundaryThreads_B(iBlock)% B_III(0,j,k)*No2Si_V(UnitB_)

       !\
       ! Heat flux equals PAvr * UHeat (spent for radiation) +
       ! Pi * U * (5/2) + Pe * U * (5/2) (carried away by outflow), 
       ! the pressure ratio being  Pe = Pi * AverageIonCharge
       ! 
       ! Temperature gradient equals the heat flux divided by kappa0*T**2.5
       !/
       DTeOverDsSiOut = ( PeSiOut/sqrt(AverageIonCharge) * &
            Value_V(UHeat_)                           -  & !Radiation losses
            AWHeating                                    & !AW Heating
            +USiIn * (PeSiOut/AverageIonCharge)          & !5/2*U*(Pi+Pe)
            *(inv_gm1 +1)*(1 + AverageIonCharge) )/      &
            (HeatCondParSi * TeSiIn**2.50)    

       AMajorOut = AWValue_V(APlusBC_)
       if(DoTestMe)then
          write(*,*)'Dimensionless length characteristic of AW dissipation=',&
            AWLength
          write(*,*)'AMajorOut=    ', AMajorOut
          write(*,*)'Final dT/ds=  ',DTeOverDsSiOut,&
               ', dT/ds*Length=',DTeOverDsSiOut*&
               BoundaryThreads_B(iBlock)% LengthSi_III(0,j,k),' K'
       end if
       RETURN
    end if
    !Full 1D model
    nPoint = BoundaryThreads_B(iBlock)% nPoint_II(j,k)
    if(iAction/=DoInit_)then
       !\
       ! Retrieve temperature and pressure distribution
       !/
       TeSi_I(1:nPoint) = BoundaryThreads_B(iBlock)%TSi_III(1-nPoint:0,j,k)
       PeSi_I(1:nPoint) = BoundaryThreads_B(iBlock)%PSi_III(1-nPoint:0,j,k)
       TeSi_I(nPoint)   = TeSiIn
    end if

    iTableRadCool = i_lookup_table('radcool')
    call set_initial_thread

    !\
    ! Outputs
    !/
    DTeOverDsSiOut = (TeSi_I(nPoint) - TeSi_I(nPoint-1))/&
         (BoundaryThreads_B(iBlock)% LengthSi_III(0,j,k) - &
         BoundaryThreads_B(iBlock)% LengthSi_III(-1,j,k))
    PeSiOut = PeSi_I(nPoint)
    RhoNoDimOut = RhoNoDimCoef*PeSi_I(nPoint)/TeSi_I(nPoint)
    AMajorOut = AWValue_V(APlusBC_)
    if(DoTestMe)then
       write(*,*)'Pressure 1D (SI) = ',PeSiOut
       write(*,*)'AMajorOut=    ', AMajorOut
       write(*,*)'Final dT/ds=  ',DTeOverDsSiOut,&
            ', dT/ds*Length=',DTeOverDsSiOut*&
            BoundaryThreads_B(iBlock)% LengthSi_III(0,j,k),' K'
    end if
  contains
    !=======================
    subroutine get_res_heating(nIterIn)
      integer, intent(in)::nIterIn
      integer:: iPoint
      real    :: SqrtRho, RhoNoDim
      !--------------------
      !\
      ! Now prepare to calculate Alfven wave amplitude distribution
      !/
      Xi_I(0) = 0.0
      do iPoint=1,nPoint
         !\
         ! 1. Calculate sqrt(RhoNoDim)
         !/
         RhoNoDim = RhoNoDimCoef*PeSi_I(iPoint)/TeSi_I(iPoint)
         if(RhoNoDim<=0.0)then
            write(*,*)'iPoint=',iPoint,' PeSi_I(iPoint)=',&
                 PeSi_I(iPoint),' TeSi_I(iPoint)=',TeSi_I(iPoint)
            call CON_stop('Non-positive Density in '//NameSub)
         end if
         SqrtRho = sqrt(RhoNoDim)
         !\
         ! 2. Calculate Alfven wave speed
         !/
         VaLog_I(iPoint) = &
            log(BoundaryThreads_B(iBlock)% B_III(iPoint-nPoint,j,k)/SqrtRho)
         !\
         ! 3. Calculate dimensionless length (in terms of the dissipation length
         !/
         DXi_I(iPoint) = &
              (BoundaryThreads_B(iBlock)% DXi_III(iPoint-nPoint,j,k) - &
              BoundaryThreads_B(iBlock)% DXi_III(iPoint-1-nPoint,j,k))*&
              sqrt(SqrtRho)
         Xi_I(iPoint) = Xi_I(iPoint-1) + DXi_I(iPoint) 

      end do
    
      !\
      ! 4. Calculate the reflection coefficient
      !/
      ReflCoef_I(1) = abs(VaLog_I(2) - VaLog_I(1))/(0.50*DXi_I(2) + DXi_I(1))
      ReflCoef_I(0) =  ReflCoef_I(1)
      do iPoint = 2, nPoint-2
         ReflCoef_i(iPoint) = &
              abs(VaLog_I(iPoint+1) - VaLog_I(iPoint))/&
              (0.50*(DXi_I(iPoint+1) +  DXi_I(iPoint)))
      end do
      ReflCoef_I(nPoint-1) = abs(VaLog_I(nPoint) - VaLog_I(nPoint-1))/&
           (0.50*DXi_I(nPoint-1) + DXi_I(nPoint))
      ReflCoef_I(nPoint) =     ReflCoef_I(nPoint-1)
  
      !\
      ! Solve amplitudes of the Alfven waves (arrays there have dimension)
      ! (0:nI)
      !/
      call solve_a_plus_minus(&
           nI=nPoint,                      & 
           ReflCoef_I=ReflCoef_I(0:nPoint),&
           Xi_I=Xi_I(0:nPoint),            &
           AMinusBC=AMinorIn,              &
           Heating=AWValue_V(AWHeating_),  &
           APlusBC=AWValue_V(APlusBC_),    &
           APlusOut_I=APlus_I(0:nPoint),   &
           AMinusOut_I=AMinus_I(0:nPoint), &
           nIterIn=nIterIn)

      ResHeating_I = 0.0
      ResHeating_I(1:nPoint-1) = &
           (APlus_I(0:nPoint-2)**2 - AMinus_I(0:nPoint-2)**2)&
           -(APlus_I(1:nPoint-1  )**2 - AMinus_I(1:nPoint-1)**2)
    end subroutine get_res_heating
    !=============================
    subroutine set_pressure
      integer::iPoint
      !---------------
      !\
      ! First value is now the product of the thread length in meters times
      ! a geometric mean pressure, so that
      !/
      PeSi_I(1) = Value_V(LengthPAvrSi_)*sqrt(AverageIonCharge)/&
           BoundaryThreads_B(iBlock)% LengthSi_III(1-nPoint,j,k)
      !\
      !   Hydrostatic equilibrium in an isothermal corona: 
      !    d(N_i*k_B*(Z*T_e +T_i) )/dr=G*M_sun*N_I*M_i*d(1/r)/dr
      ! => N_i*Te\propto exp(cGravPot/TeSi*(M_i[amu]/(1+Z))*\Delta(R_sun/r)) 
      !/
      !GravHydroStat = cGravPot*MassIon_I(1)/(AverageIonCharge + 1)
      do iPoint = 2, nPoint
         PeSi_I(iPoint) = PeSi_I(iPoint-1)*&
            exp( (0.5/TeSi_I(iPoint) + 0.5/TeSi_I(iPoint-1))*&
            GravHydroStat*&
            (BoundaryThreads_B(iBlock)%RInv_III(iPoint-nPoint,j,k) - &
            BoundaryThreads_B(iBlock)%RInv_III(iPoint-1-nPoint,j,k)))
      end do
    end subroutine set_pressure
    !==========================
    subroutine get_heat_cond
      !\
      ! The number of points to solve temperature is nPoint - 1
      !/
      U_I(1:nPoint-1) = -BoundaryThreads_B(iBlock)% BDsInv_III(1-nPoint:-1,j,k)
      L_I(2:nPoint-1) = U_I(1:nPoint-2)
      M_I(2:nPoint-1) = -U_I(2:nPoint-1) - L_I(2:nPoint-1)
      M_I(1) = -U_I(1) + Value_V(DHeatFluxXOverU_)*&
           BoundaryThreads_B(iBlock)% BDsInv_III(-nPoint,j,k)
      ResHeatCond_I(1:nPoint-1) = &
           ( Cons_I(1:nPoint-1) - Cons_I(2:nPoint))*U_I(1:nPoint-1)
      !\
      ! 4. Add left heat flux to the TR
      !/
      ResHeatCond_I(1) = ResHeatCond_I(1) - Value_V(HeatFluxLength_)*&
           BoundaryThreads_B(iBlock)% BDsInv_III(-nPoint,j,k)
      !\
      ! 5. Add left heat fluxes
      !/
      ResHeatCond_I(2:nPoint-1) = ResHeatCond_I(2:nPoint-1) &
           +( Cons_I(2:nPoint-1) - Cons_I(1:nPoint-2))*L_I(2:nPoint-1) 
    end subroutine get_heat_cond
    !=============
    subroutine get_rad_cooling
      integer:: iPoint
      !-------------------
      ResCooling_I = 0.0;
      do iPoint = 1, nPoint-1
         call interpolate_lookup_table(iTableRadCool,&
              TeSi_I(iPoint), 1.0e2, ValCooling)
         ResCooling_I(iPoint) = &
              -BoundaryThreads_B(iBlock)%DsOverB_III(iPoint-nPoint,j,k)&
              *ValCooling(1)*cCoolingPerPe2*&
              (PeSi_I(iPoint)/TeSi_I(iPoint))**2
      end do
    end subroutine get_rad_cooling
    !======================
    subroutine set_initial_thread
      !\
      ! Loop variable
      !/
      
      integer :: iPoint
      integer        :: iIter
      !-----------------    
      !\
      ! As a first approximation, recover Te from the analytical solution
      !/
      TeSi_I(nPoint) = TeSiIn
      do iPoint = nPoint-1, 1, -1
         call interpolate_lookup_table(&
              iTable=iTable,           &
              iVal=LengthPAvrSi_,      &
              ValIn=PeSiOut/sqrt(AverageIonCharge)*                       &
              BoundaryThreads_B(iBlock)% LengthSi_III(iPoint-nPoint,j,k), &
              Arg2In=1.0e8,            &
              Value_V=Value_V,         &
              Arg1Out=TeSi_I(iPoint),  & 
              DoExtrapolate=.false.)
      end do
      TeSi_I(1:nPoint) = max(TeSiMin, TeSi_I(1:nPoint))
      Cons_I(1:nPoint) = cTwoSevenths*HeatCondParSi*TeSi_I(1:nPoint)**3.50
      !\
      ! The analytical solution assumes constant pressure and 
      ! no heating. Calculate the pressure distribution
      !/ 
      call set_pressure
      do iIter = 1,nIter
         !\
         !Shape the source.
         !/
         call get_res_heating(nIterIn=nIter+iIter)
         call get_heat_cond
         call get_rad_cooling
         DCons_I = 0.0;  ResEnthalpy_I = 0.0
         !\
         ! Add enthalpy correction
         !/
         if(USiLtd>0)then
            FluxConst = USiLtd * (PeSi_I(nPoint)/AverageIonCharge)& !5/2*U*Pi
                 *(inv_gm1 +1)*(1 + AverageIonCharge)/&
                 (TeSiIn*PoyntingFluxPerBSi*&
                 BoundaryThreads_B(iBlock)% B_III(0,j,k)*No2Si_V(UnitB_))
            do iPoint = 1, nPoint-2
               EnthalpyCorrection = FluxConst*TeSi_I(iPoint)
               ResEnthalpy_I(iPoint)   = &
                    ResEnthalpy_I(iPoint)   - EnthalpyCorrection
               ResEnthalpy_I(iPoint+1) = &
                    ResEnthalpy_I(iPoint+1) + EnthalpyCorrection
            end do
            EnthalpyCorrection = FluxConst*TeSi_I(nPoint-1)
            ResEnthalpy_I(nPoint-1)   = &
                 ResEnthalpy_I(nPoint-1)  - EnthalpyCorrection
         elseif(USiLtd<0)then
            FluxConst = USiLtd * (PeSiIn/AverageIonCharge)  & !5/2*U*Pi
                 *(inv_gm1 +1)*(1 + AverageIonCharge)/&
                 (TeSiIn*PoyntingFluxPerBSi*&
                 BoundaryThreads_B(iBlock)% B_III(0,j,k)*No2Si_V(UnitB_))
            do iPoint = 2, nPoint-1
               EnthalpyCorrection = FluxConst*TeSi_I(iPoint)
               ResEnthalpy_I(iPoint)   = &
                    ResEnthalpy_I(iPoint)   + EnthalpyCorrection
               ResEnthalpy_I(iPoint-1) = &
                    ResEnthalpy_I(iPoint-1) - EnthalpyCorrection
            end do
            EnthalpyCorrection = FluxConst*TeSi_I(nPoint)
            ResEnthalpy_I(nPoint-1)   = &
                 ResEnthalpy_I(nPoint-1)  - EnthalpyCorrection
         end if
         
         !\
         !cGravPot = cGravitation*mSun*cAtomicMass/&
         !/   (cBoltzmann*rSun)
         !\
         !energy flux needed to raise the mass flux rho*u to the heliocentric 
         !distance r equals: rho*u*G*Msun*(1/R_sun -1/r)=
         !=k_B*N_i*M_i(amu)*u*cGravPot*(1-R_sun/r)=
         !=P_e/T_e*cGravPot*(M_sun[amu]/Z)*u*(1/R_sun -1/r)
         !/
         !GravityCoef =  cGravPot/TeSiIn*MassIon_I(1)*          & 
         !     (1 - BoundaryThreads_B(iBlock)%RInv_III(0,j,k) )
         
         Res_I(1:nPoint-1) = &
              ResHeating_I(1:nPoint-1) +  ResCooling_I(1:nPoint-1) +&
              ResEnthalpy_I(1:nPoint-1) + ResHeatCond_I(1:nPoint-1)
         !\
         ! Solve equation
         ! 
         !/
         !=====================================================================!
         ! This routine solves three-diagonal system of equations:             !
         !  ||m_1 u_1  0....        || ||w_1|| ||r_1||                         !
         !  ||l_2 m_2 u_2...        || ||w_2|| ||r_2||                         !
         !  || 0  l_3 m_3 u_3       ||.||w_3||=||r_3||                         !
         !  ||...                   || ||...|| ||...||                         !
         !  ||.............0 l_n m_n|| ||w_n|| ||r_n||                         !
         ! From: Numerical Recipes, Chapter 2.6, p.40.
         ! I do not remember, why it is misspelled tridag instead of tridiag   !
         !=====================================================================!
         call tridag(n=nPoint-1,  &
              L_I=L_I(1:nPoint-1),&
              M_I=M_I(1:nPoint-1),&
              U_I=U_I(1:nPoint-1),&
              R_I=Res_I(1:nPoint-1),&
              W_I=DCons_I(1:nPoint-1))
         
         Cons_I(1:nPoint-1) = &
              max(ConsMin,Cons_I(1:nPoint-1) + DCons_I(1:nPoint-1))
         TeSi_I(1:nPoint-1) = &
              (3.50*Cons_I(1:nPoint-1)/HeatCondParSi)**cTwoSevenths
         
         call interpolate_lookup_table(iTable, TeSi_I(1), 1.0e8, Value_V, &
              DoExtrapolate=.false.)
         call set_pressure
      end do
    end subroutine set_initial_thread
  end subroutine solve_boundary_thread
  !=========================================================================
  !============================================================================!
  ! This routine solves three-diagonal system of equations:                    !
  !  ||m_1 u_1  0....        || ||w_1|| ||r_1||                                !
  !  ||l_2 m_2 u_2...        || ||w_2|| ||r_2||                                !
  !  || 0  l_3 m_3 u_3       ||.||w_3||=||r_3||                                !
  !  ||...                   || ||...|| ||...||                                !
  !  ||.............0 l_n m_n|| ||w_n|| ||r_n||                                !
  ! From: Numerical Recipes, Chapter 2.6, p.40.                                !
  !============================================================================!
  subroutine tridag(n,L_I,M_I,U_I,R_I,W_I)
    implicit none
    !--------------------------------------------------------------------------!
    integer, intent(in):: n
    real, intent(in):: L_I(n),M_I(n),U_I(n),R_I(n)
    real, intent(out):: W_I(n)
    !--------------------------------------------------------------------------!
    integer:: j
    real:: Aux,Aux_I(2:n)
    !--------------------------------------------------------------------------!
    if (M_I(1).eq.0.0) then
       call CON_stop('Error in tridag: M_I(1)=0')
    end if
    Aux = M_I(1)
    W_I(1) = R_I(1)/Aux
    do j=2,n
       Aux_I(j) = U_I(j-1)/Aux
       Aux = M_I(j)-L_I(j)*Aux_I(j)
       if (Aux.eq.0.0) then
          write(*,*)'j, M_I(j), L_I(j), Aux_I(j) = ',j, M_I(j),L_I(j),Aux_I(j)
          call CON_stop('Tridag failed')
       end if
       W_I(j) = (R_I(j)-L_I(j)*W_I(j-1))/Aux
    end do
    do j=n-1,1,-1
       W_I(j) = W_I(j)-Aux_I(j+1)*W_I(j+1)
    end do
    !------------------------------------ DONE --------------------------------!
  end subroutine tridag
  !=========================================================================
  subroutine set_field_line_thread_bc(nGhost, iBlock, nVarState, State_VG, &
               iImplBlock, IsLinear)
    use ModAdvance,      ONLY: State_VGB
    use BATL_lib, ONLY:  MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    use ModFaceGradient, ONLY: get_face_gradient
    use ModPhysics,      ONLY: No2Si_V, Si2No_V, UnitTemperature_, &
         UnitEnergyDens_, UnitU_, UnitX_, inv_gm1
    use ModMultiFluid,   ONLY: UseMultiIon, MassIon_I, ChargeIon_I, iRhoIon_I
    use ModVarIndexes,   ONLY: Rho_, Pe_, p_, Bx_, Bz_, &
         RhoUx_, RhoUz_, EHot_
    use ModGeometry,     ONLY: Xyz_DGB 
    use ModConst,        ONLY: cTolerance
    use ModImplicit,     ONLY: iTeImpl
    use ModWaves
    use ModHeatFluxCollisionless, ONLY: UseHeatFluxCollisionless, &
         get_gamma_collisionless
    use ModMain, ONLY: BlkTest, ProcTest, jTest, kTest
    use ModProcMH, ONLY: iProc
    integer, intent(in):: nGhost
    integer, intent(in):: iBlock
    integer, intent(in):: nVarState
    real, intent(inout):: State_VG(nVarState,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
 
    ! Optional arguments when called by semi-implicit scheme
    integer, optional, intent(in):: iImplBlock
    logical, optional, intent(in):: IsLinear

    !\
    ! Determines, which action should be done with the thread 
    ! before setting the BC
    !/
    integer:: iAction
    
    logical:: IsNewBlock
    integer :: i, j, k, Major_, Minor_
    real :: FaceGrad_D(3), TeSi, PeSi, BDir_D(3), U_D(3), B_D(3), SqrtRho
    real :: PeSiOut, U, AMinor, AMajor, DTeOverDsSi, DTeOverDs, TeGhost, Gamma
    real :: RhoNoDimOut, MinusDeltaROverBR
    logical:: DoTest, DoTestMe
    character(len=*), parameter :: NameSub = 'set_thread_bc'
    !--------------------------------------------------------------------------
    if(iBlock==BLKtest.and.iProc==PROCtest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    endif

    IsNewBlock = .true.
    if(present(iImplBlock))then
       if(BoundaryThreads_B(iBlock)%iAction/=Done_)&
            call CON_stop('Algorithm error in '//NameSub)
       iAction=Impl_
    else
       iAction=BoundaryThreads_B(iBlock)%iAction
    end if

    !\
    ! Start from floating boundary values
    !/
    do k = MinK, MaxK; do j = MinJ, maxJ; do i = 1 - nGhost, 0
       State_VG(:, i,j,k) = State_VG(:,1, j, k)
    end do; end do; end do
    !\
    ! Fill in the temperature array
    !/
    if(present(iImplBlock))then
       do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
          Te_G(i, j, k) = State_VG(iTeImpl,i,j,k)
       end do; end do; end do
    else
       if(UseMultiIon)then
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             Te_G(i,j,k) = State_VG(Pe_,i,j,k) &
                  /sum(ChargeIon_I*State_VG(iRhoIon_I,i,j,k)/MassIon_I)
          end do; end do; end do
       elseif(UseIdealEos)then
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             Te_G(i,j,k) = TeFraction*State_VG(iP,i,j,k) &
                  /State_VG(Rho_,i,j,k)
          end do; end do; end do
       else
          call CON_stop('Generic EOS is not applicable with threads')
       end if
    end if
    
    do k = 1, nK; do j = 1, nJ
       B_D = State_VGB(Bx_:Bz_,1,j,k,iBlock)
       BDir_D = B_D + &
            BoundaryThreads_B(iBlock) % B0Face_DII(:, j, k)
       BDir_D = BDir_D/max(sqrt(sum(BDir_D**2)), cTolerance)
       if(BoundaryThreads_B(iBlock) % SignBr_II(j, k) < 0.0)then
          BDir_D = -BDir_D
          Major_ = WaveLast_
          Minor_ = WaveFirst_
       else
          Major_ = WaveFirst_
          Minor_ = WaveLast_
       end if
       if(present(iImplBlock))then
          !\
          ! Gradient across the boundary face
          !/
          call get_face_gradient(1, 1, j, k, iBlock, &
               IsNewBlock, Te_G, FaceGrad_D, &
               UseFirstOrderBcIn=.true.)
          if(IsLinear)then
             if(BoundaryThreads_B(iBlock)%UseLimitedDTe_II(j,k))then
                State_VG(iTeImpl,0,j,k) = 0.0  
             else
                !\
                ! Apply linearization of DTe/Ds derivation in \delta Te(1,j,k)
                !/
                State_VG(iTeImpl,0,j,k) = Te_G(0,j,k) +(Te_G(0, j, k)*&
                     BoundaryThreads_B(iBlock)%DDTeOverDsOverTeTrueSi_II(j,k) &
                     /Si2No_V(UnitX_) - sum(FaceGrad_D*BDir_D))/&
                     sum(BoundaryThreads_B(iBlock)% DGradTeOverGhostTe_DII(&
                     :,j,k)* BDir_D)
             end if
             CYCLE
          end if
       end if
       !\
       ! Calculate input parameters for solving the thread
       !/
       if(present(iImplBlock))&
            BoundaryThreads_B(iBlock)%UseLimitedDTe_II(j,k)=&
            BoundaryThreads_B(iBlock) % TMax_II(j,k)<Te_G(0, j, k).or.&
            Te_G(0, j, k)<TeMin
       Te_G(0, j, k) = max(TeMin,min(Te_G(0, j, k), &
            BoundaryThreads_B(iBlock) % TMax_II(j,k)))
       TeSi = Te_G(0, j, k)*No2Si_V(UnitTemperature_)
       SqrtRho = sqrt(State_VGB(Rho_, 1, j, k, iBlock))
       AMinor = min(sqrt(State_VGB(Minor_, 1, j, k, iBlock)/&
            ( SqrtRho* PoyntingFluxPerB)  ),1.0)
       U_D = State_VGB(RhoUx_:RhoUz_, 1, j, k, iBlock)/&
            State_VGB(Rho_, 1, j, k, iBlock)
       U = sum(U_D*BDir_D)

       PeSi = PeFraction*State_VGB(iP, 1, j, k, iBlock)*No2Si_V(UnitEnergyDens_)
       call solve_boundary_thread(j=j, k=k, iBlock=iBlock, iAction=iAction,    &
            TeSiIn=TeSi, PeSiIn=PeSi, USiIn=U*No2Si_V(UnitU_), AMinorIn=AMinor,&
            DTeOverDsSiOut=DTeOverDsSi, PeSiOut=PeSiOut,&
            RhoNoDimOut=RhoNoDimOut, AMajorOut=AMajor)
       if(present(iImplBlock))then
          DTeOverDs = DTeOverDsSi * Si2No_V(UnitTemperature_)/Si2No_V(UnitX_)
          BoundaryThreads_B(iBlock)%UseLimitedDTe_II(j,k) = &
               BoundaryThreads_B(iBlock)%UseLimitedDTe_II(j,k).or.&
               DTeOverDs < 0.0
          !Do not allow heat flux from the TR to the low corona
          DTeOverDs = max(DTeOverDs,0.0)
          !\
          ! Calculate temperature in the ghost cell by adding the difference 
          ! between the required value DTeOverDs and the temperature gradient
          ! calculated with the floating BC 
          !/ 
          !\
          ! Trasformation coefficient -Delta R/b_R, avoid the division by zero.
          !/
          MinusDeltaROverBR = 1/&
            min(sum(BoundaryThreads_B(iBlock)% DGradTeOverGhostTe_DII(:, j, k) &
            * BDir_D),-0.1*sqrt(&
            sum(BoundaryThreads_B(iBlock)% DGradTeOverGhostTe_DII(:,j,k)**2))) 
          TeGhost = Te_G(0, j, k) +(&
               DTeOverDs - sum(FaceGrad_D*BDir_D) )*MinusDeltaROverBR
          Te_G(0, j, k) = TeGhost 
          State_VG(iTeImpl, 0, j, k) = Te_G(0, j, k)
          if(BoundaryThreads_B(iBlock)%UseLimitedDTe_II(j,k))CYCLE
          !\
          ! Calculate the derivative of DTeOverDs over \delta Te_G(1,j,k)
          !/
          BoundaryThreads_B(iBlock)%DDTeOverDsOverTeTrueSi_II(j,k) =&
               -DTeOverDsSi
          call solve_boundary_thread(j=j, k=k, iBlock=iBlock, &
               iAction=iAction, TeSiIn=1.01*TeSi, PeSiIn=1.01*PeSi, &
               USiIn=U*No2Si_V(UnitU_), AMinorIn=AMinor, &
               DTeOverDsSiOut=DTeOverDsSi, PeSiOut=PeSiOut,&
               RhoNoDimOut=RhoNoDimOut, AMajorOut=AMajor) 
          BoundaryThreads_B(iBlock)%DDTeOverDsOverTeTrueSi_II(j,k) = (&
               BoundaryThreads_B(iBlock)%DDTeOverDsOverTeTrueSi_II(j,k) + &
               DTeOverDsSi)/(0.01*TeSi)
          !\
          !To achieve the diagonal-dominance property
          !/
          BoundaryThreads_B(iBlock)%DDTeOverDsOverTeTrueSi_II(j,k) =max( &
               BoundaryThreads_B(iBlock)%DDTeOverDsOverTeTrueSi_II(j,k), 0.0)
          CYCLE
       end if

       !\
       ! 
       !/
       State_VG(iP, 0, j, k) = PeSiOut*Si2No_V(UnitEnergyDens_)/PeFraction
       !\
       !Exponential extrapolation of pressure
       !/
       State_VG(iP, 1-nGhost:-1, j, k) =  &
            State_VG(iP,0,j,k)**2/State_VG(iP,1,j,k)

       !\
       ! Assign ion pressure (if separate from electron one)
       if(iP/=p_)State_VG(p_, 1-nGhost:0, j, k) = &
            State_VG(iP, 1-nGhost:0, j, k)/AverageIonCharge

       State_VG(Rho_, 0, j, k) = State_VG(iP, 0, j, k)* &
             TeFraction/Te_G(1, j, k) !RhoNoDimOut
       !\
       !Exponential extrapolation of density
       !/
       State_VG(Rho_, 1-nGhost:-1, j, k) = State_VG(Rho_, 0, j, k)**2/ &
            State_VG(Rho_, 1, j, k) 

       !\
       ! Magnetic field and velocity vector components orthogonal to 
       ! B0 filed are all reflected 
       !/
       B_D = B_D - BDir_D*sum(B_D*BDir_D)
       U_D = U_D - BDir_D*U
       !\
       ! 
       !/
       U_D = -U_D; B_D = -B_D

       !\
       ! Velocity component along B0 field is reflected if inward directed
       ! Otherwise it is extrapolated inversely proportional to density
       !/
    
       do i = 1-nGhost, 0
          State_VG(Bx_:Bz_, i, j, k) = B_D
          State_VG(RhoUx_:RhoUz_, i, j, k) = State_VG(Rho_,  i, j, k) * &
               (max(-U,U*State_VG(Rho_,  1, j, k)/State_VG(Rho_,  i, j, k)) &
               *BDir_D + U_D)
          State_VG(Major_, i, j, k) = AMajor**2 * PoyntingFluxPerB *&
               sqrt( State_VG(Rho_, i, j, k) )
       end do

       if(Ehot_ > 1)then
          if(UseHeatFluxCollisionless)then
             call get_gamma_collisionless(Xyz_DGB(:,1,j,k,iBlock), Gamma)
             State_VG(Ehot_,1-nGhost:0,j,k) = &
                  State_VG(iP,1-nGhost:0,j,k)*(1.0/(Gamma - 1) - inv_gm1)
          else
             State_VG(Ehot_,1-nGhost:0,j,k) = 0.0
          end if
       end if
    end do; end do
    BoundaryThreads_B(iBlock)%iAction = Done_
  end subroutine set_field_line_thread_bc
end module ModThreadedLC
