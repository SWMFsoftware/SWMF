!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModThreadedLC
  use ModFieldLineThread, ONLY: &
       BoundaryThreads, BoundaryThreads_B, &
       LengthPAvrSi_, UHeat_, HeatFluxLength_, DHeatFluxXOverU_, &
       RadCool2Si, DoInit_, Done_, Enthalpy_, Heat_, iStage
  !\
  !   Hydrostatic equilibrium in an isothermal corona: 
  !    d(N_i*k_B*(Z*T_e +T_i) )/dr=G*M_sun*N_I*M_i*d(1/r)/dr
  ! => N_i*Te\propto exp(cGravPot/TeSi*(M_i[amu]/(1+Z))*\Delta(R_sun/r)) 
  !/
  !\
  ! The plasma properties dependent coefficient needed to evaluate the 
  ! eefect of gravity on the hydrostatic equilibrium
  !/
  use ModFieldLineThread, ONLY:  &
       GravHydroStat != cGravPot*MassIon_I(1)/(Z + 1)
  !\
  !energy flux needed to raise the mass flux rho*u to the heliocentric 
  !distance r equals: rho*u*G*Msun*(1/R_sun -1/r)=
  !=k_B*N_i*M_i(amu)u*cGravPot*(1-R_sun/r)=
  !=P_e/T_e*cGravPot*u(M_i[amu]/Z)*(1/R_sun -1/r)
  !/
  use ModFieldLineThread, ONLY:  &
       GravHydroDyn ! = cGravPot*MassIon_I(1)/Z 

  use ModCoronalHeating, ONLY:PoyntingFluxPerBSi, PoyntingFluxPerB, &
                              LPerpTimesSqrtBSi, LPerpTimesSqrtB

  use ModAdvance,    ONLY: UseElectronPressure, UseIdealEos
  use ModPhysics,    ONLY: Z => AverageIonCharge
  use ModConst,         ONLY: rSun, mSun, cBoltzmann, cAtomicMass, cGravitation
  use ModGeometry,   ONLY: Xyz_DGB
  implicit none
  !\
  ! To expsress Te in terms of P and rho.
  !/
  real    :: TeFraction, PeFraction
  integer :: iP
  !\
  ! Temperature 3D array
  !/
  real,allocatable :: Te_G(:,:,:)
  !\
  ! Arrays for 1D distributions
  !/
  real,allocatable,dimension(:)::ReflCoef_I, APlus_I, AMinus_I, &
       TeSi_I, TeSiOld_I, PeSi_I, Xi_I, Cons_I, &
       VaLog_I, DXi_I, ResHeating_I, ResCooling_I, ResEnthalpy_I,&
       ResHeatCond_I, ResGravity_I, SpecHeat_I, DeltaEnergy_I
  !\
  ! We apply ADI to solve state vector, the components of the state
  ! being temperature and log pressure. 
  ! The heating at constant pressure is characterized by the
  ! specific heat at constant pressure. For pressure the barometric
  ! formula is applied
  !/
  integer, parameter:: Cons_ = 1, Te_ = 1, LogP_ = 2
  real, allocatable, dimension(:,:):: Res_VI, DCons_VI
  real, allocatable, dimension(:,:,:) :: M_VVI, L_VVI, U_VVI
  
  !\
  ! Table numbers needed to use lookup table
  !/ 
  integer :: iTableTR, iTableRadcool  
  !\
  !Control parameters: minimum temerature, the temperature, below which 
  !the reflection turns to zero and two logicals with self-explained names
  !/
  real, parameter:: TeSiMin = 5.0e4, TeMinReflectionSi=4.0e4 
  logical        :: UseAlignedVelocity = .true., UseGravity = .false.
  !\
  ! 1. Coef to express radiation losses in terms of pressure (not density)
  ! 2. Coefficient to express dimensionless density as RhoNoDimCoef*PeSi/TeSi
  !/
  real           :: cCoolingPerPe2, RhoNoDimCoef 
  !\
  ! Misc
  !/
  real:: TeMin, ConsMin, ConsMax, TeSiMax
  !\
  ! For transforming conservative to TeSi and back 
  !/
  real, parameter:: cTwoSevenths = 2.0/7.0
  !\
  ! Gravitation potential in K
  !/
  real, parameter:: cGravPot = cGravitation*mSun*cAtomicMass/&
       (cBoltzmann*rSun)
  !\
  ! See above for the use of the latter constant
  !/

  real :: SqrtZ
  integer, parameter:: Impl_=3

  real,parameter:: CTol=1.0e-6
contains
  !=========================================================================
  subroutine init_threaded_lc
    use BATL_lib, ONLY:  MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    use ModLookupTable,  ONLY: i_lookup_table
    use ModMultiFluid,   ONLY: MassIon_I
    use ModFieldLineThread, ONLY: check_tr_table, & 
         nPointThreadMax, HeatCondParSi
    use ModPhysics,            ONLY: &
         UnitTemperature_, Si2No_V, UnitEnergyDens_
    use ModVarIndexes,         ONLY: Pe_, p_
    !-------------------
    allocate(Te_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)); Te_G = 0.0

    allocate(ReflCoef_I(0:nPointThreadMax)); ReflCoef_I = 0.0
    allocate(   APlus_I(0:nPointThreadMax));    APlus_I = 0.0
    allocate(  AMinus_I(0:nPointThreadMax));   AMinus_I = 0.0

    allocate(   TeSi_I(nPointThreadMax));     TeSi_I = 0.0
    allocate(TeSiOld_I(nPointThreadMax));  TeSiOld_I = 0.0
    allocate(   Cons_I(nPointThreadMax));     Cons_I = 0.0
    allocate( PeSi_I(nPointThreadMax));   PeSi_I = 0.0

    allocate(   SpecHeat_I(nPointThreadMax));   SpecHeat_I = 0.0
    allocate(DeltaEnergy_I(nPointThreadMax));DeltaEnergy_I = 0.0

    allocate( ResHeating_I(nPointThreadMax)); ResHeating_I = 0.0
    allocate( ResCooling_I(nPointThreadMax)); ResCooling_I = 0.0
    allocate(ResEnthalpy_I(nPointThreadMax));ResEnthalpy_I = 0.0
    allocate(ResHeatCond_I(nPointThreadMax));ResHeatCond_I = 0.0
    allocate( ResGravity_I(nPointThreadMax)); ResGravity_I = 0.0

    allocate(     Xi_I(0:nPointThreadMax));     Xi_I = 0.0
    allocate(  VaLog_I(nPointThreadMax));    VaLog_I = 0.0
    allocate(    DXi_I(nPointThreadMax));      DXi_I = 0.0

    allocate(  Res_VI(Te_:LogP_,nPointThreadMax));      Res_VI = 0.0
    allocate(DCons_VI(Te_:LogP_,nPointThreadMax));    DCons_VI = 0.0

    allocate(U_VVI(Te_:LogP_,Te_:LogP_,nPointThreadMax)); U_VVI = 0.0
    allocate(L_VVI(Te_:LogP_,Te_:LogP_,nPointThreadMax)); L_VVI = 0.0
    allocate(M_VVI(Te_:LogP_,Te_:LogP_,nPointThreadMax)); M_VVI = 0.0

    SqrtZ = sqrt(Z)

    !\
    ! TeFraction is used for ideal EOS:
    !/
    if(UseElectronPressure)then
       ! Pe = ne*Te (dimensionless) and n=rho/ionmass
       ! so that Pe = ne/n *n*Te = (ne/n)*(rho/ionmass)*Te
       ! TeFraction is defined such that Te = Pe/rho * TeFraction
       TeFraction = MassIon_I(1)/Z
       iP = Pe_
       PeFraction = 1.0
    else
       ! p = n*T + ne*Te (dimensionless) and n=rho/ionmass
       ! so that p=rho/massion *T*(1+ne/n Te/T)
       ! TeFraction is defined such that Te = p/rho * TeFraction
       TeFraction = MassIon_I(1) &
            /(1 + Z)
       iP = p_
       PeFraction = Z/(1.0 + Z) 
    end if
    !\
    ! Therefore Te = TeFraction*State_V(iP)/State_V(Rho_)
    ! Pe = PeFraction*State_V(iP)
    !/

    call check_tr_table
    iTableTR = i_lookup_table('TR')
    if(iTableTR<=0)call CON_stop('TR table is not set')
    iTableRadCool = i_lookup_table('radcool')
    if(iTableRadCool<=0)call CON_stop('Radiative cooling table is not set')

    TeMin = TeSiMin*Si2No_V(UnitTemperature_)

    ConsMin = cTwoSevenths*HeatCondParSi*TeSiMin**3.50

    !\
    !   Hydrostatic equilibrium in an isothermal corona: 
    !    d(N_i*k_B*(Z*T_e +T_i) )/dr=G*M_sun*N_I*M_i*d(1/r)/dr
    ! => N_i*Te\propto exp(cGravPot/TeSi*(M_i[amu]/(1+Z))*\Delta(R_sun/r)) 
    !/
    GravHydroStat = cGravPot*MassIon_I(1)/(Z + 1)
    !\
    !energy flux needed to raise the mass flux rho*u to the heliocentric 
    !distance r equals: rho*u*G*Msun*(1/R_sun -1/r)=
    !=k_B*N_i*M_i(amu)u*cGravPot*(1-R_sun/r)=
    !=P_e/T_e*cGravPot*u(M_i[amu]/Z)*(1/R_sun -1/r)
    !/
    GravHydroDyn  = cGravPot*MassIon_I(1)/Z
    !\
    ! With this constant, the volumetric radiative cooling rate is
    ! the table value multiplied by cCoolingPerPe2*(PeSi/TeSi)**2
    !/
    cCoolingPerPe2 = RadCool2Si/(cBoltzmann*cBoltzmann*Z)

    !\
    ! With this constant, the dimensionless density 
    ! equals RhoNoDimCoef*PeSi/TeSi
    !/
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
    use ModFieldLineThread, ONLY: HeatCondParSi
    use ModPhysics,      ONLY: InvGammaMinus1, Gamma, GammaMinus1,&
         No2Si_V, UnitX_,Si2No_V, UnitB_, UnitTemperature_
    use ModLookupTable,  ONLY: interpolate_lookup_table
    use ModMain,         ONLY: BlkTest, ProcTest, jTest, kTest
    use ModProcMH,       ONLY: iProc
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
    !            WaveEnergyDensity = (\Pi/B)\sqrt{\rho} AMinor**2 
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
    ! Arrays needed to use lookup table
    !/ 
    real    :: Value_V(4), ValCooling(1)
    !\
    !---------Used in 1D numerical model------------------------
    !/
    !\
    ! Number of TEMPERATURE nodes (first one is on the top of TR
    ! the last one is in the physical cell of the SC model
    !/
    integer :: nPoint 
    integer        :: nIter = 25
    !\
    ! Limited USiIn
    !/
    !\
    ! Corrrect density and pressure values in the ghost cell
    !/
    real :: GhostCellCorr

    logical :: DoTest, DoTestMe

    character(len=*), parameter :: NameSub = 'solve_boundary_thread'
    !-------------------------------------------------------------------------
    !\
    ! Initialize all output parameters from 0D solution
    !/
    call interpolate_lookup_table(iTableTR, TeSiIn, 1.0e8, Value_V, &
           DoExtrapolate=.false.)
    !\
    ! First value is now the product of the thread length in meters times
    ! a geometric mean pressure, so that
    !/
    PeSiOut        = Value_V(LengthPAvrSi_)*SqrtZ/&
         BoundaryThreads_B(iBlock)% LengthSi_III(0,j,k)
    DTeOverDsSiOut = PeSiOut*Value_V(UHeat_)/(HeatCondParSi*TeSiIn**2.5)
    RhoNoDimOut    = RhoNoDimCoef*PeSiOut/TeSiIn
    AMajorOut      = 1.0
    TeSiMax        = &
          BoundaryThreads_B(iBlock) % TMax_II(j,k)*No2Si_V(UnitTemperature_)
    ConsMax = cTwoSevenths*HeatCondParSi*TeSiMax**3.50
    if(iBlock==BLKtest.and.iProc==PROCtest.and.j==jTest.and.k==kTest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
       if(DoTestMe)then
          write(*,*)'TeSiIn=       ',TeSiIn,' K '
          write(*,*)'PeSiIn = ', PeSiIn
          write(*,*)'NeSiIn = ', PeSiIn/(TeSiIn*cBoltzmann)
          write(*,*)'Dimensionless density (input)=',RhoNoDimCoef*PeSiIn/TeSiIn
          write(*,*)'AMinorIn=     ', AMinorIn
          write(*,*)'USiIn=        ',USiIn,' m/s'
          write(*,*)'Thread Length=', &
               BoundaryThreads_B(iBlock)% LengthSi_III(0,j,k) &
               ,' m = ',  Si2No_V(UnitX_)*&
               BoundaryThreads_B(iBlock)% LengthSi_III(0,j,k),' Rs'
          write(*,*)'TeSiMax=       ',TeSiMax
          write(*,*)'iAction=',iAction
          write(*,*)'0D model results'
          write(*,*)'Pressure 0D=',PeSiOut
          write(*,*)'DTeOverDsSiOut=', DTeOverDsSiOut
          write(*,*)'RhoNoDimOut=', RhoNoDimOut 
       end if
    else
       DoTest=.false.; DoTestMe=.false.
    endif

    nPoint = BoundaryThreads_B(iBlock)% nPoint_II(j,k)
    if(iAction/=DoInit_)then
       !\
       ! Retrieve temperature and pressure distribution
       !/
       TeSi_I(1:nPoint) = BoundaryThreads_B(iBlock)%TSi_III(1-nPoint:0,j,k)
       PeSi_I(1:nPoint) = BoundaryThreads_B(iBlock)%PSi_III(1-nPoint:0,j,k)
       TeSi_I(nPoint)   = TeSiIn 
    end if
    select case(iAction)
    case(Enthalpy_)
       !\
       ! For a given ResHeating_I advance the effect of heating
       ! and hydrodynamic motion through a time step (or half time step
       ! for 2-stage scheme). Store temperature and pressure
       !/
       call advance_hydro_and_heating
    case(Impl_)
       call advance_heat_cond
       !\
       ! Output for temperature gradient, all the other outputs
       ! are meaningless
       !/
       DTeOverDsSiOut = max((TeSi_I(nPoint) - TeSi_I(nPoint-1))/&
            (BoundaryThreads_B(iBlock)% LengthSi_III(0,j,k) - &
            BoundaryThreads_B(iBlock)% LengthSi_III(-1,j,k)),0.0)
       if(DoTestMe)then
          write(*,*)'Final dT/ds=  ',DTeOverDsSiOut,&
               ', dT/ds*Length=',DTeOverDsSiOut*&
               BoundaryThreads_B(iBlock)% LengthSi_III(0,j,k),' K'
       end if
       RETURN
    case(Heat_)
       call advance_heat_cond
       !\
       ! Calculate AWaves and store pressure and temperature
       !/
       call get_res_heating(nIterIn=2*nIter) 
       BoundaryThreads_B(iBlock)%TSi_III(1-nPoint:0,j,k) = TeSi_I(1:nPoint) 
       BoundaryThreads_B(iBlock)%PSi_III(1-nPoint:0,j,k) = PeSi_I(1:nPoint)
    case(DoInit_)
       call set_initial_thread
    case default
       write(*,*)'iAction=',iAction
       call CON_stop('Unknown action in '//NameSub)
    end select
    !\
    ! Outputs
    !/
    GhostCellCorr =  BoundaryThreads_B(iBlock)% DeltaR_II(j,k)/&
         (1/BoundaryThreads_B(iBlock)%RInv_III(-1,j,k) - &
         1/BoundaryThreads_B(iBlock)%RInv_III(0,j,k) )
    PeSiOut = exp(log(PeSi_I(nPoint)) + &
         (log(PeSi_I(nPoint)) - log(PeSi_I(nPoint-1)))*GhostCellCorr )
    RhoNoDimOut = RhoNoDimCoef*PeSiOut/TeSi_I(nPoint)*max(&
         exp( -log(TeSi_I(nPoint  ))*GhostCellCorr  &
         +     log(TeSi_I(nPoint-1))*GhostCellCorr ),1.0)
    if(DoTestMe)then
       write(*,*)'AMajorOut=    ', AMajorOut
       write(*,*)'Before correction:'
       write(*,*)'Pressure 1D (SI) = ',PeSi_I(nPoint)
       write(*,*)'RhoNoDimOut      = ',&
            RhoNoDimCoef*PeSi_I(nPoint)/TeSi_I(nPoint)
       write(*,*)'Corrected:'
       write(*,*)'Pressure 1D (SI) = ',PeSiOut
       write(*,*)'RhoNoDimOut      = ',RhoNoDimOut
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
              (BoundaryThreads_B(iBlock)% Xi_III(iPoint-nPoint,j,k) - &
              BoundaryThreads_B(iBlock)% Xi_III(iPoint-1-nPoint,j,k))*&
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
      where(TeSi_I(1:nPoint)<TeMinReflectionSi)ReflCoef_I(1:nPoint)=0.0
      !\
      ! Solve amplitudes of the Alfven waves (arrays there have dimension)
      ! (0:nI)
      !/
      call solve_a_plus_minus(&
           nI=nPoint,                      & 
           ReflCoef_I=ReflCoef_I(0:nPoint),&
           Xi_I=Xi_I(0:nPoint),            &
           AMinusBC=AMinorIn,              &
           APlusBC=AMajorOut,              &
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
      PeSi_I(1) = Value_V(LengthPAvrSi_)*SqrtZ/&
           BoundaryThreads_B(iBlock)% LengthSi_III(1-nPoint,j,k)
      !\
      !   Hydrostatic equilibrium in an isothermal corona: 
      !    d(N_i*k_B*(Z*T_e +T_i) )/dr=G*M_sun*N_I*M_i*d(1/r)/dr
      ! => N_i*Te\propto exp(cGravPot/TeSi*(M_i[amu]/(1+Z))*\Delta(R_sun/r)) 
      !/
      !GravHydroStat = cGravPot*MassIon_I(1)/(Z + 1)
      do iPoint = 2, nPoint
         PeSi_I(iPoint) = PeSi_I(iPoint-1)*&
            exp( -BoundaryThreads_B(iBlock)%TGrav_III(iPoint-nPoint,j,k)&
            /TeSi_I(iPoint))
      end do
    end subroutine set_pressure
    !==========================
    subroutine get_heat_cond
      M_VVI = 0.0; ResHeatCond_I = 0.0; U_VVI = 0.0; L_VVI = 0.0
      !----------------
      !\
      ! Contribution from heat conduction fluxes
      !/
      !\
      ! Flux linearizations over small dCons
      ! Dimensionless flux= dCons/ds(1/( (PoyntingFlux/B)B)
      !/
      U_VVI(Cons_,Cons_,1:nPoint-1) = &
           -BoundaryThreads_B(iBlock)% BDsInvSi_III(1-nPoint:-1,j,k)
      L_VVI(Cons_,Cons_,2:nPoint-1) = U_VVI(Te_,Te_,1:nPoint-2)
      M_VVI(Cons_,Cons_,2:nPoint-1) = &
           -U_VVI(Cons_,Cons_,2:nPoint-1) - L_VVI(Cons_,Cons_,2:nPoint-1)
      !\
      ! Right heat fluxes
      !/
      ResHeatCond_I(1:nPoint-1) = &
           (Cons_I(1:nPoint-1) - Cons_I(2:nPoint))*U_VVI(Cons_,Cons_,1:nPoint-1)
      !\
      ! Add left heat flux to the TR
      !/
      ResHeatCond_I(1) = ResHeatCond_I(1) - Value_V(HeatFluxLength_)*&
           BoundaryThreads_B(iBlock)% BDsInvSi_III(-nPoint,j,k)
      !\
      ! Linearize left heat flux to the TR
      !/
      M_VVI(Cons_,Cons_,1) = -U_VVI(Cons_,Cons_,1) + Value_V(DHeatFluxXOverU_)*&
           BoundaryThreads_B(iBlock)% BDsInvSi_III(-nPoint,j,k)
      !\
      ! Add other left heat fluxes
      !/
      ResHeatCond_I(2:nPoint-1) = ResHeatCond_I(2:nPoint-1) + &
           (Cons_I(2:nPoint-1) - Cons_I(1:nPoint-2))*&
           L_VVI(Cons_,Cons_,2:nPoint-1)
      !\
      ! LogP_ terms
      !/
      M_VVI(LogP_,LogP_,1:nPoint-1) =  1.0
      !\
      ! We satisfy the equation, d(LogP) = d(Cons)*(dP/dCons)_{TR}/P
      !/
      M_VVI(LogP_,Cons_,1) = -1/&
              (Value_V(HeatFluxLength_)*Sqrt(Z))

      !\
      ! For other points we satisfy the hydrostatic equilibrium condition
      ! LogPe^{i-1}=LogPe^i+TGrav/Te^i
      ! dLogPe^i - dCons^i(TGrav/(Te^i)^2)*dTe/dCons
      !/
      M_VVI(LogP_,Cons_,2:nPoint-1) = -&
           BoundaryThreads_B(iBlock)% TGrav_III(2-nPoint:-1,j,k)/&
           (TeSi_I(2:nPoint-1)*3.50*Cons_I(2:nPoint-1))
      L_VVI(LogP_,LogP_,2:nPoint-1) = -1.0

      call get_cooling(nLast=nPoint-1)
      M_VVI(Te_,LogP_,2:nPoint-1) = &
              -2*ResCooling_I(2:nPoint-1) !=-dCooling/dLogPe
      !\
      ! For the first cell to improve the stability the term
      ! is added to diagonal, by expressing DeltaLogP in terms of DeltaCons
      !/
      M_VVI(Cons_,Cons_,1) =  M_VVI(Cons_,Cons_,1) + &
           2*ResCooling_I(1)*M_VVI(LogP_,Cons_,1)
    end subroutine get_heat_cond
    !===========================
    subroutine get_cooling(nLast)
      integer,intent(in) ::nLast
      integer ::  iPoint
      !-----------------
      ResCooling_I = 0.0;
      do iPoint = 1, nLast
         if(TeSi_I(iPoint)>1.0e8)then
            write(*,*)'Failure in heat condusction setting'
            write(*,*)'In the point Xyz=',Xyz_DGB(:,1,j,k,iBlock)
            write(*,*)'TeSiIn, PeSiIn = ', TeSiIn, PeSiIn
            write(*,*)'TeSi_I=',TeSi_I(1:nPoint)
            call CON_stop('Stop!!!')
         end if
         call interpolate_lookup_table(iTableRadCool,&
              TeSi_I(iPoint), 1.0e2, ValCooling)
         ResCooling_I(iPoint) = &
              -BoundaryThreads_B(iBlock)%DsOverBSi_III(iPoint-nPoint,j,k)&
              *ValCooling(1)*cCoolingPerPe2*&
              (PeSi_I(iPoint)/TeSi_I(iPoint))**2
      end do
    end subroutine get_cooling
    !======================
    subroutine set_initial_thread
      !\
      ! Enthalpy correction coefficient
      !/
      real    :: EnthalpyFlux, EnthalpyCorrection
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
              iTable=iTableTR,         &
              iVal=LengthPAvrSi_,      &
              ValIn=PeSiOut/SqrtZ*     &
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
      call get_res_heating(nIterIn=nIter)
      do iIter = 1,nIter
         TeSiOld_I(1:nPoint) = TeSi_I(1:nPoint)
         !\
         !Shape the source.
         !/
         call get_heat_cond
         !\
         ! Add enthalpy correction
         !/
         if(USiIn>0)then
            EnthalpyFlux = USiIn * (PeSi_I(nPoint)/Z)& !5/2*U*Pi
                 *(InvGammaMinus1 +1)*(1 + Z)/&
                 (TeSiIn*PoyntingFluxPerBSi*&
                 BoundaryThreads_B(iBlock)% B_III(0,j,k)*No2Si_V(UnitB_))
            EnthalpyCorrection = EnthalpyFlux*TeSi_I(1)
            ResEnthalpy_I(1) = EnthalpyCorrection
            do iPoint = 1, nPoint-2
               EnthalpyCorrection = EnthalpyFlux*TeSi_I(iPoint)
               ResEnthalpy_I(iPoint)   = &
                    ResEnthalpy_I(iPoint)   - EnthalpyCorrection
               ResEnthalpy_I(iPoint+1) = &
                    ResEnthalpy_I(iPoint+1) + EnthalpyCorrection
            end do
            EnthalpyCorrection = EnthalpyFlux*TeSi_I(nPoint-1)
            ResEnthalpy_I(nPoint-1)   = &
                 ResEnthalpy_I(nPoint-1)  - EnthalpyCorrection
         elseif(USiIn<0)then
            EnthalpyFlux = USiIn * (PeSiIn/Z)  & !5/2*U*Pi
                 *(InvGammaMinus1 +1)*(1 + Z)/&
                 (TeSiIn*PoyntingFluxPerBSi*&
                 BoundaryThreads_B(iBlock)% B_III(0,j,k)*No2Si_V(UnitB_))
            EnthalpyCorrection = EnthalpyFlux*TeSi_I(1)
            ResEnthalpy_I(1) = EnthalpyCorrection
            do iPoint = 2, nPoint-1
               EnthalpyCorrection = EnthalpyFlux*TeSi_I(iPoint)
               ResEnthalpy_I(iPoint)   = &
                    ResEnthalpy_I(iPoint)   + EnthalpyCorrection
               ResEnthalpy_I(iPoint-1) = &
                    ResEnthalpy_I(iPoint-1) - EnthalpyCorrection
            end do
            EnthalpyCorrection = EnthalpyFlux*TeSi_I(nPoint)
            ResEnthalpy_I(nPoint-1)   = &
                 ResEnthalpy_I(nPoint-1)  - EnthalpyCorrection
         end if
         ResEnthalpy_I=0.0
         Res_VI(Te_,1:nPoint-1) = &
              ResHeating_I(1:nPoint-1) +  ResCooling_I(1:nPoint-1) +&
              ResEnthalpy_I(1:nPoint-1) + ResHeatCond_I(1:nPoint-1)
         !==========Add Gravity Source================================
         !\
         !cGravPot = cGravitation*mSun*cAtomicMass/&
         !/   (cBoltzmann*rSun)
         !GravHydroDyn = cGravPot*MassIon_I(1)/Z
         !\
         !energy flux needed to raise the mass flux rho*u to the heliocentric 
         !distance r equals: rho*u*G*Msun*(1/R_sun -1/r)=
         !=k_B*N_i*M_i(amu)*u*cGravPot*(1-R_sun/r)=
         !=P_e/T_e*cGravPot*(M_ion[amu]/Z)*u*(1/R_sun -1/r)
         !/
         !Res_I(2:nPoint-1)=Res_I(2:nPoint-1) + 0.5*GravHydroDyn*EnthalpyFlux*(&
         !     - BoundaryThreads_B(iBlock)%RInv_III(1-nPoint:-2,j,k)&
         !     + BoundaryThreads_B(iBlock)%RInv_III(3-nPoint: 0,j,k))
         !Res_I(1) = Res_I(1) + GravHydroDyn*EnthalpyFlux*(-1 + 0.5*(&
         !       BoundaryThreads_B(iBlock)%RInv_III(1-nPoint,j,k)&
         !     + BoundaryThreads_B(iBlock)%RInv_III(2-nPoint,j,k)))
         DCons_VI = 0.0
         call tridiag_2by2_block(n=nPoint-1,  &
              L_VVI=L_VVI(:,:,1:nPoint-1),&
              M_VVI=M_VVI(:,:,1:nPoint-1),&
              U_VVI=U_VVI(:,:,1:nPoint-1),&
              R_VI=Res_VI(:,1:nPoint-1),  &
              W_VI=DCons_VI(:,1:nPoint-1))
         
         Cons_I(1:nPoint-1) = &
              max(ConsMin,Cons_I(1:nPoint-1) + DCons_VI(Te_,1:nPoint-1))
         TeSi_I(1:nPoint-1) = &
              (3.50*Cons_I(1:nPoint-1)/HeatCondParSi)**cTwoSevenths
         TeSi_I(2:nPoint-1) = max(TeSi_I(2:nPoint-1),&
              BoundaryThreads_B(iBlock)%TGrav_III(2-nPoint:-1,j,k))
         
         call interpolate_lookup_table(iTableTR, TeSi_I(1), 1.0e8, Value_V, &
              DoExtrapolate=.false.)
         call set_pressure
         if(DoTestMe)then
            write(*,*)&
                 'TeOld abs(DCons) Cons abs(DCons)/Cons Res_V TNew at iIter=',&
                 iIter
            do iPoint=1,nPoint
               write(*,*)iPoint, TeSiOld_I(iPoint), &
              abs(DCons_VI(Cons_,iPoint)),Cons_I(iPoint),         &
              abs(DCons_VI(Cons_,iPoint))/Cons_I(iPoint),        &
              Res_VI(:,iPoint), TeSi_I(iPoint)
            end do
         end if
         if(all(abs(DCons_VI(Te_,1:nPoint-1))<cTol*Cons_I(1:nPoint-1)))EXIT
      end do
      if(any(abs(DCons_VI(Te_,1:nPoint-1))>cTol*Cons_I(1:nPoint-1)))then
         write(*,*)'TeOld Te TeMin PeSi_I'
         do iPoint=1,nPoint
            write(*,*)iPoint, TeSiOld_I(iPoint),TeSi_I(iPoint),&
                 BoundaryThreads_B(iBlock)%TGrav_III(iPoint-nPoint,j,k),&
                 PeSi_I(iPoint)
         end do
         write(*,*)'Error =',maxval(&
              abs(DCons_VI(Te_,1:nPoint-1)/Cons_I(1:nPoint-1))),&
              ' at the point Xyz=',Xyz_DGB(:,1,j,k,iBlock)
         call CON_stop('Algorithm failure in set_initial_thread')
      end if
      call get_res_heating(nIterIn=nIter)
      BoundaryThreads_B(iBlock)%TSi_III(1-nPoint:0,j,k) = TeSi_I(1:nPoint) 
      BoundaryThreads_B(iBlock)%PSi_III(1-nPoint:0,j,k) = PeSi_I(1:nPoint)
    end subroutine set_initial_thread
    !==========================
    subroutine advance_hydro_and_heating
      use ModMain,    ONLY: cfl, nStage, time_accurate, Dt
      use ModAdvance, ONLY: time_BLK
      use ModPhysics, ONLY: No2Si_V, UnitT_
      !\
      ! Time step in the physical cell from which the thread originates
      ! divided by nStage if needed
      !/
      real ::DtLocal, UDt
      !\
      ! Enthalpy correction coefficient
      !/
      real    :: EnthalpyFlux, FluxConst
      !\
      ! Loop variable
      !/
      integer :: iPoint, iIter
      !-------------

      if(time_accurate)then
         DtLocal = Dt
      else
         DtLocal = cfl*time_BLK(1,j,k,iBlock)
      end if
      DtLocal = iStage*DtLocal*No2Si_V(UnitT_)/nStage

      if(DoTestMe)write(*,*)'DtLocal=', DtLocal

      SpecHeat_I(1:nPoint-1) = InvGammaMinus1 &
           * (1 + Z)/Z* &
           BoundaryThreads_B(iBlock)%DsOverBSi_III(1-nPoint:-1,j,k)

      call interpolate_lookup_table(iTableTR, TeSi_I(1), 1.0e8, Value_V, &
           DoExtrapolate=.false.)
      DeltaEnergy_I = 0.0
      call get_res_heating(nIterIn=nIter)
      ResHeating_I(1:nPoint-1) = ResHeating_I(1:nPoint-1)*DtLocal
      do iIter = 1,nIter
         TeSiOld_I(1:nPoint) = TeSi_I(1:nPoint)
         M_VVI = 0.0; L_VVI = 0.0; U_VVI = 0.0

         M_VVI(Te_,Te_,  1:nPoint-1) =  &
              PeSi_I(1:nPoint-1)*SpecHeat_I(1:nPoint-1)/TeSi_I(1:nPoint-1)
         !\
         ! LogP_ terms
         !/ 

         M_VVI(LogP_,LogP_,1:nPoint-1) =  1.0
         !\
         ! We satisfy the equation, d(LogP) = d(Cons)*(dP/dCons)_{TR}/P=
         ! HeatCondParSi*TeSi_I(1)**2.50*(dP/dCons)_{TR}/P
         !/
         M_VVI(LogP_,Cons_,1) = -HeatCondParSi*TeSi_I(1)**2.50/&
              (Value_V(HeatFluxLength_)*Sqrt(Z))

         !\
         ! For other points we satisfy the hydrostatic equilibrium condition
         ! LogPe^{i-1}=LogPe^i+TGrav/Te^i
         ! dLogPe^i - dTe^i(TGrav/(Te^i)^2) -dLogPe^{i-1}=0
         !/
         M_VVI(LogP_,Cons_,2:nPoint-1) = -&
              BoundaryThreads_B(iBlock)% TGrav_III(2-nPoint:-1,j,k)/&
              TeSi_I(2:nPoint-1)**2
         L_VVI(LogP_,LogP_,2:nPoint-1) = -1.0

         !\
         ! Multiply the heating source term and speed by the time step
         !/
         UDt = USiIn*DtLocal 
         !\
         ! Add enthalpy correction
         !/

         Res_VI = 0.0; FluxConst = 0.0; ResEnthalpy_I = 0.0

         if(UDt>0)then
            FluxConst    = UDt * PeSi_I(nPoint)/&
                 (TeSiIn*PoyntingFluxPerBSi*&
                 BoundaryThreads_B(iBlock)% B_III(0,j,k)*No2Si_V(UnitB_))

            EnthalpyFlux = FluxConst/Z& !5/2*U*Pi
                 *(InvGammaMinus1 +1)*(1 + Z)
            !\
            ! Solve equation!
            !SpecHeat*(T^{n+1}_i-T^n_i) + EnthalpyFlux*(T^{n+1}_i-T^{n+1}_{i-1})=0
            ! SpecHeat*(T^{n+1}_i-T^n_i) + EnthalpyFlux*(T^{n+1}_i-T^{n+1}_{i-1})-
            !                              -EnthalpyFlux*(T^n_i-T^n_{i-1}) = &
            !                  ResHeating  -EnthalpyFlux*(T^n_i-T^n_{i-1}) 
            !/

            L_VVI(Te_,Te_,2:nPoint-1) =  - EnthalpyFlux
            M_VVI(Te_,Te_,2:nPoint-1) =  EnthalpyFlux + M_VVI(Te_,Te_,2:nPoint-1)
            ResEnthalpy_I(2:nPoint-1) = &
                 EnthalpyFlux*(TeSi_I(1:nPoint-2) - TeSi_I(2:nPoint-1))
            ResEnthalpy_I(1)   = 0.0
         elseif(UDt<0)then
            FluxConst    = UDt * PeSiIn/&
                 (TeSiIn*PoyntingFluxPerBSi*&
                 BoundaryThreads_B(iBlock)% B_III(0,j,k)*No2Si_V(UnitB_)) 
            EnthalpyFlux = FluxConst/Z  & !5/2*U*Pi
                 *(InvGammaMinus1 +1)*(1 + Z)           
            !\
            ! Solve equation!
            ! SpecHeat*(T^{n+1}_i-T^n_i) + EnthalpyFlux*(T^{n+1}_{i+1}-T^{n+1}_i)=0
            ! SpecHeat*(T^{n+1}_i-T^n_i) + EnthalpyFlux*(T^{n+1}_{i+1}-T^{n+1}_i)&
            !                            - EnthalpyFlux*(T^n_{i+1}-T^n_i) = &
            !                        Res - EnthalpyFlux*(T^n_{i+1}-T^n_i)
            !/ 
            U_VVI(Te_,Te_,1:nPoint-1) = EnthalpyFlux
            M_VVI(Te_,Te_,2:nPoint-1) = M_VVI(Te_,Te_,2:nPoint-1) - EnthalpyFlux
            ResEnthalpy_I(1:nPoint-1) = &
                 -EnthalpyFlux*(TeSi_I(2:nPoint) - TeSi_I(1:nPoint-1))

            M_VVI(Te_,Te_,1) = M_VVI(Te_,Te_,1) - EnthalpyFlux
         end if
  
         !==========Add Gravity Source================================
         !\
         !cGravPot = cGravitation*mSun*cAtomicMass/&
         !/   (cBoltzmann*rSun)
         !GravHydroDyn = cGravPot*MassIon_I(1)/Z
         !\
         !energy flux needed to raise the mass flux rho*u to the heliocentric 
         !distance r equals: rho*u*G*Msun*(1/R_sun -1/r)=
         !=k_B*N_i*M_i(amu)*u*cGravPot*(1-R_sun/r)=
         !=P_e/T_e*cGravPot*(M_ion[amu]/Z)*u*(1/R_sun -1/r)
         !/
         !if(UseGravity)then
         !   ResGravity_I(2:nPoint-1) = 0.5*GravHydroDyn*FluxConst*(            &
!                 - BoundaryThreads_B(iBlock)%RInv_III(1-nPoint:-2,j,k)         &
 !                + BoundaryThreads_B(iBlock)%RInv_III(3-nPoint: 0,j,k))
 !           ResGravity_I(1)          = 0.5*GravHydroDyn*FluxConst*(            &
!                 - BoundaryThreads_B(iBlock)%RInv_III(1-nPoint,j,k)            &
!                 + BoundaryThreads_B(iBlock)%RInv_III(2-nPoint,j,k))
!            ResHeating_I(1:nPoint-1) = ResHeating_I(1:nPoint-1) + &
!                 ResGravity_I(1:nPoint-1) 
!         end if

         Res_VI(Te_,1:nPoint-1) = ResHeating_I(1:nPoint-1) + &
              ResEnthalpy_I(1:nPoint-1) - DeltaEnergy_I(1:nPoint-1)

         DCons_VI = 0.0
         call tridiag_2by2_block(n=nPoint-1,  &
              L_VVI=L_VVI(:,:,1:nPoint-1),&
              M_VVI=M_VVI(:,:,1:nPoint-1),&
              U_VVI=U_VVI(:,:,1:nPoint-1),&
              R_VI=Res_VI(:,1:nPoint-1),&
              W_VI=DCons_VI(:,1:nPoint-1))
         !\
         ! limit DeltaTe
         !/
         DCons_VI(Te_,1:nPoint-1) = max(min(DCons_VI(Te_,1:nPoint-1), &
                                    TeSiMax - TeSi_I(    1:nPoint-1)),&
                                    TeSiMin - TeSi_I(    1:nPoint-1))
         DCons_VI(Te_,2:nPoint-1) = max(DCons_VI(Te_,2:nPoint-1),     &
            BoundaryThreads_B(iBlock)% TGrav_III(2-nPoint:-1,j,k) -   &
                                              TeSi_I(    2:nPoint-1))   
         TeSi_I(1)          = max(TeSi_I(1) + DCons_VI(Te_,1), TeSiMin)
         TeSi_I(2:nPoint-1) = max(TeSi_I(2:nPoint-1) + DCons_VI(Te_,2:nPoint-1),&
              TeSiMin, BoundaryThreads_B(iBlock)% TGrav_III(2-nPoint:-1,j,k))
         DeltaEnergy_I(1:nPoint-1) = DeltaEnergy_I(1:nPoint-1) + &
              PeSi_I(1:nPoint-1)*SpecHeat_I(1:nPoint-1)/TeSi_I(1:nPoint-1)*&
              DCons_VI(Te_,1:nPoint-1)
         !\
         ! Set pressure for updated temperature 
         !/
         call interpolate_lookup_table(iTableTR, TeSi_I(1), 1.0e8, Value_V, &
              DoExtrapolate=.false.)
         call set_pressure
         if(all(abs(DCons_VI(Te_,1:nPoint-1))<cTol*TeSi_I(1:nPoint-1)))EXIT
      end do
      if(any(abs(DCons_VI(Te_,1:nPoint-1))>cTol*TeSi_I(1:nPoint-1)))then
         write(*,*)'TeOld Te TeMin PeSi_I'
         do iPoint=1,nPoint
            write(*,*)iPoint, TeSiOld_I(iPoint),TeSi_I(iPoint),&
                 BoundaryThreads_B(iBlock)%TGrav_III(iPoint-nPoint,j,k),&
                 PeSi_I(iPoint)
         end do
         write(*,*)'Error =',maxval(&
              abs(DCons_VI(Te_,1:nPoint-1)/TeSi_I(1:nPoint-1))),&
              ' at the point Xyz=',Xyz_DGB(:,1,j,k,iBlock)
         call CON_stop('Algorithm failure in advance_hydro')
      end if
      !\
      ! Store pressure and temperature
      !/
      if(iStage/=nStage)RETURN
      BoundaryThreads_B(iBlock)%TSi_III(1-nPoint:0,j,k) = TeSi_I(1:nPoint) 
      BoundaryThreads_B(iBlock)%PSi_III(1-nPoint:0,j,k) = PeSi_I(1:nPoint)
    end subroutine advance_hydro_and_heating
    !========================
    subroutine advance_heat_cond
      use ModMain,     ONLY: cfl, Dt, time_accurate
      use ModAdvance,  ONLY: time_BLK
      use ModPhysics,  ONLY: UnitT_, No2Si_V 
      !\
      ! Time step in the physical cell from which the thread originates
      !/
      real    :: DtLocal
      integer :: iIter, iPoint
      !-----------
      if(time_accurate)then
         DtLocal = Dt*No2Si_V(UnitT_)
      else
         DtLocal = cfl*time_BLK(1,j,k,iBlock)*No2Si_V(UnitT_)
      end if
      call interpolate_lookup_table(iTableTR, TeSi_I(1), 1.0e8, Value_V, &
           DoExtrapolate=.false.)
      Cons_I(1:nPoint) = cTwoSevenths*HeatCondParSi*TeSi_I(1:nPoint)**3.50
      SpecHeat_I(1:nPoint-1) = InvGammaMinus1  &
           * (1 + Z)/Z*  &
           BoundaryThreads_B(iBlock)%DsOverBSi_III(1-nPoint:-1,j,k)
      DeltaEnergy_I= 0.0; Res_VI=0.0
      do iIter = 1, nIter
         TeSiOld_I(1:nPoint) = TeSi_I(1:nPoint)
         call get_heat_cond
         !\
         ! Correct the temperature derivative near the TR
         ! due to the pressure dependence of radiative cooling rate, 
         ! the pressure near the TR beibg the (tabulated) function
         ! of the temperature near TR
         ! Contribution to dPe/dCons = -dCooling/dPe*dPe/dCons
         !                   = -2Cooling/Pe*Pe/(Pe*L*UHeat(T))
         !/
         Res_VI(Te_,1:nPoint-1) = - DeltaEnergy_I(1:nPoint-1) + &
              DtLocal*ResHeatCond_I(1:nPoint-1)      + &
              DtLocal*ResCooling_I(1:nPoint-1)
         U_VVI(Te_,Te_,1:nPoint-1) = DtLocal*U_VVI(Te_,Te_,1:nPoint-1)
         L_VVI(Te_,Te_,1:nPoint-1) = DtLocal*L_VVI(Te_,Te_,1:nPoint-1)
         M_VVI(Te_,Te_,1:nPoint-1) = DtLocal*M_VVI(Te_,Te_,1:nPoint-1) + &
              SpecHeat_I(1:nPoint-1)*PeSi_I(1:nPoint-1)/&!*Gamma/
              (3.50*Cons_I(1:nPoint-1))
         M_VVI(Te_,LogP_,1:nPoint-1) = DtLocal*M_VVI(Te_,LogP_,1:nPoint-1)! &
              !-GammaMinus1*SpecHeat_I(1:nPoint-1)*PeSi_I(1:nPoint-1)
         DCons_VI = 0.0
         call tridiag_2by2_block(n=nPoint-1,  &
              L_VVI=L_VVI(:,:,1:nPoint-1),&
              M_VVI=M_VVI(:,:,1:nPoint-1),&
              U_VVI=U_VVI(:,:,1:nPoint-1),&
              R_VI=Res_VI(:,1:nPoint-1),&
              W_VI=DCons_VI(:,1:nPoint-1))
         !\
         ! limit DeltaCons
         !/
         DCons_VI(Cons_,1:nPoint-1) = max(min(DCons_VI(Cons_,1:nPoint-1), &
                                      ConsMax - Cons_I(    1:nPoint-1)),  &
                                      ConsMin - Cons_I(    1:nPoint-1))
         DCons_VI(Cons_,2:nPoint-1) = max(DCons_VI(Cons_,2:nPoint-1),     &
                                      cTwoSevenths*HeatCondParSi*         &
          BoundaryThreads_B(iBlock)%TGrav_III(2-nPoint:-1,j,k)**3.50 -    &
                                      Cons_I(    2:nPoint-1)) 
         DeltaEnergy_I(1:nPoint-1) = DeltaEnergy_I(1:nPoint-1) + &
              SpecHeat_I(1:nPoint-1)*PeSi_I(1:nPoint-1)*         &
              DCons_VI(Cons_,1:nPoint-1)/(3.50*Cons_I(1:nPoint-1))
         Cons_I(1:nPoint-1) = &
              max(ConsMin,Cons_I(1:nPoint-1) + DCons_VI(Te_,1:nPoint-1))
         TeSi_I(1:nPoint-1) = &
              (3.50*Cons_I(1:nPoint-1)/HeatCondParSi)**cTwoSevenths
         TeSi_I(2:nPoint-1) = max(TeSi_I(2:nPoint-1),&
              BoundaryThreads_B(iBlock)%TGrav_III(2-nPoint:-1,j,k))
         !\
         ! Set pressure for updated temperature 
         !/
         call interpolate_lookup_table(iTableTR, TeSi_I(1), 1.0e8, Value_V, &
              DoExtrapolate=.false.)
         call set_pressure
         if(all(abs(DCons_VI(Te_,1:nPoint-1))<cTol*Cons_I(1:nPoint-1)))EXIT
      end do
      if(any(abs(DCons_VI(Te_,1:nPoint-1))>cTol*Cons_I(1:nPoint-1)))then
         write(*,*)'TeOld Te TeMin PeSi_I'
         do iPoint=1,nPoint
            write(*,*)iPoint, TeSiOld_I(iPoint),TeSi_I(iPoint),&
                 BoundaryThreads_B(iBlock)%TGrav_III(iPoint-nPoint,j,k),&
                 PeSi_I(iPoint)
         end do
         write(*,*)'Error =',maxval(&
              abs(DCons_VI(Te_,1:nPoint-1)/Cons_I(1:nPoint-1))),&
              ' at the point Xyz=',Xyz_DGB(:,1,j,k,iBlock)
         call CON_stop('Algorithm failure in advance_heat_cond')
      end if
    end subroutine advance_heat_cond
  end subroutine solve_boundary_thread
  !=========================================================================  
  ! This routine solves three-diagonal system of equations:                    !
  !  ||m_1 u_1  0....        || ||w_1|| ||r_1||                                !
  !  ||l_2 m_2 u_2...        || ||w_2|| ||r_2||                                !
  !  || 0  l_3 m_3 u_3       ||.||w_3||=||r_3||                                !
  !  ||...                   || ||...|| ||...||                                !
  !  ||.............0 l_n m_n|| ||w_n|| ||r_n||                                !
  ! Prototype: Numerical Recipes, Chapter 2.6, p.40.
  ! Here each of the compenets w_i and r_i are 2-component states and 
  ! m_i, l_i, u_i are 2*2 matrices                                             !
  !============================================================================!
  subroutine tridiag_2by2_block(n,L_VVI,M_VVI,U_VVI,R_VI,W_VI)
    implicit none
    !--------------------------------------------------------------------------!
    integer, intent(in):: n
    real, intent(in):: L_VVI(2,2,n),M_VVI(2,2,n),U_VVI(2,2,n),R_VI(2,n)
    real, intent(out):: W_VI(2,n)
    !--------------------------------------------------------------------------!
    integer:: j
    real   :: TildeM_VV(2,2), TildeMInv_VV(2,2), TildeMInvDotU_VVI(2,2,2:n)
    !--------------------------------------------------------------------------!
    !\
    ! If tilde(M)+L.Inverted(\tilde(M))\dot.U = M, then the equation 
    !      (M+L+U)W = R
    ! may be equivalently written as 
    ! (tilde(M) +L).(I + Inverted(\tilde(M)).U).W=R
    !/
    if (det(M_VVI(:,:,1)).eq.0.0) then
       call CON_stop('Error in tridiag: M_I(1)=0')
    end if
    TildeM_VV = M_VVI(:,:,1)
    TildeMInv_VV = inverted(TildeM_VV)
    !\
    ! First 2-vector element of the vector, Inverted(tilde(M) + L).R 
    !/
    W_VI(:,1) = matmul(TildeMInv_VV,R_VI(:,1))
    do j=2,n
       !\
       ! Next 2*2 blok element of the matrix, Inverted(Tilde(M)).U
       !/
       TildeMInvDotU_VVI(:,:,j) = matmul(TildeMInv_VV,U_VVI(:,:,j-1))
       !\
       ! Next 2*2 block element of matrix tilde(M), obeying the eq.
       ! tilde(M)+L.Inverted(\tilde(M))\dot.U = M
       !/
       TildeM_VV = M_VVI(:,:,j) - &
            matmul(L_VVI(:,:,j),TildeMInvDotU_VVI(:,:,j))
       if (det(TildeM_VV).eq.0.0) then
          write(*,*)'j, M_I(j), L_I(j), TildeMInvDotU_I(j) = ',j, &
               M_VVI(:,:,j),L_VVI(:,:,j),TildeMInvDotU_VVI(:,:,j)
          call CON_stop('2*2 block Tridiag failed')
       end if
       !\
       ! Next element of inverted(Tilde(M))
       !/
       TildeMInv_VV = inverted(TildeM_VV)
       !\
       ! Next 2-vector element of the vector, Inverted(tilde(M) + L).R
       ! satisfying the eq. (tilde(M) + L).W = R
       !/
       W_VI(:,j) = matmul(TildeMInv_VV,R_VI(:,j)) - &
            matmul(TildeMInv_VV,matmul(L_VVI(:,:,j),W_VI(:,j-1)))
    end do
    do j=n-1,1,-1
       !\
       ! Finally we solve equation
       ! (I + Inverted(Tilde(M)).U).W =  Inverted(tilde(M) + L).R
       !/
       W_VI(:,j) = W_VI(:,j)-matmul(TildeMInvDotU_VVI(:,:,j+1),W_VI(:,j+1))
    end do
  end subroutine tridiag_2by2_block
  !=============
  real function det(A_II)
    real, intent(in)::A_II(2,2)
    det = A_II(1,1)*A_II(2,2) - A_II(1,2)*A_II(2,1)
  end function det
  !============
  function inverted(A_II)RESULT(B_II)
    real,intent(in)::A_II(2,2)
    real           ::B_II(2,2)
    real           ::DetInv
    !------------------
    DetInv = 1/det(A_II)
    B_II(1,1) = A_II(2,2)*DetInv
    B_II(2,2) = A_II(1,1)*DetInv
    B_II(1,2) = -A_II(1,2)*DetInv
    B_II(2,1) = -A_II(2,1)*DetInv
  end function inverted
  !====================
  subroutine solve_a_plus_minus(nI, ReflCoef_I, Xi_I, AMinusBC,&
       APlusBC, nIterIn)
    !INPUT
    integer,         intent(in):: nI
    real,            intent(in):: ReflCoef_I(0:nI), Xi_I(0:nI)
    real,            intent(in):: AMinusBC         !BC for A-
    !OUTPUT
    real,           intent(out):: APlusBC          !BC for A+
    integer,optional,intent(in):: nIterIn
    real:: DeltaXi
    integer::iStep,iIter
    integer, parameter:: nIterMax = 10
    integer:: nIter
    real::Derivative, AOld, ADiffMax, AP, AM, APMid, AMMid
    !---------------------------------------------------------------------------
    APlus_I(0:nI)  = 1.0
    AMinus_I(0:nI) = AMinusBC
    if(present(nIterIn))then
       nIter=nIterIn
    else
       nIter=nIterMax
    end if
    do iIter=1,nIter
       !Go forward, integrate APlus_I with given AMinus_I
       ADiffMax = 0.0
       do iStep=1, nI
          !Predictor
          AP = APlus_I(iStep-1); AM = AMinus_I(iStep-1)
          Derivative = -AM*(max(0.0, AP - 2*AM) - max(0.0, AM - 2*AP) )*&
               min(0.5*ReflCoef_I(iStep-1)/max(AM,AP),  1.0)- &
               AP*AM
          AOld = APlus_I(iStep)
          DeltaXi = Xi_I(iStep) - Xi_I(iStep-1)

          !Corrector
          AMMid = 0.5*(AMinus_I(iStep-1) + AMinus_I(iStep))
          APMid = AP + 0.5*Derivative*DeltaXi
          Derivative = -AMMid*&
               (max(0.0, APMid -2*AMMid) -  max(0.0,AMMid - 2*APMid))*&
               min(0.250*(ReflCoef_I(iStep-1)+ReflCoef_I(iStep))/&
               max(AMMid,APMid), 1.0) - AMMid*APMid
          APlus_I(iStep) = AP + Derivative*DeltaXi
          ADiffMax = max(ADiffMax, &
               abs(AOld - APlus_I(iStep))/max(AOld,APlus_I(iStep)))
       end do
       !Go backward, integrate APlus_I with given AMinus_I
       !We integrate equation,
       !\
       ! 2da_-/d\xi=
       !=-[ max(1-2a_-/a_+,0)-max(1-2a_+/a_-,0)]* a_+ *
       ! *min(ReflCoef,2max(a_,a_+))-
       ! -2a_-a_+
       !/
       do iStep=nI - 1, 0, -1
          !Predictor
          AP = APlus_I(iStep+1); AM = AMinus_I(iStep+1)
          Derivative = -AP*(max(0.0,AP - 2*AM) - max(0.0,AM - 2*AP)     )*&
               min(0.5*ReflCoef_I(iStep+1)/max(AM,AP), 1.0) + &
               AP*AM
          AOld = AMinus_I(iStep)
          DeltaXi = Xi_I(iStep+1) - Xi_I(iStep) 

          !Corrector
          APMid = 0.5*(APlus_I(iStep+1) + APlus_I(iStep))
          AMMid = AM - 0.5*Derivative*DeltaXi
          Derivative = -APMid*&
               ( max(0.0,APMid -2*AMMid)- max(0.0,AMMid - 2*APMid) )*&
               min(0.250*(ReflCoef_I(iStep+1) + ReflCoef_I(iStep))/&
               max(AMMid, APMid), 1.0) + AMMid*APMid
          AMinus_I(iStep) = AMinus_I(iStep+1) - Derivative*DeltaXi
          ADiffMax = max(ADiffMax,&
               abs(AOld - AMinus_I(iStep))/max(AOld, AMinus_I(iStep)))
       end do
       if(ADiffMax<cTol)EXIT
    end do
    APlusBC = APlus_I(nI)
  end subroutine solve_a_plus_minus
  !=========================================================================
  subroutine set_field_line_thread_bc(nGhost, iBlock, nVarState, State_VG, &
               iImplBlock, IsLinear)
    use ModAdvance,      ONLY: State_VGB
    use BATL_lib, ONLY:  MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    use BATL_size,ONLY:  nJ, nK
    use ModPhysics,      ONLY: No2Si_V, Si2No_V, UnitTemperature_, &
         UnitEnergyDens_, UnitU_, UnitX_, InvGammaElectronMinus1
    use ModMultiFluid,   ONLY: UseMultiIon, MassIon_I, ChargeIon_I, iRhoIon_I
    use ModVarIndexes,   ONLY: Rho_, Pe_, p_, Bx_, Bz_, &
         RhoUx_, RhoUz_, EHot_
    use ModImplicit,     ONLY: iTeImpl
    use ModB0,           ONLY: B0_DGB
    use ModWaves,        ONLY: WaveFirst_, WaveLast_
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
    
    integer :: i, j, k, Major_, Minor_
    real :: TeSi, PeSi, BDir_D(3), U_D(3), U, B1_D(3), SqrtRho, DirR_D(3)
    real :: PeSiOut, AMinor, AMajor, DTeOverDsSi, DTeOverDs, GammaHere
    real :: RhoNoDimOut
    logical:: DoTest, DoTestMe
    character(len=*), parameter :: NameSub = 'set_thread_bc'
    !--------------------------------------------------------------------------
    if(present(iImplBlock))then
       if(BoundaryThreads_B(iBlock)%iAction/=Done_)&
            call CON_stop('Algorithm error in '//NameSub)
       iAction=Impl_
    else
       iAction=BoundaryThreads_B(iBlock)%iAction
    end if
    if(iAction==Done_)RETURN
    if(present(iImplBlock))then
       if(IsLinear)then
          !\
          !Version Easter 2015
          !/
          do k = MinK, MaxK; do j = MinJ, maxJ
             State_VG(iTeImpl,0,j,k) = 0.0
          end do; end do
          RETURN
       end if
    end if
    if(iBlock==BLKtest.and.iProc==PROCtest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    endif
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
       B1_D = State_VGB(Bx_:Bz_,1,j,k,iBlock)
       BDir_D = B1_D + 0.50*(B0_DGB(:, 1, j, k, iBlock) + &
            B0_DGB(:, 0, j, k, iBlock))
       BDir_D = BDir_D/max(sqrt(sum(BDir_D**2)), 1e-30)
       DirR_D = Xyz_DGB(:,1,j,k,iBlock)  
       DirR_D = DirR_D/max(sqrt(sum(DirR_D**2)),1e-30)
     
       if(sum(BDir_D*DirR_D) <  0.0)then
          BDir_D = -BDir_D
          Major_ = WaveLast_
          Minor_ = WaveFirst_
       else
          Major_ = WaveFirst_
          Minor_ = WaveLast_
       end if
       !\
       ! Calculate input parameters for solving the thread
       !/
       if(DoTestMe.and.j==jTest.and.k==kTest)then
          write(*,*)'Direction B=', BDir_D
          write(*,*)'Direction R=', DirR_D
          write(*,*)'Magnetic fields: B0True =', B0_DGB(:, 1, j, k, iBlock) 
          write(*,*)'Magnetic fields: B0Ghost=', B0_DGB(:, 0, j, k, iBlock) 
          write(*,*)'B1_D=', B1_D
          write(*,*)NameSub//': before limiting TeGhost=',&
               Te_G(0, j, k)*No2Si_V(UnitTemperature_)
       end if
       Te_G(0, j, k) = max(TeMin,min(Te_G(0, j, k), &
            BoundaryThreads_B(iBlock) % TMax_II(j,k)))
       TeSi = Te_G(0, j, k)*No2Si_V(UnitTemperature_)
       SqrtRho = sqrt(State_VGB(Rho_, 1, j, k, iBlock))
       AMinor = min(1.0,&
            sqrt(  State_VGB(Minor_, 1, j, k, iBlock)/&
            (SqrtRho* PoyntingFluxPerB)  ))
       U_D = State_VGB(RhoUx_:RhoUz_, 1, j, k, iBlock)/&
            State_VGB(Rho_, 1, j, k, iBlock)
       U = sum(U_D*BDir_D)
       if(DoTestMe.and.j==jTest.and.k==kTest)then
          write(*,*)'U_D=',U_D
          write(*,*)'Direction U', U_D/sqrt(sum(U_D**2))
       end if

       PeSi = PeFraction*State_VGB(iP, 1, j, k, iBlock)&
            *(Te_G(0,j,k)/Te_G(1,j,k))*No2Si_V(UnitEnergyDens_)
       call solve_boundary_thread(j=j, k=k, iBlock=iBlock, iAction=iAction,    &
            TeSiIn=TeSi, PeSiIn=PeSi, USiIn=U*No2Si_V(UnitU_), AMinorIn=AMinor,&
            DTeOverDsSiOut=DTeOverDsSi, PeSiOut=PeSiOut,&
            RhoNoDimOut=RhoNoDimOut, AMajorOut=AMajor)
       if(present(iImplBlock))then
          DTeOverDs = DTeOverDsSi * Si2No_V(UnitTemperature_)/Si2No_V(UnitX_)
          !\
          ! Solve equation: -(TeGhost-TeTrue)/DeltaR = 
          ! dTe/ds*(b . DirR)
          !/
          Te_G(0, j, k) = Te_G(0, j, k) - DTeOverDs*sum(BDir_D*DirR_D)*&
               BoundaryThreads_B(iBlock)% DeltaR_II(j,k)
          !\
          ! Version Easter 2015 Limit TeGhost
          !/
          Te_G(0, j, k) = max(TeMin,min(Te_G(0, j, k), &
               BoundaryThreads_B(iBlock) % TMax_II(j,k)))
          if(DoTestMe.and.j==jTest.and.k==kTest)then
             write(*,*)NameSub//': reset TeGhost=',&
                  Te_G(0, j, k)*No2Si_V(UnitTemperature_)
          end if
          State_VG(iTeImpl, 0, j, k) = Te_G(0, j, k)
          CYCLE
       end if
      
       State_VG(iP,0,j,k) = max(PeSiOut*Si2No_V(UnitEnergyDens_)/PeFraction,&
            0.90*State_VG(iP,1,j,k))
       !\
       !Extrapolation of pressure
       !/
       State_VG(iP, 1-nGhost:-1, j, k) = State_VG(iP,0,j,k)**2/&
            State_VG(iP,1,j,k)
       !\
       ! Assign ion pressure (if separate from electron one)
       !/
       if(iP/=p_)State_VG(p_, 1-nGhost:0, j, k) = &
            State_VG(iP, 1-nGhost:0, j, k)/Z

       State_VG(Rho_, 0, j, k) = max(RhoNoDimOut, 0.90*State_VG(Rho_, 1, j, k))
       !\
       !Extrapolation of density
       !/
       State_VG(Rho_, 1-nGhost:-1, j, k) = State_VG(Rho_, 0, j, k)**2&
            /State_VG(Rho_,1,j,k) 
   
       do i = 1-nGhost, 0
          !\
          ! Ghost cell value of the magnetic field: cancel radial B1 field 
          !/ 
          B1_D = State_VG(Bx_:Bz_, 1-i, j, k)
          State_VG(Bx_:Bz_, i, j, k) = B1_D - DirR_D*sum(DirR_D*B1_D)
          !\
          !Gnost cell value of velocity: keep the velocity projection
          !onto the magnetic field, if UseAlignedVelocity=.true.
          !Reflect the other components
          !/
          U_D = State_VG(RhoUx_:RhoUz_,1-i,j,k)/State_VG(Rho_,1-i,j,k)
          if(UseAlignedVelocity)then
             U   = sum(U_D*BDir_D); U_D = U_D - U*BDir_D
          else
             U = 0
          end if     
          State_VG(RhoUx_:RhoUz_, i, j, k) = -U_D*State_VG(Rho_,i,j,k) &
                + U*BDir_D*State_VG(Rho_,1,j,k)
              
          State_VG(Major_, i, j, k) = AMajor**2 * PoyntingFluxPerB *&
               sqrt( State_VG(Rho_, i, j, k) )
       end do

       if(Ehot_ > 1)then
          if(UseHeatFluxCollisionless)then
             call get_gamma_collisionless(Xyz_DGB(:,1,j,k,iBlock), GammaHere)
             State_VG(Ehot_,1-nGhost:0,j,k) = &
                  State_VG(iP,1-nGhost:0,j,k)*&
                  (1.0/(GammaHere - 1) - InvGammaElectronMinus1)
          else
             State_VG(Ehot_,1-nGhost:0,j,k) = 0.0
          end if
       end if
    end do; end do
    BoundaryThreads_B(iBlock)%iAction = Done_
  end subroutine set_field_line_thread_bc
end module ModThreadedLC
