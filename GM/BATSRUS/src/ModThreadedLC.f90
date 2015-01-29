!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModThreadedLC
  use ModFieldLineThread, ONLY: &
       BoundaryThreads, BoundaryThreads_B
  use ModCoronalHeating, ONLY:PoyntingFluxPerBSi, PoyntingFluxPerB, &
                              LPerpTimesSqrtBSi, LPerpTimesSqrtB
  use ModAdvance,    ONLY: UseElectronPressure, UseIdealEos
  use ModPhysics,    ONLY: AverageIonCharge
  use ModConst,         ONLY: rSun, mSun, cBoltzmann, cAtomicMass, cGravitation
  implicit none

  real :: TeFraction, TiFraction, PeFraction
  real,allocatable :: Te_G(:,:,:)
  real,allocatable,dimension(:):: APlus_I, AMinus_I, SqrtRho_I,&
       TeSi_I, ConsSi_I, PAvrSi_I, Res_I, U_I, L_I, M_I, Xi_I, &
       VaLog_I, ReflCoef_I

  real, parameter:: TeSiMin = 5.0e5    ![K]
  real           :: TeMin

  !\
  ! Gravitation potential
  !/
  real, parameter:: cGravPot = cGravitation*mSun*cAtomicMass/&
       (cBoltzmann*rSun)
  !\
  !   Hydrostatic equilibrium in an isothermal corona: 
  !    d(N_i*k_B*(Z*T_e +T_i) )/dr=G*M_sun*N_I*M_i*d(1/r)/dr
  ! => N_i*Te\propto exp(cGravPot/TeSi*(M_i[amu]/(1+Z))*\Delta(R_sun/r)) 
  !/
  !\
  !energy flux needed to raise the mass flux rho*u to the heliocentric 
  !distance r equals: rho*u*G*Msun*(1/R_sun -1/r)=
  !=k_B*N_i*M_i(amu)u*cGravPot*(1-R_sun/r)=
  !=P_e/T_e*cGravPot*(M_sun[amu]/Z)*(1/R_sun -1/r)
  integer:: iP
contains
  !=========================================================================
  subroutine init_threaded_lc
    use BATL_lib, ONLY:  MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    use ModMultiFluid,   ONLY: MassIon_I
    use ModFieldLineThread, ONLY: check_tr_table, get_poynting_flux, &
         nPointInThreadMax
    use ModPhysics,            ONLY: UnitTemperature_, Si2No_V
    use ModVarIndexes,         ONLY: Pe_, p_
    !-------------------
    allocate(Te_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)); Te_G = 0.0
    allocate(ReflCoef_I(nPointInThreadMax)); ReflCoef_I = 0.0
    allocate( SqrtRho_I(nPointInThreadMax));  SqrtRho_I = 0.0

    allocate( APlus_I(nPointInThreadMax));  APlus_I = 0.0
    allocate(AMinus_I(nPointInThreadMax)); AMinus_I = 0.0
    allocate(  TeSi_I(nPointInThreadMax));   TeSi_I = 0.0
    allocate(ConsSi_I(nPointInThreadMax)); ConsSi_I = 0.0
    allocate(PAvrSi_I(nPointInThreadMax)); PAvrSi_I = 0.0
    allocate(   Res_I(nPointInThreadMax));    Res_I = 0.0
    allocate(     U_I(nPointInThreadMax));      U_I = 0.0
    allocate(     L_I(nPointInThreadMax));      L_I = 0.0
    allocate(     M_I(nPointInThreadMax));      M_I = 0.0
    allocate(    Xi_I(nPointInThreadMax));     Xi_I = 0.0
    allocate( VaLog_I(nPointInThreadMax));  VaLog_I = 0.0
    ! TeFraction is used for ideal EOS:
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
    call get_poynting_flux(PoyntingFluxPerBSi)
    call check_tr_table
  end subroutine init_threaded_lc
  !================================
  !\
  ! Main routine:
  ! solves MHD equations along thread, to link the state above the
  ! inner boundary of the global solar corona to the photosphere
  !/
  subroutine solve_boundary_thread(j, k, iBlock, &
       TeSiIn, PeSiIn, USiIn, AMinorIn,          &
       DTeOverDsSiOut, PAvrSiOut, AMajorOut)
    !\
    ! USE:
    !/
    use ModFieldLineThread, ONLY: HeatCondParSi, solve_a_plus_minus
    use ModPhysics,      ONLY: inv_gm1, No2Si_V, UnitX_,Si2No_V, &
                               UnitEnergyDens_, UnitTemperature_, UnitB_
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
    integer,intent(in):: j, k, iBlock 
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
    !PAvrSiOut: The geometric mean of electron and ion pressure (\sqrt{Pe*Pi})
    !AMajorOut: For the wave propagating outward the Sun
    !            EnergyDensity = (\Pi/B)\sqrt{\rho} AMajor**2  
    !/
    real,  intent(out):: DTeOverDsSiOut, PAvrSiOut, AMajorOut

    !\
    ! Two components arrays to use lookup table
    !/ 
    real    :: Value_V(4), AWValue_V(2), Length, RhoNoDim, Heating, GravityCoef
    integer :: iTable, iTableAW
  
    !\
    ! 
    !/
    real:: Alpha, SqrtAlphaPlus1
    
    real           :: PSiMin
    integer:: nPoint, iPoint

    logical :: DoTest, DoTestMe

    character(len=*), parameter :: NameSub = 'solve_boundary_thread'
    !-------------------------------------------------------------------------
    if(iBlock==BLKtest.and.iProc==PROCtest.and.j==jTest.and.k==kTest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    endif

    iTable = i_lookup_table('TR')
    if(iTable<=0)call CON_stop('TR table is not set')

    call interpolate_lookup_table(iTable, TeSiIn, 1.0e8, Value_V, &
         DoExtrapolate=.false.)

    !\
    ! First value is now the product of the thread length in meters times
    ! a geometric mean pressure, so that
    !/
    PAvrSiOut = Value_V(1)/( BoundaryThreads_B(iBlock)% Length_III(0,j,k) * &
         No2Si_V(UnitX_))

    RhoNoDim = (PeSiIn*Si2No_V(UnitEnergyDens_)/PeFraction)*&
          TeFraction/(TeSiIn*Si2No_V(UnitTemperature_))

    !Dimmensionless length (related to the wave dissipation length)
    Length = BoundaryThreads_B(iBlock)% Length2SqrtB_III(0,j,k)*&
         sqrt(sqrt(RhoNoDim)*PoyntingFluxPerB/LperpTimesSqrtB**2)

    !Calculate Alfven waves for the given thread length and BC for ingoing wave 
    iTable = i_lookup_table('AW_TR')
    if(iTable<=0)call CON_stop('AW_TR table is not set')
    call interpolate_lookup_table(iTable, Length, AMinorIn, AWValue_V, &
         DoExtrapolate=.false.)

    Heating = AWValue_V(1)*PoyntingFluxPerBSi*&
         BoundaryThreads_B(iBlock)% B_III(0,j,k)*No2Si_V(UnitB_)
    !\
    !cGravPot = cGravitation*mSun*cAtomicMass/&
    !/   (cBoltzmann*rSun)
    !\
    !energy flux needed to raise the mass flux rho*u to the heliocentric 
    !distance r equals: rho*u*G*Msun*(1/R_sun -1/r)=
    !=k_B*N_i*M_i(amu)*u*cGravPot*(1-R_sun/r)=
    !=P_e/T_e*cGravPot*(M_sun[amu]/Z)*u*(1/R_sun -1/r)
    !/
    GravityCoef =  cGravPot/TeSiIn*MassIon_I(1)*          & 
         (1 - BoundaryThreads_B(iBlock)%RInv_III(0,j,k) )
    !\
    ! Heat flux equals PAvr * UHeat (spent for radiation) +
    ! Pi * U * (5/2) + Pe * U * (5/2) (carried away by outflow), 
    ! the pressure ratio being  Pe = Pi * AverageIonCharge
    ! 
    ! Temperature gradient equals the heat flux divided by kappa0*T**2.5
    !/
    DTeOverDsSiOut = ( PAvrSiOut * Value_V(2) - & !Radiation losses
         Heating                                & !AW Heating
         +USiIn * (PAvrSiOut/sqrt(AverageIonCharge)) *(inv_gm1 +1) * & !5/2*U*Pi
         (1 + AverageIonCharge) ) /&
         (HeatCondParSi * TeSiIn**2.50)


    PSiMin = PeSiIn * TeSiMin/TeSiIn
    Alpha = (-USiIn/Value_V(2)) * (PeSiIn/ PAvrSiOut) * &
         (inv_gm1 +1) * (1 + AverageIonCharge)/sqrt(AverageIonCharge)
    Alpha = max(Alpha, -1 +  (PSiMin/PAvrSiOut)**2)
   
    SqrtAlphaPlus1 = sqrt(1 + Alpha)
    ! PAvrSiOut =  PAvrSiOut*SqrtAlphaPlus1

    AMajorOut = AWValue_V(2)

    !Full 1D model
    !\
    ! Reduce by 0ne the number of points in order to avoid too
    ! close points to the photosphere
    !/ 
    nPoint = BoundaryThreads_B(iBlock)% nPoint_II(j,k) -1

    !\
    ! As a first approximation, recover Te from the analytical solution
    !/
    TeSi_I(nPoint) = TeSiIn
    iTable = i_lookup_table('TR')   
    do iPoint = 1, nPoint-1
       call interpolate_lookup_table(&
            iTable=iTable,       &
            iVal=1,              &
            ValIn=PAvrSiOut*     &
            BoundaryThreads_B(iBlock)% Length_III(iPoint-nPoint,j,k) * &
            No2Si_V(UnitX_),     &
            Arg2In=1.0e8,        &
            Value_V=Value_V,     &
            Arg1Out=TeSi_I(iPoint),  & 
            DoExtrapolate=.false.)
    end do
    !\
    ! The analytical solution assumes constant pressure and 
    ! no heating. Calculate the pressure distribution
    !/ 
    PAvrSi_I(1) = PAvrSiOut
    !\
    !   Hydrostatic equilibrium in an isothermal corona: 
    !    d(N_i*k_B*(Z*T_e +T_i) )/dr=G*M_sun*N_I*M_i*d(1/r)/dr
    ! => N_i*Te\propto exp(cGravPot/TeSi*(M_i[amu]/(1+Z))*\Delta(R_sun/r)) 
    !/
    GravityCoef = cGravPot*MassIon_I(1)/(AverageIonCharge + 1)
    do iPoint = 2, nPoint
       PAvrSi_I(iPoint) = PAvrSi_I(iPoint-1)*&
            exp( (0.5/TeSi_I(iPoint) + 0.5/TeSi_I(iPoint-1))*GravityCoef*&
            (BoundaryThreads_B(iBlock)%RInv_III(iPoint-nPoint,j,k) - &
            BoundaryThreads_B(iBlock)%RInv_III(iPoint-1-nPoint,j,k)))
    end do
    !\
    ! Now prepare to calculate Alfven wave amplitude distribution
    ! 1. Calculate sqrt(RhoNoDim
    !/
    do iPoint=1,nPoint
       RhoNoDim = (PAvrSi_I(iPoint)*sqrt(AverageIonCharge)*&
            Si2No_V(UnitEnergyDens_)/PeFraction)*&
            TeFraction/(TeSi_I(iPoint)*Si2No_V(UnitTemperature_))
       SqrtRho_I(iPoint) = sqrt(RhoNoDim)
    end do
    !\
    ! 2. Calculate dimensionless length (in terms of the dissipation length
    !/
    Xi_I(1) = 0.0
    do iPoint = 2, nPoint
       Xi_I(iPoint) = Xi_I(iPoint-1) + &
            (BoundaryThreads_B(iBlock)% Length2SqrtB_III(iPoint-nPoint,j,k) - &
            BoundaryThreads_B(iBlock)% Length2SqrtB_III(iPoint-1-nPoint,j,k))*&
         sqrt(0.5*(SqrtRho_I(iPoint) + SqrtRho_I(iPoint-1))&
         *PoyntingFluxPerB/LperpTimesSqrtB**2)
    end do
    !\
    ! 3. Calculate Alfven wave speed
    !/
    do iPoint = 1, nPoint
       VaLog_I(iPoint) = log(BoundaryThreads_B(iBlock)% B_III(iPoint-nPoint,j,k)&
            /SqrtRho_I(iPoint))
    end do
    !\
    ! 4. Calculate the reflection coefficient
    !/
    ReflCoef_I(1) = abs(VaLog_I(2) - VaLog_I(1))/(Xi_I(2) - Xi_I(1))
    do iPoint = 2, nPoint-1
       ReflCoef_i(iPoint) = min(abs(VaLog_I(iPoint) - VaLog_I(iPoint-1))/&
                                      (Xi_I(iPoint) -    Xi_I(iPoint-1)),&
                                abs(VaLog_I(iPoint+1) - VaLog_I(iPoint))/&
                                      (Xi_I(iPoint+1) -    Xi_I(iPoint)))
    end do
    ReflCoef_I(nPoint) = abs(VaLog_I(nPoint) - VaLog_I(nPoint-1))/&
                               (Xi_I(nPoint) -    Xi_I(nPoint-1))
    !\
    ! Solve amplitudes of the Alfven waves (arrays there have dimension
    ! (0:nI)
    !/
    call solve_a_plus_minus(&
            nI=nPoint-1,                    & 
            ReflCoef_I=ReflCoef_I(1:nPoint),&
            Xi_I=Xi_I(1:nPoint),            &
            AMinusBC=AMinorIn,              &
            Heating=AWValue_V(1),           &
            APlusBC=AWValue_V(2),           &
            APlusOut_I=APlus_I(1:nPoint),   &
            AMinusOut_I=AMinus_I(1:nPoint)  )

    if(DoTestMe)then
       write(*,*)'TeSiIn=       ',TeSiIn,' K '
       write(*,*)'TMax=         ', BoundaryThreads_B(iBlock) % TMax_II(j,k)&
            *No2Si_V(UnitTemperature_),' K'
       write(*,*)'PeSiIn = ', PeSiIn
       write(*,*)'NeSiIn = ', PeSiIn/(TeSiIn*cBoltzmann)
       write(*,*)'AMinorIn=     ', AMinorIn
       write(*,*)'USiIn=        ',USiIn,' m/s'
       write(*,*)'Thread Length=', &
            BoundaryThreads_B(iBlock)% Length_III(0,j,k) * &
            No2Si_V(UnitX_),' m, dimless length=',&
            BoundaryThreads_B(iBlock)% Length_III(0,j,k)
       write(*,*)'Dimensionless length characteristic of AW dissipation=',Length
       write(*,*)'Pressure=     ',PAvrSiOut
       write(*,*)'Alpha = ', Alpha
       write(*,*)'Poynting Flux max=',PoyntingFluxPerBSi*&
            BoundaryThreads_B(iBlock)% B_III(0,j,k)*No2Si_V(UnitB_),' W/m2'
       write(*,*)'Pointling flux in the TR=', Heating,' W/m2'
       write(*,*)'AMajorOut=    ', AMajorOut
       write(*,*)'Contributions to dT/ds:'
       write(*,*)'Radiation losses=',PAvrSiOut*Value_V(2)/SqrtAlphaPlus1 ,&
            ' W/m2'
       write(*,*)'Gravitational energy loss=',&
             USiIn * (PAvrSiOut/sqrt(AverageIonCharge))*GravityCoef
       write(*,*)'Final dT/ds=  ',DTeOverDsSiOut,&
            ', dT/ds*Length=',DTeOverDsSiOut*&
            BoundaryThreads_B(iBlock)% Length_III(0,j,k) * &
            No2Si_V(UnitX_),' K'
    end if

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
         UnitEnergyDens_, UnitU_, UnitX_, OmegaBody, inv_gm1
    use ModMultiFluid,   ONLY: UseMultiIon, MassIon_I, ChargeIon_I, iRhoIon_I
    use ModVarIndexes,   ONLY: nVar, Rho_, Pe_, p_, Bx_, Bz_, &
         RhoUx_, RhoUz_, EHot_
    use ModGeometry,     ONLY: Xyz_DGB 
    use ModSize,         ONLY: nDim
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

    logical:: IsNewBlock
    integer :: i, j, k, Major_, Minor_
    real :: FaceGrad_D(3), TeSi, PeSi, BDir_D(3), U_D(3), B_D(3), SqrtRho
    real :: PAvrSI, U, AMinor, AMajor, DTeOverDsSi, DTeOverDs, TeGhost, Gamma
    logical:: DoTest, DoTestMe
    character(len=*), parameter :: NameSub = 'set_thread_bc'
    !--------------------------------------------------------------------------
    if(iBlock==BLKtest.and.iProc==PROCtest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    endif

    IsNewBlock = .true.

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
       !\
       ! Gradient across the boundary face
       !/
       call get_face_gradient(1, 1, j, k, iBlock, &
            IsNewBlock, Te_G, FaceGrad_D, &
            UseFirstOrderBcIn=.true.)
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
          if(DoTestMe.and.j==jTest.and.k==kTest)&
               write(*,*)'iImplBlock=',iImplBlock,' IsLinear=',IsLinear
          if(IsLinear)then
             if(DoTestMe.and.j==jTest.and.k==kTest)&
                  write(*,*)'UseLimitedDTe=',&
                  BoundaryThreads_B(iBlock)%UseLimitedDTe_II(j,k)
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
       if(DoTestMe.and.j==jTest.and.k==kTest)&
            write(*,*)'Te_G(0,j,k)=',Te_G(0,j,k),' TMax=',&
            BoundaryThreads_B(iBlock) % TMax_II(j,k)
       if(present(iImplBlock))&
            BoundaryThreads_B(iBlock)%UseLimitedDTe_II(j,k)=&
            BoundaryThreads_B(iBlock) % TMax_II(j,k)<Te_G(0, j, k)
       Te_G(0, j, k) = min(Te_G(0, j, k), &
            BoundaryThreads_B(iBlock) % TMax_II(j,k))

       TeSi = Te_G(0, j, k) * No2Si_V(UnitTemperature_)

       SqrtRho = sqrt(State_VGB(Rho_, 1, j, k, iBlock))
       AMinor = sqrt(State_VGB(Minor_, 1, j, k, iBlock)/&
            ( SqrtRho* &
            PoyntingFluxPerB)  )
       U_D = State_VGB(RhoUx_:RhoUz_, 1, j, k, iBlock)/&
            State_VGB(Rho_, 1, j, k, iBlock)
       U = sum(U_D*BDir_D)

       PeSi = PeFraction*State_VGB(iP, 1, j, k, iBlock)*No2Si_V(UnitEnergyDens_)
       call solve_boundary_thread(j=j, k=k, iBlock=iBlock, &
            TeSiIn=TeSi, PeSiIn=PeSi, USiIn=U*No2Si_V(UnitU_), AMinorIn=AMinor,&
            DTeOverDsSiOut=DTeOverDsSi, PAvrSiOut=PAvrSi, AMajorOut=AMajor)
       DTeOverDs = DTeOverDsSi * Si2No_V(UnitTemperature_)/Si2No_V(UnitX_)
       if(present(iImplBlock))&
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

       TeGhost = Te_G(0, j, k) +(&
            DTeOverDs - sum(FaceGrad_D*BDir_D) )/&
            sum(BoundaryThreads_B(iBlock)% DGradTeOverGhostTe_DII(:, j, k) &
            * BDir_D)
       if(DoTestMe.and.j==jTest.and.k==kTest)then
          write(*,*)'TeSi=',TeSi,' K'
          write(*,*)'BDir_D=',BDir_D
          write(*,*)'DTeOverDs=', DTeOverDs
          write(*,*)'FaceGrad_D=',FaceGrad_D
          write(*,*)'sum(FaceGrad_D*BDir_D)=',sum(FaceGrad_D*BDir_D)
          write(*,*)'DGradTeOverGhostTe=',&
               BoundaryThreads_B(iBlock)% DGradTeOverGhostTe_DII(:, j, k)
          write(*,*)'DGradTeOverGhostTe\cdot BDir_D=',&
               sum(BoundaryThreads_B(iBlock)% DGradTeOverGhostTe_DII(:, j, k) &
               * BDir_D)
          write(*,*)'TeGhost=',TeGhost*No2Si_V(UnitTemperature_),' K' 
          Te_G(0, j, k) = TeGhost
       
          write(*,*)'Check if gradTe=dT/ds'
       
          call get_face_gradient(1, 1, j, k, iBlock, &
                IsNewBlock, Te_G, FaceGrad_D, &
                UseFirstOrderBcIn=.true.)
          write(*,*)'sum(BDir_D*FaceGrad_D)=',sum(BDir_D*FaceGrad_D)
       end if
       Te_G(0, j, k) = TeGhost 
       if(present(iImplBlock))then
          State_VG(iTeImpl, 0, j, k) = Te_G(0, j, k)
          if(BoundaryThreads_B(iBlock)%UseLimitedDTe_II(j,k))CYCLE
          !\
          ! Calculate the derivative of DTeOverDs over \delta Te_G(1,j,k)
          !/
          BoundaryThreads_B(iBlock)%DDTeOverDsOverTeTrueSi_II(j,k) =&
               -DTeOverDsSi
          call solve_boundary_thread(j=j, k=k, iBlock=iBlock, &
               TeSiIn=1.01*TeSi, PeSiIn=1.01*PeSi, &
               USiIn=U*No2Si_V(UnitU_), AMinorIn=AMinor, &
               DTeOverDsSiOut=DTeOverDsSi, PAvrSiOut=PAvrSi, AMajorOut=AMajor) 
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
       State_VG(iP, 0, j, k) = PAvrSi*Si2No_V(UnitEnergyDens_)*&
            sqrt(AverageIonCharge)/PeFraction
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
            TeFraction/Te_G(1, j, k)
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
  end subroutine set_field_line_thread_bc
end module ModThreadedLC
