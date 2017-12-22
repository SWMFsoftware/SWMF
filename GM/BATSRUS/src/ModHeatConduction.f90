!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModHeatConduction

  use BATL_lib, ONLY: &
       test_start, test_stop

  use ModHeatFluxCollisionless, ONLY: UseHeatFluxRegion, &
       rCollisional, rCollisionless
  use BATL_size, ONLY: nDim, MaxDim

  implicit none
  save

  private ! except

  ! Public methods
  public :: read_heatconduction_param
  public :: init_heat_conduction
  public :: get_heat_flux
  public :: get_ion_heat_flux
  public :: calc_ei_heat_exchange
  public :: get_impl_heat_cond_state
  public :: get_heat_conduction_rhs
  public :: add_jacobian_heat_cond
  public :: update_impl_heat_cond

  ! Logical for adding field-aligned heat conduction
  logical, public :: IsNewBlockHeatCond = .true.
  logical, public :: IsNewBlockIonHeatCond = .true.

  ! Variables for setting the field-aligned heat conduction coefficient
  character(len=20), public :: TypeHeatConduction = 'spitzer'
  logical :: DoUserHeatConduction
  character(len=20), public :: TypeIonHeatConduction = 'spitzer'
  logical :: DoUserIonHeatConduction

  ! These variables can be used in ModUser::user_material_properties
  ! to set the fraction of the field-algined heat conduction if
  ! DoWeakFieldConduction is set to true.
  real, public:: FractionFieldAligned  = -1.0
  real, public:: ElectronCollisionRate = -1.0

  ! Parameters for heat conduction in regions of weak magnetic field
  logical :: DoWeakFieldConduction = .false.

  ! Dimensionless heat conduction coefficients
  real :: HeatCondPar, IonHeatCondPar

  ! Unit conversion factor for heat conduction coefficients
  real :: Si2NoHeatCoef

  ! Coefficient for the electron-ion collision rate formula
  real :: ElectronIonCollisionCoef

  ! electron/ion temperature used for calculating heat flux
  real, allocatable :: Te_G(:,:,:), Ti_G(:,:,:)

  ! Used for ideal EOS: p = n*T + ne*Te (dimensionless) and n=rho/ionmass
  ! so that p=rho/massion *T*(1+ne/n Te/T)
  ! TiFraction is defined such that Ti = p/rho * TiFraction
  ! TeFraction is defined such that Te = p/rho * TeFraction
  real :: TiFraction, TeFraction

  ! Array needed for second order interpolation of ghost cells
  real, allocatable :: State1_VG(:,:,:,:), State2_VG(:,:,:,:)

  ! Heat flux for operator split scheme
  real, allocatable :: FluxImpl_VFD(:,:,:,:,:)

  ! Heat conduction dyad pre-multiplied by the face area
  real, allocatable :: HeatCond_DFDB(:,:,:,:,:,:)
  ! Arrays to build the Heat conduction dyad
  real, allocatable :: HeatCoef_G(:,:,:), bb_DDG(:,:,:,:,:)
  ! Arrays needed for the heat flux limiter
  real, allocatable :: FreeStreamFlux_G(:,:,:)

  ! electron-ion energy exchange
  real, allocatable :: PointCoef_CB(:,:,:,:)
  real, allocatable :: PointImpl_VCB(:,:,:,:,:)

  real:: cTeTiExchangeRate

  ! radiative cooling
  logical :: DoRadCooling = .false.

  ! Arrays for radiative cooling
  real, allocatable :: CoolHeat_CB(:,:,:,:)
  real, allocatable :: CoolHeatDeriv_CB(:,:,:,:)

contains
  !============================================================================

  subroutine read_heatconduction_param(NameCommand)

    use ModMain,      ONLY: UseHeatConduction, UseIonHeatConduction
    use ModReadParam, ONLY: read_var

    character(len=*), intent(in) :: NameCommand

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'read_heatconduction_param'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    select case(NameCommand)
    case("#HEATCONDUCTION")
       call read_var('UseHeatConduction', UseHeatConduction)
       if(UseHeatConduction)then
          call read_var('TypeHeatConduction', TypeHeatConduction)

          select case(TypeHeatConduction)
          case('user','spitzer')
          case default
             call stop_mpi(NameSub//': unknown TypeHeatConduction = ' &
                  //TypeHeatConduction)
          end select
       end if

    case("#WEAKFIELDCONDUCTION")
       call read_var('DoWeakFieldConduction', DoWeakFieldConduction)

    case("#IONHEATCONDUCTION")
       call read_var('UseIonHeatConduction', UseIonHeatConduction)
       if(UseIonHeatConduction)then
          call read_var('TypeIonHeatConduction', TypeIonHeatConduction)

          select case(TypeIonHeatConduction)
          case('user','spitzer')
          case default
             call stop_mpi(NameSub//': unknown TypeIonHeatConduction = ' &
                  //TypeIonHeatConduction)
          end select
       end if

    case default
       call stop_mpi(NameSub//' invalid NameCommand='//NameCommand)
    end select

    call test_stop(NameSub, DoTest)
  end subroutine read_heatconduction_param
  !============================================================================

  subroutine init_heat_conduction

    use BATL_size,     ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nI, nJ, nK, &
         j0_, nJp1_, k0_, nKp1_
    use ModAdvance,    ONLY: UseElectronPressure, UseAnisoPressure, UseAnisoPe
    use ModConst,      ONLY: cBoltzmann, cElectronMass, cProtonMass, &
         cEps, cElectronCharge
    use ModImplicit,   ONLY: UseSemiImplicit, nVarSemi, iTeImpl
    use ModMain,       ONLY: MaxBlock, UseHeatConduction, UseIonHeatConduction
    use ModMultiFluid, ONLY: UseMultiIon, MassIon_I
    use ModNumConst,   ONLY: cTwoPi
    use ModRadDiffusion, ONLY: UseHeatFluxLimiter
    use ModRadiativeCooling, ONLY: UseRadCooling
    use ModResistivity,  ONLY: UseHeatExchange
    use ModPhysics,    ONLY: Si2No_V, UnitEnergyDens_, UnitTemperature_, &
         UnitU_, UnitX_, UnitT_, No2Si_V, UnitN_, &
         ElectronTemperatureRatio, AverageIonCharge, CoulombLog
    use ModVarIndexes, ONLY: nVar

    real :: HeatCondParSi, IonHeatCondParSi
    real::  cTeTiExchangeRateSi

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'init_heat_conduction'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    if(.not.UseSemiImplicit)then
       if(UseHeatConduction)then
          if(.not.allocated(Te_G))allocate(Te_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK))
       else
          if(allocated(Te_G))deallocate(Te_G)
       end if
       if(UseIonHeatConduction)then
          if(.not.allocated(Ti_G))allocate(Ti_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK))
       else
          if(allocated(Ti_G))deallocate(Ti_G)
       end if
    end if

    ! TeFraction is used for ideal EOS:
    if(UseElectronPressure)then
       ! Pe = ne*Te (dimensionless) and n=rho/ionmass
       ! so that Pe = ne/n *n*Te = (ne/n)*(rho/ionmass)*Te
       ! TeFraction is defined such that Te = Pe/rho * TeFraction
       TiFraction = MassIon_I(1)
       TeFraction = MassIon_I(1)/AverageIonCharge
    else
       ! p = n*T + ne*Te (dimensionless) and n=rho/ionmass
       ! so that p=rho/massion *T*(1+ne/n Te/T)
       ! TeFraction is defined such that Te = p/rho * TeFraction
       TiFraction = MassIon_I(1) &
            /(1 + AverageIonCharge*ElectronTemperatureRatio)
       TeFraction = TiFraction*ElectronTemperatureRatio
    end if

    ! Conversion factor for heat conduction coefficient
    Si2NoHeatCoef = Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_) &
         *Si2No_V(UnitU_)*Si2No_V(UnitX_)

    ! Electron-ion collision rate coefficient for the formula
    ! ElectronIonCollision = ElectronIonCollisionCoef*Nion*Zion**2/Te^1.5
    ElectronIonCollisionCoef = &
         CoulombLog*(cElectronCharge**2/cEps)**2/(3*(cTwoPi*cBoltzmann)**1.5)

    ! To obtain the effective electron-ion collision frequency, this
    ! coefficient still need to be multiplied by Nion*Zion**2/Te**1.5.
    ! Here, we already take care of the units.
    ElectronIonCollisionCoef = ElectronIonCollisionCoef &
         *(1/Si2No_V(UnitT_))*No2Si_V(UnitN_)/No2Si_V(UnitTemperature_)**1.5

    ! electron heat conduct coefficient for single charged ions
    ! = 9.2e-12 W/(m*K^(7/2))
    HeatCondParSi = 3.2*3.0*cTwoPi/CoulombLog &
         *sqrt(cTwoPi*cBoltzmann/cElectronMass)*cBoltzmann &
         *((cEps/cElectronCharge)*(cBoltzmann/cElectronCharge))**2

    ! unit HeatCondParSi is W/(m*K^(7/2))
    HeatCondPar = HeatCondParSi &
         *Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_)**3.5 &
         *Si2No_V(UnitU_)*Si2No_V(UnitX_)

    ! ion heat conduct coefficient
    ! = 2.6e-13 W/(m*K^(7/2)) for protons
    IonHeatCondParSi = 3.9*3.0*cTwoPi/CoulombLog &
         *sqrt(cTwoPi*cBoltzmann/(MassIon_I(1)*cProtonMass))*cBoltzmann &
         *((cEps/cElectronCharge)*(cBoltzmann/cElectronCharge))**2

    ! unit HeatCondParSi is W/(m*K^(7/2))
    IonHeatCondPar = IonHeatCondParSi &
         *Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_)**3.5 &
         *Si2No_V(UnitU_)*Si2No_V(UnitX_)

    !\
    ! In hydrogen palsma, the electron-ion heat exchange is described by
    ! the equation as follows:
    ! dTe/dt = -(Te-Ti)/(tau_{ei})
    ! dTi/dt = +(Te-Ti)/(tau_{ei})
    ! The expression for 1/tau_{ei} may be found in
    ! Lifshitz&Pitaevskii, Physical Kinetics, Eq.42.5
    ! note that in the Russian edition they denote k_B T as Te and
    ! the factor 3 is missed in the denominator:
    ! 1/tau_ei = 2* CoulombLog * sqrt{m_e} (e^2/cEps)**2* Z**2 *Ni /&
    ! ( 3 * (2\pi k_B Te)**1.5 M_p). This exchange rate scales linearly
    ! with the plasma density, therefore, we introduce its ratio to
    ! the particle concentration. We calculate the temperature exchange
    ! rate by multiplying the expression for electron-ion effective
    ! collision rate,
    ! \nu_{ei} = CoulombLog/sqrt(cElectronMass)*  &
    !            ( cElectronCharge**2 / cEps)**2 /&
    !            ( 3 *(cTwoPi*cBoltzmann)**1.50 )* Ne/Te**1.5
    !  and then multiply in by the energy exchange coefficient
    !            (2*cElectronMass/cProtonMass)
    ! The calculation of the effective electron-ion collision rate is
    ! re-usable and can be also applied to calculate the resistivity:
    ! \eta = m \nu_{ei}/(e**2 Ne)
    !/

    cTeTiExchangeRateSi = &
         CoulombLog/sqrt(cElectronMass)*  &!\
         ( cElectronCharge**2 / cEps)**2 /&! effective ei collision frequency
         ( 3 *(cTwoPi*cBoltzmann)**1.50 ) &!/
         *(2*cElectronMass/cProtonMass)    ! *energy exchange per ei collision
    !\
    ! While used, this should be divided by TeSi**1.5 and multipled by
    ! atomic density, N_i in SI. We will apply dimensionless density
    ! so that the transformation coefficient shoule be multiplied by
    ! No2Si_V(UnitN_). We will employ dimensionless Te, therefore, we should
    ! divide by No2Si_V(UnitTemperature_)**1.5. We also need to convert the
    ! exchange rate from inverse seconds to dimensionless units by dividing by
    ! Si2No_V(UnitT_)
    !/

    cTeTiExchangeRate = cTeTiExchangeRateSi * &
         (1/Si2No_V(UnitT_))*No2Si_V(UnitN_)/No2Si_V(UnitTemperature_)**1.5

    DoUserHeatConduction    = TypeHeatConduction == 'user'
    DoUserIonHeatConduction = TypeIonHeatConduction == 'user'

    if(UseSemiImplicit.and..not.allocated(HeatCond_DFDB))then
       allocate( &
            State1_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK), &
            State2_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK), &
            FluxImpl_VFD(nVarSemi,nI+1,nJ+1,nK+1,nDim), &
            HeatCoef_G(0:nI+1,j0_:nJp1_,k0_:nKp1_), &
            bb_DDG(MaxDim,MaxDim,0:nI+1,j0_:nJp1_,k0_:nKp1_), &
            HeatCond_DFDB(nDim,nI+1,nJ+1,nK+1,nDim,MaxBlock), &
            Te_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK) )

       if(UseHeatFluxLimiter) &
            allocate(FreeStreamFlux_G(0:nI+1,j0_:nJp1_,k0_:nKp1_))

       if(UseElectronPressure .and. .not.UseMultiIon)then
          allocate(PointCoef_CB(nI,nJ,nK,MaxBlock))
          if(UseAnisoPressure)then
             allocate(PointImpl_VCB(2,nI,nJ,nK,MaxBlock))
          else
             allocate(PointImpl_VCB(1,nI,nJ,nK,MaxBlock))
          end if

          if(UseAnisoPe) call stop_mpi(NameSub// &
               ' heat conduction for UseAnisoPe has not been implemented yet.')

          UseHeatExchange = .false.
       end if

       if(UseRadCooling)then
          DoRadCooling = UseRadCooling
          UseRadCooling = .false.
       end if

       if(DoRadCooling) &
            allocate( &
            CoolHeat_CB(nI,nJ,nK,MaxBlock), &
            CoolHeatDeriv_CB(nI,nJ,nK,MaxBlock))

       iTeImpl = 1
    end if

    call test_stop(NameSub, DoTest)
  end subroutine init_heat_conduction
  !============================================================================

  subroutine get_heat_flux(iDir, iFace, jFace, kFace, iBlock, &
       StateLeft_V, StateRight_V, Normal_D, HeatCondCoefNormal, HeatFlux)

    use BATL_lib,        ONLY: Xyz_DGB
    use BATL_size,       ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    use ModAdvance,      ONLY: State_VGB, UseIdealEos, UseElectronPressure
    use ModFaceGradient, ONLY: get_face_gradient
    use ModPhysics,      ONLY: Si2No_V, UnitTemperature_, &
         UnitEnergyDens_, InvGammaElectronMinus1
    use ModVarIndexes,   ONLY: nVar, Rho_, p_, Pe_, Ehot_
    use ModMultifluid,   ONLY: UseMultiIon, MassIon_I, ChargeIon_I, iRhoIon_I
    use ModUserInterface ! user_material_properties
    use ModMain,         ONLY: UseFieldLineThreads, nDim, nIJK_D
    use ModGeometry,     ONLY: far_field_BCs_BLK
    use ModParallel,     ONLY: NOBLK, NeiLev
    use ModHeatFluxCollisionless, ONLY: UseHeatFluxCollisionless, &
         get_gamma_collisionless

    integer, intent(in) :: iDir, iFace, jFace, kFace, iBlock
    real,    intent(in) :: StateLeft_V(nVar), StateRight_V(nVar), Normal_D(3)
    real,    intent(out):: HeatCondCoefNormal, HeatFlux

    integer :: i, j, k, iP, iFace_D(3)
    real :: HeatCondL_D(3), HeatCondR_D(3), HeatCond_D(3), HeatCondFactor
    real :: FaceGrad_D(3), TeSi, CvL, CvR, CvSi, NumDensL, NumDensR, GammaTmp
    logical :: UseFirstOrderBc = .false.
    logical :: UseLeftStateOnly = .false., UseRightStateOnly = .false.

    character(len=*), parameter:: NameSub = 'get_heat_flux'
    !--------------------------------------------------------------------------
    ! Use first order flux across the computational domain boundary with
    ! threaded-field-line-model

    if(UseFieldLineThreads)then
       UseFirstOrderBc = far_field_BCs_BLK(iBlock)
    else
       UseFirstOrderBc = .false.
    end if
    if(UseFirstOrderBc)then
       iFace_D = (/iFace, jFace, kFace/)
       UseRightStateOnly = any(&
            iFace_D(1:nDim)==1.and.NeiLev(1:(2*nDim-1):2,iBlock)==NOBLK)
       UseLeftStateOnly =  any(&
            iFace_D(1:nDim)==nIJK_D(1:nDim)+1&
            .and.NeiLev(2:2*nDim:2,iBlock)==NOBLK)
    else
       UseRightStateOnly = .false.
       UseLeftStateOnly  = .false.
    end if

    if(IsNewBlockHeatCond)then
       if(UseMultiIon)then
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             Te_G(i,j,k) = State_VGB(Pe_,i,j,k,iBlock) &
                  /sum(ChargeIon_I*State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I)
          end do; end do; end do
       elseif(UseIdealEos)then
          iP = p_
          if(UseElectronPressure)  iP = Pe_

          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             Te_G(i,j,k) = TeFraction*State_VGB(iP,i,j,k,iBlock) &
                  /State_VGB(Rho_,i,j,k,iBlock)
          end do; end do; end do
       else
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             call user_material_properties( &
                  State_VGB(:,i,j,k,iBlock), i, j, k, iBlock, TeOut=TeSi)
             Te_G(i,j,k) = TeSi*Si2No_V(UnitTemperature_)
          end do; end do; end do
       end if
    end if

    call get_face_gradient(iDir, iFace, jFace, kFace, iBlock, &
         IsNewBlockHeatCond, Te_G, FaceGrad_D, &
         UseFirstOrderBcIn=UseFirstOrderBc)
    if(UseLeftStateOnly)then
       call get_heat_cond_coef(iDir, iFace, jFace, kFace, iBlock, &
            StateLeft_V, Normal_D, HeatCond_D)
    elseif(UseRightStateOnly)then
       call get_heat_cond_coef(iDir, iFace, jFace, kFace, iBlock, &
            StateRight_V, Normal_D, HeatCond_D)
    else

       call get_heat_cond_coef(iDir, iFace, jFace, kFace, iBlock, &
            StateLeft_V, Normal_D, HeatCondL_D)
       call get_heat_cond_coef(iDir, iFace, jFace, kFace, iBlock, &
            StateRight_V, Normal_D, HeatCondR_D)

       HeatCond_D = 0.5*(HeatCondL_D + HeatCondR_D)

    end if

    if(UseHeatFluxRegion)then
       HeatCondFactor = heat_cond_factor(iDir, iFace, jFace, kFace, iBlock)
       HeatCond_D = HeatCond_D*HeatCondFactor
    end if

    HeatFlux = -sum(HeatCond_D*FaceGrad_D)

    ! get the heat conduction coefficient normal to the face for
    ! time step restriction
    if(UseIdealEos)then
       if(UseMultiIon)then
          NumDensL = sum(ChargeIon_I*StateLeft_V(iRhoIon_I)/MassIon_I)
          NumDensR = sum(ChargeIon_I*StateRight_V(iRhoIon_I)/MassIon_I)
       else
          NumDensL = StateLeft_V(Rho_)/TeFraction
          NumDensR = StateRight_V(Rho_)/TeFraction
       end if
       if(Ehot_ > 1 .and. UseHeatFluxCollisionless)then
          call get_gamma_collisionless(Xyz_DGB(:,i,j,k,iBlock), GammaTmp)
          CvL = NumDensL/(GammaTmp - 1)
          CvR = NumDensR/(GammaTmp - 1)
       else
          ! Heat flux is carried by electrons, for single fluid w/o the
          ! electron equation, InvGammaElectronMinus1 = InvGammaMinus1
          CvL = InvGammaElectronMinus1*NumDensL
          CvR = InvGammaElectronMinus1*NumDensR
       end if
    else
       call user_material_properties(StateLeft_V, &
            iFace, jFace, kFace, iBlock, iDir, CvOut = CvSi)
       CvL = CvSi*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_)
       call user_material_properties(StateRight_V, &
            iFace, jFace, kFace, iBlock, iDir, CvOut = CvSi)
       CvR = CvSi*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_)
    end if

    if(UseLeftStateOnly)then
       HeatCondCoefNormal = sum(HeatCond_D*Normal_D)/CvL
    elseif(UseRightStateOnly)then
       HeatCondCoefNormal = sum(HeatCond_D*Normal_D)/CvR
    else
       HeatCondCoefNormal = sum(HeatCond_D*Normal_D)/min(CvL,CvR)
    end if

  end subroutine get_heat_flux
  !============================================================================

  subroutine get_heat_cond_coef(iDir, iFace, jFace, kFace, iBlock, &
       State_V, Normal_D, HeatCond_D)

    use ModAdvance,      ONLY: UseIdealEos, UseElectronPressure
    use ModB0,           ONLY: B0_DX, B0_DY, B0_DZ
    use ModMain,         ONLY: UseB0
    use ModNumConst,     ONLY: cTolerance
    use ModPhysics,      ONLY: No2Si_V, Si2No_V, UnitTemperature_, &
         ElectronGyroFreqCoef
    use ModVarIndexes,   ONLY: nVar, Bx_, Bz_, Rho_, p_, Pe_
    use ModRadiativeCooling, ONLY: DoExtendTransitionRegion, extension_factor
    use ModMultifluid,   ONLY: UseMultiIon, MassIon_I, ChargeIon_I, iRhoIon_I
    use ModUserInterface ! user_material_properties

    integer, intent(in) :: iDir, iFace, jFace, kFace, iBlock
    real, intent(in) :: State_V(nVar), Normal_D(3)
    real, intent(out):: HeatCond_D(3)

    real :: B_D(3), Bnorm, Bunit_D(3), TeSi, Te
    real :: HeatCoefSi, HeatCoef

    character(len=*), parameter:: NameSub = 'get_heat_cond_coef'
    !--------------------------------------------------------------------------
    if(UseB0)then
       select case(iDir)
       case(1)
          B_D = State_V(Bx_:Bz_) + B0_DX(:,iFace,jFace,kFace)
       case(2)
          B_D = State_V(Bx_:Bz_) + B0_DY(:,iFace,jFace,kFace)
       case(3)
          B_D = State_V(Bx_:Bz_) + B0_DZ(:,iFace,jFace,kFace)
       end select
    else
       B_D = State_V(Bx_:Bz_)
    end if

    ! The magnetic field should nowhere be zero. The following fix will
    ! push the magnitude of Bunit_D  to zero if B_D is approaching zero.
    Bnorm = sqrt(sum(B_D**2))
    Bunit_D = B_D/max(Bnorm,cTolerance)

    ! Calculate the electron temperature in SI units
    if(UseIdealEos)then
       if(UseMultiIon)then
          Te = State_V(Pe_)/sum(ChargeIon_I*State_V(iRhoIon_I)/MassIon_I)
       else
          if(UseElectronPressure)then
             Te = TeFraction*State_V(Pe_)/State_V(Rho_)
          else
             Te = TeFraction*State_V(p_)/State_V(Rho_)
          end if
       end if
       !\
       ! To calculate the extension factor for the transition
       ! region we may need the temperature in Kelvin.
       !/
       TeSi = Te*No2Si_V(UnitTemperature_)
    else
       call user_material_properties(State_V, &
            iFace, jFace, kFace, iBlock, iDir, TeOut=TeSi)
       Te = TeSi*Si2No_V(UnitTemperature_)
    end if

    if(DoWeakFieldConduction)then
       ! Initialize these public variables. The user can change them in
       ! user_material_properties when HeatCondOut is present
       ElectronCollisionRate = 0.0
       FractionFieldAligned  = -1.0
    end if

    if(DoUserHeatConduction .or. .not.UseIdealEos)then
       call user_material_properties(State_V, &
            iFace, jFace, kFace, iBlock, iDir, &
            TeIn=TeSi, HeatCondOut=HeatCoefSi)
       if(HeatCoefSi < 0.0)then
          ! Spitzer conductivity if user sets negative HeatCoefSi
          HeatCoef = HeatCondPar*Te**2.5
       else
          HeatCoef = HeatCoefSi * Si2NoHeatCoef
       end if
    else
       ! Spitzer conductivity for collisional regime
       HeatCoef = HeatCondPar*Te**2.5
    end if

    ! Artificial modified heat conduction for a smoother transition
    ! region, Linker et al. (2001)
    if(DoExtendTransitionRegion) HeatCoef = HeatCoef*extension_factor(TeSi)

    if(DoWeakFieldConduction)then
       ! If the user did not set the field-aligned fraction, calculate it as
       ! FractionFieldAligned = 1/(1 + ElectronCollisionRate/OmegaElectron)
       ! OmegaElectron = B*q_e/m_e = B*ElectronGyroFreqCoef
       ! The collision rate is the sum of the user supplied value
       ! (can be the neutral-electron rate) and
       ! the ion-electron collision rate.
       if(FractionFieldAligned < 0.0)then
          ElectronCollisionRate = ElectronCollisionRate + &
               ElectronIonCollisionCoef* &
               sum(State_V(iRhoIon_I)/MassIon_I*ChargeIon_I**2)/Te**1.5
          FractionFieldAligned = &
               1/(1 + ElectronCollisionRate/(Bnorm*ElectronGyroFreqCoef))
       end if

       HeatCond_D = HeatCoef*( &
            FractionFieldAligned*sum(Bunit_D*Normal_D)*Bunit_D &
            + (1 - FractionFieldAligned)*Normal_D )
    else
       HeatCond_D = HeatCoef*sum(Bunit_D*Normal_D)*Bunit_D
    end if

  end subroutine get_heat_cond_coef
  !============================================================================

  subroutine get_ion_heat_flux(iDir, iFace, jFace, kFace, iBlock, &
       StateLeft_V, StateRight_V, Normal_D, HeatCondCoefNormal, HeatFlux)

    use BATL_size,       ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    use ModAdvance,      ONLY: State_VGB, UseIdealEos
    use ModFaceGradient, ONLY: get_face_gradient
    use ModPhysics,      ONLY: InvGammaMinus1
    use ModVarIndexes,   ONLY: nVar, Rho_, p_

    integer, intent(in) :: iDir, iFace, jFace, kFace, iBlock
    real,    intent(in) :: StateLeft_V(nVar), StateRight_V(nVar), Normal_D(3)
    real,    intent(out):: HeatCondCoefNormal, HeatFlux

    integer :: i, j, k
    real :: HeatCondL_D(3), HeatCondR_D(3), HeatCond_D(3), HeatCondFactor
    real :: FaceGrad_D(3), CvL, CvR

    character(len=*), parameter:: NameSub = 'get_ion_heat_flux'
    !--------------------------------------------------------------------------
    if(IsNewBlockIonHeatCond)then
       if(UseIdealEos .and. .not.DoUserIonHeatConduction)then
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             Ti_G(i,j,k) = TiFraction*State_VGB(p_,i,j,k,iBlock) &
                  /State_VGB(Rho_,i,j,k,iBlock)
          end do; end do; end do
       else
          call stop_mpi(NameSub// &
               ': no ion heat conduction yet for non-ideal eos')
       end if
    end if

    call get_face_gradient(iDir, iFace, jFace, kFace, iBlock, &
         IsNewBlockIonHeatCond, Ti_G, FaceGrad_D)

    call get_ion_heat_cond_coef(iDir, iFace, jFace, kFace, iBlock, &
         StateLeft_V, Normal_D, HeatCondL_D)
    call get_ion_heat_cond_coef(iDir, iFace, jFace, kFace, iBlock, &
         StateRight_V, Normal_D, HeatCondR_D)

    HeatCond_D = 0.5*(HeatCondL_D + HeatCondR_D)

    if(UseHeatFluxRegion)then
       HeatCondFactor = heat_cond_factor(iDir, iFace, jFace, kFace, iBlock)
       HeatCond_D = HeatCond_D*HeatCondFactor
    end if

    HeatFlux = -sum(HeatCond_D*FaceGrad_D)

    ! get the heat conduction coefficient normal to the face for
    ! time step restriction
    if(UseIdealEos .and. .not.DoUserIonHeatConduction)then
       CvL = InvGammaMinus1*StateLeft_V(Rho_)/TiFraction
       CvR = InvGammaMinus1*StateRight_V(Rho_)/TiFraction
    else
       call stop_mpi(NameSub// &
            ': no ion heat conduction yet for non-ideal eos')
    end if
    HeatCondCoefNormal = sum(HeatCond_D*Normal_D)/min(CvL,CvR)

  end subroutine get_ion_heat_flux
  !============================================================================

  subroutine get_ion_heat_cond_coef(iDir, iFace, jFace, kFace, iBlock, &
       State_V, Normal_D, HeatCond_D)

    use ModAdvance,    ONLY: UseIdealEos
    use ModB0,         ONLY: B0_DX, B0_DY, B0_DZ
    use ModMain,       ONLY: UseB0
    use ModNumConst,   ONLY: cTolerance
    use ModVarIndexes, ONLY: nVar, Bx_, Bz_, Rho_, p_
    use ModUserInterface ! user_material_properties

    integer, intent(in) :: iDir, iFace, jFace, kFace, iBlock
    real, intent(in) :: State_V(nVar), Normal_D(3)
    real, intent(out):: HeatCond_D(3)

    real :: B_D(3), Bnorm, Bunit_D(3), Ti
    real :: IonHeatCoefSi, IonHeatCoef

    character(len=*), parameter:: NameSub = 'get_ion_heat_cond_coef'
    !--------------------------------------------------------------------------
    if(UseB0)then
       select case(iDir)
       case(1)
          B_D = State_V(Bx_:Bz_) + B0_DX(:,iFace,jFace,kFace)
       case(2)
          B_D = State_V(Bx_:Bz_) + B0_DY(:,iFace,jFace,kFace)
       case(3)
          B_D = State_V(Bx_:Bz_) + B0_DZ(:,iFace,jFace,kFace)
       end select
    else
       B_D = State_V(Bx_:Bz_)
    end if

    ! The magnetic field should nowhere be zero. The following fix will
    ! turn the magnitude of the field direction to zero.
    Bnorm = sqrt(sum(B_D**2))
    Bunit_D = B_D/max(Bnorm,cTolerance)

    if(UseIdealEos .and. .not.DoUserIonHeatConduction)then
       Ti = TiFraction*State_V(p_)/State_V(Rho_)

       ! Spitzer form for collisional regime
       IonHeatCoef = IonHeatCondPar*Ti**2.5
    else
       call stop_mpi(NameSub//': no ion heat conduction yet for non-ideal eos')

       call user_material_properties(State_V, &
            iFace, jFace, kFace, iBlock, iDir, IonHeatCondOut=IonHeatCoefSi)
       IonHeatCoef = IonHeatCoefSi*Si2NoHeatCoef
    end if

    HeatCond_D = IonHeatCoef*sum(Bunit_D*Normal_D)*Bunit_D

  end subroutine get_ion_heat_cond_coef
  !============================================================================

  real function heat_cond_factor(iDir, iFace, jFace, kFace, iBlock)

    use BATL_lib,    ONLY: Xyz_DGB

    integer, intent(in) :: iDir, iFace, jFace, kFace, iBlock

    real :: x_D(MaxDim), r
    !--------------------------------------------------------------------------
    select case(iDir)
    case(1)
       x_D = 0.5*(Xyz_DGB(:,iFace-1,jFace,kFace,iBlock) &
            +     Xyz_DGB(:,iFace  ,jFace,kFace,iBlock))
    case(2)
       x_D = 0.5*(Xyz_DGB(:,iFace,jFace-1,kFace,iBlock) &
            +     Xyz_DGB(:,iFace,jFace  ,kFace,iBlock))
    case(3)
       x_D = 0.5*(Xyz_DGB(:,iFace,jFace,kFace-1,iBlock) &
            +     Xyz_DGB(:,iFace,jFace,kFace  ,iBlock))
    end select
    r = sqrt(sum(x_D**2))
    if(rCollisionless < 0.0)then
       heat_cond_factor = 1.0/((r/rCollisional)**2 + 1)
    elseif(r <= rCollisional)then
       heat_cond_factor = 1.0
    else
       heat_cond_factor = &
            exp(-((r-rCollisional)/(rCollisionless-rCollisional))**2)
    end if

  end function heat_cond_factor
  !============================================================================
  subroutine calc_ei_heat_exchange

    ! Non-split operator, (almost) explicit ei heat energy exchange

    use ModMain,       ONLY: Cfl, nBlock, Unused_B, nI, nJ, nK
    use ModGeometry,   ONLY: true_cell
    use ModPhysics,    ONLY: Si2No_V, UnitTemperature_
    use ModVarIndexes, ONLY: Rho_, p_, Pe_, Ppar_
    use ModAdvance,    ONLY: time_blk, State_VGB, UseAnisoPressure, &
                             UseIdealEos
    use ModEnergy,     ONLY: calc_energy_cell
    use ModMultifluid, ONLY: ChargeIon_I,MassIon_I, iRhoIon_I, UseMultiIon
    use ModUserInterface ! user_material_properties

    real :: DtLocal, TeSi
    real :: HeatExchange, HeatExchangePeP, HeatExchangePePpar
    integer:: i, j, k, iBlock
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'calc_ei_heat_exchange'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    HeatExchange = 0.0
    HeatExchangePeP = 0.0
    HeatExchangePePpar = 0.0

    do iBlock = 1, nBlock
       if (Unused_B(iBlock)) CYCLE
      ! For the electron flux limiter, we need Te in the ghostcells
       if(UseMultiIon)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             Te_G(i,j,k) = State_VGB(Pe_,i,j,k,iBlock)/sum( &
                  ChargeIon_I*State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I)
          end do; end do; end do
       elseif(UseIdealEos)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             Te_G(i,j,k) = TeFraction &
                  *State_VGB(Pe_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
          end do; end do; end do
       else
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  i, j, k, iBlock, TeOut = TeSi)
             Te_G(i,j,k) = TeSi*Si2No_V(UnitTemperature_)
          end do; end do; end do
       end if
       !\
       ! More work is to be done for if .not.UseIdealEos or UseMultion
       !/
       if(.not.UseIdealEos)call stop_mpi(&
            'No explicit ei heat exchange for non-idealized plasmas')

       if(UseMultiion)call stop_mpi(&
            'No explicit ei heat exchange for non-idealized plasmas')

       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(.not.true_cell(i,j,k,iBlock)) CYCLE

          DtLocal = Cfl*time_BLK(i,j,k,iBlock)
          !\
          ! We apply the energy exchange rate for temperature,
          ! Ni*cTeTiExchangeRate/Te_G(i,j,k)**1.5
          !/
          ! For a hydrogen only, for ideal EOS only
          HeatExchange = cTeTiExchangeRate * &
               State_VGB(Rho_,i,j,k,iBlock)/Te_G(i,j,k)**1.5

          ! Point-implicit correction for stability: H' = H/(1+2*dt*H)
          HeatExchange = HeatExchange / (1 + 2.0*DtLocal*HeatExchange)

          HeatExchangePeP = HeatExchange &
               *(State_VGB(P_,i,j,k,iBlock) - State_VGB(Pe_,i,j,k,iBlock))

          ! Heat exchange for parallel ion pressure
          if(UseAnisoPressure)then
             HeatExchangePePpar = HeatExchange &
                  *(State_VGB(Ppar_,i,j,k,iBlock) &
                  - State_VGB(Pe_,i,j,k,iBlock))

             State_VGB(Ppar_,i,j,k,iBlock) = State_VGB(Ppar_,i,j,k,iBlock) &
                  - DtLocal*HeatExchangePePpar
          end if

          ! Heat exchange for the ions
          State_VGB(P_,i,j,k,iBlock) = State_VGB(P_,i,j,k,iBlock) &
               - DtLocal*HeatExchangePeP

          ! Heat exchange for the electrons
          State_VGB(Pe_,i,j,k,iBlock) = State_VGB(Pe_,i,j,k,iBlock) &
               + DtLocal*HeatExchangePeP
       end do; end do; end do

       call calc_energy_cell(iBlock)
    end do

    call test_stop(NameSub, DoTest)
  end subroutine calc_ei_heat_exchange
  !============================================================================

  subroutine get_impl_heat_cond_state(SemiAll_VCB, DconsDsemiAll_VCB, &
       DeltaSemiAll_VCB, DoCalcDeltaIn)

    ! Operator split, semi-implicit subroutines

    use ModVarIndexes,   ONLY: nVar, Rho_, p_, Pe_, Ppar_, Ehot_
    use ModAdvance,      ONLY: State_VGB, UseIdealEos, UseElectronPressure, &
         UseAnisoPressure, time_BLK, Source_VCB
    use ModFaceGradient, ONLY: set_block_field2, get_face_gradient
    use ModImplicit,     ONLY: nVarSemiAll, nBlockSemi, iBlockFromSemi_B, &
         iTeImpl
    use ModMain,         ONLY: Dt, time_accurate, Cfl
    use ModMultifluid,   ONLY: UseMultiIon, MassIon_I, ChargeIon_I, iRhoIon_I
    use ModNumConst,     ONLY: i_DD
    use ModPhysics,      ONLY: Si2No_V, No2Si_V, UnitTemperature_, &
         UnitEnergyDens_, UnitN_, UnitT_, AverageIonCharge, &
         InvGammaElectronMinus1
    use ModRadDiffusion, ONLY: UseHeatFluxLimiter
    use ModRadiativeCooling, ONLY: get_radiative_cooling
    use BATL_lib,        ONLY: IsCartesian, IsRzGeometry, &
         CellFace_DB, CellFace_DFB, FaceNormal_DDFB, Xyz_DGB
    use BATL_size,       ONLY: nI, nJ, nK, j0_, nJp1_, k0_, nKp1_, &
         MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    use ModHeatFluxCollisionless, ONLY: UseHeatFluxCollisionless, &
         get_gamma_collisionless
    use ModUserInterface ! user_material_properties
    use ModMain,         ONLY: UseFieldLineThreads
    use ModGeometry,     ONLY: far_field_BCs_BLK
    use ModParallel,     ONLY: NOBLK, NeiLev

    real, intent(out)  :: SemiAll_VCB(nVarSemiAll,nI,nJ,nK,nBlockSemi)
    real, intent(inout):: DconsDsemiAll_VCB(nVarSemiAll,nI,nJ,nK,nBlockSemi)
    real, optional, intent(out):: &
         DeltaSemiAll_VCB(nVarSemiAll,nI,nJ,nK,nBlockSemi)
    logical, optional, intent(in):: DoCalcDeltaIn

    logical :: DoCalcDelta
    ! Set it as allocatable
    real :: StarSemiAll_VCB(nVarSemiAll,nI,nJ,nK,nBlockSemi)
    real :: State_V(nVar)

    integer :: iDim, iDir, i, j, k, Di, Dj, Dk, iBlock, iBlockSemi, iP
    real :: GammaTmp
    real :: DtLocal
    real :: NumDens, NatomicSi, Natomic, TeTiRelaxSi, TeTiCoef, Cvi, TeSi, CvSi
    real :: HeatCoef, FreeStreamFlux, GradTe_D(3), GradTe
    real :: TeEpsilonSi = 1.0, TeEpsilon, RadCoolEpsilonR, RadCoolEpsilonL
    real :: bb_DD(nDim,nDim)
    logical :: IsNewBlockTe

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'get_impl_heat_cond_state'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    DoCalcDelta = .false.
    if(present(DoCalcDeltaIn)) DoCalcDelta=DoCalcDeltaIn

    TeEpsilon = TeEpsilonSi*Si2No_V(UnitTemperature_)

    iP = p_
    if(UseElectronPressure) iP = Pe_

    DtLocal = Dt

    do iBlockSemi = 1, nBlockSemi
       iBlock = iBlockFromSemi_B(iBlockSemi)

       IsNewBlockTe = .true.

       if(DoCalcDelta) then
          ! For the electron flux limiter, we need Te in the ghostcells
          if(UseMultiIon)then
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                State_V = State_VGB(:,i,j,k,iBlock)
                State_V(iP) = State_V(iP) + Source_VCB(iP,i,j,k,iBlock)
                Te_G(i,j,k) = State_V(Pe_)/sum( &
                     ChargeIon_I*State_V(iRhoIon_I)/MassIon_I)
             end do; end do; end do
          elseif(UseIdealEos)then
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                State_V = State_VGB(:,i,j,k,iBlock)
                State_V(iP) = State_V(iP) + Source_VCB(iP,i,j,k,iBlock)
                Te_G(i,j,k) = TeFraction &
                     *State_V(iP)/State_V(Rho_)
             end do; end do; end do
          else
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                State_V = State_VGB(:,i,j,k,iBlock)
                State_V(iP) = State_V(iP) + Source_VCB(iP,i,j,k,iBlock)
                call user_material_properties(State_V, &
                     i, j, k, iBlock, TeOut = TeSi)
                Te_G(i,j,k) = TeSi*Si2No_V(UnitTemperature_)
             end do; end do; end do
          end if
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             StarSemiAll_VCB(iTeImpl,i,j,k,iBlockSemi) = Te_G(i,j,k)
          enddo; enddo; enddo
       endif ! DoCalcDelta

       ! For the electron flux limiter, we need Te in the ghostcells
       if(UseMultiIon)then
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             Te_G(i,j,k) = State_VGB(Pe_,i,j,k,iBlock)/sum( &
                  ChargeIon_I*State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I)
          end do; end do; end do
       elseif(UseIdealEos)then
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             Te_G(i,j,k) = TeFraction &
                  *State_VGB(iP,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
          end do; end do; end do
       else
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  i, j, k, iBlock, TeOut = TeSi)
             Te_G(i,j,k) = TeSi*Si2No_V(UnitTemperature_)
          end do; end do; end do
       end if

       ! Store the electron temperature in SemiAll_VCB and the
       ! specific heat in DconsDsemiAll_VCB
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          SemiAll_VCB(iTeImpl,i,j,k,iBlockSemi) = Te_G(i,j,k)
          if(DoCalcDelta) &
               DeltaSemiAll_VCB(:,i,j,k,iBlockSemi) = &
               StarSemiAll_VCB(:,i,j,k,iBlockSemi) - &
               SemiAll_VCB(:,i,j,k,iBlockSemi)

          TeSi = Te_G(i,j,k)*No2Si_V(UnitTemperature_)

          if(UseIdealEos)then
             if(UseMultiIon)then
                ! Electron number density
                NumDens = sum( &
                     ChargeIon_I*State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I)
             else
                ! For simulations with Pe_, NumDens is the electron number
                ! density, while for single temperature simulations NumDens
                ! is the total number density
                NumDens = State_VGB(Rho_,i,j,k,iBlock)/TeFraction
             end if
             if(Ehot_ > 1 .and. UseHeatFluxCollisionless)then
                call get_gamma_collisionless(Xyz_DGB(:,i,j,k,iBlock), GammaTmp)
                DconsDsemiAll_VCB(iTeImpl,i,j,k,iBlockSemi)=NumDens/(GammaTmp-1)
             else
                DconsDsemiAll_VCB(iTeImpl,i,j,k,iBlockSemi)= &
                     InvGammaElectronMinus1*NumDens
             end if

             if(UseElectronPressure .and. .not.UseMultiIon)then
                Natomic = State_VGB(Rho_,i,j,k,iBlock)/MassIon_I(1)
                !\
                ! We apply the energy exchange rate for temperature,
                ! Ni*cTeTiExchangeRate/Te_G(i,j,k)**1.5
                ! to the electron energy density, therefore,we multiply by
                ! Ne/(\gamma -1)
                !/
                TeTiCoef = InvGammaElectronMinus1*(AverageIonCharge*Natomic)* &
                     (cTeTiExchangeRate*Natomic/Te_G(i,j,k)**1.5)
             end if

          else
             if(UseElectronPressure)then
                call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                     i, j, k, iBlock, TeOut=TeSi, CvOut = CvSi, &
                     NatomicOut = NatomicSi, TeTiRelaxOut = TeTiRelaxSi)

                Natomic = NatomicSi*Si2No_V(UnitN_)
                TeTiCoef = Natomic*TeTiRelaxSi/Si2No_V(UnitT_)
             else
                call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                     i, j, k, iBlock, TeOut=TeSi, CvOut = CvSi)
             end if

             DconsDsemiAll_VCB(iTeImpl,i,j,k,iBlockSemi) = &
                  CvSi*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_)
          end if

          if(.not.time_accurate) DtLocal = Cfl*time_BLK(i,j,k,iBlock)

          if(UseElectronPressure .and. .not.UseMultiIon)then
             Cvi = InvGammaElectronMinus1*Natomic
             PointCoef_CB(i,j,k,iBlock) = TeTiCoef/(1.0 + DtLocal*TeTiCoef/Cvi)
             PointImpl_VCB(1,i,j,k,iBlock) = State_VGB(p_,i,j,k,iBlock)/Natomic
             if(UseAnisoPressure)then
                PointImpl_VCB(2,i,j,k,iBlock) = &
                     State_VGB(Ppar_,i,j,k,iBlock)/Natomic
             end if
          end if

          if(DoRadCooling)then
             call get_radiative_cooling(i, j, k, iBlock, TeSi, &
                  CoolHeat_CB(i,j,k,iBlock))
             call get_radiative_cooling(i, j, k, iBlock, TeSi+TeEpsilonSi, &
                  RadCoolEpsilonR)
             call get_radiative_cooling(i, j, k, iBlock, TeSi-TeEpsilonSi, &
                  RadCoolEpsilonL)
             CoolHeatDeriv_CB(i,j,k,iBlock) = min(0.0, &
                  0.5*(RadCoolEpsilonR - RadCoolEpsilonL)/TeEpsilon)
          end if

       end do; end do; end do

       ! The following is because we entered the semi-implicit solve with
       ! first order ghost cells.
       State2_VG = State_VGB(:,:,:,:,iBlock)
       call set_block_field2(iBlock, nVar, State1_VG, State2_VG)

       ! Calculate the cell centered heat conduction tensor
       do k = k0_, nKp1_; do j = j0_, nJp1_; do i = 0, nI+1
          call get_heat_cond_tensor(State2_VG(:,i,j,k), i, j, k, iBlock)
       end do; end do; end do

       if(UseFieldLineThreads.and.far_field_BCs_BLK(iBlock))then
          ! First order BC at the outer boundaries
          if(NeiLev(1,iBlock)==NOBLK)then
             HeatCoef_G(0,:,:) = HeatCoef_G(1,:,:)
             bb_DDG(:,:,0,:,:) = bb_DDG(:,:,1,:,:)
          end if
          if(NeiLev(2,iBlock)==NOBLK)then
             HeatCoef_G(nI+1,:,:) = HeatCoef_G(nI,:,:)
             bb_DDG(:,:,nI+1,:,:) = bb_DDG(:,:,nI,:,:)
          end if
          if(nDim>=2)then
             if(NeiLev(3,iBlock)==NOBLK)then
                HeatCoef_G(:,0,:) = HeatCoef_G(:,1,:)
                bb_DDG(:,:,:,0,:) = bb_DDG(:,:,:,1,:)
             end if
             if(NeiLev(4,iBlock)==NOBLK)then
                HeatCoef_G(:,nJ+1,:) = HeatCoef_G(:,nJ,:)
                bb_DDG(:,:,:,nJ+1,:) = bb_DDG(:,:,:,nJ,:)
             end if
          end if
          if(nDim==3)then
             if(NeiLev(5,iBlock)==NOBLK)then
                HeatCoef_G(:,:,0) = HeatCoef_G(:,:,1)
                bb_DDG(:,:,:,:,0) = bb_DDG(:,:,:,:,1)
             end if
             if(NeiLev(6,iBlock)==NOBLK)then
                HeatCoef_G(:,:,nK+1) = HeatCoef_G(:,:,nK)
                bb_DDG(:,:,:,:,nK+1) = bb_DDG(:,:,:,:,nK)
             end if
          end if
       end if

       ! Average the cell centered heat conduction tensor to the faces
       ! and multiply with the area
       do iDim = 1, nDim
          Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
          do k = 1, nK+Dk; do j = 1, nJ+Dj; do i = 1, nI+Di
             bb_DD = 0.5*(bb_DDG(:nDim,:nDim,i,j,k) &
                  +       bb_DDG(:nDim,:nDim,i-Di,j-Dj,k-Dk))

             HeatCoef = 0.5*(HeatCoef_G(i,j,k) + HeatCoef_G(i-Di,j-Dj,k-Dk))

             if(UseHeatFluxLimiter)then
                call get_face_gradient(iDim, i, j, k, iBlock, &
                     IsNewBlockTe, Te_G, GradTe_D)

                GradTe = 0.0
                do iDir = 1, nDim
                   GradTe = GradTe + sum(bb_DD(iDir,:)*GradTe_D(:nDim)) &
                        *GradTe_D(iDir)
                end do
                GradTe = sqrt(GradTe)

                FreeStreamFlux = 0.5*(FreeStreamFlux_G(i,j,k) &
                     + FreeStreamFlux_G(i-Di,j-Dj,k-Dk))

                ! The threshold heat flux limiter model
                if(HeatCoef*GradTe > FreeStreamFlux) &
                     HeatCoef = FreeStreamFlux/GradTe
             end if

             if(IsCartesian)then
                HeatCond_DFDB(:nDim,i,j,k,iDim,iBlock) = &
                     CellFace_DB(iDim,iBlock)*HeatCoef*bb_DD(iDim,:)
             elseif(IsRzGeometry)then
                HeatCond_DFDB(:nDim,i,j,k,iDim,iBlock) = &
                     CellFace_DFB(iDim,i,j,k,iBlock)*HeatCoef*bb_DD(iDim,:)
             else
                do iDir = 1, nDim
                   HeatCond_DFDB(iDir,i,j,k,iDim,iBlock) = HeatCoef &
                        *sum( FaceNormal_DDFB(:,iDim,i,j,k,iBlock) &
                        *bb_DD(:,iDir) )
                end do
             end if
          end do; end do; end do
       end do

    end do

    call test_stop(NameSub, DoTest)
  contains
    !==========================================================================
    subroutine get_heat_cond_tensor(State_V, i, j, k, iBlock)

      use ModB0,         ONLY: B0_DGB
      use ModConst,      ONLY: cBoltzmann, cElectronmass
      use ModGeometry,   ONLY: r_BLK
      use ModMain,       ONLY: UseB0
      use ModNumConst,   ONLY: cTolerance
      use ModPhysics,    ONLY: UnitTemperature_, AverageIonCharge, &
           UnitPoynting_, ElectronGyroFreqCoef
      use ModRadDiffusion, ONLY: HeatFluxLimiter
      use ModVarIndexes, ONLY: Bx_, Bz_
      use ModRadiativeCooling, ONLY: &
           DoExtendTransitionRegion, extension_factor
      use ModUserInterface ! user_material_properties

      real, intent(in) :: State_V(nVar)
      integer, intent(in) :: i, j, k, iBlock

      ! Set HeatCoef_G(i,j,k) and bb_DDG(:,iDim,i,j,k) tensor and optionally
      ! the free stream flux FreeStreamFlux_G(i,j,k) for the cell center
      ! indexed by i,j,k of block iBlock based on its state State_V.
      ! The actual heat conduction tensor is
      ! HeatCoef_G(i,j,k)*bb_DDG(:,iDim,i,j,k)
      ! so in general bb_DDG(:,iDim,i,j,k) is the sum of the bb tensor AND
      ! an isotropic part. The isotropic part can be set as a fraction
      ! of the field aligned heat conduction.

      real :: TeSi, Te, NatomicSi, Ne, NeSi, Zav
      real :: HeatCoefSi, HeatCoef
      real :: Factor, r
      real :: Bnorm, B_D(3), Bunit_D(3)
      integer :: iDim
      !------------------------------------------------------------------------

      ! Calculate Ne, NeSi, Te, TeSi
      if(UseIdealEos)then
         if(UseMultiIon)then
            Te = State_V(Pe_)/sum(ChargeIon_I*State_V(iRhoIon_I)/MassIon_I)
            Ne = sum(ChargeIon_I*State_V(iRhoIon_I)/MassIon_I)
         else
            Te = TeFraction*State_V(iP)/State_V(Rho_)
            Ne = AverageIonCharge*State_V(Rho_)/MassIon_I(1)
         end if
         TeSi = Te*No2Si_V(UnitTemperature_)
         NeSi = Ne*No2Si_V(UnitN_)

      else
         call user_material_properties(State_V, i, j, k, iBlock, &
              TeOut=TeSi, NatomicOut=NatomicSi, AverageIonChargeOut=Zav)

         NeSi = Zav*NatomicSi
         Ne = NeSi*Si2No_V(UnitN_)
         Te = TeSi*Si2No_V(UnitTemperature_)
      end if

      if(DoWeakFieldConduction)then
         ! Initialize these public variables. The user can change them in
         ! user_material_properties when HeatCondOut is present
         ElectronCollisionRate = 0.0
         FractionFieldAligned  = -1.0
      end if

      if(DoUserHeatConduction .or. .not.UseIdealEos)then
         call user_material_properties(State_V, i, j, k, iBlock, &
              TeIn=TeSi, HeatCondOut=HeatCoefSi)
         if(HeatCoefSi < 0.0)then
            ! Spitzer conductivity if user sets negative HeatCoefSi
            HeatCoef = HeatCondPar*Te**2.5
         else
            HeatCoef = HeatCoefSi*Si2NoHeatCoef
         end if
      else
         ! Spitzer form for collisional regime
         HeatCoef = HeatCondPar*Te**2.5
      end if

      ! Artificial modified heat conduction for a smoother transition
      ! region, Linker et al. (2001)
      if(DoExtendTransitionRegion) HeatCoef = HeatCoef*extension_factor(TeSi)

      if(UseHeatFluxRegion)then
         r = r_BLK(i,j,k,iBlock)
         if(rCollisionless < 0.0)then
            Factor = 1.0/((r/rCollisional)**2 + 1)
         elseif(r <= rCollisional)then
            Factor = 1.0
         else
            Factor = exp(-((r-rCollisional)/(rCollisionless-rCollisional))**2)
         end if
         HeatCoef = Factor*HeatCoef
      end if

      HeatCoef_G(i,j,k) = HeatCoef

      if(UseHeatFluxLimiter)then
         FreeStreamFlux_G(i,j,k) = HeatFluxLimiter &
              *NeSi*cBoltzmann*TeSi*sqrt(cBoltzmann*TeSi/cElectronMass) &
              *Si2No_V(UnitPoynting_)
      end if

      if(UseB0)then
         B_D = State_V(Bx_:Bz_) + B0_DGB(:,i,j,k,iBlock)
      else
         B_D = State_V(Bx_:Bz_)
      end if

      ! The magnetic field should nowhere be zero. The following fix will
      ! turn the magnitude of the field direction to zero.
      Bnorm = sqrt(sum(B_D**2))
      Bunit_D = B_D / max( Bnorm, cTolerance )

      if(DoWeakFieldConduction)then
         ! If the user did not set the field-aligned fraction, calculate it as
         ! FractionFieldAligned = 1/(1 + ElectronCollisionRate/OmegaElectron)
         ! OmegaElectron = B*q_e/m_e = B*ElectronGyroFreqCoef
         ! The collision rate is the sum of the user supplied value
         ! (can be the neutral-electron rate) and
         ! the ion-electron collision rate
         if(FractionFieldAligned < 0.0)then
            ElectronCollisionRate = ElectronCollisionRate + &
                 ElectronIonCollisionCoef* &
                 sum(State_V(iRhoIon_I)/MassIon_I*ChargeIon_I**2)/Te**1.5
            FractionFieldAligned = &
                 1/(1 + ElectronCollisionRate/(Bnorm*ElectronGyroFreqCoef))
         end if
         do iDim = 1, 3
            bb_DDG(:,iDim,i,j,k) = ( &
                 FractionFieldAligned*Bunit_D*Bunit_D(iDim) &
                 + (1.0 - FractionFieldAligned)*i_DD(:,iDim) )
         end do
      else
         do iDim = 1, 3
            bb_DDG(:,iDim,i,j,k) = Bunit_D*Bunit_D(iDim)
         end do
      end if

    end subroutine get_heat_cond_tensor
    !==========================================================================

  end subroutine get_impl_heat_cond_state
  !============================================================================

  subroutine get_heat_conduction_rhs(iBlock, StateImpl_VG, Rhs_VC, IsLinear)

    use BATL_lib,        ONLY: store_face_flux, CellVolume_GB
    use ModSize,         ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    use ModAdvance,      ONLY: UseElectronPressure
    use ModFaceGradient, ONLY: get_face_gradient
    use ModImplicit,     ONLY: nVarSemi, iTeImpl, &
         FluxImpl_VXB, FluxImpl_VYB, FluxImpl_VZB
    use ModMain,         ONLY: nI, nJ, nK
    use ModNumConst,     ONLY: i_DD
    use ModMultiFluid,   ONLY: UseMultiIon
    use ModMain,         ONLY: UseFieldLineThreads
    use ModGeometry,     ONLY: far_field_BCs_BLK

    integer, intent(in) :: iBlock
    real, intent(inout) :: StateImpl_VG(nVarSemi,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    real, intent(out)   :: Rhs_VC(nVarSemi,nI,nJ,nK)
    logical, intent(in) :: IsLinear

    integer :: iDim, i, j, k, Di, Dj, Dk
    real :: FaceGrad_D(MaxDim)
    logical :: IsNewBlockHeatCond, UseFirstOrderBc
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'get_heat_conduction_rhs'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    IsNewBlockHeatCond = .true.
    UseFirstOrderBc = UseFieldLineThreads.and.far_field_BCs_BLK(iBlock)

    ! Calculate the electron thermal heat flux
    do iDim = 1, nDim
       Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
       do k = 1, nK+Dk; do j = 1, nJ+Dj; do i = 1, nI+Di

          ! Second-order accurate electron temperature gradient
          call get_face_gradient(iDim, i, j, k, iBlock, &
               IsNewBlockHeatCond, StateImpl_VG, FaceGrad_D,&
               UseFirstOrderBcIn=UseFirstOrderBC)

          FluxImpl_VFD(iTeImpl,i,j,k,iDim) = &
               -sum(HeatCond_DFDB(:,i,j,k,iDim,iBlock)*FaceGrad_D(:nDim))

       end do; end do; end do
    end do

    ! Store the fluxes at resolution changes for restoring conservation
    call store_face_flux(iBlock, nVarSemi, FluxImpl_VFD, &
         FluxImpl_VXB, FluxImpl_VYB, FluxImpl_VZB)

    Rhs_VC = 0.0

    do iDim = 1, nDim
       Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          Rhs_VC(:,i,j,k) = Rhs_VC(:,i,j,k) &
               -(FluxImpl_VFD(:,i+Di,j+Dj,k+Dk,iDim) &
               - FluxImpl_VFD(:,i,j,k,iDim))/CellVolume_GB(i,j,k,iBlock)
       end do; end do; end do
    end do

    if(IsLinear)then
       if(DoRadCooling)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             Rhs_VC(1,i,j,k) = Rhs_VC(1,i,j,k) &
                  + CoolHeatDeriv_CB(i,j,k,iBlock)*StateImpl_VG(iTeImpl,i,j,k)
          end do; end do; end do
       end if
    else
       if(DoRadCooling)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             Rhs_VC(1,i,j,k) = Rhs_VC(1,i,j,k) + CoolHeat_CB(i,j,k,iBlock)
          end do; end do; end do
       end if
    end if

    ! Point implicit source terms due to electron-ion energy exchange
    if(UseElectronPressure .and. .not.UseMultiIon)then
       if(IsLinear)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             Rhs_VC(1,i,j,k) = Rhs_VC(1,i,j,k) &
                  - PointCoef_CB(i,j,k,iBlock)*StateImpl_VG(iTeImpl,i,j,k)
          end do; end do; end do
       else
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             Rhs_VC(1,i,j,k) = Rhs_VC(1,i,j,k) + PointCoef_CB(i,j,k,iBlock) &
                  *(PointImpl_VCB(1,i,j,k,iBlock) &
                  - StateImpl_VG(iTeImpl,i,j,k))
          end do; end do; end do
       end if
    end if

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine get_heat_conduction_rhs
  !============================================================================

  subroutine add_jacobian_heat_cond(iBlock, nVarImpl, Jacobian_VVCI)

    ! This code can only be called from the semi-implicit scheme
    ! since this works on temperature and not energy or pressure,

    use ModAdvance,      ONLY: UseElectronPressure
    use ModFaceGradient, ONLY: set_block_jacobian_face, DcoordDxyz_DDFD
    use ModImplicit,     ONLY: UseNoOverlap, nStencil, iTeImpl
    use ModMain,         ONLY: nI, nJ, nK
    use ModNumConst,     ONLY: i_DD
    use BATL_lib,        ONLY: IsCartesianGrid, CellSize_DB, CellVolume_GB
    use ModMultiFluid,   ONLY: UseMultiIon

    integer, intent(in):: iBlock
    integer, intent(in):: nVarImpl
    real, intent(inout):: Jacobian_VVCI(nVarImpl,nVarImpl,nI,nJ,nK,nStencil)

    integer :: i, j, k, iDim, Di, Dj, Dk
    real :: DiffLeft, DiffRight, InvDcoord_D(nDim), InvDxyzVol_D(nDim), Coeff
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'add_jacobian_heat_cond'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    ! Contributions due to electron-ion energy exchange
    if(UseElectronPressure .and. .not.UseMultiIon)then
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          Jacobian_VVCI(1,1,i,j,k,1) = Jacobian_VVCI(1,1,i,j,k,1) &
               - PointCoef_CB(i,j,k,iBlock)
       end do; end do; end do
    end if

    if(DoRadCooling)then
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          Jacobian_VVCI(1,1,i,j,k,1) = Jacobian_VVCI(1,1,i,j,k,1) &
               + CoolHeatDeriv_CB(i,j,k,iBlock)
       end do; end do; end do
    end if

    InvDcoord_D = 1/CellSize_DB(:nDim,iBlock)

    if(.not.IsCartesianGrid) call set_block_jacobian_face(iBlock)

    ! the transverse diffusion is ignored in the Jacobian
    do iDim = 1, nDim
       Di = i_DD(iDim,1); Dj = i_DD(iDim,2); Dk = i_DD(iDim,3)
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          Coeff = InvDcoord_D(iDim)/CellVolume_GB(i,j,k,iBlock)
          if(IsCartesianGrid)then
             DiffLeft = Coeff*HeatCond_DFDB(iDim,i,j,k,iDim,iBlock)
             DiffRight = Coeff*HeatCond_DFDB(iDim,i+Di,j+Dj,k+Dk,iDim,iBlock)
          else
             InvDxyzVol_D = DcoordDxyz_DDFD(iDim,:nDim,i,j,k,iDim)*Coeff
             DiffLeft = sum(HeatCond_DFDB(:,i,j,k,iDim,iBlock)*InvDxyzVol_D)

             InvDxyzVol_D = DcoordDxyz_DDFD(iDim,:nDim,i+Di,j+Dj,k+Dk,iDim) &
                  *Coeff
             DiffRight = &
                  sum(HeatCond_DFDB(:,i+Di,j+Dj,k+Dk,iDim,iBlock)*InvDxyzVol_D)
          end if

          Jacobian_VVCI(iTeImpl,iTeImpl,i,j,k,1) = &
               Jacobian_VVCI(iTeImpl,iTeImpl,i,j,k,1) - (DiffLeft + DiffRight)

          if(UseNoOverlap)then
             if(  iDim==1.and.i==1  .or. &
                  iDim==2.and.j==1  .or. &
                  iDim==3.and.k==1)        DiffLeft = 0.0
             if(  iDim==1.and.i==nI .or. &
                  iDim==2.and.j==nJ .or. &
                  iDim==3.and.k==nK)       DiffRight = 0.0
          end if

          Jacobian_VVCI(iTeImpl,iTeImpl,i,j,k,2*iDim)   = &
               Jacobian_VVCI(iTeImpl,iTeImpl,i,j,k,2*iDim) + DiffLeft
          Jacobian_VVCI(iTeImpl,iTeImpl,i,j,k,2*iDim+1) = &
               Jacobian_VVCI(iTeImpl,iTeImpl,i,j,k,2*iDim+1) + DiffRight
       end do; end do; end do
    end do

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine add_jacobian_heat_cond
  !============================================================================

  subroutine update_impl_heat_cond(iBlock, iBlockSemi, &
       NewSemiAll_VC, OldSemiAll_VC, DconsDsemiAll_VC)

    ! The use ModVarIndexes has to be next to use ModAdvance for sake
    ! of the extremely advanced PGF90 12.9 compiler

    use ModAdvance,  ONLY: State_VGB, UseIdealEos, UseElectronPressure, &
         UseAnisoPressure, time_BLK
    use ModVarIndexes, ONLY: p_, Pe_, Ppar_, ExtraEint_, Ehot_
    use ModEnergy,   ONLY: calc_energy_cell
    use ModGeometry, ONLY: true_cell
    use ModImplicit, ONLY: nVarSemiAll, iTeImpl
    use ModMain,     ONLY: nI, nJ, nK, Dt, time_accurate, Cfl
    use ModPhysics,  ONLY: InvGammaElectronMinus1, GammaElectronMinus1, &
         InvGammaMinus1, GammaMinus1, No2Si_V, Si2No_V, UnitEnergyDens_, &
         UnitP_, ExtraEintMin, pMin_I, PeMin
    use ModHeatFluxCollisionless, ONLY: UseHeatFluxCollisionless, &
         get_gamma_collisionless
    use BATL_lib,    ONLY: Xyz_DGB
    use ModUserInterface ! user_material_properties
    use ModMultiFluid, ONLY: UseMultiIon, IonFirst_

    integer, intent(in) :: iBlock, iBlockSemi
    real, intent(in) :: NewSemiAll_VC(nVarSemiAll,nI,nJ,nK)
    real, intent(in) :: OldSemiAll_VC(nVarSemiAll,nI,nJ,nK)
    real, intent(in) :: DconsDsemiAll_VC(nVarSemiAll,nI,nJ,nK)

    integer :: i, j, k, iP
    real :: DeltaEinternal, Einternal, EinternalSi, PressureSi, pMin
    real :: GammaTmp
    real :: DtLocal
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'update_impl_heat_cond'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    if(UseElectronPressure)then
       iP = Pe_
       pMin = PeMin
    else
       iP = p_
       pMin = pMin_I(IonFirst_)
    end if

    DtLocal = Dt

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(.not.true_cell(i,j,k,iBlock)) CYCLE

       DeltaEinternal = DconsDsemiAll_VC(iTeImpl,i,j,k) &
            *(NewSemiAll_VC(iTeImpl,i,j,k) - OldSemiAll_VC(iTeImpl,i,j,k))

       if(UseIdealEos)then
          if(Ehot_ > 1 .and. UseHeatFluxCollisionless)then
             call get_gamma_collisionless(Xyz_DGB(:,i,j,k,iBlock), GammaTmp)

             ! Heat conduction is carried by electrons
             ! If single fluid w/o the electron equation, then
             ! InvGammaElectronMinus1 = InvGammaMinus1

             State_VGB(iP,i,j,k,iBlock) = max(pMin, (GammaTmp - 1) &
                  *(InvGammaElectronMinus1*State_VGB(iP,i,j,k,iBlock) &
                  + State_VGB(Ehot_,i,j,k,iBlock) + DeltaEinternal))
             State_VGB(Ehot_,i,j,k,iBlock) = State_VGB(iP,i,j,k,iBlock) &
                  *(1.0/(GammaTmp - 1) - InvGammaElectronMinus1)
          else
             State_VGB(iP,i,j,k,iBlock) = &
                  max(pMin, State_VGB(iP,i,j,k,iBlock) + &
                  GammaElectronMinus1*DeltaEinternal)
          end if
       else

          Einternal = InvGammaElectronMinus1*State_VGB(iP,i,j,k,iBlock) &
               + State_VGB(ExtraEint_,i,j,k,iBlock) + DeltaEinternal

          EinternalSi = Einternal*No2Si_V(UnitEnergyDens_)

          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               i, j, k, iBlock, &
               EinternalIn = EinternalSi, PressureOut = PressureSi)

          State_VGB(iP,i,j,k,iBlock) = PressureSi*Si2No_V(UnitP_)

          State_VGB(ExtraEint_,i,j,k,iBlock) = max(ExtraEintMin, &
               Einternal - InvGammaElectronMinus1*State_VGB(iP,i,j,k,iBlock))

       end if

       ! update ion pressure for energy exchange between ions and electrons
       if(UseElectronPressure .and. .not.UseMultiIon)then
          if(.not.time_accurate) DtLocal = Cfl*time_BLK(i,j,k,iBlock)
          Einternal = InvGammaMinus1*State_VGB(p_,i,j,k,iBlock) &
               + DtLocal*PointCoef_CB(i,j,k,iBlock) &
               *(NewSemiAll_VC(iTeImpl,i,j,k) - PointImpl_VCB(1,i,j,k,iBlock))

          State_VGB(p_,i,j,k,iBlock) = max(1e-30, GammaMinus1*Einternal)

          if(UseAnisoPressure)then
             Einternal = InvGammaMinus1*State_VGB(Ppar_,i,j,k,iBlock) &
                  + DtLocal*PointCoef_CB(i,j,k,iBlock) &
                  *(NewSemiAll_VC(iTeImpl,i,j,k) &
                  - PointImpl_VCB(2,i,j,k,iBlock))

             State_VGB(Ppar_,i,j,k,iBlock) = max(1e-30, GammaMinus1*Einternal)
          end if
       end if
    end do; end do; end do

    call calc_energy_cell(iBlock)

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine update_impl_heat_cond
  !============================================================================

end module ModHeatConduction
!==============================================================================
