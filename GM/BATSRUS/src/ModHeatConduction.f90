!^CFG COPYRIGHT UM
!^CFG FILE IMPLICIT
!==============================================================================
module ModHeatConduction

  use BATL_size, ONLY: nDim, MaxDim

  implicit none
  save

  private ! except

  ! Public methods
  public :: read_heatconduction_param
  public :: init_heat_conduction
  public :: get_heat_flux
  public :: get_ion_heat_flux
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

  ! Dimensionless heat conduction coefficients
  real :: HeatCondPar, IonHeatCondPar

  ! Parameters for heat conduction in regions of weak magnetic field
  logical :: DoWeakFieldConduction = .false.
  real :: BmodifySi = 1.0e-7, DeltaBmodifySi = 1.0e-8 ! modify about 1 mG
  real :: Bmodify, DeltaBmodify

  ! Parameters for heat flux region
  logical :: UseHeatFluxRegion = .false.
  real :: rCollisional, rCollisionless

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

  ! Heat conduction tensor
  real, allocatable :: HeatCond_DDG(:,:,:,:,:)
  ! Heat conduction dyad pre-multiplied by the face area
  real, allocatable :: HeatCond_DFDB(:,:,:,:,:,:)

  ! electron-ion energy exchange
  real, allocatable :: PointCoef_CB(:,:,:,:)
  real, allocatable :: PointImpl_CB(:,:,:,:)

contains

  !============================================================================

  subroutine read_heatconduction_param(NameCommand)

    use ModMain,      ONLY: UseHeatConduction, UseIonHeatConduction
    use ModReadParam, ONLY: read_var

    character(len=*), intent(in) :: NameCommand

    character(len=*), parameter :: &
         NameSub = 'ModHeatConduction::read_heatconduction_param'
    !--------------------------------------------------------------------------

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
       if(DoWeakFieldConduction)then
          call read_var('BmodifySi', BmodifySi)
          call read_var('DeltaBmodifySi', DeltaBmodifySi)
       end if

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

    case("#HEATFLUXREGION")
       call read_var('UseHeatFluxRegion', UseHeatFluxRegion)
       if(UseHeatFluxRegion)then
          call read_var('rCollisional', rCollisional)
          call read_var('rCollisionless', rCollisionless)
       end if

    case default
       call stop_mpi(NameSub//' invalid NameCommand='//NameCommand)
    end select

  end subroutine read_heatconduction_param

  !============================================================================

  subroutine init_heat_conduction

    use BATL_size,     ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nI, nJ, nK
    use ModAdvance,    ONLY: UseElectronPressure, UseIdealEos
    use ModConst,      ONLY: cBoltzmann, cElectronMass, cProtonMass, &
         cEps, cElectronCharge
    use ModImplicit,   ONLY: UseSemiImplicit, nVarSemi, iTeImpl
    use ModMain,       ONLY: MaxBlock, UseHeatConduction, UseIonHeatConduction
    use ModMultiFluid, ONLY: MassIon_I
    use ModNumConst,   ONLY: cTwoPi
    use ModPhysics,    ONLY: Si2No_V, UnitEnergyDens_, UnitTemperature_, &
         UnitU_, UnitX_, UnitB_, ElectronTemperatureRatio, AverageIonCharge
    use ModVarIndexes, ONLY: nVar

    real, parameter:: CoulombLog = 20.0
    real :: HeatCondParSi, IonHeatCondParSi

    character(len=*), parameter :: &
         NameSub = 'ModHeatConduction::init_heat_conduction'
    !--------------------------------------------------------------------------

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

    DoUserHeatConduction = .false.
    if(TypeHeatConduction == 'user') DoUserHeatConduction = .true.
    DoUserIonHeatConduction = .false.
    if(TypeIonHeatConduction == 'user') DoUserIonHeatConduction = .true.

    if(DoWeakFieldConduction)then
       Bmodify = BmodifySi*Si2No_V(UnitB_)
       DeltaBmodify = DeltaBmodifySi*Si2No_V(UnitB_)
    end if

    if(UseSemiImplicit.and..not.allocated(HeatCond_DFDB))then
       allocate( &
            State1_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK), &
            State2_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK), &
            FluxImpl_VFD(nVarSemi,nI+1,nJ+1,nK+1,nDim), &
            HeatCond_DDG(MaxDim,MaxDim,0:nI+1,0:nJ+1,0:nK+1), &
            HeatCond_DFDB(nDim,nI+1,nJ+1,nK+1,nDim,MaxBlock) )

       if(UseElectronPressure .and. .not.UseIdealEos)then
          allocate( &
               PointImpl_CB(nI,nJ,nK,MaxBlock), &
               PointCoef_CB(nI,nJ,nK,MaxBlock) )
       end if

       iTeImpl = 1
    end if

  end subroutine init_heat_conduction

  !============================================================================

  subroutine get_heat_flux(iDir, iFace, jFace, kFace, iBlock, &
       StateLeft_V, StateRight_V, Normal_D, HeatCondCoefNormal, HeatFlux)

    use BATL_size,       ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    use ModAdvance,      ONLY: State_VGB, UseIdealEos, UseElectronPressure
    use ModFaceGradient, ONLY: get_face_gradient
    use ModPhysics,      ONLY: inv_gm1, Si2No_V, UnitTemperature_, &
         UnitEnergyDens_
    use ModUser,         ONLY: user_material_properties
    use ModVarIndexes,   ONLY: nVar, Rho_, p_, Pe_

    integer, intent(in) :: iDir, iFace, jFace, kFace, iBlock
    real,    intent(in) :: StateLeft_V(nVar), StateRight_V(nVar), Normal_D(3)
    real,    intent(out):: HeatCondCoefNormal, HeatFlux

    integer :: i, j, k, iP
    real :: HeatCondL_D(3), HeatCondR_D(3), HeatCond_D(3), HeatCondFactor
    real :: FaceGrad_D(3), TeSi, CvL, CvR, CvSi

    character(len=*), parameter :: NameSub = 'ModHeatConduction::get_heat_flux'
    !--------------------------------------------------------------------------

    if(IsNewBlockHeatCond)then
       if(UseIdealEos .and. .not.DoUserHeatConduction)then
          iP = p_
          if(UseElectronPressure) iP = Pe_

          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             Te_G(i,j,k) = TeFraction*State_VGB(iP,i,j,k,iBlock) &
                  /State_VGB(Rho_,i,j,k,iBlock)
          end do; end do; end do
       else
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             call user_material_properties( &
                  State_VGB(:,i,j,k,iBlock),i,j,k,iBlock,TeOut=TeSi)
             Te_G(i,j,k) = TeSi*Si2No_V(UnitTemperature_)
          end do; end do; end do
       end if
    end if

    call get_face_gradient(iDir, iFace, jFace, kFace, iBlock, &
         IsNewBlockHeatCond, Te_G, FaceGrad_D)

    call get_heat_cond_coef(iDir, iFace, jFace, kFace, iBlock, &
         StateLeft_V, Normal_D, HeatCondL_D)
    call get_heat_cond_coef(iDir, iFace, jFace, kFace, iBlock, &
         StateRight_V, Normal_D, HeatCondR_D)

    HeatCond_D = 0.5*(HeatCondL_D + HeatCondR_D)

    if(UseHeatFluxRegion)then
         HeatCondFactor = heat_cond_factor(iDir, iFace, jFace, kFace, iBlock)
         HeatCond_D = HeatCond_D*HeatCondFactor
    end if

    HeatFlux = -sum(HeatCond_D*FaceGrad_D)

    ! get the heat conduction coefficient normal to the face for
    ! time step restriction
    if(UseIdealEos .and. .not.DoUserHeatConduction)then
       CvL = inv_gm1*StateLeft_V(Rho_)/TeFraction
       CvR = inv_gm1*StateRight_V(Rho_)/TeFraction
    else
       call user_material_properties(StateLeft_V, CvOut = CvSi)
       CvL = CvSi*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_)
       call user_material_properties(StateRight_V, CvOut = CvSi)
       CvR = CvSi*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_)
    end if
    HeatCondCoefNormal = sum(HeatCond_D*Normal_D)/min(CvL,CvR)

  end subroutine get_heat_flux

  !============================================================================

  subroutine get_heat_cond_coef(iDim, iFace, jFace, kFace, iBlock, &
       State_V, Normal_D, HeatCond_D)

    use ModAdvance,      ONLY: UseIdealEos, UseElectronPressure
    use ModB0,           ONLY: B0_DX, B0_DY, B0_DZ
    use ModMain,         ONLY: UseB0
    use ModNumConst,     ONLY: cTolerance
    use ModPhysics,      ONLY: No2Si_V, Si2No_V, UnitTemperature_, &
         UnitEnergyDens_, UnitU_, UnitX_
    use ModUser,         ONLY: user_material_properties
    use ModVarIndexes,   ONLY: nVar, Bx_, Bz_, Rho_, p_, Pe_
    use ModRadiativeCooling, ONLY: DoExtendTransitionRegion, extension_factor

    integer, intent(in) :: iDim, iFace, jFace, kFace, iBlock
    real, intent(in) :: State_V(nVar), Normal_D(3)
    real, intent(out):: HeatCond_D(3)

    real :: B_D(3), Bnorm, Bunit_D(3), TeSi, Te
    real :: HeatCoefSi, HeatCoef, FractionFieldAligned
    !--------------------------------------------------------------------------

    if(UseB0)then
       select case(iDim)
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

    if(UseIdealEos .and. .not.DoUserHeatConduction)then
       if(UseElectronPressure)then
          Te = TeFraction*State_V(Pe_)/State_V(Rho_)
       else
          Te = TeFraction*State_V(p_)/State_V(Rho_)
       end if
       TeSi = Te*No2Si_V(UnitTemperature_)

       ! Spitzer form for collisional regime
       HeatCoef = HeatCondPar*Te**2.5
    else
       ! Note we assume that the heat conduction formula for the
       ! ideal state is still applicable for the non-ideal state
       call user_material_properties(State_V, TeOut=TeSi, &
            HeatCondOut=HeatCoefSi)

       HeatCoef = HeatCoefSi &
              *Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_) &
              *Si2No_V(UnitU_)*Si2No_V(UnitX_)
    end if

    ! Artificial modified heat conduction for a smoother transition
    ! region, Linker et al. (2001)
    if(DoExtendTransitionRegion) HeatCoef = HeatCoef*extension_factor(TeSi)

    if(DoWeakFieldConduction)then
       FractionFieldAligned = 0.5*(1.0+tanh((Bnorm-Bmodify)/DeltaBmodify))

       HeatCond_D = HeatCoef*( &
            FractionFieldAligned*sum(Bunit_D*Normal_D)*Bunit_D &
            + (1.0 - FractionFieldAligned)*Normal_D )
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
    use ModPhysics,      ONLY: inv_gm1
    use ModVarIndexes,   ONLY: nVar, Rho_, p_

    integer, intent(in) :: iDir, iFace, jFace, kFace, iBlock
    real,    intent(in) :: StateLeft_V(nVar), StateRight_V(nVar), Normal_D(3)
    real,    intent(out):: HeatCondCoefNormal, HeatFlux

    integer :: i, j, k
    real :: HeatCondL_D(3), HeatCondR_D(3), HeatCond_D(3), HeatCondFactor
    real :: FaceGrad_D(3), CvL, CvR

    character(len=*), parameter :: &
         NameSub = 'ModHeatConduction::get_ion_heat_flux'
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
       CvL = inv_gm1*StateLeft_V(Rho_)/TiFraction
       CvR = inv_gm1*StateRight_V(Rho_)/TiFraction
    else
       call stop_mpi(NameSub// &
            ': no ion heat conduction yet for non-ideal eos')
    end if
    HeatCondCoefNormal = sum(HeatCond_D*Normal_D)/min(CvL,CvR)

  end subroutine get_ion_heat_flux

  !============================================================================

  subroutine get_ion_heat_cond_coef(iDim, iFace, jFace, kFace, iBlock, &
       State_V, Normal_D, HeatCond_D)

    use ModAdvance,    ONLY: UseIdealEos
    use ModB0,         ONLY: B0_DX, B0_DY, B0_DZ
    use ModMain,       ONLY: UseB0
    use ModNumConst,   ONLY: cTolerance
    use ModPhysics,    ONLY: Si2No_V, UnitEnergyDens_, UnitTemperature_, &
         UnitU_, UnitX_
    use ModUser,       ONLY: user_material_properties
    use ModVarIndexes, ONLY: nVar, Bx_, Bz_, Rho_, p_

    integer, intent(in) :: iDim, iFace, jFace, kFace, iBlock
    real, intent(in) :: State_V(nVar), Normal_D(3)
    real, intent(out):: HeatCond_D(3)

    real :: B_D(3), Bnorm, Bunit_D(3), Ti
    real :: IonHeatCoefSi, IonHeatCoef

    character(len=*), parameter :: &
         NameSub = 'ModHeatConduction::get_ion_heat_cond_coef'
    !--------------------------------------------------------------------------

    if(UseB0)then
       select case(iDim)
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

       call user_material_properties(State_V, IonHeatCondOut=IonHeatCoefSi)
       IonHeatCoef = IonHeatCoefSi &
              *Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_) &
              *Si2No_V(UnitU_)*Si2No_V(UnitX_)
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
    if(r <= rCollisional)then
       heat_cond_factor = 1.0
    else
       heat_cond_factor = &
            exp(-((r-rCollisional)/(rCollisionless-rCollisional))**2)
    end if

  end function heat_cond_factor

  !============================================================================
  ! Operator split, semi-implicit subroutines
  !============================================================================

  subroutine get_impl_heat_cond_state(StateImpl_VGB, DconsDsemi_VCB)

    use ModAdvance,      ONLY: State_VGB, UseIdealEos, UseElectronPressure, &
         time_BLK
    use ModFaceGradient, ONLY: set_block_field2
    use ModImplicit,     ONLY: nw, nImplBLK, impl2iBlk, iTeImpl, ImplCoeff
    use ModMain,         ONLY: nI, nJ, nK, MaxImplBlk, Dt, time_accurate, Cfl
    use ModNumConst,     ONLY: i_DD
    use ModPhysics,      ONLY: Si2No_V, UnitTemperature_, UnitEnergyDens_,&
         UnitN_, UnitT_, inv_gm1
    use ModUser,         ONLY: user_material_properties
    use ModVarIndexes,   ONLY: nVar, Rho_, p_, Pe_
    use BATL_lib,        ONLY: IsCartesian, IsRzGeometry, &
         CellFace_DB, CellFace_DFB, FaceNormal_DDFB

    real, intent(out) :: StateImpl_VGB(nw,0:nI+1,0:nJ+1,0:nK+1,MaxImplBlk)
    real, intent(inout) :: DconsDsemi_VCB(nw,nI,nJ,nK,MaxImplBlk)

    integer :: iDim, iDir, i, j, k, Di, Dj, Dk, iBlock, iImplBlock, iP
    real :: DtLocal
    real :: NatomicSi, Natomic, TeTiRelaxSi, TeTiCoef, Cvi, TeSi, CvSi

    integer, parameter :: jMin1 = 1 - min(1,nJ-1), jMax1 = nJ + min(1,nJ-1)
    integer, parameter :: kMin1 = 1 - min(1,nK-1), kMax1 = nK + min(1,nK-1)
    !--------------------------------------------------------------------------

    iP = p_
    if(UseElectronPressure) iP = Pe_

    DtLocal = Dt

    do iImplBlock = 1, nImplBLK
       iBlock = impl2iBLK(iImplBlock)

       ! Store the electron temperature in StateImpl_VGB and the
       ! specific heat in DconsDsemi_VCB
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(UseIdealEos .and. .not.DoUserHeatConduction)then
             StateImpl_VGB(iTeImpl,i,j,k,iImplBlock) = TeFraction &
                  *State_VGB(iP,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)

             DconsDsemi_VCB(iTeImpl,i,j,k,iImplBlock) = &
                  inv_gm1*State_VGB(Rho_,i,j,k,iBlock)/TeFraction
          else

             if(UseElectronPressure .and. .not.DoUserHeatConduction)then
                call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                     i, j, k, iBlock, TeOut=TeSi, CvOut = CvSi, &
                     NatomicOut = NatomicSi, TeTiRelaxOut = TeTiRelaxSi)

                Natomic = NatomicSi*Si2No_V(UnitN_)
                TeTiCoef = Natomic*TeTiRelaxSi/Si2No_V(UnitT_)
             else
                call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                     i, j, k, iBlock, TeOut=TeSi, CvOut = CvSi)
             end if

             StateImpl_VGB(iTeImpl,i,j,k,iImplBlock) = &
                  TeSi*Si2No_V(UnitTemperature_)

             DconsDsemi_VCB(iTeImpl,i,j,k,iImplBlock) = &
                  CvSi*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_)
          end if

          if(UseElectronPressure .and. .not.UseIdealEos)then
             if(.not.time_accurate) DtLocal = Cfl*time_BLK(i,j,k,iBlock)
             Cvi = inv_gm1*Natomic
             PointCoef_CB(i,j,k,iBlock) = TeTiCoef &
                  /(1.0 + ImplCoeff*DtLocal*TeTiCoef/Cvi)
             PointImpl_CB(i,j,k,iBlock) = State_VGB(p_,i,j,k,iBlock)/Natomic
          end if
       end do; end do; end do

       ! The following is because we entered the semi-implicit solve with
       ! first order ghost cells.
       State2_VG = State_VGB(:,:,:,:,iBlock)
       call set_block_field2(iBlock, nVar, State1_VG, State2_VG)

       ! Calculate the cell centered heat conduction tensor
       do k = kMin1, kMax1; do j = jMin1, jMax1; do i = 0, nI+1
          call get_heat_cond_tensor(State2_VG(:,i,j,k), &
               i, j, k, iBlock, HeatCond_DDG(:,:,i,j,k))
       end do; end do; end do

       ! Average the cell centered heat conduction tensor to the faces
       ! and multiply with the area
       do iDim = 1, nDim
          Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
          do k = 1, nK+Dk; do j = 1, nJ+Dj; do i = 1, nI+Di
             if(IsCartesian)then
                HeatCond_DFDB(:nDim,i,j,k,iDim,iBlock) = &
                     CellFace_DB(iDim,iBlock) &
                     *0.5*(HeatCond_DDG(iDim,:nDim,i,j,k) &
                     +     HeatCond_DDG(iDim,:nDim,i-Di,j-Dj,k-Dk))
             elseif(IsRzGeometry)then
                HeatCond_DFDB(:nDim,i,j,k,iDim,iBlock) = &
                     CellFace_DFB(iDim,i,j,k,iBlock) &
                     *0.5*(HeatCond_DDG(iDim,:nDim,i,j,k) &
                     +     HeatCond_DDG(iDim,:nDim,i-Di,j-Dj,k-Dk))
             else
                do iDir = 1, nDim
                   HeatCond_DFDB(iDir,i,j,k,iDim,iBlock) = &
                        0.5*sum( FaceNormal_DDFB(:,iDim,i,j,k,iBlock) &
                        *(HeatCond_DDG(:nDim,iDir,i,j,k) &
                        + HeatCond_DDG(:nDim,iDir,i-Di,j-Dj,k-Dk)) )
                end do
             end if
          end do; end do; end do
       end do

    end do

  contains

    subroutine get_heat_cond_tensor(State_V, i, j, k, iBlock, HeatCond_DD)

      use ModAdvance,    ONLY: UseIdealEos
      use ModB0,         ONLY: B0_DGB
      use ModGeometry,   ONLY: r_BLK
      use ModMain,       ONLY: UseB0
      use ModNumConst,   ONLY: cTolerance
      use ModPhysics,    ONLY: No2Si_V, Si2No_V, UnitTemperature_, &
           UnitEnergyDens_, UnitU_, UnitX_
      use ModUser,       ONLY: user_material_properties
      use ModVarIndexes, ONLY: nVar, Bx_, Bz_, Rho_
      use ModRadiativeCooling, ONLY: DoExtendTransitionRegion, extension_factor

      real, intent(in) :: State_V(nVar)
      integer, intent(in) :: i, j, k, iBlock
      real, intent(out) :: HeatCond_DD(3,3)

      real :: TeSi, Te
      real :: HeatCoefSi, HeatCoef
      real :: Factor, r
      real :: Bnorm, B_D(3), Bunit_D(3)
      real :: FractionFieldAligned
      integer :: iDim
      !------------------------------------------------------------------------

      if(UseIdealEos .and. .not.DoUserHeatConduction)then
         Te = TeFraction*State_V(iP)/State_V(Rho_)
         TeSi = Te*No2Si_V(UnitTemperature_)

         ! Spitzer form for collisional regime
         HeatCoef = HeatCondPar*Te**2.5
      else
         call user_material_properties(State_V, i, j, k, iBlock, &
              TeOut=TeSi, HeatCondOut=HeatCoefSi)

         HeatCoef = HeatCoefSi &
              *Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_) &
              *Si2No_V(UnitU_)*Si2No_V(UnitX_)
      end if

      ! Artificial modified heat conduction for a smoother transition
      ! region, Linker et al. (2001)
      if(DoExtendTransitionRegion) HeatCoef = HeatCoef*extension_factor(TeSi)

      if(UseHeatFluxRegion)then
         r = r_BLK(i,j,k,iBlock)
         if(r <= rCollisional)then
            Factor = 1.0
         else
            Factor = exp(-((r-rCollisional)/(rCollisionless-rCollisional))**2)
         end if
         HeatCoef = Factor*HeatCoef
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
         FractionFieldAligned = 0.5*(1.0+tanh((Bnorm-Bmodify)/DeltaBmodify))

         do iDim = 1, 3
            HeatCond_DD(:,iDim) = HeatCoef*( &
                 FractionFieldAligned*Bunit_D*Bunit_D(iDim) &
                 + (1.0 - FractionFieldAligned)*i_DD(:,iDim) )
         end do
      else
         do iDim = 1, 3
            HeatCond_DD(:,iDim) = HeatCoef*Bunit_D*Bunit_D(iDim)
         end do
      end if

    end subroutine get_heat_cond_tensor

  end subroutine get_impl_heat_cond_state

  !============================================================================

  subroutine get_heat_conduction_rhs(iBlock, StateImpl_VG, Rhs_VC, IsLinear)

    use BATL_lib,        ONLY: store_face_flux, CellVolume_GB
    use ModSize,         ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    use ModAdvance,      ONLY: UseElectronPressure, UseIdealEos
    use ModFaceGradient, ONLY: get_face_gradient
    use ModImplicit,     ONLY: nVarSemi, iTeImpl, FluxImpl_VXB, FluxImpl_VYB, &
         FluxImpl_VZB
    use ModMain,         ONLY: nI, nJ, nK
    use ModNumConst,     ONLY: i_DD

    integer, intent(in) :: iBlock
    real, intent(inout) :: StateImpl_VG(nVarSemi,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    real, intent(out)   :: Rhs_VC(nVarSemi,nI,nJ,nK)
    logical, intent(in) :: IsLinear

    integer :: iDim, i, j, k, Di, Dj, Dk
    real :: FaceGrad_D(MaxDim)
    logical :: IsNewBlockHeatCond
    !--------------------------------------------------------------------------

    IsNewBlockHeatCond = .true.

    ! Calculate the electron thermal heat flux
    do iDim = 1, nDim
       Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
       do k = 1, nK+Dk; do j = 1, nJ+Dj; do i = 1, nI+Di

          ! Second-order accurate electron temperature gradient
          call get_face_gradient(iDim, i, j, k, iBlock, &
               IsNewBlockHeatCond, StateImpl_VG, FaceGrad_D)

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

    ! Point implicit source terms due to electron-ion energy exchange
    if(UseElectronPressure .and. .not.UseIdealEos)then
       if(IsLinear)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             Rhs_VC(1,i,j,k) = Rhs_VC(1,i,j,k) &
                  - PointCoef_CB(i,j,k,iBlock)*StateImpl_VG(iTeImpl,i,j,k)
          end do; end do; end do
       else
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             Rhs_VC(1,i,j,k) = Rhs_VC(1,i,j,k) + PointCoef_CB(i,j,k,iBlock) &
                  *(PointImpl_CB(i,j,k,iBlock) - StateImpl_VG(iTeImpl,i,j,k))
          end do; end do; end do
       end if
    end if

  end subroutine get_heat_conduction_rhs

  !============================================================================

  subroutine add_jacobian_heat_cond(iBlock, Jacobian_VVCI)

    use ModAdvance,      ONLY: UseElectronPressure, UseIdealEos
    use ModFaceGradient, ONLY: set_block_jacobian_face, DcoordDxyz_DDFD
    use ModImplicit,     ONLY: iTeImpl, nVarSemi, UseNoOverlap
    use ModMain,         ONLY: nI, nJ, nK
    use ModNumConst,     ONLY: i_DD
    use BATL_lib,        ONLY: IsCartesianGrid, CellSize_DB, CellVolume_GB

    integer, parameter:: nStencil = 2*nDim + 1

    integer, intent(in) :: iBlock
    real, intent(inout) :: Jacobian_VVCI(nVarSemi,nVarSemi,nI,nJ,nK,nStencil)

    integer :: i, j, k, iDim, Di, Dj, Dk
    real :: DiffLeft, DiffRight, InvDcoord_D(nDim), InvDxyzVol_D(nDim), Coeff
    !--------------------------------------------------------------------------

    ! Contributions due to electron-ion energy exchange
    if(UseElectronPressure .and. .not.UseIdealEos)then
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          Jacobian_VVCI(1,1,i,j,k,1) = Jacobian_VVCI(1,1,i,j,k,1) &
               - PointCoef_CB(i,j,k,iBlock)
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

  end subroutine add_jacobian_heat_cond

  !============================================================================

  subroutine update_impl_heat_cond(iBlock, iImplBlock, StateImpl_VG)

    ! The use ModVarIndexes has to be next to use ModAdvance for sake
    ! of the extremely advanced PGF90 12.9 compiler 
    use ModAdvance,  ONLY: State_VGB, UseIdealEos, UseElectronPressure, &
         time_BLK
    use ModVarIndexes, ONLY: p_, Pe_, ExtraEint_
    use ModEnergy,   ONLY: calc_energy_cell
    use ModGeometry, ONLY: true_cell
    use ModImplicit, ONLY: nw, iTeImpl, DconsDsemi_VCB, ImplOld_VCB, ImplCoeff
    use ModMain,     ONLY: nI, nJ, nK, Dt, time_accurate, Cfl
    use ModPhysics,  ONLY: inv_gm1, gm1, No2Si_V, Si2No_V, UnitEnergyDens_, &
         UnitP_, ExtraEintMin
    use ModUser,     ONLY: user_material_properties


    integer, intent(in) :: iBlock, iImplBlock
    real, intent(in) :: StateImpl_VG(nw,nI,nJ,nK)

    integer :: i, j, k, iP
    real :: DeltaEinternal, Einternal, EinternalSi, PressureSi
    real :: DtLocal
    !--------------------------------------------------------------------------

    iP = p_
    if(UseElectronPressure) iP = Pe_

    DtLocal = Dt

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(.not.true_cell(i,j,k,iBlock)) CYCLE

       DeltaEinternal = DconsDsemi_VCB(iTeImpl,i,j,k,iImplBlock) &
            *(StateImpl_VG(iTeImpl,i,j,k) - ImplOld_VCB(iTeImpl,i,j,k,iBlock))

       if(UseIdealEos)then
          State_VGB(iP,i,j,k,iBlock) = max(1e-30, State_VGB(iP,i,j,k,iBlock) &
               + gm1*DeltaEinternal)
       else

          Einternal = inv_gm1*State_VGB(iP,i,j,k,iBlock) &
               + State_VGB(ExtraEint_,i,j,k,iBlock) + DeltaEinternal

          EinternalSi = Einternal*No2Si_V(UnitEnergyDens_)

          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               i, j, k, iBlock, &
               EinternalIn = EinternalSi, PressureOut = PressureSi)

          State_VGB(iP,i,j,k,iBlock) = PressureSi*Si2No_V(UnitP_)

          State_VGB(ExtraEint_,i,j,k,iBlock) = max(ExtraEintMin, &
               Einternal - inv_gm1*State_VGB(iP,i,j,k,iBlock))

       end if

       ! update ion pressure for energy exchange between ions and electrons
       if(UseElectronPressure .and. .not.UseIdealEos)then
          if(.not.time_accurate) DtLocal = Cfl*time_BLK(i,j,k,iBlock)
          Einternal = inv_gm1*State_VGB(p_,i,j,k,iBlock) &
               + DtLocal*PointCoef_CB(i,j,k,iBlock)*( &
               + ImplCoeff*(StateImpl_VG(iTeImpl,i,j,k) &
               -            PointImpl_CB(i,j,k,iBlock)) &
               + (1.0 - ImplCoeff)*(ImplOld_VCB(iTeImpl,i,j,k,iBlock) &
               -            PointImpl_CB(i,j,k,iBlock)) )

          State_VGB(p_,i,j,k,iBlock) = max(1e-30, gm1*Einternal)
       end if
    end do; end do; end do

    call calc_energy_cell(iBlock)

  end subroutine update_impl_heat_cond

end module ModHeatConduction
