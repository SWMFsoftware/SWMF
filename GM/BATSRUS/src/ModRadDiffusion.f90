!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!This code is a copyright protected software (c) 2002- University of Michigan

!============================================================================
module ModRadDiffusion

  ! Solve for gray or multigroup radiation diffusion 
  ! and/or isotropic electron heat conduction in the 
  ! non-relativistic fluid velocity limit.
  !
  ! The unknowns are radiation energy density and black body 
  ! radiation energy density (a*Te**4) representing the material
  ! (electron or ion) energy density.
  !
  ! Both unsplit and split semi-implicit schemes are implemented.
  ! The energy exchange terms are either evaluated point-implicitly,
  ! or implicitly. 
  !
  ! Notes: 
  !
  ! 1. If present, the electron/ion energy density must be the
  !    first implicit variable! 
  !
  ! 2. In the split semi-implicit scheme all variables
  !    use the point-implicit energy exchange!

  use ModImplicit,   ONLY: UseAccurateRadiation
  use ModVarIndexes, ONLY: p_, nWave
  use BATL_size, ONLY: nDim, MaxDim

  implicit none
  save

  private !except

  ! Public methods
  public :: read_rad_diffusion_param
  public :: init_rad_diffusion
  public :: get_radiation_energy_flux
  public :: calc_source_rad_diffusion
  public :: set_rad_outflow_bc
  public :: get_impl_rad_diff_state
  public :: get_rad_diffusion_rhs
  public :: add_jacobian_rad_diff
  public :: update_impl_rad_diff
  public :: set_rad_diff_range

  ! Logical for adding radiation diffusion
  logical, public :: IsNewBlockRadDiffusion = .true.
  logical, public :: IsNewTimestepRadDiffusion = .true.

  ! Logical for using the electron heat flux limiter
  logical, public :: UseHeatFluxLimiter = .false.

  ! Fraction of free streaming flux that is used in the threshold model
  real, public :: HeatFluxLimiter = 0.06

  ! Local variables --------------

  ! Parameters for radiation flux limiter
  logical          :: UseRadFluxLimiter  = .false.
  character(len=20):: TypeRadFluxLimiter = 'larsen'

  ! Minimum threshold for radiation temperature
  real             :: TradMinSi, EradMin_W(nWave)

  ! Coefficients for n-temperature electron-ion-radiation model
  real, allocatable :: DiffCoef_VFDB(:,:,:,:,:,:)
  real, allocatable :: RelaxCoef_VCB(:,:,:,:,:)
  real, allocatable :: DiffSemiCoef_VGB(:,:,:,:,:)
  real, allocatable :: PointCoef_VCB(:,:,:,:,:)
  real, allocatable :: PointImpl_VCB(:,:,:,:,:)
  real, allocatable :: PlanckWeight_WCB(:,:,:,:,:)

  ! Index indicating which implicit variables involve diffusion
  integer, allocatable :: iDiff_I(:)
  integer :: nDiff, iDiffMin, iDiffMax

  ! which one are heat conduction, which one radiation diffusion
  integer :: iDiffHeat = 0, iDiffRadFirst = 0, iDiffRadLast = 0
  integer :: iDiffRadMin = 0, iDiffRadMax = 0

  ! Index indicating which nRelax implicit variables involve energy exchange
  ! with the electrons
  integer, allocatable :: iRelax_I(:)
  integer :: nRelax, iRelaxMin, iRelaxMax

  ! Index which nPoint variables involve point implicit energy exchange
  ! with the state variable stored in PointImpl_VCB
  integer :: nPoint, iPointSemi

  ! radiation energy used for calculating radiative energy flux
  real, allocatable :: Erad_WG(:,:,:,:)
  ! temporary radiation energy array needed by set_block_field
  real, allocatable :: Erad1_WG(:,:,:,:)

  ! The electron heat flux limiter corrects the electron heat conduction if
  ! the electron temperature length scale is only a few collisonal mean free
  ! paths of the electrons or smaller.
  ! Only the threshold model is currently implemented.
  ! If we define the free streaming flux as
  !   F_fs = n_e*k_B*T_e*v_th, where v_th = sqrt(k_B*T_e/m_e)
  ! is a characteristic thermal velocity, then the threshold model
  ! limits the heat conduction flux F = -\kappa grad(Te), with heat
  ! conduction coefficient \kappa, by
  !   F = -min(\kappa, f*F_fs / |grad(Te)|) * grad(Te)
  ! Here, f is the electron flux limiter.
  !
  ! electron temperature array needed for calculating the elctron temperature
  ! gradient in the heat flux limiter
  real, allocatable :: Te_G(:,:,:)
  ! temporary electron temperature array needed by set_block_field
  real, allocatable :: Te1_G(:,:,:)

  logical :: UseTemperature

  real, allocatable :: FluxImpl_VFD(:,:,:,:,:)

contains

  subroutine read_rad_diffusion_param(NameCommand)

    use ModMain,      ONLY: UseRadDiffusion
    use ModReadParam, ONLY: read_var

    character(len=*), intent(in) :: NameCommand

    character(len=*), parameter :: NameSub = 'read_rad_diffusion_param'
    !--------------------------------------------------------------------------

    select case(NameCommand)
    case("#RADIATION")
       call read_var('UseRadDiffusion', UseRadDiffusion)

       if(UseRadDiffusion)then
          call read_var('UseRadFluxLimiter', UseRadFluxLimiter)
          if(UseRadFluxLimiter)then
             call read_var('TypeRadFluxLimiter', TypeRadFluxLimiter, &
                  IsLowerCase=.true.)

             select case(TypeRadFluxLimiter)
             case("larsen","sum","max")
             case default
                call stop_mpi(NameSub//': unknown TypeRadFluxLimiter='&
                     //TypeRadFluxLimiter)
             end select
          end if
          call read_var('TradMinSi', TradMinSi)
       end if

    case("#HEATFLUXLIMITER")
       call read_var('UseHeatFluxLimiter', UseHeatFluxLimiter)
       call read_var('HeatFluxLimiter', HeatFluxLimiter)

    case("#ACCURATERADIATION")
       call read_var('UseAccurateRadiation', UseAccurateRadiation)

    case default
       call stop_mpi(NameSub//' invalid NameCommand='//NameCommand)
    end select


  end subroutine read_rad_diffusion_param

  !============================================================================

  subroutine init_rad_diffusion

    use ModAdvance,     ONLY: nWave, UseElectronPressure
    use ModMain,        ONLY: UseRadDiffusion
    use ModSize,        ONLY: nI, nJ, nK, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
         MaxBlock
    use ModImplicit,    ONLY: UseFullImplicit, &
         UseSemiImplicit, UseSplitSemiImplicit, TypeSemiImplicit, &
         nVarSemi, iEradImpl, iTeImpl, iTrImplFirst, iTrImplLast
    use ModPhysics,     ONLY: Si2No_V, UnitTemperature_, &
         cRadiationNo
    use ModVarIndexes,  ONLY: nWave, WaveFirst_
    use ModWaves,       ONLY: UseWavePressure, GammaWave

    integer :: iVarImpl
    real :: TradMin

    character(len=*), parameter :: NameSub = "init_rad_diffusion"
    !------------------------------------------------------------------------

    if(allocated(Erad_WG)) RETURN

    ! Make sure that Erad_ is correct
    if(UseRadDiffusion)then
       if(nWave < 1) call stop_mpi(NameSub// &
            ": the number of wave bins should be 1 or more")

       TradMin   = TradMinSi*Si2No_V(UnitTemperature_)
       EradMin_W = cRadiationNo*TradMin**4/nWave
    end if

    if(UseFullImplicit)then
       allocate(Erad_WG(1,MinI:MaxI,MinJ:MaxJ,MinK:MaxK))

       nDiff = 1
       allocate(iDiff_I(nDiff))
       iDiff_I(1) = WaveFirst_
       nRelax = 1
    end if

    if(UseSemiImplicit)then

       allocate(Erad_WG(nWave,MinI:MaxI,MinJ:MaxJ,MinK:MaxK))
       allocate(Erad1_WG(nWave,MinI:MaxI,MinJ:MaxJ,MinK:MaxK))

       allocate(Te_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK))
       allocate(Te1_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK))

       ! Default to zero, unless reset
       iTeImpl = 0; iTrImplFirst = 0; iTrImplLast = 0;

       UseTemperature = .false.

       select case(TypeSemiImplicit)
       case('radiation')
          if(UseElectronPressure)then
             call stop_mpi(NameSub//": Te/=Ti requires Te,Ti relaxation")
          end if

          nDiff = nWave
          allocate(iDiff_I(nDiff))

          if(nWave == 1 .or. UseSplitSemiImplicit)then
             ! Radiation-material exchange is done point-implicitly
             iTrImplFirst = 1
             iTrImplLast  = nWave
             nPoint       = nWave
             nRelax = 0
             iDiff_I = 1
          else
             ! Radiation-material exchange is done with relaxation terms
             iTeImpl = 1; iTrImplFirst = 2; iTrImplLast = nWave + 1
             nPoint = 0
             nRelax = nWave
             allocate(iRelax_I(nRelax))
             iRelax_I = (/ (iVarImpl,iVarImpl=iTrImplFirst,iTrImplLast) /)
             iDiff_I = (/ (iVarImpl,iVarImpl=iTrImplFirst,iTrImplLast) /)
          end if
          iDiffHeat = 0; iDiffRadFirst = 1; iDiffRadLast = nWave

       case('radcond')
          iTeImpl = 1; iTrImplFirst = 2; iTrImplLast = iTrImplFirst + nWave - 1
          nDiff = 1 + nWave
          allocate(iDiff_I(nDiff))
          iDiffHeat = 1; iDiffRadFirst = 2; iDiffRadLast = nWave + 1

          if(UseSplitSemiImplicit)then
             UseTemperature = .true.
             iDiff_I = 1
             nRelax  = 0
             nPoint  = 1 + nWave
          else
             iDiff_I = &
                  (/ iTeImpl, (iVarImpl,iVarImpl=iTrImplFirst,iTrImplLast) /)
             nRelax = nWave
             allocate(iRelax_I(nRelax))
             iRelax_I = (/ (iVarImpl,iVarImpl=iTrImplFirst,iTrImplLast) /)
             if(UseElectronPressure)then
                nPoint = 1
             else
                nPoint = 0
             end if
          end if
       case('cond')
          UseTemperature = .true.
          iTeImpl = 1
          nDiff = 1
          allocate(iDiff_I(nDiff))
          iDiff_I = iTeImpl
          iDiffHeat = 1; iDiffRadFirst = 0; iDiffRadLast = 0
          nRelax = 0
          if(UseElectronPressure)then
             nPoint = 1
          else
             nPoint = 0
          end if

       end select

       iEradImpl = iTrImplFirst

       if(nPoint>0)then
          allocate(PointCoef_VCB(nPoint,nI,nJ,nK,MaxBlock))
          PointCoef_VCB = 0.0
          allocate(PointImpl_VCB(nPoint,nI,nJ,nK,MaxBlock))
          PointImpl_VCB = 0.0
       end if

       allocate(DiffSemiCoef_VGB(nDiff,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
       DiffSemiCoef_VGB = 0.0

       if(UseAccurateRadiation) &
            allocate(FluxImpl_VFD(nVarSemi,nI+1,nJ+1,nK+1,nDim))

       if(TypeSemiImplicit=='radiation' .or. TypeSemiImplicit=='radcond')then
          allocate(PlanckWeight_WCB(nWave,nI,nJ,nK,MaxBlock))
          PlanckWeight_WCB = 1.0
       end if
    end if

    if(nRelax>0)then
       allocate(RelaxCoef_VCB(nRelax,nI,nJ,nK,MaxBlock))
       RelaxCoef_VCB = 0.0
    end if

    allocate(DiffCoef_VFDB(nDiff,1:nI+1,1:nJ+1,1:nK+1,nDim,MaxBlock))
    DiffCoef_VFDB = 0.0 ! make sure all elements are initialized

    ! Setup for wave infrastructure
    if(UseRadDiffusion)then
       UseWavePressure = .true.
       GammaWave = 4.0/3.0       ! relativistic gamma for photon field
    end if

    ! For the unsplit semi-implicit scheme these are the defaults
    iRelaxMin   = 1
    iRelaxMax   = nRelax
    iDiffMin    = 1
    iDiffMax    = nDiff
    iDiffRadMin = iDiffRadFirst
    iDiffRadMax = iDiffRadLast
    iPointSemi  = 1

  end subroutine init_rad_diffusion

  !============================================================================

  subroutine get_radiation_energy_flux(iDir, i, j, k, iBlock, &
       StateLeft_V, StateRight_V, Normal_D, RadDiffCoef, EradFlux)

    use ModAdvance,      ONLY: State_VGB, Erad_
    use ModFaceGradient, ONLY: get_face_gradient
    use ModVarIndexes,   ONLY: nVar

    integer, intent(in) :: iDir, i, j, k, iBlock
    real,    intent(in) :: StateLeft_V(nVar), StateRight_V(nVar), Normal_D(3)
    real,    intent(out):: RadDiffCoef, EradFlux

    real :: FaceGrad_D(3), DiffCoefL, DiffCoefR
    !--------------------------------------------------------------------------

    if(IsNewBlockRadDiffusion) &
         Erad_WG(1,:,:,:) = State_VGB(Erad_,:,:,:,iBlock)

    call get_face_gradient(iDir, i, j, k, iBlock, &
         IsNewBlockRadDiffusion, Erad_WG, FaceGrad_D)

    if(IsNewTimestepRadDiffusion)then
       call get_diffusion_coef(StateLeft_V, DiffCoefL)
       call get_diffusion_coef(StateRight_V, DiffCoefR)

       RadDiffCoef = 0.5*(DiffCoefL + DiffCoefR)
       DiffCoef_VFDB(1,i,j,k,iDir,iBlock) = RadDiffCoef
    else
       RadDiffCoef = DiffCoef_VFDB(1,i,j,k,iDir,iBlock)
    end if

    EradFlux = -RadDiffCoef*sum(Normal_D*FaceGrad_D)

  contains

    subroutine get_diffusion_coef(State_V, DiffCoef)

      use ModAdvance,     ONLY: nWave
      use ModPhysics,     ONLY: Si2No_V, UnitX_, Clight
      use ModUser,        ONLY: user_material_properties

      real, intent(in) :: State_V(nVar)
      real, intent(out):: DiffCoef

      real :: OpacityRosselandSi_W(nWave), OpacityRosseland, Grad2ByErad2
      !------------------------------------------------------------------------

      call user_material_properties(State_V, i, j, k, iBlock, iDir, &
           OpacityRosselandOut_W = OpacityRosselandSi_W)

      OpacityRosseland = OpacityRosselandSi_W(1)/Si2No_V(UnitX_)

      if(UseRadFluxLimiter)then
         Grad2ByErad2 = sum(FaceGrad_D(1:nDim)**2)/State_V(Erad_)**2

         select case(TypeRadFluxLimiter)
         case("sum")
            DiffCoef = Clight/(3*OpacityRosseland + sqrt(Grad2ByErad2))
         case("max")
            DiffCoef = Clight/max(3*OpacityRosseland,sqrt(Grad2ByErad2))
         case("larsen")
            DiffCoef = Clight/sqrt(9*OpacityRosseland**2 + Grad2ByErad2)
         end select
      else
         DiffCoef = Clight/(3*OpacityRosseland)
      end if

    end subroutine get_diffusion_coef

  end subroutine get_radiation_energy_flux

  !============================================================================

  subroutine calc_source_rad_diffusion(iBlock)

    use ModAdvance,    ONLY: State_VGB, Source_VC, Erad_, nWave
    use ModConst,      ONLY: cLightSpeed
    use ModPhysics,    ONLY: cRadiationNo, Si2No_V, UnitTemperature_, UnitT_
    use ModMain,       ONLY: nI, nJ, nK
    use ModUser,       ONLY: user_material_properties
    use ModGeometry,   ONLY: true_cell

    integer, intent(in) :: iBlock

    integer :: i, j, k
    real :: TeSi, Te
    real :: AbsorptionEmission, OpacityPlanckSi_W(nWave)
    character(len=*), parameter:: NameSub = "calc_source_rad_diffusion"
    !------------------------------------------------------------------------

    do k=1,nK; do j=1,nJ; do i=1,nI

       if(.not.true_cell(i,j,k,iBlock)) CYCLE

       if(IsNewTimestepRadDiffusion)then
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               i, j, k, iBlock, OpacityPlanckOut_W = OpacityPlanckSi_W)

          RelaxCoef_VCB(1,i,j,k,iBlock) = &
               OpacityPlanckSi_W(1)*cLightSpeed/Si2No_V(UnitT_)
       end if

       call user_material_properties(State_VGB(:,i,j,k,iBlock), &
            i, j, k, iBlock, TeOut = TeSi)

       Te = TeSi*Si2No_V(UnitTemperature_)

       ! Source term due to absorption and emission
       ! Sigma_a*(cRadiation*Te**4-Erad)
       AbsorptionEmission =  RelaxCoef_VCB(1,i,j,k,iBlock) &
            *(cRadiationNo*Te**4 - State_VGB(Erad_,i,j,k,iBlock))

       ! dErad/dt = + AbsorptionEmission
       Source_VC(Erad_,i,j,k) = Source_VC(Erad_,i,j,k) + AbsorptionEmission

    end do; end do; end do

  end subroutine calc_source_rad_diffusion

  !============================================================================
  ! Semi-implicit interface
  !============================================================================

  subroutine get_impl_rad_diff_state(StateImpl_VGB,DconsDsemi_VCB)

    use BATL_lib,    ONLY: message_pass_cell, IsCartesian, IsRzGeometry, &
         CellSize_DB, CellFace_DFB, CellVolume_B
    use BATL_size,   ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
         j0_, nJp1_, k0_, nKp1_
    use ModAdvance,  ONLY: State_VGB, UseElectronPressure, nWave, WaveFirst_, &
         WaveLast_
    use ModConst,    ONLY: cBoltzmann
    use ModImplicit, ONLY: nw, nImplBlk, impl2iBlk, TypeSemiImplicit, &
         UseSplitSemiImplicit, iTeImpl, iTrImplFirst, iTrImplLast, ImplCoeff
    use ModMain,     ONLY: x_, y_, z_, nI, nJ, nK, MaxImplBlk, Dt
    use ModNumConst, ONLY: i_DD
    use ModPhysics,  ONLY: inv_gm1, Clight, cRadiationNo, UnitN_, &
         Si2No_V, UnitTemperature_, UnitEnergyDens_, UnitX_, UnitU_, UnitT_, &
         No2Si_V
    use ModUser,     ONLY: user_material_properties
    use ModGeometry, ONLY: TypeGeometry
    use ModParallel, ONLY: NOBLK, NeiLev

    real, intent(out):: StateImpl_VGB(nw,0:nI+1,j0_:nJp1_,k0_:nKp1_,MaxImplBlk)
    real, intent(inout):: DconsDsemi_VCB(nw,nI,nJ,nK,MaxImplBlk)

    integer :: iImplBlock, iBlock, i, j, k
    real :: OpacityPlanckSi_W(nWave), OpacityRosselandSi_W(nWave)
    real :: OpacityPlanck_W(nWave), CvSi, Cv, TeSi, Te
    real :: HeatCondSi, HeatCond, TeTiRelaxSi
    real :: InvDx2, InvDy2, InvDz2
    real :: NatomicSi, Natomic, Zav, CveSi, Cve, Cvi, Ti
    real :: TeTiCoef, TeTiCoefPrime

    integer :: iMin, iMax, jMin, jMax, kMin, kMax
    integer :: iDim, Di, Dj, Dk, iDiff
    real :: Coeff, InvDCoord_D(nDim)

    real :: PlanckSi_W(nWave), Planck_W(nWave), Planck

    ! Variables for the electron heat flux limiter
    logical :: IsNewBlockTe = .false.

    character(len=*), parameter:: NameSub='get_impl_rad_diff_state'
    !--------------------------------------------------------------------------

    ! For the electron flux limiter, we need Te in the ghostcells
    if(UseHeatFluxLimiter)then
       iMin = MinI; iMax = MaxI
       jMin = MinJ; jMax = MaxJ
       kMin = MinK; kMax = MaxK
    else
       iMin = 1; iMax = nI
       jMin = 1; jMax = nJ
       kMin = 1; kMax = nK
    end if

    do iImplBlock = 1, nImplBLK

       iBlock = impl2iBLK(iImplBlock)

       IsNewBlockRadDiffusion = .true.
       IsNewBlockTe = .true.

       if(iTeImpl > 0)then
          ! The ghost cells in Te_G are only needed for the electron heat flux
          ! limiter.
          do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  i, j, k, iBlock, TeOut = TeSi)
             Te_G(i,j,k) = TeSi*Si2No_V(UnitTemperature_)
          end do; end do; end do

          if(UseTemperature)then
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                StateImpl_VGB(iTeImpl,i,j,k,iImplBlock) = Te_G(i,j,k)
             end do; end do; end do
          else
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                StateImpl_VGB(iTeImpl,i,j,k,iImplBlock) = &
                     cRadiationNo*Te_G(i,j,k)**4
             end do; end do; end do
          end if
       end if

       if(iTrImplFirst > 0)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             StateImpl_VGB(iTrImplFirst:iTrImplLast,i,j,k,iImplBlock) = &
                  State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock)
          end do; end do; end do
       end if

       InvDx2 = 0.5/CellSize_DB(x_,iBlock)
       InvDy2 = 0.5/CellSize_DB(y_,iBlock)
       InvDz2 = 0.5/CellSize_DB(z_,iBlock)

       ! calculate coefficients for linearized energy exchange and diffusion
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(UseElectronPressure)then
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  i, j, k, iBlock, &
                  OpacityPlanckOut_W = OpacityPlanckSi_W, &
                  OpacityRosselandOut_W = OpacityRosselandSi_W, &
                  CvOut = CveSi, TeOut = TeSi, NatomicOut = NatomicSi, &
                  AverageIonChargeOut = Zav, &
                  HeatCondOut = HeatCondSi, TeTiRelaxOut = TeTiRelaxSi, &
                  PlanckOut_W = PlanckSi_W)

             Natomic = NatomicSi*Si2No_V(UnitN_)
             Ti = State_VGB(p_,i,j,k,iBlock)/Natomic
             Cvi = inv_gm1*Natomic
             Te = TeSi*Si2No_V(UnitTemperature_)
             Cve = CveSi*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_)
             if(.not.UseTemperature)then
                Cvi = Cvi/(4.0*cRadiationNo*Ti**3)
                Cve = Cve/(4.0*cRadiationNo*Te**3)
             end if
             TeTiCoef = Natomic*TeTiRelaxSi/Si2No_V(UnitT_)
          else
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  i, j, k, iBlock, &
                  OpacityPlanckOut_W = OpacityPlanckSi_W, &
                  OpacityRosselandOut_W = OpacityRosselandSi_W, &
                  AverageIonChargeOut = Zav, &
                  CvOut = CvSi, TeOut = TeSi, NatomicOut = NatomicSi, &
                  HeatCondOut=HeatCondSi, PlanckOut_W = PlanckSi_W)

             Te = TeSi*Si2No_V(UnitTemperature_)
             Cv = CvSi*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_)
             if(.not.UseTemperature) Cv = Cv/(4.0*cRadiationNo*Te**3)
          end if

          OpacityPlanck_W = OpacityPlanckSi_W/Si2No_V(UnitX_)
          Planck_W = PlanckSi_W*Si2No_V(UnitEnergyDens_)
          Planck   = cRadiationNo*Te**4

          select case(TypeSemiImplicit)
          case('radiation')
             if(nWave == 1)then
                ! This coefficient is cR'' = cR/(1+dt*cR*dPlanck/dEint)
                PointCoef_VCB(1,i,j,k,iBlock) = Clight*OpacityPlanck_W(1) &
                     /(1 + ImplCoeff*Dt*Clight*OpacityPlanck_W(1) / Cv)

                ! This is just the Planck function at time level * saved
                PointImpl_VCB(:,i,j,k,iBlock) = Planck
             elseif(UseSplitSemiImplicit)then
                ! Store information for Opacity_g*(Erad_g - B_g) term
                PointCoef_VCB(:,i,j,k,iBlock) = Clight*OpacityPlanck_W
                PointImpl_VCB(:,i,j,k,iBlock) = Planck_W
             else
                ! Unsplit multigroup
                DconsDsemi_VCB(iTeImpl,i,j,k,iImplBlock) = Cv
                RelaxCoef_VCB(:,i,j,k,iBlock) = Clight*OpacityPlanck_W
                PlanckWeight_WCB(:,i,j,k,iBlock) = Planck_W/Planck
             end if

          case('radcond')
             if(UseElectronPressure)then
                DconsDsemi_VCB(iTeImpl,i,j,k,iImplBlock) = Cve
             else
                DconsDsemi_VCB(iTeImpl,i,j,k,iImplBlock) = Cv
             end if

             if(UseSplitSemiImplicit)then
                if(UseElectronPressure)then
                   PointCoef_VCB(1,i,j,k,iBlock) = TeTiCoef &
                        /(1.0 + ImplCoeff*Dt*TeTiCoef/Cvi)
                   PointImpl_VCB(1,i,j,k,iBlock) = Ti
                end if

                ! The radiation-material energy exchange is point-implicit
                PointCoef_VCB(iTrImplFirst:,i,j,k,iBlock) &
                     = Clight*OpacityPlanck_W
                PointImpl_VCB(iTrImplFirst:,i,j,k,iBlock) &
                     = Planck_W
             else
                if(UseElectronPressure)then
                   ! Store variables for 
                   ! TeTiCoef*(B_e - B_i) = TeTiCoef'*(Te - Ti) term
                   TeTiCoefPrime = TeTiCoef/Cvi &
                        /(cRadiationNo*(Te+Ti)*(Te**2+Ti**2))
                   PointCoef_VCB(1,i,j,k,iBlock) = &
                        TeTiCoefPrime*Cvi/(1.0 + ImplCoeff*Dt*TeTiCoefPrime)
                   PointImpl_VCB(1,i,j,k,iBlock) = cRadiationNo*Ti**4
                end if

                ! The radiation-material energy exchange is semi-implicit
                RelaxCoef_VCB(:,i,j,k,iBlock) = Clight*OpacityPlanck_W
                PlanckWeight_WCB(:,i,j,k,iBlock) = Planck_W/Planck
             end if

          case('cond')
             if(UseElectronPressure)then
                DconsDsemi_VCB(iTeImpl,i,j,k,iImplBlock) = Cve

                PointCoef_VCB(1,i,j,k,iBlock) = TeTiCoef &
                     /(1.0 + ImplCoeff*Dt*TeTiCoef/Cvi)
                PointImpl_VCB(1,i,j,k,iBlock) = Ti
             else
                DconsDsemi_VCB(iTeImpl,i,j,k,iImplBlock) = Cv
             end if
          end select

          call get_diffusion_coef

       end do; end do; end do

       if(NeiLev(1,iBlock) == NOBLK)then
          i = 0
          do k = 1, nK; do j = 1, nJ
             call get_ghostcell_diffcoef
          end do; end do
       end if
       if(NeiLev(2,iBlock) == NOBLK)then
          i = nI + 1
          do k = 1, nK; do j = 1, nJ
             call get_ghostcell_diffcoef
          end do; end do
       end if
       if(nJ > 1 .and. NeiLev(3,iBlock) == NOBLK)then
          j = 0
          do k = 1, nK; do i = 1, nI
             call get_ghostcell_diffcoef
          end do; end do
       end if
       if(nJ > 1 .and. NeiLev(4,iBlock) == NOBLK)then
          j = nJ + 1
          do k = 1, nK; do i = 1, nI
             call get_ghostcell_diffcoef
          end do; end do
       end if
       if(nK > 1 .and. NeiLev(5,iBlock) == NOBLK)then
          k = 0
          do j = 1, nJ; do i = 1, nI
             call get_ghostcell_diffcoef
          end do; end do
       end if
       if(nK > 1 .and. NeiLev(6,iBlock) == NOBLK)then
          k = nK + 1
          do j = 1, nJ; do i = 1, nI
             call get_ghostcell_diffcoef
          end do; end do
       end if

    end do

    if(UseAccurateRadiation)then
       call message_pass_cell(nDiff, DiffSemiCoef_VGB)

       do iImplBlock = 1, nImplBLK
          iBlock = impl2iBLK(iImplBlock)

          do iDim = 1, nDim
             Di = i_DD(iDim,1); Dj = i_DD(iDim,2); Dk = i_DD(iDim,3)
             do k = 1, nK+Dk; do j = 1, nJ+Dj; do i = 1, nI+Di

                DiffCoef_VFDB(:,i,j,k,iDim,iBlock) = &
                     0.5*( DiffSemiCoef_VGB(:,i-Di,j-Dj,k-Dk,iBlock) &
                     +     DiffSemiCoef_VGB(:,i,j,k,iBlock) )
             end do; end do; end do
          end do
       end do

       RETURN
    end if

    ! Message pass to fill in ghost cells 
    call message_pass_cell(nDiff, DiffSemiCoef_VGB, nWidthIn=1, &
         nProlongOrderIn=1, DoSendCornerIn=.false., DoRestrictFaceIn=.true.)

    do iImplBlock = 1, nImplBLK
       iBlock = impl2iBLK(iImplBlock)
       ! Calculate face averaged values. Include geometric factors.

       call face_equal(1,2,nI,1,nJ,1,nK)
       if(NeiLev(1,iBlock)==0.or.NeiLev(1,iBlock)==NOBLK)then
          call face_equal(1,1,1,1,nJ,1,nK)
       else if(NeiLev(1,iBlock)==-1)then
          call face_left_coarse2fine(1,1,1,1,nJ,1,nK)
       else if(NeiLev(1,iBlock)==1)then
          call face_left_fine2coarse(1,1,1,1,nJ,1,nK)
       end if
       if(NeiLev(2,iBlock)==0.or.NeiLev(2,iBlock)==NOBLK)then
          call face_equal(1,nI+1,nI+1,1,nJ,1,nK)
       else if(NeiLev(2,iBlock)==-1)then
          call face_right_coarse2fine(1,nI+1,nI+1,1,nJ,1,nK)
       else if(NeiLev(2,iBlock)==1)then
          call face_right_fine2coarse(1,nI+1,nI+1,1,nJ,1,nK)
       end if

       if(nJ > 1)then
          call face_equal(2,1,nI,2,nJ,1,nK)
          if(NeiLev(3,iBlock)==0.or.NeiLev(3,iBlock)==NOBLK)then
             call face_equal(2,1,nI,1,1,1,nK)
          else if(NeiLev(3,iBlock)==-1)then
             call face_left_coarse2fine(2,1,nI,1,1,1,nK)
          else if(NeiLev(3,iBlock)==1)then
             call face_left_fine2coarse(2,1,nI,1,1,1,nK)
          end if
          if(NeiLev(4,iBlock)==0.or.NeiLev(4,iBlock)==NOBLK)then
             call face_equal(2,1,nI,nJ+1,nJ+1,1,nK)
          else if(NeiLev(4,iBlock)==-1)then
             call face_right_coarse2fine(2,1,nI,nJ+1,nJ+1,1,nK)
          else if(NeiLev(4,iBlock)==1)then
             call face_right_fine2coarse(2,1,nI,nJ+1,nJ+1,1,nK)
          end if
       end if

       if(nK > 1)then
          call face_equal(3,1,nI,1,nJ,2,nK)
          if(NeiLev(5,iBlock)==0.or.NeiLev(5,iBlock)==NOBLK)then
             call face_equal(3,1,nI,1,nJ,1,1)
          else if(NeiLev(5,iBlock)==-1)then
             call face_left_coarse2fine(3,1,nI,1,nJ,1,1)
          else if(NeiLev(5,iBlock)==1)then
             call face_left_fine2coarse(3,1,nI,1,nJ,1,1)
          end if
          if(NeiLev(6,iBlock)==0.or.NeiLev(6,iBlock)==NOBLK)then
             call face_equal(3,1,nI,1,nJ,nK+1,nK+1)
          else if(NeiLev(6,iBlock)==-1)then
             call face_right_coarse2fine(3,1,nI,1,nJ,nK+1,nK+1)
          else if(NeiLev(6,iBlock)==1)then
             call face_right_fine2coarse(3,1,nI,1,nJ,nK+1,nK+1)
          end if
       end if

       InvDcoord_D = 1/CellSize_DB(:nDim,iBlock)

       if(IsCartesian)then
          do iDim = 1, nDim
             ! FaceYZ/dx = Volume/dx^2
             Coeff = CellVolume_B(iBlock)*InvDcoord_D(iDim)**2
             Di = i_DD(iDim,1); Dj = i_DD(iDim,2); Dk = i_DD(iDim,3)
             do k=1,nK+Dk; do j=1,nJ+Dj; do i=1,nI+Di
                do iDiff = 1, nDiff
                   DiffCoef_VFDB(iDiff,i,j,k,iDim,iBlock) = &
                        Coeff*DiffCoef_VFDB(iDiff,i,j,k,iDim,iBlock)
                end do
             enddo; enddo; enddo
          end do

       elseif(IsRzGeometry)then
          do iDim = 1, nDim
             Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
             do k = 1, nK+Dk; do j = 1, nJ+Dj; do i = 1, nI+Di
                do iDiff = 1, nDiff
                   DiffCoef_VFDB(iDiff,i,j,k,iDim,iBlock) = &
                        InvDcoord_D(iDim)*CellFace_DFB(iDim,i,j,k,iBlock) &
                        *DiffCoef_VFDB(iDiff,i,j,k,iDim,iBlock)
                end do
             end do; end do; end do
          end do

       else
          call stop_mpi(NameSub//': unimplemented TypeGeometry='//TypeGeometry)
       end if
    end do

  contains

    subroutine face_equal(iDim,iMin,iMax,jMin,jMax,kMin,kMax)

      integer, intent(in) :: iDim, iMin, iMax, jMin, jMax, kMin, kMax

      !------------------------------------------------------------------------

      Di = i_DD(iDim,1); Dj = i_DD(iDim,2); Dk = i_DD(iDim,3)
      do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
         DiffCoef_VFDB(:,i,j,k,iDim,iBlock) = 0.5*( &
              DiffSemiCoef_VGB(:,i-Di,j-Dj,k-Dk,iBlock) &
              + DiffSemiCoef_VGB(:,i,j,k,iBlock))
      enddo; enddo; enddo

    end subroutine face_equal

    !==========================================================================

    subroutine face_left_coarse2fine(iDim,iMin,iMax,jMin,jMax,kMin,kMax)

      integer, intent(in) :: iDim, iMin, iMax, jMin, jMax, kMin, kMax

      !------------------------------------------------------------------------

      Di = i_DD(iDim,1); Dj = i_DD(iDim,2); Dk = i_DD(iDim,3)
      do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
         DiffCoef_VFDB(:,i,j,k,iDim,iBlock) =  &
              (DiffSemiCoef_VGB(:,i,j,k,iBlock) &
              + 2.0*DiffSemiCoef_VGB(:,i-Di,j-Dj,k-Dk,iBlock))*4.0/9.0
      enddo; enddo; enddo

    end subroutine face_left_coarse2fine

    !==========================================================================

    subroutine face_right_coarse2fine(iDim,iMin,iMax,jMin,jMax,kMin,kMax)

      integer, intent(in) :: iDim, iMin, iMax, jMin, jMax, kMin, kMax

      !------------------------------------------------------------------------

      Di = i_DD(iDim,1); Dj = i_DD(iDim,2); Dk = i_DD(iDim,3)
      do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
         DiffCoef_VFDB(:,i,j,k,iDim,iBlock) =  &
              (DiffSemiCoef_VGB(:,i-Di,j-Dj,k-Dk,iBlock) &
              + 2.0*DiffSemiCoef_VGB(:,i,j,k,iBlock))*4.0/9.0
      enddo; enddo; enddo

    end subroutine face_right_coarse2fine

    !==========================================================================

    subroutine face_left_fine2coarse(iDim,iMin,iMax,jMin,jMax,kMin,kMax)

      integer, intent(in) :: iDim, iMin, iMax, jMin, jMax, kMin, kMax

      integer :: iShift, jShift, kShift
      !------------------------------------------------------------------------

      Di = i_DD(iDim,1); Dj = i_DD(iDim,2); Dk = i_DD(iDim,3)
      iShift = 1-Di; jShift = min(1-Dj,nJ-1); kShift = min(1-Dk,nK-1)
      do k=kMin,kMax,2-Dk; do j=jMin,jMax,2-Dj; do i=iMin,iMax,2-Di
         do iDiff = 1, nDiff
            DiffCoef_VFDB(iDiff,i:i+iShift,j:j+jShift,k:k+kShift,iDim,iBlock) &
                 = (DiffSemiCoef_VGB(iDiff,i-Di,j-Dj,k-Dk,iBlock) &
                 + (4.0/2**nDim)*sum(DiffSemiCoef_VGB( &
                 iDiff,i:i+iShift,j:j+jShift,k:k+kShift,iBlock)))*2.0/9.0
         end do
      enddo; enddo; enddo

    end subroutine face_left_fine2coarse

    !==========================================================================

    subroutine face_right_fine2coarse(iDim,iMin,iMax,jMin,jMax,kMin,kMax)

      integer, intent(in) :: iDim, iMin, iMax, jMin, jMax, kMin, kMax

      integer :: iShift, jShift, kShift, i1, j1, k1
      !------------------------------------------------------------------------

      Di = i_DD(iDim,1); Dj = i_DD(iDim,2); Dk = i_DD(iDim,3)
      iShift = 1-Di; jShift = min(1-Dj,nJ-1); kShift = min(1-Dk,nK-1)
      do k=kMin,kMax,2-Dk; do j=jMin,jMax,2-Dj; do i=iMin,iMax,2-Di
         i1=i-Di; j1=j-Dj; k1=k-Dk
         do iDiff = 1, nDiff
            DiffCoef_VFDB(iDiff,i:i+iShift,j:j+jShift,k:k+kShift,iDim,iBlock) &
                 = (DiffSemiCoef_VGB(iDiff,i,j,k,iBlock) &
                 + (4.0/2**nDim)*sum(DiffSemiCoef_VGB( &
                 iDiff,i1:i1+iShift,j1:j1+jShift,k1:k1+kShift,iBlock)))*2.0/9.0
         end do
      enddo; enddo; enddo

    end subroutine face_right_fine2coarse

    !==========================================================================

    subroutine get_diffusion_coef

      use ModAdvance,      ONLY: nWave
      use ModConst,        ONLY: cElectronMass
      use ModFaceGradient, ONLY: set_block_field3

      real :: OpacityRosseland_W(nWave), DiffRad_W(nWave)
      real :: Grad2ByErad2_W(nWave)

      ! Variables for electron heat flux limiter
      real :: FreeStreamFluxSi, GradTe, GradTeSi, NeSi
      !------------------------------------------------------------------------

      if(iDiffRadFirst > 0)then
         OpacityRosseland_W = OpacityRosselandSi_W/Si2No_V(UnitX_)

         ! Calculate the cell centered diffusion coefficients
         if(UseRadFluxLimiter)then

            if(IsNewBlockRadDiffusion)then
               Erad_WG = State_VGB(WaveFirst_:WaveLast_,:,:,:,iBlock)
               call set_block_field3(iBlock, nWave, Erad1_WG, Erad_WG)

               IsNewBlockRadDiffusion = .false.
            end if

            Grad2ByErad2_W = &
                 ((Erad_WG(:,i+1,j,k) - Erad_WG(:,i-1,j,k))*InvDx2)**2
            if(nJ > 1) Grad2ByErad2_W = Grad2ByErad2_W + &
                 ((Erad_WG(:,i,j+1,k) - Erad_WG(:,i,j-1,k))*InvDy2)**2
            if(nK > 1) Grad2ByErad2_W = Grad2ByErad2_W + &
                 ((Erad_WG(:,i,j,k+1) - Erad_WG(:,i,j,k-1))*InvDz2)**2
            Grad2ByErad2_W = Grad2ByErad2_W/ Erad_WG(:,i,j,k)**2

            ! For now only Larsen's radiation flux limiter
            DiffRad_W = Clight/sqrt(9*OpacityRosseland_W**2 + Grad2ByErad2_W)
         else
            DiffRad_W = Clight/(3*OpacityRosseland_W)
         end if

         ! Store it for message passing
         DiffSemiCoef_VGB(iDiffRadFirst:iDiffRadLast,i,j,k,iBlock) = DiffRad_W
      end if

      if(iDiffHeat > 0)then

         if(UseHeatFluxLimiter)then
            ! Correct ghost cells as needed for gradient calculation
            if(IsNewBlockTe)then
               call set_block_field3(iBlock, 1, Te1_G, Te_G)

               IsNewBlockTe = .false.
            end if

            GradTe = ((Te_G(i+1,j,k) - Te_G(i-1,j,k))*InvDx2)**2
            if(nJ > 1) GradTe = GradTe + &
                 ((Te_G(i,j+1,k) - Te_G(i,j-1,k))*InvDy2)**2
            if(nK > 1) GradTe = GradTe + &
                 ((Te_G(i,j,k+1) - Te_G(i,j,k-1))*InvDz2)**2
            GradTe = sqrt(GradTe)

            GradTeSi = GradTe*No2Si_V(UnitTemperature_)/No2Si_V(UnitX_)

            NeSi = Zav*NatomicSi
            FreeStreamFluxSi = HeatFluxLimiter &
                 *NeSi*cBoltzmann*TeSi*sqrt(cBoltzmann*TeSi/cElectronMass)

            ! The threshold heat flux limiter model
            if(HeatCondSi*GradTeSi > FreeStreamFluxSi) &
                 HeatCondSi = FreeStreamFluxSi/GradTeSi
         end if

         HeatCond = HeatCondSi &
              *Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_) &
              *Si2No_V(UnitU_)*Si2No_V(UnitX_)

         if(UseTemperature)then
            DiffSemiCoef_VGB(iDiffHeat,i,j,k,iBlock) = HeatCond
         else
            Te = TeSi*Si2No_V(UnitTemperature_)
            DiffSemiCoef_VGB(iDiffHeat,i,j,k,iBlock) = HeatCond &
                 /(4.0*cRadiationNo*Te**3)
         end if
      end if

    end subroutine get_diffusion_coef

    !==========================================================================

    subroutine get_ghostcell_diffcoef

      !------------------------------------------------------------------------

      call user_material_properties(State_VGB(:,i,j,k,iBlock), &
           i, j, k, iBlock, &
           OpacityRosselandOut_W = OpacityRosselandSi_W, &
           HeatCondOut = HeatCondSi, TeOut = TeSi, &
           NatomicOut = NatomicSi, AverageIonChargeOut = Zav)

      call get_diffusion_coef

    end subroutine get_ghostcell_diffcoef

  end subroutine get_impl_rad_diff_state

  !============================================================================

  subroutine set_rad_outflow_bc(iSide, iBlock, iImplBlock, State_VG, IsLinear)

    use BATL_lib,    ONLY: CellSize_DB
    use ModImplicit, ONLY: nVarSemi
    use ModPhysics,  ONLY: Clight
    use ModSize,     ONLY: x_, y_, z_, &
         nI, nJ, nK, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
         j0_, nJp1_, k0_, nKp1_

    integer, intent(in) :: iSide, iBlock, iImplBlock
    real, intent(inout) :: State_VG(nVarSemi,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    logical, intent(in) :: IsLinear

    integer :: iVar, i, j, k, iDiff
    real :: Coef
    character(len=*), parameter :: NameSub='set_rad_outflow_bc'
    !--------------------------------------------------------------------------

    select case(iSide)
    case(1)
       do k = 1, nK; do j = 1, nJ
          do iDiff = iDiffRadMin, iDiffRadMax
             iVar = iDiff_I(iDiff)
             Coef = 2/Clight &
                  *DiffSemiCoef_VGB(iDiff,1,j,k,iBlock)/CellSize_DB(x_,iBlock)
             State_VG(iVar,0,j,k) = State_VG(iVar,1,j,k) &
                  *(Coef - 0.5)/(Coef + 0.5)
          end do
       end do; end do
    case(2)
       do k = 1, nK; do j = 1, nJ
          do iDiff = iDiffRadMin, iDiffRadMax
             iVar = iDiff_I(iDiff)
             Coef = 2/Clight &
                  *DiffSemiCoef_VGB(iDiff,nI,j,k,iBlock)/CellSize_DB(x_,iBlock)
             State_VG(iVar,nI+1,j,k) = State_VG(iVar,nI,j,k) &
                  *(Coef - 0.5)/(Coef + 0.5)
          end do
       end do; end do
    case(3)
       do k = 1, nK; do i = 1, nI
          do iDiff = iDiffRadMin, iDiffRadMax
             iVar = iDiff_I(iDiff)
             Coef = 2/Clight &
                  *DiffSemiCoef_VGB(iDiff,i,1,k,iBlock)/CellSize_DB(y_,iBlock)
             State_VG(iVar,i,j0_,k) = State_VG(iVar,i,1,k) &
                  *(Coef - 0.5)/(Coef + 0.5)
          end do
       end do; end do
    case(4)
       do k = 1, nK; do i = 1, nI
          do iDiff = iDiffRadMin, iDiffRadMax
             iVar = iDiff_I(iDiff)
             Coef = 2/Clight &
                  *DiffSemiCoef_VGB(iDiff,i,nJ,k,iBlock)/CellSize_DB(y_,iBlock)
             State_VG(iVar,i,nJp1_,k) = State_VG(iVar,i,nJ,k) &
                  *(Coef - 0.5)/(Coef + 0.5)
          end do
       end do; end do
    case(5)
       do j = 1, nJ; do i = 1, nI
          do iDiff = iDiffRadMin, iDiffRadMax
             iVar = iDiff_I(iDiff)
             Coef = 2/Clight &
                  *DiffSemiCoef_VGB(iDiff,i,j,1,iBlock)/CellSize_DB(z_,iBlock)
             State_VG(iVar,i,j,k0_) = State_VG(iVar,i,j,1) &
                  *(Coef - 0.5)/(Coef + 0.5)
          end do
       end do; end do
    case(6)
       do j = 1, nJ; do i = 1, nI
          do iDiff = iDiffRadMin, iDiffRadMax
             iVar = iDiff_I(iDiff)
             Coef = 2/Clight &
                  *DiffSemiCoef_VGB(iDiff,i,j,nK,iBlock)/CellSize_DB(z_,iBlock)
             State_VG(iVar,i,j,nKp1_) = State_VG(iVar,i,j,nK) &
                  *(Coef - 0.5)/(Coef + 0.5)
          end do
       end do; end do
    end select

  end subroutine set_rad_outflow_bc

  !============================================================================

  subroutine get_rad_diffusion_rhs(iBlock, StateImpl_VG, Rhs_VC, IsLinear)

    use BATL_lib,        ONLY: store_face_flux, IsCartesian, CellFace_DB, &
         CellFace_DFB, CellSize_DB, CellVolume_GB
    use BATL_size,       ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    use ModFaceGradient, ONLY: set_block_field3
    use ModImplicit,     ONLY: nVarSemi, iTeImpl, FluxImpl_VXB, FluxImpl_VYB, &
         FluxImpl_VZB
    use ModLinearSolver, ONLY: pDotADotPPe, UsePDotADotP
    use ModMain,         ONLY: nI, nJ, nK
    use ModParallel,     ONLY: NeiLev, NOBLK
    use ModNumConst,     ONLY: i_DD
    use ModGeometry,     ONLY: true_cell

    integer, intent(in) :: iBlock
    real, intent(inout) :: StateImpl_VG(nVarSemi,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    real, intent(out)   :: Rhs_VC(nVarSemi,nI,nJ,nK)
    logical, intent(in) :: IsLinear

    real :: InvDcoord_D(MaxDim)
    real :: Area, EnergyExchange
    integer :: iDim, i, j, k, Di, Dj, Dk, iDiff, iRelax, iVar

    real :: StateImpl1_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    real :: StateImpl_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

    integer, parameter:: &
         jMin = 1 - 2*min(1,nJ-1), jMax = nJ + 2*min(1,nJ-1), &
         kMin = 1 - 2*min(1,nK-1), kMax = nK + 2*min(1,nK-1)

    character(len=*), parameter :: NameSub='get_rad_diffusion_rhs'
    !--------------------------------------------------------------------------

    if(UseAccurateRadiation)then

       InvDcoord_D = 1.0/CellSize_DB(:,iBlock)

       do iDiff = iDiffMin, iDiffMax
          iVar = iDiff_I(iDiff)

          do k = kMin, kMax; do j = jMin, jMax; do i = MinI, MaxI
             if(.not.true_cell(i,j,k,iBlock)) CYCLE
             StateImpl_G(i,j,k) = StateImpl_VG(iVar,i,j,k)
          end do; end do; end do

          call set_block_field3(iBlock, 1, StateImpl1_G, StateImpl_G)

          do iDim = 1, nDim
             if(IsCartesian) Area = CellFace_DB(iDim,iBlock)
             Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
             do k = 1, nK+Dk; do j = 1, nJ+Dj; do i = 1, nI+Di
                if(.not.true_cell(i,j,k,iBlock)) CYCLE
                if(.not.IsCartesian) Area = CellFace_DFB(iDim,i,j,k,iBlock)

                FluxImpl_VFD(iVar,i,j,k,iDim) = &
                     -Area*DiffCoef_VFDB(iDiff,i,j,k,iDim,iBlock) &
                     *(StateImpl_G(i,j,k) - StateImpl_G(i-Di,j-Dj,k-Dk)) &
                     *InvDcoord_D(iDim)
             end do; end do; end do
          end do
       end do

       ! Store the fluxes at resolution changes for restoring conservation
       call store_face_flux(iBlock, nVarSemi, FluxImpl_VFD, &
            FluxImpl_VXB, FluxImpl_VYB, FluxImpl_VZB)

       Rhs_VC = 0.0

       do iDim = 1, nDim
          Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.true_cell(i,j,k,iBlock)) CYCLE
             do iDiff = iDiffMin, iDiffMax
                iVar = iDiff_I(iDiff)
                Rhs_VC(iVar,i,j,k) = Rhs_VC(iVar,i,j,k) &
                     -(FluxImpl_VFD(iVar,i+Di,j+Dj,k+Dk,iDim) &
                     - FluxImpl_VFD(iVar,i,j,k,iDim)) &
                     /CellVolume_GB(i,j,k,iBlock)
             end do
          end do; end do; end do
       end do

    else

       if(NeiLev(1,iBlock)==1) call correct_left_ghostcell(1,0,0,1,nJ,1,nK)
       if(NeiLev(2,iBlock)==1) &
            call correct_right_ghostcell(1,nI+1,nI+1,1,nJ,1,nK)
       if(nJ > 1)then
          if(NeiLev(3,iBlock)==1) &
               call correct_left_ghostcell(2,1,nI,0,0,1,nK)
          if(NeiLev(4,iBlock)==1) &
               call correct_right_ghostcell(2,1,nI,nJ+1,nJ+1,1,nK)
       end if
       if(nK > 1)then
          if(NeiLev(5,iBlock)==1) &
               call correct_left_ghostcell(3,1,nI,1,nJ,0,0)
          if(NeiLev(6,iBlock)==1) &
               call correct_right_ghostcell(3,1,nI,1,nJ,nK+1,nK+1)
       end if

       Rhs_VC = 0.0

       if(nDim == 1)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.true_cell(i,j,k,iBlock)) CYCLE
             do iDiff = iDiffMin, iDiffMax
                iVar = iDiff_I(iDiff)
                Rhs_VC(iVar,i,j,k) = ( &
                     + DiffCoef_VFDB(iDiff,i+1,j,k,1,iBlock)* &
                     (   StateImpl_VG(iVar,i+1,j,k)   &
                     -   StateImpl_VG(iVar,i  ,j,k))  &
                     - DiffCoef_VFDB(iDiff,i  ,j,k,1,iBlock)* &
                     (   StateImpl_VG(iVar,i  ,j,k)   &
                     -   StateImpl_VG(iVar,i-1,j,k)) ) &
                     /CellVolume_GB(i,j,k,iBlock)
             end do
          end do; end do; end do
       elseif(nDim == 2)then
          ! No flux from Z direction
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.true_cell(i,j,k,iBlock)) CYCLE
             do iDiff = iDiffMin, iDiffMax
                iVar = iDiff_I(iDiff)
                Rhs_VC(iVar,i,j,k) = ( &
                     + DiffCoef_VFDB(iDiff,i+1,j,k,1,iBlock)* &
                     (   StateImpl_VG(iVar,i+1,j,k)   &
                     -   StateImpl_VG(iVar,i  ,j,k))  &
                     - DiffCoef_VFDB(iDiff,i  ,j,k,1,iBlock)* &
                     (   StateImpl_VG(iVar,i  ,j,k)   &
                     -   StateImpl_VG(iVar,i-1,j,k))  &
                     + DiffCoef_VFDB(iDiff,i,j+1,k,2,iBlock)* &
                     (   StateImpl_VG(iVar,i,j+1,k)   &
                     -   StateImpl_VG(iVar,i,j  ,k))  &
                     - DiffCoef_VFDB(iDiff,i,j  ,k,2,iBlock)* &
                     (   StateImpl_VG(iVar,i,j  ,k)   &
                     -   StateImpl_VG(iVar,i,j-1,k)) ) &
                     /CellVolume_GB(i,j,k,iBlock)
             end do
          end do; end do; end do
       else
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.true_cell(i,j,k,iBlock)) CYCLE
             do iDiff = iDiffMin, iDiffMax
                iVar = iDiff_I(iDiff)
                Rhs_VC(iVar,i,j,k) = ( &
                     DiffCoef_VFDB(iDiff,i+1,j,k,1,iBlock)*   &
                     (   StateImpl_VG(iVar,i+1,j,k)   &
                     -   StateImpl_VG(iVar,i  ,j,k))  &
                     - DiffCoef_VFDB(iDiff,i  ,j,k,1,iBlock)* &
                     (   StateImpl_VG(iVar,i  ,j,k)   &
                     -   StateImpl_VG(iVar,i-1,j,k))  &
                     + DiffCoef_VFDB(iDiff,i,j+1,k,2,iBlock)* &
                     (   StateImpl_VG(iVar,i,j+1,k)   &
                     -   StateImpl_VG(iVar,i,j  ,k))  &
                     - DiffCoef_VFDB(iDiff,i,j  ,k,2,iBlock)* &
                     (   StateImpl_VG(iVar,i,j  ,k)   &
                     -   StateImpl_VG(iVar,i,j-1,k))  &
                     + DiffCoef_VFDB(iDiff,i,j,k+1,3,iBlock)* &
                     (   StateImpl_VG(iVar,i,j,k+1)   &
                     -   StateImpl_VG(iVar,i,j,k  ))  &
                     - DiffCoef_VFDB(iDiff,i,j,k  ,3,iBlock)* &
                     (   StateImpl_VG(iVar,i,j,k  )   &
                     -   StateImpl_VG(iVar,i,j,k-1)) ) &
                     /CellVolume_GB(i,j,k,iBlock)
             end do
          end do; end do; end do

       end if

    end if

    ! Source terms due to energy exchange
    if(nRelax > 0)then
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(.not.true_cell(i,j,k,iBlock)) CYCLE
          do iRelax = iRelaxMin, iRelaxMax
             iVar = iRelax_I(iRelax)

             EnergyExchange = RelaxCoef_VCB(iRelax,i,j,k,iBlock) &
                  *(PlanckWeight_WCB(iRelax,i,j,k,iBlock) &
                  * StateImpl_VG(iTeImpl,i,j,k) - StateImpl_VG(iVar,i,j,k))

             ! dEvar/dt = + EnergyExchange
             Rhs_VC(iVar,i,j,k)    = Rhs_VC(iVar,i,j,k)    + EnergyExchange

             ! dEe/dt   = - EnergyExchange
             Rhs_VC(iTeImpl,i,j,k) = Rhs_VC(iTeImpl,i,j,k) - EnergyExchange
          end do
       end do; end do; end do
    end if

    ! Point implicit source terms due to energy exchange
    if(nPoint > 0)then
       if(IsLinear)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.true_cell(i,j,k,iBlock)) CYCLE
             Rhs_VC(1,i,j,k) = Rhs_VC(1,i,j,k) &
                  - PointCoef_VCB(iPointSemi,i,j,k,iBlock) &
                  *StateImpl_VG(1,i,j,k)
          end do; end do; end do
       else
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.true_cell(i,j,k,iBlock)) CYCLE
             Rhs_VC(1,i,j,k) = Rhs_VC(1,i,j,k) &
                  + PointCoef_VCB(iPointSemi,i,j,k,iBlock) &
                  *(PointImpl_VCB(iPointSemi,i,j,k,iBlock) &
                  - StateImpl_VG(1,i,j,k))
          end do; end do; end do
       end if
    end if

    if(UsePDotADotP)then
       if(nDim == 1)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.true_cell(i,j,k,iBlock)) CYCLE
             do iDiff = iDiffMin, iDiffMax
                iVar = iDiff_I(iDiff)
                pDotADotPPe = pDotADotPPe  + 0.5 *(&
                     DiffCoef_VFDB(iDiff,i+1,j,k,1,iBlock)* &
                     (   StateImpl_VG(iVar,i+1,j,k)   &
                     -   StateImpl_VG(iVar,i  ,j,k))**2  &
                     + DiffCoef_VFDB(iDiff,i  ,j,k,1,iBlock)* &
                     (   StateImpl_VG(iVar,i  ,j,k)   &
                     -   StateImpl_VG(iVar,i-1,j,k))**2 )
             end do
          end do; end do; end do
       elseif(nDim == 2)then
          ! No flux from Z direction
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.true_cell(i,j,k,iBlock)) CYCLE
             do iDiff = iDiffMin, iDiffMax
                iVar = iDiff_I(iDiff)
                pDotADotPPe = pDotADotPPe  + 0.5 *(&
                     DiffCoef_VFDB(iDiff,i+1,j,k,1,iBlock)* &
                     (   StateImpl_VG(iVar,i+1,j,k)   &
                     -   StateImpl_VG(iVar,i  ,j,k))**2  &
                     + DiffCoef_VFDB(iDiff,i  ,j,k,1,iBlock)* &
                     (   StateImpl_VG(iVar,i  ,j,k)   &
                     -   StateImpl_VG(iVar,i-1,j,k))**2  &
                     + DiffCoef_VFDB(iDiff,i,j+1,k,2,iBlock)* &
                     (   StateImpl_VG(iVar,i,j+1,k)   &
                     -   StateImpl_VG(iVar,i,j  ,k))**2  &
                     + DiffCoef_VFDB(iDiff,i,j  ,k,2,iBlock)* &
                     (   StateImpl_VG(iVar,i,j  ,k)   &
                     -   StateImpl_VG(iVar,i,j-1,k))**2 )
             end do
          end do; end do; end do
       else
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.true_cell(i,j,k,iBlock)) CYCLE
             do iDiff = iDiffMin, iDiffMax
                iVar = iDiff_I(iDiff)
                pDotADotPPe = pDotADotPPe + 0.5 *(&
                     DiffCoef_VFDB(iDiff,i+1,j,k,1,iBlock)*   &
                     (   StateImpl_VG(iVar,i+1,j,k)   &
                     -   StateImpl_VG(iVar,i  ,j,k))**2  &
                     + DiffCoef_VFDB(iDiff,i  ,j,k,1,iBlock)* &
                     (   StateImpl_VG(iVar,i  ,j,k)   &
                     -   StateImpl_VG(iVar,i-1,j,k))**2  &
                     + DiffCoef_VFDB(iDiff,i,j+1,k,2,iBlock)* &
                     (   StateImpl_VG(iVar,i,j+1,k)   &
                     -   StateImpl_VG(iVar,i,j  ,k))**2  &
                     + DiffCoef_VFDB(iDiff,i,j  ,k,2,iBlock)* &
                     (   StateImpl_VG(iVar,i,j  ,k)   &
                     -   StateImpl_VG(iVar,i,j-1,k))**2  &
                     + DiffCoef_VFDB(iDiff,i,j,k+1,3,iBlock)* &
                     (   StateImpl_VG(iVar,i,j,k+1)   &
                     -   StateImpl_VG(iVar,i,j,k  ))**2  &
                     + DiffCoef_VFDB(iDiff,i,j,k  ,3,iBlock)* &
                     (   StateImpl_VG(iVar,i,j,k  )   &
                     -   StateImpl_VG(iVar,i,j,k-1))**2 )
             end do
          end do; end do; end do

          !\
          ! Correct the contributions to the quadratic form
          ! from the boundary faces
          !/
          if(NeiLev(5,iBlock) == NOBLK)then
             k = 1
             do j = 1, nJ; do i = 1, nI; do iDiff = iDiffMin, iDiffMax
                iVar = iDiff_I(iDiff)
                pDotADotPPe = pDotADotPPe + 0.5 *&
                     DiffCoef_VFDB(iDiff,i,j,k  ,3,iBlock)* &
                     (   StateImpl_VG(iVar,i,j,k  )**2  &
                     -   StateImpl_VG(iVar,i,j,k-1)**2)
             end do; end do; end do
          end if

          if(NeiLev(6,iBlock) == NOBLK)then
             k = nK
             do j = 1, nJ; do i = 1, nI; do iDiff = iDiffMin, iDiffMax
                if(.not.true_cell(i,j,k,iBlock)) CYCLE
                iVar = iDiff_I(iDiff)
                pDotADotPPe = pDotADotPPe + 0.5 *&
                     DiffCoef_VFDB(iDiff,i,j,k+1,3,iBlock)* &
                     (  - StateImpl_VG(iVar,i,j,k+1)**2   &
                     +   StateImpl_VG(iVar,i,j,k  )**2  )
             end do; end do; end do  
          end if
       end if

       if(NeiLev(1,iBlock) == NOBLK)then
          i = 1
          do k = 1, nK; do j = 1, nJ; do iDiff = iDiffMin, iDiffMax
             if(.not.true_cell(i,j,k,iBlock)) CYCLE
             iVar = iDiff_I(iDiff)
             pDotADotPPe = pDotADotPPe + 0.5 *&
                  DiffCoef_VFDB(iDiff,i  ,j,k,1,iBlock)* &
                  (   StateImpl_VG(iVar,i  ,j,k)**2   &
                  -   StateImpl_VG(iVar,i-1,j,k)**2 )
          end do; end do; end do  
       end if

       if(NeiLev(2,iBlock) == NOBLK)then
          i = nI
          do k = 1, nK; do j = 1, nJ; do iDiff = iDiffMin, iDiffMax
             if(.not.true_cell(i,j,k,iBlock)) CYCLE
             iVar = iDiff_I(iDiff)
             pDotADotPPe = pDotADotPPe + 0.5 *&
                  DiffCoef_VFDB(iDiff,i+1,j,k,1,iBlock)*   &
                  ( - StateImpl_VG(iVar,i+1,j,k)**2   &
                  +   StateImpl_VG(iVar,i  ,j,k)**2 ) 
          end do; end do; end do  
       end if

       if(nDim > 1 .and. NeiLev(3,iBlock) == NOBLK)then
          j = 1
          do k = 1, nK; do i = 1, nI; do iDiff = iDiffMin, iDiffMax
             if(.not.true_cell(i,j,k,iBlock)) CYCLE
             iVar = iDiff_I(iDiff)
             pDotADotPPe = pDotADotPPe + 0.5 *&
                  DiffCoef_VFDB(iDiff,i,j  ,k,2,iBlock)* &
                  (   StateImpl_VG(iVar,i,j  ,k)**2   &
                  -   StateImpl_VG(iVar,i,j-1,k)**2 ) 
          end do; end do; end do  
       end if

       if(nDim > 1 .and. NeiLev(4,iBlock) == NOBLK)then
          j = nJ
          do k = 1, nK; do i = 1, nI; do iDiff = iDiffMin, iDiffMax
             if(.not.true_cell(i,j,k,iBlock)) CYCLE
             iVar = iDiff_I(iDiff)
             pDotADotPPe = pDotADotPPe + 0.5 *&
                  DiffCoef_VFDB(iDiff,i,j+1,k,2,iBlock)* &
                  ( - StateImpl_VG(iVar,i,j+1,k)**2   &
                  +   StateImpl_VG(iVar,i,j  ,k)**2 )
          end do; end do; end do 
       end if

       ! Source terms due to energy exchange
       if(nRelax > 0)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.true_cell(i,j,k,iBlock)) CYCLE
             do iRelax = iRelaxMin, iRelaxMax
                iVar = iRelax_I(iRelax)

                pDotADotPPe = pDotADotPPe + &
                     RelaxCoef_VCB(iRelax,i,j,k,iBlock) &
                     *CellVolume_GB(i,j,k,iBlock) &
                     *(StateImpl_VG(iTeImpl,i,j,k)-StateImpl_VG(iVar,i,j,k))**2
             end do
          end do; end do; end do
       end if

       ! Point implicit source terms due to energy exchange
       if(nPoint > 0 .and. IsLinear)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.true_cell(i,j,k,iBlock)) CYCLE
             pDotADotPPe = pDotADotPPe + &
                  PointCoef_VCB(iPointSemi,i,j,k,iBlock) &
                  *StateImpl_VG(1,i,j,k)**2 &
                  *CellVolume_GB(i,j,k,iBlock)
          end do; end do; end do
       end if

    else
       pDotADotPPe = 0.0
    end if

  contains

    subroutine correct_left_ghostcell(iDim,iMin,iMax,jMin,jMax,kMin,kMax)

      use ModNumConst, ONLY: i_DD

      integer, intent(in) :: iDim, iMin, iMax, jMin, jMax, kMin, kMax

      integer :: i, j, k, iShift, jShift, kShift, Di, Dj, Dk, i1, j1, k1
      integer :: iDiff, iVar
      !------------------------------------------------------------------------

      Di = i_DD(iDim,1); Dj = i_DD(iDim,2); Dk = i_DD(iDim,3)
      iShift = 1-Di; jShift = min(1-Dj,nJ-1); kShift = min(1-Dk,nK-1)
      do k=kMin,kMax,2-Dk; do j=jMin,jMax,2-Dj; do i=iMin,iMax,2-Di
         i1=i+Di; j1=j+Dj; k1=k+Dk
         do iDiff = iDiffMin, iDiffMax
            iVar = iDiff_I(iDiff)
            StateImpl_VG(iVar,i:i+iShift,j:j+jShift,k:k+kShift) = &
                 StateImpl_VG(iVar,i:i+iShift,j:j+jShift,k:k+kShift) &
                 + StateImpl_VG(iVar,i1:i1+iShift,j1:j1+jShift,k1:k1+kShift) &
                 -(2.0/2**nDim)*sum(StateImpl_VG( &
                 iVar,i1:i1+iShift,j1:j1+jShift,k1:k1+kShift))
         end do
      enddo; enddo; enddo

    end subroutine correct_left_ghostcell

    !==========================================================================

    subroutine correct_right_ghostcell(iDim,iMin,iMax,jMin,jMax,kMin,kMax)

      use ModNumConst, ONLY: i_DD

      integer, intent(in) :: iDim, iMin, iMax, jMin, jMax, kMin, kMax

      integer :: i, j, k, iShift, jShift, kShift, Di, Dj, Dk, i1, j1, k1
      integer :: iDiff, iVar
      !------------------------------------------------------------------------

      Di = i_DD(iDim,1); Dj = i_DD(iDim,2); Dk = i_DD(iDim,3)
      iShift = 1-Di; jShift = min(1-Dj,nJ-1); kShift = min(1-Dk,nK-1)
      do k=kMin,kMax,2-Dk; do j=jMin,jMax,2-Dj; do i=iMin,iMax,2-Di
         i1=i-Di; j1=j-Dj; k1=k-Dk
         do iDiff = iDiffMin, iDiffMax
            iVar = iDiff_I(iDiff)
            StateImpl_VG(iVar,i:i+iShift,j:j+jShift,k:k+kShift) = &
                 StateImpl_VG(iVar,i:i+iShift,j:j+jShift,k:k+kShift) &
                 + StateImpl_VG(iVar,i1:i1+iShift,j1:j1+jShift,k1:k1+kShift) &
                 -(2.0/2**nDim)*sum(StateImpl_VG( &
                 iVar,i1:i1+iShift,j1:j1+jShift,k1:k1+kShift))
         end do
      enddo; enddo; enddo

    end subroutine correct_right_ghostcell

  end subroutine get_rad_diffusion_rhs

  !============================================================================

  subroutine add_jacobian_rad_diff(iBlock, Jacobian_VVCI)

    use BATL_lib,    ONLY: IsCartesian, DiLevelNei_IIIB, Unset_, CellSize_DB, &
         CellFace_DB, CellFace_DFB, CellVolume_B, CellVolume_GB
    use ModImplicit, ONLY: iTeImpl, UseFullImplicit, &
         UseSemiImplicit, nVarSemi, nStencil, Stencil1_, Stencil2_, &
         Stencil3_, Stencil4_, Stencil5_, Stencil6_, Stencil7_, UseNoOverlap
    use ModMain,     ONLY: nI, nJ, nK, TypeBc_I
    use ModNumConst, ONLY: i_DD
    use ModPhysics,  ONLY: InvClight
    use ModGeometry, ONLY: true_cell

    integer, intent(in) :: iBlock
    real, intent(inout) :: Jacobian_VVCI(nVarSemi,nVarSemi,nI,nJ,nK,nStencil)

    integer :: iVar, i, j, k, iDim, Di, Dj, Dk, iDiff, iRelax
    real :: DiffLeft, DiffRight, RelaxCoef, PlanckWeight
    real :: InvDcoord_D(MaxDim), CoeffLeft, CoeffRight
    real :: Coeff0, Coeff
    !--------------------------------------------------------------------------

    if(UseSemiImplicit)then
       ! Point implicit for ions (or electrons with no heat-conduction)
       if(nPoint > 0)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.true_cell(i,j,k,iBlock)) CYCLE
             ! dSvar/dVar (diagonal)
             Jacobian_VVCI(1,1,i,j,k,1) = Jacobian_VVCI(1,1,i,j,k,1) &
                  - PointCoef_VCB(iPointSemi,i,j,k,iBlock)
          end do; end do; end do
       end if

       if(nRelax > 0)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.true_cell(i,j,k,iBlock)) CYCLE

             do iRelax = iRelaxMin, iRelaxMax
                iVar = iRelax_I(iRelax)

                RelaxCoef = RelaxCoef_VCB(iRelax,i,j,k,iBlock)
                PlanckWeight = PlanckWeight_WCB(iRelax,i,j,k,iBlock)

                ! dSvar/dVar (diagonal)
                Jacobian_VVCI(iVar,iVar,i,j,k,1) = &
                     Jacobian_VVCI(iVar,iVar,i,j,k,1) - RelaxCoef

                ! dSe/dVar (off diagonal)
                Jacobian_VVCI(iTeImpl,iVar,i,j,k,1) = &
                     Jacobian_VVCI(iTeImpl,iVar,i,j,k,1) + RelaxCoef

                ! dSe/daTe^4 (diagonal)
                Jacobian_VVCI(iTeImpl,iTeImpl,i,j,k,1) = &
                     Jacobian_VVCI(iTeImpl,iTeImpl,i,j,k,1) &
                     - RelaxCoef*PlanckWeight

                ! dSvar/daTe^4 (off diagonal)
                Jacobian_VVCI(iVar,iTeImpl,i,j,k,1) = &
                     Jacobian_VVCI(iVar,iTeImpl,i,j,k,1) &
                     + RelaxCoef*PlanckWeight
             end do
          end do; end do; end do
       end if
    end if

    ! For the fully implicit scheme:
    ! add partial derivatives of the rad diffusion term to the Jacobian that
    ! are not calculated by the general algorithm, these are for the diffusion
    ! operators the same as the semi-implicit jacobian.

    if(UseAccurateRadiation)then

       InvDcoord_D = 1/CellSize_DB(:,iBlock)
       do iDim = 1, nDim
          if(IsCartesian)then
             CoeffLeft  = CellFace_DB(iDim,iBlock) &
                  *InvDCoord_D(iDim)/CellVolume_B(iBlock)
             CoeffRight = CellFace_DB(iDim,iBlock) &
                  *InvDCoord_D(iDim)/CellVolume_B(iBlock)
          end if
          Di = i_DD(iDim,1); Dj = i_DD(iDim,2); Dk = i_DD(iDim,3)
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.true_cell(i,j,k,iBlock)) CYCLE

             if(.not.IsCartesian)then ! rz-geometry only
                CoeffLeft  = CellFace_DFB(iDim,i,j,k,iBlock) &
                     *InvDCoord_D(iDim)/CellVolume_GB(i,j,k,iBlock)
                CoeffRight = CellFace_DFB(iDim,i+Di,j+Dj,k+Dk,iBlock) &
                     *InvDCoord_D(iDim)/CellVolume_GB(i,j,k,iBlock)
             end if
             do iDiff = iDiffMin, iDiffMax
                iVar = iDiff_I(iDiff)

                DiffLeft = CoeffLeft &
                     *DiffCoef_VFDB(iDiff,i,j,k,iDim,iBlock)
                DiffRight = CoeffRight &
                     *DiffCoef_VFDB(iDiff,i+Di,j+Dj,k+Dk,iDim,iBlock)

                Jacobian_VVCI(iVar,iVar,i,j,k,1) = &
                     Jacobian_VVCI(iVar,iVar,i,j,k,1) - (DiffLeft + DiffRight)

                if(UseNoOverlap)then
                   if(  iDim==1.and.i==1  .or. &
                        iDim==2.and.j==1  .or. &
                        iDim==3.and.k==1 )       DiffLeft = 0.0
                   if(  iDim==1.and.i==nI .or. &
                        iDim==2.and.j==nJ .or. &
                        iDim==3.and.k==nK)       DiffRight = 0.0
                end if

                Jacobian_VVCI(iVar,iVar,i,j,k,2*iDim) = &
                     Jacobian_VVCI(iVar,iVar,i,j,k,2*iDim) + DiffLeft
                Jacobian_VVCI(iVar,iVar,i,j,k,2*iDim+1) = &
                     Jacobian_VVCI(iVar,iVar,i,j,k,2*iDim+1) + DiffRight
             end do
          end do; end do; end do
       end do

    else

       Coeff = 1.0
       do iDim = 1, nDim
          if(UseFullImplicit) &
               Coeff = CellFace_DB(iDim,iBlock)/CellSize_DB(iDim,iBlock)

          Di = i_DD(iDim,1); Dj = i_DD(iDim,2); Dk = i_DD(iDim,3)
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.true_cell(i,j,k,iBlock)) CYCLE
             do iDiff = iDiffMin, iDiffMax
                iVar = iDiff_I(iDiff)

                DiffLeft = Coeff/CellVolume_GB(i,j,k,iBlock) &
                     *DiffCoef_VFDB(iDiff,i,j,k,iDim,iBlock)
                DiffRight = Coeff/CellVolume_GB(i,j,k,iBlock) &
                     *DiffCoef_VFDB(iDiff,i+Di,j+Dj,k+Dk,iDim,iBlock)

                Jacobian_VVCI(iVar,iVar,i,j,k,1) = &
                     Jacobian_VVCI(iVar,iVar,i,j,k,1) - (DiffLeft + DiffRight)

                if(UseNoOverlap)then
                   if(  iDim==1.and.i==1  .or. &
                        iDim==2.and.j==1  .or. &
                        iDim==3.and.k==1 )       DiffLeft = 0.0
                   if(  iDim==1.and.i==nI .or. &
                        iDim==2.and.j==nJ .or. &
                        iDim==3.and.k==nK)       DiffRight = 0.0
                end if

                Jacobian_VVCI(iVar,iVar,i,j,k,2*iDim) = &
                     Jacobian_VVCI(iVar,iVar,i,j,k,2*iDim) + DiffLeft
                Jacobian_VVCI(iVar,iVar,i,j,k,2*iDim+1) = &
                     Jacobian_VVCI(iVar,iVar,i,j,k,2*iDim+1) + DiffRight
             end do
          end do; end do; end do
       end do

    end if

    ! For now the boundary conditions are only applied when
    ! the block edges are not neglected (e.g. for HYPRE preconditioner)
    if(UseNoOverLap) RETURN

    ! Apply boundary conditions
    if(DiLevelNei_IIIB(-1,0,0,iBlock) == Unset_)then
       select case(TypeBc_I(1))
       case('float','outflow')
          Coeff0 = 2*InvClight/CellSize_DB(1,iBlock)
          do iDiff = iDiffRadMin, iDiffRadMax
             iVar = iDiff_I(iDiff)
             do k = 1, nK; do j = 1, nJ
                if(.not.true_cell(1,j,k,iBlock)) CYCLE
                ! Taken from ModRadDiffusion::set_rad_outflow_bc
                Coeff = Coeff0*DiffSemiCoef_VGB(iDiff,1,j,k,iBlock)
                Jacobian_VVCI(iVar,iVar,1,j,k,Stencil1_) = &
                     Jacobian_VVCI(iVar,iVar,1,j,k,Stencil1_) &
                     + (Coeff - 0.5)/(Coeff + 0.5) &
                     *Jacobian_VVCI(iVar,iVar,1,j,k,Stencil2_)
             end do; end do
          end do
       end select
    end if

    if(DiLevelNei_IIIB(+1,0,0,iBlock) == Unset_)then
       select case(TypeBc_I(2))
       case('float','outflow')
          Coeff0 = 2*InvClight/CellSize_DB(1,iBlock)
          do iDiff = iDiffRadMin, iDiffRadMax
             iVar = iDiff_I(iDiff)
             do k = 1, nK; do j = 1, nJ
                if(.not.true_cell(nI,j,k,iBlock)) CYCLE
                ! Taken from ModRadDiffusion::set_rad_outflow_bc
                Coeff = Coeff0*DiffSemiCoef_VGB(iDiffRadMin,nI,j,k,iBlock)

                Jacobian_VVCI(iVar,iVar,nI,j,k,Stencil1_) = &
                     Jacobian_VVCI(iVar,iVar,nI,j,k,Stencil1_) &
                     + (Coeff - 0.5)/(Coeff + 0.5) &
                     *Jacobian_VVCI(iVar,iVar,nI,j,k,Stencil3_)
             end do; end do
          end do
       end select
    end if

    if(nJ > 1)then
       if(DiLevelNei_IIIB(0,-1,0,iBlock) == Unset_)then
          select case(TypeBc_I(3))
          case('float','outflow')
             Coeff0 = 2*InvClight/CellSize_DB(2,iBlock)
             do iDiff = iDiffRadMin, iDiffRadMax
                iVar = iDiff_I(iDiff)
                do k = 1, nK; do i = 1, nI
                   if(.not.true_cell(i,1,k,iBlock)) CYCLE
                   Coeff = Coeff0*DiffSemiCoef_VGB(iDiffRadMin,i,1,k,iBlock)
                   Jacobian_VVCI(iVar,iVar,i,1,k,Stencil1_) = &
                        Jacobian_VVCI(iVar,iVar,i,1,k,Stencil1_) &
                        + (Coeff - 0.5)/(Coeff + 0.5) &
                        *Jacobian_VVCI(iVar,iVar,i,1,k,Stencil4_)
                end do; end do
             end do
          end select
       end if

       if(DiLevelNei_IIIB(0,+1,0,iBlock) == Unset_)then
          select case(TypeBc_I(4))
          case('float','outflow')
             Coeff0 = 2*InvClight/CellSize_DB(2,iBlock)
             do iDiff = iDiffRadMin, iDiffRadMax
                iVar = iDiff_I(iDiff)
                do k = 1, nK; do i = 1, nI
                   if(.not.true_cell(i,nJ,k,iBlock)) CYCLE
                   Coeff = Coeff0*DiffSemiCoef_VGB(iDiffRadMin,i,nJ,k,iBlock)
                   Jacobian_VVCI(iVar,iVar,i,nJ,k,Stencil1_) = &
                        Jacobian_VVCI(iVar,iVar,i,nJ,k,Stencil1_) &
                        +(Coeff - 0.5)/(Coeff + 0.5) &
                        *Jacobian_VVCI(iVar,iVar,i,nJ,k,Stencil5_)
                end do; end do
             end do
          end select
       end if
    end if

    if(nK > 1)then
       if(DiLevelNei_IIIB(0,0,-1,iBlock) == Unset_)then
          select case(TypeBc_I(5))
          case('float','outflow')
             Coeff0 = 2*InvClight/CellSize_DB(3,iBlock)
             do iDiff = iDiffRadMin, iDiffRadMax
                iVar = iDiff_I(iDiff)
                do j = 1, nJ; do i = 1, nI
                   if(.not.true_cell(i,j,1,iBlock)) CYCLE
                   Coeff = Coeff0*DiffSemiCoef_VGB(iDiffRadMin,i,j,1,iBlock)
                   Jacobian_VVCI(iVar,iVar,i,j,1,Stencil1_) = &
                        Jacobian_VVCI(iVar,iVar,i,j,1,Stencil1_) &
                        + (Coeff - 0.5)/(Coeff + 0.5) &
                        *Jacobian_VVCI(iVar,iVar,i,j,1,Stencil6_)
                end do; end do
             end do
          end select
       end if

       if(DiLevelNei_IIIB(0,0,+1,iBlock) == Unset_)then
          select case(TypeBc_I(6))
          case('float','outflow')
             Coeff0 = 2*InvClight/CellSize_DB(3,iBlock)
             do iDiff = iDiffRadMin, iDiffRadMax
                iVar = iDiff_I(iDiff)
                do j = 1, nJ; do i = 1, nI
                   if(.not.true_cell(i,j,nK,iBlock)) CYCLE
                   Coeff = Coeff0*DiffSemiCoef_VGB(iDiffRadMin,i,j,nK,iBlock)
                   Jacobian_VVCI(iVar,iVar,i,j,nK,Stencil1_) = &
                        Jacobian_VVCI(iVar,iVar,i,j,nK,Stencil1_) &
                        +(Coeff - 0.5)/(Coeff + 0.5) &
                        *Jacobian_VVCI(iVar,iVar,i,j,nK,Stencil7_)
                end do; end do
             end do
          end select
       end if
    end if

  end subroutine add_jacobian_rad_diff

  !============================================================================

  subroutine update_impl_rad_diff(iBlock, iImplBlock, StateImpl_VG)

    ! The use ModVarIndexes has to be next to use ModAdvance for sake
    ! of the extremely advanced PGF90 12.9 compiler
    use ModAdvance,    ONLY: State_VGB, UseElectronPressure
    use ModVarIndexes, ONLY: Rho_, p_, ExtraEint_, Pe_, nWave, WaveFirst_
    use ModEnergy,     ONLY: calc_energy_cell
    use ModImplicit,   ONLY: nw, iTeImpl, iTrImplFirst, &
         DconsDsemi_VCB, ImplOld_VCB, ImplCoeff
    use ModMain,       ONLY: nI, nJ, nK, Dt, UseRadDiffusion
    use ModPhysics,    ONLY: inv_gm1, gm1, No2Si_V, Si2No_V, UnitEnergyDens_, &
         UnitP_, UnitRho_, UnitTemperature_, ExtraEintMin
    use ModUser,       ONLY: user_material_properties
    use ModGeometry,   ONLY: true_cell

    integer, intent(in) :: iBlock, iImplBlock
    real, intent(inout) :: StateImpl_VG(nw,nI,nJ,nK)

    integer :: i, j, k, iVarImpl, iVar, iWave
    real :: Einternal, EinternalSi, PressureSi
    real :: PeSi, Ee, EeSi
    real :: Relaxation

    character(len=*), parameter :: NameSub = 'update_impl_rad_diff'
    !--------------------------------------------------------------------------

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(.not.true_cell(i,j,k,iBlock)) CYCLE
       if(UseRadDiffusion)then
          do iWave = 1, nWave
             iVarImpl = iTrImplFirst - 1 + iWave
             iVar = WaveFirst_ - 1 + iWave

             StateImpl_VG(iVarImpl,i,j,k) = &
                  max(EradMin_W(iWave), StateImpl_VG(iVarImpl,i,j,k))

             State_VGB(iVar,i,j,k,iBlock) = StateImpl_VG(iVarImpl,i,j,k)
          end do
       end if

       if(UseElectronPressure)then
          ! electron pressure -> electron internal energy Ee
          Ee = inv_gm1*State_VGB(Pe_,i,j,k,iBlock) &
               + State_VGB(ExtraEint_,i,j,k,iBlock)

          ! electron energy update: Ee_new = Ee_old + Cv'*Delta(a*Te^4)
          Ee = Ee + DconsDsemi_VCB(iTeImpl,i,j,k,iImplBlock) &
               *(StateImpl_VG(iTeImpl,i,j,k)-ImplOld_VCB(iTeImpl,i,j,k,iBlock))

          ! ion pressure -> Einternal
          Einternal = inv_gm1*State_VGB(p_,i,j,k,iBlock)
       else
          ! ion + electron pressure -> Einternal
          Einternal = inv_gm1*State_VGB(p_,i,j,k,iBlock) &
               + State_VGB(ExtraEint_,i,j,k,iBlock)

          if(iTeImpl>0)then
             Einternal = Einternal + DconsDsemi_VCB(iTeImpl,i,j,k,iImplBlock) &
                  *(StateImpl_VG(iTeImpl,i,j,k) &
                  - ImplOld_VCB(iTeImpl,i,j,k,iBlock))
          end if
       end if

       ! Add up internal energy change from all point-implicit variables
       ! Multiple point-implicit variables can occur 
       ! in split semi-implicit scheme only.
       do iVar = 1, nPoint

          Relaxation = PointCoef_VCB(iVar,i,j,k,iBlock)*( &
               + ImplCoeff*(StateImpl_VG(iVar,i,j,k) &
               -            PointImpl_VCB(iVar,i,j,k,iBlock)) &
               + (1.0 - ImplCoeff)*(ImplOld_VCB(iVar,i,j,k,iBlock) &
               -            PointImpl_VCB(iVar,i,j,k,iBlock)) )

          if(UseElectronPressure .and. iVar > 1)then
             ! Add energy exchange between electrons and each radiation group
             ! for split semi-implicit scheme
             Ee = Ee + Dt*Relaxation
          else
             ! Add energy exchange between ions and electrons
             ! or ions+electrons and radiation (when UseElectronPressure=F)
             Einternal = Einternal + Dt*Relaxation
          end if
       end do

       if(Einternal < 0.0 .and. .not.UseElectronPressure)then
          write(*,*)NameSub,': ERROR Rho, p, TOrigSi=', &
               State_VGB(Rho_,i,j,k,iBlock)*No2Si_V(UnitRho_), &
               State_VGB(p_,i,j,k,iBlock)*No2Si_V(UnitP_), &
               State_VGB(p_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock) &
               *No2Si_V(UnitTemperature_)

          write(*,*)NameSub,': ERROR negative Eint=', Einternal
          write(*,*)NameSub,': ERROR at i,j,k,iBlock=', i, j, k, iBlock
          call stop_mpi(NameSub//' negative Eint')
       end if

       if(UseElectronPressure)then
          ! ions
          State_VGB(p_,i,j,k,iBlock) = max(1e-30, gm1*Einternal)

          ! electrons
          EeSi = Ee*No2Si_V(UnitEnergyDens_)

          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               i, j, k, iBlock, &
               EinternalIn = EeSi, PressureOut = PeSi)

          ! Set true electron pressure
          State_VGB(Pe_,i,j,k,iBlock) = PeSi*Si2No_V(UnitP_)

          ! Set ExtraEint = Electron internal energy - Pe/(gamma -1)
          State_VGB(ExtraEint_,i,j,k,iBlock) = max(ExtraEintMin, &
               Ee - inv_gm1*State_VGB(Pe_,i,j,k,iBlock))

       else
          ! ions + electrons
          EinternalSi = Einternal*No2Si_V(UnitEnergyDens_)

          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               i, j, k, iBlock, &
               EinternalIn = EinternalSi, PressureOut = PressureSi)

          ! Set true total pressure
          State_VGB(p_,i,j,k,iBlock) = PressureSi*Si2No_V(UnitP_)

          ! Set ExtraEint = Total internal energy - Ptotal/(gamma -1)
          State_VGB(ExtraEint_,i,j,k,iBlock) = max(ExtraEintMin, &
               Einternal - inv_gm1*State_VGB(p_,i,j,k,iBlock))

       end if

    end do; end do; end do

    call calc_energy_cell(iBlock)

  end subroutine update_impl_rad_diff
  !==========================================================================
  subroutine set_rad_diff_range

    use ModImplicit, ONLY: iVarSemi

    character(len=*), parameter:: NameSub = 'set_rad_diff_range'
    !-----------------------------------------------------------------------

    ! Range of relaxation (all semi-implicit variables should have relaxation)
    iRelaxMin   = iVarSemi
    iRelaxMax   = iVarSemi

    ! Range of diffusions (all semi-implicit variables should have diffusion)
    iDiffMin    = iVarSemi
    iDiffMax    = iVarSemi

    ! Range of radiation diffusion (never used for electron temperature)
    iDiffRadMin = iVarSemi
    iDiffRadMax = iVarSemi

    ! Index of point implicit energy exchange terms
    iPointSemi  = iVarSemi

  end subroutine set_rad_diff_range

end module ModRadDiffusion
