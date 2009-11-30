!^CFG COPYRIGHT UM
!^CFG FILE IMPLICIT
!============================================================================
module ModRadDiffusion

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

  ! Logical for adding radiation diffusion
  logical, public :: IsNewBlockRadDiffusion = .true.
  logical, public :: IsNewTimestepRadDiffusion = .true.

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
  real, allocatable :: PointImpl_CB(:,:,:,:)
  real, allocatable :: PlanckWeight_WCB(:,:,:,:,:)

  ! Index indicating which implicit variables involve diffusion
  integer, allocatable :: iDiff_I(:)
  integer :: nDiff
  ! which one are heat conduction, which one radiation diffusion
  integer :: iDiffHeat = 0, iDiffRadFirst = 0, iDiffRadLast = 0

  ! Index indicating which nRelax implicit variables involve energy exchange
  ! with the electrons
  integer, allocatable :: iRelax_I(:)
  integer :: nRelax

  ! Index which nPoint variables involve point implicit energy exchange
  ! with the state variable stored in PointImpl_CB
  integer, allocatable :: iPoint_I(:)
  integer :: nPoint

  ! radiation energy used for calculating radiative energy flux
  real, allocatable :: Erad_WG(:,:,:,:)
  ! temporary radiation energy array needed by set_block_field
  real, allocatable :: Erad1_WG(:,:,:,:)


  logical :: UseTemperature

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

    case default
       call stop_mpi(NameSub//' invalid NameCommand='//NameCommand)
    end select


  end subroutine read_rad_diffusion_param

  !============================================================================

  subroutine init_rad_diffusion

    use ModAdvance,     ONLY: nWave, UseElectronPressure
    use ModMain,        ONLY: UseRadDiffusion
    use ModSize,        ONLY: nI, nJ, nK, MaxBlock
    use ModImplicit,    ONLY: UseSemiImplicit, UseFullImplicit, &
         TypeSemiImplicit, iEradImpl, iTeImpl, iTrImplFirst, iTrImplLast
    use ModPhysics,     ONLY: Si2No_V, UnitTemperature_, UnitEnergyDens_, &
         cRadiationNo
    use ModVarIndexes,  ONLY: nVar, nWave, WaveFirst_
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
       allocate(Erad_WG(1,-1:nI+2,-1:nJ+2,-1:nK+2))

       nDiff = 1
       allocate(iDiff_I(nDiff))
       iDiff_I(1) = WaveFirst_
       nRelax = 1
    end if
       
    if(UseSemiImplicit)then

       allocate(Erad_WG(nWave,-1:nI+2,-1:nJ+2,-1:nK+2))
       allocate(Erad1_WG(nWave,0:nI+1,0:nJ+1,0:nK+1))

       ! Default to zero, unless reset
       iTeImpl = 0; iTrImplFirst = 0; iTrImplLast = 0;

       UseTemperature = .false.

       select case(TypeSemiImplicit)
       case('radiation')
          if(UseElectronPressure)then
             call stop_mpi(NameSub//": Te/=Ti requires Te,Ti relaxation")
          end if
          if(nWave == 1)then
             iTrImplFirst = 1; iTrImplLast = 1
             nRelax = 0
             nPoint = 1
             allocate(iPoint_I(nPoint))
             iPoint_I = iTrImplFirst
          else
             iTeImpl = 1; iTrImplFirst = 2; iTrImplLast = nWave + 1
             nPoint = 0
             nRelax = nWave
             allocate(iRelax_I(nRelax))
             iRelax_I = (/ (iVarImpl,iVarImpl=iTrImplFirst,iTrImplLast) /)
          end if
          nDiff = nWave
          allocate(iDiff_I(nDiff))
          iDiff_I = (/ (iVarImpl,iVarImpl=iTrImplFirst,iTrImplLast) /)
          iDiffHeat = 0; iDiffRadFirst = 1; iDiffRadLast = nWave

       case('radcond')
          iTeImpl = 1; iTrImplFirst = 2; iTrImplLast = iTrImplFirst + nWave - 1
          nDiff = 1 + nWave
          allocate(iDiff_I(nDiff))
          iDiff_I = (/ iTeImpl, (iVarImpl,iVarImpl=iTrImplFirst,iTrImplLast) /)
          iDiffHeat = 1; iDiffRadFirst = 2; iDiffRadLast = nWave + 1
          nRelax = nWave
          allocate(iRelax_I(nRelax))
          iRelax_I = (/ (iVarImpl,iVarImpl=iTrImplFirst,iTrImplLast) /)
          if(UseElectronPressure)then
             nPoint = 1
             allocate(iPoint_I(nPoint))
             iPoint_I = iTeImpl
          else
             nPoint = 0
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
             allocate(iPoint_I(nPoint))
             iPoint_I = iTeImpl
          else
             nPoint = 0
          end if

       end select

       iEradImpl = iTrImplFirst

       if(nPoint>0)then
          allocate(PointCoef_VCB(nPoint,nI,nJ,nK,MaxBlock))
          PointCoef_VCB = 0.0
          allocate(PointImpl_CB(nI,nJ,nK,MaxBlock))
          PointImpl_CB = 0.0
       end if

       allocate(DiffSemiCoef_VGB(nDiff,-1:nI+2,-1:nJ+2,-1:nK+2,MaxBlock))
       DiffSemiCoef_VGB = 0.0

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
         Grad2ByErad2 = sum(FaceGrad_D**2)/State_V(Erad_)**2

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
    use ModVarIndexes, ONLY: Energy_

    integer, intent(in) :: iBlock

    integer :: i, j, k
    real :: TeSi, Te
    real :: AbsorptionEmission, OpacityPlanckSi_W(nWave)
    character(len=*), parameter:: NameSub = "calc_source_rad_diffusion"
    !------------------------------------------------------------------------

    do k=1,nK; do j=1,nJ; do i=1,nI

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

       ! dE/dt    = - AbsorptionEmission
       Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) - AbsorptionEmission

    end do; end do; end do

  end subroutine calc_source_rad_diffusion

  !============================================================================
  ! Semi-implicit interface
  !============================================================================

  subroutine get_impl_rad_diff_state(StateImpl_VGB,DconsDsemi_VCB)

    use ModAdvance,  ONLY: State_VGB, UseElectronPressure, nWave, WaveFirst_, &
         WaveLast_
    use ModConst,    ONLY: cBoltzmann
    use ModImplicit, ONLY: nw, nImplBlk, impl2iBlk, TypeSemiImplicit, &
         iTeImpl, iTrImplFirst, iTrImplLast, ImplCoeff
    use ModMain,     ONLY: x_, y_, nI, nJ, nK, MaxImplBlk, Dt
    use ModNumConst, ONLY: i_DD
    use ModPhysics,  ONLY: inv_gm1, Clight, cRadiationNo, UnitN_, UnitP_, &
         Si2No_V, UnitTemperature_, UnitEnergyDens_, UnitX_, UnitU_, UnitT_, &
         No2Si_V
    use ModUser,     ONLY: user_material_properties
    use ModGeometry, ONLY: dx_BLK, dy_BLK, dz_BLK, vInv_CB, &
         UseCovariant, TypeGeometry, FaceAreaI_DFB, FaceAreaJ_DFB
    use ModParallel, ONLY: NOBLK, NeiLev
    use ModMessagePass, ONLY: message_pass_dir

    real, intent(out) :: StateImpl_VGB(nw,0:nI+1,0:nJ+1,0:nK+1,MaxImplBlk)
    real, intent(inout) :: DconsDsemi_VCB(nw,nI,nJ,nK,MaxImplBlk)

    integer :: iImplBlock, iBlock, i, j, k, iVar, iVarImpl
    real :: OpacityPlanckSi_W(nWave), OpacityRosselandSi_W(nWave)
    real :: OpacityPlanck_W(nWave), CvSi, Cv, TeSi, Te
    real :: HeatCondSi, HeatCond, TeTiRelaxSi, TeTiRelax
    real :: Grad2ByErad2, DiffRad, InvDx2, InvDy2, InvDz2
    real :: InvDx, InvDy
    real :: NatomicSi, CveSi, Cve, CviSi, Cvi, PiSi, TiSi, Ti, PeSi
    real :: TeTiCoefSi, TeTiCoef, TeTiCoefPrime

    integer :: iDim, Di, Dj, Dk, iDiff, nDimInUse
    real :: Coeff, Dxyz_D(MaxDim)

    real :: PlanckSi_W(nWave), Planck_W(nWave), Planck

    character(len=*), parameter:: NameSub='get_impl_rad_diff_state'
    !--------------------------------------------------------------------------

    do iImplBlock = 1, nImplBLK

       iBlock = impl2iBLK(iImplBlock)
       IsNewBlockRadDiffusion = .true.

       if(iTeImpl > 0)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  i, j, k, iBlock, TeOut = TeSi)
             Te = TeSi*Si2No_V(UnitTemperature_)
             if(UseTemperature)then
                StateImpl_VGB(iTeImpl,i,j,k,iImplBlock) = Te
             else
                StateImpl_VGB(iTeImpl,i,j,k,iImplBlock) = cRadiationNo*Te**4
             end if
          end do; end do; end do
       end if
       if(iTrImplFirst > 0)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             StateImpl_VGB(iTrImplFirst:iTrImplLast,i,j,k,iImplBlock) = &
                  State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock)
          end do; end do; end do
       end if

       InvDx2 = 0.5/dx_BLK(iBlock)
       InvDy2 = 0.5/dy_BLK(iBlock)
       InvDz2 = 0.5/dz_BLK(iBlock)

       ! calculate coefficients for linearized energy exchange and diffusion
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(UseElectronPressure)then
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  i, j, k, iBlock, &
                  OpacityPlanckOut_W = OpacityPlanckSi_W, &
                  OpacityRosselandOut_W = OpacityRosselandSi_W, &
                  CvOut=CveSi, TeOut = TeSi, NatomicOut=NatomicSi, &
                  HeatCondOut = HeatCondSi, TeTiRelaxOut = TeTiRelaxSi, &
                  PlanckOut_W = PlanckSi_W)

             PiSi = State_VGB(p_,i,j,k,iBlock)*No2Si_V(UnitP_)
             TiSi = PiSi/(cBoltzmann*NatomicSi)
             Ti = TiSi*Si2No_V(UnitTemperature_)
             CviSi = inv_gm1*cBoltzmann*NatomicSi
             Cvi = CviSi*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_)
             Te = TeSi*Si2No_V(UnitTemperature_)
             Cve = CveSi*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_)
             if(.not.UseTemperature)then
                Cvi = Cvi/(4.0*cRadiationNo*Ti**3)
                Cve = Cve/(4.0*cRadiationNo*Te**3)
             end if
             TeTiCoefSi = TeTiRelaxSi*cBoltzmann*NatomicSi
             TeTiCoef = TeTiCoefSi*Si2No_V(UnitEnergyDens_) &
                  /(Si2No_V(UnitTemperature_)*Si2No_V(UnitT_))
          else
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  i, j, k, iBlock, &
                  OpacityPlanckOut_W = OpacityPlanckSi_W, &
                  OpacityRosselandOut_W = OpacityRosselandSi_W, &
                  CvOut = CvSi, TeOut = TeSi, HeatCondOut=HeatCondSi, &
                  PlanckOut_W = PlanckSi_W)

             Te = TeSi*Si2No_V(UnitTemperature_)
             Cv = CvSi*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_)
             if(.not.UseTemperature) Cv = Cv/(4.0*cRadiationNo*Te**3)
          end if

          OpacityPlanck_W = OpacityPlanckSi_W/Si2No_V(UnitX_)
          if(nWave == 1)then
             Planck_W(1) = cRadiationNo*Te**4
          else
             Planck_W = PlanckSi_W*Si2No_V(UnitEnergyDens_)
          end if
          Planck = cRadiationNo*Te**4

          select case(TypeSemiImplicit)
          case('radiation')
             if(nWave == 1)then
                ! This coefficient is cR'' = cR/(1+dt*cR*dPlanck/dEint)
                PointCoef_VCB(1,i,j,k,iBlock) = Clight*OpacityPlanck_W(1) &
                     /(1 + ImplCoeff*Dt*Clight*OpacityPlanck_W(1) / Cv)

                ! This is just the Planck function at time level * saved
                PointImpl_CB(i,j,k,iBlock) = Planck
             else
                DconsDsemi_VCB(iTeImpl,i,j,k,iImplBlock) = Cv
                RelaxCoef_VCB(:,i,j,k,iBlock) = Clight*OpacityPlanck_W
                PlanckWeight_WCB(:,i,j,k,iBlock) = Planck_W/Planck
             end if

          case('radcond')
             if(UseElectronPressure)then
                TeTiCoefPrime = TeTiCoef/Cvi &
                     /(cRadiationNo*(Te+Ti)*(Te**2+Ti**2))
                PointCoef_VCB(1,i,j,k,iBlock) = &
                     TeTiCoefPrime*Cvi/(1.0 + ImplCoeff*Dt*TeTiCoefPrime)
                PointImpl_CB(i,j,k,iBlock) = cRadiationNo*Ti**4
                DconsDsemi_VCB(iTeImpl,i,j,k,iImplBlock) = Cve
             else
                DconsDsemi_VCB(iTeImpl,i,j,k,iImplBlock) = Cv
             end if

             RelaxCoef_VCB(:,i,j,k,iBlock) = Clight*OpacityPlanck_W
             PlanckWeight_WCB(:,i,j,k,iBlock) = Planck_W/Planck

          case('cond')
             if(UseElectronPressure)then
                DconsDsemi_VCB(iTeImpl,i,j,k,iImplBlock) = Cve

                PointCoef_VCB(1,i,j,k,iBlock) = TeTiCoef &
                     /(1.0 + ImplCoeff*Dt*TeTiCoef/Cvi)
                PointImpl_CB(i,j,k,iBlock) = Ti
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
       if(NeiLev(3,iBlock) == NOBLK)then
          j = 0
          do k = 1, nK; do i = 1, nI
             call get_ghostcell_diffcoef
          end do; end do
       end if
       if(NeiLev(4,iBlock) == NOBLK)then
          j = nJ + 1
          do k = 1, nK; do i = 1, nI
             call get_ghostcell_diffcoef
          end do; end do
       end if
       if(NeiLev(5,iBlock) == NOBLK)then
          k = 0
          do j = 1, nJ; do i = 1, nI
             call get_ghostcell_diffcoef
          end do; end do
       end if
       if(NeiLev(6,iBlock) == NOBLK)then
          k = nK + 1
          do j = 1, nJ; do i = 1, nI
             call get_ghostcell_diffcoef
          end do; end do
       end if

    end do

    nDimInUse = 3; if(TypeGeometry == 'rz') nDimInUse = 2

    ! Message pass to fill in ghost cells 
    call message_pass_dir(iDirMin=1,iDirMax=3,Width=1, &
         SendCorners=.false.,ProlongOrder=1,nVar=nDiff, &
         Sol_VGB=DiffSemiCoef_VGB,restrictface=.true.)

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

       if(nDimInUse==3)then
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

       if(.not.UseCovariant)then
          Dxyz_D = (/dx_BLK(iBlock), dy_BLK(iBlock), dz_Blk(iBlock)/)
          do iDim = 1, nDim
             ! FaceYZ/dx = Volume/dx^2
             Coeff = 1.0 / (Dxyz_D(iDim)**2 * vInv_CB(1,1,1,iBlock))
             Di = i_DD(iDim,1); Dj = i_DD(iDim,2); Dk = i_DD(iDim,3)
             do k=1,nK+Dk; do j=1,nJ+Dj; do i=1,nI+Di
                do iDiff = 1, nDiff
                   DiffCoef_VFDB(iDiff,i,j,k,iDim,iBlock) = &
                        Coeff*DiffCoef_VFDB(iDiff,i,j,k,iDim,iBlock)
                end do
             enddo; enddo; enddo
          end do

       elseif(TypeGeometry == 'rz')then

          InvDx = 1.0/Dx_Blk(iBlock)
          do k=1,nK; do j=1,nJ; do i=1,nI+1
             do iDiff = 1, nDiff
                DiffCoef_VFDB(iDiff,i,j,k,x_,iBlock) = &
                     InvDx*FaceAreaI_DFB(x_,i,j,k,iBlock) &
                     *DiffCoef_VFDB(iDiff,i,j,k,x_,iBlock)
             end do
          end do; end do; end do

          InvDy = 1.0/Dy_Blk(iBlock)
          do k=1,nK; do j=1,nJ+1; do i=1,nI
             do iDiff = 1, nDiff
                DiffCoef_VFDB(iDiff,i,j,k,y_,iBlock) = &
                     InvDy*FaceAreaJ_DFB(y_,i,j,k,iBlock) &
                     *DiffCoef_VFDB(iDiff,i,j,k,y_,iBlock)
             end do
          end do; end do; end do
       else
          call stop_mpi(NameSub//': unimplemented TypeGeometry=//TypeGeometry')
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
      iShift = 1-Di; jShift = 1-Dj; kShift = 1-Dk
      do k=kMin,kMax,2-Dk; do j=jMin,jMax,2-Dj; do i=iMin,iMax,2-Di
         do iDiff = 1, nDiff
            DiffCoef_VFDB(iDiff,i:i+iShift,j:j+jShift,k:k+kShift,iDim,iBlock) &
                 = (DiffSemiCoef_VGB(iDiff,i-Di,j-Dj,k-Dk,iBlock) &
                 + 0.5*sum(DiffSemiCoef_VGB( &
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
      iShift = 1-Di; jShift = 1-Dj; kShift = 1-Dk
      do k=kMin,kMax,2-Dk; do j=jMin,jMax,2-Dj; do i=iMin,iMax,2-Di
         i1=i-Di; j1=j-Dj; k1=k-Dk
         do iDiff = 1, nDiff
            DiffCoef_VFDB(iDiff,i:i+iShift,j:j+jShift,k:k+kShift,iDim,iBlock) &
                 = (DiffSemiCoef_VGB(iDiff,i,j,k,iBlock) &
                 + 0.5*sum(DiffSemiCoef_VGB( &
                 iDiff,i1:i1+iShift,j1:j1+jShift,k1:k1+kShift,iBlock)))*2.0/9.0
         end do
      enddo; enddo; enddo

    end subroutine face_right_fine2coarse

    !==========================================================================

    subroutine get_diffusion_coef

      use ModAdvance,      ONLY: nWave
      use ModFaceGradient, ONLY: set_block_field

      real :: OpacityRosseland_W(nWave), DiffRad_W(nWave)
      real :: Grad2ByErad2_W(nWave)
      !------------------------------------------------------------------------

      if(iDiffRadFirst > 0)then
         OpacityRosseland_W = OpacityRosselandSi_W/Si2No_V(UnitX_)

         ! Calculate the cell centered diffusion coefficients
         if(UseRadFluxLimiter)then

            if(IsNewBlockRadDiffusion)then
               Erad_WG = State_VGB(WaveFirst_:WaveLast_,:,:,:,iBlock)
               call set_block_field(iBlock, nWave, Erad1_WG, Erad_WG)

               IsNewBlockRadDiffusion = .false.
            end if

            Grad2ByErad2_W = &
                 ( ((Erad_WG(:,i+1,j,k) - Erad_WG(:,i-1,j,k))*InvDx2)**2 &
                 + ((Erad_WG(:,i,j+1,k) - Erad_WG(:,i,j-1,k))*InvDy2)**2 &
                 + ((Erad_WG(:,i,j,k+1) - Erad_WG(:,i,j,k-1))*InvDz2)**2 ) &
                 / Erad_WG(:,i,j,k)**2

            ! For now only Larsen's radiation flux limiter
            DiffRad_W = Clight/sqrt(9*OpacityRosseland_W**2 + Grad2ByErad2_W)
         else
            DiffRad_W = Clight/(3*OpacityRosseland_W)
         end if

         ! Store it for message passing
         DiffSemiCoef_VGB(iDiffRadFirst:iDiffRadLast,i,j,k,iBlock) = DiffRad_W
      end if

      if(iDiffHeat > 0)then
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
           HeatCondOut = HeatCondSi, TeOut = TeSi)

      call get_diffusion_coef

    end subroutine get_ghostcell_diffcoef

  end subroutine get_impl_rad_diff_state

  !============================================================================

  subroutine set_rad_outflow_bc(iSide, iBlock, iImplBlock, State_VG, IsLinear)

    use ModAdvance,  ONLY: State_VGB
    use ModImplicit, ONLY: nw
    use ModGeometry, ONLY: dx_BLK, dy_BLK, dz_BLK
    use ModMain,     ONLY: nI, nJ, nK
    use ModPhysics,  ONLY: Clight

    integer, intent(in) :: iSide, iBlock, iImplBlock
    real, intent(inout) :: State_VG(nw,-1:nI+2,-1:nJ+2,-1:nK+2)
    logical, intent(in) :: IsLinear

    integer :: iVar, i, j, k, iDiff
    real :: Coef
    character(len=*), parameter :: NameSub='set_rad_outflow_bc'
    !--------------------------------------------------------------------------

    select case(iSide)
    case(1)
       do k = 1, nK; do j = 1, nJ
          do iDiff = iDiffRadFirst, iDiffRadLast
             iVar = iDiff_I(iDiff)          
             Coef = 2/Clight &
                  *DiffSemiCoef_VGB(iDiff,1,j,k,iBlock)/dx_BLK(iBlock)
             State_VG(iVar,0,j,k) = State_VG(iVar,1,j,k) &
                  *(Coef - 0.5)/(Coef + 0.5)
          end do
       end do; end do
    case(2)
       do k = 1, nK; do j = 1, nJ
          do iDiff = iDiffRadFirst, iDiffRadLast
             iVar = iDiff_I(iDiff)
             Coef = 2/Clight &
                  *DiffSemiCoef_VGB(iDiff,nI,j,k,iBlock)/dx_BLK(iBlock)
             State_VG(iVar,nI+1,j,k) = State_VG(iVar,nI,j,k) &
                  *(Coef - 0.5)/(Coef + 0.5)
          end do
       end do; end do
    case(3)
       do k = 1, nK; do i = 1, nI
          do iDiff = iDiffRadFirst, iDiffRadLast
             iVar = iDiff_I(iDiff)
             Coef = 2/Clight &
                  *DiffSemiCoef_VGB(iDiff,i,1,k,iBlock)/dy_BLK(iBlock)
             State_VG(iVar,i,0,k) = State_VG(iVar,i,1,k) &
                  *(Coef - 0.5)/(Coef + 0.5)
          end do
       end do; end do
    case(4)
       do k = 1, nK; do i = 1, nI
          do iDiff = iDiffRadFirst, iDiffRadLast
             iVar = iDiff_I(iDiff)
             Coef = 2/Clight &
                  *DiffSemiCoef_VGB(iDiff,i,nJ,k,iBlock)/dy_BLK(iBlock)
             State_VG(iVar,i,nJ+1,k) = State_VG(iVar,i,nJ,k) &
                  *(Coef - 0.5)/(Coef + 0.5)
          end do
       end do; end do
    case(5)
       do j = 1, nJ; do i = 1, nI
          do iDiff = iDiffRadFirst, iDiffRadLast
             iVar = iDiff_I(iDiff)
             Coef = 2/Clight &
                  *DiffSemiCoef_VGB(iDiff,i,j,1,iBlock)/dz_BLK(iBlock)
             State_VG(iVar,i,j,0) = State_VG(iVar,i,j,1) &
                  *(Coef - 0.5)/(Coef + 0.5)
          end do
       end do; end do
    case(6)
       do k = j, nJ; do i = 1, nI
          do iDiff = iDiffRadFirst, iDiffRadLast
             iVar = iDiff_I(iDiff)
             Coef = 2/Clight &
                  *DiffSemiCoef_VGB(iDiff,i,j,nK,iBlock)/dz_BLK(iBlock)
             State_VG(iVar,i,j,nK+1) = State_VG(iVar,i,j,nK) &
                  *(Coef - 0.5)/(Coef + 0.5)
          end do
       end do; end do
    end select

  end subroutine set_rad_outflow_bc

  !============================================================================

  subroutine get_rad_diffusion_rhs(iBlock, StateImpl_VG, Rhs_VC, IsLinear)

    use ModGeometry, ONLY: TypeGeometry, vInv_CB
    use ModImplicit, ONLY: nw, iTeImpl
    use ModLinearSolver, ONLY: pDotADotPPe, UsePDotADotP
    use ModMain,     ONLY: nI, nJ, nK
    use ModParallel, ONLY: NeiLev, NOBLK

    integer, intent(in) :: iBlock
    real, intent(inout) :: StateImpl_VG(nw,-1:nI+2,-1:nJ+2,-1:nK+2)
    real, intent(out)   :: Rhs_VC(nw,nI,nJ,nK)
    logical, intent(in) :: IsLinear

    real :: Te, EnergyExchange
    integer :: i, j, k, iDiff, iRelax, iPoint, iVar

    character(len=*), parameter :: NameSub='get_rad_diffusion_rhs'
    !--------------------------------------------------------------------------

    if(NeiLev(1,iBlock)==1) call correct_left_ghostcell(1,0,0,1,nJ,1,nK)
    if(NeiLev(2,iBlock)==1) call correct_right_ghostcell(1,nI+1,nI+1,1,nJ,1,nK)
    if(NeiLev(3,iBlock)==1) call correct_left_ghostcell(2,1,nI,0,0,1,nK)
    if(NeiLev(4,iBlock)==1) call correct_right_ghostcell(2,1,nI,nJ+1,nJ+1,1,nK)

    if(TypeGeometry /= 'rz')then
       if(NeiLev(5,iBlock)==1) call correct_left_ghostcell(3,1,nI,1,nJ,0,0)
       if(NeiLev(6,iBlock)==1) &
            call correct_right_ghostcell(3,1,nI,1,nJ,nK+1,nK+1)
    end if

    Rhs_VC = 0.0

    if(TypeGeometry == 'rz')then
       ! No flux from Z direction
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          do iDiff = 1, nDiff
             iVar = iDiff_I(iDiff)
             Rhs_VC(iVar,i,j,k) = &
                  vInv_CB(i,j,k,iBlock) * ( &
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
                  -   StateImpl_VG(iVar,i,j-1,k)) )
          end do
       end do; end do; end do
    else
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          do iDiff = 1, nDiff
             iVar = iDiff_I(iDiff)
             Rhs_VC(iVar,i,j,k) = &
                  vInv_CB(i,j,k,iBlock) * ( &
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
                  -   StateImpl_VG(iVar,i,j,k-1)) )
          end do
       end do; end do; end do

    end if

    ! Source terms due to energy exchange
    if(nRelax > 0)then
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          do iRelax = 1, nRelax
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
             do iPoint = 1, nPoint
                iVar = iPoint_I(iPoint)

                Rhs_VC(iVar,i,j,k) = Rhs_VC(iVar,i,j,k) &
                     - PointCoef_VCB(iPoint,i,j,k,iBlock) &
                     *StateImpl_VG(iVar,i,j,k)
             end do
          end do; end do; end do
       else
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             do iPoint = 1, nPoint
                iVar = iPoint_I(iPoint)

                EnergyExchange = PointCoef_VCB(iPoint,i,j,k,iBlock) &
                     *(PointImpl_CB(i,j,k,iBlock) - StateImpl_VG(iVar,i,j,k))

                Rhs_VC(iVar,i,j,k) = Rhs_VC(iVar,i,j,k) + EnergyExchange
             end do
          end do; end do; end do
       end if
    end if

    if(UsePDotADotP)then
       if(TypeGeometry == 'rz')then
          ! No flux from Z direction
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             do iDiff = 1, nDiff
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
             do iDiff = 1, nDiff
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
             do j = 1, nJ; do i = 1, nI; do iDiff = 1, nDiff
                iVar = iDiff_I(iDiff)
                pDotADotPPe = pDotADotPPe + 0.5 *&
                     DiffCoef_VFDB(iDiff,i,j,k  ,3,iBlock)* &
                     (   StateImpl_VG(iVar,i,j,k  )**2  &
                     -   StateImpl_VG(iVar,i,j,k-1)**2)
             end do; end do; end do
          end if

          if(NeiLev(6,iBlock) == NOBLK)then
             k = nK
             do j = 1, nJ; do i = 1, nI; do iDiff = 1, nDiff
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
          do k = 1, nK; do j = 1, nJ; do iDiff = 1, nDiff
             iVar = iDiff_I(iDiff)
             pDotADotPPe = pDotADotPPe + 0.5 *&
                  DiffCoef_VFDB(iDiff,i  ,j,k,1,iBlock)* &
                  (   StateImpl_VG(iVar,i  ,j,k)**2   &
                  -   StateImpl_VG(iVar,i-1,j,k)**2 )
          end do; end do; end do  
       end if

       if(NeiLev(2,iBlock) == NOBLK)then
          i = nI
          do k = 1, nK; do j = 1, nJ; do iDiff = 1, nDiff
             iVar = iDiff_I(iDiff)
             pDotADotPPe = pDotADotPPe + 0.5 *&
                  DiffCoef_VFDB(iDiff,i+1,j,k,1,iBlock)*   &
                  ( - StateImpl_VG(iVar,i+1,j,k)**2   &
                  +   StateImpl_VG(iVar,i  ,j,k)**2 ) 
          end do; end do; end do  
       end if


       if(NeiLev(3,iBlock) == NOBLK)then
          j = 1
          do k = 1, nK; do i = 1, nI; do iDiff = 1, nDiff
             iVar = iDiff_I(iDiff)
             pDotADotPPe = pDotADotPPe + 0.5 *&
                  DiffCoef_VFDB(iDiff,i,j  ,k,2,iBlock)* &
                  (   StateImpl_VG(iVar,i,j  ,k)**2   &
                  -   StateImpl_VG(iVar,i,j-1,k)**2 ) 
          end do; end do; end do  
       end if

       if(NeiLev(4,iBlock) == NOBLK)then
          j = nJ
          do k = 1, nK; do i = 1, nI; do iDiff = 1, nDiff
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
             do iRelax = 1, nRelax
                iVar = iRelax_I(iRelax)

                pDotADotPPe = pDotADotPPe + &
                     RelaxCoef_VCB(iRelax,i,j,k,iBlock) &
                     /vInv_CB(i,j,k,iBlock)&
                     *(StateImpl_VG(iTeImpl,i,j,k)-StateImpl_VG(iVar,i,j,k))**2
             end do
          end do; end do; end do
       end if

       ! Point implicit source terms due to energy exchange
       if(nPoint > 0 .and. IsLinear)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             do iPoint = 1, nPoint
                iVar = iPoint_I(iPoint)
                pDotADotPPe = pDotADotPPe + &
                     PointCoef_VCB(iPoint,i,j,k,iBlock) &
                     *StateImpl_VG(iVar,i,j,k)**2 &
                     /vInv_CB(i,j,k,iBlock)
             end do
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
      iShift = 1-Di; jShift = 1-Dj; kShift = 1-Dk
      do k=kMin,kMax,2-Dk; do j=jMin,jMax,2-Dj; do i=iMin,iMax,2-Di
         i1=i+Di; j1=j+Dj; k1=k+Dk
         do iDiff = 1, nDiff
            iVar = iDiff_I(iDiff)
            StateImpl_VG(iVar,i:i+iShift,j:j+jShift,k:k+kShift) = &
                 StateImpl_VG(iVar,i:i+iShift,j:j+jShift,k:k+kShift) &
                 + StateImpl_VG(iVar,i1:i1+iShift,j1:j1+jShift,k1:k1+kShift) &
                 -0.25*sum(StateImpl_VG( &
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
      iShift = 1-Di; jShift = 1-Dj; kShift = 1-Dk
      do k=kMin,kMax,2-Dk; do j=jMin,jMax,2-Dj; do i=iMin,iMax,2-Di
         i1=i-Di; j1=j-Dj; k1=k-Dk
         do iDiff = 1, nDiff
            iVar = iDiff_I(iDiff)
            StateImpl_VG(iVar,i:i+iShift,j:j+jShift,k:k+kShift) = &
                 StateImpl_VG(iVar,i:i+iShift,j:j+jShift,k:k+kShift) &
                 + StateImpl_VG(iVar,i1:i1+iShift,j1:j1+jShift,k1:k1+kShift) &
                 -0.25*sum(StateImpl_VG( &
                 iVar,i1:i1+iShift,j1:j1+jShift,k1:k1+kShift))
         end do
      enddo; enddo; enddo

    end subroutine correct_right_ghostcell

  end subroutine get_rad_diffusion_rhs

  !============================================================================

  subroutine add_jacobian_rad_diff(iBlock, nVar, Jacobian_VVCI)

    use ModGeometry, ONLY: vInv_CB, dx_BLK, dy_BLK, dz_BLK, &
         fAx_BLK, fAy_BLK, fAz_BLK
    use ModImplicit, ONLY: TypeSemiImplicit, iTeImpl, UseFullImplicit, &
         UseSemiImplicit
    use ModMain,     ONLY: nI, nJ, nK
    use ModNumConst, ONLY: i_DD

    integer, parameter:: nStencil = 7 !!! 2*nDim + 1

    integer, intent(in) :: iBlock, nVar
    real, intent(inout) :: Jacobian_VVCI(nVar,nVar,nI,nJ,nK,nStencil)

    integer :: iVar, i, j, k, iDim, Di, Dj, Dk, iDiff, iRelax, iPoint
    real :: DiffLeft, DiffRight, RelaxCoef, PlanckWeight
    real :: Dxyz_D(MaxDim), Area_D(MaxDim), Coeff
    !--------------------------------------------------------------------------

    if(UseSemiImplicit)then
       if(nPoint > 0)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             do iPoint = 1, nPoint
                iVar = iPoint_I(iPoint)

                ! dSvar/dVar (diagonal)
                Jacobian_VVCI(iVar,iVar,i,j,k,1) = &
                     Jacobian_VVCI(iVar,iVar,i,j,k,1) &
                     - PointCoef_VCB(iPoint,i,j,k,iBlock)
             end do
          end do; end do; end do
       end if

       if(nRelax > 0)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             do iRelax = 1, nRelax
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

    Coeff = 1.0
    Dxyz_D = (/dx_BLK(iBlock), dy_BLK(iBlock), dz_Blk(iBlock)/)
    Area_D = (/fAx_BLK(iBlock), fAy_BLK(iBlock), fAz_BLK(iBlock)/)
    do iDim = 1, nDim
       if(UseFullImplicit) Coeff = -Area_D(iDim)/Dxyz_D(iDim)
       Di = i_DD(iDim,1); Dj = i_DD(iDim,2); Dk = i_DD(iDim,3)
       do k=1,nK; do j=1,nJ; do i=1,nI
          do iDiff = 1, nDiff
             iVar = iDiff_I(iDiff)
             DiffLeft = Coeff*vInv_CB(i,j,k,iBlock) &
                  *DiffCoef_VFDB(iDiff,i,j,k,iDim,iBlock)
             DiffRight = Coeff*vInv_CB(i,j,k,iBlock) &
                  *DiffCoef_VFDB(iDiff,i+Di,j+Dj,k+Dk,iDim,iBlock)
             Jacobian_VVCI(iVar,iVar,i,j,k,1) = &
                  Jacobian_VVCI(iVar,iVar,i,j,k,1) - (DiffLeft + DiffRight)

             if(iDim==1.and.i==1 .or. iDim==2.and.j==1 .or. iDim==3.and.k==1)&
                  DiffLeft = 0.0
             if(iDim==1.and.i==nI .or. iDim==2.and.j==nJ &
                  .or. iDim==3.and.k==nK) DiffRight = 0.0

             Jacobian_VVCI(iVar,iVar,i,j,k,2*iDim) = &
                  Jacobian_VVCI(iVar,iVar,i,j,k,2*iDim) + DiffLeft
             Jacobian_VVCI(iVar,iVar,i,j,k,2*iDim+1) = &
                  Jacobian_VVCI(iVar,iVar,i,j,k,2*iDim+1) + DiffRight
          end do
       end do; end do; end do
    end do

  end subroutine add_jacobian_rad_diff

  !============================================================================

  subroutine update_impl_rad_diff(iBlock, iImplBlock, StateImpl_VG)

    use ModAdvance,    ONLY: State_VGB, UseElectronPressure
    use ModEnergy,     ONLY: calc_energy_cell
    use ModImplicit,   ONLY: nw, iTeImpl, iTrImplFirst, iTrImplLast, &
         DconsDsemi_VCB, ImplOld_VCB, ImplCoeff
    use ModMain,       ONLY: nI, nJ, nK, Dt, UseRadDiffusion
    use ModPhysics,    ONLY: inv_gm1, g, No2Si_V, Si2No_V, UnitEnergyDens_, &
         UnitP_, UnitRho_, UnitTemperature_, PeMin
    use ModUser,       ONLY: user_material_properties
    use ModVarIndexes, ONLY: Rho_, p_, ExtraEint_, Pe_, nWave, WaveFirst_

    integer, intent(in) :: iBlock, iImplBlock
    real, intent(inout) :: StateImpl_VG(nw,nI,nJ,nK)

    integer :: i, j, k, iVarImpl, iVar, iPoint, iWave
    real :: Einternal, EinternalSi, PressureSi
    real :: PeSi, Ee, EeSi
    real :: Relaxation

    character(len=*), parameter :: NameSub = 'update_impl_rad_diff'
    !--------------------------------------------------------------------------

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
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
          ! electrons
          State_VGB(Pe_,i,j,k,iBlock) = State_VGB(Pe_,i,j,k,iBlock) &
               + (g - 1)*DconsDsemi_VCB(iTeImpl,i,j,k,iImplBlock) &
               *(StateImpl_VG(iTeImpl,i,j,k)-ImplOld_VCB(iTeImpl,i,j,k,iBlock))

          State_VGB(Pe_,i,j,k,iBlock) = max(State_VGB(Pe_,i,j,k,iBlock), PeMin)

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

       do iPoint = 1, nPoint
          iVar = iPoint_I(iPoint)

          Relaxation = PointCoef_VCB(iPoint,i,j,k,iBlock)*( &
               + ImplCoeff*(StateImpl_VG(iVar,i,j,k) &
               -            PointImpl_CB(i,j,k,iBlock)) &
               + (1.0 - ImplCoeff)*(ImplOld_VCB(iVar,i,j,k,iBlock) &
               -            PointImpl_CB(i,j,k,iBlock)) )

          Einternal = Einternal + Dt*Relaxation
       end do

       if(Einternal < 0.0)then
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
          State_VGB(p_,i,j,k,iBlock) = (g - 1)*Einternal

          ! electrons
          Ee = inv_gm1*State_VGB(Pe_,i,j,k,iBlock) &
               + State_VGB(ExtraEint_,i,j,k,iBlock)
          EeSi = Ee*No2Si_V(UnitEnergyDens_)

          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               i, j, k, iBlock, &
               EinternalIn = EeSi, PressureOut = PeSi)

          ! Set true electron pressure
          State_VGB(Pe_,i,j,k,iBlock) = PeSi*Si2No_V(UnitP_)

          ! Set ExtraEint = electron internal energy - Pe/(gamma -1)
          State_VGB(ExtraEint_,i,j,k,iBlock) = &
               Ee - inv_gm1*State_VGB(Pe_,i,j,k,iBlock)

       else
          ! ions + electrons
          EinternalSi = Einternal*No2Si_V(UnitEnergyDens_)

          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               i, j, k, iBlock, &
               EinternalIn = EinternalSi, PressureOut = PressureSi)

          State_VGB(p_,i,j,k,iBlock) = PressureSi*Si2No_V(UnitP_)

          State_VGB(ExtraEint_,i,j,k,iBlock) = &
               Einternal - inv_gm1*State_VGB(p_,i,j,k,iBlock)

       end if

    end do; end do; end do

    call calc_energy_cell(iBlock)

  end subroutine update_impl_rad_diff

end module ModRadDiffusion
