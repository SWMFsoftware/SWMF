!^CFG COPYRIGHT UM
!^CFG FILE IMPLICIT
!============================================================================
module ModGrayDiffusion

  use ModVarIndexes, ONLY: p_
  use ModProcMH, ONLY: iProc
  !use ModMain, ONLY: iTest, jTest, kTest, BlkTest

  implicit none
  save

  private !except

  ! Public methods
  public :: init_gray_diffusion
  public :: get_radiation_energy_flux
  public :: calc_source_gray_diffusion
  public :: set_gray_outflow_bc
  public :: get_impl_gray_diff_state
  public :: get_gray_diffusion_rhs
  public :: add_jacobian_gray_diff
  public :: update_impl_gray_diff

  ! Logical for adding Gray Diffusion
  logical, public :: IsNewBlockGrayDiffusion = .true.
  logical, public :: IsNewTimestepGrayDiffusion = .true.

  ! Coefficients for two-temperature electron-radiation model
  real, allocatable :: DiffCoef_VFDB(:,:,:,:,:,:)
  real, allocatable :: RelaxCoef_VCB(:,:,:,:,:)
  real, allocatable :: DiffSemiCoef_VGB(:,:,:,:,:)
  real, allocatable :: PointSemiCoef_VCB(:,:,:,:,:)

  ! Index indicating which implicit variables involve diffusion
  integer, allocatable :: iDiff_I(:)
  integer :: nDiff

  ! Index indicating which implicit variables involve energy exchange
  ! with the electrons
  integer, allocatable :: iRelax_I(:)
  integer :: nRelax

  ! Number of point implicit coefficients
  integer :: nPoint
  ! Index which variable involve point implicit calculation
  integer :: iVarPoint
  ! Named indices for point implicit coefficients
  integer, parameter :: Relax_ = 1, Planck_ = 2

  ! radiation energy used for calculating radiative energy flux
  real, allocatable :: Erad_G(:,:,:)
  ! temporary radiation energy array needed by set_block_field
  real, allocatable :: Erad1_G(:,:,:)

  ! radiation has gamma of 4/3
  real, parameter :: GammaRel = 4.0/3.0

  ! Indices indicating which semi-implicit vars are mapped to which
  ! State_VGB variables. The total gas pressure and energy are excluded
  ! from this list, since they taken care of separately.
  integer :: nImplToState
  integer, allocatable :: iFromImpl_I(:), iToState_I(:)

  ! If we solve for the temperature, then ImplStateMin indicates the minimum
  ! allowed temperature that is used in the update from the implicit state
  ! variables to State_VGB.
  ! If the implicit solver uses energies (cRadiation*T**4),
  ! then ImplStateMin is the minimum allowed energy value
  real :: ImplStateMin

contains

  !============================================================================

  subroutine init_gray_diffusion

    use ModAdvance,     ONLY: Erad_, Ee_, UseElectronEnergy
    use ModMain,        ONLY: UseGrayDiffusion
    use ModSize,        ONLY: nI, nJ, nK, MaxBlock, nDim
    use ModImplicit,    ONLY: UseSemiImplicit, UseFullImplicit, &
         TypeSemiImplicit, iEradImpl, iTeImpl
    use ModPhysics,     ONLY: Si2No_V, UnitTemperature_, cRadiationNo
    use ModTemperature, ONLY: TradMinSi
    use ModWaves,       ONLY: UseWavePressure, GammaWave

    real :: EradMin

    character(len=*), parameter :: NameSub = "init_gray_diffusion"
    !------------------------------------------------------------------------

    if(allocated(Erad_G)) RETURN

    ! Make sure that Erad_ is correct
    if(UseGrayDiffusion)then
       if(Erad_ == 1) call stop_mpi(NameSub// &
            ": incorrect index for Erad variable in ModEquation")

       EradMin = cRadiationNo*(TradMinSi*Si2No_V(UnitTemperature_))**4
    end if

    allocate(Erad_G(-1:nI+2,-1:nJ+2,-1:nK+2))
    if(UseSemiImplicit) allocate(Erad1_G(0:nI+1,0:nJ+1,0:nK+1))

    if(UseFullImplicit)then
       nDiff = 1
       allocate(iDiff_I(nDiff))
       iDiff_I(1) = Erad_
       nRelax = 1
    end if
       
    if(UseSemiImplicit)then

       ! Default to zero, unless reset
       iTeImpl = 0; iEradImpl = 0

       select case(TypeSemiImplicit)
       case('radiation')
          if(UseElectronEnergy)then
             call stop_mpi(NameSub//": Te/=Ti requires heat conduction")
          end if
          iEradImpl = 1
          nDiff = 1
          allocate(iDiff_I(nDiff))
          iDiff_I(1) = iEradImpl
          nRelax = 0
          nPoint = 2
          iVarPoint = iEradImpl
          nImplToState = 1
          allocate(iFromImpl_I(nImplToState), iToState_I(nImplToState))
          iFromImpl_I = iEradImpl
          iToState_I  = Erad_
          ImplStateMin = EradMin

       case('radcond')
          iTeImpl = 1; iEradImpl = 2
          nDiff = 2
          allocate(iDiff_I(nDiff))
          iDiff_I = (/ iTeImpl, iEradImpl /)
          nRelax = 1
          allocate(iRelax_I(nRelax))
          iRelax_I(1) = iEradImpl
          if(UseElectronEnergy)then
             nPoint = 2
             iVarPoint = iTeImpl
             nImplToState = 2
             allocate(iFromImpl_I(nImplToState), iToState_I(nImplToState))
             iFromImpl_I = (/ iTeImpl, iEradImpl /)
             iToState_I  = (/ Ee_, Erad_ /)
          else
             nPoint = 0
             nImplToState = 1
             allocate(iFromImpl_I(nImplToState), iToState_I(nImplToState))
             iFromImpl_I = iEradImpl
             iToState_I  = Erad_
          end if
          ImplStateMin = EradMin

       case('cond')
          iTeImpl = 1
          nDiff = 1
          allocate(iDiff_I(nDiff))
          iDiff_I(1) = iTeImpl
          nRelax = 0
          if(UseElectronEnergy)then
             nPoint = 2
             iVarPoint = iTeImpl
             nImplToState = 1
             allocate(iFromImpl_I(nImplToState), iToState_I(nImplToState))
             iFromImpl_I = iTeImpl
             iToState_I  = Ee_
          else
             nPoint = 0
             nImplToState = 0
          end if
          ImplStateMin = 0.0

       end select

       if(nPoint>0)then
          allocate(PointSemiCoef_VCB(nPoint,nI,nJ,nK,MaxBlock))
          PointSemiCoef_VCB = 0.0
       end if

       allocate(DiffSemiCoef_VGB(nDiff,-1:nI+2,-1:nJ+2,-1:nK+2,MaxBlock))
       DiffSemiCoef_VGB = 0.0
    end if

    if(nRelax>0)then
       allocate(RelaxCoef_VCB(nRelax,nI,nJ,nK,MaxBlock))
       RelaxCoef_VCB = 0.0
    end if

    allocate(DiffCoef_VFDB(nDiff,1:nI+1,1:nJ+1,1:nK+1,nDim,MaxBlock))
    DiffCoef_VFDB = 0.0 ! make sure all elements are initialized

    ! Setup for wave infrastructure
    if(UseGrayDiffusion)then
       UseWavePressure = .true.
       GammaWave = GammaRel
    end if

  end subroutine init_gray_diffusion

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

    if(IsNewBlockGrayDiffusion) Erad_G = State_VGB(Erad_,:,:,:,iBlock)

    call get_face_gradient(iDir, i, j, k, iBlock, &
         IsNewBlockGrayDiffusion, Erad_G, FaceGrad_D)

    if(IsNewTimestepGrayDiffusion)then
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

      use ModAdvance,     ONLY: nOpacity
      use ModPhysics,     ONLY: Si2No_V, UnitX_, Clight
      use ModTemperature, ONLY: UseRadFluxLimiter, TypeRadFluxLimiter
      use ModUser,        ONLY: user_material_properties

      real, intent(in) :: State_V(nVar)
      real, intent(out):: DiffCoef

      real :: DiffusionOpacitySi_I(nOpacity), DiffusionOpacity, Grad2ByErad2
      !------------------------------------------------------------------------

      call user_material_properties(State_V, i, j, k, iBlock, iDir, &
           DiffusionOpacitySiOut_I = DiffusionOpacitySi_I)

      DiffusionOpacity = DiffusionOpacitySi_I(1)/Si2No_V(UnitX_)

      if(UseRadFluxLimiter)then
         Grad2ByErad2 = sum(FaceGrad_D**2)/State_V(Erad_)**2

         select case(TypeRadFluxLimiter)
         case("sum")
            DiffCoef = Clight/(3*DiffusionOpacity + sqrt(Grad2ByErad2))
         case("max")
            DiffCoef = Clight/max(3*DiffusionOpacity,sqrt(Grad2ByErad2))
         case("larsen")
            DiffCoef = Clight/sqrt(9*DiffusionOpacity**2 + Grad2ByErad2)
         end select
      else
         DiffCoef = Clight/(3*DiffusionOpacity)
      end if

    end subroutine get_diffusion_coef

  end subroutine get_radiation_energy_flux

  !============================================================================

  subroutine calc_source_gray_diffusion(iBlock)

    use ModAdvance,    ONLY: State_VGB, Source_VC, Erad_, nOpacity
    use ModConst,      ONLY: cLightSpeed
    use ModPhysics,    ONLY: cRadiationNo, Si2No_V, UnitTemperature_, UnitT_
    use ModMain,       ONLY: nI, nJ, nK
    use ModUser,       ONLY: user_material_properties
    use ModVarIndexes, ONLY: Energy_

    integer, intent(in) :: iBlock

    integer :: i, j, k
    real :: TeSi, Te
    real :: AbsorptionEmission, PlanckOpacitySi_I(nOpacity)
    character(len=*), parameter:: NameSub = "calc_source_gray_diffusion"
    !------------------------------------------------------------------------

    do k=1,nK; do j=1,nJ; do i=1,nI

       if(IsNewTimestepGrayDiffusion)then
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               i, j, k, iBlock, AbsorptionOpacitySiOut_I = PlanckOpacitySi_I)

          RelaxCoef_VCB(1,i,j,k,iBlock) = &
               PlanckOpacitySi_I(1)*cLightSpeed/Si2No_V(UnitT_)
       end if

       call user_material_properties(State_VGB(:,i,j,k,iBlock), &
            i, j, k, iBlock, TeSiOut = TeSi)

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

  end subroutine calc_source_gray_diffusion

  !============================================================================
  ! Semi-implicit interface
  !============================================================================

  subroutine get_impl_gray_diff_state(StateImpl_VGB,DconsDsemi_VCB)

    use ModAdvance,  ONLY: Erad_, State_VGB, UseElectronEnergy, nOpacity
    use ModConst,    ONLY: cBoltzmann
    use ModImplicit, ONLY: nw, nImplBlk, impl2iBlk, TypeSemiImplicit, &
         iEradImpl, iTeImpl, ImplCoeff
    use ModMain,     ONLY: nDim, x_, y_, nI, nJ, nK, MaxImplBlk, Dt
    use ModNumConst, ONLY: i_DD
    use ModPhysics,  ONLY: inv_gm1, Clight, cRadiationNo, UnitN_, UnitP_, &
         Si2No_V, UnitTemperature_, UnitEnergyDens_, UnitX_, UnitU_, UnitT_, &
         No2Si_V
    use ModUser,     ONLY: user_material_properties
    use ModTemperature, ONLY: UseRadFluxLimiter, TypeRadFluxLimiter
    use ModGeometry, ONLY: dx_BLK, dy_BLK, dz_BLK, vInv_CB, &
         UseCovariant, TypeGeometry, FaceAreaI_DFB, FaceAreaJ_DFB
    use ModParallel, ONLY: NOBLK, NeiLev
    use ModMessagePass, ONLY: message_pass_dir

    real, intent(out) :: StateImpl_VGB(nw,0:nI+1,0:nJ+1,0:nK+1,MaxImplBlk)
    real, intent(inout) :: DconsDsemi_VCB(nw,nI,nJ,nK,MaxImplBlk)

    integer :: iImplBlock, iBlock, i, j, k
    real :: PlanckOpacitySi_I(nOpacity), DiffusionOpacitySi_I(nOpacity)
    real :: DiffusionOpacity, PlanckOpacity, CvSi, Cv, TeSi, Te
    real :: HeatCondSi, HeatCond, TeTiRelaxSi, TeTiRelax
    real :: Grad2ByErad2, DiffRad, InvDx2, InvDy2, InvDz2
    real :: InvDx, InvDy
    real :: NatomicSi, CveSi, Cve, CviSi, Cvi, PiSi, TiSi, Ti, PeSi
    real :: TeTiCoefSi, TeTiCoef, TeTiCoefPrime

    integer :: iDim, Di, Dj, Dk, iDiff, nDimInUse
    real :: Coeff, Dxyz_D(3)

    character(len=*), parameter:: NameSub='get_impl_gray_diff_state'
    !--------------------------------------------------------------------------

    do iImplBlock = 1, nImplBLK

       iBlock = impl2iBLK(iImplBlock)
       IsNewBlockGrayDiffusion = .true.

       if(TypeSemiImplicit=='radiation' .or. TypeSemiImplicit=='radcond')then
          do k = 0, nK+1; do j = 0, nJ+1; do i = 0, nI+1
             StateImpl_VGB(iEradImpl,i,j,k,iImplBlock) = &
                  State_VGB(Erad_,i,j,k,iBlock)
          end do; end do; end do
       end if
       if(TypeSemiImplicit=='radcond')then
          do k = 0, nK+1; do j = 0, nJ+1; do i = 0, nI+1
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  i, j, k, iBlock, TeSiOut = TeSi)
             Te = TeSi*Si2No_V(UnitTemperature_)
             StateImpl_VGB(iTeImpl,i,j,k,iImplBlock) = cRadiationNo*Te**4
          end do; end do; end do
       elseif(TypeSemiImplicit=='cond')then
          do k = 0, nK+1; do j = 0, nJ+1; do i = 0, nI+1
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  i, j, k, iBlock, TeSiOut = TeSi)
             StateImpl_VGB(iTeImpl,i,j,k,iImplBlock) = &
                  TeSi*Si2No_V(UnitTemperature_)
          end do; end do; end do
       end if

       InvDx2 = 0.5/dx_BLK(iBlock)
       InvDy2 = 0.5/dy_BLK(iBlock)
       InvDz2 = 0.5/dz_BLK(iBlock)

       ! calculate coefficients for linearized energy exchange and diffusion
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(UseElectronEnergy)then
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  i, j, k, iBlock, &
                  AbsorptionOpacitySiOut_I = PlanckOpacitySi_I, &
                  DiffusionOpacitySiOut_I = DiffusionOpacitySi_I, &
                  CvSiOut=CveSi, TeSiOut = TeSi, NatomicSiOut=NatomicSi, &
                  HeatCondSiOut = HeatCondSi, TeTiRelaxSiOut = TeTiRelaxSi)

             PlanckOpacity = PlanckOpacitySi_I(1)/Si2No_V(UnitX_)
             Te = TeSi*Si2No_V(UnitTemperature_)
             Cve = CveSi*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_)
             PiSi = State_VGB(p_,i,j,k,iBlock)*No2Si_V(UnitP_)
             TiSi = PiSi/(cBoltzmann*NatomicSi)
             Ti = TiSi*Si2No_V(UnitTemperature_)
             CviSi = inv_gm1*cBoltzmann*NatomicSi
             Cvi = CviSi*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_)
             TeTiCoefSi = TeTiRelaxSi*cBoltzmann*NatomicSi
             TeTiCoef = TeTiCoefSi*Si2No_V(UnitEnergyDens_) &
                  /(Si2No_V(UnitTemperature_)*Si2No_V(UnitT_))
          else
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  i, j, k, iBlock, &
                  AbsorptionOpacitySiOut_I = PlanckOpacitySi_I, &
                  DiffusionOpacitySiOut_I = DiffusionOpacitySi_I, &
                  CvSiOut = CvSi, TeSiOut = TeSi, HeatCondSiOut = HeatCondSi)

             PlanckOpacity = PlanckOpacitySi_I(1)/Si2No_V(UnitX_)
             Cv = CvSi*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_)
             Te = TeSi*Si2No_V(UnitTemperature_)
          end if

          select case(TypeSemiImplicit)
          case('radiation')
             ! This coefficient is cR'' = cR/(1+dt*cR*dPlanck/dEint)
             PointSemiCoef_VCB(Relax_,i,j,k,iBlock) = Clight*PlanckOpacity  &
                  /(1 + ImplCoeff*Dt*Clight*PlanckOpacity &
                  *4.0*cRadiationNo*Te**3 / Cv)

             ! This is just the Planck function at time level * saved
             PointSemiCoef_VCB(Planck_,i,j,k,iBlock) = cRadiationNo*Te**4

          case('radcond')
             if(UseElectronEnergy)then
                DconsDsemi_VCB(iTeImpl,i,j,k,iImplBlock) = &
                     Cve/(4.0*cRadiationNo*Te**3)

                RelaxCoef_VCB(1,i,j,k,iBlock) = Clight*PlanckOpacity

                TeTiCoefPrime = TeTiCoef*4.0*Ti**3/Cvi &
                     /((Te+Ti)*(Te**2+Ti**2))
                PointSemiCoef_VCB(Relax_,i,j,k,iBlock) = &
                     TeTiCoefPrime*Cvi/(4.0*cRadiationNo*Ti**3) &
                     /(1.0 + ImplCoeff*Dt*TeTiCoefPrime)
                PointSemiCoef_VCB(Planck_,i,j,k,iBlock) = cRadiationNo*Ti**4
             else
                RelaxCoef_VCB(1,i,j,k,iBlock) = Clight*PlanckOpacity

                DconsDsemi_VCB(iTeImpl,i,j,k,iImplBlock) = &
                     Cv/(4.0*cRadiationNo*Te**3)
             end if

          case('cond')
             if(UseElectronEnergy)then
                DconsDsemi_VCB(iTeImpl,i,j,k,iImplBlock) = Cve

                PointSemiCoef_VCB(Relax_,i,j,k,iBlock) = TeTiCoef &
                     /(1.0 + ImplCoeff*Dt*TeTiCoef/Cvi)
                PointSemiCoef_VCB(Planck_,i,j,k,iBlock) = Ti
             else
                DconsDsemi_VCB(iTeImpl,i,j,k,iImplBlock) = Cv
             end if
          end select

          call get_diffusion_coef

       end do; end do; end do

       if(NeiLev(1,iBlock) == NOBLK)then
          i = 0
          do k = 1, nK; do j = 1, nJ
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  i, j, k, iBlock, &
                  DiffusionOpacitySiOut_I = DiffusionOpacitySi_I, &
                  HeatCondSiOut = HeatCondSi, &
                  TeSiOut = TeSi)
             call get_diffusion_coef
          end do; end do
       end if
       if(NeiLev(2,iBlock) == NOBLK)then
          i = nI + 1
          do k = 1, nK; do j = 1, nJ
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  i, j, k, iBlock, &
                  DiffusionOpacitySiOut_I = DiffusionOpacitySi_I, &
                  HeatCondSiOut = HeatCondSi, &
                  TeSiOut = TeSi)
             call get_diffusion_coef
          end do; end do
       end if
       if(NeiLev(3,iBlock) == NOBLK)then
          j = 0
          do k = 1, nK; do i = 1, nI
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  i, j, k, iBlock, &
                  DiffusionOpacitySiOut_I = DiffusionOpacitySi_I, &
                  HeatCondSiOut = HeatCondSi, &
                  TeSiOut = TeSi)
             call get_diffusion_coef
          end do; end do
       end if
       if(NeiLev(4,iBlock) == NOBLK)then
          j = nJ + 1
          do k = 1, nK; do i = 1, nI
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  i, j, k, iBlock, &
                  DiffusionOpacitySiOut_I = DiffusionOpacitySi_I, &
                  HeatCondSiOut = HeatCondSi, &
                  TeSiOut = TeSi)
             call get_diffusion_coef
          end do; end do
       end if
       if(NeiLev(5,iBlock) == NOBLK)then
          k = 0
          do j = 1, nJ; do i = 1, nI
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  i, j, k, iBlock, &
                  DiffusionOpacitySiOut_I = DiffusionOpacitySi_I, &
                  HeatCondSiOut = HeatCondSi, &
                  TeSiOut = TeSi)
             call get_diffusion_coef
          end do; end do
       end if
       if(NeiLev(6,iBlock) == NOBLK)then
          k = nK + 1
          do j = 1, nJ; do i = 1, nI
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  i, j, k, iBlock, &
                  DiffusionOpacitySiOut_I = DiffusionOpacitySi_I, &
                  HeatCondSiOut = HeatCondSi, &
                  TeSiOut = TeSi)
             call get_diffusion_coef
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

      use ModFaceGradient, ONLY: set_block_field

      !------------------------------------------------------------------------

      if(TypeSemiImplicit=='radiation' .or. TypeSemiImplicit=='radcond')then
         DiffusionOpacity = DiffusionOpacitySi_I(1)/Si2No_V(UnitX_)

         ! Calculate the cell centered diffusion coefficients
         if(UseRadFluxLimiter)then

            if(IsNewBlockGrayDiffusion)then
               Erad_G = State_VGB(Erad_,:,:,:,iBlock)
               call set_block_field(iBlock, 1, Erad1_G, Erad_G)

               IsNewBlockGrayDiffusion = .false.
            end if
        
            Grad2ByErad2 = &
                 (((Erad_G(i+1,j,k) - Erad_G(i-1,j,k))*InvDx2)**2 &
                 +((Erad_G(i,j+1,k) - Erad_G(i,j-1,k))*InvDy2)**2 &
                 +((Erad_G(i,j,k+1) - Erad_G(i,j,k-1))*InvDz2)**2 &
                 )/ Erad_G(i,j,k)**2

            select case(TypeRadFluxLimiter)
            case("sum")
               DiffRad = Clight/(3*DiffusionOpacity + sqrt(Grad2ByErad2))
            case("max")
               DiffRad = Clight/max(3*DiffusionOpacity,sqrt(Grad2ByErad2))
            case("larsen")
               DiffRad = Clight/sqrt(9*DiffusionOpacity**2 + Grad2ByErad2)
            end select
         else
            DiffRad = Clight/(3*DiffusionOpacity)
         end if

         ! Store it for message passing
         DiffSemiCoef_VGB(iEradImpl,i,j,k,iBlock) = DiffRad
      end if

      if(TypeSemiImplicit=='cond')then
         HeatCond = HeatCondSi &
              *Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_) &
              *Si2No_V(UnitU_)*Si2No_V(UnitX_)

         DiffSemiCoef_VGB(iTeImpl,i,j,k,iBlock) = HeatCond
      elseif(TypeSemiImplicit=='radcond')then
         HeatCond = HeatCondSi &
              *Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_) &
              *Si2No_V(UnitU_)*Si2No_V(UnitX_)
         Te = TeSi*Si2No_V(UnitTemperature_)

         DiffSemiCoef_VGB(iTeImpl,i,j,k,iBlock) = &
              HeatCond/(4.0*cRadiationNo*Te**3)
      end if

    end subroutine get_diffusion_coef

  end subroutine get_impl_gray_diff_state

  !============================================================================

  subroutine set_gray_outflow_bc(iSide, iBlock, State_VG)

    use ModImplicit, ONLY: nw, iTeImpl
    use ModGeometry, ONLY: dx_BLK, dy_BLK, dz_BLK
    use ModMain,     ONLY: nI, nJ, nK
    use ModPhysics,  ONLY: Clight
    use ModUser,     ONLY: user_material_properties

    integer, intent(in) :: iSide, iBlock
    real, intent(inout) :: State_VG(nw,-1:nI+2,-1:nJ+2,-1:nK+2)

    integer :: iVar, i, j, k, iDiff
    real :: Coef
    character(len=*), parameter :: NameSub='set_gray_outflow_bc'
    !--------------------------------------------------------------------------
    select case(iSide)
    case(1)
       do k = 1, nK; do j = 1, nJ
          do iDiff = 1, nDiff
             iVar = iDiff_I(iDiff)          
             if(iVar == iTeImpl) CYCLE
             Coef = 2/Clight &
                  *DiffSemiCoef_VGB(iDiff,1,j,k,iBlock)/dx_BLK(iBlock)
             State_VG(iVar,0,j,k) = State_VG(iVar,1,j,k) &
                  *(Coef - 0.5)/(Coef + 0.5)
          end do
       end do; end do
    case(2)
       do k = 1, nK; do j = 1, nJ
          do iDiff = 1, nDiff
             iVar = iDiff_I(iDiff)
             if(iVar == iTeImpl) CYCLE
             Coef = 2/Clight &
                  *DiffSemiCoef_VGB(iDiff,nI,j,k,iBlock)/dx_BLK(iBlock)
             State_VG(iVar,nI+1,j,k) = State_VG(iVar,nI,j,k) &
                  *(Coef - 0.5)/(Coef + 0.5)
          end do
       end do; end do
    case(3)
       do k = 1, nK; do i = 1, nI
          do iDiff = 1, nDiff
             iVar = iDiff_I(iDiff)
             if(iVar == iTeImpl) CYCLE
             Coef = 2/Clight &
                  *DiffSemiCoef_VGB(iDiff,i,1,k,iBlock)/dy_BLK(iBlock)
             State_VG(iVar,i,0,k) = State_VG(iVar,i,1,k) &
                  *(Coef - 0.5)/(Coef + 0.5)
          end do
       end do; end do
    case(4)
       do k = 1, nK; do i = 1, nI
          do iDiff = 1, nDiff
             iVar = iDiff_I(iDiff)
             if(iVar == iTeImpl) CYCLE
             Coef = 2/Clight &
                  *DiffSemiCoef_VGB(iDiff,i,nJ,k,iBlock)/dy_BLK(iBlock)
             State_VG(iVar,i,nJ+1,k) = State_VG(iVar,i,nJ,k) &
                  *(Coef - 0.5)/(Coef + 0.5)
          end do
       end do; end do
    case(5)
       do j = 1, nJ; do i = 1, nI
          do iDiff = 1, nDiff
             iVar = iDiff_I(iDiff)
             if(iVar == iTeImpl) CYCLE
             Coef = 2/Clight &
                  *DiffSemiCoef_VGB(iDiff,i,j,1,iBlock)/dz_BLK(iBlock)
             State_VG(iVar,i,j,0) = State_VG(iVar,i,j,1) &
                  *(Coef - 0.5)/(Coef + 0.5)
          end do
       end do; end do
    case(6)
       do k = j, nJ; do i = 1, nI
          do iDiff = 1, nDiff
             iVar = iDiff_I(iDiff)
             if(iVar == iTeImpl) CYCLE
             Coef = 2/Clight &
                  *DiffSemiCoef_VGB(iDiff,i,j,nK,iBlock)/dz_BLK(iBlock)
             State_VG(iVar,i,j,nK+1) = State_VG(iVar,i,j,nK) &
                  *(Coef - 0.5)/(Coef + 0.5)
          end do
       end do; end do
    end select

  end subroutine set_gray_outflow_bc

  !============================================================================

  subroutine get_gray_diffusion_rhs(iBlock, StateImpl_VG, Rhs_VC, IsLinear)

    use ModGeometry, ONLY: TypeGeometry, vInv_CB
    use ModImplicit, ONLY: nw, iTeImpl
    use ModMain,     ONLY: nI, nJ, nK
    use ModParallel, ONLY: NeiLev

    integer, intent(in) :: iBlock
    real, intent(inout) :: StateImpl_VG(nw,-1:nI+2,-1:nJ+2,-1:nK+2)
    real, intent(out)   :: Rhs_VC(nw,nI,nJ,nK)
    logical, intent(in) :: IsLinear

    real :: Te, EnergyExchange
    integer :: i, j, k, iDiff, iRelax, iVar
    character(len=*), parameter :: NameSub='get_gray_diffusion_rhs'
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

    !!! Rhs_VC = 0.0

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
    if(nRelax>0)then
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          do iRelax = 1, nRelax
             iVar = iRelax_I(iRelax)

             EnergyExchange = RelaxCoef_VCB(iRelax,i,j,k,iBlock) &
                  *(StateImpl_VG(iTeImpl,i,j,k) - StateImpl_VG(iVar,i,j,k))

             ! dEvar/dt = + EnergyExchange
             Rhs_VC(iVar,i,j,k)    = Rhs_VC(iVar,i,j,k)    + EnergyExchange

             ! dEe/dt   = - EnergyExchange
             Rhs_VC(iTeImpl,i,j,k) = Rhs_VC(iTeImpl,i,j,k) - EnergyExchange
          end do
       end do; end do; end do
    end if
       
    ! Point implicit source terms due to energy exchange
    if(nPoint>0)then
       if(IsLinear)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             Rhs_VC(iVarPoint,i,j,k) = Rhs_VC(iVarPoint,i,j,k) &
                  - PointSemiCoef_VCB(Relax_,i,j,k,iBlock) &
                  *StateImpl_VG(iVarPoint,i,j,k)
          end do; end do; end do
       else
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             EnergyExchange = PointSemiCoef_VCB(Relax_,i,j,k,iBlock) &
                  * (PointSemiCoef_VCB(Planck_,i,j,k,iBlock) &
                  -  StateImpl_VG(iVarPoint,i,j,k))

             Rhs_VC(iVarPoint,i,j,k) = Rhs_VC(iVarPoint,i,j,k) + EnergyExchange
          end do; end do; end do
       end if
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

  end subroutine get_gray_diffusion_rhs

  !============================================================================

  subroutine add_jacobian_gray_diff(iBlock, nVar, Jacobian_VVCI)

    use ModGeometry, ONLY: vInv_CB, dx_BLK, dy_BLK, dz_BLK, &
         fAx_BLK, fAy_BLK, fAz_BLK
    use ModImplicit, ONLY: TypeSemiImplicit, iTeImpl, UseFullImplicit, &
         UseSemiImplicit
    use ModMain,     ONLY: nI, nJ, nK, nDim
    use ModNumConst, ONLY: i_DD

    integer, parameter:: nStencil = 2*nDim + 1

    integer, intent(in) :: iBlock, nVar
    real, intent(inout) :: Jacobian_VVCI(nVar,nVar,nI,nJ,nK,nStencil)

    integer :: iVar, i, j, k, iDim, Di, Dj, Dk, iDiff, iRelax
    real :: DiffLeft, DiffRight, RelaxCoef
    real :: Dxyz_D(nDim), Area_D(nDim), Coeff
    !--------------------------------------------------------------------------

    if(UseSemiImplicit)then
       if(nPoint>0)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             ! dSvar/dVar (diagonal)
             Jacobian_VVCI(iVarPoint,iVarPoint,i,j,k,1) = &
                  -PointSemiCoef_VCB(Relax_,i,j,k,iBlock)
          end do; end do; end do
       end if

       if(nRelax>0)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             do iRelax = 1, nRelax
                iVar = iRelax_I(iRelax)

                RelaxCoef = RelaxCoef_VCB(iRelax,i,j,k,iBlock)

                ! dSvar/dVar (diagonal)
                Jacobian_VVCI(iVar,iVar,i,j,k,1) = -RelaxCoef

                ! dSe/dVar (off diagonal)
                Jacobian_VVCI(iTeImpl,iVar,i,j,k,1) = +RelaxCoef

                ! dSe/daTe^4 (diagonal)
                Jacobian_VVCI(iTeImpl,iTeImpl,i,j,k,1) = -RelaxCoef

                ! dSvar/daTe^4 (off diagonal)
                Jacobian_VVCI(iVar,iTeImpl,i,j,k,1) = +RelaxCoef
             end do
          end do; end do; end do
       end if
    end if

    ! For the fully implicit scheme:
    ! add partial derivatives of the gray diffusion term to the Jacobian that
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

  end subroutine add_jacobian_gray_diff

  !============================================================================

  subroutine update_impl_gray_diff(iBlock, iImplBlock, StateImpl_VG)

    use ModAdvance,  ONLY: State_VGB, Rho_, p_, ExtraEint_, &
         UseElectronEnergy, Ee_
    use ModEnergy,   ONLY: calc_energy_cell
    use ModImplicit, ONLY: nw, iTeImpl, DconsDsemi_VCB, ImplOld_VCB, ImplCoeff
    use ModMain,     ONLY: nI, nJ, nK, Dt
    use ModPhysics,  ONLY: inv_gm1, g, No2Si_V, Si2No_V, UnitEnergyDens_, &
         UnitP_, UnitRho_, UnitTemperature_
    use ModUser,     ONLY: user_material_properties

    integer, intent(in) :: iBlock, iImplBlock
    real, intent(inout) :: StateImpl_VG(nw,nI,nJ,nK)

    integer :: i, j, k, iImplToState, iVarImpl, iVar
    real :: Einternal, EinternalSi, PressureSi
    real :: PeSi, Ee, EeSi
    real :: DelEnergy, Relaxation

    character(len=*), parameter :: NameSub = 'update_impl_gray_diff'
    !--------------------------------------------------------------------------

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       ! minimum allowed temperature for all implicit variables
       do iVarImpl = 1, nw
          StateImpl_VG(iVarImpl,i,j,k) = &
               max(ImplStateMin, StateImpl_VG(iVarImpl,i,j,k))
       end do

       do iImplToState = 1, nImplToState
          iVarImpl = iFromImpl_I(iImplToState)
          iVar     = iToState_I(iImplToState)

          DelEnergy = DconsDsemi_VCB(iVarImpl,i,j,k,iImplBlock) &
               *(StateImpl_VG(iVarImpl,i,j,k) &
               - ImplOld_VCB(iVarImpl,i,j,k,iBlock))

          State_VGB(iVar,i,j,k,iBlock) = &
               State_VGB(iVar,i,j,k,iBlock) + DelEnergy
       end do

       if(UseElectronEnergy)then
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

       if(nPoint>0)then
          Relaxation = PointSemiCoef_VCB(Relax_,i,j,k,iBlock)*( &
               + ImplCoeff*(StateImpl_VG(iVarPoint,i,j,k) &
               -            PointSemiCoef_VCB(Planck_,i,j,k,iBlock)) &
               + (1.0 - ImplCoeff)*(ImplOld_VCB(iVarPoint,i,j,k,iBlock) &
               -            PointSemiCoef_VCB(Planck_,i,j,k,iBlock)) )

          Einternal = Einternal + Dt*Relaxation
       end if

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

       if(UseElectronEnergy)then
          ! ions
          State_VGB(p_,i,j,k,iBlock) = (g - 1)*Einternal

          ! electrons
          Ee = State_VGB(Ee_,i,j,k,iBlock) + State_VGB(ExtraEint_,i,j,k,iBlock)
          EeSi = Ee*No2Si_V(UnitEnergyDens_)

          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               i, j, k, iBlock, &
               EinternalSiIn = EeSi, PressureSiOut = PeSi)

          ! use true electron pressure
          State_VGB(Ee_,i,j,k,iBlock) = inv_gm1*PeSi*Si2No_V(UnitP_)

          ! Set ExtraEint = electron internal energy - Pe/(gamma -1)
          State_VGB(ExtraEint_,i,j,k,iBlock) = &
               Ee - State_VGB(Ee_,i,j,k,iBlock)

          if(State_VGB(ExtraEint_,i,j,k,iBlock)<0.0)then
             write(*,*)NameSub,': ERROR extra internal energy =', &
                  State_VGB(ExtraEint_,i,j,k,iBlock)
             write(*,*)NameSub,': ERROR at i,j,k,iBlock=', i, j, k, iBlock
             call stop_mpi(NameSub//': ERROR negative extra internal energy')
          end if

       else
          ! ions + electrons
          EinternalSi = Einternal*No2Si_V(UnitEnergyDens_)

          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               i, j, k, iBlock, &
               EinternalSiIn = EinternalSi, PressureSiOut = PressureSi)

          State_VGB(p_,i,j,k,iBlock) = PressureSi*Si2No_V(UnitP_)

          State_VGB(ExtraEint_,i,j,k,iBlock) = &
               Einternal - inv_gm1*State_VGB(p_,i,j,k,iBlock)

       end if

    end do; end do; end do

    call calc_energy_cell(iBlock)

  end subroutine update_impl_gray_diff

end module ModGrayDiffusion
