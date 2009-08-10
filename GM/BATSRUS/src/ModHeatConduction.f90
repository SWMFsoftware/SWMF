!^CFG COPYRIGHT UM
!^CFG FILE IMPLICIT
!==============================================================================
module ModHeatConduction

  implicit none
  save

  private ! except

  ! Public methods
  public :: read_heatconduction_param
  public :: init_heat_conduction
  public :: get_heat_flux
  public :: get_impl_heat_cond_state
  public :: get_heat_conduction_bc
  public :: get_heat_conduction_rhs
  public :: get_heat_cond_jacobian
  public :: update_impl_heat_cond

  ! Logical for adding field-aligned heat conduction
  logical, public :: IsNewBlockHeatConduction = .true.

  ! Variables for setting the field-aligned heat conduction coefficient
  character(len=20), public :: TypeHeatConduction = 'test'
  logical :: DoModifyHeatConduction, DoTestHeatConduction
  
  real :: HeatConductionParSi = 1.23e-11  ! taken from calc_heat_flux
  real :: TmodifySi = 3.0e5, DeltaTmodifySi = 2.0e4
  real :: HeatConductionPar, Tmodify, DeltaTmodify

  logical :: DoWeakFieldConduction = .false.
  real :: BmodifySi = 1.0e-7, DeltaBmodifySi = 1.0e-8 ! modify about 1 mG
  real :: Bmodify, DeltaBmodify

  ! electron temperature used for calculating heat flux
  real, allocatable :: Te_G(:,:,:)

  ! ratio of Te Te+T used for ideal EOS
  real :: TeFraction

  ! Heat flux for operator split scheme
  real, allocatable :: FluxImpl_X(:,:,:), FluxImpl_Y(:,:,:), FluxImpl_Z(:,:,:)

  ! Heat conduction dyad pre-multiplied by the face area
  real, allocatable :: HeatCond_DFDB(:,:,:,:,:,:)

contains

  !============================================================================

  subroutine read_heatconduction_param(NameCommand)

    use ModMain,      ONLY: UseParallelConduction
    use ModReadParam, ONLY: read_var

    character(len=*), intent(in) :: NameCommand

    character(len=*), parameter :: NameSub = 'read_heatconduction_param'
    !--------------------------------------------------------------------------

    select case(NameCommand)
    case("#PARALLELCONDUCTION")
       call read_var('UseParallelConduction', UseParallelConduction)
       if(UseParallelConduction)then
          call read_var('TypeHeatConduction', TypeHeatConduction)
          call read_var('HeatConductionParSi', HeatConductionParSi)

          select case(TypeHeatConduction)
          case('test','spitzer')
          case('modified')
             call read_var('TmodifySi', TmodifySi)
             call read_var('DeltaTmodifySi', DeltaTmodifySi)
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

    case default
       call stop_mpi(NameSub//' invalid NameCommand='//NameCommand)
    end select

  end subroutine read_heatconduction_param

  !============================================================================

  subroutine init_heat_conduction

    use ModImplicit,   ONLY: UseSemiImplicit, iTeImpl
    use ModMain,       ONLY: nI, nJ, nK, MaxBlock, nDim
    use ModMultiFluid, ONLY: MassIon_I
    use ModPhysics,    ONLY: Si2No_V, UnitEnergyDens_, UnitTemperature_, &
         UnitU_, UnitX_, UnitB_, ElectronTemperatureRatio, AverageIonCharge

    character(len=*), parameter :: NameSub = 'init_heat_conduction'
    !--------------------------------------------------------------------------

    write(*,*) "The field aligned heat conduction is under construction"
    write(*,*) "It is currently guaranteed to give wrong results !!!"
    
    if(allocated(Te_G)) RETURN

    allocate(Te_G(-1:nI+2,-1:nJ+2,-1:nK+2))

    ! Get TeFraction=Te/(Te+T) from P = n*k*T + ne*k*Te for ideal EOS
    TeFraction = MassIon_I(1)*ElectronTemperatureRatio &
         /(1 + AverageIonCharge*ElectronTemperatureRatio)

    DoTestHeatConduction = .false.
    DoModifyHeatConduction = .false.

    if(TypeHeatConduction == 'test')then
       DoTestHeatConduction = .true.
    elseif(TypeHeatConduction == 'modified')then
       DoModifyHeatConduction = .true.
    end if

    ! unit HeatConductionParSi is W/(m*K^(7/2))
    HeatConductionPar = HeatConductionParSi &
         *Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_)**3.5 &
         *Si2No_V(UnitU_)*Si2No_V(UnitX_)

    if(DoModifyHeatConduction)then
       Tmodify = TmodifySi*Si2No_V(UnitTemperature_)
       DeltaTmodify = DeltaTmodifySi*Si2No_V(UnitTemperature_)
    end if

    if(DoWeakFieldConduction)then
       Bmodify = BmodifySi*Si2No_V(UnitB_)
       DeltaBmodify = DeltaBmodifySi*Si2No_V(UnitB_)
    end if

    if(UseSemiImplicit)then
       allocate( &
            FluxImpl_X(nI+1,nJ,nK), &
            FluxImpl_Y(nI,nJ+1,nK), &
            FluxImpl_Z(nI,nJ,nK+1), &
            HeatCond_DFDB(nDim,nI+1,nJ+1,nK+1,nDim,MaxBlock) )

       iTeImpl = 1
    end if

  end subroutine init_heat_conduction

  !============================================================================

  subroutine get_heat_flux(iDir, i, j, k, iBlock, State_V, Normal_D, &
       HeatCondCoefNormal, HeatFlux)

    use ModAdvance,      ONLY: State_VGB, UseIdealEos
    use ModFaceGradient, ONLY: get_face_gradient
    use ModMain,         ONLY: nI, nJ, nK
    use ModPhysics,      ONLY: inv_gm1, Si2No_V, UnitTemperature_, &
         UnitEnergyDens_
    use ModUser,         ONLY: user_material_properties
    use ModVarIndexes,   ONLY: nVar, Rho_, p_

    integer, intent(in) :: iDir, i, j, k, iBlock
    real,    intent(in) :: State_V(nVar), Normal_D(3)
    real,    intent(out):: HeatCondCoefNormal, HeatFlux

    integer :: ii, jj, kk
    real :: FaceGrad_D(3), HeatCond_D(3), TeSi, Cv, CvSi

    character(len=*), parameter :: NameSub = 'get_heat_flux'
    !--------------------------------------------------------------------------

    if(IsNewBlockHeatConduction)then
       if(UseIdealEos)then
          Te_G = State_VGB(p_,:,:,:,iBlock)/State_VGB(Rho_,:,:,:,iBlock) &
               *TeFraction
       else
          do kk = -1, nK+2; do jj = -1, nJ+2; do ii = -1, nI+2
             call user_material_properties( &
                  State_VGB(:,ii,jj,kk,iBlock), TeSiOut=TeSi)
             Te_G(ii,jj,kk) = TeSi*Si2No_V(UnitTemperature_)
          end do; end do; end do
       end if
    end if

    call get_face_gradient(iDir, i, j, k, iBlock, &
         IsNewBlockHeatConduction, Te_G, FaceGrad_D)

    call get_heat_conduction_coef(iDir, i, j, k, iBlock, State_V, &
         Normal_D, HeatCond_D)

    HeatFlux = -sum(HeatCond_D*FaceGrad_D)

    ! get the heat conduction coefficient normal to the face for
    ! time step restriction
    if(UseIdealEos)then
       Cv = inv_gm1*State_V(Rho_)/TeFraction
    else
       call user_material_properties(State_V, CvSiOut = CvSi)
       Cv = CvSi*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_)
    end if
    HeatCondCoefNormal = sum(HeatCond_D*Normal_D)/Cv

  end subroutine get_heat_flux

  !============================================================================

  subroutine get_heat_conduction_coef(iDim, i, j, k, iBlock, State_V, &
       Normal_D, HeatCond_D)

    use ModAdvance,      ONLY: State_VGB, UseIdealEos
    use ModB0,           ONLY: B0_DX, B0_DY, B0_DZ
    use ModMain,         ONLY: UseB0
    use ModNumConst,     ONLY: cTolerance
    use ModPhysics,      ONLY: Si2No_V, UnitTemperature_
    use ModUser,         ONLY: user_material_properties
    use ModVarIndexes,   ONLY: nVar, Bx_, Bz_, Rho_, p_

    integer, intent(in) :: iDim, i, j, k, iBlock
    real, intent(in) :: State_V(nVar), Normal_D(3)
    real, intent(out):: HeatCond_D(3)

    real :: B_D(3), Bnorm, Bunit_D(3), TeSi, Te
    real :: HeatCoef, FractionSpitzer, FractionFieldAligned
    !--------------------------------------------------------------------------

    if(UseB0)then
       select case(iDim)
       case(1)
          B_D = State_V(Bx_:Bz_) + B0_DX(:,i,j,k)
       case(2)
          B_D = State_V(Bx_:Bz_) + B0_DY(:,i,j,k)
       case(3)
          B_D = State_V(Bx_:Bz_) + B0_DZ(:,i,j,k)
       end select
    else
       B_D = State_V(Bx_:Bz_)
    end if

    ! The magnetic field should nowhere be zero. The following fix will
    ! turn the magnitude of the field direction to zero.
    Bnorm = sqrt(sum(B_D**2))
    Bunit_D = B_D/max(Bnorm,cTolerance)

    if(UseIdealEos)then
       Te = TeFraction*State_V(p_)/State_V(Rho_)
    else
       ! Note we assume that the heat conduction formula for the
       ! ideal state is still applicable for the non-ideal state
       call user_material_properties(State_V, TeSiOut=TeSi)
       Te = TeSi*Si2No_V(UnitTemperature_)
    end if

    if(DoTestHeatConduction)then
       HeatCoef = 1.0
    else

       if(DoModifyHeatConduction)then
          ! Artificial modified heat conduction for a smoother transition
          ! region, Linker et al. (2001)
          FractionSpitzer = 0.5*(1.0+tanh((Te-Tmodify)/DeltaTmodify))
          HeatCoef = HeatConductionPar*(FractionSpitzer*Te**2.5 &
               + (1.0 - FractionSpitzer)*Tmodify**2.5)
       else
          ! Spitzer form for collisional regime
          HeatCoef = HeatConductionPar*Te**2.5
       end if
    end if

    if(DoWeakFieldConduction)then
       FractionFieldAligned = 0.5*(1.0+tanh((Bnorm-Bmodify)/DeltaBmodify))

       HeatCond_D = HeatCoef*( &
            FractionFieldAligned*sum(Bunit_D*Normal_D)*Bunit_D &
            + (1.0 - FractionFieldAligned)*Normal_D )
    else
       HeatCond_D = HeatCoef*sum(Bunit_D*Normal_D)*Bunit_D
    end if

  end subroutine get_heat_conduction_coef

  !============================================================================
  ! Operator split, semi-implicit subroutines
  !============================================================================

  subroutine get_impl_heat_cond_state(StateImpl_VGB, DconsDsemi_VCB)

    use ModAdvance,    ONLY: State_VGB, UseIdealEos, &
       LeftState_VX,  LeftState_VY,  LeftState_VZ,  &
       RightState_VX, RightState_VY, RightState_VZ
    use ModFaceValue,  ONLY: calc_face_value
    use ModGeometry,   ONLY: UseCovariant
    use ModImplicit,   ONLY: nw, nImplBLK, impl2iBlk, iTeImpl
    use ModMain,       ONLY: nI, nJ, nK, MaxImplBlk, x_, y_, z_
    use ModPhysics,    ONLY: Si2No_V, UnitTemperature_, UnitEnergyDens_,inv_gm1
    use ModUser,       ONLY: user_material_properties
    use ModVarIndexes, ONLY: Rho_, p_

    real, intent(out) :: StateImpl_VGB(nw,0:nI+1,0:nJ+1,0:nK+1,MaxImplBlk)
    real, intent(inout) :: DconsDsemi_VCB(nw,nI,nJ,nK,MaxImplBlk)

    integer :: i, j, k, iBlock, iImplBlock
    real :: TeSi, CvSi
    real :: Normal_D(3), Area
    real :: HeatCondL_D(3), HeatCondR_D(3)
    !--------------------------------------------------------------------------

    do iImplBlock = 1, nImplBLK
       iBlock = impl2iBLK(iImplBlock)

       if(UseIdealEos)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             StateImpl_VGB(iTeImpl,i,j,k,iImplBlock) = &
                  TeFraction*State_VGB(p_,i,j,k,iBlock) &
                  /State_VGB(Rho_,i,j,k,iBlock)
             DconsDsemi_VCB(iTeImpl,i,j,k,iImplBlock) = &
                  inv_gm1*State_VGB(Rho_,i,j,k,iBlock)/TeFraction
          end do; end do; end do
       else
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             call user_material_properties( &
                  State_VGB(:,i,j,k,iBlock), TeSiOut=TeSi, CvSiOut = CvSi)
             StateImpl_VGB(iTeImpl,i,j,k,iImplBlock) = &
                  TeSi*Si2No_V(UnitTemperature_)
             DconsDsemi_VCB(iTeImpl,i,j,k,iImplBlock) = &
                  CvSi*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_)
          end do; end do; end do
       end if

       call calc_face_value(.false.,iBlock)

       if(.not.UseCovariant) call set_cartesian_cell_face(x_, iBlock)
       do k = 1, nK; do j = 1, nJ; do i = 1, nI+1
          if(UseCovariant) call set_covariant_cell_face(x_, i, j, k, iBlock)

          call get_heat_conduction_coef(x_, i, j, k, iBlock, &
               LeftState_VX(:,i,j,k), Normal_D, HeatCondL_D)
          call get_heat_conduction_coef(x_, i, j, k, iBlock, &
               RightState_VX(:,i,j,k), Normal_D, HeatCondR_D)

          HeatCond_DFDB(:,i,j,k,1,iBlock) = 0.5*(HeatCondL_D+HeatCondR_D)*Area
       end do; end do; end do

       if(.not.UseCovariant) call set_cartesian_cell_face(y_, iBlock)
       do k = 1, nK; do j = 1, nJ+1; do i = 1, nI
          if(UseCovariant) call set_covariant_cell_face(y_, i, j, k, iBlock)

          call get_heat_conduction_coef(y_, i, j, k, iBlock, &
               LeftState_VY(:,i,j,k), Normal_D, HeatCondL_D)
          call get_heat_conduction_coef(y_, i, j, k, iBlock, &
               RightState_VY(:,i,j,k), Normal_D, HeatCondR_D)

          HeatCond_DFDB(:,i,j,k,2,iBlock) = 0.5*(HeatCondL_D+HeatCondR_D)*Area
       end do; end do; end do

       if(.not.UseCovariant) call set_cartesian_cell_face(z_, iBlock)
       do k = 1, nK+1; do j = 1, nJ; do i = 1, nI
          if(UseCovariant) call set_covariant_cell_face(z_, i, j, k, iBlock)

          call get_heat_conduction_coef(z_, i, j, k, iBlock, &
               LeftState_VZ(:,i,j,k), Normal_D, HeatCondL_D)
          call get_heat_conduction_coef(z_, i, j, k, iBlock, &
               RightState_VZ(:,i,j,k), Normal_D, HeatCondR_D)

          HeatCond_DFDB(:,i,j,k,3,iBlock) = 0.5*(HeatCondL_D+HeatCondR_D)*Area
       end do; end do; end do

    end do

  contains

    subroutine set_cartesian_cell_face(iDim, iBlock)

      use ModGeometry, ONLY: fAx_BLK, fAy_BLK, fAz_BLK

      integer, intent(in) :: iDim, iBlock
      !------------------------------------------------------------------------

      Normal_D = 0.0; Normal_D(iDim) = 1.0
      select case(iDim)
      case(x_)
         Area = fAx_BLK(iBlock)
      case(y_)
         Area = fAy_BLK(iBlock)
      case(z_)
         Area = fAz_BLK(iBlock)
      end select

    end subroutine set_cartesian_cell_face

    !==========================================================================

    subroutine set_covariant_cell_face(iDim, i, j, k, iBlock)

      use ModGeometry, ONLY: FaceAreaI_DFB, FaceAreaJ_DFB, FaceAreaK_DFB, &
           FaceArea2MinI_B, FaceArea2MinJ_B, FaceArea2MinK_B, &
           x_BLK, y_BLK, z_BLK

      integer, intent(in) :: iDim, i, j, k, iBlock

      integer :: iLeft, jLeft, kLeft
      real :: Area_D(3), Area2, Area2Min
      !------------------------------------------------------------------------

      select case(iDim)
      case(x_)
         iLeft = i-1; jLeft = j; kLeft = k
         Area_D = FaceAreaI_DFB(:,i,j,k,iBlock)
         Area2Min = FaceArea2MinI_B(iBlock)
      case(y_)
         iLeft = i; jLeft = j-1; kLeft = k
         Area_D = FaceAreaJ_DFB(:,i,j,k,iBlock)
         Area2Min = FaceArea2MinJ_B(iBlock)
      case(z_)
         iLeft = i; jLeft = j; kLeft = k-1
         Area_D = FaceAreaK_DFB(:,i,j,k,iBlock)
         Area2Min = FaceArea2MinK_B(iBlock)
      end select

      Area2 = sum(Area_D**2)

      if(Area2 < 0.5*Area2Min)then
         ! The face is at the pole
         Normal_D(x_) = x_BLK(i,j,k,iBlock) - x_BLK(iLeft,jLeft,kLeft,iBlock)
         Normal_D(y_) = y_BLK(i,j,k,iBlock) - y_BLK(iLeft,jLeft,kLeft,iBlock)
         Normal_D(z_) = z_BLK(i,j,k,iBlock) - z_BLK(iLeft,jLeft,kLeft,iBlock)
         Normal_D = Normal_D/sqrt(sum(Normal_D**2))
         Area = 0.0
      else
         Area = sqrt(Area2)
         Normal_D = Area_D/Area
      end if

    end subroutine set_covariant_cell_face

  end subroutine get_impl_heat_cond_state

  !============================================================================

  subroutine get_heat_conduction_bc(iBlock, IsLinear)

    use ModImplicit, ONLY: StateSemi_VGB, nw
    use ModMain,     ONLY: nI, nJ, nK, TypeBc_I
    use ModParallel, ONLY: NOBLK, NeiLev
    use ModUser,     ONLY: user_set_outerbcs

    integer, intent(in) :: iBlock
    logical, intent(in) :: IsLinear

    logical :: IsFound
    character(len=20), parameter :: TypeUserBc = 'usersemi'
    character(len=20), parameter :: TypeUserBcLinear = 'usersemilinear'
    character(len=*),  parameter :: NameSub = 'get_heat_conduction_bc'
    !--------------------------------------------------------------------------

    if(NeiLev(1,iBlock) == NOBLK)then
       if(TypeBc_I(1) == 'outflow' .or. TypeBc_I(1) == 'float')then
          if(IsLinear)then
             StateSemi_VGB(:,0,:,:,iBlock) = 0.0
          else
             StateSemi_VGB(:,0,:,:,iBlock) = StateSemi_VGB(:,1,:,:,iBlock)
          end if
       elseif(TypeBc_I(1) == 'user')then
          if(IsLinear)then
             StateSemi_VGB(:,0,:,:,iBlock) = 0.0
             call user_set_outerbcs(iBlock,1,TypeUserBcLinear,IsFound)
          else
             IsFound = .false.
             call user_set_outerbcs(iBlock,1,TypeUserBc,IsFound)
             if(.not. IsFound) call stop_mpi(NameSub//': unknown TypeBc=' &
                  //TypeUserBc//' on iSide=1 in user_set_outerbcs')
          end if
       elseif(TypeBc_I(1) == 'reflect')then
          StateSemi_VGB(:,0,:,:,iBlock) = StateSemi_VGB(:,1,:,:,iBlock)
       else
          call stop_mpi(NameSub//': unknown TypeBc_I(1)='//TypeBc_I(1))
       end if
    end if
    if(NeiLev(2,iBlock) == NOBLK)then
       if(TypeBc_I(2) == 'outflow' .or. TypeBc_I(2) == 'float')then
          if(IsLinear)then
             StateSemi_VGB(:,nI+1,:,:,iBlock) = 0.0
          else
             StateSemi_VGB(:,nI+1,:,:,iBlock) = StateSemi_VGB(:,nI,:,:,iBlock)
          end if
       elseif(TypeBc_I(2) == 'user')then
          if(IsLinear)then
             StateSemi_VGB(:,nI+1,:,:,iBlock) = 0.0
             call user_set_outerbcs(iBlock,2,TypeUserBcLinear,IsFound)
          else
             IsFound = .false.
             call user_set_outerbcs(iBlock,2,TypeUserBc,IsFound)
             if(.not. IsFound) call stop_mpi(NameSub//': unknown TypeBc=' &
                  //TypeUserBc//' on iSide=2 in user_set_outerbcs')
          end if
       elseif(TypeBc_I(2) == 'reflect')then
          StateSemi_VGB(:,nI+1,:,:,iBlock) = StateSemi_VGB(:,nI,:,:,iBlock)
       else
          call stop_mpi(NameSub//': unknown TypeBc_I(2)='//TypeBc_I(2))
       end if
    end if
    if(NeiLev(3,iBlock) == NOBLK)then
       if(TypeBc_I(3) == 'outflow' .or. TypeBc_I(3) == 'float')then
          if(IsLinear)then
             StateSemi_VGB(:,:,0,:,iBlock) = 0.0
          else
             StateSemi_VGB(:,:,0,:,iBlock) = StateSemi_VGB(:,:,1,:,iBlock)
          end if
       elseif(TypeBc_I(3) == 'user')then
          if(IsLinear)then
             StateSemi_VGB(:,:,0,:,iBlock) =  0.0
             call user_set_outerbcs(iBlock,3,TypeUserBcLinear,IsFound)
          else
             IsFound = .false.
             call user_set_outerbcs(iBlock,3,TypeUserBc,IsFound)
             if(.not. IsFound) call stop_mpi(NameSub//': unknown TypeBc=' &
                  //TypeUserBc//' on iSide=3 in user_set_outerbcs')
          end if
       elseif(TypeBc_I(3) == 'reflect')then
          StateSemi_VGB(:,:,0,:,iBlock) = StateSemi_VGB(:,:,1,:,iBlock)
       else
          call stop_mpi(NameSub//': unknown TypeBc_I(3)='//TypeBc_I(3))
       end if
    end if
    if(NeiLev(4,iBlock) == NOBLK) then
       if(TypeBc_I(4) == 'outflow' .or. TypeBc_I(4) == 'float')then
          if(IsLinear)then
             StateSemi_VGB(:,:,nJ+1,:,iBlock) = 0.0
          else
             StateSemi_VGB(:,:,nJ+1,:,iBlock) = StateSemi_VGB(:,:,nJ,:,iBlock)
          end if
       elseif(TypeBc_I(4) == 'user')then
          if(IsLinear)then
             StateSemi_VGB(:,:,nJ+1,:,iBlock) = 0.0
             call user_set_outerbcs(iBlock,4,TypeUserBcLinear,IsFound)
          else
             IsFound = .false.
             call user_set_outerbcs(iBlock,4,TypeUserBc,IsFound)
             if(.not. IsFound) call stop_mpi(NameSub//': unknown TypeBc=' &
                  //TypeUserBc//' on iSide=4 in user_set_outerbcs')
          end if
       elseif(TypeBc_I(4) == 'reflect')then
          StateSemi_VGB(:,:,nJ+1,:,iBlock) = StateSemi_VGB(:,:,nJ,:,iBlock)
       else
          call stop_mpi(NameSub//': unknown TypeBc_I(4)='//TypeBc_I(4))
       end if
    end if
    if(NeiLev(5,iBlock) == NOBLK) then
       if(TypeBc_I(5) == 'outflow' .or. TypeBc_I(5) == 'float')then
          if(IsLinear)then
             StateSemi_VGB(:,:,:,0,iBlock) = 0.0
          else
             StateSemi_VGB(:,:,:,0,iBlock) = StateSemi_VGB(:,:,:,1,iBlock)
          end if
       elseif(TypeBc_I(5) == 'user')then
          if(IsLinear)then
             StateSemi_VGB(:,:,:,0,iBlock) = 0.0
             call user_set_outerbcs(iBlock,5,TypeUserBcLinear,IsFound)
          else
             IsFound = .false.
             call user_set_outerbcs(iBlock,5,TypeUserBc,IsFound)
             if(.not. IsFound) call stop_mpi(NameSub//': unknown TypeBc=' &
                  //TypeUserBc//' on iSide=5 in user_set_outerbcs')
          end if
       elseif(TypeBc_I(5) == 'reflect')then
          StateSemi_VGB(:,:,:,0,iBlock) = StateSemi_VGB(:,:,:,1,iBlock)
       else
          call stop_mpi(NameSub//': unknown TypeBc_I(5)='//TypeBc_I(5))
       end if
    end if
    if(NeiLev(6,iBlock) == NOBLK)then 
       if(TypeBc_I(6) == 'outflow' .or. TypeBc_I(6) == 'float')then
          if(IsLinear)then
             StateSemi_VGB(:,:,:,nK+1,iBlock) = 0.0
          else
             StateSemi_VGB(:,:,:,nK+1,iBlock) = StateSemi_VGB(:,:,:,nK,iBlock)
          end if
       elseif(TypeBc_I(6) == 'user')then
          if(IsLinear)then
             StateSemi_VGB(:,:,:,nK+1,iBlock) = 0.0
             call user_set_outerbcs(iBlock,6,TypeUserBcLinear,IsFound)
          else
             IsFound = .false.
             call user_set_outerbcs(iBlock,6,TypeUserBc,IsFound)
             if(.not. IsFound) call stop_mpi(NameSub//': unknown TypeBc=' &
                  //TypeUserBc//' on iSide=6 in user_set_outerbcs')
          end if
       elseif(TypeBc_I(6) == 'reflect')then
          StateSemi_VGB(:,:,:,nK+1,iBlock) = StateSemi_VGB(:,:,:,nK,iBlock)
       else
          call stop_mpi(NameSub//': unknown TypeBc_I(6)='//TypeBc_I(6))
       end if
    end if

  end subroutine get_heat_conduction_bc

  !============================================================================

  subroutine get_heat_conduction_rhs(iBlock, StateImpl_VG, Rhs_VC, IsLinear)

    use ModFaceGradient, ONLY: get_face_gradient
    use ModGeometry,     ONLY: vInv_CB
    use ModImplicit,     ONLY: nw, iTeImpl
    use ModMain,         ONLY: nI, nJ, nK, x_, y_, z_

    integer, intent(in) :: iBlock
    real, intent(inout) :: StateImpl_VG(nw,-1:nI+2,-1:nJ+2,-1:nK+2)
    real, intent(out)   :: Rhs_VC(nw,nI,nJ,nK)
    logical, intent(in) :: IsLinear

    integer :: i, j, k
    real :: FaceGrad_D(3)
    logical :: IsNewBlockHeatConduction
    !--------------------------------------------------------------------------

    IsNewBlockHeatConduction = .true.

    do k = 1, nK; do j = 1, nJ; do i = 1, nI+1
       call get_face_gradient(x_, i, j, k, iBlock, &
            IsNewBlockHeatConduction, StateImpl_VG, FaceGrad_D)
       FluxImpl_X(i,j,k) = -sum(HeatCond_DFDB(:,i,j,k,1,iBlock)*FaceGrad_D)
    end do; end do; end do
    do k = 1, nK; do j = 1, nJ+1; do i = 1, nI
       call get_face_gradient(y_, i, j, k, iBlock, &
            IsNewBlockHeatConduction, StateImpl_VG, FaceGrad_D)
       FluxImpl_Y(i,j,k) = -sum(HeatCond_DFDB(:,i,j,k,2,iBlock)*FaceGrad_D)
    end do; end do; end do
    do k = 1, nK+1; do j = 1, nJ; do i = 1, nI
       call get_face_gradient(z_, i, j, k, iBlock, &
            IsNewBlockHeatConduction, StateImpl_VG, FaceGrad_D)
       FluxImpl_Z(i,j,k) = -sum(HeatCond_DFDB(:,i,j,k,3,iBlock)*FaceGrad_D)
    end do; end do; end do

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       Rhs_VC(:,i,j,k) = vInv_CB(i,j,k,iBlock) &
            *(FluxImpl_X(i,j,k) - FluxImpl_X(i+1,j,k) &
            + FluxImpl_Y(i,j,k) - FluxImpl_Y(i,j+1,k) &
            + FluxImpl_Z(i,j,k) - FluxImpl_Z(i,j,k+1) )
    end do; end do; end do

  end subroutine get_heat_conduction_rhs

  !============================================================================

  subroutine get_heat_cond_jacobian(iBlock, nVar, Jacobian_VVCI)

    use ModGeometry, ONLY: dx_BLK, dy_BLK, dz_BLK, vInv_CB, UseCovariant
    use ModImplicit, ONLY: iTeImpl
    use ModMain,     ONLY: nI, nJ, nK, nDim
    use ModNumConst, ONLY: i_DD

    integer, parameter:: nStencil = 2*nDim + 1

    integer, intent(in) :: iBlock, nVar
    real, intent(out) :: Jacobian_VVCI(nVar,nVar,nI,nJ,nK,nStencil)

    integer :: i, j, k, iDim, Di, Dj, Dk
    real :: DiffLeft, DiffRight, Dxyz_D(nDim)
    !--------------------------------------------------------------------------

    ! All elements have to be set
    Jacobian_VVCI(:,:,:,:,:,:) = 0.0

    Dxyz_D = (/dx_BLK(iBlock), dy_BLK(iBlock), dz_Blk(iBlock)/)

    ! the transverse diffusion is ignored in the Jacobian

    if(UseCovariant)then
       call get_jacobian_covariant
    else

       do iDim = 1, nDim
          Di = i_DD(iDim,1); Dj = i_DD(iDim,2); Dk = i_DD(iDim,3)
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             DiffLeft = vInv_CB(i,j,k,iBlock) &
                  *HeatCond_DFDB(iDim,i,j,k,iDim,iBlock)/Dxyz_D(iDim)
             DiffRight = vInv_CB(i,j,k,iBlock) &
                  *HeatCond_DFDB(iDim,i+Di,j+Dj,k+Dk,iDim,iBlock)/Dxyz_D(iDim)

             Jacobian_VVCI(iTeImpl,iTeImpl,i,j,k,1) = &
                  Jacobian_VVCI(iTeImpl,iTeImpl,i,j,k,1) &
                  - (DiffLeft + DiffRight)

             if(iDim==1.and.i==1 .or. iDim==2.and.j==1 .or. iDim==3.and.k==1) &
                  DiffLeft = 0.0
             if(iDim==1.and.i==nI .or. iDim==2.and.j==nJ &
                  .or. iDim==3.and.k==nK) &
                  DiffRight = 0.0

             Jacobian_VVCI(iTeImpl,iTeImpl,i,j,k,2*iDim)   = DiffLeft
             Jacobian_VVCI(iTeImpl,iTeImpl,i,j,k,2*iDim+1) = DiffRight
          end do; end do; end do
       end do
    end if

  contains

    subroutine get_jacobian_covariant

      use ModFaceGradient, ONLY: set_block_jacobian_face, DcoordDxyz_DDFD

      real :: InvDxyz_D(nDim)
      !------------------------------------------------------------------------

      call set_block_jacobian_face(iBlock)

      do iDim = 1, nDim
         Di = i_DD(iDim,1); Dj = i_DD(iDim,2); Dk = i_DD(iDim,3)
         do k = 1, nK; do j = 1, nJ; do i = 1, nI
            InvDxyz_D = DcoordDxyz_DDFD(iDim,:,i,j,k,iDim)/Dxyz_D(iDim)
            DiffLeft = vInv_CB(i,j,k,iBlock) &
                 *sum(HeatCond_DFDB(:,i,j,k,iDim,iBlock)*InvDxyz_D)

            InvDxyz_D = &
                 DcoordDxyz_DDFD(iDim,:,i+Di,j+Dj,k+Dk,iDim)/Dxyz_D(iDim)
            DiffRight = vInv_CB(i,j,k,iBlock) &
                 *sum(HeatCond_DFDB(:,i+Di,j+Dj,k+Dk,iDim,iBlock)*InvDxyz_D)

            Jacobian_VVCI(iTeImpl,iTeImpl,i,j,k,1) = &
                 Jacobian_VVCI(iTeImpl,iTeImpl,i,j,k,1) &
                 - (DiffLeft + DiffRight)

            if(iDim==1.and.i==1 .or. iDim==2.and.j==1 .or. iDim==3.and.k==1) &
                 DiffLeft = 0.0
            if(iDim==1.and.i==nI .or. iDim==2.and.j==nJ &
                 .or. iDim==3.and.k==nK) DiffRight = 0.0

            Jacobian_VVCI(iTeImpl,iTeImpl,i,j,k,2*iDim)   = DiffLeft
            Jacobian_VVCI(iTeImpl,iTeImpl,i,j,k,2*iDim+1) = DiffRight
         end do; end do; end do
      end do

    end subroutine get_jacobian_covariant
    
  end subroutine get_heat_cond_jacobian

  !============================================================================

  subroutine update_impl_heat_cond(iBlock, iImplBlock, StateImpl_VG)

    use ModAdvance,  ONLY: State_VGB, Energy_GBI, UseIdealEos, p_, &
         ExtraEint_
    use ModEnergy,   ONLY: calc_energy_cell, calc_pressure_cell
    use ModImplicit, ONLY: nw, iTeImpl, DconsDsemi_VCB, ImplOld_VCB
    use ModMain,     ONLY: nI, nJ, nK
    use ModPhysics,  ONLY: inv_gm1, No2Si_V, Si2No_V, UnitEnergyDens_, UnitP_
    use ModUser,     ONLY: user_material_properties

    integer, intent(in) :: iBlock, iImplBlock
    real, intent(in) :: StateImpl_VG(nw,nI,nJ,nK)

    integer :: i, j, k
    real :: Einternal, EinternalSi, PressureSi
    !--------------------------------------------------------------------------

    if(UseIdealEos)then
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          Energy_GBI(i,j,k,iBlock,1) = Energy_GBI(i,j,k,iBlock,1) &
               + DconsDsemi_VCB(iTeImpl,i,j,k,iImplBlock) &
               *( StateImpl_VG(iTeImpl,i,j,k) &
               -  ImplOld_VCB(iTeImpl,i,j,k,iBlock) )
       end do; end do; end do

       call calc_pressure_cell(iBlock)
    else
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          Einternal = inv_gm1*State_VGB(p_,i,j,k,iBlock) &
               + State_VGB(ExtraEint_,i,j,k,iBlock) &
               + DconsDsemi_VCB(iTeImpl,i,j,k,iImplBlock) &
               *( StateImpl_VG(iTeImpl,i,j,k) &
               -  ImplOld_VCB(iTeImpl,i,j,k,iBlock) )

          EinternalSi = Einternal*No2Si_V(UnitEnergyDens_)

          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               EinternalSiIn = EinternalSi, PressureSiOut = PressureSi)

          State_VGB(p_,i,j,k,iBlock) = PressureSi*Si2No_V(UnitP_)
          State_VGB(ExtraEint_,i,j,k,iBlock) = &
               Einternal - inv_gm1*State_VGB(p_,i,j,k,iBlock)

       end do; end do; end do

       call calc_energy_cell(iBlock)
    end if

  end subroutine update_impl_heat_cond

end module ModHeatConduction
