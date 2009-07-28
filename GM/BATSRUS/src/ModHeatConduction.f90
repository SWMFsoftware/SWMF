!^CFG COPYRIGHT UM
!==============================================================================
module ModHeatConduction

  implicit none
  save

  private ! except

  ! Public methods
  public :: read_heatconduction_param
  public :: init_heat_conduction
  public :: get_heat_flux
  public :: add_jacobian_heat_conduction

  ! Logical for adding field-aligned heat conduction
  Logical, public :: UseParallelConduction = .false.
  logical, public :: IsNewBlockHeatConduction = .true.

  ! Variables for setting the field-aligned heat conduction coefficient
  character(len=20), public :: TypeHeatConduction = 'test'
  logical :: DoModifyHeatConduction, DoTestHeatConduction
  logical :: UseIdealState = .true.
  
  real :: HeatConductionParSi = 1.23e-11  ! taken from calc_heat_flux
  real :: TmodifySi = 3.0e5, DeltaTmodifySi = 2.0e4
  real :: HeatConductionPar, Tmodify, DeltaTmodify

  logical :: DoWeakFieldConduction = .false.
  real :: BmodifySi = 1.0e-7, DeltaBmodifySi = 1.0e-8 ! modify about 1 mG
  real :: Bmodify, DeltaBmodify

  ! electron temperature used for calculating heat flux
  real, allocatable :: Te_G(:,:,:)

contains

  !============================================================================

  subroutine read_heatconduction_param(NameCommand)

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
          call read_var('UseIdealState', UseIdealState)

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

    use ModPhysics, ONLY: Si2No_V, UnitEnergyDens_, UnitTemperature_, &
         UnitU_, UnitX_, UnitB_
    use ModSize, ONLY: nI, nJ, nK

    character(len=*), parameter :: NameSub = 'init_heat_conduction'
    !--------------------------------------------------------------------------
    
    if(allocated(Te_G)) RETURN

    allocate(Te_G(-1:nI+2,-1:nJ+2,-1:nK+2))

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

  end subroutine init_heat_conduction

  !============================================================================

  subroutine get_heat_flux(iDir, i, j, k, iBlock, State_V, Normal_D, &
       HeatCondCoefNormal, HeatFlux_D)

    use ModAdvance,      ONLY: State_VGB
    use ModB0,           ONLY: B0_DX, B0_DY, B0_DZ
    use ModFaceGradient, ONLY: calc_face_gradient
    use ModMain,         ONLY: nI, nJ, nK, UseB0
    use ModMultiFluid,   ONLY: MassIon_I
    use ModNumConst,     ONLY: cTolerance
    use ModPhysics,      ONLY: inv_gm1, Si2No_V, UnitTemperature_, &
         UnitEnergyDens_, ElectronTemperatureRatio, AverageIonCharge
    use ModUser,         ONLY: user_material_properties
    use ModVarIndexes,   ONLY: nVar, Bx_, Bz_, Rho_, p_

    integer, intent(in) :: iDir, i, j, k, iBlock
    real,    intent(in) :: State_V(nVar), Normal_D(3)
    real,    intent(out):: HeatCondCoefNormal, HeatFlux_D(3)

    integer :: ii, jj, kk
    real :: B_D(3), Bunit_D(3), Bnorm, Cv, CvSi
    real :: FaceGrad_D(3), HeatCoef, TeSi, Te, &
         FractionSpitzer, FractionFieldAligned

    character(len=*), parameter :: NameSub = 'get_heat_flux'
    !--------------------------------------------------------------------------

    if(UseB0)then
       select case(iDir)
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


    if(IsNewBlockHeatConduction)then
       if(UseIdealState)then
          Te_G = State_VGB(p_,:,:,:,iBlock)/State_VGB(Rho_,:,:,:,iBlock) &
               *MassIon_I(1)*ElectronTemperatureRatio &
               /(1 + AverageIonCharge*ElectronTemperatureRatio)
       else
          do kk = -1, nK+2; do jj = -1, nJ+2; do ii = -1, nI+2
             call user_material_properties( &
                  State_VGB(:,ii,jj,kk,iBlock), TeSiOut=TeSi)
             Te_G(ii,jj,kk) = TeSi*Si2No_V(UnitTemperature_)
          end do; end do; end do
       end if
    end if

    call calc_face_gradient(iDir, i, j, k, iBlock, &
         Te_G, IsNewBlockHeatConduction, FaceGrad_D)

    if(UseIdealState)then
       Te = State_V(p_)/State_V(Rho_) &
            *MassIon_I(1)*ElectronTemperatureRatio &
            /(1 + AverageIonCharge*ElectronTemperatureRatio)
       Cv = State_V(Rho_)*inv_gm1
    else
       ! Note we assume that the heat conduction formula for the
       ! ideal state is still applicable for the non-ideal state
       call user_material_properties(State_V, TeSiOut=TeSi, CvSiOut = CvSi)
       Te = TeSi*Si2No_V(UnitTemperature_)
       Cv = CvSi*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_)
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

       HeatFlux_D = -HeatCoef*( &
            FractionFieldAligned*Bunit_D*dot_product(Bunit_D,FaceGrad_D) &
            + (1.0 - FractionFieldAligned)*FaceGrad_D )

       ! get the heat conduction coefficient normal to the face for
       ! time step restriction
       HeatCondCoefNormal = HeatCoef*( &
            FractionFieldAligned*dot_product(Bunit_D,Normal_D)**2 &
            + (1.0 - FractionFieldAligned) )/Cv
    else
       HeatFlux_D = -HeatCoef*Bunit_D*dot_product(Bunit_D,FaceGrad_D)

       ! get the heat conduction coefficient normal to the face for
       ! time step restriction
       HeatCondCoefNormal = HeatCoef*dot_product(Bunit_D,Normal_D)**2/Cv
    end if

  end subroutine get_heat_flux

  !============================================================================

  subroutine add_jacobian_heat_conduction(iBlock, nVar, Jacobian_VVCI)

    use ModMain, ONLY: nI, nJ, nK, nDim

    integer, parameter :: nStencil = 2*nDim + 1

    integer, intent(in) :: iBlock, nVar
    real, intent(inout) :: Jacobian_VVCI(nVar,nVar,nI,nJ,nK,nStencil)
    !--------------------------------------------------------------------------


  end subroutine add_jacobian_heat_conduction

end module ModHeatConduction
