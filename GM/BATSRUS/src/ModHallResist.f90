!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModHallResist

  use ModSize, ONLY: nI, nJ, nK, MaxDim, j0_, nJp1_, k0_, nKp1_

  implicit none

  SAVE

  private !except

  ! Public methods
  public :: init_hall_resist
  public :: read_hall_param
  public :: set_hall_factor_face
  public :: set_hall_factor_cell
  public :: set_ion_mass_per_charge
  public :: set_ion_mass_per_charge_point

  ! Logical for adding the Biermann battery term
  logical, public :: UseBiermannBattery = .false.

  ! Logical for adding hall resistivity
  logical, public:: UseHallResist=.false.
  logical, public:: IsNewBlockCurrent=.true.

  ! Coefficient for taking whistler wave speed into account
  real, public:: HallCmaxFactor = 1.0

  ! Adjustable coefficient for the Hall term
  ! (similar effect as changing the ion mass per charge)
  real, public:: HallFactorMax = 1.0

  ! Ion mass per charge may depend on space and time for multispecies
  real, public, allocatable:: IonMassPerCharge_G(:,:,:)
  real:: IonMassPerChargeCoef

  ! Arrays for the implicit preconditioning
  real, public, allocatable :: HallJ_CD(:,:,:,:)

  ! Hall factor on the faces and in the cell centers
  real, public, allocatable:: HallFactor_DF(:,:,:,:), HallFactor_C(:,:,:)

  ! Logical is true if call set_hall_factor* sets any non-zero hall factors
  logical, public:: IsHallBlock

  ! Local variables ---------

  ! Description of the region where Hall effect is used
  character(len=200):: StringHallRegion ='none'

  ! Indexes of regions defined with the #REGION commands
  integer, allocatable:: iRegionHall_I(:)

  ! If Hall speed is less than HallSpeedTinyDim, ignore the Hall effect
  real:: HallSpeedTinyDim = -1.0

contains
  !============================================================================
  subroutine init_hall_resist

    use ModConst,   ONLY: cElectronCharge
    use ModPhysics, ONLY: IonMassPerCharge, Si2No_V, UnitX_, UnitCharge_
    use BATL_lib,   ONLY: get_region_indexes

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='init_hall_resist'
    !--------------------------------------------------------------------------

    call set_oktest(NameSub, DoTest, DoTestMe)

    if (DoTestMe) then
       write(*,*) ''
       write(*,*) '>>>>>>>>>>>>>>>>> HALL Resistivity Parameters <<<<<<<<<<'
       write(*,*)
       write(*,*) 'HallFactorMax    = ', HallFactorMax
       write(*,*) 'HallCmaxFactor   = ', HallCmaxFactor
       write(*,*) 'IonMassPerCharge = ', IonMassPerCharge
       ! Omega_Bi=B0/IonMassPerCharge'
       write(*,*)
       write(*,*) '>>>>>>>>>>>>>>>>>                       <<<<<<<<<<<<<<<<<'
       write(*,*) ''
    end if

    if(.not.allocated(HallJ_CD)) allocate(              &
         HallJ_CD(nI,nJ,nK,MaxDim),                     &
         IonMassPerCharge_G(0:nI+1,j0_:nJp1_,k0_:nKp1_) )

    HallJ_CD = 0.0

    IonMassPerCharge_G = IonMassPerCharge

    ! This is used in combination with normalized density
    ! divided by SI charge density.
    IonMassPerChargeCoef = &
         Si2No_V(UnitX_)**3 / (cElectronCharge*Si2No_V(UnitCharge_))

    ! Get signed indexes for Hall region(s)
    call get_region_indexes(StringHallRegion, iRegionHall_I)

  end subroutine init_hall_resist

  !=========================================================================
  subroutine read_hall_param(NameCommand)

    use ModReadParam, ONLY: read_var
    use ModSemiImplVar, ONLY: IsDynamicSemiImpl

    character(len=*), intent(in):: NameCommand

    character(len=*), parameter:: NameSub = 'read_hall_param'
    !---------------------------------------------------------------------
    select case(NameCommand)
    case("#HALLRESISTIVITY")
       call read_var('UseHallResist',  UseHallResist)
       call read_var('HallFactorMax',  HallFactorMax)
       call read_var('HallCmaxFactor', HallCmaxFactor)

    case("#HALLREGION")
       call read_var('StringHallRegion', StringHallRegion)
       call read_var('HallSpeedTinyDim', HallSpeedTinyDim)
       IsDynamicSemiImpl = HallSpeedTinyDim > 0

    case("#BIERMANNBATTERY")
       call read_var("UseBiermannBattery", UseBiermannBattery)

    case default
       call stop_mpi(NameSub//' unknown command='//NameCommand)
    end select

  end subroutine read_hall_param

  !=========================================================================
  subroutine set_ion_mass_per_charge(iBlock)

    use ModAdvance, ONLY: State_VGB, UseIdealEos
    use ModVarIndexes, ONLY: UseMultiSpecies
    use ModMultiFluid, ONLY: UseMultiIon

    ! Set IonMassPerCharge_G based on average mass
    integer, intent(in) :: iBlock

    integer :: i, j, k
    !-------------------------------------------------------------------------

    ! Check if IonMassPerCharge_G varies at all. Return if it is constant.
    if(.not.UseMultiIon .and. .not.UseMultiSpecies .and. UseIdealEos) RETURN

    ! Set IonMassPerCharge_G to the average ion mass = rho_total / n_total
    ! including 1 layer of ghost cells
    do k = k0_, nKp1_; do j = j0_, nJp1_; do i = 0, nI+1
       call set_ion_mass_per_charge_point(State_VGB(:,i,j,k,iBlock), &
            IonMassPerCharge_G(i,j,k))
    end do; end do; end do

  end subroutine set_ion_mass_per_charge

  !===========================================================================

  subroutine set_ion_mass_per_charge_point(State_V, IonMassPerChargeOut)

    use ModAdvance,    ONLY: UseIdealEos
    use ModVarIndexes, ONLY: nVar, Rho_, &
         UseMultiSpecies, SpeciesFirst_, SpeciesLast_, MassSpecies_V
    use ModMultiFluid, ONLY: UseMultiIon, iRhoIon_I, MassIon_I,ChargeIon_I
    use ModPhysics,    ONLY: IonMassPerCharge
    use ModUserInterface

    real, intent(in) :: State_V(nVar)
    real, intent(out):: IonMassPerChargeOut

    real :: zAverage, NatomicSi
    !--------------------------------------------------------------------------

    if(.not.UseIdealEos)then
       call user_material_properties(State_V, &
            AverageIonChargeOut=zAverage, NatomicOut=NatomicSi)

       ! Avoid using small zAverage, since then we will generate magnetic
       ! field with the Biermann Battery term based numerical errors.
       zAverage = max(zAverage, 1.0)

       IonMassPerChargeOut = IonMassPerChargeCoef*State_V(Rho_) &
            /(zAverage*NatomicSi)

    elseif(UseMultiSpecies)then
       IonMassPerChargeOut = IonMassPerCharge*State_V(Rho_) &
            / sum(State_V(SpeciesFirst_:SpeciesLast_)/MassSpecies_V)

    elseif(UseMultiIon)then
       ! Get mass density per total number denisity
       IonMassPerChargeOut = IonMassPerCharge*sum(State_V(iRhoIon_I)) &
            / sum(State_V(iRhoIon_I)*ChargeIon_I / MassIon_I)

    else
       IonMassPerChargeOut = IonMassPerCharge

    end if

  end subroutine set_ion_mass_per_charge_point

  !=========================================================================
  subroutine set_hall_factor_cell(iBlock, UseIonMassPerCharge)

    use BATL_lib, ONLY: block_inside_regions
    use ModCurrent, ONLY: get_current
    use ModAdvance, ONLY: State_VGB
    use ModMultiFluid, ONLY: iRhoIon_I
    use ModPhysics, ONLY: Io2No_V, UnitU_

    integer, intent(in):: iBlock
    logical, optional, intent(in):: UseIonMassPerCharge

    ! Set the hall factor for the cell centers of block iBlock
    ! Also set IsHallBlock if any of the cells have a non-zero factor
    ! If UseIonMassPerCharge is present and FALSE, the hall factor is
    ! NOT multiplied with the ion mass per charge averaged to the face

    integer:: i, j, k
    real:: Current_D(MaxDim)
    real:: Rho, HallSpeedTiny2, HallSpeed2
    !----------------------------------------------------------------------
    if(.not.allocated(HallFactor_C)) allocate(HallFactor_C(nI,nJ,nK))

    if(.not.allocated(iRegionHall_I))then
       IsHallBlock = .true.
       HallFactor_C = HallFactorMax
    else
       call block_inside_regions(iRegionHall_I, iBlock, &
            size(HallFactor_C), 'cells', IsHallBlock, Value_I=HallFactor_C)

       ! Nothing to do if the block does not intersect with the Hall region
       if(.not.IsHallBlock) RETURN

       ! Multiply by HallFactorMax
       if(HallFactorMax /= 1) HallFactor_C = HallFactorMax*HallFactor_C
    end if

    if(present(UseIonMassPerCharge))then
       if(.not.UseIonMassPerCharge) RETURN
    end if

    ! Multiply with ion mass per charge
    call set_ion_mass_per_charge(iBlock)
    HallFactor_C = HallFactor_C*IonMassPerCharge_G(1:nI,1:nJ,1:nK)

    if(HallSpeedTinyDim <= 0) RETURN

    HallSpeedTiny2 = (HallSpeedTinyDim*Io2No_V(UnitU_))**2
    do k=1, nK; do j=1,nJ; do i=1, nI
       call get_current(i, j, k, iBlock, Current_D)
       Rho = sum(State_VGB(iRhoIon_I,i,j,k,iBlock))
       HallSpeed2 = sum(Current_D**2)*(HallFactor_C(i,j,k)/Rho)**2
       if(HallSpeed2 < HallSpeedTiny2) HallFactor_C(i,j,k) = 0
    end do; end do; end do

  end subroutine set_hall_factor_cell
  !=========================================================================
  subroutine set_hall_factor_face(iBlock, UseIonMassPerCharge)

    use BATL_lib, ONLY: block_inside_regions, nDim, nINode, nJNode, nKNode
    use ModAdvance, ONLY: State_VGB
    use ModMultiFluid, ONLY: iRhoIon_I
    use ModPhysics, ONLY: Io2No_V, UnitU_
    use ModCurrent, ONLY: get_current

    integer, intent(in):: iBlock
    logical, optional, intent(in):: UseIonMassPerCharge

    ! Set the hall factor for the cell faces of block iBlock
    ! Also set IsHallBlock if any of the faces have a non-zero factor
    ! If UseIonMassPerCharge is present and FALSE, the hall factor is
    ! NOT multiplied with the ion mass per charge averaged to the face

    integer:: i, j, k

    real, allocatable, save:: j2_G(:,:,:)
    real:: Current_D(MaxDim)
    real:: Rho, HallSpeedTiny2, HallSpeed2
    logical:: IsNewBlock
    !----------------------------------------------------------------------
    if(.not.allocated(HallFactor_DF)) &
         allocate(HallFactor_DF(nDim,nINode,nJNode,nKNode))

    if(.not.allocated(iRegionHall_I))then
       IsHallBlock = .true.
       HallFactor_DF = HallFactorMax
    else
       call block_inside_regions(iRegionHall_I, iBlock, &
            size(HallFactor_DF), 'face', IsHallBlock, Value_I=HallFactor_DF)

       ! Nothing to do if the block does not intersect with the Hall region
       if(.not.IsHallBlock) RETURN

       ! Multiply by HallFactorMax
       if(HallFactorMax /= 1) HallFactor_DF = HallFactorMax*HallFactor_DF
    end if

    if(present(UseIonMassPerCharge))then
       if(.not.UseIonMassPerCharge) RETURN
    end if

    ! Multiply with ion mass per charge
    call set_ion_mass_per_charge(iBlock)

    ! Check for small Hall velocity
    if(HallSpeedTinyDim > 0)then
       ! Calculate Hall speed and compare with HallSpeedTiny
       HallSpeedTiny2 = (HallSpeedTinyDim*Io2No_V(UnitU_))**2
       if(.not.allocated(j2_G)) allocate(j2_G(0:nI+1,j0_:nJp1_,k0_:nKp1_))
       do k=k0_, nKp1_; do j=j0_,nJp1_; do i=0, nI+1
          call get_current(i, j, k, iBlock, Current_D)
          j2_G(i,j,k) = sum(Current_D**2)
       end do; end do; end do
    end if

    do k=1, nK; do j=1,nJ; do i=1, nI+1
       HallFactor_DF(1,i,j,k) = HallFactor_DF(1,i,j,k)*&
            0.5*sum(IonMassPerCharge_G(i-1:i,j,k))
       if(HallSpeedTinyDim <= 0) CYCLE
       Rho = 0.5*sum(State_VGB(iRhoIon_I,i-1:i,j,k,iBlock))
       HallSpeed2 = 0.5*sum(j2_G(i-1:i,j,k))*(HallFactor_DF(1,i,j,k)/Rho)**2
       if(HallSpeed2 < HallSpeedTiny2) HallFactor_DF(1,i,j,k) = 0
    end do; end do; end do

    if(nDim == 1) RETURN

    do k=1, nK; do j=1,nJ+1; do i=1, nI
       HallFactor_DF(2,i,j,k) = HallFactor_DF(2,i,j,k)*&
            0.5*sum(IonMassPerCharge_G(i,j-1:j,k))
       if(HallSpeedTinyDim <= 0) CYCLE
       Rho = 0.5*sum(State_VGB(iRhoIon_I,i,j-1:j,k,iBlock))
       HallSpeed2 = 0.5*sum(j2_G(i,j-1:j,k))*(HallFactor_DF(2,i,j,k)/Rho)**2
       if(HallSpeed2 < HallSpeedTiny2) HallFactor_DF(2,i,j,k) = 0
    end do; end do; end do

    if(nDim == 2) RETURN

    do k=1, nK+1; do j=1,nJ; do i=1, nI
       HallFactor_DF(3,i,j,k) = HallFactor_DF(3,i,j,k)*&
            0.5*sum(IonMassPerCharge_G(i,j,k-1:k))
       if(HallSpeedTinyDim <= 0) CYCLE
       Rho = 0.5*sum(State_VGB(iRhoIon_I,i,j,k-1:k,iBlock))
       HallSpeed2 = 0.5*sum(j2_G(i,j,k-1:k))*(HallFactor_DF(3,i,j,k)/Rho)**2
       if(HallSpeed2 < HallSpeedTiny2) HallFactor_DF(3,i,j,k) = 0
    end do; end do; end do

  end subroutine set_hall_factor_face
  !=========================================================================

end module ModHallResist
