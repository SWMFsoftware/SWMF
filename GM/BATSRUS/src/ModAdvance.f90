!^CFG COPYRIGHT UM
module ModAdvance

  use ModSize
  use ModVarIndexes
  use ModB0
  use ModMain,       ONLY: UseB, UseRotatingFrame, UseGravity
  use ModIO,         ONLY: iUnitOut, write_prefix
  use ModProcMH,     ONLY: iProc, nProc

  implicit none
  save

  ! Numerical flux type
  character (len=10) :: FluxType

  ! Named index for electron pressure and velocity
  integer, parameter:: eFluid_ = nFluid + 1

  ! Additional equations are activated by the named indices of the
  ! corresponding extra variables in the equation modules.
  ! The default values of these named indices are in ModExtraVariables.
  logical, parameter:: UseElectronPressure = Pe_ > 1
  logical, parameter:: UseAnisoPressure    = Ppar_ > 1
  logical, parameter:: UseIdealEos = ExtraEint_ == 1

  logical:: UseWavePressure = .false.

  logical:: DoCalcElectricField = .false.

  !\ One of the two possible ways to treat the MHD-like systems
  !  (oartially symmetrizable, following the Godunov definition).
  !  If the UseRS7=.true. then the 7 waves Riemann Solver (RS) with 
  !  continuous  normal component of the magnetic field across the face.
  !  The number of jumps in the physical variables across the face is equal
  !  to the number of waves, resulting in the the well-posed solution os
  !  the Riemann problem. This approach is alternative to the 8-wave scheme

  logical::UseRS7

  ! Update check parameters
  logical :: UseUpdateCheck
  real :: percent_max_rho(2), percent_max_p(2)

  ! The percentage limit for species to be checked in update check
  real :: SpeciesPercentCheck = 1.0

  ! Replace density with sum of species densities (in multi-species plasma)
  logical :: DoReplaceDensity = .true.

  !\
  ! Conservative/Non-conservative parameters
  !/
  logical :: UseNonConservative

  ! Number and type of criteria
  integer :: nConservCrit
  character (len=10), allocatable :: TypeConservCrit_I(:)

  ! Geometrical parameters
  real    :: rConserv, xParabolaConserv, yParabolaConserv 

  ! Physics based parameters (to locate shocks)
  real    :: pCoeffConserv, GradPCoeffConserv

  ! Cells selected to be updated with conservative equations
  logical, allocatable :: IsConserv_CB(:,:,:,:)

  !\
  ! Block cell-centered MHD solution
  !/
  real, allocatable :: State_VGB(:,:,:,:,:)
  real, allocatable :: Energy_GBI(:,:,:,:,:)

  !\
  ! Block cell-centered MHD solution old state
  !/
  real, allocatable :: StateOld_VCB(:,:,:,:,:)
  real, allocatable :: EnergyOld_CBI(:,:,:,:,:)

  !\
  ! Block cell-centered intrinsic magnetic field, time, and temporary storage
  !/
  real, allocatable :: tmp1_BLK(:,:,:,:)
  real, allocatable :: tmp2_BLK(:,:,:,:)

  real, allocatable :: time_BLK(:,:,:,:)

  ! Arrays for the total electric field
  real, allocatable :: Ex_CB(:,:,:,:)
  real, allocatable :: Ey_CB(:,:,:,:)
  real, allocatable :: Ez_CB(:,:,:,:)

  !\
  ! Block cell-centered body forces
  !/
  real, allocatable :: fbody_x_BLK(:,:,:,:)
  real, allocatable :: fbody_y_BLK(:,:,:,:)
  real, allocatable :: fbody_z_BLK(:,:,:,:)

  !\
  ! Local cell-centered source terms and divB.
  !/
  real :: Source_VC(nVar+nFluid, nI, nJ, nK)
  real :: Theat0(nI,nJ,nK)
  real, allocatable :: DivB1_GB(:,:,:,:)

  real, dimension(0:nI+1, 0:nJ+1, 0:nK+1) :: &
       gradX_Ux, gradX_Uy, gradX_Uz, gradX_Bx, gradX_By, gradX_Bz, gradX_VAR,&
       gradY_Ux, gradY_Uy, gradY_Uz, gradY_Bx, gradY_By, gradY_Bz, gradY_VAR,&
       gradZ_Ux, gradZ_Uy, gradZ_Uz, gradZ_Bx, gradZ_By, gradZ_Bz, gradZ_VAR

  !\
  ! X Face local MHD solution array definitions.
  !/
  ! These are primitive variables (velocity)
  real :: LeftState_VX(nVar,2-gcn:nI+gcn,0:nJ+1,0:nK+1)
  real :: RightState_VX(nVar,2-gcn:nI+gcn,0:nJ+1,0:nK+1)

  real :: EDotFA_X(2-gcn:nI+gcn,0:nJ+1,0:nK+1)    !^CFG IF BORISCORR
  real :: VdtFace_X(2-gcn:nI+gcn,0:nJ+1,0:nK+1)   ! V/dt Face X

  ! Fluxes are for conservative variables (momentum)
  real :: Flux_VX(nVar+nFluid,0:nI+1,2-gcn:nJ+gcn,0:nK+1)

  real :: uDotArea_XI(2-gcn:nI+gcn,0:nJ+1,0:nK+1,nFluid+1)

  real :: bCrossArea_DX(3,2-gcn:nI+gcn,0:nJ+1,0:nK+1)

  !\
  ! Y Face local MHD solution array definitions.
  !/
  real :: LeftState_VY(nVar,0:nI+1,2-gcn:nJ+gcn,0:nK+1)
  real :: RightState_VY(nVar,0:nI+1,2-gcn:nJ+gcn,0:nK+1)

  real :: EDotFA_Y(0:nI+1,2-gcn:nJ+gcn,0:nK+1)    !^CFG IF BORISCORR
  real :: VdtFace_Y(0:nI+1,2-gcn:nJ+gcn,0:nK+1)   ! V/dt Face Y

  real :: Flux_VY(nVar+nFluid,0:nI+1,2-gcn:nJ+gcn,0:nK+1)

  real :: uDotArea_YI(0:nI+1,2-gcn:nJ+gcn,0:nK+1,nFluid+1)

  real :: bCrossArea_DY(3,0:nI+1,2-gcn:nJ+gcn,0:nK+1)

  !\
  ! Z Face local MHD solution array definitions.
  !/
  real :: LeftState_VZ(nVar,0:nI+1,0:nJ+1,2-gcn:nK+gcn)
  real :: RightState_VZ(nVar,0:nI+1,0:nJ+1,2-gcn:nK+gcn)

  real :: EDotFA_Z(0:nI+1,0:nJ+1,2-gcn:nK+gcn)    !^CFG IF BORISCORR
  real :: VdtFace_z(0:nI+1,0:nJ+1,2-gcn:nK+gcn)   ! V/dt Face Z

  real :: Flux_VZ(nVar+nFluid,0:nI+1,0:nJ+1,2-gcn:nK+gcn)

  real :: uDotArea_ZI(0:nI+1,0:nJ+1,2-gcn:nK+gcn,nFluid+1)

  real :: bCrossArea_DZ(3,0:nI+1,0:nJ+1,2-gcn:nK+gcn)

  !\
  ! Merge cells around the polar axis in spherical geometry
  !/
  logical :: DoFixAxis = .false.
  real ::    rFixAxis = 0.0, r2FixAxis = 0.0 

  !\
  ! Block type information
  !/
  integer :: iTypeAdvance_B(MaxBlock)
  integer, allocatable :: iTypeAdvance_BP(:,:)

  ! Named indexes for block types
  integer, parameter :: &
       SkippedBlock_=0,     & ! Blocks which were unused originally.
       SteadyBlock_=1,      & ! Blocks which do not change
       SteadyBoundBlock_=2, & ! Blocks surrounding the evolving blocks
       ExplBlock_=3,        & ! Blocks changing with the explicit scheme
       ImplBlock_=4           ! Blocks changing with the implicit scheme

contains

  !============================================================================

  subroutine init_mod_advance

    ! These arrays may need allocation depending on the parameters
    if(DoCalcElectricField .and. .not.allocated(Ex_CB))then
       allocate(Ex_CB(nI,nJ,nK,MaxBlock))
       allocate(Ey_CB(nI,nJ,nK,MaxBlock))
       allocate(Ez_CB(nI,nJ,nK,MaxBlock))
    end if
    if((UseGravity .or. UseRotatingFrame).and. .not.allocated(fbody_x_BLK))then
       allocate(fbody_x_BLK(nI,nJ,nK,MaxBlock))
       allocate(fbody_y_BLK(nI,nJ,nK,MaxBlock))
       allocate(fbody_z_BLK(nI,nJ,nK,MaxBlock))
    end if

    if(allocated(State_VGB)) RETURN

    ! The arrays below are allocated at the beginning (if at all)
    if(UseB)then 
       allocate(DivB1_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
       DivB1_GB = 0.0
    end if
    allocate(State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    allocate(Energy_GBI(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock,nFluid))
    allocate(StateOld_VCB(nVar,nI,nJ,nK,MaxBlock))
    allocate(EnergyOld_CBI(nI,nJ,nK,MaxBlock,nFluid))
    allocate(tmp1_BLK(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    allocate(tmp2_BLK(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    allocate(time_BLK(nI,nJ,nK,MaxBlock))
    allocate(iTypeAdvance_BP(MaxBlock,0:nProc-1))
    iTypeAdvance_B  = SkippedBlock_
    iTypeAdvance_BP = SkippedBlock_

    if(iProc==0)then
       call write_prefix
       write(iUnitOut,'(a)') 'init_mod_advance allocated arrays'
    end if

  end subroutine init_mod_advance

  !============================================================================

  subroutine clean_mod_advance

    if(allocated(State_VGB))        deallocate(State_VGB)
    if(allocated(Energy_GBI))      deallocate(Energy_GBI)
    if(allocated(StateOld_VCB))    deallocate(StateOld_VCB)
    if(allocated(EnergyOld_CBI))   deallocate(EnergyOld_CBI)
    if(allocated(tmp1_BLK))        deallocate(tmp1_BLK)
    if(allocated(tmp2_BLK))        deallocate(tmp2_BLK)
    if(allocated(time_BLK))        deallocate(time_BLK)
    if(allocated(DivB1_GB))        deallocate(DivB1_GB)
    if(allocated(Ex_CB))           deallocate(Ex_CB, Ey_CB, Ez_CB)
    if(allocated(fbody_x_BLK))     deallocate(fbody_x_BLK, fbody_y_BLK, &
         fbody_z_BLK)
    if(allocated(iTypeAdvance_BP)) deallocate(iTypeAdvance_BP)

    if(iProc==0)then
       call write_prefix
       write(iUnitOut,'(a)') 'clean_mod_advance deallocated arrays'
    end if

  end subroutine clean_mod_advance

  !============================================================================

end module ModAdvance
