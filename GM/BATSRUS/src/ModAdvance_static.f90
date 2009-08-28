!^CFG COPYRIGHT UM
Module ModAdvance
  use ModSize
  use ModVarIndexes
  use ModIO,         ONLY: iUnitOut, write_prefix
  use ModProcMH,     ONLY: iProc, nProc
  use ModB0

  implicit none
  save

  ! Logical parameter indicating static vs. dynamic allocation
  logical, parameter :: IsDynamicAdvance = .false.

  ! Numerical flux type
  character (len=10) :: FluxType

  ! Named index for electron pressure and velocity
  integer, parameter:: eFluid_ = nFluid + 1

  ! Additional equations are activated by the named indices of the
  ! corresponding extra variables in the equation modules.
  ! The default values of these named indices are in ModExtraVariables.
  logical, parameter:: UseElectronPressure = Pe_ > 1
  logical, parameter:: UseParallelPressure = Ppar_ > 1
  logical, parameter:: UseIdealEos = ExtraEint_ == 1
  logical, parameter:: UseElectronEnergy = Ee_ > 1

  ! Number of Rosseland, Planck opacities used for radiation diffusion
  integer, parameter:: nOpacity = max(1, nWave)

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

  ! Use the same speeds in the solvers for all fluids, or not
  logical :: UseTotalSpeed = .true.

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
  real:: State_VGB(nVar,1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn, MaxBlock)
  real:: Energy_GBI(1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn, MaxBlock, nFluid)

  !\
  ! Block cell-centered MHD solution old state
  !/
  real :: StateOld_VCB(nVar, nI, nJ, nK, MaxBlock)
  real :: EnergyOld_CBI(nI, nJ, nK, MaxBlock, nFluid)

  !\
  ! Block cell-centered intrinsic magnetic field, time, and temporary storage
  !/
  real,  dimension(1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn, MaxBlock) :: &
       tmp1_BLK, tmp2_BLK

  real :: time_BLK(nI, nJ, nK, MaxBlock)

  ! Array for storing dB0/dt derivatives
  real, allocatable :: Db0Dt_CDB(:,:,:,:,:)

  ! Arrays for the total electric field
  real, dimension(nI, nJ, nK, MaxBlock) :: Ex_CB, Ey_CB, Ez_CB

  !\
  ! Block cell-centered body forces
  !/
  real, dimension(nI, nJ, nK, MaxBlock) :: &
       fbody_x_BLK, fbody_y_BLK, fbody_z_BLK

  !\
  ! Local cell-centered source terms and divB.
  !/
  real :: Source_VC(nVar+nFluid, nI, nJ, nK)
  real :: Theat0(nI,nJ,nK)
  real :: DivB1_GB(1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn, MaxBlock)

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
  ! The number of the face variables, which are corrected at the
  ! resolution changes
  !/

  !\
  ! Merge cells around the polar axis in spherical geometry
  !/
  logical :: DoFixAxis = .false.
  real ::    rFixAxis = 0.0, r2FixAxis = 0.0 

  !\
  ! Block type information
  !/
  integer              :: iTypeAdvance_B(MaxBlock)
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

    if(allocated(iTypeAdvance_BP)) RETURN
    allocate(iTypeAdvance_BP(MaxBlock,0:nProc-1))
    iTypeAdvance_B  = SkippedBlock_
    iTypeAdvance_BP = SkippedBlock_

    if(IsDynamicAdvance .and. iProc==0)then
       call write_prefix
       write(iUnitOut,'(a)') 'init_mod_advance allocated arrays'
    end if

  end subroutine init_mod_advance

  !============================================================================

  subroutine clean_mod_advance

    if(allocated(iTypeAdvance_BP)) deallocate(iTypeAdvance_BP)

    if(IsDynamicAdvance .and. iProc==0)then
       call write_prefix
       write(iUnitOut,'(a)') 'clean_mod_advance deallocated arrays'
    end if

  end subroutine clean_mod_advance

  !============================================================================

end Module ModAdvance
