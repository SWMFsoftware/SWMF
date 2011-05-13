!^CFG COPYRIGHT UM
Module ModAdjoint

  !DESCRIPTION: 
  !
  ! This module contains adjoint-specific variables, flags, and
  ! subroutines.  Adjoint capability is available for explicit and
  ! semi-implicit discretizations.  Not all code options are
  ! supported.
  !
  !EOP

  use ModSize
  use ModVarIndexes
  use ModMain,        ONLY: iteration_number, n_step
  use ModIO,          ONLY: dn_output, restart_, iUnitOut, write_prefix
  use ModProcMH,      ONLY: iProc, nProc

  implicit none

  public init_mod_adjoint
  public clean_mod_adjoint
  public read_adjoint_parameters
  public store_block_buffer
  public recall_block_buffer

  ! Logical flag  indicating that this is an adjoint run
  logical :: DoAdjoint = .false.

  ! How often should adjoint be saved?
  integer :: DnSaveAdjoint = -1
  real    :: DtSaveAdjoint = -1.0

  !\
  ! Forward solution buffer to avoid reading restarts at every iteration
  !/
  integer           :: nBuffer = -1
  integer           :: iBuffer = -1
  real, allocatable :: Buffer_State_VGB(:,:,:,:,:,:)

  !\
  ! Simulation time buffer
  !/
  real, allocatable :: Buffer_time_simulation(:)

  !\
  ! Adjoint solution, and temp/prev value
  !/
  real, allocatable :: Adjoint_VGB(:,:,:,:,:)
  real, allocatable :: AdjEnergy_GBI(:,:,:,:,:)
  real, allocatable :: AdjointPrev_VGB(:,:,:,:,:)
  real, allocatable :: AdjEnergyPrev_GBI(:,:,:,:,:)

  !\
  ! Local cell-centered adjoint source terms
  !/
  real :: AdjSource_VC(nVar+nFluid, nI, nJ, nK)

  !\
  ! Enumerated types for temporary block buffer storage of state
  !/
  integer, parameter :: &
       AdjPreUpdate_   = 1,       &
       AdjUserUpdate1_ = 2,       &
       AdjUserUpdate2_ = 3,       &
       AdjPreEnergyP_  = 4,       &
       AdjLast_        = 5

  ! Block buffer for temporary state storage
  real :: BlockBuf_State(nVar+nFluid, nI, nJ, nK,AdjLast_-1) = 0.0

  !\
  ! X Face local MHD adjoint solution array definitions.
  !/
  ! These are primitive variables (velocity)
  real :: AdjLeftState_VX(nVar,2-gcn:nI+gcn,0:nJ+1,0:nK+1)
  real :: AdjRightState_VX(nVar,2-gcn:nI+gcn,0:nJ+1,0:nK+1)

  ! Fluxes are for conservative variables (momentum)
  real :: AdjFlux_VX(nVar+nFluid,0:nI+1,2-gcn:nJ+gcn,0:nK+1)

  real :: AdjuDotArea_XI(2-gcn:nI+gcn,0:nJ+1,0:nK+1,nFluid+1)

  !\
  ! Y Face local MHD adjoint solution array definitions.
  !/
  real :: AdjLeftState_VY(nVar,0:nI+1,2-gcn:nJ+gcn,0:nK+1)
  real :: AdjRightState_VY(nVar,0:nI+1,2-gcn:nJ+gcn,0:nK+1)

  real :: AdjFlux_VY(nVar+nFluid,0:nI+1,2-gcn:nJ+gcn,0:nK+1)
  
  real :: AdjuDotArea_YI(0:nI+1,2-gcn:nJ+gcn,0:nK+1,nFluid+1)


  !\
  ! Z Face local MHD adjoint solution array definitions.
  !/
  real :: AdjLeftState_VZ(nVar,0:nI+1,0:nJ+1,2-gcn:nK+gcn)
  real :: AdjRightState_VZ(nVar,0:nI+1,0:nJ+1,2-gcn:nK+gcn)

  real :: AdjFlux_VZ(nVar+nFluid,0:nI+1,0:nJ+1,2-gcn:nK+gcn)

  real :: AdjuDotArea_ZI(0:nI+1,0:nJ+1,2-gcn:nK+gcn,nFluid+1)



contains

  !============================================================================

  subroutine init_mod_adjoint

    use ModMain,      ONLY: nIter, time_accurate

    character(len=*), parameter:: NameSub = 'init_mod_adjoint'

    ! The adjoint mode relies on including the last restart.H file from
    ! the forward run into the PARAM.in of the adjoint run.
    ! From the read of the restart.H file, n_step is set to the final
    ! iteration number
    iteration_number = n_step                 ! start from last iteration
    nIter            = n_step                 ! probably not necessary
    nBuffer          = dn_output(restart_)+1  ! number of buffer states

    ! check if adjoint is supported
    if(.not.time_accurate ) call stop_mpi(NameSub//': need time_accurate = true')
    if(iteration_number<=0) call stop_mpi(NameSub//': need iteration_number > 0')
    if(nBuffer         <=0) call stop_mpi(NameSub//': need DnSaveRestart > 0')

    if(.not.allocated(Buffer_State_VGB))then
       allocate(Buffer_State_VGB(nVar,1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn,MaxBlock,nBuffer))
       allocate(     Adjoint_VGB(nVar,1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn,MaxBlock))
       allocate( AdjointPrev_VGB(nVar,1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn,MaxBlock))
       allocate(     AdjEnergy_GBI(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn,MaxBlock,nFluid))
       allocate( AdjEnergyPrev_GBI(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn,MaxBlock,nFluid))

       allocate(Buffer_time_simulation(nBuffer))

       if(iProc==0)then
          call write_prefix
          write(iUnitOut,'(a)') 'init_mod_adjoint allocated arrays'
       end if

    end if

    call init_adjoint_solution

  end subroutine init_mod_adjoint

  !============================================================================

  subroutine clean_mod_adjoint

    if(.not.allocated(Buffer_State_VGB)) return
    deallocate(Buffer_State_VGB)
    deallocate(Adjoint_VGB)
    deallocate(AdjointPrev_VGB)
    deallocate(AdjEnergy_GBI)
    deallocate(AdjEnergyPrev_GBI)

  end subroutine clean_mod_adjoint

  !============================================================================


  subroutine read_adjoint_parameters(NameCommand)

    use ModReadParam, ONLY: read_var
    
    character(len=*), intent(in) :: NameCommand
    character(len=*), parameter:: NameSub = 'read_adjoint_parameters'
    !--------------------------------------------------------------------------
    
    select case(NameCommand)
    case("#ADJOINT")
       ! DoAdjoint
       ! Starting iteration stored in MaxIteration (read in elsewhere)
       ! DnSaveAdjoint, DtSaveAdjoint
       call read_var('DoAdjoint',DoAdjoint)
       if(DoAdjoint)then
          call read_var('DnSaveAdjoint',DnSaveAdjoint)
          call read_var('DtSaveAdjoint',DtSaveAdjoint)
       end if
    case default
       call stop_mpi(NameSub//' unknown NameCommand='//NameCommand)
    end select
    
  end subroutine read_adjoint_parameters

  !============================================================================

  subroutine init_adjoint_solution

    character(len=*), parameter:: NameSub = 'init_adjoint_solution'

    ! initialize adjoint at final time
    ! call stop_mpi(NameSub//': To be implemented')
    Adjoint_VGB(:,:,:,:,:) = 0.

    ! TODO: add source

  end subroutine init_adjoint_solution


  !============================================================================

  subroutine store_block_buffer(iBlock, ind)

    use ModAdvance,    ONLY: State_VGB, Energy_GBI

    integer, intent(in) :: iBlock
    integer, intent(in) :: ind
    
    integer :: i,j,k

    ! store State and Energy in a block buffer
    do k=1,nK; do j=1,nJ; do i=1,nI
       BlockBuf_State(1:nVar,i,j,k, ind) = State_VGB(1:nVar,i,j,k,iBlock)
       BlockBuf_State(nVar+1:nVar+nFluid,i,j,k, ind) = Energy_GBI(i,j,k,iBlock,1:nFluid)
    end do; end do; end do
  end subroutine store_block_buffer

  !============================================================================

  
  subroutine recall_block_buffer(iBlock, ind)

    use ModAdvance,    ONLY: State_VGB, Energy_GBI

    integer, intent(in) :: iBlock
    integer, intent(in) :: ind
    
    integer :: i,j,k
    
    ! store State and Energy in a block buffer
    do k=1,nK; do j=1,nJ; do i=1,nI
       State_VGB(1:nVar,i,j,k,iBlock) = BlockBuf_State(1:nVar,i,j,k, ind)
       Energy_GBI(i,j,k,iBlock,1:nFluid) = BlockBuf_State(nVar+1:nVar+nFluid,i,j,k, ind)
    end do; end do; end do

  end subroutine recall_block_buffer
  

end Module ModAdjoint
