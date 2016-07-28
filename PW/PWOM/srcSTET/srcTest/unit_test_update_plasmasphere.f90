program unit_test_update_plasmasphere
  use ModSeGrid
  use ModSeMpi

  use ModMPI
  use CON_planet, ONLY: init_planet_const, set_planet_defaults

  integer :: iError
  
  !-----------------------------------------------------------------------------

  !****************************************************************************
  ! Initiallize MPI and get number of processors and rank of given processor
  !****************************************************************************

  write(*,*) 'Initiallizing MPI'

  !---------------------------------------------------------------------------
  call MPI_INIT(iError)
  iComm = MPI_COMM_WORLD

  call MPI_COMM_RANK(iComm,iProc,iError)
  call MPI_COMM_SIZE(iComm,nProc,iError)

  !\
  ! Initialize the planetary constant library and set Earth
  ! as the default planet.
  !/
  write(*,*) 'Initiallizing Planet'

  call init_planet_const
  call set_planet_defaults

  write(*,*) 'starting se_grid_test'
!  call se_grid_test
  call se_update_state_test
!  call background_test
!  call se_update_state_iono_test
!  call se_update_state_iono_test_transport
end program unit_test_update_plasmasphere

!============================================================================
! UNIT test for SE update states
subroutine se_update_state_test
  use ModSeGrid, only:create_se_test_grid,nLine,nPoint,nIono,nPlas
  use ModSePlot, only: plot_state
  use ModSeState,only: allocate_state_arrays, iphiup,iphidn,phiup,phidn,&
                       initplas,check_time,delt,epsilon, update_se_state, &
                       liphiup,liphidn
  integer :: iLine=1, flag=1
  logical :: DoSavePreviousAndReset = .true.
  integer :: nStep
  logical :: IsOpen
  
  ! thermal background  
  real, allocatable :: eThermalDensity_C(:),eThermalTemp_C(:)
  
  !--------------------------------------------------------------------------
  ! Allocate the background right
  if(.not.allocated(eThermalDensity_C)) &
       allocate(eThermalDensity_C(nPoint))
  if(.not.allocated(eThermalTemp_C)) &
       allocate(eThermalTemp_C(nPoint))
  
  ! First set up the grid that we will update the state in (this is the same 
  ! as the unit test for the grid).
  write(*,*) 'creating grid'
  call create_se_test_grid
  
  ! Allocate the state arrays
  write(*,*) 'allocating state arrays'
  call allocate_state_arrays
  
  ! Fill the background arrays
  !    eThermalDensity_IC(:,:) = 0.00001
  eThermalDensity_C(:) = 1.0e2
  eThermalTemp_C(:) = 1.0
  
  ! Define the initial state in the ionosphere
  do iLine=1,nLine
     iphiup(iLine,:,:,:)=1.0e5
     iphidn(iLine,:,:,:)=1.0e5
     
     liphiup(iLine,:,:,:)=1.0e5
     liphidn(iLine,:,:,:)=1.0e5
     
     phiup(iLine,:,:,:)=0.00001
     phidn(iLine,:,:,:)=0.00001
  end do
  
  ! Set the timestep and convergence criteria
  delt=1.0e5
  epsilon = 0.4
  
  write(*,*) 'Starting Time loop'
  nStep = 0
  
  ! plot initial state
  call plot_state(1,nStep,time,iphiup,iphidn,phiup,phidn)
  
  IsOpen = .true.
  TIME_LOOP: do while (flag == 1)
     ! Initialize the plasmasphere
     write(*,*) 'Initializing plasmasphere'
     call initplas(1,DoSavePreviousAndReset)
     
     ! update the SE state
     write(*,*) 'update se state'
     call update_se_state(1, eThermalDensity_C,&
          eThermalTemp_C,IsOpen)
     
     ! check convergence
     write(*,*) 'check for convergence'
     call check_time(1,flag)
     
     ! increment step
     nStep=nStep+1
     
     ! plot output
     call plot_state(1,nStep,time,iphiup,iphidn,phiup,phidn)
     
     
  end do TIME_LOOP
  
end subroutine se_update_state_test


!============================================================================
! The following subroutines are here so that we can use SWMF library routines
! Also some features available in SWMF mode only require empty subroutines
! for compilation of the stand alone code.
!============================================================================
subroutine CON_stop(StringError)
  use ModSeMpi, ONLY : iProc,iComm
  use ModMpi
  implicit none
  character (len=*), intent(in) :: StringError

  ! Local variables:
  integer :: iError,nError
  !----------------------------------------------------------------------------

  write(*,*)'Stopping execution! me=',iProc,&
       ' with msg:'
  write(*,*)StringError
  call MPI_abort(iComm, nError, iError)
  stop

end subroutine CON_stop

subroutine CON_set_do_test(String,DoTest,DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe

  DoTest = .false.; DoTestMe = .false.

end subroutine CON_set_do_test

subroutine CON_io_unit_new(iUnit)

  use ModIoUnit, ONLY: io_unit_new
  implicit none
  integer, intent(out) :: iUnit

  iUnit = io_unit_new()

end subroutine CON_io_unit_new

