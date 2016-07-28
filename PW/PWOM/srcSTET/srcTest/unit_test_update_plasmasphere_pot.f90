program unit_test_update_plasmasphere_pot
  use ModSeGrid
  use ModSeMpi

  use ModMPI
  use CON_planet, ONLY: init_planet_const, set_planet_defaults,is_planet_init

  integer :: iError
!  character(len=5) :: NamePlanet = 'EARTH'
  character(len=7) :: NamePlanet = 'JUPITER'
  logical :: IsPlanetSet=.false.  
  
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
  write(*,*) 'Initiallizing Planet ', NamePlanet

  call init_planet_const

  if (NamePlanet == 'EARTH') then
     call set_planet_defaults
     IsPlanetSet = .true.
  else
     IsPlanetSet = is_planet_init(NamePlanet)
  endif
  
  if (.not.IsPlanetSet) then
     call CON_stop('Planet not set. Stopping STET')
  endif

  write(*,*) 'starting se_grid_test'

  call se_update_state_test

end program unit_test_update_plasmasphere_pot

!============================================================================
! UNIT test for SE update states
subroutine se_update_state_test
  use ModSeGrid, only:create_se_test_grid,set_grid_pot,nLine,nPoint,nIono,&
       nPlas,nTop,Efield_IC,FieldLineGrid_IC,nThetaAlt_IIC,nThetaAlt_II,nEnergy
  use ModSePlot, only: plot_state,plot_state_pot,plot_omni_pot
  use ModSeState,only: allocate_state_arrays, iphiup,iphidn,phiup,phidn,&
                       specup,specdn,initplas_pot,check_time_pot,delt,epsilon, &
                       update_se_state_pot, &
                       liphiup,liphidn
  integer :: iLine=1, flag=1
  logical :: DoSavePreviousAndReset = .true.
  integer :: nStep, iPoint, iEnergy
  logical :: IsOpen
  real    :: time
  
  ! thermal background  
  real, allocatable :: eThermalDensity_C(:),eThermalTemp_C(:)
  real,parameter :: cCmToM = 1.0e-2
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
  
  ! input electric field (assume constant in plasmapshere and zero in iono)
  ! choose values to have a 5 V drop from top of iono to equator
  Efield_IC(iLine,:)=0.0
  Efield_IC(iLine,nIono+1:nTop-1)=5.0&
       /((FieldLineGrid_IC(iLine,nTop)-FieldLineGrid_IC(iLine,nIono))*cCmToM)
  Efield_IC(iLine,nTop+1:nPoint-nIono)=-5.0&
       /((FieldLineGrid_IC(iLine,nTop)-FieldLineGrid_IC(iLine,nIono))*cCmToM)

  ! set the region of exisitence for the defined Efield
  call set_grid_pot(iLine)

!  do iPoint=1,nPoint
!     write(*,*)iPoint,nThetaAlt_II(iLine,iPoint),nThetaAlt_IIC(iLine,50,iPoint)
!  enddo
!  call con_stop('')


  ! Allocate the state arrays
  write(*,*) 'allocating state arrays'
  call allocate_state_arrays
  
  ! Fill the background arrays
  !eThermalDensity_C(:) = 0.00001
  eThermalDensity_C(:) = 1.0e2
  eThermalTemp_C(:) = 1.0
  
  !set initial states
  do iLine=1,nLine
     ! Define the initial state in the ionosphere
     
     iphiup(iLine,:,:,:)=0
     iphidn(iLine,:,:,:)=0
     
     liphiup(iLine,:,:,:)=0.0
     liphidn(iLine,:,:,:)=0.0
     !just fill flux in upward loss cone in ionosphere
     do iEnergy=1,nEnergy
        iphiup (iLine,0:nThetaAlt_IIC(iLine,iEnergy,nIono),1:nIono,iEnergy)=1.0e5
        liphiup(iLine,0:nThetaAlt_IIC(iLine,iEnergy,nIono),1:nIono,iEnergy)=1.0e5
     end do

     !Initial plasmasphere set to small number
     phiup(iLine,:,:,:)=0.0001
     phidn(iLine,:,:,:)=0.0001
  end do
  
  ! Set the timestep and convergence criteria
  delt=1.0e5
  epsilon = 0.1
  
  write(*,*) 'Starting Time loop'
  nStep = 0
  time=0
  ! plot initial state
  !call plot_state(1,nStep,time,iphiup,iphidn,phiup,phidn)
  
  IsOpen = .false.
  TIME_LOOP: do while (flag == 1)
     ! Initialize the plasmasphere
     write(*,*) 'Initializing plasmasphere'
     call initplas_pot(1,DoSavePreviousAndReset)
     
     ! update the SE state
     write(*,*) 'update se state_pot'
     call update_se_state_pot(1, eThermalDensity_C,&
          eThermalTemp_C,IsOpen)
     
     ! check convergence
     write(*,*) 'check for convergence'
     call check_time_pot(1,flag)
     
     ! increment step
     nStep=nStep+1
     
     ! plot output
!     call plot_state(1,nStep,time,iphiup,iphidn,phiup,phidn)
     
     call plot_state_pot(1,90,nStep,time,iphiup,iphidn,phiup,phidn)

     call plot_state_pot(1,1,nStep,time,iphiup,iphidn,phiup,phidn)
     call plot_state_pot(1,2,nStep,time,iphiup,iphidn,phiup,phidn)
     call plot_state_pot(1,3,nStep,time,iphiup,iphidn,phiup,phidn)
     call plot_state_pot(1,4,nStep,time,iphiup,iphidn,phiup,phidn)
     call plot_state_pot(1,5,nStep,time,iphiup,iphidn,phiup,phidn)
     call plot_state_pot(1,6,nStep,time,iphiup,iphidn,phiup,phidn)
     

     call plot_omni_pot(1,nStep,time,specup,specdn)
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

