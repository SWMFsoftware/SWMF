program unit_test_update_iono
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
  write(*,*) 'Initiallizing Planet'

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

  call se_update_state_iono_test
end program unit_test_update_iono

  !============================================================================
  ! UNIT test for update_se_state_iono. Tests update in case with no 
  ! plasmasphere included. Ionospheres are kept decoupled
  subroutine se_update_state_iono_test
    use ModSeBackground
    use ModSeGrid, only:create_se_test_grid,nLine,nPoint,nIono,nPlas
    use ModSeCross,only: SIGS,SIGI,SIGA
    use ModSePlot, only: plot_state,plot_omni_iono
    use ModSeState,only: allocate_state_arrays, iphiup,iphidn,phiup,phidn,&
         initiono,check_time,delt,epsilon,update_se_state_iono,liphiup,liphidn,&
         lphiup,lphidn,specup, specdn    
    integer :: iLine=1, flag=1
    logical :: DoSavePreviousAndReset = .true.
    integer :: nStep
    logical :: IsIono1
    real    :: Ap(7), F107=80, F107A=80, t=0
    !--------------------------------------------------------------------------

    ! First set up the grid that we will update the state in (this is the same 
    ! as the unit test for the grid).
    write(*,*) 'creating grid'
    call create_se_test_grid

    ! Allocate the state arrays
    write(*,*) 'allocating state arrays'
    call allocate_state_arrays
    
    ! Allocate the background right
    write(*,*) 'allocating background arrays'
    call allocate_background_arrays
    
    !align dipole and rotation
    DoAlignDipoleRot = .true.
    
    ! set location of field line
    mLat_I(iLine)=60.0 
    mLon_I(iLine)=0.0

    !set glat and glon coords
    call set_footpoint_locations(iLine)

    ! Get the neutral atmosphere and photo e production spectrum
    AP(:)=4.0
    call get_neutrals_and_pe_spectrum(iLine,F107,F107A,AP)

    ! Fill the background arrays
    write(*,*) 'filling background arrays'
    call fill_thermal_plasma_empirical(iLine,F107,F107A,t)
    
    ! Define the initial state in the ionosphere
    iphiup(iLine,:,:,:)=0.0
    iphidn(iLine,:,:,:)=0.0
    
    liphiup(iLine,:,:,:)=0.0
    liphidn(iLine,:,:,:)=0.0
    
    phiup(iLine,:,:,:)=0.0
    phidn(iLine,:,:,:)=0.0
    
    iphiup(iLine,:,:,:)=0.0
    iphidn(iLine,:,:,:)=0.0
 
    
    ! Set the timestep and convergence criteria
    delt=1.0e5
    epsilon = 0.4

    write(*,*) 'Starting Time loop'
    nStep = 0

    ! plot initial state
    call plot_state(1,nStep,time,iphiup,iphidn,phiup,phidn)

    !do test for iono1
    IsIono1 = .true.
    TIME_LOOP: do while (flag == 1)
       ! Initialize the plasmasphere
       write(*,*) 'Initializing plasmasphere'
       call initiono(1,DoSavePreviousAndReset)
       
       ! update the SE state for iono1
       write(*,*) 'update se state for iono1'
       call update_se_state_iono(iLine,IsIono1,eThermalDensity_IC(iLine,:),&
            eThermalTemp_IC(iLine,:),nNeutralSpecies,&
            NeutralDens1_IIC(iLine,:,:),SIGS,SIGI,SIGA,&
            ePhotoProdSpec1_IIC(iLine,:,:))

       ! update the SE state for iono2
       write(*,*) 'update se state for iono2'
       call update_se_state_iono(iLine,.not.IsIono1,&
            eThermalDensity_IC(iLine,:),&
            eThermalTemp_IC(iLine,:),nNeutralSpecies,&
            NeutralDens2_IIC(iLine,:,:),SIGS,SIGI,SIGA,&
            ePhotoProdSpec2_IIC(iLine,:,:))

       ! check convergence
       write(*,*) 'check for convergence'
       call check_time(1,flag)
       
       ! increment step
       nStep=nStep+1

       ! plot output
       call plot_state(1,nStep,time,iphiup,iphidn,phiup,phidn)
       
       call plot_omni_iono(iLine,nStep,time,specup,specdn,.true.)
       call plot_omni_iono(iLine,nStep,time,specup,specdn,.false.)

    end do TIME_LOOP
    
  end subroutine se_update_state_iono_test


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



