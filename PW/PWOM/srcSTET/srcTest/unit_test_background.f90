program unit_test_background
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

  call background_test
end program unit_test_background

  !============================================================================
  ! UNIT test for SE update states
  subroutine background_test
    use ModSeGrid, only:create_se_test_grid,nLine,nPoint,nIono,nPlas
    use ModSeBackground,only: allocate_background_arrays,mLat_I,mLon_I, &
         set_footpoint_locations,fill_thermal_plasma_empirical,plot_background,&
         plot_ephoto_prod,DoAlignDipoleRot,get_neutrals_and_pe_spectrum,&
         plot_ionization_rate
    
    integer :: iLine=1, flag=1, nStep=0
    real    :: time=0
    logical :: DoSavePreviousAndReset = .true.
    
    real :: Ap(7), F107=80, F107A=80, t=0
    !--------------------------------------------------------------------------

    ! First set up the grid that we will update the state in (this is the same 
    ! as the unit test for the grid).
    write(*,*) 'creating grid'
    call create_se_test_grid

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
    
    ! Fill the background arrays
    write(*,*) 'filling background arrays'
    call fill_thermal_plasma_empirical(iLine,F107,F107A,t)
    
    ! Get the neutral atmosphere and photo e production spectrum
    AP(:)=4.0
    write(*,*) 'get_neutrals_and_pe_spectrum'
    call get_neutrals_and_pe_spectrum(iLine,F107,F107A,AP)

    write(*,*) 'plot_background'
    ! plot initial state
    call plot_background(iLine,nStep,time)
    
    write(*,*) 'plot_ephoto_prod'
    ! plot ephoto production
    call plot_ephoto_prod(iLine,nStep,time)

    write(*,*) 'plot_ionization_rate'
    ! plot ionization rate
    call plot_ionization_rate(iLine,nStep,time)

  end subroutine background_test

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



