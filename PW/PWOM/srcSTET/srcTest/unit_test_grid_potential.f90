program unit_test_grid_potential
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
  write(*,*) 'starting se_grid_test'

  call grid_potential_test

end program unit_test_grid_potential

!============================================================================
! UNIT test for SE update states
subroutine grid_potential_test
  use ModSeGrid, only:create_se_test_grid,set_grid_pot,plot_grid_pot,&
       nLine,nPoint,nIono,nPoint,nTop,Efield_IC,FieldLineGrid_IC
  !use ModSePlot, only: plot_along_field

  integer :: iLine=1, flag=1
  logical :: DoSavePreviousAndReset = .true.
  integer :: nStep
  logical :: IsOpen
  real,parameter :: cCmToM = 1.0e-2
  !--------------------------------------------------------------------------
  
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

  ! plot the potential, and min total energy considered
  !call plot_along_field(iLine)

  !write output for several total energies
  call plot_grid_pot(iLine,0,90,0.0)
  call plot_grid_pot(iLine,0,80,0.0)
  call plot_grid_pot(iLine,0,70,0.0)
  call plot_grid_pot(iLine,0,60,0.0)
  call plot_grid_pot(iLine,0,50,0.0)
  call plot_grid_pot(iLine,0,40,0.0)
  call plot_grid_pot(iLine,0,30,0.0)
  call plot_grid_pot(iLine,0,20,0.0)
  call plot_grid_pot(iLine,0,10,0.0)
  call plot_grid_pot(iLine,0,10,0.0)
  call plot_grid_pot(iLine,0,9,0.0)
  call plot_grid_pot(iLine,0,8,0.0)
  call plot_grid_pot(iLine,0,7,0.0)
  call plot_grid_pot(iLine,0,6,0.0)
  call plot_grid_pot(iLine,0,5,0.0)
  call plot_grid_pot(iLine,0,4,0.0)
  call plot_grid_pot(iLine,0,3,0.0)
  call plot_grid_pot(iLine,0,2,0.0)
  call plot_grid_pot(iLine,0,1,0.0)

end subroutine grid_potential_test



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

