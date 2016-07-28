program main_stet
  use ModSeGrid
  use ModSeMpi
  use ModSeBackground,only: ZEP
  use ModMPI
  use CON_planet, ONLY: init_planet_const, set_planet_defaults,is_planet_init
  use ModNumConst,    ONLY: cRadToDeg
  implicit none
  
  integer :: iError
  
  integer :: nLineIn, iLine
  real :: L

  real :: Coord_ID(1,2) ! Lat and Lon in degrees for each line
  real :: Ap_I(7), F107, F107A
  integer,parameter :: Lat_=1 ,Lon_=2 !named parameters for Coord_ID

  real :: Time = 0.0

  character(len=5) :: NamePlanet = 'EARTH'
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
  write(*,*) 'Initiallizing Planet',NamePlanet

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
  
  !\
  ! Parameters should be read here, for now just set them by hand
  !/ 
  nLineIn=1
!  L = 4.0
  L = 8.55
  Coord_ID(1,Lat_)=acos(sqrt(1.0/L))*cRadToDeg
  Coord_ID(1,Lon_)=0.0
  AP_I(:)=4.0
  F107 = 60.0
  F107A= 60.0
  ZEP=2
  ! Include the parallel electric field?
  DoIncludePotential=.true.

  !\
  ! Initialize the stet code
  !/
  write(*,*) 'Initializing stet'
  call stet_init(nLineIn,Coord_ID,Ap_I,F107,F107A, Time)
  
  
  write(*,*) 'Running stet'
  do iLine=1,nLineIn
     call stet_run(iLine,.true.,.false.)
  end do

end program main_stet


!============================================================================
! The following subroutines are here so that we can use SWMF library routines
! Also some features available in SWMF mode only require empty subroutines
! for compilation of the stand alone code.
!============================================================================
subroutine CON_stop(StringError)
  use ModSeMpi, ONLY : iProc,iComm
  use ModSeState, ONLY : Time
  use ModMpi
  implicit none
  character (len=*), intent(in) :: StringError

  ! Local variables:
  integer :: iError,nError
  !----------------------------------------------------------------------------

  write(*,*)'Stopping execution! me=',iProc,' at time=',Time,&
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

