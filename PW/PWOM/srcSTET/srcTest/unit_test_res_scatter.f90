program unit_test_res_scatter
!  use ModSeGrid
!  use ModSeMpi
!  use ModMPI
!  use CON_planet, ONLY: init_planet_const, set_planet_defaults

!  integer :: iError
  use ModSeProduction,only: unit_test_plas_resonant_scattering
  !-----------------------------------------------------------------------------

  !****************************************************************************
  ! Initiallize MPI and get number of processors and rank of given processor
  !****************************************************************************

!  write(*,*) 'Initiallizing MPI'

  !---------------------------------------------------------------------------
!  call MPI_INIT(iError)
!  iComm = MPI_COMM_WORLD

!  call MPI_COMM_RANK(iComm,iProc,iError)
!  call MPI_COMM_SIZE(iComm,nProc,iError)

  !\
  ! Initialize the planetary constant library and set Earth
  ! as the default planet.
  !/
 ! write(*,*) 'Initiallizing Planet'

!  call init_planet_const
!  call set_planet_defaults


  call unit_test_plas_resonant_scattering
end program unit_test_res_scatter
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
