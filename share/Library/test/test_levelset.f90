program test_levelset

  use ModLevelSet, test => test_levelset

  implicit none

  integer:: iError
  !---------------------------------------------------------------------------
  call MPI_init(iError)
  call test
  call MPI_finalize(iError)

end program test_levelset

subroutine CON_stop(StringError)

  implicit none
  character (len=*), intent(in) :: StringError
  !----------------------------------------------------------------------------

  write(*,'(a)')StringError

  stop

end subroutine CON_stop

