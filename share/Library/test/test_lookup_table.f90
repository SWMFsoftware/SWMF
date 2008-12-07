program test_lookup_table

  use ModLookupTable, ONLY: test => test_lookup_table
  use ModMpi

  implicit none

  integer :: iError
  !--------------------------------------------------------------------------
  call MPI_init(iError)
  call test
  call MPI_finalize(iError)

end program test_lookup_table

!=============================================================================

subroutine CON_stop(StringError)

  implicit none

  character (len=*), intent(in) :: StringError

  write(*,'(a)') 'ERROR: '//StringError

  stop

end subroutine CON_stop
