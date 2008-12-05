program test_lookup_table

  use ModLookupTable, ONLY: test => test_lookup_table

  implicit none

  call test

end program test_lookup_table

subroutine CON_stop(StringError)

  implicit none

  character (len=*), intent(in) :: StringError

  write(*,'(a)') 'ERROR: '//StringError

  stop

end subroutine CON_stop
