program test_sort

  use ModSort, ONLY: test => sort_test

  implicit none

  call test

end program test_sort

subroutine CON_stop(StringError)

  implicit none

  character (len=*), intent(in) :: StringError

  write(*,'(a)')StringError

  stop

end subroutine CON_stop
