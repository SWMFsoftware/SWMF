program test_time_convert

  use ModTimeConvert, ONLY: test_time

  implicit none

  call test_time

end program test_time_convert

subroutine CON_stop(StringError)

  implicit none

  character (len=*), intent(in) :: StringError

  write(*,'(a)')StringError

  stop

end subroutine CON_stop
