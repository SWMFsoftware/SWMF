program test_freq

  use ModFreq, ONLY: test => test_freq

  implicit none

  call test

end program test_freq

subroutine CON_stop(StringError)

  implicit none

  character (len=*), intent(in) :: StringError

  write(*,'(a)')StringError

  stop

end subroutine CON_stop
