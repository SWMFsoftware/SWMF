program test_axes

  use CON_axes, ONLY: test => test_axes

  implicit none

  call test

end program test_axes

subroutine CON_stop(StringError)

  implicit none

  character (len=*), intent(in) :: StringError

  write(*,'(a)')StringError

  stop

end subroutine CON_stop

subroutine CON_set_do_test(String,DoTest,DoTestMe)

  implicit none

  character (len=*), intent(in) :: String
  logical, intent(out) :: DoTest,DoTestMe

  DoTest   = .false.
  DoTestMe = .false.

end subroutine CON_set_do_test
