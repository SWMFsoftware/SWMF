program test_geopack

  use CON_geopack, ONLY: test => CON_test_geopack

  implicit none

  call test

end program test_geopack

subroutine CON_stop(StringError)

  implicit none

  character (len=*), intent(in) :: StringError

  write(*,'(a)')StringError

  stop

end subroutine CON_stop
