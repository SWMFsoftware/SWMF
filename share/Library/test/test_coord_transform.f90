program test_coord_transform

  use ModCoordTransform, ONLY: test => test_coord_transform

  implicit none

  call test

end program test_coord_transform

subroutine CON_stop(StringError)

  implicit none

  character (len=*), intent(in) :: StringError

  write(*,'(a)')StringError

  stop

end subroutine CON_stop
