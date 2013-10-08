program test_interpolate_amr

  use ModInterpolateAMRGrid, test => test_interpolate_amr_grid

  implicit none

  call test

end program test_interpolate_amr

subroutine CON_stop(StringError)

  implicit none
  character (len=*), intent(in) :: StringError
  !----------------------------------------------------------------------------

  write(*,'(a)')StringError
  write(*,'(a)')'!!! SWMF_ABORT !!!'
  stop

end subroutine CON_stop

