!^CFG COPYRIGHT UM

subroutine set_oktest(str,oktest,oktest_me)

  implicit none

  character (len=*) :: str
  logical :: oktest, oktest_me

  oktest = .true.
  oktest_me = .true.

  ! This routine does nothing.

end subroutine set_oktest
