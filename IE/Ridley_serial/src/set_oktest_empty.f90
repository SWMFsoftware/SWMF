!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine set_oktest(str,oktest,oktest_me)

  implicit none

  character (len=*) :: str
  logical :: oktest, oktest_me

  oktest = .true.
  oktest_me = .true.

  ! This routine does nothing.

end subroutine set_oktest
