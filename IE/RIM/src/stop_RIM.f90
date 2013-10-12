!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine stop_RIM(str)

  implicit none

  character (len=*), intent(in) :: str

  call CON_stop("IE/RIM Error: "//str)

end subroutine stop_RIM
