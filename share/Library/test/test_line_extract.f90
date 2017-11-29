!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program test_line

  use CON_line_extract

  implicit none

  integer :: iError

  call MPI_init(iError)
  call line_test
  call MPI_finalize(iError)

end program test_line


