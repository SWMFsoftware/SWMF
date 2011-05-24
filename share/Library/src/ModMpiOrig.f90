module ModMpiOrig

  implicit none
  ! iRealPrec = 1 if the code is compiled with 8 byte reals and 0 otherwise
  integer, parameter :: iRealPrec = (1.00000000011 - 1.0)*10000000000.0
  include 'mpif.h'
  

end module ModMpiOrig
