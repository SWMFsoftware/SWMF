module ModMpiOrig
  implicit none
  integer, parameter :: iRealPrec = 0
  include 'mpif90.h'
end module ModMpiOrig

module ModMpi
  implicit none
  integer, parameter :: iRealPrec = (1.00000000011 - 1.0)*10000000000.0
  include 'mpif90.h'
end module ModMpi
