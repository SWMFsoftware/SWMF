!^CFG COPYRIGHT UM
!^CFG FILE USERFILES
!========================================================================
module ModUser
  ! This is the default user module which contains empty methods defined
  ! in ModUserEmpty.f90
  !
  ! Please see the documentation, and the files ModUserEmpty.f90 and 
  ! srcUser/ModUserExamples.f90 for information about what the different user
  ! subroutines do and how to implement them for your specific problem.
  use ModUserEmpty
  include 'user_module.h' !list of public methods

  !\
  ! Here you must define a user routine Version number and a 
  ! descriptive string.
  !/
  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: &
       NameUserModule = 'DEFAULT USER ROUTINES, K.C. Hansen'

end module ModUser
