! !  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
! !  For more information, see http://csem.engin.umich.edu/tools/swmf
module CON_variables
  use ModReadParam, ONLY: lStringLine
  use ModUtilities, ONLY: lVerbose, StringTest, iProcTest

  implicit none

  ! This module contains a number of public variables used by CON
  ! Version number for SWMF
  real, parameter :: VersionSwmf = 2.4

  ! Logical to decide if the SWMF running stand alone or as part of something
  logical :: IsStandAlone = .true.

  ! Error code for the SWMF when it runs as part of a larger application
  integer :: iErrorSwmf = 0

  ! Description string for the problem being solved
  character (len=lStringLine) :: StringDescription='Please describe me!'

  ! Shall we use timing module
  logical :: UseTiming = .true.
  integer :: DnTiming  = -2       ! Show timing at the end of run only

  ! How strict shoule we be.
  ! If true stop with error if false write warning but try to correct problem.
  logical :: UseStrict=.true.

end module CON_variables
!==============================================================================
