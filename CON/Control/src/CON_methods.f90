! !  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
! !  For more information, see http://csem.engin.umich.edu/tools/swmf
!
!
!
! Simple external subroutines which are used by CON and can be used
! by the science components too. The use of simple external subroutines
! as opposed to module procedures makes the access possible from
! source code writtn in F77 or C/C++.
! For the same reason there are no optional arguments or other advanced
! F90 features.

! revision history:
!  21Aug03 - Gabor Toth <gtoth@umich.edu> - initial prototype/prolog/code
!            Some methods are based on similar methods used in BATSRUS.
!            Some methods borrow ideas from the MPEU library.
!            Some methods are new.
!  28Aug03 - Added external access functions to the ModIoUnit class
!            for sake of F77/C codes.
!  17Mar04 - The self contained methods are moved into ModUtilities

! ROUTINE: CON_set_do_test - set logicals for testing
subroutine CON_set_do_test(String,DoTest,DoTestMe)

  use ModUtilities, ONLY: util_set_do_test => CON_set_do_test

  implicit none

  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe

  ! See ModUtilities.
  !----------------------------------------------------------------------------
  call util_set_do_test(String, DoTest, DoTestMe)

end subroutine CON_set_do_test
!==============================================================================

! ROUTINE: CON_stop - print error message and stop/abort execution
subroutine CON_stop(StringError)

  use ModUtilities, ONLY: util_stop => CON_stop
  implicit none

  character (len=*), intent(in) :: StringError

  ! This subroutine is used to abort the run with an error report.
  ! It provides an external subroutine interface to ModUtilities::CON\_stop.
  ! Open I/O units are closed and empty output files are deleted before abort.
  ! This will only be done on the aborting processor(s).

  !----------------------------------------------------------------------------
  call util_stop(StringError)
end subroutine CON_stop
!==============================================================================

! ROUTINE: CON_io_unit_new - provide an unused unit number
subroutine CON_io_unit_new(iUnit)
  ! USES
  use ModIoUnit, ONLY: io_unit_new
  implicit none
  integer, intent(out) :: iUnit
  ! This external subroutine is an access method for non-F90 source.
  ! The file should be opened right after the unit number was obtained
  ! so that the unit number gets locked. When the file is closed, the
  ! unit number is automatically released.
  !----------------------------------------------------------------------------
  iUnit = io_unit_new()
end subroutine CON_io_unit_new
!==============================================================================

! ROUTINE: CON_io_unit_tmp - provide a temporary unit number for open and close
subroutine CON_io_unit_tmp(iUnit)
  use ModIoUnit, ONLY: UnitTMP_
  implicit none
  integer, intent(out) :: iUnit
  ! This external subroutine is an access method for non-F90 source.
  ! The file using a temporary unit number should be closed before
  ! any other file could be opened.
  !----------------------------------------------------------------------------
  iUnit = UnitTMP_
end subroutine CON_io_unit_tmp
!==============================================================================

