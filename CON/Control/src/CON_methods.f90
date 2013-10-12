! !  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
! !  For more information, see http://csem.engin.umich.edu/tools/swmf
!
!BOP
!
!MODULE: CON_methods - Simple Methods for CON and Components
!
!DESCRIPTION:
! Simple external subroutines which are used by CON and can be used 
! by the science components too. The use of simple external subroutines
! as opposed to module procedures makes the access possible from
! source code writtn in F77 or C/C++. 
! For the same reason there are no optional arguments or other advanced
! F90 features.

!REVISION HISTORY:
!  21Aug03 - Gabor Toth <gtoth@umich.edu> - initial prototype/prolog/code
!            Some methods are based on similar methods used in BATSRUS.
!            Some methods borrow ideas from the MPEU library.
!            Some methods are new.
!  28Aug03 - Added external access functions to the ModIoUnit class
!            for sake of F77/C codes.
!  17Mar04 - The self contained methods are moved into ModUtilities
!EOP

!BOP
!ROUTINE: CON_set_do_test - set logicals for testing
!INTERFACE:
subroutine CON_set_do_test(String,DoTest,DoTestMe)

  !USES:
  use CON_world
  use CON_variables, ONLY: StringTest, tStartTest, nIterStartTest, iProcTest,&
       lVerbose
  use CON_time, ONLY: DoTimeAccurate, tSimulation, nIteration

  implicit none

  !INPUT ARGUMENTS:
  character (len=*), intent(in)  :: String
  !OUTPUT ARGUMENTS:
  logical          , intent(out) :: DoTest, DoTestMe

  !DESCRIPTION:
  ! Based on the input String and the test variables set 
  ! the output variables. Normally DoTest is true if the String
  ! can be found in StringTest, while DoTestMe is true if
  ! DoTest is true and the processor rank is equal to iProcTest.
  ! This behaviour maybe modified by nIterStartTest and tStartTest
  ! which can define the starting iteration and time for the test.
  ! Also the high verbosity levels will produce extra information.
  !EOP
  !----------------------------------------------------------------------------
  !BOC
  if(nIteration >= nIterStartTest .or. &
       (DoTimeAccurate .and. tSimulation>=tStartTest))then

     DoTest   = i_sub_string(' '//StringTest,' '//String//' ')>0
     DoTestMe = DoTest .and. i_proc()==iProcTest
     if(DoTestMe)then
        write(*,*)String,' at iter=',nIteration
     else if(lVerbose>=100)then
        write(*,*)String,' CALLED by me=',i_proc(),' at iter=',nIteration
     else if(i_proc()==iProcTest .and. lVerbose>=10 )then
        write(*,*)String,' CALLED at iter=',nIteration
     endif
  else
     DoTest   = .false.
     DoTestMe = .false.
  end if
  !EOC
contains
  !===========================================================================

  integer function i_sub_string(StringA,StringB)

    ! This is needed to avoid some SGI f90 compiler bug 
    ! (which results in a memory leak) if we use
    !
    ! index(' '//StringTest,' '//str//' ')
    !
    ! directly.

    implicit none

    character (len=*), intent(in) :: StringA, StringB

    i_sub_string=index(StringA, StringB)

  end function i_sub_string

end subroutine CON_set_do_test

!BOP =========================================================================
!ROUTINE: CON_stop - print error message and stop/abort execution
!INTERFACE:
subroutine CON_stop(StringError)

  !USES:
  use ModMpi
  use CON_world, ONLY: world_abort
  implicit none

  !INPUT ARGUMENTS:
  character (len=*), intent(in) :: StringError

  !DESCRIPTION:
  ! This subroutine is used to abort the run with an error report.
  ! It provides an external subroutine interface to CON\_world::world\_abort.
  ! Open IO units are closed and empty output files are deleted before abort.
  ! This will only be done on the aborting processor(s).
  !EOP
  !----------------------------------------------------------------------------
  !BOC
  call world_abort(StringError)
  !EOC
end subroutine CON_stop

!BOP ==========================================================================
!ROUTINE: CON_io_unit_new - provide an unused unit number
!INTERFACE:
subroutine CON_io_unit_new(iUnit)
  !USES
  use ModIoUnit, ONLY: io_unit_new
  implicit none
  !OUTPUT ARGUMENTS:
  integer, intent(out) :: iUnit
  !DESCRIPTION:
  ! This external subroutine is an access method for non-F90 source.
  ! The file should be opened right after the unit number was obtained
  ! so that the unit number gets locked. When the file is closed, the
  ! unit number is automatically released.
  !EOP
  !BOC
  iUnit = io_unit_new()
  !EOC
end subroutine CON_io_unit_new

!BOP ==========================================================================
!ROUTINE: CON_io_unit_tmp - provide a temporary unit number for open and close
!INTERFACE:
subroutine CON_io_unit_tmp(iUnit)
  !USES:
  use ModIoUnit, ONLY: UnitTMP_
  implicit none
  !OUTPUT ARGUMENTS:
  integer, intent(out) :: iUnit
  !DESCRIPTION:
  ! This external subroutine is an access method for non-F90 source.
  ! The file using a temporary unit number should be closed before
  ! any other file could be opened.
  !EOP
  !BOC
  iUnit = UnitTMP_
  !EOC
end subroutine CON_io_unit_tmp

