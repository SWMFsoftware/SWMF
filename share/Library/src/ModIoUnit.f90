!BOP -------------------------------------------------------------------
!
!MODULE: ModIoUnit - general utilities for Fortran I/O units.
!
!DESCRIPTION:
!
! This module provides various utilities related to Fortran I/O units.
! In particular independently developped components can use the 
! io\_unit\_new() function to obtain an unused IO unit for extended use.
! 
! The unit number in UNITTMP\_ is a safe unit number to open and close a file 
! if no other file is opened between the open and close and all programs
! use ModIoUnit to obtain unit numbers.
!
! Standard output has unit number STDOUT\_=6. This constant is easier to read.
!
! The io\_unit\_clean subroutine closes all open IO units and deletes the
! empty files.
!
! The methods in this module can be tested by running the 
! io\_unit\_test subroutine.
!
!INTERFACE:

module ModIoUnit

  implicit none

  private ! except

  !PUBLIC MEMBER FUNCTIONS:

  public :: io_unit_new    ! Return an unused unit number for extended use
  public :: io_unit_clean  ! Close open units, delete empty files
  public :: io_unit_test   ! Test the functionality of this module

  !PUBLIC DATA MEMBERS:

  public :: UNITTMP_       ! For open read/write close without intervening open
  public :: STDOUT_        ! Fortran unit number for standard output

  integer, parameter :: STDOUT_       = 6     ! Standard output
  integer, parameter :: UNITTMP_      = 9     ! Temporary unit number

  !LOCAL VARIABLES:

  integer, parameter :: MinUnitNumber = 10    ! Smallest allowed unit number
  integer, parameter :: MaxUnitNumber = 1000  ! Largest allowed unit number

  integer :: iUnitMax = UNITTMP_              ! The largest unit number used

  !REVISION HISTORY:
  ! 01Aug03  Gabor Toth <gtoth@umich.edu> - initial prototype/prolog/code
  !EOP ___________________________________________________________________

contains

  function io_unit_new()  result(iUnit)

    !  Returns a unit number of a unit that exists and is not connected

    integer :: iUnit
    logical :: IsExisting, IsOpened
    integer :: iError

    do iUnit = MinUnitNumber, MaxUnitNumber
       inquire (&
            unit   = iUnit,       &
            exist  = IsExisting,  &
            opened = IsOpened,    &
            iostat = iError)
       if (IsExisting .and. .not. IsOpened .and. iError == 0) then
          iUnitMax = max(iUnitMax, iUnit)
          return
       end if
    end do

    iUnit = -1

  end function io_unit_new
  !===========================================================================
  subroutine io_unit_clean

    ! Close all open units for this processor
    integer :: iUnit, iError
    logical :: IsOpen
    character(len=100) :: Name
    character :: String
    !------------------------------------------------------------------------
    do iUnit = UNITTMP_,iUnitMax

       inquire(iUnit,OPENED=IsOpen,NAME=Name)
       if(IsOpen)then
          ! Clos file so that output is flushed
          close(iUnit)
          ! Try to open file and read 1 character
          open(iUnit,FILE=Name,STATUS='old',IOSTAT=iError)
          if(iError/=0) CYCLE
          read(iUnit,'(a1)',IOSTAT=iError)String
          if(iError<0)then
             ! Delete empty files
             close(iUnit,STATUS='delete')
          else
             ! Close file again
             close(iUnit)
          end if
       end if
    end do

  end subroutine io_unit_clean
  !==========================================================================
  subroutine io_unit_test

    integer :: iUnit
    logical :: IsExisting
    !---------------------------------------------------------------------

    iUnit = io_unit_new()
    write(*,*)'iUnit=',iUnit
    open(iUnit,file='baba1',status='unknown',form='formatted')
    write(iUnit,*)1

    iUnit = io_unit_new()
    write(*,*)'iUnit=',iUnit
    open(iUnit,file='baba2',status='unknown',form='unformatted')
    write(iUnit)1

    iUnit = io_unit_new()
    write(*,*)'iUnit=',iUnit
    open(iUnit,file='empty1',status='unknown',form='formatted')

    iUnit = io_unit_new()
    write(*,*)'iUnit=',iUnit
    open(iUnit,file='empty2',status='unknown',form='unformatted')

    call io_unit_clean

    inquire(file='baba1',exist=IsExisting)
    if(.not.IsExisting)write(*,*)'test baba1 failed'

    inquire(file='baba2',exist=IsExisting)
    if(.not.IsExisting)write(*,*)'test baba2 failed'

    inquire(file='empty1',exist=IsExisting)
    if(IsExisting)write(*,*)'test empty1 failed'

    inquire(file='empty2',exist=IsExisting)
    if(IsExisting)write(*,*)'test empty2 failed'

    write(*,*)'io_unit_test done'

  end subroutine io_unit_test

end module ModIoUnit
