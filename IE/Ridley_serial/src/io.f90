! Copyright (c) 1994 Unicomp, Inc.  All rights reserved.
!
! Developed at Unicomp, Inc.
!
! Permission to use, copy, modify, and distribute this
! software is freely granted, provided that this notice 
! is preserved.

! NAG F90

module IO

   implicit none

!  Default input and output units:

   integer, parameter :: DEFAULT_INPUT_UNIT = 5
   integer, parameter :: DEFAULT_OUTPUT_UNIT = 6

!  Number and value of preconnected units

   integer, parameter :: NUMBER_OF_PRECONNECTED_UNITS = 3
   integer, parameter :: PRECONNECTED_UNITS (NUMBER_OF_PRECONNECTED_UNITS) = &
                         (/ 0, 5, 6 /)

!  Values returned to IOSTAT for end of record and end of file

   integer, parameter :: END_OF_RECORD = -2
   integer, parameter :: END_OF_FILE = -1

!  Largest allowed unit number (or a large number, if none)

   integer, parameter :: MAX_UNIT_NUMBER = 1000

contains

   function new_unit ()  result (result)

!  Returns a unit number of a unit that exists and is not connected
   
      integer :: result
      logical :: exists, opened
      integer :: ios
   
      do result = 1, max_unit_number
         if (result == DEFAULT_INPUT_UNIT .or. &
             result == DEFAULT_OUTPUT_UNIT) cycle
         if (any (result == PRECONNECTED_UNITS)) cycle
         inquire (unit = result,  &
                  exist = exists,  &
                  opened = opened,  &
                  iostat = ios)
         if (exists .and. .not. opened .and. ios == 0) return
      end do
   
      result = -1
   
   end function new_unit

end module IO
