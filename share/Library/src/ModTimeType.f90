!^CFG COPYRIGHT UM
!BOP =========================================================================
!
!MODULE: ModTimeType - defines a type for storing date and time information
!
!DESCRIPTION:
!
! Define a type for storing date and time information in various formats.
! One format stores, year, month, day, hour, minute, second as integers,
! and FracSecond as a double precision real.
! Another format uses one double precision real to count the number
! of seconds since some base time.
! Finally the year...second information is provided in a string padded
! with 0-s, so it can be used in file names.
!
! The conversion between various formats is handled by CON\_time.
!
!INTERFACE:
module ModTimeType
  !
  !USES:
  !
  use ModKind

  implicit none

  !PUBLIC TYPES:
  type TimeType
     integer           :: iYear
     integer           :: iMonth
     integer           :: iDay
     integer           :: iHour
     integer           :: iMinute
     integer           :: iSecond
     real(Real8_)      :: FracSecond
     real(Real8_)      :: Time         ! time in seconds since base time
     character(len=14) :: String       ! string with year...second.
  end type TimeType

  !EOP

end module ModTimeType
