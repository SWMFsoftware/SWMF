!^CFG COPYRIGHT UM
!
!QUOTE: \clearpage
!
!BOP
!
!QUOTE: \section{Time Related Variables and Methods}
!
!MODULE: CON_time - time related variables and methods 
!
!DESCRIPTION:
! The methods in this module can convert between various representations
! of time. The date can be defined as year, month, day, hour, minute, second
! and a fractional second, which is real. Alternatively an integer
! millisecond can be used. This representation is easy to read, but difficult
! to manipulate. An easy to manipulate representation is the time measured
! from some base time, namely the 00:00 GMT on January 1 of iYearBase.
! This time is measured in seconds with an 8 byte real number.
! Finally, for output and for constructing file names, a string
! representation is convenient. The string contains year, month,
! day,hour, minute, second in 14 characters. This order of the date and
! time makes alphabetic order to be the same as temporal order.
! The string is padded with 0-s so it can be used to form file names.
! The conversion routines can convert between the different fields
! of a variable of TimeType, or they can convert between integer arrays
! (year...millisecond) and an 8 byte real variable (seconds since base time).
!
! \bigskip
!
! This module also contains the time related variables of SWMF such
! as current and maximum time step and simulation time, or maximum
! CPU time allowed, etc.
!
! \bigskip
!
! The module also provides 
! a data type TypeFreq and the is\_time\_to function for an easy handling
! of actions done with some frequency. The frequency can be defined in
! terms of time steps or simulation time. The initial step or time can
! be also given, which is adjusted to the current step and time with
! the adjust\_freq subroutine.

!INTERFACE:
module CON_time

  !USES:
  use ModKind
  use ModConst
  use ModTimeType

  implicit none

  save

  !PUBLIC TYPES:

  ! Derived type to store data about the frequency of some action
  type FreqType
     logical :: DoThis     ! Do the action or not
     integer :: Dn         ! Frequency in terms of time steps
     real    :: Dt         ! Frequency in terms of seconds
     integer :: nNext      ! The next time step the action should be done
     real    :: tNext      ! The next time the action should be done
  end type FreqType

  !PUBLIC MEMBER FUNCTIONS:
  public :: init_time        ! Initialize start time
  public :: adjust_freq      ! Adjust initial step and time to current values
  public :: is_time_to       ! Returns true it is time to act and updates act

  public :: time_int_to_real ! Convert integer time info to real
  interface time_int_to_real
     module procedure time_int_to_real1, time_int_to_real2
  end interface

  public :: time_real_to_int !Convert real time info into integer
  interface time_real_to_int
     module procedure time_real_to_int1, time_real_to_int2
  end interface

  !PUBLIC DATA MEMBERS:

  ! iYearBase MUST follow : mod(iYearBase,4) == 1
  ! iYearBase MUST be AFTER 1900
  ! This particular value is required by the UA component GITM:
  integer, parameter :: iYearBase = 1965

  ! The earliest year which is already correctly handled
  integer, parameter :: iYearMin  = iYearBase

  ! The latest year which is still correctly handled
  integer, parameter :: iYearMax  = 2099

  ! The session index starting with 1
  integer :: iSession = 1

  ! Is this a time accurate run?
  logical :: DoTimeAccurate=.true.

  ! Number of time steps/iterations from the beginning
  integer :: nStep=0

  ! Number of time steps/iterations since last restart
  integer :: nIteration=0

  ! Maximum number of iterations since last restart
  integer :: MaxIteration=0

  ! Initial and current date/time
  type(TimeType) :: TimeStart, TimeCurrent

  ! Simulation time
  real :: tSimulation = 0.0

  ! Maximum simulation time
  real :: tSimulationMax = 0.0

  ! Maximum CPU time
  real(Real8_) :: CpuTimeStart    ! Initial value returned by MPI_WTIME
  real :: CpuTimeMax = -1.0       ! Maximum time allowed for the run

  ! Shall we check for the stop file?
  logical :: DoCheckStopFile = .true.

  ! How often shall we check cpu time and stop file
  ! Default is every time step
  type(FreqType):: CheckStop = FreqType(.false.,-1,-1.0,0,0.0)

  !REVISION HISTORY:
  ! 01Aug03 Aaron Ridley and G. Toth - initial implementation
  ! 22Aug03 G. Toth - added TypeFreq and is_time_to function
  ! 25Aug03 G. Toth - added adjust_freq subroutine
  !EOP

  ! February will be adjusted.....
  integer, dimension(1:12), private :: nDayInMonth_I = (/ &
       31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

  character(len=*), parameter, private :: NameMod='CON_time'

contains

  !BOP ========================================================================
  !IROUTINE: init_time - initialize start time
  !INTERFACE:
  subroutine init_time
    !EOP
    character (len=*), parameter :: NameSub = NameMod//'::init_time'

    !BOC
    TimeStart % iYear   = 2000
    TimeStart % iMonth  = 3
    TimeStart % iDay    = 21
    TimeStart % iHour   = 10
    TimeStart % iMinute = 45
    TimeStart % iSecond = 0
    TimeStart % FracSecond = 0.0

    call time_int_to_real(TimeStart)
    !EOC

  end subroutine init_time

  !BOP ========================================================================
  !IROUTINE: is_valid_int_time - check if the integer time info is valid
  !INTERFACE:
  logical function is_valid_int_time(Time)
    !INPUT/OUTPUT ARGUMENTS:
    type(TimeType), intent(inout) :: Time

    !DESCRIPTION:
    ! Check if the integer description of the time is valid.
    ! The year may be corrected, e.g. 66 --> 1966, 03 --> 2003.
    ! Return false if time is not valid.
    !EOP

    call fix_year(Time % iYear)
    is_valid_int_time = .false.
    if(Time % iYear < iYearMin .or. Time % iYear > iYearMax) RETURN
    if(Time % iMonth > 12 .or. Time % iMonth < 1) RETURN
    if(Time % iMonth == 2) call fix_february(Time % iYear)
    if(Time % iDay < 1 .or. Time % iDay > nDayInMonth_I(Time % iMonth)) RETURN
    if(Time % iHour < 0   .or. Time % iHour   > 23) RETURN
    if(Time % iMinute < 0 .or. Time % iMinute > 59) RETURN
    if(Time % iSecond < 0 .or. Time % iSecond > 59) RETURN
    if(Time % FracSecond < 0.0_Real8_ .or. Time % FracSecond > 1.0_Real8_) &
         RETURN
    is_valid_int_time = .true.

  end function is_valid_int_time

  !BOP ========================================================================
  !IROUTINE: time_int_to_string - write integer info into string
  !INTERFACE:
  subroutine time_int_to_string(Time)
    !INPUT/OUTPUT ARGUMENTS:
    type(TimeType), intent(inout) :: Time
    !DESCRIPTION:
    ! Convert integer time info into the string field of the Time variable
    !EOP
    character (len=*), parameter :: NameSub=NameMod//'::time_int_to_string'
    !-------------------------------------------------------------------------
    if(.not.is_valid_int_time(Time))then
       write(*,*)NameSub,': invalid Time = ',Time
       call CON_stop(NameSub//' SWMF_ERROR invalid time')
    end if

    write(Time % String,'(i4.4,5(i2.2))') &
         Time % iYear, Time % iMonth, Time % iDay, &
         Time % iHour, Time % iMinute, Time % iSecond

  end subroutine time_int_to_string

  !BOP =======================================================================
  !IROUTINE: time_int_to_real2: integer array into double precision time
  !INTERFACE:
  subroutine time_int_to_real2(iTime_I,Time)

    !INPUT ARGUMENTS:
    integer, intent(in) :: iTime_I(1:7)

    !OUTPUT ARGUMENTS:
    real(real8_), intent(out) :: Time
    !DESCRIPTION:
    ! Convert an integer array containing 
    ! year, month, day, hour, minute, second, millisecond
    !into the number of seconds since 00:00 January 1 of the base year.
    !EOP
    character (len=*), parameter :: NameSub = NameMod//'::time_int_to_real2'
    type(TimeType) :: TimeTmp
    !-------------------------------------------------------------------------
    TimeTmp % iYear   = iTime_I(1)
    TimeTmp % iMonth  = iTime_I(2)
    TimeTmp % iDay    = iTime_I(3)
    TimeTmp % iHour   = iTime_I(4)
    TimeTmp % iMinute = iTime_I(5)
    TimeTmp % iSecond = iTime_I(6)
    TimeTmp % FracSecond = iTime_I(7)/1000.0

    call time_int_to_real(TimeTmp)

    Time = TimeTmp % Time

  end subroutine time_int_to_real2

  !BOP ========================================================================
  !IROUTINE: time_int_to_real - convert integer time info into real time info
  !INTERFACE:
  subroutine time_int_to_real1(Time)

    !INPUT/OUTPUT ARGUMENTS:
    type(TimeType), intent(inout) :: Time

    !DESCRIPTION:
    ! Convert the integer fields containing year ... second and the
    ! fractional second into the double precision seconds counted
    ! from the beginning of the base year. Also fill in the string field
    ! of Time.
    !EOP
    character (len=*), parameter :: NameSub = NameMod//'::time_int_to_real1'

    !-------------------------------------------------------------------------
    if(.not.is_valid_int_time(Time))then
       write(*,*)NameSub,': invalid Time = ',Time
       call CON_stop(NameSub//' SWMF_ERROR invalid time')
    end if

    Time % Time = &
         ((Time%iYear-iYearBase) * 365 + n_leap_day(Time%iYear) + &
         n_day_of_year(Time%iYear, Time%iMonth, Time%iDay)-1)*cSecondPerDay + &
         Time%iHour * cSecondPerHour + &
         Time%iMinute * cSecondPerMinute + &
         Time%iSecond + &
         Time%FracSecond

    call time_int_to_string(Time)

  end subroutine time_int_to_real1

  !BOP ========================================================================
  !IROUTINE: time_real_to_int - convert the real field into the integer fields
  !INTERFACE:
  subroutine time_real_to_int1(Time)

    !INPUT/OUTPUT ARGUMENTS:
    type(TimeType), intent(inout) :: Time
    !DESCRIPTION:
    ! Convert the number of seconds counted from the beginning of the base year
    ! to the integer fields and the fractional second field. Also fill in
    ! the string field of Time.
    !EOP
    character (len=*), parameter :: NameSub = NameMod//'::time_real_to_int'

    integer :: iYear, iMonth, iDay, nLeapYear
    real(Real8_) :: TimeRemaining

    iYear = floor(Time%Time/cSecondPerYear) + iYearBase
    !write(*,*) 'iYear=',iYear
    do
       nLeapYear = n_leap_day(iYear)
       !write(*,*) 'nLeapYear=',nLeapYear
       iDay = floor((Time%Time - (iYear-iYearBase)*cSecondPerYear)/&
            cSecondPerDay) - nLeapYear
       !write(*,*)'iDay=', iDay
       if(iDay >= 0) EXIT
       iYear = iYear - 1
    end do
    !write(*,*) 'iYear, nLeapYear, is_leap_year, iDay=',&
    !     iYear,nLeapYear, is_leap_year(iYear),iDay

    TimeRemaining = Time % Time - (iYear-iYearBase) * cSecondPerYear
    TimeRemaining = TimeRemaining - (iDay+nLeapYear)*cSecondPerDay

    Time % iHour = floor(TimeRemaining/cSecondPerHour)
    TimeRemaining = TimeRemaining - Time % iHour * cSecondPerHour

    Time % iMinute = floor(TimeRemaining/cSecondPerMinute)
    TimeRemaining = TimeRemaining - Time % iMinute*cSecondPerMinute

    Time % iSecond = floor(TimeRemaining)

    Time % FracSecond = TimeRemaining - Time % iSecond

    iMonth = 1;
    call fix_february(iYear)

    do while (iDay >= nDayInMonth_I(iMonth))
       iDay = iDay - nDayInMonth_I(iMonth)
       iMonth = iMonth + 1
    end do

    Time % iYear = iYear
    Time % iMonth = iMonth
    Time % iDay = iDay + 1

    call time_int_to_string(Time)

  end subroutine time_real_to_int1

  !BOP ========================================================================
  !IROUTINE: time_real_to_int2 - double precicion seconds into int array
  !INTERFACE:
  subroutine time_real_to_int2(Time, iTime_I)
    !INPUT ARGUMENTS:
    real(real8_), intent(in) :: Time
    !OUTPUT ARGUMENTS:
    integer, intent(out) :: iTime_I(1:7)
    !DESCRIPTION:
    ! Convert the double precision number of seconds since the beginning
    ! of the base year into an integer array of year, month, day, hour,
    ! minute, second, millisecond.
    !EOP
    character (len=*), parameter :: NameSub = NameMod//'::time_real_to_int2'
    type(TimeType) :: TimeTmp
    !-------------------------------------------------------------------------
    TimeTmp % Time = Time

    call time_real_to_int(TimeTmp)

    iTime_I(1) = TimeTmp % iYear
    iTime_I(2) = TimeTmp % iMonth
    iTime_I(3) = TimeTmp % iDay
    iTime_I(4) = TimeTmp % iHour
    iTime_I(5) = TimeTmp % iMinute
    iTime_I(6) = TimeTmp % iSecond
    iTime_I(7) = TimeTmp % FracSecond * 1000.0

  end subroutine time_real_to_int2

  !BOP ========================================================================
  !IROUTINE: fix_february - fix the length of february for leap years
  !INTERFACE:
  subroutine fix_february(iYear)
    !INPUT ARGUMENTS:
    integer, intent(in) :: iYear
    !EOP
    !BOC
    if(is_leap_year(iYear))then
       nDayInMonth_I(2) = 29
    else
       nDayInMonth_I(2) = 28
    end if
    !EOC
  end subroutine fix_february

  !BOP =======================================================================
  !IROUTINE: fix_year - convert 2 digit year into 4 digit year
  !INTERFACE:
  subroutine fix_year(iYear)
    !INPUT ARGUMENTS:
    integer, intent(inout) :: iYear

    !DESCRIPTION:
    ! Attempt to fix 2 digit years. Assumption : 
    !\begin{verbatim}
    !  0-49 --> 2000-2049
    ! 50-99 --> 1950-1999
    !\end{verbatim}
    ! Using a 4 digit year is safer. You should convert before using.
    ! The 4 digit years are checked to be in the iYearMin$-$iMaxYear range
    !EOP

    character (len=*), parameter :: NameSub = NameMod//'::fix_year'

    select case(iYear)
    case(0:49)
       iYear = iYear + 2000
    case(50:99)
       iYear = iYear + 1900
    case(iYearMin:iYearMax)
       ! Fine
    case default
       write(*,*)NameSub,': iYear = ',iYear,' iYearMin = ',iYearMin,&
            'iYearMax = ',iYearMax
       call CON_stop(NameSub//' SWMF_ERROR year is out of range')
    end select

  end subroutine fix_year

  !BOP ========================================================================
  !IROUTINE: n_leap_day - number of leap days in the years since base year
  !INTERFACE:
  integer function n_leap_day(iYear) 
    !INPUT ARGUMENTS:
    integer, intent(in) :: iYear

    !DESCRIPTION:
    ! This formula ASSUMES that there are NO leap years which are
    ! the exceptions ( /100 but !/400, e.g. 2100).
    ! The leap day in iYear itself is not counted!
    !EOP
    character (len=*), parameter :: NameSub = NameMod//'::n_leap_day'
    !BOC
    n_leap_day = (iYear - iYearBase)/4
    !EOC
  end function n_leap_day

  !BOP ========================================================================
  !IROUTINE: is_leap_year - return true if the year is a leap year
  !INTERFACE:
  logical function is_leap_year(iYear)
    !INPUT ARGUMENTS:
    integer, intent(in) :: iYear
    !EOP
    !BOC
    is_leap_year = mod(iYear,4)   == 0 .and. &
         (mod(iYear,100) /= 0 .or. mod(iYear,400) == 0)
    !EOC
  end function is_leap_year

  !BOP ========================================================================
  !IROUTINE - number of days since the beginning of this year
  !INTERFACE:
  integer function n_day_of_year(iYear, iMonth, iDay)
    !INPUT ARGUMENTS:
    integer, intent(in) :: iYear, iMonth, iDay
    !DESCRIPTION:
    ! Calculate the number of days since the beginning of iYear.
    ! January 1 returns 1. Leap years are taken into account.
    !EOP
    character (len=*), parameter :: NameSub = NameMod//'::n_day_of_year'

    call fix_february(iYear)
    n_day_of_year = sum(nDayInMonth_I(1:iMonth-1)) + iDay

  end function n_day_of_year

  !BOP ========================================================================
  !IROUTINE: adjust_freq - adjust initial values to current values
  !INTERFACE:
  subroutine adjust_freq(Act,nStep,tSim)

    !INPUT/OUTPUT ARGUMENTS:
    type(FreqType), intent(inout) :: Act    ! Frequency of some action

    !INPUT ARGUMENTS:
    integer,           intent(in) :: nStep  ! Current step number
    real,              intent(in) :: tSim   ! Current simulation time

    !DESCRIPTION:
    ! Adjust the nNext and tNext fields of Act based on the
    ! current time step nStep and current simulation time tSim.
    !EOP
    !------------------------------------------------------------------------

    if(.not. Act % DoThis) RETURN

    if(Act % Dn > 0 .and. Act % nNext < nStep) &
         Act % nNext = nStep - mod(nStep - Act % nNext, Act % Dn) + Act % Dn

    if(Act % Dt > 0 .and. Act % tNext < tSim ) &
         Act % tNext = tSim  - mod(tSim  - Act % tNext, Act % Dt) + Act % Dt

  end subroutine adjust_freq
  !BOP ========================================================================
  !IROUTINE: is_time_to - is it time to do something
  !INTERFACE:
  function is_time_to(Act,nStep,tSimulation,DoCheckOnly) result(IsTimeTo)

    !INPUT/OUTPUT ARGUMENTS:
    type(FreqType), intent(inout) :: Act

    !INPUT ARGUMENTS:
    integer,           intent(in) :: nStep
    real,              intent(in) :: tSimulation
    logical, optional, intent(in) :: DoCheckOnly

    !RETURN VALUE:
    logical :: IsTimeTo

    !DESCRIPTION:
    ! Based on the next step/time info in Act and the tSimulation and 
    ! nStep values (and the DoTimeAccurate variable) decide 
    ! if Act should be done. If the answer is yes and the optional
    ! parameter DoCheckOnly is not present, then the next step/time 
    ! values in Act are increased with the step/time frequencies.
    ! If the increased values do not reach the tSimulation and
    ! nStep values then increase the next step/time relative to the
    ! values of tSimulation and nStep.
    !EOP
    !------------------------------------------------------------------------
    if(.not.Act % DoThis)then
       IsTimeTo = .false.
       RETURN
    end if

    if(Act % Dn == 0)then
       IsTimeTo = .true.
    else
       IsTimeTo = .false.

       if(DoTimeAccurate)then
          if(Act % Dt >=0) then
             IsTimeTo = tSimulation >= Act % tNext
          elseif(Act % Dn > 0)then
             IsTimeTo = nStep >= Act % nNext
          end if
       else
          if(Act % Dn > 0)then
             IsTimeTo = nStep >= Act % nNext
          end if
       end if
    end if

    if(IsTimeTo .and. .not.present(DoCheckOnly))then

       if(Act % Dt > 0)then
          Act % tNext = Act % tNext + Act % Dt
          if(Act % tNext < tSimulation) &
               Act % tNext = tSimulation + Act % Dt
       end if

       if(Act % Dn > 0) then
          Act % nNext = Act % nNext + Act % Dn
          if(Act % nNext <= nStep) &
               Act % nNext = nStep + Act % Dn
       end if
    end if

  end function is_time_to

  !BOP ========================================================================
  !IROUTINE: test_time - test the methods in this module
  !INTERFACE:
  subroutine test_time

    !LOCAL VARIABLES:
    logical, parameter :: IsVerbose = .false.
    type(FreqType)     :: Act
    integer, parameter :: nAct = 10
    logical            :: DoAct(nAct)
    logical, parameter :: F=.false., T=.true.
    !EOP
    !-------------------------------------------------------------------------
    !BOC
    write(*,*)'Testing time conversion routines'
    TimeStart % iHour = 0; TimeStart % iMinute = 0; TimeStart % iSecond = 1
    call check_all_days
    TimeStart % iHour =23; TimeStart % iMinute =59; TimeStart % iSecond =59
    call check_all_days
    write(*,'(a,i5,a,i5)')'Successfully tested all days from Jan 1',&
         iYearMin,' to Dec 31',iYearMax

    write(*,*)'Testing is_time_to function'

    DoTimeAccurate = .false.

    Act = FreqType(.false., 1,  1.0, 0, 0.0); call check_act
    if(any(DoAct))stop 'Error: DoThis=.false. should yield all false'

    Act = FreqType(.true., -1, -1.0, 0, 0.0); call check_act
    if(any(DoAct))stop 'Error: frequency -1 -1.0 should yield all false'

    Act = FreqType(.true., 0, -1.0, 2, 0.0); call check_act
    if(.not.all(DoAct))stop 'Error: frequencies 0 -1.0 should yield all true'

    Act = FreqType(.true., 2, -1.0, 2, 0.0); call check_act
    if( any(DoAct .neqv. (/F,T,F,T,F,T,F,T,F,T/)) ) &
       stop 'Error with frequencies 2 -1.0 start at 2'

    Act = FreqType(.true., 4, -1.0, -100, 0.0); call check_act
    if( any(DoAct .neqv. (/T,F,F,F,T,F,F,F,T,F/)) ) &
       stop 'Error with frequency 4 -1.0 start at -100'

    DoTimeAccurate = .true.

    Act = FreqType(.true., -1, 0.0, 0, 0.0); call check_act
    if(.not.all(DoAct))stop 'Error: frequency -1 0.0 should yield all true'

    Act = FreqType(.true., -1, 1.0, 0, 0.0); call check_act
    if(.not.all(DoAct))stop 'Error: frequency -1 1.0 should yield all true'

    Act = FreqType(.true., -1, 2.0, 0, 1.9); call check_act
    if( any(DoAct .neqv. (/F,T,F,T,F,T,F,T,F,T/)) ) &
       stop 'Error with frequency 2.0 start at 1.9'

    Act = FreqType(.true., -1, 4.0, 0, -100.0); call check_act
    if( any(DoAct .neqv. (/T,F,F,F,T,F,F,F,T,F/)) ) &
       stop 'Error with frequency 4.0 start at -100.0'

    write(*,*)'Testing adjust_freq subroutine'

    Act = FreqType(.true.,10,0.0,0,0.0)
    call adjust_freq(Act,15,0.0)
    if( Act % nNext /= 20) &
         stop 'Error with freq=10,0.0,0,0.0 adjust_freq(15,0.0)'


    Act = FreqType(.true.,10,0.0,3,0.0)
    call adjust_freq(Act,15,0.0)
    if( Act % nNext /= 23) &
         stop 'Error with freq=10,0.0,3,0.0 adjust_freq(15,0.0)'

    Act = FreqType(.true.,10,0.0,3,0.0)
    call adjust_freq(Act,0,0.0)
    if( Act % nNext /= 3) &
         stop 'Error with freq=10,0.0,0,0.0 adjust_freq(0,0.0)'

    write(*,*)'Successful'

    !EOC

  contains

    !======================================================================

    subroutine check_act
      integer :: iAct
      real    :: t
      t=0.0
      do iAct=1,nAct
         t=t+1.0           
         DoAct(iAct)=is_time_to(Act,iAct,t)
      end do
      if(IsVerbose)then
         write(*,*)'DoTimeAccurate=',DoTimeAccurate
         write(*,*)'Act  =',Act
         write(*,*)'DoAct=',DoAct
      end if
    end subroutine check_act

    !======================================================================

    subroutine check_all_days
      integer :: iYear, iMonth, iDay

      do iYear = iYearMin,iYearMax
         TimeStart % iYear = iYear
         do iMonth = 1, 12
            TimeStart % iMonth = iMonth
            call fix_february(iYear)
            do iDay = 1, nDayInMonth_I(iMonth)
               TimeStart % iDay = iDay

               ! Convert to real time
               call time_int_to_real(TimeStart)
               TimeCurrent % Time = TimeStart % Time
               ! Convert back
               call time_real_to_int(TimeCurrent)

               if(TimeCurrent % String /= TimeStart % String) then
                  write(*,*)'TimeStart  =',TimeStart % String
                  write(*,*)'TimeCurrent=', TimeCurrent % String
                  stop
               end if
            end do
         end do
      end do
    end subroutine check_all_days

  end subroutine test_time

  !============================================================================

end module CON_time
