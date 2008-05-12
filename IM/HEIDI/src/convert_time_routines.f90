
!integer function jday(year, mon, day) result(Julian_Day)
!
!  implicit none
!
!  integer :: i
!  integer, dimension(1:12) :: dayofmon
!  integer :: year, mon, day
!
!  dayofmon(1) = 31
!  dayofmon(2) = 28
!  dayofmon(3) = 31
!  dayofmon(4) = 30
!  dayofmon(5) = 31
!  dayofmon(6) = 30
!  dayofmon(7) = 31
!  dayofmon(8) = 31
!  dayofmon(9) = 30
!  dayofmon(10) = 31
!  dayofmon(11) = 30
!  dayofmon(12) = 31
!
!  if (mod(year,4).eq.0) dayofmon(2) = dayofmon(1) + 1
!  Julian_Day = 0
!  do i = 1, mon-1
!     Julian_Day = Julian_Day + dayofmon(i)
!  enddo
!  Julian_Day = Julian_Day + day
!
!end


subroutine increment_real_world_time(dt_in)

  use ModMain

  implicit none
  
  real, intent(in) :: dt_in
  integer :: carry
  logical :: notdone

  integer, dimension(1:12) :: dayofmon

  real*8, external :: Time_Of_Year

  dayofmon(1) = 31
  dayofmon(2) = 28
  dayofmon(3) = 31
  dayofmon(4) = 30
  dayofmon(5) = 31
  dayofmon(6) = 30
  dayofmon(7) = 31
  dayofmon(8) = 31
  dayofmon(9) = 30
  dayofmon(10) = 31
  dayofmon(11) = 30
  dayofmon(12) = 31

  if (mod(Start_Time_Array(1),4).eq.0) dayofmon(2) = dayofmon(1) + 1

  Time_Array = Start_Time_Array

! Milliseconds

  Time_Array(7) = Time_Array(7) + (dt_in-int(dt_in))*1000

! Seconds

  Time_Array(6) = Time_Array(6) + int(dt_in)

! Now, carry everything.

! Milliseconds

!  carry = Time_Array(7)/1000
!  Time_Array(7) = mod(Time_Array(7),1000)

!  Time_Array(6) = Time_Array(6) + carry

! Seconds

  carry = Time_Array(6)/60
  Time_Array(6) = mod(Time_Array(6),60)

  Time_Array(5) = Time_Array(5) + carry

! Minutes

  carry = Time_Array(5)/60
  Time_Array(5) = mod(Time_Array(5),60)

  Time_Array(4) = Time_Array(4) + carry

! Hours

  carry = Time_Array(4)/24
  Time_Array(4) = mod(Time_Array(4),24)

  Time_Array(3) = Time_Array(3) + carry

! Days

  notdone = .true.

  do while (notdone)

     if (Time_Array(2) > 12) then
        if (mod(Time_Array(1),4).eq.0) dayofmon(2) = dayofmon(1) - 1
        Time_Array(1) = Time_Array(1) + 1
        if (mod(Time_Array(1),4).eq.0) dayofmon(2) = dayofmon(1) + 1
        Time_Array(2) = Time_Array(2) - 12
     endif

     if (Time_Array(3) > dayofmon(Time_Array(2))) then

        Time_Array(3) = Time_Array(3) - dayofmon(Time_Array(2))
        Time_Array(2) = Time_Array(2) + 1

     else

        notdone = .false.

     endif

  enddo

  Real_Time_of_Year = Time_Of_Year(Time_Array)

  return

end subroutine increment_real_world_time


real*8 function Time_Of_Year(time_array) result(current_time)

  implicit none
  
  integer, dimension(1:7)  :: time_array
  integer :: tmp
  integer, external :: jday

  current_time = time_array(7)/1000.0
  current_time = current_time + time_array(6)
  current_time = current_time + time_array(5) * 60.0
  current_time = current_time + time_array(4) * 60.0 * 60.0

  tmp = jday(time_array(1),time_array(2),time_array(3))
  current_time = current_time + tmp * 60.0 * 60.0 * 24.0

end function Time_Of_Year
