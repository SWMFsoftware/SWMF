!The use /PATH_TO_BIN/TIME_CONV.exe <data.in>TIME.in
!Format of the data.in file (the filename is arbitrary:
!
!2047         Carrigton rotation 
!360          Carrington coordinate (0 corresponds to the end time, 360 to the start time.
program time_convert
  use ModUT
  use ModConst
  use ModTimeConvert
  use CON_geopack_internal
  implicit none
  integer::iCRIn, lLongitudeIn
  integer::iCR
  real::LongitudeR, LongitudeL, Longitude,LongitudeIn
  integer::iYearL,iMonthL,iDayL,iHourL,iMinL
  integer::iYearR,iMonthR,iDayR,iHourR,iMinR
  integer,parameter::iSec=0
  real::TimeL,TimeR, Time
  real::Frac=0.0
  real::TimeStart
  integer, dimension(7) :: iTime_I,iTimeL_I,iTimeR_I
  integer:: iTime
  !================
  read(*,*)iCRIn
  read(*,*)lLongitudeIn

  if(lLongitudeIn==360)then
     iCRIn=iCRIn-1
     lLongitudeIn=0
  end if

  LongitudeIn = lLongitudeIn * cDegToRad
  !This is an arbitrarily chosen Carrington rotation with the choice for 
  !left and right time covering the end time of the rotation
  iCR = 2008

  iTimeL_I=(/2003,10,23,13,0,0,0/)
  call time_int_to_real(iTimeL_I,TimeL)
  iTimeR_I=(/2003,10,23,14,0,0,0/)
  call time_int_to_real(iTimeR_I,TimeR)

  if(lLongitudeIn/=0)then
     TimeL=TimeL - (lLongitudeIn /360.0) * 28 * 86400
     TimeR=TimeR - (lLongitudeIn /360.0) * 27 * 86400
  end if
  call solve_time
  !\
  ! Now we have a time which corrsponds to the given Carrington Longitude in the rotation 2008
  !/
  do while(iCR/=iCRIn)
     if(iCR>iCRIn)then
        TimeR=TimeR-27*86400
        TimeL=TimeL-28*86400
        iCR = iCR-1
        call solve_time
     else
        TimeR=TimeR+28*86400
        TimeL=TimeL+27*86400
        iCR = iCR+1
        call solve_time
     end if
  end do
  write(*,'(a)')'#STARTTIME'
  do iTime=1,5
     write(*,*)iTime_I(iTime)
  end do
  write(*,'(a)')' 0'
  write(*,'(a)')' 0.0'
  write(*,*)
  write(*,'(a)')'#END'
  stop
contains
 subroutine solve_time
  call time_real_to_int(TimeL,iTimeL_I)
  call CON_recalc(iTimeL_I(1),iTimeL_I(2),iTimeL_I(3),iTimeL_I(4),iTimeL_I(5),iTimeL_I(6))
  LongitudeL=CarringtonLongitude

  call time_real_to_int(TimeR,iTimeR_I)
  call CON_recalc(iTimeR_I(1),iTimeR_I(2),iTimeR_I(3),iTimeR_I(4),iTimeR_I(5),iTimeR_I(6))
  LongitudeR=CarringtonLongitude

  if(LongitudeL<=cPi.and.LongitudeR>=cPi)then
     LongitudeR=LongitudeR-cTwoPi
     if(LongitudeIn > cPi)LongitudeIn=LongitudeIn-cTwoPi
  end if
  !write(*,*)iCR,TimeL,TimeR,LongitudeL,LongitudeR
  do while( any( iTimeL_I(1:5)/=iTimeR_I(1:5)).and.TimeR-TimeL>10.0)
     Time=0.50*(TimeL+TimeR)
     call time_real_to_int(Time,iTime_I)
     call CON_recalc(iTime_I(1),iTime_I(2),iTime_I(3),iTime_I(4),iTime_I(5),iTime_I(6))
     Longitude=CarringtonLongitude
     if(LongitudeL<=cPi.and.Longitude>=cPi)Longitude=Longitude-cTwoPi
     if(Longitude - LongitudeIn < 0.0)then
        iTimeR_I=iTime_I
        TimeR = Time
        LongitudeR = Longitude
     else
        iTimeL_I=iTime_I
        TimeL = Time
        LongitudeL = Longitude
     end if
!     write(*,*)TimeL,Time,TimeR,LongitudeL,Longitude,LongitudeR
!     write(*,*)iTimeL_I,iTimeR_I
  end do
end subroutine solve_time
end program time_convert
