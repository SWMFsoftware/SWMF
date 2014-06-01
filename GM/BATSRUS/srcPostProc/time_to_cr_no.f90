!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!The use /PATH_TO_BIN/TIME_TOCRNO.exe <TIME.in>data.in
!Format of the data.in file (the filename is arbitrary:
!
!2047         Carrigton rotation 
!360          Carrington coordinate (0 corresponds to the end time, 360 to the start time.
program time_tocrno
  use ModUT
  use ModConst
  use ModTimeConvert
  use CON_geopack_internal
  implicit none
  integer::lLongitudeIn
  integer::iCR
  real::LongitudeR, LongitudeL, Longitude,LongitudeIn
  real::TimeL,TimeR, Time, TimeIn
  integer, dimension(7) :: iTime_I,iTimeL_I,iTimeR_I,iTimeIn_I
  integer:: iTime
  character(LEN=20)::TypeRead
  !================
  read(*,*)TypeRead
  if(TypeRead(1:10)/='#STARTTIME')then
     write(*,*)'Wrong input file'
     stop
  end if
  do iTime=1,6
     read(*,*)iTimeIn_I(iTime)
  end do
  iTimeIn_I(7)=0
 
  call time_int_to_real(iTimeIn_I,TimeIn)
  call CON_recalc(iTimeIn_I(1),iTimeIn_I(2),iTimeIn_I(3),iTimeIn_I(4),iTimeIn_I(5),iTimeIn_I(6))
  LongitudeIn=CarringtonLongitude
  
  !This is an arbitrarily chosen Carrington rotation with the choice for 
  !left and right time covering the end time of the rotation
  iCR = 2008

  iTimeL_I=(/2003,10,23,13,0,0,0/)
  call time_int_to_real(iTimeL_I,TimeL)
  
  iTimeR_I=(/2003,10,23,14,0,0,0/)
  call time_int_to_real(iTimeR_I,TimeR)
 
  lLongitudeIn = nint(LongitudeIn * cRadToDeg)
  if(lLongitudeIn/=0)then
     TimeL=TimeL - (lLongitudeIn /360.0) * 28 * 86400
     TimeR=TimeR - (lLongitudeIn /360.0) * 27 * 86400
  end if
  call solve_time
  !\
  ! Now we have a time which corrsponds to the given Carrington Longitude in the rotation 2008
  !/
  do while(abs(Time-TimeIn)>1.5d2)
     if(Time>TimeIn)then
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
  write(*,'(i4.4)')iCR
  if(lLongitudeIn < 10)then
     write(*,'(i1.1)')lLongitudeIn
  elseif(lLongitudeIn < 100)then
     write(*,'(i2.2)')lLongitudeIn
  else
     write(*,'(i3.3)')lLongitudeIn
  end if
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
  end if
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
  end do
end subroutine solve_time
end program time_tocrno
