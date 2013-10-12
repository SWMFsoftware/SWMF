!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
Module ModUT
  !UT stands for the Universal Time
  implicit none
  double precision,parameter::JulianCentury=36525.000000000000000000d0,&
                    JulianDayOfEpocZeroDay=2415020.000000000000000d0
  private::JulianDay,FDAY
contains
  integer function JulianDay(iYear,iMonth,iDay)
    !Coded by A.Ridley
    !Comment: valid for 1800<iYear<2100
    integer,intent(in)::iYear,iMonth,iDay
    integer, dimension(1:12),parameter :: nDayInMonth_I = (/ &
         31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
    IF(iYear.LT.1801.OR.iYear.GT.2099)then
       write(*,*)'Does not work for the year of ',iYear
    end IF
    JulianDay=iDay
    if(iMonth>1)JulianDay=JulianDay+sum(nDayInMonth_I(1:iMonth-1))
    if(iMonth>2.and.mod(iYear,4)==0)JulianDay=JulianDay+1
  end function JulianDay
  !---------------------------------------------------------------------------
  double precision function FDAY(iHour,iMin,iSec,FracSecond)
    integer,intent(in)::iHour,iMin,iSec
    real,intent(in)::FracSecond
    optional::FracSecond,iMin,iSec

    FDAY=IHOUR*3600
    if(present(iMin))FDAY=FDAY+iMIN*60
    if(present(iSec))FDAY=FDAY+ISEC
    if(present(FracSecond))FDAY=FDAY+dble(FracSecond)
    FDAY=FDAY/86400.0000000000000D0       !Fraction of a day
  end function FDAY
  !--------------------------------------------------------------------------
  double precision function JD(iYear,iMonth,iDay,iHour,iMin,iSec,FracSecond)
    integer,intent(in)::iYear,iMonth,iDay,iHour,iMin,iSec
    real,intent(in)::FracSecond
    optional::FracSecond,iMin,iSec
    JD=day_of_epoc(iYear,iMonth,iDay,iHour,iMin,iSec,FracSecond)+&
         JulianDayOfEpocZeroDay
  end function JD
  !--------------------------------------------------------------------------
  double precision function day_of_epoc(iYear,iMonth,iDay,iHour,iMin,iSec,FracSecond)
    integer,intent(in)::iYear,iMonth,iDay,iHour,iMin,iSec
    real,intent(in)::FracSecond
    optional::FracSecond,iMin,iSec
    day_of_epoc=365*(IYear-1800)+(IYear-1801)/4+&
         JulianDay(iYear,iMonth,iDay)-&
         0.5000000000000000000000000D0+& !Offset noon-midnight
         FDAY(iHour,iMin,iSec)&
         -JulianCentury& !Offset from the previous century
         +(1-iYear/1900) !Eclude 1899,Dec,31,noon:1900,Jan,1,noon day
                         !which is added as a fake leap day 1900,Feb,29
  end function day_of_epoc
 !------------------------------------------------------------------------
  double precision function day_2000_epoc(iYear,iMonth,iDay,iHour,iMin,iSec,FracSecond)
    integer,intent(in)::iYear,iMonth,iDay,iHour,iMin,iSec
    real,intent(in)::FracSecond
    optional::FracSecond,iMin,iSec
    day_2000_epoc=day_of_epoc(iYear,iMonth,iDay,iHour,iMin,iSec,FracSecond)-&
         JulianCentury
  end function day_2000_epoc
  !------------------------------------------------------------------------
  subroutine test_ut
    write(*,'(a)')
    write(*,'(a)')'Coded by I.Sokolov.'  
    write(*,'(a)')
    write(*,'(a,F10.1,a)')'It was a long (look at it yourself)  Julian day,JD=',&
         JD(2005,3,28,12),','
    write(*,'(a,F10.1,a)')'when I mentioned,that between   Jan,1,2000,noon,JD=',&
         JD(2000,1,1,12),','
    write(*,'(a,F10.1,a)')'and                             Jan,1,1900,noon,JD=',&
         JD(1900,1,1,12),','
    write(*,'(a)')'seem to be 25 leap days, - yes, Feb,29,1900 is among them!'
    write(*,'(a,F10.1,a)')'And between Jan,1,1900,noon and Jan,1,1854,noon,JD=',&
         JD(1854,1,1,12),','
    write(*,'(a)')'seem to be 10 leap days only - among 46 years should be 11 leap ones.'
    write(*,'(a)')
    write(*,'(a)')'This is because the time Dec,31,1899,11.59.00'
    write(*,*)'is only 1 min=',&
         day_of_epoc(1899,12,31,11,59),' days prior to Jan,1,1900,noon,'
    write(*,'(a,F3.1,a)')'which is accepted to be the start of epoc 1900.',&
         day_of_epoc(1900,1,1,12),','
    write(*,'(a,F3.1,a)')'while Jan,1,2000,noon is the start of epoc 2000.',&
         day_2000_epoc(2000,1,1,12),'.'
    write(*,'(a)')
    write(*,'(a,F10.1,a)')'Carrington rotations started at Nov,9,1853,    JD=',&
         JD(1853,11,9,22),','
    write(*,'(a)')'so we need the universal time scale for this year too.'
    write(*,'(a)')'You see that we have it now.'
  end subroutine test_ut
end Module ModUT
