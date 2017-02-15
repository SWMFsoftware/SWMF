!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine logfileDGCPM(nStepIn)

  use ModMainDGCPM,   ONLY: mgridpot
  use ModTimeDGCPM,   ONLY: StartTime, CurrentTime
  use ModIoDGCPM,     ONLY: Kp, cOutputDir
  use ModTimeConvert, ONLY: TimeType, time_real_to_int
  use ModProcPS,      ONLY: iProc
  use ModIoUnit,      ONLY: UnitTmp_

  implicit none

  integer, intent(in)          ::  nStepIn

  character(len=19), save :: StringFileTime
  
  logical, save  :: IsFirstTime = .true.
  type(TimeType) :: TimeNow

  !----------------------------------------------------------------------------

  if (iProc == 0) then
     
     ! Use TimeType to get a well-formatted date/time string.
     TimeNow%Time = CurrentTime
     call time_real_to_int(TimeNow)
     
     if (IsFirstTime) then
        ! Create time stamp for file name.
        write(StringFileTime, '(i4.4,i2.2,i2.2,"_",i2.2,i2.2,i2.2,"_",i3.3)') &
             TimeNow%iYear, TimeNow%iMonth, TimeNow%iDay, &
             TimeNow%iHour, TimeNow%iMinute, TimeNow%iSecond, &
             floor(TimeNow%FracSecond*1000.0)
          
        ! Open file and write header:
        open(unit=UnitTmp_, &
             file=cOutputDir//"PS_t"//StringFileTime//".log",status="replace")
        write(UnitTmp_,'(a)')  &
             'Dynamic Global Core Plasma Model'
        write(UnitTmp_,'(a)') &
             'i year mo dy hr mn sc msc cpp[kV] Kp'
        !else
        !   open(unit=UnitTmp_, &
        !        file=trim(dir)//"PS.log",status="old", position="append")
        !endif
        IsFirstTime = .false.
     else
        open(unit=UnitTmp_, status="old", position="append",&
             file=cOutputDir//"PS_t"//StringFileTime//".log")
     endif

     write(UnitTmp_,'(i10.10,i5,5(1x,i2.2),1x,i3.3, 2(1x,f8.3))') &
          nStepIn, TimeNow%iYear, TimeNow%iMonth,  TimeNow%iDay, &
          TimeNow%iHour,     TimeNow%iMinute, TimeNow%iSecond, &
          floor(TimeNow%FracSecond*1000.0), &
          (maxval(mgridpot) - minval(mgridpot))/1000.0, Kp
     close(UnitTmp_)

  endif

end subroutine logfileDGCPM

