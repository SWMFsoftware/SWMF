!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine logfileDGCPM(dir, i3)

  use ModMainDGCPM,   ONLY: mgridpot,CrossPotential, rkph
  use ModCoupleDGCPM, ONLY: IsCoupled
  use ModTimeDGCPM,   ONLY: StartTime, CurrentTime
  use ModTimeConvert, ONLY: TimeType, time_real_to_int
  use ModProcPS,      ONLY: iProc
  use ModIoUnit,      ONLY: UnitTmp_

  implicit none

  character (len=*), intent(in) :: dir
  integer, intent(in) ::  i3
  logical :: IsFirstTime = .true.
  type(TimeType) :: TimeDgcpmNow

  real kp
  !----------------------------------------------------------------------------

  if (iProc == 0) then

     if (IsFirstTime) then
        !if (i3 == 0 .and. StartTime == CurrentTime) then
           open(unit=UnitTmp_, &
                file=trim(dir)//"PS.log",status="replace")
           write(UnitTmp_,'(a)')  &
                'Dynamic Global Core Plasma Model'
           write(UnitTmp_,'(a)') &
                'i year mo dy hr mn sc msc cpp[kV]'
        !else
        !   open(unit=UnitTmp_, &
        !        file=trim(dir)//"PS.log",status="old", position="append")
        !endif
        IsFirstTime = .false.
     else
        open(unit=UnitTmp_, &
             file=trim(dir)//"PS.log",status="old", position="append")
     endif

     ! Use TimeType to get a well-formatted date/time string.
     TimeDgcpmNow%Time = CurrentTime
     call time_real_to_int(TimeDgcpmNow)

     write(UnitTmp_,'(i10.10,i5,5i3,1x,i3.3, f8.3)') &
          i3, TimeDgcpmNow%iYear, TimeDgcpmNow%iMonth, TimeDgcpmNow%iDay, &
          TimeDgcpmNow%iHour, TimeDgcpmNow%iMinute, TimeDgcpmNow%iSecond, &
          floor(TimeDgcpmNow%FracSecond*1000.0), &
          (maxval(mgridpot) - minval(mgridpot))/1000.0
     close(UnitTmp_)

  endif

end subroutine logfileDGCPM

