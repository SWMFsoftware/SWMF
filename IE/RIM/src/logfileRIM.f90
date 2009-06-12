
subroutine logfileRIM(dir)

  use ModRIM,      ONLY : nSolve,StartTime,CurrentTime,TimeArray,&
       ThetaTilt,cpcpn,cpcps
  use ModProcIE,   ONLY : iProc
  use ModNumConst, ONLY: cRadToDeg
  use ModIoUnit,   ONLY: UnitTmp_

  implicit none

  character (len=*), intent(in) :: dir
  logical :: IsFirstTime = .true.
  !----------------------------------------------------------------------------

  if (iProc == 0) then

     if (IsFirstTime) then

        if (nSolve == 0 .and. StartTime == CurrentTime) then
           open(unit=UnitTmp_, &
                file=trim(dir)//"/IE.log",status="replace")
           write(UnitTmp_,'(a)')  &
                'Ridley Ionosphere Model, [deg] and [kV]'
           write(UnitTmp_,'(a)') &
                'nsolve t yy mm dd hh mm ss ms tilt cpcpn cpcps'
        endif
        IsFirstTime = .false.
     else
        open(unit=UnitTmp_, &
             file=trim(dir)//"/IE.log",status="old", position="append")
     endif

     write(UnitTmp_,'(i8,es13.5,i5,5i3,i4,f8.2,2es13.5)') &
          nSolve, CurrentTime-StartTime, TimeArray(1:7), &
          ThetaTilt*cRadToDeg, cpcpn, cpcps
     close(UnitTmp_)

  endif

end subroutine logfileRIM

