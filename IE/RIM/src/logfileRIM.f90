
subroutine logfileRIM(dir)

  use ModRIM
  use ModProcIE
  use ModUtilities, ONLY: flush_unit
  use ModNumConst

  implicit none

  character (len=*), intent(in) :: dir
  logical :: IsFirstTime = .true.
  integer :: iLogFileUnit_ = 57

  if (iProc == 0) then

     if (IsFirstTime) then

        call CON_io_unit_new(iLogFileUnit_)

        if (nSolve == 0 .and. StartTime == CurrentTime) then
           open(unit=iLogFileUnit_, &
                file=trim(dir)//"/logIE.dat",status="replace")
           write(iLogFileUnit_,'(a)')  &
                'Ridley Ionosphere Model, [deg] and [kV]'
           write(iLogFileUnit_,'(a)') &
                'nsolve t yy mm dd hh mm ss ms ttilt ptilt cpcpn cpcps'
        endif

     else

        open(unit=iLogFileUnit_, &
             file=trim(dir)//"/logIE.dat",status="old", position="append")

     endif

     write(iLogFileUnit_,'(i4,1p,e13.5,i5,5i3,i4,2f7.2,2e13.5)') &
          nSolve, CurrentTime-StartTime, TimeArray(1:7), &
          ThetaTilt*cRadToDeg, 0.0, cpcpn, cpcps
     close(iLogFileUnit_)

  endif

end subroutine logfileRIM

