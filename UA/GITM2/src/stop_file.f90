
!\
! -------------------------------------------------------------------
! -------------------------------------------------------------------
!/

subroutine check_stop

  use ModGITM
  use ModTime
  use ModInputs, only: CPUTimeMax, iOutputUnit_
  use ModMpi
  implicit none

  real, external :: get_timing

  real*8  :: EndTimeLocal
  logical :: IsThere
  integer :: iError

  call report("check_stop",2)

  EndTimeLocal = EndTime

  inquire(file="GITM.STOP",EXIST=IsThere)
  if (IsThere) then
     if (iProc == 0) write(*,*) "GITM.STOP file found. Exiting."
     EndTimeLocal = CurrentTime - 1.0
  endif

  if (get_timing("GITM") > CPUTimeMax) then
     if (iProc == 0) write(*,*) "CPUTimeMax Exceeded. Exiting."
     EndTimeLocal = CurrentTime - 1.0
  endif

  call MPI_AllREDUCE(EndTimeLocal, EndTime,  &
       1, MPI_DOUBLE_PRECISION, MPI_MIN, iCommGITM, iError)

end subroutine check_stop

!\
! -------------------------------------------------------------------
! -------------------------------------------------------------------
!/

subroutine delete_stop

  use ModGITM
  use ModInputs, only: iOutputUnit_

  implicit none

  logical :: IsThere

  call report("delete_stop",2)

  inquire(file='GITM.STOP',EXIST=IsThere)
  if (IsThere .and. iProc == 0) then
     open(iOutputUnit_, file = 'GITM.STOP', status = 'OLD')
     close(iOutputUnit_, status = 'DELETE')
  endif

  inquire(file='GITM.DONE',EXIST=IsThere)
  if (IsThere .and. iProc == 0) then
     open(iOutputUnit_, file = 'GITM.DONE', status = 'OLD')
     close(iOutputUnit_, status = 'DELETE')
  endif

end subroutine delete_stop
