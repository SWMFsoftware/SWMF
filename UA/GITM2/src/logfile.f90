
subroutine logfile(dir)

  use ModGITM
  use ModInputs
  use ModTime
  use ModMpi

  implicit none

  character (len=*), intent(in) :: dir
  character (len=8) :: cIter

  real    :: minTemp, maxTemp, localVar
  integer :: iError

  if (.not. IsOpenLogFile .and. iProc == 0) then

     IsOpenLogFile = .true.
     call CON_io_unit_new(iLogFileUnit_)

     write(cIter,"(i8.8)") iStep

     open(unit=iLogFileUnit_, &
          file=dir//"/log"//cIter//".dat",status="unknown")

     write(iLogFileUnit_,'(a)') &
          "   iStep yyyy mm dd hh mm ss  ms  min(t) min(t)"

  endif

  localVar = minval(temperature(1:nLons,1:nLats,1:nAlts,1:nBlocks))
  call MPI_AllREDUCE(localVar, minTemp, 1, MPI_REAL, MPI_MIN, &
       iCommGITM, iError)

  localVar = maxval(temperature(1:nLons,1:nLats,1:nAlts,1:nBlocks))
  call MPI_AllREDUCE(localVar, maxTemp, 1, MPI_REAL, MPI_MAX, &
       iCommGITM, iError)

  if (iProc == 0) then
     minTemp = minTemp*TempUnit(1,1,1)
     maxTemp = maxTemp*TempUnit(1,1,nalts)
     write(iLogFileUnit_,"(i8,i5,5i3,i4,2f7.1)") &
          iStep, iTimeArray, minTemp, maxTemp
  endif

end subroutine logfile
