
subroutine get_log_info(GlobalMinTemp, GlobalMaxTemp, &
     GlobalMinVertVel, GlobalMaxVertVel)

  use ModGITM

  real, intent(out) :: GlobalMinTemp, GlobalMaxTemp
  real, intent(out) :: GlobalMinVertVel, GlobalMaxVertVel
  integer :: iBlock

  do iBlock = 1, nBlocks

     call calc_rates(iBlock)

     if (iBlock == 1) then

        GlobalMaxTemp = 0.0
        GlobalMinTemp = 1.0e32

        GlobalMaxVertVel = 0.0
        GlobalMinVertVel = 1.0e32

     end if

     GlobalMaxTemp = max(GlobalMaxTemp, &
          maxval(Temperature(1:nLons,1:nLats,1:nAlts,iBlock)* &
          TempUnit(1:nLons,1:nLats,1:nAlts)))
     GlobalMinTemp = min(GlobalMinTemp, &
          minval(Temperature(1:nLons,1:nLats,1:nAlts,iBlock)* &
          TempUnit(1:nLons,1:nLats,1:nAlts)))

     GlobalMaxVertVel = max(GlobalMaxVertVel, &
          maxval(VerticalVelocity(1:nLons,1:nLats,1:nAlts,:,iBlock)))
     GlobalMinVertVel = min(GlobalMinVertVel, &
          minval(VerticalVelocity(1:nLons,1:nLats,1:nAlts,:,iBlock)))

  enddo

end subroutine get_log_info

subroutine logfile(dir)

  use ModGITM
  use ModInputs
  use ModTime
  use ModMpi

  implicit none

  character (len=*), intent(in) :: dir
  character (len=8) :: cIter

  real    :: minTemp, maxTemp, localVar, minVertVel, maxVertVel
  integer :: iError

  if (.not. IsOpenLogFile .and. iProc == 0) then

     IsOpenLogFile = .true.
     call CON_io_unit_new(iLogFileUnit_)

     write(cIter,"(i8.8)") iStep

     open(unit=iLogFileUnit_, &
          file=dir//"/log"//cIter//".dat",status="unknown")

     write(iLogFileUnit_,'(a)') &
          "   iStep yyyy mm dd hh mm ss  ms      dt min(t) max(t)"// &
          " min(VV) max(VV)"

  endif

  call get_log_info(MinTemp, MaxTemp, MinVertVel, MaxVertVel)

  localVar = MinTemp
  call MPI_AllREDUCE(localVar, minTemp, 1, MPI_REAL, MPI_MIN, &
       iCommGITM, iError)

  localVar = MaxTemp
  call MPI_AllREDUCE(localVar, maxTemp, 1, MPI_REAL, MPI_MAX, &
       iCommGITM, iError)

  localVar = MinVertVel
  call MPI_AllREDUCE(localVar, minVertVel, 1, MPI_REAL, MPI_MIN, &
       iCommGITM, iError)

  localVar = MaxVertVel
  call MPI_AllREDUCE(localVar, maxVertVel, 1, MPI_REAL, MPI_MAX, &
       iCommGITM, iError)

  if (iProc == 0) then
     write(iLogFileUnit_,"(i8,i5,5i3,i4,f8.4,4f13.5)") &
          iStep, iTimeArray, dt, minTemp, maxTemp, minVertVel, maxVertVel
  endif

end subroutine logfile
