
subroutine get_log_info(GlobalMinTemp, GlobalMaxTemp, &
     GlobalMinVertVel, GlobalMaxVertVel, AverageTemp, AverageVertVel)

  use ModGITM

  real, intent(out) :: GlobalMinTemp, GlobalMaxTemp
  real, intent(out) :: GlobalMinVertVel, GlobalMaxVertVel
  real, intent(out) :: AverageTemp, AverageVertVel

  integer :: iBlock
  real    :: CellNumber
  !--------------------------------------------------------------------------

  GlobalMaxTemp    = 0.0
  GlobalMinTemp    = 1.0e32
  GlobalMaxVertVel = 0.0
  GlobalMinVertVel = 1.0e32
  AverageTemp      = 0.0
  AverageVertVel   = 0.0

  do iBlock = 1, nBlocks

     call calc_rates(iBlock)

     GlobalMaxTemp = max(GlobalMaxTemp, &
          maxval(Temperature(1:nLons,1:nLats,1:nAlts,iBlock)* &
          TempUnit(1:nLons,1:nLats,1:nAlts)))
     GlobalMinTemp = min(GlobalMinTemp, &
          minval(Temperature(1:nLons,1:nLats,1:nAlts,iBlock)* &
          TempUnit(1:nLons,1:nLats,1:nAlts)))
     AverageTemp = AverageTemp + &
          sum(Temperature(1:nLons,1:nLats,1:nAlts,iBlock) &
          *   TempUnit(1:nLons,1:nLats,1:nAlts))

     GlobalMaxVertVel = max(GlobalMaxVertVel, &
          maxval(VerticalVelocity(1:nLons,1:nLats,1:nAlts,:,iBlock)))
     GlobalMinVertVel = min(GlobalMinVertVel, &
          minval(VerticalVelocity(1:nLons,1:nLats,1:nAlts,:,iBlock)))
     AverageVertVel = AverageVertVel + &
          sum(VerticalVelocity(1:nLons,1:nLats,1:nAlts,:,iBlock))

  enddo

  CellNumber = real(nLons*nLats)*real(nAlts*nBlocks)

  AverageTemp    = AverageTemp    / CellNumber
  AverageVertVel = AverageVertVel / (CellNumber * nSpecies)

end subroutine get_log_info

!==============================================================================

subroutine logfile(dir)

  use ModGITM
  use ModInputs
  use ModTime
  use ModMpi
  use ModUtilities, ONLY: flush_unit

  implicit none

  character (len=*), intent(in) :: dir
  character (len=8) :: cIter

  real    :: minTemp, maxTemp, localVar, minVertVel, maxVertVel
  real    :: AverageTemp, AverageVertVel
  integer :: iError

  if (.not. IsOpenLogFile .and. iProc == 0) then

     IsOpenLogFile = .true.
     call CON_io_unit_new(iLogFileUnit_)

     write(cIter,"(i8.8)") iStep

     open(unit=iLogFileUnit_, &
          file=dir//"/log"//cIter//".dat",status="replace")

     write(iLogFileUnit_,'(a)') "GITM2 log file"
     write(iLogFileUnit_,'(a)') &
          "   iStep yyyy mm dd hh mm ss  ms      dt "// &
          "min(T) max(T) mean(T) min(VV) max(VV) mean(VV)"

  endif

  call get_log_info(MinTemp, MaxTemp, MinVertVel, MaxVertVel, &
       AverageTemp, AverageVertVel)

  localVar = MinTemp
  call MPI_REDUCE(localVar, minTemp, 1, MPI_REAL, MPI_MIN, &
       0, iCommGITM, iError)

  localVar = MaxTemp
  call MPI_REDUCE(localVar, maxTemp, 1, MPI_REAL, MPI_MAX, &
       0, iCommGITM, iError)

  localVar = MinVertVel
  call MPI_REDUCE(localVar, minVertVel, 1, MPI_REAL, MPI_MIN, &
       0, iCommGITM, iError)

  localVar = MaxVertVel
  call MPI_REDUCE(localVar, maxVertVel, 1, MPI_REAL, MPI_MAX, &
       0, iCommGITM, iError)

  LocalVar = AverageTemp
  call MPI_REDUCE(LocalVar, AverageTemp, 1, MPI_REAL, MPI_SUM, &
       0, iCommGITM, iError)

  LocalVar = AverageVertVel
  call MPI_REDUCE(LocalVar, AverageVertVel, 1, MPI_REAL, MPI_SUM, &
       0, iCommGITM, iError)

  if (iProc == 0) then

     AverageTemp = AverageTemp / nProcs
     AverageVertVel = AverageVertVel / nProcs

     write(iLogFileUnit_,"(i8,i5,5i3,i4,f8.4,6f13.5)") &
          iStep, iTimeArray, dt, minTemp, maxTemp, AverageTemp, &
          minVertVel, maxVertVel, AverageVertVel

     call flush_unit(iLogFileUnit_)
  endif

end subroutine logfile
