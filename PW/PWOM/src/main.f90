program pw

  use ModPwom
  use ModFieldLine
  use ModMpi
  use ModReadParam
  implicit none

  !****************************************************************************
  ! Initiallize MPI and get number of processors and rank of given processor
  !****************************************************************************

  !---------------------------------------------------------------------------
  call MPI_INIT(errcode)
  iComm = MPI_COMM_WORLD

  call MPI_COMM_RANK(iComm,iProc,errcode)
  call MPI_COMM_SIZE(iComm,nProc,errcode)

  !****************************************************************************
  ! Read the input file
  !****************************************************************************
  IsStandAlone = .true.
  NameInput    = 'pw.input'

  call read_file('pw.input',iComm)
  call read_init('  ',iSessionIn=1,iLineIn=0)

  call PW_set_parameters('READ')
  ! call PW_set_parameters('CHECK')

  call PW_initialize

  !****************************************************************************
  ! Move the flux tube, solve each fieldline, and advance the time
  !****************************************************************************

  if (DoTimeAccurate) then
     TIMELOOP:do
        if (Time >= Tmax) exit TIMELOOP
        DtHorizontal = min(DtHorizontalOrig, Tmax - Time)
        if (DtHorizontal < 1.0e-6) then
           Time = Tmax
           exit TIMELOOP
        endif
        do iLine=1,nLine
           
           ! move_line moves the flux tube, then we can use the angular
           !position to get the lat and lon
           
           call move_line
           
           !  Call the flux tube to be solved
           
           call PW_advance_line
        enddo
        !Output the electrodynamics info
        if (DoPlotElectrodynamics) then
           if (floor(Time/DtPlotElectrodynamics) &
                /= floor((Time-DtHorizontal)/DtPlotElectrodynamics) ) &
                call PW_print_electrodynamics
        endif
     enddo TIMELOOP
  else
     NLOOP:do
        if (nStep >= MaxStep) exit NLOOP
        DtHorizontal = DtHorizontalOrig
        do iLine=1,nLine
           ! move_line moves the flux tube, then we can use the angular
           !position to get the lat and lon
           call move_line
           
           !  Call the flux tube to be solved
           call PW_advance_line
        enddo
        
     enddo NLOOP
  end if

  !****************************************************************************
  !  Write output, use cartesian coords for output
  !****************************************************************************

  if (nLog == -1) then
     do iLine=1,nLine
        close(iUnitOutput(iLine))
     enddo
  elseif(nLog ==0) then
     !do nothing in this case
  elseif(nLog==iLineGlobal(iLine)) then
     close(iUnitOutput(iLine))
  else
  end if

  !\
  ! Deallocate variables needed for simulation
  !/
  deallocate(r_C, State_CVI, GeoMagLat_I,GeoMagLon_I,     &
       ThetaLine_I, PhiLine_I, xLine_I, yLine_I, zLine_I, &
       xLineOld_I, yLineOld_I, zLineOld_I, UthetaLine_I,  &
       UphiLine_I, UxLine_I, UyLine_I, UzLine_I,          &
       OmegaLine_I, JrLine_I, iThetaLine_I,iPhiLine_I,    &
       NameRestartIn, NameRestart, NameGraphics,          &
       NameOutput,  iUnitRestart, iUnitRestartIn,         &
       iUnitGraphics,iUnitOutput, iLineGlobal)
  

  call MPI_FINALIZE(errcode)

end program pw

!============================================================================
! The following subroutines are here so that we can use SWMF library routines
! Also some features available in SWMF mode only require empty subroutines
! for compilation of the stand alone code.
!============================================================================
subroutine CON_stop(StringError)
  use ModPWOM, ONLY : iProc,iComm,Time
  use ModMpi
  implicit none
  character (len=*), intent(in) :: StringError

  ! Local variables:
  integer :: iError,nError
  !----------------------------------------------------------------------------

  write(*,*)'Stopping execution! me=',iProc,' at time=',Time,&
       ' with msg:'
  write(*,*)StringError
  call MPI_abort(iComm, nError, iError)
  stop

end subroutine CON_stop

subroutine CON_set_do_test(String,DoTest,DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe

  DoTest = .false.; DoTestMe = .false.

end subroutine CON_set_do_test

subroutine CON_io_unit_new(iUnit)

  use ModIoUnit, ONLY: io_unit_new
  implicit none
  integer, intent(out) :: iUnit

  iUnit = io_unit_new()

end subroutine CON_io_unit_new

