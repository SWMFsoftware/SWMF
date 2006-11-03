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

  TIMELOOP:do
     if (Time >= Tmax) exit TIMELOOP
     do iLine=1,nLine

        ! MoveFluxTube moves the flux tube, then we can use the angular
        !position to get the lat and lon

        call MoveFluxTube

        !  Call the flux tube to be solved

        call AdvancePWline
     enddo
  enddo TIMELOOP


  !****************************************************************************
  !  Write output, use cartesian coords for output
  !****************************************************************************


  do iLine=1,nLine
     CLOSE(UNIT=iUnitGraphics(iLine))
  enddo
  close(UNIT=iUnitOutput)

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

