!BOP
!
! !DESCRIPTION:
!  This is the SWMF Gridded Component, which acts as an interface to the SWMF.
!
!\begin{verbatim}

module SwmfGridCompMod

  ! ESMF Framework module
  use ESMF_Mod

  ! Named indexes for integer time arrays
  use ModTimeArray

  implicit none
  private

  public SetServices

contains

  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp) :: gcomp
    integer :: rc

    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_SETINIT, my_init, &
         ESMF_SINGLEPHASE, rc)
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_SETRUN, my_run, &
         ESMF_SINGLEPHASE, rc)
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_SETFINAL, my_final, &
         ESMF_SINGLEPHASE, rc)

  end subroutine SetServices

  !===========================================================================

  subroutine my_init(gcomp, importState, exportState, externalclock, rc)
    type(ESMF_GridComp) :: gcomp
    type(ESMF_State) :: importState
    type(ESMF_State) :: exportState
    type(ESMF_Clock) :: externalclock
    integer          :: rc

    logical          :: IsLastSession ! true if SWMF has a single session
    type(ESMF_VM)    :: vm
    integer          :: iComm, iProc
    type(ESMF_Time)  :: StartTime
    integer          :: iStartTime_I(Year_:Millisec_)
    type(ESMF_TimeInterval) :: SimTime, RunDuration
    integer(ESMF_KIND_I4)   :: iSecond, iMilliSec
    real(ESMF_KIND_R8)      :: TimeSim, TimeStop
    !------------------------------------------------------------------------
    ! Obtain the VM for the SWMF gridded component
    call ESMF_GridCompGet(gcomp,vm=vm)

    ! Obtain the MPI communicator for the VM
    call ESMF_VMGet(vm, mpiCommunicator=iComm)

    ! Obtain the start time from the clock 
    call ESMF_ClockGet(externalclock, &
         startTime=StartTime, currSimTime=SimTime, runDuration=RunDuration)
    call ESMF_TimeGet(StartTime,   &
         yy=iStartTime_I(Year_),   &
         mm=iStartTime_I(Month_),  &
         dd=iStartTime_I(Day_),    &
         h =iStartTime_I(Hour_),   &
         m =iStartTime_I(Minute_), &
         s =iStartTime_I(Second_), &
         ms=iStartTime_I(Millisec_))

    ! Obtain the simulation time from the clock
    call ESMF_TimeIntervalGet(SimTime, s=iSecond, ms=iMillisec)
    TimeSim = iSecond + iMillisec/1000.0

    ! Obtain the final simulation time from the clock
    call ESMF_TimeIntervalGet(RunDuration, s=iSecond, ms=iMillisec)
    TimeStop = iSecond + iMillisec/1000.0

    ! Initialze the SWMF with this MPI communicator and start time
    call ESMF_LogWrite("SWMF_initialize routine called", ESMF_LOG_INFO)
    call SWMF_initialize(iComm, iStartTime_I, &
         TimeSim, TimeStop, IsLastSession, rc)
    call ESMF_LogWrite("SWMF_initialize routine returned", ESMF_LOG_INFO)
    if(rc /= 0)then
       call ESMF_LogWrite("SWMF_initialize FAILED", ESMF_LOG_ERROR)
       call ESMF_VMGet(vm, localPET=iProc)
       if(iProc == 0)write(0, *) "SWMF_initialize FAILED"
       rc = ESMF_FAILURE
    endif

  end subroutine my_init

  !===========================================================================

  subroutine my_run(gcomp, importState, exportState, externalclock, rc)
    type(ESMF_GridComp) :: gcomp
    type(ESMF_State) :: importState
    type(ESMF_State) :: exportState
    type(ESMF_Clock) :: externalclock
    integer :: rc

    logical          :: DoStop ! true if SWMF requests a stop
    type(ESMF_VM)    :: vm
    integer          :: iProc
    !------------------------------------------------------------------------
    
    call ESMF_LogWrite("SWMF_run routine called", ESMF_LOG_INFO)
    call SWMF_run(DoStop, rc)
    call ESMF_LogWrite("SWMF_run routine returned", ESMF_LOG_INFO)
    if(rc /= 0)then
       call ESMF_LogWrite("SWMF_run FAILED", ESMF_LOG_ERROR)
       call ESMF_VMGet(vm, localPET=iProc)
       if(iProc == 0)write(0, *) "SWMF_run FAILED"
       rc = ESMF_FAILURE
    endif

  end subroutine my_run

  !===========================================================================

  subroutine my_final(gcomp, importState, exportState, externalclock, rc)
    type(ESMF_GridComp) :: gcomp
    type(ESMF_State) :: importState
    type(ESMF_State) :: exportState
    type(ESMF_Clock) :: externalclock
    integer :: rc

    type(ESMF_VM)    :: vm
    integer          :: iProc
    !------------------------------------------------------------------------
    call ESMF_LogWrite("SWMF_finalize routine called", ESMF_LOG_INFO)
    call SWMF_finalize(rc)
    call ESMF_LogWrite("SWMF_finalize routine returned", ESMF_LOG_INFO)
    if(rc /= 0)then
       call ESMF_LogWrite("SWMF_finalize FAILED", ESMF_LOG_ERROR)
       call ESMF_VMGet(vm, localPET=iProc)
       if(iProc == 0)write(0, *) "SWMF_finalize FAILED"
       rc = ESMF_FAILURE
    endif

  end subroutine my_final

end module SwmfGridCompMod

!\end{verbatim}

