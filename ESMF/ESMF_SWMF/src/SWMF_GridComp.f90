!BOP
!
! !DESCRIPTION:
!  This is the SWMF Gridded Component, which acts as an interface to the SWMF.
!
!\begin{verbatim}

module SwmfGridCompMod

  ! ESMF Framework module
  use ESMF_Mod

  ! Named indexes for integer time arrays and access to MHD data
  use ESMF_SWMF_Mod, only: Year_,Month_,Day_,Hour_,Minute_,Second_,MilliSec_,&
       add_mhd_fields

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
    call ESMF_LogWrite("SWMF_GridComp:init routine called", ESMF_LOG_INFO)
    rc = ESMF_FAILURE

    ! Add MHD fields to the SWMF import state
    call add_mhd_fields(gComp, importState, -1.0, rc)
    if(rc /= ESMF_SUCCESS) RETURN

    ! Obtain the VM for the SWMF gridded component
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if(rc /= ESMF_SUCCESS) RETURN

    ! Obtain the MPI communicator for the VM
    call ESMF_VMGet(vm, mpiCommunicator=iComm, rc=rc)
    if(rc /= ESMF_SUCCESS) RETURN

    ! Obtain the start time from the clock 
    call ESMF_ClockGet(externalclock, startTime=StartTime, &
         currSimTime=SimTime, runDuration=RunDuration, rc=rc)
    if(rc /= ESMF_SUCCESS) RETURN

    call ESMF_TimeGet(StartTime,   &
         yy=iStartTime_I(Year_),   &
         mm=iStartTime_I(Month_),  &
         dd=iStartTime_I(Day_),    &
         h =iStartTime_I(Hour_),   &
         m =iStartTime_I(Minute_), &
         s =iStartTime_I(Second_), &
         ms=iStartTime_I(Millisec_), &
         rc=rc)
    if(rc /= ESMF_SUCCESS) RETURN

    ! Obtain the simulation time from the clock
    call ESMF_TimeIntervalGet(SimTime, s=iSecond, ms=iMillisec, rc=rc)
    if(rc /= ESMF_SUCCESS) RETURN
    TimeSim = iSecond + iMillisec/1000.0

    ! Obtain the final simulation time from the clock
    call ESMF_TimeIntervalGet(RunDuration, s=iSecond, ms=iMillisec, rc=rc)
    if(rc /= ESMF_SUCCESS) RETURN
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
       RETURN
    endif

    rc = ESMF_SUCCESS
    call ESMF_LogWrite("SWMF_GridComp:init routine returned", ESMF_LOG_INFO)

  end subroutine my_init

  !===========================================================================

  subroutine my_run(gComp, importState, exportState, clock, rc)

    use ESMF_SWMF_Mod, ONLY: NameSwmfComp, DoBlockAllSwmf, iProcCoupleSwmf, &
         NameField_V, nVar, iMax, jMax, yMin, yMax, zMin, zMax

    type(ESMF_GridComp), intent(inout) :: gComp
    type(ESMF_State),    intent(in) :: ImportState
    type(ESMF_State),    intent(in) :: ExportState
    type(ESMF_Clock),    intent(in) :: Clock
    integer,             intent(out):: rc

    ! Access to the MHD data
    real(ESMF_KIND_R8), pointer     :: Ptr(:,:)
    real(ESMF_KIND_R8), allocatable :: Mhd_VII(:,:,:)
    integer                         :: iVar

    ! Access to time
    type(ESMF_TimeInterval) :: SimTime, TimeStep
    integer(ESMF_KIND_I4)   :: iSec, iMilliSec

    ! Parameters for the SWMF_run interface
    logical            :: DoStop          ! true if SWMF requests a stop
    real(ESMF_KIND_R8) :: tCouple         ! Coupling time
    real(ESMF_KIND_R8) :: tSimSwmf        ! SWMF Simulation time

    ! Misc variables
    type(ESMF_VM)      :: vm
    integer            :: iProc
    !------------------------------------------------------------------------
    call ESMF_LogWrite("SWMF_GridComp:run routine called", ESMF_LOG_INFO)
    rc = ESMF_FAILURE

    ! Get processor rank
    call ESMF_GridCompGet(gComp, vm=vm, rc=rc)
    if(rc /= ESMF_SUCCESS) call my_error('ESMF_GridCompGet failed')
    call ESMF_VMGet(vm, localPET=iProc, rc=rc)
    if(rc /= ESMF_SUCCESS) call my_error('ESMF_VMGet failed')

    ! Obtain pointer to the MHD data obtained from the ESMF component
    allocate(Mhd_VII(nVar, iMax, jMax), stat=rc)
    if(rc /= 0) call my_error('allocate(Mhd_VII) failed')

    if(iProc == iProcCoupleSwmf)then
       ! Copy fields into an array
       do iVar = 1, nVar
          nullify(Ptr)
          call ESMF_StateGetDataPointer(ImportState, NameField_V(iVar), Ptr, &
               rc=rc)
          if(rc/=ESMF_SUCCESS) call my_error('ESMF_StateGetDataPointer failed')
          Mhd_VII(iVar,:,:) = Ptr
       end do
       !write(*,*)'SWMF_GridComp shape of Ptr=',shape(Ptr)
       !write(*,*)'SWMF_GridComp value of Mhd=',Mhd_VII(:,1,1)
    end if

    ! Send MHD data to the GM processors in the SWMF
    !write(*,*)'!!! SWMF_GridComp SWMF_couple Mhd=',Mhd_VII(:,1,1)
    call SWMF_couple('ESMF_IH', NameSwmfComp, 'GSM', &
         nVar, iMax, jMax, yMin, yMax, zMin, zMax, Mhd_VII, rc)
    if(rc /= 0)call my_error('SWMF_couple failed')

    deallocate(Mhd_VII)

    ! Get the next coupling time from the clock
    call ESMF_ClockGet(clock, CurrSimTime=SimTime, TimeStep=TimeStep,rc=rc)
    if(rc /= ESMF_SUCCESS) call my_error('ESMF_ClockGet failed')

    ! Calculate simulation time for next coupling
    SimTime = SimTime + TimeStep
    call ESMF_TimeIntervalGet(SimTime, s=iSec, ms=iMilliSec, rc=rc)
    if(rc /= ESMF_SUCCESS) call my_error('ESMF_TimeIntervalGet failed')
    tCouple = iSec + 0.001*iMilliSec

    call ESMF_LogWrite("SWMF_run routine called!", ESMF_LOG_INFO)
    write(*,*)'SWMF_run starts  with tCouple =',tCouple
    if(DoBlockAllSwmf)then
       call SWMF_run('**', tCouple, tSimSwmf, DoStop, rc)
    else
       call SWMF_run(NameSwmfComp, tCouple, tSimSwmf, DoStop, rc)
    end if
    write(*,*)'SWMF_run returns with tSimSwmf=',tSimSwmf
    call ESMF_LogWrite("SWMF_run routine returned!", ESMF_LOG_INFO)
    if(rc /= 0)call my_error('SWMF_run failed')

    rc = ESMF_SUCCESS
    call ESMF_LogWrite("SWMF_GridComp:run routine returned", ESMF_LOG_INFO)

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

  !===========================================================================

  subroutine my_error(String)

    ! Since the error flag is not returned from my_run due to the 
    ! non-blocking flag, we have to do something drastic here

    character(len=*), intent(in) :: String

    write(*,*)'ERROR in SwmfGridCompMod:run: ',String
    
    call ESMF_Finalize

  end subroutine my_error

end module SwmfGridCompMod

!\end{verbatim}

