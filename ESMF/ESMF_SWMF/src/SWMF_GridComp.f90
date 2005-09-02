!BOP
!
! !DESCRIPTION:
!  This is the SWMF Gridded Component, which acts as an interface to the SWMF.
!
!\begin{verbatim}

module SwmfGridCompMod

  ! ESMF Framework module
  use ESMF_Mod

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
    integer :: iComm, iProc
    !------------------------------------------------------------------------
    ! Obtain the VM for the SWMF gridded component
    call ESMF_GridCompGet(gcomp,vm=vm)

    ! Obtain the MPI communicator for the VM
    call ESMF_VMGet(vm, mpiCommunicator=iComm)

    ! Initialze the SWMF with this MPI communicator
    call ESMF_LogWrite("SWMF_initialize routine called", ESMF_LOG_INFO)
    call SWMF_initialize(iComm, IsLastSession, rc)
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

