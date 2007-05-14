!BOP
!
! !DESCRIPTION:
! Code for the ESMF_SWMF Gridded Component which creates 3 child Components:
!  an ESMF and an SWMF Gridded Component which perform a computation 
!  and a Coupler component which mediates the data exchange between them.
!
!\begin{verbatim}

module ESMF_SWMF_GridCompMod

  ! ESMF Framework module
  use ESMF_Mod

  ! User Component registration routines
  use SwmfGridCompMod,    only : SWMF_SetServices    => SetServices
  use EsmfGridCompMod,    only : ESMF_SetServices    => SetServices
  use EsmfSwmfCouplerMod, only : Coupler_SetServices => SetServices

  implicit none
  private

  public ESMF_SWMF_SetServices

  type(ESMF_GridComp), save :: SwmfComp, SwmfRootComp, EsmfComp
  type(ESMF_CplComp),  save :: CouplerComp
  type(ESMF_State),    save :: SwmfImport, EsmfExport

contains
  !============================================================================
  subroutine ESMF_SWMF_SetServices(gcomp, rc)
    type(ESMF_GridComp) :: gcomp
    integer :: rc
    !-------------------------------------------------------------------------

    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_SETINIT, my_init, &
         ESMF_SINGLEPHASE, rc)
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_SETRUN, my_run, &
         ESMF_SINGLEPHASE, rc)
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_SETFINAL, my_final, &
         ESMF_SINGLEPHASE, rc)

  end subroutine ESMF_SWMF_SetServices

  !============================================================================

  subroutine my_init(gComp, ImportState, ExportState, ParentClock, rc)

    use ESMF_SWMF_Mod, only: iMax, jMax, yMin, yMax, zMin, zMax, &
         iProcRootSwmf, iProcRootEsmf, iProcLastEsmf, iProcLastSwmf

    type(ESMF_GridComp), intent(inout)  :: gComp
    type(ESMF_State),    intent(in)  :: ImportState
    type(ESMF_State),    intent(in)  :: ExportState
    type(ESMF_Clock),    intent(inout)  :: Parentclock
    integer,             intent(out) :: rc

    type(ESMF_VM)      :: ParentVM
    type(ESMF_Grid)    :: SwmfGrid, EsmfGrid
    type(ESMF_DELayout):: Layout

    integer :: i, iProc, nProc

    !-------------------------------------------------------------------------

    call ESMF_LogWrite("ESMF_SWMF Initialize routine called", ESMF_LOG_INFO)
    rc = ESMF_FAILURE

    ! Get the layout associated with this component
    call ESMF_GridCompGet(gComp, vm=parentvm, rc=rc)
    if(rc /= ESMF_SUCCESS) return

    ! Get the number of processors
    call ESMF_VMGet(parentvm, petCount=nProc, localPET=iProc, rc=rc)
    if(rc /= ESMF_SUCCESS) return

    ! Create two grids
    SwmfGrid = ESMF_GridCreateHorzXYUni(counts=(/iMax, jMax/), &
         minGlobalCoordPerDim=(/yMin, zMin/), &
         maxGlobalCoordPerDim=(/yMax, zMax/), &
         name="SWMF grid", rc=rc)
    if(rc /= ESMF_SUCCESS) return

    EsmfGrid = ESMF_GridCreateHorzXYUni(counts=(/iMax, jMax/), &
         minGlobalCoordPerDim=(/yMin, zMin/), &
         maxGlobalCoordPerDim=(/yMax, zMax/), &
         name="ESMF grid", rc=rc)
    if(rc /= ESMF_SUCCESS) return

    ! Create the SWMF Gridded component
    SwmfComp = ESMF_GridCompCreate(parentvm, name="SWMF Gridded Component", & 
         petlist = (/ (i, i=iProcRootSwmf, iProcLastSwmf) /), &
         grid=SwmfGrid, rc=rc)
    if(rc /= ESMF_SUCCESS) return

    ! Create the ESMF Gridded component(s, there could be more than one here)
    EsmfComp = ESMF_GridCompCreate(parentvm, name="ESMF Gridded Component", &
         petlist = (/ (i, i=iProcRootEsmf, iProcLastEsmf) /), &
         grid=EsmfGrid, rc=rc)
    if(rc /= ESMF_SUCCESS) return

    ! Create the Coupler component
    CouplerComp = ESMF_CplCompCreate(parentvm, &
         name="ESMF-SWMF Coupler Component", rc=rc)
    if(rc /= ESMF_SUCCESS) return

    call ESMF_LogWrite("Component Creates finished", ESMF_LOG_INFO)

    ! Call the SetServices routine for each so they can register their
    ! subroutines for Init, Run, and Finalize
    call ESMF_GridCompSetServices(SwmfComp,  SWMF_SetServices, rc)
    if(rc /= ESMF_SUCCESS) return
    call ESMF_GridCompSetServices(EsmfComp,  ESMF_SetServices, rc)
    if(rc /= ESMF_SUCCESS) return
    call ESMF_CplCompSetServices(CouplerComp, Coupler_SetServices, rc)
    if(rc /= ESMF_SUCCESS) return

    ! Create Import and Export State objects in order to pass data
    ! between the Coupler and the Gridded Components

    SwmfImport = ESMF_StateCreate("SWMF Import", ESMF_STATE_IMPORT, rc=rc)
    if(rc /= ESMF_SUCCESS) return
    EsmfExport = ESMF_StateCreate("ESMF Export", ESMF_STATE_EXPORT, rc=rc)
    if(rc /= ESMF_SUCCESS) return

    ! Each of the subcomponents initialize themselves.
    call ESMF_GridCompInitialize(EsmfComp, exportState = EsmfExport, &
         clock=parentclock, rc=rc)
    if(rc /= ESMF_SUCCESS) RETURN

    call ESMF_GridCompInitialize(SwmfComp, importState = SwmfImport, &
         clock=parentclock, rc=rc)
    if(rc /= ESMF_SUCCESS) RETURN

    call ESMF_CplCompInitialize(CouplerComp, &
         importstate = EsmfExport, exportstate = SwmfImport, &
         clock=parentclock, rc=rc)
    if(rc /= ESMF_SUCCESS) RETURN

    rc = ESMF_SUCCESS
    call ESMF_LogWrite("ESMF-SWMF Grid Component Initialize finished", &
         ESMF_LOG_INFO)

  end subroutine my_init

  !============================================================================

  subroutine my_run(gcomp, importState, exportState, parentclock, rc)
    type(ESMF_GridComp), intent(inout) :: gcomp
    type(ESMF_State),    intent(inout) :: importState
    type(ESMF_State),    intent(inout) :: exportState
    type(ESMF_Clock),    intent(inout)    :: parentclock
    integer,             intent(out)   :: rc

    ! Local variables
    type(ESMF_Clock) :: localclock
    !-------------------------------------------------------------------------

    call ESMF_LogWrite("ESMF-SWMF Run routine called", ESMF_LOG_INFO)

    rc = ESMF_FAILURE

    ! make our own local copy of the clock
    localclock = ESMF_ClockCreate(parentclock, rc)
    if(rc /= ESMF_SUCCESS) RETURN

    do
       if(ESMF_ClockIsStopTime(localclock, rc)) EXIT
       if(rc /= ESMF_SUCCESS) RETURN

       ! Couple the subcomponents first so that SWMF has the input from ESMF
       call ESMF_CplCompRun(CouplerComp, EsmfExport, SwmfImport, localclock, &
            blockingflag=ESMF_NONBLOCKING, rc=rc)
       if(rc /= ESMF_SUCCESS) RETURN

       ! Run the subcomponents concurrently if possible
       call ESMF_GridCompRun(SwmfComp, importState=SwmfImport, &
            clock=localClock, blockingflag=ESMF_NONBLOCKING, rc=rc)
       if(rc /= ESMF_SUCCESS) RETURN

       call ESMF_GridCompRun(EsmfComp, exportState=EsmfExport, &
            clock=localClock, blockingflag=ESMF_NONBLOCKING, rc=rc)
       if(rc /= ESMF_SUCCESS) RETURN

       ! Advance the time
       call ESMF_ClockAdvance(localclock, rc=rc)
       if(rc /= ESMF_SUCCESS) RETURN

    end do

    call ESMF_ClockDestroy(localclock, rc)
    if(rc /= ESMF_SUCCESS) RETURN

    rc = ESMF_SUCCESS
    call ESMF_LogWrite("ESMF-SWMF Run finished", ESMF_LOG_INFO)

  end subroutine my_run

  !============================================================================

  subroutine my_final(gcomp, importState, exportState, parentclock, rc)
    type(ESMF_GridComp), intent(in) :: gcomp
    type(ESMF_State),    intent(in)  :: importState
    type(ESMF_State),    intent(in)  :: exportState
    type(ESMF_Clock),    intent(inout)  :: parentclock
    integer,             intent(out) :: rc

    integer :: iError 
    !-------------------------------------------------------------------------

    call ESMF_LogWrite("ESMF-SWMF Finalize routine called", ESMF_LOG_INFO)

    ! If something fails, try finalizing other things, but return with failure
    ! Assume success 
    iError = ESMF_SUCCESS

    ! Give each of the subcomponents and the coupler a chance to finalize 
    call ESMF_GridCompFinalize(SwmfComp, importState=SwmfImport, &
         clock=parentclock, rc=rc)
    if(rc /= ESMF_SUCCESS) iError = rc

    call ESMF_GridCompFinalize(EsmfComp, exportState=EsmfExport, &
         clock=parentClock, rc=rc)
    if(rc /= ESMF_SUCCESS) iError = rc

    call ESMF_CplCompFinalize(CouplerComp, SwmfImport, EsmfExport, &
         parentClock, rc=rc)
    if(rc /= ESMF_SUCCESS) iError = rc

    ! Now remove the Components to free up their resources
    call ESMF_GridCompDestroy(SwmfComp, rc)
    if(rc /= ESMF_SUCCESS) iError = rc

    call ESMF_GridCompDestroy(EsmfComp, rc)
    if(rc /= ESMF_SUCCESS) iError = rc

    call ESMF_CplCompDestroy(CouplerComp, rc)
    if(rc /= ESMF_SUCCESS) iError = rc

    rc = iError
    call ESMF_LogWrite( "ESMF-SWMF Finalize routine finished", ESMF_LOG_INFO)

  end subroutine my_final

end module ESMF_SWMF_GridCompMod

!\end{verbatim}

