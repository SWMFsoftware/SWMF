module ESMF_SWMF_GridCompMod

  ! Code for the ESMF_SWMF Gridded Component which creates 3 child Components:
  !  an ESMF and an SWMF Gridded Component which perform a computation 
  !  and a Coupler component which mediates the data exchange between them.
  
  ! ESMF Framework module
  use ESMF

  ! User Component registration routines
  use SwmfGridCompMod,    only : SwmfSetServices => SetServices
  use EsmfGridCompMod,    only : EsmfSetServices => SetServices
  use EsmfSwmfCouplerMod, only : CouplerSetServices => SetServices

  implicit none
  private

  public:: ESMF_SWMF_SetServices

  type(ESMF_GridComp), save :: SwmfComp, SwmfRootComp, EsmfComp
  type(ESMF_CplComp),  save :: CouplerComp
  type(ESMF_State),    save :: SwmfImport, EsmfExport

contains
  !============================================================================
  subroutine ESMF_SWMF_SetServices(gcomp, rc)

    type(ESMF_GridComp) :: gcomp
    integer, intent(out):: rc
    !--------------------------------------------------------------------------
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=my_init, rc=rc)
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         userRoutine=my_run, rc=rc)
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_FINALIZE, &
         userRoutine=my_final, rc=rc)

  end subroutine ESMF_SWMF_SetServices
  !============================================================================
  subroutine my_init(gComp, ImportState, ExportState, ParentClock, rc)

    use ESMF_SWMF_Mod, only: nLon, nLat, LonMin, LonMax, LatMin, LatMax, &
         iProcRootSwmf, iProcRootEsmf, iProcLastEsmf, iProcLastSwmf

    type(ESMF_GridComp):: gComp
    type(ESMF_State):: ImportState, ExportState
    type(ESMF_Clock):: Parentclock
    integer, intent(out) :: rc

    type(ESMF_VM)      :: ParentVM
    type(ESMF_Grid)    :: SwmfGrid, EsmfGrid
    type(ESMF_DELayout):: Layout

    integer :: i, iProc, nProc
    !--------------------------------------------------------------------------
    call ESMF_LogWrite("ESMF_SWMF_GridComp init called", ESMF_LOGMSG_INFO)
    rc = ESMF_FAILURE

    ! Get the layout associated with this component
    call ESMF_GridCompGet(gComp, vm=parentvm, rc=rc)
    if(rc /= ESMF_SUCCESS) call my_error('ESMF_GridCompGet')

    ! Get the number of processors
    call ESMF_VMGet(parentvm, petCount=nProc, localPET=iProc, rc=rc)
    if(rc /= ESMF_SUCCESS) call	my_error('ESMF_VMGet')

    ! Create two grids
    SwmfGrid = ESMF_GridCreate1PeriDimUfrm(maxIndex=[nLon, nLat], &
         minCornerCoord=[LonMin, LatMin], maxCornerCoord=[LonMax, LatMax], &
         staggerLocList=[ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER], &
         name="SWMF grid", rc=rc)
    if(rc /= ESMF_SUCCESS)call	my_error('ESMF_GridCreate1PeriDimUfrm Swmf')

    EsmfGrid = ESMF_GridCreate1PeriDimUfrm(maxIndex=[nLon, nLat], &
         minCornerCoord=[LonMin, LatMin], maxCornerCoord=[LonMax, LatMax], &
         staggerLocList=[ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER], &
         name="ESMF grid", rc=rc)
    if(rc /= ESMF_SUCCESS)call my_error('ESMF_GridCreate1PeriDimUfrm Esmf')

    ! Create the SWMF Gridded component
    SwmfComp = ESMF_GridCompCreate(name="SWMF Gridded Component", & 
         grid=SwmfGrid, petlist = [ (i, i=iProcRootSwmf, iProcLastSwmf) ], &
         rc=rc)
    if(rc /= ESMF_SUCCESS)call my_error('ESMF_GridCompCreate Swmf')

    ! Create the ESMF Gridded component(s, there could be more than one here)
    EsmfComp = ESMF_GridCompCreate(name="ESMF Gridded Component", &
         grid=EsmfGrid, petlist = [ (i, i=iProcRootEsmf, iProcLastEsmf) ], &
         rc=rc)
    if(rc /= ESMF_SUCCESS)call my_error('ESMF_GridCompCreate Esmf')

    ! Create the Coupler component
    CouplerComp = ESMF_CplCompCreate(name="ESMF-SWMF Coupler Component", &
         petlist = [ (i, i=0, nProc-1) ], rc=rc)
    if(rc /= ESMF_SUCCESS)call my_error('ESMF_GridCompCreate Coupler')

    call ESMF_LogWrite("Component Creates finished", ESMF_LOGMSG_INFO)

    ! Call the SetServices routine for each so they can register their
    ! subroutines for Init, Run, and Finalize
    call ESMF_GridCompSetServices(SwmfComp, &
         userRoutine=SwmfSetServices, rc=rc)
    if(rc /= ESMF_SUCCESS)call my_error('ESMF_GridCompSetServices Swmf')
    call ESMF_GridCompSetServices(EsmfComp, &
         userRoutine=EsmfSetServices, rc=rc)
    if(rc /= ESMF_SUCCESS)call my_error('ESMF_GridCompSetServices Esmf')
    call ESMF_CplCompSetServices(CouplerComp, &
         userRoutine=CouplerSetServices, rc=rc)
    if(rc /= ESMF_SUCCESS)call my_error('ESMF_GridCompSetServices Coupler')

    ! Create Import and Export State objects in order to pass data
    ! between the Coupler and the Gridded Components
    SwmfImport = ESMF_StateCreate(name="SWMF Import", &
         stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
    if(rc /= ESMF_SUCCESS)call my_error('ESMF_StateCreate SwmfImport')
    EsmfExport = ESMF_StateCreate(name="ESMF Export", &
         stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
    if(rc /= ESMF_SUCCESS)call my_error('ESMF_StateCreate EsmfExport')

    ! Each of the subcomponents initialize themselves.
    call ESMF_GridCompInitialize(EsmfComp, exportState = EsmfExport, &
         clock=parentclock, rc=rc)
    if(rc /= ESMF_SUCCESS)call my_error('ESMF_GridCompInitialize Esmf')

    call ESMF_GridCompInitialize(SwmfComp, importState = SwmfImport, &
         clock=parentclock, rc=rc)
    if(rc /= ESMF_SUCCESS)call my_error('ESMF_GridCompInitialize Swmf')

    call ESMF_CplCompInitialize(CouplerComp, &
         importstate = EsmfExport, exportstate = SwmfImport, &
         clock=parentclock, rc=rc)
    if(rc /= ESMF_SUCCESS)call my_error('ESMF_CplCompInitialize')

    rc = ESMF_SUCCESS
    call ESMF_LogWrite("ESMF-SWMF GridComp init finished", &
         ESMF_LOGMSG_INFO)

  end subroutine my_init
  !============================================================================
  subroutine my_run(gcomp, importState, exportState, parentclock, rc)

    use ESMF_SWMF_Mod, ONLY: SyncFlag

    type(ESMF_GridComp):: gcomp
    type(ESMF_State)   :: importState
    type(ESMF_State)   :: exportState
    type(ESMF_Clock)   :: parentclock
    integer, intent(out):: rc

    ! Local variables
    type(ESMF_Clock) :: localclock
    !-------------------------------------------------------------------------
    call ESMF_LogWrite("ESMF-SWMF Run routine called", ESMF_LOGMSG_INFO)

    rc = ESMF_FAILURE

    ! make our own local copy of the clock
    localclock = ESMF_ClockCreate(parentclock, rc=rc)
    if(rc /= ESMF_SUCCESS)call my_error('ESMF_ClockCreate')

    do
       if(ESMF_ClockIsStopTime(localclock, rc=rc)) EXIT
       if(rc /= ESMF_SUCCESS)call my_error('ESMF_ClockIsStopTime')

       ! Couple the subcomponents first so that SWMF has the input from ESMF
       call ESMF_CplCompRun(CouplerComp, &
            importstate = EsmfExport, exportstate = SwmfImport, &
            clock=localclock, syncflag=SyncFlag, rc=rc)
       if(rc /= ESMF_SUCCESS)call my_error('ESMF_CplCompRun')

       ! Run the subcomponents concurrently if possible
       call ESMF_GridCompRun(SwmfComp, importState=SwmfImport, &
            clock=localClock, syncflag=SyncFlag, rc=rc)
       if(rc /= ESMF_SUCCESS)call my_error('ESMF_GridCompRun Swmf')

       call ESMF_GridCompRun(EsmfComp, exportState=EsmfExport, &
            clock=localClock, syncflag=SyncFlag, rc=rc)
       if(rc /= ESMF_SUCCESS)call my_error('ESMF_GridCompRun Esmf')

       ! Advance the time
       call ESMF_ClockAdvance(localclock, rc=rc)
       if(rc /= ESMF_SUCCESS)call my_error('ESMF_ClockAdvance')

    end do

    call ESMF_ClockDestroy(localclock, rc=rc)
    if(rc /= ESMF_SUCCESS)call my_error('ESMF_ClockDestroy')

    rc = ESMF_SUCCESS
    call ESMF_LogWrite("ESMF-SWMF Run finished", ESMF_LOGMSG_INFO)

  end subroutine my_run
  !============================================================================
  subroutine my_final(gcomp, importState, exportState, parentclock, rc)

    type(ESMF_GridComp):: gcomp
    type(ESMF_State):: importState
    type(ESMF_State):: exportState
    type(ESMF_Clock):: parentclock
    integer, intent(out) :: rc

    integer :: iError 
    !-------------------------------------------------------------------------
    call ESMF_LogWrite("ESMF-SWMF Finalize routine called", ESMF_LOGMSG_INFO)

    ! If something fails, try finalizing other things, but return with failure
    ! Assume success 
    iError = ESMF_SUCCESS

    ! Give each of the subcomponents and the coupler a chance to finalize 
    call ESMF_GridCompFinalize(SwmfComp, importState=SwmfImport, &
         clock=Parentclock, rc=rc)
    if(rc /= ESMF_SUCCESS) iError = rc

    call ESMF_GridCompFinalize(EsmfComp, exportState=EsmfExport, &
         clock=ParentClock, rc=rc)
    if(rc /= ESMF_SUCCESS) iError = rc

    call ESMF_CplCompFinalize(CouplerComp, &
         importstate=EsmfExport, exportstate=SwmfImport, &
         clock=ParentClock, rc=rc)
    if(rc /= ESMF_SUCCESS) iError = rc

    ! Now remove the Components to free up their resources
    call ESMF_GridCompDestroy(SwmfComp, rc=rc)
    if(rc /= ESMF_SUCCESS) iError = rc

    call ESMF_GridCompDestroy(EsmfComp, rc=rc)
    if(rc /= ESMF_SUCCESS) iError = rc

    call ESMF_CplCompDestroy(CouplerComp, rc=rc)
    if(rc /= ESMF_SUCCESS) iError = rc

    rc = iError
    call ESMF_LogWrite( "ESMF-SWMF Finalize routine finished", &
         ESMF_LOGMSG_INFO)

  end subroutine my_final
  !============================================================================
  subroutine my_error(String)
    
    character(len=*), intent(in) :: String

    write(*,*)'ERROR in ESMF_SWMF_GridCompMod: ',String
    
    call ESMF_finalize

  end subroutine my_error
  !============================================================================
end module ESMF_SWMF_GridCompMod
