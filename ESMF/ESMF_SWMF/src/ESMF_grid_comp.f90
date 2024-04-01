module ESMF_grid_comp

  ! Code for the ESMFSWMF Gridded Component which creates 4 child Components:
  ! IPE, SWMF, IE and IPE_ie_coupler.
  
  ! ESMF Framework module
  use ESMF
  use ESMFSWMF_variables, ONLY: write_log, write_error

  ! User Component registration routines
  use ESMFSWMF_variables, ONLY: &
       iProcRootSwmf, iProcRootEsmf, iProcLastEsmf, iProcLastSwmf, SyncFlag

  use SWMF_grid_comp, ONLY: SwmfSetServices => SetServices
  use IPE_grid_comp,  ONLY: EsmfSetServices => SetServices
  use IPE_ie_coupler, ONLY: CouplerSetServices => SetServices

  implicit none
  private

  public:: ESMF_set_services

  type(ESMF_GridComp), save :: SwmfComp, SwmfRootComp, EsmfComp
  type(ESMF_CplComp),  save :: CouplerComp
  type(ESMF_State),    save :: SwmfImport, EsmfExport

contains
  !============================================================================
  subroutine ESMF_set_services(gComp, rc)

    type(ESMF_GridComp) :: gComp
    integer, intent(out):: rc
    !--------------------------------------------------------------------------
    call ESMF_GridCompSetEntryPoint(gComp, ESMF_METHOD_INITIALIZE, &
         userRoutine=my_init, rc=rc)
    call ESMF_GridCompSetEntryPoint(gComp, ESMF_METHOD_RUN, &
         userRoutine=my_run, rc=rc)
    call ESMF_GridCompSetEntryPoint(gComp, ESMF_METHOD_FINALIZE, &
         userRoutine=my_final, rc=rc)

  end subroutine ESMF_set_services
  !============================================================================
  subroutine my_init(gComp, ImportState, ExportState, ParentClock, rc)

    type(ESMF_GridComp):: gComp
    type(ESMF_State):: ImportState, ExportState
    type(ESMF_Clock):: Parentclock
    integer, intent(out) :: rc

    type(ESMF_VM)      :: ParentVM
    type(ESMF_Grid)    :: SwmfGrid, EsmfGrid
    type(ESMF_DELayout):: Layout

    integer :: i, iProc, nProc
    !--------------------------------------------------------------------------
    call write_log("ESMF_gric_comp init called")
    rc = ESMF_FAILURE

    ! Get the layout associated with this component
    call ESMF_GridCompGet(gComp, vm=parentvm, rc=rc)
    if(rc /= ESMF_SUCCESS) call my_error('ESMF_GridCompGet')

    ! Get the number of processors
    call ESMF_VMGet(parentvm, petCount=nProc, localPET=iProc, rc=rc)
    if(rc /= ESMF_SUCCESS) call	my_error('ESMF_VMGet')

    ! Create the SWMF Gridded component
    SwmfComp = ESMF_GridCompCreate(name="SWMF Gridded Component", & 
         petlist = [ (i, i=iProcRootSwmf, iProcLastSwmf) ], rc=rc)
    if(rc /= ESMF_SUCCESS)call my_error('ESMF_GridCompCreate Swmf')

    ! Create the ESMF Gridded component(s, there could be more than one here)
    EsmfComp = ESMF_GridCompCreate(name="ESMF Gridded Component", &
         petlist = [ (i, i=iProcRootEsmf, iProcLastEsmf) ], rc=rc)
    if(rc /= ESMF_SUCCESS)call my_error('ESMF_GridCompCreate Esmf')

    ! Create the Coupler component
    CouplerComp = ESMF_CplCompCreate(name="ESMF-SWMF Coupler Component", &
         petlist = [ (i, i=0, nProc-1) ], rc=rc)
    if(rc /= ESMF_SUCCESS)call my_error('ESMF_GridCompCreate Coupler')

    call write_log("Component Creates finished")

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
    call write_log("ESMF_grid_comp init finished")

  end subroutine my_init
  !============================================================================
  subroutine my_run(gComp, ImportState, ExportState, ParentClock, rc)

    type(ESMF_GridComp):: gComp
    type(ESMF_State)   :: ImportState
    type(ESMF_State)   :: ExportState
    type(ESMF_Clock)   :: Parentclock
    integer, intent(out):: rc

    ! Local variables
    type(ESMF_Clock) :: localclock
    !-------------------------------------------------------------------------
    call write_log("ESMF_grid_comp run routine called")

    rc = ESMF_FAILURE

    ! make our own local copy of the clock
    LocalClock = ESMF_ClockCreate(parentclock, rc=rc)
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
    call write_log("ESMF_grid_comp run finished")

  end subroutine my_run
  !============================================================================
  subroutine my_final(gComp, ImportState, ExportState, Parentclock, rc)

    type(ESMF_GridComp):: gcomp
    type(ESMF_State):: importState
    type(ESMF_State):: exportState
    type(ESMF_Clock):: parentclock
    integer, intent(out) :: rc

    integer :: iError 
    !-------------------------------------------------------------------------
    call write_log("ESMF-SWMF Finalize routine called")

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
    call write_log( "ESMF-SWMF Finalize routine finished")

  end subroutine my_final
  !============================================================================
  subroutine my_error(String)
    
    character(len=*), intent(in) :: String

    call write_error("ESMF_grid_comp "//String)

  end subroutine my_error
  !============================================================================
end module ESMF_grid_comp
