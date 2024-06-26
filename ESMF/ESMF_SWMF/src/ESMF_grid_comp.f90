module ESMF_grid_comp

  ! Code for the ESMFSWMF Gridded Component which creates 4 child Components:
  ! IPE, SWMF, IE and IPE_ie_coupler.

  ! ESMF Framework module
  use ESMF

  use ESMFSWMF_variables, ONLY: &
       iProcRootSwmf, iProcRootEsmf, iProcLastEsmf, iProcLastSwmf, &
       iProc0SwmfComp, iProcLastSwmfComp, nProcSwmfComp, &
       SyncFlag, write_log, write_error

  ! User Component registration routines
  use SWMF_grid_comp, ONLY: swmf_set_services    => set_services
  use IPE_grid_comp,  ONLY: ipe_set_services     => set_services
  use RIM_grid_comp,  ONLY: rim_set_services     => set_services
  use IPERIM_coupler, ONLY: coupler_set_services => set_services

  implicit none
  private

  public:: ESMF_set_services

  type(ESMF_GridComp), save :: IpeComp, SwmfComp, RimComp
  type(ESMF_CplComp),  save :: CouplerComp
  type(ESMF_State),    save :: IpeExport, RimImport

contains
  !============================================================================
  subroutine ESMF_set_services(gComp, iError)

    type(ESMF_GridComp) :: gComp
    integer, intent(out):: iError
    !--------------------------------------------------------------------------
    call ESMF_GridCompSetEntryPoint(gComp, ESMF_METHOD_INITIALIZE, &
         userRoutine=my_init, rc=iError)
    call ESMF_GridCompSetEntryPoint(gComp, ESMF_METHOD_RUN, &
         userRoutine=my_run, rc=iError)
    call ESMF_GridCompSetEntryPoint(gComp, ESMF_METHOD_FINALIZE, &
         userRoutine=my_final, rc=iError)

  end subroutine ESMF_set_services
  !============================================================================
  subroutine my_init(gComp, ImportState, ExportState, ParentClock, iError)

    type(ESMF_GridComp):: gComp
    type(ESMF_State):: ImportState, ExportState
    type(ESMF_Clock):: Parentclock
    integer, intent(out) :: iError

    type(ESMF_VM)      :: ParentVM
    type(ESMF_Grid)    :: SwmfGrid, EsmfGrid
    type(ESMF_DELayout):: Layout

    integer, allocatable:: iProcCouple_I(:)
    integer :: i
    !--------------------------------------------------------------------------
    call write_log("ESMF_gric_comp init called")
    iError = ESMF_FAILURE

    ! Create the SWMF Gridded component
    SwmfComp = ESMF_GridCompCreate(name="SWMF Gridded Component", &
         petlist = [ (i, i=iProcRootSwmf, iProcLastSwmf) ], rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridCompCreate SWMF')

    ! Create the ESMF Gridded component(s, there could be more than one here)
    IpeComp = ESMF_GridCompCreate(name="IPE Gridded Component", &
         petlist = [ (i, i=iProcRootEsmf, iProcLastEsmf) ], rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridCompCreate IPE')

    RimComp = ESMF_GridCompCreate(name="RIM Gridded Component", &
         petlist = [ (i, i=iProc0SwmfComp, iProcLastSwmfComp) ], rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridCompCreate RIM')

    ! Create the Coupler component.
    ! Assume iProcRootEsmf <= iProc0SwmfComp <= iProcLastSwmfComp
    if(iProcRootEsmf == iProc0SwmfComp .and. nProcSwmfComp == 1)then
       ! All share the same PET
       iProcCouple_I = [ iProcRootEsmf ]
    elseif(iProcRootEsmf < iProc0SwmfComp .and. nProcSwmfComp == 2)then
       ! All run on different PETs
       iProcCouple_I = [ iProcRootEsmf, iProc0SwmfComp, iProcLastSwmfComp ]
    else
       ! They use 2 PETs
       iProcCouple_I = [ iProcRootEsmf, iProcLastSwmfComp ]
    end if
    CouplerComp = ESMF_CplCompCreate(name="ESMF-SWMF Coupler Component", &
         petlist = iProcCouple_I, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridCompCreate Coupler')
    deallocate(iProcCouple_I)

    call write_log("Component Creates finished")

    ! Call the SetServices routine for each so they can register their
    ! subroutines for Init, Run, and Finalize
    call ESMF_GridCompSetServices(IpeComp, &
         userRoutine=ipe_set_services, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridCompSetServices Esmf')

    call ESMF_GridCompSetServices(SwmfComp, &
         userRoutine=swmf_set_services, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridCompSetServices Swmf')

    call ESMF_GridCompSetServices(RimComp, &
         userRoutine=rim_set_services, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridCompSetServices Swmf')

    call ESMF_CplCompSetServices(CouplerComp, &
         userRoutine=coupler_set_services, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridCompSetServices Coupler')

    ! Create Import and Export State objects in order to pass data
    ! between the Coupler and the Gridded Components
    RimImport = ESMF_StateCreate(name="RIM Import", &
         stateintent=ESMF_STATEINTENT_IMPORT, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_StateCreate RimImport')
    IpeExport = ESMF_StateCreate(name="IPE Export", &
         stateintent=ESMF_STATEINTENT_EXPORT, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_StateCreate IpeExport')

    ! Each of the subcomponents initialize themselves.
    call ESMF_GridCompInitialize(IpeComp, exportState = IpeExport, &
         clock=parentclock, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridCompInitialize IPE')

    call ESMF_GridCompInitialize(SwmfComp, clock=parentclock, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridCompInitialize SWMF')

    call ESMF_GridCompInitialize(RimComp, importState=RimImport, &
         clock=parentclock, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridCompInitialize RIM')

    call ESMF_CplCompInitialize(CouplerComp, &
         importState = IpeExport, exportState=RimImport, &
         clock=ParentClock, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_CplCompInitialize')

    iError = ESMF_SUCCESS
    call write_log("ESMF_grid_comp init finished")

  end subroutine my_init
  !============================================================================
  subroutine my_run(gComp, ImportState, ExportState, ParentClock, iError)

    type(ESMF_GridComp):: gComp
    type(ESMF_State)   :: ImportState
    type(ESMF_State)   :: ExportState
    type(ESMF_Clock)   :: Parentclock
    integer, intent(out):: iError

    ! Local variables
    type(ESMF_Clock) :: LocalClock
    !--------------------------------------------------------------------------
    call write_log("ESMF_grid_comp run routine called")

    iError = ESMF_FAILURE

    ! make our own local copy of the clock
    LocalClock = ESMF_ClockCreate(ParentClock, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_ClockCreate')

    do
       if(ESMF_ClockIsStopTime(LocalClock, rc=iError)) EXIT
       if(iError /= ESMF_SUCCESS)call my_error('ESMF_ClockIsStopTime')

       ! Couple the subcomponents first so that SWMF has the input from ESMF
       call ESMF_CplCompRun(CouplerComp, &
            importstate=IpeExport, exportstate=RimImport, &
            clock=LocalClock, syncflag=SyncFlag, rc=iError)
       if(iError /= ESMF_SUCCESS)call my_error('ESMF_CplCompRun')

       ! Run RIM to copy/check the received state
       call ESMF_GridCompRun(RimComp, importState=RimImport, &
            clock=LocalClock, syncflag=SyncFlag, rc=iError)
       if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridCompRun Swmf')

       ! Run SWMF
       call ESMF_GridCompRun(SwmfComp, &
            clock=LocalClock, syncflag=SyncFlag, rc=iError)
       if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridCompRun Swmf')

       ! Run IPE
       call ESMF_GridCompRun(IpeComp, exportState=IpeExport, &
            clock=LocalClock, syncflag=SyncFlag, rc=iError)
       if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridCompRun Esmf')

       ! Advance the time
       call ESMF_ClockAdvance(localclock, rc=iError)
       if(iError /= ESMF_SUCCESS)call my_error('ESMF_ClockAdvance')

    end do

    call ESMF_ClockDestroy(localclock, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_ClockDestroy')

    iError = ESMF_SUCCESS
    call write_log("ESMF_grid_comp run finished")

  end subroutine my_run
  !============================================================================
  subroutine my_final(gComp, ImportState, ExportState, Parentclock, iError)

    type(ESMF_GridComp):: gComp
    type(ESMF_State):: ImportState
    type(ESMF_State):: ExportState
    type(ESMF_Clock):: ParentClock
    integer, intent(out) :: iError
    !--------------------------------------------------------------------------
    call write_log("ESMF-SWMF Finalize routine called")

    call ESMF_GridCompFinalize(SwmfComp, clock=Parentclock, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('ESMF_GridCompFinalize SwmfComp')

    ! Give each of the subcomponents and the coupler a chance to finalize
    call ESMF_GridCompFinalize(RimComp, importState=RimImport, &
         clock=Parentclock, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('ESMF_GridCompFinalize RimComp')

    call ESMF_GridCompFinalize(IpeComp, exportState=IpeExport, &
         clock=ParentClock, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('ESMF_GridCompFinalize IpeComp')

    call ESMF_CplCompFinalize(CouplerComp, &
         importState=IpeExport, exportState=RimImport, &
         clock=ParentClock, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('ESMF_CplCompFinalize')

    ! Now remove the Components to free up their resources
    call ESMF_GridCompDestroy(IpeComp, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('ESMF_GridCompDestroy IpeComp')

    call ESMF_GridCompDestroy(SwmfComp, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('ESMF_GridCompDestroy SwmfComp')

    call ESMF_GridCompDestroy(RimComp, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('ESMF_GridCompDestroy RimComp')

    call ESMF_CplCompDestroy(CouplerComp, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('ESMF_CplCompDestroy')

    iError = ESMF_SUCCESS
    call write_log( "ESMF-SWMF Finalize routine finished")

  end subroutine my_final
  !============================================================================
  subroutine my_error(String)

    character(len=*), intent(in) :: String
    !--------------------------------------------------------------------------
    call write_error("ESMF_grid_comp "//String)

  end subroutine my_error
  !============================================================================
end module ESMF_grid_comp
!==============================================================================
