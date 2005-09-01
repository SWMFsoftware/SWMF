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
  use SwmfGridCompMod,     only : SWMF_SetServices    => SetServices
  use EsmfGridCompMod,     only : ESMF_SetServices    => SetServices
  use EsmfSwmfCouplermod, only : Coupler_SetServices => SetServices

  implicit none
  private

  public ESMF_SWMF_SetServices

  type(ESMF_GridComp), save :: comp1Grid, comp2Grid
  type(ESMF_CplComp), save :: compCoupler
  character(len=ESMF_MAXSTR), save :: gname1, gname2, cname
  type(ESMF_State), save :: G1imp, G1exp, G2imp, G2exp
  type(ESMF_State), save :: Cplimp, Cplexp

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

  subroutine my_init(gcomp, importState, exportState, parentclock, rc)
    type(ESMF_GridComp) :: gcomp
    type(ESMF_State) :: importState
    type(ESMF_State) :: exportState
    type(ESMF_Clock) :: parentclock
    integer :: rc

    type(ESMF_VM)   :: parentvm
    type(ESMF_Grid) :: parentgrid

    !-------------------------------------------------------------------------

    call ESMF_LogWrite("ESMF_SWMF Initialize routine called", ESMF_LOG_INFO)

    ! Get the layout and grid associated with this component
    call ESMF_GridCompGet(gcomp, vm=parentvm, grid=parentgrid, rc=rc)

    ! Create the SWMF Gridded component
    gname1 = "SWMF Gridded Component"
    comp1Grid = ESMF_GridCompCreate(parentvm, name=gname1, & 
         grid=parentgrid, rc=rc)

    ! Create thw ESMF Gridded component(s, there could be more than one here)
    gname2 = "ESMF Gridded Component"
    comp2Grid = ESMF_GridCompCreate(parentvm, name=gname2, &
         grid=parentgrid, rc=rc)

    ! Create the Coupler component
    cname = "ESMF-SWMF Coupler Component"
    compCoupler = ESMF_CplCompCreate(parentvm, name=cname, rc=rc)

    call ESMF_LogWrite("Component Creates finished", ESMF_LOG_INFO)

    ! Now call the SetServices routine for each so they can register their
    ! subroutines for Init, Run, and Finalize
    call ESMF_GridCompSetServices(comp1Grid,  SWMF_SetServices, rc)
    call ESMF_GridCompSetServices(comp2Grid,  ESMF_SetServices, rc)
    call ESMF_CplCompSetServices(compCoupler, Coupler_SetServices, rc)


    ! Now create Import and Export State objects in order to pass data
    ! between the Coupler and the Gridded Components
    G1imp = ESMF_StateCreate("SWMF Import", ESMF_STATE_IMPORT)
    G1exp = ESMF_StateCreate("SWMF Export", ESMF_STATE_EXPORT)

    G2imp = ESMF_StateCreate("ESMF Import", ESMF_STATE_IMPORT)
    G2exp = ESMF_StateCreate("ESMF Export", ESMF_STATE_EXPORT)

    Cplimp = ESMF_StateCreate("ESMF-SWMF Coupler Import", ESMF_STATE_IMPORT)
    Cplexp = ESMF_StateCreate("ESMF-SWMF Coupler Export", ESMF_STATE_EXPORT)

    call ESMF_StateAddState(Cplimp, G1imp, rc)
    call ESMF_StateAddState(Cplimp, G2imp, rc)
    call ESMF_StateAddState(Cplexp, G1exp, rc)
    call ESMF_StateAddState(Cplexp, G2exp, rc)

    ! Now give each of the subcomponents a chance to initialize themselves.
    call ESMF_GridCompInitialize(comp1Grid, G1imp, G1exp, parentclock, rc=rc)
    if(rc /= ESMF_SUCCESS) RETURN

    call ESMF_GridCompInitialize(comp2Grid, G2imp, G2exp, parentclock, rc=rc)
    if(rc /= ESMF_SUCCESS) RETURN

    call ESMF_CplCompInitialize(compCoupler, Cplimp, Cplexp, parentclock, &
         rc=rc)
    if(rc /= ESMF_SUCCESS) RETURN

    call ESMF_LogWrite("ESMF-SWMF Grid Component Initialize finished", &
         ESMF_LOG_INFO)

  end subroutine my_init

  !============================================================================

  subroutine my_run(gcomp, importState, exportState, parentclock, rc)
    type(ESMF_GridComp) :: gcomp
    type(ESMF_State) :: importState
    type(ESMF_State) :: exportState
    type(ESMF_Clock) :: parentclock
    integer :: rc
    !-------------------------------------------------------------------------

    call ESMF_LogWrite("ESMF-SWMF Run routine called", ESMF_LOG_INFO)

    ! Now run and couple the subcomponents (assumes an SWMF->ESMF coupling)
    call ESMF_GridCompRun(comp1Grid, G1imp, G1exp, parentclock, rc=rc)
    if(rc /= ESMF_SUCCESS) RETURN

    call ESMF_CplCompRun(compCoupler, G1exp, G2imp, parentclock, rc=rc)
    if(rc /= ESMF_SUCCESS) RETURN

    call ESMF_GridCompRun(comp2Grid, G2imp, G2exp, parentclock, rc=rc)
    if(rc /= ESMF_SUCCESS) RETURN

    call ESMF_LogWrite("ESMF-SWMF Run finished", ESMF_LOG_INFO)

  end subroutine my_run

  !============================================================================

  subroutine my_final(gcomp, importState, exportState, parentclock, iError)
    type(ESMF_GridComp) :: gcomp
    type(ESMF_State) :: importState
    type(ESMF_State) :: exportState
    type(ESMF_Clock) :: parentclock
    integer :: iError
    integer :: rc
    !-------------------------------------------------------------------------

    call ESMF_LogWrite("ESMF-SWMF Finalize routine called", ESMF_LOG_INFO)

    ! Assume success
    iError = ESMF_SUCCESS

    ! Give each of the subcomponents and the coupler a chance to finalize 
    call ESMF_GridCompFinalize(comp1Grid, G1imp, G1exp, parentclock, rc=rc)
    if(rc /= ESMF_SUCCESS) iError = rc

    call ESMF_GridCompFinalize(comp2Grid, G2imp, G2exp, parentclock, rc=rc)
    if(rc /= ESMF_SUCCESS) iError = rc

    call ESMF_CplCompFinalize(compCoupler, G1exp, G2imp, parentclock, rc=rc)
    if(rc /= ESMF_SUCCESS) iError = rc

    ! Now remove the Components to free up their resources
    call ESMF_GridCompDestroy(comp1Grid, rc)
    if(rc /= ESMF_SUCCESS) iError = rc

    call ESMF_GridCompDestroy(comp2Grid, rc)
    if(rc /= ESMF_SUCCESS) iError = rc

    call ESMF_CplCompDestroy(compCoupler, rc)
    if(rc /= ESMF_SUCCESS) iError = rc

    call ESMF_LogWrite( "ESMF-SWMF Finalize routine finished", ESMF_LOG_INFO)

  end subroutine my_final

end module ESMF_SWMF_GridCompMod

!\end{verbatim}

