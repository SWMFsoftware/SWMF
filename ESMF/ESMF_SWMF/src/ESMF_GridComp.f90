!BOP
!
! !DESCRIPTION:
!  A template for the ESMF gridded components to be coupled with the SWMF.
!
!\begin{verbatim}

module EsmfGridCompMod

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


  subroutine my_init(gcomp, importState, exportState, externalclock, rc)
    type(ESMF_GridComp) :: gcomp
    type(ESMF_State) :: importState
    type(ESMF_State) :: exportState
    type(ESMF_Clock) :: externalclock
    integer :: rc


    call ESMF_LogWrite("User initialize routine called", ESMF_LOG_INFO)

  end subroutine my_init


  subroutine my_run(gcomp, importState, exportState, externalclock, rc)
    type(ESMF_GridComp) :: gcomp
    type(ESMF_State) :: importState
    type(ESMF_State) :: exportState
    type(ESMF_Clock) :: externalclock
    integer :: rc

    call ESMF_LogWrite("User run routine called", ESMF_LOG_INFO)

  end subroutine my_run


  subroutine my_final(gcomp, importState, exportState, externalclock, rc)
    type(ESMF_GridComp) :: gcomp
    type(ESMF_State) :: importState
    type(ESMF_State) :: exportState
    type(ESMF_Clock) :: externalclock
    integer :: rc

    call ESMF_LogWrite("User finalize routine called", ESMF_LOG_INFO)

  end subroutine my_final

end module EsmfGridCompMod

!\end{verbatim}
