!-------------------------------------------------------------------------

!BOP
!
! !DESCRIPTION:
!  A skeletal coupler component for coupling the SWMF with ESMF component(s).
!
!
!\begin{verbatim}

module EsmfSwmfCouplerMod

  ! ESMF Framework module
  use ESMF_Mod

  implicit none
  private

  public SetServices

contains

  subroutine SetServices(ccomp, rc)
    type(ESMF_CplComp) :: ccomp
    integer :: rc

    call ESMF_CplCompSetEntryPoint(ccomp, ESMF_SETINIT, my_init, &
         ESMF_SINGLEPHASE, rc)
    call ESMF_CplCompSetEntryPoint(ccomp, ESMF_SETRUN, my_run, &
         ESMF_SINGLEPHASE, rc)
    call ESMF_CplCompSetEntryPoint(ccomp, ESMF_SETFINAL, my_final, &
         ESMF_SINGLEPHASE, rc)

  end subroutine SetServices

  !===========================================================================

  subroutine my_init(ccomp, importstate, exportstate, externalclock, rc)
    type(ESMF_CplComp) :: ccomp
    type(ESMF_State) :: importstate, exportstate
    type(ESMF_Clock) :: externalclock
    integer :: rc

    type(ESMF_State) :: state1, state2

    call ESMF_LogWrite("Coupler Initialize routine called", ESMF_LOG_INFO)

    call ESMF_StateGetState(importstate,  "SWMF Import", state1, rc)
    call ESMF_StateGetState(importstate,  "ESMF Import", state2, rc)

    call ESMF_LogWrite("Coupler Initialize routine returning", ESMF_LOG_INFO)

  end subroutine my_init

  !==========================================================================

  subroutine my_run(ccomp, importstate, exportstate, externalclock, rc)
    type(ESMF_CplComp) :: ccomp
    type(ESMF_State) :: importstate, exportstate
    type(ESMF_Clock) :: externalclock
    integer :: rc

    call ESMF_LogWrite("Coupler Run routine called", ESMF_LOG_INFO)

    call ESMF_LogWrite("Coupler Run routine returning", ESMF_LOG_INFO)

  end subroutine my_run

  !==========================================================================

  subroutine my_final(ccomp, importstate, exportstate, externalclock, rc)
    type(ESMF_CplComp) :: ccomp
    type(ESMF_State) :: importstate, exportstate
    type(ESMF_Clock) :: externalclock
    integer :: rc

    call ESMF_LogWrite("Coupler Finalize routine called", ESMF_LOG_INFO)

    call ESMF_LogWrite("Coupler Finalize routine returning", ESMF_LOG_INFO)

  end subroutine my_final

end module EsmfSwmfCouplerMod

!\end{verbatim}

