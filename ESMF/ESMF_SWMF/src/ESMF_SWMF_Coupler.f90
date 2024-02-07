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

  type(ESMF_RouteHandle), save :: RouteHandle

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

  subroutine my_init(cComp, ImportState, ExportState, Clock, rc)

    use ESMF_SWMF_Mod, only: NameField_V

    type(ESMF_CplComp), intent(inout) :: cComp
    type(ESMF_State),   intent(inout) :: ImportState
    type(ESMF_State),   intent(inout) :: ExportState
    type(ESMF_Clock),   intent(in)    :: Clock
    integer, intent(out)              :: rc

    type(ESMF_Field):: Field1, FIeld2

    type(ESMF_VM)    :: vm
    !------------------------------------------------------------------------
    call ESMF_LogWrite("Coupler Initialize routine called", ESMF_LOG_INFO)
    rc = ESMF_FAILURE

    ! Get VM from coupler component to use in computing the redist
    call ESMF_CplCompGet(cComp, vm=vm, rc=rc)
    if (rc /= ESMF_SUCCESS) RETURN

    ! Since the states were created on non-overlapping PE-s, need to reconcile
    call ESMF_StateReconcile(ImportState, vm, rc=rc)
    if (rc /= ESMF_SUCCESS) RETURN
    call ESMF_StateReconcile(ExportState, vm, rc=rc)
    if (rc /= ESMF_SUCCESS) RETURN

    call ESMF_StateGetField(ImportState, NameField_V(1), Field1, rc=rc)
    if(rc /= ESMF_SUCCESS) RETURN

    call ESMF_StateGetField(ExportState, NameField_V(1), Field2, rc=rc)
    if(rc /= ESMF_SUCCESS) RETURN

    ! Precompute communication patterns: it is the same for all the fields
    call ESMF_FieldRedistStore(Field1, Field2, vm, routehandle=RouteHandle, &
         rc=rc)

    rc = ESMF_SUCCESS
    call ESMF_LogWrite("Coupler Initialize routine returning", ESMF_LOG_INFO)

  end subroutine my_init

  !==========================================================================

  subroutine my_run(cComp, ImportState, ExportState, Clock, rc)

    use ESMF_SWMF_Mod, only: nVar, NameField_V

    type(ESMF_CplComp), intent(in)    :: cComp
    type(ESMF_State),   intent(in)    :: ImportState
    type(ESMF_State),   intent(inout) :: ExportState
    type(ESMF_Clock),   intent(in)    :: Clock
    integer,            intent(out)   :: rc

    type(ESMF_Field) :: Field1, Field2
    integer          :: iVar
    !------------------------------------------------------------------------
    call ESMF_LogWrite("Coupler Run routine called", ESMF_LOG_INFO)

    rc = ESMF_FAILURE

    ! BundleRedist does not work, use FieldRedist
    do iVar=1,nVar
       call ESMF_StateGetField(ImportState, NameField_V(iVar), Field1, rc=rc)
       if (rc /= ESMF_SUCCESS) &
            call my_error('ESMF_StateGetField(ImportState) failed')

       call ESMF_StateGetField(ExportState, NameField_V(iVar), Field2, rc=rc)
       if (rc /= ESMF_SUCCESS) &
            call my_error('ESMF_StateGetField(ExportState) failed')

       call ESMF_FieldRedist(Field1, Field2, RouteHandle, rc=rc)
       if (rc /= ESMF_SUCCESS) call my_error('ESMF_FieldRedist failed')

    end do

    rc = ESMF_SUCCESS
    call ESMF_LogWrite("Coupler Run routine returning", ESMF_LOG_INFO)

  end subroutine my_run

  !==========================================================================

  subroutine my_final(cComp, ImportState, ExportState, Clock, rc)
    type(ESMF_CplComp), intent(in)  :: cComp
    type(ESMF_State),   intent(in)  :: ImportState, ExportState
    type(ESMF_Clock),   intent(in)  :: Clock
    integer,            intent(out) :: rc
    !------------------------------------------------------------------------
    call ESMF_LogWrite("Coupler Finalize routine called", ESMF_LOG_INFO)
    rc = ESMF_FAILURE
    call ESMF_FieldRedistRelease(RouteHandle, rc=rc)
    if (rc /= ESMF_SUCCESS) RETURN

    rc = ESMF_SUCCESS
    call ESMF_LogWrite("Coupler Finalize routine returning", ESMF_LOG_INFO)

  end subroutine my_final

  !===========================================================================

  subroutine my_error(String)

    ! Since the error flag is not returned from my_run due to the 
    ! non-blocking flag, we have to do something drastic here

    character(len=*), intent(in) :: String

    write(*,*)'ERROR in EsmfSwmfCouplerMod:run: ',String
    
    call ESMF_Finalize

  end subroutine my_error

end module EsmfSwmfCouplerMod

!\end{verbatim}

