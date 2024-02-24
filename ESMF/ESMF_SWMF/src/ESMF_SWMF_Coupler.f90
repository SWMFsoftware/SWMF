module EsmfSwmfCouplerMod

  !  A skeletal coupler component for coupling the SWMF with ESMF component(s).

  ! ESMF Framework module
  use ESMF

  implicit none
  private

  public SetServices

  type(ESMF_RouteHandle), save :: RouteHandle

contains
  !============================================================================
  subroutine SetServices(ccomp, rc)
    type(ESMF_CplComp) :: ccomp
    integer, intent(out):: rc

    call ESMF_CplCompSetEntryPoint(ccomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=my_init, rc=rc)
    call ESMF_CplCompSetEntryPoint(ccomp, ESMF_METHOD_RUN, &
         userRoutine=my_run, rc=rc)
    call ESMF_CplCompSetEntryPoint(ccomp, ESMF_METHOD_FINALIZE, &
         userRoutine=my_final, rc=rc)

  end subroutine SetServices
  !============================================================================
  subroutine my_init(cComp, ImportState, ExportState, Clock, rc)

    use ESMF_SWMF_Mod, ONLY: NameFieldEsmf_V

    type(ESMF_CplComp):: cComp
    type(ESMF_State):: ImportState, ExportState
    type(ESMF_Clock):: Clock
    integer, intent(out):: rc

    type(ESMF_Field):: Field1, FIeld2

    type(ESMF_VM)    :: vm
    !--------------------------------------------------------------------------
    call ESMF_LogWrite("Coupler Initialize routine called", ESMF_LOGMSG_INFO)
    rc = ESMF_FAILURE

    ! Get VM from coupler component to use in computing the redist
    call ESMF_CplCompGet(cComp, vm=vm, rc=rc)
    if (rc /= ESMF_SUCCESS) call my_error('ESMF_CplCompGet')

    ! Since the states were created on non-overlapping PE-s, need to reconcile
    call ESMF_StateReconcile(ImportState, Vm, rc=rc)
    if (rc /= ESMF_SUCCESS) call my_error('ESMF_StateReconcile Import')
    call ESMF_StateReconcile(ExportState, Vm, rc=rc)
    if (rc /= ESMF_SUCCESS) call my_error('ESMF_StateReconcile Export')

    call ESMF_StateGet(ImportState, NameFieldEsmf_V(1), Field1, rc=rc)
    if(rc /= ESMF_SUCCESS) call my_error('ESMF_StateGet Import')

    call ESMF_StateGet(ExportState, NameFieldEsmf_V(1), Field2, rc=rc)
    if(rc /= ESMF_SUCCESS) call my_error('ESMF_StateGet Export')

    ! Precompute communication patterns: it is the same for all the fields
    call ESMF_FieldRedistStore(srcField=Field1, dstField=Field2, &
         routehandle=RouteHandle, rc=rc)
    if(rc /= ESMF_SUCCESS) call my_error('ESMF_FieldRedistStor')

    rc = ESMF_SUCCESS
    call ESMF_LogWrite("Coupler Initialize routine returning", &
         ESMF_LOGMSG_INFO)

  end subroutine my_init
  !============================================================================
  subroutine my_run(cComp, ImportState, ExportState, Clock, rc)

    use ESMF_SWMF_Mod, only: nVarEsmf, NameFieldEsmf_V

    type(ESMF_CplComp):: cComp
    type(ESMF_State):: ImportState, ExportState
    type(ESMF_Clock):: Clock
    integer, intent(out):: rc

    type(ESMF_Field):: Field1, Field2
    character(len=4):: NameField
    integer         :: iVar
    !--------------------------------------------------------------------------
    call ESMF_LogWrite("Coupler Run routine called", ESMF_LOGMSG_INFO)

    rc = ESMF_FAILURE

    ! BundleRedist does not work, use FieldRedist
    do iVar = 1, nVarEsmf
       NameField = NameFieldEsmf_V(iVar)
       call ESMF_StateGet(ImportState, NameField, Field1, rc=rc)
       if (rc /= ESMF_SUCCESS) &
            call my_error('ESMF_StateGet Import for '//NameField)

       call ESMF_StateGet(ExportState, NameField, Field2, rc=rc)
       if (rc /= ESMF_SUCCESS) &
            call my_error('ESMF_StateGet Export for '//NameField)

       call ESMF_FieldRedist(Field1, Field2, RouteHandle, rc=rc)
       if (rc /= ESMF_SUCCESS) &
            call my_error('ESMF_FieldRedist for '//NameField)

    end do

    rc = ESMF_SUCCESS
    call ESMF_LogWrite("Coupler Run routine returning", ESMF_LOGMSG_INFO)

  end subroutine my_run
  !============================================================================
  subroutine my_final(cComp, ImportState, ExportState, Clock, rc)

    type(ESMF_CplComp):: cComp
    type(ESMF_State):: ImportState, ExportState
    type(ESMF_Clock):: Clock
    integer, intent(out):: rc
    !--------------------------------------------------------------------------
    call ESMF_LogWrite("Coupler Finalize routine called", ESMF_LOGMSG_INFO)
    rc = ESMF_FAILURE
    call ESMF_FieldRedistRelease(RouteHandle, rc=rc)
    if (rc /= ESMF_SUCCESS) call my_error('ESMF_FieldRedistRelease')

    rc = ESMF_SUCCESS
    call ESMF_LogWrite("Coupler Finalize routine returning", ESMF_LOGMSG_INFO)

  end subroutine my_final
  !============================================================================
  subroutine my_error(String)

    ! Since the error flag is not returned from my_run due to the 
    ! non-blocking flag, we have to do something drastic here

    character(len=*), intent(in) :: String

    write(*,*)'ERROR in EsmfSwmfCouplerMod: ',String
    
    call ESMF_Finalize

  end subroutine my_error
  !============================================================================
end module EsmfSwmfCouplerMod
