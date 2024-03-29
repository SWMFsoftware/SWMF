module EsmfSwmfCouplerMod

  !  A skeletal coupler component for coupling the SWMF with ESMF component(s).

  ! ESMF Framework module
  use ESMF
  use ESMF_SWMF_Mod, ONLY: log_write

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
    type(ESMF_VM):: Vm
    type(ESMF_State):: ImportState, ExportState
    type(ESMF_Clock):: Clock
    integer, intent(out):: rc

    type(ESMF_Field):: Field1, FIeld2
    !--------------------------------------------------------------------------
    call log_write("Coupler Initialize routine called")
    rc = ESMF_FAILURE

    call log_write("Coupler ESMF_CplCompGet")
    call ESMF_CplCompGet(cComp, vm=Vm, rc=rc)
    if (rc/=ESMF_SUCCESS) call my_error('ESMF_CplCompGet Vm')

    call log_write("Coupler ESMF_StateReconcile(ImportState)")
    call ESMF_StateReconcile(ImportState, vm=Vm, rc=rc)
    if (rc/=ESMF_SUCCESS) call my_error('ESMF_StateReconcile Import')

    call log_write("Coupler ESMF_StateReconcile(ExportState)")
    call ESMF_StateReconcile(ExportState, vm=Vm, rc=rc)
    if (rc/=ESMF_SUCCESS) call my_error('ESMF_StateReconcile Export')

    call log_write("Coupler ESMF_StateGet(ImportState)")
    call ESMF_StateGet(ImportState, NameFieldEsmf_V(1), Field1, rc=rc)
    if(rc /= ESMF_SUCCESS) call my_error('ESMF_StateGet Import')

    call log_write("Coupler ESMF_StateGet(ExportState)")
    call ESMF_StateGet(ExportState, NameFieldEsmf_V(1), Field2, rc=rc)
    if(rc /= ESMF_SUCCESS) call my_error('ESMF_StateGet Export')

    call log_write("Coupler ESMF_FieldRegridStore")
    ! This should be moved into my_run when the grids rotate !!!
    call ESMF_FieldRegridStore(srcField=Field1, dstfield=Field2, &
         routeHandle=RouteHandle, &
         regridmethod=ESMF_REGRIDMETHOD_BILINEAR, rc=rc)

    if(rc /= ESMF_SUCCESS) call my_error('ESMF_FieldRegridStore')

    rc = ESMF_SUCCESS
    call log_write("Coupler Initialize routine returning")

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
    call log_write("Coupler Run routine called")

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

       !call ESMF_FieldRedist(Field1, Field2, RouteHandle, rc=rc)
       call ESMF_FieldRegrid(Field1, Field2, RouteHandle, rc=rc)

       if (rc /= ESMF_SUCCESS) &
            call my_error('ESMF_FieldRegrid for '//NameField)

    end do

    rc = ESMF_SUCCESS
    call log_write("Coupler Run routine returning")

  end subroutine my_run
  !============================================================================
  subroutine my_final(cComp, ImportState, ExportState, Clock, rc)

    type(ESMF_CplComp):: cComp
    type(ESMF_State):: ImportState, ExportState
    type(ESMF_Clock):: Clock
    integer, intent(out):: rc
    !--------------------------------------------------------------------------
    call log_write("Coupler Finalize routine called")
    rc = ESMF_FAILURE
    call ESMF_FieldRegridRelease(RouteHandle, rc=rc)
    if (rc /= ESMF_SUCCESS) call my_error('ESMF_FieldRedistRelease')

    rc = ESMF_SUCCESS
    call log_write("Coupler Finalize routine returning")

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
