module IPERIM_coupler

  ! Couple IPE with RIM

  ! ESMF Framework module
  use ESMF
  use ESMFSWMF_variables, ONLY: &
       nVarEsmf, NameFieldEsmf_V, &
       write_log, write_error, iProc, iProcRootEsmf, iProcRootSwmf

  implicit none
  private

  public:: set_services

  type(ESMF_RouteHandle), save :: RouteHandle

contains
  !============================================================================
  subroutine set_services(cComp, iError)

    type(ESMF_CplComp) :: cComp
    integer, intent(out):: iError
    !--------------------------------------------------------------------------
    call ESMF_CplCompSetEntryPoint(cComp, ESMF_METHOD_INITIALIZE, &
         userRoutine=my_init, rc=iError)
    call ESMF_CplCompSetEntryPoint(cComp, ESMF_METHOD_RUN, &
         userRoutine=my_run, rc=iError)
    call ESMF_CplCompSetEntryPoint(cComp, ESMF_METHOD_FINALIZE, &
         userRoutine=my_final, rc=iError)

  end subroutine set_services
  !============================================================================
  subroutine my_init(cComp, ImportState, ExportState, Clock, iError)

    type(ESMF_CplComp):: cComp
    type(ESMF_VM):: Vm
    type(ESMF_State):: ImportState, ExportState
    type(ESMF_Clock):: Clock
    integer, intent(out):: iError

    type(ESMF_Field):: Field1, FIeld2
    !--------------------------------------------------------------------------
    if(iProc /= iProcRootEsmf .and. iProc < iProcRootSwmf) RETURN
    call write_log("Coupler Initialize routine called")
    iError = ESMF_FAILURE

    call write_log("Coupler ESMF_CplCompGet")
    call ESMF_CplCompGet(cComp, vm=Vm, rc=iError)
    if (iError/=ESMF_SUCCESS) call my_error('ESMF_CplCompGet Vm')

    call write_log("Coupler ESMF_StateReconcile(ImportState)")
    call ESMF_StateReconcile(ImportState, vm=Vm, rc=iError)
    if (iError/=ESMF_SUCCESS) call my_error('ESMF_StateReconcile Import')

    call write_log("Coupler ESMF_StateReconcile(ExportState)")
    call ESMF_StateReconcile(ExportState, vm=Vm, rc=iError)
    if (iError/=ESMF_SUCCESS) call my_error('ESMF_StateReconcile Export')

    call write_log("Coupler ESMF_StateGet(ImportState)")
    call ESMF_StateGet(ImportState, NameFieldEsmf_V(1), Field1, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('ESMF_StateGet Import')

    call write_log("Coupler ESMF_StateGet(ExportState)")
    call ESMF_StateGet(ExportState, NameFieldEsmf_V(1), Field2, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('ESMF_StateGet Export')

    call write_log("Coupler ESMF_FieldRegridStore")
    ! This should be moved into my_run when the grids rotate !!!
    call ESMF_FieldRegridStore(srcField=Field1, dstfield=Field2, &
         routeHandle=RouteHandle, &
         regridmethod=ESMF_REGRIDMETHOD_BILINEAR, rc=iError)

    if(iError /= ESMF_SUCCESS) call my_error('ESMF_FieldRegridStore')

    iError = ESMF_SUCCESS
    call write_log("Coupler Initialize routine returning")

  end subroutine my_init
  !============================================================================
  subroutine my_run(cComp, ImportState, ExportState, Clock, iError)

    type(ESMF_CplComp):: cComp
    type(ESMF_State):: ImportState, ExportState
    type(ESMF_Clock):: Clock
    integer, intent(out):: iError

    type(ESMF_Field):: Field1, Field2
    character(len=4):: NameField
    integer         :: iVar
    !--------------------------------------------------------------------------
    if(iProc /= iProcRootEsmf .and. iProc < iProcRootSwmf) RETURN
    call write_log("Coupler Run routine called")

    iError = ESMF_FAILURE

    do iVar = 1, nVarEsmf
       NameField = NameFieldEsmf_V(iVar)
       call ESMF_StateGet(ImportState, NameField, Field1, rc=iError)
       if (iError /= ESMF_SUCCESS) &
            call my_error('ESMF_StateGet Import for '//NameField)

       call ESMF_StateGet(ExportState, NameField, Field2, rc=iError)
       if (iError /= ESMF_SUCCESS) &
            call my_error('ESMF_StateGet Export for '//NameField)

       call ESMF_FieldRegrid(Field1, Field2, RouteHandle, rc=iError)

       if (iError /= ESMF_SUCCESS) &
            call my_error('ESMF_FieldRegrid for '//NameField)

    end do

    iError = ESMF_SUCCESS
    call write_log("Coupler Run routine returning")

  end subroutine my_run
  !============================================================================
  subroutine my_final(cComp, ImportState, ExportState, Clock, iError)

    type(ESMF_CplComp):: cComp
    type(ESMF_State):: ImportState, ExportState
    type(ESMF_Clock):: Clock
    integer, intent(out):: iError
    !--------------------------------------------------------------------------
    if(iProc /= iProcRootEsmf .and. iProc < iProcRootSwmf) RETURN
    call write_log("Coupler Finalize routine called")
    iError = ESMF_FAILURE
    call ESMF_FieldRegridRelease(RouteHandle, rc=iError)
    if (iError /= ESMF_SUCCESS) call my_error('ESMF_FieldRedistRelease')

    iError = ESMF_SUCCESS
    call write_log("Coupler Finalize routine returning")

  end subroutine my_final
  !============================================================================
  subroutine my_error(String)

    character(len=*), intent(in) :: String
    !--------------------------------------------------------------------------
    call write_error('IPERIM_coupler '//String)

  end subroutine my_error
  !============================================================================
end module IPERIM_coupler
!==============================================================================
