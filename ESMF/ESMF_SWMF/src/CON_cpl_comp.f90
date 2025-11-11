module CON_cpl_comp

  ! This is the CON (Connector) Coupler Component, which acts as an
  ! interface to transfer data between SWMF and IPE.

  ! ESMF Framework module
  use ESMF
  use NUOPC

  use NUOPC_Connector, ONLY: conSS => SetServices
  use NUOPC_Connector, ONLY: label_ComputeRouteHandle
  use NUOPC_Connector, ONLY: label_ExecuteRouteHandle
  use NUOPC_Connector, ONLY: label_ReleaseRouteHandle
  use NUOPC_Connector, ONLY: NUOPC_ConnectorGet
  use NUOPC_Connector, ONLY: NUOPC_ConnectorSet

  use ESMFSWMF_variables, ONLY: write_log, write_error

  implicit none

  private

  public :: set_services

contains
  !============================================================================
  subroutine set_services(cComp, iError)
    type(ESMF_CplComp) :: cComp
    integer, intent(out) :: iError

    ! derive from NUOPC_Connector
    call NUOPC_CompDerive(cComp, conSS, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('NUOPC_CompDerive')
    call NUOPC_CompSpecialize(cComp, specLabel=label_ComputeRouteHandle, &
         specRoutine=compute_rh, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('NUOPC_CompSpecialize')
    call NUOPC_CompSpecialize(cComp, specLabel=label_ExecuteRouteHandle, &
         specRoutine=execute_rh, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('NUOPC_CompSpecialize')
    call NUOPC_CompSpecialize(cComp, specLabel=label_ReleaseRouteHandle, &
         specRoutine=release_rh, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('NUOPC_CompSpecialize')

  end subroutine set_services
  !============================================================================
  subroutine compute_rh(cComp, iError)
    type(ESMF_CplComp) :: cComp
    integer, intent(out) :: iError

    type(ESMF_FieldBundle) :: srcFields, dstFields
    type(ESMF_State) :: state
    type(ESMF_RouteHandle) :: rh
    integer :: srcTermProcessing = 0
    logical :: ignoreDegenerate = .true.
    !--------------------------------------------------------------------------
    call write_log("CON_cpl_comp:compute_rh routine called")
    iError = ESMF_SUCCESS

    ! Get the source and destination field bundles from the Connector Component
    call NUOPC_ConnectorGet(cComp, srcFields=srcFields, &
      dstFields=dstFields, state=state, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('NUOPC_ConnectorGet failed')

    ! Compute Route Handle routine for the Connector Component    
    call ESMF_FieldBundleRegridStore(srcFields, dstFields, routehandle=rh, &
      regridmethod=ESMF_REGRIDMETHOD_BILINEAR, srcTermProcessing=srcTermProcessing, &
      ignoreDegenerate=ignoreDegenerate, unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
      rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('ESMF_FieldBundleRegridStore failed')

    ! Set routehandle in connector
    call NUOPC_ConnectorSet(cComp, rh=rh, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('NUOPC_ConnectorSet failed')

    call write_log("CON_cpl_comp:compute_rh routine returned")

  end subroutine compute_rh
  !============================================================================
  subroutine execute_rh(cComp, iError)
    type(ESMF_CplComp) :: cComp
    integer, intent(out) :: iError

    type(ESMF_FieldBundle) :: srcFields, dstFields
    type(ESMF_State) :: state
    type(ESMF_RouteHandle) :: rh
    integer :: srcTermProcessing = 0
    logical :: ignoreDegenerate = .true.
    !--------------------------------------------------------------------------
    call write_log("CON_cpl_comp:execute_rh routine called")
    iError = ESMF_SUCCESS

    ! Get the source and destination field bundles from the Connector Component
    call NUOPC_ConnectorGet(cComp, srcFields=srcFields, &
      dstFields=dstFields, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('NUOPC_ConnectorGet failed')

    ! Compute Route Handle routine for the Connector Component    
    call ESMF_FieldBundleRegridStore(srcFields, dstFields, routehandle=rh, &
      regridmethod=ESMF_REGRIDMETHOD_BILINEAR, srcTermProcessing=srcTermProcessing, &
      ignoreDegenerate=ignoreDegenerate, unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
      rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('ESMF_FieldBundleRegridStore failed')

    ! Perform interpolation
    call ESMF_FieldBundleRegrid(srcFields, dstFields, &
      routehandle=rh, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('ESMF_FieldBundleRegrid failed')

    call write_log("CON_cpl_comp:execute_rh routine returned")

  end subroutine execute_rh
  !============================================================================
  subroutine release_rh(cComp, iError)
    type(ESMF_CplComp) :: cComp
    integer, intent(out) :: iError

    type(ESMF_RouteHandle) :: rh    
    !--------------------------------------------------------------------------
    call write_log("CON_cpl_comp:release_rh routine called")

    iError = ESMF_SUCCESS

    call NUOPC_ConnectorGet(cComp, rh=rh, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('NUOPC_ConnectorGet failed')

    call ESMF_FieldBundleRegridRelease(rh, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('ESMF_FieldBundleRegridRelease failed')

    call write_log("CON_cpl_comp:release_rh routine returned")

  end subroutine release_rh
  !============================================================================
  !============================================================================
  subroutine my_error(String)
    ! Write out error message and stop
    character(len=*), intent(in) :: String
    !--------------------------------------------------------------------------
    call write_error('CON_cpl_comp '//String)
  end subroutine my_error
  !============================================================================
end module CON_cpl_comp
!==============================================================================