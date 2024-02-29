module EsmfGridCompMod

  ! ESMF Framework module
  use ESMF

  implicit none
  private

  public:: SetServices

contains
  !============================================================================
  subroutine SetServices(gcomp, rc)

    type(ESMF_GridComp) :: gcomp
    integer, intent(out):: rc

    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=my_init, rc=rc)
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         userRoutine=my_run, rc=rc)
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_FINALIZE, &
         userRoutine=my_final, rc=rc)

  end subroutine SetServices
  !============================================================================
  subroutine my_init(gComp, importState, exportState, externalClock, rc)

    use ESMF_SWMF_Mod, ONLY: add_fields, nVarEsmf, NameFieldEsmf_V, &
         nLonEsmf, nLatEsmf

    type(ESMF_GridComp):: gComp
    type(ESMF_State)   :: importState
    type(ESMF_State)   :: exportState
    type(ESMF_Clock)   :: externalClock
    integer, intent(out):: rc

    ! Access to the MHD data
    type(ESMF_Field):: Field
    real(ESMF_KIND_R8), pointer :: Ptr(:,:)
    integer                     :: iVar, i, j
    character(len=4):: NameField
    ! Units
    real, parameter :: nT=1e-9, amu=1.6726*1e-27, cc=1e-6, kms=1e3, kb=1.38E-23
    !-------------------------------------------------------------------------
    call ESMF_LogWrite("ESMFGridComp init called", ESMF_LOGMSG_INFO)
    rc = ESMF_FAILURE

    ! Add MHD fields to the export state
    call add_fields(gComp, ExportState, IsFromEsmf=.true., rc=rc)
    if(rc /= ESMF_SUCCESS) call my_error("add_fields failed")

    ! Initialize the data
    do iVar = 1, nVarEsmf
       ! Get pointers to the variables in the export state
       nullify(Ptr)
       NameField = NameFieldEsmf_V(iVar)
       call ESMF_StateGet(ExportState, itemName=NameField, field=Field, rc=rc)
       if(rc /= ESMF_SUCCESS) call my_error("ESMF_StateGet for "//NameField)
            
       call ESMF_FieldGet(Field, farrayPtr=Ptr, rc=rc) 
       if(rc /= ESMF_SUCCESS) call my_error("ESMF_FieldGet for "//NameField)

       if(rc /= ESMF_SUCCESS) RETURN
       select case(NameField)
       case('Ped')
          Ptr =  5.0*amu/cc             ! 5 amu/cc
       case('Hall')
          Ptr = -400.0*kms              ! -400km/s
       case default
          write(*,*)'ERROR in ESMF_GridComp:init: unknown NameField=',&
               NameField,' for iVar=',iVar
          rc = ESMF_FAILURE; return
       end select

       ! Add coordinate dependence (5% in Y, 10% in Z)
       if(nLonEsmf > 1 .and. nLatEsmf > 1)then
          do j = 1, nLatEsmf; do i = 1, nLonEsmf
             Ptr(i,j) = Ptr(i,j) &
                  * (0.95 + 0.1*(i - 1.0)/(nLonEsmf - 1)) &
                  * (0.90 + 0.2*(j - 1.0)/(nLatEsmf - 1))
          end do; end do
       end if

    end do

    rc = ESMF_SUCCESS
    call ESMF_LogWrite("ESMFGridComp init returned", ESMF_LOGMSG_INFO)

  end subroutine my_init
  !============================================================================
  subroutine my_run(gComp, importState, exportState, externalclock, rc)

    type(ESMF_GridComp):: gComp
    type(ESMF_State)   :: importState
    type(ESMF_State)   :: exportState
    type(ESMF_Clock)   :: externalclock
    integer, intent(out):: rc

    type(ESMF_Field):: Field
    
    ! Access to the MHD data
    real(ESMF_KIND_R8), pointer :: Ptr(:,:)
    !--------------------------------------------------------------------------
    call ESMF_LogWrite("ESMFGridComp run called", ESMF_LOGMSG_INFO)
    rc = ESMF_FAILURE

    ! Get pointers to the MHD variables in the export state
    nullify(Ptr)
    call ESMF_StateGet(ExportState, itemName='Hall', field=Field, rc=rc)
    if(rc /= ESMF_SUCCESS) call my_error("ESMF_StateGet for Hall")
    call ESMF_FieldGet(Field, farrayPtr=Ptr, rc=rc) 
    if(rc /= ESMF_SUCCESS) call my_error("ESMF_FieldGet for Hall")

    ! Update state by changing Hall conductivity
    write(*,*)'ESMFGridComp:run old Hall=',Ptr(1,1)
    Ptr = Ptr - 2.0/30   ! Change by -2 in 1 minute = 30 couplings
    write(*,*)'ESMFGridComp:run new Hall=',Ptr(1,1)

    rc = ESMF_SUCCESS
    call ESMF_LogWrite("ESMFGridComp run returned", ESMF_LOGMSG_INFO)

  end subroutine my_run
  !============================================================================
  subroutine my_final(gcomp, importState, exportState, externalclock, rc)

    type(ESMF_GridComp) :: gcomp
    type(ESMF_State) :: importState
    type(ESMF_State) :: exportState
    type(ESMF_Clock) :: externalclock
    integer, intent(out):: rc

    call ESMF_LogWrite("ESMFGridComp finalize called", ESMF_LOGMSG_INFO)
    call ESMF_LogWrite("ESMFGridComp finalize returned", ESMF_LOGMSG_INFO)

  end subroutine my_final
  !============================================================================
  subroutine my_error(String)

    character(len=*), intent(in) :: String

    write(*,*)'ERROR in EsmfGridCompMod: ',String
    
    call ESMF_finalize

  end subroutine my_error
  !============================================================================
end module EsmfGridCompMod

!\end{verbatim}
