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

    use ESMF_SWMF_Mod, ONLY: add_fields, nVar, NameField_V, nLon, nLat

    type(ESMF_GridComp):: gComp
    type(ESMF_State)   :: importState
    type(ESMF_State)   :: exportState
    type(ESMF_Clock)   :: externalClock
    integer, intent(out):: rc

    ! Access to the MHD data
    type(ESMF_Field):: Field
    real(ESMF_KIND_R8), pointer :: Ptr(:,:)
    integer                     :: iVar, i, j
    ! Units
    real, parameter :: nT=1e-9, amu=1.6726*1e-27, cc=1e-6, kms=1e3, kb=1.38E-23
    !-------------------------------------------------------------------------
    call ESMF_LogWrite("ESMFGridComp init called", ESMF_LOGMSG_INFO)
    rc = ESMF_FAILURE

    ! Add MHD fields to the export state
    call add_fields(gComp, ExportState, rc=rc)
    if(rc /= ESMF_SUCCESS) call my_error("add_fields failed")

    ! Initialize the data
    do iVar = 1, nVar
       ! Get pointers to the variables in the export state
       nullify(Ptr)
       call ESMF_StateGet(ExportState, itemName=NameField_V(iVar), &
            field=Field, rc=rc)
       if(rc /= ESMF_SUCCESS) call my_error("ESMF_StateGet failed for " &
            //trim(NameField_V(iVar)))
            
       call ESMF_FieldGet(Field, farrayPtr=Ptr, rc=rc) 
       if(rc /= ESMF_SUCCESS) call my_error("ESMF_FieldGet failed for " &
            //trim(NameField_V(iVar)))

       if(rc /= ESMF_SUCCESS) RETURN
       select case(NameField_V(iVar))
       case('Rho')
          Ptr =  5.0*amu/cc             ! 5 amu/cc
       case('Ux')
          Ptr = -400.0*kms              ! -400km/s
       case('Bz')
          Ptr = 1.0*nT                  ! 1 nT
       case('P')
          Ptr = 5.0/cc*kb*100000.0      ! 5/cc*kBoltzmann*100000 K 
       case('Uy', 'Uz', 'Bx', 'By')
          Ptr =  0.0
       case default
          write(*,*)'ERROR in ESMF_GridComp:init: unknown NameField_V=',&
               NameField_V(iVar),' for iVar=',iVar
          rc = ESMF_FAILURE; return
       end select

       ! Add coordinate dependence (5% in Y, 10% in Z)
       if(nLon > 1 .and. nLat > 1)then
          do j = 1, nLat; do i = 1, nLon
             Ptr(i,j) = Ptr(i,j) &
                  * (0.95 + 0.1*(i - 1.0)/(nLon - 1.0)) &
                  * (0.90 + 0.2*(j - 1.0)/(nLat - 1.0))
          end do; end do
       end if

    end do

    rc = ESMF_SUCCESS
    call ESMF_LogWrite("ESMFGridComp init returned", ESMF_LOGMSG_INFO)

  end subroutine my_init
  !============================================================================
  subroutine my_run(gComp, importState, exportState, externalclock, rc)

    use ESMF_SWMF_Mod, ONLY: NameField_V

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
    call ESMF_StateGet(ExportState, itemName='Bz', field=Field, rc=rc)
    if(rc /= ESMF_SUCCESS) call my_error("ESMF_StateGet failed")
    call ESMF_FieldGet(Field, farrayPtr=Ptr, rc=rc) 
    if(rc /= ESMF_SUCCESS) call my_error("ESMF_FieldGet failed")

    ! Update MHD state by changing Bz
    write(*,*)'ESMFGridComp:run old Bz=',Ptr(1,1)
    Ptr = Ptr - 2.0e-9/30   ! Change Bz by -2nT in 1 minute = 30 couplings
    write(*,*)'ESMFGridComp:run new Bz=',Ptr(1,1)

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
