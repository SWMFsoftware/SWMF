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

  !===========================================================================

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

  !===========================================================================

  subroutine my_init(gComp, importState, exportState, externalclock, rc)

    use ESMF_SWMF_Mod, ONLY: add_mhd_fields, nVar, NameField_V, iMax, jMax

    type(ESMF_GridComp), intent(inout):: gComp
    type(ESMF_State),    intent(in)   :: importState
    type(ESMF_State),    intent(inout):: exportState
    type(ESMF_Clock),    intent(in)   :: externalclock
    integer,             intent(out)  :: rc

    ! Access to the MHD data
    real(ESMF_KIND_R8), pointer :: Ptr(:,:)
    integer                     :: iVar, i, j
    ! Units
    real, parameter :: nT=1e-9, amu=1.6726*1e-27, cc=1e-6, kms=1e3, kb=1.38E-23
    !-------------------------------------------------------------------------
    call ESMF_LogWrite("ESMFGridComp init called",  ESMF_LOG_INFO)
    rc = ESMF_FAILURE

    ! Add MHD fields to the export state
    call add_mhd_fields(gComp, ExportState, 7.77, rc=rc)
    if(rc /= ESMF_SUCCESS) RETURN

    ! Initialize the MHD data
    do iVar = 1, nVar
       ! Get pointers to the MHD variables in the export state
       nullify(Ptr)
       call ESMF_StateGetDataPointer(ExportState, NameField_V(iVar), Ptr, &
            rc=rc)
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
       if(iMax > 1 .and. jMax > 1)then
          do j=1, jMax; do i = 1, iMax
             Ptr(i, j) = Ptr(i,j) &
                  * (0.95 + 0.1*(i-1.0)/(iMax-1.0)) &
                  * (0.90 + 0.2*(j-1.0)/(jMax-1.0))
          end do; end do
       end if

    end do

    rc = ESMF_SUCCESS
    call ESMF_LogWrite("ESMFGridComp init returned", ESMF_LOG_INFO)

  end subroutine my_init

  !===========================================================================

  subroutine my_run(gComp, importState, exportState, externalclock, rc)

    use ESMF_SWMF_Mod, ONLY: NameField_V

    type(ESMF_GridComp), intent(in)   :: gComp
    type(ESMF_State),    intent(in)   :: importState
    type(ESMF_State),    intent(inout):: exportState
    type(ESMF_Clock),    intent(in)   :: externalclock
    integer,             intent(out)  :: rc

    ! Access to the MHD data
    real(ESMF_KIND_R8), pointer :: Ptr(:,:)
    !-------------------------------------------------------------------------
    call ESMF_LogWrite("ESMFGridComp run called", ESMF_LOG_INFO)
    rc = ESMF_FAILURE

    ! Get pointers to the MHD variables in the export state
    nullify(Ptr)
    call ESMF_StateGetDataPointer(ExportState, 'Bz', Ptr, rc=rc)
    if(rc /= ESMF_SUCCESS) call my_error('ESMF_StateGetDataPointer failed')

    ! Update MHD state by changing Bz
    write(*,*)'ESMFGridComp:run old Bz=',Ptr(1,1)
    Ptr = Ptr - 2.0e-9/30   ! Change Bz by -2nT in 1 minute = 30 couplings
    write(*,*)'ESMFGridComp:run new Bz=',Ptr(1,1)

    rc = ESMF_SUCCESS
    call ESMF_LogWrite("ESMFGridComp run returned", ESMF_LOG_INFO)

  end subroutine my_run

  !===========================================================================

  subroutine my_final(gcomp, importState, exportState, externalclock, rc)
    type(ESMF_GridComp) :: gcomp
    type(ESMF_State) :: importState
    type(ESMF_State) :: exportState
    type(ESMF_Clock) :: externalclock
    integer :: rc

    call ESMF_LogWrite("ESMFGridComp finalize called", ESMF_LOG_INFO)
    call ESMF_LogWrite("ESMFGridComp finalize returned", ESMF_LOG_INFO)

  end subroutine my_final

  !===========================================================================

  subroutine my_error(String)

    ! Since the error flag is not returned from my_run due to the 
    ! non-blocking flag, we have to do something drastic here

    character(len=*), intent(in) :: String

    write(*,*)'ERROR in EsmfGridCompMod:run: ',String
    
    call ESMF_finalize

  end subroutine my_error

end module EsmfGridCompMod

!\end{verbatim}
