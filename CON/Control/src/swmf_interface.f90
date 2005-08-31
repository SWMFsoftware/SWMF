!^CMP COPYRIGHT UM

!==============================================================================
! The subroutines in this file provide an interface to the SWMF library, 
! when it does not run in stand alone mode.
! These are all external subroutines, so that the external application 
! does not have to compile any SWMF modules 
! (which avoids a lot of compilation problems).
! All subroutines return an error code, which is 0 on success.
!==============================================================================

subroutine initialize_swmf(iComm, iError)
  use CON_main,      ONLY: initialize
  use CON_variables, ONLY: iErrorSwmf
  implicit none
  integer, intent(in) :: iComm  ! The MPI communicator for the SWMF
  integer, intent(out):: iError
  !---------------------------------------------------------------------------
  call initialize(iComm)
  iError = iErrorSwmf
end subroutine initialize_swmf

!==============================================================================

subroutine run_swmf(iError)
  use CON_main,      ONLY: run
  use CON_variables, ONLY: iErrorSwmf
  implicit none
  integer, intent(out):: iError
  !---------------------------------------------------------------------------
  call run
  iError = iErrorSwmf
end subroutine run_swmf

!==============================================================================

subroutine finalize_swmf(iError)
  use CON_main,      ONLY: finalize
  use CON_variables, ONLY: iErrorSwmf
  implicit none
  integer, intent(out):: iError
  !---------------------------------------------------------------------------
  call finalize
  iError = iErrorSwmf
end subroutine finalize_swmf
