! Wrapper for Eruptive Event (EE) component
!==============================================================================
subroutine EE_set_param(CompInfo, TypeAction)

  use CON_comp_info

  implicit none

  type(CompInfoType), intent(inout) :: CompInfo   ! Information for this comp.
  character (len=*), intent(in)     :: TypeAction ! What to do
  !----------------------------------------------------------------------------

end subroutine EE_set_param

!==============================================================================
subroutine EE_init_session(iSession, TimeSimulation)

  implicit none

  integer, intent(in) :: iSession         ! session number (starting from 1)
  real,    intent(in) :: TimeSimulation   ! seconds from start time
  !----------------------------------------------------------------------------

end subroutine EE_init_session

!==============================================================================
subroutine EE_finalize(TimeSimulation)

  use EE_ModMain, ONLY: time_loop

  implicit none

  real, intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='EE_finalize'
  !----------------------------------------------------------------------------
  ! We are not advancing in time any longer
  time_loop = .false.

  call EE_BATS_save_files('FINAL')

  call EE_BATSRUS_finalize

end subroutine EE_finalize

!==============================================================================
subroutine EE_save_restart(TimeSimulation)

  implicit none

  real, intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='EE_save_restart'
  !----------------------------------------------------------------------------
  call EE_BATS_save_files('RESTART')

end subroutine EE_save_restart

!==============================================================================
subroutine EE_run(TimeSimulation, TimeSimulationLimit)

  implicit none

  real, intent(inout) :: TimeSimulation   ! current time of component
  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded

  character(len=*), parameter :: NameSub='EE_run'

  logical :: DoTest, DoTestMe
  !----------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)

end subroutine EE_run
