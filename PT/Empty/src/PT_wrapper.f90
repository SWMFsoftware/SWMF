!^CFG COPYRIGHT UM
! Wrapper for the empty PTOM (PT) component
!==========================================================================
subroutine PT_set_param(CompInfo, TypeAction)

  use CON_comp_info

  implicit none

  character (len=*), parameter :: NameSub='PT_set_param'

  ! Arguments
  type(CompInfoType), intent(inout) :: CompInfo   ! Information for this comp.
  character (len=*), intent(in)     :: TypeAction ! What to do
  !-------------------------------------------------------------------------
  select case(TypeAction)
  case('VERSION')
     call put(CompInfo,&
          Use        =.false., &
          NameVersion='Empty', &
          Version    =0.0)

  case default
     call CON_stop(NameSub//': PT_ERROR: empty version cannot be used!')
  end select

end subroutine PT_set_param

!==============================================================================

subroutine PT_init_session(iSession, TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  integer,  intent(in) :: iSession         ! session number (starting from 1)
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='PT_init_session'

  call CON_stop(NameSub//': PT_ERROR: empty version cannot be used!')

end subroutine PT_init_session

!==============================================================================

subroutine PT_finalize(TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='PT_finalize'

  call CON_stop(NameSub//': PT_ERROR: empty version cannot be used!')

end subroutine PT_finalize

!==============================================================================

subroutine PT_save_restart(TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='PT_save_restart'

  call CON_stop(NameSub//': PT_ERROR: empty version cannot be used!')

end subroutine PT_save_restart

!==============================================================================

subroutine PT_run(TimeSimulation,TimeSimulationLimit)

  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  real, intent(inout) :: TimeSimulation   ! current time of component

  !INPUT ARGUMENTS:
  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded

  character(len=*), parameter :: NameSub='PT_run'

  call CON_stop(NameSub//': PT_ERROR: empty version cannot be used!')

end subroutine PT_run

