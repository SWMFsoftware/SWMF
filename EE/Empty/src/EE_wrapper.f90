!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
! Wrapper for the empty Eruptive Event generator (EE) component
!==========================================================================
subroutine EE_set_param(CompInfo, TypeAction)

  use CON_comp_info

  implicit none

  character (len=*), parameter :: NameSub='EE_set_param'

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
     call CON_stop(NameSub//': EE_ERROR: empty version cannot be used!')
  end select

end subroutine EE_set_param

!==============================================================================

subroutine EE_init_session(iSession, TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  integer,  intent(in) :: iSession         ! session number (starting from 1)
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='EE_init_session'

  call CON_stop(NameSub//': EE_ERROR: empty version cannot be used!')

end subroutine EE_init_session

!==============================================================================

subroutine EE_finalize(TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='EE_finalize'

  call CON_stop(NameSub//': EE_ERROR: empty version cannot be used!')

end subroutine EE_finalize

!==============================================================================

subroutine EE_save_restart(TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='EE_save_restart'

  call CON_stop(NameSub//': EE_ERROR: empty version cannot be used!')

end subroutine EE_save_restart

!==============================================================================

subroutine EE_run(TimeSimulation,TimeSimulationLimit)

  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  real, intent(inout) :: TimeSimulation   ! current time of component

  !INPUT ARGUMENTS:
  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded

  character(len=*), parameter :: NameSub='EE_run'

  call CON_stop(NameSub//': EE_ERROR: empty version cannot be used!')

end subroutine EE_run

