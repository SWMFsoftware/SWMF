!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
! Wrapper for Upper Atmosphere (UA) component
!==========================================================================
subroutine UA_set_param(CompInfo, TypeAction)

  use CON_comp_info

  implicit none

  character (len=*), parameter :: NameSub='UA_set_param'

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
     call CON_stop(NameSub//': UA_ERROR: empty version cannot be used!')
  end select

end subroutine UA_set_param

!==============================================================================

subroutine UA_init_session(iSession, TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  integer,  intent(in) :: iSession         ! session number (starting from 1)
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='UA_init_session'

  call CON_stop(NameSub//': UA_ERROR: empty version cannot be used!')

end subroutine UA_init_session

!==============================================================================

subroutine UA_finalize(TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='UA_finalize'

  call CON_stop(NameSub//': UA_ERROR: empty version cannot be used!')

end subroutine UA_finalize

!==============================================================================

subroutine UA_save_restart(TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='UA_save_restart'

  call CON_stop(NameSub//': UA_ERROR: empty version cannot be used!')

end subroutine UA_save_restart

!==============================================================================

subroutine UA_run(TimeSimulation,TimeSimulationLimit)

  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  real, intent(inout) :: TimeSimulation   ! current time of component

  !INPUT ARGUMENTS:
  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded

  character(len=*), parameter :: NameSub='UA_run'

  call CON_stop(NameSub//': UA_ERROR: empty version cannot be used!')

end subroutine UA_run

!==============================================================================

subroutine UA_fill_electrodynamics(UAr2_fac, UAr2_ped, UAr2_hal, &
            UAr2_lats, UAr2_mlts)

  character(len=*), parameter :: NameSub='UA_fill_electrodynamics'

  call CON_stop(NameSub//': UA_ERROR: empty version cannot be used!')

end subroutine UA_fill_electrodynamics

!==============================================================================

subroutine UA_calc_electrodynamics(UAi_nMLTs, UAi_nLats)

  character(len=*), parameter :: NameSub='UA_calc_electrodynamics'

  call CON_stop(NameSub//': UA_ERROR: empty version cannot be used!')

end subroutine UA_calc_electrodynamics

