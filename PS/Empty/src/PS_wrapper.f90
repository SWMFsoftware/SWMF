!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
! Wrapper for the "empty" Plasmasphere (PS) component
!==========================================================================
subroutine PS_set_param(CompInfo, TypeAction)

  use CON_comp_info

  implicit none

  character (len=*), parameter :: NameSub='PS_set_param'

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
     call CON_stop(NameSub//': PS_ERROR: empty version cannot be used!')
  end select

end subroutine PS_set_param

!==============================================================================

subroutine PS_init_session(iSession, TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  integer,  intent(in) :: iSession         ! session number (starting from 1)
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='PS_init_session'

  call CON_stop(NameSub//': PS_ERROR: empty version cannot be used!')

end subroutine PS_init_session

!==============================================================================

subroutine PS_finalize(TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='PS_finalize'

  call CON_stop(NameSub//': PS_ERROR: empty version cannot be used!')

end subroutine PS_finalize

!==============================================================================

subroutine PS_save_restart(TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='PS_save_restart'

  call CON_stop(NameSub//': PS_ERROR: empty version cannot be used!')

end subroutine PS_save_restart

!==============================================================================

subroutine PS_run(TimeSimulation,TimeSimulationLimit)

  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  real, intent(inout) :: TimeSimulation   ! current time of component

  !INPUT ARGUMENTS:
  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded

  character(len=*), parameter :: NameSub='PS_run'

  call CON_stop(NameSub//': PS_ERROR: empty version cannot be used!')

end subroutine PS_run

!==============================================================================

subroutine PS_put_from_ie(Buffer_IIV, iSize, jSize, nVar)

  implicit none
  character (len=*),parameter :: NameSub='PS_put_from_ie'

  integer, intent(in)           :: iSize, jSize, nVar
  real, intent(out)             :: Buffer_IIV(iSize,jSize,nVar)

  call CON_stop(NameSub//': PS_ERROR: empty version cannot be used!')

end subroutine PS_put_from_ie

!==============================================================================

subroutine PS_get_grid_size(iSize,jSize)

  implicit none
  integer, intent(out) :: iSize, jSize

  character (len=*),parameter ::  NameSub='PS_get_grid_size'
  call CON_stop(NameSub//': PS_ERROR: empty version cannot be used!')

end subroutine PS_get_grid_size

!==============================================================================

subroutine PS_get_grid(MLTs,Lats)

  implicit none
  real:: MLTs(:,:), Lats(:,:)

  character (len=*),parameter ::  NameSub='PS_get_grid'
  call CON_stop(NameSub//': PS_ERROR: empty version cannot be used!')

end subroutine PS_get_grid

!==============================================================================
