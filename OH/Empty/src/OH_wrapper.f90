!^CFG COPYRIGHT UM
! Wrapper for an "empty" Outer Heliosphere (OH) component
!==========================================================================
subroutine OH_set_param(CompInfo, TypeAction)

  use CON_comp_info

  implicit none

  character (len=*), parameter :: NameSub='OH_set_param'

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
     call CON_stop(NameSub//': OH_ERROR: empty version cannot be used!')
  end select

end subroutine OH_set_param

!==============================================================================

subroutine OH_init_session(iSession, TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  integer,  intent(in) :: iSession         ! session number (starting from 1)
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='OH_init_session'

  call CON_stop(NameSub//': OH_ERROR: empty version cannot be used!')

end subroutine OH_init_session

!==============================================================================

subroutine OH_finalize(TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='OH_finalize'

  call CON_stop(NameSub//': OH_ERROR: empty version cannot be used!')

end subroutine OH_finalize

!==============================================================================

subroutine OH_save_restart(TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='OH_save_restart'

  call CON_stop(NameSub//': OH_ERROR: empty version cannot be used!')

end subroutine OH_save_restart

!==============================================================================

subroutine OH_run(TimeSimulation,TimeSimulationLimit)

  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  real, intent(inout) :: TimeSimulation   ! current time of component

  !INPUT ARGUMENTS:
  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded

  character(len=*), parameter :: NameSub='OH_run'

  call CON_stop(NameSub//': OH_ERROR: empty version cannot be used!')

end subroutine OH_run

!===============================================================

subroutine OH_synchronize_refinement(iProc0,iCommUnion)

  implicit none
  integer, intent(in) ::iProc0,iCommUnion
  character(len=*), parameter :: NameSub='OH_synchronize_refinement'

  call CON_stop(NameSub//': OH_ERROR: empty version cannot be used!')

end subroutine OH_synchronize_refinement

!===============================================================
subroutine OH_get_for_ih(&
     nPartial,iGetStart,Get,W,State_V,nVar)
! derived type parameters, it is easier not to declare them
  character(len=*), parameter :: NameSub='OH_get_for_ih'

  call CON_stop(NameSub//': OH_ERROR: empty version cannot be used!')
end subroutine OH_get_for_ih
!===================================================================!
subroutine OH_set_buffer_grid(DD,CompID_)
  character(len=*), parameter :: NameSub='OH_set_buffer_grid'

  call CON_stop(NameSub//': OH_ERROR: empty version cannot be used!')
end subroutine OH_set_buffer_grid
!===================================================================!
subroutine OH_get_a_line_point(&
     nPartial,iGetStart,Get,W,State_V,nVar)
  use CON_coupler, ONLY: IndexPtrType, WeightPtrType
  implicit none
  
  !INPUT ARGUMENTS:
  integer,intent(in)::nPartial,iGetStart,nVar
  type(IndexPtrType),intent(in)::Get
  type(WeightPtrType),intent(in)::W
  real,dimension(nVar),intent(out)::State_V
  
  character(len=*), parameter :: NameSub='OH_get_a_line_point'

  call CON_stop(NameSub//': OH_ERROR: empty version cannot be used!')
end subroutine OH_get_a_line_point
!===================================================================!
subroutine OH_match_ibc
  implicit none
  character(len=*), parameter :: NameSub='OH_match_ibc'

  call CON_stop(NameSub//': OH_ERROR: empty version cannot be used!')
end subroutine OH_match_ibc
