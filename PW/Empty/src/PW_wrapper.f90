!^CFG COPYRIGHT UM
! Wrapper for the empty PWOM (PW) component
!==========================================================================
subroutine PW_set_param(CompInfo, TypeAction)

  use CON_comp_info

  implicit none

  character (len=*), parameter :: NameSub='PW_set_param'

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
     call CON_stop(NameSub//': PW_ERROR: empty version cannot be used!')
  end select

end subroutine PW_set_param

!==============================================================================

subroutine PW_init_session(iSession, TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  integer,  intent(in) :: iSession         ! session number (starting from 1)
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='PW_init_session'

  call CON_stop(NameSub//': PW_ERROR: empty version cannot be used!')

end subroutine PW_init_session

!==============================================================================

subroutine PW_finalize(TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='PW_finalize'

  call CON_stop(NameSub//': PW_ERROR: empty version cannot be used!')

end subroutine PW_finalize

!==============================================================================

subroutine PW_save_restart(TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='PW_save_restart'

  call CON_stop(NameSub//': PW_ERROR: empty version cannot be used!')

end subroutine PW_save_restart

!==============================================================================

subroutine PW_run(TimeSimulation,TimeSimulationLimit)

  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  real, intent(inout) :: TimeSimulation   ! current time of component

  !INPUT ARGUMENTS:
  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded

  character(len=*), parameter :: NameSub='PW_run'

  call CON_stop(NameSub//': PW_ERROR: empty version cannot be used!')

end subroutine PW_run

!==============================================================================

subroutine PW_put_from_ie(Buffer_IIV, iSize, jSize, nVar, &
                 Name_V, iBlock)
  implicit none

  character(len=*), parameter :: NameSub='PW_put_from_ie'

  !INPUT ARGUMENTS:
  integer, intent(in):: iSize, jSize, nVar, iBlock
  real, intent(in) :: Buffer_IIV(iSize, jSize, nVar)
  character(len=*), intent(in) :: Name_V(nVar)

  call CON_stop(NameSub//': PW_ERROR: empty version cannot be used!')

end subroutine PW_put_from_ie
