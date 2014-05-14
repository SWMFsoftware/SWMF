!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module PC_wrapper

  ! Wrapper for an "empty" Particle in Cell (PC) component

  implicit none

  private ! except

  public:: PC_set_param
  public:: PC_init_session
  public:: PC_run
  public:: PC_save_restart
  public:: PC_finalize

contains
  !==========================================================================
  subroutine PC_set_param(CompInfo, TypeAction)

    use CON_comp_info

    character (len=*), parameter :: NameSub='PC_set_param'

    ! Arguments
    type(CompInfoType), intent(inout):: CompInfo   ! Information for this comp.
    character (len=*), intent(in)    :: TypeAction ! What to do
    !-------------------------------------------------------------------------
    select case(TypeAction)
    case('VERSION')
       call put(CompInfo,&
            Use        =.false., &
            NameVersion='Empty', &
            Version    =0.0)

    case default
       call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')
    end select

  end subroutine PC_set_param

  !============================================================================

  subroutine PC_init_session(iSession, TimeSimulation)

    !INPUT PARAMETERS:
    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='PC_init_session'

    call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')

  end subroutine PC_init_session

  !============================================================================

  subroutine PC_finalize(TimeSimulation)

    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='PC_finalize'

    call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')

  end subroutine PC_finalize

  !============================================================================

  subroutine PC_save_restart(TimeSimulation)

    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='PC_save_restart'

    call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')

  end subroutine PC_save_restart

  !============================================================================

  subroutine PC_run(TimeSimulation,TimeSimulationLimit)

    !INPUT/OUTPUT ARGUMENTS:
    real, intent(inout):: TimeSimulation   ! current time of component

    !INPUT ARGUMENTS:
    real, intent(in):: TimeSimulationLimit ! simulation time not to be exceeded

    character(len=*), parameter :: NameSub='PC_run'

    call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')

  end subroutine PC_run

end module PC_wrapper
