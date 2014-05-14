!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module EE_wrapper

  ! Wrapper for the empty Eruptive Event generator (EE) component

  use ModUtilities, ONLY: flush_unit
  use ModIoUnit, ONLY: io_unit_new

  implicit none

  private ! except

  public:: EE_set_param
  public:: EE_init_session
  public:: EE_run
  public:: EE_save_restart
  public:: EE_finalize

contains

  !==========================================================================
  subroutine EE_set_param(CompInfo, TypeAction)

    use CON_comp_info

    character (len=*), parameter :: NameSub='EE_set_param'

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
       call CON_stop(NameSub//': EE_ERROR: empty version cannot be used!')
    end select

  end subroutine EE_set_param

  !============================================================================

  subroutine EE_init_session(iSession, TimeSimulation)

    !INPUT PARAMETERS:
    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='EE_init_session'

    call CON_stop(NameSub//': EE_ERROR: empty version cannot be used!')

  end subroutine EE_init_session

  !============================================================================

  subroutine EE_finalize(TimeSimulation)

    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='EE_finalize'

    call CON_stop(NameSub//': EE_ERROR: empty version cannot be used!')

  end subroutine EE_finalize

  !============================================================================

  subroutine EE_save_restart(TimeSimulation)

    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='EE_save_restart'

    call CON_stop(NameSub//': EE_ERROR: empty version cannot be used!')

  end subroutine EE_save_restart

  !============================================================================

  subroutine EE_run(TimeSimulation,TimeSimulationLimit)

    !INPUT/OUTPUT ARGUMENTS:
    real, intent(inout):: TimeSimulation   ! current time of component

    !INPUT ARGUMENTS:
    real, intent(in):: TimeSimulationLimit ! simulation time not to be exceeded

    character(len=*), parameter :: NameSub='EE_run'

    call CON_stop(NameSub//': EE_ERROR: empty version cannot be used!')

  end subroutine EE_run

end module EE_wrapper
