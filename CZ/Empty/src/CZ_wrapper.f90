!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module CZ_wrapper

  ! Wrapper for the empty Convection Zone (CZ) component

  implicit none

  private ! except

  public:: CZ_set_param
  public:: CZ_init_session
  public:: CZ_run
  public:: CZ_save_restart
  public:: CZ_finalize

contains

  !==========================================================================
  subroutine CZ_set_param(CompInfo, TypeAction)

    use CON_comp_info

    character (len=*), parameter :: NameSub='CZ_set_param'

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
       call CON_stop(NameSub//': CZ_ERROR: empty version cannot be used!')
    end select

  end subroutine CZ_set_param

  !============================================================================

  subroutine CZ_init_session(iSession, TimeSimulation)

    !INPUT PARAMETERS:
    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='CZ_init_session'

    call CON_stop(NameSub//': CZ_ERROR: empty version cannot be used!')

  end subroutine CZ_init_session

  !============================================================================

  subroutine CZ_finalize(TimeSimulation)

    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='CZ_finalize'

    call CON_stop(NameSub//': CZ_ERROR: empty version cannot be used!')

  end subroutine CZ_finalize

  !============================================================================

  subroutine CZ_save_restart(TimeSimulation)

    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='CZ_save_restart'

    call CON_stop(NameSub//': CZ_ERROR: empty version cannot be used!')

  end subroutine CZ_save_restart

  !============================================================================

  subroutine CZ_run(TimeSimulation,TimeSimulationLimit)

    !INPUT/OUTPUT ARGUMENTS:
    real, intent(inout):: TimeSimulation   ! current time of component

    !INPUT ARGUMENTS:
    real, intent(in):: TimeSimulationLimit ! simulation time not to be exceeded

    character(len=*), parameter :: NameSub='CZ_run'

    call CON_stop(NameSub//': CZ_ERROR: empty version cannot be used!')

  end subroutine CZ_run

end module CZ_wrapper
