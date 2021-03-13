!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_wrapper

  implicit none

  save

  private ! except

  public:: SP_set_param
  public:: SP_init_session
  public:: SP_run
  public:: SP_save_restart
  public:: SP_finalize

  ! coupling with MHD components
  public:: SP_adjust_lines
  public:: SP_put_coupling_param
  public:: SP_do_extract_lines

contains
  !============================================================================

  subroutine SP_do_extract_lines(DoExtract)
    logical, intent(out):: DoExtract

    character(len=*), parameter:: NameSub = 'SP_do_extract_lines'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': cannot call the empty version')
  end subroutine SP_do_extract_lines
  !============================================================================
  subroutine SP_run(TimeSimulation,TimeSimulationLimit)

    real,intent(inout)::TimeSimulation
    real,intent(in)::TimeSimulationLimit
    !--------------------------------------------------------------------------
    call CON_stop('Can not call SP_run')
  end subroutine SP_run
  !============================================================================
  subroutine SP_init_session(iSession,TimeSimulation)

    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time
    !--------------------------------------------------------------------------
    call CON_stop('Can not call SP_init_session')
  end subroutine SP_init_session
  !============================================================================
  subroutine SP_finalize(TimeSimulation)

    real,intent(in)::TimeSimulation
    !--------------------------------------------------------------------------
    call CON_stop('Can not call SP_finalize')
  end subroutine SP_finalize
  !============================================================================
  subroutine SP_set_param(CompInfo,TypeAction)
    use CON_comp_info

    type(CompInfoType),intent(inout)       :: CompInfo
    character(len=*), intent(in)           :: TypeAction
    !--------------------------------------------------------------------------
    select case(TypeAction)
    case('VERSION')
       call put(CompInfo,&
            Use        =.false., &
            NameVersion='Empty', &
            Version    =0.0)

    case default
       call CON_stop('Can not call SP_set_param for '//trim(TypeAction))
    end select
  end subroutine SP_set_param
  !============================================================================
  subroutine SP_save_restart(TimeSimulation)

    real,     intent(in) :: TimeSimulation
    !--------------------------------------------------------------------------
    call CON_stop('Can not call SP_save restart')
  end subroutine SP_save_restart
  !============================================================================
  subroutine SP_put_coupling_param(Source_, TimeIn)
    real,           intent(in):: TimeIn
    integer,        intent(in):: Source_
    character(len=*), parameter:: NameSub = 'SP_put_coupling_param'
    !--------------------------------------------------------------------------
    call CON_stop('SP:'//NameSub//': cannot call the empty version')
  end subroutine SP_put_coupling_param
  !============================================================================
  subroutine SP_adjust_lines(Source_)
    integer, intent(in):: Source_
    character(len=*), parameter:: NameSub = 'SP_adjust_lines'
    !--------------------------------------------------------------------------
    call CON_stop('SP:'//NameSub//': cannot call the empty version')
  end subroutine SP_adjust_lines
  !============================================================================
end module SP_wrapper
