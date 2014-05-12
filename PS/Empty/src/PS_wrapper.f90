!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module PS_wrapper

  ! Wrapper for the "empty" Plasmasphere (PS) component

  implicit none

  private ! except
 
  public:: PS_set_param
  public:: PS_init_session
  public:: PS_run
  public:: PS_save_restart
  public:: PS_finalize

  ! coupling with IE
  public:: PS_put_from_ie

contains
  !==========================================================================
  subroutine PS_set_param(CompInfo, TypeAction)

    use CON_comp_info

    character (len=*), parameter :: NameSub='PS_set_param'

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
       call CON_stop(NameSub//': PS_ERROR: empty version cannot be used!')
    end select

  end subroutine PS_set_param

  !============================================================================

  subroutine PS_init_session(iSession, TimeSimulation)

    !INPUT PARAMETERS:
    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='PS_init_session'

    call CON_stop(NameSub//': PS_ERROR: empty version cannot be used!')

  end subroutine PS_init_session

  !============================================================================

  subroutine PS_finalize(TimeSimulation)

    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='PS_finalize'

    call CON_stop(NameSub//': PS_ERROR: empty version cannot be used!')

  end subroutine PS_finalize

  !============================================================================

  subroutine PS_save_restart(TimeSimulation)

    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='PS_save_restart'

    call CON_stop(NameSub//': PS_ERROR: empty version cannot be used!')

  end subroutine PS_save_restart

  !============================================================================

  subroutine PS_run(TimeSimulation,TimeSimulationLimit)

    !INPUT/OUTPUT ARGUMENTS:
    real, intent(inout) :: TimeSimulation   ! current time of component

    !INPUT ARGUMENTS:
    real, intent(in):: TimeSimulationLimit ! simulation time not to be exceeded

    character(len=*), parameter :: NameSub='PS_run'

    call CON_stop(NameSub//': PS_ERROR: empty version cannot be used!')

  end subroutine PS_run

  !============================================================================

  subroutine PS_put_from_ie(iSize, jSize, Buffer_II)

    integer, intent(in):: iSize, jSize
    real, intent(out)  :: Buffer_II(iSize,jSize)

    character (len=*),parameter :: NameSub='PS_put_from_ie'

    call CON_stop(NameSub//': PS_ERROR: empty version cannot be used!')

  end subroutine PS_put_from_ie

end module PS_wrapper
