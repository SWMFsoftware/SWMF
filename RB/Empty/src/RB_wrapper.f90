!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module RB_wrapper

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !               Space Weather Modeling Framework (SWMF)                !
  !    Center for Space Environment Modeling, The University of Michigan !
  !-----------------------------------------------------------------------

  ! !TITLE: Wrapper for the "empty" Radiation Belts Component 
  ! !AUTHORS: Ovsei Volberg
  ! !AFFILIATION: CSEM, The University of Michigan
  ! !DATE: April 13, 2004 - the inital version was written
  ! !INTRODUCTION: This wrapper provides the "empty" interface
  !                for the  Radiation Belts Component             

  implicit none

  private ! except
  
  public:: RB_set_param
  public:: RB_init_session
  public:: RB_run
  public:: RB_save_restart
  public:: RB_finalize

  ! coupling with GM
  public:: RB_put_from_gm
  public:: RB_put_sat_from_gm

  ! coupling with IE
  public:: RB_put_from_ie

contains

  subroutine RB_set_param(CompInfo, TypeAction)

    use CON_comp_info

    character (len=*), intent(in)     :: TypeAction ! which action to perform
    type(CompInfoType), intent(inout) :: CompInfo   ! component information

    character (len=*), parameter :: NameSub='RB_set_param'

    !-------------------------------------------------------------------------
    select case(TypeAction)
    case('VERSION')
       call put(CompInfo,                         &
            Use=.false.,                           &
            NameVersion='Empty', &
            Version=0.0)
    case default
       call CON_stop(NameSub//' RB_ERROR: empty version cannot be used!')
    end select

  end subroutine RB_set_param

  !============================================================================

  subroutine RB_init_session(iSession, TimeSimulation)

    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='RB_init_session'

    call CON_stop(NameSub//': RB_ERROR: empty version cannot be used!')

  end subroutine RB_init_session

  !===========================================================================

  subroutine RB_run(TimeSimulation,TimeSimulationLimit)

    real, intent(in):: TimeSimulationLimit ! simulation time not to be exceeded

    real, intent(inout) :: TimeSimulation   ! current time of component
    character(len=*), parameter :: NameSub='RB_run'

    call CON_stop(NameSub//': RB_ERROR: empty version cannot be used!')
  end subroutine RB_run

  !===========================================================================

  subroutine RB_finalize(TimeSimulation)

    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='RB_finalize'

    call CON_stop(NameSub//': RB_ERROR: empty version cannot be used!')
  end subroutine RB_finalize

  !===========================================================================

  subroutine RB_save_restart(TimeSimulation)

    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='RB_save_restart'

    call CON_stop(NameSub//': RB_ERROR: empty version cannot be used!')
  end subroutine RB_save_restart

  !============================================================================

  subroutine RB_put_from_gm(Buffer_IIV, iSizeIn, jSizeIn, nVarIn,&
       BufferLine_VI, nVarLine, nPointLine, NameVar, tSimulation)

    integer, intent(in) :: iSizeIn, jSizeIn, nVarIn
    real,    intent(in) :: Buffer_IIV(iSizeIn,jSizeIn,nVarIn)
    integer, intent(in) :: nVarLine, nPointLine
    real,    intent(in) :: BufferLine_VI(nVarLine, nPointLine)

    character (len=*),intent(in) :: NameVar
    real, intent(in) :: tSimulation

    character(len=*), parameter :: NameSub='RB_put_from_gm'

    call CON_stop(NameSub//': RB_ERROR: empty version cannot be used!')
  end subroutine RB_put_from_gm
  !============================================================================
  subroutine RB_put_from_ie(Buffer_IIV, iSize, jSize, nVarIn, &
       Name_V, iBlock)

    character(len=*), parameter :: NameSub='RB_put_from_ie'

    !INPUT ARGUMENTS:
    integer, intent(in):: iSize, jSize, nVarIn, iBlock
    real, intent(in) :: Buffer_IIV(iSize, jSize, nVarIn)
    character(len=*), intent(in) :: Name_V(nVarIn)

    call CON_stop(NameSub//': RB_ERROR: empty version cannot be used!')
  end subroutine RB_put_from_ie
  !============================================================================
  subroutine RB_put_sat_from_gm(nSats, Buffer_I, Buffer_III)

    integer, intent(in)            :: nSats
    real, intent(in)               :: Buffer_III(4,2,nSats)
    character(len=100), intent(in) :: Buffer_I(nSats)
  end subroutine RB_put_sat_from_gm

end module RB_wrapper
