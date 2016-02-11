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

  ! Point coupler interface
  public:: EE_get_grid_info
  public:: EE_find_points

  ! EE-SC coupling
  public:: EE_get_for_SC
  public:: EE_put_from_sc

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

  !============================================================================
  subroutine EE_find_points(nDimIn, nPoint, Xyz_DI, iProc_I)

    integer, intent(in) :: nDimIn                ! dimension of positions
    integer, intent(in) :: nPoint                ! number of positions
    real,    intent(in) :: Xyz_DI(nDimIn,nPoint) ! positions
    integer, intent(out):: iProc_I(nPoint)       ! processor owning position

    character(len=*), parameter:: NameSub = 'EE_find_points'

    call CON_stop(NameSub//': EE_ERROR: empty version cannot be used!')

  end subroutine EE_find_points
  !============================================================================
  subroutine EE_get_grid_info(nDimOut, iGridOut, iDecompOut)

    integer, intent(out):: nDimOut    ! grid dimensionality
    integer, intent(out):: iGridOut   ! grid index (increases with AMR)
    integer, intent(out):: iDecompOut ! decomposition index 

    character(len=*), parameter :: NameSub='EE_get_grid_info'

    call CON_stop(NameSub//': EE_ERROR: empty version cannot be used!')

  end subroutine EE_get_grid_info
  !============================================================================
  subroutine EE_get_for_sc(IsNew, NameVar, nVarIn, nDimIn, nPoint, Xyz_DI, &
       Data_VI)
    logical,          intent(in):: IsNew   ! true for new point array
    character(len=*), intent(in):: NameVar ! List of variables
    integer,          intent(in):: nVarIn  ! Number of variables in Data_VI
    integer,          intent(in):: nDimIn  ! Dimensionality of positions
    integer,          intent(in):: nPoint  ! Number of points in Xyz_DI

    real, intent(in) :: Xyz_DI(nDimIn,nPoint)  ! Position vectors
    real, intent(out):: Data_VI(nVarIn,nPoint) ! Data array

    character(len=*), parameter:: NameSub='EE_get_for_sc'

    call CON_stop(NameSub//': EE_ERROR: empty version cannot be used!')

  end subroutine EE_get_for_sc
  !============================================================================
  subroutine EE_put_from_sc( &
       NameVar, nVar, nPoint, Data_VI, iPoint_I, Pos_DI)

    character(len=*), intent(inout):: NameVar ! List of variables
    integer,          intent(inout):: nVar    ! Number of variables in Data_VI
    integer,          intent(inout):: nPoint  ! Number of points in Pos_DI

    real,    intent(in), optional:: Data_VI(:,:)    ! Recv data array
    integer, intent(in), optional:: iPoint_I(nPoint)! Order of data
    real, intent(out), allocatable, optional:: Pos_DI(:,:) ! Position vectors

    character(len=*), parameter :: NameSub='EE_put_from_sc'

    call CON_stop(NameSub//': EE_ERROR: empty version cannot be used!')

  end subroutine EE_put_from_sc

end module EE_wrapper
