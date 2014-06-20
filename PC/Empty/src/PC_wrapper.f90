!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module PC_wrapper

  ! Wrapper for an "empty" Particle in Cell (PC) component

  implicit none

  private ! except

  public:: PC_set_param
  public:: PC_init_session
  public:: PC_finilize_init_session
  public:: PC_run
  public:: PC_save_restart
  public:: PC_finalize
  public:: PC_find_points
  public:: PC_get_for_gm
  public:: PC_get_grid_info
  public:: PC_put_from_gm
  public:: PC_put_from_gm_dt
  public:: PC_put_from_gm_init

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

  subroutine PC_finilize_init_session

    !INPUT PARAMETERS:
    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='PC_finilize_init_session'
    !--------------------------------------------------------------------------

    call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')

  end subroutine PC_finilize_init_session

  !============================================================================

  subroutine PC_init_session(iSession, TimeSimulation)

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

  !==============================================================================

  subroutine PC_find_points(nDimIn, nPoint, Xyz_DI, iProc_I)

    implicit none

    integer, intent(in) :: nDimIn                ! dimension of position vectors
    integer, intent(in) :: nPoint                ! number of positions
    real,    intent(in) :: Xyz_DI(nDimIn,nPoint) ! positions
    integer, intent(out):: iProc_I(nPoint)       ! processor owning position

    character(len=*), parameter:: NameSub = 'PC_find_points'
    !--------------------------------------------------------------------------

    call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')

  end subroutine PC_find_points

  !==============================================================================

  subroutine PC_get_for_gm(IsNew, NameVar, nVarIn, nDimIn, nPoint, Xyz_DI, &
       Data_VI)

    use CON_coupler, ONLY: i_proc, PC_, n_proc

    implicit none

    logical,          intent(in):: IsNew   ! true for new point array
    character(len=*), intent(in):: NameVar ! List of variables
    integer,          intent(in):: nVarIn  ! Number of variables in Data_VI
    integer,          intent(in):: nDimIn  ! Dimensionality of positions
    integer,          intent(in):: nPoint  ! Number of points in Xyz_DI

    real, intent(in) :: Xyz_DI(nDimIn,nPoint)  ! Position vectors
    real, intent(out):: Data_VI(nVarIn,nPoint) ! Data array

    character(len=*), parameter :: NameSub='GM_get_for_pc'
    !--------------------------------------------------------------------------

    call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')

  end subroutine PC_get_for_gm

  !==============================================================================

  subroutine PC_get_grid_info(nDimOut, iGridOut, iDecompOut)

    implicit none

    integer, intent(out):: nDimOut    ! grid dimensionality
    integer, intent(out):: iGridOut   ! grid index (increases with AMR)
    integer, intent(out):: iDecompOut ! decomposition index

    character(len=*), parameter :: NameSub = 'PC_get_grid_info'
    !---------------------------------------------------------------------------

    call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')

  end subroutine PC_get_grid_info

  !==============================================================================

  subroutine PC_put_from_gm( &
       NameVar, nVar, nPoint, Data_VI, iPoint_I, Pos_DI)

    use CON_coupler, ONLY: i_proc, PC_, n_proc

    implicit none

    character(len=*), intent(inout):: NameVar ! List of variables
    integer,          intent(inout):: nVar    ! Number of variables in Data_VI
    integer,          intent(inout):: nPoint  ! Number of points in Pos_DI
    real,    intent(in), optional  :: Data_VI(:,:)    ! Recv data array
    integer, intent(in), optional  :: iPoint_I(nPoint)! Order of data
    real, intent(out), allocatable, optional:: Pos_DI(:,:)               ! Position vectors

    character(len=*), parameter :: NameSub='PC_put_from_gm'
    !--------------------------------------------------------------------------

    call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')

  end subroutine PC_put_from_gm

  !==============================================================================

  subroutine PC_put_from_gm_dt(DtSi)

    implicit none

    real,    intent(in) :: DtSi

    character(len=*), parameter :: NameSub = 'PC_put_from_gm_dt'
    !---------------------------------------------------------------------------

    call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')

  end subroutine PC_put_from_gm_dt
  !==============================================================================
  subroutine PC_put_from_gm_init(ParamInt_I, ParamReal_I, n)

    implicit none

    integer, intent(in) :: n
    integer, intent(in) :: ParamInt_I(4)
    real,    intent(in) :: ParamReal_I(n)

    character(len=*), parameter :: NameSub = 'PC_put_from_gm_init'
    !---------------------------------------------------------------------------

    call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')

  end subroutine PC_put_from_gm_init

end module PC_wrapper
