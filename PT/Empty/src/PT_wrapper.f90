!^CFG COPYRIGHT UM
! Wrapper for the empty PTOM (PT) component
!==========================================================================
subroutine PT_set_param(CompInfo, TypeAction)

  use CON_comp_info

  implicit none

  character (len=*), parameter :: NameSub='PT_set_param'

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
     call CON_stop(NameSub//': PT_ERROR: empty version cannot be used!')
  end select

end subroutine PT_set_param

!==============================================================================

subroutine PT_init_session(iSession, TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  integer,  intent(in) :: iSession         ! session number (starting from 1)
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='PT_init_session'

  call CON_stop(NameSub//': PT_ERROR: empty version cannot be used!')

end subroutine PT_init_session

!==============================================================================

subroutine PT_finalize(TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='PT_finalize'

  call CON_stop(NameSub//': PT_ERROR: empty version cannot be used!')

end subroutine PT_finalize

!==============================================================================

subroutine PT_save_restart(TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='PT_save_restart'

  call CON_stop(NameSub//': PT_ERROR: empty version cannot be used!')

end subroutine PT_save_restart

!==============================================================================

subroutine PT_run(TimeSimulation,TimeSimulationLimit)

  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  real, intent(inout) :: TimeSimulation   ! current time of component

  !INPUT ARGUMENTS:
  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded

  character(len=*), parameter :: NameSub='PT_run'

  call CON_stop(NameSub//': PT_ERROR: empty version cannot be used!')

end subroutine PT_run

!==============================================================================
subroutine PT_get_grid_info(nDimOut, iGridOut, iDecompOut)

  implicit none

  integer, intent(out):: nDimOut    ! grid dimensionality
  integer, intent(out):: iGridOut   ! grid index (increases with AMR)
  integer, intent(out):: iDecompOut ! decomposition index
  
  character(len=*), parameter :: NameSub = 'PT_get_grid_info'

  call CON_stop(NameSub//': PT_ERROR: empty version cannot be used!')

end subroutine PT_get_grid_info
!==============================================================================
subroutine PT_put_from_gm(UseData, &
     NameVar, nVar, nPoint, Pos_DI, Data_VI, iPoint_I)

  implicit none

  logical,          intent(in)   :: UseData ! true when data is transferred
                                            ! false if positions are asked
  character(len=*), intent(inout):: NameVar ! List of variables
  integer,          intent(inout):: nVar    ! Number of variables in Data_VI
  integer,          intent(inout):: nPoint  ! Number of points in Pos_DI

  real, pointer:: Pos_DI(:,:)               ! Position vectors

  real,    intent(in):: Data_VI(nVar,nPoint)! Recv data array
  integer, intent(in):: iPoint_I(nPoint)    ! Order of data

  character(len=*), parameter :: NameSub='PT_put_from_gm'
  
  call CON_stop(NameSub//': PT_ERROR: empty version cannot be used!')

end subroutine PT_put_from_gm

