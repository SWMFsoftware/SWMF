!^CFG COPYRIGHT UM
! Wrapper for the "empty" Inner Magnetosphere (IM) component
!==========================================================================
subroutine IM_set_param(CompInfo, TypeAction)

  use CON_comp_info

  implicit none

  character (len=*), parameter :: NameSub='IM_set_param'

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
     call CON_stop(NameSub//': IM_ERROR: empty version cannot be used!')
  end select

end subroutine IM_set_param

!==============================================================================

subroutine IM_init_session(iSession, TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  integer,  intent(in) :: iSession         ! session number (starting from 1)
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='IM_init_session'

  call CON_stop(NameSub//': IM_ERROR: empty version cannot be used!')

end subroutine IM_init_session

!==============================================================================

subroutine IM_finalize(TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='IM_finalize'

  call CON_stop(NameSub//': IM_ERROR: empty version cannot be used!')

end subroutine IM_finalize

!==============================================================================

subroutine IM_save_restart(TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='IM_save_restart'

  call CON_stop(NameSub//': IM_ERROR: empty version cannot be used!')

end subroutine IM_save_restart

!==============================================================================

subroutine IM_run(TimeSimulation,TimeSimulationLimit)

  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  real, intent(inout) :: TimeSimulation   ! current time of component

  !INPUT ARGUMENTS:
  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded

  character(len=*), parameter :: NameSub='IM_run'

  call CON_stop(NameSub//': IM_ERROR: empty version cannot be used!')

end subroutine IM_run

!==============================================================================

subroutine IM_print_variables(NameSource)

  implicit none
  character(len=*), parameter :: NameSub='IM_print_variables'

  character(len=*),intent(in) :: NameSource

  call CON_stop(NameSub//': IM_ERROR: empty version cannot be used!')

end subroutine IM_print_variables

!==============================================================================

subroutine IM_put_from_ie(nPoint,iPointStart,Index,Weight,DoAdd,Buff_V,nVar)

  ! Some of the arguments have complicated derived type definitions
  ! so they are left unspecified for sake of simplicity

  character(len=*), parameter   :: NameSub='IM_put_from_ie'

  call CON_stop(NameSub//': IM_ERROR: empty version cannot be used!')

end subroutine IM_put_from_ie

!==============================================================================

subroutine IM_put_from_ie_complete

  implicit none

  character(len=*), parameter   :: NameSub='IM_put_from_ie_complete'

  call CON_stop(NameSub//': IM_ERROR: empty version cannot be used!')

end subroutine IM_put_from_ie_complete

!==============================================================================

subroutine IM_put_from_gm(Buffer_IIV,iSizeIn,jSizeIn,nVarIn,NameVar)

  implicit none
  character (len=*),parameter :: NameSub='IM_put_from_gm'

  integer, intent(in) :: iSizeIn,jSizeIn,nVarIn
  real, dimension(iSizeIn,jSizeIn,nVarIn), intent(in) :: Buffer_IIV
  character (len=*),intent(in)       :: NameVar

  call CON_stop(NameSub//': IM_ERROR: empty version cannot be used!')

end subroutine IM_put_from_gm

!==============================================================================

subroutine IM_get_for_gm(Buffer_II,iSizeIn,jSizeIn,NameVar)

  implicit none
  character (len=*),parameter :: NameSub='IM_get_for_gm'

  integer, intent(in) :: iSizeIn,jSizeIn
  real, dimension(iSizeIn,jSizeIn), intent(out) :: Buffer_II
  character (len=*),intent(in)       :: NameVar

  call CON_stop(NameSub//': IM_ERROR: empty version cannot be used!')

end subroutine IM_get_for_gm

!==============================================================================

subroutine IM_put_sat_from_gm(nSats, Buffer_I, Buffer_III)

  implicit none
  character (len=*), parameter :: NameSub='IM_put_sat_from_gm'

  integer, intent(in)            :: nSats
  real, intent(in)               :: Buffer_III(3,2,nSats)
  character(len=100), intent(in) :: Buffer_I(nSats)

  call CON_stop(NameSub//': IM_ERROR: empty version cannot be used!')

end subroutine IM_put_sat_from_gm

!==============================================================================

subroutine IM_get_for_ie(nPoint,iPointStart,Index,Weight,Buff_V,nVar)

  use CON_router,   ONLY: IndexPtrType, WeightPtrType

  implicit none
  character(len=*), parameter :: NameSub='IM_get_for_ie'

  integer,intent(in)            :: nPoint, iPointStart, nVar
  real,intent(out)              :: Buff_V(nVar)
  type(IndexPtrType),intent(in) :: Index
  type(WeightPtrType),intent(in):: Weight

  call CON_stop(NameSub//': IM_ERROR: empty version cannot be used!')

end subroutine IM_get_for_ie

