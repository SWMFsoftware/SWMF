!^CFG COPYRIGHT UM
! Wrapper for the "empty" Ionosphere Electrodynamics (IE) component
!==========================================================================
subroutine IE_set_param(CompInfo, TypeAction)

  use CON_comp_info

  implicit none

  character (len=*), parameter :: NameSub='IE_set_param'

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
     call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')
  end select

end subroutine IE_set_param
!=============================================================================
subroutine IE_check_allocation_south(nPoint)
  implicit none
  integer,intent(in)::nPoint
  call CON_stop(&
       'IE_check_allocation_south: IE_ERROR: empty version cannot be used!')
end subroutine IE_check_allocation_south
!=============================================================================
subroutine IE_check_allocation_north(nPoint)
  implicit none
  integer,intent(in)::nPoint
  call CON_stop(&
       'IE_check_allocation_north: IE_ERROR: empty version cannot be used!')
end subroutine IE_check_allocation_north

!==============================================================================

subroutine IE_init_session(iSession, TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  integer,  intent(in) :: iSession         ! session number (starting from 1)
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='IE_init_session'

  call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

end subroutine IE_init_session

!==============================================================================

subroutine IE_finalize(TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='IE_finalize'

  call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

end subroutine IE_finalize

!==============================================================================

subroutine IE_save_restart(TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='IE_save_restart'

  call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

end subroutine IE_save_restart

!==============================================================================

subroutine IE_run(TimeSimulation,TimeSimulationLimit)

  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  real, intent(inout) :: TimeSimulation   ! current time of component

  !INPUT ARGUMENTS:
  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded

  character(len=*), parameter :: NameSub='IE_run'

  call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

end subroutine IE_run

!==============================================================================

subroutine IE_get_for_gm_swmf(nPartial,iGetStart,Get,W,State_V,nVar)
  use CON_router                          
  implicit none
  integer,intent(in)::nPartial,iGetStart,nVar
  type(IndexPtrType),intent(in)::Get
  type(WeightPtrType),intent(in)::W
  real,dimension(nVar),intent(out)::State_V
  character(len=*), parameter :: NameSub='IE_get_for_gm_swmf'
  
  call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

end subroutine IE_get_for_gm_swmf

!==============================================================================

subroutine IE_get_for_gm(Buffer_II,iSize,jSize,NameVar)

  implicit none
  character (len=*),parameter :: NameSub='IE_get_for_gm'

  integer, intent(in)           :: iSize,jSize
  real, intent(out)             :: Buffer_II(iSize,jSize)
  character (len=*),intent(in)  :: NameVar

  call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

end subroutine IE_get_for_gm
!==============================================================================

subroutine IE_put_from_gm_swmf(State_V,nVar,iBlock,nPoint,ColatLim)
  implicit none
  !INPUT ARGUMENTS
  integer,intent(in)::nVar,iBlock,nPoint
  real,dimension(nVar),intent(in):: State_V
  real,intent(inout)::ColatLim
  character (len=*),parameter :: NameSub='IE_put_from_gm_swmf'

  call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

end subroutine IE_put_from_gm_swmf

!==============================================================================

subroutine IE_put_from_gm(Buffer_IV,nPoint,nVar,NameVar)

  implicit none
  character (len=*),parameter :: NameSub='IE_put_from_gm'
  integer,          intent(in) :: nPoint, nVar
  real,             intent(in) :: Buffer_IV(nPoint,nVar)
  character(len=*) ,intent(in) :: NameVar


  call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

end subroutine IE_put_from_gm

!==============================================================================

subroutine IE_get_for_im(nPoint,iPointStart,Index,Weight,Buff_V,nVar)

  ! Some of the arguments are complicated derived type
  ! For sake of simplicity the arguments are not declared

  character (len=*),parameter :: NameSub='IE_get_for_im'

  call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

end subroutine IE_get_for_im

!==============================================================================

subroutine IE_interpolate(iBlock,ColatLim)

  implicit none

  character (len=*),parameter :: NameSub='IE_interpolate'

  integer, intent(in) :: iBlock     ! Block index 1 for north, 2 for south
  real,    intent(in) :: ColatLim   ! Colatitude limit for the hemisphere

  call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

end subroutine IE_interpolate

!==============================================================================

subroutine initialize_ie_ua_buffers(iOutputError)

  implicit none

  integer :: iOutputError

  character (len=*),parameter :: NameSub='initialize_ie_ua_buffers'

  call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

end subroutine initialize_ie_ua_buffers

!==============================================================================

subroutine IE_put_from_UA(Buffer_III, iBlock, nMLTs, nLats, nVarsToPass)

  implicit none

  integer, intent(in) :: nMlts, nLats, iBlock, nVarsToPass
  real, dimension(nMlts, nLats, nVarsToPass), intent(in) :: Buffer_III

  character (len=*),parameter :: NameSub='IE_put_from_UA'

  call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

end subroutine IE_put_from_UA

!==============================================================================

subroutine IE_get_for_ua(Buffer_II,iSize,jSize,NameVar,NameHem,tSimulation)

  implicit none

  integer,          intent(in)  :: iSize,jSize
  real,             intent(out) :: Buffer_II(iSize,jSize)
  character (len=*),intent(in)  :: NameVar
  character (len=*),intent(in)  :: NameHem
  real,             intent(in)  :: tSimulation

  character (len=*),parameter :: NameSub='IE_get_for_ua'

  call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

end subroutine IE_get_for_ua

!==============================================================================

subroutine SPS_put_into_ie(Buffer_II, iSize, jSize, NameVar, iBlock)

  implicit none

  integer, intent(in)           :: iSize,jSize
  real, intent(in)              :: Buffer_II(iSize,jSize)
  character (len=*),intent(in)  :: NameVar
  integer,intent(in)            :: iBlock

  character (len=*), parameter :: NameSub='SPS_put_into_ie'

  call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

end subroutine SPS_put_into_ie

!==============================================================================

subroutine ionosphere_fac(iBlock) !srcIE/iono_coupling.f90
  implicit none
  integer,intent(in)::iBlock
  character (len=*), parameter :: NameSub='ionosphere_fac'

  call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

end subroutine ionosphere_fac

!==============================================================================

