!^CFG COPYRIGHT UM
! Wrapper for an "empty" Solar Corona (PC) component
!==========================================================================
subroutine PC_set_param(CompInfo, TypeAction)

  use CON_comp_info

  implicit none

  character (len=*), parameter :: NameSub='PC_set_param'

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
     call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')
  end select

end subroutine PC_set_param

!==============================================================================

subroutine PC_init_session(iSession, TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  integer,  intent(in) :: iSession         ! session number (starting from 1)
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='PC_init_session'

  call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')

end subroutine PC_init_session

!==============================================================================

subroutine PC_finalize(TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='PC_finalize'

  call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')

end subroutine PC_finalize

!==============================================================================

subroutine PC_save_restart(TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='PC_save_restart'

  call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')

end subroutine PC_save_restart

!==============================================================================

subroutine PC_run(TimeSimulation,TimeSimulationLimit)

  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  real, intent(inout) :: TimeSimulation   ! current time of component

  !INPUT ARGUMENTS:
  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded

  character(len=*), parameter :: NameSub='PC_run'

  call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')

end subroutine PC_run

!===============================================================

subroutine PC_synchronize_refinement(iProc0,iCommUnion)

  implicit none
  integer, intent(in) ::iProc0,iCommUnion
  character(len=*), parameter :: NameSub='PC_synchronize_refinement'

  call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')

end subroutine PC_synchronize_refinement

!===============================================================

subroutine PC_get_for_gm(&
     nPartial,iGetStart,Get,W,State_V,nVar,TimeCoupling)

  ! derived type parameters, it is easier not to declare them
  character(len=*), parameter :: NameSub='PC_get_for_gm'

  call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')

end subroutine PC_get_for_gm
!===================================================================!
subroutine PC_get_for_mh(&
     nPartial,iGetStart,Get,W,State_V,nVar)
! derived type parameters, it is easier not to declare them
  character(len=*), parameter :: NameSub='PC_get_for_mh'

  call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')
end subroutine PC_get_for_mh
!===================================================================!
subroutine PC_get_for_mh_with_xyz(&
     nPartial,iGetStart,Get,W,State_V,nVar)
! derived type parameters, it is easier not to declare them
  character(len=*), parameter :: NameSub='PC_get_for_mh_with_xyz'

  call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')
end subroutine PC_get_for_mh_with_xyz
!===================================================================!
subroutine PC_set_buffer_grid(DD,CompID_)
  character(len=*), parameter :: NameSub='PC_set_buffer_grid'

  call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')
end subroutine PC_set_buffer_grid
!===================================================================!
subroutine PC_get_for_sp(&
     nPartial,iGetStart,Get,W,State_V,nVar)
  use CON_router, ONLY: IndexPtrType, WeightPtrType
  implicit none

  !INPUT ARGUMENTS:
  integer,intent(in)::nPartial,iGetStart,nVar
  type(IndexPtrType),intent(in)::Get
  type(WeightPtrType),intent(in)::W
  real,dimension(nVar),intent(out)::State_V

  character(len=*), parameter :: NameSub='PC_get_for_sp'

  call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')
end subroutine PC_get_for_sp
!===================================================================!
subroutine PC_get_a_line_point(&
     nPartial,iGetStart,Get,W,State_V,nVar)
  use CON_router, ONLY: IndexPtrType, WeightPtrType
  implicit none
  
  !INPUT ARGUMENTS:
  integer,intent(in)::nPartial,iGetStart,nVar
  type(IndexPtrType),intent(in)::Get
  type(WeightPtrType),intent(in)::W
  real,dimension(nVar),intent(out)::State_V
  
  character(len=*), parameter :: NameSub='PC_get_a_line_point'

  call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')
end subroutine PC_get_a_line_point
!===================================================================!
subroutine PC_put_from_mh(nPartial,&
     iPutStart,&
     Put,& 
     Weight,&
     DoAdd,&
     StateSI_V,&
     nVar)
  !USES:
  integer,intent(in)::nPartial,iPutStart,nVar
  logical,intent(in)::DoAdd
  real,dimension(nVar),intent(in)::StateSI_V

  character (len=*), parameter :: NameSub='PC_put_from_mh.f90'

  call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')
end subroutine PC_put_from_mh
!===================================================================!
subroutine PC_match_ibc
  implicit none
  character(len=*), parameter :: NameSub='PC_match_ibc'

  call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')

end subroutine PC_match_ibc
!===============================================================
subroutine PC_set_buffer_grid_get_info(CompID_, &
     nR, nPhi, nTheta, BufferMinMax_DI)

  implicit none

  integer, intent(in)     :: CompID_
  integer, intent(out)    :: nR, nPhi, nTheta
  real, intent(out)       :: BufferMinMax_DI(3,2)

  character(len=*), parameter :: NameSub = 'PC_set_buffer_grid_get_info'
  ! ---------------------------------------------------------------

  call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')

end subroutine PC_set_buffer_grid_get_info
!===============================================================
subroutine PC_save_global_buffer(nVar, nR, nPhi, nTheta, BufferIn_VG)

  implicit none

  integer,intent(in) :: nVar, nR, nPhi, nTheta
  real,intent(in)    :: BufferIn_VG(nVar, nR, nPhi, nTheta)

  character(len=*), parameter :: NameSub = 'PC_save_global_buffer'
  !-------------------------------------------------------------

 call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')

end subroutine PC_save_global_buffer
!===============================================================
subroutine PC_get_for_global_buffer(&
     nR, nPhi,nTheta, BufferMinMax_DI, Buffer_VG)

  character (len=*), parameter :: NameSub='PC_get_for_buffer_grid'
  !-------------------------------------------------------------

  call CON_stop(NameSub//': PC_ERROR: empty version cannot be used!')

end subroutine PC_get_for_global_buffer
