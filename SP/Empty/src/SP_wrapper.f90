subroutine SP_run(TimeSimulation,TimeSimulationLimit)
  implicit none
  real,intent(inout)::TimeSimulation
  real,intent(in)::TimeSimulationLimit
  call CON_stop('Can not call SP_run')
end subroutine SP_run
!========================================================================
!======================================================================
subroutine SP_init_session(iSession,TimeSimulation)
  implicit none
  integer,  intent(in) :: iSession         ! session number (starting from 1)
  real,     intent(in) :: TimeSimulation   ! seconds from start time
  call CON_stop('Can not call SP_init_session')
end subroutine SP_init_session
!======================================================================
subroutine SP_finalize(TimeSimulation)
  implicit none
  real,intent(in)::TimeSimulation
  call CON_stop('Can not call SP_finalize')
end subroutine SP_finalize
!=========================================================
subroutine SP_set_param(CompInfo,TypeAction)
  use CON_comp_info
  implicit none
  type(CompInfoType),intent(inout)       :: CompInfo
  character(len=*), intent(in)           :: TypeAction
  !-------------------------------------------------------------------------
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
!=========================================================
subroutine SP_save_restart(TimeSimulation) 
  implicit none
  real,     intent(in) :: TimeSimulation 
  call CON_stop('Can not call SP_save restart')
end subroutine SP_save_restart
!=========================================================
subroutine SP_put_input_time(TimeIn)
  implicit none
  real,     intent(in)::TimeIn
  call CON_stop('Can not call SP_get_input_time')
end subroutine SP_put_input_time
!===================================================================
subroutine SP_put_from_mh(nPartial,iPutStart,Put,W,DoAdd,Buff_I,nVar)
  use CON_router, ONLY: IndexPtrType, WeightPtrType
  implicit none
  integer,intent(in)::nPartial,iPutStart,nVar
  type(IndexPtrType),intent(in)::Put
  type(WeightPtrType),intent(in)::W
  logical,intent(in)::DoAdd
  real,dimension(nVar),intent(in)::Buff_I
  call CON_stop('Can not put ih data')
end subroutine SP_put_from_mh
!===================================================================
subroutine SP_get_line_param(DsOut,XyzOut_D,DSCOut,DIHOut)
  implicit none
  real,intent(out)::DsOut,XyzOut_D,DSCOut,DIHOut
  call CON_stop('Can not get line parameters from SP')
end subroutine SP_get_line_param
