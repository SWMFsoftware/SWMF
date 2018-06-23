!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!=============================================================!
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
  public:: SP_get_line_param
  public:: SP_put_from_mh
  public:: SP_put_line
  public:: SP_adjust_lines
  public:: SP_interface_point_coords
  public:: SP_get_bounds_comp
  public:: SP_put_coupling_param
  public:: SP_n_particle
  public:: SP_copy_old_state
  public:: SP_check_ready_for_mh

contains

  subroutine SP_check_ready_for_mh(IsReady)
    logical, intent(out):: IsReady
    character(len=*), parameter:: NameSub='SP_check_ready_for_mh'
    !---------------------------------------------------------------
    call CON_stop(NameSub//': cannot call the empty version')
  end subroutine SP_check_ready_for_mh
  !========================================================================  
  subroutine SP_run(TimeSimulation,TimeSimulationLimit)

    real,intent(inout)::TimeSimulation
    real,intent(in)::TimeSimulationLimit
    call CON_stop('Can not call SP_run')
  end subroutine SP_run
  !======================================================================
  subroutine SP_init_session(iSession,TimeSimulation)

    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time
    call CON_stop('Can not call SP_init_session')
  end subroutine SP_init_session
  !======================================================================
  subroutine SP_finalize(TimeSimulation)

    real,intent(in)::TimeSimulation
    call CON_stop('Can not call SP_finalize')
  end subroutine SP_finalize
  !=========================================================
  subroutine SP_set_param(CompInfo,TypeAction)
    use CON_comp_info

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

    real,     intent(in) :: TimeSimulation 
    call CON_stop('Can not call SP_save restart')
  end subroutine SP_save_restart
 
  !===================================================================
  subroutine SP_put_from_mh(nPartial,iPutStart,Put,W,DoAdd,Buff_I,nVar)
    use CON_router, ONLY: IndexPtrType, WeightPtrType

    integer,intent(in)::nPartial,iPutStart,nVar
    type(IndexPtrType),intent(in)::Put
    type(WeightPtrType),intent(in)::W
    logical,intent(in)::DoAdd
    real,dimension(nVar),intent(in)::Buff_I
    call CON_stop('Can not put mh data')
  end subroutine SP_put_from_mh
  !===================================================================
  subroutine SP_get_line_param(DsOut, XyzOut_D, DSCOut, DIHOut)

    real,intent(out):: DsOut, XyzOut_D(3), DSCOut, DIHOut
    call CON_stop('Can not get line parameters from SP')

  end subroutine SP_get_line_param
  !===================================================================

  subroutine SP_get_bounds_comp(ThisModel_, RMinOut, RMaxOut)
    ! return the values of component boundaries as set in SP component
    integer, intent(in):: ThisModel_
    real,   intent(out):: RMinOut, RMaxOut
    character(len=*), parameter:: NameSub='SP_get_solar_corona_boundary'
    !-----------------------------------------------------------------
    call CON_stop('SP: '//NameSub//' : cannot call the empty version')
  end subroutine SP_get_bounds_comp
  !===================================================================
  subroutine SP_interface_point_coords(&
       nDim, Xyz_D, nIndex, iIndex_I,&
       IsInterfacePoint)
    integer,intent(in)   :: nDim
    real,   intent(inout):: Xyz_D(nDim)
    integer,intent(in)   :: nIndex
    integer,intent(inout):: iIndex_I(nIndex)
    logical,intent(out)  :: IsInterfacePoint
    character(len=*), parameter:: NameSub='SP_interface_point_coords'
    !---------------------------------------------------------------
    call CON_stop('SP:'//NameSub//': cannot call the empty version')
  end subroutine SP_interface_point_coords
  !=======================================
  subroutine SP_put_coupling_param(iModelIn, rMinIn, rMaxIn, TimeIn,&
       rBufferLoIn, rBufferUpIn)
    real,           intent(in):: TimeIn
    integer,        intent(in):: iModelIn
    real,           intent(in):: rMinIn, rMaxIn
    real, optional, intent(in):: rBufferLoIn, rBufferUpIn
    character(len=*), parameter:: NameSub='SP_put_coupling_param'
    !---------------------------------------------------------------
    call CON_stop('SP:'//NameSub//': cannot call the empty version')
  end subroutine SP_put_coupling_param
  !===================================================================
  subroutine SP_put_line(nPartial, iPutStart, Put,&
       Weight, DoAdd, Coord_D, nVar)
    use CON_router, ONLY: IndexPtrType, WeightPtrType
    integer, intent(in) :: nPartial, iPutStart, nVar
    type(IndexPtrType), intent(in) :: Put
    type(WeightPtrType),intent(in) :: Weight
    logical,            intent(in) :: DoAdd
    real,               intent(in) :: Coord_D(nVar) !nVar=nDim
    call CON_stop('Can not put line parameters')
  end subroutine SP_put_line
  !===========================
  integer function SP_n_particle(iBlockLocal)
    integer, intent(in) :: iBlockLocal
    !----------------------------
    call CON_stop('Can not find nParticle for empty SP')
  end function SP_n_particle
  !=================================================================
  subroutine SP_copy_old_state
    character(len=*), parameter:: NameSub='SP_copy_old_state'
    !---------------------------------------------------------------
    call CON_stop('SP:'//NameSub//': cannot call the empty version')
  end subroutine SP_copy_old_state
  !===========================
  subroutine SP_adjust_lines(DoInit)
    Logical, intent(in):: DoInit
    character(len=*), parameter:: NameSub='SP_adjust_lines'
    !---------------------------------------------------------------
    call CON_stop('SP:'//NameSub//': cannot call the empty version')
  end subroutine SP_adjust_lines
end module SP_wrapper
