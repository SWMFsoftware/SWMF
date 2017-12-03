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
  public:: SP_put_input_time
  public:: SP_put_from_sc
  public:: SP_put_from_ih
  public:: SP_put_line
  public:: SP_synchronize_grid
  public:: SP_get_cell_index
  public:: SP_get_particle_index
  public:: SP_adjust_lines
  public:: SP_interface_point_coords_for_ih
  public:: SP_interface_point_coords_for_ih_extract
  public:: SP_interface_point_coords_for_sc
  public:: SP_get_domain_boundary
  public:: SP_put_r_min
  public:: SP_n_particle
  public:: SP_copy_old_state
  public:: SP_do_extract
  public:: SP_assign_lagrangian_coords
contains

  subroutine SP_do_extract(DoExtract)
    logical, intent(out):: DoExtract
    character(len=*), parameter:: NameSub='SP_do_extract'
    !---------------------------------------------------------------
    call CON_stop('SP:'//NameSub//': cannot call the empty version')
  end subroutine SP_do_extract
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
  !=========================================================
  subroutine SP_put_input_time(TimeIn)

    real,     intent(in)::TimeIn
    call CON_stop('Can not call SP_get_input_time')
  end subroutine SP_put_input_time
  !===================================================================
  subroutine SP_put_from_sc(nPartial,iPutStart,Put,W,DoAdd,Buff_I,nVar)
    use CON_router, ONLY: IndexPtrType, WeightPtrType

    integer,intent(in)::nPartial,iPutStart,nVar
    type(IndexPtrType),intent(in)::Put
    type(WeightPtrType),intent(in)::W
    logical,intent(in)::DoAdd
    real,dimension(nVar),intent(in)::Buff_I
    call CON_stop('Can not put ih data')
  end subroutine SP_put_from_sc
  !===================================================================
  subroutine SP_put_from_ih(nPartial,iPutStart,Put,W,DoAdd,Buff_I,nVar)
    use CON_router, ONLY: IndexPtrType, WeightPtrType

    integer,intent(in)::nPartial,iPutStart,nVar
    type(IndexPtrType),intent(in)::Put
    type(WeightPtrType),intent(in)::W
    logical,intent(in)::DoAdd
    real,dimension(nVar),intent(in)::Buff_I
    call CON_stop('Can not put ih data')
  end subroutine SP_put_from_ih
  !===================================================================
  subroutine SP_get_line_param(DsOut, XyzOut_D, DSCOut, DIHOut)

    real,intent(out):: DsOut, XyzOut_D(3), DSCOut, DIHOut
    call CON_stop('Can not get line parameters from SP')

  end subroutine SP_get_line_param
  !===================================================================

  subroutine SP_get_domain_boundary(RScOut, RIhOut)
    ! return the value of the solar corona boundary as set in SP component
    real, intent(out):: RScOut, RIhOut
    character(len=*), parameter:: NameSub='SP_get_solar_corona_boundary'
    !-----------------------------------------------------------------
    call CON_stop('SP: '//NameSub//' : cannot call the empty version')
  end subroutine SP_get_domain_boundary

  !===================================================================

  subroutine SP_put_r_min(R)
    real, intent(in)::R
    character(len=*), parameter:: NameSub='SP_put_r_min'
    call CON_stop('SP:'//NameSub//': cannot call the empty version')
  end subroutine SP_put_r_min

  !===================================================================
  subroutine SP_interface_point_coords_for_sc(&
       GridDescriptor, iBlockUsed, nDim, Xyz_D, nIndex,iIndex_I,&
       IsInterfacePoint)
    use CON_grid_descriptor
    type(LocalGDType),intent(in)::GridDescriptor
    integer,intent(in)   :: iBlockUsed
    integer,intent(in)   :: nDim
    real,   intent(inout):: Xyz_D(nDim)
    integer,intent(in)   :: nIndex
    integer,intent(inout):: iIndex_I(nIndex)
    logical,intent(out)  :: IsInterfacePoint
    character(len=*), parameter:: NameSub='SP_interface_point_coords_for_sc'
    !---------------------------------------------------------------
    call CON_stop('SP:'//NameSub//': cannot call the empty version')
  end subroutine SP_interface_point_coords_for_sc
  !===================================================================
  subroutine SP_interface_point_coords_for_ih(&
       GridDescriptor, iBlockUsed, nDim, Xyz_D, nIndex, iIndex_I,&
       IsInterfacePoint)
    use CON_grid_descriptor
    type(LocalGDType),intent(in)::GridDescriptor
    integer,intent(in)   :: iBlockUsed
    integer,intent(in)   :: nDim
    real,   intent(inout):: Xyz_D(nDim)
    integer,intent(in)   :: nIndex
    integer,intent(inout):: iIndex_I(nIndex)
    logical,intent(out)  :: IsInterfacePoint
    character(len=*), parameter:: NameSub='SP_interface_point_coords_for_ih'
    !---------------------------------------------------------------
    call CON_stop('SP:'//NameSub//': cannot call the empty version')
  end subroutine SP_interface_point_coords_for_ih
  !===================================================================
  subroutine SP_interface_point_coords_for_ih_extract(&
       GridDescriptor, iBlockUsed, nDim, Xyz_D, nIndex, iIndex_I,&
       IsInterfacePoint)
    use CON_grid_descriptor
    type(LocalGDType),intent(in)::GridDescriptor
    integer,intent(in)   :: iBlockUsed
    integer,intent(in)   :: nDim
    real,   intent(inout):: Xyz_D(nDim)
    integer,intent(in)   :: nIndex
    integer,intent(inout):: iIndex_I(nIndex)
    logical,intent(out)  :: IsInterfacePoint
    character(len=*), parameter:: NameSub='SP_interface_point_coords_for_ih'
    !---------------------------------------------------------------
    call CON_stop('SP:'//NameSub//': cannot call the empty version')
  end subroutine SP_interface_point_coords_for_ih_extract

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
  !========================================================================
  subroutine SP_copy_old_state
    character(len=*), parameter:: NameSub='SP_copy_old_state'
    !---------------------------------------------------------------
    call CON_stop('SP:'//NameSub//': cannot call the empty version')
  end subroutine SP_copy_old_state
  !===========================
  subroutine SP_adjust_lines(iComp)
    integer, intent(in):: iComp
    character(len=*), parameter:: NameSub='SP_adjust_lines'
    !---------------------------------------------------------------
    call CON_stop('SP:'//NameSub//': cannot call the empty version')
  end subroutine SP_adjust_lines
  !===========================
  subroutine SP_synchronize_grid(iComm)
    integer, intent(in):: iComm
    character(len=*), parameter:: NameSub='SP_synchronize_grid'
    !---------------------------------------------------------------
    call CON_stop('SP:'//NameSub//': cannot call the empty version')
  end subroutine SP_synchronize_grid
  !===========================
  subroutine SP_get_particle_index(i1,i2,i3)
    integer, intent(in) :: i1,i2
    integer, intent(out):: i3
    character(len=*), parameter:: NameSub='SP_get_particle_index'
    !---------------------------------------------------------------
    call CON_stop('SP:'//NameSub//': cannot call the empty version')
  end subroutine SP_get_particle_index
  !===========================
  subroutine SP_get_cell_index(i1,i2,i3)
    integer, intent(in) :: i1,i2
    integer, intent(out):: i3
    character(len=*), parameter:: NameSub='SP_get_cell_index'
    !---------------------------------------------------------------
    call CON_stop('SP:'//NameSub//': cannot call the empty version')
  end subroutine SP_get_cell_index
  !===========================
  subroutine SP_assign_lagrangian_coords
    character(len=*), parameter:: NameSub='assign_lagrangian_coords'
    !---------------------------------------------------------------
    call CON_stop('SP_'//NameSub//': cannot call the empty version')
  end subroutine SP_assign_lagrangian_coords
end module SP_wrapper
