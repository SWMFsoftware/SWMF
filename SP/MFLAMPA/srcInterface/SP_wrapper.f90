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
  public:: SP_put_from_mh

contains

  subroutine SP_run(TimeSimulation,TimeSimulationLimit)
    use ModMain, ONLY: run
    real,intent(inout)::TimeSimulation
    real,intent(in)::TimeSimulationLimit
    !--------------------------------------------------------------------------
    call run
  end subroutine SP_run

  !========================================================================

  subroutine SP_init_session(iSession,TimeSimulation)

    use ModMain, ONLY: initialize

    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time
    !--------------------------------------------------------------------------
    call initialize
  end subroutine SP_init_session

  !======================================================================

  subroutine SP_finalize(TimeSimulation)
    
    use ModMain, ONLY: finalize

    real,intent(in)::TimeSimulation
    !--------------------------------------------------------------------------
    call finalize
  end subroutine SP_finalize

  !=========================================================

  subroutine SP_set_param(CompInfo,TypeAction)
    use CON_comp_info
    use ModMain, ONLY: read_param, iComm, iProc, nProc

    type(CompInfoType),intent(inout):: CompInfo
    character(len=*),  intent(in)   :: TypeAction

    character(len=*), parameter :: NameSub='SP_set_param'
    !-------------------------------------------------------------------------
    select case(TypeAction)
    case('VERSION')
       call put(CompInfo,&
            Use        =.true., &
            NameVersion='Empty', &
            Version    =0.0)
    case('MPI')
       call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc)
       if(nProc/=1)call CON_stop(&
            'The present SEP version can not use more than 1 PE')
    case('STDOUT')
       ! placeholder
    case('CHECK')
       ! placeholder
    case('READ')
       call read_param(TypeAction)
    case('GRID')
       call SP_set_grid
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
  subroutine SP_put_from_mh(nPartial,iPutStart,Put,W,DoAdd,Buff_I,nVar)
    use CON_router, ONLY: IndexPtrType, WeightPtrType

    integer,intent(in)::nPartial,iPutStart,nVar
    type(IndexPtrType),intent(in)::Put
    type(WeightPtrType),intent(in)::W
    logical,intent(in)::DoAdd
    real,dimension(nVar),intent(in)::Buff_I
    call CON_stop('Can not put ih data')
  end subroutine SP_put_from_mh
  !===================================================================
  subroutine SP_get_line_param(DsOut, XyzOut_D, DSCOut, DIHOut)

    real,intent(out):: DsOut, XyzOut_D(3), DSCOut, DIHOut
    call CON_stop('Can not get line parameters from SP')

  end subroutine SP_get_line_param

  !===================================================================

  subroutine SP_set_grid

    use CON_coupler,    ONLY: &
         set_coord_system, &
         ! CON_coupler::CON_router::CON_grid_descriptor::CON_grid_storage::
         init_decomposition, get_root_decomposition, bcast_decomposition
    use CON_world,      ONLY: is_proc0
    use CON_comp_param, ONLY: SP_
    use ModConst,       ONLY: rSun
    use ModNumConst,    ONLY: cHalfPi, cTwoPi
    use ModSize,        ONLY: iIdMin, iIdMax

    ! Initialize 3D grid with NON-TREE structure
    call init_decomposition(&
         GridID_ = SP_,&
         CompID_ = SP_,&
         nDim    = 3)

    ! Construct decomposition
    if(is_proc0(SP_))&
         call get_root_decomposition(&
         SP_,&
         iRootMapDim_D = (/1, 1, 1/),&
         XyzMin_D      = (/real(iIdMin), 0.0,   -cHalfPi/),&
         XyzMax_D      = (/real(iIdMax), cTwoPi, cHalfPi/),&
         nCells_D      = (/1, 1, 1/))
    call bcast_decomposition(SP_)

    ! Coordinate system is Heliographic Inertial Coordinate System (HGI)
    ! with length measured in solar radii
    call set_coord_system(&
         GridID_   = SP_, &
         TypeCoord ='HGI',&
         UnitX     = rSun)
  end subroutine SP_set_grid

end module SP_wrapper
