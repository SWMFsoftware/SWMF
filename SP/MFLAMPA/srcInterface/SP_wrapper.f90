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
  public:: SP_put_input_time
  public:: SP_put_from_mh
  public:: SP_put_line
  public:: SP_get_line_origin
  public:: SP_get_interface

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
    call CON_stop('Can not put mh data')
  end subroutine SP_put_from_mh

  !===================================================================

  subroutine SP_set_grid

    use CON_coupler,    ONLY: &
         set_coord_system, &
         init_decomposition, get_root_decomposition, bcast_decomposition
    use CON_world,      ONLY: is_proc0
    use CON_comp_param, ONLY: SP_
    use ModConst,       ONLY: rSun
    use ModSize,        ONLY: nDim, nLat, nLon, &
         iParticleMin, iParticleMax, nParticle
    use ModMain,        ONLY: LatMin, LatMax, LonMin, LonMax, &
         iGrid_IA, Block_, Proc_

    ! Initialize 3D grid with NON-TREE structure
    call init_decomposition(&
         GridID_ = SP_,&
         CompID_ = SP_,&
         nDim    = nDim)

    ! Construct decomposition
    if(is_proc0(SP_))&
         call get_root_decomposition(&
         GridID_       = SP_,&
         iRootMapDim_D = (/1, nLat, nLon/),&
         XyzMin_D      = (/real(iParticleMin), LatMin, LonMin/),&
         XyzMax_D      = (/real(iParticleMax), LatMax, LonMax/),&
         nCells_D      = (/nParticle , 1, 1/),&
         PE_I          = iGrid_IA(Proc_,:),&
         iBlock_I      = iGrid_IA(Block_,:))
    call bcast_decomposition(SP_)

    ! Coordinate system is Heliographic Inertial Coordinate System (HGI)
    ! with length measured in solar radii
    call set_coord_system(&
         GridID_   = SP_, &
         TypeCoord ='HGI',&
         UnitX     = rSun)
  end subroutine SP_set_grid

  !===================================================================
  
  subroutine SP_get_line_origin(CoordOriginOut_DA)
    use ModSize, ONLY: nDim, nNode
    use ModMain, ONLY: CoordOrigin_DA
    real, intent(out):: CoordOriginOut_DA(nDim, nNode)
    !------------------------------------------------
    CoordOriginOut_DA = CoordOrigin_DA
  end subroutine SP_get_line_origin

  !===================================================================
  
  subroutine SP_get_interface(CoordInterfaceOut_DA)
    use ModSize, ONLY: nDim, nNode, iParticleMin, iParticleMax
    use ModMain, ONLY: State_VIB, iNode_B, nBlock
    ! get the last points of the field lines;
    ! called after recv'd lines extracted from SC,
    ! thus last points are at interface between SC and IH
    real, intent(out):: CoordInterfaceOut_DA(nDim, nNode)
    ! loop variable
    integer:: iBlock
    ! index of the last on the current line
    integer:: iParticleLast
    character(len=*), parameter:: NameSub='SP_get_interface'
    !------------------------------------------------
    ! each proc fills its own part
    CoordInterfaceOut_DA = 0.0
    do iBlock = 1, nBlock
       ! find the last particle on this line/block
       iParticleLast = iParticleMax
       do while(State_VIB(1, iParticleLast, iBlock) < 0)
          if(iParticleLast == iParticleMin)&
               call MPI_stop(NameSub//': Line is not found')
          iParticleLast = iParticleLast - 1
       end do
       ! get its HGI coordinates
       CoordInterfaceOut_DA(:, iNode_B(iBlock)) = &
            State_VIB(1:nDim, iParticleLast, iBlock)
    end do

  end subroutine SP_get_interface

  !===================================================================
  
  subroutine SP_put_line(nParticle, ParticleData_II)
    use ModSize, ONLY: nDim
    use ModMain, ONLY: iGrid_IA, State_VIB, Block_, Proc_, iProc
    ! store field lines data extracted elsewhere
    integer, intent(in):: nParticle
    real,    intent(in):: ParticleData_II(nDim+2, nParticle)
    
    ! loop variable
    integer:: iParticle
    ! indices of the particle
    integer:: iLine, iIndex
    character(len=*), parameter:: NameSub='SP_put_line'
    !----------------------------------------------------------------
    do iParticle = 1, nParticle
       iLine  = nint(ParticleData_II(1, iParticle))
       iIndex = nint(ParticleData_II(2, iParticle))
       if(iGrid_IA(Proc_, iLine) /= iProc)&
            call MPI_stop(NameSub//': Incorrect message pass')
       State_VIB(1:nDim, iIndex, iGrid_IA(Block_,iLine)) = &
            ParticleData_II(3:2+nDim, iParticle)
    end do
  end subroutine SP_put_line

end module SP_wrapper
