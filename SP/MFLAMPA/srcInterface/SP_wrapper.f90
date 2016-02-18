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
  public:: SP_request_line
  public:: SP_put_line


  ! variables requested via coupling: coordinates, 
  ! field line and particles indexes
  character(len=*), parameter:: NameVarRequest = 'xx yy zz fl id'
  integer,          parameter:: nVarRequest = 5


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

  subroutine SP_request_line(NameVar, nVar, iDirIn, CoordOut_DA)
    use ModSize, ONLY: nDim, nNode, R_, Lat_, Lon_
    use ModMain, ONLY: iGrid_IA, State_VIB, iNode_B,&
         Proc_, Block_, Begin_, End_, iProc, iComm, nBlock
    use ModMpiOrig
    ! request coordinates of field lines' beginning/origin/end
    ! as well as names variables to be imported
    !---------------------------------------------------------------
    character(len=*), intent(out):: NameVar
    integer,          intent(out):: nVar
    integer,          intent(in) :: iDirIn
    real,             intent(out):: CoordOut_DA(nDim, nNode)

    ! directions requested
    integer, parameter:: iDirBegin_ = -1, iDirOrigin_ = 0, iDirEnd_ = 1

    ! loop variables
    integer:: iParticle, iBlock, iNode
    ! indices of the particle
    integer:: iLine, iIndex
    integer:: iMin_A(nNode),iMax_A(nNode)
    integer:: iError
    character(len=*), parameter:: NameSub='SP_request_line'
    !----------------------------------------------------------------
    ! indicate variables requested
    NameVar = NameVarRequest
    nVar    = nVarRequest
    ! each processor fills only its own nodes; reset all
    CoordOut_DA = 0
    select case(iDirIn)
    case(iDirBegin_)
       ! get coordinates of the 1st points on field lines
       do iBlock = 1, nBlock
          iNode = iNode_B(iBlock); iParticle = iGrid_IA(Begin_, iNode)
          CoordOut_DA(:, iNode) = &
               State_VIB((/R_,Lat_,Lon_/), iParticle, iBlock)
       end do
    case(iDirOrigin_)
       ! get coordinates of the origin points of field lines
       do iBlock = 1, nBlock
          CoordOut_DA(:, iNode_B(iBlock)) = &
               State_VIB((/R_,Lat_,Lon_/), 0, iBlock)
       end do
    case(iDirEnd_)
       ! get coordinates of the last points on field lines
       do iBlock = 1, nBlock
          iNode = iNode_B(iBlock); iParticle = iGrid_IA(End_, iNode)
          CoordOut_DA(:, iNode_B(iBlock)) = &
               State_VIB((/R_,Lat_,Lon_/), iParticle, iBlock)
       end do
    case default
       call CON_stop(NameSub//': invalid request of field line coordinates')
    end select
    !\
    ! Collect all coords on the root
    !/
    if(iProc==0)then
       call MPI_Reduce(MPI_IN_PLACE, CoordOut_DA, nDim*nNode, MPI_REAL, &
            MPI_SUM, 0, iComm, iError)
    else
       call MPI_Reduce(CoordOut_DA, CoordOut_DA, nDim*nNode, MPI_REAL, &
            MPI_SUM, 0, iComm, iError)
    end if

  end subroutine SP_request_line

  !===================================================================
  
  subroutine SP_put_line(NameVar, nVar, nParticle, Data_VI)
    use ModSize, ONLY: nDim, nNode
    use ModMain, ONLY: iGrid_IA, State_VIB, iNode_B,&
         Proc_, Block_, Begin_, End_, iProc, iComm
    use ModMpiOrig
    ! store particle data extracted elsewhere
    !---------------------------------------------------------------
    character(len=*), intent(in):: NameVar
    integer,          intent(in):: nVar
    integer,          intent(in):: nParticle
    real,             intent(in):: Data_VI(nVar, nParticle)

    ! loop variables
    integer:: iParticle, iBlock, iNode
    ! indices of the particle
    integer:: iLine, iIndex
    integer:: iMin_A(nNode),iMax_A(nNode)
    integer:: iError
    character(len=*), parameter:: NameSub='SP_put_line'
    !----------------------------------------------------------------
    ! check correctness
    if(index(NameVar, NameVarRequest) == 0 .or. nVar /= nVarRequest)&
         call CON_stop(NameSub//': a different set variables was requested')
    ! store passed particles
    do iParticle = 1, nParticle
       iLine  = nint(Data_VI(4, iParticle))
       iIndex = nint(Data_VI(5, iParticle))
       iGrid_IA(Begin_, iLine) = MIN(iGrid_IA(Begin_,iLine), iIndex)
       iGrid_IA(Begin_, iLine) = MAX(iGrid_IA(Begin_,iLine), iIndex)
       if(iGrid_IA(Proc_, iLine) /= iProc)&
            call CON_stop(NameSub//': Incorrect message pass')
       State_VIB(1:nDim, iIndex, iGrid_IA(Block_,iLine)) = &
            Data_VI(1:nDim, iParticle)
    end do
    !\
    ! Update begin/end points on all procs
    !/
    iMin_A = iGrid_IA(Begin_, :); iMax_A = iGrid_IA(End_,:)
    call MPI_Allreduce(MPI_IN_PLACE, iMin_A, nDim*nNode, MPI_REAL, &
         MPI_MIN, iComm, iError)
    call MPI_Allreduce(MPI_IN_PLACE, iMax_A, nDim*nNode, MPI_REAL, &
         MPI_MAX, iComm, iError)
    iGrid_IA(Begin_, :) = iMin_A; iGrid_IA(End_,:) = iMax_A
  end subroutine SP_put_line

end module SP_wrapper
