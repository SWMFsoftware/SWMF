!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!=============================================================!
module SP_wrapper
 
  use ModNumConst, ONLY: cHalfPi
 
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
  public:: SP_get_grid_descriptor_param
  public:: SP_get_line_all

  ! variables requested via coupling: coordinates, 
  ! field line and particles indexes
  character(len=*), parameter:: NameVarRequest = 'xx yy zz fl id'
  integer,          parameter:: nVarRequest = 5


contains

  subroutine SP_run(TimeSimulation,TimeSimulationLimit)
    use ModMain, ONLY: advance, TimeGlobal
    real,intent(inout)::TimeSimulation
    real,intent(in)::TimeSimulationLimit
    !--------------------------------------------------------------------------
    TimeGlobal = TimeSimulation
    call advance
    TimeSimulation = TimeSimulationLimit
  end subroutine SP_run

  !========================================================================

  subroutine SP_init_session(iSession,TimeSimulation)

    use ModMain, ONLY: initialize

    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    logical, save:: IsInitialized = .false.
    !--------------------------------------------------------------------------
    if(IsInitialized)&
         RETURN
    IsInitialized = .true.
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
    use ModMain, ONLY: check, read_param, iComm, iProc, nProc

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
       call check
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
    use ModCoordTransform, ONLY: xyz_to_rlonlat
    use ModMain, ONLY: State_VIB
    use ModConst, ONLY: rSun

    integer,intent(in)::nPartial,iPutStart,nVar
    type(IndexPtrType),intent(in)::Put
    type(WeightPtrType),intent(in)::W
    logical,intent(in)::DoAdd
    real,dimension(nVar),intent(in)::Buff_I

    real:: Xyz_D(3), Coord_D(3)
    real:: B_D(3)

    integer:: i, j, k, iBlock
    !------------------------------------------------------------
    i      = Put%iCB_II(1,iPutStart)
    j      = Put%iCB_II(2,iPutStart)
    k      = Put%iCB_II(3,iPutStart)
    iBlock = Put%iCB_II(4,iPutStart)

    Xyz_D = Buff_I((nVar-2):nVar)

    B_D = Buff_I(5:7)

    ! convert from SI
    Xyz_D = Xyz_D / rSun 
    call xyz_to_rlonlat(Xyz_D, Coord_D)


    if(DoAdd)then
!       State_VIB(1:3,i,iBlock) = State_VIB(1:3,i,iBlock) + Coord_D
       State_VIB(4:6,i,iBlock) = State_VIB(4:6,i,iBlock) + B_D
    else
!       State_VIB((/1,2,3/),i,iBlock) = Coord_D
       State_VIB((/4,5,6/),i,iBlock) = B_D
    end if
  end subroutine SP_put_from_mh

  !===================================================================

  subroutine SP_set_grid

    use CON_coupler,    ONLY: &
         set_coord_system, &
         init_decomposition, get_root_decomposition, bcast_decomposition
    use CON_world,      ONLY: is_proc0
    use CON_comp_param, ONLY: SP_
    use ModConst,       ONLY: rSun
    use ModMain,        ONLY: &
         LatMin, LatMax, LonMin, LonMax, &
         iGrid_IA, Block_, Proc_, &
         nDim, nLat, nLon, &
         iParticleMin, iParticleMax, nParticle,&
         TypeCoordSystem

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
         GridID_      = SP_, &
         TypeCoord    = TypeCoordSystem, &
         TypeGeometry = 'spherical', &
         UnitX        = rSun)
  end subroutine SP_set_grid

  !===================================================================

  subroutine SP_request_line(iDirIn, nLine, CoordOut_DI, iIndexOut_II)
    use ModMain, ONLY: &
         iGrid_IA, State_VIB, iNode_B,&
         Proc_, Block_, Begin_, End_, iProc, iComm, nBlock, &
         nDim, nNode, R_, Lat_, Lon_, &
         get_node_indexes
    ! request coordinates & indices of field lines' beginning/origin/end
    ! for the current processor
    !---------------------------------------------------------------
    integer,              intent(in) :: iDirIn
    integer,              intent(out):: nLine
    real,    allocatable, intent(out):: CoordOut_DI(:, :)
    integer, allocatable, intent(out):: iIndexOut_II(:,:)

    ! directions requested
    integer, parameter:: iDirBegin_ = -1, iDirOrigin_ = 0, iDirEnd_ = 1

    ! radius-lon-lat coordinates
    real:: Coord_D(nDim)
    ! loop variables
    integer:: iParticle, iBlock, iNode
    ! indices of the particle
    integer:: iLine, iIndex
    integer:: iMin_A(nNode),iMax_A(nNode)
    integer:: iError
    character(len=*), parameter:: NameSub='SP_request_line'
    !----------------------------------------------------------------
    ! number of lines on this processor
    nLine = nBlock

    ! prepare containers to hold the request
    if(allocated(CoordOut_DI)) deallocate(CoordOut_DI)
    allocate(CoordOut_DI(nDim, nBlock))
    if(allocated(iIndexOut_II)) deallocate(iIndexOut_II)
    allocate(iIndexOut_II(nDim+1, nBlock))! 3 cell + 1 block index

    ! each processor fills only its own nodes
    select case(iDirIn)
    case(iDirBegin_)
       ! get coordinates of the 1st points on field lines
       do iBlock = 1, nBlock
          iNode = iNode_B(iBlock); iParticle = iGrid_IA(Begin_, iNode)
          CoordOut_DI(:, iBlock) = &
               State_VIB((/R_,Lon_,Lat_/), iParticle, iBlock)
          iIndexOut_II(1, iBlock) = iParticle
          call get_node_indexes(iNode, &
               iIndexOut_II(2, iBlock), iIndexOut_II(3, iBlock))
          iIndexOut_II(4, iBlock) = iBlock
       end do
    case(iDirOrigin_)
       ! get coordinates of the origin points of field lines
       do iBlock = 1, nBlock
          iNode = iNode_B(iBlock)
          CoordOut_DI(:, iBlock) = &
               State_VIB((/R_,Lon_,Lat_/), 0, iBlock)
          iIndexOut_II(1, iBlock) = 0
          call get_node_indexes(iNode, &
               iIndexOut_II(2, iBlock), iIndexOut_II(3, iBlock))
          iIndexOut_II(4, iBlock) = iBlock
       end do
    case(iDirEnd_)
       ! get coordinates of the last points on field lines
       do iBlock = 1, nBlock
          iNode = iNode_B(iBlock); iParticle = iGrid_IA(End_, iNode)
          CoordOut_DI(:, iBlock) = &
               State_VIB((/R_,Lon_,Lat_/), iParticle, iBlock)
          iIndexOut_II(1, iBlock) = 0
          call get_node_indexes(iNode, &
               iIndexOut_II(2, iBlock), iIndexOut_II(3, iBlock))
          iIndexOut_II(4, iBlock) = iBlock
       end do
    case default
       call CON_stop(NameSub//': invalid request of field line coordinates')
    end select
  end subroutine SP_request_line

  !===================================================================

  subroutine SP_put_line(nParticle, Coord_DI, iIndex_II, iDirIn)
    use ModMain, ONLY: &
         iGrid_IA, State_VIB, iNode_B,&
         Proc_, Block_, Begin_, End_, iProc, iComm, &
         nDim, nNode, iParticleMin, iParticleMax, Lat_, Lon_, R_
    use ModMpi
    ! store particle coordinates extracted elsewhere
    !---------------------------------------------------------------
    integer, intent(in):: nParticle
    real,    intent(in):: Coord_DI( nDim,   nParticle)
    integer, intent(in):: iIndex_II(nDim+1, nParticle)
    integer, intent(in):: iDirIn

    ! directions where to put particles
    integer, parameter:: iDirBegin_ = -1, iDirOrigin_ = 0, iDirEnd_ = 1

    ! cartesian coordinates
    real:: Xyz_D(nDim)
    ! radius-lon-lat coordinates
    real:: Coord_D(nDim)
    ! loop variables
    integer:: iParticle, iBlock, iNode
    ! indices of the particle
    integer:: iLine, iIndex
    integer:: iMin_A(nNode),iMax_A(nNode), iOffset_A(nNode)
    integer:: iError
    character(len=*), parameter:: NameSub='SP_put_line'
    !----------------------------------------------------------------
    ! list of min/max index of active particles is needed if iDirIn /= 0,
    ! i.e. need to add particle to the beginning/end of list
    select case(iDirIn)
    case(iDirBegin_)
       iOffset_A = iGrid_IA(Begin_,:)
    case(iDirOrigin_)
       iOffset_A = 0
    case(iDirEnd_)
       iOffset_A = iGrid_IA(End_, :)
    case default
       call CON_stop(NameSub//': invalid call')
    end select
    ! store passed particles
    do iParticle = 1, nParticle
       iBlock  = iIndex_II(4, iParticle)
       iLine = iNode_B(iBlock)
       iIndex = iIndex_II(1, iParticle) + iOffset_A(iLine)
       if(iIndex < iParticleMin)&
            call CON_stop(NameSub//': particle index is below limit')
       if(iIndex > iParticleMax)&
            call CON_stop(NameSub//': particle index is above limit')
       iGrid_IA(Begin_, iLine) = MIN(iGrid_IA(Begin_,iLine), iIndex)
       iGrid_IA(End_,   iLine) = MAX(iGrid_IA(End_,  iLine), iIndex)
       if(iGrid_IA(Proc_, iLine) /= iProc)&
            call CON_stop(NameSub//': Incorrect message pass')
       ! store data
       State_VIB((/R_, Lon_, Lat_/), iIndex, iBlock) = &
            Coord_DI(1:nDim, iParticle)
    end do
    !\
    ! Update begin/end points on all procs
    !/
    iMin_A = iGrid_IA(Begin_, :); iMax_A = iGrid_IA(End_,:)
    call MPI_Allreduce(MPI_IN_PLACE, iMin_A, nNode, MPI_INTEGER, &
         MPI_MIN, iComm, iError)
    call MPI_Allreduce(MPI_IN_PLACE, iMax_A, nNode, MPI_INTEGER, &
         MPI_MAX, iComm, iError)
    iGrid_IA(Begin_, :) = iMin_A; iGrid_IA(End_,:) = iMax_A
  end subroutine SP_put_line

  !===================================================================

  subroutine SP_get_grid_descriptor_param(&
       iGridMin_D, iGridMax_D, Displacement_D)
    use ModMain, ONLY: nDim, iParticleMin, iParticleMax
    integer, intent(out):: iGridMin_D(nDim)
    integer, intent(out):: iGridMax_D(nDim)
    real,    intent(out):: Displacement_D(nDim)
    !-----------------------------------------
    iGridMin_D = (/iParticleMin, 1, 1/)
    iGridMax_D = (/iParticleMax, 1, 1/)
    Displacement_D = 0.0
  end subroutine SP_get_grid_descriptor_param

  !===================================================================

  subroutine SP_get_line_all(Xyz_DI)
    use ModMain, ONLY: iProc, iComm, Block_, Proc_, Begin_, End_,&
         iGrid_IA, State_VIB, &
         nDim, nNode, nParticle, R_, Lat_, Lon_, iParticleMin,iParticleMax
    use ModCoordTransform, ONLY: rlonlat_to_xyz
    use ModMpi
    real, pointer:: Xyz_DI(:, :)

    integer:: iNode, iParticle, iBlock, iError
    ! radius-lon-lat coordinates
    real:: Coord_D(nDim)
    !-----------------------------------------
    Xyz_DI = 0.0
    do iNode = 1, nNode
       if(iGrid_IA(Proc_, iNode) /= iProc)&
            CYCLE
       iBlock = iGrid_IA(Block_, iNode)
       do iParticle = iParticleMin, iParticleMax
          if(  iParticle < iGrid_IA(Begin_, iNode) .or. &
               iParticle > iGrid_IA(End_,   iNode)) &
               CYCLE
          Coord_D = State_VIB((/R_,Lon_,Lat_/), iParticle, iBlock)
          call rlonlat_to_xyz(Coord_D, &
               Xyz_DI(:, (iNode-1)*nParticle+iParticle-iParticleMin+1) )
       end do
    end do
    call MPI_Allreduce(MPI_IN_PLACE, Xyz_DI, nParticle*nNode*nDim, MPI_REAL, &
         MPI_SUM, iComm, iError)

 end subroutine SP_get_line_all

end module SP_wrapper
