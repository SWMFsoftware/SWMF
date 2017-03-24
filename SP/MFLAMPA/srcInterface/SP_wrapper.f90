!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!=============================================================!
module SP_wrapper

  use ModNumConst, ONLY: cHalfPi
  use ModConst, ONLY: rSun, cProtonMass
  use ModCoordTransform, ONLY: xyz_to_rlonlat, rlonlat_to_xyz
  use SP_ModMain, ONLY: &
       run, initialize, finalize, check, read_param,&
       get_node_indexes, &
       iComm, iProc, nProc, &
       nDim, nNode, nLat, nLon, nBlock,&
       iParticleMin, iParticleMax, nParticle,&
       RMin, RSc, LatMin, LatMax, LonMin, LonMax, &
       iGridGlobal_IA, iGridLocal_IB, State_VIB, iNode_B, TypeCoordSystem,&
       CoordMin_DI, &
       Block_, Proc_, Begin_, End_, &
       R_, Lat_, Lon_, Rho_, Bx_,By_,Bz_,B_, Ux_,Uy_,Uz_, T_, RhoOld_, BOld_
  use CON_comp_info
  use CON_router, ONLY: IndexPtrType, WeightPtrType
  use CON_coupler, ONLY: &
       set_coord_system, &
       init_decomposition, get_root_decomposition, bcast_decomposition, &
       iVar_V, DoCoupleVar_V, &
       Density_, RhoCouple_, Pressure_, PCouple_, &
       Momentum_, RhoUxCouple_, RhoUzCouple_, &
       BField_, BxCouple_, BzCouple_
  use CON_world, ONLY: is_proc0
  use CON_comp_param, ONLY: SP_, SC_, IH_

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
  public:: SP_interface_point_coords_for_ih
  public:: SP_interface_point_coords_for_sc
  public:: SP_put_line
  public:: SP_get_grid_descriptor_param
  public:: SP_get_solar_corona_boundary
  public:: SP_put_r_min

  ! variables requested via coupling: coordinates, 
  ! field line and particles indexes
  character(len=*), parameter:: NameVarCouple = 'rho p mx my mz bx by bz'

  ! particles that need to be requested from MH component:
  ! 1st index - the request's index in buffer
  ! 2nd index - content of the request
  integer, pointer:: iRequestSc_II(:,:)=>null(), iRequestIh_II(:,:)=>null()
  ! indices that define a request: field line and particle indices
  integer, parameter  :: nRequestIndex = 2
  ! current size of the request
  integer:: nRequestSc, nRequestIh
  ! whether the request has been sent
  logical:: IsSentRequestSc, IsSentRequestIh
  ! estimated number of requests per field line
  integer, parameter  :: nRequestPerLine = nParticle / 100

contains

  subroutine add_to_request(iComp, iRequest_I)
    ! auxilary function to add request to the appropriate buffer
    !-------------------------------------------------------------------
    ! component to be requested from
    integer, intent(in):: iComp
    ! request to be added
    integer, intent(in):: iRequest_I(nRequestIndex)

    character(len=100):: StringError
    character(len=*), parameter:: NameSub='SP:add_to_request'
    !------------------------------------------------------------
    ! each component has its own request buffer
    select case(iComp)
    case(SC_)
       ! add the request to SC buffer
       nRequestSc = nRequestSc + 1
       ! if the current size of the buffer is exceeded -> double its size
       if(nRequestSc >  ubound(iRequestSc_II,2)) &
            call double_buffer(iRequestSc_II)
       ! put the request to buffer
       iRequestSc_II(:,nRequestSc) = iRequest_I
    case(IH_)
       ! add the request to IH buffer
       nRequestIh = nRequestIh + 1
       ! if the current size of the buffer is exceeded -> double its size
       if(nRequestIh >  ubound(iRequestIh_II,2)) &
            call double_buffer(iRequestIh_II)
       ! put the request to buffer
       iRequestIh_II(:,nRequestIh) = iRequest_I
    case default
       ! requesting from iComp isn't implemented -> ERROR
       write(StringError,'(a,i2)') &
            ": isn't implemented for requesting data from component ", iComp
       call CON_stop(NameSub//StringError)
    end select
  contains
    subroutine double_buffer(iBuffer_II)
      ! double the size of the buffer
      integer, pointer, intent(inout):: iBuffer_II(:,:)
      ! auxilary pointer
      integer, pointer:: iAux_II(:,:)
      !------------------------------------------------------
      ! allocated memeory of doubled size
      allocate(iAux_II(nRequestIndex,2*ubound(iBuffer_II,2)))
      ! transfer data to new buffer
      iAux_II(:,1:ubound(iBuffer_II,2)) = iBuffer_II
      ! deallocated the old buffer
      deallocate(iBuffer_II)
      ! reassign pointers
      iBuffer_II => iAux_II; nullify(iAux_II)
    end subroutine double_buffer
  end subroutine add_to_request

  !========================================================================

  subroutine SP_run(TimeSimulation,TimeSimulationLimit)
    real,intent(inout)::TimeSimulation
    real,intent(in)::TimeSimulationLimit
    !--------------------------------------------------------------------------
    call run(TimeSimulationLimit)
    TimeSimulation = TimeSimulationLimit
  end subroutine SP_run

  !========================================================================

  subroutine SP_init_session(iSession,TimeSimulation)
    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    logical, save:: IsInitialized = .false.
    !--------------------------------------------------------------------------
    if(IsInitialized)&
         RETURN
    IsInitialized = .true.
     call initialize(TimeSimulation)
  end subroutine SP_init_session

  !======================================================================

  subroutine SP_finalize(TimeSimulation)
    real,intent(in)::TimeSimulation
    !--------------------------------------------------------------------------
    call finalize
  end subroutine SP_finalize

  !=========================================================

  subroutine SP_set_param(CompInfo,TypeAction)
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

    integer,intent(in)::nPartial,iPutStart,nVar
    type(IndexPtrType),intent(in)::Put
    type(WeightPtrType),intent(in)::W
    logical,intent(in)::DoAdd
    real,dimension(nVar),intent(in)::Buff_I
    integer:: iRho, iP, iMx, iMz, iBx, iBz
    integer:: i, j, k, iBlock
    integer:: iPartial
    real:: Weight
    real:: Aux

    real, external:: energy_in

    character(len=*), parameter:: NameSub='SP_put_from_mh'
    !------------------------------------------------------------
    ! check consistency: momentum and pressure are needed together with density
    if(.not. DoCoupleVar_V(Density_) .and. &
         (DoCoupleVar_V(Pressure_) .or. DoCoupleVar_V(Momentum_)))&
         call CON_Stop(NameSub//': pressure or momentum is coupled,'//&
         ' but density is not')

    ! indices of variables in the buffer
    iRho= iVar_V(RhoCouple_)
    iP  = iVar_V(PCouple_)
    iMx = iVar_V(RhoUxCouple_)
    iMz = iVar_V(RhoUzCouple_)
    iBx = iVar_V(BxCouple_)
    iBz = iVar_V(BzCouple_)   
    ! auxilary factor to account for value of DoAdd
    Aux = 0.0
    if(DoAdd) Aux = 1.0

    do iPartial = 0, nPartial-1
       ! cell and block indices
       i      = Put%iCB_II(1, iPutStart + iPartial)
       j      = Put%iCB_II(2, iPutStart + iPartial)
       k      = Put%iCB_II(3, iPutStart + iPartial)
       iBlock = Put%iCB_II(4, iPutStart + iPartial)
       ! interpolation weight
       Weight = W%Weight_I(   iPutStart + iPartial)
       ! put the data
       ! NOTE: State_VIB must be reset to zero before putting coupled data
       if(DoCoupleVar_V(Density_))&
            State_VIB(Rho_,i,iBlock) = Aux * State_VIB(Rho_,i,iBlock) + &
            Buff_I(iRho)/cProtonMass * Weight
       if(DoCoupleVar_V(Pressure_))&
            State_VIB(T_,i,iBlock) = Aux * State_VIB(T_,i,iBlock) + &
            Buff_I(iP)/Buff_I(iRho)*cProtonMass/energy_in('kev') * Weight
       if(DoCoupleVar_V(Momentum_))&
            State_VIB(Ux_:Uz_,i,iBlock) = Aux * State_VIB(Ux_:Uz_,i,iBlock) + &
            Buff_I(iMx:iMz) / Buff_I(iRho) * Weight
       if(DoCoupleVar_V(BField_))&
            State_VIB(Bx_:Bz_,i,iBlock) = Aux * State_VIB(Bx_:Bz_,i,iBlock) + &
            Buff_I(iBx:iBz) * Weight
    end do
  end subroutine SP_put_from_mh

  !===================================================================

  subroutine SP_set_grid

    ! Initialize 3D grid with NON-TREE structure
    call init_decomposition(&
         GridID_ = SP_,&
         CompID_ = SP_,&
         nDim    = nDim)

    ! Construct decomposition
    if(is_proc0(SP_))&
         call get_root_decomposition(&
         GridID_       = SP_,&
         iRootMapDim_D = (/1, nLon, nLat/),&
         XyzMin_D      = (/real(iParticleMin)-0.5, LonMin, LatMin/),&
         XyzMax_D      = (/real(iParticleMax)+0.5, LonMax, LatMax/),&
         nCells_D      = (/nParticle , 1, 1/),&
         PE_I          = iGridGlobal_IA(Proc_,:),&
         iBlock_I      = iGridGlobal_IA(Block_,:))
    call bcast_decomposition(SP_)

    ! Coordinate system is Heliographic Inertial Coordinate System (HGI)
    ! with length measured in solar radii
    call set_coord_system(&
         GridID_      = SP_, &
         TypeCoord    = TypeCoordSystem, &
         TypeGeometry = 'spherical', &
         NameVar      = NameVarCouple, &
         UnitX        = rSun)
  end subroutine SP_set_grid

  !===================================================================

  subroutine SP_get_solar_corona_boundary(RScOut)
    ! return the value of the solar corona boundary as set in SP component
    real, intent(out):: RScOut
    !-----------------------------------------------------------------
    RScOut = RSc
  end subroutine SP_get_solar_corona_boundary

  !===================================================================

  subroutine SP_put_r_min(RMinIn)
    ! save the lower boundary of the domain as set in other components;
    ! compute coordinates of the footprints of field lines
    real, intent(in):: RMinIn
    
    ! exisiting particle with lowest index along line
    real:: Coord1_D(nDim), Xyz1_D(nDim) 
    ! direction of the field at Xyz1_D
    real:: Dir1_D(nDim) 
    ! dot product Xyz1 and Dir1 and its sign
    real:: Dot, S
    ! variable to compute coords of the footprints
    real:: Alpha

    ! the footprint
    real:: Coord0_D(nDim), Xyz0_D(nDim) 

    integer:: iBlock    ! loop variable
    !-----------------------------------------------------------------
    RMin = RMinIn

    ! compute coordinates of footprints for each field lines
    do iBlock = 1, nBlock
       ! get the coordinates of lower particle
       Coord1_D = &
            State_VIB((/R_, Lon_, Lat_/), iGridLocal_IB(Begin_,iBlock), iBlock)
       call rlonlat_to_xyz(Coord1_D, Xyz1_D) 

       ! get the field direction at this location
       Dir1_D = &
            State_VIB((/Bx_, By_, Bz_/), iGridLocal_IB(Begin_,iBlock), iBlock)
       Dir1_D = Dir1_D / sqrt(sum(Dir1_D**2))
       
       ! their dot product and its sign
       Dot = sum(Dir1_D*Xyz1_D)
       S   = sign(1.0, Dot)

       ! Xyz0 is distance Alpha away from Xyz1:
       ! Xyz0 = Xyz1 + Alpha * Dir1 and R0 = RMin =>
       Alpha = S * sqrt(Dot**2 - Coord1_D(R_)**2 + RMin**2) - Dot
       Xyz0_D = Xyz1_D + Alpha * Dir1_D
       call xyz_to_rlonlat(Xyz0_D, CoordMin_DI(:,iBlock)) 
    end do
  end subroutine SP_put_r_min

  !===================================================================

  subroutine SP_interface_point_coords(iComp,&
       nDim, Xyz_D, nIndex, iIndex_I, IsInterfacePoint)
    ! interface points (request), which needed to be communicated
    ! to other components to perform field line extraction and
    ! perform further coupling with SP:
    ! the framework tries to determine Xyz_D of such points,
    ! SP changes them to the correct values
    integer, intent(in)  :: iComp
    integer,intent(in)   :: nDim
    real,   intent(inout):: Xyz_D(nDim)
    integer,intent(in)   :: nIndex
    integer,intent(inout):: iIndex_I(nIndex)
    logical,intent(out)  :: IsInterfacePoint
    !---------------------------------------------------------------
    integer:: iBlock, iRequest     ! loop variables
    integer:: iRequest_I(nRequestIndex)
    logical, save:: IsFirstCall = .true.
    character(len=100):: StringError
    character(len=*), parameter:: NameSub='SP_interface_point_coords'
    !----------------------------------------------------------------
    if(IsFirstCall)then
       IsFirstCall = .false.
       ! need to initalize request buffer;
       ! initial size = number of lines:
       allocate(iRequestSc_II(nRequestIndex,nBlock))
       allocate(iRequestIh_II(nRequestIndex,nBlock))
       
       ! fill the request buffer for SC
       do iBlock = 1, nBlock
          iRequestSc_II(:,iBlock) = (/iBlock, 1/)
       end do
       
       ! reset the sizes of requests
       nRequestSc = nBlock
       nRequestIh = 0

       IsSentRequestSc = .false.
       IsSentRequestIh = .false.
    end if

    ! choose the request buffer based iComp and reset its size counter
    select case(iComp)
    case(SC_)
       if(.not.IsSentRequestSc) IsSentRequestSc = .true.
       call check_request(iRequestSc_II, nRequestSc, IsInterfacePoint)
    case(IH_)
       if(.not.IsSentRequestIh) IsSentRequestIh = .true.
       call check_request(iRequestIh_II, nRequestIh, IsInterfacePoint)
    case default
       write(StringError,'(a,i2)') &
            ": isn't implemented for interface with component ", iComp
       call CON_stop(NameSub//StringError)
    end select
  contains
    subroutine check_request(iRequest_II, nRequest, IsInterfacePointLocal)
      ! check whether the point is requested and return its correct coordinates
      integer, pointer, intent(inout):: iRequest_II(:,:)
      integer,          intent(inout):: nRequest
      logical,          intent(out)  :: IsInterfacePointLocal
      !--------------------------------------------------------------
      do iRequest = 1, nRequest
         if(all(iRequest_II(:,iRequest)==iIndex_I((/nIndex,1/))))then
            IsInterfacePointLocal = .true.
            Xyz_D = State_VIB((/R_,Lon_,Lat_/),iIndex_I(1),iIndex_I(nIndex))
            RETURN
         end if
      end do
      IsInterfacePointLocal = .false.
    end subroutine check_request
  end subroutine SP_interface_point_coords

  !===================================================================
  subroutine SP_interface_point_coords_for_sc(&
       GridDescriptor, lGlobalTreeNode, nDim, Xyz_D, nIndex,iIndex_I,&
       IsInterfacePoint)
    use CON_grid_descriptor
    type(GridDescriptorType),intent(in)::GridDescriptor
    integer,intent(in)   :: lGlobalTreeNode
    integer,intent(in)   :: nDim
    real,   intent(inout):: Xyz_D(nDim)
    integer,intent(in)   :: nIndex
    integer,intent(inout):: iIndex_I(nIndex)
    logical,intent(out)  :: IsInterfacePoint
    !---------------------------------------------------------------
    call SP_interface_point_coords(SC_,&
         nDim,Xyz_D,nIndex,iIndex_I, IsInterfacePoint)
  end subroutine SP_interface_point_coords_for_sc
  !===================================================================
  subroutine SP_interface_point_coords_for_ih(&
       GridDescriptor, lGlobalTreeNode, nDim, Xyz_D, nIndex, iIndex_I,&
       IsInterfacePoint)
    use CON_grid_descriptor
    type(GridDescriptorType),intent(in)::GridDescriptor
    integer,intent(in)   :: lGlobalTreeNode
    integer,intent(in)   :: nDim
    real,   intent(inout):: Xyz_D(nDim)
    integer,intent(in)   :: nIndex
    integer,intent(inout):: iIndex_I(nIndex)
    logical,intent(out)  :: IsInterfacePoint
    !---------------------------------------------------------------
    call SP_interface_point_coords(IH_,&
         nDim,Xyz_D,nIndex,iIndex_I, IsInterfacePoint)
  end subroutine SP_interface_point_coords_for_ih

  !===================================================================

  subroutine SP_put_line(nPut, Coord_DI, iIndex_II)
    use ModMpi
    ! store particle coordinates extracted elsewhere
    !---------------------------------------------------------------
    integer, intent(in):: nPut
    real,    intent(in):: Coord_DI( nDim,   nPut)
    integer, intent(in):: iIndex_II(nDim+1, nPut)

    ! cartesian coordinates
    real:: Xyz_D(nDim)
    ! radius-lon-lat coordinates
    real:: Coord_D(nDim)
    ! loop variables
    integer:: iPut, iBlock, iNode
    ! indices of the particle
    integer:: iLine, iParticle
    integer:: iMin_A(nNode),iMax_A(nNode)
    integer:: iError
    ! which field lines are being extracted
    logical:: WasInSc, IsInSc, WasUndef
    logical, save:: IsFirstCall = .true.
    logical, save,allocatable:: DoneExtractSolarCorona_B(:)
    integer, parameter:: nVarReset  = 8
    integer, parameter:: &
         VarReset_I(nVarReset) = (/Rho_,Bx_,By_,Bz_,T_,Ux_,Uy_,Uz_/)
    character(len=*), parameter:: NameSub='SP_put_line'
    !----------------------------------------------------------------
    if(IsFirstCall)then
       ! initialize information of which 
       IsFirstCall = .false.
       allocate(DoneExtractSolarCorona_B(nBlock))
       DoneExtractSolarCorona_B = .false.
    end if

    ! if the request has been sent to other components -> reset its size
    if(IsSentRequestSc) nRequestSc = 0
    if(IsSentRequestIh) nRequestIh = 0

    ! store passed particles
    do iPut = 1, nPut
       iBlock = iIndex_II(4, iPut)
       iLine  = iNode_B(iBlock)
       iParticle = iIndex_II(1, iPut)

       if(iParticle < iParticleMin)&
            call CON_stop(NameSub//': particle index is below limit')
       if(iParticle > iParticleMax)&
            call CON_stop(NameSub//': particle index is above limit')
       iGridLocal_IB(Begin_,iBlock)=MIN(iGridLocal_IB(Begin_,iBlock),iParticle)
       iGridLocal_IB(End_,  iBlock)=MAX(iGridLocal_IB(End_,  iBlock),iParticle)
       if(iGridGlobal_IA(Proc_, iLine) /= iProc)&
            call CON_stop(NameSub//': Incorrect message pass')
       
       ! check if the particle has crossed the solar corona boundary
       WasUndef = State_VIB(R_, iParticle, iBlock) < 0.0
       WasInSc  = State_VIB(R_, iParticle, iBlock) < RSc
       IsInSc   = Coord_DI( R_, iPut)              < RSc
       if(.not.WasUndef .and. WasInSc .and. .not. IsInSc)then
          ! particle existed and crossed SC boundary
          call add_to_request(IH_,(/iBlock, iParticle/))
       elseif(.not.DoneExtractSolarCorona_B(iBlock) .and. &
            WasUndef .and. .not. IsInSc)then
          ! particle didn't exist, is the 1st beyond SC
          call add_to_request(IH_,(/iBlock, iParticle/))
          DoneExtractSolarCorona_B(iBlock) = .true.
       end if

       ! keep some variables as "old" state
       State_VIB((/RhoOld_,BOld_/),  iParticle, iBlock) = &
            State_VIB((/Rho_,B_/),   iParticle, iBlock)
       ! reset others 
       State_VIB(VarReset_I,         iParticle, iBlock) = 0.0
       ! put coordinates
       State_VIB((/R_, Lon_, Lat_/), iParticle, iBlock) = &
            Coord_DI(1:nDim, iPut)
    end do
  end subroutine SP_put_line

  !===================================================================

  subroutine SP_get_grid_descriptor_param(&
       iGridMin_D, iGridMax_D, Displacement_D)
    integer, intent(out):: iGridMin_D(nDim)
    integer, intent(out):: iGridMax_D(nDim)
    real,    intent(out):: Displacement_D(nDim)
    !-----------------------------------------
    iGridMin_D = (/iParticleMin, 1, 1/)
    iGridMax_D = (/iParticleMax, 1, 1/)
    Displacement_D = 0.0
  end subroutine SP_get_grid_descriptor_param

end module SP_wrapper
