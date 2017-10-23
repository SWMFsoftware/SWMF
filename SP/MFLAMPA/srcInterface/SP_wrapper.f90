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
       RMin, RBufferMin, RBufferMax, RMax, LatMin, LatMax, LonMin, LonMax, &
       iGridGlobal_IA, iGridLocal_IB, State_VIB, Distribution_IIB,&
       iNode_B, TypeCoordSystem, &
       CoordMin_DI, DataInputTime, &
       Block_, Proc_, Begin_, End_, &
       X_, Y_, Z_, Rho_, Bx_,By_,Bz_,B_, Ux_,Uy_,Uz_, T_, RhoOld_, BOld_
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
  public:: SP_put_from_sc
  public:: SP_put_from_ih
  public:: SP_interface_point_coords_for_ih
  public:: SP_interface_point_coords_for_ih_extract
  public:: SP_interface_point_coords_for_sc
  public:: SP_put_line
  public:: SP_adjust_lines
  public:: SP_get_grid_descriptor_param
  public:: SP_get_domain_boundary
  public:: SP_put_r_min
  public:: SP_n_particle
  public:: SP_copy_old_state

  ! variables requested via coupling: coordinates, 
  ! field line and particles indexes
  character(len=*), parameter:: NameVarCouple = 'rho p mx my mz bx by bz'

contains
  !========================
  integer function SP_n_particle(iBlockLocal)
    integer, intent(in) :: iBlockLocal
    SP_n_particle = iGridLocal_IB(End_,  iBlockLocal)
  end function SP_n_particle

  !========================================================================

  subroutine SP_run(TimeSimulation,TimeSimulationLimit)
    real,intent(inout)::TimeSimulation
    real,intent(in)::TimeSimulationLimit
    !--------------------------------------------------------------------------
    call run(TimeSimulation,TimeSimulationLimit)
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
    DataInputTime = TimeIn
  end subroutine SP_put_input_time

  !===================================================================

  subroutine SP_put_from_mh(iComp, nPartial,iPutStart,Put,W,DoAdd,Buff_I,nVar)
    integer, intent(in):: iComp
    integer,intent(in)::nPartial,iPutStart,nVar
    type(IndexPtrType),intent(in)::Put
    type(WeightPtrType),intent(in)::W
    logical,intent(in)::DoAdd
    real,dimension(nVar),intent(in)::Buff_I
    integer:: iRho, iP, iMx, iMz, iBx, iBz
    integer:: i, j, k, iBlock
    integer:: iPartial
    real:: Weight
    real:: R
    real:: Aux

    real, external:: energy_in

    character(len=100):: StringError
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
       if(is_in_buffer(State_VIB(X_:Z_,i,iBlock)))then
          R = sqrt(sum(State_VIB(X_:Z_,i,iBlock)**2))
          select case(iComp)
          case(SC_)
             Weight = Weight * (0.50 - 0.50 * &
                  tanh(2*(2*R-RBufferMax-RBufferMin)/(RBufferMax-RBufferMin)))
          case(IH_)
             Aux    = 1.0
             Weight = Weight * (0.50 + 0.50 * &
                  tanh(2*(2*R-RBufferMax-RBufferMin)/(RBufferMax-RBufferMin)))
          case default
             write(StringError,'(a,i2)') &
                  ": isn't implemented for interface with component ", iComp
             call CON_stop(NameSub//StringError)
          end select
       end if
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
  subroutine SP_put_from_sc(nPartial,iPutStart,Put,W,DoAdd,Buff_I,nVar)
    integer,intent(in)::nPartial,iPutStart,nVar
    type(IndexPtrType),intent(in)::Put
    type(WeightPtrType),intent(in)::W
    logical,intent(in)::DoAdd
    real,dimension(nVar),intent(in)::Buff_I
    !------------------------------------------------------------
    call SP_put_from_mh(SC_,nPartial,iPutStart,Put,W,DoAdd,Buff_I,nVar)
  end subroutine SP_put_from_sc
  !===================================================================
  subroutine SP_put_from_ih(nPartial,iPutStart,Put,W,DoAdd,Buff_I,nVar)
    integer,intent(in)::nPartial,iPutStart,nVar
    type(IndexPtrType),intent(in)::Put
    type(WeightPtrType),intent(in)::W
    logical,intent(in)::DoAdd
    real,dimension(nVar),intent(in)::Buff_I
    !------------------------------------------------------------
    call SP_put_from_mh(IH_,nPartial,iPutStart,Put,W,DoAdd,Buff_I,nVar)
  end subroutine SP_put_from_ih
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
         TypeGeometry = 'cartesian', &
         NameVar      = NameVarCouple, &
         UnitX        = rSun)
  end subroutine SP_set_grid

  !===================================================================

  subroutine SP_get_domain_boundary(RScOut, RIhOut)
    ! return the value of the solar corona boundary as set in SP component
    real, intent(out):: RScOut, RIhOut
    !-----------------------------------------------------------------
    RScOut = RBufferMax
    RIhOut = RMax
  end subroutine SP_get_domain_boundary

  !===================================================================

  subroutine SP_put_r_min(RMinIn)
    ! save the lower boundary of the domain as set in other components;
    ! compute coordinates of the footprints of field lines
    real, intent(in):: RMinIn
    
    ! existing particle with lowest index along line
    real:: Xyz1_D(nDim)
    ! direction of the field at Xyz1_D
    real:: Dir1_D(nDim) 
    ! dot product Xyz1 and Dir1 and its sign
    real:: Dot, S
    ! variable to compute coords of the footprints
    real:: Alpha

    integer:: iBlock    ! loop variable
    !-----------------------------------------------------------------
    RMin = RMinIn

    ! compute coordinates of footprints for each field lines
    do iBlock = 1, nBlock
       ! get the coordinates of lower particle
       Xyz1_D = State_VIB(X_:Z_, iGridLocal_IB(Begin_,iBlock), iBlock)

       ! get the field direction at this location
       Dir1_D = &
            State_VIB((/Bx_, By_, Bz_/), iGridLocal_IB(Begin_,iBlock), iBlock)
       Dir1_D = Dir1_D / sqrt(sum(Dir1_D**2))
       
       ! their dot product and its sign
       Dot = sum(Dir1_D*Xyz1_D)
       S   = sign(1.0, Dot)

       ! Xyz0, the footprint, is distance Alpha away from Xyz1:
       ! Xyz0 = Xyz1 + Alpha * Dir1 and R0 = RMin =>
       Alpha = S * sqrt(Dot**2 - sum(Xyz1_D**2) + RMin**2) - Dot
       CoordMin_DI(:,iBlock) = Xyz1_D + Alpha * Dir1_D
    end do
  end subroutine SP_put_r_min

  !===================================================================

  subroutine SP_interface_point_coords(iComp, SendBuffer, &
       nDim, Xyz_D, nIndex, iIndex_I, IsInterfacePoint)
    ! interface points (request), which needed to be communicated
    ! to other components to perform field line extraction and
    ! perform further coupling with SP:
    ! the framework tries to determine Xyz_D of such points,
    ! SP changes them to the correct values
    integer,intent(in)   :: iComp
    logical,intent(in)   :: SendBuffer
    integer,intent(in)   :: nDim
    real,   intent(inout):: Xyz_D(nDim)
    integer,intent(in)   :: nIndex
    integer,intent(inout):: iIndex_I(nIndex)
    logical,intent(out)  :: IsInterfacePoint
    !---------------------------------------------------------------
    integer:: iParticle, iBlock
    logical:: IsSc
    character(len=100):: StringError
    character(len=*), parameter:: NameSub='SP_interface_point_coords'
    !----------------------------------------------------------------
    ! choose the request buffer based iComp and reset its size counter
    select case(iComp)
    case(SC_)
       IsSc = .true.
    case(IH_)
       IsSc = .false.
    case default
       write(StringError,'(a,i2)') &
            ": isn't implemented for interface with component ", iComp
       call CON_stop(NameSub//StringError)
    end select

    iParticle = iIndex_I(1)
    iBlock    = iIndex_I(4)
    ! first, check whether the particle is within bounds
    IsInterfacePoint = &
         iParticle >= iGridLocal_IB(Begin_,iBlock) .and. &
         iParticle <= iGridLocal_IB(End_,iBlock)
    ! second, check whether the particle is within the appropriate domain
    if(IsInterfacePoint)&
         IsInterfacePoint = &
         is_in_buffer(State_VIB(X_:Z_, iParticle, iBlock)) .and. SendBuffer .or.&
         .not.(IsSc.eqv.sum(State_VIB(X_:Z_,iParticle,iBlock)**2)>RBufferMax**2)
    ! lastly, fix coordinates
    if(IsInterfacePoint)&
         Xyz_D = State_VIB(X_:Z_, iParticle, iBlock)
  end subroutine SP_interface_point_coords

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
    !---------------------------------------------------------------
    call SP_interface_point_coords(SC_, .true., &
         nDim,Xyz_D,nIndex,iIndex_I, IsInterfacePoint)
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
    !---------------------------------------------------------------
    call SP_interface_point_coords(IH_, .true., &
         nDim,Xyz_D,nIndex,iIndex_I, IsInterfacePoint)
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
    !---------------------------------------------------------------
    call SP_interface_point_coords(IH_, .false., &
         nDim,Xyz_D,nIndex,iIndex_I, IsInterfacePoint)
  end subroutine SP_interface_point_coords_for_ih_extract

  !===================================================================

  subroutine SP_put_line(iComp, nPut, Coord_DI, iIndex_II)
    ! store particle coordinates extracted elsewhere
    !---------------------------------------------------------------
    integer, intent(in):: iComp
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
    logical:: PutInBuffer

    character(len=*), parameter:: NameSub='SP_put_line'
    !----------------------------------------------------------------
    PutInBuffer = iComp == SC_
    ! store passed particles
    do iPut = 1, nPut
       iBlock = iIndex_II(4, iPut)
       iLine  = iNode_B(iBlock)
       iParticle = iIndex_II(1, iPut)

       if(iParticle < iParticleMin)&
            call CON_stop(NameSub//': particle index is below limit')
       if(iParticle > iParticleMax)&
            call CON_stop(NameSub//': particle index is above limit')
       if(iGridGlobal_IA(Proc_, iLine) /= iProc)&
            call CON_stop(NameSub//': Incorrect message pass')

       if(.not.PutInBuffer .and. is_in_buffer(Coord_DI(1:nDim, iPut)) )&
            CYCLE

       ! put coordinates
       State_VIB(X_:Z_, iParticle, iBlock) = Coord_DI(1:nDim, iPut)
       iGridLocal_IB(Begin_,iBlock)=MIN(iGridLocal_IB(Begin_,iBlock),iParticle)
       iGridLocal_IB(End_,  iBlock)=MAX(iGridLocal_IB(End_,  iBlock),iParticle)
    end do
  end subroutine SP_put_line

  !===================================================================

  subroutine SP_adjust_lines
    ! once new geometry of lines has been put, account for some particles
    ! exiting the domain (can happen both at the beginning and the end)
    integer:: iParticle, iBlock, iBegin, iEnd, iEndNew ! loop variables
    logical:: IsMissingCurr, IsMissingPrev
    real   :: R2
    !--------------------------------------------------------------------
    do iBlock = 1, nBlock
       iBegin = iGridLocal_IB(Begin_,iBlock)
       iEnd   = iGridLocal_IB(End_,  iBlock)
       IsMissingPrev = all(State_VIB(X_:Z_,iBegin,iBlock)==0.0)
       R2 = sum(State_VIB(X_:Z_,iBegin,iBlock)**2)
       do iParticle = iBegin + 1, iEnd
          IsMissingCurr = all(State_VIB(X_:Z_,iParticle,iBlock)==0.0)
          if(IsMissingCurr .and. R2 > RBufferMin**2)then
             iGridLocal_IB(End_,  iBlock) = iParticle - 1
             EXIT
          end if
          
          if(.not.IsMissingCurr)then
             R2 = sum(State_VIB(X_:Z_,iParticle,iBlock)**2)
             if(IsMissingPrev)then
                iGridLocal_IB(Begin_,  iBlock) = iParticle
             end if
          end if
          IsMissingPrev = IsMissingCurr
       end do
       ! adjust data storage
       iEndNew= iGridLocal_IB(End_, iBlock) - &
            (iGridLocal_IB(Begin_, iBlock) - 1)
       if(  iGridLocal_IB(Begin_,iBlock) /= 1 .or.&
            iGridLocal_IB(End_,  iBlock) /= iEndNew)then
          iBegin = iGridLocal_IB(Begin_, iBlock)
          iEnd   = iGridLocal_IB(End_,   iBlock)
          iGridLocal_IB(Begin_, iBlock) = 1
          iGridLocal_IB(End_,   iBlock) = iEndNew
          State_VIB(:, 1:iEndNew, iBlock) = &
               State_VIB(:, iBegin:iEnd, iBlock)
          Distribution_IIB(:,1:iEndNew, iBlock) = &
               Distribution_IIB(:,iBegin:iEnd, iBlock)
       end if
    end do
  end subroutine SP_adjust_lines

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
  !========================================================================

  subroutine SP_copy_old_state
    ! copy current state to old state for all field lines
    integer:: iBegin, iEnd, iBlock
    integer, parameter:: nVarReset  = 11
    integer, parameter:: &
         VarReset_I(nVarReset) = (/X_,Y_,Z_,Rho_,Bx_,By_,Bz_,T_,Ux_,Uy_,Uz_/)
    !--------------------------------------------------------------------------
    do iBlock = 1, nBlock
       iBegin = iGridLocal_IB(Begin_,iBlock)
       iEnd   = iGridLocal_IB(End_,  iBlock)
       State_VIB((/RhoOld_,BOld_/), iBegin:iEnd, iBlock) = &
            State_VIB((/Rho_,B_/),  iBegin:iEnd, iBlock)
       ! reset other variables
       State_VIB(VarReset_I,iBegin:iEnd, iBlock) = 0.0
    end do
  end subroutine SP_copy_old_state

  !========================================================================

  function is_in_buffer(Xyz_D) Result(IsInBuffer)
    real,   intent(in) :: Xyz_D(nDim)
    logical:: IsInBuffer
    real:: R2
    !---------------------------------------------
    R2 = sum(Xyz_D**2)
    IsInBuffer = R2 >= RBufferMin**2 .and. R2 < RBufferMax**2
  end function is_in_buffer

end module SP_wrapper
