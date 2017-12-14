!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!=============================================================!
module SP_wrapper

  use ModConst, ONLY: rSun, cProtonMass, energy_in
  use SP_ModMain, ONLY: &
       run, initialize, check, read_param, save_restart, &
       get_node_indexes, append_particles, &
       DoRestart, &
       iComm, iProc, nProc, &
       nDim, nLat, nLon, nBlock, nParticleMax, &
       RMin, RBufferMin, RBufferMax, RMax, LatMin, LatMax, LonMin, LonMax, &
       iGridGlobal_IA, iGridLocal_IB, State_VIB, Distribution_IIB,&
       iNode_B, TypeCoordSystem, &
       ParamLocal_IB, DataInputTime, Offset_,&
       Block_, Proc_, End_, Shock_, ShockOld_, XMin_, ZMin_, Length_,&
       LagrID_,X_,Y_,Z_, Rho_, Bx_,By_,Bz_,B_, Ux_,Uy_,Uz_, T_, RhoOld_,BOld_,&
       Wave1_, Wave2_
  use CON_comp_info
  use CON_router, ONLY: IndexPtrType, WeightPtrType
  use CON_coupler, ONLY: &
       set_coord_system, SP_, is_proc0, i_comm, i_proc0, &
       init_decomposition, get_root_decomposition, bcast_decomposition, &
       iVar_V, DoCoupleVar_V, &
       Density_, RhoCouple_, Pressure_, PCouple_, &
       Momentum_, RhoUxCouple_, RhoUzCouple_, &
       BField_, BxCouple_, BzCouple_, &
       Wave_, WaveFirstCouple_, WaveLastCouple_
  use ModMpi
  use CON_world, ONLY: is_proc0, is_proc, n_proc

  implicit none

  save

  private ! except

  public:: SP_set_param
  public:: SP_init_session
  public:: SP_run
  public:: SP_save_restart
  public:: SP_finalize

  ! coupling with MHD components
  public:: SP_put_from_mh
  public:: SP_interface_point_coords
  public:: SP_put_line
  public:: SP_adjust_lines
  public:: SP_get_bounds_comp
  public:: SP_n_particle
  public:: SP_check_ready_for_mh
  public:: SP_put_coupling_param
  ! variables requested via coupling: coordinates, 
  ! field line and particles indexes
  character(len=*), parameter:: NameVarCouple =&
       'rho p mx my mz bx by bz i01 i02 pe'
  integer :: Model_ = -1
  integer, parameter:: Lower_=0, Upper_=1
  real :: rInterfaceMin, rInterfaceMax
contains
  !\Interface routines to be called from super-structure only  
  subroutine SP_check_ready_for_mh(IsReady)
    use ModMpi
    logical, intent(out):: IsReady

    integer :: iError
    character(len=*), parameter:: NameSub='SP_check_ready_for_mh'
    !--------------
    ! when restarting, line data is available, i.e. ready to couple with mh;
    ! get value at SP root and broadcast to all SWMF processors
    if(is_proc0(SP_)) &
         IsReady = DoRestart
    call MPI_Bcast(IsReady, 1, MPI_LOGICAL, i_proc0(SP_), i_comm(), iError)
  end subroutine SP_check_ready_for_mh
  !===================================================================
  subroutine SP_get_bounds_comp(ThisModel_, RMinOut, RMaxOut)
    use ModMpi
    ! return the MHD boundaries as set in SP component
    integer, intent(in )  :: ThisModel_
    real,    intent(out)  :: RMinOut, RMaxOut
    integer :: iError
    real    :: rAux_I(2)
    character(len=*), parameter :: NameSub = 'SP_get_bounds_comp'
    !-----------------------------------------------------------------
    if(is_proc0(SP_))then
       select case(ThisModel_)
       case(Lower_)
          rAux_I(1) = RMin
          rAux_I(2) = RBufferMax
       case(Upper_)
          rAux_I(1) = RBufferMin
          rAux_I(2) = RMax
       case default
          call CON_stop('Incorrect model ID in '//NameSub)
       end select
    end if
    call MPI_Bcast(rAux_I(1), 2, MPI_REAL, i_proc0(SP_), i_comm(), iError)
       RMinOut = rAux_I(1); RMaxOut = rAux_I(2)
    
  end subroutine SP_get_bounds_comp 
  ! Above routines may be called from superstructure only.
  !/ 
  !========================================================================
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
    real:: TimeAux ! to satisfy intent of arguments in run()
    !--------------------------------------------------------------------
    TimeAux = TimeSimulation
    call run(TimeAux, TimeSimulation, .true.)
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
            NameVersion='MFLAMPA', &
            Version    =0.90)
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
    real:: TimeAux ! to satisfy intent of arguments in run()
    !--------------------------------------------------------------------
    TimeAux = TimeSimulation
    call run(TimeAux, TimeSimulation)
    call save_restart
  end subroutine SP_save_restart


  !===================================================================

  subroutine SP_put_from_mh(nPartial,iPutStart,Put,W,DoAdd,Buff_I,nVar)
  !subroutine SP_put_from_mh(iComp, nPartial,iPutStart,Put,W,DoAdd,Buff_I,nVar)
  !  integer, intent(in):: iComp
    integer,intent(in)::nPartial,iPutStart,nVar
    type(IndexPtrType),intent(in)::Put
    type(WeightPtrType),intent(in)::W
    logical,intent(in)::DoAdd
    real,dimension(nVar),intent(in)::Buff_I
    integer:: iRho, iP, iMx, iMz, iBx, iBz, iWave1, iWave2
    integer:: i, j, k, iBlock
    integer:: iPartial
    real:: Weight
    real:: R, Aux

    character(len=100):: StringError
    character(len=*), parameter:: NameSub='SP_put_from_mh'
    !------------------------------------------------------------
    !check consistency of DoCoupleVar_V
    if(.not. DoCoupleVar_V(Density_) .and. &
         (DoCoupleVar_V(Pressure_) .or. DoCoupleVar_V(Momentum_)))&
         call CON_Stop(NameSub//': pressure or momentum is coupled,'//&
         ' but density is not')

    ! indices of variables in the buffer
    iRho  = iVar_V(RhoCouple_)
    iP    = iVar_V(PCouple_)
    iMx   = iVar_V(RhoUxCouple_)
    iMz   = iVar_V(RhoUzCouple_)
    iBx   = iVar_V(BxCouple_)
    iBz   = iVar_V(BzCouple_)   
    iWave1= iVar_V(WaveFirstCouple_)
    iWave2= iVar_V(WaveLastCouple_)
    Aux = 0
    if(DoAdd)Aux = 1.0
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
          select case(Model_)
          case(Lower_)
          !select case(iComp)
          !case(SC_)
             !if(Model_/=Lower_)call CON_stop('Incorrect Model_for SC_')
             Weight = Weight * (0.50 - 0.50 * &
                  tanh(2*(2*R-RBufferMax-RBufferMin)/(RBufferMax-RBufferMin)))
          case(Upper_)
          !case(IH_)
             !if(Model_/=Upper_)call CON_stop('Incorrect Model_for IH_')
             Aux = 1.0   
             Weight = Weight * (0.50 + 0.50 * &
                  tanh(2*(2*R-RBufferMax-RBufferMin)/(RBufferMax-RBufferMin)))
          case default
             write(StringError,'(a,i2)') &
                  ": isn't implemented for interface with component ", Model_
             call CON_stop(NameSub//StringError)
          end select
       end if
       ! put the data
       ! NOTE: State_VIB must be reset to zero before putting coupled data
       if(DoCoupleVar_V(Density_))&
            State_VIB(Rho_,i,iBlock) = Aux*State_VIB(Rho_,i,iBlock) + &
            Buff_I(iRho)/cProtonMass * Weight
       if(DoCoupleVar_V(Pressure_))&
            State_VIB(T_,i,iBlock) = Aux*State_VIB(T_,i,iBlock) + &
            Buff_I(iP)/Buff_I(iRho)*cProtonMass/energy_in('kev') * Weight
       if(DoCoupleVar_V(Momentum_))&
            State_VIB(Ux_:Uz_,i,iBlock) = Aux*State_VIB(Ux_:Uz_,i,iBlock) + &
            Buff_I(iMx:iMz) / Buff_I(iRho) * Weight
       if(DoCoupleVar_V(BField_))&
            State_VIB(Bx_:Bz_,i,iBlock) = Aux*State_VIB(Bx_:Bz_,i,iBlock) + &
            Buff_I(iBx:iBz) * Weight
       if(DoCoupleVar_V(Wave_))&
            State_VIB(Wave1_:Wave2_,i,iBlock) = &
            Aux*State_VIB(Wave1_:Wave2_,i,iBlock) + &
            Buff_I(iWave1:iWave2)*Weight
    end do
  end subroutine SP_put_from_mh
  !============================
  subroutine SP_set_grid
    logical, save:: IsInitialized = .false.
    !------------------------------------------------------------
    if(IsInitialized)&
         RETURN
    IsInitialized = .true.

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
         CoordMin_D    = (/0.5, LonMin, LatMin/),&
         CoordMax_D    = (/real(nParticleMax)+0.5, LonMax, LatMax/),&
         nCells_D      = (/nParticleMax, 1, 1/),&
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
  !================================
  subroutine SP_put_coupling_param(iModelIn, rMinIn, rMaxIn, TimeIn)
    integer, intent(in) :: iModelIn
    real,    intent(in) :: rMinIn, rMaxIn
    real,     intent(in)::TimeIn
    !-----------------
    call SP_put_interface_bounds(iModelIn, rMinIn, rMaxIn)
    call SP_put_input_time(TimeIn)
  end subroutine SP_put_coupling_param
  subroutine SP_put_interface_bounds(iModelIn, rMinIn, rMaxIn)
    integer, intent(in) :: iModelIn
    real,    intent(in) :: rMinIn, rMaxIn
    !-----------------
    rInterfaceMin = rMinIn; rInterfaceMax = rMaxIn; Model_ = iModelIn
  end subroutine SP_put_interface_bounds
  !=========================================================
  subroutine SP_put_input_time(TimeIn)
    real,     intent(in)::TimeIn
    if(DataInputTime >= TimeIn)RETURN
    select case(Model_)
    case(Lower_)
       call copy_old_state
    case(Upper_)
       call CON_stop(&
            "Time in coupling to IH differs from that to SC")
    case default
       call CON_stop("Wrong model name")
    end select
    DataInputTime = TimeIn
  contains
    subroutine copy_old_state
      use SP_ModGrid, ONLY: nVarRead
      ! copy current state to old state for all field lines
      integer:: iEnd, iBlock
      !-----------------------------------------------------------------------
      do iBlock = 1, nBlock
         iEnd   = iGridLocal_IB(End_,  iBlock)
         State_VIB((/RhoOld_,BOld_/), 1:iEnd, iBlock) = &
              State_VIB((/Rho_,B_/),  1:iEnd, iBlock)
         ! reset other variables
         State_VIB(1:nVarRead,1:iEnd, iBlock) = 0.0
      end do
    end subroutine copy_old_state

  end subroutine SP_put_input_time
  !===================================================================
  subroutine SP_interface_point_coords(nDim, Xyz_D, &
       nIndex, iIndex_I, IsInterfacePoint)
    ! interface points (request), which needed to be communicated
    ! to other components to perform field line extraction and
    ! perform further coupling with SP:
    ! the framework tries to determine Xyz_D of such points,
    ! SP changes them to the correct values
    integer,intent(in)   :: nDim
    real,   intent(inout):: Xyz_D(nDim)
    integer,intent(in)   :: nIndex
    integer,intent(inout):: iIndex_I(nIndex)
    logical,intent(out)  :: IsInterfacePoint
    integer:: iParticle, iBlock
    real:: R2
    character(len=*), parameter:: NameSub='SP_interface_point_coords'
    !----------------------------------------------------------------
    iParticle = iIndex_I(1)
    iBlock    = iIndex_I(4)
    ! Check whether the particle is within interface bounds
    R2 = sum(State_VIB(X_:Z_,iParticle,iBlock)**2)
    IsInterfacePoint = &
         R2 >= rInterfaceMin**2 .and. R2 < rInterfaceMax**2
    ! Fix coordinates to be used in mapping
    if(IsInterfacePoint)then
       Xyz_D = State_VIB(X_:Z_, iParticle, iBlock)
    end if
  end subroutine SP_interface_point_coords
  !============================
  subroutine SP_put_line(nPartial, iPutStart, Put,&
       Weight, DoAdd, Coord_D, nVar)
    integer, intent(in) :: nPartial, iPutStart, nVar
    type(IndexPtrType), intent(in) :: Put
    type(WeightPtrType),intent(in) :: Weight
    logical,            intent(in) :: DoAdd
    real,               intent(in) :: Coord_D(nVar) !nVar=nDim

    ! indices of the particle
    integer:: iBlock, iParticle

    character(len=*), parameter:: NameSub='SP_put_line'
    !----------------------------------------------------------------
    ! store passed particles
    iBlock    = Put%iCB_II(4,iPutStart)
    iParticle = Put%iCB_II(1,iPutStart) + iGridLocal_IB(Offset_,iBlock)
    ! put coordinates
    State_VIB(X_:Z_,  iParticle, iBlock) = Coord_D(1:nDim)
    iGridLocal_IB(End_,  iBlock)=MAX(iGridLocal_IB(End_,  iBlock),iParticle)
  end subroutine SP_put_line

  !===================================================================
  !\
  ! Called from coupler after the updated grid point location are 
  ! received from the other component (SC, IH). Determines whether some
  ! grid points should be added/deleted
  !/
  subroutine SP_adjust_lines(DoInit, DoAdjustStart, DoAdjustEnd)
    !\
    ! If DoAdjustStart, the points in the starting portion of the line are
    ! processed, if DoAdjustEnd - the same for the end points
    !/
    logical, intent(in) :: DoInit, DoAdjustStart, DoAdjustEnd
    ! once new geometry of lines has been put, account for some particles
    ! exiting the domain (can happen both at the beginning and the end)
    integer:: iParticle, iBlock, iBegin,  iEnd, iOffset ! loop variables
    logical:: IsMissingCurr, IsMissingPrev
    real   :: R2
    
    character(len=*), parameter:: NameSub = "SP_adjust_lines"
    character(len=100):: StringError
    !--------------------------------------------------------------------
    if(DoInit)then
       if(DoAdjustStart)then
          do iBlock = 1, nBlock
             call SP_set_line_foot_b(iBlock)
          end do
       end if
       RETURN
    end if
    do iBlock = 1, nBlock
       !\
       ! Called after the grid points are received from the 
       ! component, nullify offset
       !/
       if(DoAdjustStart)iGridLocal_IB(Offset_,iBlock) = 0
       iBegin = 1
       iEnd   = iGridLocal_IB(End_,  iBlock)
       IsMissingPrev = all(State_VIB(X_:Z_,1,iBlock)==0.0)
       R2 = sum(State_VIB(X_:Z_,1,iBlock)**2)
       do iParticle = 2, iEnd
          IsMissingCurr = all(State_VIB(X_:Z_,iParticle,iBlock)==0.0)

          if(IsMissingCurr .and. R2 > RBufferMin**2)then
             if(DoAdjustEnd)&
                  iGridLocal_IB(End_,  iBlock) = iParticle - 1
             EXIT
          end if

          if(.not.IsMissingCurr)then
             R2 = sum(State_VIB(X_:Z_,iParticle,iBlock)**2)
             if(IsMissingPrev)then
                iBegin = iParticle
             end if
          end if
          IsMissingPrev = IsMissingCurr
       end do
       if(DoAdjustStart.and.iBegin /= 1)then
          !\
          ! Offset particle arrays
          !/
          iEnd   = iGridLocal_IB(End_,   iBlock) 
          iOffset = 1 - iBegin
          iGridLocal_IB(End_,   iBlock) = iEnd  + iOffset
          iGridLocal_IB(Shock_, iBlock) = iGridLocal_IB(Shock_, iBlock) +&
               iOffset
          iGridLocal_IB(ShockOld_, iBlock) = iGridLocal_IB(ShockOld_, iBlock)+&
               iOffset
          iGridLocal_IB(Offset_, iBlock)   = iGridLocal_IB(Offset_, iBlock) +&
               iOffset
          State_VIB(:, 1:iGridLocal_IB(End_,iBlock), iBlock) = &
               State_VIB(:, iBegin:iEnd, iBlock)
          Distribution_IIB(:,1:iGridLocal_IB(End_,iBlock), iBlock) = &
               Distribution_IIB(:,iBegin:iEnd, iBlock)
          ! need to recalculate footpoints
          call SP_set_line_foot_b(iBlock)
       end if
    end do
    ! may need to add particles to the beginning of lines
    if(DoAdjustStart) call append_particles
    !\
    ! Called after the grid points are received from the 
    ! component, nullify offset
    !/
    if(DoAdjustEnd)iGridLocal_IB(Offset_,1:nBlock) = 0
  contains
    subroutine SP_set_line_foot_b(iBlock)
      integer, intent(in) :: iBlock

      ! existing particle with lowest index along line
      real:: Xyz1_D(nDim)
      ! direction of the field at Xyz1_D and segment vectors between particles
      real, dimension(nDim):: Dir0_D, Dist1_D, Dist2_D
      ! dot product Xyz1 and Dir1 and its sign
      real:: Dot, S
      ! distances between particles
      real:: Dist1, Dist2
      ! variable to compute coords of the footprints
      real:: Alpha
      !---------------
      ! get the coordinates of lower particle
      Xyz1_D = State_VIB(X_:Z_, 1, iBlock)

      ! generally, field direction isn't known
      ! approximate it using directions of first 2 segments of the line
      Dist1_D = &
           State_VIB((/X_, Y_, Z_/), 1, iBlock) - &
           State_VIB((/X_, Y_, Z_/), 2, iBlock)
      Dist1 = sqrt(sum(Dist1_D**2))
      Dist2_D = &
           State_VIB((/X_, Y_, Z_/), 2, iBlock) - &
           State_VIB((/X_, Y_, Z_/), 3, iBlock)
      Dist2 = sqrt(sum(Dist2_D**2))
      Dir0_D = ((2*Dist1+Dist2)*Dist1_D - Dist1*Dist2_D) / (Dist1 + Dist2)

      Dir0_D = Dir0_D / sqrt(sum(Dir0_D**2))

      ! dot product and sign: used in computation below
      Dot = sum(Dir0_D*Xyz1_D)
      S   = sign(1.0, Dot)

      ! Xyz0, the footprint, is distance Alpha away from Xyz1:
      ! Xyz0 = Xyz1 + Alpha * Dir0 and R0 = RMin =>
      Alpha = S * sqrt(Dot**2 - sum(Xyz1_D**2) + RMin**2) - Dot

      ! store new footpoint of the line
      ParamLocal_IB(XMin_:ZMin_,iBlock) = Xyz1_D + Alpha * Dir0_D
      ! length is used to decide when need to append new particles:
      ! use distance between first two particles on the line
      ParamLocal_IB(Length_,    iBlock) = Dist1
    end subroutine SP_set_line_foot_b
  end subroutine SP_adjust_lines
  !=============================
 
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
