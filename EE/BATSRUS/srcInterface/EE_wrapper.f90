!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module EE_wrapper

  ! Wrapper for the BATSRUS Eruptive Event (EE) component

  use EE_ModBatsrusMethods, ONLY: &
       BATS_init_session, BATS_setup, BATS_advance, BATS_save_files, &
       BATS_finalize
  use ModUtilities, ONLY: CON_set_do_test, CON_stop

  implicit none

  private ! except

  public:: EE_set_param
  public:: EE_init_session
  public:: EE_run
  public:: EE_save_restart
  public:: EE_finalize

  ! Point coupler interface
  public:: EE_get_grid_info
  public:: EE_find_points

  ! EE-SC coupling
  public:: EE_get_for_SC
  public:: EE_put_from_sc

  ! Pointer coupling
  public:: EE_use_pointer

contains
  !============================================================================
  subroutine EE_set_param(CompInfo, TypeAction)

    use CON_comp_info
    use CON_physics, ONLY: get_time
    use EE_BATL_lib, ONLY: iProc, nProc, iComm
    use EE_ModIO, ONLY: iUnitOut, StringPrefix, STDOUT_, NamePlotDir
    use EE_ModSetParameters, ONLY: set_parameters
    use EE_ModRestartFile, ONLY: NameRestartInDir, NameRestartOutDir
    use EE_ModMain, ONLY : NameThisComp, &
         IsTimeAccurate, tSimulation, StartTime, iStartTime_I
    use ModTimeConvert, ONLY: time_real_to_int

    type(CompInfoType), intent(inout):: CompInfo   ! Information for this comp.
    character (len=*), intent(in)    :: TypeAction ! What to do

    logical :: DoTest, DoTestMe

    character(len=*), parameter:: NameSub = 'EE_set_param'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTest)write(*,*)NameSub,' called with TypeAction, iProc=', &
         TypeAction,iProc

    select case(TypeAction)
    case('VERSION')
       call put(CompInfo,&
            Use        =.true.,                        &
            NameVersion='EE_BATSRUS (Univ. of Michigan)')
    case('MPI')
       call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc,&
            Name=NameThisComp)

       NamePlotDir(1:2)       = NameThisComp
       NameRestartInDir(1:2)  = NameThisComp
       NameRestartOutDir(1:2) = NameThisComp
    case('READ')
       call set_parameters('READ')
    case('CHECK')
       call get_time( &
            DoTimeAccurateOut = IsTimeAccurate, &
            tSimulationOut=tSimulation, &
            tStartOut         = StartTime)
       call time_real_to_int(StartTime, iStartTime_I)

       call set_parameters('CHECK')
    case('STDOUT')
       iUnitOut=STDOUT_
       if(iProc==0)then
          StringPrefix = NameThisComp//':'
       else
          write(StringPrefix,'(a,i4.4,a)')NameThisComp,iProc,':'
       end if
    case('FILEOUT')
       call get(CompInfo, iUnitOut=iUnitOut)
       StringPrefix=''
    case('GRID')
       call EE_set_grid
    case default
       call CON_stop(NameSub//' SWMF_ERROR: invalid TypeAction='//TypeAction)
    end select

  end subroutine EE_set_param
  !============================================================================

  subroutine EE_set_grid

    use CON_coupler
    use CON_comp_param, ONLY: EE_
    use EE_domain_decomposition
    use EE_ModMain, ONLY: TypeCoordSystem, nVar, NameVarCouple
    use EE_ModAdvance, ONLY: State_VGB
    use EE_ModPhysics, ONLY: No2Si_V, UnitX_
    use EE_ModGeometry, ONLY: TypeGeometry, RadiusMin, RadiusMax
    use EE_BATL_lib, ONLY: CoordMin_D, CoordMax_D

    character(len=*), parameter:: NameSub = 'EE_set_grid'
    !--------------------------------------------------------------------------

    if(done_dd_init(EE_))RETURN

    call init_decomposition( &
         GridID_ = EE_, &
         CompID_ = EE_, &
         nDim = 3,      &
         IsTreeDD = .true.)

    call set_coord_system( &
         GridID_ = EE_, &
         TypeCoord = TypeCoordSystem, &
         UnitX = No2Si_V(UnitX_), &
         nVar = nVar, &
         NameVar = NameVarCouple,   &
         TypeGeometry=TypeGeometry, &
         Coord1_I = [ RadiusMin, RadiusMax ], &
         Coord2_I = [ CoordMin_D(2), CoordMax_D(2) ], &
         Coord3_I = [ CoordMin_D(3), CoordMax_D(3) ]  )

    if(is_proc(EE_)) Grid_C(EE_)%State_VGB => State_VGB

    if(is_proc(EE_))then
       ! set the local grid
       call init_decomposition(&
            Domain = MH_Domain, &
            CompID_ = EE_, &
            nDim = 3, &
            IsTreeDD = .true.)

       ! Get the octree root array
       call MH_get_root_decomposition(MH_Domain)

       ! Get the whole octree after the initial refinement
       call MH_update_local_decomposition(MH_Domain)
       MH_Domain%IsLocal=.true.
    end if

    ! Repeat the initialization at the global grid level:
    ! Octree root array:
    if(is_proc0(EE_))call MH_get_root_decomposition(EE_)

    ! Broadcast root array:
    call bcast_decomposition(EE_)

    ! Synchronize global and local grids:
    call synchronize_refinement( &
         GridID_     = EE_, &
         LocalDomain = MH_Domain)

  end subroutine EE_set_grid
  !============================================================================

  subroutine EE_init_session(iSession, TimeSimulation)

    integer, intent(in) :: iSession         ! session number (starting from 1)
    real,    intent(in) :: TimeSimulation   ! seconds from start time

    logical :: IsUninitialized = .true.
    logical :: DoTest, DoTestMe

    character(len=*), parameter:: NameSub = 'EE_init_session'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(IsUninitialized)then
       call BATS_setup
       IsUninitialized = .false.
    end if
    call BATS_init_session

    if(DoTest)write(*,*)NameSub,' finished for session ',iSession

  end subroutine EE_init_session
  !============================================================================

  subroutine EE_finalize(TimeSimulation)

    use EE_ModMain, ONLY: IsTimeLoop

    real, intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter:: NameSub = 'EE_finalize'
    !--------------------------------------------------------------------------
    ! We are not advancing in time any longer
    IsTimeLoop = .false.

    call BATS_save_files('FINAL')

    call BATS_finalize

  end subroutine EE_finalize
  !============================================================================

  subroutine EE_save_restart(TimeSimulation)

    real, intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter:: NameSub = 'EE_save_restart'
    !--------------------------------------------------------------------------
    call BATS_save_files('RESTART')

  end subroutine EE_save_restart
  !============================================================================

  subroutine EE_run(TimeSimulation, TimeSimulationLimit)

    use EE_BATL_lib, ONLY: iProc
    use EE_ModMain,  ONLY: tSimulation

    real, intent(inout):: TimeSimulation   ! current time of component
    real, intent(in):: TimeSimulationLimit ! simulation time not to be exceeded

    logical :: DoTest, DoTestMe

    character(len=*), parameter:: NameSub = 'EE_run'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTest)write(*,*)NameSub,' called with tSim, tSimLimit, iProc=',&
         TimeSimulation, TimeSimulationLimit, iProc

    if(abs(tSimulation-TimeSimulation)>0.0001) then
       write(*,*)NameSub, &
            ' EE time=',tSimulation,' SWMF time=',TimeSimulation
       call CON_stop(NameSub//': EE and SWMF simulation times differ')
    end if

    call BATS_advance(TimeSimulationLimit)

    ! Return time after the time step
    TimeSimulation = tSimulation

  end subroutine EE_run
  !============================================================================

  subroutine EE_find_points(nDimIn, nPoint, Xyz_DI, iProc_I)

    use EE_BATL_lib,   ONLY: MaxDim, find_grid_block
    use EE_ModPhysics, ONLY: Si2No_V, UnitX_

    integer, intent(in) :: nDimIn                ! dimension of positions
    integer, intent(in) :: nPoint                ! number of positions
    real,    intent(in) :: Xyz_DI(nDimIn,nPoint) ! positions
    integer, intent(out):: iProc_I(nPoint)       ! processor owning position

    ! Find array of points and return processor indexes owning them
    ! Could be generalized to return multiple processors...

    real:: Xyz_D(MaxDim) = 0.0
    integer:: iPoint, iBlock

    character(len=*), parameter:: NameSub = 'EE_find_points'
    !--------------------------------------------------------------------------
    do iPoint = 1, nPoint
       Xyz_D(1:nDimIn) = Xyz_DI(:,iPoint)*Si2No_V(UnitX_)
       call find_grid_block(Xyz_D, iProc_I(iPoint), iBlock)
    end do

  end subroutine EE_find_points
  !============================================================================
  subroutine EE_get_grid_info(nDimOut, iGridOut, iDecompOut)

    use EE_BATL_lib, ONLY: nDim
    use EE_ModMain,  ONLY: iNewGrid, iNewDecomposition

    integer, intent(out):: nDimOut    ! grid dimensionality
    integer, intent(out):: iGridOut   ! grid index (increases with AMR)
    integer, intent(out):: iDecompOut ! decomposition index

    ! Return basic grid information useful for model coupling.
    ! The decomposition index increases with load balance and AMR.
    character(len=*), parameter:: NameSub = 'EE_get_grid_info'
    !--------------------------------------------------------------------------

    nDimOut    = nDim
    iGridOut   = iNewGrid
    iDecompOut = iNewDecomposition

  end subroutine EE_get_grid_info
  !============================================================================
  subroutine EE_get_for_sc(IsNew, NameVar, nVarIn, nDimIn, nPoint, Xyz_DI, &
       Data_VI)

    ! Interpolate Data_VI from EE at the list of positions Xyz_DI
    ! required by SC

    use EE_ModPhysics, ONLY: Si2No_V, UnitX_, No2Si_V, iUnitCons_V
    use EE_ModAdvance, ONLY: State_VGB, Bx_, Bz_, nVar
    use EE_ModVarIndexes, ONLY: nVar
    use EE_ModB0,      ONLY: UseB0, get_b0
    use EE_BATL_lib,   ONLY: iProc, nDim, MaxDim, MinIJK_D, MaxIJK_D, &
         find_grid_block
    use EE_ModIO, ONLY: iUnitOut
    use CON_coupler, ONLY: nVarBuffer, iVarSource_V
    use ModInterpolate, ONLY: interpolate_vector

    logical,          intent(in):: IsNew   ! true for new point array
    character(len=*), intent(in):: NameVar ! List of variables
    integer,          intent(in):: nVarIn  ! Number of variables in Data_VI
    integer,          intent(in):: nDimIn  ! Dimensionality of positions
    integer,          intent(in):: nPoint  ! Number of points in Xyz_DI

    real, intent(in) :: Xyz_DI(nDimIn,nPoint)  ! Position vectors
    real, intent(out):: Data_VI(nVarIn,nPoint) ! Data array

    real:: Xyz_D(MaxDim), B0_D(MaxDim)
    real:: Dist_D(MaxDim), State_V(nVar)
    integer:: iCell_D(MaxDim)

    integer, allocatable, save:: iBlockCell_DI(:,:)
    real,    allocatable, save:: Dist_DI(:,:)

    integer:: iPoint, iBlock, iProcFound, iVarBuffer, iVar

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'EE_get_for_sc'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! If nDim < MaxDim, make sure that all elements are initialized
    Dist_D = -1.0
    Xyz_D  =  0.0

    if(IsNew)then
       if(DoTest)write(iUnitOut,*) NameSub,': iProc, nPoint=', iProc, nPoint

       if(allocated(iBlockCell_DI)) deallocate(iBlockCell_DI, Dist_DI)
       allocate(iBlockCell_DI(0:nDim,nPoint), Dist_DI(nDim,nPoint))

       do iPoint = 1, nPoint

          Xyz_D(1:nDim) = Xyz_DI(:,iPoint)*Si2No_V(UnitX_)
          call find_grid_block(Xyz_D, iProcFound, iBlock, iCell_D, Dist_D, &
               UseGhostCell = .true.)

          if(iProcFound /= iProc)then
             write(*,*) NameSub,' ERROR: Xyz_D, iProcFound=', Xyz_D, iProcFound
             call CON_stop(NameSub//' could not find position on this proc')
          end if

          ! Store block and cell indexes and distances for interpolation
          iBlockCell_DI(0,iPoint)      = iBlock
          iBlockCell_DI(1:nDim,iPoint) = iCell_D(1:nDim)
          Dist_DI(:,iPoint)            = Dist_D(1:nDim)

       end do
    end if

    do iPoint = 1, nPoint

       Xyz_D(1:nDim) = Xyz_DI(:,iPoint)*Si2No_V(UnitX_)

       ! Use stored block and cell indexes and distances
       iBlock          = iBlockCell_DI(0,iPoint)
       iCell_D(1:nDim) = iBlockCell_DI(1:nDim,iPoint)
       Dist_D(1:nDim)  = Dist_DI(:,iPoint)

       State_V = interpolate_vector(State_VGB(:,:,:,:,iBlock), nVar, nDim, &
            MinIJK_D, MaxIJK_D, iCell_D=iCell_D, Dist_D=Dist_D)

       if(UseB0)then
          call get_b0(Xyz_D, B0_D)
          State_V(Bx_:Bz_) = State_V(Bx_:Bz_) + B0_D
       end if

       do iVarBuffer = 1, nVarBuffer
          iVar = iVarSource_V(iVarBuffer)
          Data_VI(iVarBuffer,iPoint) = State_V(iVar)*No2Si_V(iUnitCons_V(iVar))
       end do

    end do

  end subroutine EE_get_for_sc
  !============================================================================
  subroutine EE_get_sc_region(NameVar, nVar, nPoint, Pos_DI, Data_VI, iPoint_I)

    ! This function will be called 3 times :
    !
    ! 1) Count outer boundary ghost cells to be obtained from SC
    !
    ! 2) Return the Xyz_DGB coordinates of these cells
    !
    ! 3) Recieve Data_VI from SC and put them into State_VGB.
    !    The indexing array iPoint_I is needed to maintain the same order as
    !    the original position array Pos_DI was given in 2)

    use EE_BATL_lib,     ONLY: Xyz_DGB, nBlock, Unused_B, &
         MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
         CoordMin_D, CoordMax_D, CoordMin_DB, CellSize_DB
    use EE_ModGeometry,  ONLY: IsBoundary_B, r_GB
    use EE_ModPhysics,   ONLY: No2Si_V, UnitX_, Si2No_V, iUnitCons_V
    use EE_ModMain,      ONLY: UseB0
    use EE_ModB0,        ONLY: B0_DGB
    use EE_ModAdvance,   ONLY: State_VGB, Bx_, Bz_
    use EE_ModMultiFluid, ONLY: nIonFluid
    use CON_coupler,     ONLY: Grid_C, SC_, nVarBuffer, iVarTarget_V
    character(len=*), intent(inout):: NameVar ! List of variables
    integer,          intent(inout):: nVar    ! Number of variables in Data_VI
    integer,          intent(inout):: nPoint  ! Number of points in Pos_DI
    real, intent(inout), allocatable, optional :: Pos_DI(:,:)  ! Positions

    real,    intent(in), optional:: Data_VI(:,:)    ! Recv data array
    integer, intent(in), optional:: iPoint_I(nPoint)! Order of data

    logical :: DoCountOnly
    integer :: i, j, k, iBlock, iPoint, iVarBuffer, iVar
    real    :: Coord_D(3), Xyz_D(3)
    real    :: rMinSc

    character(len=*), parameter:: NameSub = 'EE_get_sc_region'
    !--------------------------------------------------------------------------
    DoCountOnly = nPoint < 1

    rMinSc = Grid_C(SC_)%Coord2_I(1)

    ! Find ghost cells in the SC domain
    iPoint = 0
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       if(.not.IsBoundary_B(iBlock)) CYCLE
       do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
          ! Set generalized coordinate
          Coord_D = &
               CoordMin_DB(:,iBlock) + ([i,j,k]-0.5)*CellSize_DB(:,iBlock)

          ! Exclude points that are inside the domain
          if(all(Coord_D > CoordMin_D) .and. all(Coord_D < CoordMax_D)) CYCLE

          ! Exclude points below the SC bottom boundary
          if(r_GB(i,j,k,iBlock) < rMinSc) CYCLE

          ! Found a point to be set by SC
          iPoint = iPoint + 1
          if(DoCountOnly) CYCLE

          if(present(Data_VI))then
             ! Put Data_VI obtained from SC into State_VGB
             do iVarBuffer = 1, nVarBuffer
                iVar = iVarTarget_V(iVarBuffer)
                State_VGB(iVar,i,j,k,iBlock) = &
                     Data_VI(iVarBuffer,iPoint_I(iPoint)) &
                     *Si2No_V(iUnitCons_V(iVar))
             end do
             ! Set variables not defined by coupling
             ! if(UseElectronPressure .and. &
             !     .not. DoCoupleVar_V(ElectronPressure_)then
             if(UseB0) State_VGB(Bx_:Bz_,i,j,k,iBlock) = &
                  State_VGB(Bx_:Bz_,i,j,k,iBlock) - B0_DGB(:,i,j,k,iBlock)
          else
             ! Provide position to SC
             Pos_DI(:,iPoint) = Xyz_DGB(:,i,j,k,iBlock)*No2Si_V(UnitX_)
          end if

       end do; end do; end do
    end do

    if(DoCountOnly) nPoint = iPoint

  end subroutine EE_get_sc_region
  !============================================================================

  subroutine EE_put_from_sc( &
       NameVar, nVar, nPoint, Data_VI, iPoint_I, Pos_DI)

    use EE_BATL_lib,    ONLY: nDim

    character(len=*), intent(inout):: NameVar ! List of variables
    integer,          intent(inout):: nVar    ! Number of variables in Data_VI
    integer,          intent(inout):: nPoint  ! Number of points in Pos_DI

    real,    intent(in), optional:: Data_VI(:,:)    ! Recv data array
    integer, intent(in), optional:: iPoint_I(nPoint)! Order of data
    real, intent(out), allocatable, optional:: Pos_DI(:,:) ! Position vectors

    character(len=*), parameter:: NameSub = 'EE_put_from_sc'
    !--------------------------------------------------------------------------

    if(.not. present(Data_VI))then
       nPoint=0;
       ! get nPoint
       call EE_get_sc_region(NameVar, nVar, nPoint, Pos_DI)

       if(allocated(Pos_DI)) deallocate(Pos_DI)
       allocate(Pos_DI(nDim,nPoint))

       ! get Pos_DI
       call EE_get_sc_region(NameVar, nVar, nPoint, Pos_DI)

       RETURN
    end if

    ! set State variables
    call EE_get_sc_region(NameVar, nVar, nPoint, Pos_DI, Data_VI, iPoint_I)

  end subroutine EE_put_from_sc
  !============================================================================

  subroutine EE_use_pointer(iComp, tSimulation)

    use CON_coupler, ONLY: NameComp_I, Grid_C
    use EE_ModMain, ONLY: nVarComp2, NameVarComp2, StateComp2_VGB

    integer, intent(in):: iComp
    real,    intent(in):: tSimulation

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'EE_use_pointer'
    !--------------------------------------------------------------------------
    nVarComp2      =  Grid_C(iComp)%nVar
    NameVarComp2   =  Grid_C(iComp)%NameVar
    StateComp2_VGB => Grid_C(iComp)%State_VGB

    if(DoTestMe)then
       write(*,*) NameSub,' called from component    =', NameComp_I(iComp)
       write(*,*) NameSub,' nVarComp2, NameVarComp2  =',  nVarComp2, trim(NameVarComp2)
!!!       write(*,*) NameSub,' StateComp2_VGB(:,1,1,1,1)=', StateComp2_VGB(:,1,1,1,1)
    end if

    call EE_user_action('POINTERCOUPLING_'//NameComp_I(iComp))

  end subroutine EE_use_pointer
  !============================================================================

end module EE_wrapper
!==============================================================================
