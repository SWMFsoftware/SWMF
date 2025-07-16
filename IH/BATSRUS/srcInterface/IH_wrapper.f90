!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module IH_wrapper

  ! Wrapper for IH_BATSRUS Inner Heliosphere (IH) component
  ! This file gets copied into Solar Corona and Outer Heliosphere.

  use CON_coupler, ONLY: &
       CON_set_do_test, CON_stop, IH_, GridType, LocalGridType
  use IH_ModBatsrusMethods, ONLY: &
       BATS_init_session, BATS_setup, BATS_advance, BATS_save_files, &
       BATS_finalize
  use IH_ModBatsrusUtility, ONLY: stop_mpi, error_report
  use IH_ModUpdateStateFast, ONLY: sync_cpu_gpu

  implicit none
  SAVE
  private ! except

  ! CON wrapper
  public:: IH_set_param
  public:: IH_init_session
  public:: IH_run
  public:: IH_save_restart
  public:: IH_finalize

  ! Global buffer coupling
  public:: IH_get_for_global_buffer
  public:: IH_xyz_to_coord, IH_coord_to_xyz

  ! Coupling toolkit
  public:: IH_synchronize_refinement
  public:: IH_get_for_mh
  public:: IH_put_from_mh
  public:: IH_is_coupled_block
  public:: IH_interface_point_coords
  public:: IH_n_particle

  ! Public variables to be set/reset by a coupler. Needed to transform
  ! vector state variables obtained via the coupler.
  Character(len=3),    public :: TypeCoordSource    ! Coords of coupled model
  real,                public :: SourceToIH_DD(3,3) ! Transformation matrix
  real,                public :: TimeMhToIH = -1.0  ! Time of coupling

  type(GridType),      public :: IH_Grid            ! Grid (MHD data)
  type(GridType),      public :: IH_LineGrid        ! Global GD for lines
  type(LocalGridType), public :: IH_LocalGrid       ! Local GD (MHD data)
  type(LocalGridType), public :: IH_LocalLineGrid   ! Local GD for lines

  ! Coupling with SC
  public:: IH_set_buffer_grid_get_info
  public:: IH_save_global_buffer
  public:: IH_match_ibc

  ! Point coupling
  public:: IH_get_grid_info
  public:: IH_find_points

  ! Coupling with SP
  public:: IH_check_ready_for_sp
  public:: IH_extract_line
  public:: IH_put_particles
  public:: IH_get_particle_indexes
  public:: IH_get_particle_coords

  ! Coupling with GM
  public:: IH_get_for_gm

  ! Coupling with PT
  public:: IH_get_for_pt
  public:: IH_put_from_pt
  public:: IH_get_for_pt_dt

  ! Coupling with EE (for SC)
  public:: IH_get_for_ee
  public:: IH_put_from_ee

contains
  !============================================================================
  subroutine IH_init_session(iSession, TimeSimulation)

    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    logical :: IsUninitialized = .true.
    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'IH_init_session'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest, DoTestMe)

    if(IsUninitialized)then
       call BATS_setup
       IsUninitialized = .false.
    end if
    call BATS_init_session

    if(DoTest)write(*,*)NameSub,' finished for session ',iSession

  end subroutine IH_init_session
  !============================================================================
  subroutine IH_set_param(CompInfo, TypeAction)

    use CON_comp_info
    use IH_BATL_lib, ONLY: iProc, nProc, iComm
    use IH_ModIO, ONLY: iUnitOut, StringPrefix, STDOUT_, NamePlotDir
    use IH_ModSetParameters, ONLY: set_parameters
    use IH_ModRestartFile, ONLY: NameRestartInDir, NameRestartOutDir
    use IH_ModMain, ONLY : NameThisComp, &
         IsTimeAccurate, tSimulation, StartTime, iStartTime_I
    use CON_physics, ONLY: get_time
    use ModTimeConvert, ONLY: time_real_to_int

    type(CompInfoType), intent(inout):: CompInfo   ! Information for this comp.
    character(len=*), intent(in)     :: TypeAction ! What to do

    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'IH_set_param'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    if(DoTest)write(*,*)NameSub,' called with TypeAction,iProc=', &
         TypeAction,iProc

    select case(TypeAction)
    case('VERSION')
       call put(CompInfo, Use=.true., &
            NameVersion='IH_BATSRUS (Univ. of Michigan)')
    case('MPI')
       call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc, &
            Name=NameThisComp)

       NamePlotDir(1:2)       = NameThisComp
       NameRestartInDir(1:2)  = NameThisComp
       NameRestartOutDir(1:2) = NameThisComp
    case('READ','CHECK')
       call get_time(DoTimeAccurateOut=IsTimeAccurate, &
            tSimulationOut=tSimulation, tStartOut=StartTime)
       call time_real_to_int(StartTime, iStartTime_I)

       call set_parameters(TypeAction)
    case('STDOUT')
       iUnitOut=STDOUT_
       if(iProc == 0)then
          StringPrefix = NameThisComp//':'
       else
          write(StringPrefix,'(a,i4.4,a)')NameThisComp, iProc, ':'
       end if
    case('FILEOUT')
       call get(CompInfo,iUnitOut=iUnitOut)
       StringPrefix=''
    case('GRID')
       call IH_set_grid
    case default
       call CON_stop(NameSub//' SWMF_ERROR: invalid TypeAction='//TypeAction)
    end select

  end subroutine IH_set_param
  !============================================================================
  subroutine IH_finalize(TimeSimulation)

    use IH_ModMain, ONLY: IsTimeLoop

    real, intent(in):: TimeSimulation   ! seconds from start time

    integer :: iError
    character(len=*), parameter:: NameSub = 'IH_finalize'
    !--------------------------------------------------------------------------
    ! We are not advancing in time any longer
    IsTimeLoop = .false.

    call BATS_save_files('FINAL')

    call error_report('PRINT',0.,iError,.true.)

  end subroutine IH_finalize
  !============================================================================
  subroutine IH_save_restart(TimeSimulation)

    real, intent(in):: TimeSimulation   ! seconds from start time

    character(len=*), parameter:: NameSub = 'IH_save_restart'
    !--------------------------------------------------------------------------
    call BATS_save_files('RESTART')

  end subroutine IH_save_restart
  !============================================================================
  subroutine IH_run(TimeSimulation, TimeSimulationLimit)

    use IH_BATL_lib, ONLY: iProc
    use IH_ModMain, ONLY: tSimulation

    real, intent(inout):: TimeSimulation   ! current time of component

    real, intent(in):: TimeSimulationLimit ! simulation time not to be exceeded

    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'IH_run'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    if(DoTest)write(*,*)NameSub,' called with tSim, tSimLimit, iProc=',&
         TimeSimulation, TimeSimulationLimit, iProc

    if(abs(tSimulation-TimeSimulation)>0.0001) then
       write(*,*)NameSub, &
            ' IH time=',tSimulation,' SWMF time=',TimeSimulation
       call CON_stop(NameSub//': IH and SWMF simulation times differ')
    end if

    call BATS_advance(TimeSimulationLimit)

    ! Return time after the time step
    TimeSimulation = tSimulation

  end subroutine IH_run
  !============================================================================
  subroutine IH_get_grid_info(nDimOut, iGridOut, iDecompOut)

    use IH_BATL_lib, ONLY: nDim
    use IH_ModMain, ONLY: iNewGrid, iNewDecomposition

    integer, intent(out):: nDimOut    ! grid dimensionality
    integer, intent(out):: iGridOut   ! grid index (increases with AMR)
    integer, intent(out):: iDecompOut ! decomposition index

    ! Return basic grid information useful for model coupling.
    ! The decomposition index increases with load balance and AMR.
    character(len=*), parameter:: NameSub = 'IH_get_grid_info'
    !--------------------------------------------------------------------------

    nDimOut    = nDim
    iGridOut   = iNewGrid
    iDecompOut = iNewDecomposition

  end subroutine IH_get_grid_info
  !============================================================================
  subroutine IH_find_points(nDimIn, nPoint, Xyz_DI, iProc_I)

    use IH_BATL_lib, ONLY: MaxDim, find_grid_block
    use IH_ModPhysics, ONLY: Si2No_V, UnitX_

    integer, intent(in) :: nDimIn                ! dimension of positions
    integer, intent(in) :: nPoint                ! number of positions
    real,    intent(in) :: Xyz_DI(nDimIn,nPoint) ! positions
    integer, intent(out):: iProc_I(nPoint)       ! processor owning position

    ! Find array of points and return processor indexes owning them
    ! Could be generalized to return multiple processors...

    real:: Xyz_D(MaxDim) = 0.0
    integer:: iPoint, iBlock

    character(len=*), parameter:: NameSub = 'IH_find_points'
    !--------------------------------------------------------------------------
    do iPoint = 1, nPoint
       Xyz_D(1:nDimIn) = Xyz_DI(:,iPoint)*Si2No_V(UnitX_)
       call find_grid_block(Xyz_D, iProc_I(iPoint), iBlock)
    end do

  end subroutine IH_find_points
  !============================================================================
  subroutine IH_get_point_data(iBlockCell_DI, Dist_DI, IsNew, &
       NameVar, nVarIn, nDimIn, nPoint, Xyz_DI, Data_VI, DoSendAllVar)

    ! Generic routine for providing point data to another component
    ! If DoSendAllVar is true, send all variables in State_VGB
    ! Otherwise send the variables defined by iVarSource_V

    use IH_ModPhysics, ONLY: Si2No_V, UnitX_, UnitU_, No2Si_V, iUnitCons_V,UnitT_,No2Io_V
    use IH_ModAdvance, ONLY: State_VGB, Rho_, RhoUx_, RhoUz_, Bx_, Bz_
    use IH_ModVarIndexes, ONLY: nVar
    use IH_ModB0, ONLY: UseB0, get_b0, B0_DGB
    use IH_BATL_lib, ONLY: nDim, MaxDim, MinIJK_D, MaxIJK_D, iProc, &
         nBlock, MaxBlock, Unused_B, find_grid_block, &
         nI, nJ, nK, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nG, Xyz_DNB
    use IH_ModIO, ONLY: iUnitOut
    use CON_coupler, ONLY: iVarSource_V
    use ModInterpolate, ONLY: interpolate_vector, interpolate_scalar
    use IH_ModCellGradient, ONLY: calc_divergence
    use IH_BATL_pass_cell, ONLY: message_pass_cell
    use ModUtilities, ONLY: split_string

    logical,          intent(in):: IsNew   ! true for new point array
    integer, allocatable, intent(inout):: iBlockCell_DI(:,:) ! interp. index
    real,    allocatable, intent(inout):: Dist_DI(:,:)       ! interp. weight

    character(len=*), intent(in):: NameVar ! List of variables
    integer,          intent(in):: nVarIn  ! Number of variables in Data_VI
    integer,          intent(in):: nDimIn  ! Dimensionality of positions
    integer,          intent(in):: nPoint  ! Number of points in Xyz_DI

    real, intent(in) :: Xyz_DI(nDimIn,nPoint)  ! Position vectors
    real, intent(out):: Data_VI(nVarIn,nPoint) ! Data array

    logical, intent(in), optional:: DoSendAllVar

    logical:: DoSendAll

    real:: Xyz_D(MaxDim), B0_D(MaxDim)
    real:: Dist_D(MaxDim), State_V(nVar)
    integer:: iCell_D(MaxDim)

    ! calculation of divu
    logical:: UseDivU, UseDivUDx
    real:: DivU, DivUDx
    real, allocatable:: u_DG(:,:,:,:), DivU_GB(:,:,:,:), DivUDx_GB(:,:,:,:)
    integer:: i, j, k
    integer:: iPoint, iBlock, iProcFound, iVarBuffer, iVar

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'IH_get_point_data'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! Get updated State and B0 onto the CPU for the coupler
    call sync_cpu_gpu('update on CPU', NameSub, State_VGB, B0_DGB)

    DoSendAll = .false.
    if(present(DoSendAllVar)) DoSendAll = DoSendAllVar
    UseDivU   = index(NameVar//' ',' divu ') > 0
    UseDivUDx = index(NameVar,' divudx ') > 0

    if(DoTestMe)then
       write(*,*)NameSub,': DoSendAll, nVar, nVarIn, NameVar=', &
            DoSendAll, nVar, nVarIn, trim(NameVar)
       if(.not.DoSendAll) &
            write(*,*)NameSub,': iVarSource_V=', iVarSource_V(1:nVarIn)
       write(*,*)NameSub,': UseDivU=', UseDivU,' UseDivUDx=', UseDivUDx
    end if

    ! calculate divu if needed
    if (UseDivU .or. UseDivUDx) then
       allocate( &
            u_DG(3,MinI:MaxI,MinJ:MaxJ,MinK:MaxK), &
            DivU_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
       if(UseDivUDx) allocate( &
            DivUDx_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))

       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          ! Calculate velocity
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             u_DG(:,i,j,k) = State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock)/ &
                  State_VGB(Rho_,i,j,k,iBlock)
          end do; end do; end do

          ! Calculate div(u) in physical cell centers
          call calc_divergence(iBlock, u_DG, nG, DivU_GB(:,:,:,iBlock), &
               UseBodyCellIn=.true.)

          if (UseDivUdX) then
             ! Calculate DivU*CellSize
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                DivUdX_GB(i,j,k,iBlock) = DivU_GB(i,j,k,iBlock) * norm2( &
                     Xyz_DNB(:,i+1,j+1,k+1,iBlock)-Xyz_DNB(:,i,j,k,iBlock))
             end do; end do; end do
          end if

       end do
       deallocate(u_DG)
       ! Fill in ghost cells
       if(UseDivU)   call message_pass_cell(DivU_GB)
       if(UseDivUDx) call message_pass_cell(DivUDx_GB)
    end if

    ! If nDim < MaxDim, make sure that all elements are initialized
    Dist_D = -1.0
    Xyz_D  =  0.0

    if(IsNew)then
       ! Find points and store cell indexes and weights

       if(DoTest)write(iUnitOut,*) NameSub,': iProc, nPoint=', iProc, nPoint

       if(allocated(iBlockCell_DI)) deallocate(iBlockCell_DI, Dist_DI)
       allocate(iBlockCell_DI(0:nDim,nPoint), Dist_DI(nDim,nPoint))

       do iPoint = 1, nPoint

          Xyz_D(1:nDim) = Xyz_DI(:,iPoint)*Si2No_V(UnitX_)
          call find_grid_block(Xyz_D, iProcFound, iBlock, iCell_D, Dist_D, &
               UseGhostCell = .true.)

          if(iProcFound /= iProc)then
             write(*,*) NameSub,' ERROR: Xyz_D, iProcFound=', Xyz_D, iProcFound
             call stop_mpi(NameSub//' could not find position on this proc')
          end if

          ! Store block and cell indexes and distances for interpolation
          iBlockCell_DI(0,iPoint)      = iBlock
          iBlockCell_DI(1:nDim,iPoint) = iCell_D(1:nDim)
          Dist_DI(:,iPoint)            = Dist_D(1:nDim)

       end do
    end if

    ! Interpolate coupled variables to the point positions
    do iPoint = 1, nPoint
       Xyz_D(1:nDim) = Xyz_DI(:,iPoint)*Si2No_V(UnitX_)

       ! Use stored block and cell indexes and distances
       iBlock          = iBlockCell_DI(0,iPoint)
       iCell_D(1:nDim) = iBlockCell_DI(1:nDim,iPoint)
       Dist_D(1:nDim)  = Dist_DI(:,iPoint)

       ! Interpolate
       State_V = interpolate_vector(State_VGB(:,:,:,:,iBlock), nVar, &
            nDim, MinIJK_D, MaxIJK_D, iCell_D=iCell_D, Dist_D=Dist_D)
       if(UseDivU) DivU = interpolate_scalar(DivU_GB(:,:,:,iBlock), &
            nDim, MinIJK_D, MaxIJK_D, iCell_D=iCell_D, Dist_D=Dist_D)
       if(UseDivUdX) DivUdX = interpolate_scalar(DivUDx_GB(:,:,:,iBlock), &
            nDim, MinIJK_D, MaxIJK_D, iCell_D=iCell_D, Dist_D=Dist_D)

       ! Provide full B
       if(UseB0)then
          call get_b0(Xyz_D, B0_D)
          State_V(Bx_:Bz_) = State_V(Bx_:Bz_) + B0_D
       end if

       if(DoSendAll)then
          ! Fill buffer with interpolated values converted to SI units
          Data_VI(1:nVar,iPoint) = State_V*No2Si_V(iUnitCons_V)
          ! DivU is right after the nVar variables
          if(UseDivU) Data_VI(nVar+1,iPoint) = DivU/No2Io_V(UnitT_)
          ! DivUDx is the last one in the buffer
          if(UseDivUdX) Data_VI(nVarIn,iPoint) = DivUdX/No2Io_V(UnitT_)
       else
          do iVarBuffer = 1, nVarIn
             iVar = iVarSource_V(iVarBuffer)
             Data_VI(iVarBuffer,iPoint) = &
                  State_V(iVar)*No2Si_V(iUnitCons_V(iVar))
          end do
       end if
    end do ! iPoint

    if(allocated(DivU_GB))   deallocate(DivU_GB)
    if(allocated(DivUDx_GB)) deallocate(DivUDx_GB)

  end subroutine IH_get_point_data
  !============================================================================
  subroutine IH_set_grid

    use CON_comp_param
    use IH_domain_decomposition
    use CON_coupler
    use IH_ModMain, ONLY: TypeCoordSystem, nVar, NameVarCouple, nG
    use IH_ModPhysics, ONLY:No2Si_V, UnitX_
    use IH_ModGeometry, ONLY: RadiusMin, RadiusMax
    use IH_BATL_geometry, ONLY: TypeGeometry, IsGenRadius, LogRGen_I
    use IH_BATL_lib, ONLY: CoordMin_D, CoordMax_D, Particle_I,  &
         rRound0, rRound1
    use IH_ModParticleFieldLine, ONLY: iKindReg, UseParticles

    logical:: DoTest,DoTestMe
    logical:: UseParticleLine = .false.
    integer:: nParticle = 0, iError = 0
    integer:: IH_nG = 0
    character(len=*), parameter:: NameSub = 'IH_set_grid'
    !--------------------------------------------------------------------------
    DoTest=.false.;DoTestMe=.false.
    ! Here we should set the IH (MH) grid descriptor
    if(done_dd_init(IH_))then
       ! The coord system may be reset
       if(is_proc0(IH_)) Grid_C(IH_)%TypeCoord = TypeCoordSystem
       call MPI_bcast(Grid_C(IH_)%TypeCoord, 3, MPI_CHARACTER, &
            i_proc0(IH_), i_comm(), iError)
       RETURN
    end if

    call init_decomposition(&
         GridID_=IH_,&
         CompID_=IH_,&
         nDim=3,     &
         IsTreeDD = .true.)

    if(IsGenRadius)then
       call set_coord_system(            &
            GridID_   = IH_,             &
            TypeCoord = TypeCoordSystem, &
            UnitX     = No2Si_V(UnitX_), &
            nVar      = nVar,            &
            NameVar   = NameVarCouple,   &
            TypeGeometry = TypeGeometry, &
            Coord1_I  = LogRGen_I,       &
            Coord2_I  = [RadiusMin, RadiusMax])
    elseif(TypeGeometry=='roundcube')then
       call set_coord_system(&
            GridID_   = IH_,             &
            TypeCoord = TypeCoordSystem, &
            UnitX     = No2Si_V(UnitX_), &
            nVar      = nVar,            &
            NameVar   = NameVarCouple,   &
            TypeGeometry = TypeGeometry, &
            Coord2_I  = [rRound0, rRound1])
    else
       call set_coord_system(&
            GridID_   = IH_,             &
            TypeCoord = TypeCoordSystem, &
            UnitX     = No2Si_V(UnitX_), &
            nVar      = nVar,            &
            NameVar   = NameVarCouple,   &
            TypeGeometry = TypeGeometry, &
            Coord2_I  = [RadiusMin, RadiusMax])
    end if

    if(is_Proc(IH_))then
       ! Initialize the local grid

       call init_decomposition(Domain=MH_Domain,&
            CompID_=IH_, nDim=3, IsTreeDD=.true.)

       ! Get the octree root array
       call MH_get_root_decomposition(MH_Domain)

       ! Get the whole octree after the initial refinement
       call MH_update_local_decomposition(MH_Domain)

       MH_Domain%IsLocal=.true.
    end if

    call CON_set_do_test('test_grids', DoTest, DoTestMe)
    ! Repeat the initialization at the global grid level:
    ! Octree root array:
    if(is_proc0(IH_))call MH_get_root_decomposition(IH_)

    ! Broadcast root array:
    call bcast_decomposition(IH_)

    ! Synchronize global and local grids:
    call synchronize_refinement(GridID_=IH_, LocalDomain=MH_Domain)
    if(is_proc0(IH_))IH_nG = nG
    call MPI_bcast(IH_nG, 1, MPI_INTEGER, i_proc0(IH_), i_comm(), iError)
    call set_standard_grid_descriptor(IH_,  & ! CompID_
         IH_nG,                             & ! Gcn
         CellCentered_,                     & ! Grid type
         IH_Grid)
    if(is_proc(IH_))call set_local_gd(iProc = i_proc(),   &
         Grid = IH_Grid, LocalGrid = IH_LocalGrid)
    if(is_proc0(IH_))UseParticleLine = UseParticles
    call MPI_bcast(UseParticleLine, 1, MPI_LOGICAL,&
         i_proc0(IH_), i_comm(), iError)
    if(UseParticleLine)then
       call init_decomposition_dd(MH_LineDecomposition, IH_, nDim=1)
       if(is_proc0(IH_))then
          nParticle = Particle_I(iKindReg)%nParticleMax
          call get_root_decomposition_dd(MH_LineDecomposition, &
               [n_proc(IH_)],        &  ! One "block" per processor
               [0.50],               &  ! "Coordinate" is a global point number
                                ! factors are converted separately to prevent integer overflow
               [real(n_proc(IH_))*real(nParticle) + 0.50], &
               [nParticle])             ! nParticle cells per "block" (proc)
       end if
       call bcast_decomposition_dd(MH_LineDecomposition)
       if(DoTest.and.is_proc0(IH_))call show_domain_decomp(    &
            MH_LineDecomposition)
       ! Set local GD on the Particle_I structure
       call set_standard_grid_descriptor(MH_LineDecomposition, &
            Grid=IH_LineGrid)
       if(is_proc(IH_))call set_local_gd(iProc=i_proc(),     &
            Grid=IH_LineGrid, LocalGrid=IH_LocalLineGrid)
    end if

  end subroutine IH_set_grid
  !============================================================================
  subroutine IH_xyz_to_coord(XyzIn_D, CoordOut_D)

    use CON_coupler
    use ModCoordTransform, ONLY: &
         atan2_check, xyz_to_sph, xyz_to_rlonlat
    real, intent(in ) :: XyzIn_D(3)
    real, intent(out) :: CoordOut_D(3)

    real:: x, y, Coord2_I(2)
    integer, parameter:: x_=1, y_=2, z_=3, r_=1
    integer:: Phi_
    character(len=20)  :: TypeGeometry
    character(len=*), parameter:: NameSub = 'IH_xyz_to_coord'
    !--------------------------------------------------------------------------
    TypeGeometry = Grid_C(IH_)%TypeGeometry
    if(TypeGeometry(1:9)  == 'cartesian')then
       CoordOut_D = XyzIn_D
       RETURN
    elseif(TypeGeometry(1:3)  == 'cyl')then
       Phi_ = 2
       x = XyzIn_D(x_); y = XyzIn_D(y_)
       CoordOut_D(r_)   = sqrt(x**2 + y**2)
       CoordOut_D(Phi_) = atan2_check(y, x)
       CoordOut_D(z_)   = XyzIn_D(z_)
    elseif(TypeGeometry(1:3)  == 'sph')then
       call xyz_to_sph(XyzIn_D, CoordOut_D)
    elseif(TypeGeometry(1:3)  == 'rlo')then
       call xyz_to_rlonlat(XyzIn_D, CoordOut_D)
    elseif(TypeGeometry(1:9) == 'roundcube')then
       Coord2_I = Grid_C(IH_)%Coord2_I
       call xyz_to_roundcube_coord(&
           XyzIn_D, Coord2_I(1), Coord2_I(2), CoordOut_D)
    else
       call CON_stop(NameSub// &
            ' not yet implemented for TypeGeometry='//TypeGeometry)
    end if
    if(index(TypeGeometry,'lnr')  > 0)then
       CoordOut_D(1) = log(max(CoordOut_D(1), 1e-30))
    elseif(index(TypeGeometry,'genr') > 0)then
       call radius_to_gen(CoordOut_D(1))
    end if

  contains
    !==========================================================================
    subroutine radius_to_gen(r)

      use ModInterpolate, ONLY: find_cell
      real, intent(inout) :: r

      integer       :: nRgen
      integer       :: i
      real          :: dCoord
      real, pointer :: LogRgen_I(:)
      !------------------------------------------------------------------------
      LogRgen_I => Grid_C(IH_)%Coord1_I
      nRgen    = Grid_C(IH_)%nCoord_D(1)
      call find_cell(0, nRgen-1, alog(r), &
           i, dCoord, LogRgen_I, DoExtrapolate=.true.)
      r = (i + dCoord)/(nRgen - 1)

    end subroutine radius_to_gen
    !==========================================================================
    subroutine xyz_to_roundcube_coord(&
         XyzIn_D, rRound0, rRound1, CoordOut_D)

      real, intent(in ) :: XyzIn_D(3), rRound0, rRound1
      real, intent(out) :: CoordOut_D(3)
      real, parameter   :: SqrtNDim = sqrt(3.0)
      real:: r2, Dist1, Dist2, Coef1, Coef2
      !------------------------------------------------------------------------
      r2 = sum(XyzIn_D**2)
      if (r2 > 0.0) then
         ! L1 and L2 distance
         Dist1 = maxval(abs(XyzIn_D))
         Dist2 = sqrt(r2)
         if (rRound1 > rRound0 ) then
            ! The rounded grid is outside of the non-distorted part
            if (Dist1 > rRound0) then
               ! Outside the undistorted region
               ! Assume Coord = w * Xyz and Replace Xyz in coord_to_xyz
               ! We have a quadratic equation of w.
               ! w^2 - Coef1*w
               ! - 4*Dist1/(rRound1-rRound0)*(SqrtNDim*Dist1/Dist2 - 1) = 0

               Coef1 = -1 &
                    + rRound0/(rRound1-rRound0)*(dist1*SqrtNDim/Dist2 - 1)
               Coef2 = Coef1**2 &
                    + 4*Dist1/(rRound1-rRound0)*(SqrtNDim*Dist1/Dist2 - 1)
               CoordOut_D = XyzIn_D/(-Coef1 + sqrt(Coef2))*2
            else
               ! No distortion
               CoordOut_D = XyzIn_D
            end if

         else
            ! The rounded (distorted) grid is inside of the non-distorted part
            if (Dist2 < rRound1) then
               ! Solving w^2-w+Coef1 = 0
               Coef1 = Dist1/rRound1*(1 - Dist1/Dist2)
               CoordOut_D = XyzIn_D / (1+sqrt(1-4*Coef1))*2
            else
               ! Solving w^2+Coef1*w+Coef2 = 0
               Coef1 = -1 + (1 - Dist1/Dist2)/(rRound0-rRound1)*rRound0
               Coef2 = -(1 - Dist1/Dist2)/(rRound0 - rRound1)*Dist1
               Coef2 = (-Coef1 + sqrt(Coef1**2 - 4*Coef2))*0.5
               CoordOut_D = XyzIn_D / Coef2
            end if
         end if

      else
         CoordOut_D = 0.0
      end if
    end subroutine xyz_to_roundcube_coord
    !==========================================================================
  end subroutine IH_xyz_to_coord
  !============================================================================
  subroutine IH_coord_to_xyz(CoordIn_D, XyzOut_D)

    use CON_coupler
    use ModCoordTransform, ONLY: sph_to_xyz, rlonlat_to_xyz
    real, intent(in) :: CoordIn_D(3)
    real, intent(out):: XyzOut_D( 3)

    real               :: Coord_D(3), r, Phi, Coord2_I(2)
    integer, parameter :: x_=1, r_=1
    integer            :: Phi_
    character(len=20)  :: TypeGeometry
    character(len=*), parameter:: NameSub = 'IH_coord_to_xyz'
    !--------------------------------------------------------------------------
    TypeGeometry = Grid_C(IH_)%TypeGeometry
    if(TypeGeometry(1:9)  == 'cartesian')then
       XyzOut_D = CoordIn_D
       RETURN
    endif

    Coord_D = CoordIn_D
    if(index(TypeGeometry,'lnr')  > 0)then
       Coord_D(1) = exp(Coord_D(1))
    elseif(index(TypeGeometry,'genr') > 0)then
       call gen_to_radius(Coord_D(1))
    end if

    if(TypeGeometry(1:3)  == 'cyl')then
       Phi_ = 2
       r = Coord_D(r_); Phi = Coord_D(Phi_)
       XyzOut_D(1) = r*cos(Phi)
       XyzOut_D(2) = r*sin(Phi)
       XyzOut_D(3) = Coord_D(3)
    elseif(TypeGeometry(1:3)  == 'sph')then
       call sph_to_xyz(Coord_D, XyzOut_D)
    elseif(TypeGeometry(1:3)  == 'rlo')then
       call rlonlat_to_xyz(Coord_D, XyzOut_D)
    elseif(TypeGeometry(1:9) == 'roundcube')then
       Coord2_I = Grid_C(IH_)%Coord2_I
       call roundcube_coord_to_xyz(&
            Coord_D, Coord2_I(1), Coord2_I(2), XyzOut_D)
    else
       call CON_stop(NameSub// &
            ' not yet implemented for TypeGeometry='//TypeGeometry)
    end if
  contains
    !==========================================================================
    subroutine gen_to_radius(r)

      use ModInterpolate, ONLY: linear

      ! Convert generalized radial coordinate to true radial coordinate
      real, intent(inout):: r
      integer       :: nRgen
      real, pointer :: LogRgen_I(:)
      !------------------------------------------------------------------------
      LogRgen_I => Grid_C(IH_)%Coord1_I
      nRgen = Grid_C(IH_)%nCoord_D(1)

      ! interpolate the LogRgen_I array for the general coordinate
      r = exp(linear(LogRgen_I, 0, nRgen-1, r*(nRgen-1), DoExtrapolate=.true.))

    end subroutine gen_to_radius
    !==========================================================================
    subroutine roundcube_coord_to_xyz(&
         CoordIn_D, rRound0, rRound1, XyzOut_D)

      real, intent(in) :: CoordIn_D(3), rRound0, rRound1
      real, intent(out):: XyzOut_D( 3)
      real:: r2, Dist1, Dist2, Weight
      real, parameter   :: SqrtNDim = sqrt(3.0)
      !------------------------------------------------------------------------
      r2 = sum(CoordIn_D**2)
      ! L1 and L2 distances from origin
      ! L1 distance is constant on the surface of a cube
      ! L2 distance is constant on the surface of a sphere
      Dist1 = maxval(abs(CoordIn_D))
      Dist2 = sqrt(r2)

      if (r2 > 0.0) then
         if (rRound0 < rRound1) then
            ! Non-distorted grid inside, round grid outside
            Weight = (Dist1 - rRound0)/(rRound1 - rRound0)
         elseif (Dist1 < rRound1) then
            ! the rounded grid is inside and the point is inside rRound1
            ! The distortion is 0 at the origin and maximum at Dist1=rRound1.
            Weight = Dist1/rRound1
         else
            ! the rounded grid is inside and the point is outside rRound1
            Weight = (rRound0 - Dist1)/(rRound0 - rRound1)
         endif
         ! Limit weight to be in the [0,1] interval
         Weight = min(1., max(0.0, Weight))

         if (rRound0 < rRound1) then
            ! Expand coordinate outward
            ! For a fully rounded grid we expand the generalized coordinate
            ! by Dist1*SqrtNDim/Dist2, so along the main diagonals there is
            ! no stretch and along the axes the expansion is SqrtNDim.
            ! For the partially rounded grid the expansion factor is reduced.
            ! The minimum expansion factor is 1 in the non-distorted region.
            XyzOut_D = (1 + Weight*(Dist1*SqrtNDim/Dist2 - 1)) * CoordIn_D
         else
            ! Contract coordinate inward
            ! In this case the grid is contracted along
            ! the main diagonals by a factor up to sqrt(nDim)
            ! and there is no contraction along the axes
            XyzOut_D = (1 + Weight*(Dist1/Dist2 - 1)) * CoordIn_D
         end if

      else
         XyzOut_D = 0.0
      end if
    end subroutine roundcube_coord_to_xyz
    !=========================================================================
  end subroutine IH_coord_to_xyz
  !============================================================================
  subroutine IH_synchronize_refinement(iProc0,iCommUnion)

    use IH_domain_decomposition
    use CON_comp_param

    integer, intent(in) ::iProc0,iCommUnion

    ! Synchronize the local grid decomposition to accomodate the
    ! grid change.
    !--------------------------------------------------------------------------
    if(is_proc(IH_)) &
         call MH_update_local_decomposition(MH_Domain)

    call synchronize_refinement(GridID_=IH_, &
         LocalDomain=MH_Domain, iProcUnion=iProc0, iCommUnion=iCommUnion)

  end subroutine IH_synchronize_refinement
  !============================================================================
  subroutine IH_set_buffer_grid_get_info(nR, nLon, nLat, BufferMinMax_DI)

    use IH_domain_decomposition, ONLY: is_proc
    use IH_ModBuffer,            ONLY: BuffR_, nRBuff, nLonBuff, nLatBuff,&
         BufferMin_D, BufferMax_D
    ! spherical buffer coupling
    use IH_ModBuffer,      ONLY: TypeCoordSource
    use CON_coupler,       ONLY: iCompSourceCouple,  Grid_C

    integer, intent(out)    :: nR, nLon, nLat
    real, intent(out)       :: BufferMinMax_DI(3,2)

    logical  :: DoTest, DoTestMe

    character(len=*), parameter:: NameSub = 'IH_set_buffer_grid_get_info'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! Return buffer size and limits to SWMF calling routine
    BufferMinMax_DI(:,1) = BufferMin_D
    BufferMinMax_DI(:,2) = BufferMax_D

    nR   = nRBuff
    nLon = nLonBuff
    nLat = nLatBuff

    if(DoTest) then
       write(*,*) NameSub,': with nR, nLon, nLat = ',nR, nLon, nLat
       write(*,*) 'BufferMin_D: ', BufferMin_D
       write(*,*) 'BufferMax_D: ', BufferMax_D
    end if
    TypeCoordSource = Grid_C(iCompSourceCouple) % TypeCoord

  end subroutine IH_set_buffer_grid_get_info
  !============================================================================
  subroutine IH_save_global_buffer(nVarCouple, nR, nLon, nLat, BufferIn_VC)

    use IH_ModBuffer, ONLY: BufferState_VG, fill_in_buffer_grid_gc
    use IH_ModMessagePass, ONLY: exchange_messages
    ! spherical buffer coupling
    use CON_coupler, ONLY: iVar_V, DoCoupleVar_V
    use IH_ModAdvance, ONLY: nVar, &
         UseElectronPressure, UseAnisoPressure, UseMultiSpecies
    use IH_ModSaMhd, ONLY: UseSaMhd
    use IH_ModVarIndexes, ONLY: &
         Rho_, Ux_, Uz_, RhoUx_, RhoUz_, Bx_, Bz_, p_,             &
         WaveFirst_, WaveLast_, Pe_, Ppar_, nFluid, BperU_, Ehot_, &
         ChargeStateFirst_, ChargeStateLast_, WDiff_, Lperp_, nIonFluid
    use CON_coupler,   ONLY: &
         Bfield_, ElectronPressure_, AnisoPressure_, Wave_, &
         MultiFluid_, MultiSpecie_, CollisionlessHeatFlux_, SaMhd_, &
         DoLperp_, DoWDiff_, LperpCouple_, WDiffCouple_, &
         RhoCouple_, RhoUxCouple_, RhoUzCouple_, PCouple_, &
         BxCouple_, BzCouple_, PeCouple_, PparCouple_, SaMhdCouple_, &
         WaveFirstCouple_, WaveLastCouple_, EhotCouple_, &
         ChargeState_, ChargeStateFirstCouple_, ChargeStateLastCouple_
    use IH_ModMultiFluid, ONLY: IsFullyCoupledFluid, iRho_I, iP_I
    use IH_ModPhysics, ONLY: No2Si_V, Si2No_V, UnitRho_, UnitB_, UnitX_
    use IH_ModPhysics, ONLY: UnitRhoU_, UnitEnergyDens_, UnitP_, UnitU_, &
         BodyRho_I, BodyP_I
    use IH_ModMain, ONLY: UseOuterHelio

    integer, intent(in) :: nVarCouple, nR, nLon, nLat
    real, intent(in)    :: BufferIn_VC(nVarCouple,nR,nLon,nLat)

    integer :: iFluid

    ! Convert from SI units to normalized units
    character(len=*), parameter:: NameSub = 'IH_save_global_buffer'
    !--------------------------------------------------------------------------
    BufferState_VG(Rho_,:,1:nLon,1:nLat) = &
         BufferIn_VC(iVar_V(RhoCouple_),:,:,:) &
         *Si2No_V(UnitRho_)
    ! Transform to primitive variables
    BufferState_VG(RhoUx_:RhoUz_,:,1:nLon,1:nLat) = &
         BufferIn_VC(iVar_V(RhoUxCouple_):iVar_V(RhoUzCouple_),:,:,:) &
         *Si2No_V(UnitRhoU_)
    if(DoCoupleVar_V(Bfield_)) &
         BufferState_VG(Bx_:Bz_,:,1:nLon,1:nLat) = &
         BufferIn_VC(iVar_V(BxCouple_):iVar_V(BzCouple_),:,:,:) &
         *Si2No_V(UnitB_)

    if(DoCoupleVar_V(Wave_)) &
         BufferState_VG(WaveFirst_:WaveLast_,:,1:nLon,1:nLat) = BufferIn_VC( &
         iVar_V(WaveFirstCouple_):iVar_V(WaveLastCouple_),:,:,:) &
         *Si2No_V(UnitEnergyDens_)

    if(DoCoupleVar_V(ChargeState_)) &
         BufferState_VG(ChargeStateFirst_:ChargeStateLast_,:,1:nLon,1:nLat) = &
         BufferIn_VC(iVar_V(ChargeStateFirstCouple_) &
         :iVar_V(ChargeStateLastCouple_),:,:,:)&
         *Si2No_V(UnitRho_)
    if(DoCoupleVar_V(SaMhd_))then
       BufferState_VG(BperU_,:,1:nLon,1:nLat) = &
            BufferIn_VC(iVar_V(SaMhdCouple_),:,:,:)* &
            Si2No_V(UnitB_)/Si2No_V(UnitU_)
    elseif(UseSaMhd .and. BperU_ > 1)then
       BufferState_VG(BperU_,:,1:nLon,1:nLat) = 0.0
    end if
    if(DoCoupleVar_V(DoLperp_))then
       BufferState_VG(Lperp_,:,1:nLon,1:nLat) = &
            BufferIn_VC(iVar_V(LperpCouple_),:,:,:)*&
            Si2No_V(UnitX_)*sqrt(Si2No_V(UnitB_))
    end if
    if(DoCoupleVar_V(DoWDiff_))then
       BufferState_VG(WDiff_,:,1:nLon,1:nLat) = &
            BufferIn_VC(iVar_V(WDiffCouple_),:,:,:)*&
            Si2No_V(UnitEnergyDens_)
    elseif(WDiff_>1)then
       BufferState_VG(WDiff_,:,1:nLon,1:nLat) = 0.0
    end if

    BufferState_VG(p_,:,1:nLon,1:nLat) = &
         BufferIn_VC(iVar_V(PCouple_),:,:,:)&
         *Si2No_V(UnitP_)
    if(DoCoupleVar_V(ElectronPressure_))then
       BufferState_VG(Pe_,:,1:nLon,1:nLat) = &
            BufferIn_VC(iVar_V(PeCouple_),:,:,:)*Si2No_V(UnitP_)
    else if(UseElectronPressure)then
       ! Split pressure equally between ions and electrons
       BufferState_VG(Pe_,:,1:nLon,1:nLat) = &
            0.5*BufferState_VG(p_,:,1:nLon,1:nLat)
       BufferState_VG(p_ ,:,1:nLon,1:nLat) = &
            BufferState_VG(Pe_,:,1:nLon,1:nLat)
    end if

    if(DoCoupleVar_V(AnisoPressure_))then
       BufferState_VG(Ppar_,:,1:nLon,1:nLat) = &
            BufferIn_VC(iVar_V(PparCouple_),:,:,:)*Si2No_V(UnitP_)
    else if(UseAnisoPressure)then
       BufferState_VG(Ppar_,:,1:nLon,1:nLat) = &
            BufferIn_VC(iVar_V(PCouple_),:,:,:)*Si2No_V(UnitP_)
    end if

    if(DoCoupleVar_V(CollisionlessHeatFlux_))then
       BufferState_VG(Ehot_,:,1:nLon,1:nLat) = &
            BufferIn_VC(iVar_V(EhotCouple_),:,:,:)&
            *Si2No_V(UnitEnergyDens_)
    endif

    if(UseOuterHelio)then
       do iFluid = nIonFluid + 1, nFluid
          BufferState_VG(iRho_I(iFluid),:,1:nLon,1:nLat) = BodyRho_I(iFluid)
          BufferState_VG(iP_I(iFluid),:,1:nLon,1:nLat)   = BodyP_I(iFluid)
       end do
    end if

    if( .not. DoCoupleVar_V(MultiFluid_) .and. nFluid > 1 .or. &
         .not. DoCoupleVar_V(MultiSpecie_) .and. UseMultiSpecies)then
       ! Values for neutrals / ions should be prescribed in set_BCs.f90
       IsFullyCoupledFluid = .false.
    else
       IsFullyCoupledFluid = .true.
    end if
    ! Make sure that ghost cells get filled after
    call fill_in_buffer_grid_gc
    ! Fill in the cells, covered by the buffer grid, including ghost cells
    call exchange_messages(UseBufferIn = .true.)

  end subroutine IH_save_global_buffer
  !============================================================================
  subroutine IH_match_ibc

    use IH_ModMessagePass, ONLY: exchange_messages
    use IH_ModIO, ONLY: IsRestartCoupler
    use IH_ModBuffer, ONLY: match_ibc
    use IH_ModAdvance, ONLY: State_VGB

    character(len=*), parameter :: StringTest ='IH_fill_buffer_only'
    logical:: DoTest, DoTestMe

    character(len=*), parameter:: NameSub = 'IH_match_ibc'
    !--------------------------------------------------------------------------
    if(IsRestartCoupler) RETURN
    call CON_set_do_test(StringTest, DoTest, DoTestMe)

    ! Fill in the physical cells, which are outside the buffer grid
    if(.not. DoTest) call match_ibc

    ! Fill in the cells, covered by the buffer grid, including ghost cells
    call exchange_messages(UseBufferIn=.true.)

  end subroutine IH_match_ibc
  !============================================================================
  subroutine IH_get_for_global_buffer( &
       nR, nLon, nLat, BufferMinMax_DI, Buffer_VG)

    ! This subroutines fills a buffer grid by interpolating from a source
    ! IH_BATSRUS grid using second-order trilinear interpolation.

    ! The buffer grid can be a spherical shell, or a segment of such a shell.

    ! All state variables in the source grid are interpolated, but only those
    ! needed for coupling (as determined by CON_coupler) are actually passed.

    ! The filled buffer state vector is converted to SI units and vector
    ! quantities are rotated to the target component coordinate system.

    ! INPUT:

    ! nR, nLon, nLat: grid spacing for the buffer grid
    ! BufferMinMAx_DI : Buffer grid minimum and maximum coordinates, in all
    ! dimensions.

    ! OUTPUT:

    ! Buffer_VG : defined for all coupling variables and all buffer grid points
    ! (including buffer ghost cells).

    ! REVISION HISTORY
    ! 30Dec2011 R. Oran   - initial version

    use IH_ModSize, ONLY: nI, nJ, nK, MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    use IH_ModMain, ONLY: UseB0, rUpperModel
    use IH_ModAdvance, ONLY: State_VGB, UseElectronPressure
    use IH_ModB0, ONLY: B0_DGB
    use IH_ModPhysics, ONLY: UnitX_,&
         No2Si_V, UnitRho_, UnitP_, UnitRhoU_, UnitB_, UnitEnergyDens_, UnitU_
    use IH_ModVarIndexes, ONLY: &
         Rho_, RhoUx_, RhoUz_, Bx_, Bz_, P_, Pe_, &
         Ppar_, WaveFirst_, WaveLast_, Ehot_, nVar, &
         ChargeStateFirst_, ChargeStateLast_, BperU_, WDiff_, Lperp_
    use CON_coupler, ONLY: &
         RhoCouple_, RhoUxCouple_,&
         RhoUzCouple_, PCouple_, BxCouple_, BzCouple_,  &
         PeCouple_, PparCouple_, WaveFirstCouple_,  &
         WaveLastCouple_, Bfield_, Wave_, EhotCouple_, &
         AnisoPressure_, ElectronPressure_,&
         CollisionlessHeatFlux_, ChargeStateFirstCouple_, &
         ChargeStateLastCouple_, ChargeState_, iVar_V, &
         DoCoupleVar_V, nVarCouple, SaMhd_, SaMhdCouple_ ,&
         DoLperp_, DoWDiff_, LperpCouple_, WDiffCouple_
    use ModCoordTransform, ONLY: rlonlat_to_xyz
    use ModInterpolate, ONLY: trilinear
    use IH_BATL_lib, ONLY: &
         iProc, find_grid_block, xyz_to_coord, CoordMin_DB, CellSize_DB

    ! Buffer size and limits
    integer,intent(in) :: nR, nLon, nLat
    real, intent(in)   :: BufferMinMax_DI(3,2)

    ! State variables to be fiiled in all buffer grid points
    real, intent(out):: Buffer_VG(nVarCouple,nR,nLon,nLat)

    ! variables for defining the buffer grid
    integer :: nCell_D(3)
    real    :: SphMin_D(3), SphMax_D(3), dSph_D(3), Sph_D(3)

    ! Variables for interpolating from a grid block to a buffer grid point

    ! Store complete interpolated state vector
    real :: StateInPoint_V(nVar)

    ! Store interpolated state variables needed for coupling
    real :: Buffer_V(nVarCouple), B0_D(3)

    ! Buffer grid cell center coordinates
    real :: CoordBuffer_D(3), XyzBuffer_D(3)

    ! Buffer grid cell center position  normalized by grid spacing
    ! (in IH_BATSRUS grid generalized coordinates)
    real :: BufferNorm_D(3)

    ! variable indices in buffer
    integer   :: &
         iRhoCouple,              &
         iRhoUxCouple,            &
         iRhoUzCouple,            &
         iPCouple,                &
         iPeCouple,               &
         iPparCouple,             &
         iBxCouple,               &
         iBzCouple,               &
         iWaveFirstCouple,        &
         iWaveLastCouple,         &
         iChargeStateFirstCouple, &
         iChargeStateLastCouple,  &
         iEhotCouple

    integer   :: iBlock, iPe, iR, iLon, iLat
    logical   :: DoTest, DoTestMe

    character(len=*), parameter:: NameSub = 'IH_get_for_global_buffer'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    Buffer_VG = 0.0

    ! get variable indices in buffer
    iRhoCouple              = iVar_V(RhoCouple_)
    iRhoUxCouple            = iVar_V(RhoUxCouple_)
    iRhoUzCouple            = iVar_V(RhoUzCouple_)
    iPCouple                = iVar_V(PCouple_)
    iPeCouple               = iVar_V(PeCouple_)
    iPparCouple             = iVar_V(PparCouple_)
    iBxCouple               = iVar_V(BxCouple_)
    iBzCouple               = iVar_V(BzCouple_)
    iWaveFirstCouple        = iVar_V(WaveFirstCouple_)
    iWaveLastCouple         = iVar_V(WaveLastCouple_)
    iEhotCouple             = iVar_V(EhotCouple_)
    iChargeStateFirstCouple = iVar_V(ChargeStateFirstCouple_)
    iChargeStateLastCouple  = iVar_V(ChargeStateLastCouple_)

    ! Calculate buffer grid spacing
    nCell_D  = [nR, nLon, nLat]
    SphMin_D = BufferMinMax_DI(:,1)
    SphMax_D = BufferMinMax_DI(:,2)

    ! Save the upper boundary radius as the limit for LOS integration span
    rUpperModel = SphMax_D(1)

    dSph_D     = (SphMax_D - SphMin_D)/real(nCell_D)
    dSph_D(1) = (SphMax_D(1) - SphMin_D(1))/(nCell_D(1) - 1)

    ! Get updated State and B0 onto the CPU for the coupler
    call sync_cpu_gpu('update on CPU', NameSub, State_VGB, B0_DGB)

    ! Loop over buffer grid points
    do iLat = 1, nLat ; do iLon = 1, nLon ; do iR = 1, nR

       ! Find the coordinates of the current buffer grid point,
       Sph_D = SphMin_D + [real(iR - 1), real(iLon)-0.5, real(iLat)-0.5]*dSph_D
       ! Find Xyz coordinates of the grid point
       call rlonlat_to_xyz(Sph_D, XyzBuffer_D)

       ! Find the block and PE in the IH_BATSRUS grid
       call find_grid_block(XyzBuffer_D, iPe, iBlock)

       ! Check if this block belongs to this processor
       if (iProc /= iPe) CYCLE

       ! Convert buffer grid point Xyz to IH_BATSRUS generalized coords
       call xyz_to_coord(XyzBuffer_D, CoordBuffer_D)

       ! Buffer grid point gen coords normalized by the block grid spacing
       BufferNorm_D = (CoordBuffer_D - CoordMin_DB(:,iBlock)) &
            /CellSize_DB(:,iBlock) + 0.5

       ! Interpolate from the true solution block to the buffer grid point
       StateInPoint_V = trilinear(State_VGB(:,:,:,:,iBlock),      &
            nVar, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, BufferNorm_D)

       ! Fill in the coupled state variables, convert to SI units

       Buffer_V(iRhoCouple)= StateInPoint_V(rho_)*No2Si_V(UnitRho_)
       Buffer_V(iRhoUxCouple:iRhoUzCouple) = &
            StateInPoint_V(rhoUx_:rhoUz_)*No2Si_V(UnitRhoU_)

       if(DoCoupleVar_V(Bfield_)) then
          if(UseB0)then
             B0_D = &
                  trilinear(B0_DGB(:,:,:,:,iBlock), &
                  3, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
                  BufferNorm_D, DoExtrapolate = .TRUE.)
             Buffer_V(iBxCouple:iBzCouple) = &
                  (StateInPoint_V(Bx_:Bz_) + B0_D)*No2Si_V(UnitB_)
          else
             Buffer_V(iBxCouple:iBzCouple) = &
                  StateInPoint_V(Bx_:Bz_)*No2Si_V(UnitB_)
          end if
       end if

       Buffer_V(iPCouple)  = StateInPoint_V(p_)*No2Si_V(UnitP_)

       if(DoCoupleVar_V(Wave_)) &
            Buffer_V(iWaveFirstCouple:iWaveLastCouple) = &
            StateInPoint_V(WaveFirst_:WaveLast_)&
            * No2Si_V(UnitEnergyDens_)

       if(DoCoupleVar_V(ChargeState_)) &
            Buffer_V(iChargeStateFirstCouple:iChargeStateLastCouple) = &
            StateInPoint_V(ChargeStateFirst_:ChargeStateLast_)&
            * No2Si_V(UnitRho_)

       if(DoCoupleVar_V(ElectronPressure_))then
          Buffer_V(iPeCouple) = StateInPoint_V(Pe_)*No2Si_V(UnitP_)
       else if(UseElectronPressure)then
          Buffer_V(iPCouple) = Buffer_V(iPCouple) + StateInPoint_V(Pe_)&
               *No2Si_V(UnitP_)
       end if

       if(DoCoupleVar_V(AnisoPressure_)) Buffer_V(iPparCouple) = &
            StateInPoint_V(Ppar_)*No2Si_V(UnitP_)

       if(DoCoupleVar_V(CollisionlessHeatFlux_)) Buffer_V(iEhotCouple) = &
            StateInPoint_V(Ehot_)*No2Si_V(UnitEnergyDens_)
       if(DoCoupleVar_V(SaMhd_))Buffer_V(iVar_V(SaMhdCouple_)) = &
            StateInPoint_V(BperU_)*No2Si_V(UnitB_)/No2Si_V(UnitU_)
       if(DoCoupleVar_V(DoLperp_))Buffer_V(iVar_V(LperpCouple_)) = &
            StateInPoint_V(Lperp_)*No2Si_V(UnitX_)*sqrt(No2Si_V(UnitB_))
       if(DoCoupleVar_V(DoWDiff_))Buffer_V(iVar_V(WDiffCouple_)) = &
            StateInPoint_V(WDiff_)*No2Si_V(UnitEnergydens_)

       ! DONE - fill the buffer grid
       Buffer_VG(:,iR, iLon,iLat) = Buffer_V

    end do; end do; end do

  end subroutine IH_get_for_global_buffer
  !============================================================================
  logical function IH_is_coupled_block(iBlock)

    use IH_ModMain, ONLY: iTypeCellBc_I
    use IH_ModParallel, ONLY: Unset_, DiLevel_EB
    use IH_ModGeometry, ONLY: IsBoundary_B
    use IH_ModCellBoundary, ONLY: CoupledBC_

    integer, intent(in) :: iBlock
    character(len=*), parameter:: NameSub = 'IH_is_coupled_block'
    !--------------------------------------------------------------------------
    if(.not.IsBoundary_B(iBlock)) then
       IH_is_coupled_block = .false.
       RETURN
    end if
    ! If block is near external boundary at which the BC type is 'Coupled'
    IH_is_coupled_block = any(DiLevel_EB(1:6, iBlock) == Unset_ &
         .and. iTypeCellBc_I(1:6) == CoupledBC_)

  end function IH_is_coupled_block
  !============================================================================
  subroutine IH_interface_point_coords(nDim, Xyz_D, nIndex, iIndex_I, &
       IsInterfacePoint)

    ! XYZ for the points beyond the boundary at which the BC 'coupled' is set

    use IH_ModMain, ONLY: iTypeCellBc_I
    use IH_ModParallel, ONLY: Unset_, DiLevel_EB
    use IH_BATL_lib, ONLY: nIJK_D, coord_to_xyz
    use IH_ModCellBoundary, ONLY: CoupledBC_

    integer, intent(in)   :: nDim
    real,    intent(inout):: Xyz_D(nDim)
    integer, intent(in)   :: nIndex
    integer, intent(inout):: iIndex_I(nIndex)
    logical, intent(out)  :: IsInterfacePoint

    logical :: IsRightCoupledBoundary_D(nDim)
    logical :: IsLeftCoupledBoundary_D(nDim)
    real    :: Coord_D(nDim)
    integer :: iBlock, iCell_D(nDim)

    character(len=*), parameter:: NameSub = 'IH_interface_point_coords'
    !--------------------------------------------------------------------------
    iCell_D = iIndex_I(1:nDim); iBlock    = iIndex_I(4)
    ! Check whether the point is in the ghost cell
    IsInterfacePoint = .not.(any(iCell_D < 1).or.any(iCell_D > nIJK_D))
    if(.not.IsInterfacePoint)RETURN
    ! Check, which boundaries are coupled
    IsLeftCoupledBoundary_D = DiLevel_EB(1:5:2, iBlock) == Unset_ &
         .and. iTypeCellBc_I(1:5:2) == CoupledBC_
    IsRightCoupledBoundary_D = DiLevel_EB(2:6:2, iBlock) == Unset_ &
         .and. iTypeCellBc_I(2:6:2) == CoupledBC_
    ! Check if the point is beyond one of these boundaries
    IsInterfacePoint = any(IsLeftCoupledBoundary_D.and.iCell_D < 1).or.&
         any(IsRightCoupledBoundary_D.and.iCell_D > nIJK_D)
    if(.not.IsInterfacePoint)RETURN
    ! Fix coordinates to be used in mapping
    Coord_D =  Xyz_D; call coord_to_xyz(CoorD_D, Xyz_D)

  end subroutine IH_interface_point_coords
  !============================================================================
  subroutine IH_get_for_mh(nPartial, iGetStart, Get, W, Buff_V, nVarIn)

    ! Put state variables to buffer for coupling toolkit,
    ! to send to other component via the coupling toolkit

    use IH_ModAdvance, ONLY: State_VGB, UseElectronPressure, nVar
    use IH_ModB0, ONLY: B0_DGB
    use IH_ModPhysics, ONLY: No2Si_V, UnitRho_, UnitP_, UnitRhoU_, UnitB_,&
         UnitU_, UnitX_
    use IH_ModPhysics, ONLY: UnitEnergyDens_
    use IH_ModAdvance, ONLY: Rho_, RhoUx_, RhoUz_, Bx_, Bz_, P_, WaveFirst_, &
         WaveLast_, Pe_, Ppar_, Ehot_, ChargeStateFirst_, ChargeStateLast_, &
         BperU_, Lperp_, WDiff_
    use IH_ModMain, ONLY: UseB0

    use CON_router, ONLY: IndexPtrType, WeightPtrType
    use CON_coupler, ONLY: iVar_V, DoCoupleVar_V, &
         RhoCouple_, RhoUxCouple_, &
         RhoUzCouple_, PCouple_, BxCouple_, BzCouple_, PeCouple_, PparCouple_,&
         WaveFirstCouple_, WaveLastCouple_, Bfield_, Wave_, AnisoPressure_, &
         ElectronPressure_, EhotCouple_, Momentum_, &
         CollisionlessHeatFlux_, ChargeStateFirstCouple_, &
         ChargeStateLastCouple_, ChargeState_, SaMhd_, SaMhdCouple_,&
         DoLperp_, DoWDiff_, LperpCouple_, WDiffCouple_

    integer,             intent(in) :: nPartial, iGetStart, nVarIn
    type(IndexPtrType),  intent(in) :: Get
    type(WeightPtrType), intent(in) :: W
    real,                intent(out):: Buff_V(nVarIn)

    integer   :: iGet, i, j, k, iBlock
    real      :: Weight, State_V(nVar)

    character(len=*), parameter:: NameSub = 'IH_get_for_mh'
    !--------------------------------------------------------------------------
    i      = Get%iCB_II(1,iGetStart)
    j      = Get%iCB_II(2,iGetStart)
    k      = Get%iCB_II(3,iGetStart)
    iBlock = Get%iCB_II(4,iGetStart)
    Weight = W%Weight_I(iGetStart)
    State_V = State_VGB(:,i,j,k,iBlock)*Weight
    if(UseB0)State_V(Bx_:Bz_) = State_V(Bx_:Bz_) + &
         B0_DGB(:,i,j,k,iBlock)*Weight

    do iGet=iGetStart+1,iGetStart+nPartial-1
       i      = Get%iCB_II(1,iGet)
       j      = Get%iCB_II(2,iGet)
       k      = Get%iCB_II(3,iGet)
       iBlock = Get%iCB_II(4,iGet)
       Weight = W%Weight_I(iGet)
       State_V = State_V + State_VGB(:,i,j,k,iBlock)*Weight
       if(UseB0)State_V(Bx_:Bz_) = State_V(Bx_:Bz_) + &
            B0_DGB(:,i,j,k,iBlock)*Weight
    end do
    Buff_V = 0.0

    ! Put variables into a buffer, convert to SI units
    Buff_V(iVar_V(RhoCouple_)) = State_V(Rho_)*No2Si_V(UnitRho_)
    if(DoCoupleVar_V(Momentum_)) Buff_V(iVar_V(RhoUxCouple_):                &
         iVar_V(RhoUzCouple_)) = State_V(rhoUx_:rhoUz_)*No2Si_V(UnitRhoU_)
    if(DoCoupleVar_V(Bfield_)) Buff_V(iVar_V(BxCouple_):iVar_V(BzCouple_)) = &
         State_V(Bx_:Bz_)*No2Si_V(UnitB_)
    if(DoCoupleVar_V(Wave_)) Buff_V(iVar_V(WaveFirstCouple_):                &
         iVar_V(WaveLastCouple_)) = State_V(WaveFirst_:WaveLast_)*           &
         No2Si_V(UnitEnergyDens_)
    if(DoCoupleVar_V(ChargeState_)) Buff_V(iVar_V(ChargeStateFirstCouple_):  &
         iVar_V(ChargeStateLastCouple_)) = &
         State_V(ChargeStateFirst_:ChargeStateLast_)*No2Si_V(UnitRho_)
    if(DoCoupleVar_V(AnisoPressure_)) Buff_V(iVar_V(PparCouple_)) =          &
         State_V(Ppar_)*No2Si_V(UnitP_)
    if(DoCoupleVar_V(CollisionlessHeatFlux_)) Buff_V(iVar_V(EhotCouple_)) =  &
         State_V(Ehot_)*No2Si_V(UnitEnergyDens_)
    if(DoCoupleVar_V(SaMhd_))Buff_V(iVar_V(SaMhdCouple_)) = State_V(BperU_)*&
         No2Si_V(UnitB_)/No2Si_V(UnitU_)
    if(DoCoupleVar_V(DoLperp_))Buff_V(iVar_V(LperpCouple_)) = State_V(Lperp_)*&
         No2Si_V(UnitX_)*sqrt(No2Si_V(UnitB_))
    if(DoCoupleVar_V(DoWDiff_))Buff_V(iVar_V(WDiffCouple_)) = &
         State_V(WDiff_)*No2Si_V(UnitEnergydens_)
    if(DoCoupleVar_V(ElectronPressure_))then
       Buff_V(iVar_V(PeCouple_)) = State_V(Pe_)*No2Si_V(UnitP_)
       Buff_V(iVar_V(PCouple_ )) = State_V(P_ )*No2Si_V(UnitP_)
    elseif(UseElectronPressure)then
       Buff_V(iVar_V(PCouple_ )) = &
            (State_V(Pe_) + State_V(P_))*No2Si_V(UnitP_)
    else
       Buff_V(iVar_V(PCouple_ )) = State_V(P_)*No2Si_V(UnitP_)
    end if

  end subroutine IH_get_for_mh
  !============================================================================
  subroutine IH_extract_line(Xyz_DI, iTraceMode, iIndex_II, RSoftBoundary)

    use IH_ModParticleFieldLine, &
         ONLY: extract_particle_line, RSoftBoundaryBats=>RSoftBoundary
    real,             intent(in) :: Xyz_DI(:,:)
    integer,          intent(in) :: iTraceMode
    integer,          intent(in) :: iIndex_II(:,:)
    real,             intent(in) :: RSoftBoundary
    ! set the soft boundary
    character(len=*), parameter:: NameSub = 'IH_extract_line'
    !--------------------------------------------------------------------------
    RSoftBoundaryBats = RSoftBoundary
    ! extract field lines starting at input points
    call extract_particle_line(Xyz_DI, iTraceMode, iIndex_II, &
         UseInputInGenCoord=.true.)

  end subroutine IH_extract_line
  !============================================================================
  subroutine IH_put_particles(Xyz_DI, iIndex_II)

    use IH_BATL_lib, ONLY: nDim, put_particles
    use IH_ModParticleFieldLine, ONLY: iKindReg
    ! add particles with specified coordinates to the already existing lines
    real,    intent(in):: Xyz_DI(:,:)
    integer, intent(in):: iIndex_II(:,:)
    !--------------------------------------------------------------------------
    call put_particles( &
         iKindParticle      = iKindReg , &
         StateIn_VI         = Xyz_DI   , &
         iIndexIn_II        = iIndex_II, &
         UseInputInGenCoord = .true.   , &
         DoReplace          = .true.     )

  end subroutine IH_put_particles
  !============================================================================
  subroutine IH_get_particle_indexes(iParticle, iIndex_I)

    use IH_ModParticleFieldLine, ONLY: fl_, id_, iKindReg
    use IH_BATL_particles, ONLY: Particle_I

    integer, intent(in) :: iParticle
    integer, intent(out):: iIndex_I(2)
    character(len=*), parameter:: NameSub = 'IH_get_particle_indexes'
    !--------------------------------------------------------------------------
    iIndex_I(1) = Particle_I(iKindReg)%iIndex_II(fl_, iParticle)
    iIndex_I(2) = Particle_I(iKindReg)%iIndex_II(id_, iParticle)

  end subroutine IH_get_particle_indexes
  !============================================================================
  subroutine IH_get_particle_coords(iParticle, Xyz_D)

    use IH_BATL_lib, ONLY: nDim
    use IH_ModParticleFieldLine, ONLY: iKindReg
    use IH_BATL_particles, ONLY: Particle_I

    integer, intent(in) :: iParticle
    real,    intent(out):: Xyz_D(nDim)
    !--------------------------------------------------------------------------
    Xyz_D = Particle_I(iKindReg)%State_VI(1:nDim, iParticle)

  end subroutine IH_get_particle_coords
  !============================================================================
  subroutine IH_put_from_mh( &
       nPartial, iPutStart, Put, Weight, DoAdd, Buff_V, nVarIn)

    ! transform and put the data got from MH

    use CON_axes, ONLY: transform_velocity
    use CON_router, ONLY: IndexPtrType, WeightPtrType
    use CON_coupler, ONLY: iVar_V, DoCoupleVar_V, &
         RhoCouple_, RhoUxCouple_, RhoUzCouple_, &
         PCouple_, BxCouple_, BzCouple_, PeCouple_, EhotCouple_,   &
         PparCouple_, WaveFirstCouple_, WaveLastCouple_, Momentum_,&
         Bfield_, Wave_, ElectronPressure_, AnisoPressure_, &
         CollisionlessHeatFlux_, ChargeStateFirstCouple_, &
         ChargeStateLastCouple_, ChargeState_, SaMhd_, SaMhdCouple_,&
         DoLperp_, DoWDiff_, WDiffCouple_, LperpCouple_
    use IH_ModAdvance, ONLY: State_VGB, UseElectronPressure, &
         UseAnisoPressure
    use IH_ModB0, ONLY: B0_DGB
    use IH_ModPhysics, ONLY: Si2No_V, UnitRho_, UnitP_, UnitRhoU_, UnitB_, &
         UnitEnergyDens_, UnitX_, No2Si_V, UnitU_
    use IH_ModMain, ONLY: UseB0, TypeCoordSystem
    use IH_ModSaMhd, ONLY: UseSaMhd, get_samhd_state
    use IH_ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUz_, Bx_, Bz_, P_, &
         WaveFirst_, WaveLast_, Pe_, Ppar_, Ehot_, ChargeStateFirst_, &
         ChargeStateLast_, nVar, BperU_, Ux_, Uz_, Lperp_, WDiff_
    use IH_ModGeometry,   ONLY: Xyz_DGB

    integer,             intent(in) :: nPartial, iPutStart, nVarIn
    type(IndexPtrType),  intent(in) :: Put
    type(WeightPtrType), intent(in) :: Weight
    logical,             intent(in) :: DoAdd
    real,                intent(in) :: Buff_V(nVarIn)

    ! revision history:
    ! 18JUL03     I.Sokolov <igorsok@umich.edu> - intial prototype/code
    ! 23AUG03                                     prolog
    ! 03SEP03     G.Toth    <gtoth@umich.edu>   - simplified
    ! 05APR11     R. Oran   <oran@umich.edu>    - Use coupling indices
    !                                          derived by the coupler according
    !                                           to actual variable names
    !                                           (see use CON_coupler).
    !                                           Handle anisotropic pressure.
    !

    real:: State_V(nVar), Xyz_D(3)
    integer:: i, j, k, iBlock, iRho

    character(len=*), parameter:: NameSub = 'IH_put_from_mh'
    !--------------------------------------------------------------------------
    i      = Put%iCB_II(1,iPutStart)
    j      = Put%iCB_II(2,iPutStart)
    k      = Put%iCB_II(3,iPutStart)
    iBlock = Put%iCB_II(4,iPutStart)
    ! Location:
    Xyz_D  = Xyz_DGB(:,i,j,k,iBlock)
    iRho   = iVar_V(RhoCouple_)    ! Reusable
    ! Copy from buffer in a proper order
    ! Convert state variable in buffer to normalized units.
    State_V = 0.0 ; State_V(Rho_) = Buff_V(iRho)*Si2No_V(UnitRho_)
    ! perform vector transformation from the source model to the IH one
    if(DoCoupleVar_V(BField_))State_V(Bx_:Bz_) =  matmul(SourceToIH_DD,     &
         Buff_V( iVar_V(BxCouple_):iVar_V(BzCouple_)))*Si2No_V(UnitB_)
    if(DoCoupleVar_V(Momentum_))State_V(rhoUx_:rhoUz_) = transform_velocity(&
         TimeMhToIH, Buff_V(iVar_V(RhoUxCouple_):iVar_V(RhoUzCouple_))/     &
         Buff_V(iRho), Xyz_D*No2Si_V(UnitX_),                               &
         TypeCoordSource, TypeCoordSystem)*Buff_V(iRho)*Si2No_V(UnitRhoU_)
    if(DoCoupleVar_V(Wave_)) State_V(WaveFirst_:WaveLast_) =      &
         Buff_V(iVar_V(WaveFirstCouple_):iVar_V(WaveLastCouple_)) &
         *Si2No_V(UnitEnergyDens_)
    if(DoCoupleVar_V(CollisionlessHeatFlux_)) State_V(Ehot_) =    &
         Buff_V(iVar_V(EhotCouple_))*Si2No_V(UnitEnergyDens_)
    if(DoCoupleVar_V(ChargeState_)) &
         State_V(ChargeStateFirst_:ChargeStateLast_) = Buff_V(           &
         iVar_V(ChargeStateFirstCouple_):iVar_V(ChargeStateLastCouple_)) &
         *Si2No_V(UnitRho_)
    if(DoCoupleVar_V(ElectronPressure_))then
       State_V(Pe_) = Buff_V(iVar_V(PeCouple_))*Si2No_V(UnitP_)
       State_V(p_ ) = Buff_V(iVar_V(PCouple_ ))*Si2No_V(UnitP_)
    elseif(UseElectronPressure)then
       State_V(p_ ) = 0.5*Buff_V(iVar_V(PCouple_))*Si2No_V(UnitP_)
       State_V(Pe_) = State_V(p_)
    else
       State_V(p_)  = Buff_V(iVar_V(PCouple_))*Si2No_V(UnitP_)
    end if
    if(DoCoupleVar_V(AnisoPressure_))then
       State_V(Ppar_) = Buff_V(iVar_V(PparCouple_))*Si2No_V(UnitP_)
    else if(UseAnisoPressure)then
       State_V(Ppar_) = State_V(p_)
    end if
    if(DoCoupleVar_V(SaMhd_))then
       State_V(BperU_) = Buff_V(iVar_V(SaMhdCouple_))*&
            Si2No_V(UnitB_)/Si2No_V(UnitU_)
    elseif(UseSaMhd)then
       ! Convert to primitive variables
       State_V(Ux_:Uz_) = State_V(RhoUx_:RhoUz_)/State_V(Rho_)
       call get_samhd_state(Xyz_D, State_V)
       ! Convert back to conservative variables
       State_V(RhoUx_:RhoUz_) = State_V(Ux_:Uz_)*State_V(Rho_)
    end if
    if(DoCoupleVar_V(DoLperp_))State_V(Lperp_) = Buff_V(iVar_V(LperpCouple_))*&
         Si2No_V(UnitX_)*sqrt(Si2No_V(UnitB_))
    if(DoCoupleVar_V(DoWDiff_))State_V(WDiff_) = &
         Buff_V(iVar_V(WDiffCouple_))*Si2No_V(UnitEnergydens_)

    ! ASSIGN the local state vector
    if(DoAdd)then
       State_VGB(:,i,j,k,iBlock) = State_VGB(:,i,j,k,iBlock) + State_V
    else
       State_VGB(:,i,j,k,iBlock) = State_V
       ! Full magnetic field in the buffer. Subtract B0, if used
       if(UseB0) State_VGB(Bx_:Bz_,i,j,k,iBlock) = &
            State_VGB(Bx_:Bz_,i,j,k,iBlock)  - B0_DGB(:,i,j,k,iBlock)
    end if

  end subroutine IH_put_from_mh
  !============================================================================
  subroutine IH_check_ready_for_sp(IsReady)

    use ModMpi
    use CON_coupler, ONLY: is_proc0, i_proc0, i_comm
    use IH_ModParticleFieldLine, ONLY: UseParticles
    logical, intent(out):: IsReady

    integer :: iError
    ! get value at IH root and broadcast to all SWMF processors
    !--------------------------------------------------------------------------
    if(is_proc0(IH_)) &
         IsReady = UseParticles
    call MPI_Bcast(IsReady, 1, MPI_LOGICAL, i_proc0(IH_), i_comm(), iError)

  end subroutine IH_check_ready_for_sp
  !============================================================================
  subroutine IH_get_for_gm( &
       nPartial,iGetStart,Get,W,State_V,nVar,TimeCoupling)

    use IH_ModAdvance, ONLY: State_VGB, Rho_, RhoUx_, RhoUz_, Bx_, Bz_,P_
    use IH_ModB0, ONLY: B0_DGB
    use IH_ModPhysics, ONLY: No2Si_V, UnitRho_, UnitP_, UnitRhoU_, UnitB_
    use IH_ModMain, ONLY: UseRotatingFrame,UseB0
    use CON_router

    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType), intent(in):: Get
    type(WeightPtrType), intent(in):: W
    real, intent(out):: State_V(nVar)
    real, intent(in):: TimeCoupling

    integer::iGet, i, j, k, iBlock
    real :: Weight, Momentum_D(3),Density

    ! The meaning of state intdex in buffer and in model can be
    ! different. Below are the conventions for buffer:
    integer, parameter:: &
         BuffRho_=1, BuffRhoUx_=2, BuffRhoUz_=4, BuffBx_=5, BuffBz_=7, BuffP_=8

    character(len=*), parameter:: NameSub = 'IH_get_for_gm'
    !--------------------------------------------------------------------------
    i      = Get%iCB_II(1,iGetStart)
    j      = Get%iCB_II(2,iGetStart)
    k      = Get%iCB_II(3,iGetStart)
    iBlock = Get%iCB_II(4,iGetStart)
    Weight = W%Weight_I(iGetStart)

    Density= State_VGB(rho_,         i,j,k,iBlock)
    State_V(BuffRho_)          = Density*Weight

    Momentum_D=State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock)
    if(UseRotatingFrame)call add_density_omega_cross_r

    State_V(BuffRhoUx_:BuffRhoUz_) = Momentum_D*Weight
    if(UseB0)then
       State_V(BuffBx_:BuffBz_) = &
            (State_VGB(Bx_:Bz_,i,j,k,iBlock) + B0_DGB(:,i,j,k,iBlock))*Weight
    else
       State_V(BuffBx_:BuffBz_) = &
            State_VGB(Bx_:Bz_,i,j,k,iBlock)*Weight
    end if

    State_V(BuffP_)  = State_VGB(P_,i,j,k,iBlock)*Weight

    do iGet=iGetStart+1,iGetStart+nPartial-1
       i      = Get%iCB_II(1,iGet)
       j      = Get%iCB_II(2,iGet)
       k      = Get%iCB_II(3,iGet)
       iBlock = Get%iCB_II(4,iGet)
       Weight = W%Weight_I(iGet)

       Density = State_VGB(rho_,i,j,k,iBlock)
       State_V(BuffRho_)=State_V(BuffRho_) + Density*Weight

       Momentum_D = State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock)

       if(UseRotatingFrame)call add_density_omega_cross_r

       State_V(BuffRhoUx_:BuffRhoUz_) = State_V(BuffRhoUx_:BuffRhoUz_) &
            + Momentum_D*Weight
       if(UseB0)then
          State_V(BuffBx_:BuffBz_) = State_V(BuffBx_:BuffBz_) &
               + (State_VGB(Bx_:Bz_,i,j,k,iBlock) &
               + B0_DGB(:,i,j,k,iBlock))*Weight
       else
          State_V(BuffBx_:BuffBz_) = State_V(BuffBx_:BuffBz_) &
               + State_VGB(Bx_:Bz_,i,j,k,iBlock)*Weight
       end if
       State_V(BuffP_) = State_V(BuffP_) &
            + State_VGB(P_,i,j,k,iBlock)*Weight
    end do

    ! Convert to SI units
    State_V(BuffRho_) = State_V(BuffRho_)*No2Si_V(UnitRho_)
    State_V(BuffRhoUx_:BuffRhoUz_)= State_V(BuffRhoUx_:BuffRhoUz_) &
         *No2Si_V(UnitRhoU_)
    State_V(BuffBx_:BuffBz_)= State_V(BuffBx_:BuffBz_)*No2Si_V(UnitB_)
    State_V(BuffP_) = State_V(BuffP_)*No2Si_V(UnitP_)

  contains
    !==========================================================================
    subroutine add_density_omega_cross_r

      ! Add Omega x R term. For IH Omega_D = (0,0,OmegaBody)
      use IH_BATL_lib, ONLY: Xyz_DGB, x_, y_
      use IH_ModPhysics, ONLY: OmegaBody
      !------------------------------------------------------------------------
      Momentum_D(x_) = Momentum_D(x_) &
           - Density*OmegaBody*Xyz_DGB(y_,i,j,k,iBlock)
      Momentum_D(y_)= Momentum_D(y_) &
           + Density*OmegaBody*Xyz_DGB(x_,i,j,k,iBlock)

    end subroutine add_density_omega_cross_r
    !==========================================================================
  end subroutine IH_get_for_gm
  !============================================================================
  subroutine IH_get_for_pt( &
       IsNew, NameVar, nVarIn, nDimIn, nPoint, Xyz_DI, Data_VI)

    ! Interpolate Data_VI from IH at the list of positions Xyz_DI
    ! required by PT

    logical,          intent(in):: IsNew   ! true for new point array
    character(len=*), intent(in):: NameVar ! List of variables
    integer,          intent(in):: nVarIn  ! Number of variables in Data_VI
    integer,          intent(in):: nDimIn  ! Dimensionality of positions
    integer,          intent(in):: nPoint  ! Number of points in Xyz_DI

    real, intent(in) :: Xyz_DI(nDimIn,nPoint)  ! Position vectors
    real, intent(out):: Data_VI(nVarIn,nPoint) ! Data array

    ! Optimize search by storing indexes and distances
    integer, allocatable, save:: iBlockCell_DI(:,:)
    real,    allocatable, save:: Dist_DI(:,:)
    !--------------------------------------------------------------------------
    call IH_get_point_data(iBlockCell_DI, Dist_DI, IsNew, &
         NameVar, nVarIn, nDimIn, nPoint, Xyz_DI, Data_VI, &
         DoSendAllVar=.true.)

  end subroutine IH_get_for_pt
  !============================================================================
  subroutine IH_put_from_pt( &
       NameVar, nVarData, nPoint, Data_VI, iPoint_I, Pos_DI)

    use IH_BATL_lib, ONLY: &
         nDim, nBlock, MaxBlock, Unused_B, nI, nJ, nK, Xyz_DGB, &
         iTest, jTest, kTest, iBlockTest
    use IH_ModPhysics, ONLY: &
         No2Si_V, Si2No_V, UnitX_, UnitRho_, UnitN_, UnitRhoU_, &
         UnitEnergyDens_, UnitT_
    use IH_ModGeometry, ONLY: Used_GB
    use IH_ModAdvance, ONLY: ExtraSource_ICB

    character(len=*), intent(inout):: NameVar  ! List of variables
    integer,          intent(inout):: nVarData ! Number of variables in Data_VI
    integer,          intent(inout):: nPoint   ! Number of points in Pos_DI

    real,    intent(in), optional:: Data_VI(:,:)    ! Recv data array
    integer, intent(in), optional:: iPoint_I(nPoint)! Order of data
    real, intent(out), allocatable, optional:: Pos_DI(:,:)  ! Position vectors

    ! For unit conversion
    real, allocatable, save:: Si2No_I(:)

    integer:: i, j, k, iBlock, iPoint, iVar

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'IH_put_from_pt'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTestMe)write(*,*) NameSub,' starting with present(Data_VI)=', &
         present(Data_VI)

    if(.not. present(Data_VI))then
       ! Provide BATSRUS grid points to PT

       nPoint = 0
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.Used_GB(i,j,k,iBlock)) CYCLE
             nPoint = nPoint + 1
          end do; end do; end do
       end do

       if(allocated(Pos_DI)) deallocate(Pos_DI)
       allocate(Pos_DI(nDim,nPoint))

       iPoint = 0
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE

          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.Used_GB(i,j,k,iBlock)) CYCLE
             iPoint = iPoint + 1
             Pos_DI(1:nDim,iPoint) = &
                  Xyz_DGB(1:nDim,i,j,k,iBlock)*No2Si_V(UnitX_)
          end do; end do; end do
       end do

       if(DoTestMe)write(*,*) NameSub,' finished setting positions'

       RETURN
    end if

    ! set source terms due to neutral charge exchange
    if(.not.allocated(Si2No_I))then
       ! Set units for density, momentum and energy source terms
       allocate(Si2No_I(nVarData))
       do iVar = 1, nVarData, 5
          Si2No_I(iVar) = Si2No_V(UnitRho_)/Si2No_V(UnitT_)/Si2No_V(UnitN_)
          Si2No_I(iVar+1:iVar+3) &
               = Si2No_V(UnitRhoU_)/Si2No_V(UnitT_)/Si2No_V(UnitN_)
          Si2No_I(iVar+4) &
               = Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_)/Si2No_V(UnitN_)
       end do
    end if

    if(.not.allocated(ExtraSource_ICB)) &
         allocate(ExtraSource_ICB(nVarData,nI,nJ,nK,MaxBlock))

    iPoint = 0
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE

       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(.not.Used_GB(i,j,k,iBlock)) CYCLE
          iPoint = iPoint + 1
          if(iPoint_I(iPoint) < 0) &
               call CON_stop(NameSub//': IH point is outside of PT domain')
          ExtraSource_ICB(:,i,j,k,iBlock) = Data_VI(:,iPoint_I(iPoint))*Si2No_I
       end do; end do; end do
    end do

    if(DoTestMe)write(*,*) NameSub,' finished with source=', &
         ExtraSource_ICB(:,iTest,jTest,kTest,iBlockTest)

  end subroutine IH_put_from_pt
  !============================================================================
  subroutine IH_get_for_pt_dt(DtSi)

    ! Calculate the global time step for PC

    use IH_ModMain, ONLY: Dt
    use IH_ModPhysics, ONLY: No2Si_V, UnitT_
    use IH_ModTimeStepControl, ONLY: set_global_timestep

    real, intent(out) ::  DtSi
    !--------------------------------------------------------------------------
    ! use -1.0 so that no limit is applied on Dt
    call set_global_timestep(TimeSimulationLimit=-1.0)
    DtSi = Dt*No2Si_V(UnitT_)

  end subroutine IH_get_for_pt_dt
  !============================================================================
  subroutine IH_get_for_sc( &
       IsNew, NameVar, nVarIn, nDimIn, nPoint, Xyz_DI, Data_VI)

    ! Interpolate Data_VI from EE at the list of positions Xyz_DI
    ! required by SC

    use IH_ModPhysics, ONLY: Si2No_V, UnitX_, No2Si_V, iUnitCons_V
    use IH_ModAdvance, ONLY: State_VGB, Bx_, Bz_
    use IH_ModVarIndexes, ONLY: nVar
    use IH_ModB0, ONLY: UseB0, get_b0
    use IH_BATL_lib, ONLY: iProc, nDim, MaxDim, MinIJK_D, MaxIJK_D, &
         find_grid_block
    use IH_ModIO, ONLY: iUnitOut
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

    integer:: iPoint, iBlock, iProcFound

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'IH_get_for_sc'
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
             call stop_mpi(NameSub//' could not find position on this proc')
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

       Data_VI(1:nVar,iPoint) = State_V*No2Si_V(iUnitCons_V)

    end do

  end subroutine IH_get_for_sc
  !============================================================================
  subroutine IH_get_ee_region( &
       NameVar, nVarData, nPoint, Pos_DI, Data_VI, iPoint_I)

    ! This routine is actually for EE-SC coupling

    ! This function will be called 3 times :
    !
    ! 1) Count grid cells to be overwritten by EE (except for extra variables)
    !
    ! 2) Return the Xyz_DGB coordinates of these cells
    !
    ! 3) Recieve Data_VI from SC and put them into State_VGB.
    !    The indexing array iPoint_I is needed to maintain the same order as
    !    the original position array Pos_DI was given in 2)

    use IH_BATL_lib, ONLY: Xyz_DGB, nBlock, Unused_B, &
         IsRLonLat, nI, nJ, nK, CoordMin_DB, CellSize_DB
    use IH_ModGeometry, ONLY: r_GB
    use IH_ModPhysics, ONLY: No2Si_V, UnitX_, Si2No_V, iUnitCons_V
    use IH_ModMain, ONLY: UseB0
    use IH_ModB0, ONLY: B0_DGB
    use IH_ModAdvance, ONLY: State_VGB, Bx_, Bz_
    use IH_ModMultiFluid, ONLY: nIonFluid
    use CON_coupler, ONLY: Grid_C, EE_, iVarTarget_V
    use ModNumConst, ONLY: cPi, cTwoPi

    character(len=*), intent(inout):: NameVar ! List of variables
    integer, intent(inout):: nVarData ! Number of variables in Data_VI
    integer, intent(inout):: nPoint   ! Number of points in Pos_DI
    real, intent(inout), allocatable, optional :: Pos_DI(:,:)  ! Positions

    real,    intent(in), optional:: Data_VI(:,:)    ! Recv data array
    integer, intent(in), optional:: iPoint_I(nPoint)! Order of data

    logical :: DoCountOnly
    integer :: i, j, k, iBlock, iPoint, iVarBuffer, iVar
    real    :: CoordMinEe_D(3), CoordMaxEe_D(3), Coord_D(3)

    character(len=*), parameter:: NameSub = 'IH_get_ee_region'
    !--------------------------------------------------------------------------
    if(.not.IsRLonLat) &
         call CON_stop(NameSub//' works for spherical grid only')

    DoCountOnly = nPoint < 1

    CoordMinEe_D(1) = Grid_C(EE_)%Coord1_I(1)
    CoordMinEe_D(2) = Grid_C(EE_)%Coord2_I(1)
    CoordMinEe_D(3) = Grid_C(EE_)%Coord3_I(1)
    CoordMaxEe_D(1) = Grid_C(EE_)%Coord1_I(2)
    CoordMaxEe_D(2) = Grid_C(EE_)%Coord2_I(2)
    CoordMaxEe_D(3) = Grid_C(EE_)%Coord3_I(2)

    ! Find ghost cells in the SC domain
    iPoint = 0
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       do k = 1, nK; do j = 1, nJ; do i = 1, nI

          ! Set generalized coordinates (longitude and latitude)
          Coord_D = &
               CoordMin_DB(:,iBlock) + ([i,j,k]-0.5)*CellSize_DB(:,iBlock)

          ! Overwrite first coordinate with true radius
          Coord_D(1) = r_GB(i,j,k,iBlock)

          ! Fix longitude if min longitude of EE is negative
          if(CoordMinEe_D(2) < 0.0 .and. Coord_D(2) > cPi) &
               Coord_D(2) = Coord_D(2) - cTwoPi

          ! Check if cell is inside EE domain
          if(any(Coord_D < CoordMinEe_D)) CYCLE
          if(any(Coord_D > CoordMaxEe_D)) CYCLE

          ! Found a point to be set by EE
          iPoint = iPoint + 1
          if(DoCountOnly) CYCLE

          if(present(Data_VI))then
             ! Put Data_VI obtained from EE into State_VGB
             ! Only a subset of variables are defined by EE
             do iVarBuffer = 1, nVarData
                iVar = iVarTarget_V(iVarBuffer)
                State_VGB(iVar,i,j,k,iBlock) = &
                     Data_VI(iVarBuffer,iPoint_I(iPoint)) &
                     *Si2No_V(iUnitCons_V(iVar))
             end do
             if(UseB0) State_VGB(Bx_:Bz_,i,j,k,iBlock) = &
                  State_VGB(Bx_:Bz_,i,j,k,iBlock) - B0_DGB(:,i,j,k,iBlock)
          else
             ! Provide position to EE
             Pos_DI(:,iPoint) = Xyz_DGB(:,i,j,k,iBlock)*No2Si_V(UnitX_)
          end if

       end do; end do; end do
    end do

    if(DoCountOnly) nPoint = iPoint

  end subroutine IH_get_ee_region
  !============================================================================
  subroutine IH_put_from_ee( &
       NameVar, nVarData, nPoint, Data_VI, iPoint_I, Pos_DI)

    ! This routine is actually for EE-SC coupling

    use IH_BATL_lib, ONLY: nDim

    character(len=*), intent(inout):: NameVar ! List of variables
    integer, intent(inout):: nVarData! Number of variables in Data_VI
    integer, intent(inout):: nPoint  ! Number of points in Pos_DI

    real, intent(in), optional:: Data_VI(:,:)    ! Recv data array
    integer, intent(in), optional:: iPoint_I(nPoint)! Order of data
    real, intent(out), allocatable, optional:: Pos_DI(:,:) ! Position vectors

    character(len=*), parameter:: NameSub = 'IH_put_from_ee'
    !--------------------------------------------------------------------------

    if(.not. present(Data_VI))then
       nPoint=0;
       ! get nPoint
       call IH_get_ee_region(NameVar, nVarData, nPoint, Pos_DI)

       if(allocated(Pos_DI)) deallocate(Pos_DI)
       allocate(Pos_DI(nDim,nPoint))

       ! get Pos_DI
       call IH_get_ee_region(NameVar, nVarData, nPoint, Pos_DI)

       RETURN
    end if

    ! set State variables
    call IH_get_ee_region(NameVar, nVarData, nPoint, Pos_DI, Data_VI, iPoint_I)

  end subroutine IH_put_from_ee
  !============================================================================
  subroutine IH_get_for_ee( &
       IsNew, NameVar, nVarIn, nDimIn, nPoint, Xyz_DI, Data_VI)

    ! This routine is actually for SC-EE coupling

    ! Interpolate Data_VI from SC at the list of positions Xyz_DI
    ! required by EE

    logical,          intent(in):: IsNew   ! true for new point array
    character(len=*), intent(in):: NameVar ! List of variables
    integer,          intent(in):: nVarIn  ! Number of variables in Data_VI
    integer,          intent(in):: nDimIn  ! Dimensionality of positions
    integer,          intent(in):: nPoint  ! Number of points in Xyz_DI

    real, intent(in) :: Xyz_DI(nDimIn,nPoint)  ! Position vectors
    real, intent(out):: Data_VI(nVarIn,nPoint) ! Data array

    ! Optimize search by storing indexes and distances
    integer, allocatable, save:: iBlockCell_DI(:,:)
    real,    allocatable, save:: Dist_DI(:,:)
    !--------------------------------------------------------------------------
    call IH_get_point_data(iBlockCell_DI, Dist_DI, IsNew, &
         NameVar, nVarIn, nDimIn, nPoint, Xyz_DI, Data_VI)

  end subroutine IH_get_for_ee
  !============================================================================
  integer function IH_n_particle(iBlockLocal)

    use IH_ModParticleFieldLine, ONLY: iKindReg
    use IH_BATL_lib, ONLY: Particle_I

    integer, intent(in) :: iBlockLocal
    !--------------------------------------------------------------------------
    IH_n_particle = Particle_I(iKindReg)%nParticle

  end function IH_n_particle
  !============================================================================
end module IH_wrapper
!==============================================================================
