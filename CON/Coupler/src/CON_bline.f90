!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module CON_bline
  ! allocate the grid used in this model
#ifdef OPENACC
  use ModUtilities, ONLY: norm2
#endif
  use ModUtilities,      ONLY: check_allocate, CON_stop
  use ModConst,          ONLY: cBoltzmann
  use CON_coupler,       ONLY: MaxComp, is_proc0, i_proc, &
       GridType, LocalGridType, set_standard_grid_descriptor, set_local_gd
  implicit none

  SAVE

  PRIVATE ! Except

  !=====The following public members are available at all PEs ================
  type(GridType),      public :: BL_Grid          ! Grid descriptor (global)
  type(LocalGridType), public :: BL_LocalGrid     ! Grid descriptor (local)
  integer, public :: BL_ =-1         ! ID of the target model (SP_, PT_)
  character(len=2):: NameCompBl = '' ! Name of the target model ('SP', 'PT'...)
  character(len=3)         :: TypeCoordBl  = '' ! Coord system of the Bl Model
  Character(len=3), public :: TypeCoordMh  = '' ! Coord system of the Mh Model
  logical, public :: UseBLine_C(MaxComp)=.false. ! To switch coupler for PT/SP
  ! Logical to determine, if a particular MH component is coupled to BL
  logical, public :: IsSource4BL_C(MaxComp) = .false.
  ! Boundaries of coupled domains in SC
  real,    public :: RScMin = 1000.0, RScMax = 0.0
  ! Boundaries of coupled domains in IH
  real,    public :: RIhMin = 1000.0, RIhMax = 0.0
  ! Boundaries of coupled domains in OH
  real,    public :: ROhMin = 1000.0, ROhMax = 0.0
  ! Transformation matrix
  real,    public :: BlMh_DD(3,3)
  !===== The rest is available on the BL_ processors =========================
  public :: BL_read_param
  public :: BL_set_grid
  public :: BL_init
  public :: BL_get_origin_points
  public :: BL_get_bounds
  public :: BL_adjust_lines
  public :: BL_update_r
  public :: BL_init_foot_points       ! initialize foot points
  public :: BL_n_particle             ! return number of points at the MF line
  public :: BL_put_from_mh            ! put MHD values from MH to SP
  public :: BL_interface_point_coords ! points rMinInterface<R<rMaxInterface
  public :: BL_put_line               ! points rMinBl < R < rMaxBl
  public :: BL_is_interface_block     ! to mark unusable lines
  public :: BL_set_line_foot_b
  public :: save_mhd

  real,    public  :: TimeBl = -1.0   ! Time of the model! Time of the model
  real             :: UnitXBl         ! Unit of scale (usually RSun)
  ! Store IDs of the lower and upper model
  integer, public  :: Lower_=0, Upper_=-1
  ! Total number of lines on given PE
  integer, public  :: nLine = -1 ! = ((iProc +1)*nLineAll)/nProc - iLineAll0
  ! Offset of B lines (per line)
  integer, allocatable, public :: iOffset_B(:)
  ! Number of variables in the state vector and the identifications
  integer, public, parameter ::  nDim = 3, nMHData = 13,      &
       LagrID_     = 0, & ! Lagrangian id           ^saved/   ^set to 0
       X_          = 1, & !                         |read in  |in copy_
       Y_          = 2, & ! Cartesian coordinates   |restart  |old_stat
       Z_          = 3, & !                         v/        |saved to
       Rho_        = 4, & ! Background plasma density         |mhd1
       T_          = 5, & ! Background temperature            |
       Ux_         = 6, & !                                   |may be
       Uy_         = 7, & ! Background plasma bulk velocity   |read from
       Uz_         = 8, & !                                   |mhd1
       Bx_         = 9, & !                                   |or
       By_         =10, & ! Background magnetic field         |received
       Bz_         =11, & !                                   |from
       Wave1_      =12, & !                                   |coupler
       Wave2_      =13, & ! Alfven wave turbulence            v
       R_          =14    ! Heliocentric distance
  ! variable names
  character(len=10), public, parameter:: NameVar_V(LagrID_:nMHData)&
        = ['LagrID    ', &
       'X         ', &
       'Y         ', &
       'Z         ', &
       'Rho       ', &
       'T         ', &
       'Ux        ', &
       'Uy        ', &
       'Uz        ', &
       'Bx        ', &
       'By        ', &
       'Bz        ', &
       'Wave1     ', &
       'Wave2     ']
  character(len=61) :: NamePlotVar = ''
  character(len=16) :: NameMHData = ''
  character(len=*), parameter:: NameExtension='.out'
  ! Local variables
  integer :: iProc = -1  ! Rank    of PE in communicator of BL_
  integer :: nProc = -1  ! Number of PEs in communicator of BL_
  !
  ! Grid integer parameters:
  ! Maximum number of vertexes per line
  integer :: nVertexMax

  ! angular grid at origin surface
  integer  :: nLon  = -1
  integer  :: nLat  = -1

  ! Total number of lines on all preceeding PEs
  integer  :: iLineAll0 = -1 ! = (iProc*nLineAll)/nProc

  ! Total number of lines on all PEs
  integer  :: nLineAll = -1

  ! MHD state vector is a pointer to be joined to a target array
  ! MHData_VIB(LagrID_:nMHData, 1:nVertexMax, 1:nLine)
  real, allocatable, target :: MHData_VIB(:,:,:)

  ! Aux state vector is a pointer to be joined to a target array
  ! State_VIB(R_:nVar, 1:nVertexMax, 1:nLine)
  ! nVar may be optionally read in BL_set_grid
  integer  :: nVar = R_
  real, allocatable, target :: State_VIB(:,:,:)

  ! nVertex_B is a pointer, which is joined to a target array
  ! For stand alone version the target array is allocated here
  ! nVertex_B(1:nLine)
  integer, allocatable, target :: nVertex_B(:)

  integer, parameter :: Length_=4
  real,    allocatable, target :: FootPoint_VB(:, :)

  ! To mark unusable lines (too short or those with negative VDF)
  logical, allocatable, target :: Used_B(:)

  ! coupling parameters:
  ! domain boundaries
  real :: rInterfaceMin, rInterfaceMax

  ! buffer boundaries located near lower (Lo) or upper (Up) boudanry of domain
  real :: rBufferLo, rBufferUp

  ! Allowed range of helocentric distances
  real :: rMinBl, rMaxBl

  ! Coefficient to transform energy
  real         :: EnergySi2Io = 1.0/cBoltzmann

contains
  !============================================================================
  subroutine BL_read_param(NameCommand,iError)
    use CON_coupler, ONLY: use_comp, i_comp, SP_, PT_, SC_, IH_, OH_, &
         i_proc, n_proc, is_proc, Couple_CC
    use ModReadParam, ONLY: read_var
    character(len=*), intent(in) :: NameCommand
    integer, intent(inout)       :: iError
    integer :: nSource, iSource, iComp, DnCouple
    real    :: DtCouple
    logical :: DoThis
    character(len=2) :: NameComp = ''
    character(len=*), parameter:: NameSub = 'BL_read_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case("#FIELDLINE")
       call read_var('NameTarget', NameCompBl)
       BL_ = i_comp(NameCompBl)
       if(.not.use_comp(BL_)) then
          if(is_proc0()) write(*,'(a)') NameSub//&
               ' SWMF_ERROR for NameMaster: '//&
               NameComp//' is OFF or not registered, not MFLAMPA target'
          iError = 34
          RETURN
       end if
       if(.not.(BL_ == PT_.or.BL_==SP_))then
          if(is_proc0()) write(*,'(a)') NameSub//&
               ' SWMF_ERROR for NameMaster: '//&
               'SP or PT can be BLine target, but '//NameComp//' cannot'
          iError = 34
          RETURN
       end if
       UseBLine_C = .false.; UseBLine_C(BL_) = .true.
       if(is_proc(BL_))then
          iProc = i_proc(BL_); nProc = n_proc(BL_)
       end if
       call read_var('nSource', nSource)
       if(nSource==0)RETURN ! Uncoupled SP
       if(nSource>3)then
          if(is_proc0()) write(*,'(a)') NameSub//&
               ' SWMF_ERROR for NameMaster: '//&
               'there can be only not more than 3 source components in MFLAMPA'
          iError = 34
          RETURN
       end if
       IsSource4BL_C = .false.
       do iSource = 1, nSource
          call read_var('NameSource', NameComp)
          iComp = i_comp(NameComp)
          select case(iComp)
          case(SC_)
             if(.not.use_comp(SC_)) then
                if(is_proc0()) write(*,'(a)') NameSub//&
                     ' SWMF_ERROR for NameMaster: '//&
                     ' SC  needed in MFLAMPA is OFF or not registered'
                iError = 34
                RETURN
             end if
             call read_var('RScMin', RScMin)
             call read_var('RScMax', RScMax)
             IsSource4BL_C(SC_) = .true.
          case(IH_)
             if(.not.use_comp(IH_)) then
                if(is_proc0()) write(*,'(a)') NameSub//&
                     ' SWMF_ERROR for NameMaster: '//&
                     ' IH  needed in MFLAMPA is OFF or not registered'
                iError = 34
                RETURN
             end if
             call read_var('RIhMin', RIhMin)
             call read_var('RIhMax', RIhMax)
             IsSource4BL_C(IH_) = .true.
          case(OH_)
             if(.not.use_comp(OH_)) then
                if(is_proc0()) write(*,'(a)') NameSub//&
                     ' SWMF_ERROR for NameMaster: '//&
                     ' OH  needed in MFLAMPA is OFF or not registered'
                iError = 34
                RETURN
             end if
             call read_var('ROhMin', ROhMin)
             call read_var('ROhMax', ROhMax)
             IsSource4BL_C(OH_) = .true.
          case default
             if(is_proc0()) write(*,*) NameSub//&
                  ' SWMF_ERROR for NameMaster: '//&
                  'SC, IH and OH are only allowed as MFLAMPA sources'
             iError = 34
             RETURN
          end select
       end do
       if(.not.is_proc(BL_))RETURN
       if(IsSource4BL_C(SC_))then
          Lower_ = SC_
          rMinBl   = RScMin
       elseif(IsSource4BL_C(IH_))then
          Lower_ = IH_
          rMinBl   = RIhMin
       else
          Lower_ = OH_
          rMinBl   = ROhMin
       end if
       if(IsSource4BL_C(OH_))then
          Upper_ = OH_
          rMaxBl   = ROhMax
       elseif(IsSource4BL_C(IH_))then
          Upper_ = IH_
          rMaxBl   = RIhMax
       else
          Upper_ = SC_
          rMaxBl   = RScMax
       end if
       select  case(BL_)
       case(PT_)
          NameMHData = 'PT/plots/MH_data'
       case(SP_)
          NameMHData = 'SP/IO2/MH_data'
       end select
    case("#COUPLEFIELDLINE")
       if(BL_ < 1)then
          if(is_proc0())write(*,'(a)') NameSub//&
            ' SWMF_ERROR for NameMaster: '//&
            '#FIELDLINE command must preceed '//NameCommand
          iError = 34
          RETURN
       end if
       call read_var('DnCouple', DnCouple)
       call read_var('DtCouple', DtCouple)
       DoThis = DnCouple > 0 .or. DtCouple > 0.0
       where(IsSource4BL_C (:))
          Couple_CC(:, BL_)%Dn     = DnCouple
          Couple_CC(:, BL_)%Dt     = DtCouple
          Couple_CC(:, BL_)%DoThis = DoThis
       end where
    case default
       if(is_proc0()) write(*,'(a)') NameSub//&
            ' SWMF_ERROR for NameMaster: '//&
            'Unknown command '//NameCommand
       iError = 34
       RETURN
    end select
  end subroutine BL_read_param
  !============================================================================
  subroutine BL_init(nParticleIn, nLonIn, nLatIn, &
       rPointer_VIB,  nPointer_B,  nVarIn, StateIO_VIB, FootPointIn_VB, &
       lPointer_B)
    use ModUtilities, ONLY: join_string
    real,    intent(inout), pointer :: rPointer_VIB(:,:,:)
    integer, intent(inout), pointer :: nPointer_B(:)
    integer, intent(in)             :: nParticleIn, nLonIn, nLatIn
    integer, optional, intent(in)   :: nVarIn
    real,    intent(inout), optional, pointer :: &
         StateIO_VIB(:,:,:), FootPointIn_VB(:,:)
    logical, intent(inout), optional, pointer :: lPointer_B(:)
    !
    ! Loop variable

    integer :: iVertex
    integer:: iError
    character(len=*), parameter:: NameSub = 'BL_init'
    !--------------------------------------------------------------------------
    nVertexMax            = nParticleIn
    nLon = nLonIn; nLat = nLatIn; nLineAll = nLon*nLat
    iLineAll0 = (iProc*nLineAll)/nProc
    nLine = ((iProc+1)*nLineAll) / nProc - iLineAll0

    allocate(MHData_VIB(LagrID_:nMHData, 1:nVertexMax, 1:nLine), &
         stat=iError)
    call check_allocate(iError, NameSub//'_MHData_VIB')
    rPointer_VIB => MHData_VIB
    MHData_VIB = 0.0
    !
    ! reset lagrangian ids
    !
    do iVertex = 1, nVertexMax
       MHData_VIB(LagrID_, iVertex, 1:nLine) = real(iVertex)
    end do

    ! Set nVertex_B array

    allocate(nVertex_B(1:nLine), &
         stat=iError)
    call check_allocate(iError, NameSub//'_nVertex_B')
    nPointer_B => nVertex_B
    nVertex_B = 0

    ! Set auxiliary State_VIB array

    if(present(nVarIn))nVar = nVarIn
    allocate(State_VIB(R_:nVar, 1:nVertexMax, 1:nLine), &
         stat=iError)
    call check_allocate(iError, NameSub//'_State_VIB')
    if(present(StateIO_VIB))StateIO_VIB => State_VIB
    State_VIB = -1

    ! Set FootPoint_VB array

    allocate(FootPoint_VB(LagrID_:Length_, nLine), stat=iError)
    call check_allocate(iError, NameSub//'_FootPoint_VB')
    if(present(FootPointIn_VB))FootPointIn_VB => FootPoint_VB
    FootPoint_VB = -1

    ! Set iOffset_B array

    allocate(iOffset_B(nLine)); iOffset_B = 0

    ! set array to mark unusable line
    allocate(Used_B(nLine)); Used_B = .true.
    if(present(lPointer_B))lPointer_B => Used_B

    ! Set plot variable names

    call join_string(NameVar_V(LagrID_:nMHData), NamePlotVar)
    NamePlotVar = trim(NamePlotVar)//' LagrID x y z'
    if(iProc==0)write(*,'(a)')NamePlotVar
  end subroutine BL_init
  !============================================================================
  subroutine BL_get_origin_points(ROrigin, LonMin, LonMax, LatMin, LatMax)
    use ModCoordTransform, ONLY: rlonlat_to_xyz
    !
    ! Parameters of the grid of the bline origin points at the origin surface
    ! at R = ROrigin
    real, intent(in)  :: ROrigin, LonMin, LonMax, LatMin, LatMax

    ! line (line) number  and corresponing lon-lat indexes
    integer           :: iLat, iLon, iLine
    !
    ! Sell size on the origin surface, per line
    real              ::  DLon, DLat
    !--------------------------------------------------------------------------
    !
    ! angular grid's steps
    DLon = (LonMax - LonMin)/nLon;  DLat = (LatMax - LatMin)/nLat
    nVertex_B(1:nLine) = 1   ! One origin point per line
    do iLine = 1, nLine
       call BL_iblock_to_lon_lat(iLine, iLon, iLat)
       call rlonlat_to_xyz([ROrigin, LonMin + (iLon - 0.5)*DLon, &
            LatMin + (iLat - 0.5)*DLat], MHData_VIB(X_:Z_,1,iLine))
    end do
  end subroutine BL_get_origin_points
  !============================================================================
  subroutine BL_iblock_to_lon_lat(iBlockIn, iLonOut, iLatOut)
    ! return angular grid's indexes corresponding to this line
    integer, intent(in) :: iBlockIn
    integer, intent(out):: iLonOut
    integer, intent(out):: iLatOut
    integer :: iLineAll
    !--------------------------------------------------------------------------
    ! Get node number from line number
    iLineAll = iBlockIn + iLineAll0
    iLatOut = 1 + (iLineAll - 1)/nLon
    iLonOut = iLineAll - nLon*(iLatOut - 1)
  end subroutine BL_iblock_to_lon_lat
  !============================================================================
  subroutine BL_set_grid(TypeCoordSystem, UnitX, EnergySi2IoIn)
    use CON_coupler, ONLY: set_coord_system, is_proc, &
         init_decomposition, get_root_decomposition, bcast_decomposition
    character(len=3), intent(in)    :: TypeCoordSystem
    real,    intent(in)             :: UnitX
    real,    optional,   intent(in) :: EnergySi2IoIn
    !
    ! Proc_ and Block_ number, for a given node:
    !
    integer, parameter:: &
         Proc_  = 1, & ! Processor that has this line/node
         Block_ = 2    ! line that has this line/node
    !
    ! They is the first index values for
    ! the following array
    integer, allocatable, target :: iGridGlobal_IA(:,:)

    character(len=*), parameter:: NameVarCouple =&
         'rho p mx my mz bx by bz i01 i02 pe'
    logical, save:: IsInitialized = .false.
    !--------------------------------------------------------------------------
    if(IsInitialized)RETURN
    IsInitialized = .true.
    ! Initialize 3D grid with NON-TREE structure
    call init_decomposition(GridID_ = BL_, CompID_ = BL_, nDim = nDim)
    ! Construct decomposition

    if(is_proc0(BL_))then
       call init_indexes
       call get_root_decomposition(&
            GridID_       = BL_,&
            iRootMapDim_D = [1, nLon, nLat],&
            CoordMin_D    = [0.50, 0.50, 0.50],&
            CoordMax_D    = [nVertexMax, nLon, nLat] + 0.50,&
            nCell_D       = [nVertexMax, 1, 1],&
            PE_I          = iGridGlobal_IA(Proc_,:),&
            iBlock_I      = iGridGlobal_IA(Block_,:))
       deallocate(iGridGlobal_IA)
    end if
    call bcast_decomposition(BL_)
    ! Coordinate system is Heliographic Inertial Coordinate System (HGI)
    ! with length measured in solar radii
    call set_coord_system(&
         GridID_      = BL_            , &
         TypeCoord    = TypeCoordSystem, &
         TypeGeometry = 'cartesian'    , &
         NameVar      = NameVarCouple  , &
         UnitX        = UnitX)
    call set_standard_grid_descriptor(BL_, Grid=BL_Grid)
    if(.not.is_proc(BL_))RETURN
    UnitXBl = UnitX
    TypeCoordBl = TypeCoordSystem
    call set_local_gd(iProc = i_proc(), Grid=BL_Grid, LocalGrid=BL_LocalGrid)
    if(present(EnergySi2IoIn))EnergySi2Io = EnergySi2IoIn
  contains
    !==========================================================================
    subroutine init_indexes
      integer:: iLat, iLon, iLineAll, iLine, iProcNode

      !
      ! fill grid containers
      !
      integer:: iError
      character(len=*), parameter:: NameSub = 'init_indexes'
      !------------------------------------------------------------------------
      allocate(iGridGlobal_IA(Proc_:Block_, nLineAll), &
           stat=iError)
      call check_allocate(iError, NameSub//'_iGridGlobal_IA')
      iLine = 1; iProcNode = 0
      do iLat = 1, nLat
         do iLon = 1, nLon
            iLineAll = iLon + nLon * (iLat-1)
            !
            ! iProcNode = ceiling(real(iLineAll*nProc)/nLineAll) - 1
            !
            iGridGlobal_IA(Proc_,   iLineAll)  = iProcNode
            iGridGlobal_IA(Block_,  iLineAll)  = iLine
            if(iLineAll == ((iProcNode+1)*nLineAll)/nProc)then
               !
               ! This was the last node on the iProcNode
               !
               iLine = 1; iProcNode = iProcNode + 1
            else
               iLine = iLine + 1
            end if
         end do
      end do
    end subroutine init_indexes
    !==========================================================================
  end subroutine BL_set_grid
  !============================================================================
  subroutine BL_get_bounds(rMinIn, rMaxIn, &
       rBufferLoIn, rBufferUpIn)
    real,           intent(in) :: rMinIn, rMaxIn
    real, optional, intent(in) :: rBufferLoIn, rBufferUpIn
    ! set domain boundaries
    !--------------------------------------------------------------------------
    rInterfaceMin = rMinIn; rInterfaceMax = rMaxIn
    ! set buffer boundaries
    if(present(rBufferLoIn))then
       rBufferLo = rBufferLoIn
    else
       rBufferLo = rMinIn
    end if
    if(present(rBufferUpIn))then
       rBufferUp = rBufferUpIn
    else
       rBufferUp = rMaxIn
    end if
  end subroutine BL_get_bounds
  !============================================================================
  integer function BL_n_particle(iBlockLocal)
    integer, intent(in) :: iBlockLocal
    !--------------------------------------------------------------------------
    BL_n_particle = nVertex_B(  iBlockLocal)
  end function BL_n_particle
  !============================================================================
  ! Put MHD data in a single point
  subroutine BL_put_from_mh(nPartial, iPutStart, Put, W, DoAdd, Buff_I, nVar)
    use CON_axes,    ONLY: transform_velocity
    use CON_coupler, ONLY: iVar_V, DoCoupleVar_V   ,&
         Density_, RhoCouple_, Pressure_, PCouple_ ,&
         Momentum_, RhoUxCouple_, RhoUzCouple_     ,&
         BField_, BxCouple_, BzCouple_             ,&
         Wave_, WaveFirstCouple_, WaveLastCouple_
    use CON_router,  ONLY: IndexPtrType, WeightPtrType
    use ModConst,    ONLY: cProtonMass
    integer,intent(in)            :: nPartial, iPutStart, nVar
    type(IndexPtrType), intent(in):: Put
    type(WeightPtrType),intent(in):: W
    logical,            intent(in):: DoAdd
    real,               intent(in):: Buff_I(nVar)
    integer:: i, iLine, iRho
    real   :: Weight,  R, Aux, State_V(Rho_:Wave2_), XyzBl_D(X_:Z_)
    ! cell and line indices
    character(len=*), parameter:: NameSub = 'BL_put_from_mh'
    !--------------------------------------------------------------------------
    i       = Put%iCB_II(1, iPutStart)
    iLine   = Put%iCB_II(4, iPutStart)
    ! Location:
    XyzBl_D = MHData_VIB(X_:Z_,i,iLine)
    iRho    = iVar_V(RhoCouple_)
    ! Copy state variables from buffer. Convert  units if needed.
    State_V = 0.0 ; State_V(Rho_) = Buff_I(iRho)/cProtonMass ! a. m. u. per m3
    ! perform vector transformation from the source model to the BL one
    if(DoCoupleVar_V(BField_))State_V(Bx_:Bz_) =     &
         matmul(BlMh_DD, Buff_I(iVar_V(BxCouple_):iVar_V(BzCouple_))) ! Tesla
    if(DoCoupleVar_V(Momentum_))State_V(Ux_:Uz_) =               &! in BL model
         transform_velocity(TimeSim = TimeBl                    ,&! s
         v1_D = Buff_I(iVar_V(RhoUxCouple_):iVar_V(RhoUzCouple_))&! in MH Model
         /Buff_I(iRho)                                          ,&! in m/s
         Xyz1_D = matmul(XyzBl_D, BlMh_DD)*UnitXBl              ,&! in MH Model
         NameCoord1 = TypeCoordMh                               ,&! in MH Model
         NameCoord2 = TypeCoordBl)                                ! in BL model
    if(DoCoupleVar_V(Pressure_))State_V(T_) =        &
         Buff_I(iVar_V(PCouple_))*(cProtonMass/Buff_I(iRho))*EnergySi2Io ! K
    if(DoCoupleVar_V(Wave_))State_V(Wave1_:Wave2_) = &
         Buff_I(iVar_V(WaveFirstCouple_):iVar_V(WaveLastCouple_))        ! J/m3
    ! interpolation weight: if the point is within the buffer, the state vector
    ! is interpolated between those in the components
    Aux = 0; if(DoAdd)Aux = 1.0; Weight = 1.0; R = norm2(XyzBl_D)
    if(R >= rInterfaceMin .and. R < rBufferLo)then
       Aux = 1.0
       Weight = Weight*(0.50 + 0.50*tanh(2*(2*R - &
            RBufferLo - RInterfaceMin)/(RBufferLo - RInterfaceMin)))
    elseif(R >= rBufferUp .and. R < rInterfaceMax)then
       Weight = Weight*(0.50 - 0.50*tanh(2*(2*R - &
            RInterfaceMax - RBufferUp)/(RInterfaceMax - RBufferUp)))
    end if
    ! put the data
    MHData_VIB(Rho_:Wave2_,i,iLine) = Aux*MHData_VIB(Rho_:Wave2_,i,iLine) &
         + State_V*Weight
  end subroutine BL_put_from_mh
  !============================================================================
  subroutine BL_interface_point_coords(nDim, Xyz_D, &
       nIndex, iIndex_I, IsInterfacePoint)
    ! interface points, which need to be communicated to other components to
    ! perform field line  extraction and further coupling with SP
    integer,intent(in)   :: nDim
    real,   intent(inout):: Xyz_D(nDim)
    integer,intent(in)   :: nIndex
    integer,intent(inout):: iIndex_I(nIndex)
    logical,intent(out)  :: IsInterfacePoint
    integer  :: iVertex, iLine
    real   :: R
    character(len=*), parameter:: NameSub = 'BL_interface_point_coords'
    !--------------------------------------------------------------------------
    iVertex = iIndex_I(1); iLine    = iIndex_I(4)
    ! Fix coordinates to be used in mapping
    Xyz_D = MHData_VIB(X_:Z_, iVertex, iLine)
    ! Check whether the particle is within interface bounds
    R = norm2(Xyz_D)
    IsInterfacePoint = R >= rInterfaceMin .and. R < rInterfaceMax
  end subroutine BL_interface_point_coords
  !============================================================================
  subroutine BL_put_line(nPartial, iPutStart, Put,&
       Weight, DoAdd, XyzBl_D, nDimBl)
    use CON_router, ONLY: IndexPtrType, WeightPtrType
    integer, intent(in) :: nPartial, iPutStart, nDimBl
    type(IndexPtrType), intent(in) :: Put
    type(WeightPtrType),intent(in) :: Weight
    logical,            intent(in) :: DoAdd
    real,               intent(in) :: XyzBl_D(nDimBl)

    ! indices of the particle
    integer :: iLine, iVertex
    ! Misc
    real :: R
    character(len=*), parameter:: NameSub = 'BL_put_line'
    !--------------------------------------------------------------------------
    iLine    = Put%iCB_II(4,iPutStart)
    if(.not.Used_B(iLine))RETURN
    iVertex  = Put%iCB_II(1,iPutStart) + iOffset_B(iLine)
    if(iVertex > nVertexMax)RETURN
    R = norm2(XyzBl_D)
    if(R < rMinBl .or. R >= rMaxBl)then
       ! Sort out particles left the SP domain
       MHData_VIB(X_:Z_,iVertex, iLine) = 0.0
    else
       ! store passed particles coordinates
       MHData_VIB(X_:Z_,iVertex,iLine) = XyzBl_D
       nVertex_B(iLine) = MAX(nVertex_B(iLine), iVertex)
    end if
  end subroutine BL_put_line
  !============================================================================
  subroutine BL_init_foot_points
    integer:: iOffset, iLine, iBegin,  iEnd ! loop variables
    character(len=*), parameter:: NameSub = 'BL_init_foot_points'
    !--------------------------------------------------------------------------
    do iLine = 1, nLine
       iBegin = 1
       do while (all(MHData_VIB(X_:Z_,iBegin,iLine)==0.0))
          iBegin = iBegin + 1
       end do
       if(iBegin /= 1)then
          ! Offset particle arrays
          iEnd   = nVertex_B(iLine)
          iOffset = 1 - iBegin
          MHData_VIB(LagrID_:Z_,1:iEnd+iOffset,iLine) = &
               MHData_VIB(LagrID_:Z_,iBegin:iEnd,iLine)
          nVertex_B(iLine) = nVertex_B(iLine) + iOffset
       end if
       call BL_set_line_foot_b(iLine)
    end do
  end subroutine BL_init_foot_points
  !============================================================================
  ! Called from coupler after the updated grid point lo<cation are
  ! received from the other component (SC, IH). Determines whether some
  ! grid points should be added/deleted
  subroutine BL_adjust_lines(Source_)
    integer, intent(in) :: Source_
    ! once new geometry of lines has been put, account for some particles
    ! exiting the domain (can happen both at the beginning and the end)
    integer:: iVertex, iLine, iBegin,  iEnd ! loop variables

    logical:: DoAdjustLo, DoAdjustUp
    logical:: IsMissing

    real              :: R

    integer, parameter:: Lo_ = 1, Up_ = 2
    integer, parameter:: iIncrement_II(2,2) =reshape([1,0,0,-1],[2,2])
    integer:: iParticle_I(2), iLoop
    ! once new geometry of lines has been put, account for some particles
    ! exiting the domain (can happen both at the beginning and the end)

    ! Called after the grid points are received from the
    ! component, nullify offset
    character(len=100):: StringError

    character(len=*), parameter:: NameSub = 'BL_adjust_lines'
    !--------------------------------------------------------------------------
    DoAdjustLo = Source_ == Lower_
    DoAdjustUp = Source_ == Upper_
    line:do iLine = 1, nLine
       if(.not. Used_B(iLine))CYCLE line
       if(nVertex_B(iLine) < 10)then
          write(*,*)NameSub//':iProc in BL=', iProc, ' iLine=',iLine
          write(*,*)NameSub//':RInterface Min, Max=', RInterfaceMin, &
               RInterfaceMax
          write(*,*)NameSub//':nVertex_B', nVertex_B(iLine)
          do iVertex = 1, nVertex_B(iLine)
             write(*,*)iVertex, MHData_VIB(LagrID_:Z_,iVertex,iLine), &
                  State_VIB(R_,iVertex,iLine)
          end do
          write(*,*)NameSub//':Too short line deleted'
          ! Remove too short lines
          Used_B(iLine) = .false.
          nVertex_B(iLine)=0
          CYCLE line
       end if
       iBegin = 1
       if(DoAdjustLo)then
          iOffset_B(iLine) = 0
          do while(all(MHData_VIB(X_:Z_,iBegin,iLine)==0.0))
             iBegin = iBegin + 1
          end do
       end if
       if(DoAdjustUp)then
          do while(all(MHData_VIB(X_:Z_,nVertex_B(iLine),iLine)==0.0))
             nVertex_B(iLine) = nVertex_B(iLine) - 1
          end do
       end if
       iEnd   = nVertex_B(iLine)
       iParticle_I(:) = [iBegin, iEnd]
       if(DoAdjustUp) then
          iLoop = Up_     ! Vertexes are checked from the upper end
       else
          iLoop = Lo_     ! Vertexes are checked from the lower end
       end if
       PARTICLE: do while(iParticle_I(1) < iParticle_I(2))
          iVertex = iParticle_I(iLoop) ! iBegin for SC, iEnd for IH
          R = State_VIB(R_,iVertex,iLine)
          ! account for all missing partiles along the line;
          ! --------------------------------------------------------
          ! particle import MUST be performed in order from lower to upper
          ! components (w/respect to radius), e.g. SC, then IH, then OH
          ! adjustment are made after each import;
          ! --------------------------------------------------------
          ! lines are assumed to start in lowest model
          ! lowest model should NOT have the Lo buffer,
          ! while the upper most should NOT have Up buffer,
          ! buffer are REQUIRED between models
          ! --------------------------------------------------------
          ! adjust as follows:
          ! LOWEST model may loose particles at the start of lines
          ! and (if line ends in the model) in the tail;
          ! loop over particles Lo-2-Up and mark losses (increase iBegin)
          ! until particles reach UP buffer;
          ! if line ends in the model, loop over particles Up-2-Lo
          ! and mark losses (decrease nVertex_B) until UP buffer is reached;
          ! after that, if line reenters the model,
          ! particles within the model may not be lost,
          !
          ! MIDDLE models are not allowed to loose particles
          !
          ! HIGHEST model may only loose particles in the tail;
          ! loop Up-2-Lo and mark losses  (decrease nVertex_B)
          ! until the bottom of Lo buffer is reaced
          ! --------------------------------------------------------
          ! whenever a particle is lost in lower models -> ERROR

          ! when looping Up-2-Lo and particle is in other model ->
          ! adjustments are no longer allowed
          ! if(iLoop == Up_ .and. (&
          !     R <  RInterfaceMin&
          !     .or.&
          !     R >= RInterfaceMax)&
          !     )then
          !   DoAdjustLo = .false.
          !   DoAdjustUp = .false.
          ! end if

          ! determine whether particle is missing
          IsMissing = all(MHData_VIB(X_:Z_,iVertex,iLine)==0.0)
          if(.not.IsMissing) then
             iParticle_I = iParticle_I + iIncrement_II(:,iLoop)
             CYCLE PARTICLE
          end if
          ! if need to adjust lower, but not upper boundary -> ADJUST
          if(DoAdjustLo .and. R < 1.10*rMinBl)then
             ! push iBegin in front of current particle;
             ! it will be pushed until it finds a non-missing particle
             iBegin = iVertex + 1
             iParticle_I = iParticle_I + iIncrement_II(:,Lo_)
             CYCLE PARTICLE
          end if
          ! if need to adjust upper, but not lower boundary -> ADJUST
          if(DoAdjustUp .and. R > 0.90*rMaxBl)then
             ! push nVertex_B() below current particle;
             ! it will be pushed until it finds a non-missing particle
             nVertex_B(iLine) = iVertex - 1
             iParticle_I = iParticle_I + iIncrement_II(:,Up_)
             CYCLE PARTICLE
          end if
          ! missing point in the lower part of the domain -> ERROR
          if(R < RInterfaceMin)then
             if(Source_==Lower_)write(*,*)'Is Lower Model'
             if(Source_==Upper_)write(*,*)'Is Upper Model'
             write(*,*)'iProc in BL=', iProc
             write(*,*)'RInterface Min, Max=', RInterfaceMin, RInterfaceMax
             write(*,*)'iVertex, R=', iVertex, R
             write(*,*)'iBegin, iEnd',  iBegin, iEnd
             do iVertex = iBegin,iEnd
                write(*,*)iVertex, MHData_VIB(LagrID_:Z_,iVertex,iLine), &
                     State_VIB(R_,iVertex,iLine)
             end do
             write(*,'(a)')NameSub//": particle has been lost"
             Used_B(iLine)  = .false.
             nVertex_B(iLine) = 0
             CYCLE line
          end if
          ! missing point in the upper part of the domain -> IGNORE;
          ! if needed to adjust beginning, then it is done,
          ! switch left -> right end of range and start adjusting
          ! tail of the line, if it has reentered current part of the domain
          if(R >= RInterfaceMax)then
             ! if(iLoop == Lo_)&
             !     iLoop = Up_
             ! iParticle_I = iParticle_I + iIncrement_II(:,iLoop)
             ! if(DoAdjustLo)then
             !   DoAdjustLo = .false.
             !   DoAdjustUp = .true.
             ! end if
             iParticle_I = iParticle_I + iIncrement_II(:,iLoop)
             CYCLE PARTICLE
          end if

          ! if point used to be in a upper buffer -> IGNORE
          ! if(R >= rBufferUp .and. R < rInterfaceMax)then
          !   iParticle_I = iParticle_I + iIncrement_II(:,iLoop)
          !   CYCLE PARTICLE
          ! end if

          ! If this is the upper model and the missing point is not
          ! near the upper boundary, nothing can be done
          if(DoAdjustUp)then
             write(*,*)'In the Upper Model'
             write(*,*)'iProc in BL=', iProc
             write(*,*)'RInterface Min, Max=', RInterfaceMin, RInterfaceMax
             write(*,*)'iVertex, R=', iVertex, R
             write(*,*)'iBegin, iEnd',  iBegin, iEnd
             do iVertex = iBegin,iEnd
                write(*,*)iVertex, MHData_VIB(LagrID_:Z_,iVertex,iLine), &
                     State_VIB(R_,iVertex,iLine)
             end do
             write(*,'(a)')NameSub//": particle has been lost"
             Used_B(iLine)  = .false.
             nVertex_B(iLine) = 0
             CYCLE line
          end if

          iParticle_I = iParticle_I + iIncrement_II(:,iLoop)
       end do PARTICLE

       ! DoAdjustLo = Source_ == Lower_
       ! DoAdjustUp = Source_ == Upper_
       if(iBegin/=1)then
          ! Offset particle arrays
          iEnd   = nVertex_B(iLine)
          iOffset_B(iLine) = 1 - iBegin
          MHData_VIB(LagrID_:Z_, 1:iEnd+iOffset_B(iLine), iLine) = &
               MHData_VIB(LagrID_:Z_, iBegin:iEnd, iLine)
          State_VIB(R_         , 1:iEnd+iOffset_B(iLine), iLine) = &
               State_VIB(R_         , iBegin:iEnd, iLine)
          nVertex_B(iLine) = nVertex_B(iLine) + iOffset_B(iLine)
          ! need to recalculate footpoints
          call BL_set_line_foot_b(iLine)
       end if
    end do line
    ! may need to add particles to the beginning of lines
    if(DoAdjustLo) call append_particles
  contains
    !==========================================================================
    subroutine append_particles
      ! appends a new particle at the beginning of lines if necessary
      integer:: iLine
      real:: DistanceToMin
      real, parameter:: cTol = 1E-06

      character(len=*), parameter:: NameSub = 'append_particles'
      !------------------------------------------------------------------------
      line:do iLine = 1, nLine
         if(.not.Used_B(iLine))CYCLE line
         ! check current value of offset: if not zero, adjustments have just
         ! been made, no need to append new particles
         if(iOffset_B(iLine) /= 0 )CYCLE line
         ! check if the beginning of the line moved far enough from its
         ! footprint on the solar surface
         DistanceToMin = norm2(&
              MHData_VIB(X_:Z_,1,iLine) - FootPoint_VB(X_:Z_,iLine))
         ! skip the line if it's still close to the Sun
         if(DistanceToMin*(1.0 + cTol) < FootPoint_VB(Length_, iLine))&
              CYCLE line
         ! append a new particle
         ! check if have enough space
         if(nVertexMax == nVertex_B( iLine))call CON_Stop(NameSub//&
              ': not enough memory allocated to append a new particle')
         ! Particles ID as handled by other components keep unchanged
         ! while their order numbers in SP are increased by 1. Therefore,
         iOffset_B(iLine)  = 1
         MHData_VIB(       LagrID_:Z_,2:nVertex_B(iLine) + 1, iLine)&
              = MHData_VIB(LagrID_:Z_,1:nVertex_B(iLine),     iLine)
         State_VIB(       R_        ,2:nVertex_B(iLine) + 1, iLine)&
              = State_VIB(R_        ,1:nVertex_B(iLine),     iLine)
         nVertex_B(iLine) = nVertex_B(iLine) + 1
         ! put the new particle just above the lower boundary
         MHData_VIB(X_:Z_,  1, iLine) = &
              FootPoint_VB(X_:Z_, iLine)*(1.0 + cTol)
         State_VIB(R_,          1, iLine) = &
              sqrt(sum((MHData_VIB(X_:Z_,  1, iLine))**2))
         MHData_VIB(LagrID_,1, iLine) = MHData_VIB(LagrID_, 2, iLine) - 1.0
         FootPoint_VB(LagrID_,iLine) = MHData_VIB(LagrID_, 1, iLine) - 1.0
      end do line
    end subroutine append_particles
    !==========================================================================
  end subroutine BL_adjust_lines
  !============================================================================
  subroutine BL_update_r
    ! Loop  variables:
    integer :: iVertex, iLine
    !--------------------------------------------------------------------------
    do iLine = 1, nLine
       do iVertex = 1, nVertex_B(iLine)
          State_VIB(R_, iVertex, iLine) = &
               norm2(MHData_VIB(X_:Z_, iVertex, iLine))
       end do
    end do
  end subroutine BL_update_r
  !============================================================================
  subroutine BL_set_line_foot_b(iLine)
    integer, intent(in) :: iLine

    ! existing particle with lowest index along line
    real:: Xyz1_D(X_:Z_)
    ! direction of the field at Xyz1_D and segment vectors between particles
    real, dimension(X_:Z_):: Dir0_D, Dist1_D, Dist2_D
    ! dot product Xyz1 and Dir1 and its sign
    real:: Dot, S
    ! distances between particles
    real:: Dist1, Dist2
    ! variable to compute coords of the footprints
    real:: Alpha
    !--------------------------------------------------------------------------
    ! get the coordinates of lower particle
    Xyz1_D = MHData_VIB(X_:Z_, 1, iLine)

    ! generally, field direction isn't known
    ! approximate it using directions of first 2 segments of the line
    Dist1_D = MHData_VIB(X_:Z_, 1, iLine) - &
         MHData_VIB(X_:Z_, 2, iLine)
    Dist1 = sqrt(sum(Dist1_D**2))
    Dist2_D = MHData_VIB(X_:Z_, 2, iLine) - &
         MHData_VIB(X_:Z_, 3, iLine)
    Dist2 = sqrt(sum(Dist2_D**2))
    Dir0_D = ((2*Dist1 + Dist2)*Dist1_D - Dist1*Dist2_D)/(Dist1 + Dist2)

    Dir0_D = Dir0_D/sqrt(sum(Dir0_D**2))

    ! dot product and sign: used in computation below
    Dot = sum(Dir0_D*Xyz1_D)
    S   = sign(1.0, Dot)

    ! there are 2 possible failures of the algorithm:
    ! Failure (1):
    ! no intersection of smoothly extended line with the sphere R = rMinBl
    if(Dot**2 - sum(Xyz1_D**2) + rMinBl**2 < 0)then
       ! project first particle for new footpoint
       FootPoint_VB(X_:Z_,iLine) = Xyz1_D * rMinBl / sqrt(sum(Xyz1_D**2))
    else
       ! Xyz0, the footprint, is distance Alpha away from Xyz1:
       ! Xyz0 = Xyz1 + Alpha * Dir0 and R0 = rMinBl =>
       Alpha = S * sqrt(Dot**2 - sum(Xyz1_D**2) + rMinBl**2) - Dot
       ! Failure (2):
       ! intersection is too far from the current beginning of the line,
       ! use distance between 2nd and 3rd particles on the line as measure
       if(abs(Alpha) > Dist2)then
          ! project first particle for new footpoint
          FootPoint_VB(X_:Z_,iLine) = Xyz1_D * rMinBl / sqrt(sum(Xyz1_D**2))
       else
          ! store newly found footpoint of the line
          FootPoint_VB(X_:Z_,iLine) = Xyz1_D + Alpha * Dir0_D
       end if
    end if

    ! length is used to decide when need to append new particles:
    ! use distance between 2nd and 3rd particles on the line
    FootPoint_VB(Length_,    iLine) = Dist2
    FootPoint_VB(LagrID_,    iLine) = MHData_VIB(LagrID_,1,iLine) - 1.0
  end subroutine BL_set_line_foot_b
  !============================================================================
  subroutine make_file_name(Time, iLine, NameOut)
    ! creates a string with file name and stores in NameOut;
    ! result is as follows:
    !   StringBase_Lon=?.?_Lat=?.?][_t?].NameExtension
    real,                 intent(in) :: Time
    integer,              intent(in) :: iLine
    character(len=100),   intent(out):: NameOut

    ! timetag
    character(len=8):: StringTime
    ! lon, lat indexes corresponding to iLineAll
    integer:: iLon, iLat
    !--------------------------------------------------------------------------
    write(NameOut,'(a)')trim(NameMHData)

    call BL_iblock_to_lon_lat(iLine, iLon, iLat)
    write(NameOut,'(a,i3.3,a,i3.3)') &
         trim(NameOut)//'_',iLon,'_',iLat

    call get_time_string(Time, StringTime)
    write(NameOut,'(a,i6.6,a)')  &
         trim(NameOut)//'_t'//StringTime
    write(NameOut,'(a)') trim(NameOut)//trim(NameExtension)
  end subroutine make_file_name
  !============================================================================
  subroutine get_time_string(Time, StringTime)
    ! the subroutine converts real variable Time into a string,
    ! the structure of the string is 'ddhhmmss',
    ! i.e shows number of days, hours, minutes and seconds
    ! after the beginning of the simulation
    real,             intent(in) :: Time
    character(len=8), intent(out):: StringTime

    ! This is the value if the time is too large

    !--------------------------------------------------------------------------
    StringTime = '99999999'
    if(Time < 100.0*86400) &
         write(StringTime,'(i2.2,i2.2,i2.2,i2.2)') &
         int(                  Time          /86400.), & ! # days
         int((Time-(86400.*int(Time/86400.)))/ 3600.), & ! # hours
         int((Time-( 3600.*int(Time/ 3600.)))/   60.), & ! # minutes
         int( Time-(   60.*int(Time/   60.)))            ! # seconds
  end subroutine get_time_string
  !============================================================================
  subroutine save_mhd(Time)
    use ModPlotFile, ONLY: save_plot_file
    ! write the output data
    ! separate file is created for each field line,
    ! name format is:
    ! MH_data_<iLon>_<iLat>_n<ddhhmmss>_n<iIter>.{out/dat}
    ! name of the output file
    real, intent(in) :: Time
    character(len=100):: NameFile
    ! header for the file
    character(len=500):: StringHeader
    ! loop variables
    integer:: iLine
    ! index of last particle on the field line
    integer:: iLast
    character(len=*), parameter:: NameSub = 'save_mhd'
    !--------------------------------------------------------------------------
    ! Write ouput files themselves
    StringHeader = &
         'MFLAMPA: data along a field line; '//&
         'Coordindate system: '//trim(TypeCoordBl)//'; '
    do iLine = 1, nLine
       call make_file_name(Time,   &
            iLine         = iLine, &
            NameOut       = NameFile)
         ! get min and max particle indexes on this field line
         iLast  = nVertex_B(   iLine)
         ! print data to file
         call save_plot_file(&
              NameFile      = NameFile,     &
              StringHeaderIn= StringHeader, &
              TypeFileIn    = 'ascii',      &
              nDimIn        = 1,            &
              TimeIn        = Time,         &
              CoordMinIn_D  = [MHData_VIB(LagrID_,1,iLine)], &
              CoordMaxIn_D  = [MHData_VIB(LagrID_,iLast,iLine)], &
              NameVarIn     = NamePlotVar,  &
              VarIn_VI      = MHData_VIB(1:nMHData, 1:iLast, iLine),&
              ParamIn_I     = FootPoint_VB(LagrID_:Z_,iLine))
      end do
  end subroutine save_mhd
  !============================================================================
  logical function BL_is_interface_block(iBlockLocal)
    integer, intent(in) :: iBlockLocal
    !--------------------------------------------------------------------------
    BL_is_interface_block = Used_B(iBlockLocal)
  end function BL_is_interface_block
  !============================================================================
end module CON_bline
!==============================================================================

