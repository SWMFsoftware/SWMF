!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module CON_bline
  ! allocate the grid used in this model
#ifdef OPENACC
  use ModUtilities, ONLY: norm2
#endif
  use ModUtilities,      ONLY: check_allocate
  use ModConst,          ONLY: cBoltzmann
  use CON_coupler,       ONLY: MaxComp, is_proc0
  implicit none
  SAVE
  PRIVATE ! Except
  !
  ! public members
  !-----------The following public members are available at all PEs------------
  integer, public :: BL_ =-1                ! ID of the target model (SP_, PT_)
  logical, public :: UseBLine_C(MaxComp)  ! To switch coupler for PT
  !
  ! Boundaries of coupled domains in SC and IH
  !
  real,    public  :: RScMin = 0.0, RScMax = 0.0, RIhMin = 0.0, RIhMax = 0.0
  integer, public  :: Lower_=0, Upper_=-1
  public :: read_param
  public :: BL_set_grid
  !------------The rest of the module is accessed from BL_ model---------------
  integer         :: iProc = -1  ! Rank    of PE in communicator of BL_
  integer         :: nProc = -1  ! Number of PEs in communicator of BL_
  public :: BL_init
  public :: BL_get_origin_points
  public :: BL_get_bounds
  public :: BL_adjust_lines
  public :: BL_n_particle             ! return number of points at the MF line
  public :: BL_put_from_mh            ! put MHD values from MH to SP
  public :: BL_interface_point_coords ! points rMinInterface<R<rMaxInterface
  public :: BL_put_line               ! points rMin < R < rMax
  public :: BL_set_line_foot_b
  !
  ! Grid integer parameters:
  ! Maximum number of vertexes per line

  integer :: nVertexMax

  ! angular grid at origin surface

  integer  :: nLon  = -1
  integer  :: nLat  = -1

  ! Total number of lines on all preceeding PEs

  integer  :: iLineAll0 = -1 ! = (iProc*nLineAll)/nProc

  ! Total number of lines on given PE

  integer, public  :: nLine = -1 ! = ((iProc +1)*nLineAll)/nProc - iLineAll0

  ! Total number of lines on all PEs

  integer  :: nLineAll = -1

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
  !
  ! MHD state vector is a pointer to be joined to a target array
  ! MHData_VIB(LagrID_:nMHData, 1:nVertexMax, 1:nLine)

  real,    allocatable, target :: MHData_VIB(:,:,:)

  ! Aux state vector is a pointer to be joined to a target array
  ! State_VIB(R_:nVar, 1:nVertexMax, 1:nLine)
  ! nVar may be optionally read in BL_set_grid

  integer  :: nVar = R_
  real,    allocatable, target :: State_VIB(:,:,:)
  !
  ! nVertex_B is a pointer, which is joined to a target array
  ! For stand alone version the target array is allocated here
  ! nVertex_B(1:nLine)

  integer, allocatable, target :: nVertex_B(:)

  integer, parameter :: Length_=4
  real,    allocatable, target :: FootPoint_VB(:, :)

  !
  integer, allocatable, public :: iOffset_B(:)
  !
  ! Misc:
  integer:: iError
  ! coupling parameters:
  ! domain boundaries
  real :: rInterfaceMin, rInterfaceMax
  ! buffer boundaries located near lower (Lo) or upper (Up) boudanry of domain
  real :: rBufferLo, rBufferUp
  !
  ! Allowed range of helocentric distances
  real :: rMin, rMax
  ! Coefficient to transform energy
  real         :: EnergySi2Io = 1.0/cBoltzmann
  character(len=*), parameter:: NameMod = 'CON_mflampa::'
contains
  !============================================================================
  subroutine read_param(iError)
    use CON_coupler, ONLY: use_comp, i_comp, SP_, PT_, SC_, IH_, &
         i_proc, n_proc, is_proc
    use ModReadParam, ONLY: read_var
    integer, intent(inout) :: iError
    integer :: nSource, iSource, iComp
    character(len=2) :: NameComp = ''
    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    call read_var('NameTarget', NameComp)
    BL_ = i_comp(NameComp)
    if(.not.use_comp(BL_)) then
       if(is_proc0()) write(*,*) NameSub//&
            ' SWMF_ERROR for NameMaster: '// &
            NameComp//' is OFF or not registered, not MFLAMPA target'
       iError = 34
       RETURN
    end if
    if(.not.(BL_ == PT_.or.BL_==SP_))then
       if(is_proc0()) write(*,*) NameSub//&
            ' SWMF_ERROR for NameMaster: '// &
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
    if(nSource/=2)then
       if(is_proc0()) write(*,*) NameSub//&
            ' SWMF_ERROR for NameMaster: '// &
            'there can be only two source components in MFLAMPA'
       iError = 34
       RETURN
    end if
    do iSource = 1, nSource
       call read_var('NameSource', NameComp)
       iComp = i_comp(NameComp)
       select case(iComp)
       case(SC_)
          if(.not.use_comp(SC_)) then
             if(is_proc0()) write(*,*) NameSub//&
                  ' SWMF_ERROR for NameMaster: '// &
                  ' SC  needed in MFLAMPA is OFF or not registered'
             iError = 34
             RETURN
          end if
          call read_var('RScMin', RScMin)
          call read_var('RScMax', RScMax)
       case(IH_)
          if(.not.use_comp(IH_)) then
             if(is_proc0()) write(*,*) NameSub//&
                  ' SWMF_ERROR for NameMaster: '// &
                  ' IH  needed in MFLAMPA is OFF or not registered'
             iError = 34
             RETURN
          end if
          call read_var('RIhMin', RIhMin)
          call read_var('RIhMax', RIhMax)
       case default
          if(is_proc0()) write(*,*) NameSub//&
               ' SWMF_ERROR for NameMaster: '// &
               'SC or IH are only allowed as MFLAMPA sources'
          iError = 34
          RETURN
       end select
    end do
    if(.not.is_proc(BL_))RETURN
    RMin = RScMin; RMax = RIhMax
    Lower_ = SC_; Upper_ = IH_
  end subroutine read_param
  !============================================================================
  subroutine BL_init(nParticleIn, nLonIn, nLatIn, &
       rPointer_VIB,  nPointer_B,  nVarIn, StateIO_VIB, FootPointIn_VB)

    real,    intent(inout), pointer :: rPointer_VIB(:,:,:)
    integer, intent(inout), pointer :: nPointer_B(:)
    integer, intent(in)             :: nParticleIn, nLonIn, nLatIn
    integer, optional, intent(in)   :: nVarIn
    real,    intent(inout), optional, pointer :: &
         StateIO_VIB(:,:,:), FootPointIn_VB(:,:)

    !
    ! Loop variable

    integer :: iVertex

    character(len=*), parameter:: NameSub = 'BL_init'
    !--------------------------------------------------------------------------
    nVertexMax            = nParticleIn
    nLon = nLonIn; nLat = nLatIn; nLineAll = nLon*nLat
    iLineAll0 = (iProc*nLineAll)/nProc
    nLine = ((iProc+1)*nLineAll) / nProc - iLineAll0

    ! Set MHD array
    
    allocate(MHData_VIB(LagrID_:nMHData, 1:nVertexMax, 1:nLine), &
         stat=iError)
    call check_allocate(iError, NameSub//'MHData_VIB')
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
    call check_allocate(iError, NameSub//'nVertex_B')
    nPointer_B => nVertex_B
    nVertex_B = 0

    ! Set auxiliary State_VIB array
    
    if(present(nVarIn))nVar = nVarIn
    allocate(State_VIB(R_:nVar, 1:nVertexMax, 1:nLine), &
         stat=iError)
    call check_allocate(iError, NameSub//'State_VIB')
    if(present(StateIO_VIB))StateIO_VIB => State_VIB
    State_VIB = -1

    ! Set FootPoint_VB array
    
    allocate(FootPoint_VB(LagrID_:Length_, nLine), stat=iError)
    call check_allocate(iError, NameSub//'FootPoint_VB')
    if(present(FootPointIn_VB))FootPointIn_VB => FootPoint_VB
    FootPoint_VB = -1

    ! Set iOffset_B array

    allocate(iOffset_B(nLine)); iOffset_B = 0
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
    ! angular grid's step
    !
    DLon = (LonMax - LonMin)/nLon
    !
    ! angular grid's step
    !
    DLat = (LatMax - LatMin)/nLat
    do iLine = 1, nLine
       call BL_iblock_to_lon_lat(iLine, iLon, iLat)
       nVertex_B(iLine) = 1
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
    !
    ! Get node number from line number
    !
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
    call init_decomposition(&
         GridID_ = BL_,&
         CompID_ = BL_,&
         nDim    = nDim)
    ! Construct decomposition

    if(is_proc0(BL_))then
       call init_indexes
       call get_root_decomposition(&
            GridID_       = BL_,&
            iRootMapDim_D = [1, nLon, nLat],&
            CoordMin_D    = [0.50, 0.50, 0.50],&
            CoordMax_D    = [nVertexMax, nLon, nLat] + 0.50,&
            nCells_D      = [nVertexMax, 1, 1],&
            PE_I          = iGridGlobal_IA(Proc_,:),&
            iBlock_I      = iGridGlobal_IA(Block_,:))
       deallocate(iGridGlobal_IA)
    end if
    call bcast_decomposition(BL_)
    ! Coordinate system is Heliographic Inertial Coordinate System (HGI)
    ! with length measured in solar radii
    call set_coord_system(&
         GridID_      = BL_, &
         TypeCoord    = TypeCoordSystem, &
         TypeGeometry = 'cartesian', &
         NameVar      = NameVarCouple, &
         UnitX        = UnitX)
    if(.not.is_proc(BL_))RETURN
    if(present(EnergySi2IoIn))EnergySi2Io = EnergySi2IoIn
  contains
    !==========================================================================
    subroutine init_indexes
      integer:: iLat, iLon, iLineAll, iLine, iProcNode

      !
      ! fill grid containers
      !

      character(len=*), parameter:: NameSub = 'init_indexes'
      !------------------------------------------------------------------------
      allocate(iGridGlobal_IA(Proc_:Block_, nLineAll), &
           stat=iError)
      call check_allocate(iError, NameSub//'iGridGlobal_IA')
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
  subroutine BL_put_from_mh(nPartial,iPutStart,Put,W,DoAdd,Buff_I,nVar)
    use CON_coupler, ONLY: &
         iVar_V, DoCoupleVar_V, &
         Density_, RhoCouple_, Pressure_, PCouple_, &
         Momentum_, RhoUxCouple_, RhoUzCouple_, &
         BField_, BxCouple_, BzCouple_, &
         Wave_, WaveFirstCouple_, WaveLastCouple_
    use CON_router, ONLY: IndexPtrType, WeightPtrType
    use ModConst,   ONLY: cProtonMass
    integer,intent(in)::nPartial,iPutStart,nVar
    type(IndexPtrType),intent(in)::Put
    type(WeightPtrType),intent(in)::W
    logical,intent(in)::DoAdd
    real,dimension(nVar),intent(in)::Buff_I
    integer:: iRho, iP, iMx, iMz, iBx, iBz, iWave1, iWave2
    integer:: i, iLine
    integer:: iPartial
    real:: Weight
    real:: R, Aux

    character(len=100):: StringError

    ! check consistency of DoCoupleVar_V

    character(len=*), parameter:: NameSub = 'BL_put_from_mh'
    !--------------------------------------------------------------------------
    if(.not. DoCoupleVar_V(Density_) .and. &
         (DoCoupleVar_V(Pressure_) .or. DoCoupleVar_V(Momentum_)))&
         call CON_Stop(NameSub//': pressure or momentum is coupled,'&
         //' but density is not')
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
       ! cell and line indices
       i      = Put%iCB_II(1, iPutStart + iPartial)
       iLine = Put%iCB_II(4, iPutStart + iPartial)
       ! interpolation weight
       Weight = W%Weight_I(   iPutStart + iPartial)
       R = norm2(MHData_VIB(X_:Z_,i,iLine))
       if(R >= rInterfaceMin .and. R < rBufferLo)then
          Aux = 1.0
          Weight = Weight * (0.50 + 0.50*tanh(2*(2*R - &
               RBufferLo - RInterfaceMin)/(RBufferLo - RInterfaceMin)))
       end if
       if(R >= rBufferUp .and. R < rInterfaceMax)then
          Weight = Weight * (0.50 - 0.50*tanh(2*(2*R - &
               RInterfaceMax - RBufferUp)/(RInterfaceMax - RBufferUp)))
       end if
       ! put the data
       if(DoCoupleVar_V(Density_))&
            MHData_VIB(Rho_,i,iLine) = Aux*MHData_VIB(Rho_,i,iLine) &
            + Buff_I(iRho)/cProtonMass*Weight
       if(DoCoupleVar_V(Pressure_))&
            MHData_VIB(T_,i,iLine) = Aux*MHData_VIB(T_,i,iLine) + &
            Buff_I(iP)*(cProtonMass/Buff_I(iRho))*EnergySi2Io*Weight
       if(DoCoupleVar_V(Momentum_))&
            MHData_VIB(Ux_:Uz_,i,iLine) = Aux*MHData_VIB(Ux_:Uz_,i,iLine) + &
            Buff_I(iMx:iMz) / Buff_I(iRho) * Weight
       if(DoCoupleVar_V(BField_))&
            MHData_VIB(Bx_:Bz_,i,iLine) = Aux*MHData_VIB(Bx_:Bz_,i,iLine) + &
            Buff_I(iBx:iBz) * Weight
       if(DoCoupleVar_V(Wave_))&
            MHData_VIB(Wave1_:Wave2_,i,iLine) = &
            Aux*MHData_VIB(Wave1_:Wave2_,i,iLine) + &
            Buff_I(iWave1:iWave2)*Weight
    end do
  end subroutine BL_put_from_mh
  !============================================================================
  subroutine BL_interface_point_coords(nDim, Xyz_D, &
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
    integer:: iVertex, iLine
    real:: R2
    character(len=*), parameter:: NameSub = 'BL_interface_point_coords'
    !--------------------------------------------------------------------------
    iVertex = iIndex_I(1); iLine    = iIndex_I(4)
    ! Check whether the particle is within interface bounds
    R2 = sum(MHData_VIB(X_:Z_,iVertex,iLine)**2)
    IsInterfacePoint = &
         R2 >= rInterfaceMin**2 .and. R2 < rInterfaceMax**2
    ! Fix coordinates to be used in mapping
    if(IsInterfacePoint)&
         Xyz_D = MHData_VIB(X_:Z_, iVertex, iLine)
  end subroutine BL_interface_point_coords
  !============================================================================
  subroutine BL_put_line(nPartial, iPutStart, Put,&
       Weight, DoAdd, Coord_D, nVar)
    use CON_router, ONLY: IndexPtrType, WeightPtrType
    integer, intent(in) :: nPartial, iPutStart, nVar
    type(IndexPtrType), intent(in) :: Put
    type(WeightPtrType),intent(in) :: Weight
    logical,            intent(in) :: DoAdd
    real,               intent(in) :: Coord_D(nVar) ! nVar=nDim

    ! indices of the particle
    integer:: iLine, iVertex
    ! Misc
    real :: R2
    character(len=*), parameter:: NameSub = 'BL_put_line'
    !--------------------------------------------------------------------------
    iLine    = Put%iCB_II(4,iPutStart)
    iVertex = Put%iCB_II(1,iPutStart) + iOffset_B(iLine)
    R2 = sum(Coord_D(1:nDim)**2)
    if(R2<RMin**2.or.R2>=RMax**2)then
       ! Sort out particles left the SP domain
       MHData_VIB(X_:Z_,iVertex, iLine) = 0.0
    else
       ! store passed particles
       ! put coordinates
       MHData_VIB(X_:Z_,iVertex, iLine) = Coord_D(1:nDim)
       nVertex_B(iLine) = MAX(nVertex_B(iLine), iVertex)
    end if
  end subroutine BL_put_line
  !============================================================================
  ! Called from coupler after the updated grid point lo<cation are
  ! received from the other component (SC, IH). Determines whether some
  ! grid points should be added/deleted
  subroutine BL_adjust_lines(DoInit, Source_)
    logical, intent(in) :: DoInit
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
       iBegin = 1
       if(DoInit.and.DoAdjustLo)then
          do while (all(MHData_VIB(X_:Z_,iBegin,iLine)==0.0))
             iBegin = iBegin + 1
          end do
          if(iBegin == 1)call BL_set_line_foot_b(iLine)
       end if
       if(DoAdjustLo)then
          iOffset_B(iLine) = 0
       end if
       iEnd   = nVertex_B(  iLine)
       iParticle_I(:) = [iBegin, iEnd]
       if(DoAdjustUp) then
          iLoop = Up_
       else
          iLoop = Lo_
       end if
       PARTICLE: do while(iParticle_I(1) < iParticle_I(2))
          iVertex = iParticle_I(iLoop)
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
          if(iLoop == Up_ .and. (&
               R <  RInterfaceMin&
               .or.&
               R >= RInterfaceMax)&
               )then
             DoAdjustLo = .false.
             DoAdjustUp = .false.
          end if

          ! determine whether particle is missing
          IsMissing = all(MHData_VIB(X_:Z_,iVertex,iLine)==0.0)
          if(.not.IsMissing) then
             iParticle_I = iParticle_I + iIncrement_II(:,iLoop)
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
                write(*,*)MHData_VIB(X_:Z_,iVertex,iLine)
             end do
             call CON_stop(NameSub//": particle has been lost")
          end if
          ! missing point in the upper part of the domain -> IGNORE;
          ! if needed to adjust beginning, then it is done,
          ! switch left -> right end of range and start adjusting
          ! tail of the line, if it has reentered current part of the domain
          if(R >= RInterfaceMax)then
             if(iLoop == Lo_)&
                  iLoop = Up_
             iParticle_I = iParticle_I + iIncrement_II(:,iLoop)
             if(DoAdjustLo)then
                DoAdjustLo = .false.
                DoAdjustUp = .true.
             end if
             CYCLE PARTICLE
          end if

          ! if point used to be in a upper buffer -> IGNORE
          if(R >= rBufferUp .and. R < rInterfaceMax)then
             iParticle_I = iParticle_I + iIncrement_II(:,iLoop)
             CYCLE PARTICLE
          end if

          ! if need to adjust lower, but not upper boundary -> ADJUST
          if(DoAdjustLo .and. .not.DoAdjustUp)then
             ! push iBegin in front of current particle;
             ! it will be pushed until it finds a non-missing particle
             iBegin = iVertex + 1
             iParticle_I = iParticle_I + iIncrement_II(:,iLoop)
             CYCLE PARTICLE
          end if

          ! if need to adjust upper, but not lower boundary -> ADJUST
          if(DoAdjustUp .and. .not.DoAdjustLo)then
             ! push nVertex_B() below current particle;
             ! it will be pushed until it finds a non-missing particle
             nVertex_B(iLine) = iVertex - 1
             iParticle_I = iParticle_I + iIncrement_II(:,iLoop)
             CYCLE PARTICLE
          end if

          ! remaining case:
          ! need to adjust both boudnaries -> ADJUST,but keep longest range
          if(iVertex - iBegin > nVertex_B(iLine) - iVertex)then
             nVertex_B(iLine) = iVertex - 1
             EXIT PARTICLE
          else
             iBegin = iVertex + 1
          end if
          iParticle_I = iParticle_I + iIncrement_II(:,iLoop)
       end do PARTICLE

       DoAdjustLo = Source_ == Lower_
       DoAdjustUp = Source_ == Upper_
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
    subroutine append_particles
      ! appends a new particle at the beginning of lines if necessary
      integer:: iLine
      real:: DistanceToMin
      real, parameter:: cTol = 1E-06

      character(len=*), parameter:: NameSub = 'append_particles'
      !------------------------------------------------------------------------
      line:do iLine = 1, nLine
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
         MHData_VIB(LagrID_:Z_,  1, iLine) = &
              FootPoint_VB(LagrID_:Z_, iLine)*(1.0 + cTol)
         State_VIB(R_,          1, iLine) = &
              sqrt(sum((MHData_VIB(X_:Z_,  1, iLine))**2))
         MHData_VIB(LagrID_,1, iLine) = MHData_VIB(LagrID_, 2, iLine) - 1.0
         FootPoint_VB(LagrID_,iLine) = MHData_VIB(LagrID_, 1, iLine) - 1.0
      end do line
    end subroutine append_particles
    !==========================================================================
  end subroutine BL_adjust_lines
  subroutine adjust_line(iLine, iBegin, DoAdjustLo, DoAdjustUp)
    integer, intent(in)    :: iLine
    integer, intent(out)   :: iBegin
    logical, intent(inout) :: DoAdjustLo, DoAdjustUp

    character(len=*), parameter:: NameSub = 'adjust_line'
    !--------------------------------------------------------------------------
    
  end subroutine adjust_line
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
    ! no intersection of smoothly extended line with the sphere R = RMin
    if(Dot**2 - sum(Xyz1_D**2) + RMin**2 < 0)then
       ! project first particle for new footpoint
       FootPoint_VB(X_:Z_,iLine) = Xyz1_D * RMin / sqrt(sum(Xyz1_D**2))
    else
       ! Xyz0, the footprint, is distance Alpha away from Xyz1:
       ! Xyz0 = Xyz1 + Alpha * Dir0 and R0 = RMin =>
       Alpha = S * sqrt(Dot**2 - sum(Xyz1_D**2) + RMin**2) - Dot
       ! Failure (2):
       ! intersection is too far from the current beginning of the line,
       ! use distance between 2nd and 3rd particles on the line as measure
       if(abs(Alpha) > Dist2)then
          ! project first particle for new footpoint
          FootPoint_VB(X_:Z_,iLine) = Xyz1_D * RMin / sqrt(sum(Xyz1_D**2))
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
end module CON_bline

