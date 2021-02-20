!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module CON_bline
  ! allocate the grid used in this model
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
  logical, public :: UseMflampa_C(MaxComp)  ! To switch coupler for PT
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
  public :: get_bounds
  public :: adjust_line
  public :: BL_n_particle             ! return number of points at the MF line
  public :: BL_put_from_mh            ! put MHD values from MH to SP
  public :: BL_interface_point_coords ! points rMinInterface<R<rMaxInterface
  public :: BL_put_line               ! points rMin < R < rMax
  public :: SP_set_line_foot_b
  !
  ! Grid integer parameters:
  ! Maximum number of vertexes per line

  integer :: nParticleMax

  ! angular grid at origin surface

  integer  :: nLon  = -1
  integer  :: nLat  = -1

  ! Total number of lines on all preceeding PEs

  integer  :: iNode0 = -1 ! = (iProc*nNode)/nProc

  ! Total number of lines on given PE

  integer  :: nBlock = -1 ! = ((iProc +1)*nNode)/nProc - iNode0

  ! Total number of lines on all PEs

  integer  :: nNode = -1

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
  ! MHData_VIB(LagrID_:nMHData, 1:nParticleMax, 1:nBlock)

  real,    allocatable, target :: MHData_VIB(:,:,:)

  ! Aux state vector is a pointer to be joined to a target array
  ! State_VIB(R_:nVar, 1:nParticleMax, 1:nBlock)
  ! nVar may be optionally read in BL_set_grid

  integer  :: nVar = R_
  real,    allocatable, target :: State_VIB(:,:,:)
  !
  ! nParticle_B is a pointer, which is joined to a target array
  ! For stand alone version the target array is allocated here
  ! nParticle_B(1:nBlock)

  integer, allocatable, target :: nParticle_B(:)

  integer, parameter :: Length_=4
  real,    allocatable, target :: FootPoint_VB(:, :)

  !
  integer, allocatable, public :: iOffset_B(:)
  !
  ! Misc:
  integer:: iError
  ! coupling parameters:
  ! domain boundaries
  real, public :: rInterfaceMin, rInterfaceMax
  ! buffer boundaries located near lower (Lo) or upper (Up) boudanry of domain
  real, public :: rBufferLo, rBufferUp
  !
  ! Allowed range of helocentric distances
  real :: rMin, rMax
  public :: rMin
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
            'SP or PT can be MFLAMPA target, but '//NameComp//' cannot'
       iError = 34
       RETURN
    end if
    UseMflampa_C = .false.; UseMflampa_C(BL_) = .true.
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
    real,    intent(inout), pointer :: StateIO_VIB(:,:,:)
    integer, intent(inout), pointer :: nPointer_B(:)
    integer, intent(in)             :: nParticleIn, nVarIn
    integer, intent(in)             :: nLonIn, nLatIn
    real,    intent(inout), pointer :: FootPointIn_VB(:,:)

    !
    ! Loop variable

    integer :: iParticle

    character(len=*), parameter:: NameSub = 'BL_init'
    !--------------------------------------------------------------------------
    nParticleMax = nParticleIn
    nVar         = nVarIn
    nLon = nLonIn; nLat = nLatIn; nNode = nLon*nLat
    iNode0 = (iProc*nNode)/nProc
    nBlock = ((iProc+1)*nNode) / nProc - iNode0
    allocate(MHData_VIB(LagrID_:nMHData, 1:nParticleMax, 1:nBlock), &
         stat=iError)
    call check_allocate(iError, NameSub//'MHData_VIB')
    rPointer_VIB => MHData_VIB
    !
    MHData_VIB = 0.0
    allocate(State_VIB(R_:nVar, 1:nParticleMax, 1:nBlock), &
         stat=iError)
    call check_allocate(iError, NameSub//'State_VIB')
    StateIO_VIB => State_VIB
    !
    State_VIB = -1
    !
    ! reset lagrangian ids
    !
    do iParticle = 1, nParticleMax
       MHData_VIB(LagrID_, iParticle, 1:nBlock) = real(iParticle)
    end do
    !
    allocate(nParticle_B(1:nBlock), &
         stat=iError)
    call check_allocate(iError, NameSub//'nParticle_B')
    nPointer_B => nParticle_B
    !
    nParticle_B = 0
    allocate(FootPoint_VB(LagrID_:Length_, nBlock), stat=iError)
    call check_allocate(iError, NameSub//'FootPoint_VB')
    FootPointIn_VB=>FootPoint_VB
    FootPoint_VB = -1

    allocate(iOffset_B(nBlock)); iOffset_B = 0
  end subroutine BL_init
  !============================================================================
  subroutine BL_get_origin_points(ROrigin, LonMin, LonMax, LatMin, LatMax)
    use ModCoordTransform, ONLY: rlonlat_to_xyz
    !
    ! Parameters of the grid of the bline origin points at the origin surface
    ! at R = ROrigin

    real, intent(in)  :: ROrigin, LonMin, LonMax, LatMin, LatMax

    ! Block (line) number  and corresponing lon-lat indexes
    integer           :: iLat, iLon, iBlock
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
    do iBlock = 1, nBlock
       call iblock_to_lon_lat(iBlock, iLon, iLat)
       nParticle_B(iBlock) = 1
       call rlonlat_to_xyz([ROrigin, LonMin + (iLon - 0.5)*DLon, &
            LatMin + (iLat - 0.5)*DLat], MHData_VIB(X_:Z_,1,iBlock))
    end do
  contains
    !==========================================================================
    subroutine iblock_to_lon_lat(iBlockIn, iLonOut, iLatOut)
      ! return angular grid's indexes corresponding to this block
      integer, intent(in) :: iBlockIn
      integer, intent(out):: iLonOut
      integer, intent(out):: iLatOut

      integer :: iNode
      !------------------------------------------------------------------------
      !
      ! Get node number from block number
      !
      iNode = iBlockIn + iNode0
      iLatOut = 1 + (iNode - 1)/nLon
      iLonOut = iNode - nLon*(iLatOut - 1)
    end subroutine iblock_to_lon_lat
    !==========================================================================
  end subroutine BL_get_origin_points
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
         Block_ = 2    ! Block that has this line/node
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
            CoordMax_D    = [nParticleMax, nLon, nLat] + 0.50,&
            nCells_D      = [nParticleMax, 1, 1],&
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
      integer:: iLat, iLon, iNode, iBlock, iProcNode

      !
      ! fill grid containers
      !

      character(len=*), parameter:: NameSub = 'init_indexes'
      !------------------------------------------------------------------------
      allocate(iGridGlobal_IA(Proc_:Block_, nNode), &
           stat=iError)
      call check_allocate(iError, NameSub//'iGridGlobal_IA')
      iBlock = 1; iProcNode = 0
      do iLat = 1, nLat
         do iLon = 1, nLon
            iNode = iLon + nLon * (iLat-1)
            !
            ! iProcNode = ceiling(real(iNode*nProc)/nNode) - 1
            !
            iGridGlobal_IA(Proc_,   iNode)  = iProcNode
            iGridGlobal_IA(Block_,  iNode)  = iBlock
            if(iNode == ((iProcNode+1)*nNode)/nProc)then
               !
               ! This was the last node on the iProcNode
               !
               iBlock = 1; iProcNode = iProcNode + 1
            else
               iBlock = iBlock + 1
            end if
         end do
      end do
    end subroutine init_indexes
    !==========================================================================
  end subroutine BL_set_grid
  !============================================================================
  subroutine get_bounds(iModelIn, rMinIn, rMaxIn, &
       rBufferLoIn, rBufferUpIn)
    integer,        intent(in) :: iModelIn
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
  end subroutine get_bounds
  !============================================================================
  integer function BL_n_particle(iBlockLocal)
    integer, intent(in) :: iBlockLocal
    !--------------------------------------------------------------------------
    BL_n_particle = nParticle_B(  iBlockLocal)
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
    integer:: i, iBlock
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
       ! cell and block indices
       i      = Put%iCB_II(1, iPutStart + iPartial)
       iBlock = Put%iCB_II(4, iPutStart + iPartial)
       ! interpolation weight
       Weight = W%Weight_I(   iPutStart + iPartial)
       R = sqrt(sum(MHData_VIB(X_:Z_,i,iBlock)**2))
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
            MHData_VIB(Rho_,i,iBlock) = Aux*MHData_VIB(Rho_,i,iBlock) &
            + Buff_I(iRho)/cProtonMass*Weight
       if(DoCoupleVar_V(Pressure_))&
            MHData_VIB(T_,i,iBlock) = Aux*MHData_VIB(T_,i,iBlock) + &
            Buff_I(iP)*(cProtonMass/Buff_I(iRho))*EnergySi2Io*Weight
       if(DoCoupleVar_V(Momentum_))&
            MHData_VIB(Ux_:Uz_,i,iBlock) = Aux*MHData_VIB(Ux_:Uz_,i,iBlock) + &
            Buff_I(iMx:iMz) / Buff_I(iRho) * Weight
       if(DoCoupleVar_V(BField_))&
            MHData_VIB(Bx_:Bz_,i,iBlock) = Aux*MHData_VIB(Bx_:Bz_,i,iBlock) + &
            Buff_I(iBx:iBz) * Weight
       if(DoCoupleVar_V(Wave_))&
            MHData_VIB(Wave1_:Wave2_,i,iBlock) = &
            Aux*MHData_VIB(Wave1_:Wave2_,i,iBlock) + &
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
    integer:: iParticle, iBlock
    real:: R2
    character(len=*), parameter:: NameSub = 'BL_interface_point_coords'
    !--------------------------------------------------------------------------
    iParticle = iIndex_I(1); iBlock    = iIndex_I(4)
    ! Check whether the particle is within interface bounds
    R2 = sum(MHData_VIB(X_:Z_,iParticle,iBlock)**2)
    IsInterfacePoint = &
         R2 >= rInterfaceMin**2 .and. R2 < rInterfaceMax**2
    ! Fix coordinates to be used in mapping
    if(IsInterfacePoint)&
         Xyz_D = MHData_VIB(X_:Z_, iParticle, iBlock)
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
    integer:: iBlock, iParticle
    ! Misc
    real :: R2
    character(len=*), parameter:: NameSub = 'BL_put_line'
    !--------------------------------------------------------------------------
    R2 = sum(Coord_D(1:nDim)**2)
    ! Sort out particles left the SP domain
    if(R2<RMin**2.or.R2>=RMax**2)RETURN
    ! store passed particles
    iBlock    = Put%iCB_II(4,iPutStart)
    iParticle = Put%iCB_II(1,iPutStart) + iOffset_B(iBlock)
    ! put coordinates
    MHData_VIB(X_:Z_,iParticle, iBlock) = Coord_D(1:nDim)
    nParticle_B(iBlock) = MAX(nParticle_B(iBlock), iParticle)
  end subroutine BL_put_line
  !============================================================================
  subroutine adjust_line(iBlock, iBegin, DoAdjustLo, DoAdjustUp)
    integer, intent(in)    :: iBlock
    integer, intent(out)   :: iBegin
    logical, intent(inout) :: DoAdjustLo, DoAdjustUp
    logical:: IsMissing

    real              :: R

    integer, parameter:: Lo_ = 1, Up_ = 2
    integer, parameter:: iIncrement_II(2,2) =reshape([1,0,0,-1],[2,2])
    integer:: iParticle_I(2), iLoop
    ! once new geometry of lines has been put, account for some particles
    ! exiting the domain (can happen both at the beginning and the end)
    integer:: iParticle, iEnd, iOffset ! loop variables
    ! Called after the grid points are received from the
    ! component, nullify offset
    character(len=*), parameter:: NameSub = 'adjust_line'
    !--------------------------------------------------------------------------
    if(DoAdjustLo)iOffset_B(iBlock) = 0
    iBegin = 1
    iEnd   = nParticle_B(  iBlock)
    iParticle_I(:) = [iBegin, iEnd]
    if(DoAdjustUp) then
       iLoop = Up_
    else
       iLoop = Lo_
    end if
    PARTICLE: do while(iParticle_I(1) < iParticle_I(2))
       iParticle = iParticle_I(iLoop)
       R = State_VIB(R_,iParticle,iBlock)
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
       ! and mark losses (decrease nParticle_B) until UP buffer is reached;
       ! after that, if line reenters the model,
       ! particles within the model may not be lost,
       !
       ! MIDDLE models are not allowed to loose particles
       !
       ! HIGHEST model may only loose particles in the tail;
       ! loop Up-2-Lo and mark losses  (decrease nParticle_B)
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
       IsMissing = all(MHData_VIB(X_:Z_,iParticle,iBlock)==0.0)
       if(.not.IsMissing) then
          iParticle_I = iParticle_I + iIncrement_II(:,iLoop)
          CYCLE PARTICLE
       end if

       ! missing point in the lower part of the domain -> ERROR
       if(R < RInterfaceMin)&
            call CON_stop(NameSub//": particle has been lost")

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
          iBegin = iParticle + 1
          iParticle_I = iParticle_I + iIncrement_II(:,iLoop)
          CYCLE PARTICLE
       end if

       ! if need to adjust upper, but not lower boundary -> ADJUST
       if(DoAdjustUp .and. .not.DoAdjustLo)then
          ! push nParticle_B() below current particle;
          ! it will be pushed until it finds a non-missing particle
          nParticle_B(iBlock) = iParticle - 1
          iParticle_I = iParticle_I + iIncrement_II(:,iLoop)
          CYCLE PARTICLE
       end if

       ! remaining case:
       ! need to adjust both boudnaries -> ADJUST,but keep longest range
       if(iParticle - iBegin > nParticle_B(iBlock) - iParticle)then
          nParticle_B(iBlock) = iParticle - 1
          EXIT PARTICLE
       else
          iBegin = iParticle + 1
       end if
       iParticle_I = iParticle_I + iIncrement_II(:,iLoop)
    end do PARTICLE

    DoAdjustLo = RBufferLo == RInterfaceMin
    DoAdjustUp = RBufferUp == RInterfaceMax
  end subroutine adjust_line
  !============================================================================
  subroutine SP_set_line_foot_b(iBlock)
    integer, intent(in) :: iBlock

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
    Xyz1_D = MHData_VIB(X_:Z_, 1, iBlock)

    ! generally, field direction isn't known
    ! approximate it using directions of first 2 segments of the line
    Dist1_D = MHData_VIB(X_:Z_, 1, iBlock) - &
         MHData_VIB(X_:Z_, 2, iBlock)
    Dist1 = sqrt(sum(Dist1_D**2))
    Dist2_D = MHData_VIB(X_:Z_, 2, iBlock) - &
         MHData_VIB(X_:Z_, 3, iBlock)
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
       FootPoint_VB(X_:Z_,iBlock) = Xyz1_D * RMin / sqrt(sum(Xyz1_D**2))
    else
       ! Xyz0, the footprint, is distance Alpha away from Xyz1:
       ! Xyz0 = Xyz1 + Alpha * Dir0 and R0 = RMin =>
       Alpha = S * sqrt(Dot**2 - sum(Xyz1_D**2) + RMin**2) - Dot
       ! Failure (2):
       ! intersection is too far from the current beginning of the line,
       ! use distance between 2nd and 3rd particles on the line as measure
       if(abs(Alpha) > Dist2)then
          ! project first particle for new footpoint
          FootPoint_VB(X_:Z_,iBlock) = Xyz1_D * RMin / sqrt(sum(Xyz1_D**2))
       else
          ! store newly found footpoint of the line
          FootPoint_VB(X_:Z_,iBlock) = Xyz1_D + Alpha * Dir0_D
       end if
    end if

    ! length is used to decide when need to append new particles:
    ! use distance between 2nd and 3rd particles on the line
    FootPoint_VB(Length_,    iBlock) = Dist2
    FootPoint_VB(LagrID_,    iBlock) = MHData_VIB(LagrID_,1,iBlock) - 1.0
  end subroutine SP_set_line_foot_b
  !============================================================================
end module CON_bline

