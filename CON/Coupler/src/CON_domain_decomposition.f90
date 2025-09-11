!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module CON_domain_decomposition
  ! This file presents the class of the domain decompositions which
  ! icludes both the uniformly spaced ones (uniformly spaced with
  ! respect to some generalized coordinates) and Octree or Quadric
  ! tree for adaptive block decompositions.

  use ModUtilities, ONLY: check_allocate, CON_stop
  use CON_world
  use ModNumConst
  use ModMpi
  use ModKind, ONLY: nByteReal

  ! revision history:
  ! 6/18/03-7/11/03 Sokolov I.V. <igorsok@umich.edu> phone(734)647-4705

  implicit none

  SAVE

  real, parameter :: cThird = 1.0/3.0
  integer, parameter:: &
       Parent_           = -1, &
       MyNumberAsAChild_ =  0, &
       FirstChild_       =  1, &
       BLK_              =  2, &
       PE_               =  3, &
       None_             = -777

  type DomainType
     ! The type is to describe the domain decomposition for a
     ! component. The component should be properly registered with
     ! registry procedures and should have a unique ID
     integer :: CompID_
     ! This type describes the general decomposition for logically
     ! rectangular box domain. Describe the property of the  domain first.
     ! nDim := the dimensionality of the generalized coordinate space
     integer ::nDim
     ! CoordMin_D:=Coordinates of the left corner of the domain
     ! CoordMax_D:=Coordinates of the right corner of the domai
     real, pointer :: CoordMin_D(:), CoordMax_D(:)
     ! Root decomposition:=  decomposition of the domain for
     ! iRootMapDim_D(1)*iRootMamDim_D(2)*.....iRootMapDim_D(nDim)
     ! blocks, which are spatially uniform in geralized coordinates
     !
     ! The principle "block is not shared" is assumed, that is if the
     ! physical data at any point, geometrically belonging to the block,
     ! are allocated at the given PE, then all the points, which geometrically
     ! belong to this block, should be also allocated at the same PE.
     !
     ! (For the uniformly spaced grids the root decomposition is the
     ! map of blocks, for tree grid - the map of the roots)
     integer, pointer :: iRootMapDim_D(:)
     ! Each of  block is split for
     ! nCell_D(1)*nCell_D(2)*.....nCell_D(nDim) cells,
     ! which are spatially uniform in the geralized coordinates
     integer, pointer :: nCell_D(:)
     ! tree decomposition:= root decomposition .or.&
     !     refinement[tree decomposition]
     ! refinement:=decompose any block for 2**nDim
     ! blocks by splitting the original block for two equal parts
     ! along each dimension
     !
     ! The mathematical space of tree decompositions, for a given
     ! iRootMapDim_D can be one-to-one mapped to the space of
     ! structures, consisting of product(iRootMapDim_D) tree graphs.
     ! Each node of this tree corresponds to either a root map
     ! element, or the block, which was refined, or a block.
     !
     ! If the tree decomposition is allowed for a given domain, set
     ! IsTreeGrid=.true.
     logical :: IsTreeDD
     ! nChildren=2**nDim for the tree grid
     integer :: nChildren
     ! For tree decomposition this is a full tree structure saved as
     ! described in the figure (PE is a CPU rank in a i_comm(Global_)
     integer, pointer::iDD_II(:,:)
     !---------------------------------------------------------------
     !
     !------Root  node of tree ----------             ...
     !   RootNumber,MyNumberAsAChild=0--ChildAdresses()
     !                                                \ (1:nChildren)
     !                                                 ...
     !
     !------Non-end node of tree---------
     !                   ...             ...
     !                  /                 /
     !   ParentAddress--MyNumberAsAChild--ChildAdresses()
     !                  \                 \(1:nChildren)
     !                   ...              ...
     !
     !------End node of tree-------------
     !                 ...
     !                  /            information for used blocks
     !   ParentAddress--MyNumberAsAChild,    (PE,BLK),
     !                  \            "address of the first
     !                  ...           child"=None_
     ! nTreeNode - the number of blocks, for non-tree grid, or the
     ! number of all nodes in the tree.
     integer::nTreeNode, nAllocatedNodes
     !---------------------------------------------------------------
     ! The following is not allocated if the grid is a tree grid.
     ! The iTreeNode of the tree roots in the
     ! iDD_II array (the values for the second index)
     integer, pointer :: iRoot_I(:)
     ! In index I  passing from 1 to nChildren, all the possible
     ! combinations of 0 and 1 of dimension nDim should be
     ! unambigously related to the values of I In the grid captured by
     ! the control module, this array should be organized exaclty as
     ! it done in the "component".
     integer, pointer :: iShift_DI(:, :)
     ! First index inumerates the spatial component of the
     ! generalized coordinates, second is lGlobalNode
     ! Generalized coordinates of the left corner of the block,
     real, pointer :: CoordBlock_DI(:, :)
     ! DCoordBlock_D is the BLOCK size
     real, pointer :: DCoordBlock_DI(:,:)
     ! and DCoordCell_D is the cell size.
     real, pointer :: DCoordCell_DI(:,:)
     ! For the domain with periodic boundary conditions:
     logical, pointer::IsPeriodic_D(:)
     ! For more complicated cases (spherical and analogous geometries
     ! having a pole singularity) a special procedure
     ! glue_margin is invoked when required. To do this
     ! for a particular grid, it should have DoGlueMargins=.true.
     logical :: DoGlueMargins
     integer :: iDirMinusGlue, iDirPlusGlue, iDirCycle
     ! To search in the tree decompositions  it is optimal to start
     ! from the previously found node (from the previous search). To
     ! do this we have this previuosly found value saved in lSearch
     integer :: lSearch
     ! For the grids which can be refined, the realization number for
     ! the domain decomposition should be stored
     integer :: iRealization
     ! In the local decomposition the PE ranks are defined to be in a
     ! local communicator, for a particular component
     logical :: IsLocal
  end type DomainType

  type DomainPointerType
     type(DomainType),pointer::Ptr
  end type DomainPointerType

  ! Needed for searching and interpolating algorithms
  real, parameter :: cTol  = 0.00000010, cTol2 = cTol**(nByteReal/4)
  real, parameter :: cAlmostOne   = 1.0 - cTol2
contains
  !============================================================================
  subroutine init_decomposition_dd(Domain, CompID_, nDim, IsTreeDD, IsLocal)
    type(DomainType),intent(inout) :: Domain
    integer,           intent(in)  :: CompID_, nDim
    logical, optional, intent(in)  :: IsTreeDD
    logical, optional, intent(in)  :: IsLocal

    integer::iError,iDim,iChildFirst,iChildLast
    !--------------------------------------------------------------------------
    if(.not.use_comp(CompID_))call CON_stop('The component is not used')
    Domain%CompID_   = CompID_
    Domain%nDim      = nDim
    ! Allocate arrays with the the size of nDim
    nullify(Domain%iRootMapDim_D)
    allocate(Domain%iRootMapDim_D(nDim), stat=iError)
    call check_allocate(iError, "iRootMapDim_D")
    Domain%iRootMapDim_D = 1

    nullify(Domain%CoordMin_D)
    allocate(Domain%CoordMin_D(nDim),stat=iError)
    call check_allocate(iError,"CoordMin_D")
    Domain%CoordMin_D = 0.0

    nullify(Domain%CoordMax_D)
    allocate(Domain%CoordMax_D(nDim),stat=iError)
    call check_allocate(iError,"CoordMax_D")
    Domain%CoordMax_D = 1.0

    nullify(Domain%nCell_D)
    allocate(Domain%nCell_D(nDim),stat=iError)
    call check_allocate(iError,"nCell_D")
    Domain%nCell_D = 1

    nullify(Domain%IsPeriodic_D)
    allocate(Domain%IsPeriodic_D(nDim), stat=iError)
    call check_allocate(iError,"IsPeriodic_D")
    Domain%IsPeriodic_D  = .false.

    Domain%DoGlueMargins = .false.
    Domain%iDirMinusGlue = 0
    Domain%iDirPlusGlue  = 0
    Domain%iDirCycle     = 0

    Domain%IsTreeDD  = .false.
    Domain%nChildren = 1
    if(present(IsTreeDD))Domain%IsTreeDD = IsTreeDD
    if(Domain%IsTreeDD)then
       Domain%nChildren = 2**Domain%nDim
       nullify(Domain%iRoot_I)
       call allocate_iroot(Domain) ! Single element in array

       nullify(Domain%iShift_DI)
       allocate(Domain%iShift_DI(Domain%nDim, Domain%nChildren), stat=iError)
       call check_allocate(iError,"iShift_DI")
       Domain%iShift_DI = 0
       ! A "binary" order for iShift_DI is set as a default:
       !(0,0,0,&
       ! 1,0,0,&
       ! 0,1,0,&
       ! and so on
       !
       iChildLast=1
       do iDim=1, Domain%nDim
          iChildFirst = iChildLast + 1
          iChildLast  = min(Domain%nChildren, 2*iChildLast)
          if(iChildFirst < iChildLast)EXIT
          Domain%iShift_DI(:, iChildFirst:iChildLast) = &
               Domain%iShift_DI(:, 1:1 + iChildLast - iChildFirst)
          Domain%iShift_DI(iDim, iChildFirst:iChildLast) = 1
       end do
    end if
    Domain%IsLocal = .false.; if(present(IsLocal)) Domain%IsLocal = IsLocal

    Domain%nTreeNode      =  1; Domain%nAllocatedNodes = -1

    nullify(Domain%iDD_II)
    nullify(Domain%CoordBlock_DI)
    nullify(Domain%DCoordBlock_DI)
    nullify(Domain%DCoordCell_DI)
    call check_octree_allocation(Domain)
    Domain%lSearch      = 1
    Domain%iRealization = 0
  end subroutine init_decomposition_dd
  !============================================================================
  subroutine check_octree_allocation(Domain)
    ! Checks the array allocation for a given grid desciptor and
    ! extends its dimension if required

    type(DomainType), intent(inout) :: Domain
    integer::iError,nUbound
    !--------------------------------------------------------------------------
    nUbound = PE_
    if(Domain%IsTreeDD)nUbound = max(nUbound, Domain%nChildren)
    if(Domain%nAllocatedNodes >= Domain%nTreeNode)RETURN
    if(associated(Domain%iDD_II))deallocate( Domain%iDD_II, &
       Domain%CoordBlock_DI, Domain%DCoordBlock_DI, Domain%DCoordCell_DI)

    nullify(Domain%iDD_II)
    allocate(Domain%iDD_II(Parent_:nUBound, Domain%nTreeNode), stat=iError)
    call check_allocate(iError, 'iDD_II')
    Domain%iDD_II  = None_; Domain%iDD_II(MyNumberAsAChild_,:) = 0

    nullify(Domain%CoordBlock_DI)
    allocate(Domain%CoordBlock_DI(Domain%nDim, Domain%nTreeNode), stat=iError)
    call  check_allocate(iError,'CoordBlock_DI')
    Domain% CoordBlock_DI = 0.0

    nullify(Domain%DCoordBlock_DI)
    allocate(Domain%DCoordBlock_DI(Domain%nDim,Domain%nTreeNode),stat=iError)
    call  check_allocate(iError,'DCoordBlock_DI')
    Domain%DCoordBlock_DI = 1.0

    nullify(Domain%DCoordCell_DI)
    allocate(Domain%DCoordCell_DI(Domain%nDim, Domain%nTreeNode), stat=iError)
    call  check_allocate(iError,'DCoordCell_DI')
    Domain%DCoordCell_DI = 1.0

    Domain%nAllocatedNodes = Domain%nTreeNode
  end subroutine check_octree_allocation
  !============================================================================
  subroutine get_root_decomposition_dd(&
       Domain,&       ! Decomposition to be constructed
       iRootMapDim_D,&! As in DomainType
       CoordMin_D,&   ! As in DomainType
       CoordMax_D,&   ! As in DomainType
       nCell_D,&      ! As in DomainType
       PE_I,&         ! PE layout
       iBlock_I,&     ! Local Block Number layout
       IsPeriodic_D,& ! As in DomainType
       iShift_DI,   & ! As in DomainType
       DoGlueMargins,&! As in DomainType
       iDirMinusGlue,&! As in DomainType
       iDirPlusGlue,& ! As in DomainType
       iDirCycle)     ! As in DomainType
    ! To get a decomposition domain, even the tree one, the root
    ! decomposition should be first constructed
    !                        ATTENTION
    ! PE here are the ranks in the LOCAL communicator for the
    ! component

    type(DomainType), intent(inout) :: Domain
    integer, intent(in) :: iRootMapDim_D(Domain%nDim)
    real,    intent(in) :: CoordMin_D(Domain%nDim)
    real,    intent(in) :: CoordMax_D(Domain%nDim)
    integer, intent(in) :: nCell_D(Domain%nDim)
    integer, optional, intent(in) :: PE_I(:), iBlock_I(:)
    logical, optional, intent(in) :: IsPeriodic_D(Domain%nDim)
    integer, optional, intent(in) :: iShift_DI(Domain%nDim, Domain%nChildren)
    logical, optional, intent(in) :: DoGlueMargins
    integer, optional, intent(in) :: iDirMinusGlue
    integer, optional, intent(in) :: iDirPlusGlue
    integer, optional, intent(in) :: iDirCycle
    integer :: lBlock, MaxBlock
    ! Assign the grid parameters                                     !
    !--------------------------------------------------------------------------
    Domain%iRootMapDim_D = iRootMapDim_D
    Domain%CoordMin_D    = CoordMin_D
    Domain%CoordMax_D    = CoordMax_D
    Domain%nCell_D       = nCell_D

    if(present(DoGlueMargins)) Domain%DoGlueMargins = DoGlueMargins
    if(present(iDirMinusGlue)) Domain%iDirMinusGlue = iDirMinusGlue
    if(present(iDirPlusGlue))  Domain%iDirPlusGlue  = iDirPlusGlue
    if(present(iDirCycle))     Domain%iDirCycle     = iDirCycle
    if(present(IsPeriodic_D))  Domain%IsPeriodic_D  = IsPeriodic_D

    Domain%nTreeNode = product(iRootMapDim_D)
    call check_octree_allocation(Domain)

    Domain%iDD_II(MyNumberAsAChild_, 1:Domain%nTreeNode) = 0
    Domain%iDD_II(FirstChild_,       1:Domain%nTreeNode) = None_
    do lBlock = 1, Domain%nTreeNode
       Domain%iDD_II(Parent_, lBlock) = lBlock
    end do

    if(Domain%IsTreeDD)then
       if(ubound(Domain%iRoot_I, 1) < product(Domain%iRootMapDim_D))&
            call allocate_iroot(Domain)
       if(present(iShift_DI))Domain%iShift_DI = iShift_DI
    end if

    ! If neither PE_I nor iBlock_I is present, the root
    ! decomposition blocks are balanced for an optimal load
    if(.not.present(PE_I).and. .not.present(iBlock_I))then
       MaxBlock=(Domain%nTreeNode - 1)/n_proc(Domain%CompID_) + 1
       do lBlock=1, Domain%nTreeNode
          Domain%iDD_II(BLK_, lBlock) = mod(lBlock - 1, MaxBlock) + 1
          Domain%iDD_II(PE_, lBlock)  = (lBlock - 1)/MaxBlock
       end do
       RETURN
    end if
    if(present(PE_I))then
       ! Check size of PE_I                       !
       if(ubound(PE_I,1)/=product(iRootMapDim_D))&
            call CON_stop('Size of PE_I is wrong')
       Domain%iDD_II(PE_,1:Domain%nTreeNode) = PE_I
    else
       ! If the PE\_I is not given, it is assumed that all the blocks
       ! are at the root PE
       Domain%iDD_II(PE_,1:Domain%nTreeNode) = 0
    end if

    if(present(iBlock_I))then
       ! Check size of iBlock_I
       if(ubound(iBlock_I,1)/=product(iRootMapDim_D))&
            call CON_stop('Size of iBlock_I is wrong')
       Domain%iDD_II(BLK_, 1:Domain%nTreeNode) = iBlock_I
    else
       do lBlock=1, Domain%nTreeNode
          Domain%iDD_II(BLK_, lBlock) = lBlock
       end do
    end if
  end subroutine get_root_decomposition_dd
  !============================================================================
  subroutine allocate_iroot(Domain)
    ! Allocates the part of the grid descriptor which is only used
    ! with octree grids. Checks the allocation for the arrayes needed
    ! for octree grids, allocate if required.
    ! Checks the dimension and extends the arrays if required.
    type(DomainType), intent(inout) :: Domain
    integer::iError
    !--------------------------------------------------------------------------
    if(associated(Domain%iRoot_I))deallocate(Domain%iRoot_I)
    allocate(Domain%iRoot_I(product(Domain%iRootMapDim_D)), stat=iError)
    call check_allocate(iError,'iRoot_I'); Domain%iRoot_I = None_
  end subroutine allocate_iroot
  !============================================================================
  subroutine bcast_decomposition_dd(Domain)
    ! Broadcasts a given domain decomposition from the root PE
    ! via the global communicator iComm

    type (DomainType), intent(inout) :: Domain
    integer :: iComm, iProc0, iError
    !--------------------------------------------------------------------------
    if(Domain%IsLocal)then
       iComm = i_comm(Domain%CompID_); iProc0 = 0
    else
       iComm = i_comm(); iProc0 = i_proc0(Domain%CompID_)
    end if
    call MPI_Bcast(Domain%IsPeriodic_D(1), &
         Domain%nDim, MPI_LOGICAL, iProc0, iComm, iError)
    call MPI_Bcast(Domain%CoordMin_D(1),   &
         Domain%nDim, MPI_REAL, iProc0, iComm, iError)
    call MPI_Bcast(Domain%CoordMax_D(1),   &
         Domain%nDim, MPI_REAL, iProc0, iComm, iError)
    call MPI_Bcast(Domain%iRootMapDim_D(1),&
         Domain%nDim, MPI_INTEGER, iProc0, iComm, iError)
    call MPI_Bcast(Domain%nCell_D(1),     &
         Domain%nDim, MPI_INTEGER, iProc0, iComm, iError)
    if(Domain%IsTreeDD)then
       if(ubound(Domain%iRoot_I, 1) < product(Domain%iRootMapDim_D))&
         call allocate_iroot(Domain)
       call MPI_Bcast(Domain%iShift_DI(1,1),&
            Domain%nDim*Domain%nChildren, MPI_INTEGER, iProc0, iComm, iError)
    end if
    call bcast_indexes(Domain)
    call MPI_Bcast(Domain%DoGlueMargins, 1, MPI_LOGICAL, iProc0, iComm, iError)
    call MPI_Bcast(Domain%iDirMinusGlue, 1, MPI_INTEGER, iProc0, iComm, iError)
    call MPI_Bcast(Domain%iDirPlusGlue,  1, MPI_INTEGER, iProc0, iComm, iError)
    call MPI_Bcast(Domain%iDirCycle,     1, MPI_INTEGER, iProc0, iComm, iError)
  end subroutine bcast_decomposition_dd
  !============================================================================
  subroutine bcast_indexes( Domain, iProcUnion, iCommUnion)
    ! Broadcasts only DD\_II. This is a part of the
    ! synchronize\_refinement procedure, but it can be used separately
    ! too. If any of the optional parameters is not present, the
    ! DD\_II array is sent from the root processor of the
    ! components to all PEs in the global communicator, otherwise
    ! it itis sent FROM the PE having the rank iProcUnion in the
    ! communicator iCommUnion TO all PEs of this communicator.
    ! Recalculate local PE ranks to their values in the global
    ! communicator while broadcasting.

    type (DomainType), intent(inout) :: Domain
    integer, optional, intent(in)    :: iProcUnion,iCommUnion

    integer::iComm, iProc0, iError
    !--------------------------------------------------------------------------
    if(present(iProcUnion).and.present(iCommUnion))then
       iComm = iCommUnion; iProc0 = iProcUnion
    elseif(Domain%IsLocal)then
       iComm = i_comm(Domain%CompID_); iProc0 = 0
    else
       iComm = i_comm(); iProc0 = i_proc0(Domain%CompID_)
    end if

    call MPI_Bcast(Domain%nTreeNode, 1, MPI_INTEGER, iProc0, iComm,iError)
    call check_octree_allocation(Domain)

    if(.not.Domain%IsLocal .and. is_proc0(Domain%CompID_))then
       ! Recalculate local PE ranks to their values in the global
       ! communicator (at the root pe only)
       where(Domain%iDD_II(FirstChild_,1:Domain%nTreeNode) == None_&
            .and.Domain%iDD_II(PE_, 1:Domain%nTreeNode) /= None_)&
            Domain%iDD_II(PE_, 1:Domain%nTreeNode) = &
            i_proc0() + i_proc0(Domain%CompID_) + &
            Domain%iDD_II(PE_, 1:Domain%nTreeNode)*&
            i_proc_stride(Domain%CompID_)
    end if
    call MPI_Bcast(Domain%iDD_II(-1,1),(2+ubound(Domain%iDD_II,1))*&
         Domain%nTreeNode, MPI_INTEGER, iProc0, iComm, iError)
    call MPI_Bcast(Domain%iRealization, 1, MPI_INTEGER, iProc0, iComm, iError)
    call complete(Domain)
  end subroutine bcast_indexes
  !============================================================================
  subroutine complete(Domain)
    ! complete recovers the geometric variables in situ

    type (DomainType), intent(inout) :: Domain
    real    :: DCoordRoot_D(Domain%nDim)
    integer ::iTreeNode, lRoot, i, iDim, nDim
    integer ::iRootCounter, lParent, iChildNumber
    !--------------------------------------------------------------------------
    DCoordRoot_D = (Domain%CoordMax_D - Domain%CoordMin_D)/Domain%iRootMapDim_D
    nDim=Domain%nDim
    iRootCounter = 0
    do iTreeNode = 1, Domain%nTreeNode
       if(Domain%iDD_II(MyNumberAsAChild_, iTreeNode)==0)then
          iRootCounter = iRootCounter + 1
          Domain%DCoordBlock_DI(:, iTreeNode) = DCoordRoot_D
          Domain%DCoordCell_DI(:, iTreeNode)  = DCoordRoot_D/&
               Domain%nCell_D
          Domain%iDD_II(Parent_,iTreeNode)=iRootCounter
          lRoot=iRootCounter-1
          do iDim = 1, nDim
             i = mod(lRoot, Domain%iRootMapDim_D(iDim))
             Domain%CoordBlock_DI(iDim, iTreeNode) = &
                  Domain%CoordMin_D(iDim) + i*DCoordRoot_D(iDim)
             lRoot = (lRoot-i)/Domain%iRootMapDim_D(iDim)
          end do
          if(Domain%IsTreeDD)Domain%iRoot_I(iRootCounter) = iTreeNode
       else
          lParent = Domain%iDD_II(Parent_, iTreeNode)
          Domain%DCoordBlock_DI(:,iTreeNode)=&
               Domain%DCoordBlock_DI(:,lParent)
          Domain%DCoordBlock_DI(1:Domain%nDim, iTreeNode) = &
               0.5*Domain%DCoordBlock_DI(1:Domain%nDim, iTreeNode)
          Domain%DCoordCell_DI(:,iTreeNode)=&
               Domain%DCoordBlock_DI(:,iTreeNode)/Domain%nCell_D
          Domain%CoordBlock_DI(:,iTreeNode) = &
               Domain%CoordBlock_DI(:,lParent)
          iChildNumber = Domain%iDD_II(MyNumberAsAChild_,iTreeNode)
          Domain%CoordBlock_DI(1:Domain%nDim, iTreeNode) =           &
               Domain%CoordBlock_DI(1:Domain%nDim,iTreeNode) +       &
               Domain%DCoordBlock_DI(1:Domain%nDim, iTreeNode)*      &
               Domain%iShift_DI(:, iChildNumber)
       end if
    end do
  end subroutine complete
  !============================================================================
  subroutine synchronize_refinement_dd(&
       GlobalDomain, LocalDomain, iProcUnion, iCommUnion)
    !           SYNCHRONIZE LOCAL AND GLOBAL GRID
    !           NOTE: IF GridID\_ is used for global grid, then
    !           synchronize\_refinement is the only way to properly
    !           account for the refinement
    !
    ! If any of the optional parameters is not present, the global
    ! decomposition at all the PEs of the global communicator is
    ! synchronized with the local one at root processor of the
    ! component. Otherwise the global
    ! decomposition at all the PEs of the communicator iProcUnion is
    ! synchronized with the local one at the PE having the rank
    ! iProcUnion in the communicator iCommUnion.
    ! Recalculate local PE ranks (of the local grid) to their values
    ! in the global communicator (i\_comm()).

    type(DomainType),  intent(inout) :: GlobalDomain
    type(DomainType),  intent(in)    :: LocalDomain
    integer, optional, intent(in)    :: iProcUnion, iCommUnion
    integer :: iProc0, iComm, LocalIRealization, iError
    logical :: IsSynchronized
    !--------------------------------------------------------------------------
    if(present(iProcUnion).and.present(iCommUnion))then
       iProc0 = iProcUnion; iComm = iCommUnion
    else
       iProc0 = i_proc0(GlobalDomain%CompID_); iComm = i_comm()
    end if
    if(is_proc0(GlobalDomain%CompID_))&
         LocalIRealization = LocalDomain%iRealization

    call MPI_Bcast(LocalIRealization, 1, MPI_INTEGER, iProc0, iComm, iError)
    call MPI_Allreduce(LocalIRealization==GlobalDomain%iRealization,&
         IsSynchronized, 1, MPI_LOGICAL, MPI_LAND, iComm, iError)
    if(IsSynchronized)RETURN
    if(is_proc0(GlobalDomain%CompID_))then
       GlobalDomain%iRealization = LocalDomain%iRealization
       GlobalDomain%nTreeNode   = LocalDomain%nTreeNode
       call check_octree_allocation(GlobalDomain)
       GlobalDomain%iDD_II(:, 1:GlobalDomain%nTreeNode) = &
            LocalDomain%iDD_II(:, 1:LocalDomain%nTreeNode)
    end if
    call bcast_indexes(GlobalDomain, iProcUnion, iCommUnion)
  end subroutine synchronize_refinement_dd
  !============================================================================
  function is_left_boundary_d(Domain, iTreeNode)
    ! Returns the nDim vector, whose iDim's component is true, if
    ! along the iDim's direction there is a domain Boundary to the
    ! Left from the block with the number iTreeNode

    integer,           intent(in) :: iTreeNode
    type(DomainType),  intent(in) :: Domain
    logical  :: is_left_boundary_d(Domain%nDim)
    !--------------------------------------------------------------------------
    is_left_boundary_d = .not.(Domain%IsPeriodic_D).and.&
         Domain%CoordBlock_DI(:, iTreeNode)   < &
         cThird*Domain%DCoordBlock_DI(:,iTreeNode) + Domain%CoordMin_D
  end function is_left_boundary_d
  !============================================================================
  function is_right_boundary_d(Domain, iTreeNode)
    ! Returns the nDim vector, whose iDim's component is true, if
    ! along the iDim's direction there is a Tree Boundary to the
    ! Right from the block with the number iTreeNode

    integer,         intent(in) :: iTreeNode
    type(DomainType),intent(in) :: Domain
    logical :: is_right_boundary_d(Domain%nDim)
    !--------------------------------------------------------------------------
    is_right_boundary_d=.not.(Domain%IsPeriodic_D).and.              &
         Domain%CoordBlock_DI(:,iTreeNode) +                 &
         (1.0 + cThird)*Domain%DCoordBlock_DI(:,iTreeNode) > &
         Domain%CoordMax_D
  end function is_right_boundary_d
  !============================================================================
  integer function l_neighbor(Domain, iTreeNode, iCells_D)
    ! Tree neighbor returns the number of octree node, to which the center of
    ! "ghost cell" belongs, which is marked with iCells\_D cell index vector.

    type(DomainType), intent(inout) :: Domain
    integer,          intent(in)    :: iTreeNode
    integer,          intent(in)    :: iCells_D(Domain%nDim)

    real   :: Coord_D(Domain%nDim)
    !--------------------------------------------------------------------------
    if(any(is_left_boundary_d(Domain,iTreeNode).and.iCells_D<1).or.&
         any(is_right_boundary_d(Domain, iTreeNode).and.iCells_D > &
         Domain%nCell_D))then
       l_neighbor=None_
    else
       ! Cell center coordinates
       Coord_D = Domain%CoordBlock_DI(:,iTreeNode) + &
            Domain%DCoordCell_DI(:,iTreeNode)*(real(iCells_D) - 0.5)
       call search_in(Domain, Coord_D, l_neighbor)
    end if
  end function l_neighbor
  !============================================================================
  subroutine glue_margin(Domain, Coord_D)
    type(DomainType), intent(in)    :: Domain
    real,             intent(inout) :: Coord_D(Domain%nDim)
    ! Loop variable:
    integer :: iDir, nDim
    ! Directions to be glued
    integer :: iDirMinusGlue, iDirPlusGlue, iDirCycle
    ! Domain boundaries:
    real, pointer :: CoordMin_D(:), CoordMax_D(:)
    logical :: DoCycle
    !--------------------------------------------------------------------------
    nDim          = Domain%nDim
    CoordMin_D    => Domain%CoordMin_D
    CoordMax_D    => Domain%CoordMax_D
    iDirCycle     = Domain%iDirCycle
    iDirMinusGlue = Domain%iDirMinusGlue
    iDirPlusGlue  = Domain%iDirPlusGlue
    DoCycle = iDirCycle > 0.and.iDirCycle <= nDim
    do iDir = 1, nDim
       if(iDir == iDirCycle)CYCLE
       if(Coord_D(iDir) < CoordMin_D(iDir).and.iDir==iDirMinusGlue)then
          ! Ruturn the coordinate along the glued direction to the domain
          Coord_D(iDir) = 2*CoordMin_D(iDir) - Coord_D(iDir)
          ! Add, if needed a half cycle over cycled direction
          if(DoCycle)Coord_D(iDirCycle) = CoordMin_D(iDirCycle) + &
               modulo(Coord_D(iDirCycle) + 0.50*&
               CoordMax_D(iDirCycle) - 1.50*CoordMin_D(iDirCycle),&
               CoordMax_D(iDirCycle) - CoordMin_D(iDirCycle))
       elseif(Coord_D(iDir) > CoordMax_D(iDir).and.iDir==iDirPlusGlue)then
          ! Return the coordinate along the glued direction to the domain
          Coord_D(iDir) = 2*CoordMax_D(iDir) - Coord_D(iDir)
          ! Add, if needed a half cycle over cycled direction
          if(DoCycle)Coord_D(iDirCycle) = CoordMin_D(iDirCycle) + &
               modulo(Coord_D(iDirCycle) + 0.50*&
               CoordMax_D(iDirCycle) - 1.50*CoordMin_D(iDirCycle),&
               CoordMax_D(iDirCycle) - CoordMin_D(iDirCycle))
       end if
    end do
  end subroutine glue_margin
  !============================================================================
  subroutine search_in(Domain, Coord_D, iTreeNode)
    ! The searching tools start from here which allow to find the
    ! location in the domain decomposition, using the generalized
    !  coordinates of the point.
    !
    ! Returns iTreeNode, which is the  global number of the
    ! block which includes the point Coord_D.
    ! The starting values for Coord_D should be the generalized
    ! coordinates for the point, finally in this array there are the
    ! coordinates with respect to the left corner of the found block.
    ! The values of Coord_D are allowed to be outside of the domain,
    ! the nearest block is found in this case.

    type(DomainType), intent(inout) :: Domain
    real,             intent(inout) :: Coord_D(Domain%nDim)
    integer,          intent(out)   :: iTreeNode
    real    :: CoordTrunc_D(Domain%nDim), Discr_D(Domain%nDim)
    integer :: iChild, iDim, lFound, iShift_D(Domain%nDim), &
         iRootMinusOne_D(Domain%nDim), iRootMinusOneStart_D(Domain%nDim)
    !--------------------------------------------------------------------------
    ! Start from the result of the previous search
    lFound=Domain%lSearch

    ! Put the point inside the domain
    do iDim=1,Domain%nDim
       if(.not.Domain%IsPeriodic_D(iDim))CYCLE
       Coord_D(iDim) = Domain%CoordMin_D(iDim) +           &
            modulo(Coord_D(iDim) - Domain%CoordMin_D(iDim),&
            Domain%CoordMax_D(iDim) - Domain%CoordMin_D(iDim))
    end do
    CoordTrunc_D = Domain%CoordMin_D + max(min(Coord_D - Domain%CoordMin_D, &
         cAlmostOne*(Domain%CoordMax_D - Domain%CoordMin_D)), 0.0)

    Discr_D = (CoordTrunc_D - Domain%CoordBlock_DI(:,lFound))&
         /Domain%DCoordBlock_DI(:, lFound)
    ! Recursive search starts
    do
       ! correct Discr_D to avoid having points exactly
       ! on the boundary of a block; potentiall introduces error,
       ! but it would matter -log_2(cTol)~23 resolution level up the tree
       Discr_D = real(floor(Discr_D)) + &
            min(cAlmostOne, max(cTol2, Discr_D - real(floor(Discr_D)) ) )

       if(any(Discr_D < 0.0) .or. any(Discr_D >= 1.0))then
          iChild = Domain%iDD_II(MyNumberAsAChild_, lFound)
          if(iChild==0)then  ! This is a root
             iRootMinusOneStart_D = nint((Domain%CoordBlock_DI(:, lFound) - &
                  Domain%CoordMin_D)/Domain%DCoordBlock_DI(:,lFound))
             ! Calculate iRoot-1, jRoot-1....
             iRootMinusOne_D = floor(Discr_D) + iRootMinusOneStart_D
             Discr_D = Discr_D - floor(Discr_D)
             ! Calculate the root number, using the formula
             ! lRoot-1=iRoot-1+&
             !(jRoot-1)*RootMapDim(1)+&
             !(kRoot-1)*RootMapDim(1)*RootMapDim(2)
             lFound = 1 + iRootMinusOne_D(1)
             do iDim=2, Domain%nDim
                lFound = lFound + iRootMinusOne_D(iDim)*   &
                     product(Domain%iRootMapDim_D(1:iDim-1))
             end do
             if(Domain%IsTreeDD)lFound = Domain%iRoot_I(lFound)
          else     ! End of computations for root
             ! Descend the octree
             lFound = Domain%iDD_II(Parent_, lFound)
             Discr_D(1:Domain%nDim) = (Discr_D(1:Domain%nDim) + &
                  real(Domain%iShift_DI(:, iChild)))*0.50
          end if
       elseif(Domain%iDD_II(FirstChild_, lFound) == None_)then
          EXIT ! Octree node is found
       else
          ! Ascend the octree: calculate the shift
          Discr_D(1:Domain%nDim) = Discr_D(1:Domain%nDim)*2.0
          iShift_D           = int(Discr_D(1:Domain%nDim))
          Discr_D(1:Domain%nDim) = Discr_D(1:Domain%nDim) - real(iShift_D)
          ! Choose the child to ascend to
          iChild = 1
          do while(any(iShift_D/=Domain%iShift_DI(:,iChild)))
             iChild = iChild + 1
          end do
          lFound = Domain%iDD_II(iChild, lFound)
       end if
    end do
    ! End of recursive search
    Coord_D        = Coord_D - Domain%CoordBlock_DI(:,lFound)
    iTreeNode      = lFound
    Domain%lSearch = lFound
  end subroutine search_in
  !============================================================================
  subroutine search_cell(Domain, iTreeNode, Coord_D, iCells_D)
    type(DomainType), intent(in) :: Domain
    real,          intent(inout) :: Coord_D(Domain%nDim)
    integer,          intent(in) :: iTreeNode
    integer,         intent(out) :: iCells_D(Domain%nDim)
    real ::DCoordCells_D(Domain%nDim)
    !--------------------------------------------------------------------------
    DCoordCells_D = Domain%DCoordCell_DI(:, iTreeNode)
    iCells_D      = floor(Coord_D/DCoordCells_D)
    Coord_D       = Coord_D - DCoordCells_D*iCells_D
    iCells_D      = iCells_D + 1
  end subroutine search_cell
  !============================================================================
  subroutine pe_and_blk(Domain, iTreeNode, iPEOut, iBlockOut)
    type(DomainType),intent(in) :: Domain
    integer,         intent(in) :: iTreeNode
    integer,        intent(out) :: iPEOut,iBlockOut
    !--------------------------------------------------------------------------
    iPEOut    = Domain%iDD_II(PE_,  iTreeNode)
    iBlockOut = Domain%iDD_II(BLK_, iTreeNode)
  end subroutine pe_and_blk
  !============================================================================
  subroutine associate_dd_pointer_dd(Domain, DomainPointer)
    Type(DomainType), target, intent(in)  :: Domain
    type(DomainPointerType),  intent(out) :: DomainPointer
    !--------------------------------------------------------------------------
    nullify(DomainPointer%Ptr); DomainPointer%Ptr=>Domain
  end subroutine associate_dd_pointer_dd
  !============================================================================

end module CON_domain_decomposition
!==============================================================================
