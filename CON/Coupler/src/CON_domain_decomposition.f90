!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!
!
!
!
module CON_domain_decomposition

  ! This file presents the class of the domain decompositions which
  ! icludes both the uniformly spaced ones (uniformly spaced with
  ! respect to some generalized coordinates) and Octree or Quadric
  ! tree for adaptive block decompositions.
  
  use ModUtilities, ONLY: check_allocate
  use CON_world
  use ModNumConst
  use ModMpi
  use ModKind, ONLY: nByteReal

  ! revision history:
  ! 6/18/03-7/11/03 Sokolov I.V. <igorsok@umich.edu> phone(734)647-4705

  implicit none
  SAVE
  real, parameter :: cThird = 1.0/3.0

  integer,parameter::Parent_=-1,&
       MyNumberAsAChild_ =0,&
       FirstChild_  =1,&
       BLK_         =2,&
       PE_          =3,&
       GlobalBlock_ =4,&
       None_        =-777
  !\begin{verbatim}
  !=====================DERIVED TYPE==============================
  type DomainType
     ! The type is to describe the domain decomposition for a
     ! component. The component should be properly registered with
     ! registry procedures and should have a unique ID

     integer :: CompID_

     ! This type describes the general decomposition for rectangular
     ! box domain.
     ! Below ":=" means "according to definition".
     ! Describe the property of the ! domain first.
     !=============================DOMAIN============================
     ! genaralized coordinates := the coordinates with respect to
     !    which the domain decomposition chunks are uniform.
     ! nDim := the dimensionality of the generalized coordinate space

     integer ::nDim

     !---------------------------------------------------------------
     ! domain := (nDim)-dimensional rectangular box, such as the
     ! data used by the component CompID_ can be associated with the
     ! points inside of this domain
     !
     ! CoordMin_D:=Coordinates of the left corner of the domain
     ! CoordMax_D:=Coordinates of the right corner of the computational
     ! left corner:=the common point of all left boundary faces,
     ! each of the nDim left boundary faces being defined with
     ! respect to the correspondent direction

     real, pointer::CoordMin_D(:), CoordMax_D(:)

     !================ROOT DECOMPOSITION=============================
     ! Root decomposition:=  decomposition of the domain for
     ! iRootMapDim_D(1)*iRootMamDim_D(2)*.....iRootMapDim_D(nDim)
     ! blocks, which are spatially uniform in geralized coordinates
     ! block:=rectangular box.
     !
     !
     ! The principle "block is not shared" is assumed, that is if the
     ! physical data at any point geometrically belonging to the block
     ! are allocated at the given PE, then all the points in the
     ! Component, which geometrically belong to this block, should be
     ! also allocated at the same PE.
     !
     !(For the uniformly spaced grids the root decomposition is the
     ! map of blocks, for tree grid - the map of the roots)

     integer, pointer :: iRootMapDim_D(:)

     !================CELL DECOMPOSITION=============================
     ! Cell decomposition:=decomposition of each of the blocks for
     ! nCells_D(1)*nCells_D(2)*.....nCells_D(nDim) cells,
     ! which are spatially uniform in the geralized coordinates
     !
     ! So, it is assumed that the domain is fully decomposed for
     !"blocks" and just in an analogous manner we assume that the
     ! region covered by blocks can be decomposed for "cells".
     ! Searching programs given here allow to find the number of block
     ! point and the vector of the cell numbers, for the cell to which
     ! with coordinates x,(y,z) belongs. The availbability of these
     ! programs still DOES NOT require to use these cells for the
     ! control volume method. The actual grid points can be chosen to
     ! be  at the vertexes, faces, edges and so on. See the
     ! explanations for grid descriptors in CON_router
     ! Nevertheless usually the indexes of cell the point belongs to
     ! allow to find readily the number of the nearest node, nearest
     ! edge and so on.

     integer, pointer :: nCells_D(:)

     !============TREE DECOMPOSITION=================================
     !RECURSIVE DEFINITION:
     ! tree decomposition:= root decomposition .or.&
     !     refinement[tree decomposition]
     ! refinement:=decompose any block for 2*nDimTree (nDimTree<=nDim)
     ! blocks by splitting the original block for two equal parts
     ! along each of the first nDimTree dimensions
     !
     ! The mathematical space of tree decompositions, for a given
     ! iRootMapDim_D can be one-to-one mapped to the space of
     ! structures, consisting of product(iRootMapDim_D) tree graphs.
     ! Each node of this tree corresponds to either a root map
     ! element, or the block, which was refined, or a block.
     !
     ! If the tree decomposition is allowed for a given domain, set
     ! IsTreeGrid=.true.

     logical :: IsTreeDecomposition

     ! Octree or Quadric tree. Is only used for the decomposition with
     ! IsTreeDecomposition=.true.

     integer :: nDimTree

     ! nChildren=2**nDimTree for the tree grid

     integer :: nChildren

     !====================DESCRIPTOR FOR TREE DECOMPOSITIONS=========
     ! For tree decomposition this is a full tree structure saved as
     ! described in the figure (PE is a CPU rank in a i_comm(Global_)

     integer, pointer::iDecomposition_II(:,:)

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
     !
     !---------------------------------------------------------------
     ! For non-tree grid: each block is both a root node and an end
     ! node.
     !
     ! For each lGlobalTreeNumber (second index) the vector
     ! iDecomposition_II(:,lGlobalTreeNumber)
     ! is
     !(lGlobalTreeNumber(Root Number),&
     ! 0 (MyNumberAsAChild),&
     ! None_(FirstChild_),&
     ! PE_(lGlobalTreeNumber),&
     ! BLK_(lGlobalTreeNumber).
     !
     ! In this case iDecomposition_II is a root decomposition,
     ! which is accompanied by some markers, allowing to apply the
     ! same searching programs for tree- and non-tree- grids
     !
     !---------------------------------------------------------------
     ! nTreeNodes - the number of blocks,
     ! for non-tree grid, or the number of all nodes
     ! in the tree.

     integer::nTreeNodes, nAllocatedNodes

     !---------------------------------------------------------------
     ! The following is not allocated if the grid is a tree grid.
     ! The GlobalTree of the tree roots in the
     ! iDecomposition_II array (the values for the second index
     integer, pointer :: iRoot_I(:)
     !
     ! Introduce the array nDimTree*nChildren, which defines the
     ! spacial shift vector of the Child's zero corner with respect to
     ! the parent zero corner, measured in the Child block size units.
     ! This shifht is unambigous function of the child number
     !
     ! In index I  passing from 1 to nChildren, all the possible
     ! combinations of 0 and 1 of dimension nDimTree should be
     ! unambigously related to the values of I In the grid captured by
     ! the control module, this array should be organized exaclty as
     ! it done in the "component".

     integer, pointer :: iShift_DI(:, :)

     !===============GEOMETRY CHARACTERISTICS========================
     ! In principle, the array introduced above contains all the
     ! geometric characteristics, however it is not convenient to
     ! calculate it each time. Tree more arrays are stored for each
     ! domain decomposition
     !
     ! CoordBlock_DI are the generalized coordinates of the left corner
     ! of the block,

     real, pointer :: CoordBlock_DI(:, :)

     ! DCoordBlock_D is the BLOCK size

     real, pointer :: DCoordBlock_DI(:,:)

     ! and DCoordCell_D is the cell size.

     real, pointer :: DCoordCell_DI(:,:)

     ! First index inumerates the spatial component of the
     ! generalized coordinates, second is lGlobalNode

     !==========SPECIAL BOUNDARIES===================================
     ! Sometimes the coordinate values which are arithmetically beyond
     ! the limit values CoordMin_D, CoordMax_D, geometrically appear to be
     ! inside the domain. The simplest example is given by periodic
     ! boundary conditions. For the domain with periodic boundary
     ! conditions the following array is introduced

     logical, pointer::IsPeriodic_D(:)

     ! For more complicated cases (spherical and analogous geometries
     ! having a pole singularity) a special procedure
     ! glue_margin is invoked when required. To do this
     ! for a particular grid, it should have DoGlueMargins=.true.

     logical :: DoGlueMargins
     integer :: iDirMinusGlue, iDirPlusGlue, iDirCycle
     !==========SOME INTEGERS========================================
     ! To search in the tree decompositions  it is optimal to start
     ! from the previously found node (from the previous search). To
     ! do this we have this previuosly found value saved in lSearch

     integer :: lSearch

     ! For the grids which can be refined, the realization number for
     ! the domain decomposition should be stored

     integer :: iRealization

     !===================LOCAL AND GLOBAL DECOMPOSITIONS=============
     ! In the local decomposition the PE ranks are defined to be in a
     ! local communicator, for a particular component

     logical :: IsLocal

     ! The derived type described below consists of decompositions,
     ! for which IsLocal=.false., but if the decomposition allows
     ! refinement, there must be a local decompositions in the
     ! correspondent model which accounts for any refinements as soon
     ! as they are done. The changes in iDecomoposition_II array then
     ! can come to the array DD_I only through the procedure
     ! synchronize
     !
     ! The inverse mapping (iPE,iBlock)->global node number is mantained
     ! using the following components:

     integer :: MinBlock, MaxBlock
     integer, pointer::iGlobal_BP(:,:)

     ! The inverse mapping GlobalBlock number->global node number
     ! is mantained using the following components:

     integer::nBlockAll
     integer, pointer :: iGlobal_A(:)

  end type DomainType
  !================DERIVED TYPE===================================
  type DomainPointerType
     type(DomainType),pointer::Ptr
  end type DomainPointerType

  ! Needed for searching and interpolating algorithms
  real, parameter:: &
       cTol  = 0.00000010, &
       cTol2 = cTol**(nByteReal/4)
  real, parameter:: &
       cAlmostOne   = 1.0 - cTol2,&
       cAlmostTwo   = 2*cAlmostOne,&
       cAlmostHalf  = 0.5*cAlmostOne
  private:: is_used_block_dd
contains
  !============================================================================
  subroutine init_decomposition_dd(&
       Domain,&
       CompID_,&            ! As in DomainType
       nDim,&               ! As in DomainType
       IsTreeDecomposition,&! As in DomainType
       nDimTree, IsLocal)

    type(DomainType),intent(inout) :: Domain
    integer, intent(in)            :: CompID_, nDim
    logical, intent(in), optional  :: IsTreeDecomposition
    integer, intent(in), optional  :: nDimTree
    logical, intent(in), optional  :: IsLocal
    integer::iError,iDim,iChildFirst,iChildLast

    !--------------------------------------------------------------------------
    if(use_comp(CompID_))then
       Domain%CompID_=CompID_
    else
       call CON_stop(&
            'Unauthorized use of the component')
    end if
    Domain%IsTreeDecomposition=.false.
    if(present(IsTreeDecomposition))&
         Domain%IsTreeDecomposition=IsTreeDecomposition
       
    Domain%nDim=nDim
    if(Domain%IsTreeDecomposition)then
       if(.not.present(nDimTree))then
          Domain%nDimTree=nDim
       else
          Domain%nDimTree=nDimTree
       end if
       Domain%nChildren=&
            2**Domain%nDimTree
    else
       Domain%nChildren=1
    end if

    nullify(Domain%iRootMapDim_D)
    allocate(Domain%iRootMapDim_D(nDim),stat=iError)
    call check_allocate(iError,"iRootMapDim_D")
    Domain%iRootMapDim_D = 1

    nullify(Domain%CoordMin_D)
    allocate(Domain%CoordMin_D(nDim),stat=iError)
    call check_allocate(iError,"CoordMin_D")
    Domain%CoordMin_D = 0.0

    nullify(Domain%CoordMax_D)
    allocate(Domain%CoordMax_D(nDim),stat=iError)
    call check_allocate(iError,"CoordMax_D")
    Domain%CoordMax_D = 1.0

    nullify(Domain%nCells_D)
    allocate(Domain%nCells_D(nDim),stat=iError)
    call check_allocate(iError,"nCells_D")
    Domain%nCells_D = 1

    nullify(Domain%IsPeriodic_D)
    allocate(Domain%IsPeriodic_D(nDim),stat=iError)
    call check_allocate(iError,"IsPeriodic_D")
    Domain%IsPeriodic_D  = .false.
    Domain%DoGlueMargins = .false.
    Domain%iDirMinusGlue = 0
    Domain%iDirPlusGlue  = 0
    Domain%iDirCycle     = 0

    if(Domain%IsTreeDecomposition)then
       nullify(Domain%iRoot_I)
       call allocate_iroot(Domain)

       nullify(Domain%iShift_DI)
       allocate(Domain%iShift_DI&
            (Domain%nDimTree,&
            Domain%nChildren),stat=iError)
       call check_allocate(iError,"iShift_DI")
       Domain%iShift_DI=0
       ! A "binary" order for iShift_DI is set as a default:
       !(0,0,0,&
       ! 1,0,0,&
       ! 0,1,0,&
       ! and so on
       !
       iChildLast=1
       do iDim=1,Domain%nDimTree
          iChildFirst=iChildLast+1
          iChildLast=min(&
               Domain%nChildren,2*iChildLast)
          if(iChildFirst<iChildLast)EXIT
          Domain%iShift_DI(:,iChildFirst:iChildLast)=&
               Domain%iShift_DI(:,1:1+iChildLast-iChildFirst)
          Domain%iShift_DI(iDim,iChildFirst:iChildLast) = 1
       end do
    end if

    Domain%IsLocal=.false.
    if(present(IsLocal)) Domain%IsLocal = IsLocal
   
    Domain%nTreeNodes=1
    Domain%nAllocatedNodes=-1

    nullify(Domain%iDecomposition_II)
    nullify(Domain%CoordBlock_DI)
    nullify(Domain%DCoordBlock_DI)
    nullify(Domain%DCoordCell_DI)
    call check_octree_allocation(Domain)
    Domain%lSearch = 1
    Domain%iRealization = 0
    nullify(Domain%iGlobal_BP)
    Domain%MinBlock=1
    Domain%MaxBlock=1
    if(Domain%IsLocal)then
       allocate(Domain%iGlobal_BP(1:1, 0:-1+n_proc(Domain%CompID_)),&
            stat=iError)
    else
       allocate(Domain%iGlobal_BP(1:1, 0:i_proc_last(CompID_)),&
            stat=iError)
    end if
    call check_allocate(iError,'iGlobal_BP,first allocation')
    Domain%nBlockAll=1
    nullify(Domain%iGlobal_A)
    allocate(Domain%iGlobal_A(1:Domain%nBlockAll),stat=iError)
    call check_allocate(iError,'iGlobal_A,first allocation')
  end subroutine init_decomposition_dd
  !============================================================================

  ! Checks the array allocation for a given grid desciptor and
  ! extends its dimension if required

  subroutine check_octree_allocation(Domain)
    type(DomainType),intent(inout)::&
         Domain
    integer::iError,nUbound

    !--------------------------------------------------------------------------
    nUbound= GlobalBlock_
    if(Domain%IsTreeDecomposition)&
         nUbound=max(nUbound, Domain%nChildren)
    if(Domain%nAllocatedNodes>=Domain%nTreeNodes)RETURN
    if(associated(Domain%iDecomposition_II))then
       deallocate(Domain%iDecomposition_II,stat=iError)
       deallocate(Domain%CoordBlock_DI,stat=iError)
       deallocate(Domain%DCoordBlock_DI,stat=iError)
       deallocate(Domain%DCoordCell_DI,stat=iError)
    end if

    nullify(Domain%iDecomposition_II)
    allocate(Domain%iDecomposition_II&
         (Parent_:nUBound, Domain%nTreeNodes), stat=iError)
    call check_allocate(iError, 'iDecomposition_II')
    Domain%iDecomposition_II=None_
    Domain%iDecomposition_II(MyNumberAsAChild_,:)=0

    nullify(Domain%CoordBlock_DI)
    allocate(Domain%CoordBlock_DI(Domain%nDim, Domain%nTreeNodes)&
         ,stat=iError)
    call  check_allocate(iError,'CoordBlock_DI')
    Domain% CoordBlock_DI=0.0

    nullify(Domain%DCoordBlock_DI)
    allocate(Domain%DCoordBlock_DI&
         (Domain%nDim,Domain%nTreeNodes),stat=iError)
    call  check_allocate(iError,'DCoordBlock_DI')
    Domain%DCoordBlock_DI = 1.0

    nullify(Domain%DCoordCell_DI)
    allocate(Domain%DCoordCell_DI(Domain%nDim, Domain%nTreeNodes)&
         ,stat=iError)
    call  check_allocate(iError,'DCoordCell_DI')
    Domain%DCoordCell_DI=1.0

    Domain%nAllocatedNodes=&
         Domain%nTreeNodes
  end subroutine check_octree_allocation
  !============================================================================
  ! To get a decomposition domain, even the tree one, the root     
  ! decomposition should be first constructed                      
  !                        ATTENTION
  ! PE here are the ranks in the LOCAL communicator for the        
  ! component
  

  subroutine get_root_decomposition_dd(&
       Domain,&! Decomposition to be constructed
       iRootMapDim_D,&! As in DomainType
       CoordMin_D,&     ! As in DomainType
       CoordMax_D,&     ! As in DomainType
       nCells_D,&     ! As in DomainType
       PE_I,&         ! PE layout
       iBlock_I,&     ! Local Block Number layout
       IsPeriodic_D,& ! As in DomainType
       iShift_DI,   & ! As in DomainType
       DoGlueMargins,& ! As in DomainType
       iDirMinusGlue,&! As in DomainType
       iDirPlusGlue,& ! As in DomainType
       iDirCycle)     ! As in DomainType
    type(DomainType),intent(inout)::&
         Domain
    !--------------------------------------------------------------------------
    integer,dimension(Domain%nDim),&
         intent(in)::iRootMapDim_D
    real,dimension(Domain%nDim),&
         intent(in)::CoordMin_D
    real,dimension(Domain%nDim),&
         intent(in)::CoordMax_D
    integer,dimension(Domain%nDim),&
         intent(in)::nCells_D
    integer,dimension(:),&
         intent(in),optional:: PE_I,iBlock_I
    logical,dimension(Domain%nDim),&
         intent(in),optional::IsPeriodic_D
    integer,optional,&
         dimension(Domain%nDim,&
         Domain%nChildren),&
         intent(in)::iShift_DI
    logical, intent(in), optional :: DoGlueMargins
    integer, intent(in), optional :: iDirMinusGlue
    integer, intent(in), optional :: iDirPlusGlue
    integer, intent(in), optional :: iDirCycle
    !-------
    integer::lBlock,MaxBlock
    !---------------------------------------------------------------!
    ! Check the dimension of PE_I and iBlock_I                       !
    if(present(PE_I))then
       if(ubound(PE_I,1)/=product(iRootMapDim_D))then
          call CON_stop('The dimension of PE_I is wrong')
       end if
    end if

    if(present(iBlock_I))then
       if(ubound(iBlock_I,1)/=product(iRootMapDim_D))&
            call CON_stop('The dimension of iBlock_I is wrong')
    end if

    !---------------------------------------------------------------!
    ! Assign the grid parameters                                     !

    Domain%iRootMapDim_D=iRootMapDim_D

    Domain%CoordMin_D=CoordMin_D
    Domain%CoordMax_D=CoordMax_D
    Domain%nCells_D=nCells_D
    if(present(DoGlueMargins)) &
         Domain%DoGlueMargins = DoGlueMargins
    if(present(iDirMinusGlue))&
         Domain%iDirMinusGlue = iDirMinusGlue
    if(present(iDirPlusGlue)) &
         Domain%iDirPlusGlue  = iDirPlusGlue
    if(present(iDirCycle))    &
         Domain%iDirCycle     = iDirCycle

    if(present(IsPeriodic_D))&
         Domain%IsPeriodic_D=IsPeriodic_D

    Domain%nTreeNodes=product(iRootMapDim_D)
    call check_octree_allocation(Domain)

    Domain%iDecomposition_II&
         (MyNumberAsAChild_,1:Domain%nTreeNodes)=0
    Domain%iDecomposition_II&
         (FirstChild_,1:Domain%nTreeNodes)=None_
    do lBlock=1,Domain%nTreeNodes
       Domain%iDecomposition_II(Parent_,lBlock)=lBlock
    end do

    if(Domain%IsTreeDecomposition)then
       call check_iroot_allocation(Domain)
       if(present(iShift_DI))Domain%iShift_DI=iShift_DI
    end if

    ! If neither PE\_I nor iBlock\_I is present, the root
    ! decomposition blocks are balanced for an optimal load
    if((.not.present(PE_I)).and.(.not.present(iBlock_I)))then
       MaxBlock=(Domain%nTreeNodes - 1)/&
            n_proc(Domain%CompID_) + 1
       do lBlock=1, Domain%nTreeNodes
          Domain%iDecomposition_II(BLK_, lBlock) = &
               mod(lBlock - 1, MaxBlock) + 1
          Domain%iDecomposition_II(PE_, lBlock) =  &
               (lBlock - 1)/MaxBlock
       end do
       call set_iglobal_and_bp_dd(Domain)
       RETURN
    end if
    ! If the PE\_I is not given, it is assumed that all the blocks
    ! are at the root PE
    if(present(PE_I))then
       Domain%iDecomposition_II&
            (PE_,1:Domain%nTreeNodes)=PE_I
    else
       Domain%iDecomposition_II&
            (PE_,1:Domain%nTreeNodes)=0
    end if

    if(present(iBlock_I))then
       Domain%iDecomposition_II(BLK_, 1:Domain%nTreeNodes) = iBlock_I
    else
       do lBlock=1, Domain%nTreeNodes
          Domain%iDecomposition_II(BLK_, lBlock) = lBlock
       end do
    end if
    call set_iglobal_and_bp_dd(Domain)
  end subroutine get_root_decomposition_dd
  !============================================================================
  subroutine set_iglobal_and_bp_dd(Domain)
    type(DomainType),intent(inout)::&
         Domain
    integer::iUpper_I(2),iLower_I(2),iError,iPE,iBlock
    integer::iGlobalNode,iGlobalBlock
    !--------------------------------------------------------------------------
    iLower_I=lbound(Domain%iGlobal_BP)
    iUpper_I=ubound(Domain%iGlobal_BP)
    Domain%MinBlock=minval(&
         Domain%iDecomposition_II(&
         BLK_,1:Domain%nTreeNodes),MASK=&
         Domain%iDecomposition_II(&
         FirstChild_,1:Domain%nTreeNodes)==None_&
         .and.Domain%iDecomposition_II(&
         PE_,1:Domain%nTreeNodes)/=None_)
    Domain%MaxBlock=maxval(&
         Domain%iDecomposition_II(&
         BLK_,1:Domain%nTreeNodes),MASK=&
         Domain%iDecomposition_II(&
         FirstChild_,1:Domain%nTreeNodes)==None_&
         .and.Domain%iDecomposition_II(&
         PE_,1:Domain%nTreeNodes)/=None_)
    Domain%nBlockAll=count(&
         Domain%iDecomposition_II(&
         FirstChild_,1:Domain%nTreeNodes)==None_&
         .and.Domain%iDecomposition_II(&
         PE_,1:Domain%nTreeNodes)/=None_)
    if(iUpper_I(1)<Domain%MaxBlock.or.&
         iLower_I(1)>Domain%MinBlock)then
       deallocate(Domain%iGlobal_BP)
       allocate(Domain%iGlobal_BP(&
            Domain%MinBlock:&
            Domain%MaxBlock,&
            iLower_I(2):iUpper_I(2)),&
            stat=iError)
       call check_allocate(iError,'iGlobal_BP - reallocate')
    end if
    if(ubound(Domain%iGlobal_A,1)<&
         Domain%nBlockAll)then
       deallocate(Domain%iGlobal_A)
       allocate(Domain%iGlobal_A(&
            1:Domain%nBlockAll),&
            stat=iError)
       call check_allocate(iError,'iGlobal_BP - reallocate')
    end if
    Domain%iGlobal_BP=None_
    Domain%iGlobal_A=None_
    iGlobalBlock=0
    do iGlobalNode=1,Domain%nTreeNodes
       if(.not.is_used_block_dd(Domain,iGlobalNode))&
            CYCLE
       call pe_and_blk_dd(Domain,iGlobalNode,&
            iPE,iBlock)
       Domain%iGlobal_BP(iBlock,iPE)=iGlobalNode
       iGlobalBlock=iGlobalBlock+1
       Domain%iDecomposition_II(GlobalBlock_,&
            iGlobalNode)=iGlobalBlock
       Domain%iGlobal_A(iGlobalBlock)=iGlobalNode
    end do

  end subroutine set_iglobal_and_bp_dd
  !============================================================================
 

  ! Allocates the part of the grid descriptor which is only used
  ! with octree grids. Checks the allocation for the arrayes needed
  ! for octree grids, allocate if required.
  ! Checks the dimension and extends the arrays if required.
  !---------------------------------------------------------------!
  subroutine allocate_iroot(Domain)
    type(DomainType),intent(inout)::&
         Domain
    integer::iError
    !--------------------------------------------------------------------------
    if(associated(Domain%iRoot_I))then
       deallocate(Domain%iRoot_I)
    end if
    allocate(&
         Domain%iRoot_I(product(&
         Domain%iRootMapDim_D)),&
         stat=iError)
    call check_allocate(iError,'iRoot_I')
    Domain%iRoot_I=None_
  end subroutine allocate_iroot
  !============================================================================
  subroutine check_iroot_allocation(Domain)
    type(DomainType),intent(inout)::&
         Domain
    !--------------------------------------------------------------------------
    if(ubound(Domain%iRoot_I,1)<product(&
         Domain%iRootMapDim_D))&
         call allocate_iroot(Domain)
  end subroutine check_iroot_allocation
  !============================================================================

  ! Broadcasts a given domain decomposition from the root PE
  ! via the global communicator iComm

  subroutine bcast_decomposition_dd(Domain)

    type (DomainType),intent(inout)::&
         Domain
    integer::iComm
    integer::iProc0
    integer::iError

    !--------------------------------------------------------------------------
    if(Domain%IsLocal)then
       iComm=i_comm(Domain%CompID_)
       iProc0=0
    else
       iComm=i_comm()
       iProc0=i_proc0(Domain%CompID_)
    end if
    call MPI_Bcast(Domain%IsPeriodic_D(1),&
         Domain%nDim,MPI_LOGICAL,&
         iProc0,iComm,iError)

    call MPI_Bcast(Domain%CoordMin_D(1),&
         Domain%nDim,MPI_REAL,&
         iProc0,iComm,iError)

    call MPI_Bcast(Domain%CoordMax_D(1),&
         Domain%nDim,MPI_REAL,&
         iProc0,iComm,iError)

    call MPI_Bcast(Domain%iRootMapDim_D(1),&
         Domain%nDim,MPI_INTEGER,&
         iProc0,iComm,iError)

    call MPI_Bcast(Domain%nCells_D(1),&
         Domain%nDim,MPI_INTEGER,&
         iProc0,iComm,iError)

    if(Domain%IsTreeDecomposition)then
       call check_iroot_allocation(Domain)
       call MPI_Bcast(Domain%iShift_DI(1,1),&
            Domain%nDim*&
            Domain%nChildren,&
            MPI_INTEGER,&
            iProc0, iComm,iError)
    end if
    call bcast_indexes(&
         Domain)

    call MPI_Bcast(Domain%DoGlueMargins,&
         1,MPI_LOGICAL,iProc0,iComm,iError)

    call MPI_Bcast(Domain%iDirMinusGlue,&
         1,MPI_INTEGER,iProc0,iComm,iError)

    call MPI_Bcast(Domain%iDirPlusGlue,&
         1,MPI_INTEGER,iProc0,iComm,iError)

    call MPI_Bcast(Domain%iDirCycle,&
         1,MPI_INTEGER,iProc0,iComm,iError)

  end subroutine bcast_decomposition_dd
  !============================================================================
  ! Broadcasts only Decomposition\_II. This is a part of the
  ! synchronize\_refinement procedure, but it can be used separately
  ! too. If any of the optional parameters is not present, the
  ! Decomposition\_II array is sent from the root processor of the
  ! components to all PEs in the global communicator, otherwise
  ! it itis sent FROM the PE having the rank iProcUnion in the
  ! communicator iCommUnion TO all PEs of this communicator.
  ! Recalculate local PE ranks to their values in the global
  ! communicator while broadcasting.

  subroutine bcast_indexes(&
       Domain, iProcUnion, iCommUnion)
    type (DomainType),intent(inout)::&
         Domain
    integer,optional,intent(in)::iProcUnion,iCommUnion

    integer::iComm
    integer::iProc0
    integer::iError
    !--------------------------------------------------------------------------
    if(present(iProcUnion).and.present(iCommUnion))then
       iComm=iCommUnion
       iProc0=iProcUnion
    elseif(Domain%IsLocal)then
       iComm=i_comm(Domain%CompID_)
       iProc0=0
    else
       iComm=i_comm()
       iProc0=i_proc0(Domain%CompID_)
    end if

    call MPI_Bcast(Domain%nTreeNodes,&
         1,MPI_INTEGER,&
         iProc0, iComm,iError)
    call check_octree_allocation(Domain)

    if(.not.Domain%IsLocal.and.&
         is_proc0(Domain%CompID_))then
       ! Recalculate local PE ranks to their values in the global      !
       ! communicator (at the root pe only)                            !
       where(Domain%iDecomposition_II(&
            FirstChild_,1:Domain%nTreeNodes)&
            ==None_&
            .and.Domain%iDecomposition_II(&
            PE_,1:Domain%nTreeNodes)/=None_)&
            Domain%iDecomposition_II(&
            PE_,1:Domain%nTreeNodes)=&
            i_proc0()+i_proc0(Domain%CompID_)+&
            Domain%iDecomposition_II(&
            PE_,1:Domain%nTreeNodes)*&
            i_proc_stride(Domain%CompID_)
    end if

    call MPI_Bcast(Domain%iDecomposition_II(-1,1),&
         (2+ubound(Domain%iDecomposition_II,1))*&
         Domain%nTreeNodes,&
         MPI_INTEGER,iProc0,iComm,iError)

    call MPI_Bcast(Domain%iRealization,&
         1,MPI_INTEGER,&
         iProc0, iComm,iError)

    call complete(Domain)
  end subroutine bcast_indexes
  !============================================================================
  ! complete recovers the geometric variables in situ         

  subroutine complete(Domain)
    type (DomainType),intent(inout)::&
         Domain
    real,dimension(Domain%nDim)::DCoordRoot_D
    integer::lGlobalTreeNumber,lRoot,i,iDim,nDim
    integer::iRootCounter,lParent,iChildNumber

    !--------------------------------------------------------------------------
    DCoordRoot_D=(Domain%CoordMax_D-&
         Domain%CoordMin_D)/&
         Domain%iRootMapDim_D
    nDim=Domain%nDim

    iRootCounter=0
    do lGlobalTreeNumber=1,Domain%nTreeNodes
       if(Domain%iDecomposition_II(&
            MyNumberAsAChild_,lGlobalTreeNumber)==0)then
          iRootCounter=iRootCounter+1
          Domain%DCoordBlock_DI(&
               :,lGlobalTreeNumber)=DCoordRoot_D
          Domain%DCoordCell_DI(&
               :,lGlobalTreeNumber)=DCoordRoot_D/&
               Domain%nCells_D
          Domain%iDecomposition_II(&
               Parent_,lGlobalTreeNumber)=iRootCounter
          lRoot=iRootCounter-1
          do iDim=1,nDim
             i=mod(lRoot,&
                  Domain%iRootMapDim_D(iDim))
             Domain%CoordBlock_DI(&
                  iDim,lGlobalTreeNumber)=&
                  Domain%CoordMin_D(iDim)+&
                  i*DCoordRoot_D(iDim)
             lRoot=(lRoot-i)/&
                  Domain%iRootMapDim_D(iDim)
          end do
          if(Domain%IsTreeDecomposition)&
               Domain%iRoot_I(&
               iRootCounter)=lGlobalTreeNumber
       else
          lParent=Domain%iDecomposition_II(&
               Parent_,lGlobalTreeNumber)
          Domain%DCoordBlock_DI(&
               :,lGlobalTreeNumber)=&
               Domain%DCoordBlock_DI(:,lParent)
          Domain%DCoordBlock_DI(&
               1:Domain%nDimTree,lGlobalTreeNumber&
               )=0.5*Domain%DCoordBlock_DI(&
               1:Domain%nDimTree,lGlobalTreeNumber)
          Domain%DCoordCell_DI(&
               :,lGlobalTreeNumber)=&
               Domain%DCoordBlock_DI(&
               :,lGlobalTreeNumber)/&
               Domain%nCells_D
          Domain%CoordBlock_DI(:,lGlobalTreeNumber)&
               =Domain%CoordBlock_DI(:,lParent)
          iChildNumber=&
               Domain%iDecomposition_II(&
               MyNumberAsAChild_,lGlobalTreeNumber)
          Domain%CoordBlock_DI(&
               1:Domain%nDimTree,lGlobalTreeNumber)&
               =Domain%CoordBlock_DI(&
               1:Domain%nDimTree,lGlobalTreeNumber)&
               +Domain%DCoordBlock_DI(&
               1:Domain%nDimTree,lGlobalTreeNumber)&
               *Domain%iShift_DI(:,iChildNumber)
       end if
    end do
    call set_iglobal_and_bp_dd(Domain)
  end subroutine complete
  !============================================================================

  !\begin{verbatim}
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
  !\end{verbatim}

  subroutine synchronize_refinement_dd(&
       GlobalDomain,LocalDomain,iProcUnion,iCommUnion)

    type(DomainType),intent(inout)::GlobalDomain
    !--------------------------------------------------------------------------
    type(DomainType),&
         intent(in)::LocalDomain
    integer,intent(in),optional::iProcUnion,iCommUnion
    integer::iProc0,iComm,LocalIRealization,iError
    logical::IsSynchronized
    if(present(iProcUnion).and.present(iCommUnion))then
       iProc0=iProcUnion
       iComm=iCommUnion
    else
       iProc0=i_proc0(GlobalDomain%CompID_)
       iComm=i_comm()
    end if
    if(is_proc0(GlobalDomain%CompID_))&
         LocalIRealization=LocalDomain%iRealization

    call MPI_Bcast(LocalIRealization,&
         1,MPI_INTEGER,&
         iProc0, iComm,iError)
    call MPI_Allreduce(LocalIRealization==GlobalDomain%iRealization,&
         IsSynchronized,1,MPI_LOGICAL,MPI_LAND,iComm,iError)
    if(IsSynchronized)RETURN
    if(is_proc0(GlobalDomain%CompID_))then
       GlobalDomain%iRealization=LocalDomain%iRealization
       GlobalDomain%nTreeNodes=LocalDomain%nTreeNodes
       call check_octree_allocation(GlobalDomain)
       GlobalDomain%iDecomposition_II(:,1:GlobalDomain%nTreeNodes)=&
            LocalDomain%iDecomposition_II(:,1:LocalDomain%nTreeNodes)
    end if
    call bcast_indexes(GlobalDomain,&
         iProcUnion,iCommUnion)
  end subroutine synchronize_refinement_dd
  !============================================================================
  ! The "methods" below show how to use the information available
  ! from grid descriptors.They involve some elements of the
  ! connectivity list for blocks and searching tools
  !
  !================CONNECTIVITY LIST==============================
  !
  ! Returns the nDim vector, whose iDim's component is true, if
  ! along the iDim's direction there is a domain Boundary to the
  ! Left from the block with the number lGlobalTreeNumber

  function is_left_boundary_d(&
       Domain,lGlobalTreeNumber)
    integer,intent(in)::lGlobalTreeNumber
    type(DomainType),intent(in)::&
         Domain
    logical,dimension(Domain%nDim)::&
         is_left_boundary_d
    !--------------------------------------------------------------------------
    is_left_boundary_d=.not.(Domain%IsPeriodic_D)&
         .and.Domain%CoordBlock_DI(&
         :,lGlobalTreeNumber)<&
         cThird*Domain%DCoordBlock_DI(&
         :,lGlobalTreeNumber)+&
         Domain%CoordMin_D
  end function is_left_boundary_d
  !============================================================================
  ! Returns the nDim vector, whose iDim's component is true, if
  ! along the iDim's direction there is a Tree Boundary to the
  ! Right from the block with the number lGlobalTreeNumber

  function is_right_boundary_d(&
       Domain,lGlobalTreeNumber)
    integer,intent(in)::lGlobalTreeNumber
    type(DomainType),intent(in)::&
         Domain
    logical,dimension(Domain%nDim)::&
         is_right_boundary_d
    !--------------------------------------------------------------------------
    is_right_boundary_d=.not.&
         (Domain%IsPeriodic_D).and.&
         Domain%CoordBlock_DI(:,lGlobalTreeNumber)+&
         (1.0+cThird)*Domain%DCoordBlock_DI(&
         :,lGlobalTreeNumber)>&
         Domain%CoordMax_D
  end function is_right_boundary_d
  !============================================================================
  ! Tree neighbor
  ! returns the number of octree node, to which the center of the
  !"ghost cell" belongs, which is marked with iCells\_D cell index
  ! vector.

  integer function l_neighbor(&
       Domain,lGlobalTreeNumber,iCells_D)
    integer,intent(in)::lGlobalTreeNumber
    type(DomainType),intent(inout)::&
         Domain
    !--------------------------------------------------------------------------
    integer,dimension(Domain%nDim),&
         intent(in)::iCells_D
    real,dimension(Domain%nDim)::Coord_D
    if(any(is_left_boundary_d( Domain,&
         lGlobalTreeNumber).and.iCells_D<1).or.&
         any(is_right_boundary_d(Domain,&
         lGlobalTreeNumber).and.iCells_D>&
         Domain%nCells_D))then
       l_neighbor=None_
    else

       ! Cell center coordinates
       Coord_D = &
            Domain%CoordBlock_DI(:,lGlobalTreeNumber) + &
            Domain%DCoordCell_DI(:,lGlobalTreeNumber)*&
            (real(iCells_D) - 0.5)
       call search_in_dd(Domain, Coord_D, l_neighbor)
    end if
  end function l_neighbor
  !============================================================================
  subroutine glue_margin_dd(&
       Domain,Coord_D)
    type(DomainType),intent(in)::&
         Domain
    !--------------------------------------------------------------------------
    real,dimension(Domain%nDim),&
         intent(inout)::Coord_D
    !-----------

    ! Loop variable:
    integer :: iDir, nDim

    ! Directions to be glued
    integer :: iDirMinusGlue, iDirPlusGlue, iDirCycle

    ! Domain boundaries:
    real,dimension(:),pointer::CoordMin_D,CoordMax_D
    logical :: DoCycle
    !-----------------
    nDim = Domain%nDim
    CoordMin_D => Domain%CoordMin_D
    CoordMax_D => Domain%CoordMax_D
    iDirCycle = Domain%iDirCycle
    iDirMinusGlue = Domain%iDirMinusGlue
    iDirPlusGlue = Domain%iDirPlusGlue
    DoCycle = iDirCycle > 0.and.iDirCycle <= nDim
    do iDir = 1, nDim
       if(iDir == iDirCycle)CYCLE
       if(Coord_D(iDir) < CoordMin_D(iDir).and.iDir==iDirMinusGlue)then
          ! Ruturn the coordinate along the glued direction to the domain
          Coord_D(iDir) = 2*CoordMin_D(iDir) - Coord_D(iDir)
          ! Add, if needed a half cycle over cycled direction
          if(DoCycle)&
               Coord_D(iDirCycle) = CoordMin_D(iDirCycle) + &
               modulo(Coord_D(iDirCycle) + 0.50*&
                CoordMax_D(iDirCycle) - 1.50*CoordMin_D(iDirCycle),&
                CoordMax_D(iDirCycle) -  CoordMin_D(iDirCycle))
       elseif(Coord_D(iDir) > CoordMax_D(iDir).and.iDir==iDirPlusGlue)then
          ! Ruturn the coordinate along the glued direction to the domain
          Coord_D(iDir) = 2*CoordMax_D(iDir) - Coord_D(iDir)
          ! Add, if needed a half cycle over cycled direction
          if(DoCycle)&
               Coord_D(iDirCycle) = CoordMin_D(iDirCycle) + &
               modulo(Coord_D(iDirCycle) + 0.50*&
               CoordMax_D(iDirCycle) - 1.50*CoordMin_D(iDirCycle),&
               CoordMax_D(iDirCycle) -  CoordMin_D(iDirCycle))
       end if
    end do
  end subroutine glue_margin_dd
  !============================================================================
  !\begin{verbatim}
  !=====================SEARCH====================================
  ! The searching tools start from here which allow to find the
  ! location in the domain decomposition, using the generalized
  !  coordinates of the point.
  !
  ! Returns lGlobalTreeNumber, which is the  global number of the
  ! block which includes the point Coord_D.
  ! The starting values for Coord_D should be the generalized
  ! coordinates for the point, finally in this array there are the
  ! coordinates with respect to the left corner of the found block.
  ! The values of Coord_D are allowed to be outside of the domain,
  ! the nearest block is found in this case.
  !\end{verbatim}

  subroutine search_in_dd(&
       Domain,Coord_D,lGlobalTreeNumber)
    type(DomainType),intent(inout)::&
         Domain
    !--------------------------------------------------------------------------
    real,dimension(Domain%nDim),&
         intent(inout)::Coord_D
    integer,intent(out)::lGlobalTreeNumber
    real,dimension(Domain%nDim)::&
         CoordTrunc_D,Discr_D
    integer,dimension(Domain%nDim)::&
         iRootMinusOne_D,iRootMinusOneStart_D
    integer,dimension(Domain%nDimTree)::iShift_D
    integer::iChild,iDim,lFound

    ! Start from the result of the previous search
    lFound=Domain%lSearch

    ! Put the point inside the domain
    do iDim=1,Domain%nDim
       if(.not.Domain%IsPeriodic_D(iDim))CYCLE
       Coord_D(iDim)=Domain%CoordMin_D(iDim)+&
            modulo(Coord_D(iDim)-Domain%CoordMin_D(iDim),&
         Domain%CoordMax_D(iDim)-&
         Domain%CoordMin_D(iDim))
    end do
    CoordTrunc_D=Domain%CoordMin_D+&
         max(0.0,min(Coord_D-Domain%CoordMin_D,&
         cAlmostOne*(Domain%CoordMax_D-&
         Domain%CoordMin_D)))

    Discr_D=(CoordTrunc_D-Domain%CoordBlock_DI(&
         :,lFound))/Domain%DCoordBlock_DI(:,lFound)

    ! Recursive search starts
    do

       ! correct Discr_D to avoid having points exactly
       ! on the boundary of a block; potentiall introduces error,
       ! but it would matter -log_2(cTol)~23 resolution level up the tree
       Discr_D = real(floor(Discr_D)) + &
            min( cAlmostOne, max( cTol2, Discr_D-real(floor(Discr_D)) ) )

       if(any(Discr_D<0.0).or.any(Discr_D>=1.0))then
          iChild=Domain%iDecomposition_II(&
               MyNumberAsAChild_,lFound)
          if(iChild==0)then  ! This is a root

             iRootMinusOneStart_D=&
                  nint((Domain%CoordBlock_DI(&
                  :,lFound)-Domain%CoordMin_D)&
                  /Domain%DCoordBlock_DI(:,lFound))
             ! Calculate iRoot-1,jRoot-1....
             iRootMinusOne_D=floor(Discr_D)+&
                  iRootMinusOneStart_D
             Discr_D=Discr_D-floor(Discr_D)
             ! Calculate the root number, using the formula
             ! lRoot-1=iRoot-1+&
             !(jRoot-1)*RootMapDim(1)+&
             !(kRoot-1)*RootMapDim(1)*RootMapDim(2)
             lFound=1+iRootMinusOne_D(1)
             do iDim=2,Domain%nDim
                lFound=lFound+iRootMinusOne_D(iDim)*&
                     product(&
                     Domain%iRootMapDim_D(1:iDim-1))
             end do
             if(Domain%IsTreeDecomposition)&
                  lFound=Domain%iRoot_I(lFound)
          else     ! End of computations for root
             !---------------------------------------------------------------!
             ! Descend the octree                                             !
             lFound=Domain%iDecomposition_II(&
                  Parent_,lFound)
             Discr_D(1:Domain%nDimTree)=&
                  (Discr_D(1:Domain%nDimTree)+&
                  real(Domain%iShift_DI(&
                  :,iChild)))* 0.50
          end if
       elseif&
            (is_used_block_dd(Domain,lFound))then
          EXIT ! Octree node is found
       else
          !---------------------------------------------------------------!
          ! Ascend the octree: calculate the shift                         !
          Discr_D(1:Domain%nDimTree)=&
               Discr_D(1:Domain%nDimTree)&
               *2.0
          iShift_D=int(Discr_D(&
               1:Domain%nDimTree))
          Discr_D(1:Domain%nDimTree)=&
               Discr_D(1:Domain%nDimTree)-&
               real(iShift_D)
          !---------------------------------------------------------------!
          ! Choose the child to ascend to                                  !
          iChild=1
          do while (any(iShift_D/=&
               Domain%iShift_DI(:,iChild)))
             iChild=iChild+1
          end do
          lFound=Domain%iDecomposition_II(&
               iChild,lFound)
       end if
    end do
    !---------------------------------------------------------------!
    ! End of recursive search                                        !
    Coord_D=Coord_D-Domain%CoordBlock_DI(:,lFound)
    lGlobalTreeNumber=lFound
    Domain%lSearch=lFound
  end subroutine search_in_dd
  !============================================================================
  ! In searching the cell the original values of the generalized
  ! coordinates must be defined with respect to the left corner and
  ! a global tree node number should be known
  ! search\_cell ALWAYS returns the position of the Coord point with
  ! respect to the found cell left corner.
  !
  ! The best order to find cell number for a position of point
  ! If the block number is needed
  !\begin{verbatim}
  ! call search_in(Domain,Coord_D,lFound),
  ! iBlock=blk_decomposition(Domain,lFound)
  ! call search_cell(Domain,Coord_D,lFound,Cells_D)
  !\end{verbatim}
  ! the second step can be missed if the block number is not needed
  ! The input Coord\_D can be outside the domain,
  ! covered by grid, in this case the values of the cell
  ! can be less than unity or greater that nCells\_D. This use is
  ! not restricted, because the cell to find may be the ghostcell.
  !---------------------------------------------------------------

  subroutine search_cell_dd(&
       Domain,lGlobalTreeNumber,Coord_D,iCells_D)
    type(DomainType),intent(in)::&
         Domain
    !--------------------------------------------------------------------------
    real,dimension(Domain%nDim),&
         intent(inout)::Coord_D
    integer,intent(in)::lGlobalTreeNumber
    integer,dimension(Domain%nDim),&
         intent(out)::iCells_D
    real,dimension(Domain%nDim)::DCoordCells_D
    DCoordCells_D=Domain%DCoordCell_DI(&
         :,lGlobalTreeNumber)
    iCells_D=floor(Coord_D/DCoordCells_D)
    Coord_D=Coord_D-DCoordCells_D*iCells_D
    iCells_D=iCells_D+1
  end subroutine search_cell_dd
  !============================================================================
  subroutine pe_and_blk_dd(&
       Domain,lGlobalTreeNumber,iPEOut,iBlockOut)
    type(DomainType),intent(in)::&
         Domain
    integer,intent(in)::lGlobalTreeNumber
    ! OUTPUT ARGUMENTS
    integer,intent(out)::iPEOut,iBlockOut
    !--------------------------------------------------------------------------
    iPEOut=Domain%iDecomposition_II(&
         PE_,lGlobalTreeNumber)
    iBlockOut=Domain%iDecomposition_II(&
         BLK_,lGlobalTreeNumber)
  end subroutine pe_and_blk_dd
  !============================================================================

  integer function n_block_dd(&
       Domain,iPE)
    type(DomainType),intent(in)::&
         Domain
    integer,intent(in)::iPE
    !--------------------------------------------------------------------------
    n_block_dd=count(&
         Domain%iDecomposition_II(&
         FirstChild_,1:Domain%nTreeNodes)&
         ==None_.and.&
         Domain%iDecomposition_II(&
         PE_,1:Domain%nTreeNodes)==iPE)
  end function n_block_dd
  !============================================================================
  integer function iglobal_bp_dd(Domain,iBlock,iPE)
    type(DomainType),intent(in)::Domain
    integer,intent(in)::iBlock,iPE
    !--------------------------------------------------------------------------
    iglobal_bp_dd=Domain%iGlobal_BP(iBlock,iPE)
  end function iglobal_bp_dd
  !============================================================================
  !---------------------------------------------------------------!
  integer function iglobal_node_dd(Domain,iBlockAll)
    type(DomainType),intent(in)::Domain
    integer,intent(in)::iBlockAll
    !--------------------------------------------------------------------------
    iglobal_node_dd=Domain%iGlobal_A(iBlockAll)
  end function iglobal_node_dd
  !============================================================================
  !---------------------------------------------------------------!
  integer function iglobal_block_dd(&
       Domain,iGlobalTreeNode)
    type(DomainType),intent(in)::Domain
    integer,intent(in)::iGlobalTreeNode
    !--------------------------------------------------------------------------
    iglobal_block_dd=Domain%iDecomposition_II(&
         GlobalBlock_,iGlobalTreeNode)
  end function iglobal_block_dd
  !============================================================================
  !---------------------------------------------------------------!
  !\begin{verbatim}
  !             Access to the elements of the structures          !
  !\end{verbatim}

  integer function compid_dd(Domain)
    ! INPUT ARGUMENTS
    type(DomainType),intent(in)::&
    !--------------------------------------------------------------------------
    Domain
    compid_dd=Domain%CompID_
  end function compid_dd
  !============================================================================
  integer function ndim_dd(Domain)
    type(DomainType),intent(in)::&
         Domain
    !--------------------------------------------------------------------------
    nDim_dd=Domain%nDim
  end function ndim_dd
  !============================================================================

  function coord_block_dd(&
       Domain,lGlobalTreeNumber)
    type(DomainType),intent(in)::&
         Domain
    integer,intent(in)::lGlobalTreeNumber
    real,dimension(Domain%nDim)::coord_block_dd
    !--------------------------------------------------------------------------
    coord_block_dd=Domain%CoordBlock_DI(&
         :,lGlobalTreeNumber)
  end function coord_block_dd
  !============================================================================

  function d_coord_block_dd(&
       Domain,lGlobalTreeNumber)
    type(DomainType),intent(in)::&
         Domain
    integer,intent(in)::lGlobalTreeNumber
    real,dimension(Domain%nDim)::d_coord_block_dd
    !--------------------------------------------------------------------------
    d_coord_block_dd=Domain%DCoordBlock_DI(&
         :,lGlobalTreeNumber)
  end function d_coord_block_dd
  !============================================================================

  function d_coord_cell_dd(&
       Domain,lGlobalTreeNumber)
    type(DomainType),intent(in)::&
         Domain
    integer,intent(in)::lGlobalTreeNumber
    real,dimension(Domain%nDim)::d_coord_cell_dd
    !--------------------------------------------------------------------------
    d_coord_cell_dd=Domain%DCoordCell_DI(&
         :,lGlobalTreeNumber)
  end function d_coord_cell_dd
  !============================================================================

  logical function is_used_block_dd(&
       Domain,lGlobalTreeNumber)
    type(DomainType),intent(in)::&
         Domain
    integer,intent(in)::lGlobalTreeNumber
    !--------------------------------------------------------------------------
    is_used_block_dd=&
         Domain%iDecomposition_II(&
         FirstChild_,lGlobalTreeNumber)==None_
  end function is_used_block_dd
  !============================================================================

  integer function irealization_dd(Domain)
    type(DomainType),intent(in)::&
         Domain
    !--------------------------------------------------------------------------
    irealization_dd=Domain%iRealization
  end function irealization_dd
  !============================================================================
  subroutine associate_dd_pointer_dd(&
       Domain,DomainPointer)
    Type(DomainType),target,intent(in)::&
         Domain
    type(DomainPointerType),intent(out)::DomainPointer
    !--------------------------------------------------------------------------
    nullify(DomainPointer%Ptr)
    DomainPointer%Ptr=>Domain
  end subroutine associate_dd_pointer_dd
  
end module CON_domain_decomposition
!==============================================================================

