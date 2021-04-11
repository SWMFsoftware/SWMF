!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module CON_grid_storage
  use CON_domain_decomposition
  use ModInitGridStorage, ONLY:init_grid_storage, MaxGrid
  implicit none

  ! The resulted domain decompositions should be properly
  ! registered and should obtain the unique GridID
  type(DomainPointerType), private :: Domain_I(MaxGrid)
  logical   :: DoneDomainInit_C(MaxGrid)=.false.
  ! Introduced to bypass the HALEM compiler restrictions
  integer :: nDim_C(MaxGrid) = 0
  interface init_decomposition
     module procedure init_decomposition_id
     module procedure init_decomposition_dd
  end interface

  interface get_root_decomposition
     module procedure get_root_decomposition_id
     module procedure get_root_decomposition_dd
  end interface

  interface bcast_decomposition
     module procedure bcast_decomposition_id
     module procedure bcast_decomposition_dd
  end interface

  interface synchronize_refinement
     module procedure synchronize_refinement_id
     module procedure synchronize_refinement_dd
  end interface

  interface glue_margin
     module procedure glue_margin_id
     module procedure glue_margin_dd
  end interface glue_margin

  interface search_in
     module procedure search_in_id
     module procedure search_in_dd
  end interface

  interface search_cell
     module procedure search_cell_id
     module procedure search_cell_dd
  end interface

  interface i_global_node_a
     module procedure iglobal_node_dd
     module procedure iglobal_node_id
  end interface

  interface i_global_block
     module procedure iglobal_block_dd
     module procedure iglobal_block_id
  end interface

  interface associate_dd_pointer
     module procedure associate_dd_pointer_dd
     module procedure associate_dd_pointer_id
  end interface associate_dd_pointer
contains
  !============================================================================
  logical function done_dd_init(GridID_)
    integer,intent(in)::GridID_
    ! Returns .true. if the domain decomposition is initialaized
    !--------------------------------------------------------------------------
    done_dd_init=.false.
    if(GridID_>0 .and. GridID_<=MaxGrid)&
         done_dd_init = DoneDomainInit_C(GridID_)
  end function done_dd_init
  !============================================================================
  !---------------------------------------------------------------!
  ! Initialization for the decomposition.
  ! Note that if the decomposition is a tree decomposition and
  ! another order of children is accepted in the component,
  ! different from that assumed to be standard one
  ! the array iShift\_DI should be reassigned in set\_root.
  subroutine init_decomposition_id(&
       GridID_,&
       CompID_,&
       nDim,   &
       IsTreeDD)

    integer,intent(in)::GridID_
    integer,intent(in)::CompID_,nDim
    logical,intent(in),optional::IsTreeDD
    !--------------------------------------------------------------------------
    if(GridID_>MaxGrid.or.GridID_<=0)call CON_stop(&
         'GridID_ is not allowed to be equal to ')
    if(done_dd_init(GridID_))call CON_stop(&
         'An attempt to reinitialize the gridID_')
    call init_grid_storage(Domain_I, GridID_)
    call init_decomposition_dd(Domain_I(GridID_)%Ptr, CompID_, nDim, IsTreeDD)
    DoneDomainInit_C(GridID_) = .true.
    nDim_C(GridID_)           = nDim
  end subroutine init_decomposition_id
  !============================================================================
  ! To get a decomposition domain, even the tree one, the root     !
  ! decomposition should be first constructed                      !
  !====================ATTENTION==================================!
  ! PE here are the ranks in the LOCAL communicator for the        !
  ! component
  subroutine get_root_decomposition_id(&
       GridID_,&!     ! ID for Decomposition to be constructed
       iRootMapDim_D,&! As in DomainType
       CoordMin_D,&   ! As in DomainType
       CoordMax_D,&   ! As in DomainType
       nCell_D,&     ! As in DomainType
       PE_I,&         ! PE layout
       iBlock_I,&     ! Local Block Number layout
       IsPeriodic_D,& ! As in DomainType
       iShift_DI, &   ! As in DomainType
       DoGlueMargins,&! As in DomainType
       iDirMinusGlue,&! As in DomainType
       iDirPlusGlue,& ! As in DomainType
       iDirCycle)     ! As in DomainType)

    integer, intent(in) :: GridID_
    integer, intent(in) :: iRootMapDim_D(:)
    real,    intent(in) :: CoordMin_D(:)
    real,    intent(in) :: CoordMax_D(:)
    integer, intent(in) :: nCell_D(:)
    integer, optional, intent(in) :: PE_I(:),iBlock_I(:)
    logical, optional, intent(in) :: IsPeriodic_D(:)
    integer, optional, intent(in) :: iShift_DI(:,:)
    logical, optional, intent(in) :: DoGlueMargins
    integer, optional, intent(in) :: iDirMinusGlue
    integer, optional, intent(in) :: iDirPlusGlue
    integer, optional, intent(in) :: iDirCycle
    integer::lBlock,MaxBlock
    
    ! Check the dimension of PE_I and iBlock_I                       
    if(present(PE_I))then
       if(ubound(PE_I,1)/=product(iRootMapDim_D))&
            call CON_stop('The dimension of PE_I should be equal to '//&
            'product(iRootMapDim_D)')
    end if

    if(present(iBlock_I))then
       if(ubound(iBlock_I,1)/=product(iRootMapDim_D))&
            call CON_stop('The dimension of iBlock_I should be equal to '//&
            'product(iRootMapDim_D)')
    end if

    ! Assign the grid parameters                                     !

    Domain_I(GridID_)%Ptr%iRootMapDim_D = iRootMapDim_D
    Domain_I(GridID_)%Ptr%CoordMin_D    = CoordMin_D
    Domain_I(GridID_)%Ptr%CoordMax_D    = CoordMax_D
    Domain_I(GridID_)%Ptr%nCell_D      = nCell_D
    if(present(DoGlueMargins)) &
         Domain_I(GridID_)%Ptr%DoGlueMargins = DoGlueMargins
    if(present(iDirMinusGlue))&
         Domain_I(GridID_)%Ptr%iDirMinusGlue = iDirMinusGlue
    if(present(iDirPlusGlue)) &
         Domain_I(GridID_)%Ptr%iDirPlusGlue  = iDirPlusGlue
    if(present(iDirCycle))    &
         Domain_I(GridID_)%Ptr%iDirCycle     = iDirCycle
    if(present(IsPeriodic_D))&
         Domain_I(GridID_)%Ptr%IsPeriodic_D=IsPeriodic_D

    Domain_I(GridID_)%Ptr%nTreeNodes=product(iRootMapDim_D)
    call check_octree_allocation(Domain_I(GridID_)%Ptr)

    Domain_I(GridID_)%Ptr%iDD_II&
         (MyNumberAsAChild_,1:Domain_I(GridID_)%Ptr%nTreeNodes)=0
    Domain_I(GridID_)%Ptr%iDD_II&
         (FirstChild_,1:Domain_I(GridID_)%Ptr%nTreeNodes)=None_
    do lBlock=1,Domain_I(GridID_)%Ptr%nTreeNodes
       Domain_I(GridID_)%Ptr%iDD_II&
            (Parent_,1:Domain_I(GridID_)%Ptr%nTreeNodes)=lBlock
    end do
    if(Domain_I(GridID_)%Ptr%IsTreeDD)then
       call check_iroot_allocation(Domain_I(GridID_)%Ptr)
       if(present(iShift_DI))&
            Domain_I(GridID_)%Ptr%iShift_DI=iShift_DI
    end if

    ! If neither PE_I nor iBlock_I is present, the root              !
    ! decomposition blocks are balanced for an optimal load          !
    if((.not.present(PE_I)).and.(.not.present(iBlock_I)))then
       MaxBlock=(Domain_I(GridID_)%Ptr%nTreeNodes-1)/&
            n_proc(Domain_I(GridID_)%Ptr%CompID_)+1
       do lBlock=1,Domain_I(GridID_)%Ptr%nTreeNodes
          Domain_I(GridID_)%Ptr%iDD_II&
               (BLK_,lBlock)=mod(lBlock,MaxBlock)+1
          Domain_I(GridID_)%Ptr%iDD_II&
               (PE_,lBlock)=(lBlock-1)/MaxBlock
       end do
       call set_iglobal_and_bp_dd(Domain_I(GridID_)%Ptr)
       RETURN
    end if

    ! If the PE\_I is not given, it is assumed that all the blocks
    ! are at the root PE

    if(present(PE_I))then
       Domain_I(GridID_)%Ptr%iDD_II&
            (PE_,1:Domain_I(GridID_)%Ptr%nTreeNodes)=PE_I
    else
       Domain_I(GridID_)%Ptr%iDD_II&
            (PE_,1:Domain_I(GridID_)%Ptr%nTreeNodes)=0
    end if

    if(present(iBlock_I))then
       Domain_I(GridID_)%Ptr%iDD_II&
            (BLK_,1:Domain_I(GridID_)%Ptr%nTreeNodes)=iBlock_I
    else
       do lBlock=1, Domain_I(GridID_)%Ptr%nTreeNodes
          Domain_I(GridID_)%Ptr%iDD_II&
               (BLK_,lBlock)=lBlock
       end do
    end if
    call set_iglobal_and_bp_dd(Domain_I(GridID_)%Ptr)
  end subroutine get_root_decomposition_id
  !============================================================================
  subroutine bcast_decomposition_id(&
       GridID_)

    integer,intent(in)::GridID_
    !--------------------------------------------------------------------------
    call bcast_decomposition_dd(Domain_I(GridID_)%Ptr)
  end subroutine bcast_decomposition_id
  !============================================================================
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
  subroutine synchronize_refinement_id(&
       GridID_,LocalDomain,iProcUnion,iCommUnion)

    integer,intent(in)::GridID_
    !--------------------------------------------------------------------------
    type(DomainType),&   ! Local Grid with which the
         intent(in)::LocalDomain          ! global one is synchronized
    integer,intent(in),optional::iProcUnion,iCommUnion
    if(present(iProcUnion).and.present(iCommUnion))then
       call synchronize_refinement_dd(&
            Domain_I(GridID_)%Ptr,LocalDomain,iProcUnion,iCommUnion)
    else
       call synchronize_refinement_dd(&
            Domain_I(GridID_)%Ptr,LocalDomain)
    end if
  end subroutine synchronize_refinement_id
  !============================================================================
  subroutine glue_margin_id(GridID_,Coord_D)
    integer,intent(in)::GridID_
    !--------------------------------------------------------------------------
    real,dimension(nDim_C(GridID_)),&
         intent(inout)::Coord_D
    !--------------------------------------------------------------------------
    call glue_margin_dd(Domain_I(GridID_)%Ptr,Coord_D)
  end subroutine glue_margin_id
  !============================================================================
  subroutine search_in_id(GridID_,Coord_D,lGlobalTreeNumber)
    integer,intent(in)::GridID_
    !--------------------------------------------------------------------------
    real,dimension(nDim_C(GridID_)),&
         intent(inout)::Coord_D
    integer,intent(out)::lGlobalTreeNumber
    !--------------------------------------------------------------------------
    call search_in_dd(Domain_I(GridID_)%Ptr,Coord_D,lGlobalTreeNumber)
  end subroutine search_in_id
  !============================================================================
  subroutine search_cell_id(&
       GridID_,lGlobalTreeNumber,Coord_D,iCells_D)
    integer,intent(in)::GridID_,lGlobalTreeNumber
    !--------------------------------------------------------------------------
    real,dimension(nDim_C(GridID_)),&
         intent(inout)::Coord_D
    integer,dimension(nDim_C(GridID_)),&
         intent(out)::iCells_D
    !--------------------------------------------------------------------------
    call search_cell_dd(Domain_I(GridID_)%Ptr,&
         lGlobalTreeNumber,Coord_D,iCells_D)
  end subroutine search_cell_id
  !============================================================================
  integer function iglobal_node_id(GridID_,iBlockAll)
    integer,intent(in)::GridID_
    integer,intent(in)::iBlockAll
    !--------------------------------------------------------------------------
    iglobal_node_id=Domain_I(GridID_)%Ptr%iGlobal_A(iBlockAll)
  end function iglobal_node_id
  !============================================================================
  integer function iglobal_block_id(&
       GridID_,iGlobalTreeNode)
    integer,intent(in)::GridID_
    integer,intent(in)::iGlobalTreeNode
    !--------------------------------------------------------------------------
    iglobal_block_id=Domain_I(GridID_)%Ptr%iDD_II(&
         GlobalBlock_,iGlobalTreeNode)
  end function iglobal_block_id
  !============================================================================
  integer function i_realization(GridID_)
    integer,intent(in)::GridID_
    !--------------------------------------------------------------------------
    i_realization = Domain_I(GridID_)%Ptr%iRealization
  end function i_realization
  !============================================================================
  subroutine associate_dd_pointer_id(GridID_,DomainPointer)
    integer,intent(in)::GridID_
    type(DomainPointerType),intent(out)::DomainPointer
    !--------------------------------------------------------------------------
    nullify(DomainPointer%Ptr)
    DomainPointer%Ptr=>Domain_I(GridID_)%Ptr
  end subroutine associate_dd_pointer_id
  !============================================================================
  integer function ndim_id(GridID_)
    integer,intent(in)::GridID_
    !--------------------------------------------------------------------------
    ndim_id = Domain_I(GridID_)%Ptr%nDim
  end function ndim_id
  !============================================================================
  function ncell_id(GridID_)
    integer,intent(in)::GridID_
    integer,dimension(nDim_C(GridID_))::ncell_id
    !--------------------------------------------------------------------------
    ncell_id=Domain_I(GridID_)%Ptr%nCell_D
  end function ncell_id
  !============================================================================
  integer function compid(GridID_)
    integer,intent(in)::GridID_
    !--------------------------------------------------------------------------
    compid = Domain_I(GridID_)%Ptr%CompID_
  end function compid
end module CON_grid_storage
!==============================================================================

