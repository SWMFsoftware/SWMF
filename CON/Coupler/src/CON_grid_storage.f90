!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module CON_grid_storage
  use CON_domain_decomposition
  use ModInitGridStorage, ONLY:init_grid_storage, MaxGrid
  implicit none

  ! The resulted domain decompositions should be properly
  ! registered and should obtain the unique GridID\_
  type(DomainPointerType),dimension(MaxGrid),&
       private::&
       Domain_I
  logical :: DoneDomainInit_C(MaxGrid)=.false.

  ! Introduced to bypass the HALEM compiler restrictions
  integer,dimension(MaxGrid)::nDim_C=0
  ! INTERFACES
  ! Each of the module procedures differ with the FIRST dummy
  ! parameter, which is either GridID\_, or DomainType
  ! structure. The procedures differ with suffix id or dd from the
  ! generic name

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

  interface pe_and_blk
     module procedure pe_and_blk_id
     module procedure pe_and_blk_dd
  end interface

  interface n_block
     module procedure n_block_id
     module procedure n_block_dd
  end interface

  interface search_cell
     module procedure search_cell_id
     module procedure search_cell_dd
  end interface

  interface i_global_node_bp
     module procedure iglobal_bp_dd
     module procedure iglobal_bp_id
  end interface

  interface i_global_node_a
     module procedure iglobal_node_dd
     module procedure iglobal_node_id
  end interface

  interface i_global_block
     module procedure iglobal_block_dd
     module procedure iglobal_block_id
  end interface

  interface compid
     module procedure compid_id
     module procedure compid_dd
  end interface

  interface ndim_grid
     module procedure ndim_id
  end interface

  interface coord_block_d
     module procedure coord_block_id
     module procedure coord_block_dd
  end interface

  interface d_coord_block_d
     module procedure d_coord_block_id
     module procedure d_coord_block_dd
  end interface

  interface d_coord_cell_d
     module procedure d_coord_cell_id
     module procedure d_coord_cell_dd
  end interface

  interface i_realization
     module procedure irealization_id
  end interface

  interface associate_dd_pointer
     module procedure associate_dd_pointer_id
     module procedure associate_dd_pointer_dd
  end interface
contains
  !============================================================================
  logical function done_dd_init(GridID_)
    integer,intent(in)::GridID_
    ! Returns .true. if the domain decomposition is initialaized
    !--------------------------------------------------------------------------
    done_dd_init=.false.
    if(GridID_>0.and.GridID_<=MaxGrid)&
         done_dd_init=DoneDomainInit_C(GridID_)
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
       IsTreeDecomposition,&
       nDimTree)

    integer,intent(in)::GridID_
    integer,intent(in)::CompID_,nDim
    logical,intent(in),optional::IsTreeDecomposition
    integer,intent(in),optional::nDimTree
    !--------------------------------------------------------------------------
    if(GridID_>MaxGrid.or.GridID_<=0)call CON_stop(&
         'GridID_ is not allowed to be equal to ')
    if(done_dd_init(GridID_))call CON_stop(&
         'An attempt to reinitialize the gridID_')
    call init_grid_storage(Domain_I,GridID_)
    if(present(nDimTree))then
       call init_decomposition_dd(&
            Domain_I(GridID_)%Ptr,&
            CompID_,&
            nDim,   &
            IsTreeDecomposition,&
            nDimTree)
    elseif(present(IsTreeDecomposition))then
       call init_decomposition_dd(&
            Domain_I(GridID_)%Ptr,&
            CompID_,&
            nDim,   &
            IsTreeDecomposition)
    else
       call init_decomposition_dd(&
            Domain_I(GridID_)%Ptr,&
            CompID_,&
            nDim)
    end if
    DoneDomainInit_C(GridID_)=.true.
    nDim_C(GridID_)=nDim
  end subroutine init_decomposition_id
  !============================================================================
  !---------------------------------------------------------------!
  !==========================WITH INTERFACE=======================!
  !\begin{verbatim}
  ! To get a decomposition domain, even the tree one, the root     !
  ! decomposition should be first constructed                      !
  !====================ATTENTION==================================!
  ! PE here are the ranks in the LOCAL communicator for the        !
  ! component
  !\end{verbatim}
  subroutine get_root_decomposition_id(&
       GridID_,&!     ! ID for Decomposition to be constructed
       iRootMapDim_D,&! As in DomainType
       CoordMin_D,&     ! As in DomainType
       CoordMax_D,&     ! As in DomainType
       nCells_D,&     ! As in DomainType
       PE_I,&         ! PE layout
       iBlock_I,&     ! Local Block Number layout
       IsPeriodic_D,& ! As in DomainType
       iShift_DI, &   ! As in DomainType
       DoGlueMargins,&! As in DomainType
       iDirMinusGlue,&! As in DomainType
       iDirPlusGlue,& ! As in DomainType
       iDirCycle)     ! As in DomainType)

    integer,intent(in)::GridID_
    !--------------------------------------------------------------------------
    integer,dimension(:),&
         intent(in)::iRootMapDim_D
    real,dimension(:),&
         intent(in)::CoordMin_D
    real,dimension(:),&
         intent(in)::CoordMax_D
    integer,dimension(:),&
         intent(in)::nCells_D
    integer,dimension(:),&
         intent(in),optional:: PE_I,iBlock_I
    logical,dimension(:),&
         intent(in),optional::IsPeriodic_D
    integer,optional,&
         dimension(:,&
         :),&
         intent(in)::iShift_DI
    logical, optional, intent(in) ::  DoGlueMargins
    integer, optional, intent(in) ::  iDirMinusGlue
    integer, optional, intent(in) ::  iDirPlusGlue
    integer, optional, intent(in) ::  iDirCycle
    integer::lBlock,MaxBlock
    !---------------------------------------------------------------!
    ! Check the dimension of PE_I and iBlock_I                       !
    if(present(PE_I))then
       if(ubound(PE_I,1)/=product(iRootMapDim_D))&
            call CON_stop(&
            'The dimension of PE_I should be equal to '//&
            'product(iRootMapDim_D)')
    end if

    if(present(iBlock_I))then
       if(ubound(iBlock_I,1)/=product(iRootMapDim_D))&
            call CON_stop(&
            'The dimension of iBlock_I should be equal to '//&
            'product(iRootMapDim_D)')
    end if

    ! Assign the grid parameters                                     !

    Domain_I(GridID_)%Ptr%iRootMapDim_D=iRootMapDim_D

    Domain_I(GridID_)%Ptr%CoordMin_D=CoordMin_D
    Domain_I(GridID_)%Ptr%CoordMax_D=CoordMax_D
    Domain_I(GridID_)%Ptr%nCells_D=nCells_D
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

    Domain_I(GridID_)%Ptr%iDecomposition_II&
         (MyNumberAsAChild_,1:Domain_I(GridID_)%Ptr%nTreeNodes)=0
    Domain_I(GridID_)%Ptr%iDecomposition_II&
         (FirstChild_,1:Domain_I(GridID_)%Ptr%nTreeNodes)=None_
    do lBlock=1,Domain_I(GridID_)%Ptr%nTreeNodes
       Domain_I(GridID_)%Ptr%iDecomposition_II&
            (Parent_,1:Domain_I(GridID_)%Ptr%nTreeNodes)=lBlock
    end do
    if(Domain_I(GridID_)%Ptr%IsTreeDecomposition)then
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
          Domain_I(GridID_)%Ptr%iDecomposition_II&
               (BLK_,lBlock)=mod(lBlock,MaxBlock)+1
          Domain_I(GridID_)%Ptr%iDecomposition_II&
               (PE_,lBlock)=(lBlock-1)/MaxBlock
       end do
       call set_iglobal_and_bp_dd(Domain_I(GridID_)%Ptr)
       RETURN
    end if

    ! If the PE\_I is not given, it is assumed that all the blocks
    ! are at the root PE

    if(present(PE_I))then
       Domain_I(GridID_)%Ptr%iDecomposition_II&
            (PE_,1:Domain_I(GridID_)%Ptr%nTreeNodes)=PE_I
    else
       Domain_I(GridID_)%Ptr%iDecomposition_II&
            (PE_,1:Domain_I(GridID_)%Ptr%nTreeNodes)=0
    end if

    ! If the iBlock_I is not given, they are assumed to be enumerated!
    ! like in the following loop:                                    !
    !      lBlock=0                                                 !
    !      do kRoot=1,iRootMapDim_D(3)                              !
    !          do jRoot=1,iRootMapDim_D(2)                          !
    !             do iRoot=1,iRootMapDim_D(1)                       !
    !               lBlock=lBlock+1                                 !
    !               BLK(iRoot,jRoot,kRoot)=lBlock                   !
    !             end do                                            !
    !          end do                                               !
    !      end do                                                   !

    if(present(iBlock_I))then
       Domain_I(GridID_)%Ptr%iDecomposition_II&
            (BLK_,1:Domain_I(GridID_)%Ptr%nTreeNodes)=iBlock_I
    else
       do lBlock=1, Domain_I(GridID_)%Ptr%nTreeNodes
          Domain_I(GridID_)%Ptr%iDecomposition_II&
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
    real,dimension(nDim_C(GridID_)),&
         intent(inout)::Coord_D
    !--------------------------------------------------------------------------
    call glue_margin_dd(Domain_I(GridID_)%Ptr,Coord_D)
  end subroutine glue_margin_id
  !============================================================================
  subroutine search_in_id(GridID_,Coord_D,lGlobalTreeNumber)
    integer,intent(in)::GridID_
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
    real,dimension(nDim_C(GridID_)),&
         intent(inout)::Coord_D
    integer,dimension(nDim_C(GridID_)),&
         intent(out)::iCells_D
    !--------------------------------------------------------------------------
    call search_cell_dd(Domain_I(GridID_)%Ptr,&
         lGlobalTreeNumber,Coord_D,iCells_D)
  end subroutine search_cell_id
  !============================================================================
  subroutine pe_and_blk_id(&
       GridID_,lGlobalTreeNumber,iPEOut,iBlockOut)
    integer,intent(in)::GridID_,lGlobalTreeNumber
    integer,intent(out)::iPEOut,iBlockOut
    !--------------------------------------------------------------------------
    call pe_and_blk_dd(Domain_I(GridID_)%Ptr,&
         lGlobalTreeNumber,iPEOut,iBlockOut)
  end subroutine pe_and_blk_id
  !============================================================================
  integer function n_block_id(GridID_,iPE)
    integer,intent(in)::GridID_,iPE
    !--------------------------------------------------------------------------
    n_block_id=n_block_dd(Domain_I(GridID_)%Ptr,iPE)
  end function n_block_id
  !============================================================================
  integer function iglobal_bp_id(GridID_,iBlock,iPE)
    integer,intent(in)::GridID_
    integer,intent(in)::iBlock,iPE
    !--------------------------------------------------------------------------
    iglobal_bp_id=Domain_I(GridID_)%Ptr%iGlobal_BP(iBlock,iPE)
  end function iglobal_bp_id
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
    iglobal_block_id=Domain_I(GridID_)%Ptr%iDecomposition_II(&
         GlobalBlock_,iGlobalTreeNode)
  end function iglobal_block_id
  !============================================================================
  integer function compid_id(GridID_)
    integer,intent(in)::GridID_
    !--------------------------------------------------------------------------
    compid_id=compid_dd(Domain_I(GridID_)%Ptr)
  end function compid_id
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
    ncell_id=Domain_I(GridID_)%Ptr%nCells_D
  end function ncell_id
  !============================================================================
  function coord_min_d(GridID_)
    integer,intent(in)::GridID_
    real,dimension(nDim_C(GridID_)):: coord_min_d
    !--------------------------------------------------------------------------
    coord_min_d = Domain_I(GridID_)%Ptr%CoordMin_D
  end function coord_min_d
  !============================================================================
  function coord_max_d(GridID_)
    integer,intent(in)::GridID_
    real,dimension(nDim_C(GridID_)):: coord_max_d
    !--------------------------------------------------------------------------
    coord_max_d = Domain_I(GridID_)%Ptr%CoordMax_D
  end function coord_max_d
  !============================================================================
  function coord_block_id(GridID_,lGlobalTreeNumber)
    integer,intent(in)::GridID_,lGlobalTreeNumber
    real,dimension(nDim_C(GridID_)):: coord_block_id
    !--------------------------------------------------------------------------
    coord_block_id=coord_block_dd(&
         Domain_I(GridID_)%Ptr,lGlobalTreeNumber)
  end function coord_block_id
  !============================================================================
  function d_coord_block_id(GridID_,lGlobalTreeNumber)
    integer,intent(in)::GridID_,lGlobalTreeNumber
    real,dimension(nDim_C(GridID_)):: d_coord_block_id
    !--------------------------------------------------------------------------
    d_coord_block_id=d_coord_block_dd(&
         Domain_I(GridID_)%Ptr,lGlobalTreeNumber)
  end function d_coord_block_id
  !============================================================================
  function d_coord_cell_id(GridID_,lGlobalTreeNumber)
    integer,intent(in)::GridID_,lGlobalTreeNumber
    real,dimension(nDim_C(GridID_)):: d_coord_cell_id
    !--------------------------------------------------------------------------
    d_coord_cell_id=d_coord_cell_dd(&
         Domain_I(GridID_)%Ptr,lGlobalTreeNumber)
  end function d_coord_cell_id
  !============================================================================
  integer function irealization_id(GridID_)
    integer,intent(in)::GridID_
    !--------------------------------------------------------------------------
    irealization_id = Domain_I(GridID_)%Ptr%iRealization
  end function irealization_id
  !============================================================================
  subroutine associate_dd_pointer_id(GridID_,DomainPointer)
    integer,intent(in)::GridID_
    type(DomainPointerType),intent(out)::DomainPointer
    !--------------------------------------------------------------------------
    nullify(DomainPointer%Ptr)
    DomainPointer%Ptr=>Domain_I(GridID_)%Ptr
  end subroutine associate_dd_pointer_id
  !============================================================================
 
end module CON_grid_storage
!==============================================================================

