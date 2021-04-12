!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module CON_grid_storage
  use CON_world
  use CON_domain_decomposition, ONLY: DomainType, DomainPointerType, &
       init_decomposition_dd, get_root_decomposition_dd,  &
       bcast_decomposition_dd, synchronize_refinement_dd, &
       search_cell, search_in, &
       associate_dd_pointer_dd, glue_margin_dd, iglobal_block_dd, &
       iglobal_node_dd
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

    integer,           intent(in) :: GridID_
    integer,           intent(in) :: CompID_,nDim
    logical, optional, intent(in) :: IsTreeDD
    !--------------------------------------------------------------------------
    if(GridID_ > MaxGrid .or. GridID_ <=0 )call CON_stop(&
         'Prrohibited value for GridID_')
    if(done_dd_init(GridID_))call CON_stop(&
         'An attempt to reinitialize the GridID_')
    call init_grid_storage(Domain_I, GridID_)
    call init_decomposition_dd(Domain_I(GridID_)%Ptr, CompID_, nDim, IsTreeDD)
    DoneDomainInit_C(GridID_) = .true.;  nDim_C(GridID_) = nDim
  end subroutine init_decomposition_id
  !============================================================================
  ! To get a decomposition domain, even the tree one, the root     
  ! decomposition should be first constructed. PE here are the ranks in the
  ! LOCAL communicator for the component
  subroutine get_root_decomposition_id(&
       GridID_,      &! ID for Decomposition to be constructed
       iRootMapDim_D,&! As in DomainType
       CoordMin_D,   &! As in DomainType
       CoordMax_D,   &! As in DomainType
       nCell_D,      &! As in DomainType
       PE_I,         &! PE layout
       iBlock_I,     &! Local Block Number layout
       IsPeriodic_D, &! As in DomainType
       iShift_DI,    &! As in DomainType
       DoGlueMargins,&! As in DomainType
       iDirMinusGlue,&! As in DomainType
       iDirPlusGlue,& ! As in DomainType
       iDirCycle)     ! As in DomainType)

    integer, intent(in) :: GridID_
    integer, intent(in) :: iRootMapDim_D(nDim_C(GridID_))
    real,    intent(in) :: CoordMin_D(nDim_C(GridID_))
    real,    intent(in) :: CoordMax_D(nDim_C(GridID_))
    integer, intent(in) :: nCell_D(nDim_C(GridID_))
    integer, optional, intent(in) :: PE_I(:)
    integer, optional, intent(in) :: iBlock_I(:)
    logical, optional, intent(in) :: IsPeriodic_D(nDim_C(GridID_))
    integer, optional, intent(in) :: iShift_DI(:,:)
    logical, optional, intent(in) :: DoGlueMargins
    integer, optional, intent(in) :: iDirMinusGlue
    integer, optional, intent(in) :: iDirPlusGlue
    integer, optional, intent(in) :: iDirCycle
    !--------------------------------------------------------------------------
    call get_root_decomposition_dd(Domain_I(GridID_)%Ptr, &
       iRootMapDim_D,&! As in DomainType
       CoordMin_D,   &! As in DomainType
       CoordMax_D,   &! As in DomainType
       nCell_D,      &! As in DomainType
       PE_I,         &! PE layout
       iBlock_I,     &! Local Block Number layout
       IsPeriodic_D, &! As in DomainType
       iShift_DI,    &! As in DomainType
       DoGlueMargins,&! As in DomainType
       iDirMinusGlue,&! As in DomainType
       iDirPlusGlue,& ! As in DomainType
       iDirCycle)     ! As in DomainType)
  end subroutine get_root_decomposition_id
  !============================================================================
  subroutine bcast_decomposition_id(GridID_)

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
       GridID_, LocalDomain, iProcUnion, iCommUnion)

    integer,           intent(in) :: GridID_      ! Global grid ID
    type(DomainType),  intent(in) :: LocalDomain  ! Local Grid with which the
    integer, optional, intent(in) :: iProcUnion,iCommUnion
    !-------------------------------------------------------------------------- 
    call synchronize_refinement_dd(&
            Domain_I(GridID_)%Ptr, LocalDomain, iProcUnion, iCommUnion)
  end subroutine synchronize_refinement_id
  !============================================================================
  subroutine glue_margin_id(GridID_, Coord_D)
    integer,         intent(in) :: GridID_
    real,         intent(inout) :: Coord_D(nDim_C(GridID_))
    !--------------------------------------------------------------------------
    call glue_margin_dd(Domain_I(GridID_)%Ptr,Coord_D)
  end subroutine glue_margin_id
  !============================================================================
  integer function iglobal_node_id(GridID_,iBlockAll)
    integer, intent(in) :: GridID_
    integer, intent(in) :: iBlockAll
    !--------------------------------------------------------------------------
    iglobal_node_id = Domain_I(GridID_)%Ptr%iGlobal_A(iBlockAll)
  end function iglobal_node_id
  !============================================================================
  integer function iglobal_block_id(GridID_,  iTreeNode)
    use CON_domain_decomposition, ONLY: GlobalBlock_
    integer, intent(in) :: GridID_
    integer, intent(in) :: iTreeNode
    !--------------------------------------------------------------------------
    iglobal_block_id=Domain_I(GridID_)%Ptr%iDD_II(GlobalBlock_,iTreeNode)
  end function iglobal_block_id
  !============================================================================
  integer function i_realization(GridID_)
    integer, intent(in) :: GridID_
    !--------------------------------------------------------------------------
    i_realization = Domain_I(GridID_)%Ptr%iRealization
  end function i_realization
  !============================================================================
  subroutine associate_dd_pointer_id(GridID_,DomainPointer)
    integer, intent(in) :: GridID_
    type(DomainPointerType),intent(out)::DomainPointer
    !--------------------------------------------------------------------------
    nullify(DomainPointer%Ptr)
    DomainPointer%Ptr=>Domain_I(GridID_)%Ptr
  end subroutine associate_dd_pointer_id
  !============================================================================
  integer function ndim_id(GridID_)
    integer, intent(in) :: GridID_
    !--------------------------------------------------------------------------
    ndim_id = Domain_I(GridID_)%Ptr%nDim
  end function ndim_id
  !============================================================================
  function ncell_id(GridID_)
    integer,intent(in) :: GridID_
    integer :: ncell_id(nDim_C(GridID_))
    !--------------------------------------------------------------------------
    ncell_id = Domain_I(GridID_)%Ptr%nCell_D
  end function ncell_id
  !============================================================================
  integer function compid(GridID_)
    integer,intent(in)::GridID_
    !--------------------------------------------------------------------------
    compid = Domain_I(GridID_)%Ptr%CompID_
  end function compid
  !============================================================================
end module CON_grid_storage
!==============================================================================

