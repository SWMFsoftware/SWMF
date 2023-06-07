!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module CON_grid_descriptor

  use ModUtilities, ONLY: check_allocate, CON_stop
  use ModKind
  use ModNumConst, ONLY: cTiny

  ! This file presents the class of the grid descriptors which
  ! includes both the uniformly spaced grids (uniformly
  ! spaced with respect to some generalized coordinates) and Octree
  ! or Quadric tree for adaptive block grid.
  !
  ! The methods include the grid descriptor allocation, coordinate
  ! computations and the interpolation procedures
  use CON_domain_decomposition, ONLY: DomainType, DomainPointerType,    &
       is_left_boundary_d, is_right_boundary_d,                         &
       glue_margin, search_cell, search_in,  l_neighbor,                &
       None_,  pe_and_blk, cAlmostOne, PE_, BLK_,  FirstChild_
  use CON_grid_storage, ONLY: ndim_id, associate_dd_pointer, ncell_id

  implicit none

  PRIVATE ! Except

  public :: GridType
  type GridType
     ! The concept of the grid descriptor is very close to the domain
     ! decomposition, the pointer for the DomainType
     ! structure is the most important element of the grid descriptor.

     type(DomainPointerType)::Domain

     ! On the other hand, there is a difference too between them.
     ! The domain decomposition consists of domains, that is, volumes.
     ! These volumes can intersect only with with their faces, edges,
     ! corners, but there can be no volume overlapping between the
     ! cells or between the blocks within the domain decomposition. To
     ! the contrary, the grid, first of all, is the set of points,
     ! rather than volumes. Even the grid for the control voulme
     ! numerical method, which is widely used for solving the
     ! conservation laws hyperbolic system,should be thought of as the
     ! set of the cell centered points in the context of the methods
     ! used here. Introduce the concept of a "basic grid".
     !
     ! basic grid[domain decomposition]:=the set of the cell centered
     !                                  points for the given domain
     !                                  decomposition.
     !
     ! Again, herewithin ":=" means "according to the definition"
     !
     ! The basic grid points can be enumerated using iCB (i stands
     ! for index, CB stands   for Cell+Block).
     !
     ! iCB:=(Cell Index, Block Index)
     ! Block Index:= TreeNode when all the grid points are
     !              enumerated.or.Local block number, when only those
     !              grid points are enumerated which belongs to the
     !              piece of grid, which is associated with one of
     !              those, blocks which are allocated at the given PE
     ! Cell Index = nDim vector index (like i,j,k for 3D grid) which
     !              enumerates the grid points throughout the piece
     !              of grid which belongs to a given block
     ! Thus, for a basic grid the  block index is
     ! (1:nTreeNodes(Domain%Ptr), MASK=IsUsedBlock(1:nTreeNodes)
     !
     ! The Cell Index for a basic grid is (1:n_cells_d(Domain%Ptr))
     !
     ! Typically the data sets in the basic grid are shaped as arrays
     ! like
     ! real,dimension(n_cells_D(1),n_cells_d(2),...,&
     !     min_block(Domain%Ptr):max_block(Domain%Ptr))::DataSet_CB
     !
     ! The space of basic grids can be one-to-one mapped to the space
     ! of domain decompositions. To use THE OTHER grids, the grid
     ! descriptor concept should be extended.
     ! Recursive definition:
     ! extended grid:= basic grid.or.cell index
     !                extension[extended grid]
     ! where the cell index extension means the arbitrary decrease in
     ! the value of the minimal range of the cell index for any of the
     ! spatial direction and/or an arbitrary increase in the maximal
     ! one  The extension value can be different for different
     ! spatial directions but should be the same for all the blocks.

     integer::nDim           ! Difenes the grid dimensionality

     integer, pointer:: iPointMin_D(:)
     integer, pointer:: iPointMax_D(:)

     ! The block index IS NOT AFFECTED by this transformation and this
     ! is of crucial importance. Although the grid points of the
     ! extended grid at some value of block index= nTreeNode
     ! can belong to the domain with differing block number, we still
     ! relate these points to the block index value= nTreeNode.
     ! This can be of sense only if the complete data sets at the
     ! extended point set is asessible at the PE, to which the block
     ! nTreeNode is assigned.
     ! Typically the data sets in the extended grid are shaped as
     ! arrays, like as follows:
     ! real,dimension(iPointMin_D(1):iPointMax_D(1),&
     !               iPointMin_D(2):iPointMax_D(2),&
     !         1    ...,&
     !     min_block(Domain%Ptr):max_block(Domain%Ptr))::DataSet_GB
     !
     ! Now introduce the grids which are not cell centered ones.
     !
     ! All grids described by the present grid descriptor:= &
     !                           extended grid.or.&
     !                           displaced[extended grid]
     ! where the displacement is governed by the vector as follows:
     real, pointer::Displacement_D(:)
     ! For an arbitrary grid, for a given block the actual geometric
     ! displacement of the grid block fragment is equal to the product
     ! of Displacement_D by the grid size
  end type GridType

  ! which can be introduced for some standard grids, like
  integer,parameter, public ::CellCentered_=1,Nodes_=2
  public :: set_standard_grid_descriptor
  interface set_standard_grid_descriptor
     module procedure set_standard_grid_descriptor_id
     module procedure set_standard_grid_descriptor_dd
  end interface set_standard_grid_descriptor

  ! For a given global grid point number (in in4 or int8 format)
  ! which enumerates all points in the described grid and for a
  ! given global grid descriptor, the routine returns
  ! a global node number and nDim grid point numbers in a block

  ! With a local grid descriptor, for a given local grid point
  ! number in a format int4, which enumerates all grid points in
  ! blocks allocated at the current PE, the procedure returns the
  ! order number of used block, which enumerates all used blocks at
  ! the given PE in the same order the all blocks in the whole domain
  ! decomposition are enumerated. The nDim grid point indexes are
  ! returned too. For a global grid point number in int8 format,
  ! it returns the same output.
  public :: global_i_grid_point_to_icb
  interface global_i_grid_point_to_icb
     module procedure global_i_grid_point_to_icb4_l
     module procedure global_i_grid_point_to_icb8_l
  end interface global_i_grid_point_to_icb

  ! Returns coordinate of the grid point, in the global or
  ! local grid descriptor, for given global tree node or used
  ! block number accordingly and for given point indexes in the
  ! block
  public :: coord_grid_d
  interface coord_grid_d
     module procedure coord_grid_d_global
     module procedure coord_grid_d_local
  end interface coord_grid_d

  ! local storage for a grid descriptor passed to interpolation_amr_gc;
  ! it is shared by interpolate_amr_gc and find_amr
  type(GridType) :: GridAMR

  ! Local Grid Descriptor
  ! The local Grid describes in a compact way the parts of grid allocated
  ! at a given PE

  ! Order of indexes in the array
  ! Declared in CON_domain_decomposition:
  ! PE_ = 3, BLK_=2,
  public :: BLK_
  integer, parameter, public :: GlobalBlock_ = 4
  integer, parameter, public :: TreeNode_ = 1, GridPointFirst_ = 3

  type LocalGridType

     ! These four members are just copied from the global Grid
     integer           :: nDim           ! Defines the grid dimensionality
     integer, pointer  :: iPointMin_D(:)
     integer, pointer  :: iPointMax_D(:)
     real,    pointer  :: Displacement_D(:)

     ! The following terms are new and characterize the part of
     ! the global grid allocated at the given process
     integer:: nBlock, nPointPerBlock
     integer, pointer ::  iIndex_IB(:,:)
     real,    pointer ::  CoordBlock_DB(:,:), DCoord_DB(:,:)
  end type LocalGridType
  public :: set_local_gd, LocalGridType
  public :: clean_gd
  interface clean_gd
     module procedure clean_grid_descriptor
     module procedure clean_grid_descriptor_l
  end interface clean_gd
  public :: n_grid_points_per_block

  ! Interpolation procedures
  public :: nearest_grid_points    ! First order interpolation
  public :: bilinear_interpolation ! For uniform or node-based grid
  public :: interpolation_amr_gc   ! Employs ghostcells

contains
  !============================================================================
  function coord_grid_d_global(Grid, iTreeNode, iPoints_D)

    ! the coordintes of the grid point

    type(GridType),intent(in)::Grid
    integer, intent(in)::iTreeNode
    integer, intent(in):: iPoints_D(Grid%nDim)
    real:: coord_grid_d_global(Grid%nDim) ! return value
    !--------------------------------------------------------------------------
    coord_grid_d_global = Grid%Domain%Ptr%CoordBlock_DI(:,iTreeNode) + &
         Grid%Domain%Ptr%dCoordCell_DI(:,iTreeNode)*&
         (Grid%Displacement_D - 0.50 + real(iPoints_D))

  end function coord_grid_d_global
  !============================================================================
  function coord_grid_d_local(LocalGrid, iBlockUsed, iPoints_D)

    type(LocalGridType), intent(in) :: LocalGrid
    integer,             intent(in) :: iBlockUsed
    integer,             intent(in) :: iPoints_D(LocalGrid%nDim)
    real   :: coord_grid_d_local(LocalGrid%nDim) ! return value

    integer:: nDim
    !--------------------------------------------------------------------------
    nDim = LocalGrid%nDim
    coord_grid_d_local = LocalGrid%CoordBlock_DB(1:nDim,iBlockUsed) +&
         LocalGrid%DCoord_DB(1:nDim,iBlockUsed)*&
         (LocalGrid%Displacement_D - 0.50 + real(iPoints_D))

  end function coord_grid_d_local
  !============================================================================
  subroutine set_standard_grid_descriptor_id(&
       iGridID, nGhostGridPoints, iStandard, Grid)

    ! Allow to set the standard grid descriptor:
    ! iStandard\_=CellCentered\_ or iStandard\_=Nodes\_,
    ! with or without halo points ("ghost points")

    integer,intent(in)::iGridID
    integer,intent(in),optional::iStandard
    integer,intent(in),optional::nGhostGridPoints
    type(GridType),intent(out)::Grid
    integer::iError, iMyStandard, nGhostGridPointsMy
    !--------------------------------------------------------------------------
    call associate_dd_pointer(iGridID, Grid%Domain)
    Grid%nDim=ndim_id(iGridID)
    allocate(Grid%Displacement_D(1:Grid%nDim),stat=iError)
    call check_allocate(iError,"Displacement_D")
    allocate(Grid%iPointMin_D(1:Grid%nDim),stat=iError)
    call check_allocate(iError,"iPointMin_D")
    allocate(Grid%iPointMax_D(1:Grid%nDim),stat=iError)
    call check_allocate(iError,"iPointMax_D")
    if(present(nGhostGridPoints))then
       nGhostGridPointsMy=nGhostGridPoints
    else
       nGhostGridPointsMy=0
    end if
    if(present(iStandard))then
       iMyStandard=iStandard
    else
       iMyStandard=CellCentered_
    end if
    select case(iMyStandard)
    case(CellCentered_)
       Grid%Displacement_D= 0.0
       Grid%iPointMin_D = 1 - nGhostGridPointsMy
       Grid%iPointMax_D = ncell_id(iGridID) + &
            nGhostGridPointsMy
    case(Nodes_)
       Grid%Displacement_D=-0.50
       Grid%iPointMin_D = min(1,2-nGhostGridPointsMy)
       Grid%iPointMax_D = ncell_id(iGridID) + &
            nGhostGridPointsMy
    case default
       call CON_stop('Unknown standard for Grid Descriptor')
    end select

  end subroutine set_standard_grid_descriptor_id
  !============================================================================
  subroutine set_standard_grid_descriptor_dd(&
       Domain, nGhostGridPoints, iStandard, Grid)

    type(DomainType), target,intent(in)::Domain
    integer, intent(in), optional::iStandard
    integer, intent(in), optional::nGhostGridPoints
    type(GridType), intent(out)::Grid

    integer:: iError, iMyStandard, nGhostGridPointsMy
    !--------------------------------------------------------------------------
    nullify(Grid%Domain%Ptr)
    call associate_dd_pointer(Domain, Grid%Domain)
    Grid%nDim=Domain%nDim
    allocate(Grid%Displacement_D(&
         1:Grid%nDim),stat=iError)
    call check_allocate(iError,"Displacement_D")
    allocate(Grid%iPointMin_D(&
         1:Grid%nDim),stat=iError)
    call check_allocate(iError,"iPointMins_D")
    allocate(Grid%iPointMax_D(&
         1:Grid%nDim),stat=iError)
    call check_allocate(iError,"iPointMaxs_D")
    if(present(nGhostGridPoints))then
       nGhostGridPointsMy=nGhostGridPoints
    else
       nGhostGridPointsMy=0
    end if
    if(present(iStandard))then
       iMyStandard=iStandard
    else
       iMyStandard=CellCentered_
    end if
    select case(iMyStandard)
    case(CellCentered_)
       Grid%Displacement_D= 0.0
       Grid%iPointMin_D = 1 - nGhostGridPointsMy
       Grid%iPointMax_D = Domain%nCell_D + nGhostGridPointsMy
    case(Nodes_)
       Grid%Displacement_D = -0.50
       Grid%iPointMin_D = min(1,2-nGhostGridPointsMy)
       Grid%iPointMax_D = Domain%nCell_D + nGhostGridPointsMy
    case default
       call CON_stop('Unknown standard for Grid Descriptor')
    end select

  end subroutine set_standard_grid_descriptor_dd
  !============================================================================
  subroutine set_local_gd(iProc, Grid, LocalGrid)

    ! PE rank (in a local group for a local model, or in a
    ! global group for a global model)

    integer,                 intent(in)  :: iProc
    type(GridType),          intent(in)  :: Grid
    type(LocalGridType),     intent(out) :: LocalGrid

    integer :: nDim, nPointPerBlock, iBlock, nBlock, iBlockAll
    integer :: iTreeNode, nTreeNode, iError
    integer, pointer:: iDomain_II(:,:)
    !--------------------------------------------------------------------------
    nDim                     = Grid%nDim
    LocalGrid%nDim           = nDim
    LocalGrid%Displacement_D => Grid%Displacement_D
    LocalGrid%iPointMin_D => Grid%iPointMin_D
    LocalGrid%iPointMax_D => Grid%iPointMax_D
    ! For better readability
    iDomain_II => Grid%Domain%Ptr%iDD_II
    nTreeNode  = Grid%Domain%Ptr%nTreeNode
    ! Number of used blocks on a given process
    nBlock =  count(iDomain_II(FirstChild_, 1:nTreeNode)==None_ &
         .and. iDomain_II(PE_, 1:nTreeNode)==iProc)
    LocalGrid%nBlock           = nBlock
    nPointPerBlock =    n_grid_points_per_block(Grid)
    LocalGrid%nPointPerBlock   = nPointPerBlock

    allocate(LocalGrid%iIndex_IB(TreeNode_:GlobalBlock_,&
         1:nBlock),stat=iError)
    call check_allocate(iError,'LocalGrid%iIndex_IB')
    allocate(LocalGrid%DCoord_DB(1:nDim,1:nBlock), stat=iError)
    call check_allocate(iError,'LocalGrid%DCoord_DB')
    allocate(LocalGrid%CoordBlock_DB(1:nDim,1:nBlock),stat=iError)
    call check_allocate(iError,'LocalGrid%CoordBlock_DB')
    ! Constrain the domain decomposition arrays for a given PE
    iBlock = 0; iBlockAll = 0
    ! Loop over all tree nodes
    do  iTreeNode = 1, nTreeNode
       ! Skip unused blocks
       if(iDomain_II(FirstChild_, iTreeNode) /= None_)CYCLE
       !  iBlockAll enumerates all used blocks
       iBlockAll = iBlockAll + 1
       ! Skip all used blocks allocated at procs other than iProc
       if(iDomain_II(PE_, iTreeNode) /=iProc)CYCLE
       ! iBlock enumerates used blocks on the process iProc
       iBlock = iBlock + 1
       ! Store iTreeNode
       LocalGrid%iIndex_IB(TreeNode_, iBlock) = iTreeNode
       ! Comy and save local block ID as stored in the global DD
       LocalGrid%iIndex_IB(BLK_, iBlock)  = iDomain_II(BLK_, iTreeNode)
       ! Store iBlockAll
       LocalGrid%iIndex_IB(GlobalBlock_, iBlock) = iBlockAll
       ! Store the point index for the first point on this block
       LocalGrid%iIndex_IB(GridPointFirst_,iBlock) = (iBlock - 1)*&
            nPointPerBlock + 1
       ! Store the grid resolution for this block
       LocalGrid%DCoord_DB(:,iBlock) = &
            Grid%Domain%Ptr%DCoordCell_DI(:, iTreeNode)
       ! Store coords of the block left corner
       LocalGrid%CoordBlock_DB(:,iBlock) = &
            Grid%Domain%Ptr%CoordBlock_DI(:, iTreeNode)
    end do

  end subroutine set_local_gd
  !============================================================================
  subroutine clean_grid_descriptor(Grid)

    type(GridType),intent(inout):: Grid
    !--------------------------------------------------------------------------
    deallocate(Grid%iPointMin_D)
    deallocate(Grid%iPointMax_D)
    deallocate(Grid%Displacement_D)
    nullify(Grid%Domain%Ptr)

  end subroutine clean_grid_descriptor
  !============================================================================
  subroutine clean_grid_descriptor_l(Grid)

    type(LocalGridType), intent(inout):: Grid
    !--------------------------------------------------------------------------
    nullify(Grid%iPointMin_D)
    nullify(Grid%iPointMax_D)
    nullify(Grid%Displacement_D)
    deallocate(Grid%iIndex_IB)
    deallocate(Grid%CoordBlock_DB)
    deallocate(Grid%DCoord_DB)

  end subroutine clean_grid_descriptor_l
  !============================================================================
  integer function n_grid_points_per_block(Grid)

    ! Number of cells per block
    type(GridType)::Grid
    !--------------------------------------------------------------------------
    n_grid_points_per_block = product(&
         Grid%iPointMax_D + 1 -Grid%iPointMin_D)

  end function n_grid_points_per_block
  !============================================================================
  subroutine global_i_grid_point_to_icb4_l(&
       Grid, iPointLocal, iBlockUsed, iPoint_D)

    type(LocalGridType),    intent(in) :: Grid
    integer,              intent(in) :: iPointLocal
    integer,             intent(out) :: iBlockUsed
    integer,             intent(out) :: iPoint_D(Grid%nDim)

    ! Local variables

    ! Number of grid points in the block, for each direction
    integer  :: nGridPoints_D(Grid%nDim)
    integer  :: iDim,iMisc,iMisc1

    !--------------------------------------------------------------------------
    nGridPoints_D=1+Grid%iPointMax_D - Grid%iPointMin_D

    iMisc = iPointLocal - 1
    do idim = 1, Grid%nDim
       iMisc1 = mod(iMisc,nGridPoints_D(iDim))
       iPoint_D(iDim) = Grid%iPointMin_D(iDim)+iMisc1
       iMisc=(iMisc-iMisc1)/nGridPoints_D(iDim)
    end do
    iBlockUsed = iMisc + 1
  end subroutine global_i_grid_point_to_icb4_l
  !============================================================================
  subroutine global_i_grid_point_to_icb8_l(&
       Grid,  &
       iPointGlobal,&
       iBlockUsed,      &
       iPoint_D)
    type(LocalGridType), intent(in) :: Grid
    integer(Int8_),    intent(in) :: iPointGlobal
    integer,          intent(out) :: iBlockUsed
    integer,          intent(out) :: &
         iPoint_D(Grid%nDim)

    ! Local variables

    ! Number of grid points in the block, for each direction
    integer(Int8_)   :: nGridPoints_D(Grid%nDim)
    integer          :: iDim, iMisc1, loc(1)
    integer(Int8_)   :: i8Misc
    !--------------------------------------------------------------------------
    nGridPoints_D = 1_Int8_ + Grid%iPointMax_D -&
         Grid%iPointMin_D
    i8Misc =  iPointGlobal - 1
    do iDim = 1, Grid%nDim
       iMisc1 = int(mod(i8Misc, nGridPoints_D(iDim)))
       iPoint_D(iDim) = Grid%iPointMin_D(iDim)+iMisc1
       i8Misc = (i8Misc - iMisc1)/nGridPoints_D(iDim)
    end do
    loc = maxloc(Grid%iIndex_IB(GlobalBlock_,:),&
         MASK = Grid%iIndex_IB(GlobalBlock_,:) <= &
         (int(i8Misc)+1))
    iBlockUsed = loc(1)
  end subroutine global_i_grid_point_to_icb8_l
  !============================================================================
  subroutine nearest_grid_points( &
       nDim,Coord_D, Grid, nIndex, iIndex_II, nImage, Weight_I)

    ! The general principle of the interpoaltions with the grids
    ! defined by the present grid descriptor is as follows.
    ! Using coordinates of the point at which the interpalated value
    ! should be found, first of all, the block is found this point
    ! belongs to. Then, among the grid points belonging to the grid
    ! fragment with the found value of the block index, the nearest
    ! point is found. For the second order of interpolation, the
    ! cubic stenchil is constructed involving the nearest point too
    ! and the interpolation weights are found. If the values of the
    ! cell index for the stenchil points appear to be out of the
    ! block, the correspondent blocks are found

    type(GridType):: Grid
    integer,      intent(in):: nDim
    real,      intent(inout):: Coord_D(nDim)
    integer,      intent(in):: nIndex

    integer,     intent(out):: iIndex_II(0:nIndex,2**nDim)
    integer,     intent(out):: nImage
    real,        intent(out):: Weight_I(2**nDim)

    real:: CoordStored_D(nDim), DCoordCells_D(nDim), DCoordTolerance_D(nDim)

    integer :: iPoints_D(nDim)
    logical :: IsAtFace_D(nDim)
    real    :: Coord_DI(nDim,2**nDim)
    real, parameter:: Tolerance = 0.001
    integer:: iTreeNode, iDim, iImages

    !--------------------------------------------------------------------------
    iIndex_II = 0
    Weight_I  = 0.0
    CoordStored_D = Coord_D
    iIndex_II(1:nDim,1) = Grid%iPointMax_D+1

    do while(any(iIndex_II(1:nDim,1) > Grid%iPointMax_D))
       Coord_D = CoordStored_D

       ! Find a block to which the point belongs
       call search_in(Grid%Domain%Ptr, Coord_D, iTreeNode)
       DCoordCells_D = Grid%Domain%Ptr%dCoordCell_DI(:,iTreeNode)

       ! Check if the point is out of the computational domain
       where(is_right_boundary_d(Grid%Domain%Ptr,iTreeNode))
          Coord_D = min(Coord_D,cAlmostOne*&
               Grid%Domain%Ptr%dCoordBlock_DI(:,iTreeNode))
          CoordStored_D = Coord_D + &
               Grid%Domain%Ptr%CoordBlock_DI(:,iTreeNode)
       end where
       where(is_left_boundary_d(Grid%Domain%Ptr,iTreeNode))
          Coord_D = max(Coord_D, 0.0)
          CoordStored_D = Coord_D + &
               Grid%Domain%Ptr%CoordBlock_DI(:,iTreeNode)
       end where

       ! Store the processor number and block number. If the block
       ! number is meaningless use nIndex=nDim and this meaningless
       ! value will be overwritten later
       call pe_and_blk(Grid%Domain%Ptr,iTreeNode,&
            iIndex_II(0,1),iIndex_II(nIndex,1))

       ! Find the nearest grid point in the block
       Coord_D = Coord_D-DCoordCells_D*Grid%Displacement_D
       iPoints_D = floor(Coord_D/DCoordCells_D)
       iIndex_II(1:nDim,1) = iPoints_D+1
       DCoordTolerance_D = Tolerance*DCoordCells_D

       ! The vector Coord - CoordGrid + 0.5 DCoord
       Coord_D = Coord_D-DCoordCells_D*iPoints_D

       ! If the point coordintates are such that Coord-CoordGrid=+/-0.5*DCoord
       ! there are several equidistant grid points

       IsAtFace_D = abs(Coord_D)<DCoordTolerance_D&
            .or.abs(Coord_D-DCoordCells_D )<DCoordTolerance_D
       where(is_right_boundary_d(Grid%Domain%Ptr,iTreeNode))&
            iIndex_II(1:nDim,1) = &
            min(iIndex_II(1:nDim,1),&
            Grid%iPointMax_D)

       ! The nearest grid points may be are in the upper block
       where(iIndex_II(1:nDim,1)>&
            Grid%iPointMax_D)&
            CoordStored_D = CoordStored_D+DCoordTolerance_D*0.25
    end do
    nImage = 1;Weight_I(1) =  1.0
    if(.not.(any(IsAtFace_D)))RETURN
    Coord_DI(:,1) = CoordStored_D
    do iDim = 1,nDim
       if(IsAtFace_D(iDim))then
          Coord_DI(:,1+nImage:nImage+nImage) = Coord_DI(:,1:nImage)
          Coord_DI(iDim,1:nImage) = Coord_DI(iDim,1:nImage)-&
               DCoordTolerance_D(iDim)
          Coord_DI(iDim,1+nImage:nImage+nImage) = Coord_DI(&
               iDim,1+nImage:nImage+nImage)+&
               DCoordTolerance_D(iDim)
          nImage = nImage+nImage
       end if
    end do
    do iImages = 1,nImage
       call search_in(Grid%Domain%Ptr,&
            Coord_DI(:,iImages),&
            iTreeNode)
       call pe_and_blk(Grid%Domain%Ptr,&
            iTreeNode,&
            iIndex_II(0,iImages),&
            iIndex_II(nIndex,iImages))
       Coord_DI(:,iImages) = Coord_DI(:,iImages)-&
            Grid%Domain%Ptr%dCoordCell_DI(:,&
            iTreeNode)*Grid%Displacement_D
       call search_cell(Grid%Domain%Ptr,&
            iTreeNode,Coord_DI(:,iImages),&
            iIndex_II(1:nDim,iImages))
    end do
    iImages = 1
    do while(iImages<=nImage)
       ! Exclude the stencil nodes which are out of the
       ! computational domain
       do while(any(Grid%iPointMin_D>&
            iIndex_II(1:nDim,iImages).or.&
            Grid%iPointMax_D<&
            iIndex_II(1:nDim,iImages)))
          if(iImages == nImage)then
             nImage = nImage-1
             EXIT
          end if
          iIndex_II(:,iImages) = iIndex_II(:,nImage)
          nImage = nImage-1
       end do
       iImages = iImages+1
    end do
    Weight_I(1:nImage) = 1.0/real(nImage)

  end subroutine nearest_grid_points
  !============================================================================
  subroutine bilinear_interpolation(&
       nDim, Coord_D, Grid, nIndex, iIndex_II, nImage, Weight_I)

    ! SECOND ORDER INTERPOLATION
    ! This is a bilinear interpolation using the grid points,
    ! described with the grid descriptor.

    integer, intent(in)   :: nDim
    type(GridType)        :: Grid
    real,    intent(inout):: Coord_D(nDim)
    integer, intent(in)   :: nIndex
    integer, intent(out)  :: iIndex_II(0:nIndex,2**nDim)
    integer, intent(out)  :: nImage

    real, intent(out):: Weight_I(2**nDim)

    real   :: CoordResid_D(nDim), CoordStored_D(nDim)
    integer:: iPoints_D(nDim)
    real   :: CoordMisc_D(nDim)
    integer:: lNode_I(2**nDim)
    integer:: iTreeNode, iDim, iImages, iNewStart, nImageNew
    real   :: WeightLeft
    logical:: IsUp_DI(nDim,2**nDim), IsDown_DI(nDim,2**nDim)
    real   :: DCoord_DI(nDim,0:2**nDim)
    logical:: IsDomainBoundaryUp_D(nDim), IsDomainBoundaryDown_D(nDim)
    !--------------------------------------------------------------------------
    iIndex_II = 0
    Weight_I  = 0.0
    CoordStored_D = Coord_D
    IsUp_DI(:,1) = .true.
    call search_in(Grid%Domain%Ptr, Coord_D, iTreeNode)
    ! Find global node number, PE and number which involved the point
    do while(any(IsUp_DI(:,1)))
       ! This do loop works more than once only in case of node
       ! grid and only in case of Coord_D point belonging to a block
       ! boundary. At this case the routine needs help in deciding,
       ! to which block this point should be assigned. Otherwise,
       ! automatically IsUp_D(:,1)=.false. after first loop pass.
       Coord_D = CoordStored_D - &
            Grid%Domain%Ptr%CoordBlock_DI(:,iTreeNode)
       IsDomainBoundaryUp_D =   &
            is_right_boundary_d(Grid%Domain%Ptr,iTreeNode)
       IsDomainBoundaryDown_D = &
            is_left_boundary_d(Grid%Domain%Ptr,iTreeNode)
       where(IsDomainBoundaryUp_D)
          Coord_D = min(Coord_D, cAlmostOne*&
               Grid%Domain%Ptr%dCoordBlock_DI(:,iTreeNode))
          CoordStored_D = Coord_D + &
               Grid%Domain%Ptr%CoordBlock_DI(:,iTreeNode)
       end where
       where(IsDomainBoundaryDown_D)
          Coord_D = max(Coord_D, 0.0)
          CoordStored_D = Coord_D + &
               Grid%Domain%Ptr%CoordBlock_DI(:,iTreeNode)
       end where

       call pe_and_blk(Grid%Domain%Ptr,iTreeNode,&
            iIndex_II(0,1),iIndex_II(nIndex,1))

       ! Find DCoord for this block
       DCoord_DI(:,0) = Grid%Domain%Ptr%dCoordCell_DI(:,iTreeNode)

       CoordResid_D = Coord_D/DCoord_DI(:,0) - Grid%Displacement_D + 0.50
       iPoints_D = floor(CoordResid_D)
       CoordResid_D = CoordResid_D - real(iPoints_D)
       ! Thus calculated CoordResid_D satisfies the inequalities as       !
       ! follow:CoordResid_D>=0 and CoordResid_D<1. It is used to calculte  !
       ! weights for the eight grid points, among them the iPoints_D!
       ! being the left one with respect to all the spatial directions  !
       iIndex_II(1:nDim,1) = iPoints_D
       IsUp_DI(:,1) = iIndex_II(1:nDim,1) > Grid%iPointMax_D
       IsDown_DI(:,1) = iIndex_II(1:nDim,1) < Grid%iPointMin_D
       if(any(IsUp_DI(:,1)))iTreeNode = l_neighbor(&
            Grid%Domain%Ptr, iTreeNode, iIndex_II(1:nDim,1))
    end do

    nImage = 1; Weight_I(1) =  1.0

    do iDim = 1, nDim
       ! Exclude the stencil nodes which are out of the
       ! computational domain
       if(IsDown_DI(iDim,1).and.IsDomainBoundaryDown_D(iDim))then
          IsDown_DI(iDim,1:nImage) = .false.
          iIndex_II(iDim,1:nImage) = iIndex_II(iDim,1:nImage) + 1
          CYCLE
       end if

       if(iIndex_II(iDim,1)==Grid%iPointMax_D(iDim)&
            .and.IsDomainBoundaryUp_D(iDim))CYCLE
       if(CoordResid_D(iDim)<cTiny)CYCLE
       iNewStart = nImage + 1; nImageNew = nImage + nImage
       iIndex_II(:,iNewStart:nImageNew) = iIndex_II(:,1:nImage)
       iIndex_II(iDim,iNewStart:nImageNew) = &
            iIndex_II(iDim,iNewStart:nImageNew) + 1
       IsUp_DI(:,iNewStart:nImageNew) = IsUp_DI(:,1:nImage)
       IsUp_DI(iDim,iNewStart:nImageNew) = &
            iIndex_II(iDim,iNewStart) > Grid%iPointMax_D(iDim)
       IsDown_DI(:,iNewStart:nImageNew) = IsDown_DI(:,1:nImage)
       IsDown_DI(iDim,iNewStart:nImageNew) = .false.
       Weight_I(iNewStart:nImageNew) =  Weight_I(1:nImage)*CoordResid_D(iDim)
       WeightLeft =  1.0 -CoordResid_D(iDim)
       Weight_I(1:nImage) =  WeightLeft* Weight_I(1:nImage)
       nImage = nImageNew
    end do

    ! Check if the grid point index is within the index limits
    do iImages = 1,nImage
       if(.not.(any(IsUp_DI(:,iImages)).or.&
            any(IsDown_DI(:,iImages))))then
          DCoord_DI(:,iImages) = DCoord_DI(:,0)
          lNode_I(iImages) = iTreeNode
       else
          lNode_I(iImages) = l_neighbor(Grid%Domain%Ptr,&
               iTreeNode,iIndex_II(1:nDim,iImages))
          CoordMisc_D = coord_grid_d(Grid,&
               iTreeNode,iIndex_II(&
               1:nDim,iImages)) - &
               Grid%Domain%Ptr%CoordBlock_DI(:,lNode_I(iImages))
          call pe_and_blk(Grid%Domain%Ptr,lNode_I(iImages),&
               iIndex_II(0,iImages),iIndex_II(nIndex,iImages))
          DCoord_DI(:,iImages) = &
               Grid%Domain%Ptr%dCoordCell_DI(:,lNode_I(iImages))
          CoordMisc_D = CoordMisc_D - DCoord_DI(:,iImages)&
               *Grid%Displacement_D
          call search_cell(Grid%Domain%Ptr, lNode_I(iImages), CoordMisc_D,&
               iIndex_II(1:nDim,iImages))
       end if
    end do

  end subroutine bilinear_interpolation
  !============================================================================
  subroutine interpolation_amr_gc(&
       nDim, Coord_D, Grid, nIndex, iIndex_II, nImage, Weight_I)

    ! utilizing ghost cells
    ! This is a continuous amr interpolation using the grid points,
    ! described with the grid descriptor. It utilizes ghost cells.

    use ModInterpolateAMR, ONLY: interpolate_amr

    integer,   intent(in)   :: nDim
    type(GridType):: Grid
    real,      intent(inout):: Coord_D(nDim)
    integer,   intent(in)   :: nIndex
    ! OUTPUT ARGUMENTS
    integer,   intent(out)  :: iIndex_II(0:nIndex,2**nDim)
    integer,   intent(out)  :: nImage
    real,      intent(out)  :: Weight_I(2**nDim)

    ! memorize grid descriptor by storing it in modular variable;
    ! it is used by find_amr subroutine passed to shared interpolation routine
    !--------------------------------------------------------------------------
    GridAMR = Grid
    ! call shared interpolation subroutine
    call interpolate_amr(&
         nDim        = nDim,&
         XyzIn_D     = Coord_D,&
         nIndexes    = nIndex,&
         find        = find_amr,&
         nCell_D     = GridAMR%Domain%Ptr%nCell_D,&
         nGridOut    = nImage,&
         Weight_I    = Weight_I,&
         iIndexes_II = iIndex_II,&
         UseGhostCell= .true.)

  end subroutine interpolation_amr_gc
  !============================================================================
  subroutine find_amr(nDim, Coord_D, &
       iProc, iBlock, CoordCorner_D, DCoord_D, IsOut)

    use ModNumConst

    integer, intent(in)   :: nDim
    real,    intent(inout):: Coord_D(nDim)
    integer, intent(out)  :: iProc, iBlock
    real,    intent(out)  :: CoordCorner_D(nDim), DCoord_D(nDim)
    logical, intent(out)  :: IsOut

    integer:: iNode
    logical:: IsPeriodic_D(nDim)
    !--------------------------------------------------------------------------
    ! check if point's inside the domain
    IsPeriodic_D = GridAMR%Domain%Ptr%IsPeriodic_D
    IsOut = &
         any(GridAMR%Domain%Ptr%CoordMin_D >  Coord_D&
         .and..not.IsPeriodic_D) .or. &
         any(GridAMR%Domain%Ptr%CoordMax_D <= Coord_D&
         .and..not.IsPeriodic_D)
    if(IsOut) then
       if(GridAMR%Domain%Ptr%DoGlueMargins)then
          ! Recalculate coordinates for spherical or similar
          ! grids
          call glue_margin(GridAMR%Domain%Ptr, Coord_D)
          IsOut = &
               any(GridAMR%Domain%Ptr%CoordMin_D >  Coord_D&
               .and..not.IsPeriodic_D) .or. &
               any(GridAMR%Domain%Ptr%CoordMax_D <= Coord_D&
               .and..not.IsPeriodic_D)
          if(IsOut)RETURN
       else
          RETURN
       end if
    end if
    ! find a node containing Coord_D
    ! this call changes Coord_D to coords relative to block corner as needed
    call search_in(GridAMR%Domain%Ptr, Coord_D, iNode)

    ! now extract the block parameters
    CoordCorner_D = GridAMR%Domain%Ptr%CoordBlock_DI(:, iNode)
    DCoord_D      = GridAMR%Domain%Ptr%DCoordCell_DI(:, iNode)
    iProc         = GridAMR%Domain%Ptr%iDD_II(PE_, iNode)
    iBlock        = GridAMR%Domain%Ptr%iDD_II(BLK_, iNode)

  end subroutine find_amr
  !============================================================================
end module CON_grid_descriptor
!==============================================================================
