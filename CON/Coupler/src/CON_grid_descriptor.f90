!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!BOP
!MODULE: CON_grid_descriptor - for uniform or octree grids
!INTERFACE:
Module CON_grid_descriptor
  use ModKind
  !This file presents the class of the grid descriptors which     
  !includes both the uniformly spaced grids (uniformly    
  !spaced with respect to some generalized coordinates) and Octree
  ! or Quadric tree for adaptive block grid.                      
  !                                                               
  !The methods include the grid descriptor allocation, coordinate  
  !computations and the interpolation procedures                  

  !USES:
  use CON_grid_storage
  !REVISION HISTORY:
  ! Sokolov I.V.                                                  
  ! 6.18.03-7.11.03                                               
  ! igorsok@umich.edu                                             
  ! phone(734)647-4705                                            
  !===============================================================
  implicit none
  PRIVATE !Except                                 
  !=====================DERIVED TYPE==============================
  public :: GridDescriptorType 
  type GridDescriptorType

     !The concept of the grid descriptor is very close to the domain 
     !decomposition, the pointer for the DomainDecompositionType     
     !structure is the most important element of the grid descriptor.

     type(DDPointerType)::DD

     !On the other hand, there is a difference too between them.     
     !The domain decomposition consists of domains, that is, volumes.
     !These volumes can intersect only with with their faces, edges, 
     !corners, but there can be no volume overlapping between the    
     !cells or between the blocks within the domain decomposition. To
     !the contrary, the grid, first of all, is the set of points,    
     !rather than volumes. Even the grid for the control voulme      
     !numerical method, which is widely used for solving the         
     !conservation laws hyperbolic system,should be thought of as the  
     !set of the cell centered points in the context of the methods  
     !used here. Introduce the concept of a "basic grid".            
     !                                                               
     !basic grid[domain decomposition]:=the set of the cell centered 
     !                                  points for the given domain  
     !                                  decomposition.               
     !                                                               
     !Again, herewithin ":=" means "according to the definition"     
     !                                                               
     !The basic grid points can be enumerated using iCB (i stands    
     !for index, CB stands   for Cell+Block).                        
     !                                                               
     !iCB:=(Cell Index, Block Index)                                 
     !Block Index:= GlobalTreeNode when all the grid points are      
     !              enumerated.or.Local block number, when only those
     !              grid points are enumerated which belongs to the  
     !              piece of grid, which is associated with one of   
     !              those, blocks which are allocated at the given PE
     !Cell Index = nDim vector index (like i,j,k for 3D grid) which  
     !              enumerates the grid points throughout the piece  
     !              of grid which belongs to a given block           
     !Thus, for a basic grid the  block index is                     
     ! (1:nTreeNodes(DD%Ptr), MASK=IsUsedBlock(1:nTreeNodes)         
     !                                                               
     !The Cell Index for a basic grid is (1:n_cells_d(DD%Ptr))       
     !                                                               
     !Typically the data sets in the basic grid are shaped as arrays 
     !like                                                           
     !real,dimension(1:n_cells_D(1),1:n_cells_d(2),...,&             
     !     min_block(DD%Ptr):max_block(DD%Ptr))::DataSet_CB          
     !                                                               
     !The space of basic grids can be one-to-one mapped to the space 
     !of domain decompositions. To use THE OTHER grids, the grid     
     !descriptor concept should be extended.                         
     !Recursive definition:                                          
     !extended grid:= basic grid.or.cell index                        
     !                extension[extended grid]                       
     !where the cell index extension means the arbitrary decrease in 
     !the value of the minimal range of the cell index for any of the
     !spatial direction and/or an arbitrary increase in the maximal  
     ! one  The extension value can be different for different       
     !spatial directions but should be the same for all the blocks.  

     integer::nDim           !Difenes the grid dimensionality   

     integer,dimension(:),pointer::iGridPointMin_D
     integer,dimension(:),pointer::iGridPointMax_D

     !The block index IS NOT AFFECTED by this transformation and this
     !is of crucial importance. Although the grid points of the      
     !extended grid at some value of block index= nGlobalTreeNode    
     !can belong to the domain with differing block number, we still 
     !relate these points to the block index value= nGlobalTreeNode. 
     !This can be of sense only if the complete data sets at the     
     !extended point set is asessible at the PE, to which the block  
     !nGlobalTreeNode is assigned.                                   
     !Typically the data sets in the extended grid are shaped as     
     !arrays, like as follows:                                       
     !real,dimension(iGridPointMin_D(1):iGridPointMax_D(1),&         
     !               iGridPointMin_D(2):iGridPointMax_D(2),&         
     !         1    ...,&                                            
     !     min_block(DD%Ptr):max_block(DD%Ptr))::DataSet_GB          
     !                                                               
     !Now introduce the grids which are not cell centered ones.      
     !                                                               
     !All grids described by the present grid descriptor:= &         
     !                           extended grid.or.&                  
     !                           displaced[extended grid]            
     !where the displacement is governed by the vector as follows:   
     real,dimension(:),pointer::Displacement_D
     !For an arbitrary grid, for a given block the actual geometric   
     !displacement of the grid block fragment is equal to the product
     !of Displacement_D by the grid size
  end type GridDescriptorType

  !which can be introduced for some standard grids, like          
  integer,parameter, public ::CellCentered_=1,Nodes_=2
  public :: set_standard_grid_descriptor
  interface set_standard_grid_descriptor
     module procedure set_standard_grid_descriptor_id
     module procedure set_standard_grid_descriptor_dd
  end interface set_standard_grid_descriptor
  !\
  ! For a given global grid point number (in in4 or int8 format)
  ! which enumerates all points in the described grid and for a 
  ! given global grid descriptor, the routine returns
  ! a global node number and nDim grid point numbers in a block
  !/
  !\
  ! With a local grid descriptor, for a given local grid point
  ! number in a format int4, which enumerates all grid points in 
  ! blocks allocated at the current PE, the procedure returns the 
  ! order number of used block, which enumerates all used blocks at 
  ! the given PE in the same order the all blocks in the whole domain 
  ! decomposition are enumerated. The nDim grid point indexes are 
  ! returned too. For a global grid point number in int8 format,
  ! it returns the same output.  
  !/
  public :: global_i_grid_point_to_icb
  interface global_i_grid_point_to_icb
     module procedure global_i_grid_point_to_icb4
     module procedure global_i_grid_point_to_icb8
     module procedure global_i_grid_point_to_icb4_l
     module procedure global_i_grid_point_to_icb8_l
  end interface global_i_grid_point_to_icb
  !\
  ! Returns coordinate of the grid point, in the global or
  ! local grid descriptor, for given global tree node or used
  ! block number accordingly and for given point indexes in the 
  ! block
  !/
  public :: xyz_grid_d
  interface xyz_grid_d
     module procedure xyz_grid_d_global
     module procedure xyz_grid_d_local
  end interface xyz_grid_d
  
  public :: i_grid_point_global
  interface i_grid_point_global
     module procedure i_grid_point_global_g
     module procedure i_grid_point_global_l
  end interface i_grid_point_global
  public :: i8_grid_point_global
  interface i8_grid_point_global
     module procedure i8_grid_point_global_g
     module procedure i8_grid_point_global_l
  end interface i8_grid_point_global
  
  !local storage for a grid descriptor passed to interpolation_amr_gc;
  !it is shared by interpolate_amr_gc and find_amr
  type(GridDescriptorType) :: GDAMR
  !===========================Local Grid Descriptor===================!
  !The local GD describes in a compact way the parts of grid allocated 
  !at a given PE
  !\
  ! Order of indexes in the array
  ! Declared in CON_domain_decomposition:
  ! PE_ = 3, BLK_=2, GlobalBlock_ =4
  !/
  public :: BLK_, GlobalBlock_
  integer, parameter, public:: GlobalTreeNode_ = 1, GridPointFirst_ = 3 
  type LocalGDType
     !\
     ! These four members are just copied from the global GD
     !/
     integer           :: nDim           !Defines the grid dimensionality   
     integer, pointer  :: iGridPointMin_D(:)
     integer, pointer  :: iGridPointMax_D(:)
     real,    pointer  :: Displacement_D(:)
     !\
     ! The following terms are new and characterize the part of
     ! the global grid allocated at the given process
     !/
     integer:: nBlock, nPointPerBlock
     integer, pointer ::  iIndex_IB(:,:)
     real,    pointer ::  XyzBlock_DB(:,:), DXyz_DB(:,:)
  end type LocalGDType
  public :: set_local_gd, LocalGDType
  public :: n_grid_points_per_block
  !\
  ! Interpolation procedures
  !/
  public :: nearest_grid_points     !First order interpolation
  public :: bilinear_interpolation !For uniform or node-based grid
  public :: interpolation_amr_gc   !Employs ghostcells 
contains
  !ROUTINE: xyz_grid_d - the coordintes of the grid point
  !INTERFACE:
  function xyz_grid_d_global(GD,&
       lGlobalTreeNode,iGridPoints_D)
    !INPUT ARGUMENTS:     
    type(GridDescriptorType),intent(in)::GD                 
    integer,intent(in)::lGlobalTreeNode
    integer,dimension(GD%nDim),intent(in)::&
         iGridPoints_D
    !OUTPUT ARGUMENTS:
    real,dimension(GD%nDim)::&
         xyz_grid_d_global
    !EOP
    xyz_grid_d_global = xyz_block_d(GD%DD%Ptr,&
         lGlobalTreeNode)+& 
         d_xyz_cell_d(GD%DD%Ptr,&
         lGlobalTreeNode)*&
         (GD%Displacement_D - cHalf+real(iGridPoints_D))
  end function xyz_grid_d_global
  !===============================
  function xyz_grid_d_local(LocalGD,&
       iBlockUsed, iGridPoints_D)
    !INPUT ARGUMENTS:     
    type(LocalGDType), intent(in) :: LocalGD                 
    integer,           intent(in) :: iBlockUsed
    integer,           intent(in) :: iGridPoints_D(LocalGD%nDim)
    !OUTPUT ARGUMENTS:
    real                          :: xyz_grid_d_local(LocalGD%nDim)
    integer :: nDim
    !EOP
    nDim = LocalGD%nDim
    xyz_grid_d_local = LocalGD%XyzBlock_DB(1:nDim,iBlockUsed) +&
         LocalGD%DXyz_DB(1:nDim,iBlockUsed)*&
         (LocalGD%Displacement_D - 0.50 +real(iGridPoints_D))
  end function xyz_grid_d_local
  !===============================================================!
  !BOP
  !IROUTINE: set_standard_grid_descriptor - cell centered or node grids.
  !DESCRIPTION:
  !Allow to set the standard grid descriptor: 
  !iStandard\_=CellCentered\_ or iStandard\_=Nodes\_,
  !with or without halo points ("ghost points")
  !EOP                  
  !---------------------------------------------------------------!
  !BOP
  !INTERFACE:
  subroutine set_standard_grid_descriptor_id(&
       iGridID,&
       nGhostGridPoints,&
       iStandard,&
       GD)
    !INPUT ARGUMENTS:
    integer,intent(in)::iGridID
    integer,intent(in),optional::iStandard
    integer,intent(in),optional::nGhostGridPoints
    !OUTPUT ARGUMENTS:
    type(GridDescriptorType),intent(out)::GD
    !EOP
    integer::iError,iMyStandard,nGhostGridPointsMy
    call associate_dd_pointer(&
         iGridID,&
         GD%DD)
    GD%nDim=ndim_grid(iGridID)
    allocate(GD%Displacement_D(&
         1:GD%nDim),stat=iError)
    call check_allocate(iError,"Displacement_D")
    allocate(GD%iGridPointMin_D(&
         1:GD%nDim),stat=iError)
    call check_allocate(iError,"iGridPointMins_D")
    allocate(GD%iGridPointMax_D(&
         1:GD%nDim),stat=iError)
    call check_allocate(iError,"iGridPointMaxs_D")
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
       GD%Displacement_D= 0.0
       GD%iGridPointMin_D=&
            1-nGhostGridPointsMy
       GD%iGridPointMax_D=&
            ncells_decomposition_D(&
            iGridID)+nGhostGridPointsMy
    case(Nodes_)
       GD%Displacement_D=-cHalf
       GD%iGridPointMin_D=&
            min(1,2-nGhostGridPointsMy)
       GD%iGridPointMax_D=&
            ncells_decomposition_D(&
            iGridID)+nGhostGridPointsMy
    case default
       call CON_stop('Unknown standard for Grid Descriptor')
    end select
  end subroutine set_standard_grid_descriptor_id
  !--------------------------------------------------------!
  !BOP
  !INTERFACE:
  subroutine set_standard_grid_descriptor_dd(&
       DomainDecomposition,&
       nGhostGridPoints,&
       iStandard,&
       GD)
    !INPUT ARGUMENTS:
    type(DomainDecompositionType),target,intent(in)::&
         DomainDecomposition
    integer,intent(in),optional::iStandard
    integer,intent(in),optional::nGhostGridPoints
    !OUTPUT ARGUMENTS:
    type(GridDescriptorType),intent(out)::GD
    !EOP
    integer::iError,iMyStandard,nGhostGridPointsMy
    nullify(GD%DD%Ptr)
    call associate_dd_pointer(&
         DomainDecomposition,&
         GD%DD)
    GD%nDim=ndim_grid(&
         DomainDecomposition)
    allocate(GD%Displacement_D(&
         1:GD%nDim),stat=iError)
    call check_allocate(iError,"Displacement_D")
    allocate(GD%iGridPointMin_D(&
         1:GD%nDim),stat=iError)
    call check_allocate(iError,"iGridPointMins_D")
    allocate(GD%iGridPointMax_D(&
         1:GD%nDim),stat=iError)
    call check_allocate(iError,"iGridPointMaxs_D")
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
       GD%Displacement_D= 0.0
       GD%iGridPointMin_D=1-nGhostGridPointsMy
       GD%iGridPointMax_D=&
            ncells_decomposition_d(&
            DomainDecomposition)+nGhostGridPointsMy
    case(Nodes_)
       GD%Displacement_D=-cHalf
       GD%iGridPointMin_D=&
            min(1,2-nGhostGridPointsMy)
       GD%iGridPointMax_D=&
            ncells_decomposition_d(&
            DomainDecomposition)+nGhostGridPointsMy
    case default
       call CON_stop('Unknown standard for Grid Descriptor')
    end select
  end subroutine set_standard_grid_descriptor_dd
  !===============================================================!
  !More general grid descriptor                                   !
  !---------------------------------------------------------------!
  !BOP
  !IROUTINE: set_grid_descriptor - set more general grid
  !INTERFACE:
  subroutine set_grid_descriptor_id(&
       iGridID,&
       nDim,&
       iGridPointMin_D,&
       iGridPointMax_D,&
       Displacement_D,&
       GD)
    !INPUT ARGUMENTS:
    integer,intent(in)::iGridID
    integer,intent(in)::nDim
    integer,intent(in),dimension(nDim)::iGridPointMin_D
    integer,intent(in),dimension(nDim)::iGridPointMax_D
    real,intent(in),dimension(nDim)::Displacement_D
    !OUTPUT ARGUMENTS:
    type(GridDescriptorType),intent(out)::GD
    !EOP
    integer::iError,iMyStandard,nGhostGridPointsMy
    call associate_dd_pointer(&
         iGridID,&
         GD%DD)
    GD%nDim=nDim
    allocate(GD%Displacement_D(&
         1:GD%nDim),stat=iError)
    call check_allocate(iError,"Displacement_D")
    GD%Displacement_D=Displacement_D
    allocate(GD%iGridPointMin_D(&
         1:GD%nDim),stat=iError)
    call check_allocate(iError,"iGridPointMins_D")
    GD%iGridPointMin_D=iGridPointMin_D
    allocate(GD%iGridPointMax_D(&
         1:GD%nDim),stat=iError)
    call check_allocate(iError,"iGridPointMaxs_D")
    GD%iGridPointMax_D=iGridPointMax_D
  end subroutine set_grid_descriptor_id
  !---------------------------------------------------------------!
  !BOP
  !INTERFACE:        
  subroutine set_grid_descriptor_dd(&
       DomainDecomposition,&
       nDim,&
       iGridPointMin_D,&
       iGridPointMax_D,&
       Displacement_D,&
       GD)
    !INPUT ARGUMENTS:
    type(DomainDecompositionType),target,intent(in)::&
         DomainDecomposition
    integer,intent(in)::nDim
    integer,intent(in),dimension(nDim)::iGridPointMin_D
    integer,intent(in),dimension(nDim)::iGridPointMax_D
    real,intent(in),dimension(nDim)::Displacement_D
    !OUTPUT ARGUMENTS:
    type(GridDescriptorType),intent(out)::GD
    !EOP
    integer::iError,iMyStandard,nGhostGridPointsMy
    call associate_dd_pointer(&
         DomainDecomposition,&
         GD%DD)
    GD%nDim=nDim
    allocate(GD%Displacement_D(&
         1:GD%nDim),stat=iError)
    call check_allocate(iError,"Displacement_D")
    GD%Displacement_D=Displacement_D
    allocate(GD%iGridPointMin_D(&
         1:GD%nDim),stat=iError)
    call check_allocate(iError,"iGridPointMins_D")
    GD%iGridPointMin_D=iGridPointMin_D
    allocate(GD%iGridPointMax_D(&
         1:GD%nDim),stat=iError)
    call check_allocate(iError,"iGridPointMaxs_D")
    GD%iGridPointMax_D=iGridPointMax_D
  end subroutine set_grid_descriptor_dd
  !===============================================================!
  subroutine set_local_gd(iProc, GD, LocalGD)
    !\
    ! PE rank (in a local group for a local model, or in a 
    ! global group for a global model)
    !/
    integer,                 intent(in) :: iProc 
    type(GridDescriptorType),intent(in) :: GD
    type(LocalGDType),       intent(out):: LocalGD
    !Misc
    integer :: nDim, nBlock, nPointPerBlock
    integer :: iBlock, iError
    integer, pointer:: iDD_II(:,:), iGlobal_A(:)
    !----------------------------------
    nDim                     = GD%nDim
    LocalGD%nDim             = nDim
    LocalGD%Displacement_D  => GD%Displacement_D
    LocalGD%iGridPointMin_D => GD%iGridPointMin_D
    LocalGD%iGridPointMax_D => GD%iGridPointMax_D
    nBlock = n_block(GD%DD%Ptr, iProc)
    LocalGD%nBlock           = nBlock
    nPointPerBlock =    n_grid_points_per_block(GD)
    LocalGD%nPointPerBlock   = nPointPerBlock

    allocate(LocalGD%iIndex_IB(GlobalTreeNode_:GlobalBlock_,&
         1:nBlock),stat=iError)
    call check_allocate(iError,'LocalGD%iIndex_IB')
    !\
    ! For better readability
    !/
    iDD_II=>GD%DD%Ptr%iDecomposition_II
    iGlobal_A=>GD%DD%Ptr%iGlobal_A
    LocalGD%iIndex_IB(GlobalTreeNode_,:) = pack(iGlobal_A,&
         MASK = iDD_II(PE_,iGlobal_A)==iProc)
    LocalGD%iIndex_IB(BLK_,:) = iDD_II(BLK_,&
         LocalGD%iIndex_IB(GlobalTreeNode_,:))
    LocalGD%iIndex_IB(GlobalBlock_,:) = iDD_II(GlobalBlock_,&
         LocalGD%iIndex_IB(GlobalTreeNode_,:))
    allocate(LocalGD%DXyz_DB(1:nDim,1:nBlock),stat=iError)
    call check_allocate(iError,'LocalGD%DXyz_DB')
    allocate(LocalGD%XyzBlock_DB(1:nDim,1:nBlock),stat=iError)
    call check_allocate(iError,'LocalGD%XyzBlock_DB')
    do iBlock = 1, nBlock
       LocalGD%iIndex_IB(GridPointFirst_,iBlock) = (iBlock - 1)*&
            nPointPerBlock + 1
       LocalGD%DXyz_DB(1:nDim,iBlock)     = &
            GD%DD%Ptr%DXyzCell_DI(1:nDim, &
            LocalGD%iIndex_IB(GlobalTreeNode_,iBlock))
       LocalGD%XyzBlock_DB(1:nDim,iBlock) = &
            GD%DD%Ptr%XyzBlock_DI(1:nDim, &
            LocalGD%iIndex_IB(GlobalTreeNode_,iBlock))
    end do
  end subroutine set_local_gd
  !===============================================================!
  !BOP
  !IROUTINE: clean_grid_descriptor
  !INTERFACE:
  subroutine clean_grid_descriptor(GD)
    !INPUT ARGUMENTS:
    type(GridDescriptorType),intent(inout)::&
         GD
    !EOP
    deallocate(GD%iGridPointMin_D)
    deallocate(GD%iGridPointMax_D)
    deallocate(GD%Displacement_D)
    nullify(GD%DD%Ptr)
  end subroutine clean_grid_descriptor
  !===============================================================!
  subroutine clean_grid_descriptor_l(GD)
    !INPUT ARGUMENTS:
    type(LocalGDType),intent(inout)::&
         GD
    !EOP
    deallocate(GD%iGridPointMin_D)
    deallocate(GD%iGridPointMax_D)
    deallocate(GD%Displacement_D)
    deallocate(GD%iIndex_IB)
    deallocate(GD%XyzBlock_DB) 
    deallocate(GD%DXyz_DB)
  end subroutine clean_grid_descriptor_l
  !===============================================================!
  !Number of cells per block
  integer function n_grid_points_per_block(GD)
    type(GriddescriptorType)::GD
    n_grid_points_per_block = product(&
         GD%iGridPointMax_D + 1 -GD%iGridPointMin_D)
  end function n_grid_points_per_block
  !===============================================================!
  !          ENUMERATION                                          !
  !All the points of the grid of the type considered here can be  !
  !easily enumerated, but with the simplest way to do this there  !
  !is some amounts of unused global cell numbers, for             !
  !TreeDecompositions. The following procedures are used the      !
  !CON_router and will be used then for coupling via the MCT, to  !
  !construct the global segmentation map                          !
  !===============================================================!
  !The following procedure transforms the global grid point number!
  ! to the iCB index                                              !
  !To use this procedure, the grid descriptor should be           !
  !constructed  first                                             !
  subroutine global_i_grid_point_to_icb4(&
       GD,&
       iGridPointGlobal,&
       lGlobalTreeNode,&
       iGridPoint_D)
    type(GridDescriptorType),intent(in)::GD
    integer,intent(in)::iGridPointGlobal
    integer,intent(out):: lGlobalTreeNode
    integer,dimension(GD%nDim),intent(out)::&
         iGridPoint_D
    integer,dimension(GD%nDim)::nGridPoints_D 
    integer::iDim,iMisc,iMisc1         
    nGridPoints_D=1+GD%iGridPointMax_D-&
         GD%iGridPointMin_D
    iMisc=iGridPointGlobal-1
    do idim=1,GD%nDim
       iMisc1=mod(iMisc,nGridPoints_D(iDim))
       iGridPoint_D(iDim)=GD%iGridPointMin_D(iDim)+iMisc1
       iMisc=(iMisc-iMisc1)/nGridPoints_D(iDim)
    end do
    lGlobalTreeNode=i_global_node_a(GD%DD%Ptr,iMisc+1)
  end subroutine global_i_grid_point_to_icb4
  !===============================================================
  subroutine global_i_grid_point_to_icb8(&
       GD,&
       iGridPointGlobal,&
       lGlobalTreeNode,&
       iGridPoint_D)
    type(GridDescriptorType),intent(in)::GD
    integer(Int8_),intent(in)::iGridPointGlobal
    integer,intent(out):: lGlobalTreeNode
    integer,dimension(GD%nDim),intent(out)::&
         iGridPoint_D
    integer(Int8_),dimension(GD%nDim)::nGridPoints_D 
    integer::iDim, iMisc1
    integer(Int8_):: i8Misc
    nGridPoints_D = 1_Int8_ + GD%iGridPointMax_D -&
         GD%iGridPointMin_D
    i8Misc = iGridPointGlobal - 1_Int8_
    do idim=1,GD%nDim
       iMisc1=int(mod(i8Misc, nGridPoints_D(iDim)))
       iGridPoint_D(iDim) = GD%iGridPointMin_D(iDim) + iMisc1
       i8Misc = (i8Misc-iMisc1)/nGridPoints_D(iDim)
    end do
    lGlobalTreeNode = i_global_node_a(GD%DD%Ptr,int(i8Misc)+1)
  end subroutine global_i_grid_point_to_icb8
  !=========================================
  subroutine global_i_grid_point_to_icb4_l(&
       GD, &
       iGridPointLocal,&
       iBlockUsed,     &
       iGridPoint_D)
    type(LocalGDType),    intent(in) :: GD
    integer,              intent(in) :: iGridPointLocal
    integer,             intent(out) :: iBlockUsed
    integer,             intent(out) :: &
         iGridPoint_D(GD%nDim)
    !\
    !Local variables
    !/
    !\
    ! Number of grid points in the block, for each direction
    !/
    integer  :: nGridPoints_D(GD%nDim) 
    integer  :: iDim,iMisc,iMisc1
    !-----------------------------         
    nGridPoints_D=1+GD%iGridPointMax_D-&
         GD%iGridPointMin_D

    iMisc = iGridPointLocal - 1
    do idim = 1, GD%nDim
       iMisc1 = mod(iMisc,nGridPoints_D(iDim))
       iGridPoint_D(iDim) = GD%iGridPointMin_D(iDim)+iMisc1
       iMisc=(iMisc-iMisc1)/nGridPoints_D(iDim)
    end do
    iBlockUsed = iMisc + 1
  end subroutine global_i_grid_point_to_icb4_l
  !===============================================================
  subroutine global_i_grid_point_to_icb8_l(&
       GD,  &
       iGridPointGlobal,&
       iBlockUsed,      &
       iGridPoint_D)
    type(LocalGDType), intent(in) :: GD
    integer(Int8_),    intent(in) :: iGridPointGlobal
    integer,          intent(out) :: iBlockUsed
    integer,          intent(out) :: &
         iGridPoint_D(GD%nDim)
    !\
    !Local variables
    !/
    !\
    ! Number of grid points in the block, for each direction
    !/
    integer(Int8_)   :: nGridPoints_D(GD%nDim) 
    integer          :: iDim, iMisc1, loc(1)
    integer(Int8_)   :: i8Misc
    !-------------------------
    nGridPoints_D = 1_Int8_ + GD%iGridPointMax_D -&
         GD%iGridPointMin_D
    i8Misc =  iGridPointGlobal - 1
    do iDim = 1, GD%nDim
       iMisc1 = int(mod(i8Misc, nGridPoints_D(iDim)))
       iGridPoint_D(iDim) = GD%iGridPointMin_D(iDim)+iMisc1
       i8Misc = (i8Misc - iMisc1)/nGridPoints_D(iDim)
    end do
    loc = maxloc(GD%iIndex_IB(GlobalBlock_,:),&
         MASK = GD%iIndex_IB(GlobalBlock_,:).le.&
         (int(i8Misc)+1))
    iBlockUsed = loc(1)
  end subroutine global_i_grid_point_to_icb8_l
  !===============================================================!
  !The inverse procedures                                         !

  integer function i_grid_point_global_g(&
       GD,&
       lGlobalTreeNode,&
       iGridPoint_D)
    type(GridDescriptorType), intent(in) :: GD
    integer,                  intent(in) :: lGlobalTreeNode
    integer,                  intent(in) :: &
         iGridPoint_D(GD%nDim)
    integer :: nGridPoints_D(GD%nDim)
    integer :: iDim, iMisc, iMisc1
    !---------------------------
    nGridPoints_D(1:GD%nDim) = 1 +  &
         GD%iGridPointMax_D -       &
         GD%iGridPointMin_D
    i_grid_point_global_g = i_global_block(GD%DD%Ptr,&
         lGlobalTreeNode) - 1
    do idim = GD%nDim, 1, -1
       i_grid_point_global_g = i_grid_point_global_g*     &
            nGridPoints_D(iDim) + iGridPoint_D(iDim) -    &
            GD%iGridPointMin_D(iDim)
    end do
    i_grid_point_global_g = i_grid_point_global_g + 1
  end function i_grid_point_global_g
  !====================================
  integer(Int8_) function i8_grid_point_global_g(&
       GD,&
       lGlobalTreeNode,&
       iGridPoint_D)
    type(GridDescriptorType), intent(in) :: GD
    integer,                  intent(in) :: lGlobalTreeNode
    integer,                  intent(in) :: &
         iGridPoint_D(GD%nDim)
    integer :: nGridPoints_D(GD%nDim)
    integer :: iDim, iMisc, iMisc1
    nGridPoints_D(1:GD%nDim) = 1 + &
         GD%iGridPointMax_D - &
         GD%iGridPointMin_D
    i8_grid_point_global_g = i_global_block(GD%DD%Ptr,&
         lGlobalTreeNode) - 1
    do idim = GD%nDim, 1, -1
       i8_grid_point_global_g = i8_grid_point_global_g*&
            nGridPoints_D(iDim) + iGridPoint_D(iDim) -&
            GD%iGridPointMin_D(iDim)
    end do
    i8_grid_point_global_g = i8_grid_point_global_g + 1
  end function i8_grid_point_global_g
  !==================================
  integer function i_grid_point_global_l(&
       GD,&
       iBlockUsed,    &
       iGridPoint_D)
    type(LocalGDType), intent(in) :: GD
    integer,           intent(in) :: iBlockUsed
    integer,           intent(in) :: &
         iGridPoint_D(GD%nDim)
    integer :: nGridPoints_D(GD%nDim)
    integer :: iDim, iMisc, iMisc1
    !---------------------------
    nGridPoints_D(1:GD%nDim) = 1 +  &
         GD%iGridPointMax_D -       &
         GD%iGridPointMin_D
    i_grid_point_global_l = iBlockUsed - 1
    do idim = GD%nDim, 1, -1
       i_grid_point_global_l = i_grid_point_global_l*     &
            nGridPoints_D(iDim) + iGridPoint_D(iDim) -    &
            GD%iGridPointMin_D(iDim)
    end do
    i_grid_point_global_l = i_grid_point_global_l + 1
  end function i_grid_point_global_l
  !====================================
  integer(Int8_) function i8_grid_point_global_l(&
       GD,&
       iBlockUsed,    &
       iGridPoint_D)
    type(LocalGDType), intent(in) :: GD
    integer,           intent(in) :: iBlockUsed
    integer,           intent(in) :: &
         iGridPoint_D(GD%nDim)
    integer :: nGridPoints_D(GD%nDim)
    integer :: iDim, iMisc, iMisc1
    !----------------
    nGridPoints_D(1:GD%nDim) = 1 + &
         GD%iGridPointMax_D - &
         GD%iGridPointMin_D
    i8_grid_point_global_l = GD%iIndex_IB(&
         GlobalBlock_,iBlockUsed) - 1
    do idim = GD%nDim, 1, -1
       i8_grid_point_global_l = i8_grid_point_global_l*&
            nGridPoints_D(iDim) + iGridPoint_D(iDim) -&
            GD%iGridPointMin_D(iDim)
    end do
    i8_grid_point_global_l = i8_grid_point_global_l + 1
  end function i8_grid_point_global_l
  !==================================
  !BOP
  !IROUTINE: nearest_grid_points - first order interpolation
  !EOP
  !===============================================================!
  !======================INTERPOLATION============================!
  !BOP
  !DESCRIPTION:
  !The general principle of the interpoaltions with the grids     
  !defined by the present grid descriptor is as follows.          
  !Using coordinates of the point at which the interpalated value 
  !should be found, first of all, the block is found this point   
  !belongs to. Then, among the grid points belonging to the grid  
  !fragment with the found value of the block index, the nearest  
  !point is found. For the second order of interpolation, the     
  !cubic stenchil is constructed involving the nearest point too  
  !and the interpolation weights are found. If the values of the  
  !cell index for the stenchil points appear to be out of the     
  !block, the correspondent blocks are found                      
  !EOP 
  !BOP
  !INTERFACE:
  subroutine nearest_grid_points(nDim,Xyz_D,&
       GD,&
       nIndex,&
       iIndex_II,&
       nImage, Weight_I)
    !INPUT ARGUMENTS:                       
    type(GridDescriptorType):: GD
    integer,      intent(in):: nDim
    real,      intent(inout):: Xyz_D(nDim) 
    integer,      intent(in):: nIndex

    !OUTPUT ARGUMENTS:
    integer,     intent(out):: iIndex_II(0:nIndex,2**nDim)
    integer,     intent(out):: nImage
    real,        intent(out):: Weight_I(2**nDim)
    !EOP

    real,    dimension(nDim)::&
         XyzStored_D, DxyzCells_D,DxyzTolerance_D

    integer, dimension(nDim):: iGridPoints_D
    logical, dimension(nDim):: IsAtFace_D
    real                    :: Xyz_DI(nDim,2**nDim)
    real,          parameter:: Tolerance = 0.001
    integer::lGlobalTreeNode,iDim,iImages
    !------------------
    !\
    !Initialize arrays
    !/
    iIndex_II=0;Weight_I= 0.0
    XyzStored_D=Xyz_D
    iIndex_II(1:nDim,1)=&
         GD%iGridPointMax_D+1
    do while(any(iIndex_II(1:nDim,1)>&
         GD%iGridPointMax_D))
       Xyz_D=XyzStored_D
       !\
       ! Find a block to which the point belongs
       !/
       call search_in(GD%DD%Ptr,Xyz_D,lGlobalTreeNode)
       DxyzCells_D=d_xyz_cell_d(&
            GD%DD%Ptr,lGlobalTreeNode)
       !\
       ! Check if the point is out of the computational domain
       !/
       where(is_right_boundary_d(GD%DD%Ptr,lGlobalTreeNode))
          Xyz_D=min(Xyz_D,cAlmostOne*d_xyz_block_d(&         
               GD%DD%Ptr,lGlobalTreeNode))
          XyzStored_D=Xyz_D+xyz_block_d(&         
               GD%DD%Ptr,lGlobalTreeNode)
       end where
       where(is_left_boundary_d(GD%DD%Ptr,lGlobalTreeNode))
          Xyz_D=max(Xyz_D, 0.0)
          XyzStored_D=Xyz_D+xyz_block_d(&         
               GD%DD%Ptr,lGlobalTreeNode)
       end where
       !\
       ! Store the processor number and block number. If the block
       ! number is meaningless use nIndex=nDim and this meaningless
       ! value will be overwritten later
       !/
       call pe_and_blk(GD%DD%Ptr,lGlobalTreeNode,&
            iIndex_II(0,1),iIndex_II(nIndex,1))
       !\
       ! Find the nearest grid point in the block
       !/
       Xyz_D=Xyz_D-DxyzCells_D*GD%Displacement_D
       iGridPoints_D=floor(Xyz_D/DxyzCells_D)
       iIndex_II(1:nDim,1)=iGridPoints_D+1
       DxyzTolerance_D=Tolerance*DxyzCells_D
       !\
       ! The vector Xyz - XyzGrid + 0.5 Dxyz
       !/
       Xyz_D=Xyz_D-DxyzCells_D*iGridPoints_D
       !\
       ! If the point coordintates are such that Xyz-XyzGrid=+/-0.5*Dxyz
       ! there are several equidistant grid points
       !/

       IsAtFace_D=abs(Xyz_D)<DxyzTolerance_D&
            .or.abs(Xyz_D-DxyzCells_D )<DxyzTolerance_D
       where(is_right_boundary_d(GD%DD%Ptr,lGlobalTreeNode))&
            iIndex_II(1:nDim,1)=&
            min(iIndex_II(1:nDim,1),&
            GD%iGridPointMax_D)
       !\
       !The nearest grid points may be are in the upper block
       !/
       where(iIndex_II(1:nDim,1)>&
            GD%iGridPointMax_D)&
            XyzStored_D=XyzStored_D+DxyzTolerance_D*0.25
    end do
    nImage=1;Weight_I(1)= 1.0
    if(.not.(any(IsAtFace_D)))return
    Xyz_DI(:,1)=XyzStored_D
    do iDim=1,nDim
       if(IsAtFace_D(iDim))then
          Xyz_DI(:,1+nImage:nImage+nImage)=Xyz_DI(:,1:nImage)
          Xyz_DI(iDim,1:nImage)=Xyz_DI(iDim,1:nImage)-&
               DxyzTolerance_D(iDim)
          Xyz_DI(iDim,1+nImage:nImage+nImage)=Xyz_DI(&
               iDim,1+nImage:nImage+nImage)+&
               DxyzTolerance_D(iDim)
          nImage=nImage+nImage
       end if
    end do
    do iImages=1,nImage
       call search_in(GD%DD%Ptr,&
            Xyz_DI(:,iImages),&
            lGlobalTreeNode)
       call pe_and_blk(GD%DD%Ptr,&
            lGlobalTreeNode,&
            iIndex_II(0,iImages),&
            iIndex_II(nIndex,iImages))
       Xyz_DI(:,iImages)=Xyz_DI(:,iImages)-&
            d_xyz_cell_d(GD%DD%Ptr,&
            lGlobalTreeNode)*GD%Displacement_D
       call search_cell(GD%DD%Ptr,&
            lGlobalTreeNode,Xyz_DI(:,iImages),&
            iIndex_II(1:nDim,iImages))
    end do
    iImages=1
    do while(iImages<=nImage)
       !Exclude the stencil nodes which are out of the 
       !computational domain
       do while(any(GD%iGridPointMin_D>&
            iIndex_II(1:nDim,iImages).or.&
            GD%iGridPointMax_D<&
            iIndex_II(1:nDim,iImages)))
          if(iImages==nImage)then
             nImage=nImage-1
             EXIT
          end if
          iIndex_II(:,iImages)=iIndex_II(:,nImage)
          nImage=nImage-1
       end do
       iImages=iImages+1
    end do
    Weight_I(1:nImage)= 1.0/real(nImage)
  end subroutine nearest_grid_points
  !===============================================================!
  !=================SECOND ORDER INTERPOLATION====================!
  !BOP
  !IROUTINE: bilinear_interpolation - second order interpolation
  !EOP
  !BOP
  !DESCRIPTION:
  !This is a bilinear interpolation using the grid points,        
  !described with the grid descriptor.                            
  !EOP
  !INTERFACE:
  subroutine bilinear_interpolation(&
       nDim,&
       Xyz_D,&
       GD,&
       nIndex,& 
       iIndex_II,&
       nImage,Weight_I)
    !INPUT ARGUMENTS:
    integer,intent(in)      ::nDim
    type(GridDescriptorType)::GD     
    real,intent(inout)::Xyz_D(nDim)
    integer,intent(in)::nIndex
    !OUTPUT ARGUMENTS
    integer, intent(out) :: iIndex_II(0:nIndex,2**nDim)

    integer,intent(out)::nImage
    real,dimension(2**nDim),intent(out)::Weight_I
    !EOP
    real,dimension(nDim)::&
         XyzResid_D,XyzStored_D
    integer,dimension(nDim)::iGridPoints_D
    real,dimension(nDim):: XyzMisc_D
    integer,dimension(2**nDim)::lNode_I
    integer::lGlobalTreeNode,iDim,iImages,iNewStart,nImageNew
    real::WeightLeft
    logical,dimension(nDim,2**nDim)::&
         IsUp_DI,IsDown_DI
    real,dimension(nDim,0:2**nDim)::&
         Dxyz_DI
    logical,dimension(nDim)::&
         IsDomainBoundaryUp_D,IsDomainBoundaryDown_D

    iIndex_II=0
    Weight_I= 0.0
    XyzStored_D=Xyz_D
    IsUp_DI(:,1)=.true.
    call search_in(GD%DD%Ptr,Xyz_D,lGlobalTreeNode)
    !Find global node number, PE and number which involved the point
    do while(any(IsUp_DI(:,1)))
       !This do loop works more than once only in case of node
       !grid and only in case of Xyz_D point belonging to a block
       !boundary. At this case the routine needs help in deciding,
       !to which block this point should be assigned. Otherwise,
       !automatically IsUp_D(:,1)=.false. after first loop pass.
       Xyz_D=XyzStored_D-xyz_block_d(&
            GD%DD%Ptr,lGlobalTreeNode)
       IsDomainBoundaryUp_D=&
            is_right_boundary_d(GD%DD%Ptr,lGlobalTreeNode)
       IsDomainBoundaryDown_D=&
            is_left_boundary_d(GD%DD%Ptr,lGlobalTreeNode)
       where(IsDomainBoundaryUp_D)
          Xyz_D=min(Xyz_D,cAlmostOne*d_xyz_block_d(&         
               GD%DD%Ptr,lGlobalTreeNode))
          XyzStored_D=Xyz_D+xyz_block_d(&         
               GD%DD%Ptr,lGlobalTreeNode)
       end where
       where(IsDomainBoundaryDown_D)
          Xyz_D=max(Xyz_D, 0.0)
          XyzStored_D=Xyz_D+xyz_block_d(&         
               GD%DD%Ptr,lGlobalTreeNode)
       end where

       call pe_and_blk(GD%DD%Ptr,lGlobalTreeNode,&
            iIndex_II(0,1),iIndex_II(nIndex,1))

       !Find Dxyz for this block
       Dxyz_DI(:,0)=d_xyz_cell_d(&         
            GD%DD%Ptr,lGlobalTreeNode)

       !\
       !/

       XyzResid_D=Xyz_D/Dxyz_DI(:,0)-&
            GD%Displacement_D+cHalf
       iGridPoints_D=floor(XyzResid_D)
       XyzResid_D=XyzResid_D-real(iGridPoints_D)
       !Thus calculated XyzResid_D satisfies the inequalities as       !
       !follow:XyzResid_D>=0 and XyzResid_D<1. It is used to calculte  !
       !weights for the eight grid points, among them the iGridPoints_D!
       !being the left one with respect to all the spatial directions  !
       iIndex_II(1:nDim,1)=iGridPoints_D
       IsUp_DI(:,1)=iIndex_II(1:nDim,1)>&
            GD%iGridPointMax_D
       IsDown_DI(:,1)=iIndex_II(1:nDim,1)<&
            GD%iGridPointMin_D
       if(any(IsUp_DI(:,1)))&
            lGlobalTreeNode=l_neighbor(&
            GD%DD%Ptr,lGlobalTreeNode,&
            iIndex_II(1:nDim,1))
    end do

    nImage=1;Weight_I(1)= 1.0

    do iDim=1,nDim
       !Exclude the stencil nodes which are out of the 
       !computational domain
       if(IsDown_DI(iDim,1).and.IsDomainBoundaryDown_D(iDim))then
          IsDown_DI(iDim,1:nImage)=.false.
          iIndex_II(iDim,1:nImage)=iIndex_II(iDim,1:nImage)+1
          CYCLE
       end if

       if(iIndex_II(iDim,1)==GD%iGridPointMax_D(iDim)&
            .and.IsDomainBoundaryUp_D(iDim))CYCLE
       if(XyzResid_D(iDim)<cTiny)CYCLE
       iNewStart=nImage+1;nImageNew=nImage+nImage
       iIndex_II(:,iNewStart:nImageNew)=&
            iIndex_II(:,1:nImage)
       iIndex_II(iDim,iNewStart:nImageNew)=&
            iIndex_II(iDim,iNewStart:nImageNew)+1
       IsUp_DI(:,iNewStart:nImageNew)=IsUp_DI(:,1:nImage)
       IsUp_DI(iDim,iNewStart:nImageNew)=&
            iIndex_II(iDim,iNewStart)>&
            GD%iGridPointMax_D(iDim)
       IsDown_DI(:,iNewStart:nImageNew)=IsDown_DI(:,1:nImage)
       IsDown_DI(iDim,iNewStart:nImageNew)=.false.
       Weight_I(iNewStart:nImageNew)= &
            Weight_I(1:nImage)*XyzResid_D(iDim)
       WeightLeft= 1.0 -XyzResid_D(iDim)
       Weight_I(1:nImage)= WeightLeft* Weight_I(1:nImage)
       nImage=nImageNew
    end do
    !---------------------------------------------------------------!
    !Check if the grid point index is within the index limits       !
    do iImages=1,nImage
       if(.not.(any(IsUp_DI(:,iImages)).or.&
            any(IsDown_DI(:,iImages))))then
          Dxyz_DI(:,iImages)=Dxyz_DI(:,0)
          lNode_I(iImages)=lGlobalTreeNode
       else
          lNode_I(iImages)=l_neighbor(GD%DD%Ptr,&
               lGlobalTreeNode,iIndex_II(&
               1:nDim,iImages))
          XyzMisc_D=xyz_grid_d(GD,&
               lGlobalTreeNode,iIndex_II(&
               1:nDim,iImages))-&
               xyz_block_d(GD%DD%Ptr,lNode_I(iImages))
          call pe_and_blk(GD%DD%Ptr,lNode_I(iImages),&
               iIndex_II(0,iImages),iIndex_II(nIndex,iImages))
          Dxyz_DI(:,iImages)=d_xyz_cell_d(GD%DD%Ptr,&
               lNode_I(iImages))
          XyzMisc_D=XyzMisc_D- Dxyz_DI(:,iImages)&
               *GD%Displacement_D
          call search_cell(GD%DD%Ptr,&
               lNode_I(iImages),XyzMisc_D,&
               iIndex_II(1:nDim,iImages))
       end if
    end do
  end subroutine bilinear_interpolation
  !===============================================================!
  !BOP
  !IROUTINE: interpolation_amr_gc - continuous 2nd order interpolation on AMR
  !utilizing ghost cells
  !EOP
  !BOP
  !DESCRIPTION:
  !This is a continuous amr interpolation using the grid points,        
  !described with the grid descriptor. It utilizes ghost cells.
  !EOP
  !INTERFACE:
  subroutine interpolation_amr_gc(&
       nDim,&
       Xyz_D,&
       GD,&
       nIndex,& 
       iIndex_II,&
       nImage,Weight_I)
    use ModInterpolateAMR, ONLY: interpolate_amr
    !INPUT ARGUMENTS:
    integer,   intent(in)   :: nDim
    type(GridDescriptorType):: GD     
    real,      intent(inout):: Xyz_D(nDim)
    integer,   intent(in)   :: nIndex
    !OUTPUT ARGUMENTS
    integer,   intent(out)  :: iIndex_II(0:nIndex,2**nDim)
    integer,   intent(out)  :: nImage
    real,      intent(out)  :: Weight_I(2**nDim)
    !EOP
    !--------------------------------------------------------------------------
    ! memorize grid descriptor by storing it in modular variable;
    ! it is used by find_amr subroutine passed to shared interpolation routine
    GDAMR = GD
    ! call shared interpolation subroutine
    call interpolate_amr(&
         nDim        = nDim,&
         XyzIn_D     = Xyz_D,&
         nIndexes    = nIndex,&
         find        = find_amr,&
         nCell_D     = GDAMR%DD%Ptr%nCells_D,&
         nGridOut    = nImage,&
         Weight_I    = Weight_I,&
         iIndexes_II = iIndex_II,&
         UseGhostCell= .true.)
  end subroutine interpolation_amr_gc

  !===============================================================!

  subroutine find_amr(nDim, Xyz_D, &
       iProc, iBlock, XyzCorner_D, Dxyz_D, IsOut)
    use ModNumConst
    integer, intent(in)   :: nDim
    real,    intent(inout):: Xyz_D(nDim)
    integer, intent(out)  :: iProc, iBlock 
    real,    intent(out)  :: XyzCorner_D(nDim), Dxyz_D(nDim)
    logical, intent(out)  :: IsOut

    integer:: iNode
    logical:: IsPeriodic_D(nDim)
    !--------------------------------------------------------------------------
    ! check if point's inside the domain
    IsPeriodic_D = GDAMR%DD%Ptr%IsPeriodic_D 
    IsOut = &
         any(GDAMR%DD%Ptr%XyzMin_D >  Xyz_D&
         .and..not.IsPeriodic_D) .or. &
         any(GDAMR%DD%Ptr%XyzMax_D <= Xyz_D&
         .and..not.IsPeriodic_D)
    if(IsOut) then
       if(GDAMR%DD%Ptr%DoGlueMargins)then
          !Recalculate coordinates for spherical or similar
          !grids
          call glue_margin(GDAMR%DD%Ptr, Xyz_D)
          IsOut = &
               any(GDAMR%DD%Ptr%XyzMin_D >  Xyz_D&
               .and..not.IsPeriodic_D) .or. &
               any(GDAMR%DD%Ptr%XyzMax_D <= Xyz_D&
               .and..not.IsPeriodic_D)
          if(IsOut)RETURN
       else
          RETURN
       end if
    end if
    ! call subroutine to find a node containing Xyz_D
    ! this call changes Xyz_D to coords relative to block's corner as needed
    call search_in(GDAMR%DD%Ptr, Xyz_D, iNode)

    ! now extract block's parameters
    XyzCorner_D = GDAMR%DD%Ptr%XyzBlock_DI(:, iNode)
    Dxyz_D      = GDAMR%DD%Ptr%DxyzCell_DI(:, iNode)
    iProc       = GDAMR%DD%Ptr%iDecomposition_II(PE_, iNode) 
    iBlock      = GDAMR%DD%Ptr%iDecomposition_II(BLK_, iNode) 
  end subroutine find_amr

  !==============================END==============================!
end Module CON_grid_descriptor
!\end{verbatim}                     !^CFG UNCOMMENT IF PRINTPS  !
!\end{document}                     !^CFG UNCOMMENT IF PRINTPS  !
!=============================LINE 1257===========================!
!BOP
!IROUTINE: more methods
!DESCRIPTION:
!More methods are available, see the source file.
!EOP
