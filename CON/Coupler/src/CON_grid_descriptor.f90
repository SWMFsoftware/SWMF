!\documentstyle[12pt]{letter}      !^CFG UNCOMMENT IF PRINTPS   !
!\setlength{\textwidth}{6.5 in}    !^CFG UNCOMMENT IF PRINTPS   !
!\setlength{\oddsidemargin}{0.5 in}!^CFG UNCOMMENT IF PRINTPS   !
!\begin{document}                  !^CFG UNCOMMENT IF PRINTPS   !
!\begin{verbatim}                  !^CFG UNCOMMENT IF PRINTPS   !
!^CFG COPYRIGHT UM                                              !
!BOP
!MODULE: CON_grid_descriptor - for uniform or octree grids
!INTERFACE:
Module CON_grid_descriptor

  !DESCRIPTION:
  !The toolkit for coupling the different data sets within a      
  !single component code, or within a framework as it is now (see 
  !the date below) consists of: 
  !                                  
  !ModMpi (iclude for mpif90.h)
  !                                   
  !ModNumConst
  !                                                    
  !CON\_world
  !                                                                   
  !CON\_domain\_decomposition
  !                                     
  !CON\_grid\_descriptor
  !                                          
  !CON\_router
  !                                                    
  !CON\_global\_message\_pass                                     


  !This file presents the class of the grid descriptors which     
  !includes both the uniformly spaced grids (uniformly    
  !spaced with respect to some generalized coordinates) and Octree
  ! or Quadric tree for adaptive block grid.                      
  !                                                               
  !The methods include the grid descriptor allocation, coordinate  
  !computations and the interpolation procedures                  

  !USES:

  use CON_grid_storage
  use ModUtilities, ONLY: check_allocate
  !REVISION HISTORY:
  ! Sokolov I.V.                                                  
  ! 6.18.03-7.11.03                                               
  ! igorsok@umich.edu                                             
  ! phone(734)647-4705                                            
  !EOP
  !===============================================================

  implicit none    
  !BOP
  !DESCRIPTION:
  !\begin{verbatim}                                
  !=====================DERIVED TYPE============================== 
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
  end type GridDescriptorType

  !which can be introduced for some standard grids, like          

  integer,parameter::CellCentered_=1,Nodes_=2
  interface set_standard_grid_descriptor
     module procedure set_standard_grid_descriptor_id
     module procedure set_standard_grid_descriptor_dd
  end interface

  !or for an arbitrary grid,                                       
  interface set_grid_descriptor_unused
     module procedure set_grid_descriptor_id
     module procedure set_grid_descriptor_dd
  end interface
  !in such a manner that for a given block the actual geometric   
  !displacement of the grid block fragment is equal to the product
  !the product of the Displacement_D by the grid size:
  !\end{verbatim} 
  !EOP             
contains
  !BOP
  !IROUTINE: xyz_grid_d - the coordintes of the grid point
  !INTERFACE:
  function xyz_grid_d(GridDescriptor,&
       lGlobalTreeNode,iGridPoints_D)
    !INPUT ARGUMENTS:     
    type(GridDescriptorType),intent(in)::GridDescriptor                 
    integer,intent(in)::lGlobalTreeNode
    integer,dimension(GridDescriptor%nDim),intent(in)::&
         iGridPoints_D
    !OUTPUT ARGUMENTS:
    real,dimension(GridDescriptor%nDim)::&
         xyz_grid_d
    !EOP
    xyz_grid_d=xyz_block_d(GridDescriptor%DD%Ptr,&
         lGlobalTreeNode)+& 
         d_xyz_cell_d(GridDescriptor%DD%Ptr,&
         lGlobalTreeNode)*&
         (GridDescriptor%Displacement_D-cHalf+real(iGridPoints_D))
  end function xyz_grid_d
  !===============================================================!
  !BOP
  !IROUTINE: set_standard_grid_descriptor - cell centered or node grids.
  !DESCRIPTION:
  !Allow to set the standard grid descriptor: 
  !Standard\_=CellCentered\_ or Standard\_=Nodes\_,
  !with or without halo points ("ghost points")
  !EOP                  
  !---------------------------------------------------------------!
  !BOP
  !INTERFACE:
  subroutine set_standard_grid_descriptor_id(&
       GridID_,&
       nGhostGridPoints,&
       Standard_,&
       GridDescriptor)
    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_
    integer,intent(in),optional::Standard_
    integer,intent(in),optional::nGhostGridPoints
    !OUTPUT ARGUMENTS:
    type(GridDescriptorType),intent(out)::GridDescriptor
    !EOP
    integer::iError,MyStandard_,nGhostGridPointsMy
    call associate_dd_pointer(&
         GridID_,&
         GridDescriptor%DD)
    GridDescriptor%nDim=ndim_grid(GridID_)
    allocate(GridDescriptor%Displacement_D(&
         1:GridDescriptor%nDim),stat=iError)
    call check_allocate(iError,"Displacement_D")
    allocate(GridDescriptor%iGridPointMin_D(&
         1:GridDescriptor%nDim),stat=iError)
    call check_allocate(iError,"iGridPointMins_D")
    allocate(GridDescriptor%iGridPointMax_D(&
         1:GridDescriptor%nDim),stat=iError)
    call check_allocate(iError,"iGridPointMaxs_D")
    if(present(nGhostGridPoints))then
       nGhostGridPointsMy=nGhostGridPoints
    else
       nGhostGridPointsMy=0
    end if
    if(present(Standard_))then
       MyStandard_=Standard_
    else
       MyStandard_=CellCentered_
    end if
    select case(MyStandard_)
    case(CellCentered_)
       GridDescriptor%Displacement_D=cZero
       GridDescriptor%iGridPointMin_D=&
            1-nGhostGridPointsMy
       GridDescriptor%iGridPointMax_D=&
            ncells_decomposition_D(&
            GridID_)+nGhostGridPointsMy
    case(Nodes_)
       GridDescriptor%Displacement_D=-cHalf
       GridDescriptor%iGridPointMin_D=&
            min(1,2-nGhostGridPointsMy)
       GridDescriptor%iGridPointMax_D=&
            ncells_decomposition_D(&
            GridID_)+nGhostGridPointsMy
    case default
       call CON_stop(&
            'Unknown standard for Grid Descriptor',Standard_)
    end select
  end subroutine set_standard_grid_descriptor_id
  !---------------------------------------------------------------!
  !BOP
  !INTERFACE:
  subroutine set_standard_grid_descriptor_dd(&
       DomainDecomposition,&
       nGhostGridPoints,&
       Standard_,&
       GridDescriptor)
    !INPUT ARGUMENTS:
    type(DomainDecompositionType),target,intent(in)::&
         DomainDecomposition
    integer,intent(in),optional::Standard_
    integer,intent(in),optional::nGhostGridPoints
    !OUTPUT ARGUMENTS:
    type(GridDescriptorType),intent(out)::GridDescriptor
    !EOP
    integer::iError,MyStandard_,nGhostGridPointsMy
    nullify(GridDescriptor%DD%Ptr)
    call associate_dd_pointer(&
         DomainDecomposition,&
         GridDescriptor%DD)
    GridDescriptor%nDim=ndim_grid(&
         DomainDecomposition)
    allocate(GridDescriptor%Displacement_D(&
         1:GridDescriptor%nDim),stat=iError)
    call check_allocate(iError,"Displacement_D")
    allocate(GridDescriptor%iGridPointMin_D(&
         1:GridDescriptor%nDim),stat=iError)
    call check_allocate(iError,"iGridPointMins_D")
    allocate(GridDescriptor%iGridPointMax_D(&
         1:GridDescriptor%nDim),stat=iError)
    call check_allocate(iError,"iGridPointMaxs_D")
    if(present(nGhostGridPoints))then
       nGhostGridPointsMy=nGhostGridPoints
    else
       nGhostGridPointsMy=0
    end if
    if(present(Standard_))then
       MyStandard_=Standard_
    else
       MyStandard_=CellCentered_
    end if
    select case(MyStandard_)
    case(CellCentered_)
       GridDescriptor%Displacement_D=cZero
       GridDescriptor%iGridPointMin_D=1-nGhostGridPointsMy
       GridDescriptor%iGridPointMax_D=&
            ncells_decomposition_d(&
            DomainDecomposition)+nGhostGridPointsMy
    case(Nodes_)
       GridDescriptor%Displacement_D=-cHalf
       GridDescriptor%iGridPointMin_D=&
            min(1,2-nGhostGridPointsMy)
       GridDescriptor%iGridPointMax_D=&
            ncells_decomposition_d(&
            DomainDecomposition)+nGhostGridPointsMy
    case default
       call CON_stop(&
            'Unknown standard for Grid Descriptor',Standard_)
    end select
  end subroutine set_standard_grid_descriptor_dd
  !===============================================================!
  !More general grid descriptor                                   !
  !---------------------------------------------------------------!
  !BOP
  !IROUTINE: set_grid_descriptor - set more general grid
  !INTERFACE:
  subroutine set_grid_descriptor_id(&
       GridID_,&
       nDim,&
       iGridPointMin_D,&
       iGridPointMax_D,&
       Displacement_D,&
       GridDescriptor)
    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_
    integer,intent(in)::nDim
    integer,intent(in),dimension(nDim)::iGridPointMin_D
    integer,intent(in),dimension(nDim)::iGridPointMax_D
    real,intent(in),dimension(nDim)::Displacement_D
    !OUTPUT ARGUMENTS:
    type(GridDescriptorType),intent(out)::GridDescriptor
    !EOP
    integer::iError,MyStandard_,nGhostGridPointsMy
    call associate_dd_pointer(&
         GridID_,&
         GridDescriptor%DD)
    GridDescriptor%nDim=nDim
    allocate(GridDescriptor%Displacement_D(&
         1:GridDescriptor%nDim),stat=iError)
    call check_allocate(iError,"Displacement_D")
    GridDescriptor%Displacement_D=Displacement_D
    allocate(GridDescriptor%iGridPointMin_D(&
         1:GridDescriptor%nDim),stat=iError)
    call check_allocate(iError,"iGridPointMins_D")
    GridDescriptor%iGridPointMin_D=iGridPointMin_D
    allocate(GridDescriptor%iGridPointMax_D(&
         1:GridDescriptor%nDim),stat=iError)
    call check_allocate(iError,"iGridPointMaxs_D")
    GridDescriptor%iGridPointMax_D=iGridPointMax_D
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
       GridDescriptor)
    !INPUT ARGUMENTS:
    type(DomainDecompositionType),target,intent(in)::&
         DomainDecomposition
    integer,intent(in)::nDim
    integer,intent(in),dimension(nDim)::iGridPointMin_D
    integer,intent(in),dimension(nDim)::iGridPointMax_D
    real,intent(in),dimension(nDim)::Displacement_D
    !OUTPUT ARGUMENTS:
    type(GridDescriptorType),intent(out)::GridDescriptor
    !EOP
    integer::iError,MyStandard_,nGhostGridPointsMy
    call associate_dd_pointer(&
         DomainDecomposition,&
         GridDescriptor%DD)
    GridDescriptor%nDim=nDim
    allocate(GridDescriptor%Displacement_D(&
         1:GridDescriptor%nDim),stat=iError)
    call check_allocate(iError,"Displacement_D")
    GridDescriptor%Displacement_D=Displacement_D
    allocate(GridDescriptor%iGridPointMin_D(&
         1:GridDescriptor%nDim),stat=iError)
    call check_allocate(iError,"iGridPointMins_D")
    GridDescriptor%iGridPointMin_D=iGridPointMin_D
    allocate(GridDescriptor%iGridPointMax_D(&
         1:GridDescriptor%nDim),stat=iError)
    call check_allocate(iError,"iGridPointMaxs_D")
    GridDescriptor%iGridPointMax_D=iGridPointMax_D
  end subroutine set_grid_descriptor_dd
  !===============================================================!
  !BOP
  !IROUTINE: clean_grid_descriptor
  !INTERFACE:
  subroutine clean_grid_descriptor(GridDescriptor)
    !INPUT ARGUMENTS:
    type(GridDescriptorType),intent(inout)::&
         GridDescriptor
    !EOP
    deallocate(GridDescriptor%iGridPointMin_D)
    deallocate(GridDescriptor%iGridPointMax_D)
    deallocate(GridDescriptor%Displacement_D)
    nullify(GridDescriptor%DD%Ptr)
  end subroutine clean_grid_descriptor
  !===============================================================!
  !Number of cells per block

  integer function n_grid_points_per_block(GridDescriptor)
    type(GriddescriptorType)::GridDescriptor
    n_grid_points_per_block=product(&
         Griddescriptor%iGridPointMax_D+1-&
         GridDescriptor%iGridPointMin_D)
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

  subroutine global_i_grid_point_to_icb(&
       GridDescriptor,&
       iGridPointGlobal,&
       lGlobalTreeNode,&
       iGridPoint_D)
    type(GridDescriptorType),intent(in)::GridDescriptor
    integer,intent(in)::iGridPointGlobal
    integer,intent(out):: lGlobalTreeNode
    integer,dimension(GridDescriptor%nDim),intent(out)::&
         iGridPoint_D
    integer,dimension(GridDescriptor%nDim)::nGridPoints_D 
    integer::iDim,iMisc,iMisc1         
    nGridPoints_D=1+GridDescriptor%iGridPointMax_D-&
         GridDescriptor%iGridPointMin_D
    iMisc=iGridPointGlobal-1
    do idim=1,GridDescriptor%nDim
       iMisc1=mod(iMisc,nGridPoints_D(iDim))
       iGridPoint_D(iDim)=GridDescriptor%iGridPointMin_D(iDim)+iMisc1
       iMisc=(iMisc-iMisc1)/nGridPoints_D(iDim)
    end do
    lGlobalTreeNode=iMisc+1
  end subroutine global_i_grid_point_to_icb
  !===============================================================!
  !The inverse procedure                                          !

  integer function i_grid_point_global(&
       GridDescriptor,&
       lGlobalTreeNode,&
       iGridPoint_D)
    type(GridDescriptorType),intent(in)::GridDescriptor
    integer,intent(in):: lGlobalTreeNode
    integer,dimension(GridDescriptor%nDim),intent(in)::&
         iGridPoint_D
    integer,dimension(0:GridDescriptor%nDim)::nGridPoints_D
    integer::iDim,iMisc,iMisc1
    nGridPoints_D(1:GridDescriptor%nDim)=1+&
         GridDescriptor%iGridPointMax_D-&
         GridDescriptor%iGridPointMin_D
    nGridPoints_D(0)=1
    i_grid_point_global=(lGlobalTreeNode-1)*nGridPoints_D(&
         GridDescriptor%nDim)
    do idim=GridDescriptor%nDim,1,-1
       i_grid_point_global=i_grid_point_global*&
            nGridPoints_D(iDim-1)+iGridPoint_D(iDim)-&
            GridDescriptor%iGridPointMin_D(iDim)
    end do
    i_grid_point_global=i_grid_point_global+1
  end function i_grid_point_global
  !^CFG IF USEMCT BEGIN
  !BOP
  !IROUTINE: grid_descriptor_to_gsmap - transformation to MCT global SegMap
  !DESCRIPTION:
  !TRANSFORMATION OF A GRID DESCRIPTOR TO A GLOBAL
  !                 SEGMENTATION MAP FOR MCT
  !
  !At the time the transformation to the MCT SegMap can only be   
  !done for a grid with no ghost points.   
  !EOP
  !===============================================================!
  !
  !REVISION HISTORY:
  !I.Sokolov<igorsok@umich.edu> & J.W. Larson <larson@mcs.anl.gov>
  ! 7.18.2003                                                     

  subroutine grid_descriptor_to_gsmap(GridDescriptor&
       !^CFG END USEMCT 
       !       ,GSMap&                 !^CFG UNCOMMENT IF USEMCT
    )                        !^CFG IF USEMCT
    !       use m_GlobalSegMap,&    !^CFG UNCOMMENT IF USEMCT
    !           ONLY:GlobalSegMap   !^CFG UNCOMMENT IF USEMCT
    !^CFG IF USEMCT BEGIN
    !---------------------------------------------------------------!
    type(GridDescriptorType),intent(in)::GridDescriptor
    !^CFG END USEMCT 
    !    type(GlobalSegMap),&               !^CFG UNCOMMENT IF USEMCT
    !         intent(out)::GSMap            !^CFG UNCOMMENT IF USEMCT
    integer::iError,iGlobalTreeNode,iCounter      !^CFG IF USEMCT
    !---------------------------------------------------------------!
    !    GSMap%comp_id=MCTCompID_           !^CFG UNCOMMENT IF USEMCT
    !    GSMap%ngseg=n_block_total(&        !^CFG UNCOMMENT IF USEMCT
    !         GridDescriptor%DD%Ptr)        !^CFG UNCOMMENT IF USEMCT
    !    GSMap%gsize= GSMap%ngseg*&         !^CFG UNCOMMENT IF USEMCT
    !         n_grid_points_per_block(&     !^CFG UNCOMMENT IF USEMCT
    !         GridDescriptor)               !^CFG UNCOMMENT IF USEMCT
    !    allocate(GSMap%start(GSMap%ngseg),&!^CFG UNCOMMENT IF USEMCT
    !         GSMap%length(GSMap%ngseg), &  !^CFG UNCOMMENT IF USEMCT
    !         GSMap%pe_loc(GSMap%ngseg), &  !^CFG UNCOMMENT IF USEMCT
    !         stat = iError)                !^CFG UNCOMMENT IF USEMCT
    !    GSMap%length=n_grid_points_per_block(&    !^CFG UNCOMMENT IF USEMCT
    !    GridDescriptor)                    !^CFG UNCOMMENT IF USEMCT
    !^CFG IF USEMCT BEGIN
    call check_allocate(iError,'GSMap arrays')
    iCounter=0
    do iGlobalTreeNode=1,ntree_nodes(GridDescriptor%DD%Ptr)
       if(.not.is_used_block(&
            GridDescriptor%DD%Ptr,iGlobalTreeNode))CYCLE
       iCounter=iCounter+1
       !^CFG END USEMCT
       !       GSMap%start(iCounter)=&         !^CFG UNCOMMENT IF USEMCT
       !            (iGlobalTreeNode-1)*&      !^CFG UNCOMMENT IF USEMCT
       !            n_grid_points_per_block(&  !^CFG UNCOMMENT IF USEMCT
       !            GridDescriptor)+1          !^CFG UNCOMMENT IF USEMCT
       !       GSMap%pe_loc(iCounter)=&        !^CFG UNCOMMENT IF USEMCT
       !            pe_decomposition(&         !^CFG UNCOMMENT IF USEMCT
       !            GridDescriptor%DD%Ptr,&    !^CFG UNCOMMENT IF USEMCT
       !            iGlobalTreeNode)           !^CFG UNCOMMENT IF USEMCT
       !^CFG IF USEMCT BEGIN
    end do
  end subroutine grid_descriptor_to_gsmap
  !^CFG END USEMCT
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
  subroutine nearest_grid_points(Xyz_D,&
       GridDescriptor,&
       nIndexes,&
       Index_II,&
       nImages,Weight_I)
    !INPUT ARGUMENTS:                       
    type(GridDescriptorType)::GridDescriptor     
    real,dimension(GridDescriptor%nDim),intent(inout)::Xyz_D 
    integer,intent(in)::nIndexes
    !OUTPUT ARGUMENTS:
    integer,dimension(0:nIndexes,2**GridDescriptor%nDim)::&
         Index_II
    integer,intent(out)::nImages
    real,dimension(2**GridDescriptor%nDim),intent(out)::&
         Weight_I
    !EOP
    real,dimension(GridDescriptor%nDim)::&
         XyzStored_D, DXyzCells_D,DXyzTolerance_D
    integer,dimension(GridDescriptor%nDim)::iGridPoints_D
    logical,dimension(GridDescriptor%nDim)::IsAtFace_D
    real,dimension(GridDescriptor%nDim,2**GridDescriptor%nDim)&
         :: Xyz_DI
    real,parameter::Tolerance=cOne/cE3
    integer::lGlobalTreeNode,iDim,iImages

    Index_II=0;Weight_I=cZero
    XyzStored_D=Xyz_D
    call search_in(GridDescriptor%DD%Ptr,Xyz_D,lGlobalTreeNode)
    call pe_and_blk(GridDescriptor%DD%Ptr,lGlobalTreeNode,&
         Index_II(0,1),Index_II(nIndexes,1))
    DXyzCells_D=d_xyz_cell_d(&
         GridDescriptor%DD%Ptr,lGlobalTreeNode)
    Xyz_D=Xyz_D-DXyzCells_D*GridDescriptor%Displacement_D
    iGridPoints_D=floor(Xyz_D/DXyzCells_D)
    Xyz_D=Xyz_D-DXyzCells_D*iGridPoints_D
    Index_II(1:GridDescriptor%nDim,1)=iGridPoints_D+1
    nImages=1;Weight_I(1)=cOne

    DXyzTolerance_D=Tolerance*DXyzCells_D
    IsAtFace_D=abs(Xyz_D)<DXyzTolerance_D&
         .or.abs(Xyz_D-DXyzCells_D )<DXyzTolerance_D
    if(.not.(any(IsAtFace_D)))return

    Xyz_DI(:,1)=XyzStored_D
    do iDim=1,GridDescriptor%nDim
       if(IsAtFace_D(iDim))then
          Xyz_DI(:,1+nImages:nImages+nImages)=Xyz_DI(:,1:nImages)
          Xyz_DI(iDim,1:nImages)=Xyz_DI(iDim,1:nImages)-&
               DXyzTolerance_D(iDim)
          Xyz_DI(iDim,1+nImages:nImages+nImages)=Xyz_DI(&
               iDim,1+nImages:nImages+nImages)+&
               DXyzTolerance_D(iDim)
          nImages=nImages+nImages
       end if
    end do
    do iImages=1,nImages
       call search_in(GridDescriptor%DD%Ptr,&
            Xyz_DI(:,iImages),&
            lGlobalTreeNode)
       call pe_and_blk(GridDescriptor%DD%Ptr,&
            lGlobalTreeNode,&
            Index_II(0,iImages),&
            Index_II(nIndexes,iImages))
       Xyz_DI(:,iImages)=Xyz_DI(:,iImages)-&
            d_xyz_cell_d(GridDescriptor%DD%Ptr,&
            lGlobalTreeNode)*GridDescriptor%Displacement_D
       call search_cell(GridDescriptor%DD%Ptr,&
            lGlobalTreeNode,Xyz_DI(:,iImages),&
            Index_II(1:GridDescriptor%nDim,iImages))
    end do
    Weight_I(1:nImages)=cOne/real(nImages)
  end subroutine nearest_grid_points
  !===============================================================!
  !=================SECOND ORDER INTERPOLATION====================!
  !BOP
  !IROUTINE: bilinear_interpolation - second order interpolation
  !EOP
  !BOP
  !DESCRIPTION:
  !This is a bilinear interpoaltion using the grid points,        
  !described with the grid descriptor.                            
  !EOP
  !INTERFACE:
  subroutine bilinear_interpolation(Xyz_D,&
       GridDescriptor,&
       nIndexes,& 
       Index_II,&
       nImages,Weight_I)
    !INPUT ARGUMENTS:                  
    type(GridDescriptorType)::GridDescriptor     
    real,dimension(GridDescriptor%nDim),intent(inout)::Xyz_D 
    integer,intent(in)::nIndexes
    !OUTPUT ARGUMENTS
    integer,dimension(0:nIndexes,2**GridDescriptor%nDim)::&
         Index_II
    integer,intent(out)::nImages
    real,dimension(2**GridDescriptor%nDim),intent(out)::Weight_I
    !EOP
    real,dimension(GridDescriptor%nDim)::&
         XyzResid_D
    integer,dimension(GridDescriptor%nDim)::iGridPoints_D
    real,dimension(GridDescriptor%nDim):: XyzMisc_D
    integer::lNodeMisc
    integer::lGlobalTreeNode,iDim,iImages,iNewStart,nImagesNew
    real::WeightLeft

    Index_II=0;Weight_I=cZero

    call search_in(GridDescriptor%DD%Ptr,Xyz_D,lGlobalTreeNode)
    call pe_and_blk(GridDescriptor%DD%Ptr,lGlobalTreeNode,&
         Index_II(0,1),Index_II(nIndexes,1))

    XyzResid_D=Xyz_D/d_xyz_cell_d(&
         GridDescriptor%DD%Ptr,lGlobalTreeNode)-&
         GridDescriptor%Displacement_D+cHalf
    iGridPoints_D=floor(XyzResid_D)
    XyzResid_D=XyzResid_D-real(iGridPoints_D)

    Index_II(1:GridDescriptor%nDim,1)=iGridPoints_D
    nImages=1;Weight_I(1)=cOne

    !Thus calculated XyzResid_D satisfies the inequilities as       !
    !follow:XyzResid_D>=0 and XyzResid_D<1. It is used to calculte  !
    !weights for the eight grid points, among them the iGridPoints_D!
    !being the left one with respect to all the spatial directions  !
    do iDim=1,GridDescriptor%nDim
       if(XyzResid_D(iDim)<cTiny)CYCLE
       iNewStart=nImages+1;nImagesNew=nImages+nImages
       Index_II(:,iNewStart:nImagesNew)=&
            Index_II(:,1:nImages)
       Index_II(iDim,iNewStart:nImagesNew)=&
            Index_II(iDim,iNewStart:nImagesNew)+1
       Weight_I(iNewStart:nImagesNew)= &
            Weight_I(1:nImages)*XyzResid_D(iDim)
       WeightLeft=cOne-XyzResid_D(iDim)
       Weight_I(1:nImages)= WeightLeft* Weight_I(1:nImages)
       nImages=nImagesNew
    end do
    !---------------------------------------------------------------!
    !Check if the grid point index is within the index limits       !

    do iImages=1,nImages
       if(all(Index_II(1:GridDescriptor%nDim,iImages)>=&
            GridDescriptor%iGridPointMin_D)&
            .and.all(Index_II(1:GridDescriptor%nDim,iImages)<=&
            GridDescriptor%iGridPointMax_D))CYCLE
       XyzMisc_D=xyz_grid_d(GridDescriptor,&
            lGlobalTreeNode,Index_II(&
            1:GridDescriptor%nDim,iImages)) 
       call search_in(GridDescriptor%DD%Ptr,XyzMisc_D,lNodeMisc)
       call pe_and_blk(GridDescriptor%DD%Ptr,lNodeMisc,&
            Index_II(0,iImages),Index_II(nIndexes,iImages))

       XyzMisc_D=XyzMisc_D-&
            d_xyz_cell_d(GridDescriptor%DD%Ptr,&
            lNodeMisc)*GridDescriptor%Displacement_D
       call search_cell(GridDescriptor%DD%Ptr,&
            lNodeMisc,XyzMisc_D,&
            Index_II(1:GridDescriptor%nDim,iImages))
    end do
  end subroutine bilinear_interpolation
  !==============================END==============================!
end Module CON_grid_descriptor
!\end{verbatim}                     !^CFG UNCOMMENT IF PRINTPS  !
!\end{document}                     !^CFG UNCOMMENT IF PRINTPS  !
!=============================LINE 643==========================!
!BOP
!IROUTINE: more methods
!DESCRIPTION:
!More methods are available, see the source file.
!EOP
