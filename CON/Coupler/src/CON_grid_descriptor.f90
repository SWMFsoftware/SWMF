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
       call CON_stop('Unknown standard for Grid Descriptor')
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
    lGlobalTreeNode=i_global_node_a(GridDescriptor%DD%Ptr,iMisc+1)
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
    integer,dimension(GridDescriptor%nDim)::nGridPoints_D
    integer::iDim,iMisc,iMisc1
    nGridPoints_D(1:GridDescriptor%nDim)=1+&
         GridDescriptor%iGridPointMax_D-&
         GridDescriptor%iGridPointMin_D
    i_grid_point_global=i_global_block(GridDescriptor%DD%Ptr,&
         lGlobalTreeNode)-1
    do idim=GridDescriptor%nDim,1,-1
       i_grid_point_global=i_grid_point_global*&
            nGridPoints_D(iDim)+iGridPoint_D(iDim)-&
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
    integer::iError,iGlobalTreeNode,iBlockAll      !^CFG IF USEMCT
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
    do iBlockAll=1,n_block_total(GridDescriptor%DD%Ptr)
       iGlobalTreeNode=i_global_node_a(GridDescriptor%DD%Ptr,iBlockAll)
       !^CFG END USEMCT
       !       GSMap%start(iBlockAll)=&        !^CFG UNCOMMENT IF USEMCT
       !            (iGlobalTreeNode-1)*&      !^CFG UNCOMMENT IF USEMCT
       !            n_grid_points_per_block(&  !^CFG UNCOMMENT IF USEMCT
       !            GridDescriptor)+1          !^CFG UNCOMMENT IF USEMCT
       !       GSMap%pe_loc(iBlockAll)=&       !^CFG UNCOMMENT IF USEMCT
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
  subroutine nearest_grid_points(nDim,Xyz_D,&
       GridDescriptor,&
       nIndexes,&
       Index_II,&
       nImages,Weight_I)
    !INPUT ARGUMENTS:                       
    type(GridDescriptorType)::GridDescriptor
    integer,intent(in)::nDim
    real,dimension(nDim),intent(inout)::Xyz_D 
    integer,intent(in)::nIndexes
    !OUTPUT ARGUMENTS:
    integer,dimension(0:nIndexes,2**nDim),intent(out)::&
         Index_II
    integer,intent(out)::nImages
    real,dimension(2**nDim),intent(out)::&
         Weight_I
    !EOP
    real,dimension(nDim)::&
         XyzStored_D, DXyzCells_D,DXyzTolerance_D
    integer,dimension(nDim)::iGridPoints_D
    logical,dimension(nDim)::IsAtFace_D
    real,dimension(nDim,2**nDim)&
         :: Xyz_DI
    real,parameter::Tolerance = 0.001
    integer::lGlobalTreeNode,iDim,iImages
    logical,dimension(nDim,2**nDim)::&
         Up_DI,Down_DI
    logical,dimension(nDim)::&
         IsDomainBoundaryUp_D,IsDomainBoundaryDown_D

    Index_II=0;Weight_I=cZero
    XyzStored_D=Xyz_D
    Index_II(1:nDim,1)=&
         GridDescriptor%iGridPointMax_D+1
    do while(any(Index_II(1:nDim,1)>&
         GridDescriptor%iGridPointMax_D))
       Xyz_D=XyzStored_D
       call search_in(GridDescriptor%DD%Ptr,Xyz_D,lGlobalTreeNode)
       DXyzCells_D=d_xyz_cell_d(&
            GridDescriptor%DD%Ptr,lGlobalTreeNode)
       where(is_right_boundary_d(GridDescriptor%DD%Ptr,lGlobalTreeNode))
          Xyz_D=min(Xyz_D,cAlmostOne*d_xyz_block_d(&         
               GridDescriptor%DD%Ptr,lGlobalTreeNode))
          XyzStored_D=Xyz_D+xyz_block_d(&         
               GridDescriptor%DD%Ptr,lGlobalTreeNode)
       end where
       where(is_left_boundary_d(GridDescriptor%DD%Ptr,lGlobalTreeNode))
          Xyz_D=max(Xyz_D,cZero)
          XyzStored_D=Xyz_D+xyz_block_d(&         
               GridDescriptor%DD%Ptr,lGlobalTreeNode)
       end where
       call pe_and_blk(GridDescriptor%DD%Ptr,lGlobalTreeNode,&
            Index_II(0,1),Index_II(nIndexes,1))
       Xyz_D=Xyz_D-DXyzCells_D*GridDescriptor%Displacement_D
       iGridPoints_D=floor(Xyz_D/DXyzCells_D)
       Xyz_D=Xyz_D-DXyzCells_D*iGridPoints_D
       Index_II(1:nDim,1)=iGridPoints_D+1
       DXyzTolerance_D=Tolerance*DXyzCells_D
       IsAtFace_D=abs(Xyz_D)<DXyzTolerance_D&
            .or.abs(Xyz_D-DXyzCells_D )<DXyzTolerance_D
       where(is_right_boundary_d(GridDescriptor%DD%Ptr,lGlobalTreeNode))&
            Index_II(1:nDim,1)=&
            min(Index_II(1:nDim,1),&
            GridDescriptor%iGridPointMax_D)
       where(Index_II(1:nDim,1)>&
            GridDescriptor%iGridPointMax_D)&
            XyzStored_D=XyzStored_D+DXyzTolerance_D*0.25
    end do
    nImages=1;Weight_I(1)=cOne

   

    if(.not.(any(IsAtFace_D)))return
    Xyz_DI(:,1)=XyzStored_D
    do iDim=1,nDim
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
            Index_II(1:nDim,iImages))
    end do
    iImages=1
    do while(iImages<=nImages)
       !Exclude the stencil nodes which are out of the 
       !computational domain
       do while(any(GridDescriptor%iGridPointMin_D>&
            Index_II(1:nDim,iImages).or.&
            GridDescriptor%iGridPointMax_D<&
            Index_II(1:nDim,iImages)))
          if(iImages==nImages)then
             nImages=nImages-1
             EXIT
          end if
          Index_II(:,iImages)=Index_II(:,nImages)
          nImages=nImages-1
       end do
       iImages=iImages+1
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
  !This is a bilinear interpolation using the grid points,        
  !described with the grid descriptor.                            
  !EOP
  !INTERFACE:
  subroutine bilinear_interpolation(&
       nDim,&
       Xyz_D,&
       GridDescriptor,&
       nIndexes,& 
       Index_II,&
       nImages,Weight_I)
    !INPUT ARGUMENTS:
    integer,intent(in)      ::nDim
    type(GridDescriptorType)::GridDescriptor     
    real,intent(inout)::Xyz_D(nDim)
    integer,intent(in)::nIndexes
    !OUTPUT ARGUMENTS
    integer           ::Index_II(0:nIndexes,2**nDim)
         
    integer,intent(out)::nImages
    real,dimension(2**nDim),intent(out)::Weight_I
    !EOP
    real,dimension(nDim)::&
         XyzResid_D,XyzStored_D
    integer,dimension(nDim)::iGridPoints_D
    real,dimension(nDim):: XyzMisc_D
    integer,dimension(2**nDim)::lNode_I
    integer::lGlobalTreeNode,iDim,iImages,iNewStart,nImagesNew
    real::WeightLeft
    logical,dimension(nDim,2**nDim)::&
         Up_DI,Down_DI
    real,dimension(nDim,0:2**nDim)::&
         DXyz_DI
    logical,dimension(nDim)::&
         IsDomainBoundaryUp_D,IsDomainBoundaryDown_D

    Index_II=0
    Weight_I=cZero
    XyzStored_D=Xyz_D
    Up_DI(:,1)=.true.
    call search_in(GridDescriptor%DD%Ptr,Xyz_D,lGlobalTreeNode)
    !Find global node number, PE and number which involved the point
    do while(any(Up_DI(:,1)))
       !This do loop works more than once only in case of node
       !grid and only in case of Xyz_D point belonging to a block
       !boundary. At this case the routine needs help in deciding,
       !to which block this point should be assigned. Otherwise,
       !automatically Up_D(:,1)=.false. after first loop pass.
       Xyz_D=XyzStored_D-xyz_block_d(&
            GridDescriptor%DD%Ptr,lGlobalTreeNode)
       IsDomainBoundaryUp_D=&
            is_right_boundary_d(GridDescriptor%DD%Ptr,lGlobalTreeNode)
       IsDomainBoundaryDown_D=&
            is_left_boundary_d(GridDescriptor%DD%Ptr,lGlobalTreeNode)
       where(IsDomainBoundaryUp_D)
          Xyz_D=min(Xyz_D,cAlmostOne*d_xyz_block_d(&         
               GridDescriptor%DD%Ptr,lGlobalTreeNode))
          XyzStored_D=Xyz_D+xyz_block_d(&         
               GridDescriptor%DD%Ptr,lGlobalTreeNode)
       end where
       where(IsDomainBoundaryDown_D)
          Xyz_D=max(Xyz_D,cZero)
          XyzStored_D=Xyz_D+xyz_block_d(&         
               GridDescriptor%DD%Ptr,lGlobalTreeNode)
       end where

       call pe_and_blk(GridDescriptor%DD%Ptr,lGlobalTreeNode,&
            Index_II(0,1),Index_II(nIndexes,1))
       
       !Find DXyz for this block
       DXyz_DI(:,0)=d_xyz_cell_d(&         
            GridDescriptor%DD%Ptr,lGlobalTreeNode)
       
       !\
       !/
       
       XyzResid_D=Xyz_D/DXyz_DI(:,0)-&
            GridDescriptor%Displacement_D+cHalf
       iGridPoints_D=floor(XyzResid_D)
       XyzResid_D=XyzResid_D-real(iGridPoints_D)
       !Thus calculated XyzResid_D satisfies the inequalities as       !
       !follow:XyzResid_D>=0 and XyzResid_D<1. It is used to calculte  !
       !weights for the eight grid points, among them the iGridPoints_D!
       !being the left one with respect to all the spatial directions  !
       Index_II(1:nDim,1)=iGridPoints_D
       Up_DI(:,1)=Index_II(1:nDim,1)>&
            GridDescriptor%iGridPointMax_D
       Down_DI(:,1)=Index_II(1:nDim,1)<&
            GridDescriptor%iGridPointMin_D
       if(any(Up_DI(:,1)))&
            lGlobalTreeNode=l_neighbor(&
            GridDescriptor%DD%Ptr,lGlobalTreeNode,&
            Index_II(1:nDim,1))
    end do

    nImages=1;Weight_I(1)=cOne

    do iDim=1,nDim
       !Exclude the stencil nodes which are out of the 
       !computational domain
       if(Down_DI(iDim,1).and.IsDomainBoundaryDown_D(iDim))then
          Down_DI(iDim,1:nImages)=.false.
          Index_II(iDim,1:nImages)=Index_II(iDim,1:nImages)+1
          CYCLE
       end if
       
       if(Index_II(iDim,1)==GridDescriptor%iGridPointMax_D(iDim)&
            .and.IsDomainBoundaryUp_D(iDim))CYCLE
       if(XyzResid_D(iDim)<cTiny)CYCLE
       iNewStart=nImages+1;nImagesNew=nImages+nImages
       Index_II(:,iNewStart:nImagesNew)=&
            Index_II(:,1:nImages)
       Index_II(iDim,iNewStart:nImagesNew)=&
            Index_II(iDim,iNewStart:nImagesNew)+1
       Up_DI(:,iNewStart:nImagesNew)=Up_DI(:,1:nImages)
       Up_DI(iDim,iNewStart:nImagesNew)=&
            Index_II(iDim,iNewStart)>&
            GridDescriptor%iGridPointMax_D(iDim)
       Down_DI(:,iNewStart:nImagesNew)=Down_DI(:,1:nImages)
       Down_DI(iDim,iNewStart:nImagesNew)=.false.
       Weight_I(iNewStart:nImagesNew)= &
            Weight_I(1:nImages)*XyzResid_D(iDim)
       WeightLeft=cOne-XyzResid_D(iDim)
       Weight_I(1:nImages)= WeightLeft* Weight_I(1:nImages)
       nImages=nImagesNew
    end do
    !---------------------------------------------------------------!
    !Check if the grid point index is within the index limits       !
    do iImages=1,nImages
       if(.not.(any(Up_DI(:,iImages)).or.&
            any(Down_DI(:,iImages))))then
          DXyz_DI(:,iImages)=DXyz_DI(:,0)
          lNode_I(iImages)=lGlobalTreeNode
       else
          lNode_I(iImages)=l_neighbor(GridDescriptor%DD%Ptr,&
               lGlobalTreeNode,Index_II(&
               1:nDim,iImages))
          XyzMisc_D=xyz_grid_d(GridDescriptor,&
               lGlobalTreeNode,Index_II(&
               1:nDim,iImages))-&
               xyz_block_d(GridDescriptor%DD%Ptr,lNode_I(iImages))
          call pe_and_blk(GridDescriptor%DD%Ptr,lNode_I(iImages),&
               Index_II(0,iImages),Index_II(nIndexes,iImages))
          DXyz_DI(:,iImages)=d_xyz_cell_d(GridDescriptor%DD%Ptr,&
               lNode_I(iImages))
          XyzMisc_D=XyzMisc_D- DXyz_DI(:,iImages)&
               *GridDescriptor%Displacement_D
          call search_cell(GridDescriptor%DD%Ptr,&
               lNode_I(iImages),XyzMisc_D,&
               Index_II(1:nDim,iImages))
       end if
    end do
  end subroutine bilinear_interpolation
  !BOP
  !IROUTINE: interpolation_fix_reschange - second order at resolution change
  !EOP
  !BOP
  !DESCRIPTION:
  !This is a bilinear interpoaltion using the grid points,        
  !described with the grid descriptor.                            
  !EOP
  !INTERFACE:
  subroutine interpolation_fix_reschange(&
       nDim, &
       Xyz_D,&
       GridDescriptor,&
       nIndexes,& 
       Index_II,&
       nImages,Weight_I)
    !INPUT ARGUMENTS:
    integer,intent(in):: nDim
    real,intent(inout):: Xyz_D(nDim)
    type(GridDescriptorType)::GridDescriptor     
    integer,intent(in)::nIndexes
    !OUTPUT ARGUMENTS
    integer,dimension(0:nIndexes,2**nDim)::&
         Index_II
    integer,intent(out)::nImages
    real,dimension(2**nDim),intent(out)::Weight_I
    !EOP
    real,dimension(nDim)::&
         XyzResid_D,XyzStored_D
    integer,dimension(nDim)::iGridPoints_D
    integer,dimension(nDim)::iShift_D
    real,dimension(nDim):: XyzMisc_D
    integer,dimension(2**nDim)::&
         lNeighbor_I,lLevel_I
    integer::lGlobalTreeNode,iDim,iImages,iNewStart,nImagesNew
    real::WeightLeft
    logical,dimension(nDim,2**nDim)::&
         Up_DI,Down_DI
    real,dimension(nDim,0:2**nDim)::&
         DXyz_DI
    logical,dimension(nDim)::&
         IsDomainBoundaryUp_D,IsDomainBoundaryDown_D
    integer,dimension(0:nIndexes,2**nDim)::&
         IndexAux_II
    integer::nImagesAux,nImagesTemp
    real,dimension(2**nDim)::WeightAux_I
    logical::IsResChangeUp,IsResChangeDown
    integer::lNodeResChange,iLoopInternal
    real,dimension(nDim)::OrigCoord_D
    real,dimension(nDim)::ResChangeCoord_D
    integer::iIteration
    integer,parameter::iIterationMax=5
    integer,parameter::Higher_=1,Lower_=-1
    integer,dimension(nDim,2**nDim)::&
         iBlockPosition_DI
    real,dimension(nDim,2**nDim)::&
         UpperCoord_DI,LowerCoord_DI
    real,dimension(nDim)::&
         UpperCoordMax_D,LowerCoordMin_D

    Index_II=0
    Weight_I=cZero
    XyzStored_D=Xyz_D

    call search_in(GridDescriptor%DD%Ptr,Xyz_D,lGlobalTreeNode)
    !Find global node, which involved the point
    lLevel_I(1)=1
    nImages=1
    iIteration=0
    COARSEN:do while(any(lLevel_I(1:nImages)>0))
       !This loop works more than once if the stencil includes
       !blocks which are coarser than lGlobalTreeNode. In this case
       !the coarser block is used as the base for a coraser stencil
       Up_DI(:,1)=.true.

       do while(any(Up_DI(:,1)))
          iIteration=iIteration+1
          if(iIteration==iIterationMax)&
               call CON_stop(&
               'Algorithmic error in interpolation_fix_reschange') 
          !This do loop works more than once only in case of node
          !grid and only in case of Xyz_D point belonging to a block
          !boundary. At this case the routine needs help in deciding,
          !to which block this point should be assigned. Otherwise,
          !automatically Up_D(:,1)=.false. after first loop pass.
          Xyz_D=XyzStored_D-xyz_block_d(&
               GridDescriptor%DD%Ptr,lGlobalTreeNode)
          IsDomainBoundaryUp_D=&
               is_right_boundary_d(GridDescriptor%DD%Ptr,lGlobalTreeNode)
          IsDomainBoundaryDown_D=&
               is_left_boundary_d(GridDescriptor%DD%Ptr,lGlobalTreeNode)
          where(IsDomainBoundaryUp_D)
             Xyz_D=min(Xyz_D,cAlmostOne*d_xyz_block_d(&         
                  GridDescriptor%DD%Ptr,lGlobalTreeNode))
             XyzStored_D=Xyz_D+xyz_block_d(&         
                  GridDescriptor%DD%Ptr,lGlobalTreeNode)
          end where
          where(IsDomainBoundaryDown_D)
             Xyz_D=max(Xyz_D,cZero)
             XyzStored_D=Xyz_D+xyz_block_d(&         
                  GridDescriptor%DD%Ptr,lGlobalTreeNode)
          end where
          call pe_and_blk(GridDescriptor%DD%Ptr,lGlobalTreeNode,&
               Index_II(0,1),Index_II(nIndexes,1))

          !Find DXyz for this block
          DXyz_DI(:,0)=d_xyz_cell_d(&         
               GridDescriptor%DD%Ptr,lGlobalTreeNode)

          !\
          !/

          XyzResid_D=Xyz_D/DXyz_DI(:,0)-&
               GridDescriptor%Displacement_D+cHalf
          iGridPoints_D=floor(XyzResid_D)
          XyzResid_D=XyzResid_D-real(iGridPoints_D)
          !Thus calculated XyzResid_D satisfies the inequalities as       !
          !follow:XyzResid_D>=0 and XyzResid_D<1. It is used to calculte  !
          !weights for the eight grid points, among them the iGridPoints_D!
          !being the left one with respect to all the spatial directions  !
          Index_II(1:nDim,1)=iGridPoints_D
          Up_DI(:,1)=Index_II(1:nDim,1)>&
               GridDescriptor%iGridPointMax_D
          Down_DI(:,1)=Index_II(1:nDim,1)<&
               GridDescriptor%iGridPointMin_D
          if(any(Up_DI(:,1)))&
               lGlobalTreeNode=l_neighbor(&
               GridDescriptor%DD%Ptr,lGlobalTreeNode,&
               Index_II(1:nDim,1))
       end do
       nImages=1;Weight_I(1)=cOne
       iShift_D=0
       do iDim=1,nDim
          !Exclude the stencil nodes which are out of the 
          !computational domain
          if(Down_DI(iDim,1).and.IsDomainBoundaryDown_D(iDim))then
             Down_DI(iDim,1:nImages)=.false.
             Index_II(iDim,1:nImages)=Index_II(iDim,1:nImages)+1
             CYCLE
          end if

          if(Index_II(iDim,1)==GridDescriptor%iGridPointMax_D(iDim)&
               .and.IsDomainBoundaryUp_D(iDim))CYCLE
          if(XyzResid_D(iDim)<cTiny)CYCLE
          iShift_D(iDim)=nImages
          iNewStart=nImages+1;nImagesNew=nImages+nImages
          Index_II(:,iNewStart:nImagesNew)=&
               Index_II(:,1:nImages)
          Index_II(iDim,iNewStart:nImagesNew)=&
               Index_II(iDim,iNewStart:nImagesNew)+1
          Up_DI(:,iNewStart:nImagesNew)=Up_DI(:,1:nImages)
          Up_DI(iDim,iNewStart:nImagesNew)=&
               Index_II(iDim,iNewStart)>&
               GridDescriptor%iGridPointMax_D(iDim)
          Down_DI(:,iNewStart:nImagesNew)=Down_DI(:,1:nImages)
          Down_DI(iDim,iNewStart:nImagesNew)=.false.
          Weight_I(iNewStart:nImagesNew)= &
               Weight_I(1:nImages)*XyzResid_D(iDim)
          WeightLeft=cOne-XyzResid_D(iDim)
          Weight_I(1:nImages)= WeightLeft* Weight_I(1:nImages)
          nImages=nImagesNew
       end do
       !---------------------------------------------------------------!
       !Check if the grid point index is within the index limits       !
       do iImages=1,nImages
          if(.not.(any(Up_DI(:,iImages)).or.&
               any(Down_DI(:,iImages))))then
             DXyz_DI(:,iImages)=DXyz_DI(:,0)
             lNeighbor_I(iImages)=lGlobalTreeNode
             lLevel_I(iImages)=0
          else
             XyzMisc_D=xyz_grid_d(GridDescriptor,&
                  lGlobalTreeNode,Index_II(&
                  1:nDim,iImages))
             call search_in(GridDescriptor%DD%Ptr,&
                  XyzMisc_D,lNeighbor_I(iImages))
             DXyz_DI(:,iImages)=&
                  d_xyz_cell_d(GridDescriptor%DD%Ptr,&
                  lNeighbor_I(iImages))
             lLevel_I(iImages)=int(&
                  2*DXyz_DI(1,iImages)&
                  /DXyz_DI(1,0) -3+cTiny)
             if(lLevel_I(iImages)>0)then
                nImages=iImages
                lGlobalTreeNode=lNeighbor_I(iImages)
                CYCLE COARSEN
             elseif(lLevel_I(iImages)==0)then
                call pe_and_blk(&
                     GridDescriptor%DD%Ptr,lNeighbor_I(iImages),&
                     Index_II(0,iImages),Index_II(nIndexes,iImages))
                Index_II(1:nDim,iImages)=&
                     nint((xyz_grid_d(GridDescriptor,&
                     lGlobalTreeNode,Index_II(&
                     1:nDim,iImages))-&
                     xyz_block_d(&
                     GridDescriptor%DD%Ptr,lNeighbor_I(iImages)))/ &
                     DXyz_DI(:,iImages)+cHalf-&
                     GridDescriptor%Displacement_D)
             end if
          end if
       end do
    end do COARSEN
    iDim=i_dir_of_only_reschange(&
         nDim,&
         iShift_D,&
         nImages,&
         lLevel_I(1:nImages))
    if(iDim==0)return
    if(iDim>0)then
       iImages=1
       IMAGES:do while(iImages<=nImages)
          do while(lLevel_I(iImages)/=0)
             IsResChangeUp=Up_DI(iDim,iImages)
             IsResChangeDown=Down_DI(iDim,iImages)
             lNodeResChange=lNeighbor_I(iImages)
             if(iImages==nImages)then
                nImages=nImages-1
                EXIT IMAGES
             end if
             Index_II(:,iImages)=Index_II(:,nImages)
             Weight_I(iImages)=Weight_I(nImages)
             lNeighbor_I(iImages)=lNeighbor_I(nImages)
             lLevel_I(iImages)=lLevel_I(nImages)
             Up_DI(iDim,iImages)=Up_DI(iDim,nImages)
             Down_DI(iDim,iImages)=Down_DI(iDim,nImages)
             nImages=nImages-1
          end do
          iImages=iImages+1
       end do IMAGES
       Weight_I(1:nImages)=Weight_I(1:nImages)/&
            sum(Weight_I(1:nImages))
       !Construct the image points at the grid layers with
       !the original resolution and changed resolution
       if(IsResChangeUp)then
          XyzMisc_D=xyz_grid_d(GridDescriptor,&
               lGlobalTreeNode,&
               GridDescriptor%iGridPointMax_D)
          OrigCoord_D=XyzStored_D
          OrigCoord_D(iDim)=XyzMisc_D(iDim)
          XyzMisc_D=xyz_grid_d(GridDescriptor,&
               lNodeResChange,&
               GridDescriptor%iGridPointMin_D)
          ResChangeCoord_D=XyzStored_D
          ResChangeCoord_D(iDim)=XyzMisc_D(iDim)
       elseif(IsResChangeDown)then
          XyzMisc_D=xyz_grid_d(GridDescriptor,&
               lGlobalTreeNode,&
               GridDescriptor%iGridPointMin_D)
          OrigCoord_D=XyzStored_D
          OrigCoord_D(iDim)=XyzMisc_D(iDim)
          XyzMisc_D=xyz_grid_d(GridDescriptor,&
               lNodeResChange,&
               GridDescriptor%iGridPointMax_D)
          ResChangeCoord_D=XyzStored_D
          ResChangeCoord_D(iDim)=XyzMisc_D(iDim)
       else
          call CON_stop(&
               'Algorithmic error in interpolation_fix_reschange')
       end if
       WeightLeft=(XyzStored_D(iDim)-ResChangeCoord_D(iDim))/&
            (OrigCoord_D(iDim)-ResChangeCoord_D(iDim))
       Weight_I(1:nImages)=Weight_I(1:nImages)*WeightLeft
       call bilinear_interpolation(&
            nDim,&
            ResChangeCoord_D,&
            GridDescriptor,&
            nIndexes,& 
            IndexAux_II,&
            nImagesAux,WeightAux_I)
       if(nImages+nImagesAux>ubound(Index_II,2))&
            call CON_stop(&
            'Algorithmic error in interpolation_fix_reschange')
       Index_II(:,nImages+1:nImages+nImagesAux)=&
            IndexAux_II(:,1:nImagesAux)
       Weight_I(nImages+1:nImages+nImagesAux)=&
            WeightAux_I(1:nImagesAux)*(cOne-WeightLeft)
       nImages=nImages+nImagesAux
    else
       !Generate iBlockPosition_DI:
       iBlockPosition_DI=0
       do iDim=1,nDim
          if(iShift_D(iDim)==0)CYCLE
          iBlockPosition_DI(&
               :,1+iShift_D(iDim):2*iShift_D(iDim))=&
               iBlockPosition_DI(:,1:iShift_D(iDim))
          iBlockPosition_DI(iDim,1:iShift_D(iDim))=Lower_
          iBlockPosition_DI(&
               iDim,1+iShift_D(iDim):2*iShift_D(iDim))=&
               Higher_
       end do
       !Remove redundant images
       iImages=1
       do while(iImages<nImages)
          nImagesTemp=count(&
               lNeighbor_I(iImages:nImages)==&
               lNeighbor_I(iImages))
          if(nImagesTemp>1)then
             iLoopInternal=iImages+1
             do while(iLoopInternal<=nImages)
                if(lNeighbor_I(iLoopInternal)==&
                     lNeighbor_I(iImages))then
                   iBlockPosition_DI(:,iImages)=&
                        iBlockPosition_DI(:,iImages)+&
                        iBlockPosition_DI(:,iLoopInternal)
                   iBlockPosition_DI(:,iLoopInternal)=&
                        iBlockPosition_DI(:,nImages)
                   lNeighbor_I(iLoopInternal)=&
                        lNeighbor_I(nImages)
                   nImages=nImages-1
                else
                   iLoopInternal=iLoopInternal+1
                end if
             end do
             iBlockPosition_DI(:,iImages)=&
                  iBlockPosition_DI(:,iImages)/nImagesTemp
          end if
          iImages=iImages+1
       end do
       !Generate  UpperCoord_DI,LowerCoord_DI
       do iImages=1,nImages
          UpperCoord_DI(:,iImages)=xyz_grid_d(GridDescriptor,&
               lNeighbor_I(iImages),&
               GridDescriptor%iGridPointMax_D)
          LowerCoord_DI(:,iImages)=xyz_grid_d(GridDescriptor,&
               lNeighbor_I(iImages),&
               GridDescriptor%iGridPointMin_D)
       end do
       !Generate  UpperCoordMax_D,LowerCoordMin_D
       do iDim=1,nDim
          if(all(iBlockPosition_DI(iDim,1:nImages)==0))then
             UpperCoordMax_D(iDim)=cZero
             LowerCoordMin_D(iDim)=cZero
             CYCLE
          end if
          UpperCoordMax_D(iDim)=maxval(&
               UpperCoord_DI(iDim,1:nImages),&
               MASK=iBlockPosition_DI(iDim,1:nImages)==Lower_)
          LowerCoordMin_D(iDim)=minval(&
               LowerCoord_DI(iDim,1:nImages),&
               MASK=iBlockPosition_DI(iDim,1:nImages)==Higher_)
       end do
       !We truncate the formfactor from Lower_ blocks at 
       !LowerCoordMin from Higher_ block and vice versa
       !Again remove redundant images
       nImagesTemp=nImages
       nImages=0
       IMAGE:do iImages=1,nImagesTemp
          XyzMisc_D=XyzStored_D
          WeightLeft=cOne
          DIM: do iDim=1,nDim
             select case(iBlockPosition_DI(iDim,iImages))
             case(0)
                CYCLE DIM
             case(Lower_)
                if(XyzMisc_D(iDim)>=LowerCoordMin_D(iDim))&
                     CYCLE IMAGE
                WeightLeft=WeightLeft*&
                     (LowerCoordMin_D(iDim)-XyzMisc_D(iDim))/&
                     (LowerCoordMin_D(iDim)-&
                     UpperCoord_DI(iDim,iImages))
                XyzMisc_D(iDim)=UpperCoord_DI(iDim,iImages)
             case(Higher_)
                if(XyzMisc_D(iDim)<=UpperCoordMax_D(iDim))&
                     CYCLE IMAGE
                WeightLeft=WeightLeft*&
                     (XyzMisc_D(iDim)-UpperCoordMax_D(iDim))/&
                     (LowerCoord_DI(iDim,iImages)-&
                     UpperCoordMax_D(iDim))
                XyzMisc_D(iDim)=LowerCoord_DI(iDim,iImages)
             case default
                call CON_stop(&
                     'Error in interpolation_fix_reschange')
             end select
          end do DIM
          call bilinear_interpolation(&
               nDim,&
               XyzMisc_D,&
               GridDescriptor,&
               nIndexes,& 
               IndexAux_II,&
               nImagesAux,WeightAux_I)
          if(nImages+nImagesAux>ubound(Index_II,2))&
               call CON_stop(&
               'Algorithmic error in interpolation_fix_reschange')
          Index_II(:,nImages+1:nImages+nImagesAux)=&
               IndexAux_II(:,1:nImagesAux)
          Weight_I(nImages+1:nImages+nImagesAux)=&
               WeightAux_I(1:nImagesAux)*WeightLeft
          nImages=nImages+nImagesAux
       end do IMAGE
       Weight_I(1:nImages)= Weight_I(1:nImages)/&
            sum(Weight_I(1:nImages))
    end if
  end subroutine interpolation_fix_reschange
  !==============================================================!
  integer function i_dir_of_only_reschange(&
         nDim,&
         iShift_D,&
         nImages,&
         lLevel_I)
    integer,intent(in):: nDim,nImages
    integer,intent(in),dimension(nDim)::iShift_D
    integer,intent(in),dimension(nImages)::lLevel_I
    integer::iDir,iPattern,nPatternLength

    i_dir_of_only_reschange=0
    if(nImages==1)return
    if(all(lLevel_I(1:nImages)==0))return
    i_dir_of_only_reschange=None_

    DIRLOOP:do iDir=1,nDim
       if(iShift_D(iDir)==0)CYCLE DIRLOOP
       if(.not.all(lLevel_I(1:iShift_D(iDir))==lLevel_I(1)))&
            CYCLE DIRLOOP
       nPatternLength=iShift_D(iDir)*2
       if(.not.all(lLevel_I(iShift_D(iDir)+1:nPatternLength)&
            ==lLevel_I(iShift_D(iDir)+1)))CYCLE DIRLOOP
       do iPattern=1,nImages/nPatternLength-1
          if(.not.all(lLevel_I(1:nPatternLength)==&
               lLevel_I(iPattern*nPatternLength+1:&
               (iPattern+1)*nPatternLength)))CYCLE DIRLOOP
       end do
       i_dir_of_only_reschange=iDir
       return
    end do DIRLOOP
  end function i_dir_of_only_reschange

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
