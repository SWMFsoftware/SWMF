!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module CON_grid_storage
  use CON_domain_decomposition
  use ModInitGridStorage,ONLY:init_grid_storage, MaxGrid
  implicit none

  !BOP
  !DESCRIPTION:
  ! The resulted domain decompositions should be properly         
  ! registered and should obtain the unique GridID\_               
  !EOP
  type(DDPointerType),dimension(MaxGrid),&
       private::&
       DD_I
  logical,dimension(MaxGrid)::DoneDDInit_C=.false.

  !Introduced to bypass the HALEM compiler restrictions
  integer,dimension(MaxGrid)::nDim_C=0
  !BOP
  !DESCRIPTION:
  !INTERFACES
  !Each of the module procedures differ with the FIRST dummy     
  !parameter, which is either GridID\_, or DomainDecompositionType 
  !structure. The procedures differ with suffix id or dd from the 
  !generic name 
  !EOP
  interface clean_decomposition
     module procedure clean_decomposition_id
     module procedure clean_decomposition_dd
  end interface

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

  interface bcast_dd_indexes
     module procedure bcast_indexes_id
     module procedure bcast_indexes_dd
  end interface

  interface synchronize_refinement
     module procedure synchronize_refinement_id
     module procedure synchronize_refinement_dd
  end interface

  interface is_left_boundary_d
     module procedure is_left_boundary_id
     module procedure is_left_boundary_dd
  end interface

  interface is_right_boundary_d
     module procedure is_right_boundary_id
     module procedure is_right_boundary_dd
  end interface

  interface xyz_cell_d
     module procedure xyz_cell_id
     module procedure xyz_cell_dd
  end interface

  interface l_neighbor
     module procedure l_neighbor_id
     module procedure l_neighbor_dd
  end interface

  interface l_level_neighbor 
     module procedure l_level_neighbor_id
     module procedure l_level_neighbor_dd
  end interface

  interface search_in
     module procedure search_in_id
     module procedure search_in_dd
  end interface

  interface pe_decomposition
     module procedure pe_id
     module procedure pe_dd
  end interface

  interface blk_decomposition
     module procedure blk_id
     module procedure blk_dd
  end interface

  interface pe_and_blk
     module procedure pe_and_blk_id
     module procedure pe_and_blk_dd
  end interface

  interface n_block
     module procedure n_block_id
     module procedure n_block_dd
  end interface

  interface n_block_total
     module procedure n_block_total_id
     module procedure n_block_total_dd
  end interface

  interface min_block_pe
     module procedure min_block_pe_id
     module procedure min_block_pe_dd
  end interface

  interface max_block_pe
     module procedure max_block_pe_id
     module procedure max_block_pe_dd
  end interface

  interface min_block
     module procedure min_block_id
     module procedure min_block_dd
  end interface

  interface max_block
     module procedure max_block_id
     module procedure max_block_dd
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

  interface used_bp
     module procedure used_bp_dd
     module procedure used_bp_id
  end interface

  interface compid_grid
     module procedure compid_grid_id
     module procedure compid_grid_dd
  end interface

  interface ndim_grid
     module procedure ndim_grid_id
     module procedure ndim_grid_dd
  end interface

  interface xyz_min_d
     module procedure xyz_min_grid_id
     module procedure xyz_min_grid_dd
  end interface

  interface xyz_max_d
     module procedure xyz_max_grid_id
     module procedure xyz_max_grid_dd
  end interface

  interface root_map_d
     module procedure root_map_id
     module procedure root_map_dd
  end interface

  interface ncells_decomposition_d
     module procedure ncells_grid_id
     module procedure ncells_grid_dd
  end interface

  interface ntree_nodes
     module procedure ntree_nodes_id
     module procedure ntree_nodes_dd
  end interface

  interface xyz_block_d
     module procedure xyz_block_id
     module procedure xyz_block_dd
  end interface

  interface d_xyz_block_d
     module procedure d_xyz_block_id
     module procedure d_xyz_block_dd
  end interface

  interface d_xyz_cell_d
     module procedure d_xyz_cell_id
     module procedure d_xyz_cell_dd
  end interface

  interface used_node
     module procedure is_used_block_id
     module procedure is_used_block_dd
  end interface

  interface is_root_node
     module procedure is_root_node_id
     module procedure is_root_node_dd
  end interface

  interface i_realization
     module procedure irealization_id
     module procedure irealization_dd
  end interface

  interface is_local_grid
     module procedure is_local_grid_id
     module procedure is_local_grid_dd
  end interface

  interface is_tree
     module procedure is_tree_dd
     module procedure is_tree_id
  end interface

  interface associate_dd_pointer
     module procedure associate_dd_pointer_id
     module procedure associate_dd_pointer_dd
  end interface
contains
  !BOP
  !IROUTINE: done_dd_init - if the decomposition is initialized
  !INTERFACE:
  logical function done_dd_init(GridID_)
    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_
    !DESCRIPTION:
    !Returns .true. if the domain decomposition is initialaized
    !EOP
    done_dd_init=.false.
    if(GridID_>0.and.GridID_<=MaxGrid)&
         done_dd_init=DoneDDInit_C(GridID_)
  end function done_dd_init
  !===============================================================!
  !===============================================================!
  !Cleaner for domain decomposition                               !
  !BOP
  !IROUTINE: clean_decomposition 
  !INTERFACE:
  subroutine clean_decomposition_id(GridID_)

    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_
    !EOP
    if(done_dd_init(GridID_))then
       call clean_decomposition_dd(DD_I(GridID_)%Ptr)
       DoneDDInit_C(GridID_)=.false.
    end if
  end subroutine clean_decomposition_id
  !---------------------------------------------------------------!
  !===============================================================
  !BOP
  !IROUTINE: init_decomposition - initialize domain decomposition
  !DESCRIPTION:
  !Initialization for the decomposition.       
  !Note that if the decomposition is a tree decomposition and     
  !another order of children is accepted in the component,        
  !different from that assumed to be standard one                  
  !the array iShift\_DI should be reassigned in set\_root.
  !EOP
  !BOP          
  !INTERFACE:
  subroutine init_decomposition_id(&
       GridID_,&
       CompID_,&
       nDim,   &
       IsTreeDecomposition,&
       nDimTree)

    !INPUT ARGUMENTS:         
    integer,intent(in)::GridID_           
    integer,intent(in)::CompID_,nDim
    logical,intent(in),optional::IsTreeDecomposition
    integer,intent(in),optional::nDimTree 
    !EOP
    !---------------------------------------------------------------!
    if(GridID_>MaxGrid.or.GridID_<=0)call CON_stop&
         ('GridID_ is not allowed to be equal to',GridID_)
    if(done_dd_init(GridID_))call CON_stop(&
         'An attempt to reinitialize the gridID_',GridID_)
    call init_grid_storage(DD_I,GridID_) 
    if(present(nDimTree))then
       call init_decomposition_dd(&
            DD_I(GridID_)%Ptr,&
            CompID_,&
            nDim,   &
            IsTreeDecomposition,&
            nDimTree)
    elseif(present(IsTreeDecomposition))then
       call init_decomposition_dd(&
            DD_I(GridID_)%Ptr,&
            CompID_,&
            nDim,   &
            IsTreeDecomposition)
    else     
       call init_decomposition_dd(&
            DD_I(GridID_)%Ptr,&
            CompID_,&
            nDim)
    end if
    DoneDDInit_C(GridID_)=.true.
    nDim_C(GridID_)=nDim
  end subroutine init_decomposition_id
  !---------------------------------------------------------------!
  !==========================WITH INTERFACE=======================!
  !BOP
  !IROUTINE: get_root_decomposition - get the domain decomposition
  !DESCRIPTION:
  !\begin{verbatim}
  !===============================================================!
  !To get a decomposition domain, even the tree one, the root     !
  !decomposition should be first constructed                      !
  !====================ATTENTION==================================!
  !PE here are the ranks in the LOCAL communicator for the        !
  !component  
  !\end{verbatim}
  !EOP
  !BOP
  !INTERFACE:                                                    !
  subroutine get_root_decomposition_id(&
       GridID_,&!     !ID for Decomposition to be constructed   
       iRootMapDim_D,&!As in DomainDecompositionType            
       XyzMin_D,&     !As in DomainDecompositionType
       XyzMax_D,&     !As in DomainDecompositionType
       nCells_D,&     !As in DomainDecompositionType
       PE_I,&         !PE layout
       iBlock_I,&     !Local Block Number layout
       IsPeriodic_D,& !As in DomainDecompositionType
       iShift_DI)     !As in DomainDecompositionType

    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_
    integer,dimension(:),&
         intent(in)::iRootMapDim_D
    real,dimension(:),&
         intent(in)::XyzMin_D   
    real,dimension(:),&
         intent(in)::XyzMax_D
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
    !EOP
    integer::lBlock,MaxBlock
    !---------------------------------------------------------------!
    !Check the dimension of PE_I and iBlock_I                       !
    if(present(PE_I))then
       if(ubound(PE_I,1)/=product(iRootMapDim_D))&
            call CON_stop(&
            'The dimension of PE_I should be equal to ',&
            product(iRootMapDim_D))
    end if

    if(present(iBlock_I))then
       if(ubound(iBlock_I,1)/=product(iRootMapDim_D))&
            call CON_stop(&
            'The dimension of iBlock_I should be equal to ',&
            product(iRootMapDim_D))
    end if

    !---------------------------------------------------------------!
    !Assign the grid parameters                                     !

    DD_I(GridID_)%Ptr%iRootMapDim_D=iRootMapDim_D

    DD_I(GridID_)%Ptr%XyzMin_D=XyzMin_D
    DD_I(GridID_)%Ptr%XyzMax_D=XyzMax_D
    DD_I(GridID_)%Ptr%nCells_D=nCells_D

    if(present(IsPeriodic_D))&
         DD_I(GridID_)%Ptr%IsPeriodic_D=IsPeriodic_D

    DD_I(GridID_)%Ptr%nTreeNodes=product(iRootMapDim_D)
    call check_octree_grid_allocation(DD_I(GridID_)%Ptr)


    DD_I(GridID_)%Ptr%iDecomposition_II&
         (MyNumberAsAChild_,1:DD_I(GridID_)%Ptr%nTreeNodes)=0    
    DD_I(GridID_)%Ptr%iDecomposition_II&
         (FirstChild_,1:DD_I(GridID_)%Ptr%nTreeNodes)=None_
    do lBlock=1,DD_I(GridID_)%Ptr%nTreeNodes
       DD_I(GridID_)%Ptr%iDecomposition_II&
            (Parent_,1:DD_I(GridID_)%Ptr%nTreeNodes)=lBlock
    end do
    if(DD_I(GridID_)%Ptr%IsTreeDecomposition)then
       call check_iroot_allocation(DD_I(GridID_)%Ptr)
       if(present(iShift_DI))&
            DD_I(GridID_)%Ptr%iShift_DI=iShift_DI
    end if

    !If neither PE_I nor iBlock_I is present, the root              !
    !decomposition blocks are balanced for an optimal load          !
    if((.not.present(PE_I)).and.(.not.present(iBlock_I)))then
       MaxBlock=(DD_I(GridID_)%Ptr%nTreeNodes-1)/&
            n_proc(DD_I(GridID_)%Ptr%CompID_)+1
       do lBlock=1,DD_I(GridID_)%Ptr%nTreeNodes
          DD_I(GridID_)%Ptr%iDecomposition_II&
               (BLK_,lBlock)=mod(lBlock,MaxBlock)+1
          DD_I(GridID_)%Ptr%iDecomposition_II&
               (PE_,lBlock)=(lBlock-1)/MaxBlock
       end do
       call set_iglobal_and_bp_dd(DD_I(GridID_)%Ptr)
       return
    end if

    !If the PE\_I is not given, it is assumed that all the blocks    
    !are at the root PE                                             

    if(present(PE_I))then
       DD_I(GridID_)%Ptr%iDecomposition_II&
            (PE_,1:DD_I(GridID_)%Ptr%nTreeNodes)=PE_I
    else
       DD_I(GridID_)%Ptr%iDecomposition_II&
            (PE_,1:DD_I(GridID_)%Ptr%nTreeNodes)=0
    end if

    !If the iBlock_I is not given, they are assumed to be enumerated!
    !like in the following loop:                                    !
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
       DD_I(GridID_)%Ptr%iDecomposition_II&
            (BLK_,1:DD_I(GridID_)%Ptr%nTreeNodes)=iBlock_I
    else
       do lBlock=1, DD_I(GridID_)%Ptr%nTreeNodes
          DD_I(GridID_)%Ptr%iDecomposition_II&
               (BLK_,lBlock)=lBlock
       end do
    end if
    call set_iglobal_and_bp_dd(DD_I(GridID_)%Ptr)
  end subroutine get_root_decomposition_id
  !---------------------------------------------------------------!
  !BOP
  !INTERFACE:
  subroutine bcast_decomposition_id(&
       GridID_)

    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_
    !EOP
    call bcast_decomposition_dd(DD_I(GridID_)%Ptr)
  end subroutine bcast_decomposition_id
  !---------------------------------------------------------------!
  !IROUTINE: bcast_indexes - send the indexes of the decomposition
  !DESCRIPTION:
  !Broadcasts only Decomposition\_II. This is a part of the
  !synchronize\_refinement procedure, but it can be used separately 
  !too. If any of the optional parameters is not present, the 
  !Decomposition\_II array is sent from the root processor of the 
  !components to all PEs in the global communicator, otherwise
  !it itis sent FROM the PE having the rank iProcUnion in the 
  !communicator iCommUnion TO all PEs of this communicator.  
  !Recalculate local PE ranks to their values in the global       
  !communicator while broadcasting.
  !INTERFACE:
  subroutine bcast_indexes_id(&
       GridID_,iProcUnion,iCommUnion)

    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_
    integer,intent(in),optional::iProcUnion,iCommUnion 

    if(present(iProcUnion).and.present(iCommUnion))then
       call bcast_indexes_dd(&
            DD_I(GridID_)%Ptr,iProcUnion,iCommUnion)
    else
       call bcast_indexes_dd(&
            DD_I(GridID_)%Ptr)
    end if
  end subroutine bcast_indexes_id
  !---------------------------------------------------------------!
  !===============================================================!
  !BOP
  !IROUTINE: synchronize_refinement - update decomposition after AMR
  !DESCRIPTION:
  !\begin{verbatim}
  !           SYNCHRONIZE LOCAL AND GLOBAL GRID                  
  !           NOTE: IF GridID\_ is used for global grid, then      
  !           synchronize\_refinement is the only way to properly  
  !           account for the refinement
  !
  !If any of the optional parameters is not present, the global 
  !decomposition at all the PEs of the global communicator is 
  !synchronized with the local one at root processor of the 
  !component. Otherwise the global 
  !decomposition at all the PEs of the communicator iProcUnion is 
  !synchronized with the local one at the PE having the rank 
  !iProcUnion in the communicator iCommUnion.
  !Recalculate local PE ranks (of the local grid) to their values 
  !in the global communicator (i\_comm()).
  !\end{verbatim}                          
  !EOP
  !BOP
  !INTERFACE:
  subroutine synchronize_refinement_id(&
       GridID_,LocalDD,iProcUnion,iCommUnion)

    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_
    type(DomainDecompositionType),&   !Local Grid with which the   
         intent(in)::LocalDD          !global one is synchronized
    integer,intent(in),optional::iProcUnion,iCommUnion
    !EOP
    if(present(iProcUnion).and.present(iCommUnion))then
       call synchronize_refinement_dd(&
            DD_I(GridID_)%Ptr,LocalDD,iProcUnion,iCommUnion)
    else
       call synchronize_refinement_dd(&
            DD_I(GridID_)%Ptr,LocalDD)
    end if
  end subroutine synchronize_refinement_id
  !---------------------------------------------------------------!
  !BOP
  !INTERFACE:
  function is_left_boundary_id(&
       GridID_,lGlobalTreeNumber)
    !INPUT ARGUMENTS:
    integer,intent(in)::lGlobalTreeNumber
    integer,intent(in)::GridID_
    !OUTPUT ARGUMENTS:
    logical,dimension(nDim_C(GridID_))::&
         is_left_boundary_id
    !EOP
    is_left_boundary_id=is_left_boundary_dd(&
         DD_I(GridID_)%Ptr,lGlobalTreeNumber)
  end function is_left_boundary_id
  !---------------------------------------------------------------!
  !BOP
  !INTERFACE:
  function is_right_boundary_id(&
       GridID_,lGlobalTreeNumber)
    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_,lGlobalTreeNumber
    logical,dimension(nDim_C(GridID_))::&
         is_right_boundary_id
    !EOP
    is_right_boundary_id=&
         is_right_boundary_dd(&
         DD_I(GridID_)%Ptr,lGlobalTreeNumber)
  end function is_right_boundary_id
  !---------------------------------------------------------------!
  !BOP
  !INTERFACE:
  function xyz_cell_id(&
       GridID_,lGlobalTreeNumber,iCells_D)
    !INPUT ARGUMENTS:
    integer,intent(in):: GridID_,lGlobalTreeNumber
    integer,dimension(nDim_C(GridID_)),&
         intent(in)::iCells_D
    !EOP
    real,dimension(nDim_C(GridID_))::xyz_cell_id
    xyz_cell_id= xyz_cell_dd(&
         DD_I(GridID_)%Ptr,lGlobalTreeNumber,iCells_D)
  end function xyz_cell_id
  !---------------------------------------------------------------!
  !BOP
  !INTERFACE:                                                       
  integer function l_neighbor_id(&
       GridID_,lGlobalTreeNumber,iCells_D)
    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_,lGlobalTreeNumber
    integer,dimension(nDim_C(GridID_)),&
         intent(in)::iCells_D
    !EOP
    l_neighbor_id=l_neighbor_dd(&
         DD_I(GridID_)%Ptr,lGlobalTreeNumber,iCells_D)
  end function l_neighbor_id
  !---------------------------------------------------------------!
  !INTERFACE:                                      
  integer function l_level_neighbor_id(GridID_,&
       lGlobalTreeNumber,iCells_D)
    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_,lGlobalTreeNumber
    integer,dimension(nDim_C(GridID_)),&
         intent(in)::iCells_D
    !EOP
    l_level_neighbor_id=l_level_neighbor_DD(&
         DD_I(GridID_)%Ptr,lGlobalTreeNumber,iCells_D)
  end function l_level_neighbor_id
  !---------------------------------------------------------------!
  !BOP
  !INTERFACE:
  subroutine search_in_id(GridID_,Xyz_D,lGlobalTreeNumber)
    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_
    real,dimension(nDim_C(GridID_)),&
         intent(inout)::Xyz_D
    !OUTPUT ARGUMENTS:
    integer,intent(out)::lGlobalTreeNumber
    !EOP
    call search_in_dd(DD_I(GridID_)%Ptr,Xyz_D,lGlobalTreeNumber)
  end subroutine search_in_id
  !---------------------------------------------------------------!
  !BOP
  !INTERFACE:
  subroutine search_cell_id(&
       GridID_,lGlobalTreeNumber,Xyz_D,iCells_D)
    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_,lGlobalTreeNumber
    real,dimension(nDim_C(GridID_)),&
         intent(inout)::Xyz_D
    !OUTPUT ARGUMENTS:
    integer,dimension(nDim_C(GridID_)),&
         intent(out)::iCells_D
    !EOP
    call search_cell_dd(DD_I(GridID_)%Ptr,&
         lGlobalTreeNumber,Xyz_D,iCells_D)
  end subroutine search_cell_id
  !---------------------------------------------------------------!
  !BOP
  !INTERFACE:
  integer function pe_id(&
       GridID_,lGlobalTreeNumber)
    !INPUT ARGUMENTS: 
    integer,intent(in)::GridID_,lGlobalTreeNumber
    !EOP
    pe_id=DD_I(GridID_)%Ptr%iDecomposition_II(&
         PE_,lGlobalTreeNumber)
  end function  pe_id
  !---------------------------------------------------------------!

  !BOP
  !INTERFACE:
  integer function blk_id(&
       GridID_,lGlobalTreeNumber)
    !INPUT ARGUMENTS: 
    integer,intent(in)::GridID_,lGlobalTreeNumber
    !EOP
    blk_id=DD_I(GridID_)%Ptr%iDecomposition_II(&
         BLK_,lGlobalTreeNumber)
  end function blk_id
  !---------------------------------------------------------------!
  !BOP
  !INTERFACE:
  subroutine pe_and_blk_id(&
       GridID_,lGlobalTreeNumber,iPEOut,iBlockOut) 
    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_,lGlobalTreeNumber
    !OUTPUT ARGUMENTS:
    integer,intent(out)::iPEOut,iBlockOut
    !EOP
    call pe_and_blk_dd(DD_I(GridID_)%Ptr,&
         lGlobalTreeNumber,iPEOut,iBlockOut) 
  end subroutine pe_and_blk_id
  !---------------------------------------------------------------!
  !INTERFACE:
  integer function n_block_id(GridID_,iPE)
    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_,iPE
    !EOP
    n_block_id=n_block_dd(DD_I(GridID_)%Ptr,iPE)
  end function n_block_id
  !---------------------------------------------------------------!
  !BOP
  !INTERFACE:
  integer function n_block_total_id(GridID_)
    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_
    !EOP
    n_block_total_id=n_block_total_dd(DD_I(GridID_)%Ptr)
  end function n_block_total_id
  !---------------------------------------------------------------!
  !BOP
  !INTERFACE:
  integer function min_block_pe_id(GridID_,iPE) 
    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_,iPE
    !EOP
    min_block_pe_id=min_block_pe_dd(DD_I(GridID_)%Ptr,iPE)
  end function min_block_pe_id
  !---------------------------------------------------------------!
  !BOP
  !INTERFACE:
  integer function max_block_pe_id(GridID_,iPE) 
    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_,iPE
    !EOP
    max_block_pe_id=max_block_pe_dd(DD_I(GridID_)%Ptr,iPE)
  end function max_block_pe_id
  !---------------------------------------------------------------!
  !BOP
  !INTERFACE:
  integer function min_block_id(GridID_) 
    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_
    !EOP
    min_block_id=min_block_dd(DD_I(GridID_)%Ptr) 
  end function min_block_id
  !---------------------------------------------------------------!
  !BOP
  !INTERFACE:
  integer function max_block_id(GridID_)
    !INPUT ARGUMENTS: 
    integer,intent(in)::GridID_
    !EOP
    max_block_id=max_block_dd(DD_I(GridID_)%Ptr)
  end function max_block_id
!---------------------------------------------------------------!
  integer function iglobal_bp_id(GridID_,iBLK,iPE)
    integer,intent(in)::GridID_
    integer,intent(in)::iBLK,iPE
    iglobal_bp_id=DD_I(GridID_)%Ptr%iGlobal_BP(iBLK,iPE)
  end function iglobal_bp_id
  !---------------------------------------------------------------!
  integer function iglobal_node_id(GridID_,iBlockAll)
    integer,intent(in)::GridID_
    integer,intent(in)::iBlockAll
    iglobal_node_id=DD_I(GridID_)%Ptr%iGlobal_A(iBlockAll)
  end function iglobal_node_id
  !---------------------------------------------------------------!
  integer function iglobal_block_id(&
       GridID_,iGlobalTreeNode)
    integer,intent(in)::GridID_
    integer,intent(in)::iGlobalTreeNode
    iglobal_block_id=DD_I(GridID_)%Ptr%iDecomposition_II(&
         GlobalBlock_,iGlobalTreeNode)
  end function iglobal_block_id
  !---------------------------------------------------------------!
  logical function used_bp_id(GridID_,iBLK,iPE)
    integer,intent(in)::GridID_
    integer,intent(in)::iBLK,iPE
    used_bp_id=used_bp_dd(DD_I(GridID_)%Ptr,iBLK,iPE)
  end function used_bp_id
  !---------------------------------------------------------------!
  !BOP
  !INTERFACE:
  integer function compid_grid_id(GridID_)
    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_
    !EOP
    compid_grid_id=compid_grid_dd(DD_I(GridID_)%Ptr)
  end function compid_grid_id
  !---------------------------------------------------------------!
  !===============================================================!
  !BOP
  !IROUTINE: access to decomposition elements: ndim_grid 
  !INTERFACE:
  integer function ndim_grid_id(GridID_)
    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_
    !EOP
    ndim_grid_id=ndim_grid_dd(DD_I(GridID_)%Ptr)
  end function ndim_grid_id
  !---------------------------------------------------------------!
  !INTERFACE:
  function xyz_min_grid_id(GridID_)
    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_
    !OUTPUT ARGUMENTS:
    real,dimension(nDim_C(GridID_))::xyz_min_grid_id
    !EOP
    xyz_min_grid_id=xyz_min_grid_dd(DD_I(GridID_)%Ptr)
  end function xyz_min_grid_id
  !---------------------------------------------------------------!
  !INTERFACE:
  function xyz_max_grid_id(GridID_)
    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_
    !OUTPUT ARGUMENTS:
    real,dimension(nDim_C(GridID_))::xyz_max_grid_id
    !EOP
    xyz_max_grid_id=xyz_max_grid_dd(DD_I(GridID_)%Ptr)
  end function xyz_max_grid_id
  !---------------------------------------------------------------!
  !INTERFACE:
  function root_map_id(GridID_)
    !INPUT AGRUMENTS:
    integer,intent(in)::GridID_
    !EOP
    integer,dimension(nDim_C(GridID_))::root_map_id
    root_map_id=root_map_dd(DD_I(GridID_)%Ptr)
  end function root_map_id
  !---------------------------------------------------------------!
  !INTERFACE:
  function ncells_grid_id(GridID_)
    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_
    !OUTPUT ARGUMENTS:
    integer,dimension(nDim_C(GridID_))::ncells_grid_id
    !EOP
    ncells_grid_id=ncells_grid_dd(DD_I(GridID_)%Ptr)
  end function ncells_grid_id
  !---------------------------------------------------------------!
  !INTERFACE:
  integer function ntree_nodes_id(GridID_)
    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_
    !EOP
    ntree_nodes_id=ntree_nodes_dd(DD_I(GridID_)%Ptr)
  end function ntree_nodes_id
  !---------------------------------------------------------------!
  !INTERFACE:
  function xyz_block_id(GridID_,lGlobalTreeNumber)
    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_,lGlobalTreeNumber
    !OUTPUT ARGUMENTS:
    real,dimension(nDim_C(GridID_)):: xyz_block_id
    !EOP
    xyz_block_id=xyz_block_dd(&
         DD_I(GridID_)%Ptr,lGlobalTreeNumber)
  end function xyz_block_id
  !---------------------------------------------------------------!
  !INTERFACE:
  function d_xyz_block_id(GridID_,lGlobalTreeNumber)
    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_,lGlobalTreeNumber
    !OUTPUT ARGUMENTS:
    real,dimension(nDim_C(GridID_)):: d_xyz_block_id
    !EOP
    d_xyz_block_id=d_xyz_block_dd(&
         DD_I(GridID_)%Ptr,lGlobalTreeNumber)
  end function d_xyz_block_id
  !---------------------------------------------------------------!
  !INTERFACE:
  function d_xyz_cell_id(GridID_,lGlobalTreeNumber)
    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_,lGlobalTreeNumber
    !OUTPUT ARGUMENTS:
    real,dimension(nDim_C(GridID_)):: d_xyz_cell_id
    !EOP
    d_xyz_cell_id=d_xyz_cell_dd(&
         DD_I(GridID_)%Ptr,lGlobalTreeNumber)
  end function d_xyz_cell_id
  !---------------------------------------------------------------!
  !INTERFACE:
  logical function is_used_block_id(GridID_,lGlobalTreeNumber)
    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_,lGlobalTreeNumber
    !EOP
    is_used_block_id=is_used_block_dd(DD_I(GridID_)%Ptr,&
         lGlobalTreeNumber)
  end function is_used_block_id
  !---------------------------------------------------------------!
  logical function is_root_node_id(GridID_,lGlobalTreeNumber)
    integer,intent(in)::GridID_,lGlobalTreeNumber
    is_root_node_id=is_root_node_dd(DD_I(GridID_)%Ptr,&
         lGlobalTreeNumber)
  end function is_root_node_id
  !---------------------------------------------------------------!
  !INTERFACE:
  integer function irealization_id(GridID_)
    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_
    !EOP
    irealization_id=irealization_dd(DD_I(GridID_)%Ptr)
  end function irealization_id
  !---------------------------------------------------------------!
  logical function is_local_grid_id(GridID_)
    integer,intent(in)::GridID_
    is_local_grid_id=.false.
  end function is_local_grid_id
!---------------------------------------------------------------!
  logical function is_tree_id(GridID_)
    integer,intent(in)::GridID_
    is_tree_id=is_tree_dd(DD_I(GridID_)%Ptr)
  end function is_tree_id
  !---------------------------------------------------------------!
  subroutine associate_dd_pointer_id(GridID_,DDPointer)
    integer,intent(in)::GridID_
    type(DDPointerType),intent(out)::DDPointer
    nullify(DDPointer%Ptr)
    DDPointer%Ptr=>DD_I(GridID_)%Ptr
  end subroutine associate_dd_pointer_id
  !---------------------------------------------------------------!
end module CON_grid_storage

