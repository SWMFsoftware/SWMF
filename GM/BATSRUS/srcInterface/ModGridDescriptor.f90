!MHD grid in BATSRUS 
Module MH_domain_decomposition
  use CON_grid_storage
  use ModOctree,ONLY:adaptive_block_ptr,adaptive_block,iChildOrder_II
  use ModCube,ONLY:iShiftChild_DI
  implicit none
  logical::UseMHGridDescriptor=.true. 
  integer,private::iLastGrid=-1,iLastDecomposition=-1
  
  integer,parameter::&
       PELast_      =5,&
       LEV_         =6,&
       LEVmin_      =7,&
       LEVmax_      =8
  type(DomainDecompositionType),save::MH_DomainDecomposition
  interface MH_get_root_decomposition
     module procedure MH_get_roots_id
     module procedure MH_get_roots_dd
  end interface
  private::pack_octree,pack_soln_block
contains
  !==============================================================
  subroutine pack_octree(DomainDecomposition)
    use ModOctree,ONLY:octree_roots                                          
    use ModParallel,ONLY:proc_dims 
    type(DomainDecompositionType),intent(inout)::DomainDecomposition
    logical::DoCountOnly
    integer :: i, j, k, l,iRoot
    type (adaptive_block_ptr) :: Octree
    !---------------------------------------------------------------
    DoCountOnly=.true. !To enter the loop
    !Check allocation 
    do while (DoCountOnly)
       call check_octree_grid_allocation(DomainDecomposition)
       DoCountOnly=.false.
       l=0;iRoot=0
       DomainDecomposition%nTreeNodes=0
       do k = 1, proc_dims(3)
          do j = 1, proc_dims(2)
             do i = 1, proc_dims(1)
                Octree % ptr => octree_roots(i, j, k) % ptr
                iRoot=iRoot+1
                call pack_soln_block&
                     (octree,l,iRoot,0,DoCountOnly,DomainDecomposition)
             end do
          end do
       end do
    end do
  end subroutine pack_octree
  !==================================================================  
  recursive subroutine pack_soln_block&
       (octree,l,lParent,iChildNumber,DoCountOnly,DomainDecomposition)
    implicit none
    type (adaptive_block_ptr),intent(inout) :: octree
    integer,intent(inout) :: l
    integer,intent(in) :: lParent,iChildNumber
    logical,intent(inout) :: DoCountOnly
    type(DomainDecompositionType),intent(inout)::DomainDecomposition
    logical IsUsedBlock    
    integer ::lOctree,iCube

    type (adaptive_block_ptr) :: child

    if (associated(octree % ptr)) then
       l=l+1
       lOctree=l
       DomainDecomposition%nTreeNodes=l
       DoCountOnly=DoCountOnly.or.(lOctree>DomainDecomposition%nAllocatedNodes)
       if(.not.DoCountOnly)then
          DomainDecomposition%iDecomposition_II(Parent_,lOctree)=lParent
          DomainDecomposition%iDecomposition_II(MyNumberAsAChild_,lOctree)=iChildNumber
       end if
       IsUsedBlock = octree % ptr % used
       if (.not.IsUsedBlock) then
          do iCube = 1, 8
             if(.not.DoCountOnly)DomainDecomposition%iDecomposition_II(&
                  iChildOrder_II(iCube,iChildNumber), lOctree)=l+1
             child % ptr => &
                  octree % ptr % child(iChildOrder_II(iCube,iChildNumber))%ptr
             call pack_soln_block(child,l,lOctree,&
                  iChildOrder_II(iCube,iChildNumber),DoCountOnly,DomainDecomposition)
          end do
       else
          if(.not.DoCountOnly)then
             DomainDecomposition% iDecomposition_II(FirstChild_,lOctree)=None_
             DomainDecomposition%iDecomposition_II(GlobalBlock_,lOctree) = &
                  octree % ptr % number
             DomainDecomposition%iDecomposition_II(PE_,lOctree) = &
                  octree % ptr % PE
             DomainDecomposition%iDecomposition_II(PELast_,lOctree)=&
                  DomainDecomposition%iDecomposition_II(PE_,lOctree)
             DomainDecomposition%iDecomposition_II(BLK_,lOctree) = &
                  octree % ptr % BLK
             DomainDecomposition%iDecomposition_II(LEV_,lOctree) = &
                  octree % ptr % LEV
             DomainDecomposition%iDecomposition_II(LEVmin_,lOctree) = &
                  octree % ptr % LEVmin
             DomainDecomposition%iDecomposition_II(LEVmax_,lOctree) = &
                  octree % ptr % LEVmax 
          end if
       end if
    end if
  end subroutine pack_soln_block
  !=====================================================================!
  subroutine MH_get_roots_dd(DomainDecomposition)                         
    use ModSize,ONLY:nDim,nIJK_D                                             
    use ModParallel,ONLY:periodic3d,proc_dims                                
    use ModGeometry,ONLY:XyzMin_D,XyzMax_D                                     
    type(DomainDecompositionType),intent(inout)::DomainDecomposition  
    call get_root_decomposition_dd(&
         DomainDecomposition,&!Decomposition to be constructed
         proc_dims  ,&!As in DomainDecompositionType
         XyzMin_D,&   !As in DomainDecompositionType
         XyzMax_D,&   !As in DomainDecompositionType
         nIJK_D,&     !As in DomainDecompositionType
         IsPeriodic_D=periodic3d,& 
         iShift_DI=iShiftChild_DI)       !As in DomainDecompositionType
  end subroutine MH_get_roots_dd                                             
!----------------------------------------------------------------------------!
  subroutine MH_get_roots_id(GridID_)                         
    use ModSize,ONLY:nDim,nIJK_D                                             
    use ModParallel,ONLY:periodic3d,proc_dims                                
    use ModGeometry,ONLY:XyzMin_D,XyzMax_D                                     
    integer,intent(in)::GridID_  
    call get_root_decomposition_id(&
         GridID_,&!Decomposition to be constructed
         proc_dims  ,&!As in DomainDecompositionType
         XyzMin_D,&   !As in DomainDecompositionType
         XyzMax_D,&   !As in DomainDecompositionType
         nIJK_D,&     !As in DomainDecompositionType
         IsPeriodic_D=periodic3d,& 
         iShift_DI=iShiftChild_DI)       !As in DomainDecompositionType
  end subroutine MH_get_roots_id
!----------------------------------------------------------------------------!
  !==========================================================================!
  subroutine MH_update_local_decomposition(DomainDecomposition)                             
    use ModMain,ONLY:iNewGrid,iNewDecomposition 
    type(DomainDecompositionType),intent(inout)::&
         DomainDecomposition
    if(iNewGrid==iLastGrid.and.&
         iNewDecomposition==iLastDecomposition&
    .and.DomainDecomposition%iRealization/=0)return
    call pack_octree(DomainDecomposition)                  
    DomainDecomposition%iRealization=&
         mod(DomainDecomposition%iRealization+1,1000)   
    iLastDecomposition=iNewDecomposition
    iLastGrid=iNewGrid
    call complete_grid(DomainDecomposition)
  end subroutine MH_update_local_decomposition                                        
!----------------------------------------------------------------------------!
end Module MH_domain_decomposition
