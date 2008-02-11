!MHD grid in BATSRUS 
Module MH_domain_decomposition
  use CON_grid_storage
  use ModOctree,ONLY:adaptive_block_ptr,adaptive_block
  implicit none
  logical::UseMHGridDescriptor=.true. 
  integer,private::iLastGrid=-1,iLastDecomposition=-1
  
  integer,parameter::&
       PELast_      =5,&
       LEV_         =6,&
       LEVmin_      =7,&
       LEVmax_      =8

  ! Decide order based on child number
  integer, parameter, dimension(8,0:8):: iChildOrder_II = reshape(&
       (/4,1,8,5,6,7,2,3,&    !0
       4,1,2,3,6,7,8,5,&    !1
       6,7,8,5,4,1,2,3,&    !2
       2,1,8,7,6,5,4,3,&    !3
       4,3,6,5,8,7,2,1,&    !4
       2,1,4,3,6,5,8,7,&    !5
       8,7,6,5,4,3,2,1,&    !6
       4,1,8,5,6,7,2,3,&    !7
       4,1,8,5,6,7,2,3/),&  !8
       (/8,9/))
  integer,parameter,dimension(3,8)::iShift3_DI= reshape(&
                                      !Coords /  1,2,3   /Children 
                                               (/0,0,1,&!1
                                                 1,0,1,&!2
                                                 1,0,0,&!3
                                                 0,0,0,&!4
                                                 0,1,0,&!5
                                                 1,1,0,&!6
                                                 1,1,1,&!7
                                                 0,1,1/),&!8
                                                 (/3,8/))
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
             select case (iChildOrder_II(iCube,iChildNumber))
             case (1)
                child % ptr => octree % ptr % child1
             case (2)
                child % ptr => octree % ptr % child2
             case (3)
                child % ptr => octree % ptr % child3
             case (4)
                child % ptr => octree % ptr % child4
             case (5)
                child % ptr => octree % ptr % child5
             case (6)
                child % ptr => octree % ptr % child6
             case (7)
                child % ptr => octree % ptr % child7
             case (8)
                child % ptr => octree % ptr % child8
             end select
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
  !=========================================================================
  recursive subroutine unpack_node(octree,l,DomainDecomposition)
    use ModOctree,ONLY:global_block_ptrs
    implicit none
    type (adaptive_block_ptr),intent(inout) :: octree
    integer,intent(inout) :: l
    type(DomainDecompositionType),intent(in)::DomainDecomposition
    logical :: IsUsedBlock

    integer ::lOctree,iPE,iBLK

    type (adaptive_block_ptr) :: child

    integer :: iChild,iChildNumber,iLEV,iError

    !---------------------------------------------------------------------------

    if (associated(octree % ptr)) then
       l=l+1
       lOctree=l
       IsUsedBlock=DomainDecomposition%&
            iDecomposition_II(FirstChild_,lOctree)==None_
       octree % ptr % used    = IsUsedBlock
       iChildNumber=DomainDecomposition%iDecomposition_II(MyNumberAsAChild_,lOctree)
       octree % ptr % child_number =iChildNumber
       if(IsUsedBlock)then
          octree % ptr % number  = DomainDecomposition%&
               iDecomposition_II(GlobalBlock_,lOctree)
          octree % ptr % PE      = DomainDecomposition%&
               iDecomposition_II(PE_,lOctree)
          octree % ptr % BLK     = DomainDecomposition%&
               iDecomposition_II(BLK_,lOctree)
          octree % ptr % LEVmin  = DomainDecomposition%&
               iDecomposition_II(LEVmin_,lOctree)
          octree % ptr % LEVmax  = DomainDecomposition%&
               iDecomposition_II(LEVmax_,lOctree)
          global_block_ptrs(octree % ptr % BLK, octree % ptr % PE +1) % ptr &
               =>octree % ptr
       else
          !Choose the PE and BLK set to be used in refining the octree node.
          iLEV=octree % ptr % LEV +1
          octree % ptr % used = .false.
          do iChild = 1, 8
             nullify(child % ptr)
             allocate ( child % ptr, stat=ierror )
             if (ierror > 0) write(*,*) "initialize_octree_block: ", &
                  & " allocation error for octree"
             nullify (child % ptr % parent)
             nullify (child % ptr % child1)
             nullify (child % ptr % child2)
             nullify (child % ptr % child3)
             nullify (child % ptr % child4)
             nullify (child % ptr % child5)
             nullify (child % ptr % child6)
             nullify (child % ptr % child7)
             nullify (child % ptr % child8)

             child % ptr % number  = 0
             child % ptr % child_number = iChildOrder_II(iChild,iChildNumber)
             child % ptr % PE      = -1
             child % ptr % BLK     = -1
             child % ptr % LEV     = iLEV
             child % ptr % LEVmin  =  0
             child % ptr % LEVmax  =  99
             child % ptr % used    = .true.
             child % ptr % refine  = .false.
             child % ptr % coarsen = .false.
             child % ptr % body    = .false.
             child % ptr % IsExtraBoundaryOrPole = .false.
             child % ptr % IsOuterBoundary = .false.

             child % ptr % parent =>octree % ptr
             select case (iChildOrder_II(iChild,iChildNumber))
             case (1)
                octree % ptr % child1 =>  child % ptr
             case (2)
                octree % ptr % child2 =>  child % ptr
             case (3)
                octree % ptr % child3 =>  child % ptr
             case (4)
                octree % ptr % child4 =>  child % ptr
             case (5)
                octree % ptr % child5 =>  child % ptr
             case (6)
                octree % ptr % child6 =>  child % ptr
             case (7)
                octree % ptr % child7 =>  child % ptr
             case (8)
                octree % ptr % child8 =>  child % ptr
             end select
             call unpack_node(child, l, DomainDecomposition)
          end do
       end if
    end if
  end subroutine unpack_node
!============================================================================!

  subroutine MH_get_roots_dd(DomainDecomposition)                         
    use ModSize,ONLY:nDim,nCells                                             
    use ModParallel,ONLY:periodic3d,proc_dims                                
    use ModGeometry,ONLY:XyzMin_D,XyzMax_D                                     
    type(DomainDecompositionType),intent(inout)::DomainDecomposition  
    call get_root_decomposition_dd(&
         DomainDecomposition,&!Decomposition to be constructed
         proc_dims  ,&!As in DomainDecompositionType
         XyzMin_D,&   !As in DomainDecompositionType
         XyzMax_D,&   !As in DomainDecompositionType
         nCells,&     !As in DomainDecompositionType
         IsPeriodic_D=periodic3d,& 
         iShift_DI=iShift3_DI)       !As in DomainDecompositionType
  end subroutine MH_get_roots_dd                                             
!----------------------------------------------------------------------------!
  subroutine MH_get_roots_id(GridID_)                         
    use ModSize,ONLY:nDim,nCells                                             
    use ModParallel,ONLY:periodic3d,proc_dims                                
    use ModGeometry,ONLY:XyzMin_D,XyzMax_D                                     
    integer,intent(in)::GridID_  
    call get_root_decomposition_id(&
         GridID_,&!Decomposition to be constructed
         proc_dims  ,&!As in DomainDecompositionType
         XyzMin_D,&   !As in DomainDecompositionType
         XyzMax_D,&   !As in DomainDecompositionType
         nCells,&     !As in DomainDecompositionType
         IsPeriodic_D=periodic3d,& 
         iShift_DI=iShift3_DI)       !As in DomainDecompositionType
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
