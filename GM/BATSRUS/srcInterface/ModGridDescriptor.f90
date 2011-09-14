!MHD grid in BATSRUS 
module MH_domain_decomposition

  use CON_grid_storage
  use ModMain, ONLY: UseBatl

  implicit none

  SAVE

  interface MH_get_root_decomposition
     module procedure MH_get_roots_id
     module procedure MH_get_roots_dd
  end interface
  public:: MH_get_root_decomposition
  public:: MH_update_local_decomposition

  type(DomainDecompositionType), public:: MH_DomainDecomposition

  ! Local variables and constants
  logical:: UseMHGridDescriptor=.true. 
  integer, private:: iLastGrid = -1, iLastDecomposition = -1
  
  integer, parameter, private::    &
       PELast_      = 5, &
       LEV_         = 6, &
       LEVmin_      = 7, &
       LEVmax_      = 8

  ! Position of children relative to the parent block
  ! in the Morton ordering
  integer, parameter, private:: iShiftMorton_DI(3,8)= reshape( (/ &
       0,0,0, &
       1,0,0, &
       0,1,0, &
       1,1,0, &
       0,0,1, &
       1,0,1, &
       0,1,1, &
       1,1,1 /), (/3,8/))

  private:: show_domain_decomp, get_batl_tree, pack_octree, &
       pack_soln_block, MH_get_roots_dd, MH_get_roots_id

contains

  !====================================================================
  subroutine show_domain_decomp(Dd)
    use ModProcMH, ONLY: iProc

    type(DomainDecompositionType),intent(in):: Dd
    integer:: iNode, iChild
    !-----------------------------------------------------------------
    if(iProc /= 0) RETURN

    write(*,*)'!!! Starting show_domain_decomp'
    write(*,*)'!!! CompID_            =', Dd%CompID_
    write(*,*)'!!! nDim               =', Dd%nDim
    write(*,*)'!!! XyzMin_D           =', Dd%XyzMin_D
    write(*,*)'!!! XyzMax_D           =', Dd%XyzMax_D
    write(*,*)'!!! iRootMapDim_D      =', Dd%iRootMapDim_D
    write(*,*)'!!! IsTreeDecomposition=', Dd%IsTreeDecomposition
    write(*,*)'!!! nDimTree           =', Dd%nDimTree
    write(*,*)'!!! nChildren          =', Dd%nChildren
    write(*,*)'!!! nDim               =', Dd%nDim
    write(*,*)'!!! IsTreeDecomposition=', Dd%IsTreeDecomposition
    write(*,*)'!!! nTreeNodes         =', Dd%nTreeNodes
    write(*,*)'!!! nAllocatedNodes    =', Dd%nAllocatedNodes
    write(*,*)'!!! IsPeriodic_D       =', Dd%IsPeriodic_D
    write(*,*)'!!! DoGlueMargins      =', Dd%DoGlueMargins
    write(*,*)'!!! iRealization       =', Dd%iRealization
    write(*,*)'!!! IsLocal            =', Dd%IsLocal
    write(*,*)'!!! MinBlock           =', Dd%MinBlock
    write(*,*)'!!! MaxBlock           =', Dd%MaxBlock
    write(*,*)'!!! nBlockAll          =', Dd%nBlockAll

    write(*,*)'!!! iChild, iShift_DI'
    do iChild = 1, Dd%nChildren
       write(*,*) iChild, Dd%iShift_DI(:,iChild)
    end do

    write(*,*)'!!! iNode, iDecomposition_II'
    do iNode = 1, Dd%nTreeNodes
       write(*,*) iNode, Dd%iDecomposition_II(:,iNode)
    end do

    write(*,*)'!!! iNode, XyzBlock_DI'
    do iNode = 1, Dd%nTreeNodes
       write(*,*) iNode, Dd%XyzBlock_DI(:,iNode), Dd%DXyzCell_DI(:,iNode)
    enddo

    write(*,*)'!!! Done with show_domain_decomp'

  end subroutine show_domain_decomp
  !===================================================================
  subroutine get_batl_tree(DomainDecomposition)

    ! Avoid name conflict with Parent_ in the SWMF coupling toolkit
    use BATL_tree, ParentBatl_ => Parent_

    type(DomainDecompositionType),intent(inout)::DomainDecomposition

    integer:: iNode, iNodeParent, iChild
    !------------------------------------------------------------------

    ! Allocate arrays for nNode sized tree
    DomainDecomposition%nTreeNodes = nNode
    call check_octree_grid_allocation(DomainDecomposition)

    ! Here we assume that there are no holes in the BATL tree
    do iNode = 1, nNode

       iNodeParent = iTree_IA(ParentBatl_,iNode)
       if(iNodeParent == Unset_)then
          ! For root blocks coupling toolkit seems to set parent to itself
          DomainDecomposition%iDecomposition_II(Parent_,iNode) = iNode
          ! For root blocks coupling toolkit seems to set child index to 0
          DomainDecomposition%iDecomposition_II(MyNumberAsAChild_,iNode) = 0
       else
          DomainDecomposition%iDecomposition_II(Parent_,iNode) = iNodeParent
          ! Find child index
          do iChild = 1, nChild
             if(iTree_IA(Child0_+iChild,iNodeParent) == iNode) then
                DomainDecomposition%iDecomposition_II(MyNumberAsAChild_,iNode)&
                     =iChild
                EXIT
             end if
          end do
       end if

       if(iTree_IA(Status_,iNode) == Unused_)then 
          do iChild = 1, nChild
             ! iChildOrder_II may be required here !!!
             DomainDecomposition%iDecomposition_II(iChild,iNode) = &
                  iTree_IA(Child0_+iChild,iNode) 
          end do
       else
          DomainDecomposition%iDecomposition_II(FirstChild_,iNode) = &
               None_
          DomainDecomposition%iDecomposition_II(GlobalBlock_,iNode) = &
               iMortonNode_A(iNode)
          DomainDecomposition%iDecomposition_II(PE_,iNode) = &
               iTree_IA(Proc_,iNode)
          DomainDecomposition%iDecomposition_II(PELast_,iNode) = &
               iTree_IA(Proc_,iNode)
          DomainDecomposition%iDecomposition_II(BLK_,iNode) = &
               iTree_IA(Block_,iNode)
          DomainDecomposition%iDecomposition_II(LEV_,iNode) = &
               iTree_IA(Level_,iNode)
          DomainDecomposition%iDecomposition_II(LEVmin_,iNode) = &
               iTree_IA(MinLevel_,iNode)
          DomainDecomposition%iDecomposition_II(LEVmax_,iNode) = &
               iTree_IA(MaxLevel_,iNode)
       end if
    end do

    ! call show_domain_decomp(DomainDecomposition)

  end subroutine get_batl_tree
  !===========================================================================
  subroutine pack_octree(DomainDecomposition)

    use ModOctree,ONLY: octree_roots, adaptive_block_ptr
    use ModParallel,ONLY: proc_dims 

    type(DomainDecompositionType),intent(inout)::DomainDecomposition

    logical::DoCountOnly
    integer :: i, j, k, l,iRoot
    type (adaptive_block_ptr) :: Octree
    !-----------------------------------------------------------------------
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

    ! call show_domain_decomp(DomainDecomposition)

  end subroutine pack_octree
  !============================================================================
  recursive subroutine pack_soln_block( &
       octree,l,lParent,iChildNumber,DoCountOnly,DomainDecomposition)

    use ModOctree, ONLY:adaptive_block_ptr, iChildOrder_II

    type (adaptive_block_ptr),intent(inout) :: octree
    integer,intent(inout) :: l
    integer,intent(in) :: lParent,iChildNumber
    logical,intent(inout) :: DoCountOnly
    type(DomainDecompositionType), intent(inout):: DomainDecomposition

    logical:: IsUsedBlock
    integer:: lOctree, iCube

    type (adaptive_block_ptr) :: child
    !---------------------------------------------------------------

    if (associated(octree % ptr)) then
       l=l+1
       lOctree=l
       DomainDecomposition%nTreeNodes=l
       DoCountOnly=DoCountOnly.or.(lOctree>DomainDecomposition%nAllocatedNodes)
       if(.not.DoCountOnly)then
          DomainDecomposition%iDecomposition_II(Parent_,lOctree)=lParent
          DomainDecomposition%iDecomposition_II(MyNumberAsAChild_,lOctree) &
               =iChildNumber
       end if
       IsUsedBlock = octree % ptr % used
       if (.not.IsUsedBlock) then
          do iCube = 1, 8
             if(.not.DoCountOnly)DomainDecomposition%iDecomposition_II(&
                  iChildOrder_II(iCube,iChildNumber), lOctree)=l+1
             child % ptr => &
                  octree % ptr % child(iChildOrder_II(iCube,iChildNumber))%ptr
             call pack_soln_block(child,l,lOctree,&
                  iChildOrder_II(iCube,iChildNumber),DoCountOnly,&
                  DomainDecomposition)
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
    use ModCube, ONLY: iShiftChild_DI

    type(DomainDecompositionType),intent(inout)::DomainDecomposition  
    !--------------------------------------------------------------------------

    if(UseBatl)then
       call get_root_decomposition_dd(&
            DomainDecomposition,&!Decomposition to be constructed
            proc_dims  ,&!As in DomainDecompositionType
            XyzMin_D,&   !As in DomainDecompositionType
            XyzMax_D,&   !As in DomainDecompositionType
            nIJK_D,&     !As in DomainDecompositionType
            IsPeriodic_D=periodic3d,&
            iShift_DI=iShiftMorton_DI)
    else
       call get_root_decomposition_dd(&
            DomainDecomposition,&!Decomposition to be constructed
            proc_dims  ,&!As in DomainDecompositionType
            XyzMin_D,&   !As in DomainDecompositionType
            XyzMax_D,&   !As in DomainDecompositionType
            nIJK_D,&     !As in DomainDecompositionType
            IsPeriodic_D=periodic3d,& 
            iShift_DI=iShiftChild_DI)       !As in DomainDecompositionType
    end if
  end subroutine MH_get_roots_dd
  !=====================================================================!
  subroutine MH_get_roots_id(GridID_)                         

    use ModSize,ONLY:nDim,nIJK_D
    use ModParallel,ONLY:periodic3d,proc_dims
    use ModGeometry,ONLY:XyzMin_D,XyzMax_D
    use ModCube, ONLY: iShiftChild_DI

    integer, intent(in):: GridID_  
    !-------------------------------------------------------------------------
    if(UseBatl)then
       call get_root_decomposition_id(&
            GridID_,&!Decomposition to be constructed
            proc_dims  ,&!As in DomainDecompositionType
            XyzMin_D,&   !As in DomainDecompositionType
            XyzMax_D,&   !As in DomainDecompositionType
            nIJK_D,&     !As in DomainDecompositionType
            IsPeriodic_D=periodic3d,&
            iShift_DI=iShiftMorton_DI)
    else
       call get_root_decomposition_id(&
            GridID_,&!Decomposition to be constructed
            proc_dims  ,&!As in DomainDecompositionType
            XyzMin_D,&   !As in DomainDecompositionType
            XyzMax_D,&   !As in DomainDecompositionType
            nIJK_D,&     !As in DomainDecompositionType
            IsPeriodic_D=periodic3d,& 
            iShift_DI=iShiftChild_DI)       !As in DomainDecompositionType
    end if
  end subroutine MH_get_roots_id

  !==========================================================================
  subroutine MH_update_local_decomposition(DomainDecomposition)
    use ModMain, ONLY: iNewGrid, iNewDecomposition

    type(DomainDecompositionType), intent(inout):: DomainDecomposition
    !-----------------------------------------------------------------------

    if(iNewGrid==iLastGrid.and.&
         iNewDecomposition==iLastDecomposition&
    .and.DomainDecomposition%iRealization/=0)return

    if(UseBatl)then
       call get_batl_tree(DomainDecomposition)
    else
       call pack_octree(DomainDecomposition)                  
    end if
    DomainDecomposition%iRealization=&
         mod(DomainDecomposition%iRealization+1,1000)   
    iLastDecomposition=iNewDecomposition
    iLastGrid=iNewGrid
    call complete_grid(DomainDecomposition)

  end subroutine MH_update_local_decomposition

end module MH_domain_decomposition
