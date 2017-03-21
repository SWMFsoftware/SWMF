!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!MHD grid in BATSRUS 
module MH_domain_decomposition

  use CON_grid_storage, ProcToolkit_ => Pe_

  implicit none

  SAVE

  interface MH_get_root_decomposition
     module procedure MH_get_roots_id
     module procedure MH_get_roots_dd
  end interface

  public:: MH_get_root_decomposition
  public:: MH_update_local_decomposition

  type(DomainDecompositionType), public:: MH_DomainDecomposition
  type(DomainDecompositionType), public:: MH_LineDecomposition


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

  private:: get_batl_tree, &
       MH_get_roots_dd, MH_get_roots_id

contains

  !===========================================================================
  subroutine show_domain_decomp(Dd)
    use ModProcMH, ONLY: iProc

    type(DomainDecompositionType),intent(in):: Dd
    integer:: iNode, iChild
    !-------------------------------------------------------------------------
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

    write(*,*)'!!! iNode, XyzBlock_DI,  DXyz_DI'
    do iNode = 1, Dd%nTreeNodes
       write(*,*) iNode, Dd%XyzBlock_DI(:,iNode), Dd%DXyzCell_DI(:,iNode)
    enddo

    write(*,*)'!!! Done with show_domain_decomp'

  end subroutine show_domain_decomp
  !===========================================================================
  subroutine get_batl_tree(DomainDecomposition)

    ! Avoid name conflict with Parent_ in the SWMF coupling toolkit
    use BATL_tree, ParentBatl_ => Parent_

    type(DomainDecompositionType),intent(inout)::DomainDecomposition

    integer:: iNode, iNodeParent, iChild
    !-------------------------------------------------------------------------

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
          DomainDecomposition%iDecomposition_II(ProcToolkit_,iNode) = &
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
  subroutine MH_get_roots_dd(DomainDecomposition)                         

    use BATL_lib, ONLY: nIJK_D, IsPeriodic_D, nRoot_D, CoordMin_D, CoordMax_D

    type(DomainDecompositionType),intent(inout)::DomainDecomposition  
    !-------------------------------------------------------------------------

    call get_root_decomposition_dd(&
         DomainDecomposition,       & ! Decomposition to be constructed
         nRoot_D,                   & ! As in DomainDecompositionType
         CoordMin_D,                & ! As in DomainDecompositionType
         CoordMax_D,                & ! As in DomainDecompositionType
         nIJK_D,                    & ! As in DomainDecompositionType
         IsPeriodic_D=IsPeriodic_D, &
         iShift_DI=iShiftMorton_DI)

  end subroutine MH_get_roots_dd
  !===========================================================================
  subroutine MH_get_roots_id(GridID_)                         

    use BATL_lib, ONLY: nIJK_D, IsPeriodic_D, nRoot_D, CoordMin_D, CoordMax_D

    integer, intent(in):: GridID_  
    !-------------------------------------------------------------------------
    call get_root_decomposition_id(&
         GridID_,                   & ! Decomposition to be constructed
         nRoot_D,                   & ! As in DomainDecompositionType
         CoordMin_D,                & ! As in DomainDecompositionType
         CoordMax_D,                & ! As in DomainDecompositionType
         nIJK_D,                    & ! As in DomainDecompositionType
         IsPeriodic_D=IsPeriodic_D, &
         iShift_DI=iShiftMorton_DI)

  end subroutine MH_get_roots_id

  !==========================================================================
  subroutine MH_update_local_decomposition(DomainDecomposition)
    use ModMain, ONLY: iNewGrid, iNewDecomposition

    type(DomainDecompositionType), intent(inout):: DomainDecomposition
    !-----------------------------------------------------------------------

    if(iNewGrid==iLastGrid .and. iNewDecomposition == iLastDecomposition &
         .and. DomainDecomposition%iRealization /= 0) &
         RETURN

    call get_batl_tree(DomainDecomposition)

    DomainDecomposition%iRealization = &
         mod(DomainDecomposition%iRealization+1, 1000)
    iLastDecomposition = iNewDecomposition
    iLastGrid          = iNewGrid
    call complete_grid(DomainDecomposition)

  end subroutine MH_update_local_decomposition

end module MH_domain_decomposition
