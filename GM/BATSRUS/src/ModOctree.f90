!^CFG COPYRIGHT UM
module ModOctree
  use ModSize,ONLY:nDim
  implicit none
  save

  type adaptive_block_ptr
     type (adaptive_block), pointer :: ptr
  end type adaptive_block_ptr


  type adaptive_block
     logical :: &
          IsExtraBoundaryOrPole,IsOuterBoundary,&         
          used, refine, coarsen, body
     integer :: number, child_number, PE, BLK, LEV, LEVmin,LEVmax, &
          iRoot,jRoot, kRoot
     type (adaptive_block_ptr) :: parent
     type (adaptive_block_ptr) :: child(1:2**nDim)
  end type adaptive_block

 

  type (adaptive_block_ptr), dimension(:,:,:), allocatable :: octree_roots
  type (adaptive_block_ptr), dimension(:,:), allocatable :: global_block_ptrs
  type (adaptive_block_ptr), dimension(:), allocatable :: blocknumber_ptrs
  
  !\ 
  ! Ordering blocks along the Peano curve
  !/
  ! Decides order based on child number
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


contains
!---------------------------------------------------------------------
  logical function all_children_used(OctreeBlock)
    type(adaptive_block_ptr),intent(in)::OctreeBlock
    integer::iChild
    all_children_used=.true.
    do iChild=1,8
       all_children_used=all_children_used.and.&
            OctreeBlock%ptr%child(iChild)%ptr%used
    end do
  end function all_children_used
  !---------------------------------------------------------------------
  logical function all_children_coarsen(OctreeBlock)
    type(adaptive_block_ptr),intent(in)::OctreeBlock
    integer::iChild
    all_children_coarsen=.true.
    do iChild=1,8
       all_children_coarsen=all_children_coarsen.and.&
            OctreeBlock%ptr%child(iChild)%ptr%coarsen
    end do
  end function all_children_coarsen

end module ModOctree
