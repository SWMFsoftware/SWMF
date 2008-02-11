!^CFG COPYRIGHT UM
module ModOctree
  implicit none
  save

  type adaptive_block
     logical :: &
          IsExtraBoundaryOrPole,IsOuterBoundary,&         
          used, refine, coarsen, body
     integer :: number, child_number, PE, BLK, LEV, LEVmin,LEVmax
     type (adaptive_block), pointer :: parent
     type (adaptive_block), pointer :: child1, child2, child3, child4, &
          child5, child6, child7, child8
  end type adaptive_block

  type adaptive_block_ptr
     type (adaptive_block), pointer :: ptr
  end type adaptive_block_ptr

  type (adaptive_block_ptr), dimension(:,:,:), allocatable :: octree_roots
  type (adaptive_block_ptr), dimension(:,:), allocatable :: global_block_ptrs
  type (adaptive_block_ptr), dimension(:), allocatable :: blocknumber_ptrs

end module ModOctree
