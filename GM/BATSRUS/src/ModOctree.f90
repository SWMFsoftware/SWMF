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
     integer :: number, child_number, PE, BLK, LEV, LEVmin,LEVmax
     type (adaptive_block_ptr) :: parent
     type (adaptive_block_ptr) :: child(1:2**nDim)
  end type adaptive_block

 

  type (adaptive_block_ptr), dimension(:,:,:), allocatable :: octree_roots
  type (adaptive_block_ptr), dimension(:,:), allocatable :: global_block_ptrs
  type (adaptive_block_ptr), dimension(:), allocatable :: blocknumber_ptrs

end module ModOctree
