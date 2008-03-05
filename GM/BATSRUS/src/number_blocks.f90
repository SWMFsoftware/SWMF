!^CFG COPYRIGHT UM
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
subroutine number_soln_blocks
  use ModProcMH
  use ModMain, ONLY : nBLK,nBlockALL,lVerbose
  use ModParallel, ONLY : proc_dims
  use ModOctree
  use ModIO, ONLY : iUnitOut, write_prefix
  implicit none

  integer :: i,j,k
  type (adaptive_block_ptr) :: octree

  do i = 1, nProc*nBLK
     nullify (blocknumber_ptrs(i) % ptr)
  end do

  nBlockALL = 0
  do k=1,proc_dims(3); do j=1,proc_dims(2); do i=1,proc_dims(1)
     octree % ptr => octree_roots(i, j, k) % ptr
     call renumber_octree_blocks(octree, nBlockALL)
  end do; end do; end do
  if(iProc==0.and.lVerbose>0)then
     call write_prefix; write(iUnitOut,*) &
          'renumber_octree: finished renumbering',nBlockALL,' blocks.'
  end if

end subroutine number_soln_blocks

recursive subroutine renumber_octree_blocks(octree, iBlockALL)
  use ModProcMH
  use ModMain, ONLY : global_block_number
  use ModParallel, ONLY : iBlock_A,iProc_A,iBlockRestartALL_A
  use ModOctree
  implicit none

  integer, intent(inout) :: iBlockALL
  type (adaptive_block_ptr) :: octree, child

  integer :: iChild, ChildOrder_a(8)
  !---------------------------------------------------------------------------

  if (associated(octree % ptr)) then
     if (octree % ptr % used) then
        iBlockALL = iBlockALL +1
        ! Store number from octree for restart filename
        iBlockRestartALL_A(iBlockALL)=octree % ptr % number
        ! Overwrite number using the Peano-Hilbert ordering
        octree % ptr % number = iBlockALL
        iBlock_A(iBlockALL) = octree % ptr % BLK
        iProc_A(iBlockALL) = octree % ptr % PE
        blocknumber_ptrs(iBlockALL) % ptr => octree % ptr
        if(iProc==octree % ptr % PE)&
             global_block_number(octree % ptr % BLK)=iBlockALL
     else
        ! Decide order based on child number
        select case(octree % ptr % child_number)
        case(4)
           ChildOrder_a=(/4,3,6,5,8,7,2,1/)
        case(1)
           ChildOrder_a=(/4,1,2,3,6,7,8,5/)
        case(0,7,8)
           ChildOrder_a=(/4,1,8,5,6,7,2,3/)
        case(5)
           ChildOrder_a=(/2,1,4,3,6,5,8,7/)
        case(6)
           ChildOrder_a=(/8,7,6,5,4,3,2,1/)
        case(2)
           ChildOrder_a=(/6,7,8,5,4,1,2,3/)
        case(3)
           ChildOrder_a=(/2,1,8,7,6,5,4,3/)
        case default
           write(*,*)'iBlockALL, child_number=',&
                iBlockALL,octree % ptr % child_number
           call stop_mpi('impossible child_number')
        end select

        ! Do children one by one in the selected order
        do iChild=1,8
           child % ptr => octree % ptr % child(ChildOrder_a(iChild))%ptr
           call renumber_octree_blocks(child, iBlockALL)
        end do
     end if
  end if
end subroutine renumber_octree_blocks
