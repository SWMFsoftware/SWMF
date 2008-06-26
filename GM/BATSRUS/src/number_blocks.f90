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

  integer :: i,j,k,iBlock
  type (adaptive_block_ptr) :: octree

  do i = 1, nProc*nBLK
     nullify (blocknumber_ptrs(i) % ptr)
  end do

  iBlock = 0
  do k=1,proc_dims(3); do j=1,proc_dims(2); do i=1,proc_dims(1)
     octree % ptr => octree_roots(i, j, k) % ptr
     call renumber_octree_blocks(octree, iBlock)
  end do; end do; end do
  if(iBlock/=nBlockAll)then
     write(*,*)iBlock,nBlockAll
     call stop_mpi('Wrong count for nBlockAll')
  end if
  if(iProc==0.and.lVerbose>0)then
     call write_prefix; write(iUnitOut,*) &
          'renumber_octree: finished renumbering',nBlockALL,' blocks.'
  end if

end subroutine number_soln_blocks
!========================================================================!
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
        ChildOrder_a=iChildOrder_II(:,octree % ptr % child_number)
        ! Do children one by one in the selected order
        do iChild=1,8
           child % ptr => octree % ptr % child(ChildOrder_a(iChild))%ptr
           call renumber_octree_blocks(child, iBlockALL)
        end do
     end if
  end if
end subroutine renumber_octree_blocks
