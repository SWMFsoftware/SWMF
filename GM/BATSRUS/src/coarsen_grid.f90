!^CFG COPYRIGHT UM
!------------------------------------------------------------------------
! Front end routine to collect coarsened blocks from each processor into
!   global list before calling parallel_coarsen
!------------------------------------------------------------------------
subroutine coarsen_grid(coarsenBLK)
  use ModProcMH
  use ModMain, ONLY : nBLK
  use ModAMR, ONLY : coarsen_list
  use ModMpi
  implicit none

  logical, intent(inout) :: coarsenBLK(nBLK)
  integer :: iError

  call MPI_ALLGATHER(coarsenBLK, nBLK, MPI_LOGICAL, &
       coarsen_list, nBLK, MPI_LOGICAL, iComm, iError)

  call parallel_coarsen

end subroutine coarsen_grid

!------------------------------------------------------------------------
! parallel_coarsen assumes that "coarsen_list" is filled
!   Execute coarsening in parallel
!------------------------------------------------------------------------
subroutine parallel_coarsen
  use ModProcMH
  use ModMain, ONLY : nBLK,iNewGrid
  use ModAMR
  use ModOctree
  use ModMpi
  implicit none

  integer :: iPE, iBLK, jBLK, icube
  integer :: nPEsCrseBlk, PEsCrseBlk(8), newPE, ii
  real    :: xf1, yf1, zf1, zoff
  
  logical :: oktest, oktest_me
  
  type (adaptive_block_ptr) :: coarse_block_ptr, parent_block_ptr

  !---------------------------------------------------------------------------

  call set_oktest('coarsen',oktest,oktest_me)
  
  if (.not. any(coarsen_list)) RETURN

  ! Change gird identification number
  iNewGrid = mod( iNewGrid+1, 10000)

  nullify(coarse_block_ptr % ptr)
  nullify(parent_block_ptr % ptr)

  do iPE=1,nProc; do jBLK=1,nBLK
     ! check each solution block for coarsening
     if(.not.coarsen_list(jBLK, iPE)) CYCLE
        
     coarse_block_ptr % ptr => global_block_ptrs(jBLK, iPE) % ptr
     
     ! ensure flagged block is not top level block
     if (associated(coarse_block_ptr % ptr)) then
        parent_block_ptr % ptr => coarse_block_ptr % ptr % parent%ptr
     else
        nullify(parent_block_ptr % ptr)
     end if
     if(.not.associated(parent_block_ptr % ptr)) then 
        ! flagged block is top level block, turn off coarsen flag
        coarsen_list(jBLK, iPE) = .false.
        CYCLE
     end if
     
     ! ensure all blocks to be coarsened are in use
     if (parent_block_ptr % ptr % child(1)%ptr % used .and. &
         parent_block_ptr % ptr % child(2)%ptr % used .and. &
         parent_block_ptr % ptr % child(3)%ptr % used .and. &
         parent_block_ptr % ptr % child(4)%ptr % used .and. &
         parent_block_ptr % ptr % child(5)%ptr % used .and. &
         parent_block_ptr % ptr % child(6)%ptr % used .and. &
         parent_block_ptr % ptr % child(7)%ptr % used .and. &
         parent_block_ptr % ptr % child(8)%ptr % used ) then 
        
        local_cube(1) = parent_block_ptr % ptr % child(1)%ptr % PE
        local_cube(2) = parent_block_ptr % ptr % child(2)%ptr % PE
        local_cube(3) = parent_block_ptr % ptr % child(3)%ptr % PE
        local_cube(4) = parent_block_ptr % ptr % child(4)%ptr % PE
        local_cube(5) = parent_block_ptr % ptr % child(5)%ptr % PE
        local_cube(6) = parent_block_ptr % ptr % child(6)%ptr % PE
        local_cube(7) = parent_block_ptr % ptr % child(7)%ptr % PE
        local_cube(8) = parent_block_ptr % ptr % child(8)%ptr % PE
        local_cubeBLK(1) = parent_block_ptr % ptr % child(1)%ptr % BLK
        local_cubeBLK(2) = parent_block_ptr % ptr % child(2)%ptr % BLK
        local_cubeBLK(3) = parent_block_ptr % ptr % child(3)%ptr % BLK
        local_cubeBLK(4) = parent_block_ptr % ptr % child(4)%ptr % BLK
        local_cubeBLK(5) = parent_block_ptr % ptr % child(5)%ptr % BLK
        local_cubeBLK(6) = parent_block_ptr % ptr % child(6)%ptr % BLK
        local_cubeBLK(7) = parent_block_ptr % ptr % child(7)%ptr % BLK
        local_cubeBLK(8) = parent_block_ptr % ptr % child(8)%ptr % BLK
        
        nPEsCrseBlk = 1
        PEsCrseBlk = -1
        PEsCrseBlk(1) = local_cube(1)
        do icube = 2, 8
           newPE = 1
           do ii = 1,icube-1
              if (local_cube(icube) == local_cube(ii)) newPE = 0
           end do
           nPEsCrseBlk = nPEsCrseBlk + newPE
           if (newPE == 1) PEsCrseBlk(nPEsCrseBlk) = local_cube(icube)
        end do
        
        if(oktest) then
           if (iProc == local_cube(1)) then
              write(*,'(a,i4,/,a,i1,/,3(a,8i4,/),a,8l4,/)') &
                   "C--> coarsen_grid: PE = ",iProc, &
                   "C-     nPEsCrseBlk   = ",nPEsCrseBlk, &
                   "C-     PEsCrseBlk    = ",PEsCrseBlk, &
                   "C-     local_cube    = ",local_cube, &
                   "C-     local_cubeBLK = ",local_cubeBLK, &
                   "C-     coarsen_list  = ", &
                   coarsen_list(local_cubeBLK(1),local_cube(1)+1), &
                   coarsen_list(local_cubeBLK(2),local_cube(2)+1), &
                   coarsen_list(local_cubeBLK(3),local_cube(3)+1), &
                   coarsen_list(local_cubeBLK(4),local_cube(4)+1), &
                   coarsen_list(local_cubeBLK(5),local_cube(5)+1), &
                   coarsen_list(local_cubeBLK(6),local_cube(6)+1), &
                   coarsen_list(local_cubeBLK(7),local_cube(7)+1), &
                   coarsen_list(local_cubeBLK(8),local_cube(8)+1)
           end if
        end if
        
        ! coarsen group of 8 blocks when permitted
        if (coarsen_list(local_cubeBLK(1), local_cube(1)+1) .and. &
            coarsen_list(local_cubeBLK(2), local_cube(2)+1) .and. &
            coarsen_list(local_cubeBLK(3), local_cube(3)+1) .and. &
            coarsen_list(local_cubeBLK(4), local_cube(4)+1) .and. &
            coarsen_list(local_cubeBLK(5), local_cube(5)+1) .and. &
            coarsen_list(local_cubeBLK(6), local_cube(6)+1) .and. &
            coarsen_list(local_cubeBLK(7), local_cube(7)+1) .and. &
            coarsen_list(local_cubeBLK(8), local_cube(8)+1) .and. &
            .not. parent_block_ptr % ptr % used ) then 
           
           coarsen_list(local_cubeBLK(1), local_cube(1)+1) = .false.
           coarsen_list(local_cubeBLK(2), local_cube(2)+1) = .false.
           coarsen_list(local_cubeBLK(3), local_cube(3)+1) = .false.
           coarsen_list(local_cubeBLK(4), local_cube(4)+1) = .false.
           coarsen_list(local_cubeBLK(5), local_cube(5)+1) = .false.
           coarsen_list(local_cubeBLK(6), local_cube(6)+1) = .false.
           coarsen_list(local_cubeBLK(7), local_cube(7)+1) = .false.
           coarsen_list(local_cubeBLK(8), local_cube(8)+1) = .false.
           
           call coarsen_octree_block(parent_block_ptr, &
                                     local_cube, local_cubeBLK)
           
           do icube = 2, 8
              availableBLKs(0,local_cube(icube)) = availableBLKs(0,local_cube(icube))-1
              availableBLKs(availableBLKs(0,local_cube(icube)),local_cube(icube)) = &
                   local_cubeBLK(icube)
           end do
           
           call create_coarse_soln_block(nPEsCrseBlk, PEsCrseBlk)
           
        else
           ! coarsen not permitted, turn off coarsen flag
           coarsen_list(local_cubeBLK(1), local_cube(1)+1) = .false.
           coarsen_list(local_cubeBLK(2), local_cube(2)+1) = .false.
           coarsen_list(local_cubeBLK(3), local_cube(3)+1) = .false.
           coarsen_list(local_cubeBLK(4), local_cube(4)+1) = .false.
           coarsen_list(local_cubeBLK(5), local_cube(5)+1) = .false.
           coarsen_list(local_cubeBLK(6), local_cube(6)+1) = .false.
           coarsen_list(local_cubeBLK(7), local_cube(7)+1) = .false.
           coarsen_list(local_cubeBLK(8), local_cube(8)+1) = .false.
        end if
        
     else
        ! some blocks to be coarsened are not in use, turn off coarsen flag
        if               (parent_block_ptr%ptr%child(1)%ptr%used) &
             coarsen_list(parent_block_ptr%ptr%child(1)%ptr%BLK, &
                          parent_block_ptr%ptr%child(1)%ptr%PE+1) = .false.
        
        if               (parent_block_ptr%ptr%child(2)%ptr%used) &
             coarsen_list(parent_block_ptr%ptr%child(2)%ptr%BLK, &
                          parent_block_ptr%ptr%child(2)%ptr%PE+1) = .false.
        
        if               (parent_block_ptr%ptr%child(3)%ptr%used) &
             coarsen_list(parent_block_ptr%ptr%child(3)%ptr%BLK, &
                          parent_block_ptr%ptr%child(3)%ptr%PE+1) = .false.
        
        if               (parent_block_ptr%ptr%child(4)%ptr%used) &
             coarsen_list(parent_block_ptr%ptr%child(4)%ptr%BLK, &
                          parent_block_ptr%ptr%child(4)%ptr%PE+1) = .false.
        
        if               (parent_block_ptr%ptr%child(5)%ptr%used) &
             coarsen_list(parent_block_ptr%ptr%child(5)%ptr%BLK, &
                          parent_block_ptr%ptr%child(5)%ptr%PE+1) = .false.
        
        if               (parent_block_ptr%ptr%child(6)%ptr%used) &
             coarsen_list(parent_block_ptr%ptr%child(6)%ptr%BLK, &
                          parent_block_ptr%ptr%child(6)%ptr%PE+1) = .false.
        
        if               (parent_block_ptr%ptr%child(7)%ptr%used) &
             coarsen_list(parent_block_ptr%ptr%child(7)%ptr%BLK, &
                          parent_block_ptr%ptr%child(7)%ptr%PE+1) = .false.
        
        if               (parent_block_ptr%ptr%child(8)%ptr%used) &
             coarsen_list(parent_block_ptr%ptr%child(8)%ptr%BLK, &
                          parent_block_ptr%ptr%child(8)%ptr%PE+1) = .false.
     end if
  end do; end do
end subroutine parallel_coarsen
