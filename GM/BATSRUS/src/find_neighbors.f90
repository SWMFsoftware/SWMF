 !^CFG COPYRIGHT UM
!//////////////////////////////////////////////////////////////////////////////
!
!    North                      8--------------7
!    y,j    Top               / |            / |
!     ^     z,k             /   |          /   |
!     |   /               /     |        /     |
!     | /     West      5-------+------6       |
!     |-----> x,i       |       |      |       |
!                       |       1------+-------2
!                       |     /        |     /
!                       |   /          |   /
!                       | /            | /
!                       4--------------3 --(front, right, and bottom 3D box corner)
!
!  Point 7 is:  +x, +y, +z, West, North, Top
!
!  Point 4 is:  -x, -y, -z, East, South, Bottom
!

subroutine find_neighbors
  use ModProcMH, ONLY: iProc
  use ModMain,   ONLY: nBlock, UnusedBLK, UseBatl
  use ModParallel
  use ModPolarNeighbor
  implicit none

  integer :: i, j, k, iBlock, iLevelOut
  integer, dimension(4) :: iProcOut_I, iBlockOut_I, iChildOut_I
  logical::DoTest=.false.,DoTestMe=.false.,DoCallOKTest=.true.
  character(LEN=*),parameter::NameSub='find_neighbors'
  !--------------------------------------------------------!

  if(UseBatl) RETURN


  if(DoCallOKTest)then
     call set_oktest(NameSub,DoTest,DoTestMe)
     DoCallOKTest=.false.
  end if
  !---------------------------------------------------------------------------
  !\
  ! Loop and fill variables with neighbor information for all blocks on this PE
  !/
  do iBlock = 1, nBlock
     if (UnusedBLK(iBlock)) CYCLE

     do k=-1,1; do j=-1,1; do i=-1,1
        call treeNeighbor(iProc, iBlock, i, j, k, &
             iProcOut_I, iBlockOut_I, iChildOut_I, iLevelOut)
        BLKneighborPE   (i,j,k,:,iBlock) = iProcOut_I
        BLKneighborBLK  (i,j,k,:,iBlock) = iBlockOut_I
        BLKneighborCHILD(i,j,k,:,iBlock) = iChildOut_I
        BLKneighborLEV  (i,j,k,  iBlock) = iLevelOut
     end do; end do; end do
  end do

  ! Fill up arrays which contain information about face neighbors only

  neiLeast(1:nBlock)    = BLKneighborLEV(-1,0,0,1:nBlock)
  neiLwest(1:nBlock)    = BLKneighborLEV(+1,0,0,1:nBlock)
  neiLsouth(1:nBlock)   = BLKneighborLEV(0,-1,0,1:nBlock)
  neiLnorth(1:nBlock)   = BLKneighborLEV(0,+1,0,1:nBlock)
  neiLbot(1:nBlock)     = BLKneighborLEV(0,0,-1,1:nBlock)
  neiLtop(1:nBlock)     = BLKneighborLEV(0,0,+1,1:nBlock)

  neiPeast(:,1:nBlock)  = BLKneighborPE(-1,0,0,:,1:nBlock)
  neiPwest(:,1:nBlock)  = BLKneighborPE(+1,0,0,:,1:nBlock)
  neiPsouth(:,1:nBlock) = BLKneighborPE(0,-1,0,:,1:nBlock)
  neiPnorth(:,1:nBlock) = BLKneighborPE(0,+1,0,:,1:nBlock)
  neiPbot(:,1:nBlock)   = BLKneighborPE(0,0,-1,:,1:nBlock)
  neiPtop(:,1:nBlock)   = BLKneighborPE(0,0,+1,:,1:nBlock)

  neiBeast(:,1:nBlock)  = BLKneighborBLK(-1,0,0,:,1:nBlock)
  neiBwest(:,1:nBlock)  = BLKneighborBLK(+1,0,0,:,1:nBlock)
  neiBsouth(:,1:nBlock) = BLKneighborBLK(0,-1,0,:,1:nBlock)
  neiBnorth(:,1:nBlock) = BLKneighborBLK(0,+1,0,:,1:nBlock)
  neiBbot(:,1:nBlock)   = BLKneighborBLK(0,0,-1,:,1:nBlock)
  neiBtop(:,1:nBlock)   = BLKneighborBLK(0,0,+1,:,1:nBlock)

  neiLEV(1,:) = neiLeast
  neiLEV(2,:) = neiLwest
  neiLEV(3,:) = neiLsouth
  neiLEV(4,:) = neiLnorth
  neiLEV(5,:) = neiLbot
  neiLEV(6,:) = neiLtop

  neiPE(:,1,:) = neiPeast
  neiPE(:,2,:) = neiPwest
  neiPE(:,3,:) = neiPsouth
  neiPE(:,4,:) = neiPnorth
  neiPE(:,5,:) = neiPbot
  neiPE(:,6,:) = neiPtop

  neiBLK(:,2,:) = neiBwest
  neiBLK(:,1,:) = neiBeast
  neiBLK(:,3,:) = neiBsouth
  neiBLK(:,4,:) = neiBnorth
  neiBLK(:,5,:) = neiBbot
  neiBLK(:,6,:) = neiBtop

end subroutine find_neighbors
!=============================================================================
subroutine treeNeighbor(inPE, inBLK, ix, iy, iz, &
     outPE, outBLK, outCHILD, outLEV)
  use ModParallel, ONLY : NOBLK
  use ModOctree
  use ModCube
  implicit none

  integer, intent(in) :: inPE,inBLK,ix,iy,iz
  integer, intent(out) :: outLEV
  integer, intent(out), dimension(4) :: outPE, outBLK, outCHILD
  integer :: i,iSubFace
  logical :: noNeighbor
  type (adaptive_block_ptr) :: inBlockPtr,outBlockPtr,tmpBlockPtr

  if ( ix .lt. -1 .or. ix .gt. 1 .or. &
       iy .lt. -1 .or. iy .gt. 1 .or. &
       iz .lt. -1 .or. iz .gt. 1 ) then
     write (*,*) 'ERROR: treeNeighbor: 0: invalid direction given  ',ix,iy,iz
     !     stop
  end if
  outLEV   = NOBLK
  outPE    = NOBLK
  outBLK   = NOBLK
  outCHILD = NOBLK

  inBlockPtr%ptr => global_block_ptrs(inBLK,inPE+1)%ptr
  if (associated(inBlockPtr%ptr)) then ! ensure block is allocated
     !While applied near the pole the present subroutine is often called
     !with ix=0,  iy=0,  iz=0. To improve an efficiency, consider
     !this case now.
     if(iX*iX+iY*iY+iZ*iZ==0)then 
        OutChild(1)=InBlockPtr%ptr%child_number
        OutPe(1)=inPE
        OutBlk(1)=inBLK
        OutLev=0
        return
     end if
     call find_tree_neighbor(inBlockPtr,outBlockPtr,iX,iY,iZ,noNeighbor)
     if (noNeighbor) then
        outLEV   = NOBLK
        outPE    = NOBLK
        outBLK   = NOBLK
        outCHILD = NOBLK
     else
        if (associated(outBlockPtr%ptr)) then ! ensure block is allocated
           if (.not. associated(outBlockPtr%ptr%child(1)%ptr)) then
              if (inBlockPtr%ptr%LEV == outBlockPtr%ptr%LEV) then
                 outLEV = 0
              else
                 outLEV = 1
              end if
              outPE(1)    = outBlockPtr%ptr%PE
              outBLK(1)   = outBlockPtr%ptr%BLK
              outCHILD(1) = outBlockPtr%ptr%child_number
           else
              outLEV = -1
              !Set values for refined neighbor
              outPE  = NOBLK
              outBLK = NOBLK
              call get_children_list(iX,iY,iZ,OutCHILD)
              do iSubFace=1,4
                 if(outCHILD(iSubface)==NOBLK)EXIT
                 outPE (iSubFace)   = outBlockPtr%ptr%child(&
                      outCHILD(iSubFace))%ptr%PE
                 outBLK(iSubFace)   = outBlockPtr%ptr%child(&
                      outCHILD(iSubFace))%ptr%BLK
              end do
           end if
           !Error check
           do i=1,4
              if (outPE(i) .ge. 0) then
                 tmpBlockPtr%ptr => global_block_ptrs(outBLK(i),outPE(i)+1)%ptr
                 if (associated(tmpBlockPtr%ptr)) then ! ensure block is allocated
                    if (.not. tmpBlockPtr%ptr%used) then
                       write (*,*) 'ERROR: treeNeighbor: 4: block not used ',inPE,inBLK
                       !                       stop
                    end if
                 else
                    write (*,*) 'ERROR: treeNeighbor: 3: not associated ',inPE,inBLK
                    !                    stop
                 end if
              end if
           end do
        else
           write (*,*) 'ERROR: treeNeighbor: 2: not associated ',inPE,inBLK
           !           stop
        end if
     end if
  else
     write (*,*) 'ERROR: treeNeighbor: 1: not associated ',inPE,inBLK
     !     stop
  end if
end subroutine treeNeighbor
!============================================================================!
subroutine find_tree_neighbor(TreeNodeIn,TreeNodeOut,iX,iY,iZ,IsNoNeighbor)
  use ModParallel, ONLY : proc_dims,periodic3D
  use ModOctree
  use ModCube
  implicit none

  integer, intent(in) :: iX,iY,iZ
  type (adaptive_block_ptr),intent(in)  ::TreeNodeIn 
  type (adaptive_block_ptr),intent(out) ::TreeNodeOut
  logical, intent(inout) :: IsNoNeighbor
  !--------------------------------------------------------------------------!
  integer::iDXyz,iBin,i,j,k
  integer,dimension(3)::iXyzFromCorner_D,iRoot_D,iShift_D
  
  IsNoNeighbor=.false.

  nullify(TreeNodeOut%ptr)

  TreeNodeOut%ptr=>TreeNodeIn%ptr

  if(iX*iX+iY*iY+iZ*iZ==0)return
  
  !Assume the In block size along each coordinate to be 2

  iDXyz=2

  !The center of the neighboring block CENTER with respect to
  !the original block CORNER are
  iXyzFromCorner_D = 1 + (/iX,iY,iZ/)*iDXyz
  

  !Descend and find a common predator
  do while(TreeNodeOut%ptr%LEV>0)
     iXyzFromCorner_D = iXyzFromCorner_D + &
          iShiftChild_DI(:,TreeNodeOut%ptr%child_number)*iDXyz
     iDXyz=iDXyz*2
     TreeNodeOut%ptr=>TreeNodeOut%ptr%parent%ptr
     if(maxval(iXyzFromCorner_D)<iDXyz.and.minval(iXyzFromCorner_D)>0)then
        !Found common predator
        call ascend_tree
        return
     end if
  end do
  
  !Found root cell, find neighbor root cell
  iRoot_D(1) = TreeNodeOut%ptr%iRoot 
  iRoot_D(2) = TreeNodeOut%ptr%jRoot
  iRoot_D(3) = TreeNodeOut%ptr%kRoot
  !Check
  if (.not.associated(TreeNodeOut%ptr,&
       octree_roots(iRoot_D(1),iRoot_D(2),iRoot_D(3))%ptr))then
     write(*,*)'iRoot_D=', iRoot_D
     if(associated(octree_roots(iRoot_D(1),iRoot_D(2),iRoot_D(3))%ptr))then
        write(*,*)octree_roots(iRoot_D(1),iRoot_D(2),iRoot_D(3))%ptr%iRoot
        write(*,*)octree_roots(iRoot_D(1),iRoot_D(2),iRoot_D(3))%ptr%jRoot
        write(*,*)octree_roots(iRoot_D(1),iRoot_D(2),iRoot_D(3))%ptr%kRoot
     else
        write(*,*)'Octree with these root numbers is not associated'
     end if
     call stop_mpi('Failure in the algorithm for finding octree root')
  end if
  !Pass to the neighboring tree, modify iXyzFromCorner_D accordingly
  where(iXyzFromCorner_D<0)
     iXyzFromCorner_D = iXyzFromCorner_D + iDXyz
     iRoot_D          = iRoot_D          - 1
  end where
  where(iXyzFromCorner_D>iDXyz)
     iXyzFromCorner_D = iXyzFromCorner_D - iDXyz
     iRoot_D          = iRoot_D          + 1
  end where
  where (periodic3D.and.iRoot_D==0)iRoot_D=Proc_Dims
  where (periodic3D.and.iRoot_D>Proc_Dims)iRoot_D=iRoot_D-Proc_Dims
  if(any(iRoot_D==0.or.iRoot_D>Proc_Dims))then
     IsNoNeighbor=.true.
     return
  end if
  TreeNodeOut%ptr => octree_roots(iRoot_D(1),iRoot_D(2),iRoot_D(3))%ptr
  call ascend_tree
contains
  subroutine ascend_tree
    !Ascend and find block which contains the ppoint Xyz
    do while(iDXyz>2)
       if(.not.associated(TreeNodeOut%ptr%child(1)%ptr))EXIT
       iDXyz=iDXyz/2
       iShift_D=iXyzFromCorner_D/iDXyz
       iXyzFromCorner_D=iXyzFromCorner_D-iShift_D*iDXyz
       iBin = 4 * iShift_D(1) + 2 * iShift_D(2) + iShift_D(3) +1
       TreeNodeOut%ptr=>TreeNodeOut%ptr%child(iBin2Child_I(iBin))%ptr
    end do
    if (.not. associated(TreeNodeOut%ptr)) then
       write (*,*) 'ERROR: find_tree_neighbor: 1: not associated '
       !           stop
    end if
  end subroutine ascend_tree
end subroutine find_tree_neighbor
!===========================================================================!
subroutine fix_refinement_and_refine
  use ModProcMH
  use ModMain, ONLY : nBLK,lVerbose,nBlockAll
  use ModAMR, ONLY : refine_list
  use ModGeometry,ONLY: is_axial_geometry                
  use ModOctree
  implicit none

  integer :: i,j,k, idir, inPE,inBLK, iLev1, iLEV2, iCount
  integer :: iCountOld,iLoop
  logical :: noNeighbor
  type (adaptive_block_ptr) :: inBlockPtr,outBlockPtr,tmpBlockPtr
  
  do iLoop=1,2
     iCount = count(refine_list)
     FIX: do
        iCountOld=iCount
        do inPE=1,nProc; do inBLK=1,nBLK
           inBlockPtr%ptr => global_block_ptrs(inBLK,inPE)%ptr
           if (.not. associated(inBlockPtr%ptr)) CYCLE
           iLEV1 = inBlockPtr%ptr%LEV
           if (refine_list(inBlockPtr%ptr%BLK,inBlockPtr%ptr%PE+1)) &
                iLev1=iLev1+1
           do i=-1,1; do j=-1,1; do k=-1,1
              call find_tree_neighbor(inBlockPtr,outBlockPtr,i,j,k,noNeighbor)
              if (noNeighbor) CYCLE
              if ( .not. associated(outBlockPtr%ptr)) CYCLE
              if ( associated(outBlockPtr%ptr%child(1)%ptr)) CYCLE

              iLEV2 = outBlockPtr%ptr%LEV
                    
              if (refine_list(outBlockPtr%ptr%BLK,outBlockPtr%ptr%PE+1))&
                   iLEV2 = iLEV2 + 1
              if (( iLEV1-iLEV2 > 1 .or. &
                   (iLEV1-iLEV2 > 0 .and. &
                   ((inBlockPtr%ptr%body .and. outBlockPtr%ptr%body)&
                   .or.(inBlockPtr%ptr%IsOuterBoundary .and. &    
                   outBlockPtr%ptr%IsOuterBoundary)&
                   .or.(inBlockPtr%ptr%IsExtraBoundaryOrPole .and.& 
                   outBlockPtr%ptr%IsExtraBoundaryOrPole &
                   .and.( (i==0.and.k==0).or.(.not.is_axial_geometry()))&
                   )))) &
                   .and. .not. &
                   refine_list(outBlockPtr%ptr%BLK, outBlockPtr%ptr%PE+1) )then
                 iCount = iCount +1
                 refine_list(outBlockPtr%ptr%BLK,outBlockPtr%ptr%PE+1) = .true.
              end if
           end do; end do; end do
        end do;end do
        if(iCount==iCountOld)EXIT FIX
     end do FIX

     if (iCount == 0) return

     if(nBlockAll + 7*iCount > nProc*nBLK)then
        if(iLoop == 1)then
           if(iProc == 0)then
              write(*,*)'WARNING: Found ', iCount, ' blocks to refine'
              write(*,*)'  Insufficient number of blocks, ', &
                   'skipping refinement!'
           end if
           refine_list=.false.
           CYCLE
        else
           call stop_mpi('Insufficient number of blocks, cannot fix grid')
        end if
     end if
     call parallel_refine
     call set_body_flag
     refine_list = .false.
  end do
     
end subroutine fix_refinement_and_refine
!=============================================================================!


subroutine get_children_list_orig(iX,iY,iZ,OutCHILD)
  use ModParallel,ONLY:NOBLK
  implicit none
  integer,intent(in)::iX,iY,iZ
  integer,intent(out)::OutCHILD(4)
  integer::iDirHere
  outCHILD=NOBLK
  idirHere = 1+ (iz+1) + 3*(iy+1) + 9*(ix+1)
  select case (idirhere)
  case (1)  ! -X -Y -Z    ESB
     outCHILD(1) = 7 
  case (2)  ! -X -Y       ES 
     outCHILD(1:2) =(/ 6,7/)
  case (3)  ! -X -Y +Z    EST
     outCHILD(1) = 6
  case (4)  ! -X    -Z    E B
     outCHILD(1:2) =(/ 2,7/)
  case (5)  ! -X          E  
     outCHILD(1:4) =(/ 3,2,6,7/)
  case (6)  ! -X    +Z    E T
     outCHILD(1:2) =(/ 3,6/)
  case (7)  ! -X +Y -Z    ENB
     outCHILD(1) = 2
  case (8)  ! -X +Y       EN
     outCHILD(1:2) =(/ 3,2/)
  case (9)  ! -X +Y +Z    ENT
     outCHILD(1) = 3
  case (10) !    -Y -Z     SB
     outCHILD(1:2) =(/8,7/)
  case (11) !    -Y        S 
     outCHILD(1:4) =(/5,8,6,7/)
  case (12) !    -Y +Z     ST
     outCHILD(1:2) = (/5,6/)
  case (13) !       -Z      B
     outCHILD(1:4) =(/1,2,8,7/)
  case (14) !         
  case (15) !       +Z      T
     outCHILD(1:4) =(/4,3,5,6/)
  case (16) !    +Y -Z     NB
     outCHILD(1:2) =(/1,2/)
  case (17) !    +Y        N 
     outCHILD(1:4) =(/4,1,3,2/)
  case (18) !    +Y +Z     NT
     outCHILD(1:2) =(/4,3/)
  case (19) ! +X -Y -Z    WSB
     outCHILD(1) = 8
  case (20) ! +X -Y       WS 
     outCHILD(1:2) = (/5,8/)
  case (21) ! +X -Y +Z    WST
     outCHILD(1) = 5
  case (22) ! +X    -Z    W B 
     outCHILD(1:2) = (/1,8/)
  case (23) ! +X          W  
     outCHILD(1:4) =(/4,1,5,8/)
  case (24) ! +X    +Z    W T
     outCHILD(1:2) =(/4,5/)
  case (25) ! +X +Y -Z     WNB
     outCHILD(1) = 1
  case (26) ! +X +Y        WN 
     outCHILD(1:2) = (/4,1/)
  case (27) ! +X +Y +Z     WNT
     outCHILD(1) = 4
  end select
end subroutine get_children_list_orig
subroutine find_tree_neighbor_orig(inblkptr,outblkptr,iX,iY,iZ,noNeighbor)
  use ModOctree
  implicit none
  integer,intent(in)::iX,iY,iZ
  type (adaptive_block_ptr) :: inblkptr,outblkptr,tmp1blkptr
  logical, intent(inout) :: noNeighbor
  integer::iDir
  iDir=9*(iX+1)+3*(iY+1)+(iZ+1)+1
  call findTreeNeighbor_orig(inblkptr,outblkptr,idir,noNeighbor)
end subroutine find_tree_neighbor_orig

recursive subroutine findTreeNeighbor_orig(inblkptr,outblkptr,idir,noNeighbor)
  use ModParallel, ONLY : proc_dims,periodic3D
  use ModOctree
  implicit none

  integer, intent(in) :: idir
  type (adaptive_block_ptr) :: inblkptr,outblkptr,tmp1blkptr,tmp2blkptr
  logical, intent(inout) :: noNeighbor

  integer :: i,j,k,iFoundIt,idirNext,idirDown,iRoot,jRoot,kRoot,ii,ix,iy,iz
  integer, dimension(27,8) :: isFound,isNext,isDown

  !idir=
  !    1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
  !
  !   -X -X -X -X -X -X -X -X -X                            +X +X +X +X +X +X +X +X +X
  !   -Y -Y -Y          +Y +Y +Y -Y -Y -Y          +Y +Y +Y -Y -Y -Y          +Y +Y +Y
  !   -Z    +Z -Z    +Z -Z    +Z -Z    +Z -Z    +Z -Z    +Z -Z    +Z -Z    +Z -Z    +Z
  !
  data isFound / &
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 1, 0, 5, 8, 0, 0, 0, 0, 3, 2, 0, 6, 7, 0, &
       0, 0, 0, 4, 1, 0, 5, 8, 0, 0, 0, 0, 3, 2, 0, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
       0, 0, 0, 0, 4, 1, 0, 5, 8, 0, 0, 0, 0, 3, 2, 0, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 1, 0, 5, 8, 0, 0, 0, 0, 3, 2, 0, 6, 7, &
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 1, 0, 5, 8, 0, 0, 0, 0, 3, 2, 0, 6, 7, 0, 0, 0, &
       0, 4, 1, 0, 5, 8, 0, 0, 0, 0, 3, 2, 0, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
       4, 1, 0, 5, 8, 0, 0, 0, 0, 3, 2, 0, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
       0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 1, 0, 5, 8, 0, 0, 0, 0, 3, 2, 0, 6, 7, 0, 0, 0, 0   /
  data isNext / &
       2, 2, 3, 5, 5, 6, 5, 5, 6,11,11,12, 0, 0,15, 0, 0,15,11,11,12, 0, 0,15, 0, 0,15, &
      11,11,12, 0, 0,15, 0, 0,15,11,11,12, 0, 0,15, 0, 0,15,20,20,21,23,23,24,23,23,24, &
      10,11,11,13, 0, 0,13, 0, 0,10,11,11,13, 0, 0,13, 0, 0,19,20,20,22,23,23,22,23,23, &
       1, 2, 2, 4, 5, 5, 4, 5, 5,10,11,11,13, 0, 0,13, 0, 0,10,11,11,13, 0, 0,13, 0, 0, &
       4, 5, 5, 4, 5, 5, 7, 8, 8,13, 0, 0,13, 0, 0,16,17,17,13, 0, 0,13, 0, 0,16,17,17, &
      13, 0, 0,13, 0, 0,16,17,17,13, 0, 0,13, 0, 0,16,17,17,22,23,23,22,23,23,25,26,26, &
       0, 0,15, 0, 0,15,17,17,18, 0, 0,15, 0, 0,15,17,17,18,23,23,24,23,23,24,26,26,27, &
       5, 5, 6, 5, 5, 6, 8, 8, 9, 0, 0,15, 0, 0,15,17,17,18, 0, 0,15, 0, 0,15,17,17,18   /
  data isDown / &
       6, 7, 6, 3, 2, 3, 6, 7, 6, 5, 8, 5, 0, 0, 4, 0, 0, 5, 6, 7, 6, 0, 0, 3, 0, 0, 6, &
       5, 8, 5, 0, 0, 4, 0, 0, 5, 6, 7, 6, 0, 0, 3, 0, 0, 6, 5, 8, 5, 4, 1, 4, 5, 8, 5, &
       8, 5, 8, 1, 0, 0, 8, 0, 0, 7, 6, 7, 2, 0, 0, 7, 0, 0, 8, 5, 8, 1, 4, 1, 8, 5, 8, &
       7, 6, 7, 2, 3, 2, 7, 6, 7, 8, 5, 8, 1, 0, 0, 8, 0, 0, 7, 6, 7, 2, 0, 0, 7, 0, 0, &
       2, 3, 2, 7, 6, 7, 2, 3, 2, 1, 0, 0, 8, 0, 0, 1, 4, 1, 2, 0, 0, 7, 0, 0, 2, 3, 2, &
       1, 0, 0, 8, 0, 0, 1, 4, 1, 2, 0, 0, 7, 0, 0, 2, 3, 2, 1, 4, 1, 8, 5, 8, 1, 4, 1, &
       0, 0, 4, 0, 0, 5, 4, 1, 4, 0, 0, 3, 0, 0, 6, 3, 2, 3, 4, 1, 4, 5, 8, 5, 4, 1, 4, &
       3, 2, 3, 6, 7, 6, 3, 2, 3, 0, 0, 4, 0, 0, 5, 4, 1, 4, 0, 0, 3, 0, 0, 6, 3, 2, 3   /

  noNeighbor = .false.

  if (inblkptr%ptr%LEV .eq.0) then
     !Found root cell, find neighbor root cell
     do i=1,proc_dims(1)
        do j=1,proc_dims(2)
           do k=1,proc_dims(3)
              if (associated(inblkptr%ptr,octree_roots(i,j,k)%ptr)) then
                 iRoot=i
                 jRoot=j
                 kRoot=k
              end if
           end do
        end do
     end do
     !Break idir into components
     ii=idir-1
     ix=-1
     iy=-1
     iz=-1
     do i=1,2
        if (ii .ge. 9) then
           ix = ix + 1
           ii = ii - 9
        end if
     end do
     do i=1,2
        if (ii .ge. 3) then
           iy = iy + 1
           ii = ii - 3
        end if
     end do
     do i=1,2
        if (ii .ge. 1) then
           iz = iz + 1
           ii = ii - 1
        end if
     end do
     ix = ix + iRoot
     iy = iy + jRoot
     iz = iz + kRoot

     !\
     ! periodic boudary conditions
     !/
     if(periodic3D(1).and.ix==0)              ix=proc_dims(1)
     if(periodic3D(1).and.ix==proc_dims(1)+1) ix=1
     if(periodic3D(2).and.iy==0)              iy=proc_dims(2)
     if(periodic3D(2).and.iy==proc_dims(2)+1) iy=1
     if(periodic3D(3).and.iz==0)              iz=proc_dims(3)
     if(periodic3D(3).and.iz==proc_dims(3)+1) iz=1

     if ( ix .ge. 1 .and. ix .le.proc_dims(1) .and. &
          iy .ge. 1 .and. iy .le.proc_dims(2) .and. &
          iz .ge. 1 .and. iz .le.proc_dims(3) ) then
        outblkptr%ptr => octree_roots(ix,iy,iz)%ptr
     else
        noNeighbor = .true.
     end if
  else
     iFoundIt = isFound(idir,inblkptr%ptr%child_number)
     select case (iFoundIt)
     case (1,2,3,4,5,6,7,8)
        tmp1blkptr%ptr => inblkptr%ptr%parent%ptr
        
        outblkptr%ptr => tmp1blkptr%ptr%child(iFoundIt)%ptr
     case default
        idirNext = isNext(idir,inblkptr%ptr%child_number)
        idirDown = isDown(idir,inblkptr%ptr%child_number)
        tmp1blkptr%ptr => inblkptr%ptr%parent%ptr
        if (associated(tmp1blkptr%ptr)) then ! ensure block is allocated
           call findTreeNeighbor_orig(tmp1blkptr,tmp2blkptr,idirNext,noNeighbor)
           if (.not. noNeighbor) then
              if (associated(tmp2blkptr%ptr)) then ! ensure block is allocated
                 select case (idirDown)
                 case (1,2,3,4,5,6,7,8)
                    
                    outblkptr%ptr => tmp2blkptr%ptr%child(iDirDown)%ptr
                    if (.not. associated(outblkptr%ptr)) then
                       !neighbor must be at coarser level keep previous block
                       outblkptr%ptr => tmp2blkptr%ptr
                    end if
                 case default
                    write (*,*) 'ERROR: findTreeNeighbor: 4: no direction down tree '
                    !                    stop
                 end select
              else
                 write (*,*) 'ERROR: findTreeNeighbor: 3: not associated '
                 !                 stop
              end if
           end if
        else
           write (*,*) 'ERROR: findTreeNeighbor: 2: not associated - parent%ptr'
           !           stop
        end if
     end select
     if (.not. noNeighbor) then
        if (.not. associated(outblkptr%ptr)) then
           write (*,*) 'ERROR: findTreeNeighbor: 1: not associated '
           !           stop
        end if
     end if
  end if
end subroutine findTreeNeighbor_orig

