!^CFG COPYRIGHT UM
!//////////////////////////////////////////////////////////////////////////////////
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
  use ModProcMH
  use ModMain, ONLY : nBlockMax,unusedBLK
  use ModAMR
  use ModParallel
  use ModMpi
  implicit none

  integer :: i,j,k, inBLK,inPE, returnedLEV, iError
  integer, dimension(4) :: returnedPE,returnedBLK,returnedCHILD

  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  ! This block of code is done only for backward compatability with the
  !   old find_neighbor routine and the old neighbor variables needed.
  !/

  neiLtop = NOBLK
  neiPtop = NOBLK
  neiBtop = NOBLK
  neiLbot = NOBLK
  neiPbot = NOBLK
  neiBbot = NOBLK
  neiLeast = NOBLK
  neiPeast = NOBLK
  neiBeast = NOBLK
  neiLwest = NOBLK
  neiPwest = NOBLK
  neiBwest = NOBLK
  neiLnorth = NOBLK
  neiPnorth = NOBLK
  neiBnorth = NOBLK
  neiLsouth = NOBLK
  neiPsouth = NOBLK
  neiBsouth = NOBLK

  inPE = iProc
  do inBLK = 1, nBlockMax
     if (.not.unusedBLK(inBLK)) then
        !east
        call treeNeighbor(inPE,inBLK,-1, 0, 0,returnedPE,returnedBLK,returnedCHILD,returnedLEV)
        neiLeast(inBLK) = returnedLEV
        neiPeast(:,inBLK) = returnedPE
        neiBeast(:,inBLK) = returnedBLK

        !west
        call treeNeighbor(inPE,inBLK, 1, 0, 0,returnedPE,returnedBLK,returnedCHILD,returnedLEV)
        neiLwest(inBLK) = returnedLEV
        neiPwest(:,inBLK) = returnedPE
        neiBwest(:,inBLK) = returnedBLK

        !north
        call treeNeighbor(inPE,inBLK, 0, 1, 0,returnedPE,returnedBLK,returnedCHILD,returnedLEV)
        neiLnorth(inBLK) = returnedLEV
        neiPnorth(:,inBLK) = returnedPE
        neiBnorth(:,inBLK) = returnedBLK

        !south
        call treeNeighbor(inPE,inBLK, 0,-1, 0,returnedPE,returnedBLK,returnedCHILD,returnedLEV)
        neiLsouth(inBLK) = returnedLEV
        neiPsouth(:,inBLK) = returnedPE
        neiBsouth(:,inBLK) = returnedBLK

        !top
        call treeNeighbor(inPE,inBLK, 0, 0, 1,returnedPE,returnedBLK,returnedCHILD,returnedLEV)
        neiLtop(inBLK) = returnedLEV
        neiPtop(:,inBLK) = returnedPE
        neiBtop(:,inBLK) = returnedBLK

        !bottom
        call treeNeighbor(inPE,inBLK, 0, 0,-1,returnedPE,returnedBLK,returnedCHILD,returnedLEV)
        neiLbot(inBLK) = returnedLEV
        neiPbot(:,inBLK) = returnedPE
        neiBbot(:,inBLK) = returnedBLK

     end if
  end do

  neiLEV(east_,:) = neiLeast
  neiLEV(west_,:) = neiLwest
  neiLEV(south_,:)= neiLsouth
  neiLEV(north_,:)= neiLnorth
  neiLEV(bot_,:)  = neiLbot
  neiLEV(top_,:)  = neiLtop

  neiPE(:,1,:) = neiPeast
  neiBLK(:,1,:) = neiBeast
  neiPE(:,2,:) = neiPwest
  neiBLK(:,2,:) = neiBwest
  neiPE(:,3,:) = neiPsouth
  neiBLK(:,3,:) = neiBsouth
  neiPE(:,4,:) = neiPnorth
  neiBLK(:,4,:) = neiBnorth
  neiPE(:,5,:) = neiPbot
  neiBLK(:,5,:) = neiBbot
  neiPE(:,6,:) = neiPtop
  neiBLK(:,6,:) = neiBtop

  
  call MPI_ALLGATHER(unusedBLK,      nBLK, MPI_LOGICAL, &
       unusedBlock_BP, nBLK, MPI_LOGICAL, iComm, iError)

  call MPI_BARRIER(iComm, iError)

  !\
  ! This block of code is done only for backward compatability with the
  !   old find_neighbor routine and the old neighbor variables needed.
  !///////////////////////////////////////////////////////////////////////

  !\
  ! Loop and fill variables with neighbor information for all blocks on this PE
  !/
  inPE = iProc
  do inBLK = 1,nBlockMax
     if (.not.unusedBLK(inBLK)) then
        do i=-1,1
           do j=-1,1
              do k=-1,1
                 call treeNeighbor(inPE,inBLK,i,j,k, &
                      returnedPE,returnedBLK,returnedCHILD,returnedLEV)
                 BLKneighborPE   (i,j,k,:,inBLK) = returnedPE
                 BLKneighborBLK  (i,j,k,:,inBLK) = returnedBLK
                 BLKneighborCHILD(i,j,k,:,inBLK) = returnedCHILD
                 BLKneighborLEV  (i,j,k,  inBLK) = returnedLEV
              end do
           end do
        end do
     end if
  end do

end subroutine find_neighbors

subroutine treeNeighbor(inPE, inBLK, ix, iy, iz, outPE, outBLK, outCHILD, outLEV)
  use ModParallel, ONLY : NOBLK
  use ModOctree
  implicit none

  integer, intent(in) :: inPE,inBLK,ix,iy,iz
  integer, intent(out) :: outLEV
  integer, intent(out), dimension(4) :: outPE, outBLK, outCHILD
  integer :: i,idir
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

  idir = 1+ (iz+1) + 3*(iy+1) + 9*(ix+1)
  inBlockPtr%ptr => global_block_ptrs(inBLK,inPE+1)%ptr
  if (associated(inBlockPtr%ptr)) then ! ensure block is allocated
     call findTreeNeighbor(inBlockPtr,outBlockPtr,idir,noNeighbor)
     if (noNeighbor) then
        outLEV   = NOBLK
        outPE    = NOBLK
        outBLK   = NOBLK
        outCHILD = NOBLK
     else
        if (associated(outBlockPtr%ptr)) then ! ensure block is allocated
           if (.not. associated(outBlockPtr%ptr%child1)) then
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
              select case (idir)
              case (1)  ! -X -Y -Z    ESB
                 outPE (1)   = outBlockPtr%ptr%child7%PE
                 outBLK(1)   = outBlockPtr%ptr%child7%BLK
                 outCHILD(1) = outBlockPtr%ptr%child7%child_number
              case (2)  ! -X -Y       ES 
                 outPE (1)   = outBlockPtr%ptr%child6%PE
                 outBLK(1)   = outBlockPtr%ptr%child6%BLK
                 outCHILD(1) = outBlockPtr%ptr%child6%child_number

                 outPE (2)   = outBlockPtr%ptr%child7%PE
                 outBLK(2)   = outBlockPtr%ptr%child7%BLK
                 outCHILD(2) = outBlockPtr%ptr%child7%child_number
              case (3)  ! -X -Y +Z    EST
                 outPE (1)   = outBlockPtr%ptr%child6%PE
                 outBLK(1)   = outBlockPtr%ptr%child6%BLK
                 outCHILD(1) = outBlockPtr%ptr%child6%child_number
              case (4)  ! -X    -Z    E B
                 outPE (1)   = outBlockPtr%ptr%child2%PE
                 outBLK(1)   = outBlockPtr%ptr%child2%BLK
                 outCHILD(1) = outBlockPtr%ptr%child2%child_number

                 outPE (2)   = outBlockPtr%ptr%child7%PE
                 outBLK(2)   = outBlockPtr%ptr%child7%BLK
                 outCHILD(2) = outBlockPtr%ptr%child7%child_number
              case (5)  ! -X          E  
                 outPE (1)   = outBlockPtr%ptr%child3%PE
                 outBLK(1)   = outBlockPtr%ptr%child3%BLK
                 outCHILD(1) = outBlockPtr%ptr%child3%child_number

                 outPE (2)   = outBlockPtr%ptr%child2%PE
                 outBLK(2)   = outBlockPtr%ptr%child2%BLK
                 outCHILD(2) = outBlockPtr%ptr%child2%child_number

                 outPE (3)   = outBlockPtr%ptr%child6%PE
                 outBLK(3)   = outBlockPtr%ptr%child6%BLK
                 outCHILD(3) = outBlockPtr%ptr%child6%child_number

                 outPE (4)   = outBlockPtr%ptr%child7%PE
                 outBLK(4)   = outBlockPtr%ptr%child7%BLK
                 outCHILD(4) = outBlockPtr%ptr%child7%child_number
              case (6)  ! -X    +Z    E T
                 outPE (1)   = outBlockPtr%ptr%child3%PE
                 outBLK(1)   = outBlockPtr%ptr%child3%BLK
                 outCHILD(1) = outBlockPtr%ptr%child3%child_number

                 outPE (2)   = outBlockPtr%ptr%child6%PE
                 outBLK(2)   = outBlockPtr%ptr%child6%BLK
                 outCHILD(2) = outBlockPtr%ptr%child6%child_number
              case (7)  ! -X +Y -Z    ENB
                 outPE (1)   = outBlockPtr%ptr%child2%PE
                 outBLK(1)   = outBlockPtr%ptr%child2%BLK
                 outCHILD(1) = outBlockPtr%ptr%child2%child_number
              case (8)  ! -X +Y       EN
                 outPE (1)   = outBlockPtr%ptr%child3%PE
                 outBLK(1)   = outBlockPtr%ptr%child3%BLK
                 outCHILD(1) = outBlockPtr%ptr%child3%child_number

                 outPE (2)   = outBlockPtr%ptr%child2%PE
                 outBLK(2)   = outBlockPtr%ptr%child2%BLK
                 outCHILD(2) = outBlockPtr%ptr%child2%child_number
              case (9)  ! -X +Y +Z    ENT
                 outPE (1)   = outBlockPtr%ptr%child3%PE
                 outBLK(1)   = outBlockPtr%ptr%child3%BLK
                 outCHILD(1) = outBlockPtr%ptr%child3%child_number
              case (10) !    -Y -Z     SB
                 outPE (1)   = outBlockPtr%ptr%child8%PE
                 outBLK(1)   = outBlockPtr%ptr%child8%BLK
                 outCHILD(1) = outBlockPtr%ptr%child8%child_number

                 outPE (2)   = outBlockPtr%ptr%child7%PE
                 outBLK(2)   = outBlockPtr%ptr%child7%BLK
                 outCHILD(2) = outBlockPtr%ptr%child7%child_number
              case (11) !    -Y        S 
                 outPE (1)   = outBlockPtr%ptr%child5%PE
                 outBLK(1)   = outBlockPtr%ptr%child5%BLK
                 outCHILD(1) = outBlockPtr%ptr%child5%child_number

                 outPE (2)   = outBlockPtr%ptr%child8%PE
                 outBLK(2)   = outBlockPtr%ptr%child8%BLK
                 outCHILD(2) = outBlockPtr%ptr%child8%child_number

                 outPE (3)   = outBlockPtr%ptr%child6%PE
                 outBLK(3)   = outBlockPtr%ptr%child6%BLK
                 outCHILD(3) = outBlockPtr%ptr%child6%child_number

                 outPE (4)   = outBlockPtr%ptr%child7%PE
                 outBLK(4)   = outBlockPtr%ptr%child7%BLK
                 outCHILD(4) = outBlockPtr%ptr%child7%child_number
              case (12) !    -Y +Z     ST
                 outPE (1)   = outBlockPtr%ptr%child5%PE
                 outBLK(1)   = outBlockPtr%ptr%child5%BLK
                 outCHILD(1) = outBlockPtr%ptr%child5%child_number

                 outPE (2)   = outBlockPtr%ptr%child6%PE
                 outBLK(2)   = outBlockPtr%ptr%child6%BLK
                 outCHILD(2) = outBlockPtr%ptr%child6%child_number
              case (13) !       -Z      B
                 outPE (1)   = outBlockPtr%ptr%child1%PE
                 outBLK(1)   = outBlockPtr%ptr%child1%BLK
                 outCHILD(1) = outBlockPtr%ptr%child1%child_number

                 outPE (2)   = outBlockPtr%ptr%child2%PE
                 outBLK(2)   = outBlockPtr%ptr%child2%BLK
                 outCHILD(2) = outBlockPtr%ptr%child2%child_number

                 outPE (3)   = outBlockPtr%ptr%child8%PE
                 outBLK(3)   = outBlockPtr%ptr%child8%BLK
                 outCHILD(3) = outBlockPtr%ptr%child8%child_number

                 outPE (4)   = outBlockPtr%ptr%child7%PE
                 outBLK(4)   = outBlockPtr%ptr%child7%BLK
                 outCHILD(4) = outBlockPtr%ptr%child7%child_number
              case (14) !         
              case (15) !       +Z      T
                 outPE (1)   = outBlockPtr%ptr%child4%PE
                 outBLK(1)   = outBlockPtr%ptr%child4%BLK
                 outCHILD(1) = outBlockPtr%ptr%child4%child_number

                 outPE (2)   = outBlockPtr%ptr%child3%PE
                 outBLK(2)   = outBlockPtr%ptr%child3%BLK
                 outCHILD(2) = outBlockPtr%ptr%child3%child_number

                 outPE (3)   = outBlockPtr%ptr%child5%PE
                 outBLK(3)   = outBlockPtr%ptr%child5%BLK
                 outCHILD(3) = outBlockPtr%ptr%child5%child_number

                 outPE (4)   = outBlockPtr%ptr%child6%PE
                 outBLK(4)   = outBlockPtr%ptr%child6%BLK
                 outCHILD(4) = outBlockPtr%ptr%child6%child_number
              case (16) !    +Y -Z     NB
                 outPE (1)   = outBlockPtr%ptr%child1%PE
                 outBLK(1)   = outBlockPtr%ptr%child1%BLK
                 outCHILD(1) = outBlockPtr%ptr%child1%child_number

                 outPE (2)   = outBlockPtr%ptr%child2%PE
                 outBLK(2)   = outBlockPtr%ptr%child2%BLK
                 outCHILD(2) = outBlockPtr%ptr%child2%child_number
              case (17) !    +Y        N 
                 outPE (1)   = outBlockPtr%ptr%child4%PE
                 outBLK(1)   = outBlockPtr%ptr%child4%BLK
                 outCHILD(1) = outBlockPtr%ptr%child4%child_number

                 outPE (2)   = outBlockPtr%ptr%child1%PE
                 outBLK(2)   = outBlockPtr%ptr%child1%BLK
                 outCHILD(2) = outBlockPtr%ptr%child1%child_number

                 outPE (3)   = outBlockPtr%ptr%child3%PE
                 outBLK(3)   = outBlockPtr%ptr%child3%BLK
                 outCHILD(3) = outBlockPtr%ptr%child3%child_number

                 outPE (4)   = outBlockPtr%ptr%child2%PE
                 outBLK(4)   = outBlockPtr%ptr%child2%BLK
                 outCHILD(4) = outBlockPtr%ptr%child2%child_number
              case (18) !    +Y +Z     NT
                 outPE (1)   = outBlockPtr%ptr%child4%PE
                 outBLK(1)   = outBlockPtr%ptr%child4%BLK
                 outCHILD(1) = outBlockPtr%ptr%child4%child_number

                 outPE (2)   = outBlockPtr%ptr%child3%PE
                 outBLK(2)   = outBlockPtr%ptr%child3%BLK
                 outCHILD(2) = outBlockPtr%ptr%child3%child_number
              case (19) ! +X -Y -Z    WSB
                 outPE (1)   = outBlockPtr%ptr%child8%PE
                 outBLK(1)   = outBlockPtr%ptr%child8%BLK
                 outCHILD(1) = outBlockPtr%ptr%child8%child_number
              case (20) ! +X -Y       WS 
                 outPE (1)   = outBlockPtr%ptr%child5%PE
                 outBLK(1)   = outBlockPtr%ptr%child5%BLK
                 outCHILD(1) = outBlockPtr%ptr%child5%child_number

                 outPE (2)   = outBlockPtr%ptr%child8%PE
                 outBLK(2)   = outBlockPtr%ptr%child8%BLK
                 outCHILD(2) = outBlockPtr%ptr%child8%child_number
              case (21) ! +X -Y +Z    WST
                 outPE (1)   = outBlockPtr%ptr%child5%PE
                 outBLK(1)   = outBlockPtr%ptr%child5%BLK
                 outCHILD(1) = outBlockPtr%ptr%child5%child_number
              case (22) ! +X    -Z    W B 
                 outPE (1)   = outBlockPtr%ptr%child1%PE
                 outBLK(1)   = outBlockPtr%ptr%child1%BLK
                 outCHILD(1) = outBlockPtr%ptr%child1%child_number

                 outPE (2)   = outBlockPtr%ptr%child8%PE
                 outBLK(2)   = outBlockPtr%ptr%child8%BLK
                 outCHILD(2) = outBlockPtr%ptr%child8%child_number
              case (23) ! +X          W  
                 outPE (1)   = outBlockPtr%ptr%child4%PE
                 outBLK(1)   = outBlockPtr%ptr%child4%BLK
                 outCHILD(1) = outBlockPtr%ptr%child4%child_number

                 outPE (2)   = outBlockPtr%ptr%child1%PE
                 outBLK(2)   = outBlockPtr%ptr%child1%BLK
                 outCHILD(2) = outBlockPtr%ptr%child1%child_number

                 outPE (3)   = outBlockPtr%ptr%child5%PE
                 outBLK(3)   = outBlockPtr%ptr%child5%BLK
                 outCHILD(3) = outBlockPtr%ptr%child5%child_number

                 outPE (4)   = outBlockPtr%ptr%child8%PE
                 outBLK(4)   = outBlockPtr%ptr%child8%BLK
                 outCHILD(4) = outBlockPtr%ptr%child8%child_number
              case (24) ! +X    +Z    W T
                 outPE (1)   = outBlockPtr%ptr%child4%PE
                 outBLK(1)   = outBlockPtr%ptr%child4%BLK
                 outCHILD(1) = outBlockPtr%ptr%child4%child_number

                 outPE (2)   = outBlockPtr%ptr%child5%PE
                 outBLK(2)   = outBlockPtr%ptr%child5%BLK
                 outCHILD(2) = outBlockPtr%ptr%child5%child_number
              case (25) ! +X +Y -Z     WNB
                 outPE (1)   = outBlockPtr%ptr%child1%PE
                 outBLK(1)   = outBlockPtr%ptr%child1%BLK
                 outCHILD(1) = outBlockPtr%ptr%child1%child_number
              case (26) ! +X +Y        WN 
                 outPE (1)   = outBlockPtr%ptr%child4%PE
                 outBLK(1)   = outBlockPtr%ptr%child4%BLK
                 outCHILD(1) = outBlockPtr%ptr%child4%child_number

                 outPE (2)   = outBlockPtr%ptr%child1%PE
                 outBLK(2)   = outBlockPtr%ptr%child1%BLK
                 outCHILD(2) = outBlockPtr%ptr%child1%child_number
              case (27) ! +X +Y +Z     WNT
                 outPE (1)   = outBlockPtr%ptr%child4%PE
                 outBLK(1)   = outBlockPtr%ptr%child4%BLK
                 outCHILD(1) = outBlockPtr%ptr%child4%child_number
              end select
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

recursive subroutine findTreeNeighbor(inblkptr,outblkptr,idir,noNeighbor)
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
        tmp1blkptr%ptr => inblkptr%ptr%parent
        call neighborAssignment(tmp1blkptr,outblkptr,iFoundIt)
     case default
        idirNext = isNext(idir,inblkptr%ptr%child_number)
        idirDown = isDown(idir,inblkptr%ptr%child_number)
        tmp1blkptr%ptr => inblkptr%ptr%parent
        if (associated(tmp1blkptr%ptr)) then ! ensure block is allocated
           call findTreeNeighbor(tmp1blkptr,tmp2blkptr,idirNext,noNeighbor)
           if (.not. noNeighbor) then
              if (associated(tmp2blkptr%ptr)) then ! ensure block is allocated
                 select case (idirDown)
                 case (1,2,3,4,5,6,7,8)
                    call neighborAssignment(tmp2blkptr,outblkptr,idirDown)
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
           write (*,*) 'ERROR: findTreeNeighbor: 2: not associated - parent'
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
end subroutine findTreeNeighbor

subroutine neighborAssignment(inblkptr,outblkptr,idir)
  use ModOctree
  implicit none

  integer, intent(in) :: idir
  type (adaptive_block_ptr) :: inblkptr,outblkptr

  select case (idir)
  case (1)
     outblkptr%ptr => inblkptr%ptr%child1
  case (2)
     outblkptr%ptr => inblkptr%ptr%child2
  case (3)
     outblkptr%ptr => inblkptr%ptr%child3
  case (4)
     outblkptr%ptr => inblkptr%ptr%child4
  case (5)
     outblkptr%ptr => inblkptr%ptr%child5
  case (6)
     outblkptr%ptr => inblkptr%ptr%child6
  case (7)
     outblkptr%ptr => inblkptr%ptr%child7
  case (8)
     outblkptr%ptr => inblkptr%ptr%child8
  end select
end subroutine neighborAssignment

subroutine fixRefinementLevels
  use ModProcMH
  use ModMain, ONLY : nBLK,lVerbose
  use ModAMR, ONLY : refine_list
  use ModGeometry,ONLY: TypeGeometry                !^CFG IF NOT CARTESIAN
  use ModOctree
  implicit none

  integer :: i,j,k, idir, inPE,inBLK, iLEV1,iLEV2, iCount, maxLev, curLev
  logical :: noNeighbor
  type (adaptive_block_ptr) :: inBlockPtr,outBlockPtr,tmpBlockPtr

  call set_body_flag
  do
     iCount = 0
     refine_list = .false.
     !find max level
     maxLev=0
     do inPE = 1,nProc
        do inBLK = 1,nBLK
           inBlockPtr%ptr => global_block_ptrs(inBLK,inPE)%ptr
           if (associated(inBlockPtr%ptr)) &
                maxLev = max(maxLev, inBlockPtr%ptr%LEV)
        end do
     end do
     !loop from max level down to zero
     do curLev = maxLev,1,-1
        do inPE=1,nProc; do inBLK=1,nBLK
           inBlockPtr%ptr => global_block_ptrs(inBLK,inPE)%ptr
           if (associated(inBlockPtr%ptr)) then
              if (inBlockPtr%ptr%LEV == curLev) then
                 iLEV1 = inBlockPtr%ptr%LEV
                 if (refine_list(inBlockPtr%ptr%BLK,inBlockPtr%ptr%PE+1)) &
                      iLEV1 = iLEV1 + 1
                 do i=-1,1; do j=-1,1; do k=-1,1
                    idir = 1+ (k+1) + 3*(j+1) + 9*(i+1)
                    call findTreeNeighbor(inBlockPtr,outBlockPtr,idir,noNeighbor)
                    if (.not. noNeighbor)then
                       if(associated(outBlockPtr%ptr)) then
                          if (.not. associated(outBlockPtr%ptr%child1)) then
                             iLEV2 = outBlockPtr%ptr%LEV
                             if (refine_list(outBlockPtr%ptr%BLK,outBlockPtr%ptr%PE+1)) &
                                  iLEV2 = iLEV2 + 1
                             if ( (iLEV1-iLEV2 > 1) .or. &
                                  (iLEV1-iLEV2 > 0 .and. &
                                  ((inBlockPtr%ptr%body .and. outBlockPtr%ptr%body)&
                                  .or.(inBlockPtr%ptr%IsOuterBoundary .and. &    !^CFG IF FACEOUTERBC BEGIN
                                  outBlockPtr%ptr%IsOuterBoundary)&
                                  .or.(inBlockPtr%ptr%IsExtraBoundary .and.& 
                                  outBlockPtr%ptr%IsExtraBoundary&
                                  .and.( i==0.and.k==0.or.TypeGeometry=='cartesian')& !^CFG IF NOT CARTESIAN
                                  )&           !^CFG END FACEOUTERBC
                                  ))) then
                                refine_list(outBlockPtr%ptr%BLK,outBlockPtr%ptr%PE+1) = .true.
                                iCount = iCount +1
                             end if
                          end if
                       end if
                    end if
                 end do; end do; end do
              end if
           end if
        end do; end do
     end do
     if (iCount == 0) EXIT
     if (iProc == 0.and.lVerbose>0)&
          write (*,*) '    FixRefinementLevels found ',iCount,' blocks to refine'
     call parallel_refine
     call set_body_flag
  end do

end subroutine fixRefinementLevels
!^CFG IF NOT CARTESIAN BEGIN
subroutine find_axial_neighbor(iPEIn,iBLKIn,iPEOut,iBLKout)
  use ModParallel, ONLY : proc_dims
  use ModOctree
  implicit none
  integer,intent(in)::iPEIn,iBLKIn
  integer,intent(out)::iPEOut,iBLKOut
  integer,dimension(99)::iChild_I
  integer::i,j,k,iRoot,jRoot,kRoot,iLevel,iLevelIn
  type (adaptive_block_ptr) :: InBlkPtr,OutBlkPtr
 
  nullify (InBlkPtr % ptr)
  InBlkPtr%ptr=> global_block_ptrs(iBlkIn,iPEIn+1)%ptr
  iLevelIn=InBlkPtr%ptr%LEV

  do iLevel=iLevelIn,1,-1
     iChild_I(iLevel)=InBlkPtr%ptr%child_number
     InBlkPtr%ptr=>InBlkPtr%ptr%parent
  end do
     !Found root cell, find neighbor root cell
  do i=1,proc_dims(1)
     do j=1,proc_dims(2)
        do k=1,proc_dims(3)
           if (associated(InBlkPtr%ptr,octree_roots(i,j,k)%ptr)) then
              iRoot=i
              jRoot=j
              kRoot=k
           end if
        end do
     end do
  end do
  jRoot=1+mod(jRoot-1+proc_dims(2)/2,proc_dims(2))
  
  nullify (OutBlkPtr%ptr)
  
  OutBlkPtr%ptr=>octree_roots(iRoot,jRoot,kRoot)%ptr
  do iLevel=1,iLevelIn
     call  neighborAssignment(OutBlkPtr,OutBlkPtr,iChild_I(iLevel))
     if (.not.associated(OutBlkPtr%ptr))&
          call stop_mpi('Wrong Axial Neighbor is found') ! ensure block is allocated
  end do
     iPEOut  = outBlkPtr%ptr%PE
     iBlkOut = outBlkPtr%ptr%BLK
end subroutine find_axial_neighbor
!^CFG END CARTESIAN
