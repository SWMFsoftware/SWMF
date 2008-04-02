module ModOctreeNew

  implicit none
  save

  private ! except

  public:: init_mod_octree
  public:: set_root_block
  public:: refine_block
  public:: coarsen_block
  public:: test_octree

  integer, parameter :: nDim = 3
  integer, parameter :: nChild = 2**nDim

  integer, allocatable :: iOctree_IA(:,:)

  integer, parameter :: &
       Status_   =  1, &
       Parent_   =  2, &
       Child0_   =  2, &
       Child1_   =  3, &
       Child2_   =  4, &
       Child3_   =  5, &
       Child4_   =  6, &
       Child5_   =  7, &
       Child6_   =  8, &
       Child7_   =  9, &
       Child8_   = 10, &
       Proc_     = 11, &
       Level_    = 12, &
       LevelMin_ = 13, &
       LevelMax_ = 14, &
       iCoord_   = 15, &
       jCoord_   = 16, &
       kCoord_   = 17, &
       iLeft_    = 18, &
       iRight_   = 19, &
       jLeft_    = 20, &
       jRight_   = 21, &
       kLeft_    = 22, &
       kRight_   = 23

  integer, parameter :: nInfo = 23

  ! Possible values for the status variable
  integer, parameter :: Skipped_=0, Unused_=1, Used_=2, Refine_=3, Coarsen_=4

  ! Index for non-existing neighbors
  integer, parameter :: NoBlock_=0

  ! Opposite directions
  integer, parameter :: jNeighbor_I(iLeft_:kRight_) = &
       (/ iRight_, iLeft_, jRight_, jLeft_, kRight_, kLeft_ /)

  integer :: MaxBlockAll
  integer :: nRoot_D(nDim)
  logical :: IsPeriodic_D(nDim)
  integer, allocatable :: iBlockRoot_III(:,:,:)

contains

  subroutine init_mod_octree(nBlock)

    ! Initialize the octree array with nBlock blocks

    integer, intent(in) :: nBlock
    !----------------------------------------------------------------------
    if(allocated(iOctree_IA)) RETURN

    MaxBlockAll = nBlock
    allocate(iOctree_IA(nInfo, MaxBlockAll))

    ! Initialize all elements and make blocks skipped
    iOctree_IA = Skipped_

  end subroutine init_mod_octree

  !==========================================================================

  integer function i_block_new()

    ! Find a skipped element in the iOctree_IA array

    integer :: iBlock

    do iBlock = 1, MaxBlockAll
       if(iOctree_IA(Status_, iBlock) == Skipped_)then
          i_block_new = iBlock
          return
       end if
    end do
    ! Could not find any skipped block
    i_block_new = -1

  end function i_block_new

  !==========================================================================

  subroutine set_root_block(nRootIn_D, IsPeriodicIn_D)

    integer, intent(in) :: nRootIn_D(nDim)
    logical, intent(in) :: IsPeriodicIn_D(nDim)

    integer :: iRoot, jRoot, kRoot, iBlock
    integer :: iRootMax, jRootMax, kRootMax
    !-----------------------------------------------------------------------

    if(allocated(iBlockRoot_III)) RETURN

    nRoot_D      = nRootIn_D
    IsPeriodic_D = IsPeriodicIn_D
    iRootMax = nRoot_D(1); jRootMax = nRoot_D(2); kRootMax = nRoot_D(3)

    allocate(iBlockRoot_III(iRootMax, jRootMax, kRootMax))

    do kRoot = 1, kRootMax
       do jRoot = 1, jRootMax
          do iRoot = 1, iRootMax
             iBlock = i_block_new()

             iBlockRoot_III(iRoot, jRoot, kRoot) = iBlock

             iOctree_IA(Status_, iBlock)         = Used_
             iOctree_IA(Parent_, iBlock)         = NoBlock_
             iOctree_IA(Child1_:Child8_, iBlock) = NoBlock_
             iOctree_IA(Level_ , iBlock)         = 1
             iOctree_IA(iCoord_, iBlock)         = iRoot
             iOctree_IA(jCoord_, iBlock)         = jRoot
             iOctree_IA(kCoord_, iBlock)         = kRoot

          end do
       end do
    end do

    ! Set neighbor info
    do kRoot = 1, kRootMax; do jRoot = 1, jRootMax; do iRoot = 1, iRootMax
       iBlock = iBlockRoot_III(iRoot, jRoot, kRoot)

       if(iRoot > 1)then
          iOctree_IA(iLeft_, iBlock) = iBlockRoot_III(iRoot-1, jRoot, kRoot)
       elseif(IsPeriodic_D(1))then
          iOctree_IA(iLeft_, iBlock) = iBlockRoot_III(iRootMax, jRoot, kRoot)
       else
          iOctree_IA(iLeft_, iBlock) = NoBlock_
       end if

       if(iRoot < iRootMax)then
          iOctree_IA(iRight_, iBlock) = iBlockRoot_III(iRoot+1, jRoot, kRoot)
       elseif(IsPeriodic_D(1))then
          iOctree_IA(iRight_, iBlock) = iBlockRoot_III(1, jRoot, kRoot)
       else
          iOctree_IA(iRight_, iBlock) = NoBlock_
       end if

       if(jRoot > 1)then
          iOctree_IA(jLeft_, iBlock) = iBlockRoot_III(iRoot, jRoot-1, kRoot)
       elseif(IsPeriodic_D(2))then
          iOctree_IA(jLeft_, iBlock) = iBlockRoot_III(iRoot, jRootMax, kRoot)
       else
          iOctree_IA(jLeft_, iBlock) = NoBlock_
       end if

       if(jRoot < jRootMax)then
          iOctree_IA(jRight_, iBlock) = iBlockRoot_III(iRoot, jRoot+1, kRoot)
       elseif(IsPeriodic_D(2))then
          iOctree_IA(jRight_, iBlock) = iBlockRoot_III(iRoot, 1, kRoot)
       else
          iOctree_IA(jRight_, iBlock) = NoBlock_
       end if

       if(kRoot > 1)then
          iOctree_IA(kLeft_, iBlock) = iBlockRoot_III(iRoot, jRoot, kRoot-1)
       elseif(IsPeriodic_D(3))then
          iOctree_IA(kLeft_, iBlock) = iBlockRoot_III(iRoot, jRoot, kRootMax)
       else
          iOctree_IA(kLeft_, iBlock) = NoBlock_
       end if

       if(kRoot < kRootMax)then
          iOctree_IA(kRight_, iBlock) = iBlockRoot_III(iRoot, jRoot, kRoot+1)
       elseif(IsPeriodic_D(3))then
          iOctree_IA(kRight_, iBlock) = iBlockRoot_III(iRoot, jRoot, 1)
       else
          iOctree_IA(kRight_, iBlock) = NoBlock_
       end if

    end do; end do; end do

  end subroutine set_root_block

  !==========================================================================
  subroutine refine_block(iBlock)

    integer, intent(in) :: iBlock

    integer :: iShift, jShift, kShift, iChild
    integer :: iBlockChild, iBlockChild_III(2,2,2)
    integer :: iSide, jSide, jBlock, jChild, jBlockChild
    !----------------------------------------------------------------------

    iOctree_IA(Status_, iBlock) = Unused_

    iChild = Child0_
    do kShift = 1,2; do jShift=1,2; do iShift=1,2
       iChild = iChild+1

       iBlockChild = i_block_new()
       iBlockChild_III(iShift, jShift, kShift) = iBlockChild

       iOctree_IA(iChild, iBlock)               = iBlockChild

       iOctree_IA(Status_, iBlockChild)         = Used_

       iOctree_IA(Parent_, iBlockChild)         = iBlock

       iOctree_IA(Child1_:Child8_, iBlockChild) = NoBlock_

       iOctree_IA(Level_, iBlockChild)          = &
            iOctree_IA(Level_, iBlock) + 1

       iOctree_IA(LevelMin_:LevelMax_, iBlockChild) = &
            iOctree_IA(LevelMin_:LevelMax_, iBlock)

       ! Calculate the coordinates of the child block
       iOctree_IA(iCoord_:kCoord_, iBlockChild) = &
            2*(iOctree_IA(iCoord_:kCoord_, iBlock) - 1) &
            + (/ iShift, jShift, kShift /)
          
    end do; end do; end do

    ! Figure out neighbors of children
    do kShift = 1,2; do jShift=1,2; do iShift=1,2

       iBlockChild = iBlockChild_III(iShift, jShift, kShift)

       ! Do i-directional neighbors first
       iSide = iLeft_ + iShift - 1; jSide = iRight_ + 1 - iShift

       ! One of the left/right neighbors is another child of the same parent
       jBlockChild = iBlockChild_III(3-iShift, jShift, kShift)
       iOctree_IA(jSide, iBlockChild) = jBlockChild
       iOctree_IA(iSide, jBlockChild) = iBlockChild

       ! The other neighbor is the appropriate child of parent's neighbor
       jBlock = iOctree_IA(iSide, iBlock)
       if(jBlock /= NoBlock_)then

          jChild      = Child1_ + (2-iShift) + (jShift-1)*2 + (kShift-1)*4
          jBlockChild = iOctree_IA(jChild, jBlock)

          iOctree_IA(iSide, iBlockChild) = jBlockChild
          if(jBlockChild /= NoBlock_) &
               iOctree_IA(jSide, jBlockChild) = iBlockChild
       else
          iOctree_IA(iSide, iBlockChild) = NoBlock_
       end if

       ! Do j-directional neighbors
       iSide = jLeft_ + jShift - 1; jSide = jRight_ + 1 - jShift

       ! One of the left/right neighbors is another child of the same parent
       jBlockChild = iBlockChild_III(iShift, 3-jShift, kShift)
       iOctree_IA(jSide, iBlockChild) = jBlockChild
       iOctree_IA(iSide, jBlockChild) = iBlockChild

       ! The other neighbor is the appropriate child of parent's neighbor
       jBlock = iOctree_IA(iSide, iBlock)
       if(jBlock /= NoBlock_)then
          jChild      = Child1_ + iShift-1 + (2-jShift)*2 + (kShift-1)*4
          jBlockChild = iOctree_IA(jChild, jBlock)
          iOctree_IA(iSide, iBlockChild) = jBlockChild
          if(jBlockChild /= NoBlock_) &
               iOctree_IA(jSide, jBlockChild) = iBlockChild
       else
          iOctree_IA(iSide, iBlockChild) = NoBlock_
       end if

       ! Do k-directional neighbors
       iSide = kLeft_ + kShift - 1; jSide = kRight_ + 1 - kShift

       ! One of the left/right neighbors is another child of the same parent
       jBlockChild = iBlockChild_III(iShift, jShift, 3-kShift)
       iOctree_IA(jSide, iBlockChild) = jBlockChild
       iOctree_IA(iSide, jBlockChild) = iBlockChild

       ! The other neighbor is the appropriate child of parent's neighbor
       jBlock = iOctree_IA(iSide, iBlock)
       if(jBlock /= NoBlock_)then
          jChild      = Child1_ + iShift-1 + (jShift-1)*2 + (2-kShift)*4
          jBlockChild = iOctree_IA(jChild, jBlock)
          iOctree_IA(iSide, iBlockChild) = jBlockChild
          if(jBlockChild /= NoBlock_) &
               iOctree_IA(jSide, jBlockChild) = iBlockChild
       else
          iOctree_IA(iSide, iBlockChild) = NoBlock_
       end if

    end do; end do; end do

  end subroutine refine_block

  !==========================================================================
  subroutine coarsen_block(iBlock)
    integer, intent(in) :: iBlock

    integer :: iChild, iBlockChild, iNeighbor, jBlockChild
    !-----------------------------------------------------------------------

    do iChild = Child1_, Child8_
       iBlockChild = iOctree_IA(iChild, iBlock)

       ! Tell neighbors that this block does not exist any longer
       do iNeighbor = iLeft_, kRight_
          jBlockChild = iOctree_IA(iNeighbor, iBlockChild)
          if(jBlockChild /= NoBlock_) &
               iOctree_IA(jNeighbor_I(iNeighbor), jBlockChild) = NoBlock_
       end do
       ! Wipe out the child block
       iOctree_IA(:, Status_) = Skipped_
    end do

    ! Make block used with no children
    iOctree_IA(Status_, iBlock) = Used_
    iOctree_IA(Child1_:Child8_, iBlock) = NoBlock_

  end subroutine coarsen_block

  !==========================================================================

  subroutine test_octree

    integer :: iBlock
    character(len=*), parameter :: NameSub = 'test_octree'
    !-----------------------------------------------------------------------
    write(*,*)'Testing init_mod_octree'
    call init_mod_octree(100)
    if(MaxBlockAll /= 100) &
         write(*,*)'init_mod_octtree faild, MaxBlockAll=',&
         MaxBlockAll, ' should be 100'

    write(*,*)'Testing i_block_new()'
    iBlock = i_block_new()
    if(iBlock /= 1) &
         write(*,*)'i_block_new() failed, iBlock=',iBlock,' should be 1'

    write(*,*)'Testing set_root_block'
    call set_root_block( (/1,2,3/), (/.true., .true., .false./) )

    if(any( nRoot_D /= (/1,2,3/) )) &
         write(*,*) 'set_root_block failed, nRoot_D=',nRoot_D,&
         ' should be 1,2,3'

    if(any( iBlockRoot_III /= reshape( (/1,2,3,4,5,6/), (/1,2,3/) ) )) &
         write(*,*) 'set_root_block failed, iBlockRoot_III=',iBlockRoot_III,&
         ' should be 1,2,3,4,5,6'

    if(any( iOctree_IA(iCoord_:kCoord_,4) /= (/1,2,2/) )) &
         write(*,*) 'set_root_block failed, coordinates of block four=',&
         iOctree_IA(iCoord_:kCoord_,4), 'should be 1,2,2'

    if(any( iOctree_IA(iLeft_:kRight_,5) /= (/5,5,6,6,3,0/) )) &
         write(*,*) 'set_root_block failed, neighbors of block four=',&
         iOctree_IA(iLeft_:kRight_,4), 'should be 5,5,6,6,3,0'

    ! Coordinates of root blocks
    !do iBlock = 1,6
    !   write(*,*)iBlock,iOctree_IA(iCoord_:kCoord_,iBlock)
    !end do
    ! Neighbors of root blocks
    !do iBlock = 1,6
    !   write(*,*)iBlock,iOctree_IA(iLeft_:kRight_,iBlock)
    !end do

    write(*,*)'Testing refine_block'
    call refine_block(4)

    write(*,*)'iOctree_IA(:,4)=',iOctree_IA(:,4)
    write(*,*)'iOctree_IA(:,7)=',iOctree_IA(:,7)

    call refine_block(3)

    write(*,*)'iOctree_IA(:,3)=',iOctree_IA(:,3)
    write(*,*)'iOctree_IA(:,7)=',iOctree_IA(:,7)

    write(*,*)'Testing coarsen_block'
    call coarsen_block(3)

    write(*,*)'iOctree_IA(:,3)=',iOctree_IA(:,3)
    write(*,*)'iOctree_IA(:,7)=',iOctree_IA(:,7)

  end subroutine test_octree

end module ModOctreeNew
