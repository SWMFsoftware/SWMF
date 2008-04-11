module ModOctreeNew

  implicit none
  save

  private ! except

  public:: init_mod_octree
  public:: set_root_block
  public:: refine_block
  public:: coarsen_block
  public:: find_point
  public:: test_octree

  integer, parameter :: MaxDim = 3
  integer, parameter :: nDim   = 3
  integer, parameter :: nChild = 2**nDim

  integer, allocatable :: iOctree_IA(:,:)

  integer, parameter :: &
       Status_   = 1, &
       Level_    = 2, &
       LevelMin_ = 3, &
       LevelMax_ = 4, &
       Parent_   = 5, &
       Child0_   = 6, &
       Child1_   = Child0_ + 1,      &
       ChildLast_= Child0_ + nChild, &
       Proc_     = Child0_ + 1,      & ! Overlaps with child 1
       Block_    = Child0_ + 2,      & ! Overlaps with child 2
       Coord0_   = ChildLast_,       &
       Coord1_   = Coord0_ + 1,      &
       CoordLast_= Coord0_ + nDim

  integer, parameter :: nInfo = CoordLast_

  ! Possible values for the status variable
  integer, parameter :: Skipped_=0, Unused_=1, Used_=2, Refine_=3, Coarsen_=4

  ! Index for non-existing block and level differences
  integer, parameter :: NoBlock_ = -100

  ! Neighbor information
  integer, allocatable :: DiLevelNei_IIIB(:,:,:,:), iBlockNei_IIIB(:,:,:,:)

  ! Deepest AMR level relative to root blocks (limited by 32 bit integers)
  integer, parameter :: MaxLevel = 30

  ! The maximum integer coordinate for a given level below root blocks
  ! The loop variable has to be declared to work-around NAG f95 bug
  integer :: L__
  integer, parameter :: MaxCoord_I(0:MaxLevel) = (/ (2**L__, L__=0,MaxLevel) /)

  ! Maximum number of blocks including unused and skipped ones
  integer :: MaxBlockAll = 0

  ! Maximum number of blocks per processor
  integer :: MaxBlock = 0

  ! Number of levels below root in level (that has occured at any time)
  integer :: nLevel = 0

  ! The number of root blocks in all dimensions, and altogether
  integer :: nRoot_D(MaxDim) = 0, nRoot = 0

  ! Periodicity of the domain per dimension
  logical :: IsPeriodic_D(MaxDim) = .false.

  ! Cylindrical or spherical coordinates
  logical :: IsSpherical = .false., IsCylindrical = .false.

contains

  subroutine init_mod_octree(nBlockProc, nBlockAll)

    ! Initialize the octree array with nBlock blocks

    integer, intent(in) :: nBlockProc ! Max number of blocks per processor
    integer, intent(in) :: nBlockAll  ! Max number of blocks altogether
    !----------------------------------------------------------------------
    if(allocated(iOctree_IA)) RETURN

    MaxBlockAll = nBlockAll
    allocate(iOctree_IA(nInfo, MaxBlockAll))

    ! Initialize all elements and make blocks skipped
    iOctree_IA = Skipped_

    MaxBlock = nBlockProc
    allocate(iBlockNei_IIIB(0:3,0:3,0:3,MaxBlock))
    allocate(DiLevelNei_IIIB(-1:1,-1:1,-1:1,MaxBlock))

    ! Initialize all elements and make neighbors unknown
    iBlockNei_IIIB  = NoBlock_
    DiLevelNei_IIIB = NoBlock_

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

    integer, intent(in) :: nRootIn_D(MaxDim)
    logical, intent(in) :: IsPeriodicIn_D(MaxDim)

    integer :: iRoot, jRoot, kRoot, iBlock, Ijk_D(MaxDim)
    !-----------------------------------------------------------------------

    nRoot_D      = nRootIn_D
    nRoot        = product(nRoot_D)
    IsPeriodic_D = IsPeriodicIn_D

    ! Use the first product(nRoot_D) blocks as root blocks in the octree
    iBlock = 0
    do kRoot = 1, nRoot_D(3)
       do jRoot = 1, nRoot_D(2)
          do iRoot = 1, nRoot_D(1)

             Ijk_D = (/ iRoot, jRoot, kRoot /)

             iBlock = iBlock + 1
             iOctree_IA(Status_, iBlock)            = Used_
             iOctree_IA(Parent_, iBlock)            = NoBlock_
             iOctree_IA(Child1_:ChildLast_, iBlock) = NoBlock_
             iOctree_IA(Level_ , iBlock)            = 0
             iOctree_IA(Coord1_:CoordLast_, iBlock) = Ijk_D(1:nDim)

          end do
       end do
    end do

    ! Set neighbor info
    !do iBlock = 1, nRoot
    !   call find_neighbors(iBlock)
    !end do

  end subroutine set_root_block

  !==========================================================================
  subroutine refine_block(iBlock)

    integer, intent(in) :: iBlock

    integer :: iChild, DiChild, iLevelChild, iProc, iBlockProc, iCoord_D(nDim)
    integer :: iDim, iBlockChild
    !----------------------------------------------------------------------

    iOctree_IA(Status_, iBlock) = Unused_

    iLevelChild = iOctree_IA(Level_, iBlock) + 1
    iProc       = iOctree_IA(Proc_,  iBlock)
    iBlockProc  = iOctree_IA(Block_, iBlock)

    ! Keep track of number of levels
    nLevel = max(nLevel, iLevelChild)
    if(nLevel > MaxLevel) &
         call CON_stop('Error in refine_block: too many levels')

    iCoord_D = 2*iOctree_IA(Coord1_:CoordLast_, iBlock) - 1

    do iChild = Child1_, ChildLast_

       iBlockChild = i_block_new()

       iOctree_IA(iChild, iBlock) = iBlockChild

       iOctree_IA(Status_,   iBlockChild) = Used_
       iOctree_IA(Level_,    iBlockChild) = iLevelChild
       iOctree_IA(LevelMin_, iBlockChild) = iOctree_IA(LevelMin_, iBlock)
       iOctree_IA(LevelMax_, iBlockChild) = iOctree_IA(LevelMax_, iBlock)
       iOctree_IA(Parent_,   iBlockChild) = iBlock
       iOctree_IA(Child1_:ChildLast_, iBlockChild) = NoBlock_

       ! This overwrites the two first children (saves memory)
       iOctree_IA(Proc_,     iBlockChild) = iProc
       iOctree_IA(Block_,    iBlockChild) = iBlockProc

       ! Calculate the coordinates of the child block
       DiChild = iChild - Child1_
       do iDim = 1, nDim
          iOctree_IA(Coord0_+iDim, iBlockChild) = &
               iCoord_D(iDim) + ibits(DiChild, iDim-1, 1)
       end do

    end do

    ! Find neighbors of children
    !do iChild = Child1_, ChildLast_
    !
    !   iBlockChild = iOctree_IA(iChild, iBlock)
    !   call find_neighbors(iBlockChild)
    !
    !end do

    ! Should also redo neighbors of the parent block

  end subroutine refine_block

  !==========================================================================
  subroutine coarsen_block(iBlock)

    integer, intent(in) :: iBlock

    integer :: iChild, iBlockChild1, iBlockChild
    !-----------------------------------------------------------------------

    do iChild = Child1_, ChildLast_
       iBlockChild = iOctree_IA(iChild, iBlock)

       ! Wipe out the child block
       iOctree_IA(Status_, iBlockChild) = Skipped_
    end do

    ! Make this block used with no children
    iOctree_IA(Status_, iBlock) = Used_

    iBlockChild1 = iOctree_IA(Child1_, iBlock)
    iOctree_IA(Child1_:ChildLast_, iBlock) = NoBlock_

    ! set proc and block info from child1
    iOctree_IA(Proc_,  iBlock) = iOctree_IA(Proc_,  iBlockChild1)
    iOctree_IA(Block_, iBlock) = iOctree_IA(Block_, iBlockChild1)

  end subroutine coarsen_block

  !==========================================================================
  subroutine find_point(XyzIn_D, iBlock)

    ! Find the block that contains a point. The point coordinates should
    ! be given in generalized coordinates normalized to the domain size:
    ! XyzIn_D = (XyzOrig_D - XyzMin_D)/(XyzMax_D-XyzMin_D)

    real, intent(in):: XyzIn_D(MaxDim)
    integer, intent(out):: iBlock

    real :: Xyz_D(MaxDim)
    integer :: iLevel, iChild
    integer :: Ijk_D(MaxDim), iCoord_D(nDim), iBit_D(nDim)
    !----------------------------------------------------------------------
    ! Scale coordinates so that 1 <= Xyz_D <= nRoot_D+1
    Xyz_D = 1.0 + nRoot_D*max(0.0, min(1.0, XyzIn_D))

    ! Get root block index
    Ijk_D = min(int(Xyz_D), nRoot_D)

    ! Root block indexes are ordered
    iBlock = Ijk_D(1) + nRoot_D(1)*((Ijk_D(2)-1) + nRoot_D(2)*(Ijk_D(3)-1))

    if(iOctree_IA(Status_, iBlock) == Used_) RETURN

    ! Get normalized coordinates within root block and scale it up
    ! to the largest resolution
    iCoord_D = (Xyz_D(1:nDim) - Ijk_D(1:nDim))*MaxCoord_I(nLevel)

    ! Go down the tree using bit information
    do iLevel = nLevel-1,0,-1
       iBit_D = ibits(iCoord_D, iLevel, 1)
       iChild = sum(iBit_D*MaxCoord_I(0:nDim-1)) + Child1_
       iBlock = iOctree_IA(iChild, iBlock)

       if(iOctree_IA(Status_, iBlock) == Used_) RETURN
    end do


  end subroutine find_point

  !==========================================================================
  logical function is_point_inside_block(Xyz_D, iBlock)

    real,    intent(in):: Xyz_D(nDim)
    integer, intent(in):: iBlock

    integer :: iLevel
    real    :: XyzStart_D(nDim), XyzEnd_D(nDim)
    !-------------------------------------------------------------------------
    iLevel = iOctree_IA(Level_, iBlock)
    XyzStart_D = (iOctree_IA(Coord1_:CoordLast_,iBlock)-1.0) &
         /MaxCoord_I(iLevel)/nRoot_D(1:nDim)
    XyzEnd_D   = (iOctree_IA(Coord1_:CoordLast_,iBlock)+0.0) &
         /MaxCoord_I(iLevel)/nRoot_D(1:nDim)

    is_point_inside_block= all(Xyz_D >= XyzStart_D) .and. all(Xyz_D < XyzEnd_D)

    if(  any(Xyz_D > XyzEnd_D) .or. any(Xyz_D < XyzStart_D) ) then
       write(*,*)'Error in is_point_inside_block'
       write(*,*)'iBlock, iLevel=',iBlock, iLevel
       write(*,*)'Block start coord=',XyzStart_D
       write(*,*)'Block end   coord=',XyzEnd_D
       write(*,*)'Point coordinates=',Xyz_D
    end if

  end function is_point_inside_block

  !===========================================================================

  subroutine find_neighbors(iBlock)

    integer, intent(in):: iBlock

    integer :: iLevel, i, j, k, Di, Dj, Dk, jBlock
    real :: Scale_D(MaxDim), x, y, z
    !-----------------------------------------------------------------------

    ! We should convert local block into global block index or vice-versa

    iLevel  = iOctree_IA(Level_, iBlock)
    Scale_D = (1.0/MaxCoord_I(iLevel))/nRoot_D
    do k=0,3
       Dk = nint((k - 1.5)/1.5)
       if(nDim < 3)then
          if(k/=1) CYCLE
          z = 0.3
       else
          z = (iOctree_IA(CoordLast_, iBlock) + 0.4*k - 1.1)*Scale_D(3)
          if(z > 1.0 .or. z < 0.0)then
             if(IsPeriodic_D(3))then
                z = modulo(z, 1.0)
             else
                iBlockNei_IIIB(:,:,k,iBlock) = NoBlock_
                DiLevelNei_IIIB(:,:,Dk,iBlock) = NoBlock_
                CYCLE
             end if
          end if
       end if
       do j=0,3
          Dj = nint((j - 1.5)/1.5)
          if(nDim < 2)then
             if(j/=1) CYCLE
             y = 0.3
          else
             y = (iOctree_IA(Coord0_+2, iBlock) + 0.4*j - 1.1)*Scale_D(2)
             if(y > 1.0 .or. y < 0.0)then
                if(IsPeriodic_D(2))then
                   y = modulo(y, 1.0)
                elseif(IsSpherical)then
                   ! Push back theta and go around half way in phi
                   y = max(0.0, min(1.0, y))
                   z = modulo( z+0.5, 1.0)
                else
                   iBlockNei_IIIB(:,j,k,iBlock) = NoBlock_
                   DiLevelNei_IIIB(:,Dj,Dk,iBlock) = NoBlock_
                   CYCLE
                end if
             end if
          end if
          do i=0,3
             ! Exclude inner points
             if(0<i.and.i<3.and.0<j.and.j<3.and.0<k.and.k<3) CYCLE

             Di = nint((i - 1.5)/1.5)

             ! If neighbor is not finer, fill in the i=2 or j=2 or k=2 elements
             if(DiLevelNei_IIIB(Di,Dj,Dk,iBlock) >= 0)then
                if(i==2)then
                   iBlockNei_IIIB(i,j,k,iBlock) = iBlockNei_IIIB(1,j,k,iBlock)
                   CYCLE
                end if
                if(j==2)then
                   iBlockNei_IIIB(i,j,k,iBlock) = iBlockNei_IIIB(i,1,k,iBlock)
                   CYCLE
                end if
                if(k==2)then
                   iBlockNei_IIIB(i,j,k,iBlock) = iBlockNei_IIIB(i,j,1,iBlock)
                   CYCLE
                end if
             end if

             x = (iOctree_IA(Coord1_, iBlock) + 0.4*i - 1.1)*Scale_D(1)
             if(x > 1.0 .or. x < 0.0)then
                if(IsPeriodic_D(1))then
                   x = modulo(x, 1.0)
                elseif(IsCylindrical .and. x < 0.0)then
                   ! Push back radius and go around half way in phi direction
                   x = 0.0
                   z = modulo( z+0.5, 1.0)
                else
                   iBlockNei_IIIB(i,j,k,iBlock) = NoBlock_
                   DiLevelNei_IIIB(Di,Dj,Dk,iBlock) = NoBlock_
                   CYCLE
                end if
             end if

             call find_point( (/x, y, z/), jBlock)
             iBlockNei_IIIB(i,j,k,iBlock) = jBlock
             DiLevelNei_IIIB(Di,Dj,Dk,iBlock) = &
                  iLevel - iOctree_IA(Level_, jBlock)
          end do
       end do
    end do

  end subroutine find_neighbors
  !==========================================================================

  subroutine test_octree

    integer :: iBlock
    real:: XyzTest_D(MaxDim)
    integer:: Int_D(MaxDim)

    character(len=*), parameter :: NameSub = 'test_octree'
    !-----------------------------------------------------------------------

    write(*,*)'Testing init_mod_octree'
    call init_mod_octree(50, 100)
    if(MaxBlock /= 50) &
         write(*,*)'init_mod_octtree faild, MaxBlock=',&
         MaxBlock, ' should be 50'

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

    Int_D = (/1,2,2/)

    if(any( iOctree_IA(Coord1_:CoordLast_,4) /= Int_D(1:nDim) )) &
         write(*,*) 'set_root_block failed, coordinates of block four=',&
         iOctree_IA(Coord1_:CoordLast_,4), ' should be ',Int_D(1:nDim)

    write(*,*)'Testing find_point'
    XyzTest_D = (/0.99,0.99,0.9/)
    call find_point(XyzTest_D, iBlock)
    if(iBlock /= nRoot)write(*,*)'ERROR: Test find point failed, iBlock=',&
         iBlock,' instead of',nRoot

    if(.not.is_point_inside_block(XyzTest_D(1:nDim), iBlock)) &
         write(*,*)'ERROR: Test find point failed'
    
    write(*,*)'Testing refine_block'
    ! Refine the block where the point was found and find it again
    call refine_block(iBlock)

    call find_point(XyzTest_D,iBlock)
    if(.not.is_point_inside_block(XyzTest_D(1:nDim), iBlock)) &
         write(*,*)'ERROR: Test find point failed'

    write(*,*)'Testing find_neighbors'
    call find_neighbors(5)
    write(*,*)'DiLevelNei_IIIB(:,:,:,5)=',DiLevelNei_IIIB(:,:,:,5)
    write(*,*)'iBlockNei_IIIB(:,:,:,5)=',iBlockNei_IIIB(:,:,:,5)

    write(*,*)'Testing coarsen_block'
    ! Coarsen back the last root block and find point again
    call coarsen_block(nRoot)
    call find_point(XyzTest_D,iBlock)
    if(iBlock /= nRoot)write(*,*)'ERROR: coarsen_block faild, iBlock=',&
         iBlock,' instead of',nRoot
    if(.not.is_point_inside_block(XyzTest_D(1:nDim), iBlock)) &
         write(*,*)'ERROR: is_point_inside_block failed'

  end subroutine test_octree

end module ModOctreeNew
