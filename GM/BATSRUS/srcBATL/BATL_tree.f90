module BATL_tree

  use BATL_size, ONLY: MaxBlock, nBlock, MaxDim, nDim, nDimAmr

  implicit none
  save

  private ! except

  public:: init_tree
  public:: clean_tree
  public:: set_tree_root
  public:: refine_tree_node
  public:: coarsen_tree_node
  public:: adapt_tree
  public:: get_tree_position
  public:: find_tree_node
  public:: distribute_tree
  public:: move_tree
  public:: write_tree_file
  public:: read_tree_file
  public:: show_tree
  public:: test_tree

  integer, public, parameter :: nChild = 2**nDimAmr

  ! Global tree information
  integer, public, allocatable :: iTree_IA(:,:)

  integer, public, parameter :: &
       Status_   =  1, &
       Level_    =  2, &
       Proc_     =  3, & ! processor index
       Block_    =  4, & ! block index
       MinLevel_ =  5, & ! minimum level allowed
       MaxLevel_ =  6, & ! maximum level allowed
       Coord0_   =  6, &
       Coord1_   =  7, &
       Coord2_   =  8, &
       Coord3_   =  9, &
       CoordLast_=  9, &
       Parent_   = 10, & ! Parent_ must be 
       Child0_   = 10, & ! equal to Child0_
       Child1_   = Child0_ + 1,      &
       ChildLast_= Child0_ + nChild

  ! Number of items stored in iTree_IA
  integer, parameter :: nInfo = ChildLast_

  character(len=10), parameter:: NameTreeInfo_I(Child0_+8) = (/ &
       'Status   ', &
       'Level    ', &
       'Proc     ', &
       'Block    ', &
       'MinLevel ', &
       'MaxLevel ', &
       'Coord1   ', &
       'Coord2   ', &
       'Coord3   ', &
       'Parent   ', &
       'Child1   ', &
       'Child2   ', &
       'Child3   ', &
       'Child4   ', &
       'Child5   ', &
       'Child6   ', &
       'Child7   ', &
       'Child8   ' /)

  ! New status (refine, coarsen etc) requested for nodes
  integer, public, allocatable :: iStatusNew_A(:)
  integer,         allocatable :: iStatusAll_A(:) ! needed for MPI_allreduce

  ! New processor index of a given node after next load balance
  integer, public, allocatable :: iProcNew_A(:)

  ! Mapping from local block index to global node index
  integer, public, allocatable :: iNode_B(:)

  logical, public, allocatable :: &
       Unused_B(:), Unused_BP(:,:)   ! Unused blocks on local/all processors

  integer, public, allocatable :: &
       DiLevelNei_IIIB(:,:,:,:),  &  ! Level difference relative to neighbors 
       iNodeNei_IIIB(:,:,:,:)        ! Node index of neighboring blocks

  ! Index for unset values (that are otherwise larger)
  integer, public, parameter :: Unset_ = -100

  ! Possible values for the status variable
  integer, public, parameter :: &
       Unused_=-1, Refine_=-2, DontCoarsen_=-3, Coarsen_=-4, &  !(to be) unused
       Used_=1, RefineNew_=2, CoarsenNew_=4                     !(to be) used

  ! Number of total and used nodes (leaves of the node tree)
  integer, public :: nNode = 0, nNodeUsed = 0

  ! Ordering along the Morton-Hilbert space filling curve
  integer, public, allocatable :: iNodeMorton_I(:)

  ! Number of levels below root in level (that has occured at any time)
  integer, public :: nLevel = 0

  ! Local variables -----------------------------------------------
  character(len=*), parameter:: NameMod = "BATL_tree"

  ! Deepest AMR level relative to root nodes (limited by 32 bit integers)
  integer, parameter :: MaxLevel = 30

  ! The maximum integer coordinate for a given level below root nodes
  ! Implied do loop was not understooed by the pgf90 compiler, so list them
  integer, parameter :: MaxCoord_I(0:MaxLevel) = &
       (/ 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, &
       16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, &
       4194304, 8388608, 16777216, 33554432, 67108864, 134217728, &
       268435456, 536870912, 1073741824 /)

  ! The number of root nodes in all dimensions, and altogether
  integer :: nRoot_D(MaxDim) = 0, nRoot = 0

  ! Maximum number of nodes including unused and skipped ones
  integer :: MaxNode = 0

  ! The index along the Morton curve is global so that it can be used by the 
  ! recursive subroutine order_children 
  integer :: iMorton

  ! Needed for compact_tree
  integer, allocatable:: iNodeNew_A(:)

contains

  subroutine init_tree(MaxBlockIn)

    use BATL_mpi, ONLY: nProc

    ! Initialize the tree assuming MaxBlockIn blocks per processor

    integer, intent(in) :: MaxBlockIn ! Max number of blocks per processor
    !----------------------------------------------------------------------
    if(allocated(iTree_IA)) RETURN

    ! Store tree size and maximum number of blocks/processor
    MaxBlock = MaxBlockIn
    MaxNode  = ceiling(nProc*MaxBlock*(1 + 1.0/(nChild - 1)))

    ! Allocate and initialize all elements of tree as unset
    allocate(iTree_IA(nInfo, MaxNode));                 iTree_IA       = Unset_
    allocate(iNodeMorton_I(MaxNode));                   iNodeMorton_I   = Unset_
    allocate(iStatusNew_A(MaxNode));                    iStatusNew_A   = Unset_
    allocate(iStatusAll_A(MaxNode));                    iStatusAll_A   = Unset_
    allocate(iProcNew_A(MaxNode));                      iProcNew_A     = Unset_
    allocate(iNodeNew_A(MaxNode));                      iNodeNew_A     = Unset_
    allocate(iNode_B(MaxBlock));                        iNode_B        = Unset_
    allocate(Unused_B(MaxBlock));                       Unused_B       = .true.
    allocate(Unused_BP(MaxBlock,0:nProc-1));            Unused_BP      = .true.
    allocate(iNodeNei_IIIB(0:3,0:3,0:3,MaxBlock));      iNodeNei_IIIB  = Unset_
    allocate(DiLevelNei_IIIB(-1:1,-1:1,-1:1,MaxBlock)); DiLevelNei_IIIB= Unset_

    ! Initialize minimum and maximum levels of refinement
    iTree_IA(MinLevel_,:) = 0;
    iTree_IA(MaxLevel_,:) = MaxLevel

  end subroutine init_tree

  !==========================================================================
  subroutine clean_tree

    if(.not.allocated(iTree_IA)) RETURN
    deallocate(iTree_IA, iNodeMorton_I, iStatusNew_A, iStatusAll_A, &
         iProcNew_A, iNodeNew_A, &
         iNode_B, Unused_B, Unused_BP, &
         iNodeNei_IIIB, DiLevelNei_IIIB)

    MaxNode = 0

  end subroutine clean_tree
  !==========================================================================

  integer function i_node_new()

    ! Find a skipped element in the iTree_IA array
    
    integer :: iNode
    !-----------------------------------------------------------------------
    do iNode = 1, MaxNode
       if(iTree_IA(Status_, iNode) == Unset_)then
          i_node_new = iNode
          RETURN
       end if
    end do
    ! Could not find any skipped node
    call CON_stop('i_node_new: ran out of nodes')

  end function i_node_new

  !==========================================================================

  subroutine set_tree_root(nRootIn_D)

    integer, optional, intent(in) :: nRootIn_D(nDim)

    integer :: iRoot, jRoot, kRoot, iNode, Ijk_D(MaxDim)
    !-----------------------------------------------------------------------

    ! Set number of root blocks: default or input arguments
    nRoot_D = 1
    if(present(nRootIn_D)) nRoot_D(1:nDim) = nRootIn_D
    nRoot   = product(nRoot_D)

    ! Use the first product(nRoot_D) nodes as root nodes in the tree
    iNode = 0
    do kRoot = 1, nRoot_D(3)
       do jRoot = 1, nRoot_D(2)
          do iRoot = 1, nRoot_D(1)

             Ijk_D = (/ iRoot, jRoot, kRoot /)

             iNode = iNode + 1
             iTree_IA(Status_, iNode)            = Used_
             iTree_IA(Parent_, iNode)            = Unset_
             iTree_IA(Child1_:ChildLast_, iNode) = Unset_
             iTree_IA(Level_ , iNode)            = 0
             iTree_IA(Coord1_:CoordLast_, iNode) = Ijk_D

          end do
       end do
    end do

    nNodeUsed = nRoot
    nNode     = nRoot

  end subroutine set_tree_root

  !==========================================================================
  subroutine refine_tree_node(iNode)

    integer, intent(in) :: iNode

    integer :: iChild, DiChild, iLevelChild, iProc, iBlock
    integer :: iCoord_D(nDimAmr)
    integer :: iDim, iNodeChild

    character(len=*), parameter:: NameSub='refine_tree_node'
    !----------------------------------------------------------------------

    if(iTree_IA(Status_, iNode) == Unused_) &
         call CON_stop(NameSub//' trying to refine and unused block')

    iTree_IA(Status_, iNode) = Refine_

    iLevelChild = iTree_IA(Level_, iNode) + 1
    iProc       = iTree_IA(Proc_,  iNode)
    iBlock      = iTree_IA(Block_, iNode)

    ! Keep track of number of levels
    nLevel = max(nLevel, iLevelChild)
    if(nLevel > MaxLevel) &
         call CON_stop('Error in refine_tree_node: too many levels')

    iCoord_D = 2*iTree_IA(Coord1_:Coord0_+nDimAmr, iNode) - 1

    do iChild = Child1_, ChildLast_

       iNodeChild = i_node_new()

       iTree_IA(iChild, iNode) = iNodeChild

       iTree_IA(Status_,   iNodeChild) = RefineNew_
       iTree_IA(Level_,    iNodeChild) = iLevelChild
       iTree_IA(MinLevel_, iNodeChild) = iTree_IA(MinLevel_,iNode)
       iTree_IA(MaxLevel_, iNodeChild) = iTree_IA(MaxLevel_,iNode)
       iTree_IA(Parent_,   iNodeChild) = iNode
       iTree_IA(Child1_:ChildLast_, iNodeChild) = Unset_

       ! Data will come from the parent's proc/block
       iTree_IA(Proc_,     iNodeChild) = iProc
       iTree_IA(Block_,    iNodeChild) = iBlock

       ! Calculate the coordinates of the child node
       DiChild = iChild - Child1_
       do iDim = 1, nDimAmr
          iTree_IA(Coord0_+iDim, iNodeChild) = &
               iCoord_D(iDim) + ibits(DiChild, iDim-1, 1)
       end do

       ! The non-AMR coordinates remain the same as for the parent node
       do iDim = nDimAmr+1, MaxDim
          iTree_IA(Coord0_+iDim, iNodeChild) = iTree_IA(Coord0_+iDim, iNode)
       end do

    end do

    ! Keep track of used nodes in the future tree 
    nNodeUsed = nNodeUsed + nChild - 1

  end subroutine refine_tree_node

  !==========================================================================
  subroutine coarsen_tree_node(iNode)

    integer, intent(in) :: iNode

    integer :: iChild, iNodeChild
    !-----------------------------------------------------------------------

    do iChild = Child1_, ChildLast_
       iNodeChild = iTree_IA(iChild, iNode)

       ! Set the status of the child node
       iTree_IA(Status_, iNodeChild) = Coarsen_
    end do

    ! Make this node used with no children
    iTree_IA(Status_, iNode) = CoarsenNew_

    ! Keep track of used nodes in the future tree
    nNodeUsed = nNodeUsed - nChild + 1

  end subroutine coarsen_tree_node

  !===========================================================================

  subroutine adapt_tree

    use BATL_size, ONLY: iRatio, jRatio, kRatio
    use BATL_mpi, ONLY: iComm, nProc, iProc, barrier_mpi
    use ModMpi, ONLY: MPI_allreduce, MPI_INTEGER, MPI_MAX

    ! All processors can request some status changes in iStatusNew_A.
    ! Here we collect requests, check for proper nesting, 
    ! limitations on level, number of blocks, etc, 
    ! modify iStatusNew_A and set iTree_IA.

    integer:: nNodeUsedNow, iMorton, iBlock, iStatus, iError
    integer:: iNode, iNodeParent, iNodeChild_I(nChild)
    integer:: jNode, jNodeParent, jNodeChild_I(nChild)
    integer:: iLevel, iLevelNew, iLevelMax, iLevelMin
    integer:: jLevel, jLevelNew
    integer:: iSide, iSideMin, iSideMax
    integer:: jSide, jSideMin, jSideMax
    integer:: kSide, kSideMin, kSideMax

    logical, parameter :: DoTest = .false., DoTestNei = .false.
    !------------------------------------------------------------------------

    ! Collect the local status requests into a global request
    if(nProc > 1)then
       call MPI_allreduce(iStatusNew_A, iStatusAll_A, nNode, MPI_INTEGER, &
            MPI_MAX, iComm, iError)
       iStatusNew_A = iStatusAll_A
    end if

    ! Check max and min levels and coarsening of all siblings
    iLevelMin = nLevel
    iLevelMax = 1
    do iMorton = 1, nNodeUsed
       iNode   = iNodeMorton_I(iMorton)

       iLevel    = iTree_IA(Level_,iNode)
       iLevelMin = min(iLevelMin, iLevel)

       ! Check MaxLevel_ of node to see if it can be refined 
       if(iStatusNew_A(iNode) == Refine_)then
          if(iLevel >= iTree_IA(MaxLevel_,iNode))then
             iStatusNew_A(iNode) = Unset_
             CYCLE
          end if
          iLevelMax = max(iLevelMax, iLevel)
       end if

       ! Only nodes to be coarsened need further checking
       if(iStatusNew_A(iNode) /= Coarsen_) CYCLE

       ! Check MinLevel_ of node to see if it can be coarsened
       if(iLevel <= iTree_IA(MinLevel_,iNode))then
          iStatusNew_A(iNode) = Unset_
          CYCLE
       end if
       iLevelMax = max(iLevelMax, iLevel)

       ! Check if all siblings want to be coarsened
       iNodeParent = iTree_IA(Parent_,iNode)
       iNodeChild_I = iTree_IA(Child1_:ChildLast_,iNodeParent)
       if(.not.all(iStatusNew_A(iNodeChild_I) == Coarsen_))then
          ! Cancel coarsening requests for all siblings
          where(iStatusNew_A(iNodeChild_I) == Coarsen_) &
               iStatusNew_A(iNodeChild_I) = Unset_
       end if

    end do

    ! Check proper nesting. Go down level by level. No need to check base level
    ! Changes in the requests will be applied to all siblings immediately
    do iLevel = iLevelMax, max(iLevelMin, 1), -1

       ! Parallel processing of nodes (blocks)
       BLOCKLOOP: do iBlock = 1, nBlock

          if(Unused_B(iBlock)) CYCLE BLOCKLOOP
          iNode = iNode_B(iBlock)
          if(iTree_IA(Level_,iNode) /= iLevel) CYCLE BLOCKLOOP

          ! Calculate requested level
          if(iStatusNew_A(iNode) == Refine_)then
             iLevelNew = iLevel + 1
          elseif(iStatusNew_A(iNode) == Coarsen_)then
             iLevelNew = iLevel - 1
          else
             CYCLE
          end if

          ! Check neighbors around this corner of the parent block

          ! Loop from 0 to 1 or from 2 to 3 in the side index 
          ! depending on which side this node
          ! is relative to its parent. If there is no refinement in some 
          ! direction, then loop from 0 to 3 (and skip 2).

          if(iRatio==1)then
             iSideMin = 0; iSideMax = 3
          else
             iSideMin = 2*modulo(iTree_IA(Coord1_,iNode)-1, 2)
             iSideMax = iSideMin + 1
          end if
          if(nDim < 2)then
             ! 1D, no neighbors in j direction
             jSideMin = 1; jSideMax = 1
          elseif(jRatio == 1)then
             ! 2D or 3D but no AMR, check neighbors in both directions
             jSideMin = 0; jSideMax = 3
          else
             ! 2D or 3D and AMR: check only the directions 
             ! corresponding to the corner occupied by this child
             jSideMin = 2*modulo(iTree_IA(Coord2_,iNode)-1, 2)
             jSideMax = jSideMin + 1
          end if
          if(nDim < 3)then
             ! 1 or 2D
             kSideMin = 1; kSideMax = 1
          elseif(kRatio == 1)then
             ! 3D but but no AMR, check neighbors in both directions
             kSideMin = 0; kSideMax = 3
          else
             ! 3D and AMR, check only the directions
             ! corresponding to the corner occupied by this child
             kSideMin = 2*modulo(iTree_IA(Coord3_,iNode)-1, 2)
             kSideMax = kSideMin + 1
          end if

          ! Loop through the at most seven neighbors
          do kSide = kSideMin, kSideMax
             if(kRatio == 1 .and. kSide == 2) CYCLE
             do jSide = jSideMin, jSideMax
                if(jRatio == 1 .and. jSide == 2) CYCLE
                do iSide = iSideMin, iSideMax
                   if(iRatio == 1 .and. iSide == 2) CYCLE

                   jNode = iNodeNei_IIIB(iSide,jSide,kSide,iBlock)

                   ! Don't check the node itself
                   if(iNode == jNode) CYCLE

                   ! Don't check if neighbor is outside the domain
                   if(jNode == Unset_) CYCLE

                   ! Get the current and requested level for the neighbor node
                   jLevel = iTree_IA(Level_,jNode)

                   jLevelNew = jLevel
                   if(iStatusNew_A(jNode) == Refine_)  jLevelNew = jLevel + 1
                   if(iStatusNew_A(jNode) == Coarsen_) jLevelNew = jLevel - 1

                   ! Fix levels if difference is too much
                   if(iLevelNew >= jLevelNew + 2)then

                      if(jLevel > 0)then
                         jNodeParent = iTree_IA(Parent_,jNode) 
                         jNodeChild_I= iTree_IA(Child1_:ChildLast_,jNodeParent)

                         ! Neighbor and its siblings cannot be coarsened
                         iStatusNew_A(jNodeChild_I) = &
                              max(iStatusNew_A(jNodeChild_I), DontCoarsen_)
                      endif

                      ! If neighbow was coarser it has to be refined
                      if(jLevel < iLevel) iStatusNew_A(jNode) = Refine_

                   elseif(iLevelNew <= jLevelNew - 2)then
                      ! Cannot coarsen this node
                      iNodeParent = iTree_IA(Parent_,iNode)
                      iStatusNew_A(iTree_IA(Child1_:ChildLast_,iNodeParent)) &
                           = DontCoarsen_
                      
                      CYCLE BLOCKLOOP
                   end if

                end do ! iSide
             end do ! jSide
          end do ! kSide
       end do BLOCKLOOP

       ! Collect the local status requests into a global request
       if(nProc > 1)then
          call MPI_allreduce(iStatusNew_A, iStatusAll_A, nNode, MPI_INTEGER, &
               MPI_MAX, iComm, iError)
          iStatusNew_A = iStatusAll_A
       end if

    end do ! levels

    nNodeUsedNow = nNodeUsed

    ! Coarsen first to reduce number of nodes and used blocks
    do iMorton = 1, nNodeUsedNow
       iNode   = iNodeMorton_I(iMorton)
       iStatus = iStatusNew_A(iNode)

       if(iStatus /= Coarsen_) CYCLE

       iNodeParent = iTree_IA(Parent_,iNode)

       ! Coarsen the parent node based on the request stored in the first child
       if(iTree_IA(Child1_,iNodeParent) /= iNode) CYCLE

       call coarsen_tree_node(iNodeParent)
    end do

    ! Refine next
    do iMorton = 1, nNodeUsedNow
       iNode   = iNodeMorton_I(iMorton)
       iStatus = iStatusNew_A(iNode)

       if(iStatus /= Refine_) CYCLE

       ! Refine tree node
       call refine_tree_node(iNode)

       if(nNodeUsed > MaxBlock*nProc) EXIT

    end do

    iStatusNew_A = Unset_

  end subroutine adapt_tree
  !==========================================================================
  subroutine get_tree_position(iNode, PositionMin_D, PositionMax_D)

    integer, intent(in) :: iNode
    real,    intent(out):: PositionMin_D(MaxDim), PositionMax_D(MaxDim)

    ! Calculate normalized position of the edges of node inode.
    ! Zero is at the minimum boundary of the grid, one is at the max boundary

    integer :: iLevel
    integer :: MaxIndex_D(MaxDim)
    !------------------------------------------------------------------------
    iLevel = iTree_IA(Level_, iNode)

    MaxIndex_D = nRoot_D
    MaxIndex_D(1:nDimAmr) = MaxCoord_I(iLevel)*MaxIndex_D(1:nDimAmr)

    ! Convert to real by adding -1.0 or 0.0 for the two edges, respectively
    PositionMin_D = (iTree_IA(Coord1_:CoordLast_,iNode) - 1.0)/MaxIndex_D
    PositionMax_D = (iTree_IA(Coord1_:CoordLast_,iNode) + 0.0)/MaxIndex_D

  end subroutine get_tree_position

  !==========================================================================
  subroutine find_tree_node(CoordIn_D, iNode)

    ! Find the node that contains a point. The point coordinates should
    ! be given in generalized coordinates normalized to the domain size:
    ! CoordIn_D = (CoordOrig_D - CoordMin_D)/(CoordMax_D-CoordMin_D)

    real, intent(in):: CoordIn_D(MaxDim)
    integer, intent(out):: iNode

    real :: Coord_D(MaxDim)
    integer :: iLevel, iChild
    integer :: Ijk_D(MaxDim), iCoord_D(nDimAmr), iBit_D(nDimAmr)
    !----------------------------------------------------------------------
    ! Scale coordinates so that 1 <= Coord_D <= nRoot_D+1
    Coord_D = 1.0 + nRoot_D*max(0.0, min(1.0, CoordIn_D))

    ! Get root node index
    Ijk_D = min(int(Coord_D), nRoot_D)

    ! Root node indexes are ordered
    iNode = Ijk_D(1) + nRoot_D(1)*((Ijk_D(2)-1) + nRoot_D(2)*(Ijk_D(3)-1))

    if(iTree_IA(Status_, iNode) == Used_) RETURN

    ! Get normalized coordinates within root node and scale it up
    ! to the largest resolution
    iCoord_D = (Coord_D(1:nDimAmr) - Ijk_D(1:nDimAmr))*MaxCoord_I(nLevel)

    ! Go down the tree using bit information
    do iLevel = nLevel-1,0,-1
       iBit_D = ibits(iCoord_D, iLevel, 1)
       iChild = sum(iBit_D*MaxCoord_I(0:nDimAmr-1)) + Child1_
       iNode = iTree_IA(iChild, iNode)

       if(iTree_IA(Status_, iNode) == Used_) RETURN
    end do

  end subroutine find_tree_node

  !==========================================================================
  logical function is_point_inside_node(Position_D, iNode)

    ! Check if position is inside node or not

    real,    intent(in):: Position_D(MaxDim)
    integer, intent(in):: iNode

    real    :: PositionMin_D(MaxDim), PositionMax_D(MaxDim)
    !-------------------------------------------------------------------------
    call get_tree_position(iNode, PositionMin_D, PositionMax_D)

    ! Include min edge but exclude max edge for sake of uniqueness
    is_point_inside_node = &
         all(Position_D >= PositionMin_D) .and. &
         all(Position_D <  PositionMax_D)

  end function is_point_inside_node

  !===========================================================================

  subroutine find_neighbor(iBlock)

    use BATL_geometry, ONLY: IsCylindrical, IsSpherical, IsPeriodic_D

    integer, intent(in):: iBlock

    integer :: iNode, iLevel, i, j, k, Di, Dj, Dk, jNode
    real :: Scale_D(MaxDim), x, y, z

    logical, parameter :: DoTestMe = .false.
    !-----------------------------------------------------------------------
    iNode = iNode_B(iBlock)
    if(DoTestMe)write(*,*)'Starting find neighbors for node ',iNode

    ! Get AMR level of the node
    iLevel = iTree_IA(Level_,iNode)

    ! Calculate scaling factor from integer index to 0<x,y,z<1 real coordinates
    Scale_D = 1.0/nRoot_D
    Scale_D(1:nDimAmr) = Scale_D(1:nDimAmr)/MaxCoord_I(iLevel)

    if(DoTestMe)then
       write(*,*)'iNode, iLevel, Scale_D=', iNode, iLevel, Scale_D
       write(*,*)'scaled coordinates=', &
            iTree_IA(Coord1_:CoordLast_, iNode)*Scale_D
    end if

    ! Fill in self-referring info
    iNodeNei_IIIB(1:2,1:2,1:2,iBlock) = iNode
    DiLevelNei_IIIB(0,0,0,iBlock)     = 0

    ! Loop through neighbors
    do k=0,3
       Dk = nint((k - 1.5)/1.5)
       if(nDim < 3)then
          if(k/=1) CYCLE
          z = 0.3
       else
          z = (iTree_IA(Coord3_, iNode) + 0.4*k - 1.1)*Scale_D(3)
          if(z > 1.0 .or. z < 0.0)then
             if(IsPeriodic_D(3))then
                z = modulo(z, 1.0)
             else
                iNodeNei_IIIB(:,:,k,iBlock) = Unset_
                DiLevelNei_IIIB(:,:,Dk,iBlock) = Unset_
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
             y = (iTree_IA(Coord2_, iNode) + 0.4*j - 1.1)*Scale_D(2)
             if(y > 1.0 .or. y < 0.0)then
                if(IsPeriodic_D(2))then
                   y = modulo(y, 1.0)
                elseif(IsSpherical)then
                   ! Push back theta and go around half way in phi
                   y = max(0.0, min(1.0, y))
                   z = modulo( z+0.5, 1.0)
                else
                   iNodeNei_IIIB(:,j,k,iBlock) = Unset_
                   DiLevelNei_IIIB(:,Dj,Dk,iBlock) = Unset_
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
                   iNodeNei_IIIB(i,j,k,iBlock) = iNodeNei_IIIB(1,j,k,iBlock)
                   CYCLE
                end if
                if(j==2)then
                   iNodeNei_IIIB(i,j,k,iBlock) = iNodeNei_IIIB(i,1,k,iBlock)
                   CYCLE
                end if
                if(k==2)then
                   iNodeNei_IIIB(i,j,k,iBlock) = iNodeNei_IIIB(i,j,1,iBlock)
                   CYCLE
                end if
             end if

             x = (iTree_IA(Coord1_, iNode) + 0.4*i - 1.1)*Scale_D(1)
             if(x > 1.0 .or. x < 0.0)then
                if(IsPeriodic_D(1))then
                   x = modulo(x, 1.0)
                elseif(IsCylindrical .and. x < 0.0)then
                   ! Push back radius and go around half way in phi direction
                   x = 0.0
                   z = modulo( z+0.5, 1.0)
                else
                   iNodeNei_IIIB(i,j,k,iBlock) = Unset_
                   DiLevelNei_IIIB(Di,Dj,Dk,iBlock) = Unset_
                   CYCLE
                end if
             end if

             call find_tree_node( (/x, y, z/), jNode)

             iNodeNei_IIIB(i,j,k,iBlock) = jNode
             DiLevelNei_IIIB(Di,Dj,Dk,iBlock) = &
                  iLevel - iTree_IA(Level_, jNode)

             if(DoTestMe)write(*,'(a,3i2,3f6.3,i4)') &
                  'i,j,k,x,y,z,jNode=',i,j,k,x,y,z,jNode

          end do
       end do
    end do
    
  end subroutine find_neighbor

  !==========================================================================

  subroutine compact_tree

    ! Eliminate holes from the tree

    ! Amount of shift for each node
    integer :: iNode, iNodeSkipped, iNodeOld, i, iBlock
    !-------------------------------------------------------------------------
    ! Set impossible initial value
    iNodeSkipped = MaxNode + 1

    do iNode = 1, MaxNode !!! nNode ???

       if(iTree_IA(Status_, iNode) == Unset_)then
          ! Store the first skipped position
          iNodeSkipped = min(iNodeSkipped, iNode)
       elseif(iNodeSkipped < iNode)then
          ! Move node to the first skipped position
          iTree_IA(:,iNodeSkipped) = iTree_IA(:,iNode)
          iTree_IA(Status_, iNode) = Unset_
          ! Store new node index
          iNodeNew_A(iNode) = iNodeSkipped
          ! Advance iNodeSkipped
          iNodeSkipped = iNodeSkipped + 1
       else
          ! The node did not move
          iNodeNew_A(iNode) = iNode
       endif
    end do

    ! Apply shifts
    do iNode = 1, MaxNode
       if(iTree_IA(Status_, iNode) == Unset_) EXIT
       do i = Parent_, ChildLast_
          iNodeOld = iTree_IA(i, iNode)
          if(iNodeOld /= Unset_) &
               iTree_IA(i, iNode) = iNodeNew_A(iNodeOld)
       end do
    end do

    ! Set number of nodes in the tree (note the EXIT above)
    nNode = iNode - 1

    ! Fix the node indexes along the Morton curve
    do iMorton = 1, nNodeUsed
       iNodeOld = iNodeMorton_I(iMorton)
       iNodeMorton_I(iMorton) = iNodeNew_A(iNodeOld)
    end do
    
    ! Fix iNode_B indexes
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       iNodeOld = iNode_B(iBlock)
       iNode_B(iBlock) = iNodeNew_A(iNodeOld)
    end do

    ! Reset iNodeNew_A
    iNodeNew_A = Unset_

  end subroutine compact_tree

  !==========================================================================

  subroutine write_tree_file(NameFile)

    use ModIoUnit, ONLY: UnitTmp_
    use BATL_mpi, ONLY: iProc, barrier_mpi

    character(len=*), intent(in):: NameFile

    ! Write tree information into a file
    !-------------------------------------------------------------------------
    call compact_tree

    if(iProc == 0)then
       open(UnitTmp_, file=NameFile, status='replace', form='unformatted')

       write(UnitTmp_) nNode, nInfo
       write(UnitTmp_) nDim, nDimAmr
       write(UnitTmp_) nRoot_D(1:nDim)
       write(UnitTmp_) iTree_IA(:,1:nNode)

       close(UnitTmp_)
    end if
    call barrier_mpi

  end subroutine write_tree_file
  
  !==========================================================================

  subroutine read_tree_file(NameFile)

    use ModIoUnit, ONLY: UnitTmp_

    character(len=*), intent(in):: NameFile

    ! Read tree information from a file

    integer :: nInfoIn, nNodeIn, nDimIn, nDimAmrIn, nRootIn_D(nDim)
    character(len=*), parameter :: NameSub = 'read_tree_file'
    !----------------------------------------------------------------------

    open(UnitTmp_, file=NameFile, status='old', form='unformatted')

    read(UnitTmp_) nNodeIn, nInfoIn
    if(nNodeIn > MaxNode)then
       write(*,*) NameSub,' nNodeIn, MaxNode=',nNodeIn, MaxNode 
       call CON_stop(NameSub//' too many nodes in tree file!')
    end if
    read(UnitTmp_) nDimIn, nDimAmrIn
    if(nDimIn /= nDim)then
       write(*,*) NameSub,' nDimIn, nDim=',nDimIn, nDim
       call CON_stop(NameSub//' nDim is different in tree file!')
    end if
    if(nDimAmrIn /= nDimAmr)then
       write(*,*) NameSub,' nDimAmrIn, nDimAmr=',nDimAmrIn, nDimAmr
       call CON_stop(NameSub//' nDimAmr is different in tree file!')
    end if
    read(UnitTmp_) nRootIn_D

    call set_tree_root(nRootIn_D)

    read(UnitTmp_) iTree_IA(:,1:nNodeIn)
    close(UnitTmp_)

    ! It is probably ordered already...
    call order_tree

  end subroutine read_tree_file
  
  !==========================================================================
  subroutine distribute_tree(DoMove, iTypeNode_A)

    ! Order tree with the space filling curve then
    ! - if DoMove=T, assign tree nodes to processors and blocks immediately
    ! - if DoMove=F, set iProcNew_A only with future processor index
    ! - if iTypeNode_A is present, it contains block types 1, 2, .., nType
    !   each type is balanced separately. The total is also balanced.

    use BATL_mpi, ONLY: nProc
    
    ! Are nodes moved immediately or just assigned new processor/node
    logical, intent(in):: DoMove

    ! Optional block type. Each type is balanced separately
    integer, intent(in), optional:: iTypeNode_A(MaxNode)

    integer :: nNodePerProc, iMorton, iNode, iBlockTo, iProcTo

    integer :: iType, nType, iProcStart, iProcStop, iProcExtraBlock
    integer, allocatable :: iNodeType_I(:), nNodeType_I(:), iProcType_I(:), &
         iBlock_P(:), nBlockType_PI(:,:)
    !------------------------------------------------------------------------
    if(DoMove)Unused_BP = .true.

    ! Initialize processor and block indexes
    iProcNew_A = Unset_

    ! Create iNodeMorton_I
    call order_tree

    ! Check if there are multiple node types that need separate balancing
    if(present(iTypeNode_A))then

       ! Find number of types and allocate arrays
       nType = maxval(iTypeNode_A, MASK=iTree_IA(Status_,:)>=Used_)

       ! Allocate load balance tables: 
       ! number of blocks per type and per processor and type
       allocate(&
            iNodeType_I(nType), nNodeType_I(nType), &
            iProcType_I(nType), iBlock_P(0:nProc-1), &
            nBlockType_PI(0:nProc-1,nType))

       ! Initialize number of type, counter, processor index and block index
       ! for various node types
       nNodeType_I = 0
       iNodeType_I = 0
       iProcType_I = 0
       iBlock_P    = 0

       ! Count number of nodes for each type. 
       do iNode = 1, nNode
          if(iTree_IA(Status_,iNode)<=0) CYCLE
          iType = iTypeNode_A(iNode)
          if(iType > 0) nNodeType_I(iType) = nNodeType_I(iType) + 1
       end do

       ! write(*,*)'nType, nNodeType_I=', nType, nNodeType_I

       ! Construct load balance table for various types
       do iType = 1, nType
          ! minimum number of blocks of type iType for each processor
          nBlockType_PI(:,iType) = nNodeType_I(iType)/nProc

          ! The processors with extra blocks are filled in 
          ! from nProc-1 backwards
          iProcStart = nProc - modulo(sum(nNodeType_I(1:iType)),nProc)
          iProcStop  = iProcStart + modulo(nNodeType_I(iType),nProc) - 1
          do iProcExtraBlock = iProcStart, iProcStop
             iProcTo = modulo(iProcExtraBlock,nProc)
             nBlockType_PI(iProcTo,iType) = nBlockType_PI(iProcTo,iType) + 1
          end do

          ! convert nBlockType_PI to cummulative load table for easier use
          do iProcTo = 1, nProc-1
             nBlockType_PI(iProcTo,iType) = nBlockType_PI(iProcTo,iType) &
                  + nBlockType_PI(iProcTo-1,iType)
          end do

          ! write(*,*)'iType, nBlockType_PI=', iType, nBlockType_PI
       end do

       ! Distribute the nodes over the processors
       do iMorton = 1, nNodeUsed

          ! Get the node index and type
          iNode = iNodeMorton_I(iMorton)
          iType = iTypeNode_A(iNode)

          ! Increase the index for this node type
          iNodeType_I(iType) = iNodeType_I(iType) + 1

          ! Select target processor. 
          ! Use iProcType_I to remember last proc. used for the given type
          do iProcTo = iProcType_I(iType), nProc-1
             if(iNodeType_I(iType) <= nBlockType_PI(iProcTo,iType))then
                iProcType_I(iType) = iProcTo
                EXIT
             end if
          end do

          ! Assign future processor index for node
          iProcNew_A(iNode) = iProcTo

          if(DoMove)then
             ! Assign block index right away
             iBlock_P(iProcTo) = iBlock_P(iProcTo) + 1
             iBlockTo = iBlock_P(iProcTo)
             iTree_IA(Block_,iNode) = iBlockTo
             Unused_BP(iBlockTo,iProcTo) = .false.
          end if

       end do

       ! write(*,*)'iProcNew_A=', iProcNew_A(1:nNode)

       deallocate(iNodeType_I, nNodeType_I, iProcType_I, iBlock_P, &
            nBlockType_PI)

    else

       ! An upper estimate on the number of nodes per processor
       nNodePerProc = (nNodeUsed-1)/nProc + 1

       do iMorton = 1, nNodeUsed
          iNode = iNodeMorton_I(iMorton)

          iProcTo  = (iMorton-1)/nNodePerProc

          ! Assign future processor index for node
          iProcNew_A(iNode) = iProcTo

          if(DoMove)then
             ! Assign block index right away
             iBlockTo = modulo(iMorton-1, nNodePerProc) + 1
             iTree_IA(Block_,iNode) = iBlockTo
             Unused_BP(iBlockTo,iProcTo) = .false.
          end if

       end do
    end if

    if(DoMove) call move_tree

  end subroutine distribute_tree

  !==========================================================================

  subroutine move_tree

    ! Finish the load balancing (with or without data movement)
    ! Set status for newly used and unused/unset nodes.
    ! Then compact the tree and find all the neighbors

    use BATL_mpi, ONLY: iProc

    integer:: iMorton, iNode, iNodeChild, iNodeParent, iChild, iBlock
    !-----------------------------------------------------------------------
    ! Update local Unused_B array
    Unused_B = Unused_BP(:,iProc)

    ! Update nBlock too as we loop through the used blocks
    nBlock = 0
    do iMorton = 1, nNodeUsed
       iNode = iNodeMorton_I(iMorton)

       ! Move the node to new processor/node
       iTree_IA(Proc_,iNode) = iProcNew_A(iNode)

       if(iTree_IA(Status_,iNode) == CoarsenNew_) then
          ! Remove the children of newly coarsened blocks from the tree

          do iChild = Child1_, ChildLast_
             iNodeChild = iTree_IA(iChild, iNode)
             iTree_IA(:,iNodeChild) = Unset_
          end do
          iTree_IA(Child1_:ChildLast_, iNode) = Unset_

       elseif(iTree_IA(Status_,iNode) == RefineNew_)then
          ! Make the parent of newly refined blocks unused

          iNodeParent = iTree_IA(Parent_, iNode)
          iTree_IA(Proc_:Block_,iNodeParent) = Unset_
          iTree_IA(Status_,iNodeParent)      = Unused_
       end if

       ! Now newly formed blocks are simply used
       iTree_IA(Status_,iNode) = Used_

       ! Set local information for this processor
       if(iProc == iTree_IA(Proc_,iNode))then
          iBlock = iTree_IA(Block_,iNode)
          iNode_B(iBlock)  = iNode
          if(.not.Unused_B(iBlock)) nBlock = max(nBlock, iBlock)
       end if

    end do

    ! Now that we removed children of coarsened blocks, compact the tree
    call compact_tree

    ! Set neighbor info
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       call find_neighbor(iBlock)
    end do

  end subroutine move_tree
  !==========================================================================
  subroutine order_tree

    ! Set iNodeMorton_I indirect index array according to 
    ! 1. root node order
    ! 2. Morton ordering for each root node

    integer :: iNode, iRoot, jRoot, kRoot
    !-----------------------------------------------------------------------
    nNode = nRoot
    iNode = 0
    iMorton = 0
    iNodeMorton_I = Unset_
    do kRoot = 1, nRoot_D(3)
       do jRoot = 1, nRoot_D(2)
          do iRoot = 1, nRoot_D(1)
             ! Root nodes are the first ones
             iNode = iNode + 1

             ! All root nodes are handled as if they were first child
             call order_children(iNode)
          end do
       end do
    end do

    nNodeUsed = iMorton

  end subroutine order_tree
  !==========================================================================
  recursive subroutine order_children(iNode)

    ! Recursively apply Morton ordering for nodes below a root block.
    ! Store result into iNodeMorton_I using the global iMorton index.

    integer, intent(in) :: iNode
    integer :: iChild
    !-----------------------------------------------------------------------
    nNode = max(nNode, iNode)

    if(iTree_IA(Status_, iNode) >= Used_)then
       iMorton = iMorton + 1
       iNodeMorton_I(iMorton) = iNode
    else
       do iChild = Child1_, ChildLast_
          call order_children(iTree_IA(iChild, iNode))
       end do
    end if

  end subroutine order_children
  !==========================================================================

  subroutine show_tree(String, DoShowNei)

    use BATL_geometry, ONLY: IsPeriodic_D

    character(len=*), intent(in):: String
    logical, optional,intent(in):: DoShowNei

    ! Show complete tree information. Also write out string as an identifier.

    character(len=10) :: Name
    character(len=200):: Header
    integer:: iInfo, iNode, iBlock
    !-----------------------------------------------------------------------
    Header = 'iNode'
    do iInfo = 1, nInfo
       Name = NameTreeInfo_I(iInfo)
       Header(7*iInfo+1:7*(iInfo+1)-1) = Name(1:6)
    end do

    write(*,*) String
    write(*,*) trim(Header)
    do iNode = 1, MaxNode
       if(iTree_IA(Status_, iNode) == Unset_) CYCLE
       write(*,'(100i7)') iNode, iTree_IA(:, iNode)
    end do

    if(.not.present(DoShowNei)) RETURN
    if(.not.DoShowNei) RETURN

    write(*,*)'nNode, nNodeUsed, nBlock=',nNode, nNodeUsed, nBlock
    write(*,*)'iNodeMorton_I =', iNodeMorton_I(1:nNodeUsed)
    write(*,*)'IsPeriodic_D =', IsPeriodic_D

    iNode = iNodeMorton_I(1)
    iBlock = iTree_IA(Block_,iNode)
    write(*,*)'DiLevelNei_IIIB(:,0,0,First)=', DiLevelNei_IIIB(:,0,0,iBlock)
    write(*,*)'DiLevelNei_IIIB(0,:,0,First)=', DiLevelNei_IIIB(0,:,0,iBlock)
    write(*,*)'DiLevelNei_IIIB(0,0,:,First)=', DiLevelNei_IIIB(0,0,:,iBlock)
    write(*,*)'iNodeNei_IIIB(:,1,1,  First)=',   iNodeNei_IIIB(:,1,1,iBlock)
    write(*,*)'iNodeNei_IIIB(1,:,1,  First)=',   iNodeNei_IIIB(1,:,1,iBlock)
    write(*,*)'iNodeNei_IIIB(1,1,:,  First)=',   iNodeNei_IIIB(1,1,:,iBlock)

!    iNode = iNodeMorton_I(nNodeUsed)
!    write(*,*)'DiLevelNei_IIIA(:,0,0, Last)=', DiLevelNei_IIIA(:,0,0,iNode)
!    write(*,*)'DiLevelNei_IIIA(0,:,0, Last)=', DiLevelNei_IIIA(0,:,0,iNode)
!    write(*,*)'DiLevelNei_IIIA(0,0,:, Last)=', DiLevelNei_IIIA(0,0,:,iNode)
!    write(*,*)'iNodeNei_IIIA(:,1,1,   Last)=',   iNodeNei_IIIA(:,1,1,iNode)
!    write(*,*)'iNodeNei_IIIA(1,:,1,   Last)=',   iNodeNei_IIIA(1,:,1,iNode)
!    write(*,*)'iNodeNei_IIIA(1,1,:,   Last)=',   iNodeNei_IIIA(1,1,:,iNode)

  end subroutine show_tree

  !==========================================================================

  subroutine test_tree

    use BATL_mpi, ONLY: iProc, nProc
    use BATL_geometry, ONLY: init_geometry, IsPeriodic_D

    integer, parameter:: MaxBlockTest            = 50
    integer, parameter:: nRootTest_D(MaxDim)     = (/3,2,1/)
    logical, parameter:: IsPeriodicTest_D(MaxDim)= (/.true., .true., .false./)
    real,    parameter:: CoordTest_D(MaxDim)     = 0.99

    integer :: iNode, Int_D(MaxDim)
    integer, allocatable:: iTypeNode_I(:)

    logical :: DoTestMe
 
    character(len=*), parameter :: NameSub = 'test_tree'
    !-----------------------------------------------------------------------

    DoTestMe = iProc == 0

    if(DoTestMe)write(*,*)'Testing init_tree'
    call init_tree(MaxBlockTest)
    if(MaxBlock /= MaxBlockTest) &
         write(*,*)'init_tree failed, MaxBlock=',&
         MaxBlock, ' should be ',MaxBlockTest
    if(MaxNode /= ceiling(MaxBlockTest*nProc*(1 + 1.0/(2**nDimAmr-1)))) &
         write(*,*)'init_tree failed, MaxNode=', MaxNode, &
         ' should be', ceiling(50*nProc*(1 + 1.0/(2**nDimAmr-1)))

    if(DoTestMe)write(*,*)'Testing init_geometry'
    call init_geometry('cartesian', IsPeriodicTest_D(1:nDim))
    if(any(IsPeriodic_D(1:nDim) .neqv. IsPeriodicTest_D(1:nDim))) &
         write(*,*)'init_geometry failed, IsPeriodic_D=',&
         IsPeriodic_D(1:nDim), ' should be ', IsPeriodicTest_D(1:nDim)

    if(DoTestMe)write(*,*)'Testing i_node_new()'
    iNode = i_node_new()
    if(iNode /= 1) &
         write(*,*)'i_node_new() failed, iNode=',iNode,' should be 1'

    if(DoTestMe)write(*,*)'Testing set_tree_root'
    call set_tree_root( nRootTest_D(1:nDim))

    if(DoTestMe)call show_tree('after set_tree_root')

    if(any( nRoot_D(1:nDim) /= nRootTest_D(1:nDim) )) &
         write(*,*) 'set_tree_root failed, nRoot_D=',nRoot_D(1:nDim),&
         ' should be ',nRootTest_D(1:nDim)

    Int_D = (/3,1,1/)

    if(any( iTree_IA(Coord1_:Coord0_+nDim,3) /= Int_D(1:nDim) )) &
         write(*,*) 'set_tree_root failed, coordinates of node four=',&
         iTree_IA(Coord1_:Coord0_+nDim,3), ' should be ',Int_D(1:nDim)

    if(DoTestMe)write(*,*)'Testing find_tree_node'
    call find_tree_node(CoordTest_D, iNode)
    if(iNode /= nRoot)write(*,*)'ERROR: Test find point failed, iNode=',&
         iNode,' instead of',nRoot

    if(.not.is_point_inside_node(CoordTest_D, iNode)) &
         write(*,*)'ERROR: Test find point failed'
    
    if(DoTestMe)write(*,*)'Testing distribute_tree 1st'
    call distribute_tree(.true.)
    if(DoTestMe)call show_tree('after distribute_tree 1st', .true.)

    if(DoTestMe)write(*,*)'Testing refine_tree_node'
    ! Refine the node where the point was found and find it again
    call refine_tree_node(iNode)

    if(DoTestMe)write(*,*)'Testing distribute_tree 2nd'

    ! Set node type to level+1
    allocate(iTypeNode_I(MaxNode))
    iTypeNode_I = 1 + iTree_IA(Level_,:)
    call distribute_tree(.true.,iTypeNode_I)
    if(DoTestMe)call show_tree('after distribute_tree with type=level', .true.)

    ! Set node type to the second coordinate index
    iTypeNode_I = iTree_IA(Coord2_,:)
    call distribute_tree(.true.,iTypeNode_I)
    if(DoTestMe)call show_tree('after distribute_tree with type=Coord2', .true.)

    ! Use default (single type)
    call distribute_tree(.true.)
    if(DoTestMe)call show_tree('after distribute_tree 2nd', .true.)

    call find_tree_node(CoordTest_D,iNode)
    if(.not.is_point_inside_node(CoordTest_D, iNode)) &
         write(*,*)'ERROR: Test find point failed for iNode=',iNode

    ! Refine another node
    if(DoTestMe)write(*,*)'nRoot=',nRoot
    call refine_tree_node(2)

    if(DoTestMe)call show_tree('after another refine_tree_node')

    if(DoTestMe)write(*,*)'Testing coarsen_tree_node'

    ! Coarsen back the last root node and find point again
    call coarsen_tree_node(nRoot)
    if(DoTestMe)call show_tree('after coarsen_tree_node')

    ! Distribute the new tree
    if(DoTestMe)write(*,*)'Testing distribute_tree 3rd'
    call distribute_tree(.true.)
    if(DoTestMe)call show_tree('after distribute_tree 3rd', .true.)

    call find_tree_node(CoordTest_D,iNode)
    if(iNode /= nRoot)write(*,*) &
         'ERROR: coarsen_tree_node+compact failed, iNode=',&
         iNode,' instead of',nRoot
    if(.not.is_point_inside_node(CoordTest_D, iNode)) &
         write(*,*)'ERROR: is_point_inside_node failed'

    if(iTree_IA(Status_, nNode+1) /= Unset_) &
         write(*,*)'ERROR: compact_tree failed, nNode=', nNode, &
         ' but status of next node is', iTree_IA(Status_, nNode+1), &
         ' instead of ', Unset_
    if(any(iTree_IA(Status_, 1:nNode) == Unset_)) &
         write(*,*)'ERROR: compact_tree faild, nNode=', nNode, &
         ' but iTree_IA(Status_, 1:nNode)=', &
         iTree_IA(Status_, 1:nNode),' contains unset=', Unset_
    call find_tree_node(CoordTest_D,iNode)
    if(iNode /= nRoot)write(*,*)'ERROR: compact_tree faild, iNode=',&
         iNode,' instead of',nRoot
    if(.not.is_point_inside_node(CoordTest_D, iNode)) &
         write(*,*)'ERROR: is_point_inside_node failed'

    if(DoTestMe)write(*,*)'Testing write_tree_file'
    call write_tree_file('tree.rst')

    if(DoTestMe)write(*,*)'Testing read_tree_file'
    iTree_IA = Unset_
    nRoot_D = 0
    call read_tree_file('tree.rst')
    if(DoTestMe)call show_tree('after read_tree')

    call find_tree_node(CoordTest_D,iNode)
    if(iNode /= nRoot)write(*,*)'ERROR: compact_tree failed, iNode=',&
         iNode,' instead of',nRoot

    if(DoTestMe)write(*,*)'Testing distribute_tree 4th'
    call distribute_tree(.true.)
    if(DoTestMe)call show_tree('after distribute_tree 4th', .true.)

    if(DoTestMe)write(*,*)'Testing clean_tree'
    call clean_tree
    if(DoTestMe)write(*,*)'MaxNode=', MaxNode

  end subroutine test_tree

end module BATL_tree
