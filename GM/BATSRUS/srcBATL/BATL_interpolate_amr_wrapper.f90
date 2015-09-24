module BATL_interpolate_amr_wrapper

  use BATL_mpi,  ONLY: iProc
  use BATL_size
  use BATL_tree, ONLY: Unset_, find_tree_node, Proc_, Block_,&
       get_tree_position, iTree_IA
  use BATL_geometry, ONLY: CoordMin_D, CoordMax_D, IsPeriodic_D

  use ModInterpolateAMR, ONLY: interpolate_amr


  implicit none

  SAVE
  private ! except

  public:: interpolate_amr_wrapper

  ! vector that keeps AMR directions
  logical,parameter:: IsAmr_D(MaxDim)= (/&
       iRatio_D(1) > 1,&
       iRatio_D(2) > 1,&
       iRatio_D(3) > 1/)

  ! non-AMR direction: 
  ! only 1 such direction, if 2 or more => interpolate_amr is not called
  ! (must be handled outside of wrapper
  ! MIN and MAX are added in order to keep value in range 1 to nDim
  integer, parameter:: iDimNoAmr = & 
       MAX(1,MIN(1*(2-iRatio) + 2*(2-jRatio) + 3*(2-kRatio), nDim))

  ! vector that keeps point's coordinate in non-AMR dimensions
  real   :: CoordNoAmr

  !\
  ! variables to keep track of nodes in the case of non-AMR direction
  ! IMPORTANT:
  !   order in which nodes are found is used as block id in interpolate_amr,
  !   the actual id is accessed as iNode_II(:,iNodeIdLocal):
  !   - iNode_II(1,:) are ids directly passed to interpolate_amr
  !   - iNode_II(2,:) are needed in the case of non-AMR directions:
  !     it stores ids of neighbors of corresponding iNode_II(1,:) if point
  !     falls close to the blocks boundary along non-AMR direction,
  !     otherwise it is the same as iNode_II(1,:)
  !-------------------------------------------------------------------------
  ! number of nodes found so far
  integer:: nNodeFound

  ! the actual node ids are stored here
  integer:: iNode_II(2,2**MaxDim) 

  ! cell indexes along non-AMR direction:
  ! iCellNoAmr_I(1) - for the cell which center is left  to the given point
  ! iCellNoAmr_I(2) - for the cell which center is right to the given point
  integer   :: iCellNoAmr_I(2)

  ! interpolation weight along non-AMR direction
  real      :: WeightNoAmr

contains

  subroutine interpolate_amr_wrapper(XyzIn_D, &
       nCell, iCell_II, Weight_I, IsSecondOrder)
    ! Find the grid cells surrounding the point Xyz_D.
    ! nCell returns the number of cells found on the processor.
    ! iCell_II returns the block+cell indexes for each cell.
    ! Weight_I returns the interpolation weights calculated 
    !                                 using AMR interpolation procedure
    ! IsSecondOrder returns whether the result is 2nd order interpolation
    real,    intent(in) :: XyzIn_D(MaxDim)
    integer, intent(out):: nCell
    integer, intent(out):: iCell_II(0:nDim,2**nDim)
    real,    intent(out):: Weight_I(2**nDim)
    logical, intent(out):: IsSecondOrder

    ! arguments for interpolate_amr subroutine
    integer:: nCell_D(MaxDim) = (/nI, nJ, nK/)
    integer:: iIndexes_II(0:nDimAmr+1,2**nDimAmr)
    integer:: nGridOut

    integer:: iProc_I(2**nDim)
    integer:: iNode !for code readability
    integer:: iCell, iDim, iDimAmr ! loop variables

    ! order of dimensions to correctly place AMR and non-AMR directions
    integer:: iDimOrder_I(nDim)
    !--------------------------------------------------------------------    
    ! coords along non-AMR direction
    if(nDimAmr < nDim) then
       ! this case is valid only for nDim=MaxDim=3
       CoordNoAmr = XyzIn_D(iDimNoAmr)
       iDimOrder_I = UNPACK((/1,2,3/), IsAmr_D(1:nDim), SPREAD(3,1,nDim))
    end if

    ! mark the beginning of the interpolation
    nNodeFound = 0
    iProc_I    = Unset_
    iNode_II   = Unset_

    ! apply general interpolate_amr routine but exclude non-AMR dimensions
    ! need to account for these dimensions separately,
    ! NOTE: routine doesn't sort out blocks located on other processors,
    !       this is done in the end of the current subroutine
    call interpolate_amr(&
         nDim          = nDimAmr, &
         XyzIn_D       = PACK(XyzIn_D, IsAmr_D), &
         nIndexes      = nDimAmr+1, &
         find          = find, &
         nCell_D       = PACK(nCell_D, IsAmr_D), &
         nGridOut      = nGridOut, &
         Weight_I      = Weight_I, &
         iIndexes_II   = iIndexes_II, &
         IsSecondOrder = IsSecondOrder, &
         UseGhostCell  = .false.)

    ! copy results indices
    iProc_I(            1:nGridOut) = iIndexes_II(0,         1:nGridOut)
    iCell_II(1:nDimAmr, 1:nGridOut) = iIndexes_II(1:nDimAmr, 1:nGridOut)
    iCell_II(0,         1:nGridOut) = iIndexes_II(nDimAmr+1, 1:nGridOut)

    ! check if there is a non-AMR direction, handle it if present
    if(nDimAmr < nDim) then
       ! local node ids & cell indexes along AMR directions
       ! are the same for cells iCell and nGridOut+iCell
       iCell_II(0:nDimAmr,nGridOut+1:2*nGridOut)=iCell_II(0:nDimAmr,1:nGridOut)

       ! store indexes and restore shape along non-AMR direction
       iCell_II(  nDim,         1:  nGridOut) = iCellNoAmr_I(1)
       iCell_II(  nDim,nGridOut+1:2*nGridOut) = iCellNoAmr_I(2)
       iCell_II(1:nDim,1:2*nGridOut) = iCell_II(iDimOrder_I,1:2*nGridOut)

       ! correct weights
       Weight_I(nGridOut+1:2*nGridOut)= Weight_I(1:nGridOut) *    WeightNoAmr
       Weight_I(         1:  nGridOut)= Weight_I(1:nGridOut) * (1-WeightNoAmr)

       ! get correct block and processor ids
       do iCell = 1, nGridOut
          ! proper global node, block and proc ids for iCell
          iNode = iNode_II(1,iCell_II(0,iCell))
          if(iNode > 0)then
             iCell_II(0, iCell) = iTree_IA(Block_,iNode)
             iProc_I(    iCell) = iTree_IA(Proc_, iNode)
          end if
          ! proper global node, block and proc ids for nGridOut+iCell
          iNode = iNode_II(2,iCell_II(0,nGridOut+iCell))
          if(iNode > 0)then
             iCell_II(0, nGridOut+iCell) = iTree_IA(Block_,iNode)
             iProc_I(    nGridOut+iCell) = iTree_IA(Proc_, iNode)
          end if
       end do

       ! number for cells has doubled
       nGridOut = 2*nGridOut
    end if


    ! sort out cells located on other processors
    nCell = 0
    do iCell = 1, nGridOut
       if(iProc_I(iCell) == iProc)then
          nCell = nCell + 1
          iCell_II(:, nCell) = iCell_II(:, iCell)
          Weight_I(nCell) = Weight_I(iCell)
       end if
    end do

  end subroutine interpolate_amr_wrapper

  !============================================================================

  subroutine find(nDimIn, Coord_D, &
       iProc, iBlock, CoordCorner_D, DCoord_D, IsOut)
    integer, intent(in) :: nDimIn
    !\
    ! "In"- the coordinates of the point, "out" the coordinates of the
    ! point with respect to the block corner. In the most cases
    ! XyzOut_D = XyzIn_D - XyzCorner_D, the important distinction,
    ! however, is the periodic boundary, near which the jump in the
    ! stencil coordinates might occur. To handle the latter problem,
    ! we added the "out" intent. The coordinates for the stencil
    ! and input point are calculated and recalculated below with
    ! respect to the block corner.
    !/
    real,  intent(inout):: Coord_D(nDimIn)
    integer, intent(out):: iProc, iBlock !processor and block number
    !\
    ! Block left corner coordinates and the grid size:
    !
    real,    intent(out):: CoordCorner_D(nDimIn), DCoord_D(nDimIn)
    logical, intent(out):: IsOut !Point is out of the domain.

    real   :: CoordFull_D(MaxDim)       ! Full Gen coords of point
    real   :: CoordCornerFull_D(MaxDim) ! Full Gen coords of corner
    real   :: DCoordFull_D(MaxDim)      ! Full mesh sizes in gen coords
    real   :: CoordTree_D(MaxDim)       ! Normalized gen coords
    integer   :: iDim, iDimAmr ! loop variable
    integer   :: iNode         ! tree node number
    real,dimension(MaxDim):: PositionMin_D, PositionMax_D
    !\
    ! Used in normalization
    !/
    real   :: CoordTreeSize_D(MaxDim) 
    !\
    ! In the case of a non-AMR direction shows, where in iNodeFound_II
    ! to store a newly found node and its neighbor
    !/
    integer, save:: iFirst = 1, iSecond = 2
    !\
    ! in the case of non-AMR direction: displacement towards neighbor
    !/
    real:: Displace = 0.0
    !--------------------------------------------------------------------
    CoordTreeSize_D = CoordMax_D - CoordMin_D
    ! Coord_D will be used to find block, copy input data into it
    ! restore non-AMR directions if necessary
    CoordFull_D = 0
    if(nDimAmr == nDim)then !NOTE: nDimIn = nDimAmr
       CoordFull_D(1:nDim) = Coord_D(1:nDim)
    else
       CoordFull_D = UNPACK(Coord_D, IsAmr_D, SPREAD(CoordNoAmr,1,MaxDim))
    end if
    !\
    ! Calculate normalized per (CoordMax_D-CoordMin_D) coordinates 
    ! for tree search
    !/
    CoordTree_D = (CoordFull_D - CoordMin_D)/CoordTreeSize_D
    !\
    ! For periodic boundary conditions fix the input coordinate if
    ! beyond the tree bounadaries
    !/
    where(IsPeriodic_D)CoordTree_D = modulo(CoordTree_D, 1.0)
    !\
    ! call internal BATL find subroutine
    !/
    call find_tree_node(CoordIn_D=CoordTree_D, iNode=iNode)
  
    !\
    ! Check if the block is found
    !/
    if(iNode<=0)then
       IsOut = .true.
       RETURN
    end if
    !\
    !position has been found
    !/
    IsOut = .false.
    iBlock = iTree_IA(Block_,iNode)
    iProc  = iTree_IA(Proc_, iNode)
    call get_tree_position(iNode=iNode,                &
                           PositionMin_D=PositionMin_D,&
                           PositionMax_D=PositionMax_D )
    CoordCornerFull_D =  CoordMin_D + PositionMin_D * CoordTreeSize_D
    DCoordFull_D      =  (PositionMax_D-PositionMin_D)*CoordTreeSize_D/nIjk_D
    ! corner coordinates of the found block
    CoordCorner_D = PACK(CoordCornerFull_D, IsAmr_D)
    ! cell size of the found block
    DCoord_D = PACK(DCoordFull_D, IsAmr_D)

    ! subtract coordinates of the corner from point's coordinates
    Coord_D = PACK((CoordTree_D - PositionMin_D)*CoordTreeSize_D, IsAmr_D)

    ! if there is no non-AMR direction => no need to do extra work
    if(nDimAmr == nDim) RETURN !NOTE: nDimIn = nDimAmr

    !\
    ! take care of the non-AMR direction:
    ! store in global variables weight along the direction and cell indexes
    !/   
    if(nNodeFound == 0)then
       ! compute both once in the beginning
       WeightNoAmr = 0.5 + nIJK_D(iDimNoAmr)*&
            (CoordTree_D(  iDimNoAmr) - PositionMin_D(iDimNoAmr)) /&
            (PositionMax_D(iDimNoAmr) - PositionMin_D(iDimNoAmr))
       iCellNoAmr_I(1) = floor(WeightNoAmr)
       WeightNoAmr     = WeightNoAmr - iCellNoAmr_I(1)      
       iCellNoAmr_I(2) = iCellNoAmr_I(1) + 1  

       ! iFirst corresponds to node found above, to be stored in iNode_II(1,:),
       ! iSecond will be searched for below, to be stored in iNode_II(2,:)
       iFirst = 1; iSecond = 2

       ! take care of the cases when point is close to the boundary
       if    (iCellNoAmr_I(1) == 0)then
          iCellNoAmr_I(1) = nIJK_D(iDimNoAmr)
          ! inverse order in which nodes in iNode_II are found
          iFirst = 2; iSecond = 1 
          ! displacement towards neighbor: down
          Displace = - DCoordFull_D(iDimNoAmr) / CoordTreeSize_D(iDimNoAmr)
       elseif(iCellNoAmr_I(2) == nIJK_D(iDimNoAmr)+1)then
          iCellNoAmr_I(2) = 1
          ! displacement towards neighbor: up
          Displace = + DCoordFull_D(iDimNoAmr) / CoordTreeSize_D(iDimNoAmr)
       end if
    end if

    ! store parameters of the node and return local node id as block id
    nNodeFound = nNodeFound + 1
    iNode_II(iFirst,nNodeFound) = iNode
    iBlock = nNodeFound

    ! check if point is close to the boundary of the block
    if(iCellNoAmr_I(1) > iCellNoAmr_I(2))then
       ! need to find a neighboring node
       CoordTree_D(iDimNoAmr) = CoordTree_D(iDimNoAmr) + Displace
       call find_tree_node(CoordIn_D=CoordTree_D, iNode=iNode)   
    end if

    ! store the second node
    iNode_II(iSecond, nNodeFound) = iNode

  end subroutine find
  !============
end module BATL_interpolate_amr_wrapper
