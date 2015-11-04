module BATL_interpolate_amr

  use BATL_mpi,  ONLY: iProc
  use BATL_size, ONLY: nIJK_D, nDim, nDimAmr, MaxDim, &
       iRatio, jRatio, kRatio, iRatio_D
  use BATL_tree, ONLY: Unset_, find_tree_node, Proc_, Block_,&
       get_tree_position, iTree_IA
  use BATL_geometry, ONLY: CoordMin_D, CoordMax_D, IsPeriodic_D

  use ModInterpolateAMR, ONLY: interpolate_amr_shared=>interpolate_amr, cTol2


  implicit none

  SAVE
  private ! except

  public:: interpolate_amr, find_block_to_interpolate_gc

  ! non-AMR direction: 
  ! only 1 such direction, if 2 or more => interpolate_amr is not called
  ! (must be handled outside this module)
  ! MIN and MAX are added in order to keep value in range 1 to nDim
  integer, parameter:: iDimNoAmr = & 
       MAX(1,MIN(1*(2 - iRatio) + 2*(2 - jRatio) + 3*(2 - kRatio), nDim))

  ! order of dimensions to correctly place AMR and non-AMR directions
  ! it depends on iDimNoAmr:
  ! iOrder_I = (/2, 3, 1/) for iDimNoAmr = 1
  ! iOrder_I = (/1, 3, 2/) for iDimNoAmr = 2
  ! iOrder_I = (/1, 2, 3/) for iDimNoAmr = 3
  ! iOrder_I = (/1, 2, 3/) for nDimAmr = nDim
  ! MIN and MAX are added in order to keep value in range 1 to nDim
  integer,parameter:: iOrder_I(MaxDim) = (/&
       3 - iRatio, & 
       MIN(6 - iRatio - jRatio, nDim),  &
       MAX(-3 + 2*iRatio + jRatio, 1) /)

  ! point's coordinate in non-AMR dimensions
  real   :: CoordNoAmr, DisplacedCoordNoAmr

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
  logical   :: IsSecondOrderNoAmr

contains

  subroutine interpolate_amr(XyzIn_D, &
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

    integer:: iIndexes_II(0:nDimAmr+1,2**nDimAmr)
    integer:: nGridOut

    integer:: iProc_I(2**nDim)
    integer:: iNode !for code readability
    integer:: iCell ! loop variables
    !--------------------------------------------------------------------    
    ! coords along non-AMR direction
    if(nDimAmr < nDim) then
       ! this case is valid only for nDim=MaxDim=3, nDimAmr=2
       CoordNoAmr = XyzIn_D(iDimNoAmr)
    end if

    ! mark the beginning of the interpolation
    nNodeFound = 0
    iProc_I    = Unset_
    iNode_II   = Unset_

    ! apply general interpolate_amr routine but exclude non-AMR dimensions
    ! need to account for these dimensions separately,
    ! NOTE: routine doesn't sort out blocks located on other processors,
    !       this is done in the end of the current subroutine
    call interpolate_amr_shared(&
         nDim          = nDimAmr, &
         XyzIn_D       = XyzIn_D(iOrder_I(1:nDimAmr)), &
         nIndexes      = nDimAmr+1, &
         find          = find, &
         nCell_D       = nIJK_D(iOrder_I(1:nDimAmr)), &
         nGridOut      = nGridOut, &
         Weight_I      = Weight_I, &
         iIndexes_II   = iIndexes_II, &
         IsSecondOrder = IsSecondOrder, &
         UseGhostCell  = .false.)
    if(nGridOut==0)RETURN

    ! copy results indices
    iProc_I(            1:nGridOut) = iIndexes_II(0,         1:nGridOut)
    iCell_II(1:nDimAmr, 1:nGridOut) = iIndexes_II(1:nDimAmr, 1:nGridOut)
    iCell_II(0,         1:nGridOut) = iIndexes_II(nDimAmr+1, 1:nGridOut)

    ! check if there is a non-AMR direction, 
    if(nDimAmr == nDim)then
       call sort_out_other_procs
       RETURN
    end if

    !Handle the non-AMR direction
    IsSecondOrder = IsSecondOrder .and. IsSecondOrderNoAmr
    
    ! store indexes and restore shape along non-AMR direction
    iCell_II(  nDim,         1:  nGridOut) = iCellNoAmr_I(1)
    if(WeightNoAmr < cTol2)then
       iCell_II(iOrder_I(1:nDim),1:nGridOut) = iCell_II(1:nDim,1:nGridOut)
       call sort_out_other_procs
       RETURN
    end if
    ! local node ids & cell indexes along AMR directions
    ! are the same for cells iCell and nGridOut+iCell
    iCell_II(0:nDimAmr,nGridOut+1:2*nGridOut) = iCell_II(0:nDimAmr,1:nGridOut)
    iCell_II(     nDim,nGridOut+1:2*nGridOut) = iCellNoAmr_I(2)
    iCell_II(iOrder_I(1:nDim),1:2*nGridOut)   = iCell_II(1:nDim,1:2*nGridOut)

    ! correct weights
    Weight_I(nGridOut+1:2*nGridOut)= Weight_I(1:nGridOut) *    WeightNoAmr
    Weight_I(         1:  nGridOut)= Weight_I(1:nGridOut)*(1 - WeightNoAmr)
    if(iCellNoAmr_I(1) > iCellNoAmr_I(2))then    
       ! get correct block and processor ids
       do iCell = 1, nGridOut
          ! proper global node, block and proc ids for iCell
          iNode = iNode_II(1,iCell_II(0,iCell))
          iCell_II(0, iCell) = iTree_IA(Block_,iNode)
          iProc_I(    iCell) = iTree_IA(Proc_, iNode)
          ! proper global node, block and proc ids for nGridOut+iCell
          iNode = iNode_II(2,iCell_II(0,nGridOut+iCell))
          iCell_II(0, nGridOut+iCell) = iTree_IA(Block_,iNode)
          iProc_I(    nGridOut+iCell) = iTree_IA(Proc_, iNode)
       end do
    end if
    ! number of cells has doubled
    nGridOut = 2*nGridOut
    call sort_out_other_procs
  contains
    subroutine sort_out_other_procs
      ! sort out cells located on other processors
      nCell = 0
      do iCell = 1, nGridOut
         if(iProc_I(iCell) == iProc)then
            nCell = nCell + 1
            iCell_II(:, nCell) = iCell_II(:, iCell)
            Weight_I(nCell) = Weight_I(iCell)
         end if
      end do
    end subroutine sort_out_other_procs

  end subroutine interpolate_amr
  !============================
  subroutine fix_tree_coord(CoordTree_D)
    use BATL_geometry, ONLY: IsLatitudeAxis, IsSphericalAxis, &
         IsCylindricalAxis
    real,intent(inout) :: CoordTree_D(MaxDim)  !Normalized gen coords
    !-----------
    !\
    ! For periodic boundary conditions fix the input coordinate if
    ! beyond the tree bounadaries
    !/
    where(IsPeriodic_D)CoordTree_D = modulo(CoordTree_D, 1.0)
    !\
    ! Check specific boundary conditions for particular geometries
    !/
    if(IsLatitudeAxis)then
       !\
       !spherical: r, lon, lat coordinates
       !/
       if(CoordTree_D(3) > 1.0)then
          !\
          ! reflect third coordinate, 
          ! add half of full range to the second one.
          !/
          CoordTree_D(3) = 2.0 - CoordTree_D(3) 
          CoordTree_D(2) = modulo(CoordTree_D(2) + 0.50, 1.0)
       elseif(CoordTree_D(3) < 0.0)then
          !\
          ! reflect third coordinate, 
          ! add half of full range to the second one.
          !/
          CoordTree_D(3) = -CoordTree_D(3) 
          CoordTree_D(2) = modulo(CoordTree_D(2) + 0.50, 1.0)
       end if
    elseif(IsSphericalAxis)then
       !\
       ! spherical: r, theta, phi
       !/
       if(CoordTree_D(2) > 1.0)then
          !\
          ! reflect second coordinate, 
          ! add half of full range to the third one.
          !/
          CoordTree_D(2) = 2.0 - CoordTree_D(2) 
          CoordTree_D(3) = modulo(CoordTree_D(3) + 0.50, 1.0)
       elseif(CoordTree_D(2) < 0.0)then
          !\
          ! reflect second coordinate, 
          ! add half of full range to the third one.
          !/
          CoordTree_D(2) = -CoordTree_D(2) 
          CoordTree_D(3) = modulo(CoordTree_D(3) + 0.50, 1.0)
       end if
    elseif(IsCylindricalAxis)then
       !\
       ! cylindrical: r, phi, z
       !/
       if(CoordTree_D(1) < 0.0)then
          !\
          ! reflect first coordinate, 
          ! add half of full range to the second one.
          !/
          CoordTree_D(1) = -CoordTree_D(1) 
          CoordTree_D(2) = modulo(CoordTree_D(2) + 0.50, 1.0)
       end if
    end if
  end subroutine fix_tree_coord
  !=======================
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
    integer, save:: First_ = 1, Second_ = 2
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
       CoordFull_D(iOrder_I(1:nDimAmr)) = Coord_D
       CoordFull_D(iDimNoAmr) = CoordNoAmr
    end if
    !\
    ! Calculate normalized per (CoordMax_D-CoordMin_D) coordinates 
    ! for tree search
    !/
    CoordTree_D = (CoordFull_D - CoordMin_D)/CoordTreeSize_D
    call fix_tree_coord(CoordTree_D)
    IsOut = any(CoordTree_D < 0.0 .or. CoordTree_D >= 1.0)
    if(IsOut)RETURN
    !\
    ! call internal BATL find subroutine
    !/
    call find_tree_node(CoordIn_D=CoordTree_D, iNode=iNode)
  
    !\
    ! Check if the block is found
    !/
    if(iNode<=0)then
       call CON_stop('Failure in BATL_interpolate_amr:find')
    end if
    !\
    !position has been found
    !/
    iBlock = iTree_IA(Block_,iNode)
    iProc  = iTree_IA(Proc_, iNode)
    call get_tree_position(iNode=iNode,                &
                           PositionMin_D=PositionMin_D,&
                           PositionMax_D=PositionMax_D )
    ! corner coordinates for the found block
    CoordCornerFull_D =  CoordMin_D + PositionMin_D * CoordTreeSize_D
    CoordCorner_D = CoordCornerFull_D(iOrder_I(1:nDimAmr))

    ! cell size for the found block
    DCoordFull_D      =  (PositionMax_D - PositionMin_D)*CoordTreeSize_D/nIjk_D
    DCoord_D = DCoordFull_D(iOrder_I(1:nDimAmr))

    ! subtract coordinates of the corner from point's coordinates
    Coord_D = (CoordTree_D(iOrder_I(1:nDimAmr)) - &
              PositionMin_D(iOrder_I(1:nDimAmr)))*&
              CoordTreeSize_D(iOrder_I(1:nDimAmr))

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

       ! First_ corresponds to node found above, to be stored in iNode_II(1,:),
       ! Second_ will be searched for below, to be stored in iNode_II(2,:)
       First_ = 1; Second_ = 2
       IsSecondOrderNoAmr = .true.
       ! take care of the cases when point is close to the boundary
       if    (iCellNoAmr_I(1) == 0)then
          ! displacement towards neighbor: down
          DisplacedCoordNoAmr = CoordTree_D(iDimNoAmr) - &
               DCoordFull_D(iDimNoAmr) / CoordTreeSize_D(iDimNoAmr)
          if(IsPeriodic_D(iDimNoAmr))&
               DisplacedCoordNoAmr = &
               modulo(DisplacedCoordNoAmr, 1.0)
          if(DisplacedCoordNoAmr < 0.0)then
             IsSecondOrderNoAmr = .false.
             WeightNoAmr = 0.0
             iCellNoAmr_I(1) = iCellNoAmr_I(2)
             RETURN
          end if
          iCellNoAmr_I(1) = nIJK_D(iDimNoAmr)
          ! inverse order in which nodes in iNode_II are found
          First_ = 2; Second_ = 1 
       elseif(iCellNoAmr_I(2) == nIJK_D(iDimNoAmr)+1)then
        ! displacement towards neighbor: up
          DisplacedCoordNoAmr = CoordTree_D(iDimNoAmr) + &
               DCoordFull_D(iDimNoAmr) / CoordTreeSize_D(iDimNoAmr)
          if(IsPeriodic_D(iDimNoAmr))&
               DisplacedCoordNoAmr = &
               modulo(DisplacedCoordNoAmr, 1.0)
          if(DisplacedCoordNoAmr >= 1.0)then
             IsSecondOrderNoAmr = .false.
             WeightNoAmr = 0.0
             iCellNoAmr_I(2) = iCellNoAmr_I(1)
             RETURN
          end if  
          iCellNoAmr_I(2) = 1
       end if
    end if

    nNodeFound = nNodeFound + 1

    ! check if point is close to the boundary of the block
    if(iCellNoAmr_I(1) > iCellNoAmr_I(2))then
       ! store parameters of the node and return local node id as block id
       iNode_II(First_,nNodeFound) = iNode
       iBlock = nNodeFound
       ! need to find a neighboring node
       CoordTree_D(iDimNoAmr) =  DisplacedCoordNoAmr
       call find_tree_node(CoordIn_D=CoordTree_D, iNode=iNode)
       !\
       ! Check if the block is found
       !/
       if(iNode<=0)then
          call CON_stop('Failure in BATL_interpolate_amr:find')
       end if
       ! store the second node
       iNode_II(Second_, nNodeFound) = iNode   
    end if
  end subroutine find
 !======================================
  subroutine find_block_to_interpolate_gc(Coord_D, iPeOut, iBlockOut)
    real,    intent(in)  :: Coord_D(MaxDim)
    integer, intent(out) :: iPeOut, iBlockOut
    character(LEN=*),parameter:: NameSub = 'find_block_to_interpolate_gc'
    real      :: CoordTree_D(MaxDim)    ! Normalized gen coords
    integer   :: iNode                  ! Tree node number
    integer   :: iDim                   ! Loop 
    real,dimension(MaxDim):: PositionMin_D, PositionMax_D
    real      :: DCoord_D(nDim)         ! Cell size 
    real      :: DCoordInv_D(nDim)      ! Cell size inv
    real      :: CoordCorner_D(nDim)    ! Coords of the block left corner
    !\
    ! Used in normalization
    !/
    real   :: CoordTreeSize_D(MaxDim) 
    !\
    !Direction along which the point goes out of the block inner part
    !/
    integer :: iShift_D(nDim)
    !\
    ! Grid, consisting of Coord_D, displaced by +/- DCoord_D 
    !/
    real    :: GridCoord_DI(MaxDim,2**nDimAmr)
    !\
    ! The grid point belonging to the initially found block
    !/
    integer   :: iGridBlock
    !\
    ! Number of grid points:
    !/
    integer   :: nGrid
    !\
    ! Loop variable
    !/
    integer   :: iGrid
    !\
    ! Refinement level of the neighboring block
    !/
    integer   :: iLevel
    !\
    ! Dimension, along which there is AMR
    !/
    integer:: iDimAmr
    !-------
    CoordTreeSize_D = CoordMax_D - CoordMin_D
    !\
    ! Calculate normalized per (CoordMax_D-CoordMin_D) coordinates 
    ! for tree search
    !/
    CoordTree_D = (Coord_D - CoordMin_D)/CoordTreeSize_D
        !\
    ! call internal BATL find subroutine
    !/
    call find_tree_node(CoordIn_D=CoordTree_D, iNode=iNode)
  
    !\
    ! Check if the block is found
    !/
    if(iNode<=0)then
       call CON_stop('Failure in '//NameSub)
    end if
    !\
    !position has been found
    !/
    iBlockOut = iTree_IA(Block_,iNode)
    iPeOut    = iTree_IA(Proc_, iNode)
    call get_tree_position(iNode=iNode,                &
                           PositionMin_D=PositionMin_D,&
                           PositionMax_D=PositionMax_D )
    ! corner coordinates for the found block
    CoordCorner_D =  CoordMin_D(1:nDim) + &
         PositionMin_D(1:nDim) * CoordTreeSize_D(1:nDim)

    ! cell size for the found block
    DCoord_D      =  (PositionMax_D(1:nDim) - PositionMin_D(1:nDim))&
         *CoordTreeSize_D(1:nDim)/nIjk_D(1:nDim)
    DCoordInv_D = 1/DCoord_D 
    !\
    ! Check if the block is suitable to interpolate with ghost cells
    !/
    
    iShift_D = floor((Coord_D(1:nDim) - CoordCorner_D  - 0.50*DCoord_D)*&
         DCoordInv_D/(nIJK_D(1:nDim) - 1))
    if(all(iShift_D(1:nDim)==0))RETURN
    !\
    ! Fill in grid point coordinates
    !/
    iGridBlock = 1 ; nGrid = 1; GridCoord_DI(:,1) = Coord_D
    do iDim = 1, nDim
       if(iRatio_D(iDim) < 2 & !No amr
            .or. iShift_D(iDim) ==0)CYCLE !No grid shift in this direction
       iDimAmr = iDim
       GridCoord_DI(:,nGrid+1:2*nGrid) = GridCoord_DI(:,1:nGrid)
       if(iShift_D(iDim)==-1)then
          GridCoord_DI(iDim,1:nGrid) =  &
               GridCoord_DI(iDim,1:nGrid) - DCoord_D(iDim)
          iGridBlock = iGridBlock + nGrid
       else !iShift_D(iDim)==+1
          GridCoord_DI(iDim,nGrid+1:2*nGrid) =  &
               GridCoord_DI(iDim,nGrid+1:2*nGrid) + DCoord_D(iDim)
       end if
       nGrid = 2*nGrid
    end do
    do iGrid = 1, nGrid
       if(iGrid==iGridBlock)CYCLE
       !\
       !Find tree node into which the displaced grid point falls
       !/
       !\
       ! Calculate normalized per (CoordMax_D-CoordMin_D) coordinates 
       ! for tree search
       !/
       CoordTree_D = (Coord_D - CoordMin_D)/CoordTreeSize_D
       call fix_tree_coord(CoordTree_D)
       !\
       ! Check if the grid point is out of domain
       !/
       if(any(CoordTree_D < 0.0 .or. CoordTree_D >= 1.0))CYCLE
       !\
       ! call internal BATL find subroutine
       !/
       call find_tree_node(CoordIn_D=CoordTree_D, iNode=iNode)
       !\
       ! Check if the block is found
       !/
       if(iNode<=0)then
          call CON_stop('Failure in '//NameSub)
       end if
       call get_tree_position(iNode=iNode,                &
                           PositionMin_D=PositionMin_D,&
                           PositionMax_D=PositionMax_D )

       ! cell size for the found block
       DCoord_D      =  (PositionMax_D(1:nDim) - PositionMin_D(1:nDim))&
            *CoordTreeSize_D(1:nDim)/nIjk_D(1:nDim)
       iLevel = 1 - floor(DCoord_D(iDimAmr)*DCoordInv_D(iDimAmr)+ 0.001)
       !\                     ^
       ! For expression above | equal to 2 , 1, 0.5 correspondingly
       ! iLevel = -1, 0, 1, meaning that the neighboring block
       ! is coarser, at the same resolution or finer than the basic 
       ! one.
       !/
       if(iLevel==-1)RETURN !Neighbor is coarser
       if(iLevel==1)then    !Finer neighbor is chosen
          iBlockOut = iTree_IA(Block_,iNode)
          iPeOut    = iTree_IA(Proc_, iNode)
          RETURN
       end if
    end do
  end subroutine find_block_to_interpolate_gc
end module BATL_interpolate_amr
