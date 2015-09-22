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
  logical:: IsAmr_D(MaxDim)

  ! vector that keeps point's coordinate in non-AMR dimensions
  real   :: XyzNonAmr_D(MaxDim), XyzNei_D(MaxDim)

  ! variables to keep track of blocks and their parameters
  ! IMPORTANT:
  !   routine uses order in which blocks are found as block id,
  !   the actual id is accessed as iBlock_I(iBlockIdLocal)
  integer:: nBlockFound
  integer:: iBlock_I(2*2**MaxDim) ! the actual block ids are stored here
  real   :: XyzCorner_DI(MaxDim, 2*2**MaxDim)
  real   :: Dxyz_DI(MaxDim, 2*2**MaxDim)

  !\
  ! Index for the cell which center is left to the given point
  !/
  integer   :: iCellFloor_D(MaxDim)     
  !\
  ! Normalized distance from the point 
  ! to the cell with the index iCellsFloor_D
  !/
  real      :: Dist_D(MaxDim)



  ! this block is needed to avoid circular dependency:
  ! - interpolate_amr_wrapper, find_block use routines from module BATL_grid;
  ! - since find_block is called by shared AMR interpolation routine,
  !   its interface can't be changed and it can use routines from BATL_grid
  !   only if they are stored in pointer kept as module variables;
  ! - these pointers are set in interpolate_amr_wrapper
  interface
     !\
     ! interface for find_grid_block procedure from BATL_grid module:
     ! detailed description can be found in BATL_grid.f90
     !/
     subroutine find_grid_block_interface(XyzIn_D, &
          iProcOut, iBlockOut, iCellOut_D, DistOut_D, iNodeOut, &
          CoordMinBlockOut_D, CoordMaxBlockOut_D, CellSizeOut_D, & 
          UseGhostCell)

       use BATL_size

       implicit none
       real,    intent(in) :: XyzIn_D(MaxDim)       ! Cartesian coords of point
       integer, intent(out):: iBlockOut, iProcOut   ! Block and proc indexes
       integer, intent(out), optional:: iCellOut_D(MaxDim)!Closest cell indexes
       real,    intent(out), optional:: DistOut_D(MaxDim) ! Normalized distance
       integer, intent(out), optional:: iNodeOut     ! Tree node index
       real,    intent(out), optional:: CoordMinBlockOut_D(MaxDim)!block corner
       real,    intent(out), optional:: CoordMaxBlockOut_D(MaxDim)!block corner
       real,    intent(out), optional:: CellSizeOut_D(MaxDim)!cellsize in block
       logical, intent(in),  optional:: UseGhostCell ! Use ghost cells or not
     end subroutine find_grid_block_interface
     !=========================================================================
     !\
     ! interface for interpolate_grid procedure from BATL_grid module
     ! detailed description can be found in BATL_grid.f90
     subroutine interpolate_grid_interface(Xyz_D, nCell, iCell_II, Weight_I)

       use BATL_size

       implicit none
       real,    intent(in) :: Xyz_D(MaxDim)
       integer, intent(out):: nCell  
       integer, intent(out):: iCell_II(0:nDim,2**nDim)
       real,    intent(out):: Weight_I(2**nDim)
     end subroutine interpolate_grid_interface
  end interface

  ! pointers keep routines from BATL_grid passed to interpolate_amr_wrapper
  procedure(find_grid_block_interface), pointer::find_grid_block_ptr  => NULL()
  procedure(interpolate_grid_interface),pointer::interpolate_grid_ptr => NULL()

contains

  subroutine interpolate_amr_wrapper(XyzIn_D, &
       nCell, iCell_II, Weight_I, &
       find_grid_block, &
       interpolate_grid)
    ! Find the grid cells surrounding the point Xyz_D.
    ! nCell returns the number of cells found on the processor.
    ! iCell_II returns the block+cell indexes for each cell.
    ! Weight_I returns the interpolation weights calculated 
    !                                 using AMR interpolation procedure
    real,    intent(in) :: XyzIn_D(MaxDim)
    integer, intent(out):: nCell
    integer, intent(out):: iCell_II(0:nDim,2**nDim)
    real,    intent(out):: Weight_I(2**nDim)
    procedure(find_grid_block_interface) :: find_grid_block
    procedure(interpolate_grid_interface):: interpolate_grid

    ! arguments for interpolate_amr subroutine
    integer:: nCell_D(MaxDim) = (/nI, nJ, nK/)
    integer:: iIndexes_II(0:nDimAmr+1,2**nDimAmr)
    integer:: nGridOut
    logical:: IsSecondOrder

    integer:: iProc_I(2**nDim)
    integer:: iBlockNei, iProcNei  
    integer:: iCell, iDim, iDimAmr ! loop variables

    real   :: Misc
    integer:: i_D(MaxDim)
    !--------------------------------------------------------------------    
    ! set pointers to the passed routines
    find_grid_block_ptr  => find_grid_block
    interpolate_grid_ptr => interpolate_grid

    ! the grid is simple => no need to use sophisticated AMR routine
    if(nDimAmr==1)then
       call interpolate_grid_ptr(XyzIn_D, nCell, iCell_II, Weight_I)
       RETURN
    end if

    ! store coords in Non-AMR dimenstions
    IsAmr_D = iRatio_D > 1
    where(IsAmr_D)
       XyzNonAmr_D = 0
    elsewhere
       XyzNonAmr_D = XyzIn_D
    end where

    ! mark the beginning of the interpolation
    nBlockFound = 0
    iProc_I     = Unset_
    iBlock_I    = Unset_

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

    ! check if point is close to the domain's boundary:
    ! in this case the result isn't 2nd order accurate, handle it separately
    if(.not. IsSecondOrder)then
       call interpolate_grid_ptr(XyzIn_D, nCell, iCell_II, Weight_I)
       RETURN
    end if

    ! copy results indices
    iProc_I(            1:nGridOut) = iIndexes_II(0,         1:nGridOut)
    iCell_II(1:nDimAmr, 1:nGridOut) = iIndexes_II(1:nDimAmr, 1:nGridOut)
    iCell_II(0,         1:nGridOut) = iIndexes_II(nDimAmr+1, 1:nGridOut)


    ! check if there are any non-AMR directions, handle them if present
    if(nDimAmr < nDim) then
       ! restore correct shape and account for non-AMR dimensions if necessary
       iDimAmr = nDimAmr
       do iDim = nDim, 1, -1
          if(.not. IsAmr_D(iDim)) CYCLE
          iCell_II(iDim,:) = iCell_II(iDimAmr,:)
          iDimAmr = iDimAmr - 1
       end do

       ! find weights along non-AMR dimensions
       do iDim = 1, nDim
          if(IsAmr_D(iDim)) CYCLE

          ! copy cell indices for the other half of points
          iCell_II(:,nGridOut+1:2*nGridOut) = iCell_II(:,1:nGridOut)

          ! normalized coords with respect to block's corner
          ! if 1 <= Misc <= nIJK_D(iDim) then point is far enough from boundary
          ! otherwise need to find neighbors
          Misc = 0.5 + &
               ( XyzIn_D(iDim) - XyzCorner_DI(iDim, iCell_II(0,1)) ) / &
               Dxyz_DI(iDim, iCell_II(0,1))

          ! indices of the cell near the point
          iCell_II(iDim, :) = floor(Misc)

          ! interpolation weight along non-AMR direction
          Misc = Misc - iCell_II(iDim,1)

          if    ( iCell_II(iDim,1) == 0 )then
             ! the point is close to the bottom boundary of block
             ! correct cell indices along non-AMR dimension
             iCell_II(iDim,         1:  nGridOut) = 1
             iCell_II(iDim,nGridOut+1:2*nGridOut) = nIJK_D(iDim)

             ! correct weights
             Weight_I(nGridOut+1:2*nGridOut) = Weight_I(1:nGridOut)*(1. - Misc)
             Weight_I(         1:  nGridOut) = Weight_I(1:nGridOut)*      Misc

             ! cycle through already found cells inside current block
             ! find their neighbors
             do iCell = 1, nGridOut
                ! location of the neighbor: current cells + displacement
                i_D = 1
                i_D(1:nDim) = iCell_II(1:nDim, iCell)
                XyzNei_D = XyzCorner_DI(:,iCell_II(0,iCell)) + &
                     ((/i_D(1), i_D(2),i_D(3)/) - 0.5) * &
                     Dxyz_DI(:,iCell_II(0,iCell))
                XyzNei_D(iDim) = XyzNei_D(iDim)-Dxyz_DI(iDim,iCell_II(0,iCell))

                call find_grid_block_ptr(XyzNei_D, iProcNei, iBlockNei)
                ! write correct processors and blocks
                iProc_I(    nGridOut+iCell) = iProcNei
                nBlockFound = nBlockFound + 1
                iBlock_I(nBlockFound) = iBlockNei
                iCell_II(0, nGridOut+iCell) = nBlockFound
             end do

          elseif( iCell_II(iDim,1) == nIJK_D(iDim) )then
             ! the point is close to the top boundary of block
             ! correct cell indices along non-AMR dimension
             iCell_II(iDim,nGridOut+1:2*nGridOut) = 1
             ! correct weights
             Weight_I(nGridOut+1:2*nGridOut) = Weight_I(1:nGridOut)*      Misc
             Weight_I(         1:  nGridOut) = Weight_I(1:nGridOut)*(1. - Misc)

             ! cycle through already found cells inside current block
             ! find their neighbors
             do iCell = 1, nGridOut
                ! location of the neighbor: current cells + displacement
                i_D = 1
                i_D(1:nDim) = iCell_II(1:nDim, iCell)
                XyzNei_D = XyzCorner_DI(:,iCell_II(0,iCell)) + &
                     ((/i_D(1), i_D(2),i_D(3)/) - 0.5) * &
                     Dxyz_DI(:,iCell_II(0,iCell))
                XyzNei_D(iDim) = XyzNei_D(iDim)+Dxyz_DI(iDim,iCell_II(0,iCell))

                call find_grid_block_ptr(XyzNei_D, iProcNei, iBlockNei)
                ! write correct processors and blocks
                iProc_I(    nGridOut+iCell) = iProcNei
                nBlockFound = nBlockFound + 1
                iBlock_I(nBlockFound) = iBlockNei
                iCell_II(0, nGridOut+iCell) = nBlockFound
             end do

          else
             ! the point is far enough from boundary

             ! correct cell indices along non-AMR direction
             iCell_II(iDim,nGridOut+1:2*nGridOut) =iCell_II(iDim,1:nGridOut) +1
             ! correct weights
             Weight_I(nGridOut+1:2*nGridOut) = Weight_I(1:nGridOut)*      Misc
             Weight_I(         1:  nGridOut) = Weight_I(1:nGridOut)*(1. - Misc)
             ! write correct processors
             iProc_I( nGridOut+1:2*nGridOut) = iProc_I(1:nGridOut) 
          end if

          ! each non-AMR direction doubles number of points in the stencil
          nGridOut = 2*nGridOut

       end do
    end if

    ! sort out cells located on other processors
    nCell = 0
    do iCell = 1, nGridOut
       if(iProc_I(iCell) == iProc)then
          nCell = nCell + 1
          iCell_II(       0, nCell) = iBlock_I(iCell_II(0,iCell))
          iCell_II(1:nDim, nCell) = iCell_II(1:nDim, iCell)
          Weight_I(nCell) = Weight_I(iCell)
       end if
    end do

  end subroutine interpolate_amr_wrapper

  !============================================================================

  subroutine find_block(nDimIn, XyzIn_D, &
       iProcOut, iBlockOut, XyzCornerOut_D, DxyzOut_D, IsOut)
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
    real,  intent(inout):: XyzIn_D(nDimIn)
    integer, intent(out):: iProcOut, iBlockOut !processor and block number
    !\
    ! Block left corner coordinates and the grid size:
    !
    real,    intent(out):: XyzCornerOut_D(nDimIn), DxyzOut_D(nDimIn)
    logical, intent(out):: IsOut !Point is out of the domain.

    real   :: Xyz_D(MaxDim) ! full Cartesian coords of point
    real   :: XyzCorner_D(MaxDim), Dxyz_D(MaxDim)! full Cartesian coords
    integer:: iDim, iDimAmr ! loop variable
    !--------------------------------------------------------------------
    ! check correctness
    if(nDimIn /= nDimAmr) &
         call CON_stop("Number of dimensions is not correct")

    ! Xyz_D will be used to find block, copy input data into it
    ! restore non-AMR directions if necessary
    Xyz_D = 0
    if(nDimIn == nDim)then
       Xyz_D(1:nDim) = XyzIn_D(1:nDim)
    else
       iDimAmr = nDimIn
       do iDim = nDim, 1, -1
          if(IsAmr_D(iDim))then
             Xyz_D(iDim) = XyzIn_D(iDimAmr)
             iDimAmr = iDimAmr - 1
          else
             Xyz_D(iDim) = XyzNonAmr_D(iDim)
          end if
       end do
    end if

    ! call internal BATL find subroutine
    call find_grid_block_ptr(Xyz_D, iProcOut, iBlockOut,&
         CoordMinBlockOut_D = XyzCorner_D, &
         CellSizeOut_D      = Dxyz_D)

    ! check if position has been found
    if(iProcOut==Unset_ .OR. iBlockOut==Unset_)then
       IsOut = .true.
       RETURN
    end if

    !position has been found
    IsOut = .false.

    ! corner coordinates of the found block
    XyzCornerOut_D = PACK(XyzCorner_D, IsAmr_D)

    ! cell size of the found block
    DxyzOut_D = PACK(Dxyz_D, IsAmr_D)

    ! subtract coordinates of the corner from point's coordinates
    XyzIn_D = XyzIn_D - XyzCornerOut_D

    ! store parameters of the block
    nBlockFound = nBlockFound + 1
    XyzCorner_DI(:, nBlockFound) = XyzCorner_D
    Dxyz_DI(     :, nBlockFound) = Dxyz_D
    iBlock_I(       nBlockFound) = iBlockOut
    iBlockOut                    = nBlockFound

  end subroutine find_block
  !============
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
    !--------------------------------------------------------------------
    CoordTreeSize_D = CoordMax_D - CoordMin_D
    ! Coord_D will be used to find block, copy input data into it
    ! restore non-AMR directions if necessary
    CoordFull_D = 0
    if(nDimIn == nDim)then
       CoordFull_D(1:nDim) = Coord_D(1:nDim)
    else
       iDimAmr = nDimIn
       do iDim = nDim, 1, -1
          if(IsAmr_D(iDim))then
             CoordFull_D(iDim) = Coord_D(iDimAmr)
             iDimAmr = iDimAmr - 1
          else
             CoordFull_D(iDim) = XyzNonAmr_D(iDim) !Should be renamed
          end if
       end do
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

    Dist_D = 0.5 + nIJK_D*(CoordFull_D - PositionMin_D)/&
         (PositionMax_D - PositionMin_D)
    !\
    ! In the following two global arrays at iDimNoAmr position
    ! there are the cell index and interpolation weight
    !/   
    iCellFloor_D = floor(Dist_D) 
    Dist_D = Dist_D - iCellFloor_D
    !\
    ! The stuff below is only needed to handle the non-amr direction?
    !/
    ! store parameters of the block
    nBlockFound = nBlockFound + 1
    XyzCorner_DI(:, nBlockFound) = CoordCornerFull_D
    Dxyz_DI(     :, nBlockFound) = DCoordFull_D
    iBlock_I(       nBlockFound) = iBlock
    iBlock                       = nBlockFound
    !\
    ! And how about iProc?
    !/
  end subroutine find
  !============
end module BATL_interpolate_amr_wrapper
