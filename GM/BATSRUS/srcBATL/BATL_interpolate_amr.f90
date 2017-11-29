module BATL_interpolate_amr

  use BATL_mpi,  ONLY: iProc
  use BATL_size, ONLY: nIJK_D, nDim, nDimAmr, MaxDim, &
       iRatio, jRatio, kRatio, iRatio_D, Dim2_
  use BATL_tree, ONLY: Unset_, find_tree_node, Proc_, Block_,&
       get_tree_position, iTree_IA
  use BATL_geometry, ONLY: CoordMin_D, DomainSize_D, IsPeriodic_D

  use ModInterpolateAMR, ONLY: cTol2, &
       interpolate_amr_shared             =>interpolate_amr, &
       interpolate_extended_stencil_shared=>interpolate_extended_stencil, &
       get_reference_block_shared         =>get_reference_block

  use ModUtilities, ONLY: CON_stop

  implicit none

  SAVE
  private ! except

  public:: interpolate_amr

  ! non-AMR direction:
  ! only 1 such direction, if 2 or more => interpolate_amr is not called
  ! (must be handled outside this module)
  ! MIN and MAX are added in order to keep value in range 1 to nDim
  integer, parameter:: iDimNoAmr = &
       max(1,min(1*(2 - iRatio) + 2*(2 - jRatio) + 3*(2 - kRatio), nDim))

  ! order of dimensions to correctly place AMR and non-AMR directions
  ! it depends on iDimNoAmr:
  ! iOrder_I = (/2, 3, 1/) for iDimNoAmr = 1
  ! iOrder_I = (/1, 3, 2/) for iDimNoAmr = 2
  ! iOrder_I = (/1, 2, 3/) for iDimNoAmr = 3
  ! iOrder_I = (/1, 2, 3/) for nDimAmr = nDim
  ! MIN and MAX are added in order to keep value in range 1 to nDim
  integer, parameter:: iOrder_I(MaxDim) = (/&
       3 - iRatio, &
       min(6 - iRatio - jRatio, nDim),  &
       max(-3 + 2*iRatio + jRatio, 1) /)

  ! point's coordinate in non-AMR dimensions
  real   :: CoordNoAmr, DisplacedCoordTreeNoAmr
  logical:: IsNearBlockBoundaryNoAmr

  !\
  ! variables to keep track of nodes in the case of non-AMR direction
  !-------------------------------------------------------------------------
  ! number of nodes found so far
  integer:: nNodeFound

  ! cell indexes along non-AMR direction:
  ! iCellNoAmr_I(1) - for the cell which center is left  to the given point
  ! iCellNoAmr_I(2) - for the cell which center is right to the given point
  integer   :: iCellNoAmr_I(2)

  ! interpolation weight along non-AMR direction
  real      :: WeightNoAmr
  logical   :: IsSecondOrderNoAmr

contains
  !============================================================================

  subroutine interpolate_amr(CoordIn_D, &
       nCell, iCell_II, Weight_I, IsSecondOrder)

    ! Find the grid cells surrounding the point Coord_D.
    ! nCell returns the number of cells found on the processor.
    ! iCell_II returns the block+cell indexes for each cell.
    ! Weight_I returns the interpolation weights calculated
    !                                 using AMR interpolation procedure
    ! IsSecondOrder returns whether the result is 2nd order interpolation

    real,    intent(in) :: CoordIn_D(MaxDim)
    integer, intent(out):: nCell
    integer, intent(out):: iCell_II(0:nDim,2**nDim)
    real,    intent(out):: Weight_I(2**nDim)
    logical, optional, intent(out):: IsSecondOrder

    integer:: iIndexes_II(0:nDimAmr+1,2**nDimAmr)
    integer:: nGridOut

    ! used in the case of non-AMR direction when point falls close to block's
    ! boundary and shared interpolation procedure is called 2nd time
    integer:: nGridOutAux

    integer:: iProc_I(2**nDim)
    integer:: iCell, iDim ! loop variables

    ! coords along non-AMR direction
    !--------------------------------------------------------------------------
    if(nDimAmr < nDim) then
       ! this case is valid only for nDim=MaxDim=3, nDimAmr=2
       CoordNoAmr = CoordIn_D(iDimNoAmr)
    end if

    ! mark the beginning of the interpolation
    nNodeFound = 0
    iProc_I    = Unset_

    ! apply general interpolate_amr routine but exclude non-AMR dimensions
    ! need to account for these dimensions separately,
    ! NOTE: routine doesn't sort out blocks located on other processors,
    !       this is done in the end of the current subroutine
    call interpolate_amr_shared(&
         nDim          = nDimAmr, &
         XyzIn_D       = CoordIn_D(iOrder_I(1:nDimAmr)), &
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

    !\
    ! Handle the non-AMR direction
    !/
    ! store indexes and restore shape along non-AMR direction
    iCell_II(  nDim,         1:  nGridOut) = iCellNoAmr_I(1)

    if(WeightNoAmr < cTol2)then
       ! this includes the case IsSecondOrderNoAmr = .false.:
       ! in this case WeightNoAmr is set to 0
       do iDim = 1, nDim
          iCell_II(iOrder_I(iDim),1:nGridOut) = iCell_II(iDim,1:nGridOut)
       end do
       call sort_out_other_procs
       ! return whether interpolation is of 2nd order
       if(present(IsSecondOrder))&
            IsSecondOrder = IsSecondOrder .and. IsSecondOrderNoAmr
       RETURN
    end if

    if(.not.IsNearBlockBoundaryNoAmr)then
       !\
       ! the other half of the stencil is in the same block:
       !/
       ! local node ids & cell indexes along AMR directions
       ! are the same for cells iCell and nGridOut+iCell
       iCell_II(0:nDimAmr,nGridOut+1:2*nGridOut)=iCell_II(0:nDimAmr,1:nGridOut)
       iCell_II(     nDim,nGridOut+1:2*nGridOut)=iCellNoAmr_I(2)
       ! copy processor indices
       iProc_I(nGridOut+1:2*nGridOut) = iProc_I(1:nGridOut)
       ! restore shape
       do iDim = 1, nDim
          iCell_II(iOrder_I(iDim),1:2*nGridOut) =iCell_II(iDim,1:2*nGridOut)
       end do
       ! correct weights
       Weight_I(nGridOut+1:2*nGridOut)= Weight_I(1:nGridOut) *    WeightNoAmr
       Weight_I(         1:  nGridOut)= Weight_I(1:nGridOut)*(1 - WeightNoAmr)
       ! number of cells has doubled
       nGridOut = 2*nGridOut
    else
       !\
       ! the point is far from domain's boundary =>
       ! the other part of the stencil exists and is in a different block(s)
       !/
       ! change coordinate along non-AMR direction to the displaced one:
       ! now restored full coordinates will fall in the appropriate block(s)
       CoordNoAmr = CoordMin_D(iDimNoAmr) + &
            DisplacedCoordTreeNoAmr * DomainSize_D(iDimNoAmr)

       ! call shared interpolation procedure 2nd time
       call interpolate_amr_shared(&
            nDim          = nDimAmr, &
            XyzIn_D       = CoordIn_D(iOrder_I(1:nDimAmr)), &
            nIndexes      = nDimAmr+1, &
            find          = find, &
            nCell_D       = nIJK_D(iOrder_I(1:nDimAmr)), &
            nGridOut      = nGridOutAux, &
            Weight_I      = Weight_I(nGridOut+1:nGridOut+2**nDimAmr), &
            iIndexes_II   = iIndexes_II, &
            IsSecondOrder = IsSecondOrder, &
            UseGhostCell  = .false.)
       if(nGridOutAux==0)&
            call CON_stop("Failure in BATL_interpolate_amr:interpolate_amr; Part of interpolation stencil can't be found")

       ! copy results indices
       iProc_I(            nGridOut+1:nGridOut+nGridOutAux) = &
            iIndexes_II(0,         1:nGridOutAux)
       iCell_II(1:nDimAmr, nGridOut+1:nGridOut+nGridOutAux) = &
            iIndexes_II(1:nDimAmr, 1:nGridOutAux)
       iCell_II(0,         nGridOut+1:nGridOut+nGridOutAux) = &
            iIndexes_II(nDimAmr+1, 1:nGridOutAux)

       ! store indexes and restore shape along non-AMR direction
       iCell_II(  nDim, nGridOut+1:nGridOut+nGridOutAux) = iCellNoAmr_I(2)
       do iDim = 1, nDim
          iCell_II(iOrder_I(iDim),1:nGridOut+nGridOutAux) = &
               iCell_II(iDim,1:nGridOut+nGridOutAux)
       end do

       ! correct weights
       Weight_I(         1:nGridOut) = Weight_I(1:nGridOut) * (1 - WeightNoAmr)
       Weight_I(nGridOut+1:nGridOut+nGridOutAux) = &
            Weight_I(nGridOut+1:nGridOut+nGridOutAux) * WeightNoAmr

       ! number of cells has changed
       nGridOut = nGridOut + nGridOutAux
    end if

    call sort_out_other_procs

  contains
    !==========================================================================
    subroutine sort_out_other_procs
      ! sort out cells located on other processors
      !------------------------------------------------------------------------
      nCell = 0
      do iCell = 1, nGridOut
         if(iProc_I(iCell) == iProc)then
            nCell = nCell + 1
            iCell_II(:, nCell) = iCell_II(:, iCell)
            Weight_I(nCell) = Weight_I(iCell)
         end if
      end do
    end subroutine sort_out_other_procs
    !==========================================================================

  end subroutine interpolate_amr
  !============================================================================
  subroutine fix_tree_coord(CoordTree_D)
    use BATL_geometry, ONLY: IsLatitudeAxis, IsSphericalAxis, &
         IsCylindricalAxis
    real,intent(inout) :: CoordTree_D(MaxDim)  ! Normalized gen coords
    !--------------------------------------------------------------------------
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
       ! spherical: r, lon, lat coordinates
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
  !============================================================================
  subroutine find(nDimIn, Coord_D, &
       iProc, iBlock, CoordCorner_D, dCoord_D, IsOut)
    integer, intent(in) :: nDimIn
    !\
    ! "In"- the coordinates of the point, "out" the coordinates of the
    ! point with respect to the block corner. In the most cases
    ! CoordOut_D = CoordIn_D - CoordCorner_D, the important distinction,
    ! however, is the periodic boundary, near which the jump in the
    ! stencil coordinates might occur. To handle the latter problem,
    ! we added the "out" intent. The coordinates for the stencil
    ! and input point are calculated and recalculated below with
    ! respect to the block corner.
    !/
    real,  intent(inout):: Coord_D(nDimIn)
    integer, intent(out):: iProc, iBlock ! processor and block number
    !\
    ! Block left corner coordinates and the grid size:
    !
    real,    intent(out):: CoordCorner_D(nDimIn), dCoord_D(nDimIn)
    logical, intent(out):: IsOut ! Point is out of the domain.

    real   :: CoordFull_D(MaxDim)       ! Full Gen coords of point
    real   :: CoordCornerFull_D(MaxDim) ! Full Gen coords of corner
    real   :: dCoordFull_D(MaxDim)      ! Full mesh sizes in gen coords
    real   :: CoordTree_D(MaxDim)       ! Normalized gen coords
    integer:: iNode                     ! tree node number
    real   :: PositionMin_D(MaxDim), PositionMax_D(MaxDim)
    !--------------------------------------------------------------------------
    ! Coord_D will be used to find block, copy input data into it
    ! restore non-AMR directions if necessary
    CoordFull_D = 0
    if(nDimAmr == nDim)then ! NOTE: nDimIn = nDimAmr
       CoordFull_D(1:nDim) = Coord_D(1:nDim)
    else
       CoordFull_D(iOrder_I(1:nDimAmr)) = Coord_D
       CoordFull_D(iDimNoAmr) = CoordNoAmr
    end if
    !\
    ! Calculate normalized (to DomainSize_D) coordinates for tree search
    !/
    CoordTree_D = (CoordFull_D - CoordMin_D)/DomainSize_D
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
    ! position has been found
    !/
    iBlock = iTree_IA(Block_,iNode)
    iProc  = iTree_IA(Proc_, iNode)
    call get_tree_position(iNode=iNode,                &
         PositionMin_D=PositionMin_D,&
         PositionMax_D=PositionMax_D )
    ! corner coordinates for the found block
    CoordCornerFull_D =  CoordMin_D + PositionMin_D*DomainSize_D
    CoordCorner_D = CoordCornerFull_D(iOrder_I(1:nDimAmr))

    ! cell size for the found block
    dCoordFull_D = (PositionMax_D - PositionMin_D)*DomainSize_D/nIjk_D
    dCoord_D = dCoordFull_D(iOrder_I(1:nDimAmr))

    ! subtract coordinates of the corner from point's coordinates
    Coord_D = (CoordTree_D(iOrder_I(1:nDimAmr)) - &
         PositionMin_D(iOrder_I(1:nDimAmr)))*&
         DomainSize_D(iOrder_I(1:nDimAmr))

    ! if there is no non-AMR direction => no need to do extra work
    if(nDimAmr == nDim) RETURN ! NOTE: nDimIn = nDimAmr

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

       ! set displaced tree coordinates to -1:
       !   if unchanged => point if far from block boundary along iDimNoAmr
       DisplacedCoordTreeNoAmr = -1.0
       IsSecondOrderNoAmr = .true.
       IsNearBlockBoundaryNoAmr = .false.
       ! take care of the cases when point is close to the boundary
       if    (iCellNoAmr_I(1) == 0)then
          ! displacement towards neighbor: down
          IsNearBlockBoundaryNoAmr = .true.
          DisplacedCoordTreeNoAmr = CoordTree_D(iDimNoAmr) - &
               dCoordFull_D(iDimNoAmr) / DomainSize_D(iDimNoAmr)
          if(IsPeriodic_D(iDimNoAmr))&
               DisplacedCoordTreeNoAmr = &
               modulo(DisplacedCoordTreeNoAmr, 1.0)
          if(DisplacedCoordTreeNoAmr < 0.0)then
             IsSecondOrderNoAmr = .false.
             WeightNoAmr = 0.0
             iCellNoAmr_I(1) = iCellNoAmr_I(2)
             RETURN
          end if
          iCellNoAmr_I(2) = nIJK_D(iDimNoAmr)
          iCellNoAmr_I(1) = 1
          ! inverse weight along iDimNoAmr
          WeightNoAmr = 1 - WeightNoAmr
       elseif(iCellNoAmr_I(2) == nIJK_D(iDimNoAmr)+1)then
          ! displacement towards neighbor: up
          IsNearBlockBoundaryNoAmr = .true.
          DisplacedCoordTreeNoAmr = CoordTree_D(iDimNoAmr) + &
               dCoordFull_D(iDimNoAmr) / DomainSize_D(iDimNoAmr)
          if(IsPeriodic_D(iDimNoAmr))&
               DisplacedCoordTreeNoAmr = &
               modulo(DisplacedCoordTreeNoAmr, 1.0)
          if(DisplacedCoordTreeNoAmr >= 1.0)then
             IsSecondOrderNoAmr = .false.
             WeightNoAmr = 0.0
             iCellNoAmr_I(2) = iCellNoAmr_I(1)
             RETURN
          end if
          iCellNoAmr_I(2) = 1
       end if
    end if

    nNodeFound = nNodeFound + 1
  end subroutine find
  !============================================================================

end module BATL_interpolate_amr
!==============================================================================
