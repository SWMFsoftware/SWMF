!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModInterpolateAMR
  !\
  !Generalize bilinear and trilinear interpolation for AMR grids
  !The data are given at the cell-centered grid which consists of AMR blocks.
  !The interpolation is free of any jumps at the resolution interfaces
  !including edges and corners
  !/
  !\
  !USE
  !/
  !=================ARRAYS FOR A CUBE====================================!
  !For a cubic stencil enumerated as follows:
  !
  !       ^
  ! z-axis|
  !       |  7----------8
  !       | /|         /|
  !       |/ |        / |
  !       5----------6  |
  !       |  |   _   |  |
  !       |  |   /|y-axi$
  !       |  |  /    |  |
  !       |  | /     |  |
  !       |  |/      |  |
  !       |  3----------4
  !       | /        | /
  !       |/         |/
  !       1----------2------> x-axis
  !       
  ! we provide several functions characterizing its geometry in terms of 
  ! the grid point numbers. 

  use ModCubeGeometry, ONLY: iFace_IDI
  !Vertexes (enumerated by the first undex), which form 
  !the face of direction iDir (second index) including the 
  !given vertex (the third index). 
  !When the first index equals 1,2,3,4, the vertex, accordingly: 
  !(1) coincides with the given one
  !(2) is connected with the given one by the edge of direction iDir+1
  !(3) is connected with the given one by the edge of direction iDir+2
  !(4) is connected to the given one by the face diagonal of direction iDir

  use ModCubeGeometry, ONLY: iOppositeFace_IDI
  !Vertexes (enumerated by the first undex), which form 
  !the face of direction iDir (second index) and does not 
  !include the given vertex (the third index). 
  !When the first index equals 1,2,3,4, the vertex, accordingly: 
  !(1) is connected to the given one by the edge of direction iDir
  !(2) is connected to (1) by the edge of direction iDir+1
  !(3) is connected to (1) by the edge of direction iDir+2
  !(4) is connected to  by the face diagonal of direction iDir

  use ModCubeGeometry, ONLY: iEdge_ID  
  !Number of the vertex connected by 
  !the edge of direction iDir (second index) 
  !with the given vertex (first index)

  !=================ARRAYS FOR A RECTANGLE==============!
  !For a rectangular stencil enumerated as follows:
  !
  !       ^
  ! y-axis|
  !       |   
  !       3----------4  
  !       |          |  
  !       |          !
  !       |          |  
  !       |          |  
  !       |          |  
  !       |          !
  !       |          |
  !       |          |
  !       1----------2------> x-axis
  !       
  ! we provide several functions characterizing its geometry in terms of 
  ! the grid point numbers. 

  use ModCubeGeometry, ONLY: iSide_IDI
  !Vertexes (enumerated by the first undex), which form 
  !the side of direction iDir (second index) including the 
  !given vertex (the third index). 
  !When the first index equals 1,2,3,4, the vertex, accordingly: 
  !(1) coincides with the given one
  !(2) is connected with the given one by the edge of direction iDir

  use ModCubeGeometry, ONLY: iOppositeSide_IDI
  !Vertexes (enumerated by the first undex), which form 
  !the side of direction iDir (second index) and does not include the
  !given vertex (the third index). 
  !When the first index equals 1,2 the vertex, accordingly: 
  !(1) is connected to the given one by the side of direction 1+mod(iDir,2)
  !(2) is connected to (1) by the edge of direction iDir
  !---------------------REFINEMENT CHARACTERSTICS------------
  use ModCubeGeometry, ONLY: Coarse_, Fine_

  !\
  ! Different cases of 3D stencil 
  !/
  use ModCubeGeometry, ONLY: Uniform_, Face_, Edge_, OneFine_, OneCoarse_,& 
        Transition2Edge_, Transition2Corner_, TransitionJunction_
  !\
  ! Analogous for 2D
  !/
  use ModCubeGeometry, ONLY: Trapezoid_, Rhombus_ 
  !\
  ! Array to sort stencil
  !/
  use ModCubeGeometry, ONLY: iSortStencil3_II, Case_, Dir_, Grid_
  !\
  !Array used to sort stencil: 
  !(1) assign a type to it (Case_) 
  !(2) assign a grid point of the stencil chosen as the basic one.
  !(Grid_) example, if a a single point in the stencil is coarse and all
  !other are fine, the 'basic point' is the coarse point.
  !(3) assign a direction, if the orientation if this sort of stencil
  !has a characteristic direction (Dir_). For example, for the 'resolution 
  !edge' of y direction iDir equals 1.
  !/
  implicit none
  PRIVATE !Except
  SAVE
  integer, parameter:: &
       x_       = 1,               &
       y_       = 2,               &
       z_       = 3
  !============================================================================
  !Interpolation on the block AMR grid
  !\
  !Calculates interpolation weights
  !/
  !\
  !Example of application for SERIAL calculation of the
  !interpolated value of the state vector sampled in the 
  !grid points as
  !State_VGB(nVar,nI,nJ,nK,nBlock) array, where nI=nJ=nK=4
  !in point Xyz_D looks as follows:
  !
  !call interpolate_amr(&
  !  nDim=3,              &!number of dimensions
  !  XyzIn_D=,Xyz_D,      &!Point in which to interpolate
  !  nIndexes = 4,        &!Three cell indexes plus one block index
  !  nCell_D  = (/4,4,4/) &!
  !  find = find_subroutine , &! Search in grid
  !  nGridOut=nGridOut,   &! Number of points in the output stencil
  !  Weight_I=Weight_I,   &! Weight coefficients
  !  iIndexes_II=iIndexes_II) !Cell+block indexes to be used in interpolation
  !
  !  if(nGridOut < 1) call CON_stop('Interpolation failed')
  !  Value_V(1:nVar) = 0
  !  do iGrid = 1, nGridOut
  !     Value_V = Value_V + &
  !       State_VGB(:,iIndexes_II(1,iGrid), iIndexes_II(2,iGrid), &
  !             iIndexes_II(3,iGrid), iIndexes_II(4,iGrid))*&
  !                               Weight_I(iGrid)
  !  end do
  ! For PARALLEL code the processor number at which the corresponding part
  ! of  State_VGB is allocated is provided in iIndex_II(0,:) components
  ! of the index output array
  !/
  public interpolate_amr 
contains
  !====================================================================
  subroutine interpolate_amr(nDim, XyzIn_D, nIndexes, nCell_D, find, &
       nGridOut, Weight_I, iIndexes_II, IsSecondOrder)
    use ModInterpolateBody, ONLY: resolution_corner
    !\
    !USE: The tools to fill in sort stencil array
    !/
    use ModCubeGeometry, ONLY: DoInit, init_sort_stencil
    use ModKind, ONLY: nByteReal
    !\
    ! INPUT PARAMETERS
    !/
    !\
    ! Number of dimensions
    !/
    integer, intent(in) :: nDim 
    !\
    ! Number of indexes. Usually, nIndexes = nDim + 1, for three cell indexes 
    ! and a block number. If the latter information is not needed, 
    ! nIndexes=nDim should be added
    !/
    integer, intent(in) :: nIndexes
    !\
    ! Point coordinates
    !/
    real,    intent(in) :: XyzIn_D(nDim)
    !\
    ! Block AMR Grid characteristic: number of cells
    !/ 
    integer, intent(in) :: nCell_D(nDim)
    !\
    ! Yet another Block AMR Grid characteristic:
    ! the search routine, which returns, for a given point,
    ! the block to which this points belong and the processor, at which
    ! this block is allocated, as well as the block parameters:
    ! coordinates of the left corner (the point in the block with the 
    ! minvalue of coordinates) and (/Dx, Dy, Dz/)
    interface 
       subroutine find(nDim, Xyz_D, &
            iProc, iBlock, XyzCorner_D, Dxyz_D, IsOut)
         implicit none
         integer, intent(in) :: nDim
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
         real,  intent(inout):: Xyz_D(nDim)
         integer, intent(out):: iProc, iBlock !processor and block number
         !\
         ! Block left corner coordinates and the grid size:
         !/
         real,    intent(out):: XyzCorner_D(nDim), Dxyz_D(nDim)
         logical, intent(out):: IsOut !Point is out of the domain.
       end subroutine find
    end interface
    !\
    !OUTPUT PARAMETERS
    !/
    !\
    !Number of grid points involved into interpolation
    !/
    integer, intent(out):: nGridOut      
    !\
    ! Interpolation weights (only the first nGridOut values are meaningful
    !/
    real,    intent(out):: Weight_I(2**nDim) 
    !\
    ! Cell(+block) indexes and processor number for grid points to be
    ! invilved into interpolation. iProc numbers are stored in 
    ! iIndexes_II(0,:) 
    !/
    integer, intent(out):: iIndexes_II(0:nIndexes,2**nDim)
    !\
    !The following is true if stencil does not employ 
    !the out-of-grid points
    !/
    logical, intent(out), optional:: IsSecondOrder  
    !\
    ! Local variables
    !/
    !\
    ! The output array in generate_basic_stencil routine
    ! This is the set of numbers of the elements of the extended stencil
    ! to be included into the basic stencil
    !/
    integer, dimension(2**nDim):: iOrderExtended_I
    !\
    ! The output array in interpolate_amr2,3  routine
    ! This is the set of numbers of the elements of the basic stencil
    ! to be involved in interpolation
    !/
    integer, dimension(2**nDim):: iOrder_I
    !\
    ! Output parameter of interpolate_amr3
    ! If .true., the basic stencil should be re-evaluated
    !/
    logical                    :: IsCorner
    !\
    ! Basic stencil; points:
    !/
    real                       :: XyzGrid_DI(nDim,2**nDim)
    !\
    ! Basic stencil; refinement level in the grid points:
    !/
    integer                    :: iLevel_I(2**nDim)
    !\
    ! Coordinates of the input point may be recalculated within
    ! routine; inverse of (/Dx, Dy, Dz/) is reused
    !/
    real, dimension(nDim)      :: Xyz_D, DxyzInv_D, XyzStencil_D

    !\
    ! Extended stencil, the generous estimate for its size is 2**(2*nDim) 
    !/
    real, dimension(nDim,64) :: XyzExtended_DI 
    integer :: iLevelExtended_I(64)
    integer :: iIndexesExtended_II(0:nIndexes,64)
    integer :: nExtendedStencil
    logical :: IsOutExtended_I(64)

    !\
    ! The same extended stencil in a structured form:
    ! A cubic 2*2*2 grid with 2*2*2 subgrids covering each vertex
    !/
    real       :: XyzGrid_DII(nDim,0:2**nDim,2**nDim)
    integer, dimension(2**nDim):: iBlock_I  , iProc_I 
    integer                    :: iCellIndexes_DII(nDim,2**nDim,2**nDim)
    integer, dimension(2**nDim):: nSubGrid_I, iLevelSubgrid_I
    logical, dimension(2**nDim):: IsOut_I

    real, dimension(nDim) :: XyzBasicBlock_D, DxyzBasicBlock_D
    logical               :: IsOutOfDomain
    !\
    ! The sort of stencil as derived from EXTENDED stencil
    !/
    integer:: iCaseExtended

    !\
    ! UseGhostPoint = .true. if in the extended stencil there are
    ! out-of-grid points 
    logical:: UseGhostPoint, IsExtended
    !/
    !\
    !Just 2**nDim
    integer:: nGrid
    !/
    !\
    ! To improve the algorithm stability against roundoff errors
    !/
    real, parameter:: cTol = 0.00000010
    real           :: cTol2
    !\
    ! Shift of the iGrid point in the stencil with respect to the
    ! first one
    !/
    integer, dimension(3,8),parameter :: iShift_DI = reshape((/&
         0, 0, 0,   1, 0, 0,   0, 1, 0,   1, 1, 0, &
         0, 0, 1,   1, 0, 1,   0, 1, 1,   1, 1, 1/),(/3,8/)) 

    integer:: iGridOutOfBlock
    !------------------------
    cTol2 = cTol**(nByteReal/4)
    IsExtended = .false.; UseGhostPoint = .false.

    nGrid = 2**nDim !Number of points in a basic stencil

    !\
    ! Initialize 
    !/
    iOrder_I    = 0; iIndexes_II = 0; Weight_I    = 0
    nGridOut = -1; Xyz_D = XyzIn_D ; IsOut_I = .false.
    if(present(IsSecondOrder))IsSecondOrder = .false.
    !\
    ! Find block to which the point belong
    !/ 
    call find(nDim, Xyz_D, iProc_I(1), iBlock_I(1), &
         XyzBasicBlock_D, DxyzBasicBlock_D, IsOutOfDomain)
    if(IsOutOfDomain)then
       !\
       ! The algorithm does not work for a point out of the computation
       ! domain. It could, but in this case too much information
       ! about grid should be brought to the table - the domain size,
       ! periodicity etc
       !/
       nGridOut = -1
       RETURN
    end if
    !\
    ! Now Xyz_D is given  with respect to the main block corner
    !/
    call get_main_block(iGridOutOfBlock,IsFirstCall = .true.)
    !\
    !The interpolation is done, if the stencil is within a single block
    !/
    if(nGridOut > 0) then
       if(present(IsSecondOrder))IsSecondOrder = .true.
       RETURN
    end if
    call get_other_blocks(iGridOutOfBlock)
    !\
    ! Prolong behind the boundary, if needed
    !/
    if(UseGhostPoint)then
       call prolong_beyond_boundary
       !\
       !The interpolation may be done inside generate_extended_stencil 
       !/
       if(nGridOut > 0) then
          if(present(IsSecondOrder))IsSecondOrder = .false.
          RETURN
       end if
    end if
    call generate_extended_stencil
    !\
    !The interpolation may be done inside generate_extended_stencil 
    !/
    if(nGridOut > 0) then
       if(present(IsSecondOrder))IsSecondOrder = .not. UseGhostPoint
       RETURN
    end if
    if(DoInit)call init_sort_stencil
    select case(nDim)
    case(2)
       call fix_basic_stencil2(&
            Xyz_D, DxyzInv_D, XyzGrid_DII(:,0,:), iLevelSubgrid_I, &
            XyzStencil_D, iCaseExtended)
       call generate_basic_stencil(&
            nDim, XyzStencil_D, nExtendedStencil,                      &
            XyzExtended_DI(:,1:nExtendedStencil), DxyzInv_D, iOrderExtended_I)
       XyzGrid_DI = XyzExtended_DI(:,iOrderExtended_I)
       iLevel_I = iLevelExtended_I(iOrderExtended_I)
       IsOut_I = IsOutExtended_I(iOrderExtended_I)
       iIndexes_II = &
            iIndexesExtended_II(:,iOrderExtended_I)
       call interpolate_amr2(&
            Xyz_D , XyzGrid_DI, iLevel_I, IsOut_I, iCaseExtended, &
            nGridOut, Weight_I, iOrder_I)
       if(nGridOut < 1)&
            call CON_stop('Failure in interpolate amr grid2') 
    case(3)
       !\
       ! For edges and corners we need to find, if point Xyz is
       ! really within the corner, or it fals into some transition 
       ! region (from edge to corner, from resolution interface to
       ! edge. iCaseExtended for these cases takes the values 
       ! Transition2Edge_, Transition2Corner_, TransitionJunction_
       !/ 
       call fix_basic_stencil3(&
            Xyz_D, DxyzInv_D, XyzGrid_DII(:,0,:), iLevelSubgrid_I, &
            XyzStencil_D, iCaseExtended, IsCorner)
       if(.not.IsCorner)then
          call generate_basic_stencil(&
               nDim, XyzStencil_D, nExtendedStencil,              &
               XyzExtended_DI(:,1:nExtendedStencil),              &
               DxyzInv_D, iOrderExtended_I)
          XyzGrid_DI = XyzExtended_DI(:,iOrderExtended_I)
          iLevel_I = iLevelExtended_I(iOrderExtended_I)
          iIndexes_II = &
               iIndexesExtended_II(:,iOrderExtended_I)
          IsOut_I = IsOutExtended_I(iOrderExtended_I)
          call interpolate_amr3(&
               Xyz_D , XyzGrid_DI, iLevel_I, IsOut_I, iCaseExtended,&
               nGridOut, Weight_I, iOrder_I, IsCorner)
       end if
       if(IsCorner)then
          !\
          ! This is not necessarily esleif,
          ! in case of sophisticated transition from edge to corner 
          ! amr3 routins finally decides if the point should be 
          ! interpolated with the corner stencil (if the interpolation 
          ! algorithm for tramsition region fails)
          !/
          call generate_corner_stencil
          call resolution_corner(Xyz_D , XyzGrid_DI, iLevel_I,&
               nGridOut, Weight_I, iOrder_I)
       end if
    end select
    !\
    ! Eliminate repeating grid points and points with zero weight
    ! which may be behind the domain boundary
    !/
    call sort_out

    iIndexes_II(:, 1:nGridOut) = iIndexes_II(:,iOrder_I(1:nGridOut))

    if(present(IsSecondOrder))IsSecondOrder = .not.any(IsOut_I)
  contains
    subroutine get_main_block(iGridOutOfBlock,IsFirstCall)
      integer, intent(out):: iGridOutOfBlock
      logical,optional,intent(in):: IsFirstCall
      !\
      ! Fills in XyzGrid_DII(:,0,:) - coarse grid
      ! Fills in indexes for the points of the stencil 
      ! belonging to this block. Returns the maximum
      ! number of the grid poit which is out of the main
      ! block
      !/
      !\
      ! Loop variable
      !/
      integer:: iGrid
      !\
      !Displacement measured in grid sizes or in their halfs
      !/ 
      integer, dimension(nDim) :: iShift_D
      !\
      ! Misc
      !/
      real :: XyzMisc_D(nDim)
      !------------------------------------
      iLevelSubgrid_I =  0; nSubGrid_I      = 1
      iCellIndexes_DII = 0; XyzGrid_DII     = 0
      DxyzInv_D = 1/DxyzBasicBlock_D
      XyzMisc_D = Xyz_D*DxyzInv_D + 0.50
      !
      !       y ^
      !         |    i,j+1         i+1,j+1
      !         |     +              +
      !         |          o (x,y)           (x,y,z)= [(x,y,z)-
      !         |     +<------dx---->+        -(x_ijk,y_ijk,z_ijk)]
      !         |  x_ij=(i-0.5)dx   i+1,j    +(dx,dy,dz)(i,j,k)-0.5*
      !         |  y_ij=(j-0.5)dy             (dx,dy,dz)
      !         |----------------------------->
      !       block           x
      !       corner, xyz=0

      iCellIndexes_DII(:,1,1) = floor(XyzMisc_D)
      !\
      !Calculate coordinates of the left corner of a stencil
      !/
      XyzMisc_D = XyzMisc_D - iCellIndexes_DII(:,1,1)
      if(present(IsFirstCall))then

         if(all(iCellIndexes_DII(:,1,1) > 0).and.&
              all(  iCellIndexes_DII(:,1,1) < nCell_D))then
            !\
            ! The whole interpolation stencil is within this block
            !/
            call interpolate_uniform( XyzMisc_D)
            !\
            ! Form index array and sort out zero weights
            !/
            nGridOut = 0 ; cTol2 = 2*cTol2
            do iGrid = 1, nGrid
               if(Weight_I(iGrid) < cTol2)CYCLE
               nGridOut = nGridOut + 1
               iIndexes_II(0,       nGridOut) = iProc_I(1)
               iIndexes_II(nIndexes,nGridOut) = iBlock_I(1)
               iIndexes_II(1:nDim,  nGridOut) = iCellIndexes_DII(:,1,1) + &
                    iShift_DI(1:nDim,iGrid)
               Weight_I(nGridOut) = Weight_I(iGrid)
            end do
            RETURN  !All interpolation is done, ready to exit
            ! elseif(all(abs(XyzMisc_D) < cTol2))then
            !    !\
            !    ! Xyz coincides with the grid point
            !    ! Commented out to satisfy iFort
            !    !/
            !    nGridOut = 1; Weight_I =0; Weight_I(1) = 1
            !    iIndexes_II(0,       1) = iProc_I(1)
            !    iIndexes_II(nIndexes,1) = iBlock_I(1)
            !    iIndexes_II(1:nDim,  1) = iCellIndexes_DII(:,1,1)
            !    RETURN
         end if
      end if
      XyzGrid_DII(:,0,1) = DxyzBasicBlock_D*(iCellIndexes_DII(:,1,1) - 0.50)
      !\ 
      ! now XyzMisc_D = (Xyz_D-XyzGrid_DII(:,0,1))/Dxyz satisfies 
      ! inequalities: XyzMisc_D >= 0 and XyzMisc_D < 1. Strengthen 
      ! these inequalities 
      !/
      XyzMisc_D = min(1 - cTol2,&
           max(XyzMisc_D, cTol2 ))
      Xyz_D = XyzGrid_DII(:,0,1) + XyzMisc_D*DxyzBasicBlock_D
      !\
      !Calculate other grid points, check if all points belong to 
      !the found block
      !/
      iGridOutOfBlock = -1
      do iGrid = nGrid, 1, -1
         iShift_D = iShift_DI(1:nDim,iGrid)
         iBlock_I(iGrid) = iBlock_I(1)
         iCellIndexes_DII(:,1,iGrid) = &
              iCellIndexes_DII(:,1,1) + iShift_D
         XyzGrid_DII(:,0,iGrid) = &
              XyzGrid_DII(:,0,1) + iShift_D*DxyzBasicBlock_D
         if(any(iCellIndexes_DII(:,1,iGrid) < 1).or.&
              any(iCellIndexes_DII(:,1,iGrid) > nCell_D))then
            !\
            !This grid point is out of block, mark it
            !/
            iProc_I(iGrid) = -1 
            iGridOutOfBlock = max(iGridOutOfBlock,iGrid)
         else
            XyzGrid_DII(:,1,iGrid) = XyzGrid_DII(:,0,iGrid)
            nSubGrid_I(iGrid) = 1
            iProc_I(iGrid) = iProc_I(1)
         end if
      end do
    end subroutine get_main_block
    !=====================
    subroutine interpolate_uniform(Dimless_D)
      !\
      !The displacement of the point from the stencil left corner
      !normalized by the grid size. Not exceeding zero and is less than 1.
      !/
      real, dimension(nDim), intent(in)::Dimless_D
      !\
      !Loop variables
      !/
      !------------------
      Weight_I(1) = (1 - Dimless_D(1))*(1 - Dimless_D(2))
      Weight_I(2) =      Dimless_D(1) *(1 - Dimless_D(2))
      Weight_I(3) = (1 - Dimless_D(1))*     Dimless_D(2)
      Weight_I(4) =      Dimless_D(1) *     Dimless_D(2)
      if(nDim==3)then
         Weight_I(5:8) = Weight_I(1:4)*     Dimless_D(3)
         Weight_I(1:4) = Weight_I(1:4)*(1 - Dimless_D(3))
      end if
    end subroutine interpolate_uniform
    !=====================
    recursive subroutine get_other_blocks(iGridOutOfBlock)
      integer, intent(inout)::iGridOutOfBlock
      !\
      ! Loop variable
      !/
      integer:: iGrid 

      integer:: iGridStored
      !\
      ! For using find routine with inout argument
      !/
      real  :: XyzMisc_D(nDim)
      !\
      ! Output parameters of find routine
      !/
      real    :: XyzCorner_D(nDim), Dxyz_D(nDim)
      logical :: IsOut
      !----------------------------
      if(iGridOutOfBlock == -1)RETURN
      iGridStored = iGridOutOfBlock
      !\
      ! For the grid point not belonging to the block
      ! find the block they belong to
      !/ 

      !\
      !Recalculate absolute coordinates for
      !the grid point which is out of the block
      !/
      XyzMisc_D = XyzBasicBlock_D + XyzGrid_DII(:,0,iGridStored)
      !\
      ! Find neighboring block
      !/
      call find(nDim, XyzMisc_D, &
           iProc_I(iGridStored), iBlock_I(iGridStored), &
           XyzCorner_D, Dxyz_D, IsOut)
      if(IsOut)then
         UseGhostPoint = .true.
         iProc_I(iGridStored) = 0 !For not processing this point again
         IsOut_I(iGridStored) = .true.
         XyzGrid_DII(:,1,iGridStored) = XyzGrid_DII(:,0,iGridStored)
         !\
         ! Find the next out-of-block point
         !/
         do iGrid = iGridStored - 1, 1, -1
            if(iProc_I(iGrid)==-1)then
               iGridOutOfBlock = iGrid
               call get_other_blocks(iGridOutOfBlock)
               RETURN
            end if
         end do
         RETURN
      end if
      iLevelSubgrid_I(iGridStored) = 1 - &
           floor(Dxyz_D(1)*DXyzInv_D(1)+ cTol)
      !\                     ^
      ! For expression above | equal to 2 , 1, 0.5 correspondingly
      ! iLevel = -1, 0, Fine_, meaning that the neighboring block
      ! is coarser, at the same resolution or finer than the basic 
      ! one.
      !/
      select case(iLevelSubgrid_I(iGridStored))
      case(-1)
         !The neighboring block is coarser and should be used as the 
         !stencil base. Now Xyz_D and XyzGrid_DII(:,0,iGridStored) are
         !defined with respect to the original block, the latter point
         !in the coarser neighboring block has a coordinates XyzMisc_D
         !So that
         Xyz_D = Xyz_D - XyzGrid_DII(:,0,iGridStored) + XyzMisc_D
         !\
         ! Return to the beginning of algorithm. Before the COARSEN
         ! loop the block and processor number for the basic block
         ! were stored in iProc_I(1), iBlock_I(1)
         !/
         iProc_I( 1) = iProc_I(iGridStored)
         iBlock_I(1) = iBlock_I(iGridStored)
         XyzBasicBlock_D = XyzCorner_D
         DxyzBasicBlock_D = Dxyz_D
         call get_main_block(iGridOutOfBlock)
      case(0)  ! (New Dxyz_D)*Stored DXyzInv =1
         call get_block(iGridOutOfBlock, XyzMisc_D, Dxyz_D)
      case(Fine_  )  !1, (New Dxyz_D)*Stored DXyzInv =0.5
         IsExtended = .true. 
         call get_fine_block(iGridOutOfBlock, XyzMisc_D, Dxyz_D)
      end select
      call get_other_blocks(iGridOutOfBlock)
    end subroutine get_other_blocks
    !========================
    subroutine get_block(iGridOutOfBlock, Xyz_D, Dxyz_D)
      integer, intent(inout)::iGridOutOfBlock
      real, dimension(nDim), intent(in):: Xyz_D, Dxyz_D
      !\
      ! Fills in the indexes for the grid boints belonging
      ! to the block, which it at the same resolution as the
      ! main block Returns the maximum number of the grid poit
      ! which is out of all blocks  found so far
      !/
      !\
      ! Loop variable
      !/
      integer:: iGrid
      !\
      !Displacement measured in grid sizes or in their halfs
      !/ 
      integer, dimension(nDim) :: iShift_D
      !/
      integer :: iGridStored
      !---------------------
      iGridStored = iGridOutOfBlock
      XyzGrid_DII(:,1,iGridStored) = XyzGrid_DII(:,0,iGridStored)
      !\
      ! Calculate cell indexes as we did before. Use nint
      ! instead of int as long as XyzMisc_D is very close to 
      ! the grid point 
      !/
      iCellIndexes_DII(:,1,iGridStored) = &
           nint(Xyz_D*DxyzInv_D + 0.50)
      !\
      ! A single point in the subgrid
      !/
      nSubGrid_I(iGridStored)      = 1
      !\
      ! Check if there are more grid points belonging to the
      ! newly found block
      !/
      !\
      ! Store shift of iGridStored point
      !/
      iShift_D = iShift_DI(1:nDim,iGridStored)
      iGridOutOfBlock = -1
      do iGrid = iGridStored - 1, 1, -1
         if(iProc_I(iGrid)/=-1)CYCLE !This point is done earlier
         iCellIndexes_DII(:,1,iGrid) = &
              iCellIndexes_DII(:,1,iGridStored) +&
              iShift_DI(1:nDim,iGrid) - iShift_D 
         !\
         ! Check if the point is in the newly found block
         !/ 
         if(any(iCellIndexes_DII(:,1,iGrid) < 1).or.&
              any(iCellIndexes_DII(:,1,iGrid) > nCell_D))then
            !\
            !This grid point is out of block, mark it for further work
            !/
            iGridOutOfBlock = max(iGridOutOfBlock, iGrid)
         else
            iProc_I( iGrid) = iProc_I( iGridStored)
            iBlock_I(iGrid) = iBlock_I(iGridStored)
            XyzGrid_DII(:,1,iGrid) = XyzGrid_DII(:,0,iGrid)
            nSubGrid_I(     iGrid) = 1
            iLevelSubGrid_I(iGrid) = 0
         end if
      end do
    end subroutine get_block
    !=====================
    subroutine get_fine_block(iGridOutOfBlock, Xyz_D, Dxyz_D)
      integer, intent(inout)::iGridOutOfBlock
      real, dimension(nDim), intent(in):: Xyz_D, Dxyz_D
      !\
      ! Fills in the indexes for the grid boints belonging
      ! to the block, which it at the same resolution as the
      ! main block Returns the maximum number of the grid poit
      ! which is out of all blocks  found so far
      !/
      !\
      ! Loop variables
      !/
      integer:: iGrid, iSubGrid
      !\
      !Displacement measured in grid sizes or in their halfs
      !/ 
      integer, dimension(nDim) :: iShift_D
      !/
      integer :: iGridStored
      iGridStored = iGridOutOfBlock
      !\
      ! Fine subgrid is displaced by half of finer subgrid size
      !/
      XyzGrid_DII(:,1,iGridStored) = XyzGrid_DII(:,0,iGridStored) -&
           0.50*Dxyz_D
      !\
      ! Calculate cell indexes as we did above. Note that DxyzInv_D
      ! is twice less than needed, because it is calculated for the
      ! whole stencil, not for the finer subgrid
      !/
      iCellIndexes_DII(:,1,iGridStored) = &
           floor(2*Xyz_D*DxyzInv_D + 0.50)
      !\
      ! All points in the 2*2*2 finer subgrid are involved
      !/
      nSubGrid_I(iGridStored)      = nGrid
      do iSubGrid = 2, nGrid
         iShift_D = iShift_DI(1:nDim,iSubGrid)
         XyzGrid_DII(:,iSubGrid,iGridStored) = &
              XyzGrid_DII(:,1,iGridStored) + &
              Dxyz_D*iShift_D
         iCellIndexes_DII(:,iSubGrid,iGridStored) = &
              iCellIndexes_DII(:,1,iGridStored) + &
              iShift_D
      end do
      !\
      ! Check if there are more grid points belonging to the
      ! newly found block
      !/
      iGridOutOfBlock = -1
      do iGrid = iGridStored -1, 1, -1
         if(iProc_I(iGrid)/=-1)CYCLE !This point is done earlier
         iCellIndexes_DII(:,1,iGrid) = &
              iCellIndexes_DII(:,1,iGridStored) + 2*(&
              iShift_DI(1:nDim,iGrid) - iShift_DI(1:nDim,iGridStored))
         if(any(iCellIndexes_DII(:,1,iGrid) < 1).or.&
              any(iCellIndexes_DII(:,1,iGrid) > nCell_D))then
            !\
            !This grid point is out of block
            !/
            iGridOutOfBlock = max(iGridOutOfBlock, iGrid)
         else
            iProc_I(iGrid ) = iProc_I( iGridStored)
            iBlock_I(iGrid) = iBlock_I(iGridStored)
            XyzGrid_DII(:,1,iGrid) = XyzGrid_DII(:,0,iGrid) -&
                 0.50*Dxyz_D
            nSubGrid_I(iGrid)      = nGrid
            iLevelSubGrid_I(iGrid) = Fine_
            do iSubGrid = 2, nGrid
               iShift_D = iShift_DI(1:nDim,iSubGrid)
               XyzGrid_DII(:,iSubGrid,iGrid) = &
                    XyzGrid_DII(:,1,iGrid) + &
                    Dxyz_D*iShift_D
               iCellIndexes_DII(:,iSubGrid,iGrid) = &
                    iCellIndexes_DII(:,1,iGrid) + &
                    iShift_D
            end do
         end if
      end do
    end subroutine get_fine_block
    !=================
    subroutine prolong_beyond_boundary
      !\
      ! Handle points behind the boundary
      !/
      !\
      ! Loop variables
      !/
      integer:: iGrid, iSubgrid, iDir
      integer:: iLoc
      !---------------- 

      select case(count(IsOut_I))
      case(3, 7) !3,7 are nGrid -1 for nDim=2,3 
         !\
         !Find the only physical point in the stencil 
         !/
         iLoc = maxloc(iLevelSubGrid_I,MASK=.not.IsOut_I, DIM=1)
         nGridOut = 1
         Weight_I = 0; Weight_I(1) = 1
         iIndexes_II(0,       1) = iProc_I( iLoc) 
         iIndexes_II(nIndexes,1) = iBlock_I(iLoc)
         iIndexes_II(1:nDim,1  ) = iCellIndexes_DII(:,1,iLoc)
         RETURN
      case(4)
         ! nDim = 3
         !\
         ! One of the faces is fully out of the domain
         !/
         !\
         !Find point in the domain
         !/
         iLoc = maxloc(iLevelSubGrid_I,MASK=.not.IsOut_I,DIM=1)
         do iDir = 1, nDim
            if(all(&
                 IsOut_I(iOppositeFace_IDI(:,iDir,iLoc))))then
               if(IsExtended)then
                  !\
                  ! Prolong grid from the physical face to the ghost one
                  ! accounting for the difference in resolution
                  !/
                  do iGrid = 1,4
                     call prolong(iFace_IDI(iGrid,iDir,iLoc), &
                          iOppositeFace_IDI(iGrid,iDir,iLoc))
                  end do
               else
                  !\
                  !Grid near the boundary is uniform. Put the point to the
                  !physical face to nullify the ghost cell contributions
                  !/
                  Xyz_D(iDir) = XyzGrid_DII(iDir,1,iLoc)
               end if
               EXIT
            end if
         end do
      case(2)
         ! Analogous to the previous case, but nDim = 2 
         !\
         !Find point in the domain
         !/
         iLoc = maxloc(iLevelSubGrid_I,MASK=.not.IsOut_I, DIM=1)
         do iDir = 1, 2
            if(all(&
                 IsOut_I(iOppositeSide_IDI(:,iDir,iLoc))))then
               if(IsExtended)then
                  !\
                  ! Prolong grid from the physical face to the ghost 
                  ! accounting for the difference in resolution
                  !/
                  do iGrid = 1,2
                     call prolong(iSide_IDI(iGrid,iDir,iLoc), &
                          iOppositeSide_IDI(iGrid,iDir,iLoc))
                  end do
               else
                  !\
                  !Put the point to the physical face, which assign them
                  !zero weight
                  !/
                  Xyz_D(3 - iDir) = XyzGrid_DII(3 - iDir,1,iLoc)
               end if
               EXIT
            end if
         end do
      case(6)
         !\
         ! Three-dimensional case, only two grid vertexes are
         ! inside the domain. From each of the two physical
         ! points the grid should be prolonged to three
         ! ghost points along the plane
         !/
         iLoc = maxloc(iLevelSubGrid_I,MASK=.not.IsOut_I, DIM=1)
         do iDir  = 1, nDim
            if(.not.IsOut_I(iEdge_ID(iLoc,iDir)))then
               if(IsExtended)then
                  do iGrid = 2,4
                     call prolong(iLoc,&
                          iFace_IDI(iGrid,iDir,iLoc))
                     call prolong(iEdge_ID(iLoc,iDir),&
                          iOppositeFace_IDI(iGrid,iDir,iLoc))
                  end do
                  EXIT
               else
                  !\
                  ! Put the point to the physical edge
                  !/
                  Xyz_D(1 + mod(iDir,3)) = &
                       XyzGrid_DII(1 + mod(iDir,3),1,iLoc)
                  Xyz_D(1 + mod(iDir+1,3)) = &
                       XyzGrid_DII(1 + mod(iDir+1,3),1,iLoc)
               end if
            end if
         end do
      end select
    end subroutine prolong_beyond_boundary
    !======================
    !\
    ! Used to prolong grid behind the boundary
    !/
    subroutine prolong(iGridPhys, iGridGhost)          
      integer, intent(in) :: iGridPhys, iGridGhost
      !Loop variable
      integer :: iSubGrid
      !--------------------
      nSubGrid_I(iGridGhost) = nSubGrid_I(iGridPhys)
      iLevelSubgrid_I(iGridGhost) = iLevelSubgrid_I(iGridPhys) 
      do iSubGrid = 1, nSubGrid_I(iGridGhost)
         XyzGrid_DII(:,iSubGrid,iGridGhost) = &
              XyzGrid_DII(:,iSubGrid,iGridPhys) +&
              XyzGrid_DII(:,0,iGridGhost)       -&
              XyzGrid_DII(:,0,iGridPhys )
      end do
    end subroutine prolong
    !======================
    subroutine generate_extended_stencil
      integer:: iGrid, iSubGrid
      if(IsExtended)then
         !\ 
         ! calculate all arrays for extended stencil
         !/
         nExtendedStencil    = 0
         XyzExtended_DI      = 0
         iIndexesExtended_II = 0
         iLevelExtended_I    = 0
         IsOutExtended_I = .false.
         do iGrid = 1, nGrid
            do iSubgrid = 1, nSubgrid_I(iGrid)
               nExtendedStencil = nExtendedStencil + 1
               XyzExtended_DI(:,nExtendedStencil) = &
                    XyzGrid_DII(:,iSubGrid,iGrid) 
               iIndexesExtended_II(0,nExtendedStencil) = &
                    iProc_I(iGrid)
               iIndexesExtended_II(nIndexes,nExtendedStencil) = &
                    iBlock_I(iGrid)
               iIndexesExtended_II(1:nDim,nExtendedStencil)   = &
                    iCellIndexes_DII(:,iSubGrid,iGrid)
               iLevelExtended_I(nExtendedStencil) = iLevelSubgrid_I(iGrid)
               IsOutExtended_I(nExtendedStencil) = IsOut_I(iGrid)
            end do
         end do
      else
         !\
         ! Calculate nDim-linear interpolation on uniform grid
         !/
         call interpolate_uniform((Xyz_D - XyzGrid_DII(:,1,1))*DxyzInv_D)
         !\
         ! Form index array and sort out zero weights
         !/
         nGridOut = 0 ; cTol2 = 2*cTol2
         do iGrid = 1, nGrid
            if(Weight_I(iGrid) < cTol2)CYCLE
            nGridOut = nGridOut + 1
            iIndexes_II(0,       nGridOut) = iProc_I( iGrid)
            iIndexes_II(nIndexes,nGridOut) = iBlock_I(iGrid)
            iIndexes_II(1:nDim,  nGridOut) = iCellIndexes_DII(:,1,iGrid)
            Weight_I(nGridOut) = Weight_I(iGrid)
         end do
      end if
    end subroutine generate_extended_stencil
    !==================
    subroutine generate_corner_stencil
      integer:: iGrid
      iLevel_I = iLevelSubgrid_I
      ! IsOut_I = IsOut_I
      do iGrid = 1, nGrid
         iIndexes_II(0,        iGrid) =  iProc_I(iGrid)
         iIndexes_II(nIndexes, iGrid) = iBlock_I(iGrid)
         if(iLevel_I(iGrid)==Coarse_)then
            XyzGrid_DI(:,iGrid) = XyzGrid_DII(:,1,iGrid)
            iIndexes_II(1:nDim,iGrid) = iCellIndexes_DII(:,1,iGrid)
         else
            XyzGrid_DI(:,iGrid) = XyzGrid_DII(:,9 - iGrid,iGrid)
            iIndexes_II(1:nDim,iGrid) = &
                 iCellIndexes_DII(:,9 - iGrid, iGrid)
         end if
      end do
    end subroutine generate_corner_stencil
    !==================
    subroutine sort_out
      !\
      ! Sorts out zero weights and repeating points
      !/
      integer:: nStored !To store starting nGridOut
      integer:: iLoc    !To find location of repeating index
      integer:: iGrid   !Loop variable
      !\
      ! Form index array and sort out zero weights
      !/
      nStored =  nGridOut 
      nGridOut = 0 
      cTol2 = 2*cTol2
      iLoc = 0
      ALL:do iGrid = 1, nStored
         if(Weight_I(iGrid) < cTol2)CYCLE
         do iLoc = 1, nGridOut
            if(iOrder_I(iLoc)==iOrder_I(iGrid))then
               Weight_I(iLoc) = Weight_I(iLoc) + Weight_I(iGrid)
               CYCLE ALL
            end if
         end do
         nGridOut = nGridOut + 1
         iOrder_I(nGridOut) = iOrder_I(iGrid)
         Weight_I(nGridOut) = Weight_I(iGrid)
      end do ALL
    end subroutine sort_out
  end subroutine interpolate_amr
  !=========================================================================
  subroutine  interpolate_amr2(&
       Xyz_D, XyzGrid_DI, iLevel_I, IsOut_I, iCaseExtended,&
       nGridOut, Weight_I, iOrder_I)
    use ModCubeGeometry, ONLY: iSortStencil2_II
    integer,parameter :: nGrid = 4, nDim = 2

    character(LEN=*),parameter:: NameSub='interpolate_amr2'

    !\
    !Input parameters
    !/
    !\
    !The location at which to interpolate the data
    !/
    real, intent(in):: Xyz_D(nDim) 

    !Grid point coordinates !2 coordinate, 4 points
    real  ,   intent(in) :: XyzGrid_DI(nDim,nGrid) 

    !The refinement level at each grid point. By one higher level of refinement
    !assumes the cell size reduced by a factor of 0.5
    integer,  intent(in) :: iLevel_I(nGrid)

    !\
    ! Logical which marks "ghost" points, not belonging to the computational 
    ! domain
    !/
    logical, intent(in) :: IsOut_I(nGrid)
    integer, intent(in) :: iCaseExtended
    !\
    !Output parameters
    !/
    !The number of grid points to be ultimately included into 
    !the interpolation stencil. If nGridOut < nGridIn, only the first 
    !nGridOut lines are meaningful in the output
    integer, intent(out) :: nGridOut

    !The weight coefficients array.
    real   , intent(out) :: Weight_I(nGrid)

    !Order(numbers) of grid points used for the interpolation
    integer, intent(out) :: iOrder_I(nGrid)

    integer ::  iDim !Loop variable
    !\
    !Parameters fully characterizing a pattern of refinement for a stencil
    !Calculated using iSortStencil_III array
    !/
    integer :: iCase, iGrid, iDir  

    !\
    ! 3 - iDir: 
    ! for coordinate x_ the perpendicular cordinate is y_ and vice versa 
    !/
    integer ::  iDirPerp              

    !\
    !Number of grid points returned from 1D routine
    !/
    integer :: nGridOut1    
    !\
    ! To find lines at "constant" cordinate
    !/
    real    :: dXyzSmall_D(nDim)       
    !\       
    ! for calculating on a uniform grid 
    !/
    real    :: dXyzInv_D(  nDim)  
    !\
    ! Misc
    !/
    integer:: iAux             
    real, dimension(nDim):: &
         Aux_D, X1_D = 0, X2_D = 0, X3_D = 0, X4_D = 0               

    !-------------
    !\
    ! Check the "ghost" grid points 
    !/
    if(count(IsOut_I)>=2)then
       do iDir=1, nDim
          if(all(IsOut_I(iOppositeSide_IDI(:,iDir,1))))then
             !\
             ! OppositeSide should be removed
             !/
             iOrder_I(1:2) = iSide_IDI(:,iDir,1)
             iOrder_I(3:4) = iOppositeSide_IDI(:,iDir,1)
          elseif(all(IsOut_I(iSide_IDI(:,iDir,1))))then
             !\
             ! Side should be removed
             !/
             iOrder_I(1:2) = iOppositeSide_IDI(:,iDir,1)
             iOrder_I(3:4) = iSide_IDI(:,iDir,1)
          else
             CYCLE
          end if
          call one_d_interpolate(iPoint1=1, iPoint2=2)
          nGridOut = nGridOut1
          RETURN
       end do
    end if
    !\
    ! Calculate the stencil size
    !/
    iAux = 2
    do iDim = 1, nDim
       Aux_D(iDim) = XyzGrid_DI(iDim,iAux) - XyzGrid_DI(iDim,1)
       dXyzSmall_D(iDim) = 0.10*abs(Aux_D(iDim))
       dXyzInv_D(iDim)   = 1/Aux_D(iDim)
       iAux = 2*iAux - 1
    end do
    !\
    !/
    iCase = i_case(iLevel_I)
    iGrid = iSortStencil2_II(Grid_,iCase)
    iDir  = iSortStencil2_II(Dir_, iCase)
    iCase = iSortStencil2_II(Case_,iCase)
    if(iCaseExtended == Transition2Edge_.and.iCase==OneCoarse_)&
         iCase = Transition2Edge_
    iOrder_I(1:2) = iSide_IDI(:,iDir,iGrid)
    iOrder_I(3:4) = iOppositeSide_IDI(:,iDir,iGrid)
    !\
    ! Frist conbinations of the resolution interfaces: 
    ! no resolution interface
    !/
    select case(iCase)
    case(Uniform_)
       !\                       C  C
       ! No refinement
       !/                       C  C
       Aux_D = (Xyz_D - XyzGrid_DI(:,1))*dXyzInv_D
       nGridOut = 4
       Weight_I(1) = (1 - Aux_D(1))*(1 - Aux_D(2))
       Weight_I(2) =      Aux_D(1) *(1 - Aux_D(2))
       Weight_I(3) = (1 - Aux_D(1))*     Aux_D(2)
       Weight_I(4) =      Aux_D(1) *     Aux_D(2)
       !\
       !Find resolution interfaces
       !/
    case(Face_)
       iDirPerp = 3 - iDir
       !\                             C  C
       ! Edges going along iDir-axis    F F   or   F F   ---> iDir
       !/                                           C  C
       Aux_D(iDirPerp) = &
            dXyzInv_D(iDirPerp)*(Xyz_D(iDirPerp) - &
            XyzGrid_DI(iDirPerp,iOrder_I(1) )  )       
       !\
       ! Interpolate along edge X1X2
       !/
       call one_d_interpolate(iPoint1=1,iPoint2=2)
       !\
       ! Apply weight for interpolation along another axis
       !/
       Weight_I(1:nGridOut1) =  &
            Weight_I(1:nGridOut1)*(1 - Aux_D(iDirPerp))
       !\
       ! May need to remove a grid point once behind the boundary
       !/
       if(nGridOut1==1)then
          iOrder_I(2:4) = iOrder_I((/3,4,2/))
       end if
       nGridOut = nGridOut1
       !\
       ! Interpolate along edge X3X4
       !/
       call one_d_interpolate(iPoint1=nGridOut +1, &
            iPoint2=nGridOut + 2)
       !\
       ! Apply weight for interpolation along another axis
       !/
       Weight_I(nGridOut + 1:nGridOut + nGridOut1) =  &
            Weight_I(nGridOut + 1:nGridOut + nGridOut1)*&
            Aux_D(iDirPerp)
       nGridOut = nGridOut + nGridOut1
    case(OneCoarse_)
       !     F-F
       !    / \|
       !   C---F
       call triangles(iTriangle1_I=(/2, 3, 1/), iTriangle2_I=(/2, 3, 4/))
    case(OneFine_,Rhombus_,Transition2Edge_)
       !   C---C
       !    \ /|
       !     F-C
       !   C-F
       !   |/\
       !   F--C
       !       |   |  
       !     3F--4F|   
       !------------    2F here
       !      \ /  !
       !       1C  | 
       call triangles(iTriangle1_I=(/1, 4, 2/), iTriangle2_I=(/1, 4, 3/))
    end select
  contains
    !========================
    subroutine one_d_interpolate(iPoint1, iPoint2)
      integer, intent(in):: iPoint1, iPoint2 !First and second point
      !\
      ! Check if the grid points are out of the grid boundary
      !/
      X1_D(iDir) = XyzGrid_DI(iDir,iOrder_I(iPoint1))
      X2_D(iDir) = XyzGrid_DI(iDir,iOrder_I(iPoint2))
      if(IsOut_I(iOrder_I(iPoint2)))then
         nGridOut1 = 1
         Weight_I(iPoint1) = 1
      elseif(IsOut_I(iOrder_I(iPoint1)))then
         iOrder_I((/iPoint1,iPoint2/)) = iOrder_I((/iPoint2,iPoint1/))
         nGridOut1 = 1
         Weight_I(iPoint1) = 1
      else
         nGridOut1 = 2
         Weight_I(iPoint2) = &
              (Xyz_D(iDir) - X1_D(iDir)) / (X2_D(iDir) - X1_D(iDir))
         Weight_I(iPoint1) = 1 - Weight_I(iPoint2) 
      end if
    end subroutine one_d_interpolate
    !=========================
    real function cross_product(a_D, b_D)
      real,dimension(nDim),intent(in) :: a_D, b_D
      !----------
      cross_product = a_D(x_)* b_D(y_) - a_D(y_)*b_D(x_)
    end function cross_product
    !=======
    subroutine triangles(iTriangle1_I,iTriangle2_I)
      integer,          intent(in) :: iTriangle1_I(3)
      integer,optional, intent(in) :: iTriangle2_I(3)
      real, dimension(nDim) :: X1_D, X2_D, X3_D
      !-------
      X1_D=XyzGrid_DI(:,iOrder_I(iTriangle1_I(1)))
      X2_D=XyzGrid_DI(:,iOrder_I(iTriangle1_I(2)))
      X3_D=XyzGrid_DI(:,iOrder_I(iTriangle1_I(3)))

      call triangle(X1_D, X2_D, X3_D)
      if(all(Weight_I(1:3)>=0.0))then
         iOrder_I = (/iOrder_I(iTriangle1_I(1)),iOrder_I(iTriangle1_I(2)),&
              iOrder_I(iTriangle1_I(3)), iOrder_I(10 - sum(iTriangle1_I))/) 
      else
         if(.not.present(iTriangle2_I))RETURN
         X1_D=XyzGrid_DI(:,iOrder_I(iTriangle2_I(1)))
         X2_D=XyzGrid_DI(:,iOrder_I(iTriangle2_I(2)))
         X3_D=XyzGrid_DI(:,iOrder_I(iTriangle2_I(3)))

         call triangle(X1_D, X2_D, X3_D)
         iOrder_I = (/iOrder_I(iTriangle2_I(1)),iOrder_I(iTriangle2_I(2)),&
              iOrder_I(iTriangle2_I(3)), iOrder_I(10 - sum(iTriangle2_I))/) 
      end if
    end subroutine triangles
    !============
    subroutine triangle(X1_D, X2_D, X3_D)
      real, dimension(nDim), intent(in):: X1_D, X2_D, X3_D
      !-------------
      Weight_I = 0; nGridOut = -1
      Weight_I(3) = cross_product(Xyz_D - X1_D,X2_D - X1_D)/&
           cross_product(X3_D - X1_D,X2_D - X1_D)
      if(Weight_I(3) >= 0.0)then
         nGridOut = 3
         Weight_I(2) = cross_product(X3_D-X1_D,Xyz_D - X1_D)/&
              cross_product(X3_D - X1_D,X2_D - X1_D)
         Weight_I(1) = 1 - Weight_I(2) - Weight_I(3)
      end if
    end subroutine triangle
  end subroutine interpolate_amr2
  !===========================
  subroutine  interpolate_amr3(&
       Xyz_D, XyzGrid_DI, iLevel_I, IsOut_I, iCaseExtended,&
       nGridOut, Weight_I, iOrder_I, IsCorner)
    use ModInterpolateBody, ONLY: resolution_corner
    integer,parameter :: nGrid = 8, nDim = 3
    character(LEN=*),parameter:: NameSub='interpolate_amr3'
    !\
    !Input parameters
    !/
    !\
    !The location at which to interpolate the data
    !/
    real   , intent(in) :: Xyz_D(nDim) 
    !\
    !Grid point coordinates !3 coordinate, 8 points
    !/
    real  ,   intent(in) :: XyzGrid_DI(nDim,nGrid)

    !The refinement level at each grid point. By one higher level
    ! of refinement assumes the cell size reduced by a factor of 0.5
    integer,  intent(in) :: iLevel_I(nGrid)
    !\
    ! Logical which marks "ghost" points, not belonging to the computational 
    ! domain
    !/
    logical, intent(in) :: IsOut_I(nGrid)
    integer, intent(in) :: iCaseExtended
    !\
    !Output parameters
    !/
    !The number of grid points to be ultimately included into 
    !the interpolation stencil If nGridOut < nGrid only the
    !first nGridOut lines in the output arrays are meaningful.
    !/
    integer, intent(out) :: nGridOut
    !\
    !The weight coefficients array.
    real   , intent(out) :: Weight_I(nGrid)
    !\
    !Order(numbers) of grid points used for the interpolation
    integer, intent(out) :: iOrder_I(nGrid)
    !\
    ! The logical is true, if (for some resolution combinations) the 
    ! basic stencil for the point at which to inpertpolate is not 
    ! applicable
    !/
    logical, intent(out):: IsCorner
    !\
    !Loop variables
    !/
    integer :: iDim, iLoop
    !\
    ! To find location of fine or coarse points
    !/
    !Number of grid points returned from 1D interpolation routine
    integer :: nGridOut1     
    !\
    ! Minimum and maximum values of coordinates
    !/
    real    :: XyzMin_D(nDim), XyzMax_D(nDim)   
    !\
    ! To find points with  "equal" cordinates
    !/
    real    :: dXyzSmall_D(nDim)      
    !\
    ! for calculating on a uniform grid 
    !/          
    real    :: dXyzInv_D(nDim)                
    real    :: Aux_D(nDim), AuxCoarse, AuxFine   ! Misc
    !\
    ! To call two-dimensional interpolation
    !/
    real    :: Xyz2_D(x_:y_), XyzGrid2_DI(x_:y_,4)
    integer :: nGridOut2, iOrder2_I(4),iLevel2_I(4)
    logical :: IsOut2_I(4)

    !\
    ! To find the sort of corner stencil
    !/
    integer:: iCase, iDir, iGrid
    !-------------
    Weight_I = 0
    !\
    ! We store in advance the 'basic' grid point
    ! and orientation of all possible stencil configurations
    ! and now we extracted this information
    !/ 
    iCase = i_case(iLevel_I)
    !\
    ! Check junction and transitions
    !/
    if(iSortStencil3_II(Case_,iCase) > Edge_&
         .and.iCaseExtended >= Transition2Corner_)&
         iCase = iCaseExtended
    if(iSortStencil3_II(Case_,iCase) == Edge_&
         .and.iCaseExtended == Transition2Edge_)&
         iCase = iCaseExtended 
    iGrid = iSortStencil3_II(Grid_,iCase)
    iDir  = iSortStencil3_II(Dir_,iCase)
    iCase = iSortStencil3_II(Case_,iCase)
    IsCorner = .false. 
    !\
    ! Check if the grid points are out of the grid boundary
    !/
    if(count(IsOut_I)>=4)then
       do iDim = 1, nDim
          if(all(IsOut_I(iOppositeFace_IDI(:,iDim,1))))then
             !\
             ! Upper face iDir should be removed
             !/
             iOrder_I(1:4) = iFace_IDI(        :,iDim,1)
             iOrder_I(5:8) = iOppositeFace_IDI(:,iDim,1)
          elseif(all(IsOut_I(iFace_IDI(:,iDim,1))))then
             !\
             ! Lower Face iDir should be removed
             !/
             iOrder_I(1:4) = iOppositeFace_IDI(:,iDim,1)
             iOrder_I(5:8) = iFace_IDI(:,iDim,1)
          else
             CYCLE
          end if
          !\
          ! iDir,1 + mod(iDir,3),1+mod(iDir+1,3) is 1,2,3 or 2,3,1 or 3,1,2
          !/
          Xyz2_D(x_) = Xyz_D(1 + mod(iDim,3))
          Xyz2_D(y_) = Xyz_D(1 + mod(iDim + 1, 3))
          XyzGrid2_DI(x_,1:4) = XyzGrid_DI(1 + mod(iDim    ,3),iOrder_I(1:4))
          XyzGrid2_DI(y_,1:4) = XyzGrid_DI(1 + mod(iDim + 1,3),iOrder_I(1:4))

          iLevel2_I(1:4) = iLevel_I(iOrder_I(1:4))
          IsOut2_I(1:4)     = IsOut_I(iOrder_I(1:4))

          call interpolate_amr2(&
               Xyz2_D,  XyzGrid2_DI, iLevel2_I, IsOut2_I,  iCase, &
               nGridOut, Weight_I(1:4), iOrder2_I)

          iOrder_I(1:4) = iOrder_I(iOrder2_I)
          RETURN
       end do
    end if
    !\ 
    ! Calculate the stencil size
    !/
    do iDim = 1, nDim
       dXyzSmall_D(iDim) = 0.10*(&
            maxval(XyzGrid_DI(iDim,:)) -  minval(XyzGrid_DI(iDim,:)))
       dXyzInv_D(iDim)   = 0.10/dXyzSmall_D(iDim)
    end do


    iOrder_I(1:4) = iFace_IDI(:,iDir,iGrid)
    iOrder_I(5:8) = iOppositeFace_IDI(:,iDir,iGrid)
    !\
    ! Reassign iCase for edges appearing from transitions, which require 
    ! stencil fix
    !/
    if(iCase==Edge_.and.&
         (iCaseExtended == TransitionJunction_.or.&
         (iCaseExtended == Transition2Corner_.and.&
         iDir/=iSortStencil3_II(Dir_,Transition2Corner_))))&
         iCase = Transition2Edge_
    select case(iCase)
    case(Face_)
       !\
       ! Faces going along plane iDir (which is the direction of normal)    
       !                               
       !/
       Xyz2_D(x_) = Xyz_D(1 + mod(iDir    ,3))
       Xyz2_D(y_) = Xyz_D(1 + mod(iDir + 1,3))
       !\
       ! Interpolate along lower face
       !/
       XyzGrid2_DI(x_,1:4) = XyzGrid_DI(1 + mod(iDir    ,3),iOrder_I(1:4))
       XyzGrid2_DI(y_,1:4) = XyzGrid_DI(1 + mod(iDir + 1,3),iOrder_I(1:4))
       iLevel2_I(        1:4) = iLevel_I(        iOrder_I(1:4))
       IsOut2_I = IsOut_I(iOrder_I(1:4))
       call interpolate_amr2(&
            Xyz2_D, XyzGrid2_DI, iLevel2_I, IsOut2_I, iCase, &
            nGridOut2, Weight_I(1:4), iOrder2_I)
       iOrder_I(1:4) = iOrder_I(iOrder2_I)
       !\
       ! Apply weight for interpolation along z-axis
       !/
       Aux_D(iDir) = dXyzInv_D(iDir)*&
            (Xyz_D(iDir) - XyzGrid_DI(iDir,iOrder_I(1)))
       Weight_I(1:nGridOut2) = Weight_I(1:nGridOut2)*(1 - Aux_D(iDir))
       !\
       ! May need to remove a grid points once behind the boundary
       !/
       if(    nGridOut2==2)then
          iOrder_I(3:8) = iOrder_I((/5,6,7,8,3,4/))
       elseif(nGridOut2==3)then
          iOrder_I(4:8) = iOrder_I((/5,6,7,8,4/))
       elseif(nGridOut2==1)then
          iOrder_I(2:8) = iOrder_I((/5,6,7,8,2,3,4/))
       end if
       nGridOut = nGridOut2
       !\
       ! Interpolate along upper face
       !/
       XyzGrid2_DI(x_,1:4) = XyzGrid_DI(1 + mod(iDir    ,3),&
            iOrder_I(nGridOut2 + 1:nGridOut2+4))
       XyzGrid2_DI(y_,1:4) = XyzGrid_DI(1 + mod(iDir + 1,3),&
            iOrder_I(nGridOut2 + 1:nGridOut2 + 4))
       iLevel2_I(1:4) = iLevel_I(iOrder_I(nGridOut2 + 1:nGridOut2 + 4))
       IsOut2_I = IsOut_I(iOrder_I(nGridOut2 + 1:nGridOut2 + 4))
       call interpolate_amr2(&
            Xyz2_D, XyzGrid2_DI, iLevel2_I, IsOut2_I, iCase, &
            nGridOut2, Weight_I(nGridOut+1:nGridOut+4), iOrder2_I)
       do iGrid=1,nGridOut2
          iOrder_I(nGridOut+iGrid)=iOrder_I(nGridOut+iOrder2_I(iGrid))
       end do
       !\
       ! Apply weight for interolation along kAxis
       !/
       Weight_I(nGridOut + 1:nGridOut + nGridOut2) = &
            Weight_I(nGridOut + 1:nGridOut + nGridOut2)*Aux_D(iDir)
       !\
       ! May need to remove a grid points once behind the boundary
       !/
       nGridOut = nGridOut + nGridOut2
    case(Edge_,Transition2Edge_)
       !\
       ! Edges going along resolution edge
       ! Opposite "faces" have the same geometry in projection
       ! term "face" is used to refer to 4 upper-, left-, 
       ! front-most etc. grid points they do not belong to plane
       !/      
       !\
       ! "Faces" have the same geometry in projection on iDir-plane
       !/
       Xyz2_D(x_:y_) = Xyz_D((/1 + mod(iDir,3),1 + mod(1 + iDir,3)/))
       !\
       ! Interpolation weight along iDir axis are calculated separately 
       ! for coarse and fine points
       !/
       XyzMin_D(iDir)   = minval(XyzGrid_DI(iDir,:),iLevel_I/=Fine_)
       AuxCoarse  = (Xyz_D(iDir) - XyzMin_D(iDir))/(&
            maxval(XyzGrid_DI(iDir,:),iLevel_I/=Fine_) - XyzMin_D(iDir))

       XyzMin_D(iDir)   = minval(XyzGrid_DI(iDir,:),iLevel_I==Fine_)
       AuxFine  = (Xyz_D(iDir) - XyzMin_D(iDir))/(&
            maxval(XyzGrid_DI(iDir,:),iLevel_I==Fine_) - XyzMin_D(iDir)) 
       !\
       ! Interpolate along lower face
       !/
       XyzGrid2_DI(x_:y_,1:4) = &
            XyzGrid_DI((/1 + mod(iDir,3),1 + mod(iDir+1,3)/),iOrder_I(1:4))
       iLevel2_I(        1:4) = iLevel_I( iOrder_I(1:4))
       IsOut2_I = .false.
       call interpolate_amr2(&
            Xyz2_D, XyzGrid2_DI, iLevel2_I, IsOut2_I, iCase, &
            nGridOut2, Weight_I(1:4), iOrder2_I)

       iOrder_I(1:4) = iOrder_I(iOrder2_I  )
       iOrder_I(5:8) = iOrder_I(iOrder2_I+4)

       !\
       ! May need to remove grid points from the stencil
       !/
       nGridOut = 2 * nGridOut2
       if(    nGridOut2==2)then
          iOrder_I(  3:4) = iOrder_I((/5,6/))
       elseif(nGridOut2==3)then
          iOrder_I(  4:6) = iOrder_I((/5,6,7/))
       end if
       !\
       ! Apply weight for interpolation along iDir
       !/
       do iGrid = 1,nGridOut2
          if(iLevel_I(iOrder_I(iGrid))==Fine_)then
             Weight_I(nGridOut2+iGrid) = Weight_I(iGrid)*     AuxFine
             Weight_I(          iGrid) = Weight_I(iGrid)*(1-  AuxFine)
          else
             if(IsOut_I(iOrder_I(iGrid)))then
                Weight_I(nGridOut2+iGrid) = Weight_I(iGrid)
                Weight_I(          iGrid) = 0
             elseif(&
                  IsOut_I(iOrder_I(nGridOut2 + iGrid)))then
                Weight_I(nGridOut2+iGrid) = 0
                Weight_I(          iGrid) = Weight_I(iGrid)
             else
                Weight_I(nGridOut2+iGrid) = Weight_I(iGrid)* AuxCoarse
                Weight_I(          iGrid) = Weight_I(iGrid)*(1 - AuxCoarse)
             end if
          end if
       end do
    case(TransitionJunction_)
       !\
       ! Check corner transition junction
       !/

       !\
       ! Junction around direction iDir
       !/
       ! Face with all Fine points is numbered 5:8
       ! Face with one to three Coarse points is numbered 1:4
       !
       !
       !  F5--F7
       !  |\ /|
       !  | C1|-----3
       !  |/|\|
       !  F6--F8
       !    |      4
       !    |
       !    2
       !
       !----------------------------------------------------------------
       !\
       ! Stencil part C1F5F6F7F8 is symmetric with respect to the plane 158
       ! Find if point Xyz is from the same side of the plane as point 2 is
       ! or as point 3 is
       !/
       if(sum((Xyz_D - XyzGrid_DI(:,iOrder_I(1)))*&
            (XyzGrid_DI(:,iOrder_I(6)) - XyzGrid_DI(:,iOrder_I(7)))) > 0)then
          iGrid = 2
       else
          iGrid = 3
       end if
       if(iLevel_I(iOrder_I(iGrid))==Coarse_)then
          call parallel_rays(Dir_D=0.50*(&
               XyzGrid_DI(:,iOrder_I(8)) + XyzGrid_DI(:,iOrder_I(5))) - &
               XyzGrid_DI(:,iOrder_I(1)),&
               iURectangle1_I=(/5, 6, 7, 8/),&
               iDTriangle1_I=(/1, iGrid, 8/))
       else
          call pyramids(iRectangular1_I = (/5, 6, 7, 8, 1/),&
               iRectangular2_I= (/iGrid, 4, iGrid +4, 8, 1/))
       end if
       IsCorner = nGridOut <= 3
    case(Transition2Corner_)
       call interpolate_corner_transition(iDir)
    case default
       call resolution_corner(&
            Xyz_D, XyzGrid_DI, iLevel_I, nGridOut, Weight_I, iOrder_I) 
    end select
  contains
    !==========================
    subroutine pyramids(&
         iTetrahedron1_I,iTetrahedron2_I,iTetrahedron3_I,&
         iTetrahedron4_I,iTetrahedron5_I,iTetrahedron6_I,&
         iRectangular1_I,iRectangular2_I,iRectangular3_I,&
         iTrapezoidal1_I,iTrapezoidal2_I)
      use ModInterpolateBody, ONLY:gen_pyramids=>pyramids
      integer,intent(in),optional,dimension(4)::&
           iTetrahedron1_I, iTetrahedron2_I, iTetrahedron3_I,&
           iTetrahedron4_I, iTetrahedron5_I, iTetrahedron6_I
      integer,intent(in),optional,dimension(5)::&
           iRectangular1_I, iRectangular2_I, iRectangular3_I,&
           iTrapezoidal1_I, iTrapezoidal2_I
      integer:: iOrderHere_I(1:5)
      real:: XyzGridHere_DI(nDim,nGrid)
      !---------------------------
      XyzGridHere_DI = XyzGrid_DI(:,iOrder_I)
      call gen_pyramids(&
           XyzGridHere_DI, Xyz_D,iOrderHere_I, Weight_I, nGridOut,&
           iTetrahedron1_I,iTetrahedron2_I,iTetrahedron3_I,&
           iTetrahedron4_I,iTetrahedron5_I,iTetrahedron6_I,&
           iRectangular1_I,iRectangular2_I,iRectangular3_I,&
           iTrapezoidal1_I,iTrapezoidal2_I)
      if(nGridOut > -1)iOrder_I(1:nGridOut) = &
           iOrder_I(iOrderHere_I(1:nGridOut))
    end subroutine pyramids
    !======================
    subroutine parallel_rays(Dir_D, &
         iDRectangle1_I, iDTrapezoid1_I,   &
         iDTriangle1_I, iDTriangle2_I, iDTriangle3_I,&
         iDTriangle4_I, iDTriangle5_I, iDTriangle6_I,&
         iURectangle1_I, iURectangle2_I, iURectangle3_I,&
         iUTrapezoid1_I,                             &
         iUTriangle1_I, iUTriangle2_I, iUTriangle3_I,&
         iUTriangle4_I)
      use ModInterpolateBody,ONLY: gen_parallel_rays=>parallel_rays
      !Direction of parallel rays
      real, intent(in) :: Dir_D(nDim)
      !\
      !Up subfaces ure in the direction of  Dir_D from Xyz point
      !Down faces are in the direction of -Dir_D
      !/
      integer, intent(in), optional, dimension(3)::&
           iDTriangle1_I, iDTriangle2_I, iDTriangle3_I,&
           iDTriangle4_I, iDTriangle5_I, iDTriangle6_I,&
           iUTriangle1_I, iUTriangle2_I, iUTriangle3_I,&
           iUTriangle4_I
      integer, intent(in), optional, dimension(4)::&
           iDRectangle1_I, iDTrapezoid1_I, iUTrapezoid1_I,&
           iURectangle1_I, iURectangle2_I, iURectangle3_I
      !\
      !The point belonging to one of the up and down subfaces
      !/ 
      real, dimension(nDim) :: XyzUp_D, XyzDown_D, X1_D, X2_D, X3_D, X4_D
      !Misc
      real:: AlphaUp, AlphaDown, XyzGridHere_DI(nDim,nGrid)
      integer:: nGridOutUp, nGridOutDown, iOrderHere_I(nGrid)
      !\
      ! Consider a ray passing through Xyz
      ! Solve equation:
      ! Xyz + Dir*AlphaUp = XyzUp 
      ! its solution is 
      ! AlphaUp = (UX1-Xyz)\cdot[(UX2-UX1)\times(UX3-UX1)]/&
      !          Dir\cdot[(UX2-UX1)\times(UX3-UX1)]
      ! XyzUp belongs to the up subface
      !/
      !\
      ! Solve equation:
      ! Xyz - Dir*AlphaDown = XyzDown 
      ! its solution is 
      ! AlphaDown = (Xyz-DX1)\cdot[(DX2-DX1)\times(DX3-DX1)]/&
      ! Dir\cdot[(DX2-DX1)\times(DX3-DX1)]
      ! XyzUp belongs to the Down subface
      !/
      !\
      ! As long as Xyz=(XyzDown*AlphaUp + XyzUp*AlphaDown)/(AlphaUp+AlphaDown)
      ! the weights for the upper face should be multiplied by 
      ! AlphaDown/(AlphaUp+AlphaDown)
      ! for the down face - by AlphaUp/(AlphaUp+AlphaDown)
      !/
      XyzGridHere_DI = XyzGrid_DI(:,iOrder_I)
      call gen_parallel_rays(&
           XyzGridHere_DI, Xyz_D, iOrderHere_I, Weight_I, nGridOut, Dir_D, &
           iDRectangle1_I, iDTrapezoid1_I,   &
           iDTriangle1_I, iDTriangle2_I, iDTriangle3_I,&
           iDTriangle4_I, iDTriangle5_I, iDTriangle6_I,&
           iURectangle1_I, iURectangle2_I, iURectangle3_I,&
           iUTrapezoid1_I,                             &
           iUTriangle1_I, iUTriangle2_I, iUTriangle3_I,&
           iUTriangle4_I)
      if(nGridOut > 0)iOrder_I(1:nGridOut) =  &
           iOrder_I(iOrderHere_I(1:nGridOut))
    end subroutine parallel_rays
    !===========================
    subroutine interpolate_corner_transition(iEdgeDir)
      integer, intent(in) :: iEdgeDir
      integer             :: iAxis,jAxis,kAxis
      !\
      ! Auxilary variables to determine vertical weight
      !/
      real :: AlphaKAxisCoarse, AlphaKAxisFine
      real :: dXyzUp, dXyzDown 
      integer :: nFine, iFine_I(3)
      integer :: nCoarse, iCoarse_I(3)
      !\
      ! Subroutine intepolates in transitional near a resolution corner
      ! Points facing a corner a numbered 5:8
      !
      !-------------------------------------------------------------------
      kAxis = iEdgeDir
      iAxis = 1 + mod(iEdgeDir    ,3)
      jAxis = 1 + mod(iEdgeDir + 1,3)
      Xyz2_D(x_:y_) = Xyz_D((/iAxis,jAxis/))

      !\
      ! Interpolate along face 1:4
      !/
      XyzGrid2_DI(x_:y_,1:4) = XyzGrid_DI((/iAxis,jAxis/),iOrder_I(1:4))
      iLevel2_I(        1:4) = iLevel_I(                  iOrder_I(1:4))
      IsOut2_I = IsOut_I(iOrder_I(1:4))
      call interpolate_amr2(&
           Xyz2_D, XyzGrid2_DI, iLevel2_I, IsOut2_I, Edge_,&
           nGridOut2, Weight_I(1:4), iOrder2_I)

      !\
      ! Rearrange points according to 2D interpolation output
      !/
      iOrder_I(1:4) = iOrder_I(iOrder2_I  )
      iOrder_I(5:8) = iOrder_I(iOrder2_I+4)
      !\
      ! Find number of Fine points in the output
      ! Move them to the begining of iOrder_I
      !/
      nFine   = 0
      do iGrid = 1, nGridOut2
         if(iLevel_I(iOrder_I(iGrid))==Fine_)then
            nFine   = nFine   + 1
            iOrder_I((/nFine,  iGrid  /)) = iOrder_I((/iGrid,  nFine  /))
            Weight_I((/nFine,  iGrid  /)) = Weight_I((/iGrid,  nFine  /))
            iOrder_I((/nFine+4,iGrid+4/)) = iOrder_I((/iGrid+4,nFine+4/))
         end if
      end do
      nCoarse = nGridOut2 - nFine
      !\
      ! Find vertical weights for Coarse
      ! if above each Coarse point in the output there is another Coarse point
      ! then that another Coarse point also should be used in the interpolation
      ! if not substitute corresponding point with Coarse point from range 1:4
      !/
      select case(nCoarse)
      case(0)
         !\
         ! no Coarse points
         ! give AlphaKAxisCoarse value 1 for calculating AlphaKAxisFine
         !/
         AlphaKAxisCoarse = 1
      case(1)
         if(iLevel_I(iOrder_I(nFine+5))==Coarse_)then
            !\
            !     F-----F
            !    /|    /|
            !   C-----C |   ---->kAxis
            !    \|    \|
            !     F-----F
            !/
            AlphaKAxisCoarse = &
                 (Xyz_D(kAxis) - XyzGrid_DI(kAxis,iOrder_I(nFine + 5)))/&
                 (XyzGrid_DI(kAxis,iOrder_I(nFine + 1)) - &
                 XyzGrid_DI(kAxis,iOrder_I(nFine +5)))
         else
            !\
            !     F-----F
            !     |\   /|
            !     |  C  |   ---->kAxis
            !     |/   \|
            !     F-----F
            !/
            AlphaKAxisCoarse = 1
            iOrder_I(nFine+5) = iOrder_I(nFine+1)
         end if
      case(2)
         if(iLevel_I(iOrder_I(nFine+5))==Coarse_.and.&
              iLevel_I(iOrder_I(nFine+6))==Coarse_)then
            !\
            !     C---C
            !    /|\  |
            !   F---F |      ---->kAxis
            !    \|/  |
            !     C---C
            !/
            call parallel_rays(&
                 Dir_D = &
                 XyzGrid_DI(:,iOrder_I(4))-XyzGrid_DI(:,iOrder_I(2)),&
                 iDRectangle1_I = &
                 (/2,3,6,7/),&
                 iUTriangle1_I = &
                 (/1,2, 5/))
            IsCorner = .false.
            if(nGridOut>0)RETURN
            AlphaKAxisCoarse = 0.5
         else
            !\
            !     F-----F       kAxis
            !    /|    /|       ^
            !   C-----C |       |
            !    \|    \|       |
            !     F-----F
            !
            !           OR
            !
            !     C
            !    /|\
            !   F---F      ---->kAxis
            !    \|/
            !     C
            !/
            AlphaKAxisCoarse = 1
            iOrder_I(nFine+5) = iOrder_I(nFine+1)
            iOrder_I(nFine+6) = iOrder_I(nFine+2)

         end if
      end select

      !\
      ! Find vertical weights for Fine
      !/
      if(nGridOut2==2)Weight_I(3:4)=0
      if(AlphaKAxisCoarse==1)then
         dXyzDown = Xyz_D(kAxis) - &
              sum( XyzGrid_DI(kAxis,iOrder_I(1:nGridOut2))*&
              Weight_I(1:nGridOut2) )
         dXyzUp   = Xyz_D(kAxis) - sum(&
              XyzGrid_DI(kAxis,iOrder_I(5:4 + nGridOut2))*&
              Weight_I(1:nGridOut2))
      else
         if(nFine==2)then
            dXyzDown = Xyz_D(kAxis) - XyzGrid_DI(kAxis,iOrder_I(1))
            dXyzUp   = Xyz_D(kAxis) - XyzGrid_DI(kAxis,iOrder_I(5))
         else
            dXyzDown = 1
            dXyzUp   = 1
         end if
      end if

      IsCorner = dXyzUp * dXyzDown > 0.0
      if(IsCorner)RETURN

      if(dXyzUp==0)then
         AlphaKAxisFine = 0.0
      else
         AlphaKAxisFine = abs(dXyzUp)/(abs(dXyzUp) + abs(dXyzDown))
      end if

      !\
      ! Remove points
      ! Apply vertical weight
      !/
      if(AlphaKAxisCoarse == 1)then
         nGridOut = nGridOut2 + nFine
         iOrder_I(nGridOut2 + 1:nGridOut) = iOrder_I(5:4+nFine)
         Weight_I(nGridOut2 + 1:nGridOut) = &
              Weight_I(1:  nFine)*(1 - AlphaKAxisFine)
         Weight_I(          1:nFine   ) = Weight_I(1:  nFine)*AlphaKAxisFine
      else
         nGridOut = 2*nGridOut2
         iOrder_I(nGridOut2 + 1:nGridOut       ) = iOrder_I(5:4 + nGridOut2)
         Weight_I(nGridOut2  +1:nGridOut2+nFine) = &
              Weight_I(1:nFine)*(1 - AlphaKAxisFine)
         Weight_I(nGridOut2 + nFine+1:nGridOut ) = &
              Weight_I(nFine + 1:nGridOut2)*(1 - AlphaKAxisCoarse)
         Weight_I(1:nFine) = Weight_I(1:nFine)*AlphaKAxisFine
         Weight_I(nFine + 1:nGridOut2) = &
              Weight_I(nFine + 1:nGridOut2)*AlphaKAxisCoarse
      end if
    end subroutine interpolate_corner_transition
  end subroutine interpolate_amr3
  !==================================
  integer function i_case(iLevel_I)
    integer, intent(in):: iLevel_I(:)
    integer:: iGrid, nGrid
    !--------------------
    nGrid = size(iLevel_I,DIM=1)
    !\
    ! Generate iCase from 'binary' iLevel_I
    !/
    i_case = iLevel_I(nGrid)
    do iGrid = nGrid - 1, 1, -1
       i_case = i_case + i_case + iLevel_I(iGrid)
    end do
  end function i_case

  !=======================================================================
  subroutine fix_basic_stencil3( Xyz_D, dXyzInv_D, XyzGrid_DI, iLevel_I,&
       XyzStencil_D, iCase, IsCorner)
    integer, parameter:: nDim = 3, nGrid = 8
    !\
    ! Point where to interpolate; inverse of (/Dx,Dy,Dz/)
    !/
    real,    intent(in) :: Xyz_D(nDim), dXyzInv_D(nDim)
    !\
    ! The rectangular grid combining centers of the refined subgrids and
    ! coarse vertexes
    !/
    real, intent(in):: XyzGrid_DI(nDim,nGrid)
    !\
    ! The refinement level pattern of the extended stencil
    !/
    integer, intent(in)::iLevel_I(nGrid)
    !\
    ! Point about which to construct basic stencil.
    !/

    real, intent(out) :: XyzStencil_D(nDim)

    !\
    ! For three-dimensional corner within this routine
    ! it is more convenient to figure out if Xyz point is
    ! inside the domain of transition from an edge to the central
    ! corner part or it belongs to a junction of such two domains
    !/
    integer, intent(out):: iCase

    !\
    ! Is true if the point is inside the resolution corner
    !/
    logical, intent(out):: IsCorner
    !\
    !Dimensionless displacement from the first grid point to Xyz_D
    !/
    real :: Dimless_D(nDim)
    !\
    !Three components of Discr_D equal to -1,0 or 1 each. The first component 
    !equals -1, if Xyz point is close to refined face x=0, +1 if it is close
    !to refined face x=1, 0 otherwise
    !/
    integer:: iDiscr_D(nDim), iDir, iDim, iGrid, jGrid
    !\ 
    !Dir = 0 - corner, positive iDirs - face, negative iDir - edges, which
    !occur at the intersection of the refined faces, which are close to
    !Xyz
    !/
    integer, parameter, dimension(2,-1:1,-1:1,-1:1):: &
         iGridDir_IIII = reshape((/&
         1, 0,    1,-1,  2, 0, 1,-2,1,3,2,-2,  3, 0, 3,-1, 4, 0,   &!z=0!
         1,-3,    1, 2,  2,-3, 1, 1,0,0,2, 1,  3,-3, 3, 2, 4,-3,   &!   !
         5, 0,    5,-1,  6, 0, 5,-2,5,3,6,-2,  7, 0, 7,-1, 8, 0 /),&!z=1!
         (/2,3,3,3/))
    !x=0,y=0!y=0 !x=1,y=0!x=0 !    !x=1 !x=0,y=1!y=1 !x=1,y=1!
    !\
    !Transitions junction's signatures
    !/
    real:: XyMin, z
    !\
    ! Average coordinates for the center of the resolution corner,
    ! for centers of faces or edges and distance from them to Xyz
    !/
    real   :: XyzAvr_D(nDim), XyzAvr_DD(nDim, nDim), Distance_D(nDim) 
    !\
    ! To treat Edges
    !/
    logical:: IsEdge
    integer:: iDirEdge
    !----------------------
    XyzStencil_D = Xyz_D
    IsCorner = .false.
    !For 2 dimensions the search of a basic stencil based on the distance
    !from the point to the point coordinates works well.
    iCase = i_case(iLevel_I)
    if(iSortStencil3_II(Case_,iCase) <= Face_)RETURN
    IsEdge = iSortStencil3_II(Case_,iCase) == Edge_
    if(IsEdge) iDirEdge = iSortStencil3_II(Dir_,iCase)
    !\
    !There is also no need to look for transitions in faces
    !/
    Dimless_D = (Xyz_D - XyzGrid_DI(:,1))*DxyzInv_D 
    iDiscr_D = 0
    do iDim = 1, nDim
       if(Dimless_D(iDim) <  0.250.and.any(iLevel_I(&
            iFace_IDI(:,iDim,1))==Fine_))iDiscr_D(iDim) = -1
       if(Dimless_D(iDim) >=  0.750.and.any(iLevel_I(&
            iOppositeFace_IDI(:,iDim,1))==Fine_))iDiscr_D(iDim) = 1
    end do
    if(IsEdge) iDiscr_D(iDirEdge) = 0
    iGrid = iGridDir_IIII(1, iDiscr_D(1), iDiscr_D(2), iDiscr_D(nDim))
    if(iGrid==0)then
       IsCorner = .not.IsEdge
       RETURN
    end if
    !\
    !Xyz point is in the transition or transition junction domain
    !/
    XyzAvr_D = 0.50*(XyzGrid_DI(:,1) + XyzGrid_DI(:,8))
    iDir = iGridDir_IIII(2, iDiscr_D(1), iDiscr_D(2), iDiscr_D(nDim))
    if(iDir > 0)then
       !Xyz point os close to a single refined face - transition!
       iCase = Transition2Corner_
       iSortStencil3_II(Grid_,Transition2Corner_) = iGrid
       iSortStencil3_II(Dir_ ,Transition2Corner_) = iDir
       !\
       ! Displace the stencil center toward the face center
       !/
       XyzStencil_D(1 + mod(iDir,3)) = XyzAvr_D(1 + mod(iDir,3))
       XyzStencil_D(1 + mod(iDir + 1,3)) = &
            XyzAvr_D(1 + mod(iDir + 1,3))
    elseif(iDir < 0)then
       !\
       !The point is close to two faces intersescting at the
       !edge of direction -iDir
       !/
       iDir = -iDir
       jGrid = iEdge_ID(iGrid,iDir)
       !\
       ! Edge consists of the coarse points or finer subgrid 
       ! iGrid,jGrid
       !/
       if(iLevel_I(iGrid)==Coarse_.and.iLevel_I(jGrid)==Fine_)then
          !\
          !Jucntion of two transition regions which goes along
          !the edge of direction iDir and closed with fine
          !subgrid at jGrid end
          iCase = TransitionJunction_
          iSortStencil3_II(Grid_,TransitionJunction_) = iGrid
          iSortStencil3_II(Dir_ ,TransitionJunction_) = iDir    
          XyzStencil_D(iDir) = XyzAvr_D(iDir)
       elseif(iLevel_I(jGrid)==Coarse_.and.iLevel_I(iGrid)==Fine_)then
          !\
          !Jucntion of two transition regions which goes along
          !the edge of direction iDir and closed with fine
          !subgrid at iGrid end
          iCase = TransitionJunction_
          iSortStencil3_II(Grid_,TransitionJunction_) = jGrid
          iSortStencil3_II(Dir_ ,TransitionJunction_) = iDir   
          XyzStencil_D(iDir) = XyzAvr_D(iDir)
       else
          !\
          ! We need to judge to which of the two transition regions
          ! intersecting along the coarse edge, that is without 
          ! forming a specific junction region as closed with the  
          ! fine subface, point Xyz belongs 
          !/
          do iDim = 0,1
             XyzAvr_DD(:,1+iDim) = 0.50*(&
                  XyzGrid_DI(:,&
                  iFace_IDI(1,1 + mod(iDir + iDim,3),iGrid)) +&
                  XyzGrid_DI(:,&
                  iFace_IDI(4,1 + mod(iDir + iDim,3),iGrid)) )
             Distance_D(1 + iDim) = &
                  sum(((Xyz_D - XyzAvr_DD(:,1+iDim))*DxyzInv_D)**2) 
          end do
          iDim = minloc(Distance_D(1:2),DIM=1)
          iDir = 1 + mod(iDir + iDim - 1,3)
          iCase = Transition2Corner_
          iSortStencil3_II(Grid_,Transition2Corner_) = iGrid
          iSortStencil3_II(Dir_ ,Transition2Corner_) = iDir
          !\
          ! Displace the stencil center toward the face center
          !/
          XyzStencil_D(1 + mod(iDir,3)) = XyzAvr_D(1 + mod(iDir,3))
          XyzStencil_D(1 + mod(iDir + 1,3)) = &
               XyzAvr_D(1 + mod(iDir + 1,3))
       end if
    else
       !\
       ! Xyz point is close to three refined planes
       !/
       select case(count(iLevel_I(iEdge_ID(iGrid,:))==Fine_))
       case(0)
          !\
          ! We need to judge to which of the three transition regions
          ! intersecting at the corner without forming a specific
          ! junction region, point Xyz belongs 
          !/
          do iDim = 1,nDim
             XyzAvr_DD(:,iDim) = 0.50*(&
                  XyzGrid_DI(:,iFace_IDI(1,iDim,iGrid)) +&
                  XyzGrid_DI(:,iFace_IDI(4,iDim,iGrid)) )
             Distance_D(iDim) = &
                  sum(((Xyz_D - XyzAvr_DD(:,iDim))*DxyzInv_D)**2) 
          end do
          iDir = minloc(Distance_D,DIM=1)
          iCase = Transition2Corner_
          iSortStencil3_II(Grid_,Transition2Corner_) = iGrid
          iSortStencil3_II(Dir_ ,Transition2Corner_) = iDir
          !\
          ! Displace the stencil center toward the face center
          !/
          XyzStencil_D(1 + mod(iDir,3)) = XyzAvr_D(1 + mod(iDir,3))
          XyzStencil_D(1 + mod(iDir + 1,3)) = &
               XyzAvr_D(1 + mod(iDir + 1,3))
       case(2,3)
          !\
          ! We need to judge to which of the two or three transition 
          ! junctions point Xyz belongs 
          !/
          do iDim = 1,nDim
             !\
             !If there is no junction region along this direction,
             !that is iLevel_I(iEdge_ID(iGrid,iDim))=0, move the
             !center of this edge far apart, otherwise take the 
             !center of the edge  
             !/
             XyzAvr_DD(:,iDim) = (&
                  XyzGrid_DI(:,iGrid)*iLevel_I(iEdge_ID(iGrid,iDim))&
                  + XyzGrid_DI(:,iEdge_ID(iGrid,iDim)) )/&
                  (iLevel_I(iEdge_ID(iGrid,iDim)) + 1)
             Distance_D(iDim) = &
                  sum(((Xyz_D - XyzAvr_DD(:,iDim))*DxyzInv_D)**2) 
          end do
          iDir = minloc(Distance_D,DIM=1)
          iCase = TransitionJunction_
          iSortStencil3_II(Grid_,TransitionJunction_) = iGrid
          iSortStencil3_II(Dir_ ,TransitionJunction_) = iDir
          XyzStencil_D(iDir) = XyzAvr_D(iDir)
       case(1)
          !\
          ! We need to find the only transition junction and then  
          ! judge if point Xyz belongs to this transition junction  
          ! or to the transition region near the lower face of the
          ! same direction. 
          !/
          iDir = minloc(iLevel_I(iEdge_ID(iGrid,:)), MASK=&
               iLevel_I(iEdge_ID(iGrid,:))==Fine_, DIM=1)
          !
          !        ^     7F----8F
          !   iDir |     /    /
          !            5F----6F
          !                  /3C
          !                 /     4F
          !                /________
          !               1C       2C
          !\
          ! As a separator between the transition junction 1C5F6F7F8F
          ! and a transition region near face 1C2C3C4F we use a
          !  surface, z = min(x,y)*Alpha, which passes through
          !  point 4F at Alpha = 1/3 and through point 8F at Alpha=3. 
          !/
          z = (Xyz_D(iDir) - XyzGrid_DI(iDir,iGrid))/  &
               (XyzGrid_DI(iDir,iEdge_ID(iGrid,iDir)) -&
               XyzGrid_DI(iDir,iGrid))
          XyMin = min( (Xyz_D(1 + mod(iDir,3))        -&
               XyzGrid_DI(1 + mod(iDir,3),iGrid))/     &
               (XyzGrid_DI(1 + mod(iDir,3),            &
               iEdge_ID(iGrid,1 + mod(iDir,3)))       -&
               XyzGrid_DI(1 + mod(iDir,3),iGrid)),     &
               (Xyz_D(1 + mod(iDir + 1,3))            -&
               XyzGrid_DI(1 + mod(iDir + 1,3),iGrid))/ &
               (XyzGrid_DI(1 + mod(iDir + 1,3),        &
               iEdge_ID(iGrid,1 + mod(iDir + 1,3)))   -&
               XyzGrid_DI(1 + mod(iDir + 1,3),iGrid)) )
          if(z  >= 3*XyMin)then
             iCase = TransitionJunction_
             iSortStencil3_II(Grid_,TransitionJunction_) = iGrid
             iSortStencil3_II(Dir_ ,TransitionJunction_) = iDir
             XyzStencil_D(iDir) = XyzAvr_D(iDir)
          elseif(z <= XyMin/3)then
             iCase = Transition2Corner_
             iSortStencil3_II(Grid_,Transition2Corner_) = iGrid
             iSortStencil3_II(Dir_ ,Transition2Corner_) = iDir
             !\
             ! Displace the stencil center toward the face center
             !/
             XyzStencil_D(1 + mod(iDir,3)) = XyzAvr_D(1 + mod(iDir,3))
             XyzStencil_D(1 + mod(iDir + 1,3)) = &
                  XyzAvr_D(1 + mod(iDir + 1,3))
          else
             IsCorner = .true.
             XyzStencil_D = XyzAvr_D 
          end if
       end select
    end if
    if(IsEdge)then
       iCase = Transition2Edge_
       iSortStencil3_II(Dir_, iCase) = iDirEdge
       XyzStencil_D(iDirEdge) = Xyz_D(iDirEdge)
    end if
  end subroutine fix_basic_stencil3
  !=======================================================================
  subroutine fix_basic_stencil2( Xyz_D, dXyzInv_D, XyzGrid_DI, iLevel_I,&
       XyzStencil_D, iCase)
    integer, parameter:: nDim = 2, nGrid = 4
    !\
    ! Point where to interpolate; inverse of (/Dx,Dy,Dz/)
    !/
    real,    intent(in) :: Xyz_D(nDim), dXyzInv_D(nDim)
    !\
    ! The rectangular grid combining centers of the refined subgrids and
    ! coarse vertexes
    !/
    real, intent(in):: XyzGrid_DI(nDim,nGrid)
    !\
    ! The refinement level pattern of the extended stencil
    !/
    integer, intent(in)::iLevel_I(nGrid)
    !\
    ! Point about which to construct basic stencil.
    !/

    real, intent(out) :: XyzStencil_D(nDim)

    !\
    ! For three-dimensional corner within this routine
    ! it is more convenient to figure out if Xyz point is
    ! inside the domain of transition from an edge to the central
    ! corner part or it belongs to a junction of such two domains
    !/
    integer, intent(out):: iCase
    !\
    !Dimensionless displacement from the first grid point to Xyz_D
    !/
    real :: Dimless_D(nDim)
    !\
    !Three components of Discr_D equal to -1,0 or 1 each. The first component 
    !equals -1, if Xyz point is close to refined face x=0, +1 if it is close
    !to refined face x=1, 0 otherwise
    !/
    integer:: iDiscr_D(nDim), iDir, iDim, iGrid, jGrid
    !\ 
    !Dir = 0 - corner, positive iDirs - face, negative iDir - edges, which
    !occur at the intersection of the refined faces, which are close to
    !Xyz
    !/
    integer, parameter, dimension(2,-1:1,-1:1):: &
         iGridDir_III = reshape((/&
         1,0,     1,2,  2,0,  1,1,  0,0, 2,1,  3,0,   3,2, 4,0/),(/2,3,3/))
    !x=0,y=0!y=0 !x=1,y=0!x=0 !    !x=1 !x=0,y=1!y=1 !x=1,y=1!
    real:: XyMin, z
    !\
    ! Average coordinates for the center of the resolution corner,
    ! for centers of faces or edges and distance from them to Xyz
    !/
    real   :: XyzAvr_D(nDim), XyzAvr_DD(nDim, nDim), Distance_D(nDim) 
    !----------------------
    XyzStencil_D = Xyz_D
    iCase = i_case(iLevel_I)
    !\                                    F F
    !We do not care about trapezoids like  C   C
    !/
    if(iSortStencil3_II(Case_,iCase) <= Face_)RETURN
    Dimless_D = (Xyz_D - XyzGrid_DI(:,1))*DxyzInv_D 
    iDiscr_D = 0
    do iDim = 1, nDim
       if(Dimless_D(iDim) <  0.250.and.any(iLevel_I(&
            iSide_IDI(:,3 - iDim,1))==Fine_))iDiscr_D(iDim) = -1
       if(Dimless_D(iDim) >=  0.750.and.any(iLevel_I(&
            iOppositeSide_IDI(:,3 - iDim,1))==Fine_))iDiscr_D(iDim) = 1
    end do
    iGrid = iGridDir_III(1, iDiscr_D(1), iDiscr_D(2))
    XyzAvr_D = 0.50*(XyzGrid_DI(:,1) + XyzGrid_DI(:,4))
    if(iGrid==0)then 
       !\
       ! The point is in the close proximity of the resolution edge
       !/
       XyzStencil_D = XyzAvr_D 
       RETURN
    end if
    !\
    !Xyz point is in the transition domain
    !/   
    iDir = iGridDir_III(2, iDiscr_D(1), iDiscr_D(2))
    if(iDir== 0)then
       !\
       !The point is close to two faces
       !/
       do iDim = 1,nDim
          XyzAvr_DD(:,iDim) = 0.50*(&
               XyzGrid_DI(:,&
               iSide_IDI(1,3 - iDim,iGrid)) +&
               XyzGrid_DI(:,&
               iSide_IDI(2,3 - iDim,iGrid)) )
          Distance_D(iDim) = &
               sum(((Xyz_D - XyzAvr_DD(:,iDim))*DxyzInv_D)**2) 
       end do
       iDir = minloc(Distance_D(1:2),DIM=1)
       !\
       !else Xyz point is close to a single refined face 
       !/
    end if
    iCase = Transition2Edge_
    !\
    ! Displace the stencil center toward the face center
    !/
    XyzStencil_D(3 - iDir) = XyzAvr_D(3 - iDir)
  end subroutine fix_basic_stencil2
  !=====================
  subroutine generate_basic_stencil(&
       nDim, XyzStencil_D, nExtendedStencil, XyzExtended_DI, dXyzInv_D, &
       iOrderExtended_I)
    !\
    ! Dimensionality; number of points in the extended stencil
    !/
    integer, intent(in) :: nDim, nExtendedStencil
    !\
    ! Point about which to construct basic stencil; inverse of (/Dx,Dy,Dz/)
    !/
    real,    intent(in) :: XyzStencil_D(nDim), dXyzInv_D(nDim)
    !\
    ! Points of the extended stencil
    !/
    real,    intent(in) :: XyzExtended_DI(nDim, nExtendedStencil)
    !\
    ! Basic stencil - point numbers in XyzExtended array
    !/
    integer,    intent(out):: iOrderExtended_I(2**nDim)
    integer:: iGrid, nGrid, iPoint !Loop variables
    logical:: IsMask_I(nExtendedStencil) 
    logical:: IsBelowExtended_DI(nDim,nExtendedStencil)
    real   :: Distance_I(nExtendedStencil)
    !\
    ! Arrays of logical values Xyz_D < XyzGrid_DI(:,iGrid) for iGrid'th grid
    ! point of the basic stencil for point Xyz
    !/
    logical, parameter, dimension(3,8):: IsBelow_DI=reshape((/&
         .false., .false., .false., &
         .true. , .false., .false., &
         .false., .true. , .false., &
         .true. , .true. , .false., &
         .false., .false., .true. , &
         .true. , .false., .true. , &
         .false., .true. , .true. , &
         .true. , .true. , .true. /)&
         , (/3,8/) )
    !-------------------------------------
    nGrid = 2**nDim
    iOrderExtended_I = -1

    do iPoint = 1,nExtendedStencil
       IsBelowExtended_DI(:,iPoint) = XyzStencil_D < XyzExtended_DI(:,iPoint)
       Distance_I(iPoint) = sum( &
            ((XyzStencil_D - XyzExtended_DI(:,iPoint))*dXyzInv_D)**2)
    end do
    do iGrid = 1, nGrid
       do iPoint = 1, nExtendedStencil
          !\
          ! Mark all the candidates for the role of iGrid point of  
          ! the basic stencil
          !/
          IsMask_I(iPoint) = all(&
               IsBelowExtended_DI(:,iPoint).eqv.IsBelow_DI(1:nDim,iGrid))
       end do
       !\
       !Among the candidates chose the closest one
       !/
       iOrderExtended_I(iGrid) = minloc(Distance_I, MASK = IsMask_I, DIM=1)
    end do
  end subroutine generate_basic_stencil
end module ModInterpolateAMR
