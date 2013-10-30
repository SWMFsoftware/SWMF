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
  use ModCubeGeometry, ONLY: Uniform_, Face_, Edge_, FiveTetrahedra_, & 
       OneFine_, OneCoarse_, FineMainDiag_, FineFaceDiag_,            &
       CoarseMainDiag_, CoarseFaceDiag_, FineEdgePlusOne_,            &
       ThreeFineOnFace_, CoarseEdgePlusOne_, ThreeCoarseOnFace_,      &
       ThreeCoarseOnFacePlusOne_, CoarseChain_, Transition2Edge_,     &
       Transition2Corner_, TransitionJunction_
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

  integer, parameter:: &
       x_       = 1,               &
       y_       = 2,               &
       z_       = 3
  
contains
 !====================================================================
  subroutine interpolate_amr(nDim, XyzIn_D, nIndexes, nCell_D, find, &
       nGridOut, Weight_I, iIndexes_II, IsSecondOrder)
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
    ! Output parameter of interpolate_amr2,3
    ! If .true., the basic stencil should be re-evaluated
    !/
    logical                    :: DoStencilFix
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
    ! routine; inverse of (/Dx, Dy, Dz/) is reused; if .true.
    ! DoFixStencil the input coordinate to reevaluate the basic stencil
    ! should be provided: 
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
    !------------------------
    cTol2 = cTol**(nByteReal/4)
    IsExtended = .false.; UseGhostPoint = .false.

    nGrid = 2**nDim !Number of points in a basic stencil
    !\
    ! Initialize 
    !/
    iOrder_I    = 0; iIndexes_II = 0; Weight_I    = 0
    nGridOut = -1; Xyz_D = XyzIn_D ; IsOut_I = .false.

    call generate_extended_stencil
    !\
    !The interpolation may be done inside generate_extended_stencil 
    !/
    if(nGridOut > 0) then
       if(present(IsSecondOrder))IsSecondOrder = .not. UseGhostPoint
       RETURN
    end if
    if(DoInit)call init_sort_stencil
    !\
    ! Failure in constructing the extended stencil may  occur if the 
    ! input point is out of AMR grid.
    !/
    if(nExtendedStencil < 2**nDim)then

       nGridOut = -1
       RETURN
    end if
    select case(nDim)
    case(2)
       call fix_basic_stencil2(&
            Xyz_D, DxyzInv_D, XyzGrid_DII(:,0,:), iLevelSubgrid_I, &
            XyzStencil_D, iCaseExtended)
       call generate_basic_stencil(&
            nDim, XyzStencil_D, nExtendedStencil,                      &
            XyzExtended_DI(:,1:nExtendedStencil), DxyzInv_D, iOrderExtended_I)
       if(any(iOrderExtended_I < 1))&
            call CON_stop('Failure in constructing basic stencil') 
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
       DoStencilFix = .true.  
       call fix_basic_stencil3(&
            Xyz_D, DxyzInv_D, XyzGrid_DII(:,0,:), iLevelSubgrid_I, &
            XyzStencil_D, iCaseExtended)
       do while(DoStencilFix)
          call generate_basic_stencil(&
               nDim, XyzStencil_D, nExtendedStencil,              &
               XyzExtended_DI(:,1:nExtendedStencil),              &
               DxyzInv_D, iOrderExtended_I)
          if(any(iOrderExtended_I < 1))&
               call CON_stop('Failure in constructing basic stencil') 
          XyzGrid_DI = XyzExtended_DI(:,iOrderExtended_I)
          iLevel_I = iLevelExtended_I(iOrderExtended_I)
          iIndexes_II = &
               iIndexesExtended_II(:,iOrderExtended_I)
          IsOut_I = IsOutExtended_I(iOrderExtended_I)
          call interpolate_amr3(&
               Xyz_D , XyzGrid_DI, iLevel_I, IsOut_I, iCaseExtended,&
               nGridOut, Weight_I, iOrder_I, DoStencilFix)
          if(DoStencilFix)then
             iCaseExtended = i_case(iLevelSubGrid_I)
             XyzStencil_D = (XyzGrid_DII(:,0,1) + XyzGrid_DII(:,0,8))/2
          end if
       end do
       if(nGridOut < 1)then
          write(*,*)'Xyz_D=',Xyz_D
          write(*,*)'XyzGrid_DI=',XyzGrid_DI
          write(*,*)'iCaseExtended=',iCaseExtended
          write(*,*)'iSortStencil=',iSortStencil3_II(:,iCaseExtended)
          call CON_stop('Failure in interpolate amr grid3') 
       end if
    case default
       call CON_stop('Only 2D and 3D AMR grids are implemented')
    end select
    call sort_out
    iIndexes_II(:, 1:nGridOut) = iIndexes_II(:,iOrder_I(1:nGridOut))
    if(present(IsSecondOrder))&
         IsSecondOrder = .not.any(IsOut_I)
  contains
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
    !=====================
    subroutine generate_extended_stencil
      !\
      !This routine extracts the grid points which in principle 
      !can be involved into interpolation
      !/
      !\
      ! Loop variables
      !/
      integer:: iGrid, iSubgrid, iDir 
      !\
      ! To evaluate grid points not belonging to the basic block
      !/
      integer:: iGridStored, iGridCheck
      !\
      ! For using find routine with inout argument
      !/
      real  :: XyzMisc_D(nDim), XyzStored_D(nDim)
      !\
      ! Output parameters of find routine
      !/
      real    :: XyzCorner_D(nDim), Dxyz_D(nDim)
      logical :: IsOut
      !\
      !Displacement measured in grid sizes or in their halfs
      !/ 
      integer, dimension(nDim) :: iShift_D
      !\
      ! Shift of the iGrid point in the stencil with respect to the
      ! first one
      !/
      integer, dimension(3,8),parameter :: iShift_DI = reshape((/&
           0, 0, 0,   1, 0, 0,   0, 1, 0,   1, 1, 0, &
           0, 0, 1,   1, 0, 1,   0, 1, 1,   1, 1, 1/),(/3,8/))
      !\
      ! Find block to which the point belong
      !/ 
      call find(nDim, Xyz_D, &
           iProc_I(1), iBlock_I(1),XyzCorner_D, Dxyz_D, IsOut)
      !\
      ! Now Xyz_D is given  with respect to the block corner
      !/
      if(IsOut)then
         !\
         ! The algorithm does not work for a point out of the computation
         ! domain. It could, but in this case too much information
         ! about grid should be brought to the table - the domain size,
         ! periodicity etc
         !/
         nExtendedStencil = -1
         RETURN
      end if
      !\
      ! Initialize iLevelSubgrid_I as all coarser, to enter the loop
      !/
      iLevelSubgrid_I = -1

      COARSEN:do while(any(iLevelSubgrid_I==-1) )
         !\
         !This loop works more than once if the stencil includes
         !blocks which are coarser than the first one found. In this 
         !case the coarser block is used as the base for a stencil
         !/
         iLevelSubgrid_I = 0
         nSubGrid_I      = 1

         iCellIndexes_DII = 0
         XyzGrid_DII      = 0
         DxyzInv_D = 1/Dxyz_D
         XyzStored_D = XyzCorner_D
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
         !if(all(abs(XyzMisc_D) < cTol2))then
            !\
            ! Xyz coincides with the grid point
            ! Commented out to satisfy iFort
            !/
            !$ nGridOut = 1; Weight_I =0; Weight_I(1) = 1
            !$ iIndexes_II(0,       1) = iProc_I(1)
            !$ iIndexes_II(nIndexes,1) = iBlock_I(1)
            !$ iIndexes_II(1:nDim,  1) = iCellIndexes_DII(:,1,1)
            !$ RETURN
         !elseif(all(iCellIndexes_DII(:,1,1) > 0).and.&
         !     all(  iCellIndexes_DII(:,1,1) < nCell_D))then
            !\
            ! The whole interpolation stencil is within this block
            !/
        !    call interpolate_uniform( XyzMisc_D)
            !\
            ! Form index array and sort out zero weights
            !/
        !   nGridOut = 0 ; cTol2 = 2*cTol2
        !     do iGrid = 1, nGrid
        !    if(Weight_I(iGrid) < cTol2)CYCLE
            !   nGridOut = nGridOut + 1
             !  iIndexes_II(0,       nGridOut) = iProc_I(1)
              ! iIndexes_II(nIndexes,nGridOut) = iBlock_I(1)
               !iIndexes_II(1:nDim,  nGridOut) = iCellIndexes_DII(:,1,1) + &
       !            iShift_DI(1:nDim,iGrid)
       !        Weight_I(nGridOut) = Weight_I(iGrid)
       !     end do
       !     RETURN
       !  end if
         
         XyzGrid_DII(:,0,1) = Dxyz_D*(iCellIndexes_DII(:,1,1) - 0.50)
         !\ 
         ! now XyzMisc_D = (Xyz_D-XyzGrid_DII(:,0,1))/Dxyz satisfies 
         ! inequalities: XyzMisc_D >= 0 and XyzMisc_D < 1. Strengthen 
         ! these inequalities 
         !/
         XyzMisc_D = min(1 - cTol2,&
              max(XyzMisc_D, cTol2 ))
         Xyz_D = XyzGrid_DII(:,0,1) + XyzMisc_D*Dxyz_D
         !\
         !Calculate other grid points, check if all points belong to 
         !the found block
         !/
         iGridCheck = -1
         do iGrid = nGrid, 1, -1
            iShift_D = iShift_DI(1:nDim,iGrid)
            iBlock_I(iGrid) = iBlock_I(1)
            iCellIndexes_DII(:,1,iGrid) = &
                 iCellIndexes_DII(:,1,1) + iShift_D
            XyzGrid_DII(:,0,iGrid) = &
                 XyzGrid_DII(:,0,1) + iShift_D*Dxyz_D
            if(any(iCellIndexes_DII(:,1,iGrid) < 1).or.&
                 any(iCellIndexes_DII(:,1,iGrid) > nCell_D))then
               !\
               !This grid point is out of block, mark it
               !/
               iProc_I(iGrid) = -1 
               if(iGridCheck==-1)iGridCheck = iGrid
            else
               XyzGrid_DII(:,1,iGrid) = XyzGrid_DII(:,0,iGrid)
               nSubGrid_I(iGrid) = 1
               iProc_I(iGrid) = iProc_I(1)
            end if
         end do
         iGridStored = iGridCheck
         !\
         ! For the grid points not belonging to the block
         ! find the block they belong to
         !/ 
         NEIBLOCK:do while(iGridStored/= -1)
            !\
            !Recalculate absolute coordinates for
            !the grid point which is out of the block
            !/
            XyzMisc_D = XyzStored_D + XyzGrid_DII(:,0,iGridStored)
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
               nSubGrid_I(iGridStored)      = 1
               !\
               ! Find the next out-of-block point
               !/
               do iGrid = iGridStored - 1, 1, -1
                  if(iProc_I(iGrid)==-1)then
                     iGridStored = iGrid
                     CYCLE NEIBLOCK
                  end if
               end do
               EXIT COARSEN
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
               CYCLE COARSEN
            case(0)  ! (New Dxyz_D)*Stored DXyzInv =1
               XyzGrid_DII(:,1,iGridStored) = XyzGrid_DII(:,0,iGridStored)
               !\
               ! Calculate cell indexes as we did before. Use nint
               ! instead of int as long as XyzMisc_D is very close to 
               ! the grid point 
               !/
               iCellIndexes_DII(:,1,iGridStored) = &
                    nint(XyzMisc_D*DxyzInv_D + 0.50)
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
               iGridCheck = -1
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
                     if(iGridCheck==-1)iGridCheck = iGrid
                  else
                     iProc_I( iGrid) = iProc_I( iGridStored)
                     iBlock_I(iGrid) = iBlock_I(iGridStored)
                     XyzGrid_DII(:,1,iGrid) = XyzGrid_DII(:,0,iGrid)
                     nSubGrid_I(     iGrid) = 1
                     iLevelSubGrid_I(iGrid) = 0
                  end if
               end do
            case(Fine_  )  !1, (New Dxyz_D)*Stored DXyzInv =0.5
               IsExtended = .true. 
               !\
               !The number of points which may potentially
               !/
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
                    floor(2*XyzMisc_D*DxyzInv_D + 0.50)
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
               iGridCheck = -1
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
                     if(iGridCheck==-1)iGridCheck = iGrid
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
            end select
            iGridStored = iGridCheck
         end do NEIBLOCK
      end do COARSEN
      !\
      ! Handle points behind the boundary
      !/
      if(UseGhostPoint)then
         select case(count(IsOut_I))
         case(3, 7) !3,7 are nGrid -1 for nDim=2,3 
            !\
            !Find the only physical point in the stencil 
            !/
            iGridCheck = maxloc(iLevelSubGrid_I,MASK=.not.IsOut_I, DIM=1)
            nGridOut = 1
            Weight_I = 0; Weight_I(1) = 1
            iIndexes_II(0,       1) = iProc_I( iGridCheck) 
            iIndexes_II(nIndexes,1) = iBlock_I(iGridCheck)
            iIndexes_II(1:nDim,1  ) = iCellIndexes_DII(:,1,iGridCheck)
            RETURN
         case(4)
            ! nDim = 3
            !\
            ! One of the faces is fully out of the domain
            !/
            !\
            !Find point in the domain
            !/
            iGridCheck = maxloc(iLevelSubGrid_I,MASK=.not.IsOut_I,DIM=1)
            do iDir = 1, nDim
               if(all(&
                    IsOut_I(iOppositeFace_IDI(:,iDir,iGridCheck))))then
                  if(IsExtended)then
                     !\
                     ! Prolong grid from the physical face to the ghost one
                     ! accounting for the difference in resolution
                     !/
                     do iGrid = 1,4
                        call prolong(iFace_IDI(iGrid,iDir,iGridCheck), &
                             iOppositeFace_IDI(iGrid,iDir,iGridCheck))
                     end do
                  else
                     !\
                     !Grid near the boundary is uniform. Put the point to the
                     !physical face to nullify the ghost cell contributions
                     !/
                     Xyz_D(iDir) = XyzGrid_DII(iDir,1,iGridCheck)
                  end if
                  EXIT
               end if
            end do
         case(2)
            ! Analogous to the previous case, but nDim = 2 
            !\
            !Find point in the domain
            !/
            iGridCheck = maxloc(iLevelSubGrid_I,MASK=.not.IsOut_I, DIM=1)
            do iDir = 1, 2
               if(all(&
                    IsOut_I(iOppositeSide_IDI(:,iDir,iGridCheck))))then
                  if(IsExtended)then
                     !\
                     ! Prolong grid from the physical face to the ghost 
                     ! accounting for the difference in resolution
                     !/
                     do iGrid = 1,2
                        call prolong(iSide_IDI(iGrid,iDir,iGridCheck), &
                             iOppositeSide_IDI(iGrid,iDir,iGridCheck))
                     end do
                  else
                     !\
                     !Put the point to the physical phase, which assign them
                     !zero weight
                     !/
                     Xyz_D(3 - iDir) = XyzGrid_DII(3 - iDir,1,iGridCheck)
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
            iGridCheck = maxloc(iLevelSubGrid_I,MASK=.not.IsOut_I, DIM=1)
            do iDir  = 1, nDim
               if(.not.IsOut_I(iEdge_ID(iGridCheck,iDir)))then
                  if(IsExtended)then
                     do iGrid = 2,4
                        call prolong(iGridCheck,&
                             iFace_IDI(iGrid,iDir,iGridCheck))
                        call prolong(iEdge_ID(iGridCheck,iDir),&
                             iOppositeFace_IDI(iGrid,iDir,iGridCheck))
                     end do
                     EXIT
                  else
                     !\
                     ! Put the point to the physical edge
                     !/
                     Xyz_D(1 + mod(iDir,3)) = &
                          XyzGrid_DII(1 + mod(iDir,3),1,iGridCheck)
                     Xyz_D(1 + mod(iDir+1,3)) = &
                          XyzGrid_DII(1 + mod(iDir+1,3),1,iGridCheck)
                  end if
               end if
            end do
         end select
      end if
      !\
      ! Finalize:
      !/
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
    !======================
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
       nGridOut, Weight_I, iOrder_I,&
       DoStencilFix)

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
    logical, intent(out):: DoStencilFix
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
    ! To call subroutine pyramid
    !/
    integer, parameter:: Rectangular_=1, Trapezoidal_=2

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
    DoStencilFix = .false. 
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
       DoStencilFix = nGridOut <= 3
    case(Transition2Corner_)
       call interpolate_corner_transition(iDir)
    case(FiveTetrahedra_)
       !\
       ! One configuration resulting in the corner split for five tetrahedra
       ! is occured in many cases (26), however, it is treated in the same  
       ! way in all these cases. If the given grid point is connected by 
       ! three edges with three coarse points, while the across the main  
       ! diagonal point is connected by three edges (or, equivalently, the  
       ! given point is connected by three face diagonals) with three coarse  
       ! points), then all faces are trangulated   
       !/
       !\
       ! View from  point iOrder
       !            xy face diag
       !   x-axis      |      y-axis
       !         \     V     /
       !          \  / F    /           
       !           C   | \C           
       !           |\  | / |            
       !           | \ |/  |            
       !           |  \/   |
       !           |   1   |
       !  xz facediag/ | \yz face diag
       !           /   |  \|
       !           F   |   F
       !             \ C /
       !               |
       !               V
       !               z-axis

       call pyramids(&
            iTetrahedron1_I=(/1, 2, 6, 4/),&
            iTetrahedron2_I=(/1, 3, 7, 4/),&
            iTetrahedron3_I=(/1, 5, 7, 6/),&
            iTetrahedron4_I=(/1, 7, 6, 4/),&
            iTetrahedron5_I=(/8, 7, 6, 4/))
    case(OneFine_)    
       !\
       ! View from the fine point
       !            xy face diag
       !   x-axis      |      y-axis
       !         \     V     /
       !          \   C4    /      Connections to opposite faces: main diag+
       !           C2/ | \C3       3478 - y edge, yz face diag, xy face daig
       !           |\  | / |       2468 - x edge, xy face diag, xz face diag
       !           | \ |/  |       5687 - z face. yz face diag, xy face diag
       !           |  \/   |
       !           |   F1  |
       !  xz facediag/ | \yz face diag
       !           /   |  \|
       !          C6   |   C7
       !             \ C5/
       !               |
       !               V
       !               z-axis
       call pyramids(&
            iRectangular1_I=(/2, 4, 6, 8, 1/) ,&
            iRectangular2_I=(/3, 4, 7, 8, 1/)  ,&
            iRectangular3_I=(/5, 6, 7, 8, 1/))
    case(OneCoarse_)                  ! 8 cases , totally 79 left 176
       !\
       ! View from the coarse point
       !            xy face diag
       !   x-axis            y-axis
       !         \          /
       !          \    F4  /           
       !           F2----- F3            
       !           |\    / |            
       !           | \  / /|            
       !           |\ \/   |            
       !           |   C1  |
       !           | \ | / |
       !           |   |   |
       !          F6  \ / F7
       !             \ F5/
       !               |
       !               V
       !               z-axis
       call pyramids(iTetrahedron1_I=(/1,2,3,5/))
       if(nGridOut>1)RETURN
       call parallel_rays(&
            Dir_D=XyzGrid_DI(:,iOrder_I(8)) - XyzGrid_DI(:,iOrder_I(1)), &
            iURectangle1_I=(/2, 4, 6, 8/),&
            iURectangle2_I=(/3, 4, 7, 8/),&
            iURectangle3_I=(/5, 6, 7, 8/),&
            iDTriangle1_I=(/2, 3, 5/),    &
            iDTriangle2_I=(/2, 3, 4/),    &
            iDTriangle3_I=(/2, 6, 5/),    &
            iDTriangle4_I=(/3, 5, 7/))
    case(FineMainDiag_)
       !Full traingulation of all faces, each coarse vertex is connected  
       !by edge with one fine point and with face diagonal to the other. 
       !\
       ! View from the fine point
       !            xy face diag
       !   x-axis      |      y-axis
       !         \     V     /
       !          \   C4    /           
       !           C2/ | \C3           The point F8 is also fine 
       !           |\  | / |           The line of sight goes along the
       !           | \ |/  |           main diagonal F1F8 which is the 
       !           |  \/   |           common side  of 6 tetrahedron
       !           |   F1  |
       !  xz facediag/ | \yz face diag
       !           /   |  \|
       !          C6   |   C7
       !             \ C5/
       !               |
       !               V
       !               z-axis
       call pyramids(&
            iTetrahedron1_I=(/1, 2, 4, 8/), &
            iTetrahedron2_I=(/1, 2, 6, 8/), &
            iTetrahedron3_I=(/1, 3, 4, 8/), &
            iTetrahedron4_I=(/1, 3, 7, 8/), &
            iTetrahedron5_I=(/1, 5, 6, 8/), &
            iTetrahedron6_I=(/1, 5, 7, 8/))
    case(FineFaceDiag_)
       !
       !       ^
       !iDir-axis
       !       |  7C---------8C
       !       | /|         /|
       !       |/ |        / |
       !       5C---------6C |
       !       |  |   _   |  |
       !       |  |   /|1+mod(iDir+1,3) axis
       !       |  |  /    |  |
       !       |  | /     |  |
       !       |  |/      |  |
       !       |  3C---------4F
       !       | /        | /
       !       |/         |/
       !       1F---------2C-----> 1+mod(iDir,3)-axis
       !
       !
       ! 1 Remove tetrahedra 1F 4F 2C 6C and 1F 4F 3C 7C
       call pyramids(&
            iTetrahedron1_I=(/1, 4, 2, 6/), &
            iTetrahedron2_I=(/1, 4, 3, 7/))
       if(nGridOut > 1)RETURN
       !
       !    7C------ 8C 
       !  5C--\---6C/|     
       !   |  _ \_  \|
       !   | /    --4F
       !  1F----/     
       !\
       ! Center of the lower face
       !/
       XyzMin_D = &
            0.50*(XyzGrid_DI(:,iOrder_I(1)) + XyzGrid_DI(:,iOrder_I(4)))
       !\
       ! Center of the upper face
       !/
       XyzMax_D = 0.50*&
            (XyzGrid_DI(:,iOrder_I(5))+XyzGrid_DI(:,iOrder_I(8)))
       call parallel_rays(Dir_D=XyzMax_D - XyzMin_D,&
            iURectangle1_I=(/5, 6, 7, 8/),&
            iDTriangle1_I=(/1, 4, 6/),&
            iDTriangle2_I=(/1, 4, 7/),&
            iDTriangle3_I=(/4, 6, 8/),&
            iDTriangle4_I=(/4, 7, 8/),&
            iDTriangle5_I=(/1, 5, 7/),&
            iDTriangle6_I=(/1, 5, 6/))
    case(CoarseMainDiag_)  
       !Two coarse pints across the main diagonal
       !\
       ! View from the coarse point
       !            xy face diag
       !   x-axis            y-axis
       !         \          /
       !          \    F4  /           
       !           F2----- F3            
       !           |\    / |            
       !           | \  / /|            
       !           |\ \/   |            
       !           |   C1  |
       !           | \ | / |
       !           |   |   |
       !          F6  \ / F7
       !             \ F5/
       !               |
       !               V
       !               z-axis
       call pyramids(&
            iTetrahedron1_I=(/1, 2, 3, 5/),&
            iTetrahedron2_I=(/8, 7, 6, 4/))
       if(nGridOut>1)RETURN
       call parallel_rays(&
            Dir_D=XyzGrid_DI(:,iOrder_I(8)) - XyzGrid_DI(:,iOrder_I(1)), &
            iUTriangle1_I=(/7, 6, 4/),&
            iUTriangle2_I=(/7, 6, 5/),&
            iUTriangle3_I=(/7, 4, 3/),&
            iUTriangle4_I=(/6, 4, 2/),&
            iDTriangle1_I=(/2, 3, 5/),&
            iDTriangle2_I=(/2, 3, 4/),&
            iDTriangle3_I=(/2, 5, 6/),&
            iDTriangle4_I=(/3, 5, 7/))
    case(CoarseFaceDiag_)
       !       ^
       !iDir-axis
       !       |  7F---------8F
       !       | /|         /|
       !       |/ |        / |
       !       5F---------6F |
       !       |  |   _   |  |
       !       |  |   /|1+mod(iDir+1,3) axis
       !       |  |  /    |  |
       !       |  | /     |  |
       !       |  |/      |  |
       !       |  3F---------4C
       !       | /        | /
       !       |/         |/
       !       1C---------2F-----> 1+mod(iDir,3)-axis
       !
       !
       ! 1 Remove tetrahedra 1C 2F 3F 5F and 4C 2F 3F 8F
       call pyramids(&
            iTetrahedron1_I=(/1, 2, 3, 5/),&
            iTetrahedron2_I=(/4, 2, 3, 8/))
       if(nGridOut < 1)then
          !
          !    7F------ 8F 
          !  5F------6F//     
          !   \|      | |
          !    3F_    |/
          !        \_ 2F
          !\
          ! Center of the lower face
          !/
          XyzMin_D = 0.50*&
               (XyzGrid_DI(:,iOrder_I(1))+XyzGrid_DI(:,iOrder_I(4)))
          !\
          ! Center of the upper face
          !/
          XyzMax_D = 0.50*&
               (XyzGrid_DI(:,iOrder_I(5)) + XyzGrid_DI(:,iOrder_I(8)))
          call parallel_rays(Dir_D=XyzMax_D - XyzMin_D,&
               iURectangle1_I=(/5, 6, 7, 8/),&
               iDTriangle1_I=(/2, 3, 5/),&
               iDTriangle2_I=(/2, 3, 8/))
       end if
    case(FineEdgePlusOne_)
       !       ^
       !iDir-axis
       !       |  7C---------8F
       !       | /|         /|
       !       |/ |        / |
       !       5C---------6C |
       !       |  |   _   |  |
       !       |  |   /|1+mod(iDir+1,3) axis
       !       |  |  /    |  |
       !       |  | /     |  |
       !       |  |/      |  |
       !       |  3C---------4F
       !       | /        | /
       !       |/         |/
       !       1F---------2C-----> 1+mod(iDir,3)-axis
       !
       !                    F4 F8
       !         Trapezoid: C3 C7
       !                
       !         Trapezoid: F4 F8
       !                    C2 C6
       !     F1 is apex of trapezoidal pyramids 
       !     Tetrahedra left F1C5F8C6 + F1C5F8C7
       call pyramids(&
            iTetrahedron1_I=(/1, 5, 8, 6/),&
            iTetrahedron2_I=(/1, 5, 8, 7/),&
            iTrapezoidal1_I=(/2, 6, 4, 8, 1/),&
            iTrapezoidal2_I=(/3, 7, 4, 8, 1/))
    case(ThreeFineOnFace_)
       !
       !       ^
       !iDir-axis
       !       |  7C---------8C
       !       | /|         /|
       !       |/ |        / |
       !       5C---------6C |
       !       |  |   _   |  |
       !       |  |   /|1+mod(iDir+1,3) axis
       !       |  |  /    |  |
       !       |  | /     |  |
       !       |  |/      |  |
       !       |  3F---------4C
       !       | /        | /
       !       |/         |/
       !       1F---------2F-----> 1+mod(iDir,3)-axis
       !
       !                  F1 F3
       !       Trapezoid: C5 C7       
       !               
       !       Trapezoid: F1 F2
       !                  C5 C6

       call pyramids(&
            iTetrahedron1_I=(/2, 3, 8, 4/))
       if(nGridOut> 1)RETURN
       !\
       !Interpolation in parralel rays
       !/
       !The view for the lower subfaces is as follows
       !  C5----------C6
       !  | \         |
       !  |   F1---F2 |
       !  |   |   /   |
       !  !   |  /  \ | 
       !  |   F3      |
       !  | /     \   |
       !  C7----------C8      
       call parallel_rays(Dir_D=&
            XyzGrid_DI(:,iOrder_I(8)) - &
            XyzGrid_DI(:,iOrder_I(4)),&
            iDTriangle1_I=(/1, 2, 3/),&
            iDTriangle2_I=(/8, 2, 3/),&
            iDTriangle3_I=(/8, 7, 3/),&
            iDTriangle4_I=(/8, 2, 6/),&
            iURectangle1_I=(/5, 6, 7, 8/))    
    case(CoarseEdgePlusOne_)
       !
       !       ^
       !iDir-axis
       !       |  7F---------8C
       !       | /|         /|
       !       |/ |        / |
       !       5F---------6F |
       !       |  |   _   |  |
       !       |  |   /|1+mod(iDir+1,3) axis
       !       |  |  /    |  |
       !       |  | /     |  |
       !       |  |/      |  |
       !       |  3F---------4C
       !       | /        | /
       !       |/         |/
       !       1C---------2F-----> 1+mod(iDir,3)-axis
       !
       !
       !                  F3 F7
       !       Trapezoid: C4 C8       
       !         
       !       Trapezoid: F2 F6
       !                  C4 C8
       !                
       call pyramids(&
            iTetrahedron1_I=(/1, 2, 3, 5/),&
            iRectangular1_I=(/2, 3, 6, 7, 5/))
       if(nGridOut> 1)RETURN
       call parallel_rays(Dir_D=&
            XyzGrid_DI(:,iOrder_I(3)) - &
            XyzGrid_DI(:,iOrder_I(2)),&
            iDTrapezoid1_I=(/4, 8, 2, 6/),&
            iUTrapezoid1_I=(/4, 8, 3, 7/)) 
    case(ThreeCoarseOnFace_)
       !
       !       ^
       !iDir-axis
       !       |  7F---------8F
       !       | /|         /|
       !       |/ |        / |
       !       5F---------6F |
       !       |  |   _   |  |
       !       |  |   /|1+mod(iDir+1,3) axis
       !       |  |  /    |  |
       !       |  | /     |  |
       !       |  |/      |  |
       !       |  3C---------4F
       !       | /        | /
       !       |/         |/
       !       1C---------2C-----> 1+mod(iDir,3)-axis
       !
       !
       !            F5 F7
       ! Trapezoid: C1 C3       
       !           
       ! Trapezoid: F5 F6
       !            C1 C2
       ! Rectangle: F5 F6 F7 F8
       !/               
       !    Common apex F4
       call pyramids(&
            iRectangular1_I=(/7, 8, 5, 6, 4/),&
            iTrapezoidal1_I=(/1, 2, 5, 6, 4/),&
            iTrapezoidal2_I=(/1, 3, 5, 7, 4/))
    case(ThreeCoarseOnFacePlusOne_)
       !\
       ! face 1:4 has only one Coarse point 1
       ! face 5:8 has three Coarse points iOrder((/6,7,8/))
       !/
       !-----------------------------------------------
       !
       !       ^
       !iDir-axis
       !       |  7C---------8C
       !       | /|         /|
       !       |/ |        / |
       !       5F---------6C |
       !       |  |   _   |  |
       !       |  |   /|1+mod(iDir+1,3) axis
       !       |  |  /    |  |
       !       |  | /     |  |
       !       |  |/      |  |
       !       |  3F---------4F
       !       | /        | /
       !       |/         |/
       !       1C---------2F-----> 1+mod(iDir,3)-axis
       !
       !                F2 F4
       !      Trapezoid C6 C8
       !                                           
       !                F3 F4
       !      Trapezoid C7 C8
       !                
       !
       !Common apex F5
       call pyramids(&
            iTrapezoidal1_I=(/6, 8, 2, 4, 5/),&
            iTrapezoidal2_I=(/7, 8, 3, 4, 5/),&
            iTetrahedron1_I=(/2, 3, 4, 5/),   &
            iTetrahedron2_I=(/1, 2, 3, 5/) )
    case(CoarseChain_)
       !\
       !Chain of four coarse points connected with mutually 
       !orthogonal edges
       !/
       !-----------------------------------------------
       !
       !       ^
       !iDir-axis
       !       |  C7---------8F
       !       | /|         /|
       !       |/ |        / |
       !       5C---------6F |
       !       |  |   _   |  |
       !       |  |   /|1+mod(iDir+1,3) axis
       !       |  |  /    |  |
       !       |  | /     |  |
       !       |  |/      |  |
       !       |  3F---------4F
       !       | /        | /
       !       |/         |/
       !       1C---------2C-----> 1+mod(iDir,3)-axis
       !
       !                  F6 F8
       !        Trapezoid C5 C7
       !                                           
       !      
       !                  F3 F4
       !        Trapezoid C1 C2                
       !
       call pyramids(&
            iTrapezoidal1_I=(/1, 2, 3, 4, 6/),& 
            iTetrahedron1_I=(/8, 3, 4, 6/),&  !see #1 below 
            iTetrahedron2_I=(/1, 6, 3, 5/),&  !see #2 below
            iTrapezoidal2_I=(/5, 7, 6, 8, 3/))!see #3 below
       ! #1 The leftover is above subfaces 136 and 364
       ! #2 The leftover is above subfaces 136 and 368
       ! #3 The leftover is above subfaces 365 and 368
    end select
  contains
    !==========================
    subroutine pyramids(&
         iTetrahedron1_I,iTetrahedron2_I,iTetrahedron3_I,&
         iTetrahedron4_I,iTetrahedron5_I,iTetrahedron6_I,&
         iRectangular1_I,iRectangular2_I,iRectangular3_I,&
         iTrapezoidal1_I,iTrapezoidal2_I)
      integer,intent(in),optional,dimension(4)::&
           iTetrahedron1_I, iTetrahedron2_I, iTetrahedron3_I,&
           iTetrahedron4_I, iTetrahedron5_I, iTetrahedron6_I
      integer,intent(in),optional,dimension(5)::&
           iRectangular1_I, iRectangular2_I, iRectangular3_I,&
           iTrapezoidal1_I, iTrapezoidal2_I
      !---------------------------
      nGridOut = -1; Weight_I = -1
      if(present(iTetrahedron1_I))then
         call tetrahedron(&
              XyzGrid_DI(:,iOrder_I(iTetrahedron1_I(1))),&
              XyzGrid_DI(:,iOrder_I(iTetrahedron1_I(2))),&
              XyzGrid_DI(:,iOrder_I(iTetrahedron1_I(3))),&
              XyzGrid_DI(:,iOrder_I(iTetrahedron1_I(4))))
         if(all(Weight_I(1:4).ge.0.0))then
            nGridOut = 4
            iOrder_I(1:4) = iOrder_I(iTetrahedron1_I)
            RETURN
         elseif(present(iTetrahedron2_I))then
            call tetrahedron(&
                 XyzGrid_DI(:,iOrder_I(iTetrahedron2_I(1))),&
                 XyzGrid_DI(:,iOrder_I(iTetrahedron2_I(2))),&
                 XyzGrid_DI(:,iOrder_I(iTetrahedron2_I(3))),&
                 XyzGrid_DI(:,iOrder_I(iTetrahedron2_I(4))))
            if(all(Weight_I(1:4).ge.0.0))then
               nGridOut = 4
               iOrder_I(1:4) = iOrder_I(iTetrahedron2_I)
               RETURN
            elseif(present(iTetrahedron3_I))then
               call tetrahedron(&
                    XyzGrid_DI(:,iOrder_I(iTetrahedron3_I(1))),&
                    XyzGrid_DI(:,iOrder_I(iTetrahedron3_I(2))),&
                    XyzGrid_DI(:,iOrder_I(iTetrahedron3_I(3))),&
                    XyzGrid_DI(:,iOrder_I(iTetrahedron3_I(4))))
               if(all(Weight_I(1:4).ge.0.0))then
                  nGridOut = 4
                  iOrder_I(1:4) = iOrder_I(iTetrahedron3_I)
                  RETURN
               elseif(present(iTetrahedron4_I))then
                  call tetrahedron(&
                       XyzGrid_DI(:,iOrder_I(iTetrahedron4_I(1))),&
                       XyzGrid_DI(:,iOrder_I(iTetrahedron4_I(2))),&
                       XyzGrid_DI(:,iOrder_I(iTetrahedron4_I(3))),&
                       XyzGrid_DI(:,iOrder_I(iTetrahedron4_I(4))))
                  if(all(Weight_I(1:4).ge.0.0))then
                     nGridOut = 4
                     iOrder_I(1:4) = iOrder_I(iTetrahedron4_I)
                     RETURN
                  elseif(present(iTetrahedron5_I))then
                     call tetrahedron(&
                          XyzGrid_DI(:,iOrder_I(iTetrahedron5_I(1))),&
                          XyzGrid_DI(:,iOrder_I(iTetrahedron5_I(2))),&
                          XyzGrid_DI(:,iOrder_I(iTetrahedron5_I(3))),&
                          XyzGrid_DI(:,iOrder_I(iTetrahedron5_I(4))))
                     if(all(Weight_I(1:4).ge.0.0))then
                        nGridOut = 4
                        iOrder_I(1:4) = iOrder_I(iTetrahedron5_I)
                        RETURN
                     elseif(present(iTetrahedron6_I))then
                        call tetrahedron(&
                             XyzGrid_DI(:,iOrder_I(iTetrahedron6_I(1))),&
                             XyzGrid_DI(:,iOrder_I(iTetrahedron6_I(2))),&
                             XyzGrid_DI(:,iOrder_I(iTetrahedron6_I(3))),&
                             XyzGrid_DI(:,iOrder_I(iTetrahedron6_I(4))))
                        if(all(Weight_I(1:4).ge.0.0))then
                           nGridOut = 4
                           iOrder_I(1:4) = iOrder_I(iTetrahedron6_I)
                           RETURN
                        end if ! 6
                     end if    ! 5
                  end if       ! 4
               end if          ! 3
            end if             ! 2
         end if                ! 1
      end if                   ! no tetrahedron
      if(present(iRectangular1_I))then
         call pyramid(iBase=Rectangular_,&
              X1_D=XyzGrid_DI(:,iOrder_I(iRectangular1_I(1))),&
              X2_D=XyzGrid_DI(:,iOrder_I(iRectangular1_I(2))),&
              X3_D=XyzGrid_DI(:,iOrder_I(iRectangular1_I(3))),&
              X4_D=XyzGrid_DI(:,iOrder_I(iRectangular1_I(4))),&
              X5_D=XyzGrid_DI(:,iOrder_I(iRectangular1_I(5))) )
         if(all(Weight_I(1:5).ge.0.0))then
            nGridOut = 5
            iOrder_I(1:5) = iOrder_I(iRectangular1_I)
            RETURN
         elseif(present(iRectangular2_I))then
            call pyramid(iBase=Rectangular_,&
                 X1_D=XyzGrid_DI(:,iOrder_I(iRectangular2_I(1))),&
                 X2_D=XyzGrid_DI(:,iOrder_I(iRectangular2_I(2))),&
                 X3_D=XyzGrid_DI(:,iOrder_I(iRectangular2_I(3))),&
                 X4_D=XyzGrid_DI(:,iOrder_I(iRectangular2_I(4))),&
                 X5_D=XyzGrid_DI(:,iOrder_I(iRectangular2_I(5))) )
            if(all(Weight_I(1:5).ge.0.0))then
               nGridOut = 5
               iOrder_I(1:5) = iOrder_I(iRectangular2_I)
               RETURN
            elseif(present(iRectangular3_I))then
               call pyramid(iBase=Rectangular_,&
                    X1_D=XyzGrid_DI(:,iOrder_I(iRectangular3_I(1))),&
                    X2_D=XyzGrid_DI(:,iOrder_I(iRectangular3_I(2))),&
                    X3_D=XyzGrid_DI(:,iOrder_I(iRectangular3_I(3))),&
                    X4_D=XyzGrid_DI(:,iOrder_I(iRectangular3_I(4))),&
                    X5_D=XyzGrid_DI(:,iOrder_I(iRectangular3_I(5))) )
               if(all(Weight_I(1:5).ge.0.0))then
                  nGridOut = 5
                  iOrder_I(1:5) = iOrder_I(iRectangular3_I)
                  RETURN
               end if          ! 3
            end if             ! 2
         end if                ! 1
      end if                   ! no rectangular
      if(present(iTrapezoidal1_I))then
         call pyramid(iBase=Trapezoidal_,&
              X1_D=XyzGrid_DI(:,iOrder_I(iTrapezoidal1_I(1))),&
              X2_D=XyzGrid_DI(:,iOrder_I(iTrapezoidal1_I(2))),&
              X3_D=XyzGrid_DI(:,iOrder_I(iTrapezoidal1_I(3))),&
              X4_D=XyzGrid_DI(:,iOrder_I(iTrapezoidal1_I(4))),&
              X5_D=XyzGrid_DI(:,iOrder_I(iTrapezoidal1_I(5))) )
         if(all(Weight_I(1:5).ge.0.0))then
            nGridOut = 5
            iOrder_I(1:5) = iOrder_I(iTrapezoidal1_I)
            RETURN
         elseif(present(iTrapezoidal2_I))then
            call pyramid(iBase=Trapezoidal_,&
                 X1_D=XyzGrid_DI(:,iOrder_I(iTrapezoidal2_I(1))),&
                 X2_D=XyzGrid_DI(:,iOrder_I(iTrapezoidal2_I(2))),&
                 X3_D=XyzGrid_DI(:,iOrder_I(iTrapezoidal2_I(3))),&
                 X4_D=XyzGrid_DI(:,iOrder_I(iTrapezoidal2_I(4))),&
                 X5_D=XyzGrid_DI(:,iOrder_I(iTrapezoidal2_I(5))) )
            if(all(Weight_I(1:5).ge.0.0))then
               nGridOut = 5
               iOrder_I(1:5) = iOrder_I(iTrapezoidal2_I)
               RETURN
            end if             ! 2
         end if                ! 1
      end if                   ! no trapezoid
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
      real:: AlphaUp, AlphaDown
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

      nGridOut = -1; nGridOutUp= -1; nGridOutDown = -1
      Weight_I = -1; iOrderHere_I = 0
      if(present(iUTriangle1_I))then
         X1_D = XyzGrid_DI(:,iOrder_I(iUTriangle1_I(1)))
         X2_D = XyzGrid_DI(:,iOrder_I(iUTriangle1_I(2)))
         X3_D = XyzGrid_DI(:,iOrder_I(iUTriangle1_I(3)))
         AlphaUp = triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
              triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
         XyzUp_D = Xyz_D + AlphaUp * Dir_D
         call triangle(X1_D, X2_D, X3_D, XyzUp_D)
         if(all(Weight_I(1:3)>=0.0))then
            !\
            !The ray for point Xyz is projected into this triangle
            !/
            if(AlphaUp==0.0)then
               !\
               !Point Xyz belongs to this traingle
               !/
               iOrder_I(1:3) = iOrder_I(iUTriangle1_I)
               nGridOut = 3
               RETURN
            elseif(AlphaUp < 0)then
               !\
               ! Point is above the upper triangle
               ! return with negative weights and nGridOut = -1
               RETURN 
            else
               iOrderHere_I(1:3) = iOrder_I(iUTriangle1_I)
               nGridOutUp = 3
               !\
               ! Exit if for triangles with positive AlphaUp and nGridOut = 3
               !/
            end if
         elseif(present(iUTriangle2_I))then
            X1_D = XyzGrid_DI(:,iOrder_I(iUTriangle2_I(1)))
            X2_D = XyzGrid_DI(:,iOrder_I(iUTriangle2_I(2)))
            X3_D = XyzGrid_DI(:,iOrder_I(iUTriangle2_I(3)))
            AlphaUp = triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
                 triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
            XyzUp_D = Xyz_D + AlphaUp * Dir_D
            call triangle(X1_D, X2_D, X3_D, XyzUp_D)
            if(all(Weight_I(1:3)>=0.0))then
               !\
               !The ray for point Xyz is projected into this triangle
               !/
               if(AlphaUp==0.0)then
                  !\
                  !Point Xyz belongs to this traingle
                  !/
                  iOrder_I(1:3) = iOrder_I(iUTriangle2_I)
                  nGridOut = 3
                  RETURN
               elseif(AlphaUp < 0)then
                  !\
                  ! Point is above the upper triangle
                  ! return with negative weights and nGridOut = -1
                  RETURN 
               else             
                  iOrderHere_I(1:3) = iOrder_I(iUTriangle2_I)
                  nGridOutUp = 3
                  !\
                  ! Exit if for triangles with positive AlphaUp
                  !/
               end if
            elseif(present(iUTriangle3_I))then
               X1_D = XyzGrid_DI(:,iOrder_I(iUTriangle3_I(1)))
               X2_D = XyzGrid_DI(:,iOrder_I(iUTriangle3_I(2)))
               X3_D = XyzGrid_DI(:,iOrder_I(iUTriangle3_I(3)))
               AlphaUp = &
                    triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
                    triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
               XyzUp_D = Xyz_D + AlphaUp * Dir_D
               call triangle(X1_D, X2_D, X3_D, XyzUp_D)
               if(all(Weight_I(1:3)>=0.0))then
                  !\
                  !The ray for point Xyz is projected into this triangle
                  !/
                  if(AlphaUp==0.0)then
                     !\
                     !Point Xyz belongs to this traingle
                     !/
                     iOrder_I(1:3) = iOrder_I(iUTriangle3_I)
                     nGridOut = 3
                     RETURN
                  elseif(AlphaUp < 0)then
                     !\
                     ! Point is above the upper triangle
                     ! return with negative weights and nGridOut = -1
                     RETURN 
                  else        
                     iOrderHere_I(1:3) = iOrder_I(iUTriangle3_I)
                     nGridOutUp = 3
                     !\
                     ! Exit if for triangles with positive AlphaUp 
                     !/

                  end if
               elseif(present(iUTriangle4_I))then
                  X1_D = XyzGrid_DI(:,iOrder_I(iUTriangle4_I(1)))
                  X2_D = XyzGrid_DI(:,iOrder_I(iUTriangle4_I(2)))
                  X3_D = XyzGrid_DI(:,iOrder_I(iUTriangle4_I(3)))
                  AlphaUp = &
                       triple_product(X1_D - Xyz_D,X2_D - X1_D,X3_D - X1_D)/&
                       triple_product(Dir_D       ,X2_D - X1_D,X3_D - X1_D)
                  XyzUp_D = Xyz_D + AlphaUp * Dir_D
                  call triangle(X1_D, X2_D, X3_D, Xyz_D)
                  if(all(Weight_I(1:3)>=0.0))then
                     !\
                     !The ray for point Xyz is projected into this triangle
                     !/
                     if(AlphaUp==0.0)then
                        !\
                        !Point Xyz belongs to this traingle
                        !/
                        iOrder_I(1:3) = iOrder_I(iUTriangle4_I)
                        nGridOut = 3
                        RETURN
                     elseif(AlphaUp < 0.0)then
                        !\
                        ! Point is above the upper triangle
                        ! return with negative weights and nGridOut = -1
                        RETURN 
                     else   
                        iOrderHere_I(1:3) = iOrder_I(iUTriangle4_I)
                        nGridOutUp = 3
                        !\
                        ! Exit if for triangles with positive AlphaUp 
                        !/
                     end if
                  end if
               end if       !4
            end if          !3
         end if             !2
      end if                !1
      if(nGridOutUp < 1)then
         if(present(iURectangle1_I))then
            X1_D = XyzGrid_DI(:,iOrder_I(iURectangle1_I(1)))
            X2_D = XyzGrid_DI(:,iOrder_I(iURectangle1_I(2)))
            X3_D = XyzGrid_DI(:,iOrder_I(iURectangle1_I(3)))
            X4_D = XyzGrid_DI(:,iOrder_I(iURectangle1_I(4)))
            AlphaUp = triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
                 triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
            XyzUp_D = Xyz_D + AlphaUp * Dir_D
            call rectangle(X1_D, X2_D, X3_D, X4_D, XyzUp_D)
            if(all(Weight_I(1:4)>=0.0))then
               !\
               !The ray for point Xyz is projected into this rectangle
               !/
               if(AlphaUp==0.0)then
                  !\
                  !Point Xyz belongs to this rectangle
                  !/
                  iOrder_I(1:4) = iOrder_I(iURectangle1_I)
                  nGridOut = 4
                  RETURN
               elseif(AlphaUp < 0.0)then
                  !\
                  ! Point is above the upper rectangle
                  ! return with negative weights and nGridOut = -1
                  RETURN 
               else 
                  iOrderHere_I(1:4) = iOrder_I(iURectangle1_I)
                  nGridOutUp = 4
                  !\
                  ! Exit if for rectangles with positive AlphaUp 
                  !/
               end if
            elseif(present(iURectangle2_I))then
               X1_D = XyzGrid_DI(:,iOrder_I(iURectangle2_I(1)))
               X2_D = XyzGrid_DI(:,iOrder_I(iURectangle2_I(2)))
               X3_D = XyzGrid_DI(:,iOrder_I(iURectangle2_I(3)))
               X4_D = XyzGrid_DI(:,iOrder_I(iURectangle2_I(4)))
               AlphaUp = triple_product(X1_D - Xyz_D,X2_D - X1_D,X3_D - X1_D)/&
                    triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
               XyzUp_D = Xyz_D + AlphaUp * Dir_D
               call rectangle(X1_D, X2_D, X3_D, X4_D, XyzUp_D)
               if(all(Weight_I(1:4)>=0.0))then
                  !\
                  !The ray for point Xyz is projected into this rectangle
                  !/
                  if(AlphaUp==0.0)then
                     !\
                     !Point Xyz belongs to this rectangle
                     !/
                     iOrder_I(1:4) = iOrder_I(iURectangle2_I)
                     nGridOut = 4
                     RETURN
                  elseif(AlphaUp < 0.0)then
                     !\
                     ! Point is above the upper rectangle
                     ! return with negative weights and nGridOut = -1
                     RETURN 
                  else 
                     !\
                     ! Exit if for rectangles with positive AlphaUp 
                     !/
                     iOrderHere_I(1:4) = iOrder_I(iURectangle2_I)
                     nGridOutUp = 4
                  end if
               elseif(present(iURectangle3_I))then
                  X1_D = XyzGrid_DI(:,iOrder_I(iURectangle3_I(1)))
                  X2_D = XyzGrid_DI(:,iOrder_I(iURectangle3_I(2)))
                  X3_D = XyzGrid_DI(:,iOrder_I(iURectangle3_I(3)))
                  X4_D = XyzGrid_DI(:,iOrder_I(iURectangle3_I(4)))
                  AlphaUp = &
                       triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
                       triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
                  XyzUp_D = Xyz_D + AlphaUp * Dir_D
                  call rectangle(X1_D, X2_D, X3_D, X4_D, XyzUp_D)
                  if(all(Weight_I(1:4)>=0.0))then
                     !\
                     !The ray for point Xyz is projected into this rectangle
                     !/
                     if(AlphaUp==0.0)then
                        !\
                        !Point Xyz belongs to this rectangle
                        !/
                        iOrder_I(1:4) = iOrder_I(iURectangle3_I)
                        nGridOut = 4
                        RETURN
                     elseif(AlphaUp < 0.0)then
                        !\
                        ! Point is above the upper rectangle
                        ! return with negative weights and nGridOut = -1
                        RETURN 
                     else 
                        !\
                        ! Exit if for rectangles with positive AlphaUp 
                        !/
                        iOrderHere_I(1:4) = iOrder_I(iURectangle3_I)
                        nGridOutUp = 4
                     end if
                  end if
               end if     !3
            end if        !2
         end if           !1
      end if
      if(nGridOutUp < 1)then
         if(present(iUTrapezoid1_I))then
            X1_D = XyzGrid_DI(:,iOrder_I(iUTrapezoid1_I(1)))
            X2_D = XyzGrid_DI(:,iOrder_I(iUTrapezoid1_I(2)))
            X3_D = XyzGrid_DI(:,iOrder_I(iUTrapezoid1_I(3)))
            X4_D = XyzGrid_DI(:,iOrder_I(iUTrapezoid1_I(4)))
            AlphaUp = triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
                 triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
            XyzUp_D = Xyz_D + AlphaUp * Dir_D
            call Trapezoid(X1_D, X2_D, X3_D, X4_D, XyzUp_D)
            if(all(Weight_I(1:4)>=0.0))then
               !\
               !The ray for point Xyz is projected into this trapezoid
               !/
               if(AlphaUp==0.0)then
                  !\
                  !Point Xyz belongs to this rectangle
                  !/
                  iOrder_I(1:4) = iOrder_I(iUTrapezoid1_I)
                  nGridOut = 4
                  RETURN
               elseif(AlphaUp < 0.0)then
                  !\
                  ! Point is above the upper trapezoid
                  ! return with negative weights and nGridOut = -1
                  RETURN 
               else 
                  !\
                  ! Exit if for rectangles with positive AlphaUp 
                  !/
                  iOrderHere_I(1:4) = iOrder_I(iUTrapezoid1_I)
                  nGridOutUp = 4
               end if
            end if
         end if
      end if
      !\
      !If no intersection point with upper boundary is found
      !/
      if(nGridOutUp == -1)RETURN
      !\
      !Calculate low face
      if(present(iDTriangle1_I))then
         Weight_I(4:3+nGridOutUp) = Weight_I(1:nGridOutUp)
         Weight_I(1:3) = 0
         iOrderHere_I(4:3+nGridOutUp) = iOrderHere_I(1:nGridOutUp)
         iOrderHere_I(1:3) = 0
         X1_D = XyzGrid_DI(:,iOrder_I(iDTriangle1_I(1)))
         X2_D = XyzGrid_DI(:,iOrder_I(iDTriangle1_I(2)))
         X3_D = XyzGrid_DI(:,iOrder_I(iDTriangle1_I(3)))
         AlphaDown = &
              triple_product(Xyz_D - X1_D, X2_D - X1_D, X3_D - X1_D)/&
              triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
         XyzDown_D = Xyz_D - AlphaDown*Dir_D
         call triangle(X1_D, X2_D, X3_D, XyzDown_D)
         if(all(Weight_I(1:3)>=0.0))then
            !\
            !The ray for point Xyz is projected into this triangle
            !/
            if(AlphaDown==0.0)then
               !\
               !Point Xyz belongs to this triangle
               !/
               iOrder_I(1:3) = iOrder_I(iDTriangle1_I)
               nGridOut = 3
               RETURN
            elseif(AlphaDown < 0.0)then
               !\
               ! Point is below the lower triangle
               ! return with negative weights and nGridOut = -1
               !/
               RETURN 
            else
               !\
               ! Exit if for triangles with positive AlphaDown 
               !/
               iOrderHere_I(1:3) = iOrder_I(iDTriangle1_I)
               nGridOutDown = 3
            end if
         elseif(present(iDTriangle2_I))then
            X1_D = XyzGrid_DI(:,iOrder_I(iDTriangle2_I(1)))
            X2_D = XyzGrid_DI(:,iOrder_I(iDTriangle2_I(2)))
            X3_D = XyzGrid_DI(:,iOrder_I(iDTriangle2_I(3)))
            AlphaDown = &
                 triple_product(Xyz_D - X1_D, X2_D - X1_D, X3_D - X1_D)/&
                 triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
            XyzDown_D = Xyz_D - AlphaDown * Dir_D
            call triangle(X1_D, X2_D, X3_D, XyzDown_D)
            if(all(Weight_I(1:3)>=0.0))then
               !\
               !The ray for point Xyz is projected into this triangle
               !/
               if(AlphaDown==0.0)then
                  !\
                  !Point Xyz belongs to this triangle
                  !/
                  iOrder_I(1:3) = iOrder_I(iDTriangle2_I)
                  nGridOut = 3
                  RETURN
               elseif(AlphaDown < 0.0)then
                  !\
                  ! Point is below the lower triangle
                  ! return with negative weights and nGridOut = -1
                  !/
                  RETURN
               else
                  !\
                  ! Exit if for triangles with positive AlphaDown 
                  !/
                  iOrderHere_I(1:3) = iOrder_I(iDTriangle2_I)
                  nGridOutDown = 3
               end if
            elseif(present(iDTriangle3_I))then
               X1_D = XyzGrid_DI(:,iOrder_I(iDTriangle3_I(1)))
               X2_D = XyzGrid_DI(:,iOrder_I(iDTriangle3_I(2)))
               X3_D = XyzGrid_DI(:,iOrder_I(iDTriangle3_I(3)))
               AlphaDown = &
                    triple_product(Xyz_D - X1_D,X2_D - X1_D,X3_D - X1_D)/&
                    triple_product(Dir_D       ,X2_D - X1_D,X3_D - X1_D)
               XyzDown_D = Xyz_D - AlphaDown * Dir_D
               call triangle(X1_D, X2_D, X3_D, XyzDown_D)
               if(all(Weight_I(1:3)>=0.0))then
                  !\
                  !The ray for point Xyz is projected into this triangle
                  !/
                  if(AlphaDown==0.0)then
                     !\
                     !Point Xyz belongs to this triangle
                     !/
                     iOrder_I(1:3) = iOrder_I(iDTriangle3_I)
                     nGridOut = 3
                     RETURN
                  elseif(AlphaDown < 0.0)then
                     !\
                     ! Point is below the lower triangle
                     ! return with negative weights and nGridOut = -1
                     !/
                     RETURN
                  else
                     !\
                     ! Exit if for triangles with positive AlphaDown 
                     !/
                     iOrderHere_I(1:3) = iOrder_I(iDTriangle3_I)
                     nGridOutDown = 3
                  end if
               elseif(present(iDTriangle4_I))then
                  X1_D = XyzGrid_DI(:,iOrder_I(iDTriangle4_I(1)))
                  X2_D = XyzGrid_DI(:,iOrder_I(iDTriangle4_I(2)))
                  X3_D = XyzGrid_DI(:,iOrder_I(iDTriangle4_I(3)))
                  AlphaDown = &
                       triple_product(Xyz_D - X1_D,X2_D - X1_D,X3_D - X1_D)/&
                       triple_product(Dir_D       ,X2_D - X1_D,X3_D - X1_D)
                  XyzDown_D = Xyz_D - AlphaDown * Dir_D
                  call triangle(X1_D, X2_D, X3_D, XyzDown_D)
                  if(all(Weight_I(1:3)>=0.0))then
                     !\
                     !The ray for point Xyz is projected into this triangle
                     !/
                     if(AlphaDown==0.0)then
                        !\
                        !Point Xyz belongs to this triangle
                        !/
                        iOrder_I(1:3) = iOrder_I(iDTriangle4_I)
                        nGridOut = 3
                        RETURN
                     elseif(AlphaDown < 0.0)then
                        !\
                        ! Point is below the lower triangle
                        ! return with negative weights and nGridOut = -1
                        !/
                        RETURN
                     else
                        !\
                        ! Exit if for triangles with positive AlphaDown 
                        !/
                        iOrderHere_I(1:3) = iOrder_I(iDTriangle4_I)
                        nGridOutDown = 3
                     end if
                  elseif(present(iDTriangle5_I))then
                     X1_D = XyzGrid_DI(:,iOrder_I(iDTriangle5_I(1)))
                     X2_D = XyzGrid_DI(:,iOrder_I(iDTriangle5_I(2)))
                     X3_D = XyzGrid_DI(:,iOrder_I(iDTriangle5_I(3)))
                     AlphaDown = &
                          triple_product(&
                          Xyz_D - X1_D,X2_D - X1_D,X3_D - X1_D)/&
                          triple_product(&
                          Dir_D       ,X2_D - X1_D,X3_D - X1_D)
                     XyzDown_D = Xyz_D - AlphaDown * Dir_D
                     call triangle(X1_D, X2_D, X3_D, XyzDown_D)
                     if(all(Weight_I(1:3)>=0.0))then
                        !\
                        !The ray for point Xyz is projected into this triangle
                        !/
                        if(AlphaDown==0.0)then
                           !\
                           !Point Xyz belongs to this triangle
                           !/
                           iOrder_I(1:3) = iOrder_I(iDTriangle5_I)
                           nGridOut = 3
                           RETURN
                        elseif(AlphaDown < 0.0)then
                           !\
                           ! Point is below the lower triangle
                           ! return with negative weights and nGridOut = -1
                           !/
                           RETURN
                        else
                           !\
                           ! Exit if for triangles with positive AlphaDown 
                           !/
                           iOrderHere_I(1:3) = iOrder_I(iDTriangle5_I)
                           nGridOutDown = 3
                        end if
                     elseif(present(iDTriangle6_I))then
                        X1_D = XyzGrid_DI(:,iOrder_I(iDTriangle6_I(1)))
                        X2_D = XyzGrid_DI(:,iOrder_I(iDTriangle6_I(2)))
                        X3_D = XyzGrid_DI(:,iOrder_I(iDTriangle6_I(3)))
                        AlphaDown = &
                             triple_product(&
                             Xyz_D - X1_D,X2_D - X1_D,X3_D - X1_D)/&
                             triple_product(&
                             Dir_D       ,X2_D - X1_D,X3_D - X1_D)
                        XyzDown_D = Xyz_D - AlphaDown * Dir_D
                        call triangle(X1_D, X2_D, X3_D, XyzDown_D)
                        if(all(Weight_I(1:3)>=0.0))then
                           !\
                           !The ray for point Xyz is projected into this 
                           !triangle
                           !/
                           if(AlphaDown==0.0)then
                              !\
                              !Point Xyz belongs to this triangle
                              !/
                              iOrder_I(1:3) = iOrder_I(iDTriangle6_I)
                              nGridOut = 3
                              RETURN
                           elseif(AlphaDown < 0.0)then
                              !\
                              ! Point is below the lower triangle
                              ! return with negative weights and nGridOut = -1
                              !/
                              RETURN
                           else
                              !\
                              ! Exit if for triangles with positive AlphaDown 
                              !/
                              iOrderHere_I(1:3) = iOrder_I(iDTriangle6_I)
                              nGridOutDown = 3
                           end if
                        end if
                     end if !6
                  end if    !5
               end if       !4
            end if          !3
         end if             !2
      end if                !1
      if(nGridOutDown <1)then
         if(present(iDRectangle1_I))then
            if(present(iDTriangle1_I))then
               Weight_I(5:4+nGridOutUp) = Weight_I(4:3+nGridOutUp)
               Weight_I(1:4) = 0
               iOrderHere_I(5:4+nGridOutUp) = iOrderHere_I(4:3+nGridOutUp)
               iOrderHere_I(1:4)=0
            else
               Weight_I(5:4+nGridOutUp) = Weight_I(1:nGridOutUp)
               Weight_I(1:4) = 0
               iOrderHere_I(5:4+nGridOutUp) = iOrderHere_I(1:nGridOutUp)
               iOrderHere_I(1:4)=0
            end if
            X1_D = XyzGrid_DI(:,iOrder_I(iDRectangle1_I(1)))
            X2_D = XyzGrid_DI(:,iOrder_I(iDRectangle1_I(2)))
            X3_D = XyzGrid_DI(:,iOrder_I(iDRectangle1_I(3)))
            X4_D = XyzGrid_DI(:,iOrder_I(iDRectangle1_I(4)))
            AlphaDown = triple_product(Xyz_D - X1_D,X2_D - X1_D,X3_D - X1_D)/&
                 triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
            XyzDown_D = Xyz_D - AlphaDown * Dir_D
            call rectangle(X1_D, X2_D, X3_D, X4_D, XyzDown_D)
            if(all(Weight_I(1:4)>=0.0))then
               !\
               !The ray for point Xyz is projected into this rectangle
               !/
               if(AlphaDown==0.0)then
                  !\
                  !Point Xyz belongs to this rectangle
                  !/
                  iOrder_I(1:4) = iOrder_I(iDRectangle1_I)
                  nGridOut = 4
                  RETURN
               elseif(AlphaDown < 0.0)then
                  !\
                  ! Point is below the lower rectangle
                  ! return with negative weights and nGridOut = -1
                  !/
                  RETURN
               else
                  !\
                  ! Exit if for rectangles with positive AlphaDown 
                  !/
                  iOrderHere_I(1:4) = iOrder_I(iDRectangle1_I)
                  nGridOutDown = 4
               end if
            end if
         end if
      end if           !1
      if(nGridOutDown < 1)then
         if(present(iDTrapezoid1_I))then
            if(.not.present(iDRectangle1_I))then
               if(present(iDTriangle1_I))then
                  Weight_I(5:4+nGridOutUp) = Weight_I(4:3+nGridOutUp)
                  Weight_I(1:4)=0
                  iOrderHere_I(5:4+nGridOutUp) = iOrderHere_I(4:3+nGridOutUp)
                  iOrderHere_I(1:4) = 0
               else
                  Weight_I(5:4+nGridOutUp) = Weight_I(1:nGridOutUp)
                  Weight_I(1:4)=0
                  iOrderHere_I(5:4+nGridOutUp) = iOrderHere_I(1:nGridOutUp)
                  iOrderHere_I(1:4) = 0
               end if
            end if
            X1_D = XyzGrid_DI(:,iOrder_I(iDTrapezoid1_I(1)))
            X2_D = XyzGrid_DI(:,iOrder_I(iDTrapezoid1_I(2)))
            X3_D = XyzGrid_DI(:,iOrder_I(iDTrapezoid1_I(3)))
            X4_D = XyzGrid_DI(:,iOrder_I(iDTrapezoid1_I(4)))
            AlphaDown = triple_product(Xyz_D - X1_D,X2_D - X1_D,X3_D - X1_D)/&
                 triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
            XyzDown_D = Xyz_D - AlphaDown * Dir_D
            call trapezoid(X1_D, X2_D, X3_D, X4_D, XyzDown_D)
            if(all(Weight_I(1:4)>=0.0))then
               !\
               !The ray for point Xyz is projected into this trapezoid
               !/
               if(AlphaDown==0.0)then
                  !\
                  !Point Xyz belongs to this rectangle
                  !/
                  iOrder_I(1:4) = iOrder_I(iDTrapezoid1_I)
                  nGridOut = 4
                  RETURN
               elseif(AlphaDown < 0.0)then
                  !\
                  ! Point is below the lower rectangle
                  ! return with negative weights and nGridOut = -1
                  !/
                  RETURN
               else
                  !\
                  ! Exit if for rectangles with positive AlphaDown 
                  !/
                  iOrderHere_I(1:4) = iOrder_I(iDTrapezoid1_I)
                  nGridOutDown = 4
               end if
            end if
         end if           !1
      end if
      !\
      !If no intersection point with the down subface is found
      !/
      if(nGridOutDown == -1)RETURN 
      !\
      !Apply the weight for interpolation along the rate
      !/
      nGridOut = nGridOutUp + nGridOutDown
      iOrder_I(1:nGridOut) = iOrderHere_I(1:nGridOut)
      Weight_I(1:nGridOutDown) = &
           Weight_I(1:nGridOutDown)*AlphaUp/(AlphaUp + AlphaDown)
      Weight_I(nGridOutDown+1:nGridOut) = &
           Weight_I(nGridOutDown+1:nGridOut)*AlphaDown/(AlphaUp + AlphaDown)
    end subroutine parallel_rays
    !========================
    function cross_product(a_D, B_d)
      real, dimension(nDim), intent(in) :: a_D, b_D
      real, dimension(nDim) :: cross_product
      !-----------------------------------------------
      cross_product(x_) = a_D(y_)*b_D(z_) - a_D(z_)*b_D(y_)
      cross_product(y_) = a_D(z_)*b_D(x_) - a_D(x_)*b_D(z_)
      cross_product(z_) = a_D(x_)*b_D(y_) - a_D(y_)*b_D(x_)
    end function cross_product
    !=========================
    real function triple_product(a_D, b_D, c_D)
      real, dimension(nDim), intent(in) :: a_D, b_D, c_D
      !-----------------------------------------------
      triple_product = sum(a_D * cross_product(b_D, c_D))
    end function triple_product
    !===================================================
    subroutine tetrahedron(X1_D, X2_D, X3_D, X4_D)
      !\
      ! Interpolate in the tetrahedron
      !/
      real, dimension(nDim), intent(in) :: X1_D, X2_D, X3_D, X4_D
      real:: Aux
      !-----------------------------------------------
      !\
      ! Need to solve an equation:
      ! Xyz = Weight_I(1)*X1+Weight_I(2)*X2+Weight_I(3)*X3+Weight_I(4)*X4
      ! Or, which is equivalent:
      ! Xyz -X1 = Weight_I(2)*(X2-X1)+Weight_I(3)*(X3-X1)+Weight_I(4)*(X4-X1)
      !/
      Aux= 1/triple_product( X4_D - X1_D,  X3_D - X1_D,  X2_D - X1_D)
      Weight_I(4) = &
           triple_product(Xyz_D - X1_D,  X3_D - X1_D,  X2_D - X1_D)*Aux
      Weight_I(3) = &
           triple_product( X4_D - X1_D, Xyz_D - X1_D,  X2_D - X1_D)*Aux
      Weight_I(2) = &
           triple_product( X4_D - X1_D,  X3_D - X1_D, Xyz_D - X1_D)*Aux
      Weight_I(1) = 1 - sum(Weight_I(2:4))
    end subroutine tetrahedron
    !=======================================================================
    subroutine pyramid(iBase,X1_D, X2_D, X3_D, X4_D, X5_D)
      !\
      ! Interpolate in the pyramid with base X1X2X3X4 and apex X5
      ! valid for case of rectangular or trapezoid base
      !/
      integer, intent(in) :: iBase
      real, dimension(nDim), intent(in) :: X1_D, X2_D, X3_D, X4_D, X5_D
      !\
      !Projection of Xyz point on the base of the pyramid X1X2X3X4 
      !along the line X5Xyz
      !/
      real, dimension(nDim) :: XyzP_D 
      real ::Alpha5
      !-----------------------------------------------
      !\
      ! Solve equation: X5 + (Xyz - X5)/Alpha5 = XyzP
      ! where XyzP belongs to the pydamid base.
      ! As long as (XyzP-X1_D)\cdot[(X2-X1)\times(X3-X1)]=0,
      ! we have Alpha5 =(&
      ! X5 -Xyz)\cdot[(X2-X1)\times(X3-X1)]/(X5-X1)\cdot[(X2-X1)\times(X3-X1)]
      !/

      Alpha5 = triple_product(X5_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
           triple_product(X5_D - X1_D , X2_D - X1_D, X3_D - X1_D)

      if(Alpha5==0.0)then
         if(all(Xyz_D(:)==X5_D(:)))then
            Weight_I(5) = 1
            Weight_I(1:4) = 0
            RETURN
         else
            Weight_I(1:5) = -1
            RETURN
         end if
      elseif(Alpha5 < 0.0 .or. Alpha5 > 1.0)then
         Weight_I(1:5) = -1
         RETURN
      end if

      XyzP_D = X5_D + (Xyz_D-X5_D)/Alpha5
      !\
      ! Now, Xyz = Alpha5*XyzP + (1-Alpha5)*X5
      ! Find weight of the apex point X5_D
      !/
      Weight_I(5) = 1 - Alpha5

      !\
      ! Find weights of the base points X1, X2, X3, X4
      !/
      select case(iBase)
      case(Rectangular_)
         call rectangle(X1_D, X2_D, X3_D, X4_D, XyzP_D)
      case(Trapezoidal_)
         call trapezoid(X1_D, X2_D, X3_D, X4_D, XyzP_D)
      end select

      if(any(Weight_I(1:4) < 0.0))RETURN
      !\
      ! Correct weights due to pyramid geometry
      !/
      Weight_I(1) = Weight_I(1) * (1 - Weight_I(5))
      Weight_I(2) = Weight_I(2) * (1 - Weight_I(5))
      Weight_I(3) = Weight_I(3) * (1 - Weight_I(5))
      Weight_I(4) = Weight_I(4) * (1 - Weight_I(5))

    end subroutine pyramid
    !============
    subroutine rectangle(X1_D, X2_D, X3_D, X4_D, XyzP_D)
      real, dimension(nDim), intent(in):: X1_D, X2_D, X3_D, X4_D, XyzP_D
      real:: x,y
      !---------
      !Calculate dimensionless coordinates with respect to vertex 1 
      y =  sum((X3_D - X1_D)*(XyzP_D - X1_D))/&
           sum((X3_D - X1_D)*(X3_D - X1_D))
      x =  sum((XyzP_D - X1_D)*(X2_D - X1_D))/&
           sum((X2_D - X1_D)*(X2_D - X1_D) )

      Weight_I(3)  =     y *(1 - x) ; Weight_I(4) =      y *x
      Weight_I(1) = (1 - y)*(1 - x) ; Weight_I(2) = (1 - y)*x

    end subroutine rectangle
    !==============================
    subroutine trapezoid(X1_D, X2_D, X3_D, X4_D, XyzP_D)
      real, dimension(nDim), intent(in):: X1_D, X2_D, X3_D, X4_D, XyzP_D
      real:: x,y, XyzMin_D(nDim), XyzMax_D(nDim)
      !We require that the lager base is X1X2
      !---------
      !Calculate dimensionless coordinates with respect to vertex 1
      XyzMin_D = 0.50*(X1_D + X2_D);  XyzMax_D = 0.50*(X3_D + X4_D)

      y =  sum((XyzMax_D - XyzMin_D)*(XyzP_D - XyzMin_D))/&
           sum((XyzMax_D - XyzMin_D)*(XyzMax_D - XyzMin_D))
      x =  sum((XyzP_D   - X1_D)*(X2_D - X1_D))/&
           sum( (X2_D    - X1_D)*(X2_D - X1_D))
      if( y < 0.0 .or. y > 1.0 .or. x < 0.0 .or. x > 1.0)then
         Weight_I(1:4) = -1
         RETURN
      end if
      if( x <= 0.250) then
         !Interpolation in triangle 132, 4th weight is zero
         Weight_I(3) = y                   ; Weight_I(4) = 0
         Weight_I(1) = 1 - 0.750*y - x     ; Weight_I(2) = x - 0.250*y
      elseif(x <= 0.750) then
         !Bilinear interpolation
         Weight_I(3) = y*(1 - 2*(x - 0.25)); Weight_I(4) = y*2*(x - 0.25)
         Weight_I(1) = (1 - y)*(1 - x)     ; Weight_I(2) = (1 - y)*x
      else
         !Interpolation in triangle 142, 3rd weight is zero
         Weight_I(3) = 0                   ; Weight_I(4) = y
         Weight_I(1) = 1 - 0.250*y - x     ; Weight_I(2) = x - 0.750*y
      end if

    end subroutine trapezoid
    !=====================================
    subroutine triangle(X1_D, X2_D, X3_D,  XyzP_D)
      real, dimension(nDim), intent(in):: X1_D, X2_D, X3_D, XyzP_D
      real:: CrossProductInv_D(nDim)
      !\
      ! In 2D case
      ! Weight_I(3) = cross_product(XyzP_D-X1_D,X2_D-X1_D)/&
      !               cross_product(X3_D-X1_D,X2_D-X1_D)
      ! Weight_I(2) = cross_product(X3_D-X1_D,Xyz_D-X1_D)/&
      !               cross_product(X3_D - X1_D, X2_D-X1_D)
      ! iN 3d case instead of /cross_product(X3_D - X1_D, X2_D-X1_D) we use 
      ! \cdot cross_product(X3_D - X1_D, X2_D-X1_D)/&
      !       sum(cross_product(X3_D - X1_D, X2_D-X1_D)**2)
      !/
      CrossProductInv_D = cross_product(X3_D - X1_D, X2_D - X1_D)
      CrossProductInv_D = CrossProductInv_D/sum(CrossProductInv_D**2)

      Weight_I(3) = sum(cross_product(XyzP_D - X1_D, X2_D  - X1_D)*&
           CrossProductInv_D)
      Weight_I(2) = sum(cross_product(X3_D  - X1_D, XyzP_D - X1_D)*&
           CrossProductInv_D)
      Weight_I(1) = 1 - sum(Weight_I(2:3))

    end subroutine triangle
    !====================================================================
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
                 XyzGrid_DI(kAxis,iOrder_I(nFine +5 )))
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
            DoStencilFix = .false.
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

      DoStencilFix = dXyzUp * dXyzDown > 0.0
      if(DoStencilFix)RETURN

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
       XyzStencil_D, iCase)
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
    if(iGrid==0)RETURN 
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
