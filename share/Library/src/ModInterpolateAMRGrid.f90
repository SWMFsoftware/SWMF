!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModInterpolateAMRGrid
  implicit none

  !Generalize bilinear and trilinear interpolation for AMR grids
  !The data are given at the cell-centered grid which consists of AMR blocks.
  !The interpolation is free of any jumps at the resolution interfaces
  !including edges and corners

  PRIVATE !Except
  SAVE

  interface interpolate_amr
     module procedure interpolate_block_amr
  end interface
  public interpolate_amr

  ! Unit test
  public:: test_interpolate_amr

  integer, parameter, public :: BehindTheBoundary_ = -7777

  integer, parameter:: &
       Coarse_  = 0,               &
       Fine_    = 1,               &
       x_       = 1,               &
       y_       = 2,               &
       z_       = 3
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
  ! the grid point numbers. The functions used in more than one place,
  ! therefore we delegate them here. 
  !Parameters to enumerate faces or face diagonal
  integer, parameter:: Xy_ = 3, Xz_ = 2, Yz_ = 1

  !Number of the vertex connected by 
  !the edge of direction iDir (second index) 
  !with the given vertex (first index)
  integer, dimension(8,3), parameter:: iEdge_ID = &
       reshape((/&   !Number of the connected vertex
       2 ,1, 4, 3, 6, 5, 8, 7,        & !Edge x
       3, 4, 1, 2, 7, 8, 5, 6,        & !Edge y
       5, 6, 7, 8, 1, 2, 3, 4         & !Edge z
       /),(/8,3/))

  !Number of the vertex connected by 
  !the face diagonal across the face of direction iDir (second index) 
  !with the given vertex (first index)
  integer, dimension(8,3), parameter:: iFaceDiag_ID = &
       reshape((/&   !Number of the connected vertex
       7 ,8, 5, 6, 3, 4, 1, 2,        & !Face yz
       6, 5, 8, 7, 2, 1, 4, 3,        & !Face xz
       4, 3, 2, 1, 8, 7, 6, 5         & !Face xy
       /),(/8,3/))

  !Number of the vertex connected by 
  !the main diagonal with the given vertex (index)
  integer, dimension(8), parameter:: iMainDiag_I = &
       (/ 8, 7, 6, 5, 4, 3, 2, 1 /)

  !Vertexes (enumerated by the first undex), which form 
  !the face of direction iDir (second index) including the 
  !given vertex (the third index). 
  !When the first index equals 1,2,3,4, the vertex, accordingly: 
  !(1) coincides with the given one
  !(2) is connected with the given one by the edge of direction iDir+1
  !(3) is connected with the given one by the edge of direction iDir+2
  !(4) is connected to the given one by the face diagonal of direction iDir
  integer, dimension(4,3,8), parameter:: iFace_IDI = &
       reshape((/&   ! yz face   ! xz face      ! xy face ! 
       1, 3, 5, 7,   1, 5, 2, 6,   1, 2, 3, 4, & 
       2, 4, 6, 8,   2, 6, 1, 5,   2, 1, 4, 3, &
       3, 1, 7, 5,   3, 7, 4, 8,   3, 4, 1, 2, &
       4, 2, 8, 6,   4, 8, 3, 7,   4, 3, 2, 1, &
       5, 7, 1, 3,   5, 1, 6, 2,   5, 6, 7, 8, &
       6, 8, 2, 4,   6, 2, 5, 1,   6, 5, 8, 7, &
       7, 5, 3, 1,   7, 3, 8, 4,   7, 8, 5, 6, &
       8, 6, 4, 2,   8, 4, 7, 3,   8, 7, 6, 5  &
       /), (/4,3,8/))

  !Vertexes (enumerated by the first undex), which form 
  !the face of direction iDir (second index) and does not include the
  !given vertex (the third index). 
  !When the first index equals 1,2,3,4, the vertex, accordingly: 
  !(1) is connected to the given one by the edge of direction iDir
  !(2) is connected to (1) by the edge of direction iDir+1
  !(3) is connected to (1) by the edge of direction iDir+2
  !(4) is connected to  by the face diagonal of direction iDir
  integer, dimension(4,3,8), parameter:: iOppositeFace_IDI = &
       reshape((/&   ! yz face   ! xz face      ! xy face ! 
       2, 4, 6, 8,   3, 7, 4, 8,   5, 6, 7, 8, & 
       1, 3, 5, 7,   4, 8, 3, 7,   6, 5, 8, 7, &
       4, 2, 8, 6,   1, 5, 2, 6,   7, 8, 5, 6, &
       3, 1, 7, 5,   2, 6, 1, 5,   8, 7, 6, 5, &
       6, 8, 2, 4,   7, 3, 8, 4,   1, 2, 3, 4, &
       5, 7, 1, 3,   8, 4, 7, 3,   2, 1, 4, 3, &
       8, 6, 4, 2,   5, 1, 6, 2,   3, 4, 1, 2, &
       7, 5, 3, 1,   6, 2, 5, 1,   4, 3, 2, 1  &
       /), (/4,3,8/))
  !\
  ! Arrays of logical values Xyz_D < XyzGrid_DI(:,iGrid) for iGrid'th grid
  ! point of the basic stencil for point Xyz
  !/
  logical, parameter:: IsBelow_DI(3,8)=reshape((/&
       .false., .false., .false., &
       .true. , .false., .false., &
       .false., .true. , .false., &
       .true. , .true. , .false., &
       .false., .false., .true. , &
       .true. , .false., .true. , &
       .false., .true. , .true. , &
       .true. , .true. , .true. /)&
       , (/3,8/) )
  !\
  ! Shift of the iGrid point in the stencil with respect to the
  ! first one
  !/
  integer, dimension(3,8),parameter :: iShift_DI = reshape((/&
       0, 0, 0, &
       1, 0, 0, &
       0, 1, 0, &
       1, 1, 0, &
       0, 0, 1, &
       1, 0, 1, &
       0, 1, 1, &
       1, 1, 1/),(/3,8/))
 
  !\
  ! For test: arrays of the refinement levels to be passed to
  ! find_test routine
  !/
  integer:: iLevelTest_I(8)
  !\
  !For one of the versions of the test: coordinate array
  !/
  real, allocatable::Xyz_DCB(:,:,:,:,:)

  integer:: iSeed = 1
contains
  !\
  !The random number generator.
  !From the Buneman's code TRISTAN
  !/
  subroutine init_rand(iSeedIn)
    integer,optional,intent(in)::iSeedIn
    if(present(iSeedIn))then
       iSeed=iSeedIn
    else
       iSeed=1
    end if
  end subroutine init_rand
  !=====================
  real function rand()
    iSeed = iSeed*48828125
    IF(iSeed < 0) iSeed=(iSeed+2147483647)+1
    if(iSeed==0) iSeed=1
    rand=FLOAT(iSeed)/2147483647
  end function rand
  !========================================================================
  subroutine  interpolate_amr_grid1(&
       Xyz_D, XyzGrid_DI, iLevel_I,&
       nGridOut, Weight_I, iOrder_I)

    integer,parameter :: nGrid = 2, nDim = 1

    character(LEN=*),parameter:: NameSub='interpolate_amr_grid1'

    !\
    !Input parameters
    !/
    !\
    !The location at which to interpolate the data
    !/
    real   , intent(in) :: Xyz_D(nDim) 

    !Grid point coordinates
    real  ,   intent(in) :: XyzGrid_DI(nDim,nGrid) !1 coordinate, 2 points

    !The refinement level at each grid point. By one higher level of 
    !refinement assumes the cell size reduced by a factor of 0.5
    integer,  intent(in) :: iLevel_I(nGrid)

    !\
    !Output parameters
    !/
    !The number of grid points to be ultimately included into the 
    !interpolation stencil. If nGridOut < nGridIn, only the first nGridOut 
    !elements  are meaningful in the output arrays   
    integer, intent(out) :: nGridOut  

    !The weight coefficients array.
    real   , intent(out) :: Weight_I(nGrid)

    !Order(numbers) of the grid points used for interpolation
    integer, intent(out) :: iOrder_I(nGrid)

    integer :: iGrid
    !-------------
    !\
    ! Check if the grid points are out of the grid boundary
    !/
    if(iLevel_I(2)== BehindTheBoundary_.or. XyzGrid_DI(x_,1)==Xyz_D(x_))then
       if(iLevel_I(1)==BehindTheBoundary_)&
            call CON_stop('Both points are out of the grid '//NameSub)
       iOrder_I = (/1,2/)
       nGridOut = 1
       Weight_I(1) = 1
    elseif(iLevel_I(1)== BehindTheBoundary_)then
       iOrder_I = (/2,1/)
       nGridOut = 1
       Weight_I(1) = 1
    else
       iOrder_I = (/1,2/)
       nGridOut = 2
       Weight_I(1) = &
            (XyzGrid_DI(1,2) - Xyz_D(1))/(XyzGrid_DI(1,2) - XyzGrid_DI(1,1))
       Weight_I(2) = 1 - Weight_I(1) 
    end if
  end subroutine interpolate_amr_grid1
  !=========================================================================
  subroutine  interpolate_amr_grid2(&
       Xyz_D, XyzGrid_DI, iLevel_I, &
       nGridOut, Weight_I, iOrder_I,&
       DoStencilFix, XyzStencil_D)

    integer,parameter :: nGrid = 4, nDim = 2

    character(LEN=*),parameter:: NameSub='interpolate_amr_grid2'

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
    integer,  intent(inout) :: iLevel_I(nGrid)

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

    !\
    ! The logical is true, if (for some resolution combinations) the
    ! basic stencil for the point at which to inpertpolate is not 
    ! applicable
    !/
    logical, intent(inout):: DoStencilFix
    !\
    ! If DoStencilFix==.true., the subroutine provides XyzStencil_D to be
    ! used to construct stencil, not to interpolate
    !/
    real, intent(out) :: XyzStencil_D(nDim)


    integer :: iGrid, iDir   !Loop variables
    integer :: iLoc          !Find minloc and maxloc
    !\
    !Number of grid points returned from 1D routine
    !/
    integer :: nGridOut1     
    integer :: iLevelMin     !Minimumum refinement level


    !\
    ! Min and max values of coordinates
    !/
    real    :: XyzMin_D(nDim), XyzMax_D(nDim)  
    !\
    ! To find lines at "constant" cordinate
    !/
    real    :: dXyzSmall_D(nDim)       
    !\       
    ! for calculating on a uniform grid 
    !/
    real    :: dXyzInv_D(  nDim)               
    real    :: Aux_D(nDim)                   ! Misc
    real    :: Xyz1_D(1), XyzGrid1_DI(1,2)   ! To call 1D interpolation
    integer :: iOrder1_I(2)                  ! To call 1D interpolation
    integer :: iLevel1_I(2)
    !-------------

    !\
    ! Make iLevel=0 for coarser grids and iLevel=1 for finer grids
    !/
    iLevelMin = minval(iLevel_I, MASK=iLevel_I/=BehindTheBoundary_)
    where(iLevel_I/=BehindTheBoundary_)iLevel_I = iLevel_I - iLevelMin

    if(DoStencilFix) goto 100
    !\
    ! Check if the grid points are out of the grid boundary
    !/
    if(all(iLevel_I(3:4)== BehindTheBoundary_).or.&
         all(XyzGrid_DI(y_,1:2)==Xyz_D(y_)) )then
       !\
       ! Edge 3:4 should be removed
       !/
       iOrder_I = (/1,2,3,4/)
       Xyz1_D(x_) = Xyz_D(x_) ; XyzGrid1_DI(x_,:) = XyzGrid_DI(x_,1:2)
       iLevel1_I = iLevel_I(1:2)

       call interpolate_amr_grid1(&
            Xyz1_D,  XyzGrid1_DI, iLevel1_I, &
            nGridOut, Weight_I(1:2), iOrder1_I)
       iOrder_I(1:2) = iOrder_I(iOrder1_I)
       RETURN

    elseif(all(iLevel_I(1:2)== BehindTheBoundary_))then
       !\
       ! Edge 1:2 should be removed
       !/
       iOrder_I = (/3,4,1,2/)
       Xyz1_D(x_) = Xyz_D(x_) ; XyzGrid1_DI(x_,:) = XyzGrid_DI(x_,3:4)
       iLevel1_I = iLevel_I(3:4)

       call interpolate_amr_grid1(&
            Xyz1_D,  XyzGrid1_DI, iLevel1_I, &
            nGridOut, Weight_I(1:2), iOrder1_I)
       iOrder_I(1:2) = iOrder_I(iOrder1_I)
       RETURN

    elseif(all(iLevel_I(2:4:2)== BehindTheBoundary_).or. &
         all(XyzGrid_DI(x_,1:3:2)==Xyz_D(x_)) )then
       !\
       ! Edge 2,4 should be removed
       !/
       iOrder_I = (/1,3,2,4/)
       Xyz1_D(x_) = Xyz_D(y_) ; XyzGrid1_DI(x_,:) = XyzGrid_DI(y_,1:3:2)
       iLevel1_I = iLevel_I(1:3:2)

       call interpolate_amr_grid1(&
            Xyz1_D,  XyzGrid1_DI, iLevel1_I, &
            nGridOut, Weight_I(1:2), iOrder1_I)
       iOrder_I(1:2) = iOrder_I(iOrder1_I)
       RETURN

    elseif(all(iLevel_I(1:3:2)== BehindTheBoundary_) )then
       !\
       ! Edge 1,3 should be removed
       !/
       iOrder_I = (/2,4,1,3/)
       Xyz1_D(x_) = Xyz_D(y_) ; XyzGrid1_DI(x_,:) = XyzGrid_DI(y_,2:4:2)
       iLevel1_I = iLevel_I(2:4:2)

       call interpolate_amr_grid1(&
            Xyz1_D,  XyzGrid1_DI, iLevel1_I, &
            nGridOut, Weight_I(1:2), iOrder1_I)
       iOrder_I(1:2) = iOrder_I(iOrder1_I)
       RETURN
    end if

    !\
    ! Calculate the stencil size
    !/
    do iDir = 1, nDim
       dXyzSmall_D(iDir) = 0.10*(&
            maxval(XyzGrid_DI(iDir,:)) - minval(XyzGrid_DI(iDir,:)))
       dXyzInv_D(iDir)   = 0.10/dXyzSmall_D(iDir)
    end do

    !\
    ! Frist conbinations of the resolution interfaces: 
    ! no resolution interface
    !/

    if(all(iLevel_I==0))then
       !\                       C  C
       ! No refinement
       !/                       C  C
       iOrder_I = (/1,2,3,4/)
       Aux_D = (Xyz_D - XyzGrid_DI(:,1))*dXyzInv_D
       nGridOut = 4
       Weight_I(1) = (1 - Aux_D(1))*(1 - Aux_D(2))
       Weight_I(2) =      Aux_D(1) *(1 - Aux_D(2))
       Weight_I(3) = (1 - Aux_D(1))*     Aux_D(2)
       Weight_I(4) =      Aux_D(1) *     Aux_D(2)
       RETURN

    elseif(abs(XyzGrid_DI(y_,1) - XyzGrid_DI(y_,2)) &
         < dXyzSmall_D(y_).and.&
         abs(XyzGrid_DI(y_,3) - XyzGrid_DI(y_,4)) &
         < dXyzSmall_D(y_))then
       !\                             C  C
       ! Edges going along x-axis    F F   or   F F
       !/                                     C  C
       iOrder_I = (/1,2,3,4/)
       Aux_D(y_) = &
            dXyzInv_D(y_)*(  Xyz_D(y_) - XyzGrid_DI( y_,iOrder_I(1) )  )
       Xyz1_D(x_) = Xyz_D(x_)

       !\
       ! Interpolate along edge X1X2
       !/
       XyzGrid1_DI(x_,1) = XyzGrid_DI(x_,iOrder_I(1))
       XyzGrid1_DI(x_,2) = XyzGrid_DI(x_,iOrder_I(2))
       iLevel1_I(     1) = iLevel_I(     iOrder_I(1))
       iLevel1_I(     2) = iLevel_I(     iOrder_I(2))
       call interpolate_amr_grid1(&
            Xyz1_D,  XyzGrid1_DI, iLevel1_I, &
            nGridOut1, Weight_I(1:2), iOrder1_I)
       iOrder_I(1:2) = iOrder_I(iOrder1_I)

       !\
       ! Apply weight for interpolation along another axis
       !/
       Weight_I(1:nGridOut1) =  Weight_I(1:nGridOut1) * (1 - Aux_D(y_))
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
       XyzGrid1_DI(x_,1) = XyzGrid_DI(x_,iOrder_I(nGridOut+1))
       XyzGrid1_DI(x_,2) = XyzGrid_DI(x_,iOrder_I(nGridOut+2))
       iLevel1_I(     1) = iLevel_I(     iOrder_I(nGridOut+1))
       iLevel1_I(     2) = iLevel_I(     iOrder_I(nGridOut+2))
       call interpolate_amr_grid1(&
            Xyz1_D,  XyzGrid1_DI, iLevel1_I, &
            nGridOut1, Weight_I(nGridOut+1:nGridOut+2), iOrder1_I)

       iOrder_I(nGridOut+1:nGridOut+2) = iOrder_I(nGridOut+iOrder1_I)
       !\
       ! Apply weight for interpolation along another axis
       !/
       Weight_I(nGridOut + 1:nGridOut + nGridOut1) =  &
            Weight_I(nGridOut + 1:nGridOut + nGridOut1)*Aux_D(y_)

       !\
       ! May need to remove a grid point once behind the boundary
       !/
       nGridOut = nGridOut + nGridOut1
       DoStencilFix = .false.
       RETURN

    elseif(abs(XyzGrid_DI(x_,1) - XyzGrid_DI(x_,3)) < dXyzSmall_D(x_).and.&
         abs(XyzGrid_DI(x_,2) - XyzGrid_DI(x_,4)) < dXyzSmall_D(x_))then
       !\                             C  F         F C
       ! Edges going along y-axis        F    or   F
       !/                             C              C
       iOrder_I = (/1,3,2,4/)
       Aux_D(x_) = dXyzInv_D(x_)*(Xyz_D(x_)-XyzGrid_DI(x_,iOrder_I(1)))
       Xyz1_D(x_) = Xyz_D(y_)

       !\
       ! Interpolate along edge X1X2
       !/
       XyzGrid1_DI(x_,1) = XyzGrid_DI(y_,iOrder_I(1))
       XyzGrid1_DI(x_,2) = XyzGrid_DI(y_,iOrder_I(2))
       iLevel1_I(     1) = iLevel_I(     iOrder_I(1))
       iLevel1_I(     2) = iLevel_I(     iOrder_I(2))
       call interpolate_amr_grid1(&
            Xyz1_D,  XyzGrid1_DI, iLevel1_I, &
            nGridOut1, Weight_I(1:2), iOrder1_I)
       iOrder_I(1:2) = iOrder_I(iOrder1_I)

       !\
       ! Apply weight for interpolation along another axis
       !/
       Weight_I(1:nGridOut1) =  Weight_I(1:nGridOut1) * (1 - Aux_D(x_))
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
       XyzGrid1_DI(x_,1) = XyzGrid_DI(y_,iOrder_I(nGridOut+1))
       XyzGrid1_DI(x_,2) = XyzGrid_DI(y_,iOrder_I(nGridOut+2))
       iLevel1_I(     1) = iLevel_I(     iOrder_I(nGridOut+1))
       iLevel1_I(     2) = iLevel_I(     iOrder_I(nGridOut+2))
       call interpolate_amr_grid1(&
            Xyz1_D,  XyzGrid1_DI, iLevel1_I, &
            nGridOut1, Weight_I(nGridOut+1:nGridOut+2), iOrder1_I)

       iOrder_I(nGridOut+1:nGridOut+2) = iOrder_I(nGridOut+iOrder1_I)
       !\
       ! Apply weight for interpolation along another axis
       !/
       Weight_I(nGridOut + 1:nGridOut + nGridOut1) =  &
            Weight_I(nGridOut + 1:nGridOut + nGridOut1)*Aux_D(x_)

       !\
       ! May need to remove a grid point once behind the boundary
       !/
       nGridOut = nGridOut + nGridOut1
       DoStencilFix = .false.
       RETURN
    end if

    !\
    ! Check the ends of the resolution interfaces. Near the resolution 
    ! interface endpoint the stencil may need to be reevaluated. Check 
    ! fine edges along x_ direction
    !/
    if(abs(XyzGrid_DI(y_,2) - XyzGrid_DI(y_,1)) < &
         dXyzSmall_D(y_).and.all(iLevel_I(1:2)==1))then

       if(abs(XyzGrid_DI(x_,3) - 0.50*(XyzGrid_DI(x_,2) + XyzGrid_DI(x_,1)))&
            < dXyzSmall_D(x_))then
          !\                    C (?)                                  C
          ! The configuration  F F   Check if Xyz belongs to triangle F F
          !/
          iOrder_I = (/2,3,1,4/)
          !\
          ! Interpolate on the triangle X1,X2,X3 (F-C-F)
          ! vector X1-X3 is parallel to iAxis and directed toward
          ! 4th point of the stencil
          !   X3(F) - X1(F)
          !      \   /
          !      X2(C)  X4(?)
          !/
          call interpolate_triangle(iAxisFF=x_)
          RETURN

       elseif(&
            abs(XyzGrid_DI(x_,4) - 0.50*(XyzGrid_DI(x_,2) + XyzGrid_DI(x_,1)))&
            < dXyzSmall_D(x_))then
          !\                (?) C                                          C
          ! The configuration  F F   Check if Xyz belongs to the triangle F F
          !/
          iOrder_I = (/1,4,2,3/)
          !\
          ! Interpolate on the triangle X1,X2,X3 (F-C-F)
          ! vector X1-X3 is parallel to iAxis and directed 
          ! toward 4th point of the stencil
          !   X3(F) - X1(F)
          !      \   /
          !      X2(C)  X4(?)
          !/
          call interpolate_triangle(iAxisFF=x_)
          RETURN
       end if
    end if

    if(abs(XyzGrid_DI(y_,4) - XyzGrid_DI(y_,3)) < &
         dXyzSmall_D(y_).and.all(iLevel_I(3:4)==1))then

       if(abs(XyzGrid_DI(x_,1) - 0.50*(XyzGrid_DI(x_,3) + XyzGrid_DI(x_,4)))&
            <dXyzSmall_D(x_))then
          !\                   F F                                        F F
          ! The configuration   C (?)    Check if Xyz belongs to triangle  C
          !/
          iOrder_I = (/4,1,3,2/)
          !\
          ! Interpolate on the triangle X1,X2,X3 (F-C-F)
          ! vector X1-X3 is parallel to iAxis and directed 
          ! toward 4th point of the stencil
          !   X3(F) - X1(F)
          !      \   /
          !      X2(C)  X4(?)
          !/
          call interpolate_triangle(iAxisFF=x_)
          RETURN

       elseif(&
            abs(XyzGrid_DI(x_,2) - 0.50*(XyzGrid_DI(x_,3) + XyzGrid_DI(x_,4)))&
            < dXyzSmall_D(x_))then
          !\                     F F                                     F F
          ! The configuration (?) C    Check if Xyz belongs to triangle   C
          !/
          iOrder_I = (/3,2,4,1/)
          !\
          ! Interpolate on the triangle X1,X2,X3 (F-C-F)
          ! vector X1-X3 is parallel to iAxis and directed 
          ! toward 4th point of the stencil
          !   X3(F) - X1(F)
          !      \   /
          !      X2(C)  X4(?)
          !/
          call interpolate_triangle(iAxisFF=x_)
          RETURN
       end if
    end if

    !\
    ! Check the ends of the resolution interfaces. Near the resolution 
    ! interface endpoint the stencil may need to be reevaluated. Check 
    ! fine edges along y_ direction
    !/
    if(abs(XyzGrid_DI(x_,3) - XyzGrid_DI(x_,1)) < dXyzSmall_D(x_)&
         .and.all(iLevel_I((/1,3/))==Fine_))then
       if(abs(XyzGrid_DI(y_,2) - &
            0.50*(XyzGrid_DI(y_,3) + XyzGrid_DI(y_,1)))&
            <dXyzSmall_D(y_))then
          !                     (?)
          !\                   F                                       F
          ! The configuration    C  Check if Xyz belongs to triangle C
          !/                   F                                       F
          iOrder_I = (/3,2,1,4/)
          !\
          ! Interpolate on the triangle X1,X2,X3 (F-C-F)
          ! vector X1-X3 is parallel to iAxis and directed toward 
          ! 4th point of the stencil
          !   X3(F) - X1(F)
          !      \   /
          !      X2(C)  X4(?)
          !/
          call interpolate_triangle(iAxisFF=y_)
          RETURN

       elseif(&
            abs(XyzGrid_DI(y_,4) - &
            0.50*(XyzGrid_DI(y_,3) + XyzGrid_DI(y_,1)))&
            < dXyzSmall_D(y_))then
          !\                   F                                         F
          ! The configuration    C  Check if Xyz belongs to triangle   C
          !/                   F                                         F
          !                     (?)
          iOrder_I = (/1,4,3,2/)
          !\
          ! Interpolate on the triangle X1,X2,X3 (F-C-F)
          ! vector X1-X3 is parallel to iAxis and directed 
          ! toward 4th point of the stencil
          !   X3(F) - X1(F)
          !      \   /
          !      X2(C)  X4(?)
          !/
          call interpolate_triangle(iAxisFF=y_)
          RETURN
       end if
    end if

    if(abs(XyzGrid_DI(x_,4) - XyzGrid_DI(x_,2)) < dXyzSmall_D(x_)&
         .and.all(iLevel_I((/2,4/))==1))then

       if(abs(XyzGrid_DI(y_,1) - 0.50*(XyzGrid_DI(y_,2) + XyzGrid_DI(y_,4)))&
            <dXyzSmall_D(y_))then
          !                    (?)
          !\                       F                                        F
          ! The configuration   C       Check if Xyz belongs to triangle  C
          !/                       F                                        F
          iOrder_I = (/4,1,2,3/)
          !\
          ! Interpolate on the triangle X1,X2,X3 (F-C-F)
          ! vector X1-X3 is parallel to iAxis and directed 
          ! toward 4th point of the stencil
          !   X3(F) - X1(F)
          !      \   /
          !      X2(C)  X4(?)
          !/
          call interpolate_triangle(iAxisFF=y_)
          RETURN

       elseif(abs(XyzGrid_DI(y_,3) - &
            0.50*(XyzGrid_DI(y_,2) + XyzGrid_DI(y_,4)))&
            <dXyzSmall_D(y_))then
          !\                       F                                        F
          ! The configuration    C    Check if Xyz belongs to the triangle C
          !/                       F                                        F
          !                     (?)
          iOrder_I = (/2,3,4,1/)
          !\
          ! Interpolate on the triangle X1,X2,X3 (F-C-F)
          ! vector X1-X3 is parallel to iAxis and directed 
          ! toward 4th point of the stencil
          !   X3(F) - X1(F)
          !      \   /
          !      X2(C)  X4(?)
          !/
          call interpolate_triangle(iAxisFF=y_)
          RETURN
       end if
    end if

100 continue


    DoStencilFix = .false.

    if(count(iLevel_I==Fine_)==1.and.count(iLevel_I==Coarse_)==3)then
       iLoc = maxloc(iLevel_I,DIM=1)
       select case(iLoc)
       case(1)
          !   C---C
          !    \ /|
          !     F-C
          iOrder_I = (/1,4,2,3/)

       case(2)
          !   C---C
          !   |\ /
          !   C-F
          iOrder_I = (/2,3,1,4/)

       case(3)
          !     F-C
          !    / \|
          !    C--C
          iOrder_I = (/2,3,1,4/)

       case(4)
          !   C-F
          !   |/\
          !   C--C
          iOrder_I = (/1,4,2,3/)

       end select

    elseif(count(iLevel_I==Fine_)==2.and.count(iLevel_I==Coarse_)==2)then
       if(iLevel_I(1)==Fine_)then
          !   C-F
          !   |/\
          !   F--C
          iOrder_I = (/1,4,2,3/)

       else
          !   F---C
          !   |\ /
          !   C-F
          iOrder_I = (/2,3,1,4/)

       end if


    elseif(count(iLevel_I==Fine_)==3.and.count(iLevel_I==Coarse_)==1)then
       iLoc = minloc(iLevel_I,DIM=1)
       select case(iLoc)
       case(1)
          !   F---F
          !    \ \|
          !     C-F
          iOrder_I = (/2,3,1,4/)

       case(2)
          !   F---F
          !   | //
          !   F-C
          iOrder_I = (/1,4,2,3/)

       case(3)
          !     C-F
          !    / /|
          !    F--F
          iOrder_I = (/1,4,2,3/)

       case(4)
          !   F-C
          !   |\\
          !   F--F
          iOrder_I = (/2,3,1,4/)

       end select
    end if

    !\
    !Points 1 and 2 are on the shared side of the triangles, 3 and 4 are off
    !/
    call triangulate

  contains
    real function cross_product(a_D, b_D)
      real,dimension(nDim),intent(in) :: a_D, b_D
      !----------
      cross_product = a_D(x_)* b_D(y_) - a_D(y_)*b_D(x_)
    end function cross_product
    !=======
    subroutine triangulate
      !\
      !Points 1 and 2 are on the shared side of the triangles, 3 and 4 are off
      !/
      real, dimension(nDim) :: X1_D, X2_D, X3_D, X4_D
      real :: Alpha2, Alpha3
      !-------
      X1_D=XyzGrid_DI(:,iOrder_I(1)); X2_D=XyzGrid_DI(:,iOrder_I(2))
      X3_D=XyzGrid_DI(:,iOrder_I(3)); X4_D=XyzGrid_DI(:,iOrder_I(4))
      Alpha3 = cross_product(Xyz_D - X1_D,X2_D - X1_D)/&
           cross_product(X3_D - X1_D,X2_D - X1_D)
      if(Alpha3==0.0)then
         nGridOut = 2
         Alpha2 = cross_product(X3_D-X1_D,Xyz_D - X1_D)/&
              cross_product(X3_D - X1_D,X2_D - X1_D)
         Weight_I(2) = Alpha2
      elseif(Alpha3 > 0.0)then
         nGridOut = 3
         Alpha2 = cross_product(X3_D-X1_D,Xyz_D - X1_D)/&
              cross_product(X3_D - X1_D,X2_D - X1_D)
         Weight_I(2) = Alpha2
         Weight_I(3) = Alpha3
      else
         nGridOut = 3
         Alpha2 = cross_product(X4_D - X1_D,Xyz_D - X1_D)/&
              cross_product(X4_D - X1_D,X2_D - X1_D)
         Alpha3 = cross_product(Xyz_D - X1_D,X2_D - X1_D)/&
              cross_product(X4_D - X1_D, X2_D-X1_D)
         iOrder_I(3:4) = iOrder_I((/4,3/))
         Weight_I(2) = Alpha2
         Weight_I(3) = Alpha3
      end if

      Weight_I(1) = 1 - sum(Weight_I(2:nGridOut))

    end subroutine triangulate
    !===============================================================
    subroutine interpolate_triangle(iAxisFF)
      !\
      ! Interpolate on the triangle X1,X2,X3 (F-C-F)
      ! vector X1-X3 is parallel to iAxis and directed 
      ! toward 4th point of the stencil
      !   X3(F) - X1(F)
      !      \   /
      !      X2(C)  X4(?)
      !/
      integer, intent(in) :: iAxisFF
      integer             :: iAxisPerp
      real, dimension(nDim) :: X1_D, X2_D, X3_D
      !-------
      !Vertexes
      X1_D=XyzGrid_DI(:,iOrder_I(1)); X2_D=XyzGrid_DI(:,iOrder_I(2))
      X3_D=XyzGrid_DI(:,iOrder_I(3))

      Weight_I(3) = cross_product(Xyz_D - X1_D, X2_D  - X1_D)/&
           cross_product(X3_D - X1_D, X2_D-X1_D)
      Weight_I(2) = cross_product(X3_D  - X1_D, Xyz_D - X1_D)/&
           cross_product(X3_D - X1_D, X2_D-X1_D)
      if(Weight_I(3)==0.0)then
         nGridOut = 2
      else
         nGridOut = 3
      end if
      Weight_I(1) = 1 - sum(Weight_I(2:3))
      DoStencilFix = any(Weight_I(1:3)<0)
      if(DoStencilFix)then
         iAxisPerp = 3 - iAxisFF
         XyzStencil_D(iAxisFF)  =  &
              X1_D(iAxisFF) + X2_D(iAxisFF) - X3_D(iAxisFF)
         XyzStencil_D(iAxisPerp) = &
              (X1_D(iAxisPerp) + X2_D(iAxisPerp) + X3_D(iAxisPerp))/3
      end if
    end subroutine interpolate_triangle
  end subroutine interpolate_amr_grid2

  !==================================================================
  subroutine  interpolate_amr_grid3(&
       Xyz_D, XyzGrid_DI, iLevel_I, &
       nGridOut, Weight_I, iOrder_I,&
       DoStencilFix, XyzStencil_D)

    integer,parameter :: nGrid = 8, nDim = 3


    character(LEN=*),parameter:: NameSub='interpolate_amr_grid3'

    !\
    !Input parameters
    !/
    !\
    !The location at which to interpolate the data
    !/
    real   , intent(in) :: Xyz_D(nDim) 

    !Grid point coordinates !3 coordinate, 8 points
    real  ,   intent(in) :: XyzGrid_DI(nDim,nGrid) 

    !The refinement level at each grid point. By one higher level
    ! of refinement assumes the cell size reduced by a factor of 0.5
    integer,  intent(inout) :: iLevel_I(nGrid)

    !\
    !Output parameters
    !/
    !The number of grid points to be ultimately included 
    !into the interpolation stencil
    !If nGridOut < nGridIn, only the first nGridOut lines in the 
    !output arrays are meaningful.
    integer, intent(out) :: nGridOut

    !The weight coefficients array.
    real   , intent(out) :: Weight_I(nGrid)

    !Order(numbers) of grid points used for the interpolation
    integer, intent(out) :: iOrder_I(nGrid)

    !\
    ! The logical is true, if (for some resolution combinations) the 
    ! basic stencil for the point at which to inpertpolate is not 
    ! applicable
    !/
    logical, intent(out):: DoStencilFix
    !\
    ! If DoStencilFix==.true., the subroutine provides the 
    ! XyzStencil_D to be used to construct stencil, not to interpolate!
    !/
    real, intent(out) :: XyzStencil_D(nDim)

    !\
    !Loop variables
    !/
    integer :: iGrid , jGrid, iDir
    !\
    ! Minimum refinement level used to set iLevel_I=0 in coarse points
    !/
    integer:: iLevelMin

    !\
    ! To find location of fine or coarse points
    !/ 
    integer :: iLoc

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
    real    :: dXyzInv_D(  nDim)                
    real    :: Aux_D(nDim), AuxCoarse, AuxFine   ! Misc

    !\
    ! To call two-dimensional interpolation
    !/
    real    :: Xyz2_D(x_:y_), XyzGrid2_DI(x_:y_,4), XyzStencil2_D(x_:y_)
    integer :: nGridOut2, iOrder2_I(4),iLevel2_I(4)
    logical :: DoStencilFix2

    !\
    ! To call subroutine pyramid
    !/
    integer, parameter:: Rectangular_=1, Trapezoidal_=2
    !-------------

    !\
    ! Make iLevel=0 for coarser grids and iLevel=1 for finer grids
    !/
    iLevelMin = minval(iLevel_I, MASK=iLevel_I/=BehindTheBoundary_)
    where(iLevel_I/=BehindTheBoundary_)iLevel_I = iLevel_I - iLevelMin


    !\
    ! Check if the grid points are out of the grid boundary
    !/
    do iDir = 1, nDim
       if(all(iLevel_I(iOppositeFace_IDI(:,iDir,1))== BehindTheBoundary_).or.&
            all(XyzGrid_DI(iDir,iFace_IDI(:,iDir,1)) == Xyz_D(iDir)) )then

          !\
          ! Upper face iDir should be removed
          !/
          iOrder_I(1:4) = iFace_IDI(        :,iDir,1)
          iOrder_I(5:8) = iOppositeFace_IDI(:,iDir,1)

       elseif(all(iLevel_I(iFace_IDI(:,iDir,1))== BehindTheBoundary_))then

          !\
          ! Lower Face iDir should be removed
          !/

          iOrder_I(1:4) = iOppositeFace_IDI(        :,iDir,1)
          iOrder_I(5:8) = iFace_IDI(:,iDir,1)
       else
          CYCLE
       end if
       !\
       ! iDir, 1 + mod(iDir,3), 1 + mod(iDir+1,3) is 1,2,3 or 2,3,1 or 3,1,2
       !/
       Xyz2_D(x_) = Xyz_D(1 + mod(iDir,3))
       Xyz2_D(y_) = Xyz_D(1 + mod(iDir + 1, 3))
       DoStencilFix2 = .false.
       XyzGrid2_DI(x_,1:4) = XyzGrid_DI(1 + mod(iDir    ,3),iOrder_I(1:4))
       XyzGrid2_DI(y_,1:4) = XyzGrid_DI(1 + mod(iDir + 1,3),iOrder_I(1:4))

       iLevel2_I(1:4) = iLevel_I(iOrder_I(1:4))
       call interpolate_amr_grid2(&
            Xyz2_D,  XyzGrid2_DI, iLevel2_I, &
            nGridOut, Weight_I(1:4), iOrder2_I, DoStencilFix, XyzStencil2_D)

       if(DoStencilFix)then
          XyzStencil_D(iDir)    = Xyz_D(iDir)
          XyzStencil_D(1 + mod(iDir  ,3)) = XyzStencil2_D(x_)
          XyzStencil_D(1 + mod(iDir+1,3)) = XyzStencil2_D(y_)
       end if

       iOrder_I(1:4) = iOrder_I(iOrder2_I)
       RETURN
    end do

    !\
    ! Calculate the stencil size
    !/
    do iDir = 1, nDim
       dXyzSmall_D(iDir) = 0.10*(&
            maxval(XyzGrid_DI(iDir,:)) -  minval(XyzGrid_DI(iDir,:)))
       dXyzInv_D(iDir)   = 0.10/dXyzSmall_D(iDir)
    end do

    !\
    ! First conbinations of the resolution interfaces: no resolution 
    ! interface.
    !/

    if(all(iLevel_I==0))then
       !\
       ! No refinement         
       !                       
       !/
       iOrder_I = (/1,2,3,4,5,6,7,8/)
       Aux_D = (Xyz_D - XyzGrid_DI(:,1))*dXyzInv_D
       nGridOut = 8
       Weight_I(1) = (1 - Aux_D(1))*(1 - Aux_D(2))*(1-Aux_D(3))
       Weight_I(2) =      Aux_D(1) *(1 - Aux_D(2))*(1-Aux_D(3))
       Weight_I(3) = (1 - Aux_D(1))*     Aux_D(2) *(1-Aux_D(3))
       Weight_I(4) =      Aux_D(1) *     Aux_D(2) *(1-Aux_D(3))
       Weight_I(5) = (1 - Aux_D(1))*(1 - Aux_D(2))*   Aux_D(3)
       Weight_I(6) =      Aux_D(1) *(1 - Aux_D(2))*   Aux_D(3)
       Weight_I(7) = (1 - Aux_D(1))*     Aux_D(2) *   Aux_D(3)
       Weight_I(8) =      Aux_D(1) *     Aux_D(2) *   Aux_D(3)
       DoStencilFix = .false.
       RETURN
    end if

    !\
    ! Opposite faces going along resolution interface
    !/
    do iDir = 1,3
       !Face and opposite face are planes
       if(all(abs(&
            XyzGrid_DI(iDir,iFace_IDI(1,iDir,1)) - &
            XyzGrid_DI(iDir,iFace_IDI(:,iDir,1))) &
            < dXyzSmall_D(iDir))&
            .and.all(abs(&
            XyzGrid_DI(iDir,iOppositeFace_IDI(1,iDir,1)) - &
            XyzGrid_DI(iDir,iOppositeFace_IDI(:,iDir,1)))&
            < dXyzSmall_D(iDir)))then
          !\
          ! Faces going along plane iDir (which is the direction of normal)    
          !                               
          !/
          iOrder_I(1:4) = iFace_IDI(:,iDir,1)
          iOrder_I(5:8) = iOppositeFace_IDI(:,iDir,1)

          Xyz2_D(x_) = Xyz_D(1 + mod(iDir    ,3))
          Xyz2_D(y_) = Xyz_D(1 + mod(iDir + 1,3))

          !\
          ! Interpolate along lower face
          !/
          DoStencilFix2 = .false.
          XyzGrid2_DI(x_,1:4) = XyzGrid_DI(1 + mod(iDir    ,3),iOrder_I(1:4))
          XyzGrid2_DI(y_,1:4) = XyzGrid_DI(1 + mod(iDir + 1,3),iOrder_I(1:4))
          iLevel2_I(        1:4) = iLevel_I(        iOrder_I(1:4))
          call interpolate_amr_grid2(&
               Xyz2_D, XyzGrid2_DI, iLevel2_I, &
               nGridOut2, Weight_I(1:4), iOrder2_I, &
               DoStencilFix2, XyzStencil2_D)
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
          DoStencilFix2 = .false.

          XyzGrid2_DI(x_,1:4) = XyzGrid_DI(1 + mod(iDir    ,3),&
               iOrder_I(nGridOut2 + 1:nGridOut2+4))
          XyzGrid2_DI(y_,1:4) = XyzGrid_DI(1 + mod(iDir + 1,3),&
               iOrder_I(nGridOut2 + 1:nGridOut2 + 4))
          iLevel2_I(1:4) = iLevel_I(iOrder_I(nGridOut2 + 1:nGridOut2 + 4))
          call interpolate_amr_grid2(&
               Xyz2_D, XyzGrid2_DI, iLevel2_I, &
               nGridOut2, Weight_I(nGridOut+1:nGridOut+4), iOrder2_I, &
               DoStencilFix2, XyzStencil2_D)
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
          DoStencilFix = .false.
          RETURN
       end if
    end do


    !\
    ! Edges going along resolution edge
    ! Opposite "faces" have the same geometry in projection
    ! term "face" is used to refer to 4 upper-, left-, 
    ! front-most etc. grid points they do not belong to plane
    !/
    do iDir = 1,nDim
       if(  all(abs(XyzGrid_DI(1 + mod(iDir    ,3),iFace_IDI(:,iDir,1)) - &
            XyzGrid_DI(1 + mod(iDir    ,3),iOppositeFace_IDI(:,iDir,1))) < &
            dXyzSmall_D(1 + mod(iDir    ,3)).and.&
            abs(XyzGrid_DI(1 + mod(iDir + 1,3),iFace_IDI(:,iDir,1)) - &
            XyzGrid_DI(1 + mod(iDir + 1,3),iOppositeFace_IDI(:,iDir,1))) < &
            dXyzSmall_D(1 + mod(iDir + 1,3)) ) )then
          !\
          ! "Faces" have the same geometry in projection on iDir-plane
          !/
          iOrder_I(1:4) = iFace_IDI(:,iDir,1)
          iOrder_I(5:8) = iOppositeFace_IDI(:,iDir,1)
          Xyz2_D(x_:y_) = Xyz_D((/1 + mod(iDir,3),1 + mod(1 + iDir,3)/))


          !\
          ! Interpolation weight along iDir axis are calculated seperately 
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
          DoStencilFix2 = DoStencilFix      !!??????!!!
          XyzGrid2_DI(x_:y_,1:4) = &
               XyzGrid_DI((/1 + mod(iDir,3),1 + mod(iDir+1,3)/),iOrder_I(1:4))
          iLevel2_I(        1:4) = iLevel_I( iOrder_I(1:4))
          where(iLevel2_I==BehindTheBoundary_)iLevel2_I=Coarse_
          call interpolate_amr_grid2(&
               Xyz2_D, XyzGrid2_DI, iLevel2_I, &
               nGridOut2, Weight_I(1:4), iOrder2_I, &
               DoStencilFix2, XyzStencil2_D)

          if(DoStencilFix2)then
             DoStencilFix = DoStencilFix2
             XyzStencil_D((/1 + mod(iDir,3),1 + mod(iDir+1,3)/)) = &
                  XyzStencil2_D(x_:y_)
             XyzStencil_D(iDir) = Xyz_D(iDir)
             RETURN
          end if

          iOrder_I(1:4) = iOrder_I(iOrder2_I  )
          iOrder_I(5:8) = iOrder_I(iOrder2_I+4)

          !\
          ! May need to remove a grid points from the stencil
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
                if(iLevel_I(iOrder_I(iGrid))==BehindTheBoundary_)then
                   Weight_I(nGridOut2+iGrid) = Weight_I(iGrid)
                   Weight_I(          iGrid) = 0
                elseif(&
                     iLevel_I(iOrder_I(nGridOut2 + iGrid))==BehindTheBoundary_&
                     )then
                   Weight_I(nGridOut2+iGrid) = 0
                   Weight_I(          iGrid) = Weight_I(iGrid)
                else
                   Weight_I(nGridOut2+iGrid) = Weight_I(iGrid)* AuxCoarse
                   Weight_I(          iGrid) = Weight_I(iGrid)*(1 - AuxCoarse)
                end if
             end if
          end do
          !\
          ! Remove points behind the boundary
          !/
          do while(any(iLevel_I(iOrder_I(1:nGridOut))==BehindTheBoundary_))
             if(iLevel_I(iOrder_I(nGridOut))==BehindTheBoundary_)then
                nGridOut = nGridOut - 1
             else
                iLoc = minloc(iLevel_I(iOrder_I(1:nGridOut-1)),DIM=1)
                Weight_I(iLoc) = Weight_I(nGridOut)
                iOrder_I(iLoc) = iOrder_I(nGridOut)
                nGridOut = nGridOut -1
             end if
          end do
          DoStencilFix = .false.
          RETURN
       end if
    end do

    !\
    ! Corner transition junction
    !/
    do iGrid = 1,nGrid;do iDir = 1,nDim
       if( abs(        XyzGrid_DI(1 + mod(iDir    ,3),iGrid) - &
            0.25*sum(XyzGrid_DI(1 + mod(iDir    ,3),&
            iOppositeFace_IDI(:,iDir,iGrid)))) < &
            dXyzSmall_D(1 + mod(iDir    ,3)).and.&
            abs(        XyzGrid_DI(1 + mod(iDir + 1,3),iGrid) - &
            0.25*sum(XyzGrid_DI(1 + mod(iDir + 1,3),&
            iOppositeFace_IDI(:,iDir,iGrid)))) < &
            dXyzSmall_D(1 + mod(iDir + 1,3)).and.&
            (count(iLevel_I==Coarse_)==3))then
          !\
          ! Junction around direction iDir
          !/
          iOrder_I(1:4) = iFace_IDI(:,iDir,iGrid)
          iOrder_I(5:8) = iOppositeFace_IDI(:,iDir,iGrid)
          call corner_transition_junction(iDir)
          RETURN
       end if
    end do;end do

    !\
    ! Corner transition
    !/
    do iGrid = 1,nGrid;do iDir = 1,nDim
       if( count(iLevel_I(iFace_IDI(:,iDir,iGrid))==Fine_)>0 .and.  &
            all(abs(XyzGrid_DI(iDir,iGrid)                         -&
            0.50*(XyzGrid_DI(iDir,iFace_IDI(        :,iDir,iGrid)) +&
            XyzGrid_DI(iDir,iOppositeFace_IDI(:,iDir,iGrid)))*      &
            (iLevel_I(iFace_IDI(:,iDir,iGrid))-Coarse_)            -&
            XyzGrid_DI(iDir,iFace_IDI(         :,iDir,iGrid))*      &
            (Fine_-iLevel_I(iFace_IDI(:,iDir,iGrid)))) < &
            dXyzSmall_D(iDir)))then
          !\
          ! Transition in direction iDir
          !/
          iOrder_I(1:4) = iFace_IDI(:,iDir,iGrid)
          iOrder_I(5:8) = iOppositeFace_IDI(:,iDir,iGrid)
          call interpolate_corner_transition(iDir)
          RETURN
       end if
    end do;end do



    DoStencilFix = .false.
    !---------------------Resolution corners----------------------
    !Among 255 legal combinations
    !(=2**nGrid - illegal "all Fine configurtion")
    !we considered above one all-coarse configuration, 
    !+ 6 resolution faces
    !+ 30 resolution edges, totally 37 left 218
    !
    !
    !\
    ! One configuration resulting in the corner split for five tetrahedra
    ! is occured in many cases (26), however, it is treated in the same way 
    ! in all these cases. If the given grid point is connected by three edges 
    ! with three coarse points, while the across the main diagonal point 
    ! is connected by three edges (or, equivalently, the given point is 
    ! connected by three face diagonals) with three coarse points), 
    ! then all faces are trangulated   
    !/
    !                                26 cases totally 63 left 192
    do iGrid = 1,nGrid
       if(  all(iLevel_I(iEdge_ID(iGrid,:))==Coarse_).and.&
            all(iLevel_I(iFaceDiag_ID(iGrid,:))==Fine_))then
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
               iTetrahedron1_I=(/ iGrid,iEdge_ID(iGrid,x_),&
               iFaceDiag_ID(iGrid,Xz_),iFaceDiag_ID(iGrid,Xy_)/),&
               iTetrahedron2_I=(/ iGrid,iEdge_ID(iGrid,y_),&
               iFaceDiag_ID(iGrid,Yz_),iFaceDiag_ID(iGrid,Xy_)/),&
               iTetrahedron3_I=(/ iGrid,iEdge_ID(iGrid,z_),&
               iFaceDiag_ID(iGrid,Yz_),iFaceDiag_ID(iGrid,Xz_)/),&
               iTetrahedron4_I=(/ iGrid,iFaceDiag_ID(iGrid,Yz_),&
               iFaceDiag_ID(iGrid,Xz_),iFaceDiag_ID(iGrid,Xy_)/),&
               iTetrahedron5_I=(/&
               iMainDiag_I(iGrid),iFaceDiag_ID(iGrid,Yz_),&
               iFaceDiag_ID(iGrid,Xz_),iFaceDiag_ID(iGrid,Xy_)/))
          RETURN
       end if
    end do


    select case( count(iLevel_I==Fine_))
    case(1)                          ! 8 cases , totally 71 left 184
       !Find the only fine grid
       iLoc = maxloc(iLevel_I,DIM=1)
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
            iRectangular1_I=(/&
            iOppositeFace_IDI(1,x_,iLoc),&
            iOppositeFace_IDI(2,x_,iLoc),&
            iOppositeFace_IDI(3,x_,iLoc),&
            iOppositeFace_IDI(4,x_,iLoc),&
            iLoc/)                      ,&
            iRectangular2_I=(/&
            iOppositeFace_IDI(1,y_,iLoc),&
            iOppositeFace_IDI(2,y_,iLoc),&
            iOppositeFace_IDI(3,y_,iLoc),&
            iOppositeFace_IDI(4,y_,iLoc),&
            iLoc/)                      ,&
            iRectangular3_I=(/&
            iOppositeFace_IDI(1,z_,iLoc),&
            iOppositeFace_IDI(2,z_,iLoc),&
            iOppositeFace_IDI(3,z_,iLoc),&
            iOppositeFace_IDI(4,z_,iLoc),&
            iLoc/))
       RETURN
    case(7)                                 ! 8 cases , totally 79 left 176
       !Find the only coarse grid
       iLoc = minloc(iLevel_I,DIM=1)
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
            iTetrahedron1_I=&
            (/iLoc,iEdge_ID(iLoc,x_),iEdge_ID(iLoc,y_),iEdge_ID(iLoc,z_)/))
       if(nGridOut>1)RETURN
       call parallel_rays(&
            Dir_D=XyzGrid_DI(:,iMainDiag_I(iLoc)) - XyzGrid_DI(:,iLoc), &
            iURectangle1_I=iFace_IDI(:,x_,iMainDiag_I(iLoc)),&
            iURectangle2_I=iFace_IDI(:,y_,iMainDiag_I(iLoc)),&
            iURectangle3_I=iFace_IDI(:,z_,iMainDiag_I(iLoc)),&
            iDTriangle1_I=&
            (/iEdge_ID(iLoc,x_),iEdge_ID(iLoc,y_),iEdge_ID(iLoc,z_)/),&
            iDTriangle2_I=&
            (/iEdge_ID(iLoc,x_),iEdge_ID(iLoc,y_),iFaceDiag_ID(iLoc,Xy_)/),&
            iDTriangle3_I=&
            (/iEdge_ID(iLoc,x_),iEdge_ID(iLoc,z_),iFaceDiag_ID(iLoc,Xz_)/),&
            iDTriangle4_I=&
            (/iEdge_ID(iLoc,y_),iEdge_ID(iLoc,z_),iFaceDiag_ID(iLoc,Yz_)/))
       RETURN
    case(2)
       !Find the fine grid
       iLoc = maxloc(iLevel_I,DIM=1)
       if(iLevel_I(iMainDiag_I(iLoc))==Fine_)then 

          ! 4 cases   totally 83 left 172

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
          iGrid = iMainDiag_I(iLoc) !Fine point across the main diagonal
          call pyramids(&
               iTetrahedron1_I=(/iLoc ,  &
               iEdge_ID(    iLoc,x_ ) ,  &
               iFaceDiag_ID(iLoc,Xy_) ,  &
               iGrid/),                    &
               iTetrahedron2_I=(/iLoc ,  &
               iEdge_ID(    iLoc,x_ ) ,  &
               iFaceDiag_ID(iLoc,Xz_) ,  &
               iGrid/),                    &
               iTetrahedron3_I=(/iLoc ,  &
               iEdge_ID(    iLoc,y_ ) ,  &
               iFaceDiag_ID(iLoc,Xy_) ,  &
               iGrid/),                    &
               iTetrahedron4_I=(/iLoc ,  &
               iEdge_ID(    iLoc,y_ ) ,  &
               iFaceDiag_ID(iLoc,Yz_) ,  &
               iGrid/),                    &
               iTetrahedron5_I=(/iLoc ,  &
               iEdge_ID(    iLoc,z_ ) ,  &
               iFaceDiag_ID(iLoc,Xz_) ,  &
               iGrid/),                    &
               iTetrahedron6_I=(/iLoc ,  &
               iEdge_ID(    iLoc,z_ ) ,  &
               iFaceDiag_ID(iLoc,Yz_) ,  &
               iGrid/))
          RETURN
       else                      ! 12 cases  totally 95  left 160
          do iDir = Yz_,Xy_
             if(iLevel_I(iFaceDiag_ID(iLoc, iDir))==Fine_)then
                iOrder_I(1:4) = iFace_IDI(:,iDir,iLoc)
                iOrder_I(5:8) = iOppositeFace_IDI(:,iDir,iLoc)
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
                     iTetrahedron1_I=(/&
                     iOrder_I(1), iOrder_I(4),iOrder_I(2), iOrder_I(6)/),&
                     iTetrahedron2_I=(/&
                     iOrder_I(1), iOrder_I(4),iOrder_I(3), iOrder_I(7)/))
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
                     0.50*(XyzGrid_DI(:,iOrder_I(1))+XyzGrid_DI(:,iOrder_I(4)))
                !\
                ! Center of the upper face
                !/
                XyzMax_D = 0.50*&
                     (XyzGrid_DI(:,iOrder_I(5))+XyzGrid_DI(:,iOrder_I(8)))
                call parallel_rays(Dir_D=XyzMax_D - XyzMin_D,&
                     iURectangle1_I=iOrder_I(5:8),&
                     iDTriangle1_I=(/iOrder_I(1),iOrder_I(4),iOrder_I(6)/),&
                     iDTriangle2_I=(/iOrder_I(1),iOrder_I(4),iOrder_I(7)/),&
                     iDTriangle3_I=(/iOrder_I(4),iOrder_I(6),iOrder_I(8)/),&
                     iDTriangle4_I=(/iOrder_I(4),iOrder_I(7),iOrder_I(8)/),&
                     iDTriangle5_I=(/iOrder_I(1),iOrder_I(5),iOrder_I(7)/),&
                     iDTriangle6_I=(/iOrder_I(1),iOrder_I(5),iOrder_I(6)/))
                RETURN
             end if
          end do
       end if
    case(6)
       !Find the coarse point
       iLoc = minloc(iLevel_I,DIM=1)
       if(iLevel_I(iMainDiag_I(iLoc))==Coarse_)then  
          !\
          ! 4 cases totally 99 left 156
          !/
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
          iGrid = iMainDiag_I(iLoc)
          call pyramids(&
               iTetrahedron1_I=(/iLoc,&
               iEdge_ID(iLoc,x_),iEdge_ID(iLoc,y_),iEdge_ID(iLoc,z_)/),&
               iTetrahedron2_I=(/iGrid,&
               iEdge_ID(iGrid,x_),iEdge_ID(iGrid,y_),iEdge_ID(iGrid,z_)/))
          if(nGridOut>1)RETURN
          call parallel_rays(&
               Dir_D=XyzGrid_DI(:,iGrid) - XyzGrid_DI(:,iLoc), &
               iUTriangle1_I=(/iEdge_ID(iGrid,x_),&
               iEdge_ID(iGrid,y_),iEdge_ID(iGrid ,z_)/),&
               iUTriangle2_I=(/iEdge_ID(iGrid,x_),&
               iEdge_ID(iGrid,y_),iFaceDiag_ID(iGrid,Xy_)/),&
               iUTriangle3_I=(/iEdge_ID(iGrid,x_),&
               iEdge_ID(iGrid,z_),iFaceDiag_ID(iGrid,Xz_)/),&
               iUTriangle4_I=(/iEdge_ID(iGrid,y_),&
               iEdge_ID(iGrid,z_),iFaceDiag_ID(iGrid,Yz_)/),&
               iDTriangle1_I=(/iEdge_ID(iLoc,x_),&
               iEdge_ID(iLoc,y_),iEdge_ID(iLoc,z_)/),&
               iDTriangle2_I=(/iEdge_ID(iLoc,x_),&
               iEdge_ID(iLoc,y_),iFaceDiag_ID(iLoc,Xy_)/),&
               iDTriangle3_I=(/iEdge_ID(iLoc,x_),&
               iEdge_ID(iLoc,z_),iFaceDiag_ID(iLoc,Xz_)/),&
               iDTriangle4_I=(/iEdge_ID(iLoc,y_),&
               iEdge_ID(iLoc,z_),iFaceDiag_ID(iLoc,Yz_)/))
          RETURN
       else                               ! 12 cases  totally 111  left 144
          do iDir = Yz_,Xy_
             if(iLevel_I(iFaceDiag_ID(iLoc, iDir))==Coarse_)then
                iOrder_I(1:4) = iFace_IDI(:,iDir,iLoc)
                iOrder_I(5:8) = iOppositeFace_IDI(:,iDir,iLoc)
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
                     iTetrahedron1_I=(/&
                     iOrder_I(1), iOrder_I(2),iOrder_I(3), iOrder_I(5)/),&
                     iTetrahedron2_I=(/&
                     iOrder_I(4), iOrder_I(2),iOrder_I(3), iOrder_I(8)/))
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
                        iURectangle1_I=iOrder_I(5:8),&
                        iDTriangle1_I=(/iOrder_I(2),iOrder_I(3),iOrder_I(5)/),&
                        iDTriangle2_I=(/iOrder_I(2),iOrder_I(3),iOrder_I(8)/))
                end if
                RETURN
             end if
          end do
       end if

    case(3)
       !                   24 cases   totally 135 left 120
       do iGrid = 1, nGrid
          if(iLevel_I(iGrid)==Fine_)then
             if(iLevel_I(iMainDiag_I(iGrid))==Fine_)then
                do iDir = Yz_,Xy_
                   if(iLevel_I(iFaceDiag_ID(iGrid, iDir))==Fine_)then
                      iOrder_I(1:4) = iFace_IDI(:,iDir,iGrid)
                      iOrder_I(5:8) = iOppositeFace_IDI(:,iDir,iGrid)
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
                           iTetrahedron1_I=(/&
                           iOrder_I(1),iOrder_I(5),iOrder_I(8),iOrder_I(6)/),&
                           iTetrahedron2_I=(/&
                           iOrder_I(1),iOrder_I(5),iOrder_I(8),iOrder_I(7)/),&
                           iTrapezoidal1_I=(/&
                           iOrder_I(2),iOrder_I(6),iOrder_I(4),iOrder_I(8),&
                           iOrder_I(1)/),&
                           iTrapezoidal2_I=(/&
                           iOrder_I(3),iOrder_I(7),iOrder_I(4),iOrder_I(8),&
                           iOrder_I(1)/))
                      RETURN
                   end if
                end do
             else
                !\
                ! Three fine cells in one plane   
                !/
                !\
                ! 24 cases, totally  159 left 96
                !/
                do iDir = 1,nDim
                   if(all(iLevel_I(iFace_IDI(1:3,iDir,iGrid))==Fine_))then
                      iOrder_I(1:4) = iFace_IDI(:,iDir,iGrid)
                      iOrder_I(5:8) = iOppositeFace_IDI(:,iDir,iGrid)
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
                           iTetrahedron1_I=&
                           (/iOrder_I(2),iOrder_I(3),iOrder_I(8),iOrder_I(4)/))
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
                           iDTriangle1_I=(/&
                           iOrder_I(1),iOrder_I(2),iOrder_I(3)/),&
                           iDTriangle2_I=(/&
                           iOrder_I(8),iOrder_I(2),iOrder_I(3)/),&
                           iDTriangle3_I=(/&
                           iOrder_I(8),iOrder_I(7),iOrder_I(3)/),&
                           iDTriangle4_I=(/&
                           iOrder_I(8),iOrder_I(2),iOrder_I(6)/),&
                           iURectangle1_I=(/&
                           iOrder_I(5),iOrder_I(6),&
                           iOrder_I(7),iOrder_I(8)/))    
                      RETURN
                   end if
                end do
             end if
          end if
       end do
    case(5)
       !               24 cases, totally  183 left 72
       do iGrid = 1, nGrid
          if(iLevel_I(iGrid)==Coarse_)then
             if(iLevel_I(iMainDiag_I(iGrid))==Coarse_)then
                do iDir = Yz_,Xy_
                   if(iLevel_I(iFaceDiag_ID(iGrid, iDir))==Coarse_)then
                      iOrder_I(1:4) = iFace_IDI(:,iDir,iGrid)
                      iOrder_I(5:8) = iOppositeFace_IDI(:,iDir,iGrid)
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
                      !
                      call pyramids(&
                           iTetrahedron1_I=(/&
                           iOrder_I(1),iOrder_I(2),iOrder_I(3),iOrder_I(5)/),&
                           iRectangular1_I=(/&
                           iOrder_I(2),iOrder_I(3),iOrder_I(6),iOrder_I(7),&
                           iOrder_I(5)/))
                      if(nGridOut> 1)RETURN

                      Xyz2_D(x_:y_) = &
                           Xyz_D((/1 + mod(iDir,3),1 + mod(1 + iDir,3)/))
                      !\
                      ! Interpolation weights along iDir axis are calculated 
                      ! separately for coarse and fine points
                      !/
                      !\
                      ! may need to rearrange iOrder
                      !/
                      if(XyzGrid_DI(iDir,iOrder_I(1))&
                           > XyzGrid_DI(iDIr,iOrder_I(5)))then
                         iOrder_I(1:8) = iOrder_I((/5,6,7,8,1,2,3,4/))
                      end if
                      XyzMin_D(iDir) = minval(&
                           XyzGrid_DI(iDir,:),iLevel_I/=Fine_)
                      AuxCoarse  = (Xyz_D(iDir) - XyzMin_D(iDir))/&
                           (maxval(&
                           XyzGrid_DI(iDir,:),iLevel_I/=Fine_) - &
                           XyzMin_D(iDir))

                      XyzMin_D(iDir) = minval(&
                           XyzGrid_DI(iDir,:),iLevel_I==Fine_)
                      AuxFine  = (Xyz_D(iDir) - XyzMin_D(iDir))/&
                           (maxval(&
                           XyzGrid_DI(iDir,:),iLevel_I==Fine_) - &
                           XyzMin_D(iDir))
                      !\
                      ! Interpolate along lower face
                      !/
                      DoStencilFix2 = .true.
                      XyzGrid2_DI(x_:y_,1:4) = &
                           XyzGrid_DI((/1 + mod(iDir,3),1 + mod(iDir+1,3)/),&
                           iOrder_I(1:4))
                      iLevel2_I(        1:4) = iLevel_I( iOrder_I(1:4))

                      call interpolate_amr_grid2(&
                           Xyz2_D, XyzGrid2_DI, iLevel2_I, &
                           nGridOut2, Weight_I(1:4), iOrder2_I, &
                           DoStencilFix2, XyzStencil2_D)

                      iOrder_I(1:4) = iOrder_I(iOrder2_I  )
                      iOrder_I(5:8) = iOrder_I(iOrder2_I+4)

                      !\
                      ! May need to remove a grid points from the stencil
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
                      do jGrid = 1,nGridOut2
                         if(iLevel_I(iOrder_I(jGrid))==Fine_)then
                            Weight_I(nGridOut2+jGrid) = &
                                 Weight_I(jGrid)*     AuxFine
                            Weight_I(          jGrid) = &
                                 Weight_I(jGrid)*(1-  AuxFine)
                         else
                            Weight_I(nGridOut2+jGrid) = &
                                 Weight_I(jGrid)*     AuxCoarse
                            Weight_I(          jGrid) = &
                                 Weight_I(jGrid)*(1 - AuxCoarse)
                         end if
                      end do
                      RETURN
                   end if
                end do
             end if
             !\
             ! Three coarse cells in one plane 24 cases, totally 207 left 48
             !/
             do iDir = 1,nDim
                if(all(iLevel_I(iFace_IDI(1:3, iDir, iGrid))==Coarse_))then
                   iOrder_I(1:4) = iFace_IDI(:,iDir,iGrid)
                   iOrder_I(5:8) = iOppositeFace_IDI(:,iDir,iGrid)
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
                        iRectangular1_I=&
                        (/iOrder_I(7),iOrder_I(8),iOrder_I(5),iOrder_I(6),&
                        iOrder_I(4)/),&
                        iTrapezoidal1_I=&
                        (/iOrder_I(1),iOrder_I(2),iOrder_I(5),iOrder_I(6),&
                        iOrder_I(4)/),&
                        iTrapezoidal2_I=&
                        (/iOrder_I(1),iOrder_I(3),iOrder_I(5),iOrder_I(7),&
                        iOrder_I(4)/))
                   RETURN
                end if
             end do
          end if
       end do
    case(4)
       do iGrid = 1,nGrid
          if(iLevel_I(iGrid)==Coarse_)then
             do iDir = 1,nDim
                !
                !   24 cases   totally 231  left 24
                !
                if(all(&
                     iLevel_I(iOppositeFace_IDI(2:4,iDir,iGrid))==Coarse_&
                     ))then
                   iOrder_I(1:4) = iFace_IDI(:,iDir,iGrid)
                   iOrder_I(5:8) = iOppositeFace_IDI(:,iDir,iGrid)
                   !\
                   ! face 1:4 has only one Coarse point iOrder_I(1)
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
                   !         
                   !                F3 F4
                   !      Trapezoid C7 C8
                   !                
                   !
                   !Common apex F5
                   call pyramids(&
                        iTrapezoidal1_I=&
                        (/iOrder_I(6),iOrder_I(8),iOrder_I(2),iOrder_I(4),&
                        iOrder_I(5)/),&
                        iTrapezoidal2_I=&
                        (/iOrder_I(7),iOrder_I(8),iOrder_I(3),iOrder_I(4),&
                        iOrder_I(5)/),&
                        iTetrahedron1_I=&
                        (/iOrder_I(2),iOrder_I(3),iOrder_I(4),iOrder_I(5)/),&
                        iTetrahedron2_I=&
                        (/iOrder_I(1),iOrder_I(2),iOrder_I(3),iOrder_I(5)/) )
                   RETURN
                end if
                !
                !\
                ! 24 cases totally 255 left 0
                !/
                if(iLevel_I(iEdge_ID(iGrid,iDir))==Coarse_.and.&
                     iLevel_I(iEdge_ID(iGrid,1 + mod(iDir,3)))==Coarse_.and.&
                     iLevel_I(iEdge_ID(iEdge_ID(iGrid,iDir),&
                     1 + mod(1 + iDir,3)))==Coarse_)then
                   iOrder_I(1:4) = iFace_IDI(:,iDir,iGrid)
                   iOrder_I(5:8) = iOppositeFace_IDI(:,iDir,iGrid)
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
                        iTrapezoidal1_I=&
                        (/iOrder_I(1),iOrder_I(2),iOrder_I(3),iOrder_I(4),&
                        iOrder_I(6)/),& 
                        iTetrahedron1_I=&            !see #1 below 
                        (/iOrder_I(8),iOrder_I(3),iOrder_I(4),iOrder_I(6)/),&
                        iTetrahedron2_I=&            !see #2 below
                        (/iOrder_I(1),iOrder_I(6),iOrder_I(3),iOrder_I(5)/),&
                        iTrapezoidal2_I=&            !see #3 below
                        (/iOrder_I(5),iOrder_I(7),iOrder_I(6),iOrder_I(8),&
                        iOrder_I(3)/))
                   ! #1 The leftover is above subfaces 136 and 364
                   ! #2 The leftover is above subfaces 136 and 368
                   ! #3 The leftover is above subfaces 365 and 368
                   RETURN
                end if
             end do
          end if
       end do

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
              XyzGrid_DI(:,iTetrahedron1_I(1)),&
              XyzGrid_DI(:,iTetrahedron1_I(2)),&
              XyzGrid_DI(:,iTetrahedron1_I(3)),&
              XyzGrid_DI(:,iTetrahedron1_I(4)))
         if(all(Weight_I(1:4).ge.0.0))then
            nGridOut = 4
            iOrder_I(1:4) = iTetrahedron1_I
            RETURN
         elseif(present(iTetrahedron2_I))then
            call tetrahedron(&
                 XyzGrid_DI(:,iTetrahedron2_I(1)),&
                 XyzGrid_DI(:,iTetrahedron2_I(2)),&
                 XyzGrid_DI(:,iTetrahedron2_I(3)),&
                 XyzGrid_DI(:,iTetrahedron2_I(4)))
            if(all(Weight_I(1:4).ge.0.0))then
               nGridOut = 4
               iOrder_I(1:4) = iTetrahedron2_I
               RETURN
            elseif(present(iTetrahedron3_I))then
               call tetrahedron(&
                    XyzGrid_DI(:,iTetrahedron3_I(1)),&
                    XyzGrid_DI(:,iTetrahedron3_I(2)),&
                    XyzGrid_DI(:,iTetrahedron3_I(3)),&
                    XyzGrid_DI(:,iTetrahedron3_I(4)))
               if(all(Weight_I(1:4).ge.0.0))then
                  nGridOut = 4
                  iOrder_I(1:4) = iTetrahedron3_I
                  RETURN
               elseif(present(iTetrahedron4_I))then
                  call tetrahedron(&
                       XyzGrid_DI(:,iTetrahedron4_I(1)),&
                       XyzGrid_DI(:,iTetrahedron4_I(2)),&
                       XyzGrid_DI(:,iTetrahedron4_I(3)),&
                       XyzGrid_DI(:,iTetrahedron4_I(4)))
                  if(all(Weight_I(1:4).ge.0.0))then
                     nGridOut = 4
                     iOrder_I(1:4) = iTetrahedron4_I
                     RETURN
                  elseif(present(iTetrahedron5_I))then
                     call tetrahedron(&
                          XyzGrid_DI(:,iTetrahedron5_I(1)),&
                          XyzGrid_DI(:,iTetrahedron5_I(2)),&
                          XyzGrid_DI(:,iTetrahedron5_I(3)),&
                          XyzGrid_DI(:,iTetrahedron5_I(4)))
                     if(all(Weight_I(1:4).ge.0.0))then
                        nGridOut = 4
                        iOrder_I(1:4) = iTetrahedron5_I
                        RETURN
                     elseif(present(iTetrahedron6_I))then
                        call tetrahedron(&
                             XyzGrid_DI(:,iTetrahedron6_I(1)),&
                             XyzGrid_DI(:,iTetrahedron6_I(2)),&
                             XyzGrid_DI(:,iTetrahedron6_I(3)),&
                             XyzGrid_DI(:,iTetrahedron6_I(4)))
                        if(all(Weight_I(1:4).ge.0.0))then
                           nGridOut = 4
                           iOrder_I(1:4) = iTetrahedron6_I
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
              X1_D=XyzGrid_DI(:,iRectangular1_I(1)),&
              X2_D=XyzGrid_DI(:,iRectangular1_I(2)),&
              X3_D=XyzGrid_DI(:,iRectangular1_I(3)),&
              X4_D=XyzGrid_DI(:,iRectangular1_I(4)),&
              X5_D=XyzGrid_DI(:,iRectangular1_I(5)) )
         if(all(Weight_I(1:5).ge.0.0))then
            nGridOut = 5
            iOrder_I(1:5) = iRectangular1_I
            RETURN
         elseif(present(iRectangular2_I))then
            call pyramid(iBase=Rectangular_,&
                 X1_D=XyzGrid_DI(:,iRectangular2_I(1)),&
                 X2_D=XyzGrid_DI(:,iRectangular2_I(2)),&
                 X3_D=XyzGrid_DI(:,iRectangular2_I(3)),&
                 X4_D=XyzGrid_DI(:,iRectangular2_I(4)),&
                 X5_D=XyzGrid_DI(:,iRectangular2_I(5)) )
            if(all(Weight_I(1:5).ge.0.0))then
               nGridOut = 5
               iOrder_I(1:5) = iRectangular2_I
               RETURN
            elseif(present(iRectangular3_I))then
               call pyramid(iBase=Rectangular_,&
                    X1_D=XyzGrid_DI(:,iRectangular3_I(1)),&
                    X2_D=XyzGrid_DI(:,iRectangular3_I(2)),&
                    X3_D=XyzGrid_DI(:,iRectangular3_I(3)),&
                    X4_D=XyzGrid_DI(:,iRectangular3_I(4)),&
                    X5_D=XyzGrid_DI(:,iRectangular3_I(5)) )
               if(all(Weight_I(1:5).ge.0.0))then
                  nGridOut = 5
                  iOrder_I(1:5) = iRectangular3_I
                  RETURN
               end if          ! 3
            end if             ! 2
         end if                ! 1
      end if                   ! no rectangular
      if(present(iTrapezoidal1_I))then
         call pyramid(iBase=Trapezoidal_,&
              X1_D=XyzGrid_DI(:,iTrapezoidal1_I(1)),&
              X2_D=XyzGrid_DI(:,iTrapezoidal1_I(2)),&
              X3_D=XyzGrid_DI(:,iTrapezoidal1_I(3)),&
              X4_D=XyzGrid_DI(:,iTrapezoidal1_I(4)),&
              X5_D=XyzGrid_DI(:,iTrapezoidal1_I(5)) )
         if(all(Weight_I(1:5).ge.0.0))then
            nGridOut = 5
            iOrder_I(1:5) = iTrapezoidal1_I
            RETURN
         elseif(present(iTrapezoidal2_I))then
            call pyramid(iBase=Trapezoidal_,&
                 X1_D=XyzGrid_DI(:,iTrapezoidal2_I(1)),&
                 X2_D=XyzGrid_DI(:,iTrapezoidal2_I(2)),&
                 X3_D=XyzGrid_DI(:,iTrapezoidal2_I(3)),&
                 X4_D=XyzGrid_DI(:,iTrapezoidal2_I(4)),&
                 X5_D=XyzGrid_DI(:,iTrapezoidal2_I(5)) )
            if(all(Weight_I(1:5).ge.0.0))then
               nGridOut = 5
               iOrder_I(1:5) = iTrapezoidal2_I
               RETURN
            end if             ! 2
         end if                ! 1
      end if                   ! no rectangular

    end subroutine pyramids
    !======================
    subroutine parallel_rays(Dir_D, &
         iDRectangle1_I, iDRectangle2_I, iDRectangle3_I,&
         iDTriangle1_I, iDTriangle2_I, iDTriangle3_I,&
         iDTriangle4_I, iDTriangle5_I, iDTriangle6_I,&
         iURectangle1_I, iURectangle2_I, iURectangle3_I,&
         iUTriangle1_I, iUTriangle2_I, iUTriangle3_I,&
         iUTriangle4_I, iUTriangle5_I, iUTriangle6_I)
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
           iUTriangle4_I, iUTriangle5_I, iUTriangle6_I
      integer, intent(in), optional, dimension(4)::&
           iDRectangle1_I,iDRectangle2_I,iDRectangle3_I,&
           iURectangle1_I,iURectangle2_I,iURectangle3_I
      !The point belonging to one of the up and down subfaces 
      real, dimension(nDim) :: XyzUp_D, XyzDown_D, X1_D, X2_D, X3_D, X4_D
      real:: AlphaUp, AlphaDown
      integer:: nGridOutUp, nGridOutDown
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
      ! the weghts for the upper face should be multiplied by 
      ! AlphaDown/(AlphaUp+AlphaDown)
      ! for the down face - by AlphaUp/(AlphaUp+AlphaDown)
      !/

      nGridOut = -1; nGridOutUp= -1; nGridOutDown = -1
      Weight_I = -1
      if(present(iUTriangle1_I))then
         X1_D = XyzGrid_DI(:,iUTriangle1_I(1))
         X2_D = XyzGrid_DI(:,iUTriangle1_I(2))
         X3_D = XyzGrid_DI(:,iUTriangle1_I(3))
         AlphaUp = triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
              triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
         if(AlphaUp==0.0)then
            call triangle(X1_D, X2_D, X3_D, Xyz_D)
            if(all(Weight_I(1:3)>=0.0))then
               iOrder_I(1:3) = iUTriangle1_I
               nGridOut = 3
               RETURN
            end if
         elseif(AlphaUp > 0.0)then
            XyzUp_D = Xyz_D + AlphaUp * Dir_D
            call triangle(X1_D, X2_D, X3_D, XyzUp_D)
            if(all(Weight_I(1:3)>=0.0))then
               iOrder_I(1:3) = iUTriangle1_I
               nGridOutUp = 3
               go to 700
            end if
         end if
         if(present(iUTriangle2_I))then
            X1_D = XyzGrid_DI(:,iUTriangle2_I(1))
            X2_D = XyzGrid_DI(:,iUTriangle2_I(2))
            X3_D = XyzGrid_DI(:,iUTriangle2_I(3))
            AlphaUp = triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
                 triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
            if(AlphaUp==0.0)then
               call triangle(X1_D, X2_D, X3_D, Xyz_D)
               if(all(Weight_I(1:3)>=0.0))then
                  iOrder_I(1:3) = iUTriangle2_I
                  nGridOut = 3
                  RETURN
               end if
            elseif(AlphaUp > 0.0)then
               XyzUp_D = Xyz_D + AlphaUp * Dir_D
               call triangle(X1_D, X2_D, X3_D, XyzUp_D)
               if(all(Weight_I(1:3)>=0.0))then
                  iOrder_I(1:3) = iUTriangle2_I
                  nGridOutUp = 3
                  go to 700
               end if
            end if
            if(present(iUTriangle3_I))then
               X1_D = XyzGrid_DI(:,iUTriangle3_I(1))
               X2_D = XyzGrid_DI(:,iUTriangle3_I(2))
               X3_D = XyzGrid_DI(:,iUTriangle3_I(3))
               AlphaUp = &
                    triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
                    triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
               if(AlphaUp==0.0)then
                  call triangle(X1_D, X2_D, X3_D, Xyz_D)
                  if(all(Weight_I(1:3)>=0.0))then
                     iOrder_I(1:3) = iUTriangle3_I
                     nGridOut = 3
                     RETURN
                  end if
               elseif(AlphaUp > 0.0)then
                  XyzUp_D = Xyz_D + AlphaUp * Dir_D
                  call triangle(X1_D, X2_D, X3_D, XyzUp_D)
                  if(all(Weight_I(1:3)>=0.0))then
                     iOrder_I(1:3) = iUTriangle3_I
                     nGridOutUp = 3
                     go to 700
                  end if
               end if
               if(present(iUTriangle4_I))then
                  X1_D = XyzGrid_DI(:,iUTriangle4_I(1))
                  X2_D = XyzGrid_DI(:,iUTriangle4_I(2))
                  X3_D = XyzGrid_DI(:,iUTriangle4_I(3))
                  AlphaUp = &
                       triple_product(X1_D - Xyz_D,X2_D - X1_D,X3_D - X1_D)/&
                       triple_product(Dir_D       ,X2_D - X1_D,X3_D - X1_D)
                  if(AlphaUp==0.0)then
                     call triangle(X1_D, X2_D, X3_D, Xyz_D)
                     if(all(Weight_I(1:3)>=0.0))then
                        iOrder_I(1:3) = iUTriangle4_I
                        nGridOut = 3
                        RETURN
                     end if
                  elseif(AlphaUp > 0.0)then
                     XyzUp_D = Xyz_D + AlphaUp * Dir_D
                     call triangle(X1_D, X2_D, X3_D, XyzUp_D)
                     if(all(Weight_I(1:3)>=0.0))then
                        iOrder_I(1:3) = iUTriangle4_I
                        nGridOutUp = 3
                        go to 700
                     end if
                  end if
                  if(present(iUTriangle5_I))then
                     X1_D = XyzGrid_DI(:,iUTriangle5_I(1))
                     X2_D = XyzGrid_DI(:,iUTriangle5_I(2))
                     X3_D = XyzGrid_DI(:,iUTriangle5_I(3))
                     AlphaUp = &
                          triple_product(X1_D - Xyz_D,X2_D - X1_D,X3_D - X1_D)&
                          /triple_product(Dir_D,X2_D - X1_D,X3_D - X1_D)
                     if(AlphaUp==0.0)then
                        call triangle(X1_D, X2_D, X3_D, Xyz_D)
                        if(all(Weight_I(1:3)>=0.0))then
                           iOrder_I(1:3) = iUTriangle5_I
                           nGridOut = 3
                           RETURN
                        end if
                     elseif(AlphaUp > 0.0)then
                        XyzUp_D = Xyz_D + AlphaUp * Dir_D
                        call triangle(X1_D, X2_D, X3_D, XyzUp_D)
                        if(all(Weight_I(1:3)>=0.0))then
                           iOrder_I(1:3) = iUTriangle5_I
                           nGridOutUp = 3
                           go to 700
                        end if
                     end if
                     if(present(iUTriangle6_I))then
                        X1_D = XyzGrid_DI(:,iUTriangle6_I(1))
                        X2_D = XyzGrid_DI(:,iUTriangle6_I(2))
                        X3_D = XyzGrid_DI(:,iUTriangle6_I(3))
                        AlphaUp = &
                             triple_product(&
                             X1_D - Xyz_D,X2_D - X1_D,X3_D - X1_D)/&
                             triple_product(&
                             Dir_D       ,X2_D - X1_D,X3_D - X1_D)
                        if(AlphaUp==0.0)then
                           call triangle(X1_D, X2_D, X3_D, Xyz_D)
                           if(all(Weight_I(1:3)>=0.0))then
                              iOrder_I(1:3) = iUTriangle6_I
                              nGridOut = 3
                              RETURN
                           end if
                        elseif(AlphaUp > 0.0)then
                           XyzUp_D = Xyz_D + AlphaUp * Dir_D
                           call triangle(X1_D, X2_D, X3_D, XyzUp_D)
                           if(all(Weight_I(1:3)>=0.0))then
                              iOrder_I(1:3) = iUTriangle6_I
                              nGridOutUp = 3
                              go to 700
                           end if
                        end if
                     end if !6
                  end if    !5
               end if       !4
            end if          !3
         end if             !2
      end if                !1

      if(present(iURectangle1_I))then
         X1_D = XyzGrid_DI(:,iURectangle1_I(1))
         X2_D = XyzGrid_DI(:,iURectangle1_I(2))
         X3_D = XyzGrid_DI(:,iURectangle1_I(3))
         X4_D = XyzGrid_DI(:,iURectangle1_I(4))
         AlphaUp = triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
              triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
         if(AlphaUp==0.0)then
            call rectangle(X1_D, X2_D, X3_D, X4_D, Xyz_D)
            if(all(Weight_I(1:4)>=0.0))then
               iOrder_I(1:4) = iURectangle1_I
               nGridOut = 4
               RETURN
            end if
         elseif(AlphaUp > 0.0)then
            XyzUp_D = Xyz_D + AlphaUp * Dir_D
            call rectangle(X1_D, X2_D, X3_D, X4_D, XyzUp_D)
            if(all(Weight_I(1:4)>=0.0))then
               iOrder_I(1:4) = iURectangle1_I
               nGridOutUp = 4
               go to 700
            end if
         end if
         if(present(iURectangle2_I))then
            X1_D = XyzGrid_DI(:,iURectangle2_I(1))
            X2_D = XyzGrid_DI(:,iURectangle2_I(2))
            X3_D = XyzGrid_DI(:,iURectangle2_I(3))
            X4_D = XyzGrid_DI(:,iURectangle2_I(4))
            AlphaUp = triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
                 triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
            if(AlphaUp==0.0)then
               call rectangle(X1_D, X2_D, X3_D, X4_D, Xyz_D)
               if(all(Weight_I(1:4)>=0.0))then
                  iOrder_I(1:4) = iURectangle2_I
                  nGridOut = 4
                  RETURN
               end if
            elseif(AlphaUp > 0.0)then
               XyzUp_D = Xyz_D + AlphaUp * Dir_D
               call rectangle(X1_D, X2_D, X3_D, X4_D, XyzUp_D)
               if(all(Weight_I(1:4)>=0.0))then
                  iOrder_I(1:4) = iURectangle2_I
                  nGridOutUp = 4
                  go to 700
               end if
            end if
            if(present(iURectangle3_I))then
               X1_D = XyzGrid_DI(:,iURectangle3_I(1))
               X2_D = XyzGrid_DI(:,iURectangle3_I(2))
               X3_D = XyzGrid_DI(:,iURectangle3_I(3))
               X4_D = XyzGrid_DI(:,iURectangle3_I(4))
               AlphaUp = &
                    triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
                    triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
               if(AlphaUp==0.0)then
                  call rectangle(X1_D, X2_D, X3_D, X4_D, Xyz_D)
                  if(all(Weight_I(1:4)>=0.0))then
                     iOrder_I(1:4) = iURectangle3_I
                     nGridOut = 4
                     RETURN
                  end if
               elseif(AlphaUp > 0.0)then
                  XyzUp_D = Xyz_D + AlphaUp * Dir_D
                  call rectangle(X1_D, X2_D, X3_D, X4_D, XyzUp_D)
                  if(all(Weight_I(1:4)>=0.0))then
                     iOrder_I(1:4) = iURectangle3_I
                     nGridOutUp = 4
                     go to 700
                  end if
               end if
            end if     !3
         end if        !2
      end if           !1
      !No intersection point with upper boundary is found
      Weight_I = -1; nGridOut = -1
      RETURN
700   continue
      !\
      !Calculate low face
      if(present(iDTriangle1_I))then
         Weight_I(4:3+nGridOutUp) = Weight_I(1:nGridOutUp)
         iOrder_I(4:3+nGridOutUp) = iOrder_I(1:nGridOutUp)
         X1_D = XyzGrid_DI(:,iDTriangle1_I(1))
         X2_D = XyzGrid_DI(:,iDTriangle1_I(2))
         X3_D = XyzGrid_DI(:,iDTriangle1_I(3))
         AlphaDown = &
              triple_product(Xyz_D - X1_D, X2_D - X1_D, X3_D - X1_D)/&
              triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
         if(AlphaDown==0.0)then
            call triangle(X1_D, X2_D, X3_D, Xyz_D)
            if(all(Weight_I(1:3)>=0.0))then
               iOrder_I(1:3) = iDTriangle1_I
               nGridOut = 3
               RETURN
            end if
         elseif(AlphaDown > 0.0)then
            XyzDown_D = Xyz_D - AlphaDown * Dir_D
            call triangle(X1_D, X2_D, X3_D, XyzDown_D)
            if(all(Weight_I(1:3)>=0.0))then
               iOrder_I(1:3) = iDTriangle1_I
               nGridOutDown = 3
               go to 800
            end if
         end if
         if(present(iDTriangle2_I))then
            X1_D = XyzGrid_DI(:,iDTriangle2_I(1))
            X2_D = XyzGrid_DI(:,iDTriangle2_I(2))
            X3_D = XyzGrid_DI(:,iDTriangle2_I(3))
            AlphaDown = &
                 triple_product(Xyz_D - X1_D, X2_D - X1_D, X3_D - X1_D)/&
                 triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
            if(AlphaDown==0.0)then
               call triangle(X1_D, X2_D, X3_D, Xyz_D)
               if(all(Weight_I(1:3)>=0.0))then
                  iOrder_I(1:3) = iDTriangle2_I
                  nGridOut = 3
                  RETURN
               end if
            elseif(AlphaDown > 0.0)then
               XyzDown_D = Xyz_D - AlphaDown * Dir_D
               call triangle(X1_D, X2_D, X3_D, XyzDown_D)
               if(all(Weight_I(1:3)>=0.0))then
                  iOrder_I(1:3) = iDTriangle2_I
                  nGridOutDown = 3
                  go to 800
               end if
            end if
            if(present(iDTriangle3_I))then
               X1_D = XyzGrid_DI(:,iDTriangle3_I(1))
               X2_D = XyzGrid_DI(:,iDTriangle3_I(2))
               X3_D = XyzGrid_DI(:,iDTriangle3_I(3))
               AlphaDown = &
                    triple_product(Xyz_D - X1_D,X2_D - X1_D,X3_D - X1_D)/&
                    triple_product(Dir_D       ,X2_D - X1_D,X3_D - X1_D)
               if(AlphaDown==0.0)then
                  call triangle(X1_D, X2_D, X3_D, Xyz_D)
                  if(all(Weight_I(1:3)>=0.0))then
                     iOrder_I(1:3) = iDTriangle3_I
                     nGridOut = 3
                     RETURN
                  end if
               elseif(AlphaDown > 0.0)then
                  XyzDown_D = Xyz_D - AlphaDown * Dir_D
                  call triangle(X1_D, X2_D, X3_D, XyzDown_D)
                  if(all(Weight_I(1:3)>=0.0))then
                     iOrder_I(1:3) = iDTriangle3_I
                     nGridOutDown = 3
                     go to 800
                  end if
               end if
               if(present(iDTriangle4_I))then
                  X1_D = XyzGrid_DI(:,iDTriangle4_I(1))
                  X2_D = XyzGrid_DI(:,iDTriangle4_I(2))
                  X3_D = XyzGrid_DI(:,iDTriangle4_I(3))
                  AlphaDown = &
                       triple_product(Xyz_D - X1_D,X2_D - X1_D,X3_D - X1_D)/&
                       triple_product(Dir_D       ,X2_D - X1_D,X3_D - X1_D)
                  if(AlphaDown==0.0)then
                     call triangle(X1_D, X2_D, X3_D, Xyz_D)
                     if(all(Weight_I(1:3)>=0.0))then
                        iOrder_I(1:3) = iDTriangle3_I
                        nGridOut = 3
                        RETURN
                     end if
                  elseif(AlphaDown > 0.0)then
                     XyzDown_D = Xyz_D - AlphaDown * Dir_D
                     call triangle(X1_D, X2_D, X3_D, XyzDown_D)
                     if(all(Weight_I(1:3)>=0.0))then
                        iOrder_I(1:3) = iDTriangle4_I
                        nGridOutDown = 3
                        go to 800
                     end if
                  end if
                  if(present(iDTriangle5_I))then
                     X1_D = XyzGrid_DI(:,iDTriangle5_I(1))
                     X2_D = XyzGrid_DI(:,iDTriangle5_I(2))
                     X3_D = XyzGrid_DI(:,iDTriangle5_I(3))
                     AlphaDown = &
                          triple_product(&
                          Xyz_D - X1_D,X2_D - X1_D,X3_D - X1_D)/&
                          triple_product(&
                          Dir_D       ,X2_D - X1_D,X3_D - X1_D)
                     if(AlphaDown==0.0)then
                        call triangle(X1_D, X2_D, X3_D, Xyz_D)
                        if(all(Weight_I(1:3)>=0.0))then
                           iOrder_I(1:3) = iDTriangle5_I
                           nGridOut = 3
                           RETURN
                        end if
                     elseif(AlphaDown > 0.0)then
                        XyzDown_D = Xyz_D - AlphaDown * Dir_D
                        call triangle(X1_D, X2_D, X3_D, XyzDown_D)
                        if(all(Weight_I(1:3)>=0.0))then
                           iOrder_I(1:3) = iDTriangle5_I
                           nGridOutDown = 3
                           go to 800
                        end if
                     end if
                     if(present(iDTriangle6_I))then
                        X1_D = XyzGrid_DI(:,iDTriangle6_I(1))
                        X2_D = XyzGrid_DI(:,iDTriangle6_I(2))
                        X3_D = XyzGrid_DI(:,iDTriangle6_I(3))
                        AlphaDown = &
                             triple_product(&
                             Xyz_D - X1_D,X2_D - X1_D,X3_D - X1_D)/&
                             triple_product(&
                             Dir_D       ,X2_D - X1_D,X3_D - X1_D)
                        if(AlphaDown==0.0)then
                           call triangle(X1_D, X2_D, X3_D, Xyz_D)
                           if(all(Weight_I(1:3)>=0.0))then
                              iOrder_I(1:3) = iDTriangle6_I
                              nGridOut = 3
                              RETURN
                           end if
                        elseif(AlphaDown > 0.0)then
                           XyzDown_D = Xyz_D - AlphaDown * Dir_D
                           call triangle(X1_D, X2_D, X3_D, XyzDown_D)
                           if(all(Weight_I(1:3)>=0.0))then
                              iOrder_I(1:3) = iDTriangle6_I
                              nGridOutDown = 3
                              go to 800
                           end if
                        end if
                     end if !6
                  end if    !5
               end if       !4
            end if          !3
         end if             !2
      end if                !1

      if(present(iDRectangle1_I))then
         if(present(iDTriangle1_I))then
            Weight_I(5:4+nGridOutUp) = Weight_I(4:3+nGridOutUp)
            iOrder_I(5:4+nGridOutUp) = iOrder_I(4:3+nGridOutUp)
         else
            Weight_I(5:4+nGridOutUp) = Weight_I(1:nGridOutUp)
            iOrder_I(5:4+nGridOutUp) = iOrder_I(1:nGridOutUp)
         end if
         X1_D = XyzGrid_DI(:,iDRectangle1_I(1))
         X2_D = XyzGrid_DI(:,iDRectangle1_I(2))
         X3_D = XyzGrid_DI(:,iDRectangle1_I(3))
         X4_D = XyzGrid_DI(:,iDRectangle1_I(4))
         AlphaDown = triple_product(Xyz_D - X1_D, X2_D - X1_D, X3_D - X1_D)/&
              triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
         if(AlphaDown==0.0)then
            call rectangle(X1_D, X2_D, X3_D, X4_D, Xyz_D)
            if(all(Weight_I(1:4)>=0.0))then
               iOrder_I(1:4) = iDRectangle1_I
               nGridOut = 4
               RETURN
            end if
         elseif(AlphaDown > 0.0)then
            XyzDown_D = Xyz_D - AlphaDown * Dir_D
            call rectangle(X1_D, X2_D, X3_D, X4_D, XyzDown_D)
            if(all(Weight_I(1:4)>=0.0))then
               iOrder_I(1:4) = iDRectangle1_I
               nGridOutDown = 4
               go to 800
            end if
         end if
         if(present(iDRectangle2_I))then

            X1_D = XyzGrid_DI(:,iDRectangle2_I(1))
            X2_D = XyzGrid_DI(:,iDRectangle2_I(2))
            X3_D = XyzGrid_DI(:,iDRectangle2_I(3))
            X4_D = XyzGrid_DI(:,iDRectangle2_I(4))
            AlphaDown = &
                 triple_product(Xyz_D - X1_D, X2_D - X1_D, X3_D - X1_D)/&
                 triple_product(Dir_D, X2_D - X1_D, X3_D - X1_D)
            if(AlphaDown==0.0)then
               call rectangle(X1_D, X2_D, X3_D, X4_D, Xyz_D)
               if(all(Weight_I(1:4)>=0.0))then
                  iOrder_I(1:4) = iDRectangle2_I
                  nGridOut = 4
                  RETURN
               end if
            elseif(AlphaDown > 0.0)then
               XyzDown_D = Xyz_D - AlphaDown * Dir_D
               call rectangle(X1_D, X2_D, X3_D, X4_D, XyzDown_D)
               if(all(Weight_I(1:4)>=0.0))then
                  iOrder_I(1:4) = iDRectangle2_I
                  nGridOutDown = 4
                  go to 800
               end if
            end if
            if(present(iDRectangle3_I))then

               X1_D = XyzGrid_DI(:,iDRectangle3_I(1))
               X2_D = XyzGrid_DI(:,iDRectangle3_I(2))
               X3_D = XyzGrid_DI(:,iDRectangle3_I(3))
               X4_D = XyzGrid_DI(:,iDRectangle3_I(4))
               AlphaDown = &
                    triple_product(Xyz_D - X1_D, X2_D - X1_D, X3_D - X1_D)/&
                    triple_product(Dir_D, X2_D - X1_D, X3_D - X1_D)
               if(AlphaDown==0.0)then
                  call rectangle(X1_D, X2_D, X3_D, X4_D, Xyz_D)
                  if(all(Weight_I(1:4)>=0.0))then
                     iOrder_I(1:4) = iDRectangle3_I
                     nGridOut = 4
                     RETURN
                  end if
               elseif(AlphaDown > 0.0)then
                  XyzDown_D = Xyz_D - AlphaDown * Dir_D
                  call rectangle(X1_D, X2_D, X3_D, X4_D, XyzDown_D)
                  if(all(Weight_I(1:4)>=0.0))then
                     iOrder_I(1:4) = iDRectangle3_I
                     nGridOutDown = 4
                     go to 800
                  end if
               end if
            end if     !3
         end if        !2
      end if           !1
      nGridOut = -1; Weight_I = -1  !No intersection with the down subface
      RETURN
800   continue
      !Apply the weight for interpolation along the rate
      nGridOut = nGridOutUp + nGridOutDown
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
      real:: x,y
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

    !=======================================================================
    subroutine corner_transition_junction(iEdgeDir)
      integer, intent(in) :: iEdgeDir
      integer             :: iAxis,jAxis,kAxis
      real :: AlphaKAxis
      !\
      !! Auxilary variables to determine vertical weight
      !/
      real :: dXyzUp, dXyzDown 
      !\
      ! Subroutine intepolates in junction of transitional regions near a &
      ! resolution corner
      ! Face with all Fine points is numbered 5:8
      ! Face with three Coarse points is numbered 1:4
      !
      !  F---F
      !  |\ /|
      !  | C-|-----C
      !  |/|\|
      !  F---F
      !    |     F
      !    |
      !    C
      !
      !---------------------------------------------------------------------
      kAxis = iEdgeDir
      iAxis = 1 + mod(iEdgeDir    ,3)
      jAxis = 1 + mod(iEdgeDir + 1,3)
      Xyz2_D(x_:y_) = Xyz_D((/iAxis,jAxis/))
    
      !\
      ! may need to rearrange iOrder
      !/
      if(XyzGrid_DI(iAxis,iOrder_I(5))>XyzGrid_DI(iAxis,iOrder_I(6)))then
         iOrder_I(1:8) = iOrder_I((/2,1,4,3,6,5,8,7/))
      end if
      if(XyzGrid_DI(jAxis,iOrder_I(5))>XyzGrid_DI(jAxis,iOrder_I(7)))then
         iOrder_I(1:8) = iOrder_I((/3,4,1,2,7,8,5,6/))
      end if

      !\
      ! Interpolate along face 5:8
      !/
      DoStencilFix2 = .false.
      XyzGrid2_DI(x_:y_,1:4) = XyzGrid_DI((/iAxis,jAxis/),iOrder_I(5:8))
      iLevel2_I(        1:4) = iLevel_I(                  iOrder_I(5:8))
      call interpolate_amr_grid2(&
           Xyz2_D, XyzGrid2_DI, iLevel2_I, &
           nGridOut2, Weight_I(5:8), iOrder2_I, DoStencilFix2, XyzStencil2_D)
      !\
      ! Rearrange points according to 2D interpolation output
      !/
      iOrder_I(1:4) = iOrder_I(iOrder2_I  )
      iOrder_I(5:8) = iOrder_I(iOrder2_I+4)
      nGridOut = nGridOut2

      !\
      ! find Fine point on 1:4 and substitute it 
      ! with corresponding point from face 5:8
      !/
      do iGrid = 1,4
         if(iLevel_I(iOrder_I(iGrid))==Fine_)then
            iOrder_I(iGrid) = iOrder_I(iGrid + 4)
         end if
      end do


      !\
      ! Interpolate along face 1:4
      !/
      DoStencilFix2 = .false.
      XyzGrid2_DI(x_:y_,1:4) = XyzGrid_DI((/iAxis,jAxis/),iOrder_I(1:4))
      iLevel2_I(        1:4) = iLevel_I(                  iOrder_I(1:4))
      call interpolate_amr_grid2(&
           Xyz2_D, XyzGrid2_DI, iLevel2_I, &
           nGridOut2, Weight_I(1:4), iOrder2_I, DoStencilFix2, XyzStencil2_D)
      !\
      ! Rearrange points according to 2D interpolation output
      !/
      iOrder_I(1:4) = iOrder_I(iOrder2_I  )
      iOrder_I(5:8) = iOrder_I(iOrder2_I+4)
      Weight_I(5:8) = Weight_I(iOrder2_I+4)

      !\
      ! move Fine point on 1:4 to position 1 and 
      ! corresponding point to position 5
      !/
      if(iLevel_I(iOrder_I(2))==Fine_)then
         iOrder_I((/1,2,5,6/)) = iOrder_I((/2,1,6,5/))
         Weight_I((/1,2,5,6/)) = Weight_I((/2,1,6,5/))
      end if

      !\
      ! find vertical weight
      !/
      dXyzUp   = Xyz_D(kAxis) - XyzGrid_DI(kAxis,iOrder_I(5))
      dXyzDown = Xyz_D(kAxis) - &
           sum(XyzGrid_DI(kAxis,iOrder_I(1:nGridOut2))*Weight_I(1:nGridOut2))

      !vertical_distance_point_plane(Xyz_D,(/1,2,3/),kAxis)

      DoStencilFix = dXyzUp * dXyzDown > 0.0
      if(DoStencilFix)then
         XyzStencil_D((/iAxis,jAxis/)) =               &
              XyzGrid_DI((/iAxis,jAxis/),iOrder_I(2)) + &
              2*(XyzGrid_DI((/iAxis,jAxis/),iOrder_I(1)) - &
              XyzGrid_DI((/iAxis,jAxis/),iOrder_I(2)))
         XyzStencil_D(kAxis) = Xyz_D(kAxis)
         RETURN
      end if

      if(dXyzUp == 0.0)then
         AlphaKAxis = 0.0
      else
         AlphaKAxis = abs(dXyzUp) / (abs(dXyzUp) + abs(dXyzDown))
      end if


      !\
      ! Apply vertical weight
      !/
      Weight_I(1           ) = Weight_I(5           ) * (1 - AlphaKAxis) + &
           Weight_I(1           ) *      AlphaKAxis
      Weight_I(2:nGridOut2 ) = Weight_I(2:nGridOut2 ) *      AlphaKAxis
      Weight_I(6:4+nGridOut) = Weight_I(6:4+nGridOut) * (1 - AlphaKAxis)
      !\
      ! Remove points
      !/
      iOrder_I(nGridOut2+1:nGridOut2+nGridOut-1) = iOrder_I(6:4+nGridOut)
      Weight_I(nGridOut2+1:nGridOut2+nGridOut-1) = Weight_I(6:4+nGridOut)
      nGridOut = nGridOut + nGridOut2 - 1

    end subroutine corner_transition_junction
    !=======================================================================
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
      logical :: IsAmbiguous
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
      ! may need to rearrange iOrder
      !/
      if(XyzGrid_DI(iAxis,iOrder_I(5))>XyzGrid_DI(iAxis,iOrder_I(6)))then
         iOrder_I(1:8) = iOrder_I((/2,1,4,3,6,5,8,7/))
      end if
      if(XyzGrid_DI(jAxis,iOrder_I(5))>XyzGrid_DI(jAxis,iOrder_I(7)))then
         iOrder_I(1:8) = iOrder_I((/3,4,1,2,7,8,5,6/))
      end if


      !\
      ! Interpolate along face 1:4
      !/
      DoStencilFix2 = DoStencilFix
      XyzGrid2_DI(x_:y_,1:4) = XyzGrid_DI((/iAxis,jAxis/),iOrder_I(1:4))
      iLevel2_I(        1:4) = iLevel_I(                  iOrder_I(1:4))
      call interpolate_amr_grid2(&
           Xyz2_D, XyzGrid2_DI, iLevel2_I, &
           nGridOut2, Weight_I(1:4), iOrder2_I, DoStencilFix2, XyzStencil2_D)

      !\
      ! Check if 2D interpolation requires a fix
      !
      !   F---F
      !    \ / \
      !     C---F
      !/
      if(DoStencilFix2)then
         DoStencilFix = DoStencilFix2
         !\
         ! Move to another transition region
         !/
         XyzStencil_D((/iAxis,jAxis/)) = XyzStencil2_D
         XyzStencil_D(kAxis) = Xyz_D(kAxis)
         RETURN
      end if

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
      if(count(iLevel_I(iOrder_I(5:8))==Coarse_)==1   .and.&
           count(iLevel_I(iOrder_I(1:4))==Fine_  )==1   .and.&
           iLevel_I(iOrder_I(5))  ==Fine_         .and.&
           iLevel_I(iOrder_I(6))  ==Fine_)then
         !\
         ! may be the wrong face
         !             C
         !      F     /|\
         !     /|\   / | \
         !    / | \ /F---F
         !   C-----C \ | /
         !    \ | /   \|/
         !     \|/     C
         !      F
         !/
         IsAmbiguous = .not.(DoStencilFix)
      else
         IsAmbiguous = .false.
      end if

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
                 (/iOrder_I(2),iOrder_I(3),iOrder_I(6),iOrder_I(7)/),&
                 iUTriangle1_I = &
                 (/iOrder_I(1),iOrder_I(2), iOrder_I(5)/))
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
      if(DoStencilFix)then
         XyzStencil_D(iAxis) = 0.25*sum(XyzGrid_DI(iAxis,iOrder_I(1:4)))
         XyzStencil_D(jAxis) = 0.25*sum(XyzGrid_DI(jAxis,iOrder_I(1:4)))
         if(IsAmbiguous)then
            XyzStencil_D((/iAxis,jAxis/)) = Xyz_D((/iAxis, jAxis/))
         end if
         XyzStencil_D(kAxis) = XyzGrid_DI(kAxis,iOrder_I(3)) + &
              XyzGrid_DI(kAxis,iOrder_I(5))-XyzGrid_DI(kAxis,iOrder_I(1))
         RETURN
      end if

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
  end subroutine interpolate_amr_grid3
  !============================================================================
  !Interpolation on the block AMR grid
  !\
  !Calculates interpolation weights
  !/
  !\
  !Example of application for SERIAL calculation of the
  !interpolated value of the state vector sampled in the 
  !grid points as
  !State_VGB(nVar,nI,nJ,nK,nBlock) array:
  !
  !call interpolate_amr(&
  ! 3, Xyz_D, 4, (/4,4,4/) ...., nGridOut, Weight_I,iIndexes_II)
  ! if(nGridOut < 1) call CON_stop('Interpolation failed')
  ! Value_V(1:nVar) = 0
  ! do iGrid = 1, nGridOut
  !    Value_V = Value_V + &
  !     State_VGB(:,iIndexes_II(1,iGrid), iIndexes_II(2,iGrid), &
  !             iIndexes_II(3,iGrid), iIndexes_II(4,iGrid))*&
  !                               Weight_I(iGrid)
  !end do 
  !/
  subroutine interpolate_block_amr(nDim, XyzIn_D, nIndexes, nCell_D, find, &
       nGridOut, Weight_I, iIndexes_II, iLevelOut_I)
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
    integer, intent(out), optional:: iLevelOut_I(2**nDim)  !For diagnostics

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
    ! The output array in interpolate_amr_grid2,3  routine
    ! This is the set of numbers of the elements of the basic stencil
    ! to be involved in interpolation
    !/
    integer, dimension(2**nDim):: iOrder_I
    !\
    ! Output parameter of interpolate_amr_grid2,3
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
    real, dimension(nDim,2**(2*nDim)) :: XyzExtended_DI 
    integer :: iLevelExtended_I(2**(2*nDim))
    integer :: iIndexesExtended_II(0:nIndexes,2**(2*nDim))
    integer :: nExtendedStencil, nIter

    !\
    ! The same extended stencil in a structured form:
    ! A cubic 2*2*2 grid with 2*2*2 subgrids covering each vertex
    !/
    real       :: XyzGrid_DII(nDim,0:2**nDim,2**nDim)
    integer, dimension(2**nDim):: iBlock_I  , iProc_I 
    integer                    :: iCellIndexes_DII(nDim,2**nDim,2**nDim)
    integer, dimension(2**nDim):: nSubGrid_I, iLevelSubgrid_I
    !--------------------
    !\
    ! Initialize output arrays
    !/

    iOrder_I    = 0
    iIndexes_II = 0
    Weight_I    = 0 
    if(present(iLevelOut_I))iLevelOut_I=0

    Xyz_D = XyzIn_D
    call generate_extended_stencil
    if(nExtendedStencil<1)call CON_stop('No extended stencil!')

    if(nExtendedStencil < 2**nDim)then
       !\
       ! Failure in constructing the extended stencil
       ! May occur if the input point is out of the
       ! AMR grid.
       !/
       nGridOut = -1
       RETURN
    end if
    call generate_basic_stencil(&
         nDim, Xyz_D, nExtendedStencil,        &
         XyzExtended_DI(:,1:nExtendedStencil), &
         DxyzInv_D, iOrderExtended_I)
    if(any(iOrderExtended_I < 1))&
         call CON_stop('Failure in constructing basic stencil') 


    XyzGrid_DI = XyzExtended_DI(:,iOrderExtended_I)
    iLevel_I = iLevelExtended_I(iOrderExtended_I)
    iIndexes_II = &
         iIndexesExtended_II(:,iOrderExtended_I)
    DoStencilFix = .false. 
    select case(nDim)
    case(2)
       call interpolate_amr_grid2(&
            Xyz_D , XyzGrid_DI, iLevel_I, &
            nGridOut, Weight_I, iOrder_I,&
            DoStencilFix, XyzStencil_D)
       nIter=0
       do while(DoStencilFix)
          nIter=nIter+1
          if(nIter==2)call CON_stop('Infinite loop in fixing stencil')
          call generate_basic_stencil(&
               nDim, XyzStencil_D, nExtendedStencil, &
               XyzExtended_DI(:,1:nExtendedStencil), &
               DxyzInv_D, iOrderExtended_I)
          if(any(iOrderExtended_I < 1))&
               call CON_stop('Failure in constructing basic stencil') 
          XyzGrid_DI = XyzExtended_DI(&
               :,iOrderExtended_I)
          iLevel_I = iLevelExtended_I(iOrderExtended_I)
          iIndexes_II = &
               iIndexesExtended_II(:,iOrderExtended_I)
          call interpolate_amr_grid2(&
               Xyz_D , XyzGrid_DI, iLevel_I, &
               nGridOut, Weight_I, iOrder_I,&
               DoStencilFix, XyzStencil_D)
       end do
       if(nGridOut < 1)&
            call CON_stop('Failure in interpolate amr grid2') 
    case(3)
       call interpolate_amr_grid3(&
            Xyz_D , XyzGrid_DI, iLevel_I, &
            nGridOut, Weight_I, iOrder_I,&
            DoStencilFix, XyzStencil_D)
       nIter = 0
       do while(DoStencilFix)
          nIter=nIter+1
          if(nIter==3)call CON_stop('Infinite loop in fixing stencil')
          call generate_basic_stencil(&
               nDim, XyzStencil_D, nExtendedStencil,        &
               XyzExtended_DI(:,1:nExtendedStencil), &
               DxyzInv_D, iOrderExtended_I)
          if(any(iOrderExtended_I < 1))&
               call CON_stop('Failure in constructing basic stencil') 
          XyzGrid_DI = XyzExtended_DI(&
               :,iOrderExtended_I)
          iLevel_I = iLevelExtended_I(iOrderExtended_I)
          iIndexes_II = &
               iIndexesExtended_II(:,iOrderExtended_I)
          call interpolate_amr_grid3(&
               Xyz_D , XyzGrid_DI, iLevel_I, &
               nGridOut, Weight_I, iOrder_I,&
               DoStencilFix, XyzStencil_D)
       end do
       if(nGridOut < 1)&
            call CON_stop('Failure in interpolate amr grid3') 
    case default
       call CON_stop('Only 2D and 3D AMR grids are implemented')
    end select
    iIndexes_II(:, 1:nGridOut) = iIndexes_II(:,iOrder_I(1:nGridOut))
    if(present(iLevelOut_I))iLevelOut_I = iLevel_I
  contains
    subroutine generate_extended_stencil
      use ModKind, ONLY: nByteReal
      !\
      !This routine extracts the grid points which in principle 
      !can be involved into interpolation
      !/
      !\
      !Just 2**nDim
      integer:: nGrid

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
      ! To inprove the algorithm stability against roundoff errors
      !/
      real, parameter:: cTol = 0.00000010
      !------------------------

      nGrid = 2**nDim !Number of points in a basic stencil

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
         ! The algorithm cannot hanldle point out of the computation
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
         XyzGrid_DII(:,0,1) = Dxyz_D*(iCellIndexes_DII(:,1,1) - 0.50)

         XyzMisc_D = XyzMisc_D - iCellIndexes_DII(:,1,1)
         !\ 
         ! now XyzMisc_D = (Xyz_D-XyzGrid_DII(:,0,1))/Dxyz satisfies 
         ! inequalities: XyzMisc_D >= 0 and XyzMisc_D < 1. Strengthen 
         ! these inequalities making them to have a gap of cTol, if the
         ! code is compiled with a single precision, the gap of 
         ! cTol**2 otherwise:
         !/
         XyzMisc_D = min(1 - cTol**(nByteReal/4),&
              max(XyzMisc_D, cTol**(nByteReal/4) ))
         Xyz_D = XyzGrid_DII(:,0,1) + XyzMisc_D*Dxyz_D
         !\
         !Calculate other grid points, check if all points belong to 
         !the found block
         !/
         iGridCheck = -1
         do iGrid = nGrid,1,-1
            iShift_D = iShift_DI(1:nDim,iGrid)
            iBlock_I(iGrid) = iBlock_I(1)
            iCellIndexes_DII(:,1,iGrid) = &
                 iCellIndexes_DII(:,1,1) + iShift_D
            XyzGrid_DII(:,0,iGrid) = &
                 XyzGrid_DII(:,0,1) + iShift_D*Dxyz_D
            !\
            ! Check if the grid point is out of the block
            !/
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
               iProc_I(iGridStored) = 0 !For not processing this point again
               iLevelSubGrid_I(iGridStored) = BehindTheBoundary_
               XyzGrid_DII(:,1,iGridStored) = XyzGrid_DII(:,0,iGridStored)
               nSubGrid_I(iGridStored)      =1
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
            iLevelSubgrid_I(iGridStored) = &
                 -int(2*(Dxyz_D(1)*DXyzInv_D(1)) - 3 + cTol)
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
               do iGrid = iGridStored -1, 1, -1
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
               !
               ! Fine subgrid is displaced by half of finer subgrid size
               !
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
                  !\
                  ! Check if the point is in the newly found block
                  !/ 
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
      select case(count(iLevelSubgrid_I==BehindTheBoundary_))
      case(0, 3, 7) !3,7 are nGrid-1 for nDim=2,3 
         !\
         !Do nothing: if there is a single physical grid point 
         !then all coordinates of the ghost points are set earlier
         !/
      case(4)
         ! nDim = 3
         !\
         ! One of the faces is fully out of the domain
         !/

         !\
         !Find point in the domain
         !/
         iGridCheck = maxloc(iLevelSubGrid_I,DIM=1)
         do iDir = 1, nDim
            if(all(&
                 iLevelSubGrid_I(iOppositeFace_IDI(:,iDir,iGridCheck))&
                 ==BehindTheBoundary_))then
               !\
               ! Prolong grid from the physical phace to the ghost face
               !/
               do iGrid = 1,4
                  call prolong(iFace_IDI(iGrid,iDir,iGridCheck), &
                       iOppositeFace_IDI(iGrid,iDir,iGridCheck))
               end do
               EXIT
            end if
         end do
      case(2)
         ! Analogous to the previos case, but nDim = 2 

         !\
         !Check if the points behind the boundary are 1,2 or 1,3
         !/
         if(iLevelSubGrid_I(1)==BehindTheBoundary_)then
            if(iLevelSubGrid_I(2)==BehindTheBoundary_)then
               !\
               !Prolong grid from 3,4 to 1,2
               !/
               do iGrid = 1,2
                  call prolong(iGrid + 2, iGrid)
               end do
            else
               !\
               !Prolong grid from 2,4 to 1,3
               !/
               do iGrid = 1, 3, 2
                  call prolong(iGrid + 1, iGrid)
               end do
            end if
         else
            !\
            !Check if the points behind the boundary are 3,4 or 2,4
            !/
            if(iLevelSubGrid_I(3)==BehindTheBoundary_)then
               !\
               !Prolong grid from 1,2 to 3,4
               !/
               do iGrid = 1,2
                  call prolong(iGrid, iGrid +2)
               end do
            else
               !\
               !Prolong grid from 1,3 to 2,4
               !/
               do iGrid = 1, 3, 2
                  call prolong(iGrid, iGrid + 1)
               end do
            end if
         end if
      case(6)
         !\
         ! Three-dimensional case, only two grid vertexes are
         ! inside the domain. From each of the two physical
         ! points the grid should be prolonged to three
         ! ghost points along the plane
         !/
         iGridCheck = maxloc(iLevelSubGrid_I,DIM=1)
         do iDir  = 1, nDim
            if(iLevelSubGrid_I(iEdge_ID(iGridCheck,iDir))>=0)then
               do iGrid = 2,4
                  call prolong(iGridCheck,&
                       iFace_IDI(iGrid,iDir,iGridCheck))
                  call prolong(iEdge_ID(iGridCheck,iDir),&
                       iOppositeFace_IDI(iGrid,iDir,iGridCheck))
               end do
               EXIT
            end if
         end do
      end select

      !\
      ! Finalize: calculate all arrays for extended stencil
      ! in the upper routine:
      !/
      nExtendedStencil    = 0
      XyzExtended_DI      = 0
      iIndexesExtended_II = 0
      iLevelExtended_I    = 0
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
         end do
      end do
    end subroutine generate_extended_stencil
    !\
    ! Used to prolong grid behind the boundary
    !/
    subroutine prolong(iGridPhys, iGridGhost)          
      integer, intent(in) :: iGridPhys, iGridGhost
      !Loop variable
      integer :: iSubGrid
      !--------------------
      nSubGrid_I(iGridGhost) = nSubGrid_I(iGridPhys)
      do iSubGrid = 1, nSubGrid_I(iGridGhost)
         XyzGrid_DII(:,iSubGrid,iGridGhost) = &
              XyzGrid_DII(:,iSubGrid,iGridPhys) +&
              XyzGrid_DII(:,0,iGridGhost)       -&
              XyzGrid_DII(:,0,iGridPhys )
      end do
    end subroutine prolong
  end subroutine interpolate_block_amr
  !=======================================================================
  subroutine generate_basic_stencil(&
       nDim, Xyz_D, nExtendedStencil, XyzExtended_DI, dXyzInv_D, &
       iOrderExtended_I)
    !\
    ! Dimensionality; number of points in the extended stencil
    !/
    integer, intent(in) :: nDim, nExtendedStencil
    !\
    ! Point about which to construct basic stencil; inverse of (/Dx,Dy,Dz/)
    !/
    real,    intent(in) :: Xyz_D(nDim), dXyzInv_D(nDim)
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
    !-------------------------------------
    nGrid = 2**nDim
    iOrderExtended_I = -1
    if(nExtendedStencil==nGrid)then
       !\
       ! The basic stencil is identical to the extended one
       ! Check the stencil before use:
       !/
       do iGrid = 1, nGrid
          if(all(&
               Xyz_D < XyzExtended_DI(:,iGrid).eqv.IsBelow_DI(1:nDim,iGrid)&
               ))then
             iOrderExtended_I(iGrid) = iGrid
          end if
       end do
    else
       do iPoint = 1,nExtendedStencil
          IsBelowExtended_DI(:,iPoint) = Xyz_D < XyzExtended_DI(:,iPoint)
          Distance_I(iPoint) = sum( &
               ((Xyz_D - XyzExtended_DI(:,iPoint))*dXyzInv_D)**2)
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
    end if
  end subroutine generate_basic_stencil
  !===========================TESTS============================================
  subroutine test_interpolate_amr(nDim,nSample)
    integer, intent(in)::nDim, nSample
  
    integer :: iIndexes_II(0:nDim+1,2**nDim), iLevelOut_I(2**nDim)
    real, dimension(nDim):: DxyzDomain_D, DxyzCoarseBlock_D, &
         DxyzFineBlock_D, DxyzCoarse_D, &
         DxyzFine_D, Xyz_D,&
         XyzInterpolated_D, XyzCorner_D
    real:: Weight_I(2**nDim)
    !Loop variables
    integer :: iCase, iSample, iGrid, iSubGrid, i, j, k, iBlock, iDir
    
    integer :: nCell_D(3)  ! Cells per block
    integer :: iCellIndex_D(3)
    integer :: nIndexes
    integer:: iMisc , nGridOut
    !--------------------
    call init_rand()
    nCell_D = 1; nCell_D(1:nDim) = 2
    DxyzDomain_D      = 4
    DxyzCoarseBlock_D = 2
    DxyzFineBlock_D   = 1
    DxyzCoarse_D      = 1
    DxyzFine_D        = 0.5
    allocate(Xyz_DCB(nDim,nCell_D(1), nCell_D(2), nCell_D(3),&
         (2**nDim)*(2**nDim+1)))
    Xyz_DCB = 0
    do iGrid = 1, 2**nDim
       iBlock = iGrid
       XyzCorner_D = DxyzCoarseBlock_D*iShift_DI(1:nDim,iGrid)
       do k = 1, nCell_D(3)
          do j = 1, nCell_D(2)
             do i = 1, nCell_D(1)
                iCellIndex_D = (/i,j,k/)
                Xyz_DCB(:,i,j,k,iBlock) = XyzCorner_D +&
                    DxyzCoarse_D*(iCellIndex_D(1:nDim) - 0.50)
             end do
          end do
       end do
       do iSubGrid = 1, 2**nDim
          iBlock = iGrid*(2**nDim)+iSubGrid
          XyzCorner_D = DxyzCoarseBlock_D*iShift_DI(1:nDim,iGrid) + &
               DxyzFineBlock_D*iShift_DI(1:nDim,iSubGrid)
          do k = 1, nCell_D(3)
             do j = 1, nCell_D(2)
                do i = 1, nCell_D(1)
                   iCellIndex_D = (/i,j,k/)
                   Xyz_DCB(:,i,j,k,iBlock) = XyzCorner_D +&
                        DxyzFine_D*(iCellIndex_D(1:nDim) - 0.50)
                end do
             end do
          end do
       end do
    end do

    nIndexes = nDim +1
    CASE:do iCase = 0, 2**(2**nDim) - 2
       iLevelTest_I = 0; iGrid = 0
       iMisc = iCase
       do while(iMisc > 0)
          iGrid = iGrid + 1
          iLevelTest_I(iGrid) = mod(iMisc, 2)
          iMisc = (iMisc - iLevelTest_I(iGrid))/2
       end do
       write(*,*)'Case=',iLevelTest_I(1:2**nDim)
       !\
       ! We generated refinement, now sample points
       !/
       SAMPLE:do iSample = 1, nSample
          do iDir = 1, nDim
             Xyz_D(iDir) = (0.01 +0.98*rand())*DxyzDomain_D(iDir)
          end do
          ! write(*,*)'Xyz_D=',Xyz_D, ' case =', iLevelTest_I(1:2**nDim
          !\
          ! call interpolate_amr
          !/
          call interpolate_block_amr(&
               nDim=nDim, &
               XyzIn_D=Xyz_D, &
               nIndexes=nDim+1,&
               nCell_D=nCell_D(1:nDim),&
               find=find_test, &
               nGridOut=nGridOut,&
               Weight_I=Weight_I,&
               iIndexes_II=iIndexes_II,&
               iLevelOut_I=iLevelOut_I)
          !\          
          !Compare with interpolated:
          !/
          XyzInterpolated_D = 0
          do iGrid = 1, nGridOut
             iCellIndex_D = 1
             iCellIndex_D(1:nDim) = iIndexes_II(1:nDim,iGrid)
             iBlock = iIndexes_II(nIndexes,iGrid)
             XyzInterpolated_D = XyzInterpolated_D + Weight_I(iGrid)*&
                  Xyz_DCB(:,iCellIndex_D(1), iCellIndex_D(2), &
                  iCellIndex_D(3), iBlock)
          end do
          if(any(abs(Xyz_D - XyzInterpolated_D) > 1.0e-6).and.&
               all(iLevelOut_I/=BehindTheBoundary_))then
             write(*,*)'Approximation test failed'
             write(*,*)'Grid:', iLevelTest_I
             write(*,*)'nGridOut=',nGridOut
             write(*,*)'Point=', Xyz_D
             write(*,*)'Cell_D  iBlock XyzGrid_D Weight_I(iGrid)'
             do iGrid = 1, nGridOut
                iCellIndex_D = 1
                iCellIndex_D(1:nDim) = iIndexes_II(1:nDim,iGrid)
                iBlock = iIndexes_II(nIndexes,iGrid)
                write(*,*)iIndexes_II(1:nDim,iGrid), iBlock ,&
                     Xyz_DCB(:,iCellIndex_D(1), iCellIndex_D(2), &
                     iCellIndex_D(3), iBlock), Weight_I(iGrid)
             end do
             write(*,*)'Xyz_D=',Xyz_D
             write(*,*)'XyzInterpolated_D=',XyzInterpolated_D
             call CON_stop('Correct code and redo test')
          end if
       end do SAMPLE
    end do CASE
    deallocate(Xyz_DCB)
  end subroutine test_interpolate_amr
  !============================
  subroutine find_test(nDim, Xyz_D, &
            iProc, iBlock, XyzCorner_D, Dxyz_D, IsOut)
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
    real, dimension(nDim):: DxyzDomain_D, DxyzCoarseBlock_D ,&
         DxyzFineBlock_D, DxyzCoarse_D, DxyzFine_D
    integer:: iShift_D(3), iGrid, iSubGrid
    integer, dimension(0:1,0:1,0:1), parameter:: iGridFromShift_III=reshape(&
         (/1, 2, 3, 4, 5, 6, 7, 8/),(/2, 2, 2/))
    logical, dimension(nDim) :: IsAboveCenter_D
    !------------------- 
    DxyzDomain_D      = 4
    DxyzCoarseBlock_D = 2 
    DxyzFineBlock_D   = 1 
    DxyzCoarse_D      = 1
    DxyzFine_D        = 0.5
    iProc = 0; iBlock=0; XyzCorner_D=0.0; Dxyz_D = 0.0
    IsOut = any(Xyz_D < 0.0 .or. Xyz_D >= DxyzDomain_D)
    if(IsOut) RETURN
    !\
    ! Find into which coarse block the point fall
    !/ 
    IsAboveCenter_D = Xyz_D >= DxyzCoarseBlock_D
    iShift_D = 0
    where(IsAboveCenter_D)iShift_D(1:nDim) = 1
    XyzCorner_D = XyzCorner_D + DxyzCoarseBlock_D*iShift_D(1:nDim)
    Xyz_D       = Xyz_D       - DxyzCoarseBlock_D*iShift_D(1:nDim)
    iGrid = iGridFromShift_III(iShift_D(1),iShift_D(2),iShift_D(3))
    !\
    ! Check if the coarse block is used
    !/
    if(iLevelTest_I(iGrid)==0)then
       iBlock = iGrid
       Dxyz_D = DxyzCoarse_D
       RETURN
    end if
    !\
    ! The coarser block is refined, find into which fine block 
    ! the point falls
    !/
    IsAboveCenter_D = Xyz_D >= DxyzFineBlock_D
    iShift_D = 0
    where(IsAboveCenter_D)iShift_D(1:nDim) = 1
    XyzCorner_D = XyzCorner_D + DxyzFineBlock_D*iShift_D(1:nDim)
    Xyz_D       = Xyz_D       - DxyzFineBlock_D*iShift_D(1:nDim)
    iSubGrid = iGridFromShift_III(iShift_D(1),iShift_D(2),iShift_D(3))
    !write(*,*)'Fine: Xyz_D, XyzCorner_D, iSubGrid=', Xyz_D, XyzCorner_D, iSubGrid
    iBlock = iGrid*(2**nDim)+iSubGrid
    Dxyz_D = DxyzFine_D
  end subroutine find_test
end module ModInterpolateAMRGrid
