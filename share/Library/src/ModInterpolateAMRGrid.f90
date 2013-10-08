module ModInterpolateAMRGrid
  implicit none

  !Generalize bilinear and trilinear interpolation for AMR grids
  !The data are given at the cell-centered grid which consists of AMR blocks.
  !The interpolation is free of any jumps at the resolution interfaces
  !including edges and corners

  PRIVATE !Except
  SAVE
  public:: interpolate_amr_grid_1 ! One-dimensional - for completeness only

  ! Two-dimensional grid, is also used to interpolate the face projection for 3D grid

  public:: interpolate_amr_grid_2,test_interpolate_amr_grid_2


  ! Three-dimensional grid, 

  public :: interpolate_amr_grid_3, test_interpolate_amr_grid_3

  integer,parameter:: BehindTheBoundary_ = -7777, &
       Coarse_  = 0,               &
       Fine_    = 1,               &
       x_       = 1,               &
       y_       = 2,               &
       z_       = 3

  public :: BehindTheBoundary_

contains
  !===================================================================================
  subroutine  interpolate_amr_grid_1(&
       Xyz_D, XyzGrid_DI, iLevel_I,&
       nGridOut, Weight_I, iOrder_I)

    integer,parameter :: nGrid = 2, nDim = 1

    character(LEN=*),parameter:: NameSub='interpolate_amr_grid_1'

    !\
    !Input parameters
    !/
    real   , intent(in) :: Xyz_D(nDim) !The location at which to interpolate the data

    !Grid point coordinates
    real  ,   intent(in) :: XyzGrid_DI(nDim,nGrid) !1 coordinate, 2 points

    !The refinement level at each grid point. By one higher level of refinement
    !assumes the cell size reduced by a factor of 0.5
    integer,  intent(in) :: iLevel_I(nGrid)

    !\
    !Output parameters
    !/
    !The number of grid points to be ultimately included into the interpolation stencil
    !If nGridOut < nGridIn, only the first nGridOut elements are meaningful in the output
    !arrays   
    integer, intent(out) :: nGridOut  

    !The weight coefficients array.
    real   , intent(out) :: Weight_I(nGrid)

    !Order(numbers) of the grid points used for interpolation
    integer, intent(out) :: iOrder_I(nGrid)

    integer :: iGrid
    !-------------
    !\
    ! Check consistency of the input parameters.
    !/
    if(.not.(&
         XyzGrid_DI(x_,1).le.Xyz_D(x_).and.&
         XyzGrid_DI(x_,2) >  Xyz_D(x_)         )   )then
       write(*,*)'Inconsistent input in '//NameSub
       write(*,*)'Location at which to interpolate ', Xyz_D
       write(*,*)'Grid points:'
       do iGrid=1,nGrid
          write(*,*)XyzGrid_DI(:,iGrid)
       end do
       call CON_stop('Stop the code in '//NameSub)
    end if
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
       Weight_I(1) = ( XyzGrid_DI(1,2) - Xyz_D(1) )/( XyzGrid_DI(1,2) - XyzGrid_DI(1,1) )
       Weight_I(2) = 1 - Weight_I(1) 
    end if
  end subroutine interpolate_amr_grid_1
  !===================================================================================================================
  subroutine  interpolate_amr_grid_2(&
       Xyz_D, XyzGrid_DI, iLevel_I, &
       nGridOut, Weight_I, iOrder_I,&
       DoStencilFix, XyzStencil_D)

    integer,parameter :: nGrid = 4, nDim = 2

    character(LEN=*),parameter:: NameSub='interpolate_amr_grid_2'

    !\
    !Input parameters
    !/
    real   , intent(in) :: Xyz_D(nDim) !The location at which to interpolate the data

    !Grid point coordinates
    real  ,   intent(in) :: XyzGrid_DI(nDim,nGrid) !2 coordinate, 4 points

    !The refinement level at each grid point. By one higher level of refinement
    !assumes the cell size reduced by a factor of 0.5
    integer,  intent(inout) :: iLevel_I(nGrid)

    !\
    !Output parameters
    !/
    !The number of grid points to be ultimately included into the interpolation stencil
    !If nGridOut < nGridIn, only the first nGridOut lines are meaningful in the output
    integer, intent(out) :: nGridOut

    !The weight coefficients array.
    real   , intent(out) :: Weight_I(nGrid)

    !Order(numbers) of grid points used for the interpolation
    integer, intent(out) :: iOrder_I(nGrid)

    !\
    ! The logical is true, if (for some resolution combinations) the basic stencil
    ! for the point at which to inpertpolate is not applicable
    !/
    logical, intent(inout):: DoStencilFix
    !\
    ! If DoStencilFix==.true., the subroutine provides the XyzStencil_D to be
    ! used to construct stencil, not to interpolate
    !/
    real, intent(out) :: XyzStencil_D(nDim)


    integer :: iGrid, iDir   !Loop variables
    integer :: Loc(1)        !Find minloc and maxloc
    integer :: nGridOut1     !Number of grid points returned from 1D interpolation routine
    integer :: iLevelMin     !Minimumum refinement level

    real    :: XyzMin_D(nDim), XyzMax_D(nDim)    ! Minimum and maximum values of coordinates
    real    :: dXyzSmall_D(nDim)                 ! To find lines with the "constant" cordinate
    real    :: dXyzInv_D(  nDim)                 ! for calculating on a uniform grid
    real    :: Aux_D(nDim)                       ! Misc
    real    :: Xyz1_D(1), XyzGrid1_DI(1,2)       ! To call one-dimensional interpolation
    integer :: iOrder1_I(2)                      ! To call one-dimensional interpolation
    integer :: iLevel1_I(2)
    !-------------
    !\
    ! Check consistency of the input parameters.
    !/

    !\
    ! Make iLevel=0 for coarser grids and iLevel=1 for finer grids
    !/
    iLevelMin = minval(iLevel_I, MASK=iLevel_I/=BehindTheBoundary_)
    where(iLevel_I/=BehindTheBoundary_)iLevel_I = iLevel_I - iLevelMin
    if(DoStencilFix) goto 100
    if(.not.(&
         XyzGrid_DI(x_,1).le.Xyz_D(x_).and.XyzGrid_DI(y_,1).le.Xyz_D(y_).and.&
         XyzGrid_DI(x_,2)  > Xyz_D(x_).and.XyzGrid_DI(y_,2).le.Xyz_D(y_).and.&
         XyzGrid_DI(x_,3).le.Xyz_D(x_).and.XyzGrid_DI(y_,3)  > Xyz_D(y_).and.&
         XyzGrid_DI(x_,4)  > Xyz_D(x_).and.XyzGrid_DI(y_,4)  > Xyz_D(y_)) )then
       write(*,*)'Inconsistent input in '//NameSub
       write(*,*)'Location at which to interpolate ', Xyz_D
       write(*,*)'Grid points:'
       do iGrid=1,nGrid
          write(*,*)XyzGrid_DI(:,iGrid)
       end do
       write(*,*)'DoStencilFix ', DoStencilFix
       call CON_stop('Stop the code in '//NameSub)
    end if

    !\
    ! Check if the grid points are out of the grid boundary
    !/
    if(all(iLevel_I(3:4)== BehindTheBoundary_).or.all(XyzGrid_DI(y_,1:2)==Xyz_D(y_)) )then
       !\
       ! Edge 3:4 should be removed
       !/
       iOrder_I = (/1,2,3,4/)
       Xyz1_D(x_) = Xyz_D(x_) ; XyzGrid1_DI(x_,:) = XyzGrid_DI(x_,1:2)
       iLevel1_I = iLevel_I(1:2)

       call interpolate_amr_grid_1(&
            Xyz1_D,  XyzGrid1_DI, iLevel1_I, &
            nGridOut, Weight_I(1:2), iOrder1_I)
       iOrder_I(1:2) = iOrder_I(iOrder1_I)
       return

    elseif(all(iLevel_I(1:2)== BehindTheBoundary_)                                    )then
       !\
       ! Edge 1:2 should be removed
       !/
       iOrder_I = (/3,4,1,2/)
       Xyz1_D(x_) = Xyz_D(x_) ; XyzGrid1_DI(x_,:) = XyzGrid_DI(x_,3:4)
       iLevel1_I = iLevel_I(3:4)

       call interpolate_amr_grid_1(&
            Xyz1_D,  XyzGrid1_DI, iLevel1_I, &
            nGridOut, Weight_I(1:2), iOrder1_I)
       iOrder_I(1:2) = iOrder_I(iOrder1_I)
       return

    elseif(all(iLevel_I(2:4:2)== BehindTheBoundary_).or. all(XyzGrid_DI(x_,1:3:2)==Xyz_D(x_)) )then
       !\
       ! Edge 2,4 should be removed
       !/
       iOrder_I = (/1,3,2,4/)
       Xyz1_D(x_) = Xyz_D(y_) ; XyzGrid1_DI(x_,:) = XyzGrid_DI(y_,1:3:2)
       iLevel1_I = iLevel_I(1:3:2)

       call interpolate_amr_grid_1(&
            Xyz1_D,  XyzGrid1_DI, iLevel1_I, &
            nGridOut, Weight_I(1:2), iOrder1_I)
       iOrder_I(1:2) = iOrder_I(iOrder1_I)
       return

    elseif(all(iLevel_I(1:3:2)== BehindTheBoundary_)                                        )then
       !\
       ! Edge 1,3 should be removed
       !/
       iOrder_I = (/2,4,1,3/)
       Xyz1_D(x_) = Xyz_D(y_) ; XyzGrid1_DI(x_,:) = XyzGrid_DI(y_,2:4:2)
       iLevel1_I = iLevel_I(2:4:2)

       call interpolate_amr_grid_1(&
            Xyz1_D,  XyzGrid1_DI, iLevel1_I, &
            nGridOut, Weight_I(1:2), iOrder1_I)
       iOrder_I(1:2) = iOrder_I(iOrder1_I)
       return
    end if

    if(any(iLevel_I==BehindTheBoundary_) ) then
       !\
       ! There still can be a single grid point behind the grid boundary
       !/
       !  Boundary|Domain
       !          |
       !        C | C
       !          |FF
       !\
       ! In this configuration there is one grid behind the boundary, two fine and
       !  one coarse grid points. Otherwise BehindTheBoundary_ points are illegal.
       !/
       if(.not.(count(iLevel_I==Coarse_)==1&
            .and.count(iLevel_I==BehindTheBoundary_)==1))then
          write(*,*)'Wrong configuration of the grid boundary in '//NameSub
          do iGrid = 1,nGrid
             write(*,*)'iGrid=',iGrid,' iLevel_I(iGrid)=',iLevel_I(iGrid)
          end do
       end if
    end if


    !\
    ! Calculate the stencil size
    !/
    do iDir = 1, nDim
       dXyzSmall_D(iDir) = 0.10*(maxval(XyzGrid_DI(iDir,:)) -  minval(XyzGrid_DI(iDir,:)))
       dXyzInv_D(iDir)   = 0.10/dXyzSmall_D(iDir)
    end do

    !\
    ! Frist conbinations of the resolution interfaces: no resolution interface or
    ! a single resolition line
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
       return

    elseif(abs(XyzGrid_DI(y_,1) - XyzGrid_DI(y_,2)) < dXyzSmall_D(y_).and.&
         abs(XyzGrid_DI(y_,3) - XyzGrid_DI(y_,4)) < dXyzSmall_D(y_))then
       !\                             C  C
       ! Edges going along x-axis    F F   or   F F
       !/                                     C  C
       iOrder_I = (/1,2,3,4/)
       Aux_D(y_) = dXyzInv_D(y_)*(Xyz_D(y_)-XyzGrid_DI(y_,iOrder_I(1)))
       Xyz1_D(x_) = Xyz_D(x_)

       !\
       ! Interpolate along edge X1X2
       !/
       XyzGrid1_DI(x_,1) = XyzGrid_DI(x_,iOrder_I(1))
       XyzGrid1_DI(x_,2) = XyzGrid_DI(x_,iOrder_I(2))
       iLevel1_I(     1) = iLevel_I(     iOrder_I(1))
       iLevel1_I(     2) = iLevel_I(     iOrder_I(2))
       call interpolate_amr_grid_1(&
            Xyz1_D,  XyzGrid1_DI, iLevel1_I, &
            nGridOut1, Weight_I(1:2), iOrder1_I)
       iOrder_I(1:2) = iOrder_I(iOrder1_I)

       !\
       ! Apply weight for interpolation along another axis
       !/
       Weight_I(1:nGridOut1) =  Weight_I(1:nGridOut1) * (1 - Aux_D(y_))
       !\
       ! May need to remove a grid point from the stencil once behind the boundary
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
       call interpolate_amr_grid_1(&
            Xyz1_D,  XyzGrid1_DI, iLevel1_I, &
            nGridOut1, Weight_I(nGridOut+1:nGridOut+2), iOrder1_I)

       iOrder_I(nGridOut+1:nGridOut+2) = iOrder_I(nGridOut+iOrder1_I)
       !\
       ! Apply weight for interpolation along another axis
       !/
       Weight_I(nGridOut+1:nGridOut+nGridOut1) =  Weight_I(nGridOut+1:nGridOut+nGridOut1) * Aux_D(y_)

       !\
       ! May need to remove a grid point from the stencil once behind the boundary
       !/
       nGridOut = nGridOut + nGridOut1
       DoStencilFix = .false.
       return

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
       call interpolate_amr_grid_1(&
            Xyz1_D,  XyzGrid1_DI, iLevel1_I, &
            nGridOut1, Weight_I(1:2), iOrder1_I)
       iOrder_I(1:2) = iOrder_I(iOrder1_I)

       !\
       ! Apply weight for interpolation along another axis
       !/
       Weight_I(1:nGridOut1) =  Weight_I(1:nGridOut1) * (1 - Aux_D(x_))
       !\
       ! May need to remove a grid point from the stencil once behind the boundary
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
       call interpolate_amr_grid_1(&
            Xyz1_D,  XyzGrid1_DI, iLevel1_I, &
            nGridOut1, Weight_I(nGridOut+1:nGridOut+2), iOrder1_I)

       iOrder_I(nGridOut+1:nGridOut+2) = iOrder_I(nGridOut+iOrder1_I)
       !\
       ! Apply weight for interpolation along another axis
       !/
       Weight_I(nGridOut+1:nGridOut+nGridOut1) =  Weight_I(nGridOut+1:nGridOut+nGridOut1) * Aux_D(x_)

       !\
       ! May need to remove a grid point from the stencil once behind the boundary
       !/
       nGridOut = nGridOut + nGridOut1
       DoStencilFix = .false.
       return
    end if

    !\
    ! Check the ends of the resolution interfaces. Near the resolution interface endpoint
    ! the stencil may need to be reevaluated. Check fine edges along x_ direction
    !/
    if(abs(XyzGrid_DI(y_,2) - XyzGrid_DI(y_,1)) < dXyzSmall_D(y_).and.all(iLevel_I(1:2)==1))then

       if(abs(XyzGrid_DI(x_,3) - 0.50*(XyzGrid_DI(x_,2) + XyzGrid_DI(x_,1)))<dXyzSmall_D(x_))then
          !\                    C   (?)                                               C
          ! The configuration  F F    Check if the Xyz point belongs to the triangle F F
          !/
          iOrder_I = (/2,3,1,4/)
          !\
          ! Interpolate on the triangle X1,X2,X3 (F-C-F)
          ! vector X1-X3 is parallel to iAxis and directed toward 4th point of the stencil
          !   X3(F) - X1(F)
          !      \   /
          !      X2(C)  X4(?)
          !/
          call interpolate_triangle(iAxisFF=x_)
          return

       elseif(abs(XyzGrid_DI(x_,4) - 0.50*(XyzGrid_DI(x_,2) + XyzGrid_DI(x_,1)))<dXyzSmall_D(x_))then
          !\                (?) C                                                     C
          ! The configuration  F F    Check if the Xyz point belongs to the triangle F F
          !/
          iOrder_I = (/1,4,2,3/)
          !\
          ! Interpolate on the triangle X1,X2,X3 (F-C-F)
          ! vector X1-X3 is parallel to iAxis and directed toward 4th point of the stencil
          !   X3(F) - X1(F)
          !      \   /
          !      X2(C)  X4(?)
          !/
          call interpolate_triangle(iAxisFF=x_)
          return
       end if
    end if

    if(abs(XyzGrid_DI(y_,4) - XyzGrid_DI(y_,3)) < dXyzSmall_D(y_).and.all(iLevel_I(3:4)==1))then

       if(abs(XyzGrid_DI(x_,1) - 0.50*(XyzGrid_DI(x_,3) + XyzGrid_DI(x_,4)))<dXyzSmall_D(x_))then
          !\                   F F                                                      F F
          ! The configuration   C (?)    Check if the Xyz point belongs to the triangle  C
          !/
          iOrder_I = (/4,1,3,2/)
          !\
          ! Interpolate on the triangle X1,X2,X3 (F-C-F)
          ! vector X1-X3 is parallel to iAxis and directed toward 4th point of the stencil
          !   X3(F) - X1(F)
          !      \   /
          !      X2(C)  X4(?)
          !/
          call interpolate_triangle(iAxisFF=x_)
          return

       elseif(abs(XyzGrid_DI(x_,2) - 0.50*(XyzGrid_DI(x_,3) + XyzGrid_DI(x_,4)))<dXyzSmall_D(x_))then
          !\                     F F                                                   F F
          ! The configuration (?) C    Check if the Xyz point belongs to the triangle   C
          !/
          iOrder_I = (/3,2,4,1/)
          !\
          ! Interpolate on the triangle X1,X2,X3 (F-C-F)
          ! vector X1-X3 is parallel to iAxis and directed toward 4th point of the stencil
          !   X3(F) - X1(F)
          !      \   /
          !      X2(C)  X4(?)
          !/
          call interpolate_triangle(iAxisFF=x_)
          return
       end if
    end if

    !\
    ! Check the ends of the resolution interfaces. Near the resolution interface endpoint
    ! the stencil may need to be reevaluated. Check fine edges along y_ direction
    !/
    if(abs(XyzGrid_DI(x_,3) - XyzGrid_DI(x_,1)) < dXyzSmall_D(x_).and.all(iLevel_I((/1,3/))==Fine_))then

       if(abs(XyzGrid_DI(y_,2) - 0.50*(XyzGrid_DI(y_,3) + XyzGrid_DI(y_,1)))<dXyzSmall_D(y_))then
          !                     (?)
          !\                   F                                                     F
          ! The configuration    C    Check if the Xyz point belongs to the triangle   C
          !/                   F                                                     F
          iOrder_I = (/3,2,1,4/)
          !\
          ! Interpolate on the triangle X1,X2,X3 (F-C-F)
          ! vector X1-X3 is parallel to iAxis and directed toward 4th point of the stencil
          !   X3(F) - X1(F)
          !      \   /
          !      X2(C)  X4(?)
          !/
          call interpolate_triangle(iAxisFF=y_)
          return

       elseif(abs(XyzGrid_DI(y_,4) - 0.50*(XyzGrid_DI(y_,3) + XyzGrid_DI(y_,1)))<dXyzSmall_D(y_))then
          !\                   F                                                       F
          ! The configuration    C      Check if the Xyz point belongs to the triangle   C
          !/                   F                                                       F
          !                     (?)
          iOrder_I = (/1,4,3,2/)
          !\
          ! Interpolate on the triangle X1,X2,X3 (F-C-F)
          ! vector X1-X3 is parallel to iAxis and directed toward 4th point of the stencil
          !   X3(F) - X1(F)
          !      \   /
          !      X2(C)  X4(?)
          !/
          call interpolate_triangle(iAxisFF=y_)
          return
       end if
    end if

    if(abs(XyzGrid_DI(x_,4) - XyzGrid_DI(x_,2)) < dXyzSmall_D(x_).and.all(iLevel_I((/2,4/))==1))then

       if(abs(XyzGrid_DI(y_,1) - 0.50*(XyzGrid_DI(y_,2) + XyzGrid_DI(y_,4)))<dXyzSmall_D(y_))then
          !                    (?)
          !\                       F                                                       F
          ! The configuration   C       Check if the Xyz point belongs to the triangle  C
          !/                       F                                                       F
          iOrder_I = (/4,1,2,3/)
          !\
          ! Interpolate on the triangle X1,X2,X3 (F-C-F)
          ! vector X1-X3 is parallel to iAxis and directed toward 4th point of the stencil
          !   X3(F) - X1(F)
          !      \   /
          !      X2(C)  X4(?)
          !/
          call interpolate_triangle(iAxisFF=y_)
          return

       elseif(abs(XyzGrid_DI(y_,3) - 0.50*(XyzGrid_DI(y_,2) + XyzGrid_DI(y_,4)))<dXyzSmall_D(y_))then
          !\                       F                                                    F
          ! The configuration    C     Check if the Xyz point belongs to the triangle C
          !/                       F                                                    F
          !                     (?)
          iOrder_I = (/2,3,4,1/)
          !\
          ! Interpolate on the triangle X1,X2,X3 (F-C-F)
          ! vector X1-X3 is parallel to iAxis and directed toward 4th point of the stencil
          !   X3(F) - X1(F)
          !      \   /
          !      X2(C)  X4(?)
          !/
          call interpolate_triangle(iAxisFF=y_)
          return
       end if
    end if

100 continue


    DoStencilFix = .false.

    if(count(iLevel_I==Fine_)==1.and.count(iLevel_I==Coarse_)==3)then
       Loc = maxloc(iLevel_I)
       select case(Loc(1))
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
       Loc = minloc(iLevel_I)
       select case(Loc(1))
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


    else
       write(*,*)'Algorithm failuire in '//NameSub
       do iGrid = 1,nGrid
          write(*,*)'iGrid=',iGrid,' iLevel_I(iGrid)=',iLevel_I(iGrid),' XyzGrid_DI=', XyzGrid_DI(:,iGrid)
       end do
       write(*,*)'Xyz_D=', Xyz_D(:)!ADDED BY DMITRY
       call CON_stop('Stop in '//NameSub)

    end if

    !\
    !Points 1 and 2 are on the shared side of the triangles, 3 and 4 are off
    !/
    call triangulate

  contains
    real function cross_product(a_D, b_D)
      real,dimension(nDim),intent(in) :: a_D, b_D
      !----------
      cross_product = a_D(x_)* b_D(y_) - a_d(y_)*b_D(x_)
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
      Alpha3 = cross_product(Xyz_D - X1_D, X2_D-X1_D)/cross_product(X3_D - X1_D, X2_D-X1_D)
      if(Alpha3==0.0)then
         nGridOut = 2
         Alpha2 = cross_product(X3_D-X1_D,Xyz_D - X1_D)/cross_product(X3_D - X1_D, X2_D-X1_D)
         Weight_I(2) = Alpha2
      elseif(Alpha3 > 0.0)then
         nGridOut = 3
         Alpha2 = cross_product(X3_D-X1_D,Xyz_D - X1_D)/cross_product(X3_D - X1_D, X2_D-X1_D)
         Weight_I(2) = Alpha2
         Weight_I(3) = Alpha3
      else
         nGridOut = 3
         Alpha2 = cross_product(X4_D - X1_D,  Xyz_D - X1_D)/cross_product(X4_D - X1_D, X2_D-X1_D)
         Alpha3 = cross_product(Xyz_D - X1_D, X2_D  - X1_D)/cross_product(X4_D - X1_D, X2_D-X1_D)
         iOrder_I(3:4) = iOrder_I((/4,3/))
         Weight_I(2) = Alpha2
         Weight_I(3) = Alpha3
      end if

      Weight_I(1) = 1 - sum(Weight_I(2:nGridOut))

    end subroutine triangulate
    !======================================================================================================================
    subroutine interpolate_triangle(iAxisFF)
      !\
      ! Interpolate on the triangle X1,X2,X3 (F-C-F)
      ! vector X1-X3 is parallel to iAxis and directed toward 4th point of the stencil
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

      Weight_I(3) = cross_product(Xyz_D - X1_D, X2_D  - X1_D)/cross_product(X3_D - X1_D, X2_D-X1_D)
      Weight_I(2) = cross_product(X3_D  - X1_D, Xyz_D - X1_D)/cross_product(X3_D - X1_D, X2_D-X1_D)
      if(Weight_I(3)==0.0)then
         nGridOut = 2
      else
         nGridOut = 3
      end if
      Weight_I(1) = 1 - sum(Weight_I(2:3))
      DoStencilFix = any(Weight_I(1:3)<0)
      if(DoStencilFix)then
         iAxisPerp = 3 - iAxisFF
         XyzStencil_D(iAxisFF) =     X1_D(iAxisFF  ) + X2_D(iAxisFF  ) - X3_D(iAxisFF  )
         XyzStencil_D(iAxisPerp) = ( X1_D(iAxisPerp) + X2_D(iAxisPerp) + X3_D(iAxisPerp) )/3
      end if
    end subroutine interpolate_triangle
  end subroutine interpolate_amr_grid_2

  !===================================================================================================
  subroutine  interpolate_amr_grid_3(&
       Xyz_D, XyzGrid_DI, iLevel_I, &
       nGridOut, Weight_I, iOrder_I,&
       DoStencilFix, XyzStencil_D)

    integer,parameter :: nGrid = 8, nDim = 3


    character(LEN=*),parameter:: NameSub='interpolate_amr_grid_3'

    !\
    !Input parameters
    !/
    real   , intent(in) :: Xyz_D(nDim) !The location at which to interpolate the data

    !Grid point coordinates
    real  ,   intent(in) :: XyzGrid_DI(nDim,nGrid) !3 coordinate, 8 points

    !The refinement level at each grid point. By one higher level of refinement
    !assumes the cell size reduced by a factor of 0.5
    integer,  intent(inout) :: iLevel_I(nGrid)

    !\
    !Output parameters
    !/
    !The number of grid points to be ultimately included into the interpolation stencil
    !If nGridOut < nGridIn, only the first nGridOut lines in the output arrays
    !are meaningful.
    integer, intent(out) :: nGridOut

    !The weight coefficients array.
    real   , intent(out) :: Weight_I(nGrid)

    !Order(numbers) of grid points used for the interpolation
    integer, intent(out) :: iOrder_I(nGrid)

    !\
    ! The logical is true, if (for some resolution combinations) the basic stencil
    ! for the point at which to inpertpolate is not applicable
    !/
    logical, intent(out):: DoStencilFix
    !\
    ! If DoStencilFix==.true., the subroutine provides the XyzStencil_D to be
    ! used to construct stencil, not to interpolate!
    !/
    real, intent(out) :: XyzStencil_D(nDim)

    !Parameters to enumerate faces or face diagonal
    integer, parameter:: xy_ = 3, xz_ = 2, yz_ = 1

    !Number of the vertex connected by 
    !the edge of direction iDir (second index) 
    !with the given vertex (first index)
    integer, dimension(nGrid,nDim), parameter:: iEdge_ID = &
         reshape((/&   !Number of the connected vertex
         2 ,1, 4, 3, 6, 5, 8, 7,        & !Edge x
         3, 4, 1, 2, 7, 8, 5, 6,        & !Edge y
         5, 6, 7, 8, 1, 2, 3, 4         & !Edge z
         /),(/nGrid,nDim/))

    !Number of the vertex connected by 
    !the face diagonal across the face of direction iDir (second index) 
    !with the given vertex (first index)
    integer, dimension(nGrid,nDim), parameter:: iFaceDiag_ID = &
         reshape((/&   !Number of the connected vertex
         7 ,8, 5, 6, 3, 4, 1, 2,        & !Face yz
         6, 5, 8, 7, 2, 1, 4, 3,        & !Face xz
         4, 3, 2, 1, 8, 7, 6, 5         & !Face xy
         /),(/nGrid,nDim/))

    !Number of the vertex connected by 
    !the main diagonal with the given vertex (index)
    integer, dimension(nGrid), parameter:: iMainDiag_I = &
         (/ 8, 7, 6, 5, 4, 3, 2, 1 /)

    !Vertexes (enumerated by the first undex), which form 
    !the face of direction iDir (second index) including the 
    !given vertex (the third index). 
    !When the first index equals 1,2,3,4, the vertex, accordingly: 
    !(1) coincides with the given one
    !(2) is connected with the given one by the edge of direction iDir+1
    !(3) is connected with the given one by the edge of direction iDir+2
    !(4) is connected to the given one by the face diagonal of direction iDir
    integer, dimension(nGrid/2,nDim,nGrid), parameter:: iFace_IDI = &
         reshape((/&   ! yz face   ! xz face      ! xy face ! 
         1, 3, 5, 7,   1, 5, 2, 6,   1, 2, 3, 4, & 
         2, 4, 6, 8,   2, 6, 1, 5,   2, 1, 4, 3, &
         3, 1, 7, 5,   3, 7, 4, 8,   3, 4, 1, 2, &
         4, 2, 8, 6,   4, 8, 3, 7,   4, 3, 2, 1, &
         5, 7, 1, 3,   5, 1, 6, 2,   5, 6, 7, 8, &
         6, 8, 2, 4,   6, 2, 5, 1,   6, 5, 8, 7, &
         7, 5, 3, 1,   7, 3, 8, 4,   7, 8, 5, 6, &
         8, 6, 4, 2,   8, 4, 7, 3,   8, 7, 6, 5  &
         /), (/nGrid/2,nDim,nGrid/))

    !Vertexes (enumerated by the first undex), which form 
    !the face of direction iDir (second index) and does not include the
    !given vertex (the third index). 
    !When the first index equals 1,2,3,4, the vertex, accordingly: 
    !(1) is connected to the given one by the edge of direction iDir
    !(2) is connected to (1) by the edge of direction iDir+1
    !(3) is connected to (1) by the edge of direction iDir+2
    !(4) is connected to  by the face diagonal of direction iDir
    integer, dimension(nGrid/2,nDim,nGrid), parameter:: iOppositeFace_IDI = &
         reshape((/&   ! yz face   ! xz face      ! xy face ! 
         2, 4, 6, 8,   3, 7, 4, 8,   5, 6, 7, 8, & 
         1, 3, 5, 7,   4, 8, 3, 7,   6, 5, 8, 7, &
         4, 2, 8, 6,   1, 5, 2, 6,   7, 8, 5, 6, &
         3, 1, 7, 5,   2, 6, 1, 5,   8, 7, 6, 5, &
         6, 8, 2, 4,   7, 3, 8, 4,   1, 2, 3, 4, &
         5, 7, 1, 3,   8, 4, 7, 3,   2, 1, 4, 3, &
         8, 6, 4, 2,   5, 1, 6, 2,   3, 4, 1, 2, &
         7, 5, 3, 1,   6, 2, 5, 1,   4, 3, 2, 1  &
         /), (/nGrid/2,nDim,nGrid/))

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
    integer :: Loc(1)

    integer :: nGridOut1     !Number of grid points returned from 1D interpolation routine


    !\
    ! Minimum and maximum values of coordinates
    !/
    real    :: XyzMin_D(nDim), XyzMax_D(nDim)   


    real    :: dXyzSmall_D(nDim)                 ! To find points with  "equal" cordinates
    real    :: dXyzInv_D(  nDim)                 ! for calculating on a uniform grid
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

    if(DoStencilFix) goto 300
    !\
    ! Check consistency of the input parameters.
    !/
    if(.not.(&
         XyzGrid_DI(x_,1).le.Xyz_D(x_).and.XyzGrid_DI(y_,1).le.Xyz_D(y_)&
         .and.XyzGrid_DI(z_,1).le.Xyz_D(z_).and.&
         XyzGrid_DI(x_,2)  > Xyz_D(x_).and.XyzGrid_DI(y_,2).le.Xyz_D(y_)&
         .and.XyzGrid_DI(z_,2).le.Xyz_D(z_).and.&
         XyzGrid_DI(x_,3).le.Xyz_D(x_).and.XyzGrid_DI(y_,3)  > Xyz_D(y_)&
         .and.XyzGrid_DI(z_,3).le.Xyz_D(z_).and.&
         XyzGrid_DI(x_,4)  > Xyz_D(x_).and.XyzGrid_DI(y_,4)  > Xyz_D(y_)&
         .and.XyzGrid_DI(z_,4).le.Xyz_D(z_).and.&
         XyzGrid_DI(x_,5).le.Xyz_D(x_).and.XyzGrid_DI(y_,5).le.Xyz_D(y_)&
         .and.XyzGrid_DI(z_,5)  > Xyz_D(z_).and.&
         XyzGrid_DI(x_,6)  > Xyz_D(x_).and.XyzGrid_DI(y_,6).le.Xyz_D(y_)&
         .and.XyzGrid_DI(z_,6)  > Xyz_D(z_).and.&
         XyzGrid_DI(x_,7).le.Xyz_D(x_).and.XyzGrid_DI(y_,7)  > Xyz_D(y_)&
         .and.XyzGrid_DI(z_,7)  > Xyz_D(z_).and.&
         XyzGrid_DI(x_,8)  > Xyz_D(x_).and.XyzGrid_DI(y_,8)  > Xyz_D(y_)&
         .and.XyzGrid_DI(z_,8)  > Xyz_D(z_))  )then
       write(*,*)'Inconsistent input in '//NameSub
       write(*,*)'Location at which to interpolate ', Xyz_D
       write(*,*)'Grid points:'
       do iGrid=1,nGrid
          write(*,*)XyzGrid_DI(:,iGrid)
       end do
       call CON_stop('Stop the code in '//NameSub)
    end if
300 continue

    !\
    ! Check if the grid points are out of the grid boundary
    !/
    do iDir = 1, nDim
       if(  all(iLevel_I(iOppositeFace_IDI(:,iDir,1))== BehindTheBoundary_).or.&
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
       call interpolate_amr_grid_2(&
            Xyz2_D,  XyzGrid2_DI, iLevel2_I, &
            nGridOut, Weight_I(1:4), iOrder2_I, DoStencilFix2, XyzStencil2_D)

       if(DoStencilFix2)then
          DoStencilFix = .true.
          XyzStencil_D(iDir)    = Xyz_D(iDir)
          XyzStencil_D(1 + mod(iDir  ,3)) = XyzStencil2_D(x_)
          XyzStencil_D(1 + mod(iDir+1,3)) = XyzStencil2_D(y_)
       end if

       iOrder_I(1:4) = iOrder_I(iOrder2_I)
       return
    end do


    if(any(iLevel_I==BehindTheBoundary_) ) then

       if(count(iLevel_I==Coarse_           ) /= &
            count(iLevel_I==BehindTheBoundary_)    &
            )then
          write(*,*)'Wrong configuration of the grid boundary in '//NameSub
          do iGrid = 1,nGrid
             write(*,*)'iGrid=',iGrid,' iLevel_I(iGrid)=',iLevel_I(iGrid)
          end do
       end if
    end if


    !\
    ! Calculate the stencil size
    !/
    do iDir = 1, nDim
       dXyzSmall_D(iDir) = 0.10*(maxval(XyzGrid_DI(iDir,:)) -  minval(XyzGrid_DI(iDir,:)))
       dXyzInv_D(iDir)   = 0.10/dXyzSmall_D(iDir)
    end do

    !\
    ! First conbinations of the resolution interfaces: no resolution interface or
    ! a single resolition line
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
       return
    end if

    !\
    ! Opposite faces going along resolution interface
    !/
    do iDir = 1,3
       !Face and opposite face are planes
       if(all(abs(&
            XyzGrid_DI(iDir,iFace_IDI(1,iDir,1)) - XyzGrid_DI(iDir,iFace_IDI(:,iDir,1))) < dXyzSmall_D(iDir))&
            .and.all(abs(&
            XyzGrid_DI(iDir,iOppositeFace_IDI(1,iDir,1)) - XyzGrid_DI(iDir,iOppositeFace_IDI(:,iDir,1)))&
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
          call interpolate_amr_grid_2(&
               Xyz2_D, XyzGrid2_DI, iLevel2_I, &
               nGridOut2, Weight_I(1:4), iOrder2_I, DoStencilFix2, XyzStencil2_D)
          iOrder_I(1:4) = iOrder_I(iOrder2_I)
          !\
          ! Apply weight for interpolation along z-axis
          !/
          Aux_D(iDir) = dXyzInv_D(iDir)*(Xyz_D(iDir) - XyzGrid_DI(iDir,iOrder_I(1)))
          Weight_I(1:nGridOut2) = Weight_I(1:nGridOut2)*(1 - Aux_D(iDir))

          !\
          ! May need to remove a grid points from the stencil once behind the boundary
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

          XyzGrid2_DI(x_,1:4) = XyzGrid_DI(1 + mod(iDir    ,3),iOrder_I(nGridOut2+1:nGridOut2+4))
          XyzGrid2_DI(y_,1:4) = XyzGrid_DI(1 + mod(iDir + 1,3),iOrder_I(nGridOut2+1:nGridOut2+4))
          iLevel2_I(        1:4) = iLevel_I(        iOrder_I(nGridOut2+1:nGridOut2+4))
          call interpolate_amr_grid_2(&
               Xyz2_D, XyzGrid2_DI, iLevel2_I, &
               nGridOut2, Weight_I(nGridOut+1:nGridOut+4), iOrder2_I, DoStencilFix2, XyzStencil2_D)
          do iGrid=1,nGridOut2
             iOrder_I(nGridOut+iGrid)=iOrder_I(nGridOut+iOrder2_I(iGrid))
          end do

          !\
          ! Apply weight for interolation along kAxis
          !/
          Weight_I(nGridOut+1:nGridOut+nGridOut2) = Weight_I(nGridOut+1:nGridOut+nGridOut2) * Aux_D(iDir)

          !\
          ! May need to remove a grid points from the stencil once behind the boundary
          !/
          nGridOut = nGridOut + nGridOut2
          DoStencilFix = .false.
          return
       end if
    end do


    !\
    ! Edges going along resolution edge
    ! Opposite "faces" have the same geometry in projection
    ! term "face" is used to refer to 4 upper-, left-, front-most etc. grid points
    ! they do not belong to plane
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
          ! Interpolation weight along iDir axis are calculated seperately for
          ! Coarse and fine points
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
          XyzGrid2_DI(x_:y_,1:4) = XyzGrid_DI((/1 + mod(iDir,3),1 + mod(iDir+1,3)/),iOrder_I(1:4))
          iLevel2_I(        1:4) = iLevel_I( iOrder_I(1:4))
          where(iLevel2_I==BehindTheBoundary_)iLevel2_I=Coarse_
          call interpolate_amr_grid_2(&
               Xyz2_D, XyzGrid2_DI, iLevel2_I, &
               nGridOut2, Weight_I(1:4), iOrder2_I, DoStencilFix2, XyzStencil2_D)

          if(DoStencilFix2)then
             DoStencilFix = DoStencilFix2
             XyzStencil_D((/1 + mod(iDir,3),1 + mod(iDir+1,3)/)) = XyzStencil2_D(x_:y_)
             XyzStencil_D(iDir) = Xyz_D(iDir)
             return
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
                if(    iLevel_I(iOrder_I(          iGrid))==BehindTheBoundary_)then
                   Weight_I(nGridOut2+iGrid) = Weight_I(iGrid)
                   Weight_I(          iGrid) = 0
                elseif(iLevel_I(iOrder_I(nGridOut2+iGrid))==BehindTheBoundary_)then
                   Weight_I(nGridOut2+iGrid) = 0
                   Weight_I(          iGrid) = Weight_I(iGrid)
                else
                   Weight_I(nGridOut2+iGrid) = Weight_I(iGrid)*     AuxCoarse
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
                Loc = minloc(iLevel_I(iOrder_I(1:nGridOut-1)))
                Weight_I(Loc(1)) = Weight_I(nGridOut)
                iOrder_I(Loc(1)) = iOrder_I(nGridOut)
                nGridOut = nGridOut -1
             end if
          end do
          DoStencilFix = .false.
          return
       end if
    end do

    !\
    ! Corner transition junction
    !/
    do iGrid = 1,nGrid;do iDir = 1,nDim
       if( abs(        XyzGrid_DI(1 + mod(iDir    ,3),iGrid) - &
            0.25*sum(XyzGrid_DI(1 + mod(iDir    ,3),iOppositeFace_IDI(:,iDir,iGrid)))) < &
            dXyzSmall_D(1 + mod(iDir    ,3)).and.&
            abs(        XyzGrid_DI(1 + mod(iDir + 1,3),iGrid) - &
            0.25*sum(XyzGrid_DI(1 + mod(iDir + 1,3),iOppositeFace_IDI(:,iDir,iGrid)))) < &
            dXyzSmall_D(1 + mod(iDir + 1,3)).and.&
            (count(iLevel_I==Coarse_)==3))then
          !\
          ! Junction around direction iDir
          !/
          iOrder_I(1:4) = iFace_IDI(:,iDir,iGrid)
          iOrder_I(5:8) = iOppositeFace_IDI(:,iDir,iGrid)
          call interpolate_corner_transition_junction(iDir)
          return
       end if
    end do;end do

    !\
    ! Corner transition
    !/
    do iGrid = 1,nGrid;do iDir = 1,nDim
       if( count(iLevel_I(iFace_IDI(:,iDir,iGrid))==Fine_)>0 .and. &
            all(abs(XyzGrid_DI(iDir,iGrid)                            - &
            0.50*(XyzGrid_DI(iDir,iFace_IDI(        :,iDir,iGrid))  + &
            XyzGrid_DI(iDir,iOppositeFace_IDI(:,iDir,iGrid))) * &
            (iLevel_I(iFace_IDI(:,iDir,iGrid))-Coarse_)-&
            XyzGrid_DI(iDir,iFace_IDI(         :,iDir,iGrid)) * &
            (Fine_-iLevel_I(iFace_IDI(:,iDir,iGrid)))) < &
            dXyzSmall_D(iDir)))then
          !\
          ! Transition in direction iDir
          !/
          iOrder_I(1:4) = iFace_IDI(:,iDir,iGrid)
          iOrder_I(5:8) = iOppositeFace_IDI(:,iDir,iGrid)
          call interpolate_corner_transition(iDir)
          return
       end if
    end do;end do



    DoStencilFix = .false.
    !---------------------Resolution corners----------------------
    !Among 255 legal combinations
    !(=2**nGrid - illegal "all Fine configurtion")
    !we considered above one all-coarse configuration, + 6 resolution faces
    !+ 30 resolution edges, totally 37 left 218
    !
    !
    !\
    ! One configuration resulting in
    ! the corner split for five tetrahedra
    ! is occured in many cases (26), however, it is treated in the same way in all these cases
    ! if the given grid point is connected by three edges with three coarse points, while the across
    ! the main diagonal point is connected by three edges (or, equivalently, the given point is connected
    ! by three face diagonals) with three coarse points), then all faces are trangulated   
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
               iTetrahedron1_I=(/ iGrid,&
               iEdge_ID(iGrid,x_),iFaceDiag_ID(iGrid,xz_),iFaceDiag_ID(iGrid,xy_)/),&
               iTetrahedron2_I=(/ iGrid,&
               iEdge_ID(iGrid,y_),iFaceDiag_ID(iGrid,yz_),iFaceDiag_ID(iGrid,xy_)/),&
               iTetrahedron3_I=(/ iGrid,&
               iEdge_ID(iGrid,z_),iFaceDiag_ID(iGrid,yz_),iFaceDiag_ID(iGrid,xz_)/),&
               iTetrahedron4_I=(/ iGrid,&
               iFaceDiag_ID(iGrid,yz_),iFaceDiag_ID(iGrid,xz_),iFaceDiag_ID(iGrid,xy_)/),&
               iTetrahedron5_I=(/ iMainDiag_I(iGrid),&
               iFaceDiag_ID(iGrid,yz_),iFaceDiag_ID(iGrid,xz_),iFaceDiag_ID(iGrid,xy_)/))
          return
       end if
    end do


    select case( count(iLevel_I==Fine_))
    case(1)                                 ! 8 cases , totally 71 left 184
       !Find the only fine grid
       Loc = maxloc(iLevel_I)
       !\
       ! View from the fine point
       !            xy face diag
       !   x-axis      |      y-axis
       !         \     V     /
       !          \   C4    /           Connections to opposite faces: main diag and 
       !           C2/ | \C3            3478 - y edge, yz face diag, xy face daig
       !           |\  | / |            2468 - x edge, xy face diag, xz face diag
       !           | \ |/  |            5687 - z face. yz face diag, xy face diag
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
            iOppositeFace_IDI(1,x_,Loc(1)),&
            iOppositeFace_IDI(2,x_,Loc(1)),&
            iOppositeFace_IDI(3,x_,Loc(1)),&
            iOppositeFace_IDI(4,x_,Loc(1)),&
            Loc(1)/)                      ,&
            iRectangular2_I=(/&
            iOppositeFace_IDI(1,y_,Loc(1)),&
            iOppositeFace_IDI(2,y_,Loc(1)),&
            iOppositeFace_IDI(3,y_,Loc(1)),&
            iOppositeFace_IDI(4,y_,Loc(1)),&
            Loc(1)/)                      ,&
            iRectangular3_I=(/&
            iOppositeFace_IDI(1,z_,Loc(1)),&
            iOppositeFace_IDI(2,z_,Loc(1)),&
            iOppositeFace_IDI(3,z_,Loc(1)),&
            iOppositeFace_IDI(4,z_,Loc(1)),&
            Loc(1)/))
       return
    case(7)                                 ! 8 cases , totally 79 left 176
       !=============ARGUEBLE CASE!!!!!!!!!!!!!!!!=====================================
       !Find the only coarse grid
       Loc = minloc(iLevel_I)
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
            iTetrahedron1_I=(/Loc(1),iEdge_ID(Loc(1),x_),iEdge_ID(Loc(1),y_),iEdge_ID(Loc(1),z_)/))
       if(nGridOut>1)return
       call parallel_rays(Dir_D=XyzGrid_DI(:,iMainDiag_I(Loc(1))) - XyzGrid_DI(:,Loc(1)), &
            iURectangle1_I=iFace_IDI(:,x_,iMainDiag_I(Loc(1))),&
            iURectangle2_I=iFace_IDI(:,y_,iMainDiag_I(Loc(1))),&
            iURectangle3_I=iFace_IDI(:,z_,iMainDiag_I(Loc(1))),&
            iDTriangle1_I=(/iEdge_ID(Loc(1),x_),iEdge_ID(Loc(1),y_),iEdge_ID(Loc(1),z_)/),&
            iDTriangle2_I=(/iEdge_ID(Loc(1),x_),iEdge_ID(Loc(1),y_),iFaceDiag_ID(Loc(1),xy_)/),&
            iDTriangle3_I=(/iEdge_ID(Loc(1),x_),iEdge_ID(Loc(1),z_),iFaceDiag_ID(Loc(1),xz_)/),&
            iDTriangle4_I=(/iEdge_ID(Loc(1),y_),iEdge_ID(Loc(1),z_),iFaceDiag_ID(Loc(1),yz_)/))
       return
    case(2)
       !Find the fine grid
       Loc = maxloc(iLevel_I)
       if(iLevel_I(iMainDiag_I(Loc(1)))==Fine_)then           ! 4 cases   totally 83 left 172

          !Full traingulation of all faces, each coarse vertex is connected by the edge 
          !with one fine point and with the face diagonal to the other. 
          !\
          ! View from the fine point
          !            xy face diag
          !   x-axis      |      y-axis
          !         \     V     /
          !          \   C4    /           
          !           C2/ | \C3            The point F8 is also fine 
          !           |\  | / |            The line of sight goes along the
          !           | \ |/  |            main diagonal F1F8 which is the common side 
          !           |  \/   |            of 6 tetrahedron
          !           |   F1  |
          !  xz facediag/ | \yz face diag
          !           /   |  \|
          !          C6   |   C7
          !             \ C5/
          !               |
          !               V
          !               z-axis
          iGrid = iMainDiag_I(Loc(1)) !Fine point across the main diagonal
          call pyramids(&
               iTetrahedron1_I=(/Loc(1) ,  &
               iEdge_ID(    Loc(1),x_ ) ,  &
               iFaceDiag_ID(Loc(1),xy_) ,  &
               iGrid/),                    &
               iTetrahedron2_I=(/Loc(1) ,  &
               iEdge_ID(    Loc(1),x_ ) ,  &
               iFaceDiag_ID(Loc(1),xz_) ,  &
               iGrid/),                    &
               iTetrahedron3_I=(/Loc(1) ,  &
               iEdge_ID(    Loc(1),y_ ) ,  &
               iFaceDiag_ID(Loc(1),xy_) ,  &
               iGrid/),                    &
               iTetrahedron4_I=(/Loc(1) ,  &
               iEdge_ID(    Loc(1),y_ ) ,  &
               iFaceDiag_ID(Loc(1),yz_) ,  &
               iGrid/),                    &
               iTetrahedron5_I=(/Loc(1) ,  &
               iEdge_ID(    Loc(1),z_ ) ,  &
               iFaceDiag_ID(Loc(1),xz_) ,  &
               iGrid/),                    &
               iTetrahedron6_I=(/Loc(1) ,  &
               iEdge_ID(    Loc(1),z_ ) ,  &
               iFaceDiag_ID(Loc(1),yz_) ,  &
               iGrid/))
          return
       else                                            ! 12 cases  totally 95  left 160
          do iDir = yz_,xy_
             if(iLevel_I(iFaceDiag_ID(Loc(1), iDir))==Fine_)then
                iOrder_I(1:4) = iFace_IDI(:,iDir,Loc(1))
                iOrder_I(5:8) = iOppositeFace_IDI(:,iDir,Loc(1))
                !
                !    7C     8C 
                !  5C     6C     
                !
                !    3C     4F
                !  1F     2C
                !
                ! 1 Remove tetrahedra 1F 4F 2C 6C and 1F 4F 3C 7C
                call pyramids(&
                     iTetrahedron1_I=(/&
                     iOrder_I(1), iOrder_I(4),iOrder_I(2), iOrder_I(6)/),&
                     iTetrahedron2_I=(/&
                     iOrder_I(1), iOrder_I(4),iOrder_I(3), iOrder_I(7)/))
                if(nGridOut > 1)return
                !
                !    7C------ 8C 
                !  5C--\---6C/|     
                !   |  _ \_  \|
                !   | /    --4F
                !  1F----/     
                !\
                ! Center of the lower face
                !/
                XyzMin_D = 0.50*(XyzGrid_DI(:,iOrder_I(1))+XyzGrid_DI(:,iOrder_I(4)))
                !\
                ! Center of the upper face
                !/
                XyzMax_D = 0.50*(XyzGrid_DI(:,iOrder_I(5))+XyzGrid_DI(:,iOrder_I(8)))
                call parallel_rays(Dir_D=XyzMax_D - XyzMin_D,&
                     iURectangle1_I=iOrder_I(5:8),&
                     iDTriangle1_I=(/iOrder_I(1),iOrder_I(4),iOrder_I(6)/),&
                     iDTriangle2_I=(/iOrder_I(1),iOrder_I(4),iOrder_I(7)/),&
                     iDTriangle3_I=(/iOrder_I(4),iOrder_I(6),iOrder_I(8)/),&
                     iDTriangle4_I=(/iOrder_I(4),iOrder_I(7),iOrder_I(8)/),&
                     iDTriangle5_I=(/iOrder_I(1),iOrder_I(5),iOrder_I(7)/),&
                     iDTriangle6_I=(/iOrder_I(1),iOrder_I(5),iOrder_I(6)/))
                !=============ARGUEABLE CASE!!!!!!!!!!!!!!!!=====================
                return
             end if
          end do
       end if
    case(6)
       !=============ARGUEABLE CASE!!!!!!!!!!!!!!!!=====================
       !Find the coarse point
       Loc = minloc(iLevel_I)
       if(iLevel_I(iMainDiag_I(Loc(1)))==Coarse_)then           ! 4 cases   totally 99 left 156
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
          iGrid = iMainDiag_I(Loc(1))
          call pyramids(&
               iTetrahedron1_I=(/Loc(1),iEdge_ID(Loc(1),x_),iEdge_ID(Loc(1),y_),iEdge_ID(Loc(1),z_)/),&
               iTetrahedron2_I=(/iGrid ,iEdge_ID(iGrid,x_ ),iEdge_ID(iGrid,y_ ),iEdge_ID(iGrid,z_)/))
          if(nGridOut>1)return
          call parallel_rays(Dir_D=XyzGrid_DI(:,iGrid) - XyzGrid_DI(:,Loc(1)), &
               iUTriangle1_I=(/iEdge_ID(iGrid ,x_),iEdge_ID(iGrid ,y_),iEdge_ID(iGrid ,z_)/),&
               iUTriangle2_I=(/iEdge_ID(iGrid ,x_),iEdge_ID(iGrid ,y_),iFaceDiag_ID(iGrid ,xy_)/),&
               iUTriangle3_I=(/iEdge_ID(iGrid ,x_),iEdge_ID(iGrid ,z_),iFaceDiag_ID(iGrid ,xz_)/),&
               iUTriangle4_I=(/iEdge_ID(iGrid ,y_),iEdge_ID(iGrid ,z_),iFaceDiag_ID(iGrid ,yz_)/),&
               iDTriangle1_I=(/iEdge_ID(Loc(1),x_),iEdge_ID(Loc(1),y_),iEdge_ID(Loc(1),z_)/),&
               iDTriangle2_I=(/iEdge_ID(Loc(1),x_),iEdge_ID(Loc(1),y_),iFaceDiag_ID(Loc(1),xy_)/),&
               iDTriangle3_I=(/iEdge_ID(Loc(1),x_),iEdge_ID(Loc(1),z_),iFaceDiag_ID(Loc(1),xz_)/),&
               iDTriangle4_I=(/iEdge_ID(Loc(1),y_),iEdge_ID(Loc(1),z_),iFaceDiag_ID(Loc(1),yz_)/))
          return
       else                                            ! 12 cases  totally 111  left 144
          do iDir = yz_,xy_
             if(iLevel_I(iFaceDiag_ID(Loc(1), iDir))==Coarse_)then
                iOrder_I(1:4) = iFace_IDI(:,iDir,Loc(1))
                iOrder_I(5:8) = iOppositeFace_IDI(:,iDir,Loc(1))
                !
                !    7F     8F 
                !  5F     6F     
                !
                !    3F     4C
                !  1C     2F
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
                   XyzMin_D = 0.50*(XyzGrid_DI(:,iOrder_I(1))+XyzGrid_DI(:,iOrder_I(4)))
                   !\
                   ! Center of the upper face
                   !/
                   XyzMax_D = 0.50*(XyzGrid_DI(:,iOrder_I(5))+XyzGrid_DI(:,iOrder_I(8)))
                   call parallel_rays(Dir_D=XyzMax_D - XyzMin_D,&
                        iURectangle1_I=iOrder_I(5:8),&
                        iDTriangle1_I=(/iOrder_I(2),iOrder_I(3),iOrder_I(5)/),&
                        iDTriangle2_I=(/iOrder_I(2),iOrder_I(3),iOrder_I(8)/))
                end if
                return
             end if
          end do
       end if

    case(3)
       !                   24 cases   totally 135 left 120
       do iGrid = 1, nGrid
          if(iLevel_I(iGrid)==Fine_)then
             if(iLevel_I(iMainDiag_I(iGrid))==Fine_)then
                do iDir = yz_,xy_
                   if(iLevel_I(iFaceDiag_ID(iGrid, iDir))==Fine_)then
                      iOrder_I(1:4) = iFace_IDI(:,iDir,iGrid)
                      iOrder_I(5:8) = iOppositeFace_IDI(:,iDir,iGrid)
                      !
                      !                                  F4 F8
                      !\       C7   F8        Trapezoid: C3 C7       
                      !     C5   C6           
                      !                       Trapezoid: F4 F8
                      !        C3   F4                   C2 C6
                      !     F1   C2           Common apex of the trapezoidal pyramids: F1
                      !                       Tetrahedra which are left F1 C5 F8 C6 and F1 C5 F8 C7 
                      !/
                      call pyramids(&
                           iTetrahedron1_I=(/iOrder_I(1),iOrder_I(5),iOrder_I(8),iOrder_I(6)/),&
                           iTetrahedron2_I=(/iOrder_I(1),iOrder_I(5),iOrder_I(8),iOrder_I(7)/),&
                           iTrapezoidal1_I=(/iOrder_I(2),iOrder_I(6),iOrder_I(4),iOrder_I(8),iOrder_I(1)/),&
                           iTrapezoidal2_I=(/iOrder_I(3),iOrder_I(7),iOrder_I(4),iOrder_I(8),iOrder_I(1)/))
                      return
                   end if
                end do
             else
                !\
                ! Three fine cells in one plane   24 cases, totally  159 left 96
                !/
                !  !!!!!!ARGUEABLE CASE!!!!!!!!!!!!!!!!!!
                do iDir = 1,nDim
                   if(all(iLevel_I(iFace_IDI(1:3,iDir,iGrid))==Fine_))then
                      iOrder_I(1:4) = iFace_IDI(:,iDir,iGrid)
                      iOrder_I(5:8) = iOppositeFace_IDI(:,iDir,iGrid)
                      !
                      !                                  F1 F3
                      !\       C7   C8        Trapezoid: C5 C7       
                      !     C5   C6           
                      !                       Trapezoid: F1 F2
                      !        F3   C4                   C5 C6
                      !     F1   F2            
                      !/
                      !=============
                      call pyramids(&
                           iTetrahedron1_I=(/iOrder_I(2),iOrder_I(3),iOrder_I(8),iOrder_I(4)/))
                      if(nGridOut> 1)return
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
                      call parallel_rays(Dir_D=XyzGrid_DI(:,iOrder_I(8)) - XyzGrid_DI(:,iOrder_I(4)),&
                           iDTriangle1_I=(/iOrder_I(1),iOrder_I(2),iOrder_I(3)/),&
                           iDTriangle2_I=(/iOrder_I(8),iOrder_I(2),iOrder_I(3)/),&
                           iDTriangle3_I=(/iOrder_I(8),iOrder_I(7),iOrder_I(3)/),&
                           iDTriangle4_I=(/iOrder_I(8),iOrder_I(2),iOrder_I(6)/),&
                           iURectangle1_I=(/iOrder_I(5),iOrder_I(6),iOrder_I(7),iOrder_I(8)/) )    

                   end if
                end do
             end if
          end if
       end do
    case(5)
       !=============ARGUEABLE CASE!!!!!!!!!!!!!!!!=====================
       !               24 cases, totally  183 left 72
       do iGrid = 1, nGrid
          if(iLevel_I(iGrid)==Coarse_)then
             if(iLevel_I(iMainDiag_I(iGrid))==Coarse_)then
                do iDir = yz_,xy_
                   if(iLevel_I(iFaceDiag_ID(iGrid, iDir))==Coarse_)then
                      iOrder_I(1:4) = iFace_IDI(:,iDir,iGrid)
                      iOrder_I(5:8) = iOppositeFace_IDI(:,iDir,iGrid)
                      !
                      !                                  F3 F7
                      !\       F7   C8        Trapezoid: C4 C8       
                      !     F5   F6           
                      !                       Trapezoid: F2 F6
                      !        F3   C4                   C4 C8
                      !     C1   F2            
                      !/
                      !=============
                      call pyramids(&
                           iTetrahedron1_I=(/iOrder_I(1),iOrder_I(2),iOrder_I(3),iOrder_I(5)/),&
                           iRectangular1_I=(/iOrder_I(2),iOrder_I(3),iOrder_I(6),iOrder_I(7),iOrder_I(5)/))
                      if(nGridOut> 1)return

                      Xyz2_D(x_:y_) = Xyz_D((/1 + mod(iDir,3),1 + mod(1 + iDir,3)/))
                      !\
                      ! Interpolation weights along iDir axis are calculated separately for
                      ! Coarse and fine points
                      !/
                      !\
                      ! may need to rearrange iOrder
                      !/
                      if(XyzGrid_DI(iDir,iOrder_I(1))>XyzGrid_DI(iDIr,iOrder_I(5)))then
                         iOrder_I(1:8) = iOrder_I((/5,6,7,8,1,2,3,4/))
                      end if
                      XyzMin_D(iDir)   = minval(XyzGrid_DI(iDir,:),iLevel_I/=Fine_)
                      AuxCoarse  = (Xyz_D(iDir) - XyzMin_D(iDir))/(&
                           maxval(XyzGrid_DI(iDir,:),iLevel_I/=Fine_) - XyzMin_D(iDir))

                      XyzMin_D(iDir)   = minval(XyzGrid_DI(iDir,:),iLevel_I==Fine_)
                      AuxFine  = (Xyz_D(iDir) - XyzMin_D(iDir))/(&
                           maxval(XyzGrid_DI(iDir,:),iLevel_I==Fine_) - XyzMin_D(iDir))
                      !\
                      ! Interpolate along lower face
                      !/
                      DoStencilFix2 = .true.
                      XyzGrid2_DI(x_:y_,1:4) = XyzGrid_DI((/1 + mod(iDir,3),1 + mod(iDir+1,3)/),iOrder_I(1:4))
                      iLevel2_I(        1:4) = iLevel_I( iOrder_I(1:4))

                      call interpolate_amr_grid_2(&
                           Xyz2_D, XyzGrid2_DI, iLevel2_I, &
                           nGridOut2, Weight_I(1:4), iOrder2_I, DoStencilFix2, XyzStencil2_D)

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
                            Weight_I(nGridOut2+jGrid) = Weight_I(jGrid)*     AuxFine
                            Weight_I(          jGrid) = Weight_I(jGrid)*(1-  AuxFine)
                         else
                            Weight_I(nGridOut2+jGrid) = Weight_I(jGrid)*     AuxCoarse
                            Weight_I(          jGrid) = Weight_I(jGrid)*(1 - AuxCoarse)
                         end if
                      end do

                      !                      call parallel_rays(Dir_D=XyzGrid_DI(:,iOrder_I(8)) - XyzGrid_DI(:,iOrder_I(4)),&
                      !                           iDTriangle1_I=(/iOrder_I(2),iOrder_I(3),iOrder_I(4)/),&
                      !                           iUTriangle2_I=(/iOrder_I(6),iOrder_I(7),iOrder_I(8)/) )
                      return
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
                   !                                  F5 F7
                   !\       F7   F8        Trapezoid: C1 C3       
                   !     F5   F6           
                   !                       Trapezoid: F5 F6
                   !        C3   F4                   C1 C2
                   !     C1   C2           Rectangle: F5 F6 F7 F8
                   !/               
                   !    Common apex F4
                   call pyramids(&
                        iRectangular1_I=(/iOrder_I(7),iOrder_I(8),iOrder_I(5),iOrder_I(6),iOrder_I(4)/),&
                        iTrapezoidal1_I=(/iOrder_I(1),iOrder_I(2),iOrder_I(5),iOrder_I(6),iOrder_I(4)/),&
                        iTrapezoidal2_I=(/iOrder_I(1),iOrder_I(3),iOrder_I(5),iOrder_I(7),iOrder_I(4)/))
                   return
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
                if(all(iLevel_I(iOppositeFace_IDI(2:4,iDir,iGrid))==Coarse_))then
                   iOrder_I(1:4) = iFace_IDI(:,iDir,iGrid)
                   iOrder_I(5:8) = iOppositeFace_IDI(:,iDir,iGrid)
                   !\
                   ! face 1:4 has only one Coarse point stored in iOrder_I(1)
                   ! face 5:8 has three Coarse points in iOrder((/6,7,8/))
                   !/
                   !-----------------------------------------------
                   !                                            F2 F4
                   !                  C7    C8        Trapezoid C6 C8
                   !                                           
                   !                F5    C6
                   !                                            F3 F4
                   !                                  Trapezoid C7 C8
                   !                  F3    F4                 
                   !
                   !                C1    F2
                   !Common apex F5
                   call pyramids(&
                        iTrapezoidal1_I=(/iOrder_I(6),iOrder_I(8),iOrder_I(2),iOrder_I(4),iOrder_I(5)/),&
                        iTrapezoidal2_I=(/iOrder_I(7),iOrder_I(8),iOrder_I(3),iOrder_I(4),iOrder_I(5)/),&
                        iTetrahedron1_I=(/            iOrder_I(2),iOrder_I(3),iOrder_I(4),iOrder_I(5)/),&
                        iTetrahedron2_I=(/iOrder_I(1),iOrder_I(2),iOrder_I(3),            iOrder_I(5)/) )
                   return
                end if
                !
                !\
                ! 24 cases totally 255 left 0
                !/
                if(iLevel_I(iEdge_ID(iGrid,iDir))==Coarse_.and.&
                     iLevel_I(iEdge_ID(iGrid,1 + mod(iDir,3)))==Coarse_.and.&
                     iLevel_I(iEdge_ID(iEdge_ID(iGrid,iDir),1 + mod(1 + iDir,3)))==Coarse_)then
                   iOrder_I(1:4) = iFace_IDI(:,iDir,iGrid)
                   iOrder_I(5:8) = iOppositeFace_IDI(:,iDir,iGrid)
                   !\
                   !Chain of four coarse points connected with mutually orthogonal edges
                   !/
                   !-----------------------------------------------
                   !                                            F6 F8
                   !                  C7    F8        Trapezoid C5 C7
                   !                                           
                   !                C5    F6
                   !                                            F3 F4
                   !                                  Trapezoid C1 C2
                   !                  F3    F4                 
                   !
                   !                C1    C2
                   call pyramids(&
                        iTrapezoidal1_I=(/iOrder_I(1),iOrder_I(2),iOrder_I(3),iOrder_I(4),iOrder_I(6)/),&  !see #1 below
                        iTetrahedron1_I=(/            iOrder_I(8),iOrder_I(3),iOrder_I(4),iOrder_I(6)/),&  !see #2 below
                        iTetrahedron2_I=(/iOrder_I(1),iOrder_I(6),iOrder_I(3),            iOrder_I(5)/),&  !see #3 below
                        iTrapezoidal2_I=(/iOrder_I(5),iOrder_I(7),iOrder_I(6),iOrder_I(8),iOrder_I(3)/) )
                   ! #1 The leftover is above subfaces 136 and 364
                   ! #2 The leftover is above subfaces 136 and 368
                   ! #3 The leftover is above subfaces 365 and 368
                   return
                end if
             end do
          end if
       end do

    end select
    write(*,*)'Algorithm failure in '//NameSub
    do iGrid = 1,nGrid
       write(*,*)'iGrid=',iGrid,'iLevel(iGrid)=',iLevel_I(iGrid),'XyzGrid_DI=', XyzGrid_DI(:, iGrid)
    end do
    call CON_stop('Stop in '//NameSub)


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
            return
         elseif(present(iTetrahedron2_I))then
            call tetrahedron(&
                 XyzGrid_DI(:,iTetrahedron2_I(1)),&
                 XyzGrid_DI(:,iTetrahedron2_I(2)),&
                 XyzGrid_DI(:,iTetrahedron2_I(3)),&
                 XyzGrid_DI(:,iTetrahedron2_I(4)))
            if(all(Weight_I(1:4).ge.0.0))then
               nGridOut = 4
               iOrder_I(1:4) = iTetrahedron2_I
               return
            elseif(present(iTetrahedron3_I))then
               call tetrahedron(&
                    XyzGrid_DI(:,iTetrahedron3_I(1)),&
                    XyzGrid_DI(:,iTetrahedron3_I(2)),&
                    XyzGrid_DI(:,iTetrahedron3_I(3)),&
                    XyzGrid_DI(:,iTetrahedron3_I(4)))
               if(all(Weight_I(1:4).ge.0.0))then
                  nGridOut = 4
                  iOrder_I(1:4) = iTetrahedron3_I
                  return
               elseif(present(iTetrahedron4_I))then
                  call tetrahedron(&
                       XyzGrid_DI(:,iTetrahedron4_I(1)),&
                       XyzGrid_DI(:,iTetrahedron4_I(2)),&
                       XyzGrid_DI(:,iTetrahedron4_I(3)),&
                       XyzGrid_DI(:,iTetrahedron4_I(4)))
                  if(all(Weight_I(1:4).ge.0.0))then
                     nGridOut = 4
                     iOrder_I(1:4) = iTetrahedron4_I
                     return
                  elseif(present(iTetrahedron5_I))then
                     call tetrahedron(&
                          XyzGrid_DI(:,iTetrahedron5_I(1)),&
                          XyzGrid_DI(:,iTetrahedron5_I(2)),&
                          XyzGrid_DI(:,iTetrahedron5_I(3)),&
                          XyzGrid_DI(:,iTetrahedron5_I(4)))
                     if(all(Weight_I(1:4).ge.0.0))then
                        nGridOut = 4
                        iOrder_I(1:4) = iTetrahedron5_I
                        return
                     elseif(present(iTetrahedron6_I))then
                        call tetrahedron(&
                             XyzGrid_DI(:,iTetrahedron6_I(1)),&
                             XyzGrid_DI(:,iTetrahedron6_I(2)),&
                             XyzGrid_DI(:,iTetrahedron6_I(3)),&
                             XyzGrid_DI(:,iTetrahedron6_I(4)))
                        if(all(Weight_I(1:4).ge.0.0))then
                           nGridOut = 4
                           iOrder_I(1:4) = iTetrahedron6_I
                           return
                        end if ! 6
                     end if    ! 5
                  end if       ! 4
               end if          ! 3
            end if             ! 2
         end if                ! 1
      end if                   ! no tetrahedron
      if(present(iRectangular1_I))then
         call pyramid(Base_=Rectangular_,&
              X1_D=XyzGrid_DI(:,iRectangular1_I(1)),&
              X2_D=XyzGrid_DI(:,iRectangular1_I(2)),&
              X3_D=XyzGrid_DI(:,iRectangular1_I(3)),&
              X4_D=XyzGrid_DI(:,iRectangular1_I(4)),&
              X5_D=XyzGrid_DI(:,iRectangular1_I(5)) )
         if(all(Weight_I(1:5).ge.0.0))then
            nGridOut = 5
            iOrder_I(1:5) = iRectangular1_I
            return
         elseif(present(iRectangular2_I))then
            call pyramid(Base_=Rectangular_,&
                 X1_D=XyzGrid_DI(:,iRectangular2_I(1)),&
                 X2_D=XyzGrid_DI(:,iRectangular2_I(2)),&
                 X3_D=XyzGrid_DI(:,iRectangular2_I(3)),&
                 X4_D=XyzGrid_DI(:,iRectangular2_I(4)),&
                 X5_D=XyzGrid_DI(:,iRectangular2_I(5)) )
            if(all(Weight_I(1:5).ge.0.0))then
               nGridOut = 5
               iOrder_I(1:5) = iRectangular2_I
               return
            elseif(present(iRectangular3_I))then
               call pyramid(Base_=Rectangular_,&
                    X1_D=XyzGrid_DI(:,iRectangular3_I(1)),&
                    X2_D=XyzGrid_DI(:,iRectangular3_I(2)),&
                    X3_D=XyzGrid_DI(:,iRectangular3_I(3)),&
                    X4_D=XyzGrid_DI(:,iRectangular3_I(4)),&
                    X5_D=XyzGrid_DI(:,iRectangular3_I(5)) )
               if(all(Weight_I(1:5).ge.0.0))then
                  nGridOut = 5
                  iOrder_I(1:5) = iRectangular3_I
                  return
               end if          ! 3
            end if             ! 2
         end if                ! 1
      end if                   ! no rectangular
      if(present(iTrapezoidal1_I))then
         call pyramid(Base_=Trapezoidal_,&
              X1_D=XyzGrid_DI(:,iTrapezoidal1_I(1)),&
              X2_D=XyzGrid_DI(:,iTrapezoidal1_I(2)),&
              X3_D=XyzGrid_DI(:,iTrapezoidal1_I(3)),&
              X4_D=XyzGrid_DI(:,iTrapezoidal1_I(4)),&
              X5_D=XyzGrid_DI(:,iTrapezoidal1_I(5)) )
         if(all(Weight_I(1:5).ge.0.0))then
            nGridOut = 5
            iOrder_I(1:5) = iTrapezoidal1_I
            return
         elseif(present(iTrapezoidal2_I))then
            call pyramid(Base_=Trapezoidal_,&
                 X1_D=XyzGrid_DI(:,iTrapezoidal2_I(1)),&
                 X2_D=XyzGrid_DI(:,iTrapezoidal2_I(2)),&
                 X3_D=XyzGrid_DI(:,iTrapezoidal2_I(3)),&
                 X4_D=XyzGrid_DI(:,iTrapezoidal2_I(4)),&
                 X5_D=XyzGrid_DI(:,iTrapezoidal2_I(5)) )
            if(all(Weight_I(1:5).ge.0.0))then
               nGridOut = 5
               iOrder_I(1:5) = iTrapezoidal2_I
               return
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
      ! AlphaUp = (UX1-Xyz)\cdot[(UX2-UX1)\times(UX3-UX1)]/Dir\cdot \cdot[(UX2-UX1)\times(UX3-UX1)]
      ! XyzUp belongs to the up subface
      !/
      !\
      ! Solve equation:
      ! Xyz - Dir*AlphaDown = XyzDown 
      ! its solution is 
      ! AlphaDown = (Xyz-DX1)\cdot[(DX2-DX1)\times(DX3-DX1)]/Dir\cdot \cdot[(DX2-DX1)\times(DX3-DX1)]
      ! XyzUp belongs to the Down subface
      !/
      !\
      ! As long as Xyz=(XyzDown*AlphaUp + XyzUp*AlphaDown)/(AlphaUp+AlphaDown)
      ! the weghts for the upper face should be multiplied by AlphaDown/(AlphaUp+AlphaDown)
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
               return
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
                  return
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
               AlphaUp = triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
                    triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
               if(AlphaUp==0.0)then
                  call triangle(X1_D, X2_D, X3_D, Xyz_D)
                  if(all(Weight_I(1:3)>=0.0))then
                     iOrder_I(1:3) = iUTriangle3_I
                     nGridOut = 3
                     return
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
                  AlphaUp = triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
                       triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
                  if(AlphaUp==0.0)then
                     call triangle(X1_D, X2_D, X3_D, Xyz_D)
                     if(all(Weight_I(1:3)>=0.0))then
                        iOrder_I(1:3) = iUTriangle4_I
                        nGridOut = 3
                        return
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
                     AlphaUp = triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
                          triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
                     if(AlphaUp==0.0)then
                        call triangle(X1_D, X2_D, X3_D, Xyz_D)
                        if(all(Weight_I(1:3)>=0.0))then
                           iOrder_I(1:3) = iUTriangle5_I
                           nGridOut = 3
                           return
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
                        AlphaUp = triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
                             triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
                        if(AlphaUp==0.0)then
                           call triangle(X1_D, X2_D, X3_D, Xyz_D)
                           if(all(Weight_I(1:3)>=0.0))then
                              iOrder_I(1:3) = iUTriangle6_I
                              nGridOut = 3
                              return
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
               return
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
                  return
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
               AlphaUp = triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
                    triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
               if(AlphaUp==0.0)then
                  call rectangle(X1_D, X2_D, X3_D, X4_D, Xyz_D)
                  if(all(Weight_I(1:4)>=0.0))then
                     iOrder_I(1:4) = iURectangle3_I
                     nGridOut = 4
                     return
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
      return
700   continue
      !\
      !Calculate low face
      if(present(iDTriangle1_I))then
         Weight_I(4:3+nGridOutUp) = Weight_I(1:nGridOutUp)
         iOrder_I(4:3+nGridOutUp) = iOrder_I(1:nGridOutUp)
         X1_D = XyzGrid_DI(:,iDTriangle1_I(1))
         X2_D = XyzGrid_DI(:,iDTriangle1_I(2))
         X3_D = XyzGrid_DI(:,iDTriangle1_I(3))
         AlphaDown = triple_product(Xyz_D - X1_D, X2_D - X1_D, X3_D - X1_D)/&
              triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
         if(AlphaDown==0.0)then
            call triangle(X1_D, X2_D, X3_D, Xyz_D)
            if(all(Weight_I(1:3)>=0.0))then
               iOrder_I(1:3) = iDTriangle1_I
               nGridOut = 3
               return
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
            AlphaDown = triple_product(Xyz_D - X1_D, X2_D - X1_D, X3_D - X1_D)/&
                 triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
            if(AlphaDown==0.0)then
               call triangle(X1_D, X2_D, X3_D, Xyz_D)
               if(all(Weight_I(1:3)>=0.0))then
                  iOrder_I(1:3) = iDTriangle2_I
                  nGridOut = 3
                  return
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
               AlphaDown = triple_product(Xyz_D - X1_D, X2_D - X1_D, X3_D - X1_D)/&
                    triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
               if(AlphaDown==0.0)then
                  call triangle(X1_D, X2_D, X3_D, Xyz_D)
                  if(all(Weight_I(1:3)>=0.0))then
                     iOrder_I(1:3) = iDTriangle3_I
                     nGridOut = 3
                     return
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
                  AlphaDown = triple_product(Xyz_D - X1_D, X2_D - X1_D, X3_D - X1_D)/&
                       triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
                  if(AlphaDown==0.0)then
                     call triangle(X1_D, X2_D, X3_D, Xyz_D)
                     if(all(Weight_I(1:3)>=0.0))then
                        iOrder_I(1:3) = iDTriangle3_I
                        nGridOut = 3
                        return
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
                     AlphaDown = triple_product(Xyz_D - X1_D, X2_D - X1_D, X3_D - X1_D)/&
                          triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
                     if(AlphaDown==0.0)then
                        call triangle(X1_D, X2_D, X3_D, Xyz_D)
                        if(all(Weight_I(1:3)>=0.0))then
                           iOrder_I(1:3) = iDTriangle5_I
                           nGridOut = 3
                           return
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
                        AlphaDown = triple_product(Xyz_D - X1_D, X2_D - X1_D, X3_D - X1_D)/&
                             triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
                        if(AlphaDown==0.0)then
                           call triangle(X1_D, X2_D, X3_D, Xyz_D)
                           if(all(Weight_I(1:3)>=0.0))then
                              iOrder_I(1:3) = iDTriangle6_I
                              nGridOut = 3
                              return
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
               return
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
            AlphaDown = triple_product(Xyz_D - X1_D, X2_D - X1_D, X3_D - X1_D)/&
                 triple_product(Dir_D, X2_D - X1_D, X3_D - X1_D)
            if(AlphaDown==0.0)then
               call rectangle(X1_D, X2_D, X3_D, X4_D, Xyz_D)
               if(all(Weight_I(1:4)>=0.0))then
                  iOrder_I(1:4) = iDRectangle2_I
                  nGridOut = 4
                  return
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
               AlphaDown = triple_product(Xyz_D - X1_D, X2_D - X1_D, X3_D - X1_D)/&
                    triple_product(Dir_D, X2_D - X1_D, X3_D - X1_D)
               if(AlphaDown==0.0)then
                  call rectangle(X1_D, X2_D, X3_D, X4_D, Xyz_D)
                  if(all(Weight_I(1:4)>=0.0))then
                     iOrder_I(1:4) = iDRectangle3_I
                     nGridOut = 4
                     return
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
      return
800   continue
      !Apply the weight for interpolation along the rate
      nGridOut = nGridOutUp + nGridOutDown
      Weight_I(1:nGridOutDown) = Weight_I(1:nGridOutDown)*AlphaUp/(AlphaUp + AlphaDown)
      Weight_I(nGridOutDown+1:nGridOut) = Weight_I(nGridOutDown+1:nGridOut)*AlphaDown/(AlphaUp + AlphaDown)
    end subroutine parallel_rays

    !========================
    function cross_product(a_D, b_d)
      real, dimension(nDim), intent(in) :: a_d, b_D
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
      ! Xyz = Weight_I(1)*X1 + Weight_I(2)*X2 + Weight_I(3)*X3 + Weight_I(4)*X4
      ! Or, which is equivalent:
      ! Xyz -X1 = Weight_I(2)*(X2-X1) + Weight_I(3)*(X3-X1) + Weight_I(4)*(X4-X1)
      !/
      Aux= 1/triple_product( X4_D - X1_D,  X3_D - X1_D,  X2_D - X1_D)
      Weight_I(4) = triple_product(Xyz_D - X1_D,  X3_D - X1_D,  X2_D - X1_D)*Aux
      Weight_I(3) = triple_product( X4_D - X1_D, Xyz_D - X1_D,  X2_D - X1_D)*Aux
      Weight_I(2) = triple_product( X4_D - X1_D,  X3_D - X1_D, Xyz_D - X1_D)*Aux
      Weight_I(1) = 1 - sum(Weight_I(2:4))
    end subroutine tetrahedron
    !=======================================================================
    subroutine pyramid(Base_,X1_D, X2_D, X3_D, X4_D, X5_D)
      !\
      ! Interpolate in the pyramid with base X1X2X3X4 and apex X5
      ! valid for case of rectangular or trapezoid base
      !/
      integer, intent(in) :: Base_
      real, dimension(nDim), intent(in) :: X1_D, X2_D, X3_D, X4_D, X5_D
      real, dimension(nDim) :: XyzP_D !Projection of Xyz point on the base of the pyramid X1X2X3X4 along the line X5Xyz
      real ::Alpha5
      !-----------------------------------------------
      !\
      ! Solve equation: X5 + (Xyz - X5)/Alpha5 = XyzP
      ! where XyzP belongs to the pydamid base.
      ! As long as (XyzP-X1_D)\cdot[(X2-X1)\times(X3-X1)]=0,
      ! we have Alpha5 =(X5 -Xyz)\cdot[(X2-X1)\times(X3-X1)]/(X5 -X1)\cdot[(X2-X1)\times(X3-X1)]
      !/

      Alpha5 = triple_product(X5_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
           triple_product(X5_D - X1_D , X2_D - X1_D, X3_D - X1_D)

      if(Alpha5==0.0)then
         if(all(Xyz_D(:)==X5_D(:)))then
            Weight_I(5) = 1
            Weight_I(1:4) = 0
            return
         else
            Weight_I(1:5) = -1
            return
         end if
      elseif(Alpha5 < 0.0 .or. Alpha5 > 1.0)then
         Weight_I(1:5) = -1
         return
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
      select case(Base_)
      case(Rectangular_)
         call rectangle(X1_D, X2_D, X3_D, X4_D, XyzP_D)
      case(Trapezoidal_)
         call trapezoid(X1_D, X2_D, X3_D, X4_D, XyzP_D)
      end select

      if(any(Weight_I(1:4) < 0.0))return
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
         return
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
      ! Weight_I(3) = cross_product(XyzP_D - X1_D, X2_D  - X1_D)/cross_product(X3_D - X1_D, X2_D-X1_D)
      ! Weight_I(2) = cross_product(X3_D  - X1_D, Xyz_D - X1_D)/cross_product(X3_D - X1_D, X2_D-X1_D)
      ! iN 3d case instead of /cross_product(X3_D - X1_D, X2_D-X1_D) we use 
      ! \cdot cross_product(X3_D - X1_D, X2_D-X1_D)/sum(cross_product(X3_D - X1_D, X2_D-X1_D)**2)
      !/
      CrossProductInv_D = cross_product(X3_D - X1_D, X2_D-X1_D)
      CrossProductInv_D = CrossProductInv_D/sum(CrossProductInv_D**2)

      Weight_I(3) = sum(cross_product(XyzP_D - X1_D, X2_D  - X1_D)*CrossProductInv_D)
      Weight_I(2) = sum(cross_product(X3_D  - X1_D, XyzP_D - X1_D)*CrossProductInv_D)
      Weight_I(1) = 1 - sum(Weight_I(2:3))

    end subroutine triangle

    !=======================================================================
    subroutine interpolate_corner_transition_junction(iEdgeDir)
      integer, intent(in) :: iEdgeDir
      integer             :: iAxis,jAxis,kAxis
      real :: AlphaKAxis
      real :: dXyzUp, dXyzDown ! Auxilary variables to determine vertical weight
      !\
      ! Subroutine intepolates in junction of transitional regions near a resolution corner\
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
      !------------------------------------------------------------------------
      kAxis = iEdgeDir
      iAxis = 1 + mod(iEdgeDir    ,3)
      jAxis = 1 + mod(iEdgeDir + 1,3)
      Xyz2_D(x_:y_) = Xyz_D((/iAxis,jAxis/))
      !write(*,*)'junction'
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
      call interpolate_amr_grid_2(&
           Xyz2_D, XyzGrid2_DI, iLevel2_I, &
           nGridOut2, Weight_I(5:8), iOrder2_I, DoStencilFix2, XyzStencil2_D)
      !\
      ! Rearrange points according to 2D interpolation output
      !/
      iOrder_I(1:4) = iOrder_I(iOrder2_I  )
      iOrder_I(5:8) = iOrder_I(iOrder2_I+4)
      nGridOut = nGridOut2

      !\
      ! find Fine point on 1:4 and substitute it with corresponding point from face 5:8
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
      call interpolate_amr_grid_2(&
           Xyz2_D, XyzGrid2_DI, iLevel2_I, &
           nGridOut2, Weight_I(1:4), iOrder2_I, DoStencilFix2, XyzStencil2_D)
      !\
      ! Rearrange points according to 2D interpolation output
      !/
      iOrder_I(1:4) = iOrder_I(iOrder2_I  )
      iOrder_I(5:8) = iOrder_I(iOrder2_I+4)
      Weight_I(5:8) = Weight_I(iOrder2_I+4)

      !\
      ! move Fine point on 1:4 to position 1 and corresponding point to position 5
      !/
      if(iLevel_I(iOrder_I(2))==Fine_)then
         iOrder_I((/1,2,5,6/)) = iOrder_I((/2,1,6,5/))
         Weight_I((/1,2,5,6/)) = Weight_I((/2,1,6,5/))
      end if

      !\
      ! find vertical weight
      !/
      dXyzUp   = Xyz_D(kAxis) - XyzGrid_DI(kAxis,iOrder_I(5))
      dXyzDown = Xyz_D(kAxis) - sum(XyzGrid_DI(kAxis,iOrder_I(1:nGridOut2))*Weight_I(1:nGridOut2))

      !vertical_distance_point_plane(Xyz_D,(/1,2,3/),kAxis)

      DoStencilFix = dXyzUp * dXyzDown > 0.0
      if(DoStencilFix)then
         XyzStencil_D((/iAxis,jAxis/)) =               &
              XyzGrid_DI((/iAxis,jAxis/),iOrder_I(2)) + &
              2*(XyzGrid_DI((/iAxis,jAxis/),iOrder_I(1))-XyzGrid_DI((/iAxis,jAxis/),iOrder_I(2)))
         !XyzStencil_D(kAxis) = 0.50*(XyzGrid_DI(kAxis,iOrder_I(2))+XyzGrid_DI(kAxis,iOrder_I(1)))
         XyzStencil_D(kAxis) = Xyz_D(kAxis)
         return
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

    end subroutine interpolate_corner_transition_junction
    !=======================================================================
    subroutine interpolate_corner_transition(iEdgeDir)
      integer, intent(in) :: iEdgeDir
      integer             :: iAxis,jAxis,kAxis
      real :: AlphaKAxisCoarse, AlphaKAxisFine
      real :: dXyzUp, dXyzDown ! Auxilary variables to determine vertical weight
      integer :: nFine, iFine_I(3)
      integer :: nCoarse, iCoarse_I(3)
      logical :: IsAmbiguous
      !\
      ! Subroutine intepolates in transitional near a resolution corner
      ! Points facing a corner a numbered 5:8
      !
      !------------------------------------------------------------------------
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
      call interpolate_amr_grid_2(&
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
         !write(*,*)'bowfix'
         DoStencilFix = DoStencilFix2
         !\
         ! Move to another transition region
         !/
         XyzStencil_D((/iAxis,jAxis/)) = XyzStencil2_D
         XyzStencil_D(kAxis) = Xyz_D(kAxis)
         return
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
            AlphaKAxisCoarse = (Xyz_D(     kAxis                  ) - XyzGrid_DI(kAxis,iOrder_I(nFine+5)))/&
                 (XyzGrid_DI(kAxis,iOrder_I(nFine+1)) - XyzGrid_DI(kAxis,iOrder_I(nFine+5)))
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
            AlphaKAxisCoarse = (Xyz_D(     kAxis                  ) - XyzGrid_DI(kAxis,iOrder_I(nFine+5)))/&
                 (XyzGrid_DI(kAxis,iOrder_I(nFine+1)) - XyzGrid_DI(kAxis,iOrder_I(nFine+5)))
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
         dXyzDown = Xyz_D(kAxis) - sum(XyzGrid_DI(kAxis,iOrder_I(1:nGridOut2))*Weight_I(1:nGridOut2))
         dXyzUp   = Xyz_D(kAxis) - sum(XyzGrid_DI(kAxis,iOrder_I(5:4+nGridOut2))*Weight_I(1:nGridOut2))
      else
         if(nFine==2)then
            dXyzDown = Xyz_D(kAxis) - XyzGrid_DI(kAxis,iOrder_I(1))
            dXyzUp   = Xyz_D(kAxis) - XyzGrid_DI(kAxis,iOrder_I(5))
         else
            dXyzDown = Xyz_D(kAxis) - sum(XyzGrid_DI(kAxis,iOrder_I(1:2))*Weight_I(1:2))-&
                 Weight_I(3)*(4*XyzGrid_DI(kAxis,iOrder_I(3))-XyzGrid_DI(kAxis,iOrder_I(7)))/3
            dXyzUp  = Xyz_D(kAxis) - sum(XyzGrid_DI(kAxis,iOrder_I((/5,2/)))*Weight_I(1:2))-&
                 Weight_I(3)*(2*XyzGrid_DI(kAxis,iOrder_I(3))+XyzGrid_DI(kAxis,iOrder_I(7)))/3
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
         return
      end if

      if(dXyzUp==0)then
         AlphaKAxisFine = 0.0
      else
         AlphaKAxisFine = abs(dXyzUp) / (abs(dXyzUp)+abs(dXyzDown))
      end if

      !\
      ! Remove points
      ! Apply vertical weight
      !/
      if(AlphaKAxisCoarse == 1)then
         nGridOut = nGridOut2 + nFine
         iOrder_I(nGridOut2+1:nGridOut) = iOrder_I(5:4+nFine)
         Weight_I(nGridOut2+1:nGridOut) = Weight_I(1:  nFine) * (1 - AlphaKAxisFine)
         Weight_I(          1:nFine   ) = Weight_I(1:  nFine) *      AlphaKAxisFine
      else
         nGridOut = 2*nGridOut2
         iOrder_I(nGridOut2      +1:nGridOut       ) = iOrder_I(      5:4+nGridOut2)
         Weight_I(nGridOut2      +1:nGridOut2+nFine) = Weight_I(      1:nFine      ) * (1 - AlphaKAxisFine)
         Weight_I(nGridOut2+nFine+1:nGridOut       ) = Weight_I(nFine+1:nGridOut2  ) * (1 - AlphaKAxisCoarse)
         Weight_I(                1:nFine          ) = Weight_I(      1:nFine      ) *      AlphaKAxisFine
         Weight_I(          nFine+1:nGridOut2      ) = Weight_I(nFine+1:nGridOut2  ) *      AlphaKAxisCoarse
      end if


    end subroutine interpolate_corner_transition







  end subroutine interpolate_amr_grid_3

  !===========================TESTS======================================================
  subroutine test_interpolate_amr_grid_2(nSample)
    !The test subtoutine nSample times generate random point in the AMR domain
    !and compares the result of the coordinate interpolation to this point with the
    !gampled coordinates
    !
    integer, intent(in) :: nSample

    integer :: iSeed=0, i, j, iBlock, loc(3), iSample, iGrid, iDir
    integer, parameter :: nDim = 2, nBlock = 44, nX=2, nY=2, nIndexes=3, nGrid = 4
    real    :: Xyz_DCB(nDim, nX, nY, nBlock)  !grid coordinates
    ! Block iOrder
    !
    ! C F C F C
    ! C F F F C
    ! C C F C C
    ! C F C F C
    !
    ! Block numbers
    !
    ! 34 35/38   39  40/43  44
    ! 20 21/24 25/28 29/32  33
    ! 12  13   14/17  18    19
    ! 1   2/5    6   7/10   11
    ! Block corner coordinates:
    real, parameter :: Xyz_DB(nDim, nBlock) = reshape((/&
         0.0, 0.0, 2.0, 0.0,  3.0, 0.0, 2.0, 1.0, 3.0, 1.0, 4.0, 0.0, 6.0, 0.0, 7.0, 0.0,  6.0, 1.0, 7.0, 1.0, 8.0, 0.0, &
         0.0, 2.0, 2.0, 2.0,  4.0, 2.0, 5.0, 2.0, 4.0, 3.0, 5.0, 3.0, 6.0, 2.0, 8.0, 2.0,  &
         0.0, 4.0, 2.0, 4.0,  3.0, 4.0, 2.0, 5.0, 3.0, 5.0, 4.0, 4.0, 5.0, 4.0, 4.0, 5.0,  5.0, 5.0, 6.0, 4.0, 7.0, 4.0, &
         6.0, 5.0, 7.0, 5.0,  8.0, 4.0,        &
         0.0, 6.0, 2.0, 6.0,  3.0, 6.0, 2.0, 7.0, 3.0, 7.0, 4.0, 6.0, 6.0, 6.0, 7.0, 6.0,  6.0, 7.0, 7.0, 7.0, 8.0, 6.0  &
         /),(/nDim, nBlock/))
    real :: dXyz_DB(nDim, nBlock) = 0.0
    integer, parameter :: iLevel_B(nBlock) = (/&
         4,5,5,5,5,4,5,5,5,4,4,4,4,5,5,5,5,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,4,4,5,5,5,5,4,5,5,5,5,4/)
    !\
    ! Parameters to call interpolate_amr_grid_2
    !/
    real :: Xyz_D(nDim)=0.0, Xyz2_D(nDim)=0.0, XyzStencil_D(nDim)=0.0, XyzGrid_DI(nDim,4)=0.0, XyzGrid2_DI(nDim,4)=0.0
    real :: Weight_I(4), Weight2_I(4)
    integer :: iOrder_I(4),iOrder2_I(4)
    real :: XyzInterpolated_D(nDim)                             ! for test of order of approximation
    real :: XyzSquaredInterpolated1_D(nDim),XyzSquaredInterpolated2_D(nDim) ! for test of continuity
    integer:: iLevel_I(4), iIndexes_II(nIndexes,4), nGridOut,iLevel_It(4), iIndexes2_II(nIndexes,4),nGridOut2
    logical :: DoStencilFix
    !---------------------
    !\
    !initialization
    !/
    do iBlock = 1, nBlock
       dXyz_DB(:,iBlock) = (/1.0,1.0/)/(iLevel_B(iBlock)-3)
    end do

    do iBlock = 1, nBlock; do j=1,nY; do i = 1,nX
       Xyz_DCB(1,i,j,iBlock) = Xyz_DB(1,iBlock)+(i-0.50)*dXyz_DB(1,iBlock)
       Xyz_DCB(2,i,j,iBlock) = Xyz_DB(2,iBlock)+(j-0.50)*dXyz_DB(2,iBlock)
    end do; end do; end do

    call init_rand()
    !\
    ! Test order of approximation
    !/
    do iSample = 1, nSample
       Xyz_D(1) = 10.0*RAND()
       Xyz_D(2) =  8.0*RAND()
       Xyz2_D(1) = Xyz_D(1) + 0.0001
       Xyz2_D(2) = Xyz_D(2) + 0.0001
       if(.not.any(Xyz_DCB(1,:,:,:).le.Xyz_D(1).and.&
            Xyz_DCB(2,:,:,:).le.Xyz_D(2)))CYCLE
       if(.not.any(Xyz_DCB(1,:,:,:) > Xyz_D(1).and.&
            Xyz_DCB(2,:,:,:).le.Xyz_D(2)))CYCLE
       if(.not.any(Xyz_DCB(1,:,:,:).le.Xyz_D(1).and.&
            Xyz_DCB(2,:,:,:)  > Xyz_D(2)))CYCLE
       if(.not.any(Xyz_DCB(1,:,:,:) > Xyz_D(1).and.&
            Xyz_DCB(2,:,:,:)  > Xyz_D(2)))CYCLE
       DoStencilFix = .false.
       call generate_basic_stencil(Xyz_D)
       if(any(XyzGrid_DI(:,:)<0))CYCLE
       call interpolate_amr_grid_2(Xyz_D, XyzGrid_DI, iLevel_I,&
            nGridOut, Weight_I, iOrder_I, DoStencilFix, XyzStencil_D)
       if(DoStencilFix)then
          if(abs(XyzStencil_D(x_)-0) .le.0.1)CYCLE
          if(abs(XyzStencil_D(x_)-10).le.0.1)CYCLE
          if(abs(XyzStencil_D(y_)-0) .le.0.1)CYCLE
          if(abs(XyzStencil_D(y_)-8) .le.0.1)CYCLE
          call generate_basic_stencil(XyzStencil_D)
          if(any(XyzGrid_DI(:,:)<0))CYCLE
          call interpolate_amr_grid_2(Xyz_D, XyzGrid_DI, iLevel_I,&
               nGridOut, Weight_I, iOrder_I, DoStencilFix, XyzStencil_D)
       end if

       XyzInterpolated_D=0.0
       XyzSquaredInterpolated1_D=0.0
       XyzSquaredInterpolated2_D=0.0

       do iGrid = 1, nGridOut
          i = iIndexes_II(     1, iOrder_I(iGrid))
          j = iIndexes_II(     2, iOrder_I(iGrid))
          iBlock = iIndexes_II(3, iOrder_I(iGrid))
          XyzInterpolated_D = XyzInterpolated_D + &
               Weight_I(iGrid) * Xyz_DCB(:, i, j, iBlock)
          XyzSquaredInterpolated1_D = XyzSquaredInterpolated1_D + &
               Weight_I(iGrid) * Xyz_DCB(:, i, j, iBlock)* Xyz_DCB(:, i, j, iBlock)
       end do

       !\
       ! Order of approximation test
       ! Check interpolated value and compare with Xyz_D
       !/
       if(sum(abs(Xyz_D(:)-XyzInterpolated_D(:))) > 1.0e-5)then
          write(*,*)'Test failed at iSample = ', iSample
          write(*,*)'Sampled Xyz_D=', Xyz_D
          write(*,*)'Interpolated value=',XyzInterpolated_D
          write(*,*)'Stencil: i j iBlock Weight'
          do iGrid = 1, nGridOut
             write(*,*)'Point #',iGrid, iIndexes_II(:,iOrder_I(iGrid)), Weight_I(iOrder_I(iGrid))
          end do
          write(*,*)'Step by step generate stencil'
          DoStencilFix = .false.
          call generate_basic_stencil(Xyz_D)
          write(*,*)'Stencil: i j iBlock Level XyzGrid_D'
          do iGrid = 1, 4
             write(*,*)'Point #',iGrid, iIndexes_II(:,iGrid), iLevel_I(iGrid), XyzGrid_DI(:,iGrid)
          end do
          call interpolate_amr_grid_2(Xyz_D, XyzGrid_DI, iLevel_I,&
               nGridOut, Weight_I, iOrder_I, DoStencilFix, XyzStencil_D)
          write(*,*)'DoStencilFix=',DoStencilFix, ' nGridOut=', nGridOut
          call CON_stop('Test failed')
       end if

       !------------------------------------------------------------------------------------------------
       DoStencilFix = .false.
       call generate_basic_stencil(Xyz2_D)
       if(any(XyzGrid_DI(:,:)<0))CYCLE
       call interpolate_amr_grid_2(Xyz2_D, XyzGrid_DI, iLevel_I,&
            nGridOut2, Weight2_I, iOrder2_I, DoStencilFix, XyzStencil_D)
       if(DoStencilFix)then
          if(abs(XyzStencil_D(x_)-0) .le.0.1)CYCLE
          if(abs(XyzStencil_D(x_)-10).le.0.1)CYCLE
          if(abs(XyzStencil_D(y_)-0) .le.0.1)CYCLE
          if(abs(XyzStencil_D(y_)-8) .le.0.1)CYCLE
          call generate_basic_stencil(XyzStencil_D)
          if(any(XyzGrid_DI(:,:)<0))CYCLE
          call interpolate_amr_grid_2(Xyz2_D, XyzGrid_DI, iLevel_I,&
               nGridOut2, Weight2_I, iOrder2_I, DoStencilFix, XyzStencil_D)
       end if

       do iGrid = 1, nGridOut2
          i = iIndexes_II(1,iOrder2_I(iGrid))
          j = iIndexes_II(2,iOrder2_I(iGrid))
          iBlock = iIndexes_II(3,iOrder2_I(iGrid))
          XyzSquaredInterpolated2_D = XyzSquaredInterpolated2_D + &
               Weight2_I(iGrid) * Xyz_DCB(:, i, j, iBlock)* Xyz_DCB(:, i, j, iBlock)
       end do
       !\
       ! Continuity test
       ! Check interpolated values of function f(x,y)=x^2+y^2 in points Xyz_D and Xyz2_D
       !/
       if(sum(abs(XyzSquaredInterpolated1_D(:)-XyzSquaredInterpolated2_D(:)-Xyz_D(:)*Xyz_D(:)+Xyz2_D(:)*Xyz2_D(:)))&
            > 5.0e-4)then
          write(*,*)'Continuity test failed at iSample = ', iSample
          write(*,*)'Sampled Xyz2_D=',  Xyz2_D
          write(*,*)'Sampled Xyz_D=', Xyz_D
          write(*,*)'Interpolated value2=',XyzSquaredInterpolated2_D
          write(*,*)'Interpolated value1=',XyzSquaredInterpolated1_D
          write(*,*)'Stencil: i j iBlock Weight'
          do iGrid = 1, nGridOut2
             write(*,*)'Point #',iGrid, iIndexes_II(:,iOrder2_I(iGrid)), Weight2_I(iOrder2_I(iGrid))
          end do
          write(*,*)'+++++++++'
          call generate_basic_stencil(Xyz_D)
          call interpolate_amr_grid_2(Xyz_D, XyzGrid_DI, iLevel_I,&
               nGridOut, Weight_I,iOrder_I, DoStencilFix, XyzStencil_D)
          if(DoStencilFix)then
             call generate_basic_stencil(XyzStencil_D)
             call interpolate_amr_grid_2(Xyz_D, XyzGrid_DI, iLevel_I,&
                  nGridOut, Weight_I,iOrder_I, DoStencilFix, XyzStencil_D)
          end if
          do iGrid = 1, nGridOut
             write(*,*)'Point #',iGrid, iIndexes_II(:,iOrder_I(iGrid)), Weight_I(iOrder_I(iGrid))
          end do
          write(*,*)'Step by step generate stencil'
          DoStencilFix = .false.
          call generate_basic_stencil(Xyz2_D)
          write(*,*)'Stencil: i j iBlock Level XyzGrid_D'
          do iGrid = 1, 4
             write(*,*)'Point #',iGrid, iIndexes_II(:,iGrid), iLevel_I(iGrid), XyzGrid_DI(:,iGrid)
          end do
          call interpolate_amr_grid_2(Xyz_D, XyzGrid_DI, iLevel_I,&
               nGridOut, Weight_I, iOrder_I, DoStencilFix, XyzStencil_D)
          write(*,*)'DoStencilFix=',DoStencilFix, ' nGridOut=', nGridOut
          call CON_stop('Test failed')
       end if

    end do

  contains
    subroutine generate_basic_stencil(XyzIn_D)
      real, intent(in):: XyzIn_D(nDim)
      integer:: nBlockNei
      integer, dimension(4):: iBlockNei_I
      real, dimension(:,:,:,:), allocatable:: XyzNei_DCB
      !-------------------------------
      !\
      ! reset XyzGrid_DI
      !/
      do iGrid = 1, nGrid
         do iDir =1, nDim
            XyzGrid_DI(iDir,iGrid) = -1.0
         end do
      end do
      nBlockNei=1
      do iBlock = 1, nBlock
         if(( (Xyz_DB(1,iBlock)+dXyz_DB(1,iBlock)-XyzIn_D(1)).le. ( 0.5+dXyz_DB(1,iBlock) )).and.&
              ( (Xyz_DB(1,iBlock)+dXyz_DB(1,iBlock)-XyzIn_D(1))  > -( 0.5+dXyz_DB(1,iBlock) )).and.&
              ( (Xyz_DB(2,iBlock)+dXyz_DB(2,iBlock)-XyzIn_D(2)).le. ( 0.5+dXyz_DB(2,iBlock) )).and.&
              ( (Xyz_DB(2,iBlock)+dXyz_DB(2,iBlock)-XyzIn_D(2))  > -( 0.5+dXyz_DB(2,iBlock) )))then
            iBlockNei_I(nBlockNei) = iBlock
            nBlockNei = nBlockNei + 1
         end if
      end do
      nBlockNei = nBlockNei - 1

      allocate(XyzNei_DCB(nDim,nX,nY,nBlockNei))

      XyzNei_DCB = Xyz_DCB(:,:,:,iBlockNei_I(1:nBlockNei))

      !-------------------------------
      !gird point 1
      Loc = minloc((XyzIn_D(2)-XyzNei_DCB(2,:,:,:))**2+&
           (XyzIn_D(1)-XyzNei_DCB(1,:,:,:))**2, &
           MASK=XyzNei_DCB(1,:,:,:).le.XyzIn_D(1).and.&
           XyzNei_DCB(2,:,:,:).le.XyzIn_D(2))
      if(.not. (any(Loc==0)))then
         iIndexes_II(:,1) = loc(:)
         iIndexes_II(3,1) = iBlockNei_I(loc(3))
         iLevel_I(   1) = iLevel_B(iBlockNei_I(loc(3)))
         XyzGrid_DI(:,1)= XyzNei_DCB(:,loc(1),loc(2),loc(3))
      end if
      !gird point 2
      Loc = minloc((XyzIn_D(2)-XyzNei_DCB(2,:,:,:))**2+&
           (XyzIn_D(1)-XyzNei_DCB(1,:,:,:))**2, &
           MASK=XyzNei_DCB(1,:,:,:) > XyzIn_D(1).and.&
           XyzNei_DCB(2,:,:,:).le.XyzIn_D(2))
      if(.not. (any(Loc==0)))then
         iIndexes_II(:,2) = loc(:)
         iIndexes_II(3,2) = iBlockNei_I(loc(3))
         iLevel_I(   2) = iLevel_B(iBlockNei_I(loc(3)))
         XyzGrid_DI(:,2)= XyzNei_DCB(:,loc(1),loc(2),loc(3))
      end if
      !gird point 3
      Loc = minloc((XyzIn_D(2)-XyzNei_DCB(2,:,:,:))**2+&
           (XyzIn_D(1)-XyzNei_DCB(1,:,:,:))**2, &
           MASK=XyzNei_DCB(1,:,:,:).le.XyzIn_D(1).and.&
           XyzNei_DCB(2,:,:,:)>XyzIn_D(2))
      if(.not. (any(Loc==0)))then
         iIndexes_II(:,3) = loc(:)
         iIndexes_II(3,3) = iBlockNei_I(loc(3))
         iLevel_I(   3) = iLevel_B(iBlockNei_I(loc(3)))
         XyzGrid_DI(:,3)= XyzNei_DCB(:,loc(1),loc(2),loc(3))
      end if
      !gird point 4
      Loc = minloc((XyzIn_D(2)-XyzNei_DCB(2,:,:,:))**2+&
           (XyzIn_D(1)-XyzNei_DCB(1,:,:,:))**2, &
           MASK=XyzNei_DCB(1,:,:,:) > XyzIn_D(1).and.&
           XyzNei_DCB(2,:,:,:) > XyzIn_D(2))
      if(.not. (any(Loc==0)))then
         iIndexes_II(:,4) = loc(:)
         iIndexes_II(3,4) = iBlockNei_I(loc(3))
         iLevel_I(   4) = iLevel_B(iBlockNei_I(loc(3)))
         XyzGrid_DI(:,4)= XyzNei_DCB(:,loc(1),loc(2),loc(3))
      end if

      deallocate(XyzNei_DCB)
    end subroutine generate_basic_stencil

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
    REAL FUNCTION RAND()
      ISEED=ISEED*48828125
      IF(ISEED < 0) ISEED=(ISEED+2147483647)+1
      if(iSeed==0) iSeed=1
      RAND=FLOAT(ISEED)/2147483647
    END FUNCTION RAND


  end subroutine test_interpolate_amr_grid_2
  !=========================================


end module ModInterpolateAMRGrid
