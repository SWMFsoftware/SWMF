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

  public :: interpolate_amr_grid_3

  integer,parameter:: BehindTheBoundary_ = -7777, &
                      Coarse_  = 0,               &
                      Fine_    = 1,               &
                      x_       = 1,               &
                      y_       = 2,               &
                      z_       = 3

  public :: BehindTheBoundary_
  
contains
  !==================================================================================================
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
            call CON_stop('The point is out of the grid '//NameSub)
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
    !If nGridOut < nGridIn, only the first nGridOut lines in  XyzGrid_DI, iIndexes_II
    !are meaningful in the output
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
    ! used to construct stencil, not to interpolate!!!
    !/
    real, intent(out) :: XyzStencil_D(nDim)


    integer :: iGrid, iDim   !Loop variables
    integer :: Loc(1)        !Find minloc and maxloc
    integer :: nGridOut1     !Number of grid points returned from 1D interpolation routine
    integer :: iLevelMin     !Minimumum refinement level

    real    :: XyzMin_D(nDim), XyzMax_D(nDim)    ! Minimum and maximum values of coordinates
    real    :: dXyzSmall_D(nDim)                 ! To find lines with the "constant" cordinate
    real    :: dXyzInv_D(  nDim)                 ! for calculating on a uniform grid
    real    :: Aux_D(nDim)                       ! Misc
    real    :: Xyz1_D(1), XyzGrid1_DI(1,2)       ! To call one-dimensional interpolation
    integer :: iOrder1_I(2)          ! To call one-dimensional interpolation
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
    do iDim = 1, nDim
        dXyzSmall_D(iDim) = 0.10*(maxval(XyzGrid_DI(iDim,:)) -  minval(XyzGrid_DI(iDim,:)))
        dXyzInv_D(iDim)   = 0.10/dXyzSmall_D(iDim)
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
    real function cross_product_2(a_D, b_D)
        real,dimension(nDim),intent(in) :: a_D, b_D
        !----------
        cross_product_2 = a_D(x_)* b_D(y_) - a_d(y_)*b_D(x_)
    end function cross_product_2
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
        Alpha3 = cross_product_2(Xyz_D - X1_D, X2_D-X1_D)/cross_product_2(X3_D - X1_D, X2_D-X1_D)
        if(Alpha3==0.0)then
            nGridOut = 2
            Alpha2 = cross_product_2(X3_D-X1_D,Xyz_D - X1_D)/cross_product_2(X3_D - X1_D, X2_D-X1_D)
            Weight_I(2) = Alpha2
        elseif(Alpha3 > 0.0)then
            nGridOut = 3
            Alpha2 = cross_product_2(X3_D-X1_D,Xyz_D - X1_D)/cross_product_2(X3_D - X1_D, X2_D-X1_D)
            Weight_I(2) = Alpha2
            Weight_I(3) = Alpha3
        else
            nGridOut = 3
            Alpha2 = cross_product_2(X4_D - X1_D,  Xyz_D - X1_D)/cross_product_2(X4_D - X1_D, X2_D-X1_D)
            Alpha3 = cross_product_2(Xyz_D - X1_D, X2_D  - X1_D)/cross_product_2(X4_D - X1_D, X2_D-X1_D)
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

        Weight_I(3) = cross_product_2(Xyz_D - X1_D, X2_D  - X1_D)/cross_product_2(X3_D - X1_D, X2_D-X1_D)
        Weight_I(2) = cross_product_2(X3_D  - X1_D, Xyz_D - X1_D)/cross_product_2(X3_D - X1_D, X2_D-X1_D)
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
    ! used to construct stencil, not to interpolate!!!
    !/
    real, intent(out) :: XyzStencil_D(nDim)

    integer :: iGrid , nGridOut2, iLevelMin, iDim
    integer :: Loc(1)
    real    :: dXyzSmall_D(nDim), dXyzInv_D(nDim)
    real    :: Aux_D(nDim)
    real    :: Xyz2_D(x_:y_), XyzGrid2_DI(x_:y_,4), XyzStencil2_D(x_:y_)
    integer :: iOrder2_I(4), iLevel2_I(4)
    logical :: DoStencilFix2
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

    if(all(iLevel_I(5:8)== BehindTheBoundary_).or. all(XyzGrid_DI(z_,1:4)==Xyz_D(z_)) )then

        !\
        ! Face 5:8 should be removed
        !/
        iOrder_I = (/1,2,3,4,5,6,7,8/)
        Xyz2_D(x_:y_) = Xyz_D(x_:y_)
        DoStencilFix2 = .false.
        XyzGrid2_DI(x_:y_,1:4) = XyzGrid_DI((/x_,y_/),iOrder_I(1:4))
        iLevel2_I(        1:4) = iLevel_I(            iOrder_I(1:4))
        call interpolate_amr_grid_2(&
                Xyz2_D,  XyzGrid2_DI, iLevel2_I, &
                nGridOut, Weight_I(1:4), iOrder2_I, DoStencilFix2, XyzStencil2_D)

       if(DoStencilFix2)then
           DoStencilFix = .true.
           XyzStencil_D(x_:y_) = XyzStencil2_D(x_:y_)
           XyzStencil_D(z_)    = Xyz_D(z_)
       end if

        iOrder_I(1:4) = iOrder_I(iOrder2_I)
       return

    elseif(all(iLevel_I(1:4)== BehindTheBoundary_)                                    )then

        !\
        ! Face 1:4 should be removed
        !/
        iOrder_I = (/5,6,7,8,1,2,3,4/)
        Xyz2_D(x_:y_) = Xyz_D(x_:y_)
        DoStencilFix2 = .false.
        XyzGrid2_DI(x_:y_,1:4) = XyzGrid_DI((/x_,y_/),iOrder_I(1:4))
        iLevel2_I(        1:4) = iLevel_I(            iOrder_I(1:4))
        call interpolate_amr_grid_2(&
            Xyz2_D,  XyzGrid2_DI, iLevel2_I, &
            nGridOut, Weight_I(1:4), iOrder2_I, DoStencilFix2, XyzStencil2_D)

       if(DoStencilFix2)then
           DoStencilFix = .true.
           XyzStencil_D(x_:y_) = XyzStencil2_D(x_:y_)
           XyzStencil_D(z_)    = Xyz_D(z_)
       end if

       iOrder_I(1:4) = iOrder_I(iOrder2_I)
       return

    elseif(all(iLevel_I(2:8:2)== BehindTheBoundary_).or. all(XyzGrid_DI(x_,1:7:2)==Xyz_D(x_)) )then

        !\
        ! Face 2,4,6,8 should be removed
        !/
        iOrder_I = (/1,3,5,7,2,4,6,8/)
        Xyz2_D(x_:y_) = Xyz_D(y_:z_)
        DoStencilFix2 = .false.
        XyzGrid2_DI(x_:y_,1:4) = XyzGrid_DI((/y_,z_/),iOrder_I(1:4))
        iLevel2_I(        1:4) = iLevel_I(            iOrder_I(1:4))
        call interpolate_amr_grid_2(&
            Xyz2_D,  XyzGrid2_DI, iLevel2_I, &
            nGridOut, Weight_I(1:4), iOrder2_I, DoStencilFix2, XyzStencil2_D)

       if(DoStencilFix2)then
           DoStencilFix = .true.
           XyzStencil_D(y_:z_) = XyzStencil2_D(x_:y_)
           XyzStencil_D(x_)    = Xyz_D(x_)
       end if

       iOrder_I(1:4) = iOrder_I(iOrder2_I)
       return

    elseif(all(iLevel_I(1:7:2)== BehindTheBoundary_)                                        )then

        !\
        ! Face 1,3,5,7 should be removed
        !/
        iOrder_I = (/2,4,6,8,1,3,5,7/)
        Xyz2_D(x_:y_) = Xyz_D(y_:z_)
        DoStencilFix2 = .false.
        XyzGrid2_DI(x_:y_,1:4) = XyzGrid_DI((/y_,z_/),iOrder_I(1:4))
        iLevel2_I(        1:4) = iLevel_I(                  iOrder_I(1:4))
        call interpolate_amr_grid_2(&
            Xyz2_D,  XyzGrid2_DI, iLevel2_I, &
            nGridOut, Weight_I(1:4), iOrder2_I, DoStencilFix2, XyzStencil2_D)

       if(DoStencilFix2)then
           DoStencilFix = .true.
           XyzStencil_D(y_:z_) = XyzStencil2_D(x_:y_)
           XyzStencil_D(x_)    = Xyz_D(x_)
       end if

       iOrder_I(1:4) = iOrder_I(iOrder2_I)
       return

    elseif(all(iLevel_I((/3,4,7,8/))== BehindTheBoundary_).or.&
           all(XyzGrid_DI(y_,(/1,2,5,6/))==Xyz_D(y_)))then

        !\
        ! Face 3,4,7,8 should be removed
        !/
        iOrder_I = (/1,2,5,6,3,4,7,8/)
        Xyz2_D(x_:y_) = Xyz_D((/x_,z_/))
        DoStencilFix2 = .false.
        XyzGrid2_DI(x_:y_,1:4) = XyzGrid_DI((/x_,z_/),iOrder_I(1:4))
        iLevel2_I(        1:4) = iLevel_I(            iOrder_I(1:4))
        call interpolate_amr_grid_2(&
            Xyz2_D,  XyzGrid2_DI, iLevel2_I, &
            nGridOut, Weight_I(1:4), iOrder2_I, DoStencilFix, XyzStencil2_D)

       if(DoStencilFix2)then
           DoStencilFix = .true.
           XyzStencil_D( (/x_,z_/) ) = XyzStencil2_D(x_:y_)
           XyzStencil_D(y_)          = Xyz_D(y_)
       end if

       iOrder_I(1:4) = iOrder_I(iOrder2_I)
       return

    elseif( all(iLevel_I((/1,2,5,6/))== BehindTheBoundary_))then

        !\
        ! Face 1,2,5,6 should be removed
        !/
       iOrder_I = (/3,4,7,8,1,2,5,6/)
        Xyz2_D(x_:y_) = Xyz_D((/x_,z_/))
        DoStencilFix2 = .false.
        XyzGrid2_DI(x_:y_,1:4) = XyzGrid_DI((/x_,z_/),iOrder_I(1:4))
        iLevel2_I(        1:4) = iLevel_I(            iOrder_I(1:4))
        call interpolate_amr_grid_2(&
            Xyz2_D,  XyzGrid2_DI, iLevel_I(1:4), &
            nGridOut, Weight_I(1:4), iOrder2_I, DoStencilFix, XyzStencil2_D)

       if(DoStencilFix2)then
           DoStencilFix = .true.
           XyzStencil_D( (/x_,z_/) ) = XyzStencil2_D(x_:y_)
           XyzStencil_D(y_)          = Xyz_D(y_)
        end if

        iOrder_I(1:4) = iOrder_I(iOrder2_I)
        return


    end if




    if(any(iLevel_I==BehindTheBoundary_) ) then
        !\
        ! There still can be three grid points behind the grid boundary - all on the same face
        !/
        !             Upper face          Lower face
        !             Domain|Boundary     Domain|Boundary
        !               * * | * *           * * | * *           
        !               * C | P *           * C | P *        + symmetry transformations
        !               * C | P *           F F | * *
        !               * * | * *           * * | * *
        !
        ! In this configuration there are three grid points behind the boundary, two fine and
        !  three coarse grid points.
        !\
        ! There still can be two grid points behind the grid boundary both on the same face
        !/
        !             Upper face          Lower face
        !             Domain|Boundary     Domain|Boundary
        !               * * | * *           * * | * *
        !               * C | P *           F F | * *        + symmetry transformations
        !               * C | P *           F F | * *
        !               * * | * *           * * | * *
        !  OR
        !             Upper face          Lower face
        !             Domain|Boundary     Domain|Boundary
        !               * * | * *           * * | * *
        !               * C | P *           F F | * *        + symmetry transformations
        !               F F | * *           * C | P *
        !               * * | * *           * * | * *
        !
        ! In this configuration there are two grid points behind the boundary, four fine and
        !  two coarse grid points.
        !\
        ! There still can be one grid point behind the grid boundary
        !/
        !             Upper face          Lower face
        !             Domain|Boundary     Domain|Boundary
        !               * * | * *           * * | * *
        !               * C | P *           F F | * *        + symmetry transformations
        !               F F | * *           F F | * *
        !               * * | * *           * * | * *
        !
        !\
        ! In this configuration there is one grid point behind the boundary, six fine and
        !  one coarse grid point.
        !/
        !\
        ! Otherwise BehindTheBoundary_ points are illegal.
        !/

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
    do iDim = 1, nDim
    dXyzSmall_D(iDim) = 0.10*(maxval(XyzGrid_DI(iDim,:)) -  minval(XyzGrid_DI(iDim,:)))
    dXyzInv_D(iDim)   = 0.10/dXyzSmall_D(iDim)
    end do

    !\
    ! First conbinations of the resolution interfaces: no resolution interface or
    ! a single resolition line
    !/

    if(all(iLevel_I==0))then
        !\
        ! No refinement         | C C |  | C C |
        !                       | C C |  | C C |
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
        return

       !\
       ! Opposite faces going along resolution interfaces
       !/

elseif(all(abs(XyzGrid_DI(z_,1) - XyzGrid_DI(z_,2:4)) < dXyzSmall_D(z_))&
           .and.&
       all(abs(XyzGrid_DI(z_,5) - XyzGrid_DI(z_,6:8)) < dXyzSmall_D(z_)))then

       !\
       ! Faces going along xy-plane    |C C|   |F F|
       !                               |C C|   |F F| 
       !/
        iOrder_I = (/1,2,3,4,5,6,7,8/)

        Xyz2_D(x_:y_) = Xyz_D(x_:y_)
        Aux_D(z_) = dXyzInv_D(z_)*(Xyz_D(z_)-XyzGrid_DI(z_,iOrder_I(1)))
        !\
        ! Interpolate along lower face
        !/
        DoStencilFix2 = .false.
        XyzGrid2_DI(x_:y_,1:4) = XyzGrid_DI(x_:y_,iOrder_I(1:4))
        iLevel2_I(        1:4) = iLevel_I(        iOrder_I(1:4))
        call interpolate_amr_grid_2(&
            Xyz2_D, XyzGrid2_DI, iLevel2_I, &
            nGridOut2, Weight_I(1:4), iOrder2_I, DoStencilFix2, XyzStencil2_D)
            iOrder_I(1:4) = iOrder_I(iOrder2_I)
        !\
        ! Apply weight for interpolation along z-axis
        !/
        Weight_I(1:nGridOut2) = Weight_I(1:nGridOut2) * (1 - Aux_D(z_))

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

        XyzGrid2_DI(x_:y_,1:4) = XyzGrid_DI(x_:y_,iOrder_I(nGridOut2+1:nGridOut2+4))
        iLevel2_I(        1:4) = iLevel_I(        iOrder_I(nGridOut2+1:nGridOut2+4))
        call interpolate_amr_grid_2(&
            Xyz2_D, XyzGrid2_DI, iLevel2_I, &
            nGridOut2, Weight_I(nGridOut+1:nGridOut+4), iOrder2_I, DoStencilFix2, XyzStencil2_D)
        do iGrid=1,nGridOut2; iOrder_I(nGridOut+iGrid)=iOrder_I(nGridOut+iOrder2_I(iGrid)); end do

        !\
        ! Apply weight for interolation along kAxis
        !/
        Weight_I(nGridOut+1:nGridOut+nGridOut2) = Weight_I(nGridOut+1:nGridOut+nGridOut2) * Aux_D(z_)

        !\
        ! May need to remove a grid points from the stencil once behind the boundary
        !/
        nGridOut = nGridOut + nGridOut2
       return

    elseif(abs(XyzGrid_DI(x_,1) - XyzGrid_DI(x_,3)) < dXyzSmall_D(x_).and.&
           abs(XyzGrid_DI(x_,3) - XyzGrid_DI(x_,5)) < dXyzSmall_D(x_).and.&
           abs(XyzGrid_DI(x_,5) - XyzGrid_DI(x_,7)) < dXyzSmall_D(x_)&
           .and.&
           abs(XyzGrid_DI(x_,2) - XyzGrid_DI(x_,4)) < dXyzSmall_D(x_).and.&
           abs(XyzGrid_DI(x_,4) - XyzGrid_DI(x_,6)) < dXyzSmall_D(x_).and.&
           abs(XyzGrid_DI(x_,6) - XyzGrid_DI(x_,8)) < dXyzSmall_D(x_))then

       !\
       ! Faces going along yz-plane    | FF |   | FF |  or   |C  C|   |C  C|
       !                               |C  C|   |C  C|       | FF |   | FF |
       !/
        iOrder_I = (/1,3,5,7,2,4,6,8/)
        Xyz2_D(x_:y_) = Xyz_D(y_:z_)
        Aux_D(x_) = dXyzInv_D(x_)*(Xyz_D(x_)-XyzGrid_DI(x_,iOrder_I(1)))
        !\
        ! Interpolate along lower face
        !/
        DoStencilFix2 = .false.
        XyzGrid2_DI(x_:y_,1:4) = XyzGrid_DI(y_:z_,iOrder_I(1:4))
        iLevel2_I(        1:4) = iLevel_I(        iOrder_I(1:4))
        call interpolate_amr_grid_2(&
            Xyz2_D, XyzGrid2_DI, iLevel2_I, &
            nGridOut2, Weight_I(1:4), iOrder2_I, DoStencilFix2, XyzStencil2_D)
        iOrder_I(1:4) = iOrder_I(iOrder2_I)
        !\
        ! Apply weight for interpolation along x_
        !/
        Weight_I(1:nGridOut2) = Weight_I(1:nGridOut2) * (1 - Aux_D(x_))

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

        XyzGrid2_DI(x_:y_,1:4) = XyzGrid_DI(y_:z_,iOrder_I(nGridOut2+1:nGridOut2+4))
        iLevel2_I(        1:4) = iLevel_I(        iOrder_I(nGridOut2+1:nGridOut2+4))
        call interpolate_amr_grid_2(&
            Xyz2_D, XyzGrid2_DI, iLevel2_I, &
            nGridOut2, Weight_I(nGridOut+1:nGridOut+4), iOrder2_I, DoStencilFix2, XyzStencil2_D)
        do iGrid=1,nGridOut2; iOrder_I(nGridOut+iGrid)=iOrder_I(nGridOut+iOrder2_I(iGrid)); end do

        !\
        ! Apply weight for interolation along x_
        !/
        Weight_I(nGridOut+1:nGridOut+nGridOut2) = Weight_I(nGridOut+1:nGridOut+nGridOut2) * Aux_D(x_)

        !\
        ! May need to remove a grid points from the stencil once behind the boundary
        !/
        nGridOut = nGridOut + nGridOut2
       return

    elseif(abs(XyzGrid_DI(y_,1) - XyzGrid_DI(y_,2)) < dXyzSmall_D(y_).and.&
           abs(XyzGrid_DI(y_,2) - XyzGrid_DI(y_,5)) < dXyzSmall_D(y_).and.&
           abs(XyzGrid_DI(y_,5) - XyzGrid_DI(y_,6)) < dXyzSmall_D(y_)&
           .and.&
           abs(XyzGrid_DI(y_,3) - XyzGrid_DI(y_,4)) < dXyzSmall_D(y_).and.&
           abs(XyzGrid_DI(y_,4) - XyzGrid_DI(y_,7)) < dXyzSmall_D(y_).and.&
           abs(XyzGrid_DI(y_,7) - XyzGrid_DI(y_,8)) < dXyzSmall_D(y_))then

       !\
       ! Faces going along xz-plane    |C F|   |C F|  or   |F C|   |F C|
       !                               |C F|   |C F|       |F C|   |F C|
       !/
        iOrder_I = (/1,2,5,6,3,4,7,8/)
        Xyz2_D(x_:y_) = Xyz_D((/x_,z_/))
        Aux_D(y_) = dXyzInv_D(y_)*(Xyz_D(y_)-XyzGrid_DI(y_,iOrder_I(1)))
        !\
        ! Interpolate along lower face
        !/
        DoStencilFix2 = .false.
        XyzGrid2_DI(x_:y_,1:4) = XyzGrid_DI((/x_,z_/),iOrder_I(1:4))
        iLevel2_I(        1:4) = iLevel_I(            iOrder_I(1:4))
        call interpolate_amr_grid_2(&
            Xyz2_D, XyzGrid2_DI, iLevel2_I, &
            nGridOut2, Weight_I(1:4), iOrder2_I, DoStencilFix2, XyzStencil2_D)
        iOrder_I(1:4) = iOrder_I(iOrder2_I)
        !\
        ! Apply weight for interpolation along kAxis
        !/
        Weight_I(1:nGridOut2) = Weight_I(1:nGridOut2) * (1 - Aux_D(y_))

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

        XyzGrid2_DI(x_:y_,1:4) = XyzGrid_DI((/x_,z_/),iOrder_I(nGridOut2+1:nGridOut2+4))
        iLevel2_I(        1:4) = iLevel_I(            iOrder_I(nGridOut2+1:nGridOut2+4))
        call interpolate_amr_grid_2(&
            Xyz2_D, XyzGrid2_DI, iLevel2_I, &
            nGridOut2, Weight_I(nGridOut+1:nGridOut+4), iOrder2_I, DoStencilFix2, XyzStencil2_D)

        do iGrid=1,nGridOut2; iOrder_I(nGridOut+iGrid)=iOrder_I(nGridOut+iOrder2_I(iGrid)); end do

        !\
        ! Apply weight for interolation along y_
        !/
        Weight_I(nGridOut+1:nGridOut+nGridOut2) = Weight_I(nGridOut+1:nGridOut+nGridOut2) * Aux_D(y_)

        !\
        ! May need to remove a grid points from the stencil once behind the boundary
        !/
        nGridOut = nGridOut + nGridOut2
       return


       !\
       ! Edges going along resolution edge
       ! Opposite "faces" may have the same geometry in projection
       ! term "face" is used to refer to 4 upper-, left-, front-most etc. grid points
       ! they do not necesserily belong to the same plane
       !/

elseif(    abs(XyzGrid_DI(x_,1) - XyzGrid_DI(x_,5)) < dXyzSmall_D(x_).and.&
           abs(XyzGrid_DI(y_,1) - XyzGrid_DI(y_,5)) < dXyzSmall_D(y_)&
           .and.&
           abs(XyzGrid_DI(x_,2) - XyzGrid_DI(x_,6)) < dXyzSmall_D(x_).and.&
           abs(XyzGrid_DI(y_,2) - XyzGrid_DI(y_,6)) < dXyzSmall_D(y_)&
           .and.&
           abs(XyzGrid_DI(x_,3) - XyzGrid_DI(x_,7)) < dXyzSmall_D(x_).and.&
           abs(XyzGrid_DI(y_,3) - XyzGrid_DI(y_,7)) < dXyzSmall_D(y_)&
           .and.&
           abs(XyzGrid_DI(x_,4) - XyzGrid_DI(x_,8)) < dXyzSmall_D(x_).and.&
           abs(XyzGrid_DI(y_,4) - XyzGrid_DI(y_,8)) < dXyzSmall_D(y_))then
       !\
       ! "Faces" have the same geometry in projection on xy-plane
       !/
        iOrder_I = (/1,2,3,4,5,6,7,8/)
        call interpolate_resolution_edge(iAxis=x_,jAxis=y_)
        return

    elseif(abs(XyzGrid_DI(z_,1) - XyzGrid_DI(z_,2)) < dXyzSmall_D(z_).and.&
           abs(XyzGrid_DI(y_,1) - XyzGrid_DI(y_,2)) < dXyzSmall_D(y_)&
           .and.&
           abs(XyzGrid_DI(z_,3) - XyzGrid_DI(z_,4)) < dXyzSmall_D(z_).and.&
           abs(XyzGrid_DI(y_,3) - XyzGrid_DI(y_,4)) < dXyzSmall_D(y_)&
           .and.&
           abs(XyzGrid_DI(z_,5) - XyzGrid_DI(z_,6)) < dXyzSmall_D(z_).and.&
           abs(XyzGrid_DI(y_,5) - XyzGrid_DI(y_,6)) < dXyzSmall_D(y_)&
           .and.&
           abs(XyzGrid_DI(z_,7) - XyzGrid_DI(z_,8)) < dXyzSmall_D(z_).and.&
           abs(XyzGrid_DI(y_,7) - XyzGrid_DI(y_,8)) < dXyzSmall_D(y_))then
       !\
       ! "Faces" have the same geometry in projection on yz-plane
       !/
        iOrder_I = (/1,3,5,7,2,4,6,8/)
        call interpolate_resolution_edge(iAxis=y_,jAxis=z_)
       return


    elseif(abs(XyzGrid_DI(x_,1) - XyzGrid_DI(x_,3)) < dXyzSmall_D(x_).and.&
           abs(XyzGrid_DI(z_,1) - XyzGrid_DI(z_,3)) < dXyzSmall_D(z_)&
           .and.&
           abs(XyzGrid_DI(x_,2) - XyzGrid_DI(x_,4)) < dXyzSmall_D(x_).and.&
           abs(XyzGrid_DI(z_,2) - XyzGrid_DI(z_,4)) < dXyzSmall_D(z_)&
           .and.&
           abs(XyzGrid_DI(x_,5) - XyzGrid_DI(x_,7)) < dXyzSmall_D(x_).and.&
           abs(XyzGrid_DI(z_,5) - XyzGrid_DI(z_,7)) < dXyzSmall_D(z_)&
           .and.&
           abs(XyzGrid_DI(x_,6) - XyzGrid_DI(x_,8)) < dXyzSmall_D(x_).and.&
           abs(XyzGrid_DI(z_,6) - XyzGrid_DI(z_,8)) < dXyzSmall_D(z_))then
       !\
       ! "Faces" have the same geometry in projection on xz-plane
       !/
        iOrder_I = (/1,2,5,6,3,4,7,8/)
        call interpolate_resolution_edge(iAxis=x_,jAxis=z_)
        return
    end if



  contains
    subroutine pyramids(iTetrahedron1_I,iTetrahedron2_I,iTetrahedron3_I,iTetrahedron4_I,iTetrahedron5_I,iTetrahedron6_I,&
                        iRectangular1_I,iRectangular2_I,iRectangular3_I,iTrapezoidal1_I,iTrapezoidal2_I,iTrapezoidal3_I)
        integer,intent(in),optional::iTetrahedron1_I(4)
        integer,intent(in),optional::iTetrahedron2_I(4)
        integer,intent(in),optional::iTetrahedron3_I(4)
        integer,intent(in),optional::iTetrahedron4_I(4)
        integer,intent(in),optional::iTetrahedron5_I(4)
        integer,intent(in),optional::iTetrahedron6_I(4)
        integer,intent(in),optional::iRectangular1_I(5)
        integer,intent(in),optional::iRectangular2_I(5)
        integer,intent(in),optional::iRectangular3_I(5)
        integer,intent(in),optional::iTrapezoidal1_I(5)
        integer,intent(in),optional::iTrapezoidal2_I(5)
        integer,intent(in),optional::iTrapezoidal3_I(5)
    end subroutine pyramids
    !========================

    real function triple_product(a_D, b_D, c_D)
    real, dimension(nDim), intent(in) :: a_D, b_D, c_D
        !-----------------------------------------------
        triple_product = a_D(x_) * ( b_D(y_)*c_D(z_) - b_D(z_)*c_D(y_) ) - &
        a_D(y_) * ( b_D(x_)*c_D(z_) - b_D(z_)*c_D(x_) ) + &
        a_D(z_) * ( b_D(x_)*c_D(y_) - b_D(y_)*c_D(x_) )
    end function triple_product
    !=======================================================================
    real function dot_product_3(a_D, b_D)
        real,dimension(nDim),intent(in) :: a_D, b_D
        !----------
        dot_product_3 = a_D(x_)* b_D(x_) + a_D(y_)* b_D(y_) + a_D(z_)* b_D(z_)
    end function dot_product_3
    !=======================================================================
    function cross_product_3(a_D, b_d)
        real, dimension(nDim), intent(in) :: a_d, b_D
        real, dimension(nDim) :: cross_product_3
        !-----------------------------------------------
        cross_product_3(x_) = a_D(y_)*b_D(z_) - a_D(z_)*b_D(y_)
        cross_product_3(y_) = a_D(z_)*b_D(x_) - a_D(x_)*b_D(z_)
        cross_product_3(z_) = a_D(x_)*b_D(y_) - a_D(y_)*b_D(x_)
    end function cross_product_3
    !=======================================================================
    function vertical_distance_point_plane(XP_D, iOrderAux_I, iAxis)
        real, dimension(nDim), intent(in) :: XP_D
        integer, dimension(3), intent(in) :: iOrderAux_I   !Points defining the plane
        integer,                  intent(in) :: iAxis              !Vertical direction
        real, dimension(nDim) ::             XAuxV_D            !Vertical vector
        real ,dimension(nDim):: X1_D, X2_D, X3_D
        real ::vertical_distance_point_plane
        !----------------------------------------------
        XAuxV_D        = (/0,0,0/)
        XAuxV_D(iAxis) = 1
        X1_D = XyzGrid_DI(:,iOrderAux_I(iOrderAux_I(1)))
        X2_D = XyzGrid_DI(:,iOrderAux_I(iOrderAux_I(2)))
        X3_D = XyzGrid_DI(:,iOrderAux_I(iOrderAux_I(3)))
        vertical_distance_point_plane = &
            triple_product(XP_D - X1_D, X2_D - X1_D, X3_D- X1_D)/&
            triple_product(XAuxV_D,     X2_D - X1_D, X3_D- X1_D)
    end function vertical_distance_point_plane
    !=======================================================================
    real function vector_length_3(a_D)
        real, dimension(nDim), intent(in) :: a_D
        !-----------------------------------------------
        vector_length_3 = sqrt(a_D(x_)*a_D(x_)+a_D(y_)*a_D(y_)+a_D(z_)*a_D(z_))
    end function vector_length_3
    !==============================================================================
    subroutine interpolate_resolution_edge(iAxis,jAxis)
        integer, intent(in) :: iAxis,jAxis
        integer             :: kAxis
        real :: AlphaKAxisCoarse, AlphaKAxisFine
        real :: dXyzUp, dXyzDown ! Auxilary variables to determine vertical weight
        !\
        ! Order of the grid points 1:4 and 5:8 should satisfy consistency for 2D interpolation
        ! in respect of axes iAxis and jAxis
        !/
        !----------------------------------------------------------------------------------------
        Xyz2_D(x_:y_) = Xyz_D((/iAxis,jAxis/))
        kAxis = 6-iAxis-jAxis
!write(*,*)'edge'

        dXyzUp   = maxval(XyzGrid_DI(kAxis,:),iLevel_I==Coarse_) - Xyz_D(kAxis)
        dXyzDown = Xyz_D(kAxis) - minval(XyzGrid_DI(kAxis,:),iLevel_I==Coarse_)
        AlphaKAxisCoarse = dXyzUp / (dXyzUp + dXyzDown)

        dXyzUp   = maxval(XyzGrid_DI(kAxis,:),iLevel_I==Fine_) - Xyz_D(kAxis)
        dXyzDown = Xyz_D(kAxis) - minval(XyzGrid_DI(kAxis,:),iLevel_I==Fine_)
        AlphaKAxisFine =  dXyzUp / (dXyzUp + dXyzDown)

        !\
        ! Interpolate along lower face
        !/
        DoStencilFix2 = DoStencilFix
        XyzGrid2_DI(x_:y_,1:4) = XyzGrid_DI((/iAxis,jAxis/),iOrder_I(1:4))
        iLevel2_I(        1:4) = iLevel_I(                  iOrder_I(1:4))
        call interpolate_amr_grid_2(&
            Xyz2_D, XyzGrid2_DI, iLevel2_I, &
            nGridOut2, Weight_I(1:4), iOrder2_I, DoStencilFix2, XyzStencil2_D)

        if(DoStencilFix2)then
            DoStencilFix = DoStencilFix2
            XyzStencil_D((/iAxis,jAxis/)) = XyzStencil2_D(x_:y_)
            XyzStencil_D(kAxis) = sum(XyzGrid_DI(kAxis,:))/8
            return
        end if

        iOrder_I(1:4) = iOrder_I(iOrder2_I  )
        iOrder_I(5:8) = iOrder_I(iOrder2_I+4)

        !\
        ! May need to remove a grid points from the stencil once behind the boundary
        !/
        nGridOut = 2 * nGridOut2
        if(    nGridOut2==2)then
            iOrder_I(  3:4) = iOrder_I((/5,6/))
        elseif(nGridOut2==3)then
            iOrder_I(  4:6) = iOrder_I((/5,6,7/))
        end if

        !\
        ! Apply weight for interpolation along kAxis
        !/
        do iGrid = 1,nGridOut2
            if(iLevel_I(iOrder_I(iGrid))==Fine_)then
                Weight_I(nGridOut2+iGrid) = Weight_I(iGrid)*(1 - AlphaKAxisFine)
                Weight_I(          iGrid) = Weight_I(iGrid)*     AlphaKAxisFine
            else
                if(    iLevel_I(iOrder_I(          iGrid))==BehindTheBoundary_)then
                    Weight_I(nGridOut2+iGrid) = Weight_I(iGrid)
                    Weight_I(          iGrid) = 0
                elseif(iLevel_I(iOrder_I(nGridOut2+iGrid))==BehindTheBoundary_)then
                    Weight_I(nGridOut2+iGrid) = 0
                    Weight_I(          iGrid) = Weight_I(iGrid)
                else
                    Weight_I(nGridOut2+iGrid) = Weight_I(iGrid)*(1 - AlphaKAxisCoarse)
                    Weight_I(          iGrid) = Weight_I(iGrid)*     AlphaKAxisCoarse
                end if
            end if
        end do
        DoStencilFix = .false.
    end subroutine interpolate_resolution_edge

       
    !=======================================================================
    subroutine interpolate_tetrahedron(X1_D, X2_D, X3_D, X4_D)
        !\
        ! Interpolate in the tetrahedron
        !/
        real, dimension(nDim), intent(in) :: X1_D, X2_D, X3_D, X4_D
        real:: Aux 
        !-----------------------------------------------
        Aux= 1/triple_product( X4_D - X1_D,  X3_D - X1_D,  X2_D - X1_D)
        Weight_I(4) = triple_product(Xyz_D - X1_D,  X3_D - X1_D,  X2_D - X1_D)*Aux
        Weight_I(3) = triple_product( X4_D - X1_D, Xyz_D - X1_D,  X2_D - X1_D)*Aux
        Weight_I(2) = triple_product( X4_D - X1_D,  X3_D - X1_D, Xyz_D - X1_D)*Aux
        Weight_I(1) = 1 - sum(Weight_I(2:4))
    end subroutine interpolate_tetrahedron
    !=======================================================================
    subroutine interpolate_pyramid(X1_D, X2_D, X3_D, X4_D, X5_D)
        !\
        ! Interpolate in the pyramid with base X1X2X3X4 and apex X5
        ! valid for case of rectangular or trapezoid base
        !/
        real, dimension(nDim), intent(in) :: X1_D, X2_D, X3_D, X4_D, X5_D
        real, dimension(nDim) :: XyzP_D !Projection of Xyz point on the base of the pyramid X1X2X3X4 along the line X5Xyz
        real ::Alpha5, Aux
        !-----------------------------------------------

        Aux = 1/triple_product( X5_D - X1_D, X2_D - X1_D, X3_D - X1_D)
        Alpha5 = triple_product(X5_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)*Aux
        if(Alpha5==0)then
            Weight_I(5) = 1
            Weight_I(1:4) = (/0,0,0,0/)
        else
            XyzP_D = X5_D + (Xyz_D-X5_D)/Alpha5
        end if
        !\
        ! Find weight of the apex point X5_D
        !/
        Weight_I(5) = triple_product(Xyz_D - X1_D, X2_D - X1_D, X3_D - X1_D)*Aux

        !\
        ! Find weights of the base points X1, X2, X3, X4
        !/

        !\
        ! Check which 2 egdes are parallel
        !/
        if(vector_length_3(cross_product_3(X2_D-X1_D,X4_D-X3_D)).le.vector_length_3(dXyzSmall_D))then
            !\
            ! X2-X1 and X4-X3 are parallel    F F    OR     C  C
            !/                              C     C         C  C
            Weight_I(3) = 0.50 * dot_product_3( X3_D-X1_D + X4_D-X2_D, XyzP_D - X1_D         )/&
                         (0.25 * dot_product_3( X3_D-X1_D + X4_D-X2_D, X3_D-X1_D + X4_D-X2_D ) )
            Weight_I(4) = Weight_I(3)

            Weight_I(1) = 1 - Weight_I(3)
            Weight_I(2) = Weight_I(1)

            Weight_I(2) = Weight_I(2) * dot_product_3(XyzP_D - X1_D, X2_D - X1_D)/&
                                        dot_product_3(  X2_D - X1_D, X2_D - X1_D)
            Weight_I(1) = Weight_I(1) - Weight_I(2)

            Weight_I(4) = Weight_I(4) * dot_product_3(XyzP_D - X3_D, X4_D - X3_D)/&
                                        dot_product_3(  X4_D - X3_D, X4_D - X3_D)
            Weight_I(3) = Weight_I(3) - Weight_I(4)

        else
            !\
            ! X3-X1 and X4-X2 are parallel  C  F
            !/                              C  F
            Weight_I(2) = 0.50 * dot_product_3( X2_D-X1_D + X4_D-X3_D, XyzP_D - X1_D         )/&
                         (0.25 * dot_product_3( X2_D-X1_D + X4_D-X3_D, X2_D-X1_D + X4_D-X3_D ) )
            Weight_I(4) = Weight_I(2)

            Weight_I(1) = 1 - Weight_I(2)
            Weight_I(3) = Weight_I(1)

            Weight_I(3) = Weight_I(3) * dot_product_3(XyzP_D - X1_D, X3_D - X1_D)/&
                                        dot_product_3(  X3_D - X1_D, X3_D - X1_D)
            Weight_I(1) = Weight_I(1) - Weight_I(3)

            Weight_I(4) = Weight_I(4) * dot_product_3(XyzP_D - X2_D, X4_D - X2_D)/&
                                        dot_product_3(  X4_D - X2_D, X4_D - X2_D)
            Weight_I(2) = Weight_I(2) - Weight_I(4)
        end if

        !\
        ! Correct weights due to pyramid geometry
        !/
        Weight_I(1) = Weight_I(1) * (1 - Weight_I(5))
        Weight_I(2) = Weight_I(2) * (1 - Weight_I(5))
        Weight_I(3) = Weight_I(3) * (1 - Weight_I(5))
        Weight_I(4) = Weight_I(4) * (1 - Weight_I(5))

    end subroutine interpolate_pyramid


  end subroutine interpolate_amr_grid_3
  !=========



  !===========================TESTS======================================================
  subroutine test_interpolate_amr_grid_2(nSample)
    !The test subtoutine nSample times generate random point in the AMR domain
    !and compares the result of the coordinate interpolation to this point with the
    !gampled coordinates
    !
    integer, intent(in) :: nSample

    integer :: iSeed=0, i, j, iBlock, loc(3), iSample, iGrid, iDim
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
    if(sum(abs(XyzSquaredInterpolated1_D(:)-XyzSquaredInterpolated2_D(:)-Xyz_D(:)*Xyz_D(:)+Xyz2_D(:)*Xyz2_D(:))) > 5.0e-4)then
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
            do iDim =1, nDim
                XyzGrid_DI(iDim,iGrid) = -1.0
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


end module ModInterpolateAMRGrid
