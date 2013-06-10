module ModInterpolateAMRGrid
  implicit none

  !Generalize bilinear and trilinear interpolation for AMR grids
  !The data are given at the cell-centered grid which consists of AMR blocks.
  !The interpolation is free of any jumps at the resolution interfaces
  !including edges and corners
  
  PRIVATE !Except
  SAVE
  public:: interpolate_on_amr_grid_1 ! One-dimensional - for completeness only

  ! Two-dimensional grid, is also used to interpolate the face projection for 3D grid

  public:: interpolate_on_amr_grid_2,test_interpolate_on_amr_grid_2 
  

  ! Three-dimensional grid, 

  public :: interpolate_on_amr_grid_3

  integer,parameter:: BehindTheBoundary_ = -7777, &
                      Coarse_ = 0,               &
                      Fine_   = 1,               &
                      x_       = 1,               &
                      y_       = 2,               &
                      z_       = 3

  public :: BehindTheBoundary_
  
contains
  !================
  subroutine  interpolate_on_amr_grid_1(&
                  nIndexes, Xyz_D, XyzGrid_DI, iIndexes_II, iLevel_I,&
                  nGridOut, Weight_I)

    integer,parameter :: nGrid = 2, nDim = 1

    character(LEN=*),parameter:: NameSub='interpolate_on_amr_grid_1'

    !\
    !Input parameters
    !/
    !The number of cell-block-Pe indexes:
    integer, intent(in) :: nIndexes    !(iCell,jCell,kCell,iBlock, Pe) - nIndexes=5
    real   , intent(in) :: Xyz_D(nDim) !The location at which to interpolate the data

    !\
    !Inout parameters
    !/
    !Grid point coordinates
    real  ,   intent(inout) :: XyzGrid_DI(nDim,nGrid) !1 coordinate, 2 points

    !Cell-Block(-Pe) indexes for grid points
    integer,  intent(inout) :: iIndexes_II(nIndexes,nGrid)  

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
       nGridOut = 1
       Weight_I(1) = 1
    elseif(iLevel_I(1)== BehindTheBoundary_)then
       nGridOut = 1
       Weight_I(1) = 1
       iLevel_I(1) = iLevel_I(2)
       iIndexes_II(:,1) =  iIndexes_II(:,2)
    else
       nGridOut = 2
       Weight_I(1) = ( XyzGrid_DI(1,2) - Xyz_D(1) )/( XyzGrid_DI(1,2) - XyzGrid_DI(1,1) )
       Weight_I(2) = 1 - Weight_I(1) 
    end if
  
  end subroutine interpolate_on_amr_grid_1    
  !================
  subroutine  interpolate_on_amr_grid_2(&
                  nIndexes, Xyz_D, XyzGrid_DI, iIndexes_II, iLevel_I,&
                  nGridOut, Weight_I, DoStencilFix, XyzStencil_D)

    integer,parameter :: nGrid = 4, nDim = 2

    character(LEN=*),parameter:: NameSub='interpolate_on_amr_grid_2'

    !\
    !Input parameters
    !/
    !The number of cell(-block(-Pe)) indexes:
    integer, intent(in) :: nIndexes    !(iCell,jCell,kCell,iBlock, Pe) - nIndexes=5
    real   , intent(in) :: Xyz_D(nDim) !The location at which to interpolate the data

    !\
    !Inout parameters
    !/
    !Grid point coordinates
    real  ,   intent(inout) :: XyzGrid_DI(nDim,nGrid) !1 coordinate, 2 points

    !Cell-Block(-Pe) indexes for grid points
    integer,  intent(inout) :: iIndexes_II(nIndexes,nGrid)  

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

    integer, dimension(nIndexes) :: i1_I, i2_I, i3_I, i4_I !Pass indexes to "triangulate"

    real    :: XyzMin_D(nDim), XyzMax_D(nDim)    ! Minimum and maximum values of coordinates
    real    :: dXyzSmall_D(nDim)                 ! To find lines with the "constant" cordinate 
    real    :: dXyzInv_D(  nDim)                 ! for calculating on a uniform grid 
    real    :: Aux_D(nDim)                       ! Misc
    real    :: Xyz1_D(1), XyzGrid1_DI(1,2)       ! To call one-dimensional interpolation 
    integer :: iIndexes1_II(nIndexes,2)          ! To call one-dimensional interpolation
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
         XyzGrid_DI(x_,2) >  Xyz_D(x_).and.XyzGrid_DI(y_,2).le.Xyz_D(y_).and.&
         XyzGrid_DI(x_,3).le.Xyz_D(x_).and.XyzGrid_DI(y_,3) >  Xyz_D(y_).and.&
         XyzGrid_DI(x_,4) >  Xyz_D(x_).and.XyzGrid_DI(y_,4) >  Xyz_D(y_)) )then
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
    if(all(iLevel_I(3:4)== BehindTheBoundary_).or.all(XyzGrid_DI(y_,1:2)==Xyz_D(y_)) )then
       !\
       ! Edge 3:4 should be removed
       !/
       Xyz1_D(x_) = Xyz_D(x_) ; XyzGrid1_DI(x_,:) = XyzGrid_DI(x_,1:2)

       call interpolate_on_amr_grid_1(&
            nIndexes, Xyz1_D,  XyzGrid1_DI, iIndexes_II(:,1:2), iLevel_I(1:2), &
            nGridOut, Weight_I(1:2) )
       return

    elseif(all(iLevel_I(1:2)== BehindTheBoundary_)                                    )then
       !\
       ! Edge 1:2 should be removed
       !/
       Xyz1_D(x_) = Xyz_D(x_) ; XyzGrid1_DI(x_,:) = XyzGrid_DI(x_,3:4)
       iIndexes_II(:,1:2) = iIndexes_II(:,3:4)
       iLevel_I(     1:2) = iLevel_I(     3:4)

       call interpolate_on_amr_grid_1(&
            nIndexes, Xyz1_D,  XyzGrid1_DI, iIndexes_II(:,1:2), iLevel_I(1:2), &
            nGridOut, Weight_I(1:2) )
       return

    elseif(all(iLevel_I(2:4:2)== BehindTheBoundary_).or. all(XyzGrid_DI(x_,1:3:2)==Xyz_D(x_)) )then
       !\
       ! Edge 2,4 should be removed
       !/
       Xyz1_D(x_) = Xyz_D(y_) ; XyzGrid1_DI(x_,:) = XyzGrid_DI(y_,1:3:2)
       iIndexes_II(:,2) = iIndexes_II(:,3)
       iLevel_I(     2) = iLevel_I(     3)

       call interpolate_on_amr_grid_1(&
            nIndexes, Xyz1_D,  XyzGrid1_DI, iIndexes_II(:,1:2), iLevel_I(1:2), &
            nGridOut, Weight_I(1:2) )
       return

    elseif(all(iLevel_I(1:3:2)== BehindTheBoundary_)                                        )then
       !\
       ! Edge 1,3 should be removed
       !/
       Xyz1_D(x_) = Xyz_D(y_) ; XyzGrid1_DI(x_,:) = XyzGrid_DI(y_,2:4:2)
       iIndexes_II(:,1) = iIndexes_II(:,2)
       iLevel_I(     1) = iLevel_I(     2)
       iIndexes_II(:,2) = iIndexes_II(:,4)
       iLevel_I(     2) = iLevel_I(     4)

       call interpolate_on_amr_grid_1(&
            nIndexes, Xyz1_D,  XyzGrid1_DI, iIndexes_II(:,1:2), iLevel_I(1:2), &
            nGridOut, Weight_I(1:2) )
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

       Aux_D(y_) = dXyzInv_D(y_)*(Xyz_D(y_)-XyzGrid_DI(y_,1))

       Xyz1_D(x_) = Xyz_D(x_)  

       !\
       ! Interpolate along lower edge
       !/
       XyzGrid1_DI(x_,:) = XyzGrid_DI(x_,1:2)
       call interpolate_on_amr_grid_1(&
            nIndexes, Xyz1_D,  XyzGrid1_DI, iIndexes_II(:,1:2), iLevel_I(1:2), &
            nGridOut1, Weight_I(1:2) )
       !\
       ! Apply weight for interpolation along y_
       !/
       Weight_I(1:nGridOut1) =  Weight_I(1:nGridOut1) * (1 - Aux_D(y_))
       !\
       ! May need to remove a grid point from the stencil once behind the boundary
       !/
       if(nGridOut1==1)then
          iIndexes_II(:,2:3) = iIndexes_II(:,3:4)
          iLevel_I(     2:3) = iLevel_I(     3:4)
       end if
       nGridOut = nGridOut1
       !\
       ! Interpolate along upper edge
       !/
       XyzGrid1_DI(x_,:) = XyzGrid_DI(x_,3:4)
       call interpolate_on_amr_grid_1(&
            nIndexes, Xyz1_D,  XyzGrid1_DI, iIndexes_II(:,nGridOut+1:nGridOut+2), &
            iLevel_I(nGridOut+1:nGridOut+2), nGridOut1, Weight_I(nGridOut+1:nGridOut+2) )
       !\
       ! Apply weight for interpolation along y_
       !/
       Weight_I(nGridOut+1:nGridOut+nGridOut1) =  Weight_I(nGridOut+1:nGridOut+nGridOut1) * Aux_D(y_)

       !\
       ! May need to remove a grid point from the stencil once behind the boundary
       !/
       nGridOut = nGridOut + nGridOut1
       return       
   
    elseif(abs(XyzGrid_DI(x_,1) - XyzGrid_DI(x_,3)) < dXyzSmall_D(x_).and.&
           abs(XyzGrid_DI(x_,2) - XyzGrid_DI(x_,4)) < dXyzSmall_D(x_))then

       !\                             C  F         F C
       ! Edges going along y-axis        F    or   F 
       !/                             C              C

       Aux_D(x_) = dXyzInv_D(x_)*(Xyz_D(x_)-XyzGrid_DI(x_,1))

       Xyz1_D(x_) = Xyz_D(y_)  

       !\
       ! Interpolate along left edge
       !/
       XyzGrid1_DI(x_,:) = XyzGrid_DI(y_,1:3:2)
       call interpolate_on_amr_grid_1(&
            nIndexes, Xyz1_D,  XyzGrid1_DI, iIndexes_II(:,1:3:2), iLevel_I(1:3:2), &
            nGridOut1, Weight_I(1:3:2) )
      
       !\
       ! May need to remove a grid point from the stencil once behind the boundary
       !/
       if(nGridOut1==1)then
          !\
          ! Apply weight for interpolation along x_
          !/
          Weight_I(1) =  Weight_I(1) * (1 - Aux_D(x_))
          iIndexes_II(:,3) = iIndexes_II(:,4)
          iLevel_I(     3) = iLevel_I(     4)
       else
          Weight_I(1:3:2) =  Weight_I(1:3:2) * (1 - Aux_D(x_))
       end if
       nGridOut = nGridOut1
       !\
       ! Interpolate along right edge
       !/
       XyzGrid1_DI(x_,:) = XyzGrid_DI(y_,2:4:2)
       call interpolate_on_amr_grid_1(&
            nIndexes, Xyz1_D,  XyzGrid1_DI, iIndexes_II(:,2:nGridOut+2:nGridOut), &
            iLevel_I(2:nGridOut+2:nGridOut), nGridOut1, Weight_I(2:nGridOut+2:nGridOut) )
       !\
       ! Apply weight for interpolation along x_
       !/
       if(nGridOut1==1)then
          Weight_I(2) = Weight_I(2) * Aux_D(x_)
       else
          Weight_I(2:nGridOut+2:nGridOut) = Weight_I(2:nGridOut+2:nGridOut) * Aux_D(x_)
       end if

       !\
       ! May need to remove a grid point from the stencil once behind the boundary
       !/

       nGridOut = nGridOut + nGridOut1
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
          call interpolate_triangle(X1_D=XyzGrid_DI(:,1),X2_D=XyzGrid_DI(:,2),X3_D=XyzGrid_DI(:,3))
          DoStencilFix = any(Weight_I(1:3)<0.0)
          if(DoStencilFix)then
             XyzStencil_D(y_) = ( XyzGrid_DI(y_,1) + XyzGrid_DI(y_,2) + XyzGrid_DI(y_,3) )/3
             XyzStencil_D(x_) =  -XyzGrid_DI(x_,1) + XyzGrid_DI(x_,2) + XyzGrid_DI(x_,3) 
          else
             nGridOut = 3
          end if
          return
          
       elseif(abs(XyzGrid_DI(x_,4) - 0.50*(XyzGrid_DI(x_,2) + XyzGrid_DI(x_,1)))<dXyzSmall_D(x_))then
          !\                (?) C                                                     C
          ! The configuration  F F    Check if the Xyz point belongs to the triangle F F
          !/
          call interpolate_triangle(X1_D=XyzGrid_DI(:,1),X2_D=XyzGrid_DI(:,2),X3_D=XyzGrid_DI(:,4))
          DoStencilFix = any(Weight_I(1:3)<0.0)
          if(DoStencilFix)then
             XyzStencil_D(y_) = ( XyzGrid_DI(y_,1) + XyzGrid_DI(y_,2) + XyzGrid_DI(y_,4) )/3
             XyzStencil_D(x_) =   XyzGrid_DI(x_,1) - XyzGrid_DI(x_,2) + XyzGrid_DI(x_,4) 
          else
             nGridOut = 3
             iIndexes_II(:,3) = iIndexes_II(:,4)
          end if
          return
       end if
    end if

    if(abs(XyzGrid_DI(y_,4) - XyzGrid_DI(y_,3)) < dXyzSmall_D(y_).and.all(iLevel_I(3:4)==1))then

       if(abs(XyzGrid_DI(x_,1) - 0.50*(XyzGrid_DI(x_,3) + XyzGrid_DI(x_,4)))<dXyzSmall_D(x_))then
          !\                   F F                                                      F F
          ! The configuration   C (?)    Check if the Xyz point belongs to the triangle  C
          !/
          call interpolate_triangle(X1_D=XyzGrid_DI(:,1),X2_D=XyzGrid_DI(:,3),X3_D=XyzGrid_DI(:,4))
          DoStencilFix = any(Weight_I(1:3)<0.0)
          if(DoStencilFix)then
             XyzStencil_D(y_) = ( XyzGrid_DI(y_,1) + XyzGrid_DI(y_,3) + XyzGrid_DI(y_,4) )/3
             XyzStencil_D(x_) =   XyzGrid_DI(x_,1) - XyzGrid_DI(x_,3) + XyzGrid_DI(x_,4) 
          else
             iIndexes_II(:,2) = iIndexes_II(:,3)
             iIndexes_II(:,3) = iIndexes_II(:,4)
             nGridOut = 3
          end if
          return

       elseif(abs(XyzGrid_DI(x_,2) - 0.50*(XyzGrid_DI(x_,3) + XyzGrid_DI(x_,4)))<dXyzSmall_D(x_))then
          !\                     F F                                                   F F
          ! The configuration (?) C    Check if the Xyz point belongs to the triangle   C
          !/
          call interpolate_triangle(X1_D=XyzGrid_DI(:,2),X2_D=XyzGrid_DI(:,3),X3_D=XyzGrid_DI(:,4))
          DoStencilFix = any(Weight_I(1:3)<0.0)
          if(DoStencilFix)then
             XyzStencil_D(y_) = ( XyzGrid_DI(y_,2) + XyzGrid_DI(y_,3) + XyzGrid_DI(y_,4) )/3
             XyzStencil_D(x_) =   XyzGrid_DI(x_,2) + XyzGrid_DI(x_,3) - XyzGrid_DI(x_,4) 
          else
             nGridOut = 3
             iIndexes_II(:,1) = iIndexes_II(:,2)
             iIndexes_II(:,2) = iIndexes_II(:,3)
             iIndexes_II(:,3) = iIndexes_II(:,4)
          end if
          return
       end if
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
          call interpolate_triangle(X1_D=XyzGrid_DI(:,1),X2_D=XyzGrid_DI(:,2),X3_D=XyzGrid_DI(:,3))
          DoStencilFix = any(Weight_I(1:3)<0.0)
          if(DoStencilFix)then
             XyzStencil_D(y_) = ( XyzGrid_DI(y_,1) + XyzGrid_DI(y_,2) + XyzGrid_DI(y_,3) )/3
             XyzStencil_D(x_) =  -XyzGrid_DI(x_,1) + XyzGrid_DI(x_,2) + XyzGrid_DI(x_,3) 
          else
             nGridOut = 3
          end if
          return
          
       elseif(abs(XyzGrid_DI(x_,4) - 0.50*(XyzGrid_DI(x_,2) + XyzGrid_DI(x_,1)))<dXyzSmall_D(x_))then
          !\                (?) C                                                     C
          ! The configuration  F F    Check if the Xyz point belongs to the triangle F F
          !/
          call interpolate_triangle(X1_D=XyzGrid_DI(:,1),X2_D=XyzGrid_DI(:,2),X3_D=XyzGrid_DI(:,4))
          DoStencilFix = any(Weight_I(1:3)<0.0)
          if(DoStencilFix)then
             XyzStencil_D(y_) = ( XyzGrid_DI(y_,1) + XyzGrid_DI(y_,2) + XyzGrid_DI(y_,4) )/3
             XyzStencil_D(x_) =   XyzGrid_DI(x_,1) - XyzGrid_DI(x_,2) + XyzGrid_DI(x_,4) 
          else
             nGridOut = 3
             iIndexes_II(:,3) = iIndexes_II(:,4)
          end if
          return
       end if
    end if

    if(abs(XyzGrid_DI(y_,4) - XyzGrid_DI(y_,3)) < dXyzSmall_D(y_).and.all(iLevel_I(3:4)==1))then

       if(abs(XyzGrid_DI(x_,1) - 0.50*(XyzGrid_DI(x_,3) + XyzGrid_DI(x_,4)))<dXyzSmall_D(x_))then
          !\                   F F                                                      F F
          ! The configuration   C (?)    Check if the Xyz point belongs to the triangle  C
          !/
          call interpolate_triangle(X1_D=XyzGrid_DI(:,1),X2_D=XyzGrid_DI(:,3),X3_D=XyzGrid_DI(:,4))
          DoStencilFix = any(Weight_I(1:3)<0.0)
          if(DoStencilFix)then
             XyzStencil_D(y_) = ( XyzGrid_DI(y_,1) + XyzGrid_DI(y_,3) + XyzGrid_DI(y_,4) )/3
             XyzStencil_D(x_) =   XyzGrid_DI(x_,1) - XyzGrid_DI(x_,3) + XyzGrid_DI(x_,4) 
          else
             iIndexes_II(:,2) = iIndexes_II(:,3)
             iIndexes_II(:,3) = iIndexes_II(:,4)
             nGridOut = 3
          end if
          return

       elseif(abs(XyzGrid_DI(x_,2) - 0.50*(XyzGrid_DI(x_,3) + XyzGrid_DI(x_,4)))<dXyzSmall_D(x_))then
          !\                     F F                                                   F F
          ! The configuration (?) C    Check if the Xyz point belongs to the triangle   C
          !/
          call interpolate_triangle(X1_D=XyzGrid_DI(:,1),X2_D=XyzGrid_DI(:,2),X3_D=XyzGrid_DI(:,4))
          DoStencilFix = any(Weight_I(1:3)<0.0)
          if(DoStencilFix)then
             XyzStencil_D(y_) = ( XyzGrid_DI(y_,2) + XyzGrid_DI(y_,3) + XyzGrid_DI(y_,4) )/3
             XyzStencil_D(x_) =   XyzGrid_DI(x_,2) + XyzGrid_DI(x_,3) - XyzGrid_DI(x_,4) 
          else
             nGridOut = 3
             iIndexes_II(:,1) = iIndexes_II(:,2)
             iIndexes_II(:,2) = iIndexes_II(:,3)
             iIndexes_II(:,3) = iIndexes_II(:,4)
          end if
          return
       end if
    end if

100 continue

    if(count(iLevel_I==Fine_)==1.and.count(iLevel_I==Coarse_)==3)then
       Loc = maxloc(iLevel_I)
       select case(Loc(1))
       case(1)
          !   C---C
          !    \ /|
          !     F-C
          i1_I=iIndexes_II(:,1)
          i2_I=iIndexes_II(:,4)
          i3_I=iIndexes_II(:,2)
          i4_I=iIndexes_II(:,3)
          call triangulate(X1_D=XyzGrid_DI(:,1),&
                           X2_D=XyzGrid_DI(:,4),&
                           X3_D=XyzGrid_DI(:,2),&
                           X4_D=XyzGrid_DI(:,3) )
          
       case(2)
          !   C---C
          !   |\ / 
          !   C-F
          i1_I=iIndexes_II(:,2)
          i2_I=iIndexes_II(:,3)
          i3_I=iIndexes_II(:,1)
          i4_I=iIndexes_II(:,4)
          call triangulate(X1_D=XyzGrid_DI(:,2), &
                           X2_D=XyzGrid_DI(:,3), &
                           X3_D=XyzGrid_DI(:,1), &
                           X4_D=XyzGrid_DI(:,4) )
          
       case(3)
          !     F-C
          !    / \|
          !    C--C
          i1_I=iIndexes_II(:,3)
          i2_I=iIndexes_II(:,2)
          i3_I=iIndexes_II(:,1)
          i4_I=iIndexes_II(:,4)
          call triangulate(X1_D=XyzGrid_DI(:,3),&
                           X2_D=XyzGrid_DI(:,2),&
                           X3_D=XyzGrid_DI(:,1),&
                           X4_D=XyzGrid_DI(:,4) )
       case(4)
          !   C-F
          !   |/\ 
          !   C--C
          i1_I=iIndexes_II(:,4)
          i2_I=iIndexes_II(:,1)
          i3_I=iIndexes_II(:,2)
          i4_I=iIndexes_II(:,3)
          call triangulate(X1_D=XyzGrid_DI(:,4),&
                           X2_D=XyzGrid_DI(:,1),&
                           X3_D=XyzGrid_DI(:,2),&
                           X4_D=XyzGrid_DI(:,3) )
 
       end select
       return
    elseif(count(iLevel_I==Fine_)==2.and.count(iLevel_I==Coarse_)==2)then
       if(iLevel_I(1)==Fine_)then
          !   C-F
          !   |/\ 
          !   F--C
          i1_I=iIndexes_II(:,1)
          i2_I=iIndexes_II(:,4)
          i3_I=iIndexes_II(:,2)
          i4_I=iIndexes_II(:,3)
          call triangulate(X1_D=XyzGrid_DI(:,1),&
                           X2_D=XyzGrid_DI(:,4),&
                           X3_D=XyzGrid_DI(:,2),&
                           X4_D=XyzGrid_DI(:,3) )

       else
          !   F---C
          !   |\ / 
          !   C-F
          i1_I=iIndexes_II(:,2)
          i2_I=iIndexes_II(:,3)
          i3_I=iIndexes_II(:,1)
          i4_I=iIndexes_II(:,4)
          call triangulate(X1_D=XyzGrid_DI(:,2),&
                           X2_D=XyzGrid_DI(:,3),&
                           X3_D=XyzGrid_DI(:,1),&
                           X4_D=XyzGrid_DI(:,4) )

       end if
       return
    elseif(count(iLevel_I==Fine_)==3.and.count(iLevel_I==Coarse_)==1)then
       Loc = minloc(iLevel_I)
       select case(Loc(1))
       case(1)
          !   F---F
          !    \ /|
          !     C-F
          i1_I=iIndexes_II(:,1)
          i2_I=iIndexes_II(:,4)
          i3_I=iIndexes_II(:,2)
          i4_I=iIndexes_II(:,3) 
          call triangulate(X1_D=XyzGrid_DI(:,1),&
                           X2_D=XyzGrid_DI(:,4),&
                           X3_D=XyzGrid_DI(:,2),&
                           X4_D=XyzGrid_DI(:,3))

       case(2)
          !   F---F
          !   |\ / 
          !   F-C
          i1_I=iIndexes_II(:,2)
          i2_I=iIndexes_II(:,3)
          i3_I=iIndexes_II(:,1)
          i4_I=iIndexes_II(:,4)
          call triangulate(X1_D=XyzGrid_DI(:,2),&
                           X2_D=XyzGrid_DI(:,3),&
                           X3_D=XyzGrid_DI(:,1),&
                           X4_D=XyzGrid_DI(:,4) )
       case(3)
          !     C-F
          !    / \|
          !    F--F
          i1_I=iIndexes_II(:,3)
          i2_I=iIndexes_II(:,2)
          i3_I=iIndexes_II(:,1)
          i4_I=iIndexes_II(:,4)
          call triangulate(X1_D=XyzGrid_DI(:,3),&
                           X2_D=XyzGrid_DI(:,2),&
                           X3_D=XyzGrid_DI(:,1),&
                           X4_D=XyzGrid_DI(:,4) )
       case(4)
          !   F-C
          !   |/\ 
          !   F--F
          i1_I=iIndexes_II(:,4)
          i2_I=iIndexes_II(:,1)
          i3_I=iIndexes_II(:,2)
          i4_I=iIndexes_II(:,3)
          call triangulate(X1_D=XyzGrid_DI(:,4),&
                           X2_D=XyzGrid_DI(:,1),&
                           X3_D=XyzGrid_DI(:,2),&
                           X4_D=XyzGrid_DI(:,3))

       end select
       return
    else
       write(*,*)'Algorithm failuire in '//NameSub
       do iGrid = 1,nGrid
          write(*,*)'iGrid=',iGrid,' iLevel_I(iGrid)=',iLevel_I(iGrid),' XyzGrid_DI=', XyzGrid_DI(:,iGrid)
       end do
       call CON_stop('Stop in '//NameSub)
       
    end if
  contains
    !========
    real function cross_product_2(a_D, b_D)
      real,dimension(nDim),intent(in) :: a_D, b_D
      !----------
      cross_product_2 = a_D(x_)* b_D(y_) - a_d(y_)*b_D(x_)
    end function cross_product_2
    !=======
    subroutine triangulate(X1_D,X2_D,X3_D,X4_D)
      !\
      !Points 1 and 2 are on the shared side of the triangles, 3 and 4 are off
      !/
      real, dimension(nDim),     intent(in) :: X1_D, X2_D, X3_D, X4_D
      real :: Alpha2, Alpha3
      !-------
      Alpha3 = cross_product_2(Xyz_D - X1_D, X2_D-X1_D)/cross_product_2(X3_D - X1_D, X2_D-X1_D)
      if(Alpha3==0.0)then
         nGridOut = 2
         iIndexes_II(:,1) = i1_I
         iIndexes_II(:,2) = i2_I
         Alpha2 = cross_product_2(X3_D-X1_D,Xyz_D - X1_D)/cross_product_2(X3_D - X1_D, X2_D-X1_D)
         Weight_I(2) = Alpha2
      elseif(Alpha3 > 0.0)then
         nGridOut = 3
         iIndexes_II(:,1) = i1_I
         iIndexes_II(:,2) = i2_I
         iIndexes_II(:,3) = i3_I
         Alpha2 = cross_product_2(X3_D-X1_D,Xyz_D - X1_D)/cross_product_2(X3_D - X1_D, X2_D-X1_D)
         Weight_I(2) = Alpha2
         Weight_I(3) = Alpha3
      else
         Alpha3 = cross_product_2(Xyz_D - X1_D, X2_D-X1_D)/cross_product_2(X4_D - X1_D, X2_D-X1_D)
         nGridOut = 3
         iIndexes_II(:,1) = i1_I
         iIndexes_II(:,2) = i2_I
         iIndexes_II(:,3) = i4_I
         Alpha2 = cross_product_2(X4_D-X1_D,Xyz_D - X1_D)/cross_product_2(X4_D - X1_D, X2_D-X1_D)
         Weight_I(2) = Alpha2
         Weight_I(3) = Alpha3
      end if
         
      Weight_I(1) = 1 - sum(Weight_I(2:nGridOut))
    end subroutine triangulate
    !=========================
    subroutine interpolate_triangle(X1_D, X2_D, X3_D)
      !\
      ! Interpolate on the triangle
      !/
      real, dimension(nDim),     intent(in) :: X1_D, X2_D, X3_D
      !-------
      Weight_I(3) = cross_product_2(Xyz_D - X1_D, X2_D-X1_D)/cross_product_2(X3_D - X1_D, X2_D-X1_D)
      Weight_I(2) = cross_product_2(X3_D-X1_D,Xyz_D - X1_D)/cross_product_2(X3_D - X1_D, X2_D-X1_D)   
      Weight_I(1) = 1 - sum(Weight_I(2:3))
    end subroutine interpolate_triangle
  end subroutine interpolate_on_amr_grid_2
  !=========
  subroutine  interpolate_on_amr_grid_3(&
                  nIndexes, Xyz_D, XyzGrid_DI, iIndexes_II, iLevel_I,&
                  nGridOut, Weight_I, DoStencilFix, XyzStencil_D)

    integer,parameter :: nGrid = 8, nDim = 3

    character(LEN=*),parameter:: NameSub='interpolate_on_amr_grid_3'

    !\
    !Input parameters
    !/
    !The number of cell-block-Pe indexes:
    integer, intent(in) :: nIndexes    !(iCell,jCell,kCell,iBlock, Pe) - nIndexes=5
    real   , intent(in) :: Xyz_D(nDim) !The location at which to interpolate the data

    !\
    !Inout parameters
    !/
    !Grid point coordinates
    real  ,   intent(inout) :: XyzGrid_DI(nDim,nGrid) !1 coordinate, 2 points

    !Cell-Block(-Pe) indexes for grid points
    integer,  intent(inout) :: iIndexes_II(nIndexes,nGrid)  

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

    integer :: iGrid , nGridOut1, iLevelMin, iDim
    real    :: dXyzSmall_D(nDim), dXyzInv_D(nDim)
    real    :: Xyz2_D(x_:y_), XyzGrid2_DI(x_:y_,4), XyzStencil2_D(x_:y_) 
    integer :: iIndexes2_II(nIndexes,4)
    !-------------
    !\
    ! Check consistency of the input parameters.
    !/
    if(.not.(&
         XyzGrid_DI(x_,1).le.Xyz_D(x_).and.XyzGrid_DI(y_,1).le.Xyz_D(y_)&
                                      .and.XyzGrid_DI(z_,1).le.Xyz_D(z_).and.&
         XyzGrid_DI(x_,2) >  Xyz_D(x_).and.XyzGrid_DI(y_,2).le.Xyz_D(y_)&
                                      .and.XyzGrid_DI(z_,2).le.Xyz_D(z_).and.&
         XyzGrid_DI(x_,3).le.Xyz_D(x_).and.XyzGrid_DI(y_,3) >  Xyz_D(y_)&
                                      .and.XyzGrid_DI(z_,3).le.Xyz_D(z_).and.&
         XyzGrid_DI(x_,4) >  Xyz_D(x_).and.XyzGrid_DI(y_,4) >  Xyz_D(y_)&
                                      .and.XyzGrid_DI(z_,4).le.Xyz_D(z_).and.&
         XyzGrid_DI(x_,5).le.Xyz_D(x_).and.XyzGrid_DI(y_,5).le.Xyz_D(y_)&
                                      .and.XyzGrid_DI(z_,5) >  Xyz_D(z_).and.&
         XyzGrid_DI(x_,6) >  Xyz_D(x_).and.XyzGrid_DI(y_,6).le.Xyz_D(y_)&
                                      .and.XyzGrid_DI(z_,6) >  Xyz_D(z_).and.&
         XyzGrid_DI(x_,7).le.Xyz_D(x_).and.XyzGrid_DI(y_,7) >  Xyz_D(y_)&
                                      .and.XyzGrid_DI(z_,7) >  Xyz_D(z_).and.&
         XyzGrid_DI(x_,8) >  Xyz_D(x_).and.XyzGrid_DI(y_,8) >  Xyz_D(y_)&
                                      .and.XyzGrid_DI(z_,8) >  Xyz_D(z_))  )then
       write(*,*)'Inconsistent input in '//NameSub
       write(*,*)'Location at which to interpolate ', Xyz_D
       write(*,*)'Grid points:'
       do iGrid=1,nGrid
          write(*,*)XyzGrid_DI(:,iGrid)
       end do
       call CON_stop('Stop the code in '//NameSub)
    end if

    if(all(iLevel_I(5:8)== BehindTheBoundary_).or. all(XyzGrid_DI(z_,1:4)==Xyz_D(z_)) )then

       Xyz2_D(x_:y_) = Xyz_D(x_:y_) ; XyzGrid2_DI(x_:y_,:) = XyzGrid_DI(x_:y_,1:4)

       call interpolate_on_amr_grid_2(&
            nIndexes, Xyz2_D,  XyzGrid2_DI, iIndexes_II(:,1:4), iLevel_I(1:4), &
            nGridOut, Weight_I(1:4), DoStencilFix, XyzStencil2_D)
       return

    elseif(all(iLevel_I(1:4)== BehindTheBoundary_)                                    )then

       Xyz2_D(x_:y_) = Xyz_D(x_:y_) ; XyzGrid2_DI(x_:y_,:) = XyzGrid_DI(x_:y_,5:8)
       iIndexes_II(:,1:4) = iIndexes_II(:,5:8)
       iLevel_I(     1:4) = iLevel_I(     5:8)

       call interpolate_on_amr_grid_2(&
            nIndexes, Xyz2_D,  XyzGrid2_DI, iIndexes_II(:,1:4), iLevel_I(1:4), &
            nGridOut, Weight_I(1:4), DoStencilFix, XyzStencil2_D )
       return

    elseif(all(iLevel_I(2:8:2)== BehindTheBoundary_).or. all(XyzGrid_DI(x_,1:7:2)==Xyz_D(x_)) )then

       Xyz2_D(x_:y_) = Xyz_D(y_:z_) ; XyzGrid2_DI(x_:y_,:) = XyzGrid_DI(y_:z_,1:7:2)
       iIndexes_II(:,2) = iIndexes_II(:,3)
       iLevel_I(     2) = iLevel_I(     3)
       iIndexes_II(:,3) = iIndexes_II(:,5)
       iLevel_I(     3) = iLevel_I(     5)
       iIndexes_II(:,4) = iIndexes_II(:,7)
       iLevel_I(     4) = iLevel_I(     7)


       call interpolate_on_amr_grid_2(&
            nIndexes, Xyz2_D,  XyzGrid2_DI, iIndexes_II(:,1:4), iLevel_I(1:4), &
            nGridOut, Weight_I(1:4), DoStencilFix, XyzStencil2_D )
       return

    elseif(all(iLevel_I(1:7:2)== BehindTheBoundary_)                                        )then

       Xyz2_D(x_:y_) = Xyz_D(y_:z_) ; XyzGrid2_DI(x_:y_,:) = XyzGrid_DI(y_:z_,2:8:2)
       iIndexes_II(:,1) = iIndexes_II(:,2)
       iLevel_I(     1) = iLevel_I(     2)
       iIndexes_II(:,2) = iIndexes_II(:,4)
       iLevel_I(     2) = iLevel_I(     4)
       iIndexes_II(:,3) = iIndexes_II(:,6)
       iLevel_I(     3) = iLevel_I(     6)
       iIndexes_II(:,4) = iIndexes_II(:,8)
       iLevel_I(     4) = iLevel_I(     8)

    elseif( (all(iLevel_I(3:4)== BehindTheBoundary_).and.all(iLevel_I(7:8)== BehindTheBoundary_))&
         .or.( all(XyzGrid_DI(y_,1:2)==Xyz_D(y_)).and.all(XyzGrid_DI(y_,5:6)==Xyz_D(y_))) )then

       Xyz2_D(x_) = Xyz_D(x_) ; Xyz2_D(y_) =  Xyz_D(z_)
       XyzGrid2_DI(x_,1:2) = XyzGrid_DI(x_,1:2)
       XyzGrid2_DI(y_,1:2) = XyzGrid_DI(z_,1:2)
       XyzGrid2_DI(x_,3:4) = XyzGrid_DI(x_,5:6)
       XyzGrid2_DI(y_,3:4) = XyzGrid_DI(z_,5:6)
       iIndexes_II(:,3:4) = iIndexes_II(:,5:6)
       iLevel_I(     3:4) = iLevel_I(     5:6)
       

       call interpolate_on_amr_grid_2(&
            nIndexes, Xyz2_D,  XyzGrid2_DI, iIndexes_II(:,1:4), iLevel_I(1:4), &
            nGridOut, Weight_I(1:4), DoStencilFix, XyzStencil2_D )
       return

    elseif( all(iLevel_I(1:2)== BehindTheBoundary_).and.all(iLevel_I(5:6)== BehindTheBoundary_)&
                                                                                   )then

       Xyz2_D(x_) = Xyz_D(x_) ; Xyz2_D(y_) =  Xyz_D(z_)
       XyzGrid2_DI(x_,1:2) = XyzGrid_DI(x_,3:4)
       XyzGrid2_DI(y_,1:2) = XyzGrid_DI(z_,3:4)
       XyzGrid2_DI(x_,3:4) = XyzGrid_DI(x_,7:8)
       XyzGrid2_DI(y_,3:4) = XyzGrid_DI(z_,7:8)
       iIndexes_II(: ,1:2) = iIndexes_II(:,3:4)
       iLevel_I(      1:2) = iLevel_I(     3:4)
       iIndexes_II(: ,3:4) = iIndexes_II(:,7:8)
       iLevel_I(      3:4) = iLevel_I(     7:8)
       

       call interpolate_on_amr_grid_2(&
            nIndexes, Xyz2_D,  XyzGrid2_DI, iIndexes_II(:,1:4), iLevel_I(1:4), &
            nGridOut, Weight_I(1:4), DoStencilFix, XyzStencil2_D )
       return

    end if

    !\
    ! Make iLevel=0 for coarser grids and iLevel=1 for finer grids
    !/
    iLevelMin = minval(iLevel_I, MASK=iLevel_I/=BehindTheBoundary_)
    where(iLevel_I/=BehindTheBoundary_)iLevel_I = iLevel_I - iLevelMin

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
  
  end subroutine interpolate_on_amr_grid_3
  !=========
  !===========================TESTS======================================================
  subroutine test_interpolate_on_amr_grid_2(nSample)
    !The test subtoutine nSample times generate random point in the AMR domain
    !and compares the result of the coordinate interpolation to this point with the 
    !gampled coordinates
    ! 
    integer, intent(in) :: nSample

    integer :: iSeed=0, i, j, iBlock, loc(3), iSample, iGrid 
    integer, parameter :: nDim = 2, nBlock = 44, nX=2, nY=2, nIndexes=3 
    real    :: Xyz_DCB(nDim, nX, nY, nBlock)  !grid coordinates
    ! Block pattern
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
    ! Parameters to call interpolate_on_amr_grid_2
    !/
    real :: Xyz_D(nDim)=0.0, XyzStencil_D(nDim)=0.0, XyzGrid_DI(nDim,4)=0.0
    real :: Weight_I(4), XyzInterpolated_D(nDim)
    integer:: iLevel_I(4), iIndexes_II(nIndexes,4), nGridOut
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

    do iSample = 1, nSample
       Xyz_D(1) = 10.0*RAND()
       Xyz_D(2) =  8.0*RAND()
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
       call interpolate_on_amr_grid_2(nIndexes, Xyz_D, XyzGrid_DI, iIndexes_II, iLevel_I,&
                  nGridOut, Weight_I, DoStencilFix, XyzStencil_D)
       if(DoStencilFix)then
          call generate_basic_stencil(XyzStencil_D)
          call interpolate_on_amr_grid_2(nIndexes, Xyz_D, XyzGrid_DI, iIndexes_II, iLevel_I,&
                  nGridOut, Weight_I, DoStencilFix, XyzStencil_D)
       end if
       XyzInterpolated_D=0.0
       do iGrid = 1, nGridOut
          i = iIndexes_II(1,iGrid); j = iIndexes_II(2,iGrid); iBlock = iIndexes_II(3,iGrid)
          XyzInterpolated_D = XyzInterpolated_D + &
               Weight_I(iGrid) * Xyz_DCB(:, i, j, iBlock)
       end do
       !\
       ! Check interpolated value and compare with Xyz_D
       !/
       if(sum(abs(Xyz_D(:)-XyzInterpolated_D(:))) > 1.0e-5)then
          write(*,*)'Test failed at iSample = ', iSample
          write(*,*)'Sampled Xyz_D=', Xyz_D
          write(*,*)'Interpolated value=',XyzInterpolated_D
          write(*,*)'Stencil: i j iBlock Weight'
          do iGrid = 1, nGridOut
             write(*,*)'Point #',iGrid, iIndexes_II(:,iGrid), Weight_I(iGrid)
          end do
          write(*,*)'Step by step generate stencil'
          DoStencilFix = .false.
          call generate_basic_stencil(Xyz_D)
          write(*,*)'Stencil: i j iBlock Level XyzGrid_D'
          do iGrid = 1, 4
             write(*,*)'Point #',iGrid, iIndexes_II(:,iGrid), iLevel_I(iGrid), XyzGrid_DI(:,iGrid)
          end do
          call interpolate_on_amr_grid_2(nIndexes, Xyz_D, XyzGrid_DI, iIndexes_II, iLevel_I,&
                  nGridOut, Weight_I, DoStencilFix, XyzStencil_D)
          write(*,*)'DoStencilFix=',DoStencilFix, ' nGridOut=', nGridOut
          call CON_stop('Test failed')
       end if
    end do
                               
  contains
    subroutine generate_basic_stencil(XyzIn_D)
      real, intent(in):: XyzIn_D(nDim)
      !-------------------------------
      !gird point 1
      Loc = minloc((XyzIn_D(2) - Xyz_DCB(2,:,:,:))**2+&
           (XyzIn_D(1) - Xyz_DCB(1,:,:,:))**2, &
           MASK=Xyz_DCB(1,:,:,:).le.XyzIn_D(1).and.&
                Xyz_DCB(2,:,:,:).le.XyzIn_D(2))
      iIndexes_II(:,1) = loc(:)
      iLevel_I(   1) = iLevel_B(loc(3))
      XyzGrid_DI(:,1)= Xyz_DCB(:,loc(1),loc(2),loc(3))
      !gird point 2
      Loc = minloc((XyzIn_D(2) - Xyz_DCB(2,:,:,:))**2+&
           (XyzIn_D(1) - Xyz_DCB(1,:,:,:))**2, &
           MASK=Xyz_DCB(1,:,:,:) > XyzIn_D(1).and.&
                Xyz_DCB(2,:,:,:).le.XyzIn_D(2))
      iIndexes_II(:,2) = loc(:)
      iLevel_I(   2) = iLevel_B(loc(3))
      XyzGrid_DI(:,2)= Xyz_DCB(:,loc(1),loc(2),loc(3))
      !gird point 3
      Loc = minloc((XyzIn_D(2) - Xyz_DCB(2,:,:,:))**2+&
           (XyzIn_D(1) - Xyz_DCB(1,:,:,:))**2, &
           MASK=Xyz_DCB(1,:,:,:).le.XyzIn_D(1).and.&
                Xyz_DCB(2,:,:,:)>XyzIn_D(2))
      iIndexes_II(:,3) = loc(:)
      iLevel_I(   3) = iLevel_B(loc(3))
      XyzGrid_DI(:,3)= Xyz_DCB(:,loc(1),loc(2),loc(3))
      !gird point 4
      Loc = minloc((XyzIn_D(2) - Xyz_DCB(2,:,:,:))**2+&
           (XyzIn_D(1) - Xyz_DCB(1,:,:,:))**2, &
           MASK=Xyz_DCB(1,:,:,:) > XyzIn_D(1).and.&
                Xyz_DCB(2,:,:,:) >XyzIn_D(2))
      iIndexes_II(:,4) = loc(:)
      iLevel_I(   4) = iLevel_B(loc(3))
      XyzGrid_DI(:,4)= Xyz_DCB(:,loc(1),loc(2),loc(3))
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
    !=====================
  end subroutine test_interpolate_on_amr_grid_2
end module ModInterpolateAMRGrid
