!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModInterpolateBody
  !\
  ! Construct interpolation weights for a point inside
  ! bodies, interolation being in terms of the values 
  ! in the vertexes of the body. Examples of the implemented
  ! bodies: tetrahedron, pyramids (rectangular and trapezoidal)
  ! and their different compositions. Interpolation in parallel 
  ! rays implements the following procedure; a ray of given 
  ! direction is passed trough a given point inside the body, 
  ! in the point of the ray intersects the body boundary, the value
  ! is interpolated over planar subface (traingle, trapezoid, 
  ! ractangle), then the point value is weighted in accordance with
  ! the distances from the given point to its intersections with the 
  ! body surface. Resolution_corner interpolates within the resolution 
  ! corner, which has "fine vertexes"  - grid points of finer 
  ! rectangular grid and "coarse vertexes" - grid points of coarser
  ! rectangular grid with twice larger grid size. Finer grid cube and 
  ! coarser grid cube have comon central point. Such grid coniguration
  ! is typical for the resolution corner of a block adaptive grid, 
  ! so that the provided routine may be used for interpolating values
  ! in the proximity of such points.  
  !/
  implicit none
  PRIVATE
  SAVE
  integer, parameter:: nDim = 3, nGrid = 8
  integer, parameter:: &
       x_       = 1,               &
       y_       = 2,               &
       z_       = 3
  
  integer, parameter:: Rectangular_=1, Trapezoidal_=2
  public:: pyramids, parallel_rays, resolution_corner
contains
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
  subroutine tetrahedron(X1_D, X2_D, X3_D, X4_D, Xyz_D, Weight_I)
    !\
    ! Interpolate in the tetrahedron
    !/
    real, dimension(nDim), intent(in) :: X1_D, X2_D, X3_D, X4_D, Xyz_D
    real, intent(out):: Weight_I(1:4)
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
  subroutine pyramid(iBase,X1_D, X2_D, X3_D, X4_D, X5_D, Xyz_D, Weight_I)
    !\
    ! Interpolate in the pyramid with base X1X2X3X4 and apex X5
    ! valid for case of rectangular or trapezoid base
    !/
    integer, intent(in) :: iBase
    real, dimension(nDim), intent(in) :: X1_D, X2_D, X3_D, X4_D, X5_D, Xyz_D
    real, intent(out):: Weight_I(1:5)
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
       call rectangle(X1_D, X2_D, X3_D, X4_D, XyzP_D, Weight_I(1:4))
    case(Trapezoidal_)
       call trapezoid(X1_D, X2_D, X3_D, X4_D, XyzP_D, Weight_I(1:4))
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
  subroutine rectangle(X1_D, X2_D, X3_D, X4_D, XyzP_D, Weight_I)
    real, dimension(nDim), intent(in):: X1_D, X2_D, X3_D, X4_D, XyzP_D
    real, intent(out):: Weight_I(1:4)
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
  subroutine trapezoid(X1_D, X2_D, X3_D, X4_D, XyzP_D, Weight_I)
    real, dimension(nDim), intent(in):: X1_D, X2_D, X3_D, X4_D, XyzP_D
    real, intent(out)::Weight_I(1:4)
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
  subroutine triangle(X1_D, X2_D, X3_D,  XyzP_D, Weight_I)
    real, dimension(nDim), intent(in):: X1_D, X2_D, X3_D, XyzP_D
    real, intent(out):: Weight_I(1:3)
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
  !=====================================
  subroutine pyramids(&
       XyzGrid_DI,Xyz_D,iOrder_I, Weight_I, nGridOut, &
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
    real, intent(in):: XyzGrid_DI(nDim,nGrid), Xyz_D(nDim)
    real, intent(out):: Weight_I(nGrid)
    integer, intent(out)::nGridOut, iOrder_I(1:5)
    !---------------------------
    nGridOut = -1; Weight_I = -1
    if(present(iTetrahedron1_I))then
       call tetrahedron(&
            XyzGrid_DI(:,iTetrahedron1_I(1)),&
            XyzGrid_DI(:,iTetrahedron1_I(2)),&
            XyzGrid_DI(:,iTetrahedron1_I(3)),&
            XyzGrid_DI(:,iTetrahedron1_I(4)),Xyz_D, Weight_I(1:4))
       if(all(Weight_I(1:4).ge.0.0))then
          nGridOut = 4
          iOrder_I(1:4) = iTetrahedron1_I
          RETURN
       elseif(present(iTetrahedron2_I))then
          call tetrahedron(&
               XyzGrid_DI(:,iTetrahedron2_I(1)),&
               XyzGrid_DI(:,iTetrahedron2_I(2)),&
               XyzGrid_DI(:,iTetrahedron2_I(3)),&
               XyzGrid_DI(:,iTetrahedron2_I(4)),Xyz_D, Weight_I(1:4))
          if(all(Weight_I(1:4).ge.0.0))then
             nGridOut = 4
             iOrder_I(1:4) = iTetrahedron2_I
             RETURN
          elseif(present(iTetrahedron3_I))then
             call tetrahedron(&
                  XyzGrid_DI(:,iTetrahedron3_I(1)),&
                  XyzGrid_DI(:,iTetrahedron3_I(2)),&
                  XyzGrid_DI(:,iTetrahedron3_I(3)),&
                  XyzGrid_DI(:,iTetrahedron3_I(4)),Xyz_D, Weight_I(1:4))
             if(all(Weight_I(1:4).ge.0.0))then
                nGridOut = 4
                iOrder_I(1:4) = iTetrahedron3_I
                RETURN
             elseif(present(iTetrahedron4_I))then
                call tetrahedron(&
                     XyzGrid_DI(:,iTetrahedron4_I(1)),&
                     XyzGrid_DI(:,iTetrahedron4_I(2)),&
                     XyzGrid_DI(:,iTetrahedron4_I(3)),&
                     XyzGrid_DI(:,iTetrahedron4_I(4)),Xyz_D, Weight_I(1:4))
                if(all(Weight_I(1:4).ge.0.0))then
                   nGridOut = 4
                   iOrder_I(1:4) = iTetrahedron4_I
                   RETURN
                elseif(present(iTetrahedron5_I))then
                   call tetrahedron(&
                        XyzGrid_DI(:,iTetrahedron5_I(1)),&
                        XyzGrid_DI(:,iTetrahedron5_I(2)),&
                        XyzGrid_DI(:,iTetrahedron5_I(3)),&
                        XyzGrid_DI(:,iTetrahedron5_I(4)),Xyz_D, Weight_I(1:4))
                   if(all(Weight_I(1:4).ge.0.0))then
                      nGridOut = 4
                      iOrder_I(1:4) = iTetrahedron5_I
                      RETURN
                   elseif(present(iTetrahedron6_I))then
                      call tetrahedron(&
                           XyzGrid_DI(:,iTetrahedron6_I(1)),&
                           XyzGrid_DI(:,iTetrahedron6_I(2)),&
                           XyzGrid_DI(:,iTetrahedron6_I(3)),&
                           XyzGrid_DI(:,iTetrahedron6_I(4)),Xyz_D, Weight_I(1:4))
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
       call pyramid(Rectangular_,&
            XyzGrid_DI(:,iRectangular1_I(1)),&
            XyzGrid_DI(:,iRectangular1_I(2)),&
            XyzGrid_DI(:,iRectangular1_I(3)),&
            XyzGrid_DI(:,iRectangular1_I(4)),&
            XyzGrid_DI(:,iRectangular1_I(5)),Xyz_D, Weight_I(1:5) )
       if(all(Weight_I(1:5).ge.0.0))then
          nGridOut = 5
          iOrder_I(1:5) = iRectangular1_I
          RETURN
       elseif(present(iRectangular2_I))then
          call pyramid(Rectangular_,&
               XyzGrid_DI(:,iRectangular2_I(1)),&
               XyzGrid_DI(:,iRectangular2_I(2)),&
               XyzGrid_DI(:,iRectangular2_I(3)),&
               XyzGrid_DI(:,iRectangular2_I(4)),&
               XyzGrid_DI(:,iRectangular2_I(5)),Xyz_D, Weight_I(1:5) )
          if(all(Weight_I(1:5).ge.0.0))then
             nGridOut = 5
             iOrder_I(1:5) = iRectangular2_I
             RETURN
          elseif(present(iRectangular3_I))then
             call pyramid(Rectangular_,&
                  XyzGrid_DI(:,iRectangular3_I(1)),&
                  XyzGrid_DI(:,iRectangular3_I(2)),&
                  XyzGrid_DI(:,iRectangular3_I(3)),&
                  XyzGrid_DI(:,iRectangular3_I(4)),&
                  XyzGrid_DI(:,iRectangular3_I(5)),Xyz_D, Weight_I(1:5) )
             if(all(Weight_I(1:5).ge.0.0))then
                nGridOut = 5
                iOrder_I(1:5) = iRectangular3_I
                RETURN
             end if          ! 3
          end if             ! 2
       end if                ! 1
    end if                   ! no rectangular
    if(present(iTrapezoidal1_I))then
       call pyramid(Trapezoidal_,&
            XyzGrid_DI(:,iTrapezoidal1_I(1)),&
            XyzGrid_DI(:,iTrapezoidal1_I(2)),&
            XyzGrid_DI(:,iTrapezoidal1_I(3)),&
            XyzGrid_DI(:,iTrapezoidal1_I(4)),&
            XyzGrid_DI(:,iTrapezoidal1_I(5)),Xyz_D, Weight_I(1:5) )
       if(all(Weight_I(1:5).ge.0.0))then
          nGridOut = 5
          iOrder_I(1:5) = iTrapezoidal1_I
          RETURN
       elseif(present(iTrapezoidal2_I))then
          call pyramid(Trapezoidal_,&
               XyzGrid_DI(:,iTrapezoidal2_I(1)),&
               XyzGrid_DI(:,iTrapezoidal2_I(2)),&
               XyzGrid_DI(:,iTrapezoidal2_I(3)),&
               XyzGrid_DI(:,iTrapezoidal2_I(4)),&
               XyzGrid_DI(:,iTrapezoidal2_I(5)),Xyz_D, Weight_I(1:5) )
          if(all(Weight_I(1:5).ge.0.0))then
             nGridOut = 5
             iOrder_I(1:5) = iTrapezoidal2_I
             RETURN
          end if             ! 2
       end if                ! 1
    end if                   ! no trapezoid
  end subroutine pyramids
  !======================
  subroutine parallel_rays(&
       XyzGrid_DI, Xyz_D, iOrder_I, Weight_I, nGridOut, Dir_D, &
       iDRectangle1_I, iDTrapezoid1_I,   &
       iDTriangle1_I, iDTriangle2_I, iDTriangle3_I,&
       iDTriangle4_I, iDTriangle5_I, iDTriangle6_I,&
       iURectangle1_I, iURectangle2_I, iURectangle3_I,&
       iUTrapezoid1_I,                             &
       iUTriangle1_I, iUTriangle2_I, iUTriangle3_I,&
       iUTriangle4_I)
    real, intent(in):: XyzGrid_DI(nDim,nGrid), Xyz_D(nDim)
    real, intent(out):: Weight_I(nGrid)
    integer, intent(out)::nGridOut, iOrder_I(8)
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
    ! the weights for the upper face should be multiplied by 
    ! AlphaDown/(AlphaUp+AlphaDown)
    ! for the down face - by AlphaUp/(AlphaUp+AlphaDown)
    !/

    nGridOut = -1; nGridOutUp= -1; nGridOutDown = -1
    Weight_I = -1; iOrder_I = 0
    if(present(iUTriangle1_I))then
       X1_D = XyzGrid_DI(:,iUTriangle1_I(1))
       X2_D = XyzGrid_DI(:,iUTriangle1_I(2))
       X3_D = XyzGrid_DI(:,iUTriangle1_I(3))
       AlphaUp = triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
            triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
       XyzUp_D = Xyz_D + AlphaUp * Dir_D
       call triangle(X1_D, X2_D, X3_D, XyzUp_D, Weight_I(1:3))
       if(all(Weight_I(1:3)>=0.0))then
          !\
          !The ray for point Xyz is projected into this triangle
          !/
          if(AlphaUp==0.0)then
             !\
             !Point Xyz belongs to this traingle
             !/
             iOrder_I(1:3) = iUTriangle1_I
             nGridOut = 3
             RETURN
          elseif(AlphaUp < 0)then
             !\
             ! Point is above the upper triangle
             ! return with negative weights and nGridOut = -1
             RETURN 
          else
             iOrder_I(1:3) = iUTriangle1_I
             nGridOutUp = 3
             !\
             ! Exit if for triangles with positive AlphaUp and nGridOut = 3
             !/
          end if
       elseif(present(iUTriangle2_I))then
          X1_D = XyzGrid_DI(:,iUTriangle2_I(1))
          X2_D = XyzGrid_DI(:,iUTriangle2_I(2))
          X3_D = XyzGrid_DI(:,iUTriangle2_I(3))
          AlphaUp = triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
               triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
          XyzUp_D = Xyz_D + AlphaUp * Dir_D
          call triangle(X1_D, X2_D, X3_D, XyzUp_D, Weight_I(1:3))
          if(all(Weight_I(1:3)>=0.0))then
             !\
             !The ray for point Xyz is projected into this triangle
             !/
             if(AlphaUp==0.0)then
                !\
                !Point Xyz belongs to this traingle
                !/
                iOrder_I(1:3) = iUTriangle2_I
                nGridOut = 3
                RETURN
             elseif(AlphaUp < 0)then
                !\
                ! Point is above the upper triangle
                ! return with negative weights and nGridOut = -1
                RETURN 
             else             
                iOrder_I(1:3) = iUTriangle2_I
                nGridOutUp = 3
                !\
                ! Exit if for triangles with positive AlphaUp
                !/
             end if
          elseif(present(iUTriangle3_I))then
             X1_D = XyzGrid_DI(:,iUTriangle3_I(1))
             X2_D = XyzGrid_DI(:,iUTriangle3_I(2))
             X3_D = XyzGrid_DI(:,iUTriangle3_I(3))
             AlphaUp = &
                  triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
                  triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
             XyzUp_D = Xyz_D + AlphaUp * Dir_D
             call triangle(X1_D, X2_D, X3_D, XyzUp_D, Weight_I(1:3))
             if(all(Weight_I(1:3)>=0.0))then
                !\
                !The ray for point Xyz is projected into this triangle
                !/
                if(AlphaUp==0.0)then
                   !\
                   !Point Xyz belongs to this traingle
                   !/
                   iOrder_I(1:3) = iUTriangle3_I
                   nGridOut = 3
                   RETURN
                elseif(AlphaUp < 0)then
                   !\
                   ! Point is above the upper triangle
                   ! return with negative weights and nGridOut = -1
                   RETURN 
                else        
                   iOrder_I(1:3) = iUTriangle3_I
                   nGridOutUp = 3
                   !\
                   ! Exit if for triangles with positive AlphaUp 
                   !/

                end if
             elseif(present(iUTriangle4_I))then
                X1_D = XyzGrid_DI(:,iUTriangle4_I(1))
                X2_D = XyzGrid_DI(:,iUTriangle4_I(2))
                X3_D = XyzGrid_DI(:,iUTriangle4_I(3))
                AlphaUp = &
                     triple_product(X1_D - Xyz_D,X2_D - X1_D,X3_D - X1_D)/&
                     triple_product(Dir_D       ,X2_D - X1_D,X3_D - X1_D)
                XyzUp_D = Xyz_D + AlphaUp * Dir_D
                call triangle(X1_D, X2_D, X3_D, Xyz_D, Weight_I(1:3))
                if(all(Weight_I(1:3)>=0.0))then
                   !\
                   !The ray for point Xyz is projected into this triangle
                   !/
                   if(AlphaUp==0.0)then
                      !\
                      !Point Xyz belongs to this traingle
                      !/
                      iOrder_I(1:3) = iUTriangle4_I
                      nGridOut = 3
                      RETURN
                   elseif(AlphaUp < 0.0)then
                      !\
                      ! Point is above the upper triangle
                      ! return with negative weights and nGridOut = -1
                      RETURN 
                   else   
                      iOrder_I(1:3) = iUTriangle4_I
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
          X1_D = XyzGrid_DI(:,iURectangle1_I(1))
          X2_D = XyzGrid_DI(:,iURectangle1_I(2))
          X3_D = XyzGrid_DI(:,iURectangle1_I(3))
          X4_D = XyzGrid_DI(:,iURectangle1_I(4))
          AlphaUp = triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
               triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
          XyzUp_D = Xyz_D + AlphaUp * Dir_D
          call rectangle(X1_D, X2_D, X3_D, X4_D, XyzUp_D, Weight_I(1:4))
          if(all(Weight_I(1:4)>=0.0))then
             !\
             !The ray for point Xyz is projected into this rectangle
             !/
             if(AlphaUp==0.0)then
                !\
                !Point Xyz belongs to this rectangle
                !/
                iOrder_I(1:4) = iURectangle1_I
                nGridOut = 4
                RETURN
             elseif(AlphaUp < 0.0)then
                !\
                ! Point is above the upper rectangle
                ! return with negative weights and nGridOut = -1
                RETURN 
             else 
                iOrder_I(1:4) = iURectangle1_I
                nGridOutUp = 4
                !\
                ! Exit if for rectangles with positive AlphaUp 
                !/
             end if
          elseif(present(iURectangle2_I))then
             X1_D = XyzGrid_DI(:,iURectangle2_I(1))
             X2_D = XyzGrid_DI(:,iURectangle2_I(2))
             X3_D = XyzGrid_DI(:,iURectangle2_I(3))
             X4_D = XyzGrid_DI(:,iURectangle2_I(4))
             AlphaUp = triple_product(X1_D - Xyz_D,X2_D - X1_D,X3_D - X1_D)/&
                  triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
             XyzUp_D = Xyz_D + AlphaUp * Dir_D
             call rectangle(X1_D, X2_D, X3_D, X4_D, XyzUp_D, Weight_I(1:4))
             if(all(Weight_I(1:4)>=0.0))then
                !\
                !The ray for point Xyz is projected into this rectangle
                !/
                if(AlphaUp==0.0)then
                   !\
                   !Point Xyz belongs to this rectangle
                   !/
                   iOrder_I(1:4) = iURectangle2_I
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
                   iOrder_I(1:4) = iURectangle2_I
                   nGridOutUp = 4
                end if
             elseif(present(iURectangle3_I))then
                X1_D = XyzGrid_DI(:,iURectangle3_I(1))
                X2_D = XyzGrid_DI(:,iURectangle3_I(2))
                X3_D = XyzGrid_DI(:,iURectangle3_I(3))
                X4_D = XyzGrid_DI(:,iURectangle3_I(4))
                AlphaUp = &
                     triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
                     triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
                XyzUp_D = Xyz_D + AlphaUp * Dir_D
                call rectangle(X1_D, X2_D, X3_D, X4_D, XyzUp_D, Weight_I(1:4))
                if(all(Weight_I(1:4)>=0.0))then
                   !\
                   !The ray for point Xyz is projected into this rectangle
                   !/
                   if(AlphaUp==0.0)then
                      !\
                      !Point Xyz belongs to this rectangle
                      !/
                      iOrder_I(1:4) = iURectangle3_I
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
                      iOrder_I(1:4) = iURectangle3_I
                      nGridOutUp = 4
                   end if
                end if
             end if     !3
          end if        !2
       end if           !1
    end if
    if(nGridOutUp < 1)then
       if(present(iUTrapezoid1_I))then
          X1_D = XyzGrid_DI(:,iUTrapezoid1_I(1))
          X2_D = XyzGrid_DI(:,iUTrapezoid1_I(2))
          X3_D = XyzGrid_DI(:,iUTrapezoid1_I(3))
          X4_D = XyzGrid_DI(:,iUTrapezoid1_I(4))
          AlphaUp = triple_product(X1_D - Xyz_D, X2_D - X1_D, X3_D - X1_D)/&
               triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
          XyzUp_D = Xyz_D + AlphaUp * Dir_D
          call Trapezoid(X1_D, X2_D, X3_D, X4_D, XyzUp_D, Weight_I(1:4))
          if(all(Weight_I(1:4)>=0.0))then
             !\
             !The ray for point Xyz is projected into this trapezoid
             !/
             if(AlphaUp==0.0)then
                !\
                !Point Xyz belongs to this rectangle
                !/
                iOrder_I(1:4) = iUTrapezoid1_I
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
                iOrder_I(1:4) = iUTrapezoid1_I
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
       iOrder_I(4:3+nGridOutUp) = iOrder_I(1:nGridOutUp)
       iOrder_I(1:3) = 0
       X1_D = XyzGrid_DI(:,iDTriangle1_I(1))
       X2_D = XyzGrid_DI(:,iDTriangle1_I(2))
       X3_D = XyzGrid_DI(:,iDTriangle1_I(3))
       AlphaDown = &
            triple_product(Xyz_D - X1_D, X2_D - X1_D, X3_D - X1_D)/&
            triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
       XyzDown_D = Xyz_D - AlphaDown*Dir_D
       call triangle(X1_D, X2_D, X3_D, XyzDown_D, Weight_I(1:3))
       if(all(Weight_I(1:3)>=0.0))then
          !\
          !The ray for point Xyz is projected into this triangle
          !/
          if(AlphaDown==0.0)then
             !\
             !Point Xyz belongs to this triangle
             !/
             iOrder_I(1:3) = iDTriangle1_I
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
             iOrder_I(1:3) = iDTriangle1_I
             nGridOutDown = 3
          end if
       elseif(present(iDTriangle2_I))then
          X1_D = XyzGrid_DI(:,iDTriangle2_I(1))
          X2_D = XyzGrid_DI(:,iDTriangle2_I(2))
          X3_D = XyzGrid_DI(:,iDTriangle2_I(3))
          AlphaDown = &
               triple_product(Xyz_D - X1_D, X2_D - X1_D, X3_D - X1_D)/&
               triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
          XyzDown_D = Xyz_D - AlphaDown * Dir_D
          call triangle(X1_D, X2_D, X3_D, XyzDown_D, Weight_I(1:3))
          if(all(Weight_I(1:3)>=0.0))then
             !\
             !The ray for point Xyz is projected into this triangle
             !/
             if(AlphaDown==0.0)then
                !\
                !Point Xyz belongs to this triangle
                !/
                iOrder_I(1:3) = iDTriangle2_I
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
                iOrder_I(1:3) = iDTriangle2_I
                nGridOutDown = 3
             end if
          elseif(present(iDTriangle3_I))then
             X1_D = XyzGrid_DI(:,iDTriangle3_I(1))
             X2_D = XyzGrid_DI(:,iDTriangle3_I(2))
             X3_D = XyzGrid_DI(:,iDTriangle3_I(3))
             AlphaDown = &
                  triple_product(Xyz_D - X1_D,X2_D - X1_D,X3_D - X1_D)/&
                  triple_product(Dir_D       ,X2_D - X1_D,X3_D - X1_D)
             XyzDown_D = Xyz_D - AlphaDown * Dir_D
             call triangle(X1_D, X2_D, X3_D, XyzDown_D, Weight_I(1:3))
             if(all(Weight_I(1:3)>=0.0))then
                !\
                !The ray for point Xyz is projected into this triangle
                !/
                if(AlphaDown==0.0)then
                   !\
                   !Point Xyz belongs to this triangle
                   !/
                   iOrder_I(1:3) = iDTriangle3_I
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
                   iOrder_I(1:3) = iDTriangle3_I
                   nGridOutDown = 3
                end if
             elseif(present(iDTriangle4_I))then
                X1_D = XyzGrid_DI(:,iDTriangle4_I(1))
                X2_D = XyzGrid_DI(:,iDTriangle4_I(2))
                X3_D = XyzGrid_DI(:,iDTriangle4_I(3))
                AlphaDown = &
                     triple_product(Xyz_D - X1_D,X2_D - X1_D,X3_D - X1_D)/&
                     triple_product(Dir_D       ,X2_D - X1_D,X3_D - X1_D)
                XyzDown_D = Xyz_D - AlphaDown * Dir_D
                call triangle(X1_D, X2_D, X3_D, XyzDown_D, Weight_I(1:3))
                if(all(Weight_I(1:3)>=0.0))then
                   !\
                   !The ray for point Xyz is projected into this triangle
                   !/
                   if(AlphaDown==0.0)then
                      !\
                      !Point Xyz belongs to this triangle
                      !/
                      iOrder_I(1:3) = iDTriangle4_I
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
                      iOrder_I(1:3) = iDTriangle4_I
                      nGridOutDown = 3
                   end if
                elseif(present(iDTriangle5_I))then
                   X1_D = XyzGrid_DI(:,iDTriangle5_I(1))
                   X2_D = XyzGrid_DI(:,iDTriangle5_I(2))
                   X3_D = XyzGrid_DI(:,iDTriangle5_I(3))
                   AlphaDown = &
                        triple_product(&
                        Xyz_D - X1_D,X2_D - X1_D,X3_D - X1_D)/&
                        triple_product(&
                        Dir_D       ,X2_D - X1_D,X3_D - X1_D)
                   XyzDown_D = Xyz_D - AlphaDown * Dir_D
                   call triangle(X1_D, X2_D, X3_D, XyzDown_D, Weight_I(1:3))
                   if(all(Weight_I(1:3)>=0.0))then
                      !\
                      !The ray for point Xyz is projected into this triangle
                      !/
                      if(AlphaDown==0.0)then
                         !\
                         !Point Xyz belongs to this triangle
                         !/
                         iOrder_I(1:3) = iDTriangle5_I
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
                         iOrder_I(1:3) = iDTriangle5_I
                         nGridOutDown = 3
                      end if
                   elseif(present(iDTriangle6_I))then
                      X1_D = XyzGrid_DI(:,iDTriangle6_I(1))
                      X2_D = XyzGrid_DI(:,iDTriangle6_I(2))
                      X3_D = XyzGrid_DI(:,iDTriangle6_I(3))
                      AlphaDown = &
                           triple_product(&
                           Xyz_D - X1_D,X2_D - X1_D,X3_D - X1_D)/&
                           triple_product(&
                           Dir_D       ,X2_D - X1_D,X3_D - X1_D)
                      XyzDown_D = Xyz_D - AlphaDown * Dir_D
                      call triangle(X1_D, X2_D, X3_D, XyzDown_D, Weight_I(1:3))
                      if(all(Weight_I(1:3)>=0.0))then
                         !\
                         !The ray for point Xyz is projected into this 
                         !triangle
                         !/
                         if(AlphaDown==0.0)then
                            !\
                            !Point Xyz belongs to this triangle
                            !/
                            iOrder_I(1:3) = iDTriangle6_I
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
                            iOrder_I(1:3) = iDTriangle6_I
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
             iOrder_I(5:4+nGridOutUp) = iOrder_I(4:3+nGridOutUp)
             iOrder_I(1:4)=0
          else
             Weight_I(5:4+nGridOutUp) = Weight_I(1:nGridOutUp)
             Weight_I(1:4) = 0
             iOrder_I(5:4+nGridOutUp) = iOrder_I(1:nGridOutUp)
             iOrder_I(1:4)=0
          end if
          X1_D = XyzGrid_DI(:,iDRectangle1_I(1))
          X2_D = XyzGrid_DI(:,iDRectangle1_I(2))
          X3_D = XyzGrid_DI(:,iDRectangle1_I(3))
          X4_D = XyzGrid_DI(:,iDRectangle1_I(4))
          AlphaDown = triple_product(Xyz_D - X1_D,X2_D - X1_D,X3_D - X1_D)/&
               triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
          XyzDown_D = Xyz_D - AlphaDown * Dir_D
          call rectangle(X1_D, X2_D, X3_D, X4_D, XyzDown_D, Weight_I(1:4))
          if(all(Weight_I(1:4)>=0.0))then
             !\
             !The ray for point Xyz is projected into this rectangle
             !/
             if(AlphaDown==0.0)then
                !\
                !Point Xyz belongs to this rectangle
                !/
                iOrder_I(1:4) = iDRectangle1_I
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
                iOrder_I(1:4) = iDRectangle1_I
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
                iOrder_I(5:4+nGridOutUp) = iOrder_I(4:3+nGridOutUp)
                iOrder_I(1:4) = 0
             else
                Weight_I(5:4+nGridOutUp) = Weight_I(1:nGridOutUp)
                Weight_I(1:4)=0
                iOrder_I(5:4+nGridOutUp) = iOrder_I(1:nGridOutUp)
                iOrder_I(1:4) = 0
             end if
          end if
          X1_D = XyzGrid_DI(:,iDTrapezoid1_I(1))
          X2_D = XyzGrid_DI(:,iDTrapezoid1_I(2))
          X3_D = XyzGrid_DI(:,iDTrapezoid1_I(3))
          X4_D = XyzGrid_DI(:,iDTrapezoid1_I(4))
          AlphaDown = triple_product(Xyz_D - X1_D,X2_D - X1_D,X3_D - X1_D)/&
               triple_product(Dir_D       , X2_D - X1_D, X3_D - X1_D)
          XyzDown_D = Xyz_D - AlphaDown * Dir_D
          call trapezoid(X1_D, X2_D, X3_D, X4_D, XyzDown_D, Weight_I(1:4))
          if(all(Weight_I(1:4)>=0.0))then
             !\
             !The ray for point Xyz is projected into this trapezoid
             !/
             if(AlphaDown==0.0)then
                !\
                !Point Xyz belongs to this rectangle
                !/
                iOrder_I(1:4) = iDTrapezoid1_I
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
                iOrder_I(1:4) = iDTrapezoid1_I
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
    Weight_I(1:nGridOutDown) = &
         Weight_I(1:nGridOutDown)*AlphaUp/(AlphaUp + AlphaDown)
    Weight_I(nGridOutDown+1:nGridOut) = &
         Weight_I(nGridOutDown+1:nGridOut)*AlphaDown/(AlphaUp + AlphaDown)
  end subroutine parallel_rays
  !===========================
  subroutine  resolution_corner(&
       Xyz_D, XyzGridIn_DI, iLevel_I, nGridOut, Weight_I, iOrder_I)
    use ModCubeGeometry, ONLY: iSortStencil3_II, iFace_IDI, iOppositeFace_IDI,&
         Case_, Grid_, Dir_, FiveTetrahedra_, & 
         OneFine_, OneCoarse_, FineMainDiag_, FineFaceDiag_,            &
         CoarseMainDiag_, CoarseFaceDiag_, FineEdgePlusOne_,            &
         ThreeFineOnFace_, CoarseEdgePlusOne_, ThreeCoarseOnFace_,      &
         ThreeCoarseOnFacePlusOne_, CoarseChain_
    character(LEN=*),parameter:: NameSub='resolution_corner'
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
    real  ,   intent(in) :: XyzGridIn_DI(nDim,nGrid)
    !\
    ! The same, but may be reordered, if needed 
    real:: XyzGrid_DI(nDim,nGrid)

    !The refinement level at each grid point. By one higher level
    ! of refinement assumes the cell size reduced by a factor of 0.5
    integer,  intent(in) :: iLevel_I(nGrid)
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

    !\
    ! To find the sort of corner stencil
    !/
    integer:: iCase, iDir, iGrid
    !-------------
    Weight_I = 0
    XyzGrid_DI = XyzGridIn_DI
    !\
    ! We store in advance the 'basic' grid point
    ! and orientation of all possible stencil configurations
    ! and now we extracted this information
    !/ 
    iCase = i_case(iLevel_I)
    iGrid = iSortStencil3_II(Grid_,iCase)
    iDir  = iSortStencil3_II(Dir_,iCase)
    iCase = iSortStencil3_II(Case_,iCase)
    iOrder_I(1:4) = iFace_IDI(:,iDir,iGrid)
    iOrder_I(5:8) = iOppositeFace_IDI(:,iDir,iGrid)
    XyzGrid_DI = XyzGrid_DI(:,iOrder_I)
    select case(iCase)
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

       call my_pyramids(&
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
       call my_pyramids(&
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
       call my_pyramids(iTetrahedron1_I=(/1,2,3,5/))
       if(nGridOut>1)RETURN
       call my_parallel_rays(&
            Dir_D=XyzGrid_DI(:,8) - XyzGrid_DI(:,1), &
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
       call my_pyramids(&
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
       call my_pyramids(&
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
            0.50*(XyzGrid_DI(:,1) + XyzGrid_DI(:,4))
       !\
       ! Center of the upper face
       !/
       XyzMax_D = 0.50*&
            (XyzGrid_DI(:,5)+XyzGrid_DI(:,8))
       call my_parallel_rays(Dir_D=XyzMax_D - XyzMin_D,&
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
       call my_pyramids(&
            iTetrahedron1_I=(/1, 2, 3, 5/),&
            iTetrahedron2_I=(/8, 7, 6, 4/))
       if(nGridOut>1)RETURN
       call my_parallel_rays(&
            Dir_D=XyzGrid_DI(:,8) - XyzGrid_DI(:,1), &
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
       call my_pyramids(&
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
               (XyzGrid_DI(:,1)+XyzGrid_DI(:,4))
          !\
          ! Center of the upper face
          !/
          XyzMax_D = 0.50*&
               (XyzGrid_DI(:,5) + XyzGrid_DI(:,8))
          call my_parallel_rays(Dir_D=XyzMax_D - XyzMin_D,&
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
       call my_pyramids(&
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

       call my_pyramids(&
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
       call my_parallel_rays(Dir_D=&
            XyzGrid_DI(:,8) - &
            XyzGrid_DI(:,4),&
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
       call my_pyramids(&
            iTetrahedron1_I=(/1, 2, 3, 5/),&
            iRectangular1_I=(/2, 3, 6, 7, 5/))
       if(nGridOut> 1)RETURN
       call my_parallel_rays(Dir_D=&
            XyzGrid_DI(:,3) - &
            XyzGrid_DI(:,2),&
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
       call my_pyramids(&
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
       call my_pyramids(&
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
       call my_pyramids(&
            iTrapezoidal1_I=(/1, 2, 3, 4, 6/),& 
            iTetrahedron1_I=(/8, 3, 4, 6/),&  !see #1 below 
            iTetrahedron2_I=(/1, 6, 3, 5/),&  !see #2 below
            iTrapezoidal2_I=(/5, 7, 6, 8, 3/))!see #3 below
       ! #1 The leftover is above subfaces 136 and 364
       ! #2 The leftover is above subfaces 136 and 368
       ! #3 The leftover is above subfaces 365 and 368
       !==========================
    end select
  contains
    subroutine my_pyramids(&
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
      integer:: iOrderHere_I(1:5)
      !---------------------------
      call pyramids(&
           XyzGrid_DI, Xyz_D,iOrderHere_I, Weight_I, nGridOut,&
           iTetrahedron1_I,iTetrahedron2_I,iTetrahedron3_I,&
           iTetrahedron4_I,iTetrahedron5_I,iTetrahedron6_I,&
           iRectangular1_I,iRectangular2_I,iRectangular3_I,&
           iTrapezoidal1_I,iTrapezoidal2_I)
      if(nGridOut > -1)iOrder_I(1:nGridOut) = &
           iOrder_I(iOrderHere_I(1:nGridOut))
    end subroutine my_pyramids
    !======================
    subroutine my_parallel_rays(Dir_D, &
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
      call parallel_rays(&
           XyzGrid_DI, Xyz_D, iOrderHere_I, Weight_I, nGridOut, Dir_D, &
           iDRectangle1_I, iDTrapezoid1_I,   &
           iDTriangle1_I, iDTriangle2_I, iDTriangle3_I,&
           iDTriangle4_I, iDTriangle5_I, iDTriangle6_I,&
           iURectangle1_I, iURectangle2_I, iURectangle3_I,&
           iUTrapezoid1_I,                             &
           iUTriangle1_I, iUTriangle2_I, iUTriangle3_I,&
           iUTriangle4_I)
      if(nGridOut > 0)iOrder_I(1:nGridOut) =  &
           iOrder_I(iOrderHere_I(1:nGridOut))
    end subroutine my_parallel_rays
    !===========================
  end subroutine resolution_corner
end module ModInterpolateBody
