module ModInterpolate

  ! Calculate second order accurate interpolation for 
  !
  ! - a uniform grid with normalized coordinates, or 
  ! - non-uniform grid with actual coordinates, or
  ! - any mixture of the two, i.e. only some of the coordinates are uniform
  !
  ! Normalized coordinates mean that the coordinates coincide with the 
  ! indexes at the grid points. For uniform grid this is a very fast algorithm.
  ! For non-uniform grid a binary search is needed. The coordinates are assumed
  ! to be either monotone increasing or monotone decreasing. 
  !
  ! One can interpolate both scalar and vector valued arrays.
  !
  ! If the coordinates are outside the allowed ranges and the DoExtrapolate
  ! argument is not present the code stops. If the DoExtrapolate argument
  ! is present and false, the last grid cell value is used. If DoExtrapolate
  ! is present and true, second order extrapolation is used.

  ! Examples of usage:
  !
  ! Cell based 2D uniform grid with ghost cells, scalar valued:
  !
  !     InterpolatedValue = bilinear(Value_II, 0, nI+1, 0, nJ+1, &
  !                         (/ (x - x0)/DeltaX, (y - y0)/DeltaY) /) )
  !
  ! Node based 2D grid with x(1)=y(1)=0.0, vector valued:
  !
  !     InterpolatedValue_V = bilinear(Value_VII, nVar, 1, nI, 1, nJ, &
  !                        (/ x/DeltaX+1, y/DeltaY+1 /) )
  !
  ! Nonuniform 3D grid with ghost cells, third coordinate is uniform, 
  ! scalar valued:
  !
  !     InterpolatedValue = trilinear(Value_III, -1, nI+2, -1, nJ+2, -1, nK+2,&
  !                       (/ x, y, (z - z0)/DeltaZ /), x_I, y_I)
  !

  implicit none

  private ! except

  public :: linear             ! 2nd order interpolation in 1D
  public :: bilinear           ! 2nd order interpolation in 2D
  public :: trilinear          ! 2nd order interpolation in 3D
  public :: find_cell          ! find cell in non-uniform grid
  public :: test_interpolation ! unit test

  character(len=*), parameter :: NameMod='ModInterpolate'

  interface linear
     module procedure linear_scalar, linear_vector
  end interface

  interface bilinear
     module procedure bilinear_scalar, bilinear_vector
  end interface

  interface trilinear
     module procedure trilinear_scalar, trilinear_vector
  end interface

contains

  !=========================================================================
  real function linear_scalar(A_I, iMin, iMax, x, x_I, DoExtrapolate)

    ! Calculate linear interpolation of A_I at position x
    ! Assume normalized coordinates unless x_I is present.
    ! If present x_I contains the coordinates in an increasing order.

    implicit none

    integer, intent(in) :: iMin, iMax
    real, intent(in)    :: A_I(iMin:iMax)
    real, intent(in)    :: x

    real,    intent(in), optional :: x_I(iMin:iMax)
    logical, intent(in), optional :: DoExtrapolate

    integer :: i1
    real    :: Dx1
    character (len=*), parameter :: NameSub=NameMod//'::linear_scalar'
    !--------------------------------------------------------------------------
    call find_cell(iMin, iMax, x, i1, Dx1, x_I, DoExtrapolate, &
         "Called from "//NameSub)

    ! Perform interpolation (or extrapolation)
    linear_scalar = (1.0 - Dx1)*A_I(i1) + Dx1*A_I(i1+1)

  end function linear_scalar

  !=========================================================================
  function linear_vector(A_VI, nVar, iMin, iMax, x, x_I, DoExtrapolate)

    ! Calculate linear interpolation of A_VI at position x
    ! Assume normalized coordinates unless x_I is present.
    ! If present x_I contains the coordinates in an increasing order.

    implicit none

    integer, intent(in) :: nVar, iMin, iMax
    real, intent(in)    :: A_VI(nVar, iMin:iMax)
    real, intent(in)    :: x

    real,    intent(in), optional :: x_I(iMin:iMax)
    logical, intent(in), optional :: DoExtrapolate

    ! return value
    real                :: linear_vector(nVar)

    integer :: i1
    real    :: Dx1
    character (len=*), parameter :: NameSub=NameMod//'::linear_vector'
    !--------------------------------------------------------------------------
    call find_cell(iMin, iMax, x, i1, Dx1, x_I, DoExtrapolate, &
         "Called from "//NameSub)

    ! Perform interpolation (or extrapolation) for multiple variables
    linear_vector = (1.0-Dx1)*A_VI(:,i1) + Dx1*A_VI(:,i1+1)

  end function linear_vector

  !=========================================================================
  real function bilinear_scalar( &
       A_II, iMin, iMax, jMin, jMax, Xy_D, x_I, y_I, DoExtrapolate)

    ! Calculate bilinear interpolation of A_II at position Xy_D
    ! Assume normalized coordinates unless x_I and/or y_I are present.
    ! If present x_I and y_I contain the coordinates in an increasing order.

    implicit none

    integer, intent(in) :: iMin, iMax, jMin, jMax
    real, intent(in)    :: A_II(iMin:iMax,jMin:jMax)
    real, intent(in)    :: Xy_D(2)

    real,    intent(in), optional :: x_I(iMin:iMax)
    real,    intent(in), optional :: y_I(jMin:jMax)
    logical, intent(in), optional :: DoExtrapolate

    integer :: i1, i2, j1, j2
    real    :: Dx1, Dx2, Dy1, Dy2
    character (len=*), parameter :: NameSub=NameMod//'::bilinear_scalar'
    !--------------------------------------------------------------------------
    call find_cell(iMin, iMax, Xy_D(1), i1, Dx1, x_I, DoExtrapolate, &
         "Called for coord1 from "//NameSub)

    call find_cell(jMin, jMax, Xy_D(2), j1, Dy1, y_I, DoExtrapolate, &
         "Called for coord2 from "//NameSub)

    ! Perform interpolation (or extrapolation)
    i2 = i1 + 1; Dx2 = 1.0 - Dx1
    j2 = j1 + 1; Dy2 = 1.0 - Dy1

    bilinear_scalar = Dy2*( Dx2*A_II(i1,j1)   &
         +                  Dx1*A_II(i2,j1))  &
         +            Dy1*( Dx2*A_II(i1,j2)   &
         +                  Dx1*A_II(i2,j2))

  end function bilinear_scalar

  !=========================================================================
  function bilinear_vector( &
       A_VII, nVar, iMin, iMax, jMin, jMax, Xy_D, x_I, y_I, DoExtrapolate)

    ! Calculate bilinear interpolation of A_VII at position Xy_D
    ! Assume normalized coordinates unless x_I and/or y_I are present.
    ! If present x_I and y_I contain the coordinates in an increasing order.

    implicit none

    integer, intent(in) :: nVar, iMin, iMax, jMin, jMax
    real, intent(in)    :: A_VII(nVar, iMin:iMax,jMin:jMax)
    real, intent(in)    :: Xy_D(2)

    real,    intent(in), optional :: x_I(iMin:iMax)
    real,    intent(in), optional :: y_I(jMin:jMax)
    logical, intent(in), optional :: DoExtrapolate

    ! return value
    real                :: bilinear_vector(nVar)

    integer :: i1, i2, j1, j2
    real :: Dx1, Dx2, Dy1, Dy2
    character (len=*), parameter :: NameSub=NameMod//'::bilinear_vector'
    !--------------------------------------------------------------------------
    call find_cell(iMin, iMax, Xy_D(1), i1, Dx1, x_I, DoExtrapolate, &
         "Called for coord1 from "//NameSub)

    call find_cell(jMin, jMax, Xy_D(2), j1, Dy1, y_I, DoExtrapolate, &
         "Called for coord2 from "//NameSub)

    ! Perform interpolation (or extrapolation) for multiple variables
    i2 = i1 + 1; Dx2 = 1.0 - Dx1
    j2 = j1 + 1; Dy2 = 1.0 - Dy1

    bilinear_vector = Dy2*( Dx2*A_VII(:,i1,j1)   &
         +                  Dx1*A_VII(:,i2,j1))  &
         +            Dy1*( Dx2*A_VII(:,i1,j2)   &
         +                  Dx1*A_VII(:,i2,j2))

  end function bilinear_vector

  !=========================================================================
  real function trilinear_scalar( &
       A_III, iMin, iMax, jMin, jMax, kMin, kMax, Xyz_D, &
       x_I, y_I, z_I, DoExtrapolate)

    ! Calculate trilinear interpolation of A_III at position Xyz_D

    implicit none
    integer, intent(in) :: iMin, iMax, jMin, jMax, kMin, kMax
    real, intent(in)    :: A_III(iMin:iMax,jMin:jMax,kMin:kMax)
    real, intent(in)    :: Xyz_D(3)

    real,    intent(in), optional :: x_I(iMin:iMax)
    real,    intent(in), optional :: y_I(jMin:jMax)
    real,    intent(in), optional :: z_I(kMin:kMax)
    logical, intent(in), optional :: DoExtrapolate

    integer :: i1, i2, j1, j2, k1, k2
    real    :: Dx1, Dx2, Dy1, Dy2, Dz1, Dz2
    character (len=*), parameter :: NameSub=NameMod//'::trilinear_scalar'
    !--------------------------------------------------------------------------
    call find_cell(iMin, iMax, Xyz_D(1), i1, Dx1, x_I, DoExtrapolate, &
         "Called for coord1 from "//NameSub)

    call find_cell(jMin, jMax, Xyz_D(2), j1, Dy1, y_I, DoExtrapolate, &
         "Called for coord2 from "//NameSub)

    call find_cell(kMin, kMax, Xyz_D(3), k1, Dz1, z_I, DoExtrapolate, &
         "Called for coord3 from "//NameSub)
    
    ! Perform interpolation (or extrapolation)
    i2 = i1 + 1; Dx2 = 1.0 - Dx1
    j2 = j1 + 1; Dy2 = 1.0 - Dy1
    k2 = k1 + 1; Dz2 = 1.0 - Dz1

    !Perform interpolation (or extrapolation)
    trilinear_scalar = Dz2*( Dy2*( Dx2*A_III(i1,j1,k1)   &
         +                         Dx1*A_III(i2,j1,k1))  &
         +                   Dy1*( Dx2*A_III(i1,j2,k1)   &
         +                         Dx1*A_III(i2,j2,k1))) &
         +             Dz1*( Dy2*( Dx2*A_III(i1,j1,k2)   &
         +                         Dx1*A_III(i2,j1,k2))  &
         +                   Dy1*( Dx2*A_III(i1,j2,k2)   &
         +                         Dx1*A_III(i2,j2,k2)))

  end function trilinear_scalar

  !===========================================================================

  function trilinear_vector( &
       A_VIII, nVar, iMin, iMax, jMin, jMax, kMin, kMax, Xyz_D, &
       x_I, y_I, z_I, DoExtrapolate)

    ! Calculate trilinear interpolation of A_III at position Xyz_D

    implicit none
    integer, intent(in) :: nVar, iMin, iMax, jMin, jMax, kMin, kMax
    real, intent(in)    :: A_VIII(nVar, iMin:iMax, jMin:jMax, kMin:kMax)
    real, intent(in)    :: Xyz_D(3)

    real,    intent(in), optional :: x_I(iMin:iMax)
    real,    intent(in), optional :: y_I(jMin:jMax)
    real,    intent(in), optional :: z_I(kMin:kMax)
    logical, intent(in), optional :: DoExtrapolate

    ! return value
    real :: trilinear_vector(nVar)

    integer :: i1, i2, j1, j2, k1, k2
    real    :: Dx1, Dx2, Dy1, Dy2, Dz1, Dz2
    character (len=*), parameter :: NameSub=NameMod//'::trilinear_vector'
    !--------------------------------------------------------------------------
    call find_cell(iMin, iMax, Xyz_D(1), i1, Dx1, x_I, DoExtrapolate, &
         "Called for coord1 from "//NameSub)

    call find_cell(jMin, jMax, Xyz_D(2), j1, Dy1, y_I, DoExtrapolate, &
         "Called for coord2 from "//NameSub)

    call find_cell(kMin, kMax, Xyz_D(3), k1, Dz1, z_I, DoExtrapolate, &
         "Called for coord3 from "//NameSub)
    
    ! Perform interpolation (or extrapolation) of multiple variables
    i2 = i1 + 1; Dx2 = 1.0 - Dx1
    j2 = j1 + 1; Dy2 = 1.0 - Dy1
    k2 = k1 + 1; Dz2 = 1.0 - Dz1

    trilinear_vector = Dz2*(Dy2*(Dx2*A_VIII(:,i1,j1,k1)   &
         +                       Dx1*A_VIII(:,i2,j1,k1))  &
         +                  Dy1*(Dx2*A_VIII(:,i1,j2,k1)   &
         +                       Dx1*A_VIII(:,i2,j2,k1))) &
         +             Dz1*(Dy2*(Dx2*A_VIII(:,i1,j1,k2)   &
         +                       Dx1*A_VIII(:,i2,j1,k2))  &
         +                  Dy1*(Dx2*A_VIII(:,i1,j2,k2)   &
         +                       Dx1*A_VIII(:,i2,j2,k2)))

  end function trilinear_vector

  !===========================================================================
  subroutine find_cell(MinCoord, MaxCoord, Coord, iCoord, dCoord, &
       Coord_I, DoExtrapolate, StringError, IsInside)

    ! Find cell index and distance from cell for either 
    ! - a uniform grid with normalized coordinate (Coord_I is NOT present)
    ! - a nonuniform grid with monotone coordinates (Coord_I is present)
    !
    ! For sake of easy usage the returned coordinate index iCoord always
    ! satisfies MinCoord <= iCoord < MaxCoord.
    !
    ! If the coordinate is out of bounds, and DoExtrapolate is not present,
    ! the code stops with an error message. If DoExtrapolate is present and
    ! false, dCoord is modified to 0 or 1 so that the last grid cell is used.
    ! If DoExtrapolate is true, iCoord and dCoord are set 
    ! corresponding to linear extrapolation from the last two grid values.
    !
    ! For interpolation the normalized distance dCoord measured from 
    ! coordinate iCoord satisfies 0.0 <= dCoord <= 1.0
    ! but for extrapolation dCoord < 0.0 or > 1.0 is also possible.
    !---
    ! FOR THE UNIFORM CASE the normalized coordinate Coord should be equal to
    ! the index at the cell centers, therefore:
    !
    ! iCoord = max(MinCoord, min(MaxCoord-1, floor(Coord)))
    ! dCoord = Coord - iCoord
    !
    ! The optional IsInside = MinCoord <= Coord <= MaxCoord
    !----
    ! IN THE NON-UNIFORM CASE the cell iCoord that is left to coordinate Coord
    ! is found with a binary search in the Coord_I coordinates.
    !
    ! The normalized distance is set to
    !    dCoord = (Coord-Coord_I(iCoord))/(Coord_I(iCoord+1)-Coord_I(iCoord))
    !
    ! The optional IsInside = Coord_I(1) <= Coord <= Coord_I(nCoord).
    !----
    ! Example for linear interpolation on a 1D uniform grid of nX cells,
    ! DeltaX cell size and the first cell center is at DeltaX/2:
    !
    !   call find_cell(1, nX, x/DeltaX+0.5, iX, d)
    !   State_V = (1.0 - d)*State_VC(:,iX) + d*State_VC(:,iX+1)
    !
    ! Example for linear interpolation on a 1D non-uniform grid 
    ! with 2 ghost cells:
    !
    !   call find_cell(-1, nI+2, x, iX, d, x_G)
    !   State_V = (1.0 - d)*State_VG(:,iX) + d*State_VG(:,iX+1)

    integer, intent(in)           :: MinCoord, MaxCoord
    real,    intent(in)           :: Coord
    integer, intent(out)          :: iCoord
    real,    intent(out), optional:: dCoord
    real,    intent(in),  optional:: Coord_I(MinCoord:MaxCoord)
    logical, intent(in),  optional:: DoExtrapolate
    character(len=*),     optional:: StringError
    logical, intent(out), optional:: IsInside

    integer:: i, Di

    character(len=*), parameter:: NameSub="ModInterpolate::find_cell"
    !------------------------------------------------------------------------

    if(.not.present(Coord_I))then
       ! Uniform grid case with normalized coordinate

       iCoord = min(MaxCoord-1, max(MinCoord, floor(Coord)))
       dCoord = Coord - iCoord

       if(Coord < MinCoord)then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) NameSub, ': ', StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord=', Coord
             call CON_stop(NameSub//': normalized coordinate is to small!')
          elseif(.not.DoExtrapolate)then
             ! Use lefttmost cell (first order accurate) 
             dCoord = 0.0
          end if
          if(present(IsInside)) IsInside = .false.
       elseif(Coord > MaxCoord)then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord=', Coord
             call CON_stop(NameSub//': normalized coordinate is too large!')
          elseif(.not.DoExtrapolate)then
             ! Use rightmost cell (first order accurate) 
             dCoord = 1.0
          endif
          if(present(IsInside)) IsInside = .false.
       else
          if(present(IsInside)) IsInside = .true.
       end if

    elseif(Coord_I(MinCoord) < Coord_I(MaxCoord))then

       ! Monotone increasing coordinates

       if(Coord < Coord_I(MinCoord))then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord, CoordMin, CoordMax=', &
                  Coord, Coord_I(MinCoord), Coord_I(MaxCoord)
             call CON_stop(NameSub//': coordinate is too small!')
          elseif(DoExtrapolate)then
             iCoord = MinCoord
             dCoord = (Coord - Coord_I(iCoord)) &
                  /   (Coord_I(iCoord+1) - Coord_I(iCoord))
          else
             iCoord = MinCoord
             dCoord = 0.0
          end if
          if(present(IsInside)) IsInside = .false.
          RETURN
       end if

       if(Coord > Coord_I(MaxCoord))then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord, CoordMin, CoordMax=', &
                  Coord, Coord_I(MinCoord), Coord_I(MaxCoord)
             call CON_stop(NameSub//': coordinate is too large!')
          elseif(DoExtrapolate)then
             iCoord = MaxCoord - 1
             dCoord = (Coord - Coord_I(iCoord))  &
                  /   (Coord_I(iCoord+1) - Coord_I(iCoord))
          else
             iCoord = MaxCoord - 1
             dCoord = 1.0
          end if
          if(present(IsInside)) IsInside = .false.
          RETURN
       end if

       if(present(IsInside)) IsInside = .true.

       ! binary search
       i  = (MinCoord + MaxCoord)/2
       Di = (MaxCoord - MinCoord)/2
       do
          Di = (Di + 1)/2
          if(Coord < Coord_I(i)) then
             i = max(MinCoord, i - Di)
          elseif(Coord > Coord_I(i+1))then
             i = min(MaxCoord-1, i + Di)
          else
             EXIT
          end if
       end do
       iCoord = i
       dCoord = (Coord             - Coord_I(iCoord)) &
            /   (Coord_I(iCoord+1) - Coord_I(iCoord))

    else

       ! Monotone decreasing coordinates

       if(Coord < Coord_I(MaxCoord))then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord, CoordMin, CoordMax=', &
                  Coord, Coord_I(MaxCoord), Coord_I(MinCoord)
             call CON_stop(NameSub//': coordinate is too small!')
          elseif(DoExtrapolate)then
             iCoord = MaxCoord - 1
             dCoord = (Coord_I(iCoord) - Coord) &
                  /   (Coord_I(iCoord) - Coord_I(iCoord+1))
          else
             iCoord = MaxCoord - 1
             dCoord = 1.0
          end if
          if(present(IsInside)) IsInside = .false.
          RETURN
       end if

       if(Coord > Coord_I(MinCoord))then
          if(.not. (present(DoExtrapolate))) then
             if(present(StringError)) write(*,*) StringError
             write(*,*) NameSub,': MinIndex, MaxIndex=', MinCoord, MaxCoord
             write(*,*) NameSub,': Coord, CoordMin, CoordMax=', &
                  Coord, Coord_I(MaxCoord), Coord_I(MinCoord)
             call CON_stop(NameSub//': coordinate is too large!')
          elseif(DoExtrapolate)then
             iCoord = MinCoord
             dCoord = (Coord_I(iCoord) - Coord)  &
                  /   (Coord_I(iCoord) - Coord_I(iCoord+1))
          else
             iCoord = MinCoord
             dCoord = 0.0
          end if
          if(present(IsInside)) IsInside = .false.
          RETURN
       end if

       if(present(IsInside)) IsInside = .true.

       ! binary search
       i  = (MinCoord + MaxCoord)/2
       Di = (MaxCoord - MinCoord)/2
       do
          Di = (Di + 1)/2
          if(Coord > Coord_I(i)) then
             i = max(MinCoord, i - Di)
          elseif(Coord < Coord_I(i+1))then
             i = min(MaxCoord-1, i + Di)
          else
             EXIT
          end if
       end do
       iCoord = i
       dCoord = (Coord_I(iCoord) - Coord  ) &
            /   (Coord_I(iCoord) - Coord_I(iCoord+1))

    end if

  end subroutine find_cell
  !===========================================================================

  subroutine test_interpolation

    real :: A_I(0:2) = (/ 10., 20., 30. /)

    real :: A_II(2,3) = reshape((/ 1., 20., 3., 40., 5., 60. /), (/2, 3/))

    real :: A_VII(2,2,3) = reshape( &
         (/1., 10., 20., 200., 3., 30., 40., 400., 5., 50., 60., 600./), &
         (/2, 2, 3/))

    real :: A_III(2,2,0:2) = reshape((/ &
         1., 20., 3., 40., &
         100., 2000., 300., 4000., &
         10000., 200000., 30000., 400000. /), (/2, 2, 3/))

    real :: A_VIII(2,2,2,0:2) = reshape((/ &
         1., -10., 20., -200., 3., -30., 40., -400., &
         100., -1000., 2000., -20000., 300., -3000., 4000., -40000.,  &
         1e4, -1e5, 2e5, -2e6, 3e4, -3e5, 4e5, -4e6 /), (/2, 2, 2, 3/))

    real :: x12_I(1:2) = (/ 1., 2./)
    real :: x13_I(1:3) = (/ 1., 2., 4./)
    real :: x02_I(0:2) = (/ 1., 2., 4./)

    integer, parameter:: MinCoord = 1, MaxCoord = 8
    real   :: Coord_I(MinCoord:MaxCoord) = &
         (/ 1.0, 2.0, 4.0, 8.0, 16.0, 17.0, 24.0, 25.0 /)
    integer:: nCoord, iSign
    real   :: Coord, dCoord, CoordMin, CoordMax
    integer:: i, iCoord
    logical:: IsInside

    real :: Result, GoodResult, Result_V(2), GoodResult_V(2)
    logical :: DoExtrapolate = .false.

    character(len=*), parameter:: NameSub=NameMod//"::test_interpolation"
    !----------------------------------------------------------------------
    ! Change sign of coordinates to test for increasing and decreasing orders
    do iSign = 1, -1, -2
       if(iSign == 1)then
          write(*,'(a)')'Testing find_cell for increasing coordinates'
       else
          write(*,'(a)')'Testing find_cell for decreasing coordinates'
       end if

       ! Change number of coordinates to test binary search
       do nCoord = MaxCoord/2, MaxCoord

          ! Search for all integer coordinates
          ! starting below and finishing above the coordinate range

          CoordMin = min(Coord_I(MinCoord), Coord_I(nCoord))
          CoordMax = max(Coord_I(MinCoord), Coord_I(nCoord))

          do i = ceiling(CoordMin) - 1, floor(CoordMax) + 1
             Coord = i
             call find_cell(MinCoord, nCoord, Coord, &
                  iCoord, dCoord, Coord_I, .false., &
                  'Called from '//NameSub, IsInside)

             if(iSign*Coord < iSign*Coord_I(MinCoord))then
                if(IsInside) write(*,*) &
                     'Test failed for nCoord, Coord=', nCoord, Coord, &
                     ', IsInside=T, should be false'
                if(iCoord /= MinCoord) write(*,*)&
                     'Test failed for nCoord, Coord=', nCoord, Coord, &
                     ', iCoord=', iCoord, ' should be ', MinCoord
                if(dCoord /= 0.0) write(*,*) &
                     'Test failed for nCoord, Coord=', nCoord, Coord, &
                     ', dCoord=', dCoord, ' should be 0.0'
                CYCLE
             end if
             if(iSign*Coord > iSign*Coord_I(nCoord))then
                if(IsInside) write(*,*) &
                     'Test failed for nCoord, Coord=', nCoord, Coord, &
                     ', IsInside=T, should be false'
                if(iCoord /= nCoord - 1) write(*,*) &
                     'Test failed for nCoord, Coord=', nCoord, Coord, &
                     ', iCoord=', iCoord, ' should be ', nCoord - 1
                if(dCoord /= 1.0) write(*,*) &
                     'Test failed for nCoord, Coord=', nCoord, Coord, &
                     ', dCoord=', dCoord, ' should be 1.0'
                CYCLE
             end if
             if(.not.IsInside) write(*,*) &
                  'Test failed for nCoord, Coord=', nCoord, Coord, &
                  ', IsInside=F, should be true'

             if(iCoord < MinCoord .or. iCoord > nCoord-1) then
                write(*,*) &
                     'Test failed for nCoord, Coord=', nCoord, Coord, &
                     ', iCoord=', iCoord, ' should be < ', MinCoord, &
                     ' and > ', nCoord - 1
                CYCLE
             end if

             if(iSign*Coord_I(iCoord) > iSign*Coord) write(*,*) &
                  'Test failed for nCoord, Coord=', nCoord, Coord, &
                  ', iSign*Coord_I(iCoord)=', iSign*Coord_I(iCoord), &
                  ' should be <= iSign*Coord'

             if(iSign*Coord_I(iCoord+1) < iSign*Coord) write(*,*)       &
                  'Test failed for nCoord, Coord=', nCoord, Coord, &
                  ', iSign*Coord_I(iCoord+1)=', iSign*Coord_I(iCoord+1), &
                  ' should be >= iSign*Coord' 
             if(abs(Coord_I(iCoord) &
                  + dCoord*(Coord_I(iCoord+1) - Coord_I(iCoord)) &
                  - Coord) > 1e-6) write(*,*) &
                  'Test failed for nCoord, Coord=', nCoord, Coord, &
                  ', Coord_I(iCoord:iCoord+1)=', Coord_I(iCoord:iCoord+1), &
                  ', but incorrect dCoord = ', dCoord
          end do
       end do
       ! Change signs of coordinates to test decreasing order
       Coord_I = -Coord_I
    end do

    !Test for normal conditions.
    write(*,'(a)')'Testing function linear for uniform grid'
    Result = linear(A_I, 0, 2, 1.1)
    GoodResult = 21.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing function linear for non-uniform grid'
    Result = linear(A_I, 0, 2, 2.2, x02_I)
    GoodResult = 21.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing function bilinear for uniform grid'
    Result = bilinear(A_II, 1, 2, 1, 3, (/1.1, 2.2/))
    GoodResult = 7.46
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing function bilinear for non-uniform grid'
    Result = bilinear(A_II, 1, 2, 1, 3, (/1.1, 2.2/), x12_I, x13_I)
    GoodResult = 7.08
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing function trilinear for uniform grid'
    Result = trilinear(A_III, 1, 2, 1, 2, 0, 2, (/1.1, 1.2, 1.3/))
    GoodResult = 11236.2
    if(abs(Result - GoodResult) > 1.e-2) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing function trilinear for nonuniform grid'
    Result = trilinear(A_III, 1, 2, 1, 2, 0, 2, (/1.1, 1.2, 1.3/), &
         x12_I, x12_I, x02_I)
    GoodResult = 112.362
    if(abs(Result - GoodResult) > 1.e-2) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    !Test out-of-bounds, no extrapolation
    write(*,'(a)')'Testing bilinear out-of-bounds: +X for uniform grid'
    Result = bilinear(A_II, 1, 2, 1, 3, (/3.,1./), DoExtrapolate=.false.)
    GoodResult = 20.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing bilinear out-of-bounds: +X for nonuniform grid'
    Result = bilinear(A_II, 1, 2, 1, 3, (/3.,1./), x12_I, x13_I, &
         DoExtrapolate=.false.)
    GoodResult = 20.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing bilinear out-of-bounds: -X for uniform grid'
    Result = bilinear(A_II, 1, 2, 1, 3, (/-3.,2./), DoExtrapolate=.false.)
    GoodResult = 3.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing bilinear out-of-bounds: -X for nonuniform grid'
    Result = bilinear(A_II, 1, 2, 1, 3, (/-3.,2./), x12_I, x13_I, &
         DoExtrapolate=.false.)
    GoodResult = 3.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing bilinear out-of-bounds: +Y for uniform grid'
    Result = bilinear(A_II, 1, 2, 1, 3, (/1.,6./), DoExtrapolate=.false.)
    GoodResult = 5.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing bilinear out-of-bounds: +Y for nonuniform grid'
    Result = bilinear(A_II, 1, 2, 1, 3, (/1.,6./), x12_I, x13_I, &
         DoExtrapolate=.false.)
    GoodResult = 5.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing bilinear out-of-bounds: -Y for uniform grid'
    Result = bilinear(A_II, 1, 2, 1, 3, (/2.,-3./), DoExtrapolate=.false.)
    GoodResult = 20.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing bilinear out-of-bounds: -Y for nonuniform grid'
    Result = bilinear(A_II, 1, 2, 1, 3, (/2.,-3./), x12_I, x13_I, &
         DoExtrapolate=.false.)
    GoodResult = 20.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing trilinear out-of-bounds: +Z for uniform grid'
    Result = trilinear(A_III, 1, 2, 1, 2, 0, 2, (/1., 1., 2.4/), &
         DoExtrapolate=.false.)
    GoodResult = 10000.0
    if(abs(Result - GoodResult) > 1.e-2) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing trilinear out-of-bounds: +Z for nonuniform grid'
    Result = trilinear(A_III, 1, 2, 1, 2, 0, 2, (/1., 1., 4.1/), &
         x12_I, x12_I, x02_I, DoExtrapolate=.false.)
    GoodResult = 10000.0
    if(abs(Result - GoodResult) > 1.e-6) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing trilinear out-of-bounds: -Z for uniform grid'
    Result = trilinear(A_III, 1, 2, 1, 2, 0, 2, (/1., 1., -0.4/), &
         DoExtrapolate=.false.)
    GoodResult = 1.0
    if(abs(Result - GoodResult) > 1.e-6) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing trilinear out-of-bounds: -Z for nonuniform grid'
    Result = trilinear(A_III, 1, 2, 1, 2, 0, 2, (/1., 1., 0.1/), &
         x12_I, x12_I, x02_I, DoExtrapolate=.false.)
    GoodResult = 1.0
    if(abs(Result - GoodResult) > 1.e-2) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    !Test extrapolation
    write(*,'(a)')'Testing bilinear extrapolation: +X uniform'
    Result = bilinear(A_II, 1, 2, 1, 3, (/2.5,1./), DoExtrapolate=.true.)
    GoodResult = 29.5
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult
    
    write(*,'(a)')'Testing bilinear extrapolation: +X nonuniform'
    Result = bilinear(A_II, 1, 2, 1, 3, (/2.5,1./), x12_I, x13_I, &
         DoExtrapolate=.true.)
    GoodResult = 29.5
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult
    
    write(*,'(a)')'Testing bilinear extrapolation: -X uniform'
    Result = bilinear(A_II, 1, 2, 1, 3, (/.5,1.5/), DoExtrapolate=.true.)
    GoodResult = -12.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing bilinear extrapolation: -X nonuniform'
    Result = bilinear(A_II, 1, 2, 1, 3, (/.5,1.5/), x12_I, x13_I, &
         DoExtrapolate=.true.)
    GoodResult = -12.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing trilinear extrapolation: +Z uniform'
    Result = trilinear(A_III, 1, 2, 1, 2, 0, 2, (/1.3, 1.9, 2.60/), &
         DoExtrapolate=.true.)
    GoodResult = 212958.38
    if(abs(Result - GoodResult) > 1.0) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing trilinear extrapolation: +Z nonuniform'
    Result = trilinear(A_III, 1, 2, 1, 2, 0, 2, (/1.3, 1.9, 5.2/), &
         x12_I, x12_I, x02_I, DoExtrapolate=.true.)
    GoodResult = 212958.38
    if(abs(Result - GoodResult) > 1.0) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing function bilinear_vector'
    Result_V = bilinear(A_VII, 2, 1, 2, 1, 3, (/1.1, 2.2/))
    GoodResult_V = (/7.46, 74.6/)
    if(any(abs(Result_V - GoodResult_V) > 1.e-5)) &
         write(*,*) 'Test failed: Result=',Result_V,&
         ' differs from ',GoodResult_V

    write(*,'(a)')'Testing function trilinear_vector'
    Result_V = trilinear(A_VIII, 2, 1, 2, 1, 2, 0, 2, (/1.1, 1.2, 1.3/))
    GoodResult_V = (/ 11236.2, -112362.0 /)
    if(any(abs(Result_V - GoodResult_V) > 1.e-2)) write(*,*) &
         'Test failed: Result=', Result_V, ' differs from ', GoodResult_V

  end subroutine test_interpolation

end module ModInterpolate
