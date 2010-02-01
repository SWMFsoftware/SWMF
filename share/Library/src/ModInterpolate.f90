module ModInterpolate

  ! Use second order interpolation in a uniform grid with normalized
  ! coordinates. The coordinates are normalized such that the
  ! the coordinates coincide with the indexes at the grid points.
  !
  ! Examples of usage:
  !
  ! Cell based grid ghost cells:
  !
  !     InterpolatedValue = bilinear(Value_II, 0, nI+1, 0, nJ+1, &
  !                         (/ (x - x(0))/DeltaX, (y - y(0))/DeltaY) /) )
  !
  ! Node based grid with x(1)=y(1)=0.0:
  !
  !     InterpolatedValue = bilinear(Value_II, 1, nI, 1, nJ, &
  !                        (/ x/DeltaX, y/DeltaY /) )

  implicit none

  private ! except

  public :: bilinear           ! 2nd order interpolation in 2D
  public :: trilinear          ! 2nd order interpolation in 3D
  public :: find_cell          ! find cell in non-uniform grid
  public :: test_interpolation ! unit test

  character(len=*), parameter :: NameMod='ModInterpolate'

  interface bilinear
     module procedure bilinear_scalar, bilinear_vector
  end interface

  interface trilinear
     module procedure trilinear_scalar, trilinear_vector
  end interface

contains

  !=========================================================================
  real function bilinear_scalar( &
       A_II, iMin, iMax, jMin, jMax, Xy_D, DoExtrapolate)

    ! Calculate bilinear interpolation of A_II at position Xy_D

    implicit none

    integer, intent(in) :: iMin, iMax, jMin, jMax
    real, intent(in)    :: A_II(iMin:iMax,jMin:jMax)
    real, intent(in)    :: Xy_D(2)

    logical, intent(in), optional :: DoExtrapolate

    integer :: i1, i2, j1, j2
    real :: Dx1, Dx2, Dy1, Dy2
    character (len=*), parameter :: NameSub=NameMod//'::bilinear_scalar'
    !--------------------------------------------------------------------------
    !Set location assuming point is inside block.
    i1 = floor(Xy_D(1))
    j1 = floor(Xy_D(2))  
    i2 = ceiling(Xy_D(1))
    j2 = ceiling(Xy_D(2))

    !If Xy_D is outside of block, change i,j,k according to DoExtrapolate.
    if(any( Xy_D < (/iMin, jMin/)) .or. any(Xy_D > (/ iMax, jMax /))) then

       !Crash if DoExtrapolate is not set.
       if(.not. (present(DoExtrapolate))) then
          write(*,*)NameSub,&
               ' ERROR: Point outside of block and DoExtrapolate is not set!'
          write(*,*)NameSub,': iMin, iMax, jMin, jMax=',iMin, iMax, jMin, jMax
          write(*,*)NameSub,': Xy_D =',Xy_D
          call CON_stop(NameSub//': normalized coordinates are out of range')
       elseif(DoExtrapolate)then
          !Extrapolate point with second order accuracy
          i1 = min(iMax-1, max(iMin, i1));   i2 = i1 + 1
          j1 = min(jMax-1, max(jMin, j1));   j2 = j1 + 1
       else
          !Move point to closest edge (first order accurate)
          i1 = min(iMax, max(iMin, i1))
          i2 = min(iMax, max(iMin, i2))
          j1 = min(jMax, max(jMin, j1))
          j2 = min(jMax, max(jMin, j2))
       endif

    endif
       
    !Set interpolation weights
    Dx1= Xy_D(1) - i1;   Dx2 = 1.0 - Dx1
    Dy1= Xy_D(2) - j1;   Dy2 = 1.0 - Dy1

    !Perform interpolation (or extrapolation)
    bilinear_scalar = Dy2*( Dx2*A_II(i1,j1)   &
         +                  Dx1*A_II(i2,j1))  &
         +            Dy1*( Dx2*A_II(i1,j2)   &
         +                  Dx1*A_II(i2,j2))

  end function bilinear_scalar

  !=========================================================================
  function bilinear_vector( &
       A_VII, nVar, iMin, iMax, jMin, jMax, Xy_D, DoExtrapolate)

    ! Calculate bilinear interpolation of A_II at position Xy_D

    implicit none

    integer, intent(in) :: nVar, iMin, iMax, jMin, jMax
    real, intent(in)    :: A_VII(nVar, iMin:iMax,jMin:jMax)
    real, intent(in)    :: Xy_D(2)

    ! return value
    real                :: bilinear_vector(nVar)

    logical, intent(in), optional :: DoExtrapolate

    integer :: i1, i2, j1, j2
    real :: Dx1, Dx2, Dy1, Dy2
    character (len=*), parameter :: NameSub=NameMod//'::bilinear_vector'
    !--------------------------------------------------------------------------
    !Set location assuming point is inside block.
    i1 = floor(Xy_D(1))
    j1 = floor(Xy_D(2))  
    i2 = ceiling(Xy_D(1))
    j2 = ceiling(Xy_D(2))

    !If Xy_D is outside of block, change i,j,k according to DoExtrapolate.
    if(any( Xy_D < (/iMin, jMin/)) .or. any(Xy_D > (/ iMax, jMax /))) then

       !Crash if DoExtrapolate is not set.
       if(.not. (PRESENT(DoExtrapolate))) then
          write(*,*)NameSub,&
               ' ERROR: Point outside of block and DoExtrapolate is not set!'
          write(*,*)NameSub,': iMin, iMax, jMin, jMax=',iMin, iMax, jMin, jMax
          write(*,*)NameSub,': Xy_D =',Xy_D
          call CON_stop(NameSub//': normalized coordinates are out of range')
       elseif(DoExtrapolate)then
          !Extrapolate point with second order accuracy
          i1 = min(iMax-1, max(iMin, i1));   i2 = i1 + 1
          j1 = min(jMax-1, max(jMin, j1));   j2 = j1 + 1
       else
          !Move point to closest edge (first order accurate)
          i1 = min(iMax, max(iMin, i1))
          i2 = min(iMax, max(iMin, i2))
          j1 = min(jMax, max(jMin, j1))
          j2 = min(jMax, max(jMin, j2))
       endif

    endif
       
    !Set interpolation weights
    Dx1= Xy_D(1) - i1;   Dx2 = 1.0 - Dx1
    Dy1= Xy_D(2) - j1;   Dy2 = 1.0 - Dy1

    !Perform interpolation (or extrapolation)
    bilinear_vector = Dy2*( Dx2*A_VII(:,i1,j1)   &
         +                  Dx1*A_VII(:,i2,j1))  &
         +            Dy1*( Dx2*A_VII(:,i1,j2)   &
         +                  Dx1*A_VII(:,i2,j2))

  end function bilinear_vector

  !=========================================================================
  real function trilinear_scalar( &
       A_III, iMin, iMax, jMin, jMax, kMin, kMax, Xyz_D, DoExtrapolate)

    ! Calculate trilinear interpolation of A_III at position Xyz_D

    implicit none
    integer, intent(in) :: iMin, iMax, jMin, jMax, kMin, kMax
    real, intent(in)    :: A_III(iMin:iMax,jMin:jMax,kMin:kMax)
    real, intent(in)    :: Xyz_D(3)
    logical, intent(in), optional :: DoExtrapolate

    integer :: i1, i2, j1, j2, k1, k2
    real    :: Dx1, Dx2, Dy1, Dy2, Dz1, Dz2
    character (len=*), parameter :: NameSub=NameMod//'::trilinear_scalar'
    !--------------------------------------------------------------------------
    !Set location assuming point is inside block.
    i1 = floor(Xyz_D(1))
    j1 = floor(Xyz_D(2))  
    k1 = floor(Xyz_D(3))
    i2 = ceiling(Xyz_D(1))
    j2 = ceiling(Xyz_D(2))
    k2 = ceiling(Xyz_D(3))

    !If Xy_D is outside of block, change i,j,k according to DoExtrapolate.
    if(any( Xyz_D < (/iMin, jMin, kMin/)) .or. &
         any(Xyz_D > (/iMax, jMax, kMax/))) then

       !Crash if DoExtrapolate is not set.
       if(.not. present(DoExtrapolate)) then
          write(*,*)NameSub,&
               ' ERROR: Point outside of block & DoExtrapolate is not set!'
          write(*,*)NameSub,' iMin, iMax, jMin, jMax, kMin, kMax=', &
               iMin, iMax, jMin, jMax, kMin, kMax
          write(*,*)NameSub,' Xyz_D =',Xyz_D
          call CON_stop(NameSub//': normalized coordinates are out of range')
       elseif(DoExtrapolate) then
          !Extrapolate point with second order accuracy
          i1 = min(iMax-1, max(iMin, i1));   i2 = i1 + 1
          j1 = min(jMax-1, max(jMin, j1));   j2 = j1 + 1
          k1 = min(kMax-1, max(kMin, k1));   k2 = k1 + 1
       else
          !Move point to closest edge (first order accurate)
          i1 = min(iMax, max(iMin, i1))
          i2 = min(iMax, max(iMin, i2))
          j1 = min(jMax, max(jMin, j1))
          j2 = min(jMax, max(jMin, j2))
          k1 = min(kMax, max(kMin, k1))
          k2 = min(kMax, max(kMin, k2))
       endif

    endif
    
    !Set interpolation weights
    Dx1= Xyz_D(1) - i1; Dx2 = 1.0 - Dx1
    Dy1= Xyz_D(2) - j1; Dy2 = 1.0 - Dy1
    Dz1= Xyz_D(3) - k1; Dz2 = 1.0 - Dz1

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
       A_VIII, nVar, iMin, iMax, jMin, jMax, kMin, kMax, Xyz_D, DoExtrapolate)

    ! Calculate trilinear interpolation of A_III at position Xyz_D

    implicit none
    integer, intent(in) :: nVar, iMin, iMax, jMin, jMax, kMin, kMax
    real, intent(in)    :: A_VIII(nVar, iMin:iMax, jMin:jMax, kMin:kMax)
    real, intent(in)    :: Xyz_D(3)
    logical, intent(in), optional :: DoExtrapolate

    ! return value
    real :: trilinear_vector(nVar)

    integer :: i1, i2, j1, j2, k1, k2
    real    :: Dx1, Dx2, Dy1, Dy2, Dz1, Dz2
    character (len=*), parameter :: NameSub=NameMod//'::trilinear_vector'
    !--------------------------------------------------------------------------
    !Set location assuming point is inside block.
    i1 = floor(Xyz_D(1))
    j1 = floor(Xyz_D(2))  
    k1 = floor(Xyz_D(3))
    i2 = ceiling(Xyz_D(1))
    j2 = ceiling(Xyz_D(2))
    k2 = ceiling(Xyz_D(3))

    !If Xy_D is outside of block, change i,j,k according to DoExtrapolate.
    if(any( Xyz_D < (/iMin, jMin, kMin/)) .or. &
         any(Xyz_D > (/iMax, jMax, kMax/))) then

       !Crash if DoExtrapolate is not set.
       if(.not. present(DoExtrapolate)) then
          write(*,*)NameSub,&
               ' ERROR: Point outside of block & DoExtrapolate is not set!'
          write(*,*)NameSub,' iMin, iMax, jMin, jMax, kMin, kMax=', &
               iMin, iMax, jMin, jMax, kMin, kMax
          write(*,*)NameSub,' Xyz_D =',Xyz_D
          call CON_stop(NameSub//': normalized coordinates are out of range')
       elseif(DoExtrapolate) then
          !Extrapolate point with second order accuracy
          i1 = min(iMax-1, max(iMin, i1));   i2 = i1 + 1
          j1 = min(jMax-1, max(jMin, j1));   j2 = j1 + 1
          k1 = min(kMax-1, max(kMin, k1));   k2 = k1 + 1
       else
          !Move point to closest edge (first order accurate)
          i1 = min(iMax, max(iMin, i1))
          i2 = min(iMax, max(iMin, i2))
          j1 = min(jMax, max(jMin, j1))
          j2 = min(jMax, max(jMin, j2))
          k1 = min(kMax, max(kMin, k1))
          k2 = min(kMax, max(kMin, k2))
       endif

    endif
    
    !Set interpolation weights
    Dx1= Xyz_D(1) - i1; Dx2 = 1.0 - Dx1
    Dy1= Xyz_D(2) - j1; Dy2 = 1.0 - Dy1
    Dz1= Xyz_D(3) - k1; Dz2 = 1.0 - Dz1

    !Perform interpolation (or extrapolation)
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
  subroutine find_cell(Coord, Coord_I, MinCoord, MaxCoord, &
       iCoord, dCoord, IsInside)

    ! Coord_I(MinCoord:MaxCoord) are cell coordinates in increasing order.
    ! The goal is to find cell index iCoord that is left to coordinate Coord.
    ! For sake of easy use MinCoord <= iCoord < MaxCoord always holds.
    ! The optional IsInside logical is set to true if 
    !    Coord_I(1) <= Coord <= Coord_I(nCoord).
    ! If Coord is inside the index iCoord is set such that 
    !    Coord_I(iCoord) <= Coord <= Coord_I(nCoord)
    ! and the normalized distance is set to
    !    dCoord = (Coord-Coord_I(iCoord))/(Coord_I(iCoord+1)-Coord_I(iCoord))
    ! If Coord is outside to the left,  then iCoord=MinCoord,   dCoord=0.0
    ! if Coord is outside to the right, then iCoord=MaxCoord-1, dCoord=1.0
    !
    ! Example for linear interpolation on a 1D non-uniform grid 
    ! with 2 ghost cells:
    !
    !   call find_cell(x, x_G, -1, nI+2, i, d)
    !   State_V = (1.0 - d)*State_VG(:,i) + d*State_VG(:,i+1)

    real,    intent(in)           :: Coord
    integer, intent(in)           :: MinCoord, MaxCoord
    real,    intent(in)           :: Coord_I(MinCoord:MaxCoord)
    integer, intent(out)          :: iCoord
    real,    intent(out), optional:: dCoord
    logical, intent(out), optional:: IsInside

    integer:: i, Di
    !------------------------------------------------------------------------

    if(Coord < Coord_I(MinCoord))then
       iCoord = MinCoord
       if(present(dCoord))   dCoord   = 0.0
       if(present(IsInside)) IsInside = .false.
       RETURN
    end if
    if(Coord > Coord_I(MaxCoord))then
       iCoord = MaxCoord - 1
       if(present(dCoord))   dCoord   = 1.0
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

    if(present(dCoord)) dCoord = &
         (Coord             - Coord_I(iCoord)) / &
         (Coord_I(iCoord+1) - Coord_I(iCoord))

  end subroutine find_cell
  !===========================================================================

  subroutine test_interpolation

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

    integer, parameter:: MinCoord = 1, MaxCoord = 8
    real   :: Coord_I(MinCoord:MaxCoord) = &
         (/ 1.0, 2.0, 4.0, 8.0, 16.0, 17.0, 24.0, 25.0 /)
    integer:: nCoord
    real   :: Coord, dCoord
    integer:: i, iCoord
    logical:: IsInside

    real :: Result, GoodResult, Result_V(2), GoodResult_V(2)
    logical :: DoExtrapolate = .false.
    !----------------------------------------------------------------------
    !Test for normal conditions.
    write(*,'(a)')'Testing function bilinear'
    Result = bilinear(A_II, 1, 2, 1, 3, (/1.1, 2.2/))
    GoodResult = 7.46
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing function trilinear'
    Result = trilinear(A_III, 1, 2, 1, 2, 0, 2, (/1.1, 1.2, 1.3/))
    GoodResult = 11236.2
    if(abs(Result - GoodResult) > 1.e-2) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    !Test out-of-bounds, no extrapolation
    write(*,'(a)')'Testing bilinear out-of-bounds: +X'
    Result = bilinear(A_II, 1, 2, 1, 3, (/3.,1./),DoExtrapolate)
    GoodResult = 20.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing bilinear out-of-bounds: -X'
    Result = bilinear(A_II, 1, 2, 1, 3, (/-3.,2./),DoExtrapolate)
    GoodResult = 3.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing bilinear out-of-bounds: +Y'
    Result = bilinear(A_II, 1, 2, 1, 3, (/1.,6./),DoExtrapolate)
    GoodResult = 5.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing bilinear out-of-bounds: -Y'
    Result = bilinear(A_II, 1, 2, 1, 3, (/2.,-3./),DoExtrapolate)
    GoodResult = 20.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing trilinear out-of-bounds: +Z'
    Result = trilinear(A_III, 1, 2, 1, 2, 0, 2, (/1., 1., 2.4/),DoExtrapolate)
    GoodResult = 10000.0
    if(abs(Result - GoodResult) > 1.e-2) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    !Test extrapolation
    DoExtrapolate = .true.
    write(*,'(a)')'Testing bilinear extrapolation out-of-bounds: +X'
    Result = bilinear(A_II, 1, 2, 1, 3, (/2.5,1./), .true.)
    GoodResult = 29.5
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult
    
    write(*,'(a)')'Testing bilinear extrapolation out-of-bounds: -X'
    Result = bilinear(A_II, 1, 2, 1, 3, (/.5,1.5/), .true.)
    GoodResult = -12.0
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing trilinear extrapolation out-of-bounds: +Z'
    Result = trilinear(A_III, 1, 2, 1, 2, 0, 2, (/1.3, 1.9, 2.60/), .true.)
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

    write(*,'(a)')'Testing find_cell'
    do nCoord = MaxCoord/2, MaxCoord
       do i = ceiling(Coord_I(MinCoord)) - 1, floor(Coord_I(nCoord)) + 1
          Coord = i
          call find_cell(Coord, Coord_I, MinCoord, nCoord, &
               iCoord, dCoord, IsInside)

          if(Coord < Coord_I(MinCoord))then
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
          if(Coord > Coord_I(nCoord))then
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

          if(Coord_I(iCoord) > Coord) write(*,*) &
               'Test failed for nCoord, Coord=', nCoord, Coord, &
               ', Coord_I(iCoord)=', Coord_I(iCoord), ' should be <= Coord'

          if(Coord_I(iCoord+1) < Coord) write(*,*)       &
               'Test failed for nCoord, Coord=', nCoord, Coord, &
               ', Coord_I(iCoord+1)=', Coord_I(iCoord+1), ' should be >= Coord' 
          if(abs(Coord_I(iCoord) + dCoord*(Coord_I(iCoord+1)-Coord_I(iCoord)) &
               - Coord) > 1e-6) write(*,*) &
               'Test failed for nCoord, Coord=', nCoord, Coord, &
               ', Coord_I(iCoord:iCoord+1)=', Coord_I(iCoord:iCoord+1), &
               ', but incorrect dCoord = ', dCoord
       end do
    end do

  end subroutine test_interpolation

end module ModInterpolate
