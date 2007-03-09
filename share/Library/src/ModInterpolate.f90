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
  public :: test_interpolation ! unit test

  character(len=*), parameter :: NameMod='ModInterpolate'

contains

  !=========================================================================
  real function bilinear(A_II, iMin, iMax, jMin, jMax, Xy_D)

    ! Calculate bilinear interpolation of A_II at position Xy_D

    implicit none
    integer, intent(in) :: iMin, iMax, jMin, jMax
    real, intent(in)    :: A_II(iMin:iMax,jMin:jMax)
    real, intent(in)    :: Xy_D(2)

    integer :: i1, i2, j1, j2
    real :: Dx1, Dx2, Dy1, Dy2
    character (len=*), parameter :: NameSub=NameMod//'::trilinear'
    !--------------------------------------------------------------------------
    if(any( Xy_D < (/iMin, jMin/)) .or. any(Xy_D > (/ iMax, jMax /))) then
       write(*,*)'iMin, iMax, jMin, jMax=',iMin, iMax, jMin, jMax
       write(*,*)'Xy_D =',Xy_D
       call CON_stop(NameSub//': normalized coordinates are out of range')
    end if

    i1 = floor(Xy_D(1))
    j1 = floor(Xy_D(2))  
    i2 = ceiling(Xy_D(1))
    j2 = ceiling(Xy_D(2))
    Dx1= Xy_D(1) - i1; Dx2 = 1.0 - Dx1
    Dy1= Xy_D(2) - j1; Dy2 = 1.0 - Dy1

    !write(*,*)'nI,nJ=',nI,nJ
    !write(*,*)'Xy_D =',Xy_D
    !write(*,*)'i1,j1,i2,j2=',i1,j1,i2,j2
    !write(*,*)'Dx1,Dx2=',Dx1,Dx2
    !write(*,*)'Dy1,Dx2=',Dx1,Dx2
    !write(*,*)'A_II(i1,*)=',A_II(i1,j1:j2)
    !write(*,*)'A_II(i2,*)=',A_II(i2,j1:j2)

    bilinear = Dy2*(   Dx2*A_II(i1,j1)   &
         +             Dx1*A_II(i2,j1))  &
         +     Dy1*(   Dx2*A_II(i1,j2)   &
         +             Dx1*A_II(i2,j2))

  end function bilinear

  !=========================================================================
  real function trilinear(A_III, iMin, iMax, jMin, jMax, kMin, kMax, Xyz_D)

    ! Calculate trilinear interpolation of A_III at position Xyz_D

    implicit none
    integer, intent(in) :: iMin, iMax, jMin, jMax, kMin, kMax
    real, intent(in)    :: A_III(iMin:iMax,jMin:jMax,kMin:kMax)
    real, intent(in)    :: Xyz_D(3)

    integer :: i1, i2, j1, j2, k1, k2
    real    :: Dx1, Dx2, Dy1, Dy2, Dz1, Dz2
    character (len=*), parameter :: NameSub=NameMod//'::trilinear'
    !--------------------------------------------------------------------------

    if(  any( Xyz_D < (/iMin, jMin, kMin/)) .or. &
         any( Xyz_D > (/iMax, jMax, kMax/))) then
       write(*,*)'iMin, iMax, jMin, jMax, kMin, kMax=',&
            iMin, iMax, jMin, jMax, kMin, kMax
       write(*,*)'Xyz_D =',Xyz_D
       call CON_stop(NameSub//': normalized coordinates are out of range')
    end if
       
    i1 = floor(Xyz_D(1))
    j1 = floor(Xyz_D(2))  
    k1 = floor(Xyz_D(3))
    i2 = ceiling(Xyz_D(1))
    j2 = ceiling(Xyz_D(2))
    k2 = ceiling(Xyz_D(3))
    Dx1= Xyz_D(1) - i1; Dx2 = 1.0 - Dx1
    Dy1= Xyz_D(2) - j1; Dy2 = 1.0 - Dy1
    Dz1= Xyz_D(3) - k1; Dz2 = 1.0 - Dz1

    trilinear = Dz2*(   Dy2*(   Dx2*A_III(i1,j1,k1)   &
         +                      Dx1*A_III(i2,j1,k1))  &
         +              Dy1*(   Dx2*A_III(i1,j2,k1)   &
         +                      Dx1*A_III(i2,j2,k1))) &
         +      Dz1*(   Dy2*(   Dx2*A_III(i1,j1,k2)   &
         +                      Dx1*A_III(i2,j1,k2))  &
         +              Dy1*(   Dx2*A_III(i1,j2,k2)   &
         +                      Dx1*A_III(i2,j2,k2)))

  end function trilinear
  !===========================================================================
  subroutine test_interpolation

    real :: A_II(2,3) = reshape((/ 1., 20., 3., 40., 5., 60. /), (/2, 3/))
    
    real :: A_III(2,2,0:2) = reshape((/ &
         1., 20., 3., 40., &
         100., 2000., 300., 4000., &
         10000., 200000., 30000., 400000. /), (/2, 2, 3/))
    
    real :: Result, GoodResult

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

  end subroutine test_interpolation

end module ModInterpolate
