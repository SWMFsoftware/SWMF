module ModInterpolate

  ! Use second order interpolation in a uniform grid with normalized
  ! coordinates. The grid indexes are assumed to go from 0 to nI-1,
  ! 0 to nJ-1 etc. The zero index can be regarded as a ghost cell
  ! or as the boundary nodal value.
  ! The coordinates are normalized such that the
  ! the coordinates coincide with the indexes at the grid points.
  ! Typical usage:
  !
  ! InterpolatedValue = bilinear(Value_II, nI+2, nJ+2, &
  !                     (/ (x - x(0))/DeltaX, (y - y(0))/DeltaY) /) )

  implicit none

  private ! except

  public :: bilinear           ! 2nd order interpolation in 2D
  public :: trilinear          ! 2nd order interpolation in 3D
  public :: test_interpolation ! unit test

  character(len=*), parameter :: NameMod='ModInterpolate'

contains

  !=========================================================================
  real function bilinear(A_II, nI, nJ, Xy_D)

    ! Calculate bilinear interpolation of A_II at position Xy_D

    implicit none
    integer, intent(in) :: nI, nJ
    real, intent(in)    :: A_II(0:nI-1,0:nJ-1)
    real, intent(in)    :: Xy_D(2)

    integer :: i1, i2, j1, j2
    real :: Dx1, Dx2, Dy1, Dy2
    character (len=*), parameter :: NameSub=NameMod//'::trilinear'
    !--------------------------------------------------------------------------
    if(any(Xy_D < 0.0)) then
       write(*,*)'Xy_D =',Xy_D
       call CON_stop(NameSub//': negative coordinate is not valid')
    end if
    if(any(Xy_D > (/ nI-1, nJ-1 /))) then
       write(*,*)'nI-1, nJ-1=',nI-1, nJ-1
       write(*,*)'Xy_D =',Xy_D
       call CON_stop(NameSub//': normalized coordinates exceed array size')
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
  real function trilinear(A_III, nI, nJ, nK, Xyz_D)

    ! Calculate trilinear interpolation of A_III at position Xyz_D

    integer, intent(in) :: nI, nJ, nK
    real, intent(in) :: A_III(0:nI-1,0:nJ-1,0:nK-1)
    real, intent(in) :: Xyz_D(3)

    integer :: i1, i2, j1, j2, k1, k2
    real    :: Dx1, Dx2, Dy1, Dy2, Dz1, Dz2
    character (len=*), parameter :: NameSub=NameMod//'::trilinear'
    !--------------------------------------------------------------------------
    if(any(Xyz_D < 0.0)) then
       write(*,*)'Xyz_D =',Xyz_D
       call CON_stop(NameSub//': negative coordinate is not valid')
    end if
    if(any(Xyz_D > (/ nI-1, nJ-1, nK-1 /))) then
       write(*,*)'nI-1, nJ-1, nK-1=',nI-1, nJ-1, nK-1
       write(*,*)'Xyz_D =',Xyz_D
       call CON_stop(NameSub//': normalized coordinates exceed array size')
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
    
    real :: A_III(2,2,3) = reshape((/ &
         1., 20., 3., 40., &
         100., 2000., 300., 4000., &
         10000., 200000., 30000., 400000. /), (/2, 2, 3/))
    
    real :: Result, GoodResult

    write(*,'(a)')'Testing function bilinear'
    Result = bilinear(A_II, 2, 3, (/0.1, 1.2/))
    GoodResult = 7.46
    if(abs(Result - GoodResult) > 1.e-5) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

    write(*,'(a)')'Testing function trilinear'
    Result = trilinear(A_III, 2, 2, 3, (/0.1, 0.2, 1.3/))
    GoodResult = 11236.2
    if(abs(Result - GoodResult) > 1.e-2) &
         write(*,*) 'Test failed: Result=',Result,' differs from ',GoodResult

  end subroutine test_interpolation

end module ModInterpolate
