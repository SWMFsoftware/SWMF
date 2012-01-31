module ModInterGen

  implicit none

  character(len=*), parameter :: NameMod='ModInterGen'

  contains

  !=========================================================================
  subroutine bilinear_general( iInSize,jInSize,iOutSize,jOutSize,&
        iIn, jIn, iOut, jOut, input, output)

    ! Calculate bilinear interpolation of A_II at position Xy_D

    implicit none

    integer,intent(in)  :: iInSize, jInSize, iOutSize, jOutSize
    real                :: iIn(iInSize), jIn(jInSize)
    real                :: iOut(iOutSize), jOut(jOutSize)
    real                :: input(iInSize, jInSize)
    real, intent(out)   :: output(iOutSize, jOutSize)

    integer :: i, j
    integer :: n=0, m=0
    real :: Dx1, Dx2, Dy1, Dy2, Xy_D(2),D_x,D_y
    character (len=*), parameter :: NameSub=NameMod//'::bilinear_general'
    !--------------------------------------------------------------------------
    
    D_x = abs(iIn(3) -iIn(2))
    D_y = abs(jIn(3) -jIn(2))

    do i=1, iOutSize
    do j=1, jOutSize 

    xy_d = (/iOut(i),jOut(j)/)

    ! Find corresponding point in input arrays
    do n=1, iInSize 
        if (xy_d(1) .LT. iIn(n)) exit
    enddo

    do m=1, jInSize 
        if (xy_d(2) .LT. jIn(m)) exit
    enddo  

    !Set interpolation weights
    Dx1 = iIn(n) - xy_d(1); Dx2 = D_x - Dx1
    Dy1 = jIn(m) - xy_d(2); Dy2 = D_Y - Dy1 

    !Perform interpolation (or extrapolation)
    output(i,j) = (Dy1/D_y)*( (Dx1/D_x)*input(n,m)   &
         +                    (Dx2/D_x)*input(n+1,m))  &
         +        (Dy2/D_y)*( (Dx1/D_x)*input(n,m+1)   &
         +                    (Dx2/D_x)*input(n+1,m+1))
    enddo
    enddo

  end subroutine bilinear_general

end module ModInterGen
