Module ModInterWrapper

contains

!=============================================================================
Subroutine bilinear_normalize(iInSize,jInSize,iOutSize,jOutSize,&
                                     iIn,jIn,iOut,jOut,A_II,Buffer_out)

    Use ModInterpolate
    Use ModKind 
    ! Calculate normalized grids for bilinear interpolation

    Implicit NONE

    integer, intent(in) :: iInSize, jInSize, iOutSize,jOutSize
    real, intent(in)    :: iIn(iInSize), jIn(jInSize),&
                           iOut(iOutSize),jOut(jOutSize)
    real, intent(in)    :: A_II(iInSize,jInSize)
    real, intent(out)   :: Buffer_out(iOutSize, jOutSize)

    integer :: i , j
    real :: Xy_D(2)
    real (Real8_) :: Dx, Dy 
    real :: iInNorm(iInSize), jInNorm(jInSize), &
            iOutNorm(iOutSize), jOutNorm(jOutSize)
    !--------------------------------------------------------------------------

    print *, "iIn"
    print *, iIn
    print *, "jIn"
    print *, jIn
    
    print *, "iOut"
    print *, iOut
    print *, "jOut"
    print *, jOut


    ! Find Grid Spacing for Normalizing
    Dx = iIn(2) - iIn(1)
    Dy = jIn(2) - jIn(1)

    ! Create Normalized Grid
    iInNorm = (iIn / Dx)
    jInNorm = (jIn / Dy)
    iOutNorm= iOut/ Dx
    jOutNorm= jOut/ Dx

!    print *, "iInNorm"
!    print *, iInNorm
!    print *, "jInNorm"
!    print *, jInNorm
    
!    print *, "iOutNorm"
!    print *, iOutNorm
!    print *, "jOutNorm"
!    print *, jOutNorm

    ! Loop over Output array size and send to Bilinear Func for Interpolation
    do i=1, iOutSize
        do j=1, jOutSize
            Xy_d = (/iOut,jOut/)
            Buffer_out(i,j) = bilinear(A_II, 1, iInSize, 1, jInSize, Xy_D,iIn, jIn, .false.)
        enddo
    enddo
end subroutine bilinear_normalize
!=============================================================================

end module ModInterWrapper
