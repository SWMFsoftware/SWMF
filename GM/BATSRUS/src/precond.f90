!^CFG COPYRIGHT UM
!^CFG FILE IMPLICIT
!=============================================================================
! GENERAL PRECONDITIONER FOR BLOCK HEPTADIAGONAL AND PENTADIAGONAL MATRICES. !
! DIRECT SOLVER FOR BLOCK TRIDIAGONAL MATRIX.                                !
!=============================================================================
!
! G. Toth 2001 
! (based on the original F77 block heptadiagonal preconditioner with
!  Eisenstat trick implemented by Auke van der Ploeg, 1997)
!
! subroutines: 
!          prehepta, Uhepta, Lhepta
!
! Usage: call prehepta with the original matrix to obtain L and U.
!        call Uhepta(.false.) to multiply a vector with U
!        call Uhepta(.true.)  to multiply a vector with U^{-1}
!        call Lhepta          to multiply a vector with L^{-1}
!
! To solve a tridiagonal system A.x=rhs use these steps
!
!        call prehepta(A)        ! A --> LU
!        x=rhs
!        call Lhepta(x)          ! x --> L^{-1}.rhs
!        call Uhepta(.true.,x)   ! x --> U^{-1}.L^{-1}.rhs = A^{-1}.rhs
!
! To solve a penta- or hepta-diagonal problem with symmetric preconditioning
!
! L^{-1}.A.U^{-1} U.x = L^{-1}.rhs
! 
! use the following steps:
!
!        call prehepta(A)                 ! A -> LU
!        call Lhepta(rhs)                 ! rhs'=L^{-1}.rhs
!        call Uhepta(.false.,x)           ! x'  = U.x (for initial guess)
!        call bicgstab(matvec_prec,x,rhs) ! solve A'.x'=rhs'
!        call Uhepta(.true.,x)            ! x = U^{-1}.x'
!
! The preconditioned matrix vector multiplication is in 
! subroutine matvec_prec(x,y,n), and it should calculate 
! y = A'.x = L^{-1}.A.U^{-1}.x as
!
!        y=x
!        call Uhepta(.true.,y)   ! multiply y with U^{-1}
!        call matvec(y)          ! multiply y with A
!        call Lhepta(y)          ! multiply y with L^{-1}
!
!
!=============================================================================

SUBROUTINE prehepta(nblock,N,M1,M2,alf_in,d,e,f,e1,f1,e2,f2)

  ! This routine constructs an incomplete block LU-decomposition 
  ! of a hepta- or penta-diagonal matrix in such a way that L+U has the 
  ! same blockstructure as A. For block tri-diagonal matrix the subroutine
  ! provides a full LU decompostion. 
  !
  ! For penta-diagonal matrix, set M2=nblock and e2 and f2 can be scalar reals.
  ! For tri-diagonal matrix set M1=M2=nblock and e1,f1,e2,f2 can be scalars.
  !
  !===================================================================
  !     Gustafsson modification
  !===================================================================
  !
  ! It is possible to encorporate the so-called Gustafsson
  ! modification for the blocks on the main diagonal.
  ! In this appoach, a splitting A=LU+R, is constructed in such a way
  ! that the block row sums of R are zero. For systems of
  ! linear equations coming from problems with an elliptic character
  ! this can strongly improve the convergence behaviour of 
  ! iterative methods. See page 22 of Phd thesis 'Preconditioning 
  ! for sparse ..... ' for an illustration of this phenomenon.

  use ModPrecond

  IMPLICIT NONE
  INTEGER, INTENT(IN)                        :: N, M1, M2, nblock
  REAL, INTENT(IN)                           :: alf_in
  REAL, INTENT(INOUT), DIMENSION(N,N,nblock) :: d,e,f,e1,f1,e2,f2

  !=======================================================================
  !     Description of arguments:
  !=======================================================================
  !
  ! nblock:  Number of diagonal blocks.
  ! N:       The size of the blocks.
  ! M1:      Distance of outer blocks to the main diagonal blocks.
  ! M2:      Distance of outer-most blocks to main diagonal blocks.
  !           1 < M1 < M2.
  !           The matrix has a blockstructure corresponding to
  !           a seven-point stencil on a three-dimensional, 
  !           rectangular grid. The blocks corresonding to the
  !           direction in which grid points are numbered
  !           first are the sub- and superdiagonal blocks.
  !           The blocks corresponding to the direction in which grid 
  !           points are numbered secondly have distance M1 from the main 
  !           diagonal blocks. Finally, the blocks corresponding to the
  !           direction in which grid points are numbered
  !           last have distance M2 from the main diagonal blocks.
  !
  ! alf_in:   The parameter for Gustafsson modification. alf>=0 means 
  !           no modification, 0> alf_in >= -1 is the valid parameter range.
  !
  ! d, e, f, e1, f1, e2, f2:
  !          on entrance: matrix A
  !          on exit: L + U - I (the diagonal of U is I)
  !           The matrix A and L+U are block heptadiagonal.
  !           The blocks are stored as follows:
  !           d(j): j=1..nblock        main diagonal
  !           e(j): j=2..nblock        sub diagonal blocks.
  !           f(j): j=1..nblock-1      super diagonal blocks.
  !           e1(j): j=M1+1..nblock    blocks in the lower-triangular
  !                                     part with distance M1 from 
  !                                     the main diagonal.
  !           f1(j): j=1..nblock-M1    blocks in the upper-triangular
  !                                     part with distance M1 from 
  !                                     the main diagonal.
  !           e2(j): j=M2+1..nblock    blocks in the lower-triangular
  !                                     part with distance M2 from 
  !                                     the main diagonal.
  !           f2(j): j=1..nblock-M2    blocks in the upper-triangular
  !                                     part with distance M2 from 
  !                                     the main diagonal.
  !          It is assumed that the
  !          blocks are not very sparse, so the sparsity pattern
  !          within the separate blocks is not exploited.
  !          For example, the (i,k)-element of the j-th block on the 
  !          main diagonal is stored in d(i,k,j).
  !
  !     Local variables:
  !
  ! dd:      a single block for manipulating diagonal block d(*,*,j)
  !
  ! pivot:   integer array which contains the sequence generated
  !          by partial pivoting in subroutine 'Lapack_getrf'.
  !
!!! Automatic arrays
!!$  REAL    :: dd(N,N)
!!$  INTEGER :: pivot(N)
  real, dimension(:,:), allocatable :: dd
  integer, dimension(:), allocatable :: pivot

  REAL    :: alf
  INTEGER :: i,j,INFO, iError
  REAL, PARAMETER :: zero=0.0, one=1.0

  ! ALF      : internal value for Gustafsson parameter
  ! INFO     : variable for LAPACK_GETRF
  ! zero, one: variables necessary for BLAS-routine dgemv. 
  !
  ! External subroutines:
  !
  ! DGEMM,   BLAS level three Matrix-Matrix Product.
  !          See 'man DGEMM' for description.
  !
  ! Lapack_getrf,  LAPACK routine, computes an LU factorization of a general 
  !          M-by-N matrix A using partial pivoting with 
  !          row interchanges. The factorization has the form
  !              A = P * L * U
  !          where P is a permutation matrix, L is lower triangular 
  !          with unit diagonal elements (lower trapezoidal if m > n),
  !          and U is upper triangular (upper trapezoidal if m < n).
  !          This is the right-looking Level 3 BLAS version of the
  !          algorithm.
  !
  ! Lapack_getrs,  LAPACK routine, solves a system of linear equations
  !          A * X = B,  A**T * X = B,  or  A**H * X = B
  !          with a general N-by-N matrix A using the LU factorization
  !          computed by LAPACK_GETRF.
  !
  !-----------------------------------------------------------------------------

  call timing_start('precond')

  ! Allocate arrays that were "Automatic"
  allocate(dd(N,N), stat=iError); call alloc_check(iError,"precond:dd")
  allocate(pivot(N), stat=iError); call alloc_check(iError,"precond:pivot")

  alf=alf_in
  IF (alf < zero) THEN
     ! (Relaxed form of) Gustafsson modification:

     IF (alf < -one) THEN
        PRINT *,'Parmeter alf replaced by -1.0'
        alf=-one
        PRINT *,' '
     END IF
  ELSE IF (alf > zero) THEN
     alf=zero
     PRINT *,' '
     PRINT *,'No Gustafsson modification.'
     PRINT *,' '
  END IF

  DO j=1,nblock

     dd = d(:,:,j)

!!! call timing_start('sweep1')

     ! D = D - E.F(j-1) - E1.F1(j-M1) - E2.F2(j-M2)

     IF (j>1 )&
          CALL BLAS_GEMM('n','n',N,N,N,-one, e(:,:,j),N, f(:,:,j- 1),N,one,dd,N)
     IF (j>M1)&
          CALL BLAS_GEMM('n','n',N,N,N,-one,e1(:,:,j),N,f1(:,:,j-M1),N,one,dd,N)
     IF (j>M2)&
          CALL BLAS_GEMM('n','n',N,N,N,-one,e2(:,:,j),N,f2(:,:,j-M2),N,one,dd,N)

!!! call timing_stop('sweep1')

     IF (alf<zero) THEN

!!! call timing_start('gustafsson')

        ! Relaxed Gustafsson modification

        ! D = D + alf*( E2.F(j-M2) + E2.F1(j-M2) + E1.F(j-M1) + E1.F2(j-M1)
        !             + E.F1(j-1)  + E.F2(j-1) )

        IF (j>M2) THEN
           CALL BLAS_GEMM('n','n',N,N,N,alf,e2(:,:,j),N, f(:,:,j-M2),N,one,dd,N)
           CALL BLAS_GEMM('n','n',N,N,N,alf,e2(:,:,j),N,f1(:,:,j-M2),N,one,dd,N)
        END IF
        IF (j>M1) THEN
           CALL BLAS_GEMM('n','n',N,N,N,alf,e1(:,:,j),N, f(:,:,j-M1),N,one,dd,N)
        END IF
        IF (j>M1.and.j-M1<=nblock-M2) THEN
           CALL BLAS_GEMM('n','n',N,N,N,alf,e1(:,:,j),N,f2(:,:,j-M1),N,one,dd,N)
        END IF
        IF (j>1.and.j-1<=nblock-M2) THEN
           CALL BLAS_GEMM('n','n',N,N,N,alf, e(:,:,j),N,f2(:,:,j-1 ),N,one,dd,N)
        END IF
        IF (j>1.and.j-1<=nblock-M1) THEN
           CALL BLAS_GEMM('n','n',N,N,N,alf, e(:,:,j),N,f1(:,:,j-1 ),N,one,dd,N)
        END IF

!!! call timing_stop('gustafsson')

     END IF

!!! call timing_start('sweep2')

     ! Compute the LU-decomposition of d(j):
     !
     !   DD = LU
     !   F2 = D^{-1}.F2, F1 = D^{-1}.F1, F = D^{-1}.F

     ! Calculate and save D^{-1} into d
     !
     ! d=I, solve dd.d = I

     CALL LAPACK_GETRF( N, N, dd, N, pivot, INFO )
     d(:,:,j)=zero
     do i=1,N
        d(i,i,j)=one
     end do
     CALL LAPACK_GETRS('n',N,N,dd,N,pivot,d(:,:,j),N,INFO)

     IF (j   < nblock)then
        dd=f(:,:,j)
        CALL BLAS_GEMM('n','n',N,N,N,one,d(:,:,j),N,dd,N,zero,f(:,:,j),N)
     end IF
     IF (j+M1<=nblock)then
        dd=f1(:,:,j)
        CALL BLAS_GEMM('n','n',N,N,N,one,d(:,:,j),N,dd,N,zero,f1(:,:,j),N)
     end IF

     IF (j+M2<=nblock)then
        dd=f2(:,:,j)
        CALL BLAS_GEMM('n','n',N,N,N,one,d(:,:,j),N,dd,N,zero,f2(:,:,j),N)
     end IF

!!! call timing_stop('sweep2')

  ENDDO

  ! DEBUG
  ! write(*,*)'D (140)  =',(( d(i,k,140),i=1,N),k=1,N)
  ! write(*,*)'E (140)  =',(( e(i,k,140),i=1,N),k=1,N)
  ! write(*,*)'E1(140)  =',((e1(i,k,140),i=1,N),k=1,N)
  ! write(*,*)'E2(140)  =',((e2(i,k,140),i=1,N),k=1,N)
  ! write(*,*)'F (140)  =',(( f(i,k,140),i=1,N),k=1,N)
  ! write(*,*)'F1(140)  =',((f1(i,k,140),i=1,N),k=1,N)
  ! write(*,*)'F2(140)  =',((f2(i,k,140),i=1,N),k=1,N)

  ! Deallocate arrays that were "Automatic"
  deallocate(dd)
  deallocate(pivot)

  call timing_stop('precond')

END SUBROUTINE prehepta

!=============================================================================
SUBROUTINE Uhepta(inverse,nblock,N,M1,M2,x,f,f1,f2)

  ! G. Toth, 2001

  ! This routine multiplies x with the upper triagonal U or U^{-1}
  ! which must have been constructed in subroutine prehepta.

  use ModPrecond

  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: inverse
  INTEGER, INTENT(IN) :: N,M1,M2,nblock
  REAL, INTENT(INOUT) :: x(N,nblock)
  REAL, INTENT(IN)    :: f(N,N,nblock), f1(N,N,nblock), f2(N,N,nblock)

  !=======================================================================
  !     Description of arguments:
  !=======================================================================
  !
  ! inverse: logical switch
  !          Multiply by U^{-1} if true, otherwise multiply by U
  !
  ! nblock:  Number of diagonal blocks.
  ! N:       the size of the blocks.
  ! M1:      distance of blocks to the main diagonal blocks.
  !          set M1=nblock for block tri-diagonal matrices!
  ! M2:      distance of outer-most blocks to main diagonal blocks.
  !          set M2=nblock for block tri- and penta-diagonal matrices.
  !
  ! x:       On input, the vector to be multiplied with U or U^{-1}.
  !          On output, the result of U.x or U^{-1}.x
  !
  ! f, f1, f2
  !           The matrix U is assumed to be in three block diagonals
  !
  !           The blocks are stored as follows:
  !           f(j): j=1..nblock-1   super diagonal blocks.
  !
  !           f1(j): j=1..nblock-M1 blocks in the upper-triangular part with
  !                                 distance M1 from the main diagonal. 
  !                                 Use scalar for block tri-diagonal matrix!
  !
  !           f2(j): j=1..nblock-M2 blocks in the upper-triangular part with 
  !                                 distance M2 from the main diagonal.
  !                                 Use scalar for block tri- and penta-diagonal
  !                                 matrix!
  !
  ! It is assumed that the blocks are not very sparse, so the sparsity pattern
  ! within the separate blocks is not exploited. For example, the (i,k)-element 
  ! of the j-th block on the super diagonal is stored in f(i,k,j).

  INTEGER :: j
  REAL, PARAMETER :: one=1.0

  ! External function
  !
  ! DGEMV,   BLAS level two Matrix-Vector Product.
  !          See 'man DGEMV' for description.
  !-----------------------------------------------------------------------
  !
  !
  call timing_start('Uhepta')

  if(N>20)then
     ! BLAS VERSION
     IF(inverse)THEN
        !  x' := U^{-1}.x = x - F.x'(j+1) - F1.x'(j+M1) - F2.x'(j+M2)
        DO j=nblock-1,1,-1
           CALL BLAS_GEMV('n',N,N,-one, f(:,:,j),N,x(:,j+1 ),1,one,x(:,j),1)
           IF (j+M1<=nblock) &
                CALL BLAS_GEMV('n',N,N,-one,f1(:,:,j),N,x(:,j+M1),1,one,x(:,j),1)
           IF (j+M2<=nblock) &
                CALL BLAS_GEMV('n',N,N,-one,f2(:,:,j),N,x(:,j+M2),1,one,x(:,j),1)
        ENDDO
     ELSE
        !  x := U.x = x + F.x(j+1) + F1.x(j+M1) + F2.x(j+M2)
        DO j=1,nblock-1
           CALL BLAS_GEMV('n',N,N,one, f(:,:,j),N,x(:,j+1 ),1,one,x(:,j),1)
           IF (j+M1<=nblock) &
                CALL BLAS_GEMV('n',N,N,one,f1(:,:,j),N,x(:,j+M1),1,one,x(:,j),1)
           IF (j+M2<=nblock) &
                CALL BLAS_GEMV('n',N,N,one,f2(:,:,j),N,x(:,j+M2),1,one,x(:,j),1)
        END DO
     END IF
  else
     ! F90 VERSION
     if(inverse)then
        !  x' := U^{-1}.x = x - F.x'(j+1) - F1.x'(j+M1) - F2.x'(j+M2)
        do j=nblock-1,1,-1
           !  x' := U^{-1}.x = x - F.x'(j+1) - F1.x'(j+M1) - F2.x'(j+M2)
           if (j+M2<=nblock) then
              x(:,j) = x(:,j) - matmul( f(:,:,j),x(:,j+1 )) &
                   - matmul(f1(:,:,j),x(:,j+M1)) &
                   - matmul(f2(:,:,j),x(:,j+M2))
           else if(j+M1<=nblock) then
              x(:,j) = x(:,j) - matmul( f(:,:,j),x(:,j+1 )) &
                   - matmul(f1(:,:,j),x(:,j+M1))
           else
              x(:,j) = x(:,j) - matmul(f(:,:,j),x(:,j+1))
           end if
        end do
     else
        !  x := U.x = x + F.x(j+1) + F1.x(j+M1) + F2.x(j+M2)
        do j=1,nblock-1
           if (j+M2<=nblock) then
              x(:,j) = x(:,j) + matmul( f(:,:,j),x(:,j+1 )) &
                   + matmul(f1(:,:,j),x(:,j+M1)) &
                   + matmul(f2(:,:,j),x(:,j+M2))
           else if (j+M1<=nblock) then
              x(:,j) = x(:,j) + matmul( f(:,:,j),x(:,j+1 )) &
                   + matmul(f1(:,:,j),x(:,j+M1))
           else
              x(:,j) = x(:,j) + matmul(f(:,:,j),x(:,j+1))
           end if
        end do
     end if
  end if

  call timing_stop('Uhepta')

END SUBROUTINE Uhepta

!=============================================================================
SUBROUTINE Lhepta(nblock,N,M1,M2,x,d,e,e1,e2)

  ! G. Toth, 2001
  !
  ! This routine multiplies x with the lower triangular matrix L^{-1},
  ! which must have been constructed in subroutine prehepta.

  use ModPrecond

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N,M1,M2,nblock
  REAL, INTENT(INOUT) :: x(N,nblock)
  REAL, INTENT(IN), DIMENSION(N,N,nblock) :: d,e,e1,e2
  !
  !=======================================================================
  !     Description of arguments:
  !=======================================================================
  !
  ! nblock:  Number of diagonal blocks.
  ! N:        the size of the blocks.
  ! M1:       distance of blocks to the main diagonal blocks.
  !           Set M1=nblock for block tri-diagonal matrix!
  ! M2:       distance of outer-most blocks to main diagonal blocks.
  !           Set M2=nblock for block tri- and penta-diagonal matrices!
  !
  ! x:        On input, the vector to be multiplied with L^{-1}.
  !           On output, the result of L^{-1}.x
  !
  ! d, e, e1, e2
  !           The matrix L is in four block diagonals.
  !
  !           The blocks are stored as follows:
  !           d(j): j=1..nblock     Contains inverse of diagonal of L
  !                                 where L is from the incomplete LU 
  !           e(j): j=2..nblock     sub diagonal blocks.
  !           e1(j): j=M1+1..nblock Blocks in the lower-triangular part with 
  !                                 distance M1 from the main diagonal.
  !                                 Use scalar for block tri-diagonal matrix!
  !           e2(j): j=M2+1..nblock Blocks in the lower-triangular part with 
  !                                 distance M2 from the main diagonal.
  !                                 Use scalar for block tri- and penta-diagonal
  !                                 matrices!
!!! automatic array
!!$  REAL    :: work(N)
  real, dimension(:), allocatable :: work

  INTEGER :: j, iError
  REAL, PARAMETER :: zero=0.0, one=1.0

  ! External subroutine
  !
  ! DGEMV,   BLAS level two Matrix-Vector Product.
  !          See 'man DGEMV' for description.
  !
  !---------------------------------------------------------------------------


  ! x' = L^{-1}.x = D^{-1}.(x - E2.x'(j-M2) - E1.x'(j-M1) - E.x'(j-1))

  call timing_start('Lhepta')

  ! Allocate arrays that were "Automatic"
  allocate(work(N), stat=iError); call alloc_check(iError,"precond:work")

  if(N>20)then
     ! BLAS VERSION
     DO j=1,nblock

        CALL BLAS_GEMV('n', N, N, one, d(:,:,j) ,N, x(:,j), 1, zero, work, 1)
        IF (j>1) CALL BLAS_GEMV('n',N,N,-one,e(:,:,j) ,N,x(:,j-1) ,1,one,work,1)
        IF (j>M1)CALL BLAS_GEMV('n',N,N,-one,e1(:,:,j),N,x(:,j-M1),1,one,work,1)
        IF (j>M2)CALL BLAS_GEMV('n',N,N,-one,e2(:,:,j),N,x(:,j-M2),1,one,work,1)
        call BLAS_COPY(N,work,1,x(:,j),1)

     ENDDO
  else
     ! F90 VERSION
     do j=1,nblock
        work = x(:,j)
        if (j>M2) then
           work = work                        &
                - matmul( e(:,:,j),x(:,j-1 )) &
                - matmul(e1(:,:,j),x(:,j-M1)) &
                - matmul(e2(:,:,j),x(:,j-M2))
        else if (j>M1) then
           work = work                        &
                - matmul( e(:,:,j),x(:,j-1 )) &
                - matmul(e1(:,:,j),x(:,j-M1))
        else if (j>1) then
           work = work                        &
                - matmul( e(:,:,j),x(:,j-1 ))
        end if
        x(:,j) = matmul( d(:,:,j),work)
     end do
  end if

  ! Deallocate arrays that were "Automatic"
  deallocate(work)

  call timing_stop('Lhepta')

END SUBROUTINE Lhepta
