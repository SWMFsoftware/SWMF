!^CFG COPYRIGHT UM
!^CFG FILE IMPLICIT
module ModPrecond
  implicit none

  interface blas_copy

     subroutine dcopy(n,dx,incx,dy,incy)
       real*8   :: dx(*),dy(*)
       integer  :: n,incx,incy
     end subroutine dcopy

     subroutine scopy(n,dx,incx,dy,incy)
       real*4   :: dx(*),dy(*)
       integer  :: n,incx,incy
     end subroutine scopy

  end interface

  interface blas_gemv

     subroutine dgemv ( TRANS, M, N, ALPHA, A, LDA, X, INCX, &
          BETA, Y, INCY )
       REAL*8           :: ALPHA, BETA
       INTEGER          :: INCX, INCY, LDA, M, N
       CHARACTER*1      :: TRANS
       REAL*8           :: A( LDA, N ), X( N ), Y( N )
     end subroutine dgemv

     subroutine sgemv ( TRANS, M, N, ALPHA, A, LDA, X, INCX, &
          BETA, Y, INCY )
       REAL*4           :: ALPHA, BETA
       INTEGER          :: INCX, INCY, LDA, M, N
       CHARACTER*1      :: TRANS
       REAL*4           :: A( LDA, N ), X( N ), Y( N )
     end subroutine sgemv

  end interface

  interface blas_gemm

     subroutine dgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, &
          BETA, C, LDC )
       CHARACTER*1 ::        TRANSA, TRANSB
       INTEGER     ::       M, N, K, LDA, LDB, LDC
       REAL*8 ::  ALPHA, BETA
       REAL*8 ::  A( LDA, N ), B( LDB, N ), C( LDC, N )
     end subroutine dgemm

     subroutine sgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, &
          BETA, C, LDC )
       CHARACTER*1 ::        TRANSA, TRANSB
       INTEGER     ::       M, N, K, LDA, LDB, LDC
       REAL*4   ::  ALPHA, BETA
       REAL*4   ::  A( LDA, N ), B( LDB, N ), C( LDC, N )
     end subroutine sgemm

  end interface

  interface lapack_getrf

     subroutine dgetrf( M, N, A, LDA, IPIV, INFO )
       integer           :: INFO, LDA, M, N
       integer           :: IPIV(N)
       real*8            :: A(LDA,N)
     end subroutine dgetrf

     subroutine sgetrf( M, N, A, LDA, IPIV, INFO )
       integer           :: INFO, LDA, M, N
       integer           :: IPIV(N)
       real*4            :: A(LDA,N)
     end subroutine sgetrf

  end interface

  interface lapack_getrs

     subroutine dgetrs( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
       CHARACTER        :: TRANS
       INTEGER          :: INFO, LDA, LDB, N, NRHS
       INTEGER          :: IPIV(N)
       real*8           :: A(LDA,N), B(LDB,N)
     end subroutine dgetrs

     subroutine sgetrs( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
       CHARACTER        :: TRANS
       INTEGER          :: INFO, LDA, LDB, N, NRHS
       INTEGER          :: IPIV(N)
       real*4           :: A(LDA,N), B(LDB,N)
     end subroutine sgetrs

  end interface

end module ModPrecond
