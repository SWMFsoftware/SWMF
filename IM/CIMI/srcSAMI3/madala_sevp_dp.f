
      SUBROUTINE madala_sevp(AX,AY,CX,CY,BB,F,phi,f11_lb,f11_ub)

!      IMPLICIT REAL*8(A-H,O-Z)                                  

      include 'param3_mpi-1.98.inc'
      
      PARAMETER (M=nnx+1, N=nnyt)
!      PARAMETER (NBLK=22, NBLK1=NBLK-1)
!      PARAMETER (NBLK=26, NBLK1=NBLK-1)
!      PARAMETER (NBLK=28, NBLK1=NBLK-1)
!      PARAMETER (NBLK=30, NBLK1=NBLK-1)
      PARAMETER (NBLK=25, NBLK1=NBLK-1)
!      PARAMETER (NBLK=36, NBLK1=NBLK-1)
      PARAMETER (M1=M-1, M2=M-2, N1=N-1, N2=N-2)
      PARAMETER (MN=M*N)

!      PARAMETER (NBLK=1+(N-5+3)/9, NBLK1=NBLK-1)


      dimension NBSIZ2(NBLK),IS(NBLK),IE(NBLK),SUMF(NBLK),IXNDX(M)
      real*8 AX(M,N),AY(M,N),CX(M,N),CY(M,N),BB(M,N)

      real*8 RINV(M2,M2,NBLK),RINV1(M2,M2,NBLK1),RCOR(M,3),RTILDA(M2),
     .       dummy(m2,m2),b1(m2),b2(m2)

      real*8 F(M,N),PHI(M,N),ERR(M,N)
      real*4 F11(M),F1N(M),F21(N),F2M(N)

      real*8 f11_lb(M)
      real*8 f11_ub(M)

!          print *,nblk

C  NBSIZ2 REPRESENTS NUBER OF INTERIOR GRID POINTS IN EACH BLOCK IN X-DIRECTION 
C N2 REPRESENTS NUMBER OF INTERIOR GRID POINTS IN Y-DIRECTION                   
C NBLK REPRESENTS NUMBER OF BLOCKS IN X-DIRECTION                               
C TOTAL NUMBER OF POINTS IN X-DIRECTION(M) SHOULD BE NBLK(NBSIZ2)+2             

C NEUMANN BOUNDARY CONDITIONS IN Y

      A11 = 1.0
      A1N = 0.

!      A11 = 0.
!      A1N = 1.0

c I think the following sets the phi derivative to 0 
c on the upper and lower boundaries

      DO 150 I = 2,M-1
         F11(I) = f11_lb(i)
         F1N(I) = f11_ub(i)
 150  CONTINUE

C CYCLIC BOUNDARY CONDITIONS IN X

      ICYC = 1
      A21 =0
      A2M =0

      do j = 2,n-1
         f21(j) = 0.
         f2m(j) = 0.
      enddo

C
C  FOR INDETERMINATE NEUMANN BOUNDARY CONDITIONS, A BOUNDARY VALUE
C  MUST BE PROVIDED AT A SINGLE BOUNDARY POINT (IX,JX) FOR THE SOLVER TO
C  OBTAIN A SOLUTION. THE BOUNDARY POINT SELECTED IN SEVP IS IX=2, JX=1.
C  FOR THE EFFICIENT CONVERGENCE OF THE SOLUTION, WE USE A ZERO VALUE
C  FOR THE SINGLE BOUNDARY POINT AND AS THE FIRST GUESS.
C
C
C TOLERANCE FOR SOLVER

       TOL = 1.0E-8
C

       IFLG = 1

       CALL SEVP(AX,AY,CX,CY,BB,RINV,RINV1,RCOR,NBSIZ2,IS,IE,IXNDX,
     1      SUMF,DUMMY,RTILDA,PHI,F,ERR,F11,F1N,F21,F2M,M,N,M1,N1,M2,N2,
     2      NBLK,NBLK1,B1,B2,A11,A1N,A21,A2M,ICYC,TOL,IFLG,M,N)
C

!       PRINT 710
! 710   FORMAT(1X,'BLOCK SIZE GIVEN BY NBSIZ2(1)+2 (1ST BLOCK) OR',
!     1    '  NBSIZ2()+3 (OTHER BLOCKS)'/,1X,' NBSIZ2:' )
!       PRINT 711, NBSIZ2
! 711   FORMAT(1X,10I4)
C

C RESIDUAL ERROR OF LAPLACIAN
!       FNORM=FMAX
!       ERMAX=0.0

!       DO 750 I=2,M1
!          DO 750 J=2,N1
!             ERMAX = AMAX1 (ERMAX,ABS(ERR(I,J)/FNORM))
! 750   CONTINUE

!        PRINT 751, ERMAX
!  751   FORMAT(1X,' MAXIMUM RESIDUAL ERROR  
!     1      [(FORCING - LAPLACIAN)/MAX',
!     1    '  FORCING ]  =  ',1PE12.4)
C

      RETURN
      END   
                                                                    
