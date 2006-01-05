    SUBROUTINE Gmresm (Matvec, N, x0, b, tol, iter_max)
    USE Rcm_variables, ONLY : iprec, rprec, n_gc
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: n
    REAL (rprec), INTENT (IN) :: b(n), tol
    REAL (rprec), INTENT (IN OUT) :: x0(n)
    INTEGER (iprec), INTENT (IN) :: iter_max

    INTERFACE
      SUBROUTINE Matvec (x,y,n)
      USE Rcm_variables, ONLY : iprec, rprec
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: n
      REAL (rprec), INTENT (IN) :: x(n)
      REAL (rprec), INTENT (OUT):: y(n)
      END SUBROUTINE Matvec
    END INTERFACE
!
!   Linear system A*X=B is solved iteratively by GMRES(m) algorithm
!   (Generalized Minimized Residuals with restarts).
!   Matrix A is only needed by routine MATVEC and is accessed through a module.
!   Since A is sparse, it is stored in the CRS format (compressed-row).
!   Preconditioner accesses matrix through the main module too.
!
!
!   Local variables:
!       nmax  : max # of grid points to treat
!       M     : max # of Krylov vectors (when exceeded, goes into restart)
!       PIVOTS: vector with inverses of diagonal of the preconditioner matrix
!       TOL   : relative error. Empirically, 1e-5 is equiv. to
!               sum of resids < 1 Volt
!       H     : Hessenberg matrix holding dot products
!       B     : holds right-hand side vector of linear system to be solved
!       CS, SN: vectors holding Givens rotations parameters
!       X0    : vector holding initial approx. on entry and solution on exit
!       RESID : vector of residuals (b-A*x)
!       AX, W : work vectors
!       Y     : vector holding coefficients of the solution expanded in the
!               Krylov vectors
!       S     : initially unit vector, then rotated by Givens rotations
!
!  This algorithm uses preconditioning based on incomplete LU factorization.
!  Namely, it uses D-ILU type.
!
    INTEGER (iprec), PARAMETER :: M = 300
    REAL (rprec), POINTER :: window (:)
    REAL (rprec), TARGET  :: h (m+1,m), s (m+1)
    REAL (rprec)          :: sn (m), cs (m), y (m), x (n,m+1)
    REAL (rprec)          :: ax (n), resid (n), w (n)
    INTEGER (iprec) :: i, jkryl, iter
    REAL (rprec) :: bnorm, rnorm, relerr, v_tmp (2)
!
!
!
    CALL Gmresm_compute_DILU ()
!
!
    bnorm = SQRT ( DOT_PRODUCT (b, b))
    IF (ABS(bnorm) < TINY(1.0_rprec)) bnorm = 1.0_rprec
!
!
    Restart_loop: DO iter = 1, iter_max   !  Begin GMRES(m) iterations (restarts):
!
!
!      ... Compute the norm of initial residuals:
!
       CALL Matvec (x0, ax, SIZE(ax))
       ax     = b - ax
       CALL Gmresm_msolve (n, ax, resid)
       rnorm  = SQRT (DOT_PRODUCT (resid, resid ))
       relerr = rnorm / bnorm
       IF (relerr < TOL ) RETURN     ! x0 is already a solution
!
!      .. Set 1st Krylov vector to R/||R||:
!
       x (:, 1) = resid / rnorm
!
!      .. Set up unit vector E1 of length RNORM:
!
       s (1) = rnorm
       s (2:m+1) = 0.0
!
!      .. Loop to generate orthonormal vectors in Krylov subspace:
!
       Iterate_loop: DO jkryl = 1, M
!
!         ... Compute A*X(Jkryl) and solve M*w=A*X(kryl) for w:
!
          CALL Matvec (x(:,jkryl), ax, SIZE(ax))
          CALL Gmresm_msolve (n, ax, w)
!
!         ... Form J-th column of H-matrix and X (Jkryl+1)
!                (modified Gramm-Schmidt process):
!
          DO i = 1, jkryl
             H (i,jkryl)   = DOT_PRODUCT ( w , x (:,i) )
             w             = w  - h (i,jkryl) * x (:,i)
          END DO
          h (jkryl+1,jkryl)  = SQRT (DOT_PRODUCT (w, w))
          x (:, jkryl+1)     = w / h(jkryl+1,jkryl)
!
!
!         .. Update QR-factorization of H. For that, 
!         .... first, apply 1, ..., (Jkryl-1)th rotations
!              to the new (Jkryl-th) column of H:
!
          DO i = 1, Jkryl-1
             v_tmp = h (i:i+1, jkryl)
             CALL Gmresm_rotate_vector (v_tmp, cs(i), sn(i), h(i:i+1,jkryl))
          END DO
!
!         .... second, compute the Jkryl-th rotation that
!              will zero H (jkryl+1,jkryl):
!
          window => h (jkryl:jkryl+1, jkryl)
          CALL Gmresm_Get_rotation ( window, cs (jkryl), sn (jkryl) )
!
!         .... third, apply Jkryl-th rotation to Jkryl-th column of H
!              and to S (rhs):
!
          window => h (jkryl:jkryl+1, jkryl)
          v_tmp = h (jkryl:jkryl+1, jkryl)
          CALL Gmresm_rotate_vector (&
               v_tmp, cs(jkryl), sn (jkryl), h(Jkryl:jkryl+1,jkryl))
          h (jkryl+1,jkryl) = 0.0
!
          window => s (jkryl : jkryl+1)
          v_tmp = s(jkryl:jkryl+1)
          CALL Gmresm_rotate_vector (&
               v_tmp, cs(jkryl), sn (jkryl), s(Jkryl:jkryl+1))
!
!
!         .. Approximate the norm of current residual:
!
          relerr = ABS (s (jkryl+1)) / bnorm
!
          IF (relerr < TOL) THEN
!
!            .. Residual is small, compute solution, exit:
!
             CALL Gmresm_u_solve (h, s(1:m), m, jkryl, y)
             DO i = 1, Jkryl
                x0  = x0 + y(i)* X(:,i)
             END DO
             EXIT restart_loop
!
          END IF
!
       END DO Iterate_loop
!
!
!      We got here because after a maximum number of Krylov vectors
!      was reached, approximated norm of residual was not small enough.
!      However, need to compute approx solution and check the actual norm
!      of residual (because the approx. norm may not be accurate due to
!      round offs):
!
       CALL Gmresm_u_solve (h, s, m, m, y)
       DO i = 1, m
          x0 = x0 + y(i)* X(:,i)
       END DO
! 
       CALL Matvec (x0, resid, SIZE(resid))
       resid = b - resid
       rnorm   = SQRT (DOT_PRODUCT (resid, resid))
       relerr = rnorm / bnorm
!
!      .. If the actual norm of residual is indeed small, exit:
! 
       IF (relerr < TOL) EXIT restart_loop
!
!      .. If not, continue by restarting...
!
    END DO restart_loop
!
!
!   Finished GMRES(m) loop. Got here either because solution was found, 
!   or because maximum number of iterations was exceeded. Check for this:
!     
    IF (relerr >= TOL) THEN
       call CON_stop('ERROR in IM/RCM2/src/gmresm.f90: '// &
            'no convergence in GMRES(m)')
    END IF
!
    RETURN
!
END SUBROUTINE Gmresm
!
!
!
!
      SUBROUTINE Gmresm_compute_DILU ()
      USE rcm_variables, ONLY : iprec, rprec, &
                                nij_pde, &
                                a_mtrx => a_mtrx_pde, &
                                row_ptr => row_ptr_pde, &
                                i_column => i_column_pde, &
                                diag_ptr => diag_ptr_pde, &
                                pivots => pivots_pde
      IMPLICIT NONE
!_____________________________________________________________________________
!     Compute the preconditioner M. Matrix A is split as A = L_a + D_a + U_a
!     (strictly-lower triangular, diagonal and strictly-upper triangular).
!     Then M = L * U = (D + L_a) * D^(-1) * (D + U_a), so only need to find
!     and store D (one diagonal). In fact, PIVOTS holds inverses of D since
!     will divide by them later. D-ILU preconditioner M is kept in PIVOTS.
!     All structures are accessed from the host subroutine.
!
!     This subroutine only modifies (computes) PIVOTS (1:nij).
!_____________________________________________________________________________
!
      INTEGER (iprec) :: irow, jcol, krow
      REAL    (rprec) :: element
      LOGICAL         :: found
!
      pivots  = a_mtrx (diag_ptr )
      DO irow = 1, nij_pde
         pivots (irow) = 1.0_rprec / pivots (irow)
         DO jcol = diag_ptr(irow)+1, row_ptr(irow+1)-1
            found = .FALSE.
            DO krow = row_ptr (i_column (jcol)), diag_ptr (i_column (jcol)) - 1
               IF (i_column (krow) == irow) THEN
                  found = .TRUE.
                  element = a_mtrx (krow)
               END IF
            END DO
            IF (found) THEN
               pivots (i_column (jcol)) = &
                  pivots (i_column (jcol)) - element*pivots(irow)*a_mtrx(jcol)
            END IF
         END DO
      END DO
      RETURN
      END SUBROUTINE Gmresm_Compute_DILU
!
!
!
!
      SUBROUTINE Gmresm_Rotate_vector(v_in, cos_theta, sin_theta, v_out)
      USE Rcm_variables, ONLY : rprec
      IMPLICIT NONE
      REAL (rprec), INTENT (IN) :: v_in (2), cos_theta, sin_theta
      REAL (rprec), INTENT (OUT) :: v_out (2)
!___________________________________________________________________________
!     Apply a plane (Givens) rotation with cos_theta, sin_theta) to a vector.
!     Rotation acts on only two elements of the vector, X and Y.
!___________________________________________________________________________
!
      REAL (rprec):: temp
!
      temp                = cos_theta * v_in (1) - sin_theta * v_in (2)
      v_out(2)    = sin_theta * v_in (1) + cos_theta * v_in (2)
      v_out(1)    = temp
      RETURN
      END SUBROUTINE Gmresm_Rotate_vector
!
!
!
      SUBROUTINE Gmresm_Msolve (n, x, y)
      USE Rcm_variables, ONLY : iprec, rprec, &
                                row_ptr => row_ptr_pde, &
                                a_mtrx => a_mtrx_pde, &
                                diag_ptr => diag_ptr_pde, &
                                i_column => i_column_pde, &
                                pivots => pivots_pde
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: n
      REAL (rprec), INTENT (IN) :: x (n)
      REAL (rprec), INTENT (OUT) :: y(n)
!_____________________________________________________________________________
!     This subroutine solves the system L *  U * y = x, where
!     M = L * U = (D + L_a) * (I + D^(-1) * U_a) is the D-ILU preconditioner.
!     Matrices L_a and U_a are strictly lower and strictly upper triangular,
!     so that A = D_a + L_a + U_a, and D comes from incomplete LU factorization
!     when computing preconditioner. Solution proceeds in the regular way by
!     forward- and then back-substition (solving L*z=x, then U*y=z).
!     A is in the compressed-row-storage format (A_MTRX, I_COLUMN, ROW_PTR).
!     Diagonal matrix D (in fact, its inverse) is stored in the PIVOTS:
!     PIVOTS(i)=1/D(i,i), and DIAPTR vector holds locations of d_i_i in amatrx.
!     Since A_MTRX and PIVOTS do not change in the potential solver once is
!     has been called, we access them by host association from GMRESM. Only
!     vector X changes from invocation to invocaton of MSOLVE, and we pass it
!     as an argument.
!**** NOTE: book by Barrett et al ("templates ...") has an error in the back-
!           substitution algorithm (p.73). Here I do it correctly.
!_____________________________________________________________________________
!
      INTEGER (iprec) :: i, j
      REAL (rprec)    :: tmp
!
      IF (SIZE(x) /= SIZE(diag_ptr) ) STOP 'ERROR IN GMRESM_MSOLVE'
!
      DO i = 1, n
         tmp = 0.0
         DO j = row_ptr (i), diag_ptr (i) - 1
            tmp = tmp + a_mtrx (j) * y (i_column (j))
         END DO
         y (i) = pivots (i) * (x(i) - tmp)
      END DO
      DO i = SIZE(x), 1, -1
         tmp = 0.0
         DO j = diag_ptr (i) + 1, row_ptr (i+1)-1
            tmp = tmp + a_mtrx(j) * y (i_column(j))
         END DO
         y (i) = y (i) - pivots (i) * tmp
      END DO
      RETURN
      END SUBROUTINE Gmresm_Msolve
!
!
!
!
      SUBROUTINE Gmresm_u_solve (a, b_rhs, nmax, n, xs) 
      USE rcm_variables, ONLY : iprec, rprec
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: n, nmax
      REAL (rprec), INTENT (IN) :: a (nmax+1,nmax), b_rhs (nmax)
      REAL (rprec), INTENT (OUT) :: xs (n)
!
!     Given an upper triangular matrix A and right-hand side vector B(n),
!     solves linear system A*x=b. A is a Hessenberg matrix (nmax+1 by nmax),
!     but only n by n section is used.
!
      INTEGER (iprec) :: j
!
      IF (n > nmax .OR. n < 1) STOP 'PROBLEM 2 IN SOLVETR'
!
      xs = b_rhs(1:n)
      DO j = N, 1, -1
         xs (j)     = xs (j) / a (j,j)
         xs (1:j-1) = xs (1:j-1) - xs (j) * a(1:j-1,j)
      END DO
      RETURN
      END SUBROUTINE Gmresm_u_solve
!
!
!
      SUBROUTINE Gmresm_Get_rotation ( vector_in, cos_theta, sin_theta )
      USE rcm_variables, ONLY : rprec
      IMPLICIT NONE
      REAL (rprec), INTENT (IN)  :: vector_in (2)
      REAL (rprec), INTENT (OUT) :: cos_theta, sin_theta
!_____________________________________________________________________________
!     Compute a Givens (plane) rotation that will act on 2 elements of a
!     vector, A and B, and will zero B. Returns cosine and sine of THETA, the
!     angle of rotation. The transformation is
!
!        A_prime = A * COS(theta) - B * SIN(theta)
!        B_prime = A * CIN(theta) - B * COS(theta)
!
!     In matrix-vector terms,
!        X = (...... A ..... B .....)^T,
!        T =
!        X_prime = T * X = (..... A_prime ...... 0 .....)^T,
!        only 2 elements of X are changed by the rotation.
!_____________________________________________________________________________
!
      REAL (rprec) :: temp
!
      IF ( ABS(vector_in (2)) < TINY(1.0_rprec) ) THEN
         cos_theta = 1.0_rprec
         sin_theta = 0.0
      ELSE IF ( ABS ( vector_in(2) ) > ABS ( vector_in (1) ) ) THEN
         temp = -vector_in (1) / vector_in (2)
         sin_theta = 1.0_rprec / SQRT( 1.0_rprec + temp**2 )
         cos_theta = temp * sin_theta
      ELSE
         temp = -vector_in (2) / vector_in (1)
         cos_theta = 1.0_rprec / SQRT( 1.0_rprec + temp**2 )
         sin_theta = temp * cos_theta
      END IF
!
      RETURN
      END SUBROUTINE Gmresm_Get_rotation
!
!
!
      SUBROUTINE Gmresm_Matvec (x, y, n)
      USE Rcm_variables, ONLY : iprec, rprec, &
                                a_mtrx_pde, i_column_pde, row_ptr_pde
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: n
      REAL (rprec), INTENT (IN)    :: x (n)
      REAL (rprec), INTENT (OUT)   :: y(n)
!____________________________________________________________________________
!     subroutine to form matrix-vector product. Matrix A of size NxN
!     is assumed to be sparse and is stored in the compressed-row (CRS) format.
!     We compute y = A*x, where y and x are both vectors of length N.
!____________________________________________________________________________
!
      INTEGER (iprec) :: i, j
!
      DO i = 1, n
         y (i) = 0.0
         DO j = row_ptr_pde (i), row_ptr_pde (i+1)-1
            y (i) = y (i) + a_mtrx_pde(j) * x (i_column_pde(j) )
         END DO
      END DO
      RETURN
      END SUBROUTINE Gmresm_Matvec
!
!
!
      SUBROUTINE Gmresm_unwrap_pde_solution ()
      USE Rcm_variables, ONLY : iprec, isize, jsize, &
                                ni_pde, v, X0_pde, n_gc
      IMPLICIT NONE
!
      INTEGER (iprec) :: i, j, krow
!
      DO j = 1, jsize
      DO i = 1, isize
         krow = ni_pde*(j-1) + i
         v (i,j) = X0_pde (krow)
      END DO
      END DO
      CALL Wrap_around_ghostcells (v, isize, jsize, n_gc)
!
      RETURN
      END SUBROUTINE Gmresm_unwrap_pde_solution
!
!
!
  SUBROUTINE Define_pde_matrix ()
  USE Rcm_variables, ONLY : iprec, rprec, isize, jsize, ncoeff, &
                            c_pde, v, imin_j, &
                            ni_pde, nj_pde, &
                            amtrx => a_mtrx_pde, &
                            bmtrx => b_mtrx_pde, &
                            iclmn => i_column_pde, &
                            rowpt => row_ptr_pde, &
                            dmtrx => diag_ptr_pde, &
                            xinit => x0_pde
  IMPLICIT NONE
!_____________________________________________________________________________
!     This subroutine is used in conjuntion with a general-purpose linear 
!     system solver to solve the elliptic PDE for the electrostatic potential.
!     Subroutine returns:
!       matrix A stored in 3 vectors AMTRX, ICLMN, ROWPT in CRS format
!         (also, DIAGPT holds locations for diagonal elements),
!       right-hand side vector BMTRX,
!       initial approximation to solution in vector XINIT.
!_____________________________________________________________________________
!
   INTEGER (iprec):: i, j, L, krow
!
!
!  This subroutine takes the coefficients of the RCM difference
!  equations approximating the MI-coupling PDE, and reformulates these
!  equations as to cast them into a general linear system A*X=B, where A is
!  an NxN square matrix, X is the unknown vector whose elements are
!  unknown values of the potential on grid points V(i,j), and B is the
!  right-hand-side vector.
!  This reformulation requires: 
!  (1) to number all grid points sequentially into a 1-dimensional sequence,
!  (2) to form A from c1-c4 and B from c5 RCM coefficients.
!  As A is going to be sparse, an additional task is to store (encode) A
!  in the Compressed-Row-Storage (CRS) format for using in the potential solver.
!
!  The boundary issues are dealt with in the following way:
!  IMIN_J array holds the first I-point inside the modeling region for given J.
!  The routine will treat the whole region (modulo ghost-cells), but outside 
!  the boundary, potential is not modified, so the linear equations for those
!  grid points are set to V(i,j) = V_in(i,j), where V_in is the initial values.
!  Matrix A is stored in 1-dim REAL array AMTRX and two INTEGER 1-dim
!     arrays ICLMN and ROWPT. We simply go along each row of A starting with
!     the 1st row, then 2nd, etc, and for each non-zero element a(p,q), we
!     write AMTRX(L)=a(p,q), ICLMN(L)=q, and L-index numbers those non-zero
!     elements sequentially. ROWPT(p) has the L-index of where p-th row
!     starts in AMTRX.
!
!  1. Numbering grid points into a 1-dim. sequence.
!     Imagine the RCM grid as extending vertically in I from I=1 (highest lat,
!     top) to I=IDIM (lowest lat., bottom) and horizontally in J from J=j1
!     (noon, left) to J=J2 (last point before noon, right). 
!     The rectangular region of the grid we treat has the size NI by NJ.
!     Each point (i,j) has number KROW (so that coefficients of the difference
!     equation on that point are on the krow-th row of A).
!
! Make sure points outside the boundary are not used:
!
  DO j = 1, jsize
     c_pde (1:4, 1:imin_j(j)-1, j) = 0.0_rprec
  END DO
  c_pde (1,isize,:) = 0.0_rprec
  c_pde (2, 1, :)   = 0.0_rprec
  DO i = 1, isize
  DO j = 1, jsize
     IF (i< imin_j(j))then
         c_pde(5,i,j) = v (i,j)
     end if
  END DO
  END DO
!
!
! RCM difference equation looks like:
!       V(i,j)=c(1,i,j)*V(i+1,j)+c(2,i,j)*V(i-1,j)+c(3,i,j)*V(i,j+1)+
!              c(4,i,j)*V(i,j-1)+c(5,i,j);
! except that if i=idim, then there is no C(1) term, and
! if i=min_j(j), there is no c(2) term.
! This is how matrix A looks like: all rows except for the first and last
! ones have 5 coefficients; 1st and last rows have 4 coefficients, but we
! define two more "fake" coefficients and set them to zero just for symmetry.
!
! The first NI matrix rows (j=j1) look like:
!
!            1 -c1   0 .............-c3........ . . .      -c4...........0
!          -c2  1  -c1 ............. 0  -c3.... . . .       0  -c4.......0
!            0  -c2  1 -c1.................-c3. . . .       ......-c4....0
!           . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!                   -c2 1 0 ..................... -c3  . .  ..........  -c4
!
! The rest of the matrix except for the last NI rows:
!
!           . . . . . -c4 . . . . . . -c2  1  -c1 . . . . . . -c3 . . . .
!
!  The last NI rows of matrix A (j=j2) look like:
!
!          -c3  0 ..... . . .     -c4............  1  -c1 ................
!           0  -c3..... . . .      0 -c4......... -c2  1  -c1 ............
!           ......-c3.. . . .      .....-c4......  0  -c2  1  -c1 ........
!           . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!           ...............-c3 0  . ............-c4 .................-c2  1
!
!  Notice: there is no c2 for i=1, j=j1,
!          there is no c1 for i=isize, j=j2
!          there is c1 for i=isize for all j but it must be 0
!          there are c2 for i=imin_j but must be zero for all j but must be 0
!- - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - - -
!
!
! Rows 1 through NI_PDE:
  j=1
  DO i = 1, isize
     krow              = ni_pde * (j - 1) + i
     L = ncoeff*(krow-1) + 1 
     bmtrx  (krow) = c_pde (5,i,j)
     Xinit  (krow) = v (i,j)
     rowpt (krow)  = L
     IF (krow == 1) THEN
        amtrx (L:L+4) = (/ +1.0_rprec,     &
                        -c_pde(1,i,j),  &
                        -c_pde(3,i,j),  &
                        -c_pde(4,i,j),  &
                         0.0 /)
        iclmn (L:L+4) = (/ krow, krow + 1, krow + ni_pde, &
                           ni_pde*(nj_pde - 1) + krow, ni_pde*nj_pde /)
        dmtrx (krow) = L 
     ELSE
        amtrx (L:L+4) = (/ -c_pde (2,i,j), &
                           +1.0_rprec,     &
                           -c_pde(1,i,j),  &
                           -c_pde(3,i,j),  &
                           -c_pde(4,i,j) /)
        iclmn (L:L+4) = (/ krow - 1, krow, krow + 1, krow + ni_pde, &
                           ni_pde*(nj_pde - 1) + krow /)
        dmtrx (krow) = L + 1
     END IF
  END DO
!
!
  DO j = 1+1, jsize-1
  DO i = 1, isize
     krow              = ni_pde * (j - 1) + i
     L                 = ncoeff*(krow-1) + 1 
     bmtrx (krow) = c_pde (5,i,j)   ! this will be RHS vector
     Xinit (krow) = v (i,j) ! initial approximation taken from prev solution
     rowpt (krow) = L
     amtrx (L:L+4)= (/ -c_pde (4,i,j), &
                       -c_pde (2,i,j), &
                       + 1.0_rprec,    &
                       -c_pde (1,i,j), &
                       -c_pde (3,i,j)    /)
     iclmn (L:L+4)= (/ krow-ni_pde, krow-1, krow, krow+1, krow+ni_pde /)
     dmtrx (krow) = L + 2
  END DO
  END DO
!
!
  j = jsize
  DO i = 1, isize
     krow              = ni_pde * (j - 1) + i
     L                 = ncoeff*(krow-1) + 1 
     rowpt (krow) = L
     IF (krow == ni_pde*nj_pde) THEN
        amtrx (L:L+4) = (/  0.0, &
                           -c_pde (3,i,j), &
                           -c_pde (4,i,j), &
                           -c_pde (2,i,j), &
                           +1.0_rprec /)
        iclmn (L:L+4) = (/ 1, i, krow - ni_pde, krow - 1, krow /)
        dmtrx (krow) = L + 4
     ELSE
        amtrx (L:L+4) = (/ -c_pde (3,i,j), &
                           -c_pde (4,i,j), &
                           -c_pde (2,i,j), &
                           +1.0_rprec,     &
                           -c_pde (1,i,j)   /)
        iclmn (L:L+4) = (/ i, krow - ni_pde, krow - 1, krow, krow + 1 /)
        dmtrx (krow) = L + 3
     END IF
     bmtrx (krow) = c_pde (5,i,j)   ! this will be RHS vector
     Xinit (krow) = v (i,j)     ! initial approximation taken from prev solution
  END DO
!
!
  rowpt (ni_pde*nj_pde+1) = ncoeff*ni_pde*nj_pde + 1    ! by definition of CRS format.
!
  RETURN
  END SUBROUTINE Define_pde_matrix
