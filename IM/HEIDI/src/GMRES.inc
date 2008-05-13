  !  SMLIB 1.0b
  !  Copyright (C) 1996 Ernst A. Meese
  !    Refer to the file copyright.doc for details and important disclaimer.
  !
    INTEGER,                              INTENT (IN)    :: m
    REAL(prec), DIMENSION(:),             INTENT (INOUT) :: x
    REAL(prec), DIMENSION(A%N),           INTENT (IN)    :: b
    LOGICAL,                    OPTIONAL, INTENT (IN)    :: init
    REAL(prec),                 OPTIONAL, INTENT (IN)    :: tol, reduction_factor
    INTEGER,                    OPTIONAL, INTENT (IN)    :: min_it, max_it, INFO
    INTEGER,                    OPTIONAL, INTENT (OUT)   :: no_of_its
    REAL(prec),                 OPTIONAL, INTENT (OUT)   :: residual_norm
    TYPE(ConvHist),             OPTIONAL, INTENT (OUT)   :: History
    !---
    INTEGER                           :: ic, min_iter, max_iter, infolevel, i, i1, j, k1, k
    REAL(prec)                        :: toler, norm, t, gam
    REAL(prec), DIMENSION (A % N,m+1) :: v
    REAL(prec), DIMENSION (A % N)     :: z
    REAL(prec), DIMENSION (m+1,m)     :: h
    REAL(prec), DIMENSION (m)         :: c, s
    REAL(prec), DIMENSION (m+1)       :: rs
    !--- 

    CALL Check_Options ()

    IF (infolevel/=NONE) PRINT F1, 'GMRES', toler, min_iter, max_iter

    IC = 0 
    v(:,1) = b - A * x
    DO
       norm = SQRT(DOT_PRODUCT(v(:,1),v(:,1)))
       IF (norm == 0.0_prec) EXIT
       v(:,1) = v(:,1) / norm
       !     ** initialize 1-st term  of rhs of hessenberg system..            
       rs(1) = norm 
       i = 0
       DO
          if (PRESENT(History).and.(i/=0.or.IC==0)) Call SetIterationHistory (history,norm)
          IF (infolevel==ALL) Print F2, ic, norm
          IF (i == m .or. ((norm <= toler .or. IC >= max_iter) .and. IC >= min_iter))  EXIT
          IC = IC + 1; i = i + 1; i1 = i + 1
          CALL solve (PC, v(:,i), z)
          v(:,i1) = A * z
          !-----------------------------------------                              
          !     modified gram - schmidt...                                        
          !-----------------------------------------                              
          DO j=1, i 
             t = DOT_PRODUCT( v(:,j), v(:,i1) )
             h(j,i) = t
             v(:,i1) = v(:,i1) - t * v(:,j)
          END DO
          t = SQRT(DOT_PRODUCT (v(:,i1), v(:,i1)))
          h(i1,i) = t 
          IF ( t /= 0.0_prec) v(:,i1) = v(:,i1) / t
          !                                                                       
          !     done with modified gram schimd and arnoldi step..                 
          !     now  update factorization of h                                   
          !                                                                       
          !--------perfrom previous transformations  on i-th column of h          
          DO k=2,i 
             k1 = k-1 
             t = h(k1,i) 
             h(k1,i) = c(k1)*t + s(k1)*h(k,i) 
             h(k,i) = -s(k1)*t + c(k1)*h(k,i) 
          END DO
          gam = SQRT(h(i,i)**2 + h(i1,i)**2) 
          !                                                                       
          !     if gamma is zero then any small value will do...                  
          !     will affect only residual estimate                                
          !                                                                       
          IF (gam .eq. 0.0_prec) gam = EPSILON(gam) 
          !                                                                       
          !     get  next plane rotation                                          
          !                                                                       
          c(i) = h(i,i)/gam 
          s(i) = h(i1,i)/gam 
          rs(i1) = -s(i)*rs(i) 
          rs(i) =  c(i)*rs(i) 
          !                                                                       
          !     determine residual norm and test for convergence-                
          !                                                                       
          h(i,i) = c(i)*h(i,i) + s(i)*h(i1,i) 
          norm = ABS(rs(i1))
       END DO
       !                                                                       
       !     now compute solution. first solve upper triangular system.        
       !                                                                       
       DO k = i, 1, -1
          rs(k) = ( rs(k) - dot_product(h(k,k+1:i), rs(k+1:i)) ) / h(k,k) 
       END DO
       !
       ! form linear combination of v(*,i)'s to get solution and call preconditioner
       !
       CALL solve (PC, MATMUL(v(:,1:i), rs(1:i)) , z)
       x = x + z
       !                                                                       
       !     restart outer loop  when necessary                                
       !                                                                       
       IF (norm <= toler .or. IC >=max_iter) EXIT 
       !                                                                       
       !     else compute residual vector and continue..                       
       !                                                                       
       DO j = i + 1, 2, -1
          rs(j-1) = -s(j-1)*rs(j) 
       END DO
       rs(2:i+1) = c(1:i) * rs(2:i+1)
       v(:,1) = MATMUL(v,rs)
    END DO

    INCLUDE "check-options.inc"

