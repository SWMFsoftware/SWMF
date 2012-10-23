!This code is a copyright protected software (c) 2002- University of Michigan


! Simple BiCGstab(\ell=1) iterative method
! Modified by G.Toth from the \ell<=2 version written
! by M.A.Botchev, Jan.'98. 
! Further modified for the BATS-R-US code (2000) by G. Toth
!
! This is the "vanilla" version of BiCGstab(\ell) as described
! in PhD thesis of D.R.Fokkema, Chapter 3.  It includes two enhancements 
! to BiCGstab(\ell) proposed by G.Sleijpen and H.van der Vorst in
! 1) G.Sleijpen and H.van der Vorst "Maintaining convergence 
!    properties of BiCGstab methods in finite precision arithmetic",
!    Numerical Algorithms, 10, 1995, pp.203-223
! 2) G.Sleijpen and H.van der Vorst "Reliable updated residuals in
!    hybrid BiCG methods", Computing, 56, 1996, 141-163
!
! {{ This code is based on:
! subroutine bistbl v1.0 1995
!
! Copyright (c) 1995 by D.R. Fokkema.
! Permission to copy all or part of this work is granted,
! provided that the copies are not made or distributed
! for resale, and that the copyright notice and this
! notice are retained.  }}
!
! This subroutine determines the solution of A.QX=RHS, where
! the matrix-vector multiplication with A is performed by
! the subroutine 'proj_matvec'. For symmetric matrix use the more efficient
! proj_cg algorithm!
!
subroutine proj_bicgstab(rhs,qx,iter,tol,typestop,info)
  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,nBLK
  use ModAdvance, ONLY : tmp1_BLK,tmp2_BLK
  use ModProject
  implicit none

  ! Arguments

  real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nBLK), intent(inout) :: rhs
  !        on input:  right-hand side vector.
  !        on output: residual vector.

  real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nBLK), intent(out):: qx
  !       on output: solution vector.

  integer, intent(inout) :: iter
  !       on input:  maximum number of iterations to be performed.
  !       on output: actual  number of iterations done.

  real, intent(inout) :: tol
  !       on input:  required (relative) 2-norm or maximum norm of residual
  !       on output: achieved (relative) 2-norm or maximum norm of residual

  character (len=3), intent(in) :: typestop
  !      Determine stopping criterion (||.|| denotes the 2-norm):
  !      typestop='rel'    -- relative stopping crit.: ||res|| <= tol*||res0||
  !      typestop='abs'    -- absolute stopping crit.: ||res|| <= tol
  !      typestop='max'    -- maximum  stopping crit.: max(abs(res)) <= tol

  ! NOTE for typestop='rel' and 'abs': 
  !            To save computational work, the value of 
  !            residual norm used to check the convergence inside the main 
  !            iterative loop is computed from 
  !            projections, i.e. it can be smaller than the true residual norm
  !            (it may happen when e.g. the 'matrix-free' approach is used).
  !            Thus, it is possible that the true residual does NOT satisfy
  !            the stopping criterion ('rel' or 'abs').
  !            The true residual norm (or residual reduction) is reported on 
  !            output in parameter TOL -- this can be changed to save 1 MATVEC
  !            (see comments at the end of the subroutine)

  integer, intent(out)   :: info
  !       Gives reason for returning:
  !     abs(info)=  0 - solution found satisfying given tolerance.
  !                 1 - iteration aborted due to division by very small value.
  !                 2 - no convergence within maximum number of iterations.
  !                 3 - initial guess satisfies the stopping criterion.
  !    sign(info)=  + - residual decreased
  !                 - - residual did not reduce


  ! Local parameters

  integer, parameter :: qz_=1,zz_=3,y0_=5,yl_=6,qy_=7

  ! Local variables (2 vectors are needed in addition to tmp1_BLK=r and tmp2_BLK=u:

  real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nBLK):: &
       bicg_r1, bicg_u1

  real :: rwork(2,7)

  logical GoOn, rcmp, xpdt
  integer nmv
  real alpha, beta, omega, rho0, rho1, sigma
  real varrho, hatgamma
  real assumedzero, rnrm0, rnrm, rnrmMax0, rnrmMax
  real mxnrmx, mxnrmr, kappa0, kappal

  logical :: oktest, oktest_me

  ! External function:

  real :: maxval_abs_BLK

  !---------------------------------------------------------------------------
  call set_oktest('proj_bicgstab',oktest, oktest_me)

  if(oktest_me)write(*,*)'Proj_BiCGSTAB tol,iter:',tol,iter

  info = 0

  if (tol<=0.0) call stop_mpi('Error in Proj_BiCGSTAB: tolerance < 0')
  if (iter<=1)  call stop_mpi('Error in Proj_BiCGSTAB: maxmatvec < 2')

  !
  !     --- Initialize first residual
  !
  assumedzero = 1.e-16
  call eq_BLK(tmp1_BLK,rhs)
  call set_BLK(qx,0.0)

  nmv = 0
  !
  !     --- Initialize iteration loop
  !

  rnrm0 = sqrt( dot_product_BLK(tmp1_BLK,tmp1_BLK))

  rnrm = rnrm0
  if(oktest_me) print *,'initial rnrm:',rnrm

  mxnrmx = rnrm0
  mxnrmr = rnrm0
  rcmp = .false.
  xpdt = .false.

  alpha = 0.0
  omega = 1.0
  sigma = 1.0
  rho0 =  1.0
  !
  !     --- Iterate
  !
  select case(typestop)
  case('rel')
     GoOn = rnrm>tol*rnrm0 .and. nmv<iter
     assumedzero = assumedzero*rnrm0
     rnrmMax = 0
     rnrmMax0 = 0

  case('abs')
     GoOn = rnrm>tol       .and. nmv<iter
     assumedzero = assumedzero*rnrm0
     rnrmMax = 0
     rnrmMax0 = 0

  case('max')
     rnrmMax0 = maxval_abs_BLK(nProc,tmp1_BLK)
     rnrmMax  = rnrmMax0
     if(oktest_me) print *,'initial rnrmMax:',rnrmMax
     GoOn = rnrmMax>tol    .and. nmv<iter
     assumedzero = assumedzero*rnrmMax
  case default
     call stop_mpi('Error in Proj_BiCGSTAB: unknown typestop value')
  end select

  if (.not.GoOn) then
     if(oktest_me) print *,'Proj_BiCGSTAB: nothing to do. info = ',info
     iter = nmv
     info = 3
     return
  end if

  do while (GoOn)
     !
     !     =====================
     !     --- The BiCG part ---
     !     =====================
     !
     rho0 = -omega*rho0

     rho1 = dot_product_BLK(rhs,tmp1_BLK)

     if (abs(rho0)<assumedzero**2) then
        info = 1
        return
     endif
     beta = alpha*(rho1/rho0)
     rho0 = rho1
     call eq_plus_times_BLK(tmp2_BLK,tmp1_BLK,-beta,tmp2_BLK)

     call proj_matvec(tmp2_BLK,bicg_u1)
     nmv = nmv+1

     !DEBUG
     !write(*,*)'u =',tmp2_BLK(1:nI,1:nJ,1:nK,BLKtest)
     !write(*,*)'u1=',bicg_u1(1:nI,1:nJ,1:nK,BLKtest)

     sigma=dot_product_BLK(rhs,bicg_u1)

     if (abs(sigma)<assumedzero**2) then
        info = 1
        return
     endif

     alpha = rho1/sigma
     call add_times_BLK(qx,alpha,tmp2_BLK)

     call add_times_BLK(tmp1_BLK,-alpha,bicg_u1)

     call proj_matvec(tmp1_BLK,bicg_r1)
     nmv = nmv+1

     rnrm = sqrt( dot_product_BLK(tmp1_BLK,tmp1_BLK) )

     mxnrmx = max (mxnrmx, rnrm)
     mxnrmr = max (mxnrmr, rnrm)

     !DEBUG
     !write(*,*)'rho0, rho1, beta, sigma, alpha, rnrm:',&
     !     rho0, rho1, beta, sigma, alpha, rnrm

     !
     !  ==================================
     !  --- The convex polynomial part ---
     !  ================================== 
     !
     !    --- Z = R'R a 2 by 2 matrix
     ! i=1,j=0
     rwork(1,1) = dot_product_BLK(tmp1_BLK,tmp1_BLK)

     ! i=1,j=1
     rwork(2,1) = dot_product_BLK(bicg_r1,tmp1_BLK)
     rwork(1,2) = rwork(2,1) 

     ! i=2,j=1
     rwork(2,2) = dot_product_BLK(bicg_r1,bicg_r1)

     !
     !   --- tilde r0 and tilde rl (small vectors)
     !
     rwork(1:2,zz_:zz_+1)   = rwork(1:2,qz_:qz_+1)
     rwork(1,y0_) = -1.0
     rwork(2,y0_) = 0.0

     rwork(1,yl_) = 0.0
     rwork(2,yl_) = -1.0
     !
     !   --- Convex combination
     !
     rwork(1:2,qy_) = rwork(1,yl_)*rwork(1:2,qz_) + &
          rwork(2,yl_)*rwork(1:2,qz_+1)

     kappal = sqrt( sum( rwork(1:2,yl_)*rwork(1:2,qy_) ) )

     rwork(1:2,qy_) = rwork(1,y0_)*rwork(1:2,qz_) + &
          rwork(2,y0_)*rwork(1:2,qz_+1)

     kappa0 = sqrt( sum( rwork(1:2,y0_)*rwork(1:2,qy_) ) )

     varrho = sum( rwork(1:2,yl_)*rwork(1:2,qy_) )  
     varrho = varrho / (kappa0*kappal)

     hatgamma = sign(1.0,varrho)*max(abs(varrho),0.7) * (kappa0/kappal)

     rwork(1:2,y0_) = -hatgamma*rwork(1:2,yl_) + rwork(1:2,y0_)

     !
     !    --- Update
     !
     omega = rwork(2,y0_)

     call add_times_BLK(tmp2_BLK,-omega,bicg_u1)

     call add_times_BLK(qx,omega,tmp1_BLK)

     call add_times_BLK(tmp1_BLK,-omega,bicg_r1)

     rwork(1:2,qy_) = rwork(1,y0_)*rwork(1:2,qz_) + &
          rwork(2,y0_)*rwork(1:2,qz_+1)

     rnrm = sqrt( sum( rwork(1:2,y0_)*rwork(1:2,qy_) ) )

     select case(typestop)
     case('rel')
        GoOn = rnrm>tol*rnrm0 .and. nmv<iter
        if(oktest_me) print *, nmv,' matvecs, ', ' ||rn||/||r0|| =',rnrm/rnrm0

     case('abs')
        GoOn = rnrm>tol       .and. nmv<iter
        if(oktest_me) print *, nmv,' matvecs, ||rn|| =',rnrm

     case('max')
        rnrmMax = maxval_abs_BLK(nProc,tmp1_BLK)
        GoOn = rnrmMax>tol    .and. nmv<iter
        if(oktest_me) print *, nmv,' matvecs, max(rn) =',rnrmMax
     end select

  end do
  !
  !     =========================
  !     --- End of iterations ---
  !     =========================

  select case(typestop)
  case('rel')
     if (rnrm>tol*rnrm0) info = 2
     tol = rnrm/rnrm0

  case('abs')
     if (rnrm>tol) info = 2
     tol = rnrm

  case('max')
     if (rnrmMax>tol) info = 2
     tol = rnrmMax
  end select

  if((typestop/='max'.and.rnrm>rnrm0).or.(typestop=='max'.and.rnrmMax&
       >rnrmMax0)) info=-info

  iter = nmv

end subroutine proj_bicgstab
