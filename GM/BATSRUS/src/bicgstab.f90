!^CFG FILE IMPLICIT
!^CFG FILE NOT SIMPLE
subroutine bicgstab(matvec,rhs,qx,nonzero,n,tol,typestop,iter,info,oktest)

  ! Simple BiCGstab(\ell=1) iterative method
  ! Modified by G.Toth from the \ell<=2 version written
  ! by M.A.Botchev, Jan.'98. 
  ! Parallelization for the BATS-R-US code (2000-2001) by G. Toth
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

  ! This subroutine determines the solution of A.QX=RHS, where
  ! the matrix-vector multiplication with A is performed by
  ! the subroutine 'matvec'. For symmetric matrix use the 
  ! more efficient conjugate gradient (CG) algorithm!

  use ModMpi
  implicit none

  ! Arguments

  interface
     subroutine matvec(a,b,n)
       ! Calculate b = M.a where M is the matrix
       integer, intent(in) :: n
       real, intent(in) ::  a(n)
       real, intent(out) :: b(n)
     end subroutine matvec
  end interface
  !        subroutine for matrix vector multiplication

  integer, intent(in) :: n
  !        on input:  number of unknowns

  real, intent(inout) :: rhs(n)
  !        on input:  right-hand side vector.
  !        on output: residual vector.

  real, intent(inout) :: qx(n)
  !        on input : initial guess if nonzero is .true.
  !        on output: solution vector.

  logical, intent(in):: nonzero
  !        on input : true  if qx contains initial guess
  !                   false if initial guess is zero

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

  integer, intent(inout) :: iter
  !       on input:  maximum number of iterations to be performed.
  !       on output: actual  number of iterations done.

  integer, intent(out)   :: info
  !       Gives reason for returning:
  !     abs(info)=  0 - solution found satisfying given tolerance.
  !                 1 - iteration aborted due to division by very small value.
  !                 2 - no convergence within maximum number of iterations.
  !                 3 - initial guess satisfies the stopping criterion.
  !    sign(info)=  + - residual decreased
  !                 - - residual did not reduce

  logical, intent(in)    :: oktest
  !       write debug info if true

  ! Local parameters

  integer, parameter :: qz_=1,zz_=3,y0_=5,yl_=6,qy_=7

  ! Local variables (only 4 big vectors are needed):

  !!! Automatic arrays for BiCGSTAB!!!
!!$  real, dimension(n):: bicg_r, bicg_u, bicg_r1, bicg_u1
  real, dimension(:), allocatable :: bicg_r, bicg_u, bicg_r1, bicg_u1

  !!! Allocatable array for initial guess !!!
  real, dimension(:), allocatable :: qx0

  real :: rwork(2,7)

  logical GoOn, rcmp, xpdt
  integer nmv
  real :: alpha, beta, omega, rho0, rho1, sigma
  real :: varrho, hatgamma
  real :: assumedzero, rnrm0, rnrm, rnrmMax0, rnrmMax
  real :: mxnrmx, mxnrmr, kappa0, kappal

  ! error codes for MPI
  integer :: iError, nError 

  !---------------------------------------------------------------------------

  ! Allocate arrays that were "Automatic"
  allocate(bicg_r(n), stat=iError); call alloc_check(iError,"bicgstab:bicg_r")
  allocate(bicg_u(n), stat=iError); call alloc_check(iError,"bicgstab:bicg_u")
  allocate(bicg_r1(n), stat=iError); call alloc_check(iError,"bicgstab:bicg_r1")
  allocate(bicg_u1(n), stat=iError); call alloc_check(iError,"bicgstab:bicg_u1")

  if(oktest)write(*,*)'BiCGSTAB tol,iter:',tol,iter

  info = 0

  if (tol<=0.0) call stop_bicgstab('Error in BiCGSTAB: tolerance < 0')
  if (iter<=1)  call stop_bicgstab('Error in BiCGSTAB: maxmatvec < 2')

  if(range(1.0)>100)then
     assumedzero=0.0000000000000001
  else
     assumedzero=0.00000001
  end if

  !
  !     --- Initialize first residual
  !
  ! Calculate initial residual
  if(nonzero)then
     ! Store initial guess into qx0
     allocate(qx0(n))
     qx0=qx
     call matvec(qx,bicg_r,n)
     bicg_r = rhs - bicg_r
  else
     bicg_r = rhs
  end if

  qx = 0.0
  bicg_u=0.0

  nmv = 0
  !
  !     --- Initialize iteration loop
  !

  rnrm0 = sqrt( dot_product_mpi(bicg_r,bicg_r,n))

  rnrm = rnrm0
  if(oktest) print *,'initial rnrm:',rnrm

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
     rnrmMax0 = maxval_abs_mpi(bicg_r,n)
     rnrmMax  = rnrmMax0
     if(oktest) print *,'initial rnrmMax:',rnrmMax
     GoOn = rnrmMax>tol    .and. nmv<iter
     assumedzero = assumedzero*rnrmMax
  case default
     call stop_bicgstab('Error in BiCGSTAB: unknown typestop value')
  end select

  if (.not.GoOn) then
     iter = nmv
     info = 3
     if(oktest) print *,'BiCGSTAB: nothing to do. info = ',info
     call deallocate_bicgstab
     return
  end if

  do while (GoOn)
     !
     !     =====================
     !     --- The BiCG part ---
     !     =====================
     !
     rho0 = -omega*rho0

     rho1 = dot_product_mpi(rhs,bicg_r,n)

     if (abs(rho0)<assumedzero**2) then
        info = 1
        call deallocate_bicgstab
        return
     endif
     beta = alpha*(rho1/rho0)
     rho0 = rho1
     bicg_u = bicg_r - beta*bicg_u

     !DEBUG
     !if(oktest)&
     !     write(*,*)'rho0, rho1, alpha, beta:',rho0, rho1, alpha, beta

     call matvec(bicg_u,bicg_u1,n)
     nmv = nmv+1

     sigma=dot_product_mpi(rhs,bicg_u1,n)

     if (abs(sigma)<assumedzero**2) then
        info = 1
        call deallocate_bicgstab
        return
     endif

     alpha = rho1/sigma

     qx     = qx     + alpha*bicg_u
     bicg_r = bicg_r - alpha*bicg_u1

     call matvec(bicg_r,bicg_r1,n)
     nmv = nmv+1

     rnrm = sqrt( dot_product_mpi(bicg_r,bicg_r,n) )

     mxnrmx = max (mxnrmx, rnrm)
     mxnrmr = max (mxnrmr, rnrm)

     !DEBUG
     !if(oktest)&
     !     write(*,*)'rho0, rho1, beta, sigma, alpha, rnrm:',&
     !     rho0, rho1, beta, sigma, alpha, rnrm

     !
     !  ==================================
     !  --- The convex polynomial part ---
     !  ================================== 
     !
     !    --- Z = R'R a 2 by 2 matrix
     ! i=1,j=0
     rwork(1,1) = dot_product_mpi(bicg_r,bicg_r,n)

     ! i=1,j=1
     rwork(2,1) = dot_product_mpi(bicg_r1,bicg_r,n)
     rwork(1,2) = rwork(2,1) 

     ! i=2,j=1
     rwork(2,2) = dot_product_mpi(bicg_r1,bicg_r1,n)

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

     bicg_u = bicg_u - omega*bicg_u1

     qx     = qx     + omega*bicg_r

     bicg_r = bicg_r - omega*bicg_r1

     rwork(1:2,qy_) = rwork(1,y0_)*rwork(1:2,qz_) + &
                      rwork(2,y0_)*rwork(1:2,qz_+1)

     rnrm = sqrt( sum( rwork(1:2,y0_)*rwork(1:2,qy_) ) )

     select case(typestop)
     case('rel')
        GoOn = rnrm>tol*rnrm0 .and. nmv<iter
        if(oktest) print *, nmv,' matvecs, ', ' ||rn||/||r0|| =',rnrm/rnrm0
     case('abs')
        GoOn = rnrm>tol       .and. nmv<iter
        if(oktest) print *, nmv,' matvecs, ||rn|| =',rnrm
     case('max')
        rnrmMax = maxval_abs_mpi(bicg_r,n)
        GoOn = rnrmMax>tol    .and. nmv<iter
        if(oktest) print *, nmv,' matvecs, max(rn) =',rnrmMax
     end select

  end do
  !
  !     =========================
  !     --- End of iterations ---
  !     =========================

  ! Add initial guess if it is non-zero
  if(nonzero)then
     qx = qx + qx0
     deallocate(qx0)
  end if

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

  call deallocate_bicgstab

contains

  !============================================================================
  real function maxval_abs_mpi(a,n)
    use ModProcMH
    implicit none
    integer, intent(in) :: n
    real, intent(in)    :: a(n)

    real :: local_max_abs, global_max_abs
    !--------------------------------------------------------------------------

    local_max_abs = maxval(abs(a))

    call MPI_allreduce(local_max_abs, global_max_abs, 1, MPI_REAL, MPI_MAX, &
         iComm, iError)

    maxval_abs_mpi = global_max_abs

  end function maxval_abs_mpi

  !============================================================================
  real function dot_product_mpi(a,b,n)
    use ModProcMH
    implicit none
    integer, intent(in) :: n
    real, intent(in)    :: a(n), b(n)

    real :: local_dot_product, global_dot_product
    !--------------------------------------------------------------------------

    local_dot_product = sum(a*b)

    call MPI_allreduce(local_dot_product, global_dot_product, 1, MPI_REAL, &
         MPI_SUM, iComm, iError)

    dot_product_mpi = global_dot_product

  end function dot_product_mpi

  !============================================================================

  subroutine stop_bicgstab(str)
    use ModProcMH
    implicit none
    character (len=*), intent(in) :: str
    !--------------------------------------------------------------------------

    write(*,*)'Stopping execution with message:'
    write(*,*)str
    call MPI_abort(iComm, nError, iError)
    stop

  end subroutine stop_bicgstab

  subroutine deallocate_bicgstab
    ! Deallocate arrays that were "Automatic"
    deallocate(bicg_r)
    deallocate(bicg_u)
    deallocate(bicg_r1)
    deallocate(bicg_u1)
  end subroutine deallocate_bicgstab

end subroutine bicgstab
