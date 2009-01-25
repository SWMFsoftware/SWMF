!^CFG COPYRIGHT UM
!BOP -------------------------------------------------------------------
!
!MODULE: ModLinearSolver - solve a system of linear equations
!
!DESCRIPTION:
!
! Contains various methods to solve linear system of equations.
! There are both serial and parallel solvers, and direct and 
! iterative solvers.
!
!INTERFACE:

module ModLinearSolver

  !USES:
  use ModMpi
  use ModBlasLapack

  implicit none
  save

  private ! except

  !PUBLIC MEMBER FUNCTIONS:
  public :: gmres           ! GMRES iterative solver
  public :: bicgstab        ! BiCGSTAB iterative solver
  public :: cg              ! CG iterative solver for symmetric positive matrix
  public :: prehepta        ! LU preconditioner for up to hepta block-diagonal 
  public :: Uhepta          ! multiply with upper block triangular matrix
  public :: Lhepta          ! multiply with lower block triangular matrix
  public :: implicit_solver ! implicit solver in 1D with 3 point stencil
  public :: test_linear_solver

  !REVISION HISTORY:
  ! 05Dec06 - Gabor Toth - initial prototype/prolog/code based on BATSRUS code
  !EOP ___________________________________________________________________

  ! Used for tests
  integer, parameter :: rho_=1, rhou_=2, p_=3
  real, parameter :: Gamma=1.6, Dx=1.0, Inv2Dx=0.5/Dx
  integer, parameter  :: UnitTmp_=9
  integer :: iStep
  real    :: Time

  !Heat conduction, linear and non-linear
  real, parameter :: cRelaxation = 10000.0
  real    :: SpecificHeat_I(  101), HeatConductionCoef_I(  2:101)
  real    :: SpecificHeat2T_I(202), HeatConductionCoef2T_I(3:202)

contains

  subroutine gmres(matvec,Rhs,Sol,IsInit,n,nKrylov,Tol,TypeStop,&
       Iter,info,DoTest,iCommIn)

    !*************************************************************
    ! This code was initially written by Youcef Saad (May 23, 1985)
    ! then revised by Henk A. van der Vorst and Mike Botchev (Oct. 1996)
    ! Rewritten into F90 and parallelized for the BATSRUS code (May 2002) 
    ! by Gabor Toth
    ! Moved into ModLinearSolver.f90 for SWMF by Gabor Toth (December 2006)
    !*************************************************************

    implicit none

    ! subroutine for matrix vector multiplication 
    interface
       subroutine matvec(a,b,n)
         ! Calculate b = M.a where M is the matrix
         integer, intent(in) :: n
         real, intent(in) ::  a(n)
         real, intent(out) :: b(n)
       end subroutine matvec
    end interface

    integer, intent(in) :: n         !  number of unknowns.
    integer, intent(in) :: nKrylov   ! size of krylov subspace
    real,    intent(in) :: Rhs(n)    ! right hand side vector
    real, intent(inout) :: Sol(n)    ! initial guess / solution vector
    logical, intent(in) :: IsInit    ! true  if Sol contains initial guess
    real, intent(inout) :: Tol       ! required / achieved residual

    !        on input : required (relative) 2-norm or maximum norm of residual
    !        on output: achieved (relative) 2-norm or maximum norm of residual
    ! eps    == tolerance for stopping criterion. process is stopped
    !           as soon as ( ||.|| is the euclidean norm):
    !           || current residual||/||initial residual|| <= eps
    !           on OUTPUT: actual achieved norm residual (if iabs.ne.0)
    !           or achieved relative residual norm reduction

    character (len=3), intent(in) :: TypeStop
    !      Determine stopping criterion (||.|| denotes the 2-norm):
    !      typestop='rel'    -- relative stopping crit.:||res|| <= Tol*||res0||
    !      typestop='abs'    -- absolute stopping crit.: ||res|| <= Tol
    !      typestop='max'    -- maximum  stopping crit.: max(abs(res)) <= Tol

    integer, intent(inout) :: Iter    ! maximum/actual number of iterations

    integer, intent(out)   :: info    ! gives reason for returning:
    !     abs(info)=  0 - solution found satisfying given tolerance.
    !                 2 - no convergence within maximum number of iterations.
    !                 3 - initial guess satisfies the stopping criterion.
    !    sign(info)=  + - residual decreased
    !                 - - residual did not reduce

    logical, intent(in)    :: DoTest  !  write debug info if true
    integer, intent(in), optional :: iCommIn   ! MPI communicator

    ! Local variables
    integer :: iComm                           ! MPI communicator
    integer :: i,i1,its,j,k,k1
    real :: coeff,Tol1,epsmac,gam,ro,ro0,t,tmp

    ! This array used to be automatic (Krylov subspace vectors)
    real, dimension(:,:), allocatable :: Krylov_II

    ! These arrays used to be automatic (Hessenberg matrix and some vectors)
    real, dimension(:,:), allocatable :: hh
    real, dimension(:),   allocatable :: c,s,rs
    !-----------------------------------------------------------------------

    if(DoTest)write(*,*)'GMRES tol,iter:',Tol,Iter

    ! Assign the MPI communicator
    iComm = MPI_COMM_SELF
    if(present(iCommIn)) iComm = iCommIn

    ! Allocate arrays that used to be automatic
    allocate(Krylov_II(n,nKrylov+2), hh(nKrylov+1,nKrylov), &
         c(nKrylov), s(nKrylov), rs(nKrylov+1))

    if(range(1.0)>100)then
       epsmac=0.0000000000000001
    else
       epsmac=0.00000001
    endif

    its = 0
    !-------------------------------------------------------------
    ! **  outer loop starts here..
    !-------------- compute initial residual vector --------------

    RESTARTLOOP: do
       !
       !           Krylov_II(1):=A*Sol
       !
       if(IsInit.or.its>0)then
          call matvec(Sol,Krylov_II,n)
          Krylov_II(:,1)=Rhs - Krylov_II(:,1)
       else
          ! Save a matvec when starting from zero initial condition
          Krylov_II(:,1)=Rhs
       endif
       !-------------------------------------------------------------
       ro = sqrt( dot_product_mpi(Krylov_II(:,1), Krylov_II(:,1), iComm ))
       if (ro == 0.0) then
          if(its == 0)then
             info=3
          else
             info = 0
          endif
          Tol = ro
          Iter = its 
          deallocate(Krylov_II, hh, c, s, rs)
          RETURN
       end if

       ! set Tol1 for stopping criterion
       if (its == 0) then
          ro0 = ro
          if(DoTest) print *,'initial rnrm:',ro0
          if (TypeStop=='abs') then
             Tol1=Tol
             if (ro <= Tol1) then ! quit if accurate enough
                info = 3
                Tol  = ro
                Iter = its
                if(DoTest) print *,'GMRES: nothing to do. info = ',info
                deallocate(Krylov_II, hh, c, s, rs)
                RETURN
             end if
          else
             Tol1=Tol*ro
          end if
       end if

       coeff = 1.0 / ro
       Krylov_II(:,1)=coeff*Krylov_II(:,1)

       ! initialize 1-st term  of rhs of hessenberg system
       rs(1) = ro
       i = 0
       KRYLOVLOOP: do
          i=i+1
          its = its + 1
          i1 = i + 1
          !
          !           Krylov_II(i1):=A*Krylov_II(i)
          !
          call matvec(Krylov_II(:,i),Krylov_II(:,i1),n) 
          !-----------------------------------------
          !  modified gram - schmidt...
          !-----------------------------------------
          do j=1, i
             t = dot_product_mpi(Krylov_II(:,j), Krylov_II(:,i1), iComm)
             hh(j,i) = t
             Krylov_II(:,i1) = Krylov_II(:,i1) - t*Krylov_II(:,j)
          end do
          t = sqrt( dot_product_mpi(Krylov_II(:,i1),Krylov_II(:,i1), iComm) )

          hh(i1,i) = t
          if (t /= 0.0)then
             t = 1.0 / t
             Krylov_II(:,i1) = t*Krylov_II(:,i1)
          endif
          !--------done with modified gram schmidt and arnoldi step

          !-------- now  update factorization of hh

          !-------- perform previous transformations  on i-th column of h
          do k=2,i
             k1 = k-1
             t = hh(k1,i)
             hh(k1,i) = c(k1)*t + s(k1)*hh(k,i)
             hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)
          end do
          gam = sqrt(hh(i,i)**2 + hh(i1,i)**2)
          if (gam == 0.0) gam = epsmac
          !-----------#  determine next plane rotation  #-------------------
          c(i) = hh(i,i)/gam
          s(i) = hh(i1,i)/gam
          rs(i1) = -s(i)*rs(i)
          rs(i) =  c(i)*rs(i)
          !---determine residual norm and test for convergence-
          hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
          ro = abs(rs(i1))
          if (DoTest) then
             select case(TypeStop)
             case('rel')
                write(*,*) its,' matvecs, ',' ||rn||/||r0|| =',ro/ro0
             case('abs')
                write(*,*) its,' matvecs, ',' ||rn|| =',ro
             end select
          end if
          if (i >= nKrylov .or. (ro <= Tol1)) exit KRYLOVLOOP
       enddo KRYLOVLOOP

       !
       ! now compute solution. first solve upper triangular system.
       !
       ! rs := hh(1:i,1:i) ^-1 * rs

       do j=i,1,-1
          if (rs(j)/=0.0) then
             rs(j)=rs(j)/hh(j,j)
             tmp = rs(j)
             do k=j-1,1,-1
                rs(k) = rs(k) - tmp*hh(k,j)
             enddo
          endif
       enddo

       ! done with back substitution..
       ! now form linear combination to get solution
       do j=1, i
          t = rs(j)
          Sol = Sol + t*Krylov_II(:,j)
       end do

       ! exit from outer loop if converged or too many iterations
       if (ro <= Tol1 .or. its >= Iter) exit RESTARTLOOP
    end do RESTARTLOOP

    Iter=its
    Tol=Tol/Tol1*ro ! (relative) tolerance achieved
    if(ro < Tol1)then
       info = 0
    elseif(ro < ro0)then
       info =  2
    else
       info = -2
    endif

    ! Deallocate arrays that used to be automatic
    deallocate(Krylov_II, hh, c, s, rs)

  end subroutine gmres
  !=========================================================================
  subroutine bicgstab(matvec,rhs,qx,nonzero,n,tol,typestop,iter,info,&
       DoTest,iCommIn)

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

    integer, intent(in) :: n       ! number of unknowns
    real, intent(inout) :: rhs(n)  ! right-hand side / residual vector.
    real, intent(inout) :: qx(n)   ! initial guess /  solution vector.

    logical, intent(in):: nonzero  ! true  if qx contains initial guess
    real, intent(inout) :: tol     ! required / achieved norm of residual
    character (len=3), intent(in) :: typestop
    !  Determine stopping criterion (||.|| denotes the 2-norm):
    !  typestop='rel'    -- relative stopping crit.: ||res|| <= tol*||res0||
    !  typestop='abs'    -- absolute stopping crit.: ||res|| <= tol
    !  typestop='max'    -- maximum  stopping crit.: max(abs(res)) <= tol

    ! NOTE for typestop='rel' and 'abs': 
    !   To save computational work, the value of 
    !   residual norm used to check the convergence inside the main 
    !   iterative loop is computed from 
    !   projections, i.e. it can be smaller than the true residual norm
    !   (it may happen when e.g. the 'matrix-free' approach is used).
    !   Thus, it is possible that the true residual does NOT satisfy
    !   the stopping criterion ('rel' or 'abs').
    !   The true residual norm (or residual reduction) is reported on 
    !   output in parameter TOL -- this can be changed to save 1 MATVEC
    !            (see comments at the end of the subroutine)

    integer, intent(inout) :: iter
    !       on input:  maximum number of iterations to be performed.
    !       on output: actual  number of iterations done.

    integer, intent(out)   :: info
    !       Gives reason for returning:
    !  abs(info)=  0 - solution found satisfying given tolerance.
    !              1 - iteration aborted due to division by very small value.
    !              2 - no convergence within maximum number of iterations.
    !              3 - initial guess satisfies the stopping criterion.
    !  sign(info)=  + - residual decreased
    !               - - residual did not reduce

    logical, intent(in)    :: DoTest        ! write debug info if true
    integer, intent(in), optional:: iCommIn ! MPI communicator

    ! Local parameters

    integer, parameter :: qz_=1,zz_=3,y0_=5,yl_=6,qy_=7

    ! Local variables (only 4 big vectors are needed):
    integer :: iComm             ! actual MPI communicator
    ! used to be automatic arrays
    real, dimension(:), allocatable :: bicg_r, bicg_u, bicg_r1, bicg_u1

    ! allocatable array for initial guess
    real, dimension(:), allocatable :: qx0

    real :: rwork(2,7)

    logical GoOn, rcmp, xpdt
    integer nmv
    real :: alpha, beta, omega, rho0, rho1, sigma
    real :: varrho, hatgamma
    real :: assumedzero, rnrm0, rnrm, rnrmMax0, rnrmMax
    real :: mxnrmx, mxnrmr, kappa0, kappal

    !--------------------------------------------------------------------------

    ! Assign the MPI communicator
    iComm = MPI_COMM_SELF
    if(present(iCommIn)) iComm = iCommIn

    ! Allocate arrays that used to be automatic
    allocate(bicg_r(n), bicg_u(n), bicg_r1(n), bicg_u1(n)); 

    if(DoTest)write(*,*)'BiCGSTAB tol,iter:',tol,iter

    info = 0

    if (tol<=0.0) call CON_stop('Error in BiCGSTAB: tolerance < 0')
    if (iter<=1)  call CON_stop('Error in BiCGSTAB: maxmatvec < 2')

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

    rnrm0 = sqrt( dot_product_mpi(bicg_r,bicg_r,iComm))

    rnrm = rnrm0
    if(DoTest) print *,'initial rnrm:',rnrm

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
       rnrmMax0 = maxval_abs_mpi(bicg_r,iComm)
       rnrmMax  = rnrmMax0
       if(DoTest) print *,'initial rnrmMax:',rnrmMax
       GoOn = rnrmMax>tol    .and. nmv<iter
       assumedzero = assumedzero*rnrmMax
    case default
       call CON_stop('Error in BiCGSTAB: unknown typestop value')
    end select

    if (.not.GoOn) then
       iter = nmv
       info = 3
       if(DoTest) print *,'BiCGSTAB: nothing to do. info = ',info
       deallocate(bicg_r, bicg_u, bicg_r1, bicg_u1)
       RETURN
    end if

    do while (GoOn)
       !
       !     =====================
       !     --- The BiCG part ---
       !     =====================
       !
       rho0 = -omega*rho0

       rho1 = dot_product_mpi(rhs,bicg_r,iComm)

       if (abs(rho0)<assumedzero**2) then
          info = 1
          deallocate(bicg_r, bicg_u, bicg_r1, bicg_u1)
          RETURN
       endif
       beta = alpha*(rho1/rho0)
       rho0 = rho1
       bicg_u = bicg_r - beta*bicg_u

       !DEBUG
       !if(DoTest)&
       !     write(*,*)'rho0, rho1, alpha, beta:',rho0, rho1, alpha, beta

       call matvec(bicg_u,bicg_u1,n)
       nmv = nmv+1

       sigma=dot_product_mpi(rhs,bicg_u1,iComm)

       if (abs(sigma)<assumedzero**2) then
          info = 1
          deallocate(bicg_r, bicg_u, bicg_r1, bicg_u1)
          RETURN
       endif

       alpha = rho1/sigma

       qx     = qx     + alpha*bicg_u
       bicg_r = bicg_r - alpha*bicg_u1

       call matvec(bicg_r,bicg_r1,n)
       nmv = nmv+1

       rnrm = sqrt( dot_product_mpi(bicg_r,bicg_r,iComm) )

       mxnrmx = max (mxnrmx, rnrm)
       mxnrmr = max (mxnrmr, rnrm)

       !DEBUG
       !if(DoTest)&
       !     write(*,*)'rho0, rho1, beta, sigma, alpha, rnrm:',&
       !     rho0, rho1, beta, sigma, alpha, rnrm

       !
       !  ==================================
       !  --- The convex polynomial part ---
       !  ================================== 
       !
       !    --- Z = R'R a 2 by 2 matrix
       ! i=1,j=0
       rwork(1,1) = dot_product_mpi(bicg_r,bicg_r,iComm)

       ! i=1,j=1
       rwork(2,1) = dot_product_mpi(bicg_r1,bicg_r,iComm)
       rwork(1,2) = rwork(2,1) 

       ! i=2,j=1
       rwork(2,2) = dot_product_mpi(bicg_r1,bicg_r1,iComm)

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
          if(DoTest) print *, nmv,' matvecs, ', ' ||rn||/||r0|| =',rnrm/rnrm0
       case('abs')
          GoOn = rnrm>tol       .and. nmv<iter
          if(DoTest) print *, nmv,' matvecs, ||rn|| =',rnrm
       case('max')
          rnrmMax = maxval_abs_mpi(bicg_r,iComm)
          GoOn = rnrmMax>tol    .and. nmv<iter
          if(DoTest) print *, nmv,' matvecs, max(rn) =',rnrmMax
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
    deallocate(bicg_r, bicg_u, bicg_r1, bicg_u1)

  end subroutine bicgstab

  !============================================================================

  subroutine cg(matvec, Rhs_I, Sol_I, IsInit, n, Tol, TypeStop, &
       nIter, iError, DoTest, iCommIn)

    !Conjugated gradients, works if and only if the matrix A
    !is definite positive

    implicit none

    ! subroutine for matrix vector multiplication 
    interface
       subroutine matvec(a,b,n)
         ! Calculate b = M.a where M is the definite positive matrix
         integer, intent(in) :: n
         real, intent(in) ::  a(n)
         real, intent(out) :: b(n)
       end subroutine matvec
    end interface

    integer, intent(in) :: n         !  number of unknowns.
    real, intent(inout) :: Rhs_I(n)    ! right hand side vector
    real, intent(inout) :: Sol_I(n)  ! initial guess / solution vector
    logical, intent(in) :: IsInit    ! true  if Sol contains initial guess
    real,    intent(inout) :: Tol       ! required / achieved residual
    integer, intent(out):: iError    ! info about convergence
    logical, intent(in) :: DoTest

    character (len=3), intent(in) :: TypeStop
    !      Determine stopping criterion (||.|| denotes the 2-norm):
    !      typestop='rel'    -- relative stopping crit.:||res|| <= Tol*||res0||
    !      typestop='abs'    -- absolute stopping crit.: ||res|| <= Tol
 

    integer, intent(inout) :: nIter  ! maximum/actual number of iterations

    integer, intent(in), optional :: iCommIn   ! MPI communicator

    ! Local variables
    integer :: iComm                           ! MPI communicator
    integer :: MaxIter
    real :: cTolerance
    real :: rDotR, pDotADotP , rDotR0, rDotRMax, Alpha, Beta

    ! These arrays used to be automatic 
    real, dimension(:), allocatable :: Vec_I, aDotVec_I
    !-----------------------------------------------------------------------

    ! Assign the MPI communicator
    iComm = MPI_COMM_SELF
    if(present(iCommIn)) iComm = iCommIn

    ! Allocate the vectors needed for CG
    allocate(Vec_I(n), aDotVec_I(n))

    MaxIter = nIter
   
    ! compute initial residual vector 

    if(IsInit)then
       call matvec(Sol_I, ADotVec_I, n)
       Rhs_I = Rhs_I - ADotVec_I
    else
       Sol_I = 0.0
    end if

    Vec_I = 0.0

    rDotR0 = dot_product_mpi(Rhs_I, Rhs_I, iComm)

    if(TypeStop=='abs')then
       rDotRMax = Tol**2
    else
       rDotRMax = Tol**2 * rDotR0
    end if

    if(DoTest)write(*,*)'CG rDotR0, rDotRMax=', rDotR0, rDotRMax

    rDotR = rDotR0

    ! Solve the problem iteratively
    do nIter = 1, MaxIter

       if( rDotR <= rDotRMax ) EXIT

       Alpha = 1.0/rDotR

       Vec_I = Vec_I + Alpha * Rhs_I
      
       call matvec(Vec_I, aDotVec_I, n)
      
       Beta = 1.0/dot_product_mpi(Vec_I, aDotVec_I, iComm)

       Rhs_I = Rhs_I - Beta * aDotVec_I
       Sol_I = Sol_I + Beta * Vec_I

       rDotR = dot_product_mpi(Rhs_I, Rhs_I, iComm)
       if(DoTest)write(*,*)'CG nIter, rDotR=',nIter, rDotR

    end do

    ! Deallocate the temporary vectors
    deallocate(Vec_I, aDotVec_I)

    ! Calculate the achieved tolerance
    if(TypeStop=='abs')then
       Tol = sqrt(rDotR)
    else
       Tol = sqrt(rDotR/rDotR0)
    end if

    ! Set the error flag
    if(rDotR0 < rDotRMax)then
       nIter  = 0
       iError = 3
    elseif(nIter <= MaxIter)then
       iError = 0 
    elseif(rDotR < rDotR0)then
       iError = 2
    else
       iError = -2
    end if

  end subroutine cg
  !============================================================================
  real function dot_product_mpi(a_I, b_I, iComm)

    real, intent(in)    :: a_I(:), b_I(:)
    integer, intent(in) :: iComm

    real :: DotProduct, DotProductMpi
    integer :: iError
    !--------------------------------------------------------------------------

    DotProduct = dot_product(a_I, b_I)

    call MPI_allreduce(DotProduct, DotProductMpi, 1, MPI_REAL, MPI_SUM, &
         iComm, iError)

    dot_product_mpi = DotProductMpi

  end function dot_product_mpi

  !============================================================================
  real function maxval_abs_mpi(a_I, iComm)

    real, intent(in)    :: a_I(:)
    integer, intent(in) :: iComm

    real :: MaxvalAbs, MaxvalAbsMpi
    integer :: iError
    !--------------------------------------------------------------------------

    MaxvalAbs = maxval(abs(a_I))

    call MPI_allreduce(MaxvalAbs, MaxvalAbsMpi, 1, MPI_REAL, MPI_MAX, &
         iComm, iError)

    maxval_abs_mpi = MaxvalAbsMpi

  end function maxval_abs_mpi

  !============================================================================
  ! GENERAL PRECONDITIONER FOR BLOCK HEPTADIAGONAL AND PENTADIAGONAL MATRICES.
  ! DIRECT SOLVER FOR BLOCK TRIDIAGONAL MATRIX.                               
  !============================================================================
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
  !============================================================================

  subroutine prehepta(nblock,N,M1,M2,alf_in,d,e,f,e1,f1,e2,f2)

    ! This routine constructs an incomplete block LU-decomposition 
    ! of a hepta- or penta-diagonal matrix in such a way that L+U has the 
    ! same blockstructure as A. For block tri-diagonal matrix the subroutine
    ! provides a full LU decompostion. 
    !
    ! For penta-diagonal matrix, set M2=nblock and e2,f2 can be omitted.
    ! For tri-diagonal matrix set M1=M2=nblock and e1,f1,e2,f2 can be omitted.
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

    INTEGER, INTENT(IN)                        :: N, M1, M2, nblock
    REAL, INTENT(IN)                           :: alf_in
    REAL, INTENT(INOUT), DIMENSION(N,N,nblock) :: d,e,f
    REAL, INTENT(INOUT), DIMENSION(N,N,nblock), optional :: e1,f1,e2,f2

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
    ! these used to be automatic arrays
    real, dimension(:,:), allocatable :: dd
    integer, dimension(:), allocatable :: pivot

    REAL    :: alf
    INTEGER :: i,j,INFO
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
    !--------------------------------------------------------------------------

    !call timing_start('precond')

    ! Allocate arrays that used to be automatic
    allocate(dd(N,N), pivot(N))

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
       ! D = D - E.F(j-1) - E1.F1(j-M1) - E2.F2(j-M2)
       IF (j>1 )&
            CALL BLAS_GEMM('n','n',N,N,N,-one, e(:,:,j),N, f(:,:,j- 1),N,one,dd,N)
       IF (j>M1)&
            CALL BLAS_GEMM('n','n',N,N,N,-one,e1(:,:,j),N,f1(:,:,j-M1),N,one,dd,N)
       IF (j>M2)&
            CALL BLAS_GEMM('n','n',N,N,N,-one,e2(:,:,j),N,f2(:,:,j-M2),N,one,dd,N)

       IF (alf<zero) THEN

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
       END IF

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
    ENDDO

    ! Deallocate arrays that used to be automatic
    deallocate(dd, pivot)

    !call timing_stop('precond')

  END SUBROUTINE prehepta
  !============================================================================
  SUBROUTINE Uhepta(inverse,nblock,N,M1,M2,x,f,f1,f2)

    ! G. Toth, 2001

    ! This routine multiplies x with the upper triagonal U or U^{-1}
    ! which must have been constructed in subroutine prehepta.
    !
    ! For penta-diagonal matrix, set M2=nblock and f2 can be omitted.
    ! For tri-diagonal matrix set M1=M2=nblock and f1,f2 can be omitted.

    LOGICAL, INTENT(IN) :: inverse
    INTEGER, INTENT(IN) :: N,M1,M2,nblock
    REAL, INTENT(INOUT) :: x(N,nblock)
    REAL, INTENT(IN)    :: f(N,N,nblock)
    REAL, INTENT(IN), OPTIONAL :: f1(N,N,nblock), f2(N,N,nblock)

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
    !call timing_start('Uhepta')

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

    !call timing_stop('Uhepta')

  end subroutine Uhepta
  !============================================================================
  subroutine Lhepta(nblock,N,M1,M2,x,d,e,e1,e2)

    ! G. Toth, 2001
    !
    ! This routine multiplies x with the lower triangular matrix L^{-1},
    ! which must have been constructed in subroutine prehepta.
    !
    ! For penta-diagonal matrix, set M2=nblock and e2 can be omitted.
    ! For tri-diagonal matrix set M1=M2=nblock and e1,e2 can be omitted.

    INTEGER, INTENT(IN) :: N,M1,M2,nblock
    REAL, INTENT(INOUT) :: x(N,nblock)
    REAL, INTENT(IN), DIMENSION(N,N,nblock) :: d,e
    REAL, INTENT(IN), DIMENSION(N,N,nblock), OPTIONAL :: e1,e2
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
    ! this used to be an automatic array
    real, dimension(:), allocatable :: work

    INTEGER :: j
    REAL, PARAMETER :: zero=0.0, one=1.0

    ! External subroutine
    !
    ! DGEMV,   BLAS level two Matrix-Vector Product.
    !          See 'man DGEMV' for description.
    !
    !--------------------------------------------------------------------------


    ! x' = L^{-1}.x = D^{-1}.(x - E2.x'(j-M2) - E1.x'(j-M1) - E.x'(j-1))

    !call timing_start('Lhepta')

    ! Allocate arrays that used to be Automatic
    allocate(work(N))

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

    ! Deallocate arrays that used to be automatic
    deallocate(work)

    !call timing_stop('Lhepta')

  end subroutine Lhepta

  !===========================================================================

  subroutine implicit_solver(ImplPar, DtImpl, DtExpl, nCell, nVar, State_GV, &
       calc_residual, update_boundary)

    real,    intent(in) :: ImplPar, DtImpl, DtExpl
    integer, intent(in) :: nCell, nVar
    real, intent(inout) :: State_GV(-1:nCell+2, nVar)

    interface
       subroutine calc_residual(nOrder, Dt, nCell, nVar, State_GV, Resid_CV)
         implicit none
         integer, intent(in) :: nOrder, nCell, nVar
         real,    intent(in) :: Dt
         real,    intent(in) :: State_GV(-1:nCell+2, nVar)
         real,    intent(out):: Resid_CV(nCell, nVar)
       end subroutine calc_residual

       subroutine update_boundary(nCell, nVar, State_GV)
         implicit none
         integer, intent(in)    :: nCell, nVar
         real,    intent(inout) :: State_GV(-1:nCell+2, nVar)
       end subroutine update_boundary
    end interface

    real, parameter :: Eps = 1.e-6
    real, allocatable, dimension(:,:) :: &
         RightHand_CV, ResidOrig_CV, StateEps_GV, ResidEps_CV

    integer, parameter :: nDiag = 3
    real, allocatable :: Matrix_VVCI(:, :, :, :)

    real, allocatable :: Norm_V(:)  ! second norm of variables
    real, allocatable :: x_I(:)     ! linear vector of right hand side/unknowns
    real :: Coeff

    integer :: iVar, jVar, i, iStencil, iDiag, iX

    logical, parameter :: DoTestMe = .false.

    !-----------------------------------------------------------------------

    allocate(StateEps_GV(-1:nCell+2, nVar), &
         RightHand_CV(nCell, nVar), ResidOrig_CV(nCell, nVar),  &
         ResidEps_CV(nCell, nVar), Norm_V(nVar), x_I(nCell*nVar), &
         Matrix_VVCI(nVar, nVar, nCell, nDiag))

    ! Make sure that ghost cells are up-to-date
    call update_boundary(nCell, nVar, State_GV)

    ! calculate right hand side
    call calc_residual(2, DtImpl, nCell, nVar, State_GV, RightHand_CV)

    ! Calculate the unperturbed residual with first order scheme
    call calc_residual(1, DtExpl, nCell, nVar, State_GV, ResidOrig_CV)

    ! Calculate the norm for the variables
    do iVar = 1, nVar
       Norm_V(iVar) = sqrt(sum(State_GV(1:nCell,iVar)**2)/nCell)
    end do

    if(DoTestMe)write(*,*)'Norm_V=',Norm_V,' Eps=',Eps

    ! Calculate the dR/dU matrix
    do jVar = 1, nVar
       do iStencil = 1, 3
          ! Get perturbed state
          StateEps_GV = State_GV
          do i = iStencil, nCell, 3
             StateEps_GV(i, jVar) = StateEps_GV(i, jVar) + Eps*Norm_V(jVar)
          end do
          call update_boundary(nCell, nVar, StateEps_GV)
          if(DoTestMe) call save_plot(iStep, Time, nCell, nVar, StateEps_GV)
          call calc_residual(1, DtExpl, nCell, nVar, StateEps_GV, ResidEps_CV)

          ! Jacobian is multiplied with -ImplPar*DtImpl
          Coeff= -Implpar*DtImpl/(Eps*Norm_V(jVar)*DtExpl)
          do i = 1, nCell
             iDiag = modulo(i-iStencil,3)+1
             do iVar = 1, nVar
                Matrix_VVCI(iVar,jVar,i,iDiag) = &
                     Coeff*(ResidEps_CV(i,iVar) - ResidOrig_CV(i,iVar))
             end do
          end do

       end do
    end do

    if(DoTestMe)then
       write(*,'(a,/,3(3f5.1,/))') 'Matrix(:,:,2,1)=',Matrix_VVCI(:,:,2,1)
       write(*,'(a,/,3(3f5.1,/))') 'Matrix(:,:,2,2)=',Matrix_VVCI(:,:,2,2)
       write(*,'(a,/,3(3f5.1,/))') 'Matrix(:,:,2,3)=',Matrix_VVCI(:,:,2,3)
    end if

    ! Add the diagonal part J = I - delta t*dR/dU
    do i=1, nCell
       do iVar = 1, nVar
          Matrix_VVCI(iVar,iVar,i,1) = Matrix_VVCI(iVar,iVar,i,1) + 1.0
       end do
    end do

    ! L-U decomposition
    call prehepta(nCell,nVar,nCell,nCell,0.0,&
         Matrix_VVCI(:,:,:,1), Matrix_VVCI(:,:,:,2), Matrix_VVCI(:,:,:,3))

    ! Put right hand side into a linear vector
    iX = 0
    do i = 1, nCell
       do iVar = 1, nVar
          iX=iX + 1
          x_I(iX) = RightHand_CV(i, iVar)
       end do
    end do
    ! x --> L^{-1}.rhs
    call Lhepta(nCell, nVar, nCell, nCell, x_I,&
         Matrix_VVCI(:,:,:,1), Matrix_VVCI(:,:,:,2))

    ! x --> U^{-1}.L^{-1}.rhs = A^{-1}.rhs
    call Uhepta(.true.,nCell,nVar,nCell,nCell, x_I, Matrix_VVCI(:,:,:,3))

    ! Update the solution (x = U^n+1 - U^n)
    iX = 0
    do i = 1, nCell
       do iVar = 1, nVar
          iX=iX + 1
          State_GV(i, iVar) = State_GV(i, iVar) + x_I(iX)
       end do
    end do

    call update_boundary(nCell, nVar, State_GV)

    deallocate(RightHand_CV, ResidOrig_CV, StateEps_GV, ResidEps_CV, &
         Norm_V, x_I, Matrix_VVCI)

  end subroutine implicit_solver

  !=======================================================================
  subroutine test_linear_solver

    integer, parameter :: nStep=20
    integer, parameter :: nCell=51, nVar=3
    real, parameter :: DtExpl = 0.1, DtImpl = 1.0, ImplPar = 1.0

    real    :: State_GV(-1:nCell+2, nVar), StateOld_GV(-1:nCell+2, nVar)
    real    :: Resid_CV(nCell, nVar)
    !---------------------------------------------------------------------
    ! initial condition
    State_GV(:,rho_) = 1.0; State_GV(5:10,rho_)=2
    State_GV(:,rhou_)= 2.0*State_GV(:,rho_);
    State_GV(:,p_)   = 3.0; State_GV(5:10,p_) = 6.0

    open(UNITTMP_, file='linear_solver.out', status='replace')
    Time = 0.0
    call save_plot(0, Time, nCell, nVar, State_GV)
    do iStep = 1, nStep

       StateOld_GV = State_GV

       Time  = Time + DtImpl

       if(.false.)then
          ! explicit
          call calc_resid_test(2, DtExpl, nCell, nVar, State_GV, Resid_CV)
          State_GV(1:nCell,:)=State_GV(1:nCell,:)+Resid_CV
          call update_bound_test(nCell, nVar, State_GV)
          call save_plot(iStep, Time, nCell, nVar, State_GV)
       else
          ! implicit
          call implicit_solver(ImplPar, DtImpl, DtExpl, nCell, nVar, &
               State_GV, calc_resid_test, update_bound_test)
          
          ! check
          call calc_resid_test(2, DtImpl, nCell, nVar, State_GV, Resid_CV)
          write(*,*)'iStep, Max error=',iStep, &
               maxval(abs(StateOld_GV(1:nCell,:)+Resid_CV-State_GV(1:nCell,:)))
          call save_plot(iStep, Time, nCell, nVar, State_GV)
       end if
    end do
    close(UNITTMP_)

  end subroutine test_linear_solver
  !=======================================================================

  subroutine save_plot(iStep, Time, nCell, nVar, State_GV)

    integer, intent(in) :: iStep
    real,    intent(in) :: Time
    integer, intent(in) :: nCell, nVar
    real,    intent(in) :: State_GV(-1:nCell+2, nVar)

    integer :: i
    !---------------------------------------------------------------------
    write(UNITTMP_,'(a79)') 'linear solver test_hd11'
    write(UNITTMP_,'(i7,1pe13.5,3i3)') iStep,Time,1,1,nVar
    write(UNITTMP_,'(3i4)')            nCell
    write(UNITTMP_,'(100es13.5)')      Gamma
    write(UNITTMP_,'(a79)') 'x rho rhou p gamma'
    do i=1,nCell
       write(UNITTMP_,'(100es18.10))') float(i), State_GV(i, rho_:p_)
    end do

  end subroutine save_plot

  !=======================================================================
  subroutine calc_resid_test(nOrder, Dt, nCell, nVar, State_GV, Resid_CV)

    integer, intent(in) :: nOrder, nCell, nVar
    real,    intent(in) :: Dt
    real,    intent(in) :: State_GV(-1:nCell+2, nVar)
    real,    intent(out):: Resid_CV(nCell, nVar)

    integer :: i
    real :: uLeft, uRight
    !-------------------------------------------------------------------------
    do i=1, nCell
       uRight = State_GV(i+1,rhou_)/State_GV(i+1,rho_)
       uLeft  = State_GV(i-1,rhou_)/State_GV(i-1,rho_)
       Resid_CV(i,rho_) = -Dt * Inv2Dx * &
            (State_GV(i+1,rhou_) - State_GV(i-1,rhou_))
       Resid_CV(i,rhou_) = -Dt * Inv2Dx * &
            ( uRight**2*State_GV(i+1,rho_) - uLeft**2*State_GV(i-1,rho_) &
            + State_GV(i+1,p_) - State_GV(i-1,p_) )
       ! dP/dt = -div(P*u) -(g-1)*P*Div(U)
       Resid_CV(i,p_) = -Dt * Inv2Dx * &
            ( uRight*State_GV(i+1,p_) - uLeft*State_GV(i-1,p_) &
            + (Gamma-1)*State_GV(i,p_)*(uRight - uLeft) )
    end do
  end subroutine calc_resid_test

  !=======================================================================
  subroutine update_bound_test(nCell, nVar, State_GV)

    integer, intent(in)    :: nCell, nVar
    real,    intent(inout) :: State_GV(-1:nCell+2, nVar)

    ! periodic
    !State_GV(-1,:)      = State_GV(nCell-1,:)
    !State_GV( 0,:)      = State_GV(nCell,:)
    !State_GV(nCell+1,:) = State_GV(1,:)
    !State_GV(nCell+2,:) = State_GV(2,:)

    ! floating
    State_GV(-1,:)      = State_GV(1,:)
    State_GV( 0,:)      = State_GV(1,:)
    State_GV(nCell+1,:) = State_GV(nCell,:)
    State_GV(nCell+2,:) = State_GV(nCell,:)

  end subroutine update_bound_test

end module ModLinearSolver


