!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!
! Wrapper to call ModLinearSolver methods from C/C++ code


!============================================================================
subroutine linear_solver_matvec(x_I, y_I, n)

  ! Fortran subroutine calling C matvec routine

  use ModLinearSolver, ONLY: iMatvecC
  use iso_c_binding

  implicit none

  interface
     subroutine linear_solver_matvec_c(x_I, y_I, n)  bind(C) 
       use iso_c_binding
       integer(c_int ), VALUE:: n
       real(c_double)        :: x_I(n)
       real(c_double)        :: y_I(n)
     end subroutine linear_solver_matvec_c
  end interface

  interface
     subroutine linear_solver_matvec_c2(x_I, y_I, n)  bind(C) 
       use iso_c_binding
       integer(c_int ), VALUE:: n
       real(c_double)        :: x_I(n)
       real(c_double)        :: y_I(n)
     end subroutine linear_solver_matvec_c2
  end interface

  integer(c_int), intent(in) :: n
  real(c_double), intent(in) :: x_I(n) 
  real(c_double), intent(out):: y_I(n)
  !--------------------------------------------------------------------------

  select case(iMatVecC)
  case(1)
     !call linear_solver_matvec_c(x_I, y_I, n)  ! MATVEC1
  case(2)
     !call linear_solver_matvec_c2(x_I, y_I, n) ! MATVEC2
  end select

end subroutine linear_solver_matvec

!============================================================================
subroutine linear_solver_gmres(iMatvec, &
     Rhs_I, x_I, lInit, n, nKrylov, &
     Tolerance, nIter, iError, lTest) bind(C)

  ! GMRES subroutine that can be called from C

  use iso_c_binding
  use ModLinearSolver, ONLY: gmres, iMatvecC

  implicit none     

  integer, intent(in):: iMatvec   ! index of C Matvec function.
  integer, intent(in):: n         ! number of unknowns.
  integer, intent(in):: nKrylov   ! size of krylov subspace
  real,    intent(in):: Rhs_I(n)  ! right hand side vector
  real, intent(inout):: x_I(n)    ! initial guess / solution vector
  integer, intent(in):: lInit     ! true if initial guess is not zero
  real, intent(inout):: Tolerance ! required / achieved residual

  integer, intent(inout):: nIter  ! maximum/actual number of iterations

  integer, intent(out)  :: iError ! gives reason for returning:
  !     abs(info)=  0 - solution found satisfying given tolerance.
  !                 2 - no convergence within maximum number of iterations.
  !                 3 - initial guess satisfies the stopping criterion.
  !    sign(info)=  + - residual decreased
  !                 - - residual did not reduce

  integer, intent(in)   :: lTest ! true if test info should be written

  ! subroutine for matrix vector multiplication 
  interface
     subroutine linear_solver_matvec(x_I, y_I, n) 
       implicit none
       ! Calculate y = M.x where M is the matrix
       integer, intent(in) :: n
       real,    intent(in) :: x_I(n)
       real,    intent(out):: y_I(n)
     end subroutine linear_solver_matvec
  end interface

  !----------------------------------------------------------------------
  iMatvecC = iMatvec

  call gmres(linear_solver_matvec, Rhs_I, x_I, lInit==1, n, nKrylov,&
       Tolerance,'rel',nIter,iError,lTest==1)

end subroutine Linear_solver_gmres
