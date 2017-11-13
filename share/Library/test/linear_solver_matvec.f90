!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!
! Wrapper to call ModLinearSolver methods from C/C++ code

!============================================================================
subroutine linear_solver_matvec(x_I, y_I, n)

  ! Fortran subroutine calling C matvec routine

  use  iso_c_binding
  implicit none
  interface
     subroutine linear_solver_matvec_c(x_I, y_I, n)  bind(C) 
       use iso_c_binding
       integer(c_int ), VALUE:: n
       real(c_double)        :: x_I(n)
       real(c_double)        :: y_I(n)
     end subroutine linear_solver_matvec_c
  end interface

  integer(c_int), intent(in) :: n
  real(c_double), intent(in) :: x_I(n) 
  real(c_double), intent(out):: y_I(n)
  !--------------------------------------------------------------------------

  call linear_solver_matvec_c(x_I, y_I, n)

end subroutine linear_solver_matvec
!============================================================================
subroutine CON_stop(String)
  implicit none
  character(len=*), intent(in):: String
  write(*,*)'ERROR: ',String
  stop
end subroutine CON_stop
