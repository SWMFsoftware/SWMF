program linear_solver_test

  use ModLinearSolver, ONLY: test_linear_solver

  implicit none

  call test_linear_solver

end program linear_solver_test

subroutine CON_stop(String)
  implicit none
  character(len=*), intent(in):: String
  write(*,*)'ERROR: ',String
  stop
end subroutine CON_stop
