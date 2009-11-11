program linear_advection_test

  use ModLinearAdvection, ONLY: test_linear_advection

  implicit none

  call test_linear_advection

end program linear_advection_test

subroutine CON_stop(String)
  implicit none
  character(len=*), intent(in):: String
  write(*,*)'ERROR: ',String
  stop
end subroutine CON_stop
