!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program test_block_data

  use ModBlockData, ONLY: test => test_block_data
  implicit none
  call test

end program test_block_data

subroutine CON_stop(String)
  character (len=*), intent(in) :: String
  write(*,*)'CON_stop called with String='
  write(*,*) String
  stop
end subroutine CON_stop

subroutine stop_mpi(String)
  character (len=*), intent(in) :: String
  write(*,*)'stop_mpi called with String='
  write(*,*) String
  stop
end subroutine Stop_mpi
