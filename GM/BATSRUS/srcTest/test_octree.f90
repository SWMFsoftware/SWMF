program test_octree

  use ModOctreeNew, ONLY: test => test_octree
  implicit none
  call test

end program test_octree

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
end subroutine stop_mpi
