program test_array
  implicit none
  integer, parameter :: iM_II(5,5)=reshape([ 3, 7, 5,21, 8,&
                                            16, 8,17,53, 7,&
                                            14, 6,35,18, 1,&
                                            13,19, 4,22,11,&
                                             9,21, 1,23,12],[5,5])
  write(*,'(4i5)')iM_II(2:5,[1,2,3,5])
end program test_array
