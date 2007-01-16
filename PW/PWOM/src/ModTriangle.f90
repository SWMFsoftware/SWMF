Module ModTriangle
  !use ModPWOM, only : maxLine, nTotalLine
  implicit none
  integer, parameter :: maxLine=500  
  integer, allocatable, dimension ( :, : )   :: triangle_node
  real , dimension ( 2, maxLine ) :: table
  integer :: node_num
  integer triangle_num


end Module ModTriangle

