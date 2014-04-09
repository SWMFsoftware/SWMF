module ModField
  use ModPar,   ONLY: in, jn, kn
  implicit none

  real, dimension(1:in,1:jn,1:kn) :: v1, v2, v3, b1, b2, b3, s, p
  
  real, dimension(1:in,1:jn,1:kn) :: v1_prev, v2_prev, v3_prev, b1_prev, b2_prev, &
       b3_prev, s_prev

end module ModField
