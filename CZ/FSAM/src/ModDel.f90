module ModDel
  use ModPar,   ONLY: in, jn, kn
  implicit none
  
  real, dimension(1:in,1:jn,1:kn) :: dels1, dels2, dels3, dels
  real, dimension(1:in,1:jn,1:kn) :: emf1s, emf2s, emf3s

end module ModDel
