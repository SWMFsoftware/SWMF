!^CFG COPYRIGHT UM
Module ModNodes
  use ModSize
  implicit none
  save

  !\
  ! Block node-centered MHD numberings
  !/
  integer :: nNodeALL
  integer, dimension(0:nI,0:nJ,0:nK,nBLK) :: NodeNumberLocal_IIIB
  integer, dimension(0:nI,0:nJ,0:nK,nBLK) :: NodeNumberGlobal_IIIB
  logical, dimension(0:nI,0:nJ,0:nK,nBLK) :: NodeUniqueGlobal_IIIB

  !\
  ! Block node-centered MHD solution and location
  !/
  real, dimension(0:nI,0:nJ,0:nK,nBLK) :: &
       NodeX_IIIB,NodeY_IIIB,NodeZ_IIIB, NodeValue_IIIB

end Module ModNodes
