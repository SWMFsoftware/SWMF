!^CFG COPYRIGHT UM
Module ModNodes
  use ModSize
  implicit none
  save

  ! Logical parameter indicating static vs. dynamic allocation
  logical, parameter :: IsDynamicNodes = .false.

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

contains
  !============================================================================
  subroutine init_mod_nodes

  end subroutine init_mod_nodes 
  !============================================================================
  subroutine clean_mod_nodes

  end subroutine clean_mod_nodes

end Module ModNodes
