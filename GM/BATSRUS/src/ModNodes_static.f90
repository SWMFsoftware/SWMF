!^CFG COPYRIGHT UM
Module ModNodes

  use ModSize
  use ModIO,         ONLY: iUnitOut, write_prefix
  use ModProcMH,     ONLY: iProc

  implicit none
  save

  ! Logical parameter indicating static vs. dynamic allocation
  logical, parameter :: IsDynamicNodes = .false.

  !\
  ! Block node-centered MHD numberings
  !/
  integer :: nNodeALL
  integer, dimension(1:nI+1,1:nJ+1,1:nK+1,nBLK) :: NodeNumberLocal_NB
  integer, dimension(1:nI+1,1:nJ+1,1:nK+1,nBLK) :: NodeNumberGlobal_NB
  logical, dimension(1:nI+1,1:nJ+1,1:nK+1,nBLK) :: NodeUniqueGlobal_NB

  !\
  ! Block node-centered MHD solution and location
  !/
  real, dimension(1:nI+1,1:nJ+1,1:nK+1,nBLK) :: &
       NodeX_NB,NodeY_NB,NodeZ_NB, NodeValue_NB

contains
  !============================================================================
  subroutine init_mod_nodes

    if(IsDynamicNodes .and. iProc==0)then
       call write_prefix
       write(iUnitOut,'(a)') 'init_mod_nodes allocated arrays'
    end if

  end subroutine init_mod_nodes 
  !============================================================================
  subroutine clean_mod_nodes

    if(IsDynamicNodes .and. iProc==0)then
       call write_prefix
       write(iUnitOut,'(a)') 'clean_mod_nodes deallocated arrays'
    end if

  end subroutine clean_mod_nodes

end Module ModNodes
