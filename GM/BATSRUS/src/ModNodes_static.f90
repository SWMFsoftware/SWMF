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
