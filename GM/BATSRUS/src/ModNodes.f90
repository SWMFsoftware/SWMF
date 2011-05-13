!^CFG COPYRIGHT UM
Module ModNodes

  use ModSize
  use ModIO,         ONLY: iUnitOut, write_prefix
  use ModProcMH,     ONLY: iProc

  implicit none
  save

  ! Logical parameter indicating static vs. dynamic allocation
  logical, parameter :: IsDynamicNodes = .true.

  !\
  ! Block node-centered MHD numberings
  !/
  integer :: nNodeALL
  integer, allocatable :: NodeNumberLocal_NB(:,:,:,:)
  integer, allocatable :: NodeNumberGlobal_NB(:,:,:,:)
  logical, allocatable :: NodeUniqueGlobal_NB(:,:,:,:)

  !\
  ! Block node-centered MHD solution and location
  !/
  real, allocatable :: NodeX_NB(:,:,:,:)
  real, allocatable :: NodeY_NB(:,:,:,:)
  real, allocatable :: NodeZ_NB(:,:,:,:)
  real, allocatable :: NodeValue_NB(:,:,:,:)

contains
  !============================================================================
  subroutine init_mod_nodes


    if(allocated(NodeNumberLocal_NB)) return
    allocate(NodeNumberLocal_NB(1:nI+1,1:nJ+1,1:nK+1,nBLK))
    allocate(NodeNumberGlobal_NB(1:nI+1,1:nJ+1,1:nK+1,nBLK))
    allocate(NodeUniqueGlobal_NB(1:nI+1,1:nJ+1,1:nK+1,nBLK))
    allocate(NodeX_NB(1:nI+1,1:nJ+1,1:nK+1,nBLK))
    allocate(NodeY_NB(1:nI+1,1:nJ+1,1:nK+1,nBLK))
    allocate(NodeZ_NB(1:nI+1,1:nJ+1,1:nK+1,nBLK))
    allocate(NodeValue_NB(1:nI+1,1:nJ+1,1:nK+1,nBLK))
    if(IsDynamicNodes .and. iProc==0)then
       call write_prefix
       write(iUnitOut,'(a)') 'init_mod_nodes allocated arrays'
    end if

  end subroutine init_mod_nodes 
  !============================================================================
  subroutine clean_mod_nodes

    if(.not.allocated(NodeNumberLocal_NB)) return
    deallocate(NodeNumberLocal_NB)
    deallocate(NodeNumberGlobal_NB)
    deallocate(NodeUniqueGlobal_NB)
    deallocate(NodeX_NB)
    deallocate(NodeY_NB)
    deallocate(NodeZ_NB)
    deallocate(NodeValue_NB)

    if(IsDynamicNodes .and. iProc==0)then
       call write_prefix
       write(iUnitOut,'(a)') 'clean_mod_nodes deallocated arrays'
    end if

  end subroutine clean_mod_nodes

end Module ModNodes
