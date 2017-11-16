!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModNodes

  use BATL_lib, ONLY: &
       test_start, test_stop

  use ModSize
  use ModIO,         ONLY: iUnitOut, write_prefix
  use ModProcMH,     ONLY: iProc

  implicit none
  save

  !\
  ! Block node-centered MHD numberings
  !/
  integer :: nNodeALL
  integer, allocatable :: NodeNumberLocal_NB(:,:,:,:)
  integer, allocatable :: NodeNumberGlobal_NB(:,:,:,:)
  logical, allocatable :: NodeUniqueGlobal_NB(:,:,:,:)

contains
  !============================================================================
  subroutine init_mod_nodes

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'init_mod_nodes'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    if(allocated(NodeNumberLocal_NB)) RETURN
    allocate(NodeNumberLocal_NB(1:nI+1,1:nJ+1,1:nK+1,MaxBlock))
    allocate(NodeNumberGlobal_NB(1:nI+1,1:nJ+1,1:nK+1,MaxBlock))
    allocate(NodeUniqueGlobal_NB(1:nI+1,1:nJ+1,1:nK+1,MaxBlock))
    if(iProc==0)then
       call write_prefix
       write(iUnitOut,'(a)') 'init_mod_nodes allocated arrays'
    end if

    call test_stop(NameSub, DoTest)
  end subroutine init_mod_nodes
  !============================================================================
  subroutine clean_mod_nodes

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'clean_mod_nodes'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    if(.not.allocated(NodeNumberLocal_NB)) RETURN
    deallocate(NodeNumberLocal_NB)
    deallocate(NodeNumberGlobal_NB)
    deallocate(NodeUniqueGlobal_NB)

    if(iProc==0)then
       call write_prefix
       write(iUnitOut,'(a)') 'clean_mod_nodes deallocated arrays'
    end if

    call test_stop(NameSub, DoTest)
  end subroutine clean_mod_nodes
  !============================================================================

end module ModNodes
!==============================================================================
