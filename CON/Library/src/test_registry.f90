!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program test_registry

  use CON_world
  use CON_comp_param
  implicit none
  logical :: UseMe
  character (len=lNameVersion) :: Name
  real :: Version
  integer :: lComp,iComp
  !---------------------------------------------------------------------------
  call world_init

  call world_setup

  UseMe = is_proc('GM')
  ! UseMe = is_proc('IX') ! Test wrong name
  ! UseMe = is_proc(53)   ! Test wrong index

  write(*,*)'is_proc(GM,IE,IH,IM), iProcWorld=',&
       is_proc(GM_),is_proc('IE'),is_proc('IH'),is_proc('IM'),i_proc()

!  if(is_proc('GM'))
  write(*,*)'GM iProc, iProcWorld=',i_proc('GM'),i_proc()

  if(is_proc('IE'))&
     call put_comp_info('IE',NameVersion='Ridley',Version=1.01)
  if(is_proc0('IE'))then
     call get_comp_info(IE_,NameVersion=Name,Version=Version)
     write(*,*)'IE ',Name,' ',Version,' iProc=',i_proc()
  end if

  if(is_proc0())then
     do lComp = 1, n_comp()
        iComp = i_comp(lComp)
        write(*,'(a,i2,a,l1,a,i4)') 'lComp=',&
             lComp,' '//NameComp_I(iComp)//' use=',use_comp(iComp),&
             ' nProc=',n_proc(iComp)
     end do
     do iComp=1,MaxComp
        write(*,'(a,l1)')'use_'//NameComp_I(iComp)//' = ',use_comp(iComp)
     end do
  end if

  call world_clean

end program test_registry

subroutine CON_stop(String)

  use ModMpi
  use CON_world, ONLY: i_proc_world
  implicit none
  character(len=*),intent(in) :: String
  integer :: iError, nError
  !--------------------------------------------------------------------------
  write(*,'(a,i3)')'!!! SWMF_ABORT !!! requested by processor ',i_proc_world()
  write(*,'(a)')String
  call MPI_abort(MPI_COMM_WORLD, nError, iError)
  stop
end subroutine CON_stop
