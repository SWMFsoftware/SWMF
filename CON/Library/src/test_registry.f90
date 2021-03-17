!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program test_registry

  use CON_world
  use CON_comp_param
  implicit none

  logical :: UseMe
  character (len=lNameVersion) :: Name
  real :: Version
  integer :: iError, lComp, iComp, iProc
  !----------------------------------------------------------------------------
  call MPI_init(iError)

  call world_init

  call world_setup

  call world_used(IsVerbose=.true.)

  call MPI_barrier(i_comm(), iError)

  UseMe = is_proc('GM')
  ! UseMe = is_proc('IX') ! Test wrong name
  ! UseMe = is_proc(53)   ! Test wrong index

  do iProc = 0, n_proc()-1
     if(iProc == i_proc()) &
          write(*,*)'iProcWorld, is_proc(GM,IE,IH,IM,CON), iProcUsed=',&
          i_proc(), &
          is_proc(GM_), is_proc('IE'), is_proc('IH'), is_proc('IM'), &
          is_proc('CON'), i_proc(CON_)
     call MPI_barrier(i_comm(), iError)
  end do

  if(is_proc('IE'))&
     call put_comp_info('IE',NameVersion='Ridley',Version=1.01)
  if(is_proc0('IE'))then
     call get_comp_info(IE_,NameVersion=Name,Version=Version)
     write(*,*)'IE ',Name,' ',Version,' iProc=',i_proc()
  end if

  ! Switch off IE
  call put_comp_info(IE_, Use=.false.)
  call world_used(IsVerbose=.false.)

  if(is_proc0())then
     write(*,'(a)') 'Switched off IE'
     write(*,'(a,i2,a,i2)') 'iProc0Used=', i_proc0(CON_), &
          ' nProcUsed = ', n_proc(CON_)
  end if
  call MPI_barrier(i_comm(), iError)

  do iProc = 0, n_proc()-1
     if(iProc == i_proc()) &
          write(*,*)'iProcWorld, is_proc(GM,IE,IH,IM,CON), iProcUsed=',&
          i_proc(), &
          is_proc(GM_), is_proc('IE'), is_proc('IH'), is_proc('IM'), &
          is_proc('CON'), i_proc(CON_)
     call MPI_barrier(i_comm(), iError)
  end do

  if(is_proc0())write(*,'(a)') 'GM threads:'
  do iProc = 0, n_proc()-1
     if(iProc == i_proc()) write(*,'(a,i3,i3,l3,i8)') &
             'iProcWorld, i_thread(GM), is_proc(CON), i_proc(GM)=', &
             i_proc(), i_thread('GM'), is_proc(CON_), i_proc(GM_)
     call MPI_barrier(i_comm(), iError)
  end do

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

  call MPI_finalize(iError)

end program test_registry
!==============================================================================

