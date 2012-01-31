!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------

subroutine stop_dgcpm(str)

  implicit none

  character (len=*), intent(in) :: str

  write(*,*) 'Stopping Execution!!!'
  write(*,*) 'Message : ',str
  stop

!!   use ModGITM
!!   use ModInputs, only: IsFramework
!!   use ModMpi
!!   implicit none
!! 
!!   character (len=*), intent(in) :: str
!!   integer :: ierror, erno
!! 
!!   if (IsFramework) then
!!      call CON_stop("UA/GITM Error: "//str)
!!   else
!!      write(*,*)'Stopping execution! iProc=',iProc,' with msg=',str
!!      call MPI_abort(iCommGITM, erno, ierror)
!!      stop
!!   endif

end subroutine stop_dgcpm

