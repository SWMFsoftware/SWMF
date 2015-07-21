! The substitute for the subroutine called by SWMF shared subroutines 
! in the case of an error
subroutine CON_stop(StringError)
  implicit none
  character (len=*), intent(in) :: StringError
  !----------------------------------------------------------------------------
  write(*,'(a)')StringError
  write(*,'(a)')'!!! AMPS_ABORT !!!'
  stop
end subroutine CON_stop
