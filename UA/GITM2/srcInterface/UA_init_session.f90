
subroutine UA_init_session(iSession, SWMFTime)

  implicit none

  real, intent(in)    :: SWMFTime
  integer, intent(in) :: iSession

  logical :: IsFirstTime = .true.

  if (IsFirstTime) then

     call initialize_gitm
     call write_output

  endif

end subroutine UA_init_session
