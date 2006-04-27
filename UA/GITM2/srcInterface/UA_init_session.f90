
subroutine UA_init_session(iSession, SWMFTime)

  use CON_physics,    ONLY: get_time
  use ModTime, only : StartTime, iTimeArray, CurrentTime

  implicit none

  real, intent(in)    :: SWMFTime
  integer, intent(in) :: iSession

  logical :: IsFirstTime = .true.

  if (IsFirstTime) then

     ! Set time related variables for UA
     call get_time(tStartOut = StartTime)

     CurrentTime = StartTime + SWMFTime
     call time_real_to_int(StartTime, iTimeArray)

     call fix_vernal_time

     call initialize_gitm(CurrentTime)
     call write_output

     IsFirstTime = .false.

  endif

end subroutine UA_init_session
