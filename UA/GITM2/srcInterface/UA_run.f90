
subroutine UA_run(SWMFTime, SWMFTimeLimit)

  use ModGITM
  use ModTime
  use ModInputs, only: iDebugLevel, Is1D
  use ModTimeConvert, ONLY: time_real_to_int, n_day_of_year

  implicit none

  save

  real, intent(in)    :: SWMFTimeLimit
  real, intent(inout) :: SWMFTime

  integer :: index,lat,long,s,i,j,a, k, jj
  logical :: done, status_ok
  integer :: ierror, CLAWiter, n
  real    :: maxi,tt
  integer :: time_array(7)

  logical :: exist, IsDone

  CurrentTime = StartTime + SWMFTime
  EndTime     = StartTime + SWMFTimeLimit

  if (iDebugLevel > 1) then
     call time_real_to_int(CurrentTime, time_array)
     write(*,"(a,i5,5i3,i3)") "> Running UA from time : ",time_array(1:7)
  endif

  call calc_pressure

  Dt = 1.e32

  call calc_timestep_vertical
  if (.not. Is1D) call calc_timestep_horizontal

  if (iDebugLevel > 1) write(*,"(a,f13.5)") "> UA_run Dt : ",Dt

  call advance

  iStep = iStep + 1

  call write_output

  SWMFTime = CurrentTime - StartTime

end subroutine UA_run
