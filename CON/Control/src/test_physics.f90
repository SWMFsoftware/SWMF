program test_physics

  use CON_physics
  use CON_planet
  use CON_planet_field
  use CON_axes
  use CON_time

  !implicit none

  call test_time
  stop

  call time_int_to_real(TimeEquinox)
  ! Equinox
  TimeStart = TimeEquinox

  ! Equinox noon
  !TimeStart % Time = TimeStart % Time + 4*3600 + 25*60

  ! Midsummer noon
  TimeStart % Time = TimeStart % Time - 3*3600 - 2*60 - 11 &
       +OrbitalPeriodEarth*0.25 

  !TimeStart % Time = TimeStart % Time + 6*3600    ! 6 pm.

  call time_real_to_int(TimeStart)

  call test_axes

  DoTimeAccurate = .true.
  call test_planet_field

end program test_physics
