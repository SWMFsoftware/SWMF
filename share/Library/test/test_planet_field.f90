!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program test_planet_field

  use CON_planet_field, ONLY: test => test_planet_field
  use CON_planet, ONLY: init_planet_const, set_planet_defaults

  implicit none

  call init_planet_const
  call set_planet_defaults
  call test

end program test_planet_field

