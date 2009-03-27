program Harmonics

  ! Transform raw magnetogram file into spherical harmonics file

  use ModMagHarmonics
  implicit none
  !----------------------------------------------------------------------------
  call read_raw_magnetogram 
  call calc_harmonics
  
end program Harmonics
