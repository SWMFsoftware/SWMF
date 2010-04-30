module CRASH_ModAcrylic
  use CRASH_ModAtomicMass
  implicit none
  SAVE
  
  !Number of elements in C_5 O_2 H_8
  integer,parameter :: nAcrylic = 3 
  
  !Relative atomic concentrations in C_5 O_2 H_8
  real,dimension(nAcrylic),parameter :: cAcrylic_I = &
         (/5.0/15, 2.0/15, 8.0/15/)

  !Atomic number for elements:
  integer,dimension(nAcrylic),parameter :: nZAcrylic_I = &
         (/6, 8, 1/)
   
  ! Averaged atomic mass
  real,parameter :: cAAcrylic = &
         (cAtomicMass_I(6) * 5.0 + &
          cAtomicMass_I(8) * 2.0 + &
          cAtomicMass_I(1) * 8.0 ) / 15
end module CRASH_ModAcrylic
