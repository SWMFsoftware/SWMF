!^CFG COPYRIGHT UM
module CRASH_ModPolyimide

  ! Constants for polyimide

  use CRASH_ModAtomicMass

  implicit none
  save
  
  !Number of elements in C_22 H_10 N_2 O_5
  integer,parameter :: nPolyimide = 4 
  
  !Relative atomic concentrations in C_22 H_10 N_2 O_5
  real,dimension(nPolyimide),parameter :: CPolyimide_I = &
         (/22.0, 10.0, 2.0, 5.0/)/(22.0 + 10.0 + 2.0 +5.0)

  !Atomic number for elements:
  integer,dimension(nPolyimide),parameter :: nZPolyimide_I = &
         (/6, 1, 7, 8/)
   
  ! Averaged atomic mass
  real,parameter :: cAPolyimide = &
         (cAtomicMass_I(6) * 22.0 + &
         cAtomicMass_I(1) * 10.0 + &
         cAtomicMass_I(7) *  2.0 + &
         cAtomicMass_I(8) *  5.0) / &
         (22.0 + 10.0 + 2.0 +5.0)

end module CRASH_ModPolyimide
