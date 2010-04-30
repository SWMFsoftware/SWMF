module CRASH_ModPolysterene
  use CRASH_ModAtomicMass
  implicit none
  SAVE
  
  !Number of elements in C_1 H_1
  integer,parameter :: nPolysterene = 2 
  
  !Relative atomic concentrations in C_1 H_1
  real,dimension(nPolysterene),parameter :: cPolysterene_I = &
         (/0.5, 0.5/)

  !Atomic number for elements:
  integer,dimension(nPolysterene),parameter :: nZPolysterene_I = &
         (/6, 1/)
   
  ! Averaged atomic mass
  real,parameter :: cAPolysterene = &
         (cAtomicMass_I(6) + &
         cAtomicMass_I(1)) / 2.0
end module CRASH_ModPolysterene
