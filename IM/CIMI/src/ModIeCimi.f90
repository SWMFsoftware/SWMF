module ModIeCimi
  use ModCimiGrid,ONLY: np, nt
  implicit none 
  
  real :: pot(np,nt)
  logical :: UseIe =.true.
  logical :: UseWeimer =.false.
contains

  subroutine init_mod_iecimi
    
    pot = 0.0

  end subroutine init_mod_iecimi

end module ModIeCimi
