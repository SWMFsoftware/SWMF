Module ModCrcm
  use ModCrcmGrid
  implicit none
  
  real    :: f2(nspec1,np1,nt1,nm,nk),phot(np,nt),Pressure_C(np,nt)
  integer :: kspec(nspec)
  real    :: dt =10., dtmax=10. ! tyical time step of crcm
  real    :: Time
  real    :: FAC_C(np,nt)

end Module ModCrcm
