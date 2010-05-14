Module ModCrcm
  use ModCrcmGrid
  use ModCrcmPlanet,ONLY: nspec
  implicit none
  
  real    :: f2(nspec,np1,nt1,nm,nk),phot(nspec,np,nt),Pressure_IC(nspec,np,nt)
  real    :: dt =1., dtmax=1. ! tyical time step of crcm
  real    :: Time
  real    :: FAC_C(np,nt)
  logical :: UseMcLimiter=.false.
  real    :: BetaLimiter = 1.5
  real    :: SDtime(np,nt,nm,nk)

end Module ModCrcm
