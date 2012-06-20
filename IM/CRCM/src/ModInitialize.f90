Module ModCrcmInitialize
  use ModCrcmGrid
  use ModCrcmPlanet,ONLY: nspec
  implicit none
  
  real :: xmm(nm),xk(nk),dphi,dmm(nm),dk(nk),delE(neng1),&
          dmu(npit1),xjac(nspec,np1,nm)

  logical :: IsEmptyInitial=.true., IsDataInitial=.false., IsGmInitial=.false.

end Module ModCrcmInitialize
