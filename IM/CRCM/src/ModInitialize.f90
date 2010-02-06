Module ModCrcmInitialize
  use ModCrcmGrid
  implicit none
  
  real :: xmm(nm),xk(nk),dphi,dmm(nm),dk(nk),delE(neng1),&
          dmu(npit1),xjac(nspec1,np1,nm),amu(nspec1)

  logical :: IsRestart =.false.
end Module ModCrcmInitialize
