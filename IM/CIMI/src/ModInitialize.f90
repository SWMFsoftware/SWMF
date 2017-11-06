Module ModCimiInitialize
  use ModCimiGrid ,ONLY: nm, nk, neng,npit1,np1
  use ModCimiPlanet,ONLY: nspec
  implicit none
  
  real :: xmm(nspec,0:nm+1),xk(0:nk+1),dphi,dmm(nspec,nm),dk(nk),delE(nspec,neng),&
          dmu(npit1),xjac(nspec,np1,nm)

  logical :: IsEmptyInitial=.false., IsDataInitial=.false., &
       IsRBSPData=.false., IsGmInitial=.true.

end Module ModCimiInitialize
