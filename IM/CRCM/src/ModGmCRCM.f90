Module ModGmCrcm
  use ModCrcmGrid,  ONLY: nLat => np, nLon => nt
  use ModCrcmPlanet,ONLY: nspec
  real, allocatable :: StateLine_VI(:,:),StateIntegral_IIV(:,:,:)
  integer :: iLineIndex_II(nLon,1:nLat),nPoint, nIntegral
  integer, parameter :: AveDens_=4, AveP_=5,AveHpRho_=7,AveOpRho_=8,&
                        AveHpP_=9,AveOpP_=10
  integer,parameter :: nVar=4
  
  real :: Den_IC(nspec,nLat,nLon),Temp_IC(nspec,nLat,nLon)
  integer :: iLatMin=22 !Minimum latitude in MHD boundary
  
  logical :: UseGm =.true.,DoneGmCoupling=.false.,DoMultiFluidGMCoupling=.false.
end Module ModGmCrcm
