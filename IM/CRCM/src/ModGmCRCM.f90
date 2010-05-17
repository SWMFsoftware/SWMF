Module ModGmCrcm

  use ModCrcmGrid,  ONLY: nLat => np, nLon => nt
  use ModCrcmPlanet,ONLY: nspec

  implicit none

  real, allocatable :: StateLine_VI(:,:),StateIntegral_IIV(:,:,:)
  integer :: iLineIndex_II(nLon,1:nLat),nPoint, nIntegral
  integer, parameter :: AveDens_=4, AveP_=5,AveHpRho_=7,AveOpRho_=8,&
                        AveHpP_=9,AveOpP_=10
  integer,parameter :: nVar=4
  
  real :: Den_IC(nspec,nLat,nLon) = 0.0, Temp_IC(nspec,nLat,nLon) = 0.0
  integer :: iLatMin=22 !Minimum latitude in MHD boundary
  
  logical :: UseGm                  = .true.
  logical :: DoneGmCoupling         = .false.
  logical :: DoMultiFluidGMCoupling = .false.

end Module ModGmCrcm
