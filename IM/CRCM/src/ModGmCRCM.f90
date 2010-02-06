Module ModGmCrcm
  use ModCrcmGrid,ONLY: nLat => np, nLon => nt
  real, allocatable :: StateLine_VI(:,:),StateIntegral_IIV(:,:,:)
  integer :: iLineIndex_II(nLon,1:nLat),nPoint, nIntegral
  integer, parameter :: AveDens_=4, AveP_=5
  integer,parameter :: nVar=4
  
  real :: rrio(nLat,nLon),ttio(nLat,nLon)
  integer :: iLatMin=22 !Minimum latitude in MHD boundary
  
  logical :: UseGm =.true.,DoneGmCoupling=.false.,DoMultiFluidGMCoupling=.false.
end Module ModGmCrcm
