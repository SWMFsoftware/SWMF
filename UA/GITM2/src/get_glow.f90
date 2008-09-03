subroutine get_glow(iLon,iLat,iBlock)

  use ModPlanet
  use ModTime
  use ModGITM
  use ModEUV, only    : WAVES,WAVEL,Num_WaveLengths_High,Flux_of_EUV

  implicit none

  integer, intent(in)    :: iLon, iLat, iBlock
  real, dimension(nAlts) :: ZRHO,ZO,ZN2,ZO2,ZNO,ZNS,ZND
  real, dimension(nAlts) :: ZTN,ZTI,ZTE,ZE,Zo4sp,Zo2p,Znop,Zn2p,Znp
  real, dimension(nEmissionWavelengths,nAlts) :: emissionrates
  real, dimension(nPhotoBins,nAlts) :: peup, pedown, pespec
  integer :: iAlt,iError,iBin

  call GL_settime(iTimeArray,utime)
  call GL_setloc(Altitude_GB(iLon,iLat,1:nAlts,iBlock),nAlts,Latitude(iLat,iBlock),&
       Longitude(iLon,iBlock),nPhotoBins,iError)
  if (iError.gt.0) then
     write(*,*) "nAlts too large for glow!.. in get_glow"
     call stop_GITM
  end if

  if (isFirstGlow) then
     call GL_init
     isFirstGlow = .False.
  endif

  if (isInitialGlow) then
     call GL_interp_flux(WaveS,WaveL,Flux_of_EUV,Num_WaveLengths_High)
     isInitialGlow = .False.
  endif

  ZRHO(1:nAlts)  = Rho(iLon,iLat,1:nAlts,iBlock)*0.000001*1000.0
  ZO(1:nAlts)    = NDensityS(iLon,iLat,1:nAlts,iO_3P_,iBlock)*0.000001
  ZN2(1:nAlts)   = NDensityS(iLon,iLat,1:nAlts,iN2_,iBlock)*0.000001
  ZO2(1:nAlts)   = NDensityS(iLon,iLat,1:nAlts,iO2_,iBlock)*0.000001
  ZNO(1:nAlts)   = NDensityS(iLon,iLat,1:nAlts,iNO_,iBlock)*0.000001
  ZNS(1:nAlts)   = NDensityS(iLon,iLat,1:nAlts,iN_4S_,iBlock)*0.000001
  ZND(1:nAlts)   = NDensityS(iLon,iLat,1:nAlts,iN_2D_,iBlock)*0.000001
 
  call GL_setND(ZRHO,ZO,ZN2,ZO2,ZNO,ZNS,ZND)

  ZTN(1:nAlts)   = Temperature(iLon,iLat,1:nAlts,iBlock) * &
             TempUnit(iLon,iLat,1:nAlts)
  ZTI(1:nAlts)   = iTemperature(iLon,iLat,1:nAlts,iBlock)
  ZTE(1:nAlts)   = eTemperature(iLon,iLat,1:nAlts,iBlock)
  
  call GL_settemp(ZTN,ZTI,ZTE)

  Zo4sp = IDensityS(iLon,iLat,1:nAlts,iO_4SP_,iBlock)*0.000001
  Zo2p =  IDensityS(iLon,iLat,1:nAlts,iO2P_,iBlock)*0.000001
  Znop = IDensityS(iLon,iLat,1:nAlts,iNOP_,iBlock)*0.000001
  ZE = IDensityS(iLon,iLat,1:nAlts,ie_,iBlock)*0.000001
  Zn2p = IDensityS(iLon,iLat,1:nAlts,iN2P_,iBlock)*0.000001
  Znp =  IDensityS(iLon,iLat,1:nAlts,iNP_,iBlock)*0.000001

  call GL_setID(Zo4sp, Zo2p, Znop, Zn2p, Znp, ZE)

  call GL_getvals(iProc,emissionrates,peup,pedown,pespec,PhotoEBins)

  do iAlt = 1, nAlts
     vEmissionRate(iLon,iLat,iAlt,:,iBlock)     = emissionrates(:,iAlt) * 100000.0
     PhotoEFluxU(iLon,iLat,iAlt,:,iBlock)       = peup(:,iAlt)   * 1000000.0
     PhotoEFluxD(iLon,iLat,iAlt,:,iBlock)       = pedown(:,iAlt) * 1000000.0
     PhotoElectronRate(iLon,iLat,iAlt,:,iBlock) = pespec(:,iAlt) * 1000000.0
  enddo


     PhotoEFluxTotal = 0

     !Sum over all Energy Bins...
     do iBin = 1, nPhotoBins
        PhotoEFluxTotal(:,:,:,iBlock,1) = PhotoEFluxTotal(:,:,:,iBlock,1) + &
             PhotoEFluxU(:,:,:,ibin,iBlock)
        PhotoEFluxTotal(:,:,:,iBlock,2) = PhotoEFluxTotal(:,:,:,iBlock,2) + &
             PhotoEFluxD(:,:,:,iBin,iBlock)
     enddo


end subroutine get_glow
