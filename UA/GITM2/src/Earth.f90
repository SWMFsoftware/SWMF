
subroutine fill_photo(photoion, photoabs)

  use ModPlanet
  use ModEUV

  implicit none

  real, intent(out) :: photoion(Num_WaveLengths_High, nIons-1)
  real, intent(out) :: photoabs(Num_WaveLengths_High, nSpecies)

  integer :: iSpecies

  PhotoAbs = 0.0
  PhotoIon = 0.0

  photoabs(:,iO_)     = PhotoAbs_O
  photoabs(:,iO2_)    = PhotoAbs_O2
  if (nSpecies > 2) then
     iSpecies = iN2_
     photoabs(:,iSpecies)    = PhotoAbs_N2
  endif
  ! This may need to be as defined below....
  photoion(:,iN2P_)   = PhotoIon_N2
  photoion(:,iO2P_)   = PhotoIon_O2
  photoion(:,iNP_)    = PhotoIon_N
  photoion(:,iO_4SP_) = PhotoIon_OPlus4S
  photoion(:,iO_2DP_) = PhotoIon_OPlus2D
  photoion(:,iO_2PP_) = PhotoIon_OPlus2P

  !     Ion_Rate_Eff_N2(:,:,:) = Ion_Rate_Eff_N2(:,:,:) +              &
  !          Intensity(:,:,:,N) * (PhotoAbs_N2(N) - PhotoIon_N2(N)) *  &
  !          NDensity(:,:,:,N_N2,index)

end subroutine fill_photo

subroutine calc_planet_sources(iBlock)

  use ModInputs
  use ModSources
  use ModGITM

  implicit none

  integer, intent(in) :: iBlock

  integer :: iAlt, iError, iDir, iLat, iLon

  real :: tmp2(nLons, nLats, nAlts)
  real :: tmp3(nLons, nLats, nAlts)
  real :: Omega(nLons, nLats, nAlts)

  !\
  ! Cooling ----------------------------------------------------------
  !/

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 4) write(*,*) "=====> NO cooling", iproc, UseNOCooling

  if (UseNOCooling) then

     !  [NO] cooling 
     ! [Reference: Kockarts,G., G.R.L.,VOL.7, PP.137-140,Feberary 1980 ]
 
     Omega = 6.5e-17 * NDensityS(1:nLons,1:nLats,1:nAlts,iO_,iBlock) /      &
          (6.5e-17 * NDensityS(1:nLons,1:nLats,1:nAlts,iO_,iBlock) + 13.3)

     ! We need to check this out. I don't like the first / sign....

     NOCooling = Planck_Constant * Speed_Light / &
          5.3e-6 * &
          Omega * 13.3 *  &
          exp(- Planck_Constant * Speed_Light / &
          (5.3e-6 * Boltzmanns_Constant * &
          Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*TempUnit)) &
          * NDensityS(1:nLons,1:nLats,1:nAlts,iNO_,iBlock)

     NOCooling = NOCooling / TempUnit / &
          (Rho(1:nLons,1:nLats,1:nAlts,iBlock)*cp(:,:,1:nAlts,iBlock))

  else

     NOCooling = 0.0

  endif

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 4) write(*,*) "=====> UseOCooling", iproc, UseOCooling

  if (UseOCooling) then 

     ! [O] cooling 
     ! Reference: Kockarts, G., Plant. Space Sci., Vol. 18, pp. 271-285, 1970
     ! We reduce the LTE 63-um cooling rate by a factor of 2 for 
     ! the non-LTE effects.[Roble,1987]         

     tmp2 = exp(-228./(Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*TempUnit))
     tmp3 = exp(-326./(Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*TempUnit))

     ! In erg/cm3/s
     OCooling = (1.69e-18*tmp2 + 4.59e-20*tmp3) * &
          (NDensityS(1:nLons,1:nLats,1:nAlts,iO_,iBlock)/1.0e6) / &
          (1.0 + 0.6*tmp2 + 0.2*tmp3)
     ! In w/m3/3
     OCooling = OCooling/10.0
     ! In our special units:
     OCooling = OCooling/ TempUnit / &
          (Rho(1:nLons,1:nLats,1:nAlts,iBlock)*cp(:,:,1:nAlts,iBlock))

  else

     OCooling = 0.0

  endif

end subroutine calc_planet_sources
