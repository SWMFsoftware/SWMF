
subroutine calc_pressure

  use ModGITM
  use ModPlanet
  use ModConstants
  use ModInputs, ONLY: iDebugLevel

  implicit none

  integer :: iSpecies

  call report("calc_pressure",2)

  Pressure    = Temperature * Rho

  NDensity = 0.0
  do iSpecies = 1, nSpecies
     NDensity(:,:,:,1:nBlocks) = &
          NDensity(:,:,:,1:nBlocks) + NDensityS(:,:,:,iSpecies,1:nBlocks)
  enddo

  IPressure = 0.0
  do iSpecies = 1, nIons-1
     IPressure(:,:,:,1:nBlocks) = IPressure(:,:,:,1:nBlocks) + &
          IDensityS(:,:,:,iSpecies,1:nBlocks) * &
          Boltzmanns_Constant * ITemperature(:,:,:,1:nBlocks)
  enddo

  ePressure(:,:,:,1:nBlocks) = &
       IDensityS(:,:,:,ie_,1:nBlocks) * Boltzmanns_Constant * &
       eTemperature(:,:,:,1:nBlocks)

  if (iDebugLevel > 2)write(*,*) &
       'calc_pressure iPres, ePres, iDens, iTemp, eTemp=',&
       sum(abs(IPressure(:,:,:,1:nBlocks))),     &
       sum(abs(ePressure(:,:,:,1:nBlocks))),     &
       sum(abs(IDensityS(:,:,:,ie_,1:nBlocks))), &
       sum(abs(ITemperature(:,:,:,1:nBlocks))),  &
       sum(abs(eTemperature(:,:,:,1:nBlocks)))

end subroutine calc_pressure
