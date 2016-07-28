! module to hold the photoelectron (well total SE) solution for a given line 
!and associated parameters
module ModPhotoElectron
  real, allocatable :: SeDens_C(:), SeFlux_C(:), SeHeat_C(:)

  !couple time to call update the SE flux
  real :: DtGetSe=120.0
  
  !minimum thermal density of electrons [/cc]
  real :: eThermalDensMin=1.0

  logical :: DoCoupleSE = .true., UseFeedbackFromSE=.true.

  !Fixed precipitation to pass to SE
  logical :: UseFixedPrecip = .false.
  real :: PrecipEnergyMin, PrecipEnergyMax, PrecipEnergyMean, PrecipEnergyFlux

  !PolarRain precipitation to pass to SE
  logical :: UsePolarRain = .false.
  real :: PolarRainEMin, PolarRainEMax, PolarRainEMean, PolarRainEFlux

  !Should SE be verbose with output?
  logical :: IsVerboseSE=.false.

end module ModPhotoElectron
