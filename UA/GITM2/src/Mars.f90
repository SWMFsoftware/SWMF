
subroutine fill_photo(photoion, photoabs)

  use ModPlanet
  use ModEUV

  implicit none

  real, intent(out) :: photoion(Num_WaveLengths_High, nIons-1)
  real, intent(out) :: photoabs(Num_WaveLengths_High, nSpecies)

  photoabs           = 0.0
  photoabs(:,iCO2_)  = PhotoAbs_CO2
  photoabs(:,iCO_)   = PhotoAbs_CO
  photoabs(:,iN2_)   = PhotoAbs_N2
  photoabs(:,iO_)    = PhotoAbs_O

  photoion           = 0.0
  photoion(:,iOP_)   = PhotoIon_OPlus4S
  photoion(:,iCO2P_) = PhotoIon_CO2

end subroutine fill_photo

subroutine calc_planet_sources(iBlock)

  use ModInputs
  use ModSources
  use ModGITM

  implicit none

  integer, intent(in) :: iBlock

  ! There are no sources specificly for Mars

end subroutine calc_planet_sources
