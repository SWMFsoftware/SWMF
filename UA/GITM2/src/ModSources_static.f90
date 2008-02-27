
module ModSources

  use ModSizeGitm
  use ModPlanet, only: nSpecies

  !\
  ! Sources for neutral temperature
  !/

  real, dimension(nLons, nLats, nAlts) :: &
       Conduction, NOCooling, OCooling, &
       AuroralHeating, JouleHeating, IonPrecipHeating

  real, dimension(nLons, nLats, nAlts,nBlocksMax) :: &
       EuvHeating,eEuvHeating, &
       RadCooling, RadCoolingRate, RadCoolingErgs, EuvHeatingErgs, &
       LowAtmosRadRate

  !\
  ! Stuff for auroral energy deposition and ionization
  !/

  real, dimension(:), allocatable :: &
       ED_grid, ED_Energies, ED_Flux, ED_Ion, ED_Heating
  integer :: ED_N_Energies, ED_N_Alts
  real, dimension(nAlts) :: ED_Interpolation_Weight
  integer, dimension(nAlts) :: ED_Interpolation_Index

  real :: AuroralIonRateS(nLons, nLats, nAlts, nSpecies, nBlocksMax)
  real :: AuroralHeatingRate(nLons, nLats, nAlts, nBlocksMax)
  real :: IonPrecipIonRateS(nLons, nLats, nAlts, nSpecies, nBlocksMax)
  real :: IonPrecipHeatingRate(nLons, nLats, nAlts, nBlocksMax)
  real :: ChemicalHeatingRate(nLons, nLats, nAlts)

  real :: VerticalTempSource(nLons, nLats, nAlts)
  real :: HorizontalTempSource(nLons, nLats, nAlts)

  real :: Diffusion(nLons, nLats, nAlts, nSpecies)
  real :: NeutralFriction(nLons, nLats, nAlts, nSpecies)
  real :: IonNeutralFriction(nLons, nLats, nAlts, nSpecies)

  real :: KappaEddyDiffusion(nLons, nLats, -1:nAlts+2, nBlocksMax)

contains
  !=========================================================================
  subroutine init_mod_sources
  end subroutine init_mod_sources
  !=========================================================================
  subroutine clean_mod_sources
  end subroutine clean_mod_sources
  !=========================================================================
end module ModSources
