
module ModSources

  use ModSize
  use ModPlanet, only: nSpecies

  !\
  ! Sources for neutral temperature
  !/

  real, dimension(nLons, nLats, nAlts) :: &
       Conduction, NOCooling, OCooling, &
       AuroralHeating, JouleHeating

  real, allocatable :: EuvHeating(:,:,:,:)
  real, allocatable :: eEuvHeating(:,:,:,:)

  !\
  ! Stuff for auroral energy deposition and ionization
  !/

  real, dimension(:), allocatable :: &
       ED_grid, ED_Energies, ED_Flux, ED_Ion, ED_Heating
  integer :: ED_N_Energies, ED_N_Alts
  real, dimension(nAlts) :: ED_Interpolation_Weight
  integer, dimension(nAlts) :: ED_Interpolation_Index

  real, allocatable :: AuroralIonRateS(:,:,:,:,:)
  real, allocatable :: AuroralHeatingRate(:,:,:,:)
  real :: ChemicalHeatingRate(nLons, nLats, nAlts)

  real :: Diffusion(nLons, nLats, nAlts, nSpecies)
  real :: NeutralFriction(nLons, nLats, nAlts, nSpecies)

contains
  !=========================================================================
  subroutine init_mod_sources

    use ModUtilities, ONLY: check_allocate
    integer :: iError
    !----------------------------------------------------------------------
    if(allocated(EuvHeating)) return
    allocate(EuvHeating(nLons, nLats, nAlts,nBlocks),stat=iError)
    call check_allocate(iError,'EuvHeating')
    allocate(eEuvHeating(nLons, nLats, nAlts,nBlocks),stat=iError)
    call check_allocate(iError,'eEuvHeating')
    allocate(AuroralIonRateS(nLons,nLats,nAlts,nSpecies,nBlocks),stat=iError)
    call check_allocate(iError,'AuroralIonRateS')
    allocate(AuroralHeatingRate(nLons,nLats,nAlts,nBlocks),stat=iError)
    call check_allocate(iError,'AuroralHeatingRate')
  end subroutine init_mod_sources
  !=========================================================================
  subroutine clean_mod_sources

    if(.not.allocated(EuvHeating)) return
    deallocate(EuvHeating)
    deallocate(eEuvHeating)
    deallocate(AuroralIonRateS)
    deallocate(AuroralHeatingRate)
  end subroutine clean_mod_sources
  !=========================================================================
end module ModSources
