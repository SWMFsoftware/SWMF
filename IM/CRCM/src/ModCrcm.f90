Module ModCrcm
  use ModCrcmGrid
  use ModCrcmPlanet,ONLY: nspec
  implicit none

  SAVE

  real    :: dt =1., dtmax=1. ! tyical time step of crcm
  real    :: Time = 0.0
  logical :: UseMcLimiter=.false.
  real    :: BetaLimiter = 1.5
  real, allocatable:: SDtime(:,:,:,:), f2(:,:,:,:,:)
  real, allocatable:: phot(:,:,:), Pressure_IC(:,:,:)
  real, allocatable:: FAC_C(:,:)

contains

  subroutine init_mod_crcm

    if(allocated(f2)) RETURN

    allocate(                      &
         SDtime(np,nt,nm,nk),      &
         f2(nspec,np1,nt1,nm,nk),  &
         phot(nspec,np,nt),        &
         Pressure_IC(nspec,np,nt), &
         FAC_C(np,nt)              &
    )

    ! Not clear why these need initialization !!!
    FAC_C = 0.0
    Pressure_IC = 0.0

  end subroutine init_mod_crcm

end Module ModCrcm
