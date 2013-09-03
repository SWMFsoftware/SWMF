Module ModCrcm
  use ModCrcmGrid
  use ModCrcmPlanet,ONLY: nspec
  implicit none

  SAVE

  real    :: dt =5., dtmax=5. ! tyical time step of crcm
  real    :: Time = 0.0
  logical :: UseMcLimiter=.false.
  real    :: BetaLimiter = 1.5
  real    :: Pmin = 1e-2
  real, allocatable:: SDtime(:,:,:,:), f2(:,:,:,:,:)
  real, allocatable:: phot(:,:,:), Ppar_IC(:,:,:), Pressure_IC(:,:,:),PressurePar_IC(:,:,:)
  real, allocatable:: FAC_C(:,:)
  real, allocatable:: Bmin_C(:,:)

  !Variable for E gain and loss
  real, allocatable :: driftin(:),driftout(:),rbsumLocal(:),rbsumGlobal(:)
  real, allocatable :: eChangeOperator_IV(:,:)
  integer, parameter :: nOperator = 5, OpBfield_=1, OpDrift_=2, OpLossCone_=3, OpChargeEx_=4, OpStrongDiff_=5
contains

  subroutine init_mod_crcm

    if(allocated(f2)) RETURN

    allocate(                      &
         SDtime(np,nt,nm,nk),      &
         f2(nspec,np1,nt1,nm,nk),  &
         phot(nspec,np,nt),        &
         Ppar_IC(nspec,np,nt),        &
         Pressure_IC(nspec,np,nt), &
         PressurePar_IC(nspec,np,nt), &
         FAC_C(np,nt),             &
         Bmin_C(np,nt)                & ! minimum B field along each field line, 
        )                            ! passed from GM to IM, now as an output of IM

    ! Not clear why these need initialization !!!
    FAC_C = 0.0
    Pressure_IC = 0.0

    ! now allocate arrays for energy tracking
    allocate(eChangeOperator_IV(nspec,nOperator), driftin(nspec), &
         driftout(nspec),rbsumLocal(nspec),rbsumGlobal(nspec))

  end subroutine init_mod_crcm

end Module ModCrcm
