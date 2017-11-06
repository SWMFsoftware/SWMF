Module ModCimi
  use ModCimiGrid
  use ModCimiPlanet,ONLY: nspec
  implicit none

  SAVE

  real    :: dt=1., dtmax=1. ! typical time step of cimi
  real    :: Time = 0.0
  logical :: UseMcLimiter=.false.,UseStrongDiff=.false.
  real    :: BetaLimiter = 1.5
  real    :: Pmin = 1e-2
  real, allocatable:: SDtime(:,:,:,:), f2(:,:,:,:,:)
  real, allocatable:: phot(:,:,:), Ppar_IC(:,:,:), &
       Pressure_IC(:,:,:), PressurePar_IC(:,:,:)
  real, allocatable:: FAC_C(:,:)
  real, allocatable:: Bmin_C(:,:)

  !Variable for E gain and loss
  real, allocatable :: driftin(:), driftout(:), &
       rbsumLocal(:), rbsumGlobal(:), rcsumLocal(:), rcsumGlobal(:)
  real, allocatable :: &
       eChangeOperator_VICI(:,:,:,:,:), pChangeOperator_VICI(:,:,:,:,:)
  real, allocatable :: eChangeLocal(:,:), eChangeGlobal(:,:)
  real, allocatable :: eTimeAccumult_ICI(:,:,:,:), pTimeAccumult_ICI(:,:,:,:)
  real, allocatable :: &
       preF(:,:,:,:), preP(:,:,:,:), Eje1(:,:,:) ! presipitation output
  integer, parameter :: &
       nOperator = 7, OpDrift_=1, OpBfield_=2, OpChargeEx_=3, &
       OpWaves_=4, OpStrongDiff_=5, OpLossCone_=6, OpLossCone0_=7
! Note order and number of operators has been changed Waves are added
! and OpLossCone0_=7 is previous in time OpLossCone (needed for
! precipitation)
                                       

! in CIMI: eChangeOperator=xle(ns,ir,ip,je+2) Here we create one array
!   for all operators (plus one dimension)
!   eTimeAccumulatv=esum(ns,ir,ip,je+2) No need for additional
!   dimention because it's total energy

  logical :: IsStandAlone=.false.
contains

  subroutine init_mod_cimi

    if(allocated(f2)) RETURN

    allocate(                         &
         SDtime(np,nt,nm,nk),         &
         f2(nspec,np1,nt1,nm,nk),     &
         phot(nspec,np,nt),           &
         Ppar_IC(nspec,np,nt),        &
         Pressure_IC(nspec,np,nt),    &
         PressurePar_IC(nspec,np,nt), &
         FAC_C(np,nt),                &
         Bmin_C(np,nt))
    ! minimum B field along each field line
    ! passed from GM to IM, now as an output of IM

    ! Not clear why these need initialization !!!
    FAC_C = 0.0
    Pressure_IC = 0.0

    ! now allocate arrays for energy tracking
    allocate(eChangeOperator_VICI(nspec,np,nt,neng+2,nOperator), &
         driftin(nspec), driftout(nspec), &
         rbsumLocal(nspec),rbsumGlobal(nspec), &
         rcsumLocal(nspec),rcsumGlobal(nspec))

    allocate(pChangeOperator_VICI(nspec,np,nt,neng+2,nOperator), &
         eChangeLocal(nspec,nOperator),eChangeGlobal(nspec,nOperator), &
         eTimeAccumult_ICI(nspec,np,nt,neng+2), &
         pTimeAccumult_ICI(nspec,np,nt,neng+2))
   
    allocate(preP(nspec,np,nt,neng+2), preF(nspec,np,nt,neng+2), &
         Eje1(nspec,np,nt))

  end subroutine init_mod_cimi

end Module ModCimi
