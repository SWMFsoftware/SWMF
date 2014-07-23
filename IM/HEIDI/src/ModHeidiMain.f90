!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModHeidiMain
  !\
  ! The main variable definition module for the HEIDI program.
  !/
  ! Include the parameter settings for MPI parallel processing
  
  use ModMpi
  use ModHeidiSize, ONLY: nR, nT, nS, nE, nPA, sLen, nPoint, dt
  use ModNumConst,  ONLY:cPi

  ! Basic physical constants: Earth radius and dipole coefficient
  real :: Re, DipoleFactor
  
  ! Define "constants" of the simulation
  real :: q,me,mp,dkp,fluxfact(NS),pi
  real :: dayr(48)
  real :: rkph(48)
  real :: f107r(48)
  real :: apr(48)
  real :: rsunr(48)
  
  ! Define independent variables and grid-related factors
  integer :: upa(NR)
  real :: DL1
  real :: DR
  real :: LZ(NR)
  real :: Z(NR)
  real :: BE(NR,Slen)
  real :: PHI(NT)
  real :: DPHI
  real :: MLT(NT)
  real    :: MAS(NS),M1(NS),WE(NE),DE(NE),EKEV(NE)
  real    :: V(NE,NS)
  real    :: VBND(NE,NS),MU(NPA),DMU(NPA),WMU(NPA),EBND(NE)
  real :: CONMU1
  real :: CONMU2
  real, allocatable :: FFACTOR(:,:,:,:)
  real, allocatable :: FACMU(:,:,:)
  real :: CONF1
  real :: CONF2
  real, allocatable :: CEDR(:,:,:,:,:)
  real, allocatable :: CIDR(:,:,:,:,:)
  
  ! Define flux variable, and a few others
  real, allocatable :: F2(:,:,:,:,:)
  real :: A, T=0.0, FGEOS(NT,NE,NPA,NS)
  
  ! Define parameters based on grid variables
  real :: ENER(NR,NS)
  real :: FACTOR(NS)
  real :: LEC(NR,NS)
  real :: ECOF(NR)
  real :: WCD(NR,NT)
  
  ! Define thermal plasma variables
  real :: xne(NR,NT)
  integer :: itherminit,ithermfirst
  
  
  ! Define parallel computing variables
  integer :: nParallelSpecies(NS)
  integer :: nSpecies, iSpecies

  real, allocatable :: funt(:,:,:)
  real, allocatable :: funi(:,:,:)
  real, allocatable :: BHeidi_III(:,:,:)
  real, allocatable :: SHeidi_III(:,:,:)
  real, allocatable :: RHeidi_III(:,:,:)
  real, allocatable :: bGradB1xHeidi_III(:,:,:)
  real, allocatable :: bGradB1yHeidi_III(:,:,:)
  real, allocatable :: bGradB1zHeidi_III(:,:,:)
  real, allocatable :: BxHeidi_III(:,:,:)
  real, allocatable :: ByHeidi_III(:,:,:)
  real, allocatable :: BzHeidi_III(:,:,:)
  real, allocatable :: pHeidi_III(:,:,:)
  real, allocatable :: rhoHeidi_III(:,:,:)
  real, allocatable :: Xyz_VIII(:,:,:,:)
  real, allocatable :: dEdt_IIII(:,:,:,:)
  real, allocatable :: VPhi_IIII(:,:,:,:)
  real, allocatable :: VR_IIII(:,:,:,:)
  real, allocatable :: dMudt_III(:,:,:,:)
  real, allocatable :: NeutralHydrogen(:,:,:)
  real, allocatable :: bFieldMagnitude_III(:,:,:)
  real, dimension(nT)             :: MhdEqPressure_I = -1.0, MhdEqDensity_I = -2.0
  logical                         :: IsBFieldNew
  real, parameter                 :: RadToDeg = 180.0/cPi

contains
  
  subroutine init_mod_heidi_main

    if(allocated(FFACTOR)) RETURN

    allocate(FFACTOR(NR,NT,NE,NPA))
    FFACTOR = 0.0
    allocate(FACMU(NPA,NR,NT))
    allocate(CEDR(NR,NT,NE,NPA,NS))
    allocate(CIDR(NR,NT,NE,NPA,NS))
    allocate(F2(NR,NT,NE,NPA,NS))
    F2 = 0.0
    allocate(funt(nPa,nR,nT))
    allocate(funi(nPa,nR,nT))
    allocate(BHeidi_III(nPoint,nR,nT))
    allocate(SHeidi_III(nPoint,nR,nT))
    allocate(RHeidi_III(nPoint,nR,nT))
    allocate(bGradB1xHeidi_III(nPoint,nR,nT))
    allocate(bGradB1yHeidi_III(nPoint,nR,nT))
    allocate(bGradB1zHeidi_III(nPoint,nR,nT))
    allocate(BxHeidi_III(nPoint,nR,nT))
    allocate(ByHeidi_III(nPoint,nR,nT))
    allocate(BzHeidi_III(nPoint,nR,nT))
    allocate(pHeidi_III(nPoint,nR,nT))
    allocate(rhoHeidi_III(nPoint,nR,nT))
    allocate(Xyz_VIII(3,nPoint,nR,nT))
    allocate(dEdt_IIII(nR,nT,nE,nPA))
    allocate(VPhi_IIII(nR,nT,nE,nPA))
    allocate(VR_IIII(nR,nT,nE,nPA))
    allocate(dMudt_III(nR,nT,nE,nPA))
    allocate(NeutralHydrogen(nR,nT,nPa))
    allocate(bFieldMagnitude_III(nPoint,nR,nT))

  end subroutine init_mod_heidi_main

  subroutine clean_mod_heidi_main

    if(.not.allocated(FFACTOR)) RETURN

    deallocate(FFACTOR)
    deallocate(FACMU)
    deallocate(CEDR)
    deallocate(CIDR)
    deallocate(F2)
    deallocate(funt)
    deallocate(funi)
    deallocate(BHeidi_III)
    deallocate(SHeidi_III)
    deallocate(RHeidi_III)
    deallocate(bGradB1xHeidi_III)
    deallocate(bGradB1yHeidi_III)
    deallocate(bGradB1zHeidi_III)
    deallocate(BxHeidi_III)
    deallocate(ByHeidi_III)
    deallocate(BzHeidi_III)
    deallocate(pHeidi_III)
    deallocate(rhoHeidi_III)
    deallocate(Xyz_VIII)
    deallocate(dEdt_IIII)
    deallocate(VPhi_IIII)
    deallocate(VR_IIII)
    deallocate(dMudt_III)
    deallocate(NeutralHydrogen)
    deallocate(bFieldMagnitude_III)

  end subroutine clean_mod_heidi_main
end module ModHeidiMain
