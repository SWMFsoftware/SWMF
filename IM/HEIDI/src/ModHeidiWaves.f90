!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
Module ModHeidiWaves
  !\
  ! Wave-particle interaction variable definition module for HEIDI
  !/
  
  use ModHeidiSize, ONLY: nR, nT, nE, nPa, Ns, eng, slen
  
  ! Define plasma anisotropy variables
  real :: BFC(NR)
  real, allocatable :: EPME(:,:,:,:)
  real, allocatable :: EPMA(:,:,:,:)
  real :: EPP(NE,NS)
  real, allocatable :: ERNM(:,:,:,:)
  real :: ERNH(NE,NS)
  real, allocatable :: ATAW(:,:,:,:)
  real, allocatable :: GTAW(:,:,:,:)

  ! Define the quasi-linear diffusion coefficients
  real :: BWAVE(NR,NT)
  real :: PAbn(NPA)
  real :: AMLA(Slen)
  real :: ENOR(ENG,3)
  real, allocatable :: ZRpabn(:,:,:)
  real :: NDAA3(NPA)
  real :: NDAAJ(ENG,NPA,3)


contains

  subroutine init_mod_heidi_waves

    if(allocated(EPME)) RETURN
    allocate(EPME(NR,NT,NE,NPA))
    allocate(EPMA(NR,NT,NE,NPA))
    allocate(ERNM(NR,NT,NE,NPA))
    allocate(ATAW(NR,NT,NE,NPA))
    allocate(GTAW(NR,NT,NE,NPA))
    allocate(ZRpabn(NR,NPA,Slen))

  end subroutine init_mod_heidi_waves

  subroutine clean_mod_heidi_waves

    if(.not.allocated(EPME)) RETURN
    deallocate(EPME)
    deallocate(EPMA)
    deallocate(ERNM)
    deallocate(ATAW)
    deallocate(GTAW)
    deallocate(ZRpabn)

  end subroutine clean_mod_heidi_waves

end Module ModHeidiWaves
