!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
Module ModHeidiWaves
  !\
  ! Wave-particle interaction variable definition module for HEIDI
  !/
  
  use ModHeidiSize, ONLY: nR, nT, nE, nPa, Ns, eng, slen
  
  ! Define plasma anisotropy variables
  real :: BFC(NR),EPME(NR,NT,NE,NPA),EPMA(NR,NT,NE,NPA),EPP(NE,NS)
  real :: ERNM(NR,NT,NE,NPA),ERNH(NE,NS),ATAW(NR,NT,NE,NPA)
  real :: GTAW(NR,NT,NE,NPA)

  ! Define the quasi-linear diffusion coefficients
  real :: BWAVE(NR,NT),PAbn(NPA),AMLA(Slen),ENOR(ENG,3)
  real :: ZRpabn(NR,NPA,Slen),NDAA3(NPA),NDAAJ(ENG,NPA,3)

end Module ModHeidiWaves
