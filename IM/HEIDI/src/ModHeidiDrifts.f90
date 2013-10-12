!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
Module ModHeidiDrifts
  !\
  ! Drift and diffusion variable definition module for the HEIDI program.
  ! Mike Liemohn, March 2006
  !/
  
  use ModHeidiSize, ONLY: nR, nT, NE, nPa, nS
  
  ! Define the atmospheric loss and charge exchange variables
  integer :: J6 = 0
  integer :: J18 = 0
  real    :: ACHARGE(NR,NT,NE,NPA,NS)=0.0 
  real    :: ATLOS(NR,NT,NE,NPA,NS) = 0.0
  
  ! Define the advection drift rate variables
  real :: VR(NR,NT,NE,NPA)  ! Radial drift
  real :: P1(NR,NT) 
  real :: P2(NR,NT,NE,NPA)= 0.0
  real :: EDOT(NR,NT,NE,NPA)
  real :: MUDOT(NR,NT,NE,NPA)
  real :: VrConv(NR,NT,NE,NPA) ! Convection velocity
  
  ! Define the pitch angle diffusion coefficient variables
  
  real :: COULE(NR,NT,NE,NPA,NS),COULI(NR,NT,NE,NPA,NS)
  real :: ATAE(NR,NT,NE,NPA,NS),BTAE(NR,NT,NE,NPA,NS),GTAE(NR,NT,NE,NPA+1,NS)
  real :: ATAI(NR,NT,NE,NPA,NS),BTAI(NR,NT,NE,NPA,NS),GTAI(NR,NT,NE,NPA+1,NS)
  
end Module ModHeidiDrifts
