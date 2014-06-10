!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
Module ModHeidiCurrents
  !\
  ! Current and potential variable definition module for the HEIDI program.
  !/
  
  use ModHeidiSize, ONLY: nS, nR, nT, sLen
  
  ! Define variables for current densities
  real, dimension(NR+3) :: Lsh, Lats
  real, dimension(NR+3,NT,NS) :: Iphi,Irad,Jion1
  real, dimension(NR+3,NT) :: BASEPOT
  integer :: Ir

  ! Define variables for FAC-driven electric potential
  real    :: Latfac(NR+3),Lonfac(NT)
  real    :: Jfac(NR+3,NT)! = 0.0
  real    :: FPOT(NR+3,NT)
  integer :: Irfac, Ilfac
  
  ! Define variables for plasma pressure and other bulk quantities
  real :: PPER(NR,NT,NS),PPAR(NR,NT,NS)
  real :: RNHT(NR,NT,NS) = 0.0
  real :: EDEN(NR,NT,NS) = 0.0
  real :: ANIS(NR,NT,NS),EPAR(NR,NT,NS),Dst(NS)
  real :: NTOT(NS),ETOT(NS),JPER(NR,NT,NS)
  real :: Nspace(NR,NT,NS),Espace(NR,NT,NS)
  
  ! Define variables for the current calculation procedure
  real :: rl(NR+3,Slen),drl(NR+3,Slen),r2(3,NR+3,NT,Slen)
  real :: dBdrB(3,NR+3,NT,Slen),Rxy(NR+3,Slen),Bz(NR+3,Slen)
  real :: Bxy(NR+3,Slen),Bf2(NR+3,Slen),ds1(NR+3,Slen),ds2(NR+3,Slen)
  real :: ds(NR+3,Slen),beta1(NR+3,Slen),beta2(NR+3,Slen)
  real :: alpha1(NR+3,Slen),delR(NR+3,Slen),gam1(NT),gam2(NT)
  real :: dRm(NR+3),dR1(NR+3,Slen),dR2(NR+3,Slen),sp(NT),cp(NT)
  real :: sr(NR+3,Slen),cr(NR+3,Slen),sr3(NR+3,Slen),BBr(NR+3,Slen)
  real :: sl(NR+3),cl(NR+3),fac1(NR+3,Slen),fac2(NR+3,Slen)
  
  ! Define variables for the trig functions in the current calculation
  real :: sg1(NT),sg2(NT),cg1(NT),cg2(NT),sb1(NR+3,Slen)
  real :: sb2(NR+3,Slen),cb1(NR+3,Slen),cb2(NR+3,Slen),sa1(NR+3,Slen)
  real :: ca1(NR+3,Slen)
  
  ! Define variables for the indexing counters in the current calculation
  integer :: Ko2,Kmax(NR+3),j1(NT),j2(NT),i1(NR+3),i2(NR+3)
  integer :: k1(NR+3,Slen),k2(NR+3,Slen),ikk1(NR+3,Slen),ikk2(NR+3,Slen)
  integer :: ik1(NR+3,Slen),ik2(NR+3,Slen)
  
end Module ModHeidiCurrents
