Module ModHeidiWaves
  !\
  ! Wave-particle interaction variable definition module for HEIDI
  ! Mike Liemohn, March 2006
  !/

	use ModHeidiSize

! Define plasma anisotropy variables
! Formerly: Common block CANIS
	real BFC(NR),EPME(NR,NE,NPA),EPMA(NR,NE,NPA),EPP(NE,NS)
	real ERNM(NR,NE,NPA),ERNH(NE,NS),ATAW(NR,NT,NE,NPA)
	real GTAW(NR,NT,NE,NPA)

! Define the quasi-linear diffusion coefficients
! Formerly: Common block CWPI
	real BWAVE(NR,NT),PAbn(NPA),AMLA(Slen),ENOR(ENG,3)
	real ZRpabn(NR,NPA,Slen),NDAA3(NPA),NDAAJ(ENG,NPA,3)

end Module ModHeidiWaves