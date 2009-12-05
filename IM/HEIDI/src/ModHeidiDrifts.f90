Module ModHeidiDrifts
  !\
  ! Drift and diffusion variable definition module for the HEIDI program.
  ! Mike Liemohn, March 2006
  !/

	use ModHeidiSize

! Define the atmospheric loss and charge exchange variables
! Formerly: Common block CE
	integer:: J6,J18
	real :: ACHAR(NR,NT,NE,NPA,NS),ATLOS(NR,NT,NE,NPA,NS)

! Define the advection drift rate variables
! Formerly: Common block CDR
	real :: VR(NR,NT,NE,NPA),P1(NR,NT),P2(NR,NT,NE,NPA),EDOT(NR,NT,NE,NPA)
	real :: MUDOT(NR,NT,NE,NPA), VrConv(NR,NT,NE,NPA)

! Define the pitch angle diffusion coefficient variables
! Formerly: Common block CCOUL
	real :: COULE(NR,NT,NE,NPA,NS),COULI(NR,NT,NE,NPA,NS)
	real :: ATAE(NR,NT,NE,NPA,NS),BTAE(NR,NT,NE,NPA,NS),GTAE(NR,NT,NE,NPA+1,NS)
	real :: ATAI(NR,NT,NE,NPA,NS),BTAI(NR,NT,NE,NPA,NS),GTAI(NR,NT,NE,NPA+1,NS)

end Module ModHeidiDrifts
