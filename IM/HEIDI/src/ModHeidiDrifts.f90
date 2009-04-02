Module ModHeidiDrifts
  !\
  ! Drift and diffusion variable definition module for the HEIDI program.
  ! Mike Liemohn, March 2006
  !/

	use ModHeidiSize

! Define the atmospheric loss and charge exchange variables
! Formerly: Common block CE
	integer:: J6,J18
	real :: ACHAR(NR,NE,NPA,NS),ATLOS(NR,NE,NPA,NS)

! Define the advection drift rate variables
! Formerly: Common block CDR
	real :: VR(NR,NT),P1(NR,NT),P2(NR,NE,NPA),EDOT(NR,NT,NE,NPA)
	real :: MUDOT(NR,NT,NPA)

! Define the pitch angle diffusion coefficient variables
! Formerly: Common block CCOUL
	real :: COULE(NR,NE,NPA,NS),COULI(NR,NE,NPA,NS)
	real :: ATAE(NE,NPA,NS),BTAE(NE,NPA,NS),GTAE(NE,NPA+1,NS)
	real :: ATAI(NE,NPA,NS),BTAI(NE,NPA,NS),GTAI(NE,NPA+1,NS)

end Module ModHeidiDrifts
