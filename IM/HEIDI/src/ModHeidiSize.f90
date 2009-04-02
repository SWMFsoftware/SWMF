Module ModHeidiSize

! Converted for HEIDI: Mike Liemohn, March 2006

! *** NOTE: This has to be used before any other ModHeidi# program ***

! Define the maximum array sizes for HEIDI
! Formerly: numv.h from the RAM version

!**  These parameters specify array sizes in HEIDI:
!**  NR		grid in the radial direction
!**  NT		grid in local time
!**  NE		grid in energy
!**  NS		number of species (e-, H+, He+, O+)
!**  NPA	grid in equatorial pitch angle
!**  NL		grid in L-shell for plane.f (Craig's thermal densities)
!**  NLT	grid in local time for plane.f
!**  ENG	grid in energy for wave coefficient variables
!**  Slen	grid in field-line distance
!**  NIJ	grid in energy for LANL boundary condition variables
!**  NSTH	number of thermal plasma species
!**  nthetacells  DGCPM output grid number in latitude
!**  nphicells	DGCPM output grid number in local time

      INTEGER NR,NT,NE,NS,NPA,NL,NLT,ENG,Slen,NIJ,NSTH
      INTEGER nthetacells,nphicells
      PARAMETER (NR=20)  ! =40 for "d0691" runs
      PARAMETER (NT=25)
      PARAMETER (NE=43)  ! =35 for old runs
      PARAMETER (NPA=72)
      PARAMETER (NS=4)
      PARAMETER (NL=35)
      PARAMETER (NLT=48)
      PARAMETER (ENG=70)
      PARAMETER (Slen=24)
      PARAMETER (NIJ=300)
      PARAMETER (NSTH=4)
      PARAMETER (nthetacells=60)
      PARAMETER (nphicells=120)

! Define array sizes as they are actually used within HEIDI
! Formerly: Common block PARAM1
	integer s,io,jo,ko,lo,iso,scalc(ns)
	real dt, DTMAX

end Module ModHeidiSize
