!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!**  These parameters specify array sizes in Vania's program:
!**  NR		grid in the radial direction
!**  NT		grid in local time
!**  NE		grid in energy
!**  NS		number of species (e-, H+, He+, O+)
!**  NPA		grid in equatorial pitch angle
!**  NL		grid in L-shell for plane.f (Craig's thermal densities)
!**  NLT		grid in local time for plane.f
!**  ENG		grid in energy for wave coefficient variables
!**  Slen	grid in field-line distance
!**  NIJ		grid in energy for LANL boundary condition variables
!**  NSTH	number of thermal plasma species
!**  nthetacells  DGCPM output grid number in latitude
!**  nphicells	DGCPM output grid number in local time

      INTEGER NR,NT,NE,NS,NPA,NL,NLT,ENG,Slen,NIJ,NSTH,nthetacells,nphicells
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
