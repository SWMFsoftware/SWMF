Module ModSizeDGCPM

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

  integer, parameter:: NR=20  ! =40 for "d0691" runs
  integer, parameter:: NT=25
  integer, parameter:: NE=43  ! =35 for old runs
  integer, parameter:: NPA=72
  integer, parameter:: NS=4
  integer, parameter:: NL=35
  integer, parameter:: NLT=48
  integer, parameter:: ENG=70
  integer, parameter:: Slen=24
  integer, parameter:: NIJ=300
  integer, parameter:: NSTH=4
  integer, parameter:: nthetacells=60
  integer, parameter:: nphicells=120

  ! Define array sizes as they are actually used within HEIDI
  ! Formerly: Common block PARAM1

  integer s,io,jo,ko,lo,iso,scalc(ns)
  real dt

end Module ModSizeDGCPM
