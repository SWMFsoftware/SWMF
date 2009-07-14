Module ModHeidiSize
  implicit none
  
  ! Converted for HEIDI: Mike Liemohn, March 2006
  ! *** NOTE: This has to be used before any other ModHeidi# program ***
  ! Define the maximum array sizes for HEIDI
  ! These parameters specify array sizes in HEIDI:
  
  integer, parameter :: NR = 20            !grid in the radial direction
  integer, parameter :: NT = 25            !grid in local time
  integer, parameter :: NE = 43            !grid in energy
  integer, parameter :: NPA = 72           !grid in equatorial pitch angle
  integer, parameter :: NS = 4             !number of species (e-, H+, He+, O+)
  integer, parameter :: NL = 35            !grid in L-shell for plane.f (Craig's thermal densities)
  integer, parameter :: NLT = 48           !grid in local time for plane.f
  integer, parameter :: ENG = 70           !grid in energy for wave coefficient variables
  integer, parameter :: Slen = 24          !grid in field-line distance
  integer, parameter :: NIJ = 300          !grid in energy for LANL boundary condition variables
  integer, parameter :: NSTH = 4           !number of thermal plasma species
  integer, parameter :: nthetacells = 60   !DGCPM output grid number in latitude
  integer, parameter :: nphicells = 120    !DGCPM output grid number in local time

  ! Define array sizes as they are actually used within HEIDI
  integer :: s,io,jo,ko,lo,iso,scalc(ns)
  real    :: dt, DTMAX,tmax

  !Define the minimum and maximum distance in the equatorial plane
  real, parameter :: RadiusMin = 1.5
  real, parameter :: RadiusMax = 7.0

end Module ModHeidiSize
