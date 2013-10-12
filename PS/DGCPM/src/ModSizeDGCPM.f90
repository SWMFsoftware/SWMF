!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
Module ModSizeDGCPM


! *** NOTE: This has to be used before any other Mod#DGCPM program ***

!**  These parameters specify array sizes in HEIDI:
!**  nthetacells    DGCPM output grid number in latitude
!**  nphicells	    DGCPM output grid number in local time
!**  nrcells        DGCPM output grid number of L shell distances
!**  thetamin       Minimum Theta (Co-Latitude) angle
!**  thatamax       Maximum Theta (Co-Latitude) angle
!**  dt             Model Timestep

  integer, parameter:: nthetacells=62
  integer, parameter:: nphicells=120
  integer, parameter:: nrcells = nthetacells
  real, parameter:: thetamin = 18.434949
  real, parameter:: thetamax = 61.385502
  real dt

end Module ModSizeDGCPM
