!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModIonoDGCPM

! Ionospheric solution array
! Testing for the RIM -> DGCPM coupler

integer, parameter :: IONO_nTheta = 91
integer, parameter :: IONO_nPsi   = 181

  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::  &
    IONO_NORTH_PHI,IONO_SOUTH_PHI,               & !Ionospheric potential
    IONO_NORTH_Theta,IONO_NORTH_Psi,             & !
    IONO_SOUTH_Theta,IONO_SOUTH_Psi,             & !
    IONO_NORTH_JR,                               & !Ionospheric current
    IONO_SOUTH_JR,                               &
    IONO_NORTH_RCM_JR,                           & !RCM current
    IONO_SOUTH_RCM_JR

end module ModIonoDGCPM
