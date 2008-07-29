
module ModIonoHeidi

  !\
  ! Ionosphere array parameters
  !/
  integer, parameter :: IONO_nTheta = 65
  integer, parameter :: IONO_nPsi   =  4*(IONO_nTheta-1)+1

  !\
  ! Ionosphere solution array definitions
  !/
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::  &
    IONO_NORTH_PHI,IONO_SOUTH_PHI,               & !Ionospheric potential
    IONO_NORTH_Theta,IONO_NORTH_Psi,             & !
    IONO_SOUTH_Theta,IONO_SOUTH_Psi,             & !
    IONO_NORTH_JR,                               & !Ionospheric current
    IONO_SOUTH_JR,                               &
    IONO_NORTH_RCM_JR,                           & !RCM current
    IONO_SOUTH_RCM_JR

  real, parameter     :: IONO_Radius = 6372.0 * 1.0e3 + 110.0 * 1.0e3

  integer, parameter  :: maxfile = 3
  character (LEN=3)   :: plot_form(maxfile)
  character (len=100) :: plot_vars(maxfile), plot_vars1

  real, dimension(1:IONO_nTheta*2-1,1:IONO_nPsi) ::  &
       IonoGmVolume, IonoGmXPoint, IonoGmYPoint, &
       IonoGmBField, &
       IonoGmDensity=-1.0, &
       IonoGmPressure=-1.0, &
       IonoGmTemperature=-1.0

end module ModIonoHeidi
