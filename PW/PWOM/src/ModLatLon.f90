!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
Module ModLatLon
    
  !DESCRIPTION:
  ! Convert PW Lat and Lon (in SMG coords) to geo and mag lat and lon

  !USES:

  implicit none

  save    ! save all variables

  private ! except

  !PUBLIC MEMBER FUNCTIONS:
  public :: convert_lat_lon         ! Convert PW Lat and Lon to geographic 
                                    ! Lat and Lon
!  public :: convert_test            ! Unit test function 

  character(len=*),  parameter :: NamePwCoord='SMG'
contains
  !============================================================================
  
  subroutine convert_lat_lon(TimeSimulation, PwLat, PwLon, GeoLat, &
       GeoLon, GeoLat2, GeoLon2, GmLat, GmLon)
    use CON_axes, ONLY: transform_matrix
    use ModNumConst, ONLY: cPi, cTwoPi, cDegToRad, cRadToDeg

    real, intent(in) :: TimeSimulation, PwLat, PwLon
    real, intent(out):: GeoLat, GeoLon, GeoLat2,GeoLon2,GmLat, GmLon
    real             :: PwGg_DD(3,3), PwGm_DD(3,3), theta, phi 
    real             :: XyzPw_D(3), XyzGg_D(3), XyzGm_D(3)
    !--------------------------------------------------------------------------
    
    ! Get polar angle (theta) and azimuthal angle (phi)
    theta = 0.5*cPi - PwLat*cDegToRad
    phi   = PwLon*cDegToRad

    ! Get xyzPw_D from PwLat and PwLon
    xyzPw_D(1) = sin(theta)*cos(phi)
    xyzPw_D(2) = sin(theta)*sin(phi)
    xyzPw_D(3) = cos(theta)

    !\
    ! get equivalent geographic coords Hemisphere 1
    !/

    ! Get transform matrix 
    PwGg_DD = &
             transform_matrix(TimeSimulation, NamePwCoord, 'GEO')
    
    ! Transform xyzPw_D to XyzGg_DI
    XyzGg_D = matmul( PwGg_DD, XyzPw_D)
    
    ! Calculate GeoLat and GeoLon 
    GeoLon = modulo(atan2(XyzGg_D(2), XyzGg_D(1)), cTwoPi) * cRadToDeg
    GeoLat = 90.0 - (acos(max(-1.0,min(1.0, XyzGg_D(3))))*cRadToDeg)

    !\
    ! get equivalent geographic coords Hemisphere 2
    !/

    ! Get transform matrix 
    PwGg_DD = &
             transform_matrix(TimeSimulation, NamePwCoord, 'GEO')
    
    ! Transform xyzPw_D to XyzGg_DI
    XyzGg_D = matmul( PwGg_DD, (/XyzPw_D(1),XyzPw_D(2),-XyzPw_D(3)/))
    
    ! Calculate GeoLat and GeoLon 
    GeoLon2 = modulo(atan2(XyzGg_D(2), XyzGg_D(1)), cTwoPi) * cRadToDeg
    GeoLat2 = 90.0 - (acos(max(-1.0,min(1.0, XyzGg_D(3))))*cRadToDeg)

    

    !\
    ! get equivalent geomagnetic coords
    !/

    ! Get transform matrix 
    PwGm_DD = &
             transform_matrix(TimeSimulation, NamePwCoord, 'MAG')
    
    ! Transform xyzPw_D to XyzGg_DI
    XyzGm_D = matmul( PwGm_DD, XyzPw_D)
    
    ! Calculate GeoLat and GeoLon 
    GmLon = modulo(atan2(XyzGm_D(2), XyzGm_D(1)), cTwoPi) * cRadToDeg
    GmLat = 90.0 - (acos(max(-1.0,min(1.0, XyzGm_D(3))))*cRadToDeg)



  end subroutine convert_lat_lon
end Module ModLatLon
