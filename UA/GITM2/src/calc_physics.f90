subroutine calc_physics(iBlock)

  use ModGITM
  use ModTime
  use ModConstants
  use ModEUV
  use ModPlanet
  use ModInputs
  implicit none

  integer, external :: jday

  integer, intent(in) :: iBlock

  integer :: itime(7)
  real*8 :: VernalDay
  real :: JulianDayEq, DaysInYear, OrbitAngle
  real :: SunDeclination, SinDec, CosDec

  integer :: iLon, iLat, iAlt

  call report("calc_physics",2)

  !\
  ! Solar Things
  !/

  OrbitAngle = 2.*pi*(CurrentTime - VernalTime)/SecondsPerYear

  SunOrbitEccentricity = SunOrbit_A                     + &
                         SunOrbit_B*cos(OrbitAngle)    + &
                         SunOrbit_C*sin(OrbitAngle)    + &
                         SunOrbit_D*cos(2.*OrbitAngle) + &
                         SunOrbit_E*sin(2.*OrbitAngle)

  SunDeclination = atan(tan(Tilt*pi/180.)*sin(OrbitAngle))

  SinDec = sin(SunDeclination)
  CosDec = cos(SunDeclination)

  LocalTime = mod((UTime/3600.0 + &
       Longitude(:,iBlock) * HoursPerDay / TwoPi), HoursPerDay)

  if (UseApex) &
       call SUBSOLR(iTimeArray(1),iJulianDay,iTimeArray(4),&
       iTimeArray(5),iTimeArray(6),SubsolarLatitude, &
       SubsolarLongitude)

  do iAlt=-1,nAlts+2

     !
     ! Compute Magnetic Local Time
     !

     if (UseApex .and. IsEarth) then
        do iLat=-1,nLats+2
           do iLon=-1,nLons+2
              call magloctm( &
                   MLongitude(iLon,iLat,iAlt,iBlock), &
                   SubsolarLatitude,   &
                   SubsolarLongitude,  &
                   MagneticPoleColat, &
                   MagneticPoleLon,   &
                   MLT(iLon,iLat,iAlt))
              if (mlt(iLon,iLat,iAlt) < 0) &
                   mlt(iLon,iLat,iAlt) = mlt(iLon,iLat,iAlt) + 24.0
           enddo
        enddo
     else

        do iLat=-1,nLats+2
           MLT(:,iLat,iAlt) = LocalTime
        enddo

        ! Since we go over the pole, 
        !we have to move the points to the proper location:

        if (Latitude(0,iBlock) < -pi/2.0) then
           MLT(:,-1,iAlt)      = MLT(:,2,iAlt) + 12.0
           MLT(:,0,iAlt)       = MLT(:,1,iAlt) + 12.0
        endif

        if (Latitude(0,iBlock) > pi/2.0) then
           MLT(:,nLats+1,iAlt) = MLT(:,nLats,iAlt) + 12.0
           MLT(:,nLats+2,iAlt) = MLT(:,nLats-1,iAlt) + 12.0
        endif

     endif

  enddo

!  write(*,*) mlt(1,1,1)

  where (MLT > 12.0)
     MLT = MLT - 24.0
  end where

  do iLon = 1, nLons
     do iLat = 1, nLats
        sza(iLon, iLat,iBlock) =  &
             acos(SinDec*sin(Latitude(iLat,iBlock)) + &
             CosDec*CosLatitude(iLat,iBlock) * &
             cos(pi*(LocalTime(iLon)-HoursPerDay/2)/(HoursPerDay/2)))

        if (DtLTERadiation < Rotation_Period) then
           call calc_avesza(iLon,iLat,iBlock, SinDec, CosDec)
        endif

     enddo
  enddo

  cosSZA = cos(SZA)
  sinSZA = sin(SZA)
  sqrtsinSZA = sqrt(abs(sinSZA))

  call calc_scaled_euv

  do iAlt = 1, nAlts
     xSolar(:,:,iAlt) = RadialDistance_GB(1:nLons,1:nLats,iAlt,iBlock) &
          * cosSZA(:,:,iBlock)
     ySolar(:,:,iAlt) = RadialDistance_GB(1:nLons,1:nLats,iAlt,iBlock) &
          * sinSZA(:,:,iBlock)
  enddo

end subroutine calc_physics


