
!\
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!/

subroutine init_b0

  use ModSize
  use ModGITM
  use ModInputs
  use ModTime

  implicit none

  real :: date, GeoLat, GeoLon, GeoAlt, ALat, ALon, LShell
  real :: xmag, ymag, zmag, r3, MagPot, bmag, MagneticPoleStrength
  integer :: iLat, iLon, iBlock, iAlt

  ! We need to fill in the Ghost Cells for the NonChanging Variables, such
  ! as the Magnetic Field:

  call report("init_b0",1)

  do iBlock=1,nBlocks
     do iAlt=-1,nAlts+2
        do iLat=-1,nLats+2
           do iLon=-1,nLons+2

              GeoLat = Latitude(iLat,iBlock)*180.0/pi
              GeoLon = Longitude(iLon,iBlock)*180.0/pi

              GeoAlt = altitude(iAlt)/1000.0

              if (GeoLat > 90.0) then
                 GeoLat = 180.0 - GeoLat
                 GeoLon = mod(GeoLon + 180.0,360.0)
              endif

              if (GeoLat < -90.0) then
                 GeoLat = -180.0 - GeoLat
                 GeoLon = mod(GeoLon + 180.0,360.0)
              endif

              alat = 0.0
              alon = 0.0
              xmag = 0.0
              ymag = 0.0
              zmag = 0.0

              call get_magfield(GeoLat,GeoLon,GeoAlt,alat,alon,xmag,ymag,zmag)

              mLatitude(iLon,iLat,iAlt,iBlock)  = alat
              mLongitude(iLon,iLat,iAlt,iBlock) = alon
              B0(iLon,iLat,iAlt,iEast_,iBlock)  = ymag
              B0(iLon,iLat,iAlt,iNorth_,iBlock) = xmag
              B0(iLon,iLat,iAlt,iUp_,iBlock)    = zmag
              B0(iLon,iLat,iAlt,iMag_,iBlock)   = &
                   sqrt(xmag*xmag + ymag*ymag + zmag*zmag)

              if (B0(iLon,iLat,iAlt,iMag_,iBlock) == 0.0) then
                 B0(iLon,iLat,iAlt,iMag_,iBlock) = 1.0e-10
              endif

              ! Calculate the magnetic dip angle, the magnetic 
              ! declination angle, and sines and cosines
              ! For now, only have positive sin of the dip angle 
              ! (and the magnetic dip angle)
              
              DipAngle(iLon,iLat,iAlt, iBlock) = &
                   atan (abs(zmag)/sqrt(xmag**2+ymag**2) )
              DecAngle(iLon,iLat,iAlt, iBlock) = atan2(ymag,xmag) * 180.0/pi

           enddo
        enddo
     enddo
  enddo

  if (UseApex .and. IsEarth) &
       call dypol( &
       MagneticPoleColat, MagneticPoleLon, MagneticPoleStrength)


  !    GyroFrequency_Electron(:,:,:,iBlock) = &
  !         Element_Charge * B0(:,:,:,Magnitude,iBlock) /  Mass_Electron

  
end subroutine init_b0

subroutine get_magfield(GeoLat,GeoLon,GeoAlt,alat,alon,xmag,ymag,zmag)

  use ModTime
  use ModPlanet
  use ModInputs
  use ModConstants, only: Pi

  implicit none

  real, intent(in) :: GeoLat, GeoLon, GeoAlt
  real, intent(out) :: alat, alon, xmag, ymag, zmag

  real :: date, bmag, LShell, r3, MagPot
  integer, external :: jday

  date = iStartTime(1) + float(iJulianDay)/float(jday(iStartTime(1),12,31))

  if (UseApex .and. IsEarth) then

     call APEX(DATE,GeoLat,GeoLon,GeoAlt,LShell, &
          alat,alon,bmag,xmag,ymag,zmag,MagPot)

     xmag =  xmag * 1.0e-9
     ymag =  ymag * 1.0e-9
     zmag = -zmag * 1.0e-9

  else

     r3 = (RBody / (RBody + GeoAlt)) ** 3

     LShell =  (RBody + GeoAlt) / RBody / (sin(pi/2 - GeoLat*pi/180.0))**2.0
     alat = acos(1.0/sqrt(LShell))*180.0/pi * sign(1.0,GeoLat)
     alon = GeoLon
     ymag = 0.0
     xmag =     - DipoleStrength * cos(GeoLat*pi/180.0) * r3
     zmag = 2.0 * DipoleStrength * sin(GeoLat*pi/180.0) * r3

  end if

end subroutine get_magfield
