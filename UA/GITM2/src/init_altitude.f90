
!--------------------------------------------------------------
!
!--------------------------------------------------------------

subroutine get_temperature(lon, lat, alt, t, h)

  use ModInputs
  use ModPlanet

  implicit none

  real, intent(in) :: lon, lat, alt
  real, intent(out) :: t, h

  real    :: tAve, tDiff, n, r, g, m
  integer :: iSpecies

  if (UseMsis) then

     call get_msis_temperature(lon, lat, alt, t, h)

  else

     tAve  = (TempMax+TempMin)/2
     tDiff = (TempMax-TempMin)/2

     if (Alt <= TempHeight) then
        t = tAve + tDiff*tanh((alt - TempHeight)/TempWidth)
     else
        t = tAve + tDiff*tanh((alt - TempHeight)/TempWidth)
     endif

     r = RBody + alt
     g = Gravitational_Constant * (RBody/r) ** 2

     h = 0.0
     n = 0.0
     m = 0.0
     do iSpecies = 1, nSpecies
        m = m + exp(LogNS0(iSpecies)) * mass(iSpecies)
        n = n + exp(LogNS0(iSpecies))
     enddo

     m = m / n
     h = Boltzmanns_Constant * t / (m*g)

  endif

end subroutine get_temperature


!-----------------------------------------------------------------
!  This is the new init_altitude version.
!  It has been determined that the optimum resolution in the
!  vertical direction is about 0.3 scale heights.  So, we force
!  the user to resolve the grid at this resolution, independent
!  of where they actually put the high altitude boundary.
!-----------------------------------------------------------------

subroutine init_altitude

  use ModGITM
  use ModInputs
  use ModTime

  implicit none

  integer :: iAlt, i, iLoop, iAltInner

  logical :: IsDone

  real :: ScaleHeights(nAlts)
  real :: OlddHFactor, dHFactor
  real :: geo_lat, geo_lst, geo_lon, geo_alt, h, t

  !
  ! Here we are creating the altitude grid.  We are
  ! assuming that we have a lower boundary and an
  ! upper boundary and we basically want to scale the
  ! grid by the scale height between the two levels.
  !
  ! This is pretty tricky, since when you change the
  ! scaling factor, you change the scale heights that
  ! you are going to use.  So, you have to keep adjusting
  ! the scaling factor until you have the altitudes and
  ! the scale height that you want.
  !

  IsDone = .false.

  dHFactor = 0.3

  do iAlt=1,nAlts

     geo_lat = 0.0
     geo_lst = 12.0
     geo_lon = mod(geo_lst*15.0 - utime/3600.0*15.0 + 360.0,360.0)
     geo_alt = AltMin
     if (iAlt > 1) geo_alt = AltMin + sum(ScaleHeights(1:iAlt-1)) * dHFactor

     geo_lon = geo_lon * pi / 180.0

     call get_temperature(geo_lon, geo_lat, geo_alt, t, h)
     ScaleHeights(iAlt) = h

  enddo

  Altitude(1)  = AltMin
  Altitude(0)  = AltMin -     dHFactor * ScaleHeights(1)
  Altitude(-1) = AltMin - 2 * dHFactor * ScaleHeights(1)

  do iAlt=2,nAlts+1
     Altitude(iAlt) = altitude(iAlt-1) + dHFactor * ScaleHeights(iAlt-1)
     if (iDebugLevel > 3) write(*,*) "Altitude, dHFactor, ScaleHeight : ", &
          Altitude(iAlt), dHFactor, ScaleHeights(iAlt-1)
  enddo

  Altitude(nAlts+2) = Altitude(nAlts+1) + dHFactor * ScaleHeights(nAlts)

  do iAlt = 0,nAlts+1
     dAlt(iAlt) = (Altitude(iAlt+1) - Altitude(iAlt-1))/2.0
  enddo
  dAlt(-1) = dAlt(0)
  dAlt(nAlts+2) = dAlt(nAlts+1)

end subroutine init_altitude

!-----------------------------------------------------------------
!  This is the old init_altitude version.
!-----------------------------------------------------------------

subroutine init_altitude_old

  use ModGITM
  use ModInputs
  use ModTime

  implicit none

  integer :: iAlt, i, iLoop, iAltInner

  logical :: IsDone

  real :: ScaleHeights(nAlts)
  real :: OlddHFactor, dHFactor
  real :: geo_lat, geo_lst, geo_lon, geo_alt, h, t

  !
  ! Here we are creating the altitude grid.  We are
  ! assuming that we have a lower boundary and an
  ! upper boundary and we basically want to scale the
  ! grid by the scale height between the two levels.
  !
  ! This is pretty tricky, since when you change the
  ! scaling factor, you change the scale heights that
  ! you are going to use.  So, you have to keep adjusting
  ! the scaling factor until you have the altitudes and
  ! the scale height that you want.
  !

  IsDone = .false.

  dHFactor = 1.0
  OlddHFactor = 0.0

  iLoop = 1

  do while (.not. IsDone)

     do iAlt=1,nAlts

        geo_lat = 0.0
        geo_lst = 12.0
        geo_lon = mod(geo_lst*15.0 - utime/3600.0*15.0 + 360.0,360.0)
        geo_alt = AltMin
        if (iAlt > 1) geo_alt = AltMin + sum(ScaleHeights(1:iAlt-1)) * dHFactor

        geo_lon = geo_lon * pi / 360.0
        call get_temperature(geo_lon, geo_lat, geo_alt, t, h)

        ScaleHeights(iAlt) = h

     enddo

     if (abs(geo_alt - AltMax) < 1.0) then
        IsDone = .true.
     else
        if (OlddHFactor == 0.0) then
           OlddHFactor = dHFactor
           dHFactor = dHFactor*(AltMax - AltMin)/(geo_alt - AltMin) 
        else
           dHFactor = &
                dHFactor*(AltMax-AltMin)/(geo_alt-AltMin)/2.0 + &
                OlddHFactor/2.0
           OlddHFactor = dHFactor
        endif
     endif

     iLoop = iLoop+1

  enddo

  Altitude(1)  = AltMin
  Altitude(0)  = AltMin -     dHFactor * ScaleHeights(1)
  Altitude(-1) = AltMin - 2 * dHFactor * ScaleHeights(1)

  do iAlt=2,nAlts+1
     Altitude(iAlt) = altitude(iAlt-1) + dHFactor * ScaleHeights(iAlt-1)
     if (iDebugLevel > 3) write(*,*) "Altitude, dHFactor, ScaleHeight : ", &
          Altitude(iAlt), dHFactor, ScaleHeights(iAlt-1)
  enddo

  Altitude(nAlts+2) = Altitude(nAlts+1) + dHFactor * ScaleHeights(nAlts)

  do iAlt = 0,nAlts+1
     dAlt(iAlt) = (Altitude(iAlt+1) - Altitude(iAlt-1))/2.0
  enddo
  dAlt(-1) = dAlt(0)
  dAlt(nAlts+2) = dAlt(nAlts+1)

end subroutine init_altitude_old
