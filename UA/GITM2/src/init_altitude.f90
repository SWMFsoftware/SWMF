
subroutine init_altitude

  use ModGITM
  use ModInputs
  use ModTime

  implicit none

  integer :: iAlt, i, iLoop, iAltInner

  logical :: IsDone

  real :: ScaleHeights(nAlts)
  real :: OlddHFactor, dHFactor
  real :: geo_lat, geo_lst, geo_lon, geo_alt, mm

  ! msis variables

  real, dimension(1:2) :: msis_temp
  real, dimension(1:8) :: msis_dens
  integer, dimension(25) :: sw
  real :: AP = 10.0

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

  sw = 1

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

        CALL GTD6(iJulianDay,utime,geo_alt/1000.0,geo_lat,geo_lon,geo_lst, &
             F107A,F107,AP,48,msis_dens,msis_temp)
        mm = (Mass(iO2_)*msis_dens(4) + &
             Mass(iO_)*msis_dens(2) + &
             Mass(iN2_)*msis_dens(3)) / &
             (msis_dens(4) + msis_dens(2) + msis_dens(3))
        ScaleHeights(iAlt) = Boltzmanns_Constant * msis_temp(2) / &
             (mm * Gravitational_Constant) 

!        geo_alt = geo_alt + dHFactor * ScaleHeights(iAlt)

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
  enddo

  Altitude(nAlts+2) = Altitude(nAlts+1) + dHFactor * ScaleHeights(nAlts)

  do iAlt = 0,nAlts+1
     dAlt(iAlt) = (Altitude(iAlt+1) - Altitude(iAlt-1))/2.0
  enddo
  dAlt(-1) = dAlt(0)
  dAlt(nAlts+2) = dAlt(nAlts+1)

end subroutine init_altitude
