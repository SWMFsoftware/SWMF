
subroutine initialize_gitm

  use ModGITM
  use ModInputs
  use ModConstants
  use ModPlanet
  use ModRates
  use ModSphereInterface
  use ModTime
  use ModEUV
  implicit none

  type (UAM_ITER) :: r_iter

  integer :: iLat, iAlt, iBlock, iSpecies, iLon

  real :: TempAve
  real :: TempDiff
  real :: InvScaleHeightS(-1:nLons+2,-1:nLats+2)

  real :: LogRho(-1:nLons+2,-1:nLats+2), NewSumRho(-1:nLons+2,-1:nLats+2)

  logical :: IsThere, IsOk, IsDone, IsFirstTime = .true.

  if (.not.IsFirstTime) return

  IsFirstTime = .false.

  call start_timing("initialize")
  call report("initialize",1)

  ! ------------------------------------------------------------------------
  ! initialize stuff
  ! ------------------------------------------------------------------------

  call init_grid

  if (DoRestart) then
     call read_inputs("UA/restartIN/header.rst")
     call set_inputs
     call read_restart("UA/restartIN")
  endif

  call set_RrTempInd

  inquire(file='GITM.STOP',EXIST=IsThere)
  if (IsThere .and. iProc == 0) then
     open(iOutputUnit_, file = 'GITM.STOP', status = 'OLD')
     close(iOutputUnit_, status = 'DELETE')
  endif

  !\
  ! Initialize the EUV wave spectrum
  !/

  call init_euv

  Gravity = 0.0

  if (.not. DoRestart) then

     if (UseStretchedAltitude) then
        call init_altitude
     else
        ! Uniform grid
        dAlt = (AltMax-AltMin)/nAlts
        do iAlt=-1,nAlts+2
           Altitude(iAlt) = AltMin + (iAlt-0.5)*dAlt(iAlt)
        enddo
     endif

  else

     do iAlt = 0,nAlts+1
        dAlt(iAlt) = (Altitude(iAlt+1) - Altitude(iAlt-1))/2.0
     enddo
     dAlt(-1) = dAlt(0)
     dAlt(nAlts+2) = dAlt(nAlts+1)

  endif

  RadialDistance = RBody + Altitude
  Gravity = &
       -Gravitational_Constant*&
       (RBody/RadialDistance) ** 2

  if (UseStretchedAltitude) then
     do iAlt=1,nAlts
        Gravity = &
             -Gravitational_Constant*&
             (RBody/((RadialDistance(iAlt)+RadialDistance(iAlt+1))/2)) ** 2
     enddo
  endif

  if (iDebugLevel > 2) then
     do iAlt=-1,nAlts+2
        write(*,*) "===>Altitude : ", &
             iAlt, Altitude(iAlt), RadialDistance(iAlt), Gravity(iAlt)
     end do
  endif

!!!!!!
!  write(*,*) "GRAVITY is set to -10!!!!!!!"
!  Gravity = -10.0
!!!  Mass = 22.5 * AMU
!!!!!!

  TempUnit = Mass(1) / Boltzmanns_Constant
!  TempUnit = 1.0

  if (Is1D) then
     Latitude(0,1)  = Latitude(1,1) - 1.0 * pi/180.0
     Latitude(-1,1) = Latitude(0,1) - 1.0 * pi/180.0
     Latitude(2,1)  = Latitude(1,1) + 1.0 * pi/180.0
     Latitude(3,1)  = Latitude(2,1) + 1.0 * pi/180.0

     Longitude(0,1)  = Longitude(1,1) - 1.0 * pi/180.0
     Longitude(-1,1) = Longitude(0,1) - 1.0 * pi/180.0
     Longitude(2,1)  = Longitude(1,1) + 1.0 * pi/180.0
     Longitude(3,1)  = Longitude(2,1) + 1.0 * pi/180.0
  endif

  ! Precalculate the physical size of cells in the Lat and Lon directions
  do iLat = 0, nLats+1
     do iAlt = 0, nAlts+1
        do iBlock = 1, nBlocks
           dLatDist_GB(iLat, iAlt, iBlock) = 0.5 * &
                (Latitude(iLat+1, iBlock) - Latitude(iLat-1, iBlock)) * &
                RadialDistance(iAlt)
           dLonDist_GB(iLat, iAlt, iBlock) = &
                (Longitude(2,iBlock) - Longitude(1,iBlock)) * &
                RadialDistance(iAlt) * abs(cos(Latitude(iLat,iBlock)))
        enddo
     enddo
  enddo

  ! Precalculate the tangent of the latitude
  TanLatitude(:,1:nBlocks) = min(abs(tan(Latitude(:,1:nBlocks))),100.0) * &
       sign(1.0,Latitude(:,1:nBlocks))

  if (IsEarth) then

     !\
     ! EARTH Specific !!!
     !/
     do iAlt=1,nAlts
        HeatingEfficiency(iAlt) = &
             max(0.6-5.56e-5*(Altitude(iAlt)/1000.-165.)**2,0.12)

        if (altitude(iAlt)/1000. > 150.) then
           eHeatingEfficiency(iAlt)= &
                MIN(0.04+0.05*(altitude(iAlt)/1000.-150.)/100., 0.4)

           eHeatingEfficiency(iAlt)= 0.04


!!$     else if (altitude(iAlt)/1000. > 200.) then
!!$        eHeatingEfficiency(iAlt)= &
!!$             MIN(0.005+0.01*(altitude(iAlt)/1000.-130.)/100., 0.4)

        else
           eHeatingEfficiency(iAlt)= &
                MAX(0.05+0.07*(altitude(iAlt)/1000.-200.)/100., 0.000001)
        endif
             
     enddo

  else

     !\
     ! Mars Specific !!!
     !/

     HeatingEfficiency = 0.22

  endif

  if (.not. DoRestart) then

     TempAve  = (TempMax+TempMin)/2/TempUnit
     TempDiff = (TempMax-TempMin)/2/TempUnit

     Velocity = 0.0
     IVelocity = 0.0
     VerticalVelocity = 0.0

     TempMax = TempMin

     if (UseMsis .and. IsEarth) then

        call init_msis

     else

        do iBlock = 1, nBlocks

           do iAlt=-1,nAlts+2
              Temperature(:,:,iAlt,iBlock)  = TempAve &
                   + TempDiff*tanh((Altitude(iAlt) - TempHeight)/TempWidth)
           enddo
           
           do iAlt=-1,nAlts+2

              if (iAlt > -1) then
                 InvScaleHeight(:,:,iAlt,iBlock)  =  &
                      -Gravity(iAlt) / &
                      Temperature(:,:,iAlt,iBlock)
              else
                 InvScaleHeight(:,:,iAlt,iBlock)  =  &
                      -Gravity(iAlt)/Temperature(:,:,iAlt,iBlock)
              endif

              Rho(:,:,iAlt,iBlock) = 0.0

              NewSumRho = 0.0

              do iSpecies = 1, nSpecies

                 if (iAlt > -1) then
                    InvScaleHeightS = -Gravity(iAlt) * &
!!!                         22.5*Mass_Proton/ &
                         Mass(iSpecies) / &
                         (Temperature(:,:,iAlt,iBlock)*TempUnit* &
                         Boltzmanns_Constant)
                 else
                    InvScaleHeightS = -Gravity(iAlt) * &
!!!                         22.5*Mass_Proton/ &
                   Mass(iSpecies) / &
                         (Temperature(:,:,iAlt,iBlock)*TempUnit* &
                         Boltzmanns_Constant)
                 endif

                 if(iAlt < 1)then
                    LogNS(:,:,iAlt,iSpecies,iBlock) = &
                         - (Altitude(iAlt)-AltMin)*InvScaleHeightS + &
                         LogNS0(iSpecies)
                 else
                    LogNS(:,:,iAlt,iSpecies,iBlock) = &
                         LogNS(:,:,iAlt-1,iSpecies,iBlock) &
                         -(Altitude(iAlt)-Altitude(iAlt-1))*InvScaleHeightS &
                         -(Temperature(:,:,iAlt,iBlock) - &
                          Temperature(:,:,iAlt-1,iBlock))/ &
                          Temperature(:,:,iAlt,iBlock)
                 endif

                 NewSumRho      = NewSumRho + &
                      Mass(iSpecies)*exp(LogNS(:,:,iAlt,iSpecies,iBlock))

              enddo

              do iSpecies=1,nSpecies

                 NDensityS(:,:,iAlt,iSpecies,iBlock) = &
                      exp(LogNS(:,:,iAlt,iSpecies,iBlock))

                 NDensity(:,:,iAlt,iBlock) = NDensity(:,:,iAlt,iBlock) + &
                      NDensityS(:,:,iAlt,iSpecies,iBlock)

                 Rho(:,:,iAlt,iBlock) = Rho(:,:,iAlt,iBlock) + &
                      Mass(iSpecies) * NDensityS(:,:,iAlt,iSpecies,iBlock)

              enddo

           enddo
           
        enddo

     endif

     if (UseIRI .and. IsEarth) then

        call init_iri

     else

        do iBlock = 1, nBlocks

           IDensityS(:,:,:,:,iBlock)    = 0.01
           IDensityS(:,:,:,ie_,iBlock)  = 0.01*(nIons-1)

        enddo

     endif

  endif

  call init_b0

  if (IsEarth) call init_energy_deposition

  if (UseApex .and. IsEarth) then
     call report("subsolr",2)
     call SUBSOLR(iTimeArray(1),iJulianDay,iTimeArray(4),&
          iTimeArray(5),iTimeArray(6),SubsolarLatitude, &
          SubsolarLongitude)
  endif

  if (.not.Is1D) call exchange_messages_sphere

  call calc_pressure

!  do iBlock = 1, nBlocks
!     call calc_rates(iBlock)
!  enddo

  call end_timing("initialize")

end subroutine initialize_gitm
