
subroutine initialize_gitm(TimeIn)

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

  real(Real8_), intent(in) :: TimeIn
  integer :: iLat, iAlt, iBlock, iSpecies, iLon

  real :: TempAve
  real :: TempDiff
  real :: InvScaleHeightS(-1:nLons+2,-1:nLats+2)

  real :: LogRho(-1:nLons+2,-1:nLats+2), NewSumRho(-1:nLons+2,-1:nLats+2)

  real :: TempUnit_const, t, h

  logical :: IsThere, IsOk, IsDone, IsFirstTime = .true.

  if (.not.IsFirstTime) return

  IsFirstTime = .false.

  call start_timing("initialize")
  call report("initialize",1)

  ! ------------------------------------------------------------------------
  ! initialize stuff
  ! ------------------------------------------------------------------------

  if (TimeIn /= CurrentTime) then
     CurrentTime = TimeIn
     call time_real_to_int(CurrentTime, iTimeArray)
     call fix_vernal_time
  endif

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
        ! Cell interface is taken to be half way between cell centers
        ! so the cell size is half of the cell center distance
        ! between cells i-1 and i+1: 
        dAlt(iAlt) = (Altitude(iAlt+1) - Altitude(iAlt-1))/2.0
     enddo
     dAlt(-1) = dAlt(0)
     dAlt(nAlts+2) = dAlt(nAlts+1)

  endif

  ! This is the cell size and its inverse
  InvDAlt   = 1.0/dAlt

  ! This is the distance between cell centers. 
  ! Note that face(i) is between cells i and i-1 (like in BATSRUS)
  dAlt_F(0:nAlts+2) = Altitude(0:nAlts+2) - Altitude(-1:nAlts+1)
  dAlt_F(-1)        = dAlt_F(0)
  InvDAlt_F         = 1.0/dAlt_F

  RadialDistance = RBody + Altitude
  InvRadialDistance = 1.0/RadialDistance

  Gravity = &
       -Gravitational_Constant*&
       (RBody/RadialDistance) ** 2

  if (UseStretchedAltitude) then
     do iAlt=1,nAlts
        Gravity(iAlt) = &
             -Gravitational_Constant*&
             (RBody/((RadialDistance(iAlt)+RadialDistance(iAlt+1))/2)) ** 2
     enddo
  endif

!   Gravity =  -Gravitational_Constant

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

  ! This is done so we don't get a /0 below.

  dLatDist_GB = 1.0
  dLatDist_FB = 1.0
  dLonDist_GB = 1.0
  dLonDist_FB = 1.0

  ! Precalculate the physical size of cells in the Lat and Lon directions
  do iLat = 0, nLats+1
     do iAlt = -1, nAlts+2
        do iBlock = 1, nBlocks

           ! This is the cell size assuming that cell interface is half way
           dLatDist_GB(iLat, iAlt, iBlock) = 0.5 * &
                (Latitude(iLat+1, iBlock) - Latitude(iLat-1, iBlock)) * &
                RadialDistance(iAlt)           
           ! This is the distance between neighboring cells
           ! Note that face(i) is between cells i and i-1 (like in BATSRUS)
           dLatDist_FB(iLat, iAlt, iBlock) = &
                (Latitude(iLat, iBlock) - Latitude(iLat-1, iBlock)) * &
                RadialDistance(iAlt)           

           ! The longitude grid is uniform in angle
           do iLon = 0, nLons+1
              dLonDist_GB(iLon, iLat, iAlt, iBlock) = 0.5 * &
                   (Longitude(iLon+1,iBlock) - Longitude(iLon-1,iBlock)) * &
                   RadialDistance(iAlt)* &
                   max(abs(cos(Latitude(iLat,iBlock))),0.01)
              dLonDist_FB(iLon, iLat, iAlt, iBlock) = &
                   (Longitude(iLon,iBlock) - Longitude(iLon-1,iBlock)) * &
                   RadialDistance(iAlt)* &
                   max(abs(cos(Latitude(iLat,iBlock))),0.01)
           enddo

           dLonDist_FB(-1, iLat, iAlt, iBlock) = &
                dLonDist_FB(0, iLat, iAlt, iBlock)
           dLonDist_FB(nLons+2, iLat, iAlt, iBlock) = &
                dLonDist_FB(nLons+1, iLat, iAlt, iBlock)

        enddo
     enddo
  enddo

  ! Fill in 2nd ghost cells
  dLatDist_GB(-1, :, 1:nBlocks)      = dLatDist_GB(0, :, 1:nBlocks)
  dLatDist_GB(nLats+2, :, 1:nBlocks) = dLatDist_GB(nLats+1, :, 1:nBlocks)

  dLatDist_FB(-1, :, 1:nBlocks)      = dLatDist_FB(0, :, 1:nBlocks)
  dLatDist_FB(nLats+2, :, 1:nBlocks) = dLatDist_FB(nLats+1, :, 1:nBlocks)

  do iLon = -1, nLons+2
     dLonDist_GB(iLon, -1, :, 1:nBlocks)      = &
          dLonDist_GB(iLon, 0, :, 1:nBlocks)
     dLonDist_GB(iLon, nLats+2, :, 1:nBlocks) = &
          dLonDist_GB(iLon, nLats+1, :, 1:nBlocks)
     dLonDist_FB(iLon, -1, :, 1:nBlocks)      = &
          dLonDist_GB(iLon, 0, :, 1:nBlocks)
     dLonDist_FB(iLon, nLats+2, :, 1:nBlocks) = &
          dLonDist_GB(iLon, nLats+1, :, 1:nBlocks)
  enddo

  InvDLatDist_GB = 1.0/dLatDist_GB
  InvDLatDist_FB = 1.0/dLatDist_FB
  InvDLonDist_GB = 1.0/dLonDist_GB
  InvDLonDist_FB = 1.0/dLonDist_FB


  ! Precalculate the tangent of the latitude
  TanLatitude(:,1:nBlocks) = min(abs(tan(Latitude(:,1:nBlocks))),100.0) * &
       sign(1.0,Latitude(:,1:nBlocks))

  call init_heating_efficiency

  if (.not. DoRestart) then

     Velocity = 0.0
     IVelocity = 0.0
     VerticalVelocity = 0.0

     if (UseMsis .and. IsEarth) then
        call init_msis
     else

        TempUnit_const = 1. * Mass(1) / Boltzmanns_Constant
        TempAve  = (TempMax+TempMin)/2/TempUnit_const
        TempDiff = (TempMax-TempMin)/2/TempUnit_const

        do iBlock = 1, nBlocks

           do iAlt=-1,nAlts+2
              call get_temperature(0.0, 0.0, Altitude(iAlt), t, h)
              Temperature(:,:,iAlt,iBlock)  = t/TempUnit_const
              eTemperature(:,:,iAlt,iBlock) = t
              iTemperature(:,:,iAlt,iBlock) = t
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
                         Mass(iSpecies) / &
                         (Temperature(:,:,iAlt,iBlock)*TempUnit_const* &
                         Boltzmanns_Constant)
                 else
                    InvScaleHeightS = -Gravity(iAlt) * &
                         Mass(iSpecies) / &
                         (Temperature(:,:,iAlt,iBlock)*TempUnit_const* &
                         Boltzmanns_Constant)
                 endif

                 if(iAlt < 2)then
                    LogNS(:,:,iAlt,iSpecies,iBlock) = &
                         - (Altitude(iAlt)-AltMin)*InvScaleHeightS + &
                         LogNS0(iSpecies)
                 endif
                 if(iAlt > 1 .and. iAlt < nAlts+2)then
                    LogNS(:,:,iAlt,iSpecies,iBlock) = &
                         LogNS(:,:,iAlt-1,iSpecies,iBlock) &
                         -(Altitude(iAlt)-Altitude(iAlt-1))*InvScaleHeightS &
                         -(Temperature(:,:,iAlt+1,iBlock) - &
                          Temperature(:,:,iAlt-1,iBlock))/ &
                          (2.0*Temperature(:,:,iAlt,iBlock))
                 endif
                 if(iAlt == nAlts+2)then
                    LogNS(:,:,iAlt,iSpecies,iBlock) = &
                         LogNS(:,:,iAlt-1,iSpecies,iBlock) &
                         -(Altitude(iAlt)-Altitude(iAlt-1))*InvScaleHeightS &
                         -(Temperature(:,:,iAlt,iBlock) - &
                          Temperature(:,:,iAlt-1,iBlock))/ &
                          (Temperature(:,:,iAlt,iBlock))
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

           IDensityS(:,:,:,:,iBlock)    = 1.00e8
           IDensityS(:,:,:,ie_,iBlock)  = 1.00e8*(nIons-1)

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
