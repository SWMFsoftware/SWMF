subroutine advance_vertical(iLon,iLat,iBlock)

  use ModGITM
  use ModPlanet, only: nSpecies, OmegaBody, nIonsAdvect
  use ModConstants, only: pi
  use ModSources, only: EUVHeating
  use ModInputs, only: UseIonAdvection, iDebugLevel
  use ModVertical, ONLY: &
       LogRho, &
       cMax1      => cMax,&
       KappaTemp1 => KappaTemp, &
       LogNS1     => LogNS, &
       Heating, &
       Vel_GD,  &
       Temp, &
       Centrifugal, Coriolis, &
       LogINS, &
       IVel, Lat, Lon, &
       VertVel
  implicit none

  integer, intent(in) :: iLon, iLat, iBlock

  integer :: iIon, iSpecies, iAlt

!  KappaTemp1 = KappaTemp(iLon,iLat,:,iBlock)
  
  if (minval(NDensityS(iLon,iLat,:,1:nSpecies,iBlock)) <= 0.0) then
     write(*,*) "negative density found!"
     call stop_gitm("Can't Continue")
  endif

  Heating     = EuvHeating(iLon,iLat,:,iBlock)
  Centrifugal = (cos(Latitude(iLat,iBlock)) * OmegaBody)**2
  Coriolis    = 2 * cos(Latitude(iLat,iBlock)) * OmegaBody
  LogRho  = log(Rho(iLon,iLat,:,iBlock))
  Vel_GD  = Velocity(iLon,iLat,:,:,iBlock)
  Temp    = Temperature(iLon,iLat,:,iBlock)
  LogNS1  = log(NDensityS(iLon,iLat,:,1:nSpecies,iBlock))
  VertVel = VerticalVelocity(iLon,iLat,:,1:nSpecies,iBlock)
  cMax1   = cMax_GDB(iLon,iLat,:,iUp_,iBlock)
  IVel    = IVelocity(iLon,iLat,:,:,iBlock)
  LogINS  = log(IDensityS(iLon,iLat,:,1:nIonsAdvect,iBlock))
!!!!  LogINS  = IDensityS(iLon,iLat,:,1:nIonsAdvect,iBlock)

  Lat = Latitude(iLat, iBlock) * 180.0/pi
  Lon = Longitude(iLon, iBlock) * 180.0/pi

!  write(*,*) "advance_vertical_1d", iLat, iLon, iBlock

  call advance_vertical_1d

  Rho(iLon,iLat,:,iBlock)                  = exp(LogRho)
  Velocity(iLon,iLat,:,:,iBlock)           = Vel_GD
  Temperature(iLon,iLat,:,iBlock)          = Temp
  LogNS(iLon,iLat,:,:,iBlock)              = LogNS1

  VerticalVelocity(iLon,iLat,:,1:nSpecies,iBlock) = VertVel

  if (minval(Temp) < 0.0) then
     write(*,*) "Temperature is negative!!!"
     do iAlt = -1,nAlts+2
        if (Temp(iAlt) < 0.0) then
           write(*,*) "iAlt : ", iAlt, Temp(iAlt)
        endif
     enddo
     call stop_gitm("Can't continue")
  endif

  if (Maxval(LogNS1) > 75.0) then
     write(*,*) "Maxval of LogNS too high!!!"
     do iAlt = -1,nAlts+2
        do iSpecies = 1, nSpecies
           if (LogNS1(iAlt,iSpecies) > 75.0) then
              write(*,*) "iSpecies, iBlock, Alt,Lon, Lat, maxval : ",&
                   iSpecies,iBlock,&
                   altitude(iAlt)/1000.0, &
                   longitude(iLon,iBlock)*180/pi, &
                   latitude(iLat,iBlock)*180/pi, LogNS1(iAlt,ispecies)
           endif
        enddo
     enddo
     call stop_gitm("Can't continue")
  endif

  nDensityS(iLon,iLat,:,1:nSpecies,iBlock) = exp(LogNS1)

  if (UseIonAdvection) then

!     if (Maxval(LogINS) > 75.0) then
!        write(*,*) "Maxval of Ion LogINS too high!!!"
!        do iAlt = -1,nAlts+2
!           do iSpecies = 1, nIonsAdvect
!              if (LogINS(iAlt,iSpecies) > 75.0) then
!                 write(*,*) "iSpecies, iBlock, Alt,Lon, Lat, maxval : ",&
!                      iSpecies,iBlock,&
!                      altitude(iAlt)/1000.0, &
!                      longitude(iLon,iBlock)*180/pi, &
!                      latitude(iLat,iBlock)*180/pi, LogINS(iAlt,ispecies)
!              endif
!           enddo
!        enddo
!        call stop_gitm("Can't continue")
!     endif
!
     IDensityS(iLon,iLat,:,1:nIonsAdvect,iBlock) = exp(LogINS)
!!!!!     IDensityS(iLon,iLat,:,1:nIonsAdvect,iBlock) = LogINS
     !\
     ! New Electron Density
     !/
     IDensityS(iLon,iLat,:,ie_,iBlock) = 0.0
     do iIon = 1, nIons-1
        IDensityS(iLon,iLat,:,ie_,iBlock) = &
             IDensityS(iLon,iLat,:,ie_,iBlock) + &
             IDensityS(iLon,iLat,:,iIon,iBlock)
     enddo
  endif

end subroutine advance_vertical
