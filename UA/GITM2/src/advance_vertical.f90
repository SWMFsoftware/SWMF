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

  integer :: iIon, iSpecies, iAlt, iVar, iDim

!  KappaTemp1 = KappaTemp(iLon,iLat,:,iBlock)
  
  if (minval(NDensityS(iLon,iLat,:,1:nSpecies,iBlock)) <= 0.0) then
     write(*,*) "negative density found!"
     call stop_gitm("Can't Continue")
  endif

  Centrifugal = (cos(Latitude(iLat,iBlock)) * OmegaBody)**2
  Coriolis    = 2 * cos(Latitude(iLat,iBlock)) * OmegaBody
  do iAlt = -1, nAlts+2
     Heating(iAlt) = EuvHeating(iLon,iLat,iAlt,iBlock)
     LogRho(iAlt)  = log(Rho(iLon,iLat,iAlt,iBlock))
     Temp(iAlt)    = Temperature(iLon,iLat,iAlt,iBlock)
     do iVar = 1, nSpecies
        LogNS1(iAlt,iVar)  = log(NDensityS(iLon,iLat,iAlt,iVar,iBlock))
        VertVel(iAlt,iVar) = VerticalVelocity(iLon,iLat,iAlt,iVar,iBlock)
     end do
     do iVar = 1, nIonsAdvect
        LogINS(iAlt,iVar) = log(IDensityS(iLon,iLat,iAlt,iVar,iBlock))
     end do
     do iDim = 1, 3
        IVel(iAlt,iDim)    = IVelocity(iLon,iLat,iAlt,iDim,iBlock)
        Vel_GD(iAlt,iDim)  = Velocity(iLon,iLat,iAlt,iDim,iBlock)
     end do
  end do
  do iAlt = 0, nAlts+1
     cMax1(iAlt)   = cMax_GDB(iLon,iLat,iAlt,iUp_,iBlock)
  end do

  Lat = Latitude(iLat, iBlock) * 180.0/pi
  Lon = Longitude(iLon, iBlock) * 180.0/pi

!  write(*,*) "advance_vertical_1d", iLat, iLon, iBlock

  call advance_vertical_1d

  do iAlt = -1, nAlts+2
     Rho(iLon,iLat,iAlt,iBlock)         = exp(LogRho(iAlt))
     Temperature(iLon,iLat,iAlt,iBlock) = Temp(iAlt)
     do iDim = 1, 3
        Velocity(iLon,iLat,iAlt,iDim,iBlock) = Vel_GD(iAlt,iDim)
     end do
     do iVar = 1, nSpecies
        VerticalVelocity(iLon,iLat,iAlt,iVar,iBlock) = VertVel(iAlt, iVar)
        LogNS(iLon,iLat,iAlt,iVar,iBlock)            = LogNS1(iAlt, iVar)
     end do
  end do

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

  do iVar = 1, nSpecies; do iAlt = -1, nAlts+2
     nDensityS(iLon,iLat,iAlt,iVar,iBlock) = exp(LogNS1(iAlt, iVar))
  end do; end do
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

     do iVar = 1, nIonsAdvect; do iAlt = -1, nAlts+2
        IDensityS(iLon,iLat,iAlt, iVar, iBlock) = exp(LogINS(iAlt, iVar))
!!!!!     IDensityS(iLon,iLat,:,1:nIonsAdvect,iBlock) = LogINS
     end do; end do
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
