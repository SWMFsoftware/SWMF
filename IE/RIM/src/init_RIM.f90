subroutine init_RIM()

  use ModProcIE
  use ModRIM
  use ModNumConst
  use ModParamRIM
  use ModPlanetConst
  use ModMpi

  implicit none
  integer :: i, iLat, iLon, iSize, iLatFrom
  character (len=100), dimension(100) :: Lines
  integer :: iError = 0
  real :: rTmp

  ! Set actual grid values
  do i=1,nLats
     Latitude(:,i)=-cHalfPi+(cPi/nLats)*(real(i)-cHalf)
  end do
  do i=0,nLons+1
     Longitude(i,:)=((cTwoPi/nLons)/real(nProc))*(iProc*nLons+real(i)-cHalf) 
  end do

  ! Center Differencing except for the edges
  dLatitude(:,1)=(Latitude(:,2)-Latitude(:,1))
  do i=2,nLats-1
     dLatitude(:,i)=(Latitude(:,i+1)-Latitude(:,i-1))
  end do
  dLatitude(:,nLats)=(Latitude(:,nLats)-Latitude(:,nLats-1))

  ! Center Differencing except for the edges
  dLongitude(0,:)=(Longitude(1,:)-Longitude(0,:))
  do i=1,nLons
     dLongitude(i,:)=(Longitude(i+1,:)-Longitude(i-1,:))
  end do
  dLongitude(i,:)=(Longitude(nLons+1,:)-Longitude(nLons,:))

  Potential = cZero
  OldPotential = cZero
  Jr = cZero

  Radius = rPlanet_I(Planet_) + IonoHeightPlanet_I(Planet_)
  rTmp = Radius/rPlanet_I(Planet_)
  SMX = rTmp * cos(Latitude) * cos(Longitude + cPi)
  SMY = rTmp * cos(Latitude) * sin(Longitude + cPi)
  SMZ = rTmp * sin(Latitude)

  ! Initialize Empirical Models

  ! This assumes a grid which can be described by 2 1-D arrays:

  nEmpiricalLats = 0

  if (.not.DoSolve) nEmpiricalLats = nLats
  if (DoSolve .and. HighLatBoundary < cHalfPi) then
     do iLat = 1, nLats
        if (abs(Latitude(1,iLat)) > HighLatBoundary) &
             nEmpiricalLats = nEmpiricalLats + 1
     enddo
  endif

  if (nEmpiricalLats > 0) then 

     Lines(1) = "#BACKGROUND"
     Lines(2) = "IE/Input/"
     Lines(3) = NameEFieldModel
     Lines(4) = NameAuroralModel
     Lines(5) = NameSolarModel
     Lines(6) = ""
     Lines(7) = "#DEBUG"
     write(Lines(8),"(i2)") iDebugLevel 
     Lines(9) = "0"
     Lines(10) = ""
     Lines(11) = "#END"

     call EIE_set_inputs(Lines)
     call EIE_Initialize(iError)

     if (iError /= 0) then
        write(*,*) "Error : ", iError
        call stop_RIM("Stopping in init_RIM, call to EIE_Initialize")
     endif

     call IO_SetnMLTs(nLons+2)
     call IO_SetnLats(nEmpiricalLats)

     allocate(EmpiricalLatitude(0:nLons+1,nEmpiricalLats))
     allocate(EmpiricalMLT(0:nLons+1,nEmpiricalLats))
     allocate(EmpiricalPotential(0:nLons+1,nEmpiricalLats))
     allocate(EmpiricalAveE(0:nLons+1,nEmpiricalLats))
     allocate(EmpiricalEFlux(0:nLons+1,nEmpiricalLats))

     write(*,*) "rLats : ", minval(Latitude), &
          maxval(Latitude)

     write(*,*) "rLon : ", minval(Longitude), &
          maxval(Longitude)

     do iLon=0,nLons+1 
        nEmpiricalLats = 1
        do iLat=1,nLats
           if (abs(Latitude(iLon,iLat)) > HighLatBoundary) then
              EmpiricalLatitude(iLon,nEmpiricalLats) = &
                   Latitude(iLon,iLat) * 180.0 / cPi
              EmpiricalMLT(iLon,nEmpiricalLats) = &
                   Longitude(iLon,iLat) * 12.0 / cPi
              nEmpiricalLats = nEmpiricalLats + 1
           endif
        enddo
     enddo

     write(*,*) "Lats : ", minval(EmpiricalLatitude), &
          maxval(EmpiricalLatitude)

     write(*,*) "MLT : ", minval(EmpiricalMLT), &
          maxval(EmpiricalMLT)

     call IO_SetGrid(EmpiricalMLT, EmpiricalLatitude, iError)
     
     if (iError /= 0) then
        write(*,*) "Error : ", iError
        call stop_RIM("Stopping in init_RIM, call to IO_SetGrid")
     endif

  endif

  ! Need to allocate memory for arrays that span all of the IE module

  allocate( &
       LatitudeAll(nLats+2,nLons*nProc+1), &
       LongitudeAll(nLats+2,nLons*nProc+1), &
       PotentialAll(nLats+2,nLons*nProc+1), &
       SigmaHAll(nLats+2,nLons*nProc+1), &
       SigmaPAll(nLats+2,nLons*nProc+1), &
       LocalVar(nLats+2,nLons*nProc+1))

  ! This is ONLY for communication to other modules

  LatitudeAll = -1.0e32
  LongitudeAll = -1.0e32

  ! Fill Arrays

  ! This is the north pole
  LatitudeAll(1,iProc*nLons+1:iProc*nLons+nLons) = 0.0
  LongitudeAll(1,iProc*nLons+1:iProc*nLons+nLons) = &
       Longitude(1:nLons,nLats)
  if (iProc == nProc-1) then
     LatitudeAll((iProc+1)*nLons+1, 1) = 0.0
     LongitudeAll((iProc+1)*nLons+1,1) = Longitude(nLons+1,nLats)
  endif

  do iLat = 1, nLats
     iLatFrom = nLats-(iLat-1)
     LatitudeAll(iLat+1, iProc*nLons+1:iProc*nLons+nLons) = &
          cPi - Latitude(1:nLons,iLatFrom)
     LongitudeAll(iLat+1,iProc*nLons+1:iProc*nLons+nLons) = &
          Longitude(1:nLons,iLatFrom)
     if (iProc == nProc-1) then
        LatitudeAll(iLat+1,(iProc+1)*nLons+1) = cPi-Latitude(nLons+1,iLatFrom)
        LongitudeAll(iLat+1,(iProc+1)*nLons+1) = Longitude(nLons+1,iLatFrom)
     endif
  enddo

  ! This is the south pole
  LatitudeAll(nLats+2,iProc*nLons+1:iProc*nLons+nLons) = cPi
  LongitudeAll(nLats+2,iProc*nLons+1:iProc*nLons+nLons) = &
       Longitude(1:nLons,1)
  if (iProc == nProc-1) then
     LatitudeAll(nLats+2,(iProc+1)*nLons+1) = cPi
     LongitudeAll(nLats+2,(iProc+1)*nLons+1) = Longitude(nLons+1,1)
  endif

  iSize = (nLats+2) * (nProc*nLons+1)

  localVar = LatitudeAll
  call MPI_REDUCE(localVar, LatitudeAll, iSize, MPI_REAL, MPI_MAX, &
       0, iComm, iError)

  localVar = LongitudeAll
  call MPI_REDUCE(localVar, LongitudeAll, iSize, MPI_REAL, MPI_MAX, &
       0, iComm, iError)

end subroutine init_RIM
