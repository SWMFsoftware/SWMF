
subroutine calc_aurora

  use ModRIM
  use ModNumConst, only: cPi
  use ModProcIE

  implicit none

  ! Variables that have an 'H' stuck on the end are only defined for a
  ! hemisphere
  real, dimension(0:nLons+1,nLats/2) :: JrH, eFluxH, AveEH
  real, allocatable :: rhoH(:,:), pH(:,:), TH(:,:), InvBH(:,:), LatH(:,:)
  real, allocatable :: OCFLBH(:)
  integer :: iLon, iLat, iLR, iLonTo, iLonFrom

  allocate( &
       rhoH(0:nLonsAll+1,nLats/2), &
       pH(0:nLonsAll+1,nLats/2), &
       TH(0:nLonsAll+1,nLats/2), &
       InvBH(0:nLonsAll+1,nLats/2), &
       LatH(0:nLonsAll+1,nLats/2), &
       OCFLBH(0:nLonsAll+1))

  ! We need to basically work from the pole outwards when figuring out the
  ! auroral oval.  Since the iLat index goes from the pole out in the Southern
  ! hemisphere and the equator to pole in the North, we will reverse the North
  ! and call a single subroutine for both hemispheres.  So, we need to move the
  ! variables into a temporary place, call the routine, then move the results
  ! back into permanent variables.

  ! South
  do iLat = 1, nLats/2
     ! These variable could be mirrored in Longitude, so we need them over the
     ! whole domain
     RhoH(:,iLat)  = OuterMagRhoAllR(:,iLat)
     PH(:,iLat)    = OuterMagPAllR(:,iLat)
     TH(:,iLat)    = OuterMagTAllR(:,iLat)
     InvBH(:,iLat) = OuterMagInvBAllR(:,iLat)
     LatH(:,iLat)  = LatitudeAllR(:,iLat)
     ! These variables we only need for local domain
     JrH(:,iLat)   = OuterMagJr(:,iLat)
  enddo

  call solve_for_aurora(RhoH, PH, TH, JrH, InvBH, LatH, OCFLBH, eFluxH, AveEH)

  ! Move results back into main variables
  do iLat = 1, nLats/2
     AveE(:,iLat)  = AveEH(:,iLat)
     eFlux(:,iLat) = eFluxH(:,iLat)
  enddo

  ! North
  do iLat = 1, nLats/2
     ! iLR = iLat Reversed (from North pole to equator)
     iLR = nLats-iLat+1
     RhoH(:,iLat)  = OuterMagRhoAllR(:,iLR)
     PH(:,iLat)    = OuterMagPAllR(:,iLR)
     TH(:,iLat)    = OuterMagTAllR(:,iLR)
     InvBH(:,iLat) = OuterMagInvBAllR(:,iLR)
     LatH(:,iLat)  = LatitudeAllR(:,iLR)
     JrH(:,iLat)   = OuterMagJr(:,iLR)
  enddo

  call solve_for_aurora(RhoH, PH, TH, JrH, InvBH, LatH, OCFLBH, eFluxH, AveEH)

  ! Move results back into main variables
  do iLat = 1, nLats/2
     ! iLR = iLat Reversed (from North pole to equator)
     iLR = nLats-iLat+1
     AveE(:,iLR)  = AveEH(:,iLat)
     eFlux(:,iLR) = eFluxH(:,iLat)
  enddo

  do iLon = 0, nLons+1
     ! Need to shift longitudes
     iLonTo   = iLon
     iLonFrom = mod(iLon + iProc*nLons, nLons*nProc)
     if (iLonFrom == 0) iLonFrom = nLons*nProc
     OCFLB(1,iLonTo) = OCFLBH(iLonFrom)
     OCFLB(2,iLonTo) = cPi - OCFLBH(iLonFrom)
  enddo

!  OCFLB(2,:) = cPi - OCFLBH

  deallocate(rhoH,pH,TH,InvBH,LatH,OCFLBH)

end subroutine calc_aurora

!------------------------------------------------------------------------

subroutine solve_for_aurora(RhoH, PH, TH, JrH, InvBH, LatH, &
     OCFLBH, eFluxH, AveEH)

  use ModRIM
  use ModParamRIM
  use ModNumConst, only: cDegToRad
  use ModProcIE

  implicit none

  ! Variables that have an 'H' stuck on the end are only defined for a
  ! hemisphere
  real, dimension(0:nLonsAll+1,nLats/2), intent(in)  :: &
       rhoH, pH, TH, InvBH, LatH
  real, dimension(0:nLons+1,nLats/2), intent(in)  :: JrH
  real, dimension(0:nLons+1,nLats/2), intent(out) :: eFluxH, AveEH
  real, dimension(0:nLonsAll+1), intent(out) :: OCFLBH

  real, allocatable :: pNorm(:,:)

  real, dimension(0:nLons+1,nLats/2) :: PolarRain_eFlux, PolarRain_AveE
  real, dimension(0:nLons+1,nLats/2) :: Discrete_eFlux, Discrete_AveE, Discrete_K
  real, dimension(0:nLons+1,nLats/2) :: Diffuse_eFlux, Diffuse_AveE
  real, allocatable :: Width(:), smooth(:), Center(:)

  integer :: iLon, iLat, nSmooth, iSubLon, l, iLonOff, iLonG, iLonM
  logical :: IsDone, IsPeakFound
  real :: MaxP
  real :: Discrete_FacAE, Discrete_FacEF
  real :: Diffuse_FacAE, Diffuse_FacEF

  allocate( &
       pNorm(0:nLonsAll+1,nLats/2), &
       Width(0:nLonsAll+1), &
       Smooth(0:nLonsAll+1), &
       Center(0:nLonsAll+1))

  iLonOff = iProc*nLons

  Discrete_FacAE = 1.0e22
  Discrete_FacEF = 1.0e22
  Diffuse_FacAE = 5.0e-11
  Diffuse_FacEF = 1.0e9
  MinPressure = 5.0e-9
  OCFLBSmoothLon = 15.0*cDegToRad

  nSmooth = OCFLBSmoothLon/(Longitude(1,1) - Longitude(0,1))

  eFluxH = 0.0001
  AveEH  = 0.5
  OCFLBH = 1.0e32

  ! One of the problems with the MHD code is that the pressure can
  ! be miserably low.  So, let's set a minimum value for the maximum pressure,
  ! bringing the maximum up to this value, but keeping the shape.

  pNorm = pH

  if (maxval(pNorm) < MinPressure) &
       pNorm = pNorm/maxval(pNorm)*MinPressure

  ! First, find open closed field-line boundary

  Width = 0.0
  do iLon = 0, nLonsAll+1

     IsDone = .false.
     IsPeakFound = .false.
     MaxP = maxval(pNorm(iLon,:))

     iLat = 1

     do while (.not. IsDone)

        if (InvBH(iLon, iLat) > 0) then

           ! Set OCFLB to latitude of first open field-line
           if (OCFLBH(iLon) == 1.0e32) &
                OCFLBH(iLon) = abs(LatH(iLon,iLat))

           ! Find the peak location of the pressure - this may be the location
           ! of the inner edge of the plasma sheet.
           if (pNorm(iLon,iLat) == MaxP) IsPeakFound = .true.

           ! Determine the width of the oval.  We want it to be greater
           ! than some minimum width
           Width(iLon) = OCFLBH(iLon) - abs(LatH(iLon,iLat))

           if (IsPeakFound .and. Width(iLon) >= MinAuroralWidth) &
                IsDone = .true.

        endif

        iLat = iLat + 1

        if (iLat == nLats/2) then
           OCFLBH(iLon) = MaxAuroralLat
           Width(iLon) = MinAuroralWidth
           IsDone = .true.
        endif

     enddo

  enddo

  do iLon = 0, nLonsAll+1
     smooth(iLon) = 0.0
     do iSubLon = iLon-nSmooth, iLon+nSmooth
        smooth(iLon) = smooth(iLon) + OCFLBH(mod(iSubLon+nLonsAll,nLonsAll))
     enddo
  enddo
  OCFLBH = smooth/(2*nSmooth+1)

  do iLon = 0, nLonsAll+1
     smooth(iLon) = 0.0
     do iSubLon = iLon-nSmooth, iLon+nSmooth
        smooth(iLon) = smooth(iLon) + Width(mod(iSubLon+nLonsAll,nLonsAll))
     enddo
  enddo
  Width = smooth/(2*nSmooth+1)

  ! We want to put the Center of the diffuse aurora about 1/2 way between
  ! the open/closed field-line boundary and the inner edge of the 
  ! plasma sheet.  I don't know why.  It just seems like a good idea.
  Center = OCFLBH - Width/2

  ! ---------------------------
  ! Polar Rain

  PolarRain_AveE  = 0.5
  PolarRain_EFlux = 0.1

  do iLon = 0, nLons+1
     iLonG = iLon + iLonOff
     do iLat = 1, nLats/2
        if (abs(LatH(iLonG,iLat)) > OCFLBH(iLonG)) then
           PolarRain_AveE(iLon,iLat) = PolarRainAveE
           PolarRain_EFlux(iLon,iLat) = PolarRainEFlux
        endif
     enddo
  enddo

  ! ---------------------------
  ! Diffuse Aurora

  Diffuse_AveE = 0.0
  Diffuse_EFlux = 0.0

  do iLon = 0, nLons+1

     ! The exponential represents the radial location of the main aurora, 
     ! or approximately the inner edge of the plasma sheet.  The cos takes
     ! into account the loss of electrons as a function of MLT due to 
     ! pitch angle scattering.  The nLons+1-iLon (i.e., mirroring the 
     ! location of the maximum pressure), is due to the electrons wanting
     ! to drift one way and the ions wanting to drift westward.

     iLonG = iLon + iLonOff
     iLonM = nLonsAll+1 - iLonG

     MaxP = maxval(pNorm(iLonM,:))

     Diffuse_EFlux(iLon,:) = &
          MaxP * Diffuse_FacEF * &
          exp(-abs(center(iLonG)-abs(LatH(iLonG,:)))/Width(iLonG))*&
          (0.375*cos(longitude(iLon,1:nLats/2))+0.625)

     where(tH(iLonM,:) > 0) &
          Diffuse_AveE(iLon,:) = (tH(iLonM,:) * Diffuse_FacAE)**0.5

     ! The average energy can sometimes be quite concentrated near
     ! the open/closed field-line boundary, which is a problem, since
     ! the eflux is spread out over significant distances. (i.e., we have
     ! the situation in which there is massive amounts of low energy
     ! electrons precipitating...)

     l = maxloc(Diffuse_AveE(iLon,:),dim=1)

     ! Let's smooth it a little bit, but keep some of the original "edge"
     Diffuse_AveE(iLon,:) = &
          0.25 * Diffuse_AveE(iLon,:) + &
          0.75 * Diffuse_AveE(iLon, l) * &
          exp(-abs(center(iLonG)-abs(LatH(iLonG,:)))/Width(iLonG))

  enddo

  ! ---------------------------
  ! Discrete Aurora

  Discrete_AveE = 0.0
  Discrete_EFlux = 0.0
  Discrete_K = 0.0

  do iLon = 0, nLons+1
     
     iLonG = iLon + iLonOff
     iLonM = nLonsAll+1 - iLonG

     do iLat = 1,nLats/2
        if (pNorm(iLonG,iLat) > 0) &
             Discrete_K(iLon,iLat) = &
             (rhoH(iLonG,iLat)**1.5) / pNorm(iLonG,iLat)
        if (JrH(iLon,iLat) > 0) &
             Discrete_EFlux(iLon,iLat) = &
             (JrH(iLon,iLat)*1e6)*Discrete_K(iLon,iLat)
     enddo
  enddo

  Discrete_AveE = Discrete_EFlux*Discrete_FacAE
  Discrete_EFlux = (JrH*1e6)*Discrete_EFlux*Discrete_FacEF

!  AveEH = Diffuse_AveE
!  EFluxH = Diffuse_EFlux

  where(Diffuse_AveE < 0.25) Diffuse_AveE = 0.25
  where(Discrete_AveE < 0.25) Discrete_AveE = 0.25

  where(Diffuse_EFlux < 0.1) Diffuse_EFlux = 0.1
  where(Discrete_EFlux < 0.1) Discrete_EFlux = 0.1

  ! Let's weight the average energy by the number flux, which is ef/av
  AveEH = &
       (Diffuse_EFlux + Discrete_EFlux + PolarRain_EFlux) / ( &
       Diffuse_EFlux/Diffuse_AveE + &
       Discrete_EFlux/Discrete_AveE + &
       PolarRain_EFlux/PolarRain_AveE) 

  EFluxH = ( &
       Diffuse_EFlux/Diffuse_AveE + &
       Discrete_EFlux/Discrete_AveE + &
       PolarRain_EFlux/PolarRain_AveE) * AveEH
       
  deallocate(pNorm, Width, Smooth, Center)

end subroutine solve_for_aurora
