
subroutine calc_aurora

  use ModRIM
  use ModNumConst, only: cPi

  implicit none

  ! Variables that have an 'H' stuck on the end are only defined for a
  ! hemisphere
  real, dimension(0:nLons+1,nLats/2) :: rhoH, pH, TH, JrH, InvBH, eFluxH, AveEH
  real, dimension(0:nLons+1) :: OCFLBH
  integer :: iLon, iLat, iLR

  ! We need to basically work from the pole outwards when figuring out the
  ! auroral oval.  Since the iLat index goes from the pole out in the Southern
  ! hemisphere and the equator to pole in the North, we will reverse the North
  ! and call a single subroutine for both hemispheres.  So, we need to move the
  ! variables into a temporary place, call the routine, then move the results
  ! back into permanent variables.

  ! South
  do iLat = 1, nLats/2
     RhoH(:,iLat) = OuterMagRho(:,iLat)
     PH(:,iLat)   = OuterMagP(:,iLat)
     TH(:,iLat)   = OuterMagT(:,iLat)
     JrH(:,iLat)   = OuterMagJr(:,iLat)
     InvBH(:,iLat)   = OuterMagInvB(:,iLat)
  enddo

  call solve_for_aurora(RhoH, PH, TH, JrH, InvBH, OCFLBH, eFluxH, AveEH)

  ! Move results back into main variables
  do iLat = 1, nLats/2
     AveE(:,iLat)  = AveEH(:,iLat)
     eFlux(:,iLat) = eFluxH(:,iLat)
  enddo
  OCFLB(1,:) = OCFLBH

  ! North
  do iLat = 1, nLats/2
     ! iLR = iLat Reversed (from North pole to equator)
     iLR = nLats-iLat+1
     RhoH(:,iLat)  = OuterMagRho(:,iLR)
     PH(:,iLat)    = OuterMagP(:,iLR)
     TH(:,iLat)    = OuterMagT(:,iLR)
     JrH(:,iLat)   = OuterMagJr(:,iLR)
     InvBH(:,iLat) = OuterMagInvB(:,iLR)
  enddo

  call solve_for_aurora(RhoH, PH, TH, JrH, InvBH, OCFLBH, eFluxH, AveEH)

  ! Move results back into main variables
  do iLat = 1, nLats/2
     ! iLR = iLat Reversed (from North pole to equator)
     iLR = nLats-iLat+1
     AveE(:,iLR)  = AveEH(:,iLat)
     eFlux(:,iLR) = eFluxH(:,iLat)
  enddo
  OCFLB(2,:) = cPi - OCFLBH

end subroutine calc_aurora

!------------------------------------------------------------------------

subroutine solve_for_aurora(RhoH, PH, TH, JrH, InvBH, OCFLBH, eFluxH, AveEH)

  use ModRIM
  use ModParamRIM
  use ModNumConst, only: cDegToRad

  implicit none

  ! Variables that have an 'H' stuck on the end are only defined for a
  ! hemisphere
  real, dimension(0:nLons+1,nLats/2), intent(in)  :: rhoH, pH, TH, JrH, InvBH
  real, dimension(0:nLons+1,nLats/2), intent(out) :: eFluxH, AveEH
  real, dimension(0:nLons+1), intent(out) :: OCFLBH

  real, dimension(0:nLons+1,nLats/2) :: PolarRain_eFlux, PolarRain_AveE
  real, dimension(0:nLons+1,nLats/2) :: Discrete_eFlux, Discrete_AveE
  real, dimension(0:nLons+1,nLats/2) :: Diffuse_eFlux, Diffuse_AveE
  real, dimension(0:nLons+1) :: Width, smooth, Center

  integer :: iLon, iLat, nSmooth, iSubLon
  logical :: IsDone, IsPeakFound
  real :: MaxP
  real :: Discrete_FacAE, Discrete_FacEF
  real :: Diffuse_FacAE, Diffuse_FacEF

  Discrete_FacAE = 2.5e22
  Discrete_FacEF = 0.2e20
  Diffuse_FacAE = 4.0e-12
  Diffuse_FacEF = 0.3e6

  nSmooth = OCFLBSmoothLon/(Longitude(1,1) - Longitude(0,1))

  eFluxH = 0.0001
  AveEH  = 0.5
  OCFLBH = 1.0e32

  ! First, find open closed field-line boundary

  Width = 0.0
  do iLon = 0, nLons+1

     IsDone = .false.
     IsPeakFound = .false.
     MaxP = maxval(pH(iLon,:))

     iLat = 1

     do while (.not. IsDone)

        if (InvBH(iLon, iLat) > 0) then

           ! Set OCFLB to latitude of first open field-line
           if (OCFLBH(iLon) == 1.0e32) OCFLBH(iLon) = abs(Latitude(iLon,iLat))

           ! Find the peak location of the pressure - this may be the location
           ! of the inner edge of the plasma sheet.
           if (pH(iLon,iLat) == MaxP) IsPeakFound = .true.

           ! Determine the width of the oval.  We want it to be greater
           ! than some minimum width
           Width(iLon) = OCFLBH(iLon) - abs(Latitude(iLon,iLat))

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

!  !!! In order to actually smooth this correctly, we need to message pass
!  !!! the OCFLB and Width to all the PEs!!!!!!
!
!  do iLon = 0, nLons+1
!     smooth(iLon) = 0.0
!     do iSubLon = iLon-nSmooth, iLon+nSmooth
!        smooth(iLon) = smooth(iLon) + OCFLBH(mod(iSubLon+nLons,nLons))
!     enddo
!  enddo
!  OCFLBH = smooth/(2*nSmooth+1)
!
!  do iLon = 0, nLons+1
!     smooth(iLon) = 0.0
!     do iSubLon = iLon-nSmooth, iLon+nSmooth
!        smooth(iLon) = smooth(iLon) + Width(mod(iSubLon+nLons,nLons))
!     enddo
!  enddo
!  Width = smooth/(2*nSmooth+1)

  ! We want to put the Center of the diffuse aurora about 1/2 way between
  ! the open/closed field-line boundary and the inner edge of the 
  ! plasma sheet.  I don't know why.  It just seems like a good idea.
  Center = OCFLBH - Width/2

  ! ---------------------------
  ! Polar Rain

  PolarRain_AveE  = 0.0
  PolarRain_EFlux = 0.0

  do iLon = 0, nLons+1
     do iLat = 1, nLats/2
        if (abs(Latitude(iLon,iLat)) > OCFLBH(iLon)) then
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
     ! or approximately the inner edge of the plasma sheet.  The cos, takes
     ! into account the loss of electrons as a function of MLT due to 
     ! pitch angle scattering.  The nLons+1-iLon (i.e., mirroring the 
     ! location of the maximum pressure), is due to the electrons wanting
     ! to drift one way and the ions wanting to drift westward.


     MaxP = maxval(pH(nLons+1-iLon,:))

     write(*,*) "lon: ",iLon, nLons+1-iLon, MaxP, center(iLon), Width(iLon)

     Diffuse_EFlux(iLon,:) = &
          MaxP * Diffuse_FacEF * &
          exp(-abs(center(iLon)-abs(Latitude(iLon,1:nLats/2)))/Width(iLon))*&
          (0.375*cos(longitude(iLon,1:nLats/2))+0.625)

     where(tH(nLons+1-iLon,:) > 0) &
          Diffuse_AveE(iLon,:) = tH(nLons+1-iLon,:) * Diffuse_FacAE

  enddo

  AveEH = Diffuse_AveE
  EFluxH = Diffuse_EFlux

end subroutine solve_for_aurora
