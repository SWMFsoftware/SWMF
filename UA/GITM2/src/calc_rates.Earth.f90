
subroutine calc_rates(iBlock)

  use ModGITM
  use ModRates
  use ModConstants
  use ModPlanet
  use ModInputs
  use ModEUV, only : SunOrbitEccentricity

  implicit none

  integer, intent(in) :: iBlock

  integer :: iAlt, iIon, iSpecies, iError, iiAlt, iLat,iLon

  real, dimension(nLons, nLats, nAlts) :: &
       Tn, Ti, TWork1, TWork2, TWork3, NO2

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2) :: &
       Ne, mnd, Te, tmp

  real :: ScaleHeight(nLons, nLats)

  real :: e2

  call report("calc_rates",2)
  call start_timing("calc_rates")

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 4) write(*,*) "=====> mean major mass", iblock

  Ne  = IDensityS(:,:,:,ie_,iBlock)

  ! We add 1 because this is in the denominator a lot, and the corners 
  ! don't have anything.
  mnd = NDensity(:,:,:,iBlock)+1.0

  MeanIonMass = 0.0
  MeanMajorMass = 0.0
  do iSpecies = 1, nSpecies
     MeanMajorMass = MeanMajorMass + &
          Mass(iSpecies) * &
          NDensityS(:,:,:,iSpecies,iBlock)/mnd
  enddo

  ! Once again, in the corners, the meanmajormass is 0.
  where (MeanMajorMass == 0) MeanMajorMass = Mass(1)

  do iIon = 1, nIons-1
     MeanIonMass = MeanIonMass + &
          MassI(iIon) * IDensityS(:,:,:,iIon,iBlock) / Ne
  enddo

  TempUnit = MeanMajorMass / Boltzmanns_Constant       

  !\
  ! These are needed for the Euv Heating and other thermodynamics:
  !/

  if (iDebugLevel > 4) write(*,*) "=====> cp and kappatemp", iblock

  do iAlt = 0, nAlts+1
     
     ! The Vibration-2 is a convertion from cp to cv.  Cv should be
     ! used when in an altitude coordinate system.  Cp is for a 
     ! pressure based system.

     cp(:,:,iAlt,iBlock) = 0.0
     Gamma(:,:,iAlt,iBlock) = 0.0  

     do iSpecies = 1, nSpecies
        cp(:,:,iAlt,iBlock) = cp(:,:,iAlt,iBlock) + &
             (Vibration(iSpecies)-2) * &
             NDensityS(1:nLons,1:nLats,iAlt,iSpecies,iBlock) * &
             (Boltzmanns_Constant / Mass(iSpecies))

!!! Gamma -1 = sum(n_i * (gamma_i -1))/n_total
!!! where gamma_i -1 = 2 / (Vibration -2) 

        Gamma(1:nLons,1:nLats,iAlt,iBlock) = &
             gamma(1:nLons,1:nLats,iAlt,iBlock) + &
             NDensityS(1:nLons,1:nLats,iAlt,iSpecies,iBlock) / & 
             (Vibration(iSpecies)-2)    

     enddo

     cp(:,:,iAlt,iBlock) = cp(:,:,iAlt,iBlock) / &
          (2.0 * NDensity(1:nLons,1:nLats,iAlt,iBlock))

     gamma(1:nLons,1:nLats,iAlt,iBlock) = &
          gamma(1:nLons,1:nLats,iAlt,iBlock) *2.0/ &   
          NDensity(1:nLons,1:nLats,iAlt,iBlock) + 1 

!     if (Is1D .and. UseKappa1DCorrection) then
!        KappaTemp(:,:,iAlt,iBlock) = KappaTemp0 * &
!             (Temperature(1:nLons,1:nLats,iAlt,iBlock) * &
!             TempUnit(1:nLons,1:nLats,iAlt) / &
!             Kappa1DCorrectionFactor)**Kappa1dCorrectionPower
!     else
         KappaTemp(:,:,iAlt,iBlock) = &
             (NDensityS(1:nLons,1:nLats,iAlt,iO2_,iBlock) / &
             NDensity(1:nLons,1:nLats,iAlt,iBlock) + &
             NDensityS(1:nLons,1:nLats,iAlt,iN2_,iBlock)/ &
             NDensity(1:nLons,1:nLats,iAlt,iBlock)) * 5.6e-4 * &
             (Temperature(1:nLons,1:nLats,ialt,iBlock)*TempUnit(1:nLons,1:nLats,iAlt)) &
             **0.75 + &
             (NDensityS(1:nLons,1:nLats,iAlt,iO_3P_,iBlock)/&
             NDensity(1:nLons,1:nLats,iAlt,iBlock)*7.59e-4) * &
             (Temperature(1:nLons,1:nLats,iAlt,iBlock) * &
             TempUnit(1:nLons,1:nLats,iAlt))**0.75
        
!        KappaTemp(:,:,iAlt,iBlock) = KappaTemp0 * &
!             (Temperature(1:nLons,1:nLats,iAlt,iBlock) * &
!             TempUnit(1:nLons,1:nLats,iAlt))**0.75
     endif

     iiAlt = iAlt
     if (iAlt == 0) iiAlt = 1
     if (iAlt == nAlts+1) iiAlt = nAlts

     ScaleHeight = &
          -Temperature(1:nLons,1:nLats,iAlt,iBlock) * &
          TempUnit(1:nLons,1:nLats,iAlt) * &
          Boltzmanns_Constant / ( &
          Gravity_GB(1:nLons,1:nLats,iAlt,iBlock) &
          * MeanMajorMass(1:nLons,1:nLats,iiAlt))

     do iLat = 1, nLats
        do iLon = 1, nLons

              KappaTemp(iLon,iLat,iAlt,iBlock) = &
                   KappaTemp(iLon,iLat,iAlt,iBlock) + &
                   KappaEddyDiffusion(iLon,iLat,iAlt,iBlock) * cp(iLon,iLat,iAlt,iBlock) * &
                   Rho(iLon,iLat,iAlt,iBlock)

        enddo
     enddo

     ViscCoef(:,:,iAlt) = 4.5e-5 * &
          (Temperature(1:nLons,1:nLats,iAlt,iBlock)*&
          TempUnit(1:nLons,1:nLats,iAlt)/ 1000.)**(-0.71)

  enddo

!  write(*,*) "mm cp:",minval(cp), maxval(cp), &
!       minval(Rho(1:nLons,1:nLats,iAlt,iBlock)*AMU/ &
!       NDensity(1:nLons,1:nLats,iAlt,iBlock)), &
!       maxval(Rho(1:nLons,1:nLats,iAlt,iBlock)*AMU/ &
!       NDensity(1:nLons,1:nLats,iAlt,iBlock))

  !\
  ! Need to get the neutral, ion, and electron temperature
  !/

  Tn = Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*TempUnit(1:nLons,1:nLats,1:nAlts)
  Ti = ITemperature(1:nLons,1:nLats,1:nAlts,iBlock)

  !\
  ! -----------------------------------------------------------
  ! Collision Frequencies
  ! -----------------------------------------------------------
  !/

  e_gyro = &
       Element_Charge * B0(:,:,:,iMag_,iBlock) /  &
       Mass_Electron

  e2 = Element_Charge * Element_Charge

!
! Ion Neutral Collision Frequency (From Kelley, 1989, pp 460):
!

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 4) write(*,*) "=====> vin",iblock

  Collisions(:,:,:,iVIN_) = 2.6e-15 * (mnd + Ne)/sqrt(MeanMajorMass/AMU)

!
! Electron Neutral Collision Frequency
!

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 4) write(*,*) "=====> ven", iblock

  Te = eTemperature(:,:,:,iBlock)
  where(te == 0.0) te = 1000.0
  Collisions(:,:,:,iVEN_) = 5.4e-16 * (mnd)*sqrt(Te)


!
! Electron Ion Collision Frequency
!

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 4) write(*,*) "=====> vei", iblock

  tmp = (34.0 + 4.18*log((TE**3.0)/(Ne*1.0e-6)))
  Collisions(:,:,:,iVEI_) = tmp*Ne*TE**(-3.0/2.0) * 1.0e-6

!!  Collisions(:,:,:,VEI) = (34.0 + 4.18*log((TE**3.0)/(Ne*1.0e-6))) &
!!       * Ne * TE**(-3.0/2.0) * 1.0e-6
!  

  i_gyro = Element_Charge * B0(:,:,:,iMag_,iBlock) / &
       MeanIonMass

  call end_timing("calc_rates")

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 4) write(*,*) "=====> Done with calc_rates"

end subroutine calc_rates
