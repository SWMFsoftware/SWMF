
subroutine calc_rates(iBlock)

  use ModGITM
  use ModRates
  use ModConstants
  use ModPlanet
  use ModInputs
  use ModEUV, only : SunOrbitEccentricity

  implicit none

  integer, intent(in) :: iBlock

  integer :: iAlt, iIon, iSpecies, iError, iiAlt

  real, dimension(nLons, nLats, nAlts) :: &
       Tn, Ti, Te, TWork1, TWork2, TWork3, NO2, &
       Ne, mnd, tmp, i_gyro, &
       Vi, Ve, MeVen, MeVei, MiVin, VeOe, ViOi, &
       Lambda_1, Lambda_2, &
       sin_dip, cos_dip, sin_dec, cos_dec

  real :: ScaleHeight(nLons, nLats)

  real :: e2

  call report("calc_rates",2)
  call start_timing("calc_rates")

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 4) write(*,*) "=====> mean major mass", iblock

  Ne  = IDensityS(1:nLons,1:nLats,1:nAlts,ie_,iBlock)
  mnd = NDensity(1:nLons,1:nLats,1:nAlts,iBlock)

  MeanIonMass = 0.0
  MeanMajorMass = 0.0
  do iSpecies = 1, nSpecies
     MeanMajorMass = MeanMajorMass + &
          Mass(iSpecies) * &
          NDensityS(1:nLons,1:nLats,1:nAlts,iSpecies,iBlock) / mnd
  enddo

  do iIon = 1, nIons-1
     MeanIonMass = MeanIonMass + &
          MassI(iIon) * IDensityS(1:nLons,1:nLats,1:nAlts,iIon,iBlock) / Ne
  enddo

  !\
  ! These are needed for the Euv Heating and other thermodynamics:
  !/

  if (iDebugLevel > 4) write(*,*) "=====> cp and kappatemp", iblock

  do iAlt = 0, nAlts+1

     cp(:,:,iAlt,iBlock) = 5.0 * NDensity(1:nLons,1:nLats,iAlt,iBlock) * &
          Boltzmanns_Constant / (2.0 * Rho(1:nLons,1:nLats,iAlt,iBlock))

     KappaTemp(:,:,iAlt,iBlock) = KappaTemp0 * &
          (Temperature(1:nLons,1:nLats,iAlt,iBlock)*TempUnit)**0.75

     iiAlt = iAlt
     if (iAlt == 0) iiAlt = 1
     if (iAlt == nAlts+1) iiAlt = nAlts

     ScaleHeight = &
          -Temperature(1:nLons,1:nLats,iAlt,iBlock) * TempUnit * &
          Boltzmanns_Constant / ( &
          Gravity(iAlt) * MeanMajorMass(1:nLons,1:nLats,iiAlt))

!     ! Add in an eddy diffusion term:
!     KappaTemp(:,:,iAlt,iBlock) = KappaTemp(:,:,iAlt,iBlock) + &
!          KappaEddy * (Pressure(1:nLons,1:nLats,iAlt,iBlock) / &
!          Pressure(1:nLons,1:nLats,1,iBlock))**2 * &
!          ScaleHeight**2

     ViscCoef(:,:,iAlt) = 4.5e-5 * &
          (Temperature(1:nLons,1:nLats,iAlt,iBlock)*TempUnit/ 1000.)**(-0.71)

  enddo

!  write(*,*) "mm cp:",minval(cp), maxval(cp), &
!       minval(Rho(1:nLons,1:nLats,iAlt,iBlock)*AMU/ &
!       NDensity(1:nLons,1:nLats,iAlt,iBlock)), &
!       maxval(Rho(1:nLons,1:nLats,iAlt,iBlock)*AMU/ &
!       NDensity(1:nLons,1:nLats,iAlt,iBlock))

  !\
  ! Need to get the neutral, ion, and electron temperature
  !/

  Tn = Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*TempUnit
  Ti = Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*TempUnit * 1.5
  Te = Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*TempUnit * 3.0

  if (UseIonChemistry) then

     write(*,*) "We have no chemistry now!!!"

  else

     RrTempDep = 1.0e-30

  endif

  !\
  ! -----------------------------------------------------------
  ! Collision Frequencies
  ! -----------------------------------------------------------
  !/

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 4) write(*,*) "=====> collisions", iblock

  e_gyro = &
       Element_Charge * B0(1:nLons,1:nLats,1:nAlts,iMag_,iBlock) /  &
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

  Collisions(:,:,:,iVEN_) = 5.4e-16 * (mnd)*sqrt(Te)

!
! Electron Ion Collision Frequency
!

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 4) write(*,*) "=====> vei", iblock

  tmp = (34.0 + 4.18*log((TE**3.0)/(Ne*1.0e-6)))
  Collisions(:,:,:,iVEI_) = tmp*Ne*TE**(-3.0/2.0) * 1.0e-6

!  Collisions(:,:,:,VEI) = (34.0 + 4.18*log((TE**3.0)/(Ne*1.0e-6))) &
!       * Ne * TE**(-3.0/2.0) * 1.0e-6
  
!
! Ion Gyrofrequency
!

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 4) write(*,*) "=====> igyro", iblock

  i_gyro = Element_Charge * B0(1:nLons,1:nLats,1:nAlts,iMag_,iBlock) / &
       MeanIonMass

  Vi = Collisions(:,:,:,iVIN_)
  Ve = Collisions(:,:,:,iVEN_) ! + Collisions(:,:,:,VEI)

  MeVen = Mass_Electron * Collisions(:,:,:,iVEN_)
  MeVei = Mass_Electron * Collisions(:,:,:,iVEI_)
  MiVin = MeanIonMass * Collisions(:,:,:,iVIN_)

  VeOe = Ve**2 + e_gyro**2
  ViOi = Vi**2 + i_gyro**2

!
! Here we calculate the conductivities. There are a number of different
! forms for there, but Song et al, 2001 show a pretty complete set which
! allows the current to be solved for (assuming that you have the electric
! field), independent of the neutral wind.
!  Conductivities are in mho/m=Siemens/m in mks units.  
! (In cgs, 1 S/m = 9e9 sec-1)

!
! This may be changed in a future version of the code.
!

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 4) write(*,*) "=====> sigmas", iblock

  Sigma_0 = e2 * Ne / MeVei

!  Alpha_Pedersen = Collisions(:,:,:,VEI) * (MiVin + MeVen) / &
!        (MeVen * Collisions(:,:,:,VEI) + MiVin * Ve)
!
!  Alpha_Hall = Collisions(:,:,:,VEI) * &
!       (MiVin - (Mass_Electron / MeanIonMass * MeVen)) / &
!       (MeVen * Collisions(:,:,:,VEI) + MiVin * Ve)
!
!  Sigma_Field_Aligned = Alpha_Pedersen * Sigma_0
!
!  Sigma_Pedersen = Sigma_Field_Aligned * &
!       (Collisions(:,:,:,VEI)*Collisions(:,:,:,VEI) / &
!       (Collisions(:,:,:,VEI)*Collisions(:,:,:,VEI) + &
!       Alpha_Hall*Alpha_Hall * e_Gyro*e_Gyro))
!
!  Sigma_Hall = Sigma_Field_Aligned * &
!       (Alpha_Hall*Collisions(:,:,:,VEI)*e_Gyro / &
!       (Collisions(:,:,:,VEI)*Collisions(:,:,:,VEI) + &
!       Alpha_Hall*Alpha_Hall * e_Gyro*e_Gyro))

  Sigma_Pedersen = ((1.0/MeVen) * (Ve*Ve/VeOe) + &
                    (1.0/MiVin) * (Vi*Vi/ViOi)) * Ne * e2

  Sigma_Hall = ((1.0/MeVen) * (Ve*e_gyro/VeOe) - &
                (1.0/MiVin) * (Vi*i_gyro/ViOi)) * Ne * e2

! Calculate the lambda's for the ion-drag term.  lambda is in sec-1.
!  In cgs units, lambda1=sigma_ped*B*B/rho, where sigma_ped is in s-1, 
!   B is in gauss,
!   and rho is in g/cm3
!  Conductivities are in mho/m=Siemens/m in mks units.  
!  (In cgs, 1 S/m = 9e9 sec-1)
!  B is in T in mks, and in gauss in cgs (1 T = 1.e4 gauss)
!  rho is in kg/m3 in mks, and in g/cm3 in cgs (1 kg/m3 = 1.e-3 g/cm3)
!   Thus, the conversion from mks to cgs is:  
!  9.e9 * 1.e4 * 1.e4 / 1.e-3 = 9.e20
!  The momentum equations are written in cgs in geographic coordinates.
!    lambda is multiplied by a velocity (cm/s in cgs or m/s in mks), so
!   no further conversion need be done after multiplying by 9.e20.
!
!  lambda_1 = sigma_Pedersen * B**2 / rho
!  lambda_2 = sigma_Hall * B**2 / rho
!  D = magnetic delination angle, I = magnetic dip angle
!  lambda_xx = lambda_1 * (1 - sinD*sinD*cosI*cosI)
!  lambda_yy = lambda_1 * (1 - cosD*cosD*cosI*cosI)
!  lambda_xy=lambda_yx = lambda_1*sinD*cosD*cosI*cosI + lambda_2*sinI

  sin_dip = B0(1:nLons,1:nLats,1:nAlts,iUp_,iBlock) / &
            B0(1:nLons,1:nLats,1:nAlts,iMag_,iBlock)
  cos_dip = sqrt(1.0 - sin_dip**2)

  tmp = sqrt(B0(1:nLons,1:nLats,1:nAlts,iEast_,iBlock)**2 + &
             B0(1:nLons,1:nLats,1:nAlts,iNorth_,iBlock)**2)
  sin_dec = B0(1:nLons,1:nLats,1:nAlts,iEast_,iBlock) / tmp
  cos_dec = B0(1:nLons,1:nLats,1:nAlts,iNorth_,iBlock) / tmp

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 4) write(*,*) "=====> lambdas", iblock

  Lambda_1 = &
       Sigma_Pedersen * &
       B0(1:nLons,1:nLats,1:nAlts,iMag_,iBlock)**2 / &
       Rho(1:nLons,1:nLats,1:nAlts,iBlock)

  Lambda_2 = &
       Sigma_Hall * B0(1:nLons,1:nLats,1:nAlts,iMag_,iBlock)**2 * &
       sin_dip / Rho(1:nLons,1:nLats,1:nAlts,iBlock)

  Lambda_xx = Lambda_1 * (1. - sin_dec**2 * cos_dip**2)

  Lambda_yy = Lambda_1 * (1. - cos_dec**2 * cos_dip**2)

  Lambda_xy = Lambda_1 * sin_dec * cos_dec + Lambda_2 * sin_dip

  call end_timing("calc_rates")

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 4) write(*,*) "=====> Done with calc_rates"

end subroutine calc_rates
