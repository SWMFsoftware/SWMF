
subroutine calc_GITM_sources(iBlock)

  use ModInputs
  use ModSources
  use ModGITM

  implicit none

  integer, intent(in) :: iBlock

  integer :: iAlt, iError, iDir, iLat, iLon, iSpecies
  integer :: iiAlt, iiLat, iiLon
  real :: tmp(nLons, nLats, nAlts+1)
  real :: tmp2(nLons, nLats, nAlts), diff(nLons, nLats, nAlts)
  real :: tmp3(nLons, nLats, -1:nAlts+2)
  real :: tmp4(-1:nLons+2, -1:nLats+2, -1:nAlts+2)
  real :: tmp5(-1:nLons+2, -1:nLats+2, -1:nAlts+2)
  real :: Omega(nLons, nLats, nAlts)
  real :: RhoI(nLons, nLats, nAlts)
  real :: ReducedMass(nLons, nLats, nAlts)
  real :: KappaEddyDiffusion(-1:nLons+2, -1:nLats+2, -1:nAlts+2)
  real :: ScaleHeight(-1:nLons+2, -1:nLats+2, -1:nAlts+2)
  real :: TmpGradient(nLons, nLats, nAlts,3)

  real :: NSoverN(-1:nLons+2, -1:nLats+2, -1:nAlts+2)
  real :: NSoverNGradient(1:nLons, 1:nLats, 1:nAlts, 3)

  real :: oldtemp(-1:nLons+2, -1:nLats+2, -1:nAlts+2)

  real :: dtsub, dttotal, dtsubmax, dtr

  call report("calc_GITM_sources",1)

  ! calc_rate is used to determine reaction rates, heating coefficients,
  ! ion-neutral collision frequency, lambdas for ion drag, etc.

  if (iDebugLevel > 4) write(*,*) "=====> going into calc_rates", iproc

  ChemicalHeatingRate = 0.0

  call calc_rates(iBlock)

  RhoI = IDensityS(1:nLons,1:nLats,1:nAlts,ie_,iBlock) * &
       MeanIonMass(1:nLons,1:nLats,1:nAlts)

  !\
  ! Gravity is a source term which is calculated in initialize.f90
  !/

  !\
  ! ---------------------------------------------------------------
  ! These terms are for Neutral Temperature
  ! ---------------------------------------------------------------
  !/

  !\
  ! Solar Heating -------------------------------------------------
  !/

  if (iDebugLevel > 4) write(*,*) "=====> solar heating", iproc

  if (UseSolarHeating .or. UseIonChemistry) then

     ! So far, calc_physics only has stuff that is needed for solar
     ! euv, such as solar zenith angles, and local time.

     call calc_physics(iBlock)
     call euv_ionization_heat(iBlock)

  endif

  if (.not. UseSolarHeating) EuvHeating = 0.0

  ! The auroral heating is specified below, after the aurora is described
  ! in get_potential

  !\
  ! Joule Heating --------------------------------------------------
  !/

  ! If you have turned off IonDrag, then Joule Heating should NOT be
  ! included either.  This is because the Neutrals can become HIGHLY
  ! seperated from the ions, causing MASSIVE Joule heating!!!

  if (UseJouleHeating .and. UseIonDrag) then

     tmp2 = Collisions(1:nLons,1:nLats,1:nAlts,iVIN_) * &
          RhoI(1:nLons,1:nLats,1:nAlts)/ &
          Rho(1:nLons,1:nLats,1:nAlts,iBlock) * &
          (MeanIonMass(1:nLons,1:nLats,1:nAlts)/AMU) /  &
          (MeanIonMass(1:nLons,1:nLats,1:nAlts)/AMU + &
          MeanMajorMass(1:nLons,1:nLats,1:nAlts)/AMU)

     do iDir = 1, 3

        JouleHeating(:,:,:) = JouleHeating(:,:,:) + tmp2 * &
             (IVelocity(1:nLons,1:nLats,1:nAlts,iDir,iBlock) - &
             Velocity(1:nLons,1:nLats,1:nAlts,iDir,iBlock))**2

!        write(*,*) "Joule Heating : ",iDir, &
!             maxval(JouleHeating) / TempUnit / minval(cp(:,:,1:nAlts,iBlock)), &
!             maxval(abs(IVelocity(1:nLons,1:nLats,1:nAlts,iDir,iBlock))), &
!             maxval(abs(Velocity(1:nLons,1:nLats,1:nAlts,iDir,iBlock)))

     enddo

     JouleHeating = JouleHeating / TempUnit / cp(:,:,1:nAlts,iBlock)

  else

     JouleHeating = 0.0

  endif

  !\
  ! Conduction ----------------------------------------------------
  !/

  if (iDebugLevel > 4) write(*,*) "=====> conduction", iproc

  if(UseConduction)then
     tmp2 = Rho(1:nLons, 1:nLats,1:nAlts, iBlock) * &
          cp(1:nLons, 1:nLats,1:nAlts, iBlock)
     call calc_conduction(iBlock, &
          Temperature(1:nLons, 1:nLats,-1:nAlts+2, iBlock), &
          KappaTemp(1:nLons, 1:nLats,0:nAlts+1, iBlock), &
          tmp2, &
          Conduction)
  else
     Conduction = 0.0
  end if

  !\
  ! Diffusion ----------------------------------------------------
  !/

  if (iDebugLevel > 4) write(*,*) "=====> diffusion ", iproc

  if (UseDiffusion)then
     do iLon = -1, nLons+2
        do iLat = -1, nLats+2
           do iAlt = -1, nAlts+2
              iiLon = min(max(iLon,1),nLons)
              iiLat = min(max(iLat,1),nLats)
              iiAlt = min(max(iAlt,1),nAlts)
              ScaleHeight(iLon,iLat,iAlt) = &
                   -Temperature(iLon,iLat,iAlt,iBlock) * TempUnit * &
                   Boltzmanns_Constant / ( &
                   Gravity(iAlt) * MeanMajorMass(iiLon,iiLat,iiAlt))
           enddo
        enddo
     enddo
     ! tmp3 is kappa
     do iAlt = -1, nAlts+2
        KappaEddyDiffusion(:,:,iAlt) = &
!             2.5e3 * (Pressure(:,:,iAlt,iBlock) / &
             1.0e4 * (Pressure(:,:,iAlt,iBlock) / &
             Pressure(:,:,1,iBlock))
     enddo

!
!     KappaEddyDiffusion = KappaEddyDiffusion/ScaleHeight
!
     do iSpecies = 1, nSpecies

        do iLon = -1, nLons+2
           do iLat = -1, nLats+2
              do iAlt = -1, nAlts+2
                 iiLon = min(max(iLon,1),nLons)
                 iiLat = min(max(iLat,1),nLats)
                 iiAlt = min(max(iAlt,1),nAlts)
                 tmp4(iLon,iLat,iAlt) = &
                      (Mass(iSpecies) / &
                      MeanMajorMass(iiLon, iiLat, iiAlt) - 1.0)
              enddo
           enddo
        enddo


!        tmp3(1:nLons, 1:nLats,1:nAlts) = KappaEddyDiffusion(1:nLons, 1:nLats,1:nAlts) * MeanMajorMass
!        tmp3(:,:,0) = KappaEddyDiffusion(1:nLons, 1:nLats,0)*&
!             MeanMajorMass(:,:,1)
!        tmp3(:,:,nAlts+1) = KappaEddyDiffusion(1:nLons, 1:nLats,nAlts+1)*&
!             MeanMajorMass(:,:,nAlts)
!

        NsOverN = NDensityS(:,:,:,iSpecies, iBlock) / &
             NDensity(:,:,:,iBlock)

        tmp2 = 1.0
!        call calc_conduction(iBlock, &
!             NDensityS(1:nLons, 1:nLats, -1:nAlts+2, iSpecies, iBlock), &
!             KappaEddyDiffusion(1:nLons, 1:nLats,0:nAlts+1), &
!             tmp2, &
!             Diff)
!        Diffusion(1:nLons,1:nLats,1:nAlts,iSpecies) = diff

        call calc_conduction(iBlock, &
             NsOverN(1:nLons, 1:nLats, -1:nAlts+2), &
             KappaEddyDiffusion(1:nLons, 1:nLats,0:nAlts+1), &
             tmp2, &
             Diff)
        Diffusion(1:nLons,1:nLats,1:nAlts,iSpecies) = &
             -diff * NDensity(1:nLons,1:nLats,1:nAlts,iBlock)

!        Diffusion(1:nLons,1:nLats,1:nAlts,iSpecies) = &
!             tmp4(1:nLons,1:nLats,1:nAlts) * &
!             diff * NDensityS(1:nLons,1:nLats,1:nAlts,iSpecies,iBlock)

!        call calc_conduction(iBlock, &
!             NDensityS(1:nLons, 1:nLats, -1:nAlts+2, iSpecies, iBlock), &
!             tmp3, &
!             tmp2, &
!             Diff)
!
!        diff = 0.0

!        Diffusion(1:nLons,1:nLats,1:nAlts,iSpecies) = 0.0

!        call UAM_Gradient(tmp5, tmpGradient, iBlock)
!

!        do iAlt = 1, 10
!           write(*,*) iSpecies, iAlt, Diffusion(1,1,iAlt,iSpecies), &
!                NDensityS(1, 1, iAlt, iSpecies, iBlock), &
!                tmp2(1,1,iAlt), KappaEddyDiffusion(1,1,iAlt)
!        enddo

!
        call UAM_Gradient(KappaEddyDiffusion, tmpGradient, iBlock)
!        
!        do iAlt = 1, 10
!           write(*,*) iSpecies, iAlt, Diffusion(1,1,iAlt,iSpecies), &
!                diff(1,1,iAlt), NDensityS(1, 1, iAlt, iSpecies, iBlock), &
!                tmp2(1,1,iAlt), KappaEddyDiffusion(1,1,iAlt), &
!                tmpGradient(1,1,iAlt, iUp_), tmp4(1,1,iAlt), tmp5(1,1,iAlt)
!        enddo
!
!        Diffusion(:,:,:,iSpecies) = Diffusion(:,:,:,iSpecies) + &
!             NDensityS(1:nLons, 1:nLats, 1:nAlts, iSpecies, iBlock) * &
!             tmpGradient(1:nLons, 1:nLats, 1:nAlts, iUp_)*&
!             tmp4(1:nLons, 1:nLats, 1:nAlts)
!
!        do iAlt = 1, 10
!           write(*,*) iSpecies, iAlt, Diffusion(1,1,iAlt,iSpecies), &
!                diff(1,1,iAlt), NDensityS(1, 1, iAlt, iSpecies, iBlock), &
!                tmp2(1,1,iAlt), KappaEddyDiffusion(1,1,iAlt), &
!                tmpGradient(1,1,iAlt, iUp_), tmp4(1,1,iAlt), tmp5(1,1,iAlt)
!        enddo

     enddo

!        NSoverN = NDensityS(:, :, :, iSpecies, iBlock) / &
!             NDensity(:, :, :, iBlock)
!        call UAM_Gradient(NSoverN, NSoverNGradient, iBlock)
!        tmp2 = 1.0e-6 / NSoverNGradient(1:nLons,1:nLats,1:nAlts, iUp_)
!        tmp2 = -1.0e-4


     do iSpecies = 1, nSpecies
        do iAlt = 1, nAlts
           do iLat = 1, nLats
              do iLon = 1, nLons
                 if (abs(Diffusion(iLon,iLat,iAlt,iSpecies)) > &
                      0.1*NDensityS(iLon,iLat,iAlt,iSpecies,iBlock)) then
write(*,*) "Limiting : ", Diffusion(iLon,iLat,iAlt,iSpecies)
                    Diffusion(iLon,iLat,iAlt,iSpecies) = &
                         max(Diffusion(iLon,iLat,iAlt,iSpecies), &
                         -0.1*NDensityS(iLon,iLat,iAlt,iSpecies,iBlock))
                    Diffusion(iLon,iLat,iAlt,iSpecies) = &
                         min(Diffusion(iLon,iLat,iAlt,iSpecies), &
                         0.1*NDensityS(iLon,iLat,iAlt,iSpecies,iBlock))
write(*,*) "         : ", Diffusion(iLon,iLat,iAlt,iSpecies)
                 endif
              enddo
           enddo
        enddo
     enddo
  else
     Diffusion = 0.0
  end if

  !\
  ! ---------------------------------------------------------------
  ! These terms are for Neutral Winds
  ! ---------------------------------------------------------------
  !/

  if (UseIonDrag) then

     tmp2 = Collisions(1:nLons,1:nLats,1:nAlts,iVIN_)*&
          RhoI/Rho(1:nLons,1:nLats,1:nAlts,iBlock)

     do iDir = 1, 3
        IonDrag(:,:,:,iDir) = tmp2 * &
             (IVelocity(1:nLons,1:nLats,1:nAlts,iDir,iBlock) - &
             Velocity(1:nLons,1:nLats,1:nAlts,iDir,iBlock))
     enddo

  else

     IonDrag = 0.0

  endif

  !\
  ! ---------------------------------------------------------------
  ! These terms are for Vertical Neutral Wind drag
  ! ---------------------------------------------------------------
  !/

  if (UseNeutralFriction) then

     call calc_neutral_friction(iBlock)

  else

     NeutralFriction = 0.0

  endif

  !\
  ! Viscosity ----------------------------------------------------
  !/

  if (iDebugLevel > 4) write(*,*) "=====> conduction", iproc

  if(UseViscosity)then
     call calc_conduction(iBlock, &
          Velocity(1:nLons, 1:nLats,-1:nAlts+2, iNorth_, iBlock), &
          ViscCoef(1:nLons, 1:nLats,0:nAlts+1), &
          Rho(1:nLons, 1:nLats,1:nAlts, iBlock), &
          Viscosity(1:nLons, 1:nLats,1:nAlts, iNorth_))
     call calc_conduction(iBlock, &
          Velocity(1:nLons, 1:nLats,-1:nAlts+2, iEast_, iBlock), &
          ViscCoef(1:nLons, 1:nLats,0:nAlts+1), &
          Rho(1:nLons, 1:nLats,1:nAlts, iBlock), &
          Viscosity(1:nLons, 1:nLats,1:nAlts, iEast_))

     ! This is totally ad hoc.....

     dtr = dt/1800.0

     if (.not.Is1D) then
        do iLon = 1, nLons
           do iLat = 1, nLats
              Viscosity(iLon,iLat,1:nAlts,iUp_) = &
                   ((Velocity(iLon-1,iLat,1:nAlts,iUp_,iBlock) - &
                   Velocity(iLon,iLat,1:nAlts,iUp_,iBlock))/ &
                   dLonDist_GB(iLat, 1:nAlts, iBlock) + &
                   (Velocity(iLon+1,iLat,1:nAlts,iUp_,iBlock) - &
                   Velocity(iLon,iLat,1:nAlts,iUp_,iBlock))/ &
                   dLonDist_GB(iLat, 1:nAlts, iBlock)) / 2.0 + &
                   ((Velocity(iLon,iLat-1,1:nAlts,iUp_,iBlock) - &
                   Velocity(iLon,iLat,1:nAlts,iUp_,iBlock))/ &
                   dLatDist_GB(iLat, 1:nAlts, iBlock) + &
                   (Velocity(iLon,iLat-1,1:nAlts,iUp_,iBlock) - &
                   Velocity(iLon,iLat,1:nAlts,iUp_,iBlock))/ &
                   dLatDist_GB(iLat, 1:nAlts, iBlock)) / 2.0 * &
                   ViscCoef(iLon, iLat,1:nAlts) / &
                   Rho(iLon, iLat, 1:nAlts, iBlock)
              do iAlt = 1, nAlts
                 if (abs(Viscosity(iLon,iLat,iAlt,iUp_)) > &
                      dtr*abs(velocity(iLon,iLat,iAlt,iUp_,iBlock))) then
                    Viscosity(iLon,iLat,iAlt,iUp_) = &
                         sign(dtr*velocity(iLon,iLat,iAlt,iUp_,iBlock), &
                         Viscosity(iLon,iLat,iAlt,iUp_))
                 endif
              enddo
           enddo
        enddo

!        write(*,*) maxval(abs(Viscosity(:,:,:,iUp_)))

     endif
  else
     Viscosity = 0.0
  end if

  !\
  ! ---------------------------------------------------------------
  ! These terms are for Ion Densities
  ! ---------------------------------------------------------------
  !/

  if (iDebugLevel > 4) write(*,*) "=====> IonChemistry", iproc, UseIonChemistry

  if (iDebugLevel > 4) write(*,*) "=====> calc_ion_v", iproc

  call get_potential(iBlock)
  call calc_efield(iBlock)
  call aurora(iBlock)
  if (UseAuroralHeating) then
     AuroralHeating = AuroralHeatingRate(:,:,:,iBlock) / &
          TempUnit / cp(:,:,1:nAlts,iBlock)
  else
     AuroralHeating = 0.0
  endif

  call calc_ion_v(iBlock)

  ! These sources are for Earth only, so let's zero them out before
  ! we call planet specific sources...
  NOCooling = 0.0
  OCooling  = 0.0

  ! This routine is located in the planet.f90 file
  call calc_planet_sources(iBlock)

  call calc_chemistry(iBlock)

!  if (UseIonChemistry) then
!     call calc_ion_density(iBlock)
!  endif
!
!  if (UseNeutralChemistry) then
!     call calc_neutral_density(iBlock)
!  endif

  ChemicalHeatingRate(:,:,:) = &
       ChemicalHeatingRate(:,:,:) * Element_Charge / &
       TempUnit / cp(1:nLons,1:nLats,1:nAlts,iBlock)

end subroutine calc_GITM_sources
