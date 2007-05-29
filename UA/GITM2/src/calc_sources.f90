
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
  real :: mmm(-1:nLons+2, -1:nLats+2, -1:nAlts+2)
  real :: mmr(-1:nLons+2, -1:nLats+2, -1:nAlts+2)
  real :: inside(-1:nLons+2, -1:nLats+2, -1:nAlts+2)
  real :: Omega(nLons, nLats, nAlts)
  real :: RhoI(nLons, nLats, nAlts)
  real :: ReducedMass(nLons, nLats, nAlts)
  real :: KappaEddyDiffusion(-1:nLons+2, -1:nLats+2, -1:nAlts+2)
  real :: ScaleHeight(-1:nLons+2, -1:nLats+2, -1:nAlts+2)
  real :: TmpGradient(nLons, nLats, nAlts,3)

  real :: NSoverN(-1:nLons+2, -1:nLats+2, -1:nAlts+2)
  real :: NSoverNGradient(1:nLons, 1:nLats, 1:nAlts, 3)
  real:: LogNum(-1:nLons+2, -1:nLats+2, -1:nAlts+2)

  real :: oldtemp(-1:nLons+2, -1:nLats+2, -1:nAlts+2)

  real :: dtsub, dttotal, dtsubmax, dtr

  real :: diffusion_velocity(nLons, nLats, nAlts,nspecies)

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
          Rho(1:nLons,1:nLats,1:nAlts,iBlock) 

     ! After reading the Thayer paper, this term needs to be taken
     ! out, since it is approximately 1/2, which is the contribution
     ! from ions heating the neutrals through conduction.
     ! * &
     !     (MeanIonMass(1:nLons,1:nLats,1:nAlts)/AMU) /  &
     !     (MeanIonMass(1:nLons,1:nLats,1:nAlts)/AMU + &
     !     MeanMajorMass(1:nLons,1:nLats,1:nAlts)/AMU)

     JouleHeating = 0.0

     do iDir = 1, 3

        JouleHeating(:,:,:) = JouleHeating(:,:,:) + tmp2 * &
             (IVelocity(1:nLons,1:nLats,1:nAlts,iDir,iBlock) - &
             Velocity(1:nLons,1:nLats,1:nAlts,iDir,iBlock))**2

     enddo

     JouleHeating = JouleHeating / &
          TempUnit(1:nLons,1:nLats,1:nAlts) / &
          cp(:,:,1:nAlts,iBlock)

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
          Temperature(1:nLons, 1:nLats,-1:nAlts+2, iBlock) * &
          TempUnit(1:nLons, 1:nLats,-1:nAlts+2), &
          KappaTemp(1:nLons, 1:nLats,0:nAlts+1, iBlock), &
          tmp2, &
          Conduction)

     Conduction = Conduction/TempUnit(1:nLons, 1:nLats,1:nAlts)
     
 else
     Conduction = 0.0
  end if

  !\
  ! Diffusion ----------------------------------------------------
  !/

  if (iDebugLevel > 4) write(*,*) "=====> diffusion ", iproc

!!---------------------------------------------
!! Yue's method
!!---------------------------------------------

  if (UseDiffusion)then
     KappaEddyDiffusion=0.
     do iAlt = -1, nAlts+2

        do iLat = 1, nLats
           do iLon = 1, nLons

              if (pressure(iLon,iLat,iAlt,iBlock) >EddyDiffusionPressure0) then
                 KappaEddyDiffusion(iLon,iLat,iAlt) = EddyDiffusionCoef

              else if (pressure(iLon,iLat,iAlt,iBlock) > &
                   EddyDiffusionPressure1) then

                 KappaEddyDiffusion(iLon,iLat,iAlt) = EddyDiffusionCoef * &
                      (pressure(iLon,iLat,iAlt,iBlock) - &
                      EddyDiffusionPressure1)/&
                      (EddyDiffusionPressure0 - EddyDiffusionPressure1)

              endif
           enddo
        enddo
     enddo

     LogNum(:,:,:) = log(NDensity(:,:,:,iblock))
     do iSpecies = 1, nSpecies
        do iAlt = 1, nAlts
           do iLat = 1, nLats
              do iLon = 1, nLons
                 Diffusion_velocity(iLon,iLat,iAlt,iSpecies) =      &
                      -KappaEddyDiffusion(iLon,iLat,iAlt)* &
                      (logNS(iLon,iLat,iAlt+1,iSpecies,iBlock)- &
                      logNS(iLon,iLat,iAlt-1,iSpecies,iBlock)- &
                      logNum(iLon,iLat,iAlt+1)+ &
                      logNum(iLon,iLat,iAlt-1))/(2*dalt(iAlt))
              enddo
           enddo
        enddo
     enddo

     do iSpecies = 1, nSpecies
        do iAlt = 2, nAlts-1
           do iLat = 1, nLats
              do iLon = 1, nLons

                 diffusion(ilon,ilat,iAlt,ispecies) =      &
                      (NdensityS(ilon,ilat,iAlt+1,ispecies,iBlock)  *&
                      diffusion_velocity(ilon,ilat,iAlt+1,ispecies)  -&
                      NdensityS(ilon,ilat,iAlt-1,ispecies,iBlock)*&
                      diffusion_velocity(ilon,ilat,iAlt-1,ispecies))/&
                      (2*dalt(iAlt))

              enddo
           enddo
        enddo
     enddo


     do iSpecies = 1, nSpecies
           do iLat = 1, nLats
              do iLon = 1, nLons

                 diffusion(ilon,ilat,1,ispecies) =      &
                      (NdensityS(ilon,ilat,2,ispecies,iBlock)  *&
                      diffusion_velocity(ilon,ilat,2,ispecies)  -&
                      NdensityS(ilon,ilat,1,ispecies,iBlock)*&
                      diffusion_velocity(ilon,ilat,1,ispecies))/&
                      (dalt(1))

                 diffusion(ilon,ilat,nAlts,ispecies) =      &
                      (NdensityS(ilon,ilat,nAlts,ispecies,iBlock)  *&
                      diffusion_velocity(ilon,ilat,nAlts,ispecies)  -&
                      NdensityS(ilon,ilat,nAlts-1,ispecies,iBlock)*&
                      diffusion_velocity(ilon,ilat,nAlts-1,ispecies))/&
                      (dalt(nAlts))
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

     ! F_iondrag = rho_i/rho * Vis * (Ui-Un)
     ! where Vis = Vin *(Ns/N)  

     do iSpecies = 1, nSpecies
        tmp2 = Collisions(1:nLons,1:nLats,1:nAlts,iVIN_)*&
             RhoI / &
             (Mass(iSpecies) * &
             NDensityS(1:nLons,1:nLats,1:nAlts,iSpecies,iBlock)) * &
            (NDensityS(1:nLons,1:nLats,1:nAlts,iSpecies,iBlock) / &
             NDensity(1:nLons,1:nLats,1:nAlts,iBlock))
             
        VerticalIonDrag(:,:,:,iSpecies) = tmp2 * &
             (IVelocity(1:nLons,1:nLats,1:nAlts,iUp_,iBlock) - &
             VerticalVelocity(1:nLons,1:nLats,1:nAlts,iSpecies,iBlock))

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
                   dLonDist_GB(iLon, iLat, 1:nAlts, iBlock) + &
                   (Velocity(iLon+1,iLat,1:nAlts,iUp_,iBlock) - &
                   Velocity(iLon,iLat,1:nAlts,iUp_,iBlock))/ &
                   dLonDist_GB(iLon,iLat, 1:nAlts, iBlock)) / 2.0 + &
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
          TempUnit(1:nLons,1:nLats,1:nAlts) / cp(:,:,1:nAlts,iBlock) / &
          rho(1:nLons,1:nLats,1:nAlts, iBlock)
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

  ChemicalHeatingRate(:,:,:) = &
       ChemicalHeatingRate(:,:,:) * Element_Charge / &
       TempUnit(1:nLons,1:nLats,1:nAlts) / cp(1:nLons,1:nLats,1:nAlts,iBlock)/&
       rho(1:nLons,1:nLats,1:nAlts,iBlock)

end subroutine calc_GITM_sources
