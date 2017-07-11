!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine calc_GITM_sources(iBlock)

  use ModInputs
  use ModSources
  use ModGITM

  implicit none

  integer, intent(in) :: iBlock

  integer :: iAlt, iError, iDir, iLat, iLon, iSpecies
  integer :: iiAlt, iiLat, iiLon
  real :: tmp(nLons, nLats, nAlts)
  real :: tmp2(nLons, nLats, 0:nAlts+1)
  real :: tmp3(nLons, nLats, 0:nAlts+1)
  real :: RhoI(nLons, nLats, nAlts)
  real :: ScaleHeight(-1:nLons+2, -1:nLats+2, -1:nAlts+2)
  real :: Prandtl(nLons,nLats,0:nalts+1)

  call report("calc_GITM_sources",1)

  ! calc_rate is used to determine reaction rates, heating coefficients,
  ! ion-neutral collision frequency, lambdas for ion drag, etc.

  if (iDebugLevel > 4) write(*,*) "=====> going into calc_rates", iproc

  ChemicalHeatingRate = 0.0
  ChemicalHeatingRateIon = 0.0
  ChemicalHeatingRateEle = 0.0

  ChemicalHeatingSpecies = 0.0
  call calc_eddy_diffusion_coefficient(iBlock)  
  call calc_rates(iBlock)
  call calc_collisions(iBlock)

  RhoI = IDensityS(1:nLons,1:nLats,1:nAlts,ie_,iBlock) * &
       MeanIonMass(1:nLons,1:nLats,1:nAlts)
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

  ! This includes Radiative Cooling....
  RadCooling = 0.0
  call calc_planet_sources(iBlock)

  ! The auroral heating is specified below, after the aurora is described
  ! in get_potential


  !\
  ! Conduction ----------------------------------------------------
  !/

  if (iDebugLevel > 4) write(*,*) "=====> conduction", iproc

  if(UseConduction)then

     tmp2 = Rho(1:nLons, 1:nLats,0:nAlts+1, iBlock) * &
          cp(1:nLons, 1:nLats,0:nAlts+1, iBlock)
     
     Prandtl = 0.0

     call calc_conduction(iBlock, &
          Temperature(1:nLons, 1:nLats,-1:nAlts+2, iBlock) * &
          TempUnit(1:nLons, 1:nLats,-1:nAlts+2), &
          KappaTemp(1:nLons, 1:nLats,0:nAlts+1, iBlock), &
          tmp2, &
          MoleConduction)
     
      Conduction = MoleConduction/TempUnit(1:nLons, 1:nLats,1:nAlts) 

      ! JMB:  07/11/2017.
      !  This accounts for turbulent thermal conduction.
      !  Currently assume that the K_turbulent is proportional to
      !  rho*cp*Ke, consistent with the TIEGCM formulation.
      !  This can be found on the NCAR TIEGCM documentation site.
      if (UseTurbulentCond) then
            Prandtl = KappaEddyDiffusion(1:nLons, 1:nLats,0:nAlts+1, iBlock) * &
                                     Rho(1:nLons, 1:nLats,0:nAlts+1, iBlock) * &
                                      Cp(1:nLons, 1:nLats,0:nAlts+1, iBlock)

            call calc_conduction(iBlock, &
                 Temperature(1:nLons, 1:nLats,-1:nAlts+2, iBlock) * &
                 TempUnit(1:nLons, 1:nLats,-1:nAlts+2), &
                 Prandtl, &
                 tmp2, &
                 EddyCond)

            Conduction = Conduction + &
                 EddyCond/TempUnit(1:nLons, 1:nLats,1:nAlts) 
      endif  ! THE USETurbulentCond Check

  else
     Conduction = 0.0
  end if

  !\
  ! ---------------------------------------------------------------
  ! These terms are for Neutral Winds
  ! ---------------------------------------------------------------
  !/

  if (UseIonDrag) then

     tmp = Collisions(1:nLons,1:nLats,1:nAlts,iVIN_)*&
          RhoI/Rho(1:nLons,1:nLats,1:nAlts,iBlock)

     do iDir = 1, 3
        IonDrag(:,:,:,iDir) = tmp * &
             (IVelocity(1:nLons,1:nLats,1:nAlts,iDir,iBlock) - &
             Velocity(1:nLons,1:nLats,1:nAlts,iDir,iBlock))
     enddo

     ! F_iondrag = rho_i/rho * Vis * (Ui-Un)
     ! where Vis = Vin *(Ns/N)  

     do iSpecies = 1, nSpecies
        tmp = Collisions(1:nLons,1:nLats,1:nAlts,iVIN_)*&
             RhoI / &
             (Mass(iSpecies) * &
             NDensityS(1:nLons,1:nLats,1:nAlts,iSpecies,iBlock)) * &
             (NDensityS(1:nLons,1:nLats,1:nAlts,iSpecies,iBlock) / &
             NDensity(1:nLons,1:nLats,1:nAlts,iBlock))

        VerticalIonDrag(:,:,:,iSpecies) = tmp * &
             (IVelocity(1:nLons,1:nLats,1:nAlts,iUp_,iBlock) - &
             VerticalVelocity(1:nLons,1:nLats,1:nAlts,iSpecies,iBlock))

     enddo

  else

     IonDrag = 0.0
     VerticalIonDrag = 0.0

  endif

  !\
  ! Viscosity ----------------------------------------------------
  !/

  if (iDebugLevel > 4) write(*,*) "=====> conduction", iproc

  if(UseViscosity)then


     call calc_viscosity(iBlock)

     call calc_conduction(iBlock, &
          Velocity(1:nLons, 1:nLats,-1:nAlts+2, iNorth_, iBlock), &
          ViscCoef(1:nLons, 1:nLats,0:nAlts+1), &
          Rho(1:nLons, 1:nLats,0:nAlts+1, iBlock), &
          Viscosity(1:nLons, 1:nLats,1:nAlts, iNorth_))

     call calc_conduction(iBlock, &
          Velocity(1:nLons, 1:nLats,-1:nAlts+2, iEast_, iBlock), &
          ViscCoef(1:nLons, 1:nLats,0:nAlts+1), &
          Rho(1:nLons, 1:nLats,0:nAlts+1, iBlock), &
          Viscosity(1:nLons, 1:nLats,1:nAlts, iEast_))

     Viscosity(:,:,:,iUp_) = 0.0

     ! JMB:  07/11/2017
     ! --- Updated to include vertical viscosity for the vertical winds.
     !     provides stability to velocity calculations.
     do iSpecies = 1, nSpecies
        call calc_conduction(iBlock, &
             VerticalVelocity(1:nLons, 1:nLats,-1:nAlts+2, iSpecies, iBlock), &
             ViscCoef(1:nLons, 1:nLats,0:nAlts+1), &
             Rho(1:nLons, 1:nLats,0:nAlts+1, iBlock), &
             Viscosity(1:nLons, 1:nLats,1:nAlts, iUp_))

        VerticalVelocity(1:nLons,1:nLats,1:nAlts,iSpecies,iBlock) = &
        VerticalVelocity(1:nLons,1:nLats,1:nAlts,iSpecies,iBlock) + &
               Viscosity(1:nLons,1:nLats,1:nAlts, iUp_)

     enddo 

  else
     Viscosity = 0.0
  end if
!
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

  ! The Emissions array was never set. Should this be here or earlier ????
  Emissions(:,:,:,:,iBlock) = 0.0
  call calc_chemistry(iBlock)

  ChemicalHeatingRate(:,:,:) = &
       ChemicalHeatingRate(:,:,:) * Element_Charge / &
       TempUnit(1:nLons,1:nLats,1:nAlts) / cp(1:nLons,1:nLats,1:nAlts,iBlock)/&
       rho(1:nLons,1:nLats,1:nAlts,iBlock)
	   
  ChemicalHeatingRateIon(:,:,:) = &
       ChemicalHeatingRateIon(:,:,:) * Element_Charge
! / &
!       cp(1:nLons,1:nLats,1:nAlts,iBlock)/&
!       rho(1:nLons,1:nLats,1:nAlts,iBlock)

  ChemicalHeatingSpecies = ChemicalHeatingSpecies * Element_Charge
end subroutine calc_GITM_sources
