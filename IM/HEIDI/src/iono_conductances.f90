!*************************************************************************
!
! CONDUCTANCE Calculation Routines
!
!*************************************************************************


!-------------------------------------------------------------------------
! FACs_to_fluxes
!
!
!
!-------------------------------------------------------------------------

subroutine FACs_to_fluxes(iModel)

  !\
  ! The goal here is to convert the ionospheric FAC pattern into a 
  ! particle precipitation pattern, which can then be turned into
  ! a conductance pattern.
  !/

  use ModIonosphere
  implicit none

  integer, intent(in) :: iModel

  real    :: aurora_peak_lat, eflux_peak_lon, ave_e_peak_lon
  real    :: aurora_peak_width, ave_e_peak
  real    :: distance, mult_fac, ef, efp, mfc
  real    :: hal_a0, hal_a1, hal_a2, ped_a0, ped_a1, ped_a2, hall, ped
  real    :: dlat, dmlt, y1, y2, x1, x2
  real    :: day_colat, dusk_colat, midnight_colat, dawn_colat
  real    :: day_fac, dusk_fac, midnight_fac, dawn_fac
  real    :: noon_mid, dusk_dawn, oval_shift, oval_tilt
  real    :: mean_colat, dev_colat, sum, night_width
  real    :: PolarCapHallConductance, PolarCap_AveE, PolarCap_EFlux
  integer :: jlat, imlt
  integer :: i,j, n, nloc
  real, dimension(1:IONO_nPsi) :: Strength_of_Oval,               &
                                  Loc_of_Oval, Width_of_Oval
  logical :: polarcap

  Hall_to_Ped_Ratio = 1.5

  if (PolarCapPedConductance > 0.0) then
     PolarCapHallConductance = Hall_to_Ped_Ratio * PolarCapPedConductance
     PolarCap_AveE = (Hall_to_Ped_Ratio/0.45)**(1.0/0.85)
     PolarCap_EFlux = ((PolarCapPedConductance*(16.0 + PolarCap_AveE**2) / &
                        (40.0*PolarCap_AveE))**2)/1000.0
  else
     PolarCap_AveE  = IONO_Min_Ave_E
     PolarCap_EFlux = IONO_Min_EFlux
  endif

  if (abs(iModel).eq.3) then

     write(6,*) "=> Using Super Smooth Oval."

!!! mult_fac is a scaling factor to convert Eflux into conductance
!!! Aaron empirically determined these coefficients based on MHD results
!!! His numbers were mult_fac = 0.25*0.33e7*2.0e-3
!!! We should increase mult_fac to account for R1 currents, perhaps 10x
!!! The first coefficient is my scaling factor (MWL, 11/22/2002)

     mult_fac   = 5.0 * 0.25*0.33e7*2.0e-3
     ave_e_peak = max(3.0, PolarCap_AveE)                        ! keV

!!! oval_shift is to offset the conductance peak from the FAC peak
!!! positive is equatorward, negative is poleward
!!! we had it set to 2.5*IONO_PI/180., meaning a slight shift equatorward
!!! Aaron uses +10 for the MHD runs to put it in between R1 and R2 FACs
!!! we should try zero as well as a poleward shift to do this too.
!!! The first coefficient is the offset in degrees (MWL, 11/22/2002)
!!! Changed to be a tilted, offset auroral oval (MWL, 10/27/2003)
!!! >> oval_shift is now the centroid of the oval offset
!!! >> oval_tilt is now the amplitude (in degrees) of the tilt
!!! >> The tilt is dawn-dusk: equatorward and poleward, respectively,
!!! >> for a negative oval_tilt.

     oval_shift = -5.00 * IONO_PI/180.0
     oval_tilt  = -0.00 * IONO_PI/180.0

     !\
     ! Do North first -----------------------------------------------------
     !/

     call Determine_Oval_Characteristics(IONO_NORTH_JR, &
                                         IONO_NORTH_Theta, &
                                         IONO_NORTH_Psi, &
                                         Loc_of_Oval, &
                                         Width_of_Oval, &
                                         Strength_of_Oval)

     if (iModel < 0) then
        do j = 1, IONO_nPsi
           if (Width_of_Oval(j) < 5.0*IONO_PI/180.0) then
              Width_of_Oval(j) = 5.0*IONO_PI/180.0
           endif
        enddo
     endif

     ave_e_peak_lon = 0.0

     do i = 1, IONO_nTheta
        do j = 1, IONO_nPsi

           ! Energy Flux

           distance = IONO_NORTH_Theta(i,j)-(Loc_of_Oval(j)+oval_shift &
                +oval_tilt*sin(IONO_NORTH_Psi(i,j)))
           IONO_NORTH_EFlux(i,j) =                                            &
                Strength_of_Oval(j) * &
                mult_fac * &
                exp(-1.0*(distance/Width_of_Oval(j))**2)

           if (distance < 0.0) then
              IONO_NORTH_EFlux(i,j) = &
                   max(IONO_NORTH_EFlux(i,j), PolarCap_EFlux)
              polarcap = .true.
           else
              IONO_NORTH_EFlux(i,j) = &
                   max(IONO_NORTH_EFlux(i,j), IONO_Min_EFlux)
              polarcap = .false.
           endif

           ! Average Energy

           IONO_NORTH_Ave_E(i,j) =                                            &
                ave_e_peak*exp(-1.0*(distance/Width_of_Oval(j))**2)

           if (PolarCap_AveE == IONO_Min_Ave_E) then
              distance = 0.25 +                                          &
                   0.75*(                                                &
                   sin(IONO_NORTH_Psi(i,j)-(ave_e_peak_lon+IONO_PI/2.0)) &
                   + 1.0)/2.0
           else
              distance = 1.0
           endif

           IONO_NORTH_Ave_E(i,j) = IONO_NORTH_Ave_E(i,j)*distance

           if (polarcap) then
              IONO_NORTH_Ave_E(i,j) = &
                   max(IONO_NORTH_Ave_E(i,j), PolarCap_AveE)
           else
              IONO_NORTH_Ave_E(i,j) = &
                   max(IONO_NORTH_Ave_E(i,j), IONO_Min_Ave_E)
           endif

        enddo
     enddo

     !\
     ! Do South next --------------------------------------------------------
     !/

     call Determine_Oval_Characteristics(IONO_SOUTH_JR, &
                                         IONO_SOUTH_Theta, &
                                         IONO_SOUTH_Psi, &
                                         Loc_of_Oval, &
                                         Width_of_Oval, &
                                         Strength_of_Oval)

     if (iModel < 0) then
        do j = 1, IONO_nPsi
           if (Width_of_Oval(j) < 5.0*IONO_PI/180.0) then
              Width_of_Oval(j) = 5.0*IONO_PI/180.0
           endif
        enddo
     endif

     Loc_of_Oval = IONO_PI - Loc_of_Oval

     ave_e_peak_lon = 0.0

     do i = 1, IONO_nTheta
        do j = 1, IONO_nPsi

           ! Energy Flux

           distance = IONO_SOUTH_Theta(i,j)-(Loc_of_Oval(j)-oval_shift &
                -oval_tilt*sin(IONO_SOUTH_Psi(i,j)))
           IONO_SOUTH_EFlux(i,j) =                                            &
                Strength_of_Oval(j) * &
                mult_fac * &
                exp(-1.0*(distance/Width_of_Oval(j))**2)

           if (distance > 0.0) then
              IONO_SOUTH_EFlux(i,j) = &
                   max(IONO_SOUTH_EFlux(i,j), PolarCap_EFlux)
              polarcap = .true.
           else
              IONO_SOUTH_EFlux(i,j) = &
                   max(IONO_SOUTH_EFlux(i,j), IONO_Min_EFlux)
              polarcap = .false.
           endif

           ! Average Energy

           IONO_SOUTH_Ave_E(i,j) =                                            &
                ave_e_peak*exp(-1.0*(distance/Width_of_Oval(j))**2)

           if (PolarCap_AveE == IONO_Min_Ave_E) then
              distance = 0.25 +                                          &
                   0.75*(                                                &
                   sin(IONO_SOUTH_Psi(i,j)-(ave_e_peak_lon+IONO_PI/2.0)) &
                   + 1.0)/2.0
           else
              distance = 1.0
           endif

           IONO_SOUTH_Ave_E(i,j) = IONO_SOUTH_Ave_E(i,j)*distance

           if (polarcap) then
              IONO_SOUTH_Ave_E(i,j) = &
                   max(IONO_SOUTH_Ave_E(i,j), PolarCap_AveE)
           else
              IONO_SOUTH_Ave_E(i,j) = &
                   max(IONO_SOUTH_Ave_E(i,j), IONO_Min_Ave_E)
           endif
        enddo
     enddo

   endif

   if (abs(iModel).eq.4) then

      write(6,*) "This Ionospheric Model is Under Construction"
      stop

   endif

   if (abs(iModel).eq.5) then 

      write(6,*) "=> Using AMIE derived relationship between FACs and Precip."

      open(unit=23,file="fac_to_cond.dat",status="old")
      read(23,*) i_cond_nmlts, i_cond_nlats

      do i = 1, i_cond_nmlts
         do j = 1, i_cond_nlats

            read(23,"(14F6.2)") cond_mlts(i),cond_lats(j), &
                 hal_a0_up(i,j),ped_a0_up(i,j),      &
                 hal_a0_do(i,j),ped_a0_do(i,j),      &
                 hal_a1_up(i,j),ped_a1_up(i,j),      &
                 hal_a1_do(i,j),ped_a1_do(i,j),      &
                 hal_a2_up(i,j),ped_a2_up(i,j),      &
                 hal_a2_do(i,j),ped_a2_do(i,j)

         enddo
      enddo

      close(23)

      do j = 1, i_cond_nlats

         hal_a0_up(i_cond_nmlts+1,j) = hal_a0_up(1,j)
         ped_a0_up(i_cond_nmlts+1,j) = ped_a0_up(1,j)
         hal_a0_do(i_cond_nmlts+1,j) = hal_a0_do(1,j)
         ped_a0_do(i_cond_nmlts+1,j) = ped_a0_do(1,j)
         hal_a1_up(i_cond_nmlts+1,j) = hal_a1_up(1,j)
         ped_a1_up(i_cond_nmlts+1,j) = ped_a1_up(1,j)
         hal_a1_do(i_cond_nmlts+1,j) = hal_a1_do(1,j)
         ped_a1_do(i_cond_nmlts+1,j) = ped_a1_do(1,j)
         hal_a2_up(i_cond_nmlts+1,j) = hal_a2_up(1,j)
         ped_a2_up(i_cond_nmlts+1,j) = ped_a2_up(1,j)
         hal_a2_do(i_cond_nmlts+1,j) = hal_a2_do(1,j)
         ped_a2_do(i_cond_nmlts+1,j) = ped_a2_do(1,j)

      enddo

      dlat = (cond_lats(1) - cond_lats(2))*IONO_PI/180.0
      dmlt = (cond_mlts(2) - cond_mlts(1))*IONO_PI/12.0

!  Do North First

      call Determine_Oval_Characteristics(IONO_NORTH_JR, &
           IONO_NORTH_Theta, &
           IONO_NORTH_Psi, &
           Loc_of_Oval, &
           Width_of_Oval, &
           Strength_of_Oval)

      Width_of_Oval = Width_of_Oval * 2.0

      do j = 1, IONO_nPsi
         do i = 1, IONO_nTheta

            y1 = IONO_NORTH_Theta(i,j)/dlat + 1.0
            if (y1 > i_cond_nlats) then
               jlat = i_cond_nlats
               y1   = 1.0
            else
               jlat = y1
               y1   = 1.0 - (y1 - jlat)
            endif
            y2 = 1.0 - y1

            x1 = mod((IONO_NORTH_Psi(i,j) + IONO_PI),2.0*IONO_PI)/dmlt + 1.0
            imlt = x1
            x1   = 1.0 - (x1 - imlt)
            x2   = 1.0 - x1

            if (iono_north_jr(i,j) > 0) then

               hal_a0 = x1*y1*hal_a0_up(imlt  ,jlat  ) + &
                        x2*y1*hal_a0_up(imlt+1,jlat  ) + &
                        x1*y2*hal_a0_up(imlt  ,jlat+1) + &
                        x2*y2*hal_a0_up(imlt+1,jlat+1)

               hal_a1 = x1*y1*hal_a1_up(imlt  ,jlat  ) + &
                        x2*y1*hal_a1_up(imlt+1,jlat  ) + &
                        x1*y2*hal_a1_up(imlt  ,jlat+1) + &
                        x2*y2*hal_a1_up(imlt+1,jlat+1)

               hal_a2 = x1*y1*hal_a2_up(imlt  ,jlat  ) + &
                        x2*y1*hal_a2_up(imlt+1,jlat  ) + &
                        x1*y2*hal_a2_up(imlt  ,jlat+1) + &
                        x2*y2*hal_a2_up(imlt+1,jlat+1)

               ped_a0 = x1*y1*ped_a0_up(imlt  ,jlat  ) + &
                        x2*y1*ped_a0_up(imlt+1,jlat  ) + &
                        x1*y2*ped_a0_up(imlt  ,jlat+1) + &
                        x2*y2*ped_a0_up(imlt+1,jlat+1)

               ped_a1 = x1*y1*ped_a1_up(imlt  ,jlat  ) + &
                        x2*y1*ped_a1_up(imlt+1,jlat  ) + &
                        x1*y2*ped_a1_up(imlt  ,jlat+1) + &
                        x2*y2*ped_a1_up(imlt+1,jlat+1)

               ped_a2 = x1*y1*ped_a2_up(imlt  ,jlat  ) + &
                        x2*y1*ped_a2_up(imlt+1,jlat  ) + &
                        x1*y2*ped_a2_up(imlt  ,jlat+1) + &
                        x2*y2*ped_a2_up(imlt+1,jlat+1)

            else

               hal_a0 = x1*y1*hal_a0_do(imlt  ,jlat  ) + &
                        x2*y1*hal_a0_do(imlt+1,jlat  ) + &
                        x1*y2*hal_a0_do(imlt  ,jlat+1) + &
                        x2*y2*hal_a0_do(imlt+1,jlat+1)

               hal_a1 = x1*y1*hal_a1_do(imlt  ,jlat  ) + &
                        x2*y1*hal_a1_do(imlt+1,jlat  ) + &
                        x1*y2*hal_a1_do(imlt  ,jlat+1) + &
                        x2*y2*hal_a1_do(imlt+1,jlat+1)

               hal_a2 = x1*y1*hal_a2_do(imlt  ,jlat  ) + &
                        x2*y1*hal_a2_do(imlt+1,jlat  ) + &
                        x1*y2*hal_a2_do(imlt  ,jlat+1) + &
                        x2*y2*hal_a2_do(imlt+1,jlat+1)

               ped_a0 = x1*y1*ped_a0_do(imlt  ,jlat  ) + &
                        x2*y1*ped_a0_do(imlt+1,jlat  ) + &
                        x1*y2*ped_a0_do(imlt  ,jlat+1) + &
                        x2*y2*ped_a0_do(imlt+1,jlat+1)

               ped_a1 = x1*y1*ped_a1_do(imlt  ,jlat  ) + &
                        x2*y1*ped_a1_do(imlt+1,jlat  ) + &
                        x1*y2*ped_a1_do(imlt  ,jlat+1) + &
                        x2*y2*ped_a1_do(imlt+1,jlat+1)

               ped_a2 = x1*y1*ped_a2_do(imlt  ,jlat  ) + &
                        x2*y1*ped_a2_do(imlt+1,jlat  ) + &
                        x1*y2*ped_a2_do(imlt  ,jlat+1) + &
                        x2*y2*ped_a2_do(imlt+1,jlat+1)

            endif

            distance = IONO_NORTH_Theta(i,j) - Loc_of_Oval(j)

            !
            ! We want minimal conductance lower than the oval
            !

            if (distance > 0.0) then
               hal_a0 = hal_a0 * exp(-1.0*(distance/(Width_of_Oval(j)*2))**2)
               ped_a0 = ped_a0 * exp(-1.0*(distance/(Width_of_Oval(j)*2))**2)
               polarcap = .false.
            else
               polarcap = .true.
            endif

            !
            ! A sort of correction on the fit
            !

            hal_a1 = hal_a0 + (hal_a1 - hal_a0)*  &
                              exp(-1.0*(distance/Width_of_Oval(j))**2)
            ped_a1 = ped_a0 + (ped_a1 - ped_a0)*  &
                              exp(-1.0*(distance/Width_of_Oval(j))**2)

            hall=hal_a0-hal_a1*exp(-abs(iono_north_jr(i,j)*1.0e9)*hal_a2**2)
            ped =ped_a0-ped_a1*exp(-abs(iono_north_jr(i,j)*1.0e9)*ped_a2**2)

            if ((hall.gt.0).and.(ped.gt.0)) then

               IONO_NORTH_Ave_E(i,j)  = ((hall/ped)/0.45)**(1.0/0.85)
               IONO_NORTH_EFlux(i,j) = (ped*(16.0+IONO_NORTH_Ave_E(i,j)**2)/ &
                    (40.0*IONO_NORTH_Ave_E(i,j)))**2/1000.0

            else

               IONO_NORTH_Ave_E(i,j) = IONO_Min_Ave_E
               IONO_NORTH_EFlux(i,j) = IONO_Min_EFlux

            endif

            if ((PolarCap_AveE > 0.0).and.(polarcap)) then
               IONO_NORTH_Ave_E(i,j) = max(IONO_NORTH_Ave_E(i,j), &
                                           PolarCap_AveE)
               IONO_NORTH_EFlux(i,j) = max(IONO_NORTH_EFlux(i,j), &
                                           PolarCap_EFlux)
            endif

         enddo
      enddo

!  Do South Next

      call Determine_Oval_Characteristics(IONO_SOUTH_JR, &
           IONO_SOUTH_Theta, &
           IONO_SOUTH_Psi, &
           Loc_of_Oval, &
           Width_of_Oval, &
           Strength_of_Oval)

      Loc_of_Oval = IONO_PI - Loc_of_Oval

      do j = 1, IONO_nPsi
         do i = 1, IONO_nTheta


            y1 = (IONO_PI-IONO_SOUTH_Theta(i,j))/dlat + 1.0
            if (y1 > i_cond_nlats) then
               jlat = i_cond_nlats
               y1   = 1.0
            else
               jlat = y1
               y1   = 1.0 - (y1 - jlat)
            endif
            y2 = 1.0 - y1

            x1 = mod((IONO_SOUTH_Psi(i,j) + IONO_PI),2.0*IONO_PI)/dmlt + 1.0
            imlt = x1
            x1   = 1.0 - (x1 - imlt)
            x2   = 1.0 - x1

            if (iono_south_jr(i,j) > 0) then

               hal_a0 = x1*y1*hal_a0_up(imlt  ,jlat  ) + &
                        x2*y1*hal_a0_up(imlt+1,jlat  ) + &
                        x1*y2*hal_a0_up(imlt  ,jlat+1) + &
                        x2*y2*hal_a0_up(imlt+1,jlat+1)

               hal_a1 = x1*y1*hal_a1_up(imlt  ,jlat  ) + &
                        x2*y1*hal_a1_up(imlt+1,jlat  ) + &
                        x1*y2*hal_a1_up(imlt  ,jlat+1) + &
                        x2*y2*hal_a1_up(imlt+1,jlat+1)

               hal_a2 = x1*y1*hal_a2_up(imlt  ,jlat  ) + &
                        x2*y1*hal_a2_up(imlt+1,jlat  ) + &
                        x1*y2*hal_a2_up(imlt  ,jlat+1) + &
                        x2*y2*hal_a2_up(imlt+1,jlat+1)

               ped_a0 = x1*y1*ped_a0_up(imlt  ,jlat  ) + &
                        x2*y1*ped_a0_up(imlt+1,jlat  ) + &
                        x1*y2*ped_a0_up(imlt  ,jlat+1) + &
                        x2*y2*ped_a0_up(imlt+1,jlat+1)

               ped_a1 = x1*y1*ped_a1_up(imlt  ,jlat  ) + &
                        x2*y1*ped_a1_up(imlt+1,jlat  ) + &
                        x1*y2*ped_a1_up(imlt  ,jlat+1) + &
                        x2*y2*ped_a1_up(imlt+1,jlat+1)

               ped_a2 = x1*y1*ped_a2_up(imlt  ,jlat  ) + &
                        x2*y1*ped_a2_up(imlt+1,jlat  ) + &
                        x1*y2*ped_a2_up(imlt  ,jlat+1) + &
                        x2*y2*ped_a2_up(imlt+1,jlat+1)

            else

               hal_a0 = x1*y1*hal_a0_do(imlt  ,jlat  ) + &
                        x2*y1*hal_a0_do(imlt+1,jlat  ) + &
                        x1*y2*hal_a0_do(imlt  ,jlat+1) + &
                        x2*y2*hal_a0_do(imlt+1,jlat+1)

               hal_a1 = x1*y1*hal_a1_do(imlt  ,jlat  ) + &
                        x2*y1*hal_a1_do(imlt+1,jlat  ) + &
                        x1*y2*hal_a1_do(imlt  ,jlat+1) + &
                        x2*y2*hal_a1_do(imlt+1,jlat+1)

               hal_a2 = x1*y1*hal_a2_do(imlt  ,jlat  ) + &
                        x2*y1*hal_a2_do(imlt+1,jlat  ) + &
                        x1*y2*hal_a2_do(imlt  ,jlat+1) + &
                        x2*y2*hal_a2_do(imlt+1,jlat+1)

               ped_a0 = x1*y1*ped_a0_do(imlt  ,jlat  ) + &
                        x2*y1*ped_a0_do(imlt+1,jlat  ) + &
                        x1*y2*ped_a0_do(imlt  ,jlat+1) + &
                        x2*y2*ped_a0_do(imlt+1,jlat+1)

               ped_a1 = x1*y1*ped_a1_do(imlt  ,jlat  ) + &
                        x2*y1*ped_a1_do(imlt+1,jlat  ) + &
                        x1*y2*ped_a1_do(imlt  ,jlat+1) + &
                        x2*y2*ped_a1_do(imlt+1,jlat+1)

               ped_a2 = x1*y1*ped_a2_do(imlt  ,jlat  ) + &
                        x2*y1*ped_a2_do(imlt+1,jlat  ) + &
                        x1*y2*ped_a2_do(imlt  ,jlat+1) + &
                        x2*y2*ped_a2_do(imlt+1,jlat+1)

            endif

            distance = IONO_SOUTH_Theta(i,j)-Loc_of_Oval(j)

            if (distance < 0.0) then
               hal_a0 = hal_a0 * exp(-1.0*(distance/(Width_of_Oval(j)*2))**2)
               ped_a0 = ped_a0 * exp(-1.0*(distance/(Width_of_Oval(j)*2))**2)
               polarcap = .false.
            else
               polarcap = .true.
            endif

            !
            ! A sort of correction on the fit
            !

            hal_a1 = hal_a0 + (hal_a1 - hal_a0)*  &
                              exp(-1.0*(distance/Width_of_Oval(j))**2)
            ped_a1 = ped_a0 + (ped_a1 - ped_a0)*  &
                              exp(-1.0*(distance/Width_of_Oval(j))**2)

            hall=hal_a0-hal_a1*exp(-abs(iono_south_jr(i,j)*1.0e9)*hal_a2**2)
            ped =ped_a0-ped_a1*exp(-abs(iono_south_jr(i,j)*1.0e9)*ped_a2**2)

            if ((hall.gt.0).and.(ped.gt.0)) then

               IONO_SOUTH_Ave_E(i,j)  = ((hall/ped)/0.45)**(1.0/0.85)
               IONO_SOUTH_EFlux(i,j) = (ped*(16.0+IONO_SOUTH_Ave_E(i,j)**2)/ &
                    (40.0*IONO_SOUTH_Ave_E(i,j)))**2/1000.0

            else

               IONO_SOUTH_Ave_E(i,j) = IONO_Min_Ave_E
               IONO_SOUTH_EFlux(i,j) = IONO_Min_EFlux

            endif

            if ((PolarCap_AveE > 0.0).and.(polarcap)) then
               IONO_SOUTH_Ave_E(i,j) = max(IONO_SOUTH_Ave_E(i,j), &
                                           PolarCap_AveE)
               IONO_SOUTH_EFlux(i,j) = max(IONO_SOUTH_EFlux(i,j),        &
                                           PolarCap_EFlux)
            endif

         enddo
      enddo

   endif

end subroutine FACs_to_fluxes


!-------------------------------------------------------------------------
! ionosphere_conductance
!
!
!
!-------------------------------------------------------------------------

subroutine ionosphere_conductance(Sigma0, SigmaH, SigmaP,               &
                                  SigmaThTh, SigmaThPs, SigmaPsPs,      &
                                  dSigmaThTh_dTheta, dSigmaThPs_dTheta, &
                                  dSigmaPsPs_dTheta,                    &
                                  dSigmaThTh_dPsi, dSigmaThPs_dPsi,     &
                                  dSigmaPsPs_dPsi,                      &
                                  Eflux, Ave_E,                         &
                                  Theta, Psi, nTheta, nPsi,             &
                                  dTheta, dPsi,                         &
                                  iModel, f107)

  !\
  ! This subroutine computes the height-integrated field-aligned and
  ! Hall and Pedersen conductances for the ionosphere at each
  ! location of the discretized solution domain.  The gradients of
  ! these quantities are also computed.
  !/

  use ModIonosphere
  use ModMain
  use ModPhysics
  implicit none

  integer :: nTheta,nPsi
  integer, intent(in) :: iModel
  real, intent(in)    :: f107
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::  &
                Sigma0, SigmaH, SigmaP, &
                SigmaThTh, SigmaThPs, SigmaPsPs, &
                dSigmaThTh_dTheta, dSigmaThPs_dTheta, dSigmaPsPs_dTheta, &
                dSigmaThTh_dPsi, dSigmaThPs_dPsi, dSigmaPsPs_dPsi, &
                Eflux, Ave_E,                                            &
                Theta, Psi, sin_clat, cos_clat, sin_lon, cos_lon, cy,    &
                oSigmaH, oSigmaP, conv_SigmaH, conv_SigmaP

  real, dimension(1:IONO_NTheta) :: dTheta
  real, dimension(1:IONO_NPsi)   :: dPsi

  integer :: i,j
  real :: f107p53, f107p49
  real :: sn, cs, sn2, cs2, cs3, cs4, C
  real :: cos_SZA, SigmaH_EUV, SigmaP_EUV, SigmaH_SCAT, SigmaP_SCAT, &
          SigmaH_STAR, SigmaP_STAR
  real :: SigmaH_Particles, SigmaP_Particles, tau
  logical :: north, old
  
  real    :: DoY, max_tilt, time_delay
  integer :: jday
  real    :: sin_doy, cos_doy, sin_tilt, cos_tilt

  if (theta(1,1) < 2.0*IONO_Theta_0) then
     north = .true.
  else
     north = .false.
  endif

  if (iModel > 1) then

     Max_Tilt = 23.5*IONO_PI/180.0

     DoY = (jday(Start_Time_Array(1),Start_Time_Array(2),Start_Time_Array(3)) &
          -80.0)*2.*IONO_PI/365.

!     if (north) write(6,*) "=> Day of Year (normalized to March 21) : ",DoY

     sin_doy = sin(DoY)
     cos_doy = cos(DoY)

     sin_tilt = sin(Max_Tilt)
     cos_tilt = cos(Max_Tilt)

     sin_clat = sin(Theta)
     cos_clat = cos(Theta)

     sin_lon  = sin(Psi - DoY + IONO_PI/2.0)
     cos_lon  = cos(Psi - DoY + IONO_PI/2.0)

     cy = cos_lon  * sin_clat * sin_doy * cos_tilt +                          &
          sin_lon  * sin_clat * cos_doy +                                     &
          cos_clat * sin_doy  * sin_tilt

! We are going to need F10.7 ^ 0.53 and F10.7 ^ 0.49 a lot,
! So, let's just store them straight away:

     f107p53 = f107**0.53
     f107p49 = f107**0.49

  endif

  if (iModel.eq.0) then

     if (north) write(6,*) "=> Using Constant Pedersen Conductance (no HALL)"

     do j = 1, nPsi
        do i = 1, nTheta

           Sigma0(i,j) = 1000.00
           SigmaH(i,j) = 0.00
           SigmaP(i,j) = 10.00

        enddo
     enddo

  endif

  if (iModel.eq.1) then

     if (north) write(6,*) "=> Using Constant Pedersen and Hall Conductance"

     do j = 1, nPsi
        do i = 1, nTheta

           Sigma0(i,j) = 1000.00
           SigmaH(i,j) = 10.00
           SigmaP(i,j) = 5.00

        enddo
     enddo

  endif

  if (iModel.eq.2) then

     if (north) write(6,*) "=> Using Solar EUV Only"

     do j = 1, nPsi
        do i = 1, nTheta

!!$        !\
!!$        ! Rasmussen and Schunk model B (JGR, Vol. 92, pp. 4491-4504, 1987).
!!$        !/

           Sigma0(i,j) = 1000.00

           cos_SZA = cy(i,j)

           if (cos_SZA > 0) then
              SigmaH_EUV = f107p53*(0.81*cos_SZA + 0.54*sqrt(cos_SZA))
              SigmaP_EUV = f107p49*(0.34*cos_SZA + 0.93*sqrt(cos_SZA))
              SigmaH_SCAT = 1.00
              SigmaP_SCAT = 0.50
           else
              SigmaH_EUV = 0.0
              SigmaP_EUV = 0.0
              SigmaH_SCAT = 1.00*(10.00**cos_SZA)
              SigmaP_SCAT = 0.50*(10.00**cos_SZA)
           end if

           SigmaH_STAR = StarLightPedConductance*2.0
           SigmaP_STAR = StarLightPedConductance

           SigmaH(i,j) = sqrt(SigmaH_EUV*SigmaH_EUV +                         &
                SigmaH_SCAT*SigmaH_SCAT +                                     &
                SigmaH_STAR*SigmaH_STAR)
           SigmaP(i,j) = sqrt(SigmaP_EUV*SigmaP_EUV +                         &
                SigmaP_SCAT*SigmaP_SCAT +                                     &
                SigmaP_STAR*SigmaP_STAR)

        enddo
     enddo

  endif

  if (iModel.ge.3) then

     do j = 1, nPsi
        do i = 1, nTheta

           Sigma0(i,j) = 1000.00

           cos_SZA = cy(i,j)

           if (cos_SZA > 0) then
              SigmaH_EUV = f107p53*(0.81*cos_SZA + 0.54*sqrt(cos_SZA))
              SigmaP_EUV = f107p49*(0.34*cos_SZA + 0.93*sqrt(cos_SZA))
              SigmaH_SCAT = 1.00
              SigmaP_SCAT = 0.50
           else
              SigmaH_EUV = 0.0
              SigmaP_EUV = 0.0
              SigmaH_SCAT = 1.00*(10.00**cos_SZA)
              SigmaP_SCAT = 0.50*(10.00**cos_SZA)
           end if

           SigmaH_STAR = StarLightPedConductance*2.0
           SigmaP_STAR = StarLightPedConductance

           !\
           ! Use Robinson's Formula to convert the Ave_E and E_Flux to 
           ! SigmaP and SigmaH
           !/

           SigmaP_Particles = 40.0 * Ave_E(i,j) /                     &
                (16.0 + Ave_E(i,j)*Ave_E(i,j))  *                     &
                sqrt(EFlux(i,j)*1000.0)

           SigmaH_Particles = 0.45 * (Ave_E(i,j)**0.85) * SigmaP_Particles

           SigmaH(i,j) = sqrt(SigmaH_EUV*SigmaH_EUV + &
                SigmaH_SCAT*SigmaH_SCAT + &
                SigmaH_STAR*SigmaH_STAR + &
                SigmaH_Particles*SigmaH_Particles)

           SigmaP_EUV = SigmaP_EUV*SigmaP_EUV + SigmaP_SCAT*SigmaP_SCAT + &
                SigmaP_STAR*SigmaP_STAR

           SigmaP_Particles = SigmaP_Particles*SigmaP_Particles

           SigmaP(i,j) = sqrt(SigmaP_EUV + SigmaP_Particles)

        enddo

     enddo

  endif

  do j = 1, nPsi
     do i = 1, nTheta

        sn = sin(Theta(i,j))
        cs = cos(Theta(i,j))
        sn2= sn*sn
        cs2 = cs*cs
        cs3 = 1.00 + 3.00*cs2
        cs4 = sqrt(cs3)
        C = 4.00*Sigma0(i,j)*cs2 + &
            SigmaP(i,j)*sn2

        SigmaThTh(i,j) = Sigma0(i,j)*SigmaP(i,j)*cs3/C
        SigmaThPs(i,j) = 2.00*Sigma0(i,j)*SigmaH(i,j)* &
                         cs*cs4/C
        SigmaPsPs(i,j) = SigmaP(i,j)+ &
                         SigmaH(i,j)*SigmaH(i,j)* &
                         sn2/C
        
        dSigmaThTh_dTheta(i,j) = 0.00
        dSigmaThTh_dPsi(i,j) = 0.00
        dSigmaThPs_dTheta(i,j) = 0.00
        dSigmaThPs_dPsi(i,j) = 0.00
        dSigmaPsPs_dTheta(i,j) = 0.00
        dSigmaPsPs_dPsi(i,j) = 0.00

     end do  
  end do

  do j = 1, nPsi
     if (j > 1 .and. j < nPsi ) then 
        do i = 2, nTheta-1
           dSigmaThTh_dTheta(i,j) = (SigmaThTh(i+1,j)-SigmaThTh(i-1,j))/ &
                                    (dTheta(i))
           dSigmaThTh_dPsi(i,j) = (SigmaThTh(i,j+1)-SigmaThTh(i,j-1))/ &
                                  (dPsi(j))
           
           dSigmaThPs_dTheta(i,j) = (SigmaThPs(i+1,j)-SigmaThPs(i-1,j))/ &
                                    (dTheta(i))
           dSigmaThPs_dPsi(i,j) = (SigmaThPs(i,j+1)-SigmaThPs(i,j-1))/ &
                                  (dPsi(j))
           
           dSigmaPsPs_dTheta(i,j) = (SigmaPsPs(i+1,j)-SigmaPsPs(i-1,j))/ &
                                    (dTheta(i))
           dSigmaPsPs_dPsi(i,j) = (SigmaPsPs(i,j+1)-SigmaPsPs(i,j-1))/ &
                                  (dPsi(j))
        end do
     else if (j == 1) then
        do i = 2, nTheta-1
           dSigmaThTh_dTheta(i,j) = (SigmaThTh(i+1,j)-SigmaThTh(i-1,j))/ &
                                    (dTheta(i))
           dSigmaThTh_dPsi(i,j) = (SigmaThTh(i,j+1)-SigmaThTh(i,nPsi-1))/ &
                                  (dPsi(j))
           
           dSigmaThPs_dTheta(i,j) = (SigmaThPs(i+1,j)-SigmaThPs(i-1,j))/ &
                                    (dTheta(i))
           dSigmaThPs_dPsi(i,j) = (SigmaThPs(i,j+1)-SigmaThPs(i,nPsi-1))/ &
                                  (dPsi(j))
           
           dSigmaPsPs_dTheta(i,j) = (SigmaPsPs(i+1,j)-SigmaPsPs(i-1,j))/ &
                                    (dTheta(i))
           dSigmaPsPs_dPsi(i,j) = (SigmaPsPs(i,j+1)-SigmaPsPs(i,nPsi-1))/ &
                                  (dPsi(j))
        end do
     else
        do i = 2, nTheta-1
           dSigmaThTh_dTheta(i,j) = (SigmaThTh(i+1,j)-SigmaThTh(i-1,j))/ &
                                    (dTheta(i))
           dSigmaThTh_dPsi(i,j) = (SigmaThTh(i,2)-SigmaThTh(i,j-1))/ &
                                  (dPsi(j))
           
           dSigmaThPs_dTheta(i,j) = (SigmaThPs(i+1,j)-SigmaThPs(i-1,j))/ &
                                    (dTheta(i))
           dSigmaThPs_dPsi(i,j) = (SigmaThPs(i,2)-SigmaThPs(i,j-1))/ &
                                  (dPsi(j))
           
           dSigmaPsPs_dTheta(i,j) = (SigmaPsPs(i+1,j)-SigmaPsPs(i-1,j))/ &
                                    (dTheta(i))
           dSigmaPsPs_dPsi(i,j) = (SigmaPsPs(i,2)-SigmaPsPs(i,j-1))/ &
                                  (dPsi(j))
        end do
     end if
  end do

end subroutine ionosphere_conductance



subroutine Determine_Oval_Characteristics(Current_in, Theta_in, Psi_in, &
                                          Loc_of_Oval, Width_of_Oval, &
                                          Strength_of_Oval)

!
! This routine calculates everything in radians away from the pole.
!

  use ModIonosphere
  implicit none

!
! inputs:
!

  real, dimension(1:IONO_nTheta, 1:IONO_nPsi), intent(in) :: &
       Current_in, Theta_in, Psi_in

!
! Outputs:
!

  real, dimension(1:IONO_nPsi), intent(out) :: &
       Loc_of_Oval, Width_of_Oval, Strength_of_Oval

!
! Working Variables:
!

  real, dimension(1:IONO_nTheta, 1:IONO_nPsi) :: &
       Current, Theta, Psi

  real, dimension(1:8) :: &
       max_fac, max_fac_colat, width, J_Save

  real    :: day_colat, dusk_colat, midnight_colat, dawn_colat
  real    :: day_fac, dusk_fac, midnight_fac, dawn_fac
  real    :: noon_mid, dusk_dawn, day_strength, night_strength
  real    :: mean_colat, dev_colat, sum, Night_Width, Day_Width

  integer :: i, j, n, nloc, dJ, J_Start, J_End

  logical :: north

!
! Reverse the Arrays for Southern Hemisphere:
!

  if (Theta_in(1,1) < IONO_PI/2.0) then
     Current = Current_in
     Theta   = Theta_in
     Psi     = Psi_in
     north   = .true.
  else
     do i = 1, IONO_nTheta
        do j = 1, IONO_nPsi
           Current(IONO_nTheta - (i-1), j) = Current_in(i,j)
           Theta(IONO_nTheta - (i-1), j)   = IONO_PI - Theta_in(i,j)
           Psi(IONO_nTheta - (i-1), j)     = Psi_in(i,j)
        enddo
     enddo
     north   = .false.
  endif

!
! Start the Oval Determination
!

  dJ = IONO_nPsi/8

  do n = 1, 8 

     !
     ! figure out location of auroral oval:
     ! Let's start a little ways away from the pole...
     !

     J_Start = max(1,(n-1)*dJ)
     J_End   = min(n*dJ, IONO_nPsi)
     max_fac(n) = abs(Current(3,J_Start))
     max_fac_colat(n) = Theta(3,J_Start)
     J_Save(n) = J_Start 
     do j = J_Start, J_End
        do i = 4, IONO_nTheta
           if (abs(Current(i,j)) > max_fac(n)) then
              max_fac(n) = abs(Current(i,j))
              max_fac_colat(n) = Theta(i,j)
              J_Save(n) = j
              nloc = i
           endif
        enddo
     enddo

     !
     ! figure out width
     !

     width(n) = 0.0
     j = J_Save(n)
     do i = nloc, IONO_nTheta
        if (abs(Current(i,j)) > max_fac(n)/5.0) then
           width(n) = abs(max_fac_colat(n) - Theta(i,j))
        endif
     enddo

     if (width(n).eq.0.0) width(n) = max_fac_colat(n)/2.0

  enddo

  night_width = 0.0
  do n=1,8
     night_width = night_width + width(n) * max_fac(n)
  enddo

  day_colat = (max_fac_colat(1) + max_fac_colat(8))/2.0
  dusk_colat = (max_fac_colat(2) + max_fac_colat(3))/2.0
  midnight_colat = (max_fac_colat(4) + max_fac_colat(5))/2.0
  dawn_colat = (max_fac_colat(6) + max_fac_colat(7))/2.0

  day_fac = (max_fac(1) + max_fac(8))/2.0
  dusk_fac = (max_fac(2) + max_fac(3))/2.0
  midnight_fac = (max_fac(4) + max_fac(5))/2.0
  dawn_fac = (max_fac(6) + max_fac(7))/2.0

  sum = 0.0
  mean_colat = 0.0
  
  do n=1,8
     mean_colat = mean_colat + max_fac_colat(n) * max_fac(n)
     sum = sum + max_fac(n)
  enddo

  mean_colat = mean_colat/sum

  Night_Width = Night_Width/sum

  if (Night_Width > 6.0*IONO_PI/180.0) Night_Width=6.0*IONO_PI/180.0
  if (Night_Width < 2.0*IONO_PI/180.0) Night_Width=2.0*IONO_PI/180.0
  Day_Width = max(Night_Width/2.0,1.0*IONO_PI/180.0)

  if (mean_colat < 15.0*IONO_PI/180.0) then
     mean_colat = 15.0*IONO_PI/180.0
     dev_colat = 0.0
  else
     dev_colat = ((day_colat - mean_colat) * day_fac - &
                  (midnight_colat - mean_colat) * midnight_fac) / &
                  (day_fac + midnight_fac)

     ! Restrict auroral location a little bit:

     if (abs(dev_colat) > mean_colat/2.0) then
        dev_colat = dev_colat*mean_colat/2.0/abs(dev_colat)
     endif
  endif

  Day_Strength = dawn_fac
  Night_Strength = sum

  do j=1,IONO_nPsi
     Loc_of_Oval(j)   = mean_colat + dev_colat*cos(Psi(1,j))
     Width_of_Oval(j) = Day_Width + &
          (Night_Width - Day_Width)*sin(Psi(1,j)/2.0)
     Strength_of_Oval(j) = Day_Strength + &
          (Night_Strength - Day_Strength)*sin(Psi(1,j)/2.0)
  enddo

end subroutine Determine_Oval_Characteristics


