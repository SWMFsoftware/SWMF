!^CFG COPYRIGHT UM
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

subroutine FACs_to_fluxes(iModel, iBlock)

  !\
  ! The goal here is to convert the ionospheric FAC pattern into a 
  ! particle precipitation pattern, which can then be turned into
  ! a conductance pattern.
  !/

  use ModIonosphere
  use IE_ModMain
  implicit none

  integer, intent(in) :: iModel  ! model number, 
                                 ! same as conductance_model in IE_ModMain

  integer, intent(in) :: iBlock  ! 1 for northern, 2 for southern hemisphere

  real    :: aurora_peak_lat, eflux_peak_lon, ave_e_peak_lon
  real    :: aurora_peak_width, ave_e_peak
  real    :: distance, mult_fac, ef, efp, mfc
  real    :: hal_a0, hal_a1, hal_a2, ped_a0, ped_a1, ped_a2, hall, ped
  real    :: dlat, dmlt, y1, y2, x1, x2
  real    :: day_colat, dusk_colat, midnight_colat, dawn_colat
  real    :: day_fac, dusk_fac, midnight_fac, dawn_fac
  real    :: noon_mid, dusk_dawn, oval_shift
  real    :: mean_colat, dev_colat, sum, night_width
  real    :: PolarCapHallConductance, PolarCap_AveE, PolarCap_EFlux
  integer :: jlat, imlt
  integer :: i,j, n, nloc
  real, dimension(1:IONO_nPsi) :: Strength_of_Oval,               &
       Loc_of_Oval, Width_of_Oval
  logical :: polarcap, IsPeakFound
  integer :: Poleward, Equatorward, Center, Width
  real :: f, MaxP, MaxT, MulFac
  real :: nden(IONO_nTheta,IONO_nPsi)
  !---------------------------------------------------------------------------
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

  if (iModel.eq.3) then

     mult_fac   = 0.25*0.33e7*2.0e-3
     ave_e_peak = max(3.0, PolarCap_AveE)                        ! keV

     oval_shift = 10.0*cDegToRad

     select case(iBlock)
     case(1)
        !\
        ! Do North first -----------------------------------------------------
        !/

        call Determine_Oval_Characteristics(IONO_NORTH_JR, &
             IONO_NORTH_Theta, &
             IONO_NORTH_Psi, &
             Loc_of_Oval, &
             Width_of_Oval, &
             Strength_of_Oval)

        ave_e_peak_lon = 0.0

        do i = 1, IONO_nTheta
           do j = 1, IONO_nPsi

              ! Energy Flux

              distance = IONO_NORTH_Theta(i,j)-(Loc_of_Oval(j)+oval_shift)
              IONO_NORTH_EFlux(i,j) =                                       &
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

              IONO_NORTH_Ave_E(i,j) =                                      &
                   ave_e_peak*exp(-1.0*(distance/Width_of_Oval(j))**2)

              if (PolarCap_AveE == IONO_Min_Ave_E) then
                 distance = 0.25 +                                          &
                      0.75*(                                                &
                      sin(IONO_NORTH_Psi(i,j)-(ave_e_peak_lon + cHalfPi)) &
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

     case(2)
        !\
        ! Do South next --------------------------------------------------------
        !/

        call Determine_Oval_Characteristics(IONO_SOUTH_JR, &
             IONO_SOUTH_Theta, &
             IONO_SOUTH_Psi, &
             Loc_of_Oval, &
             Width_of_Oval, &
             Strength_of_Oval)

        Loc_of_Oval = cPi - Loc_of_Oval

        ave_e_peak_lon = 0.0

        do i = 1, IONO_nTheta
           do j = 1, IONO_nPsi

              ! Energy Flux

              distance = IONO_SOUTH_Theta(i,j)-(Loc_of_Oval(j)-oval_shift)
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

              IONO_SOUTH_Ave_E(i,j) =                                      &
                   ave_e_peak*exp(-1.0*(distance/Width_of_Oval(j))**2)

              if (PolarCap_AveE == IONO_Min_Ave_E) then
                 distance = 0.25 +                                          &
                      0.75*(                                                &
                      sin(IONO_SOUTH_Psi(i,j)-(ave_e_peak_lon + cHalfPi))   &
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
     end select
  endif

  if (iModel.eq.4 .or. iModel.eq.5) then

     i_cond_nlats = 16
     i_cond_nmlts = 24

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

     dlat = (cond_lats(1) - cond_lats(2))*cDegToRad
     dmlt = (cond_mlts(2) - cond_mlts(1))*cPi/12.0

     select case(iBlock)
     case(1)
        !  Do North First

        call Determine_Oval_Characteristics(IONO_NORTH_JR, &
             IONO_NORTH_Theta, &
             IONO_NORTH_Psi, &
             Loc_of_Oval, &
             Width_of_Oval, &
             Strength_of_Oval)

        do j = 1, IONO_nPsi
           do i = 1, IONO_nTheta

              y1 = IONO_NORTH_Theta(i,j)/dlat + 1.0
              if (y1 > i_cond_nlats-1) then
                 jlat = i_cond_nlats-1
                 y1   = 1.0
              else
                 jlat = y1
                 y1   = 1.0 - (y1 - jlat)
              endif
              y2 = 1.0 - y1

              x1 = mod((IONO_NORTH_Psi(i,j) + cPi), cTwoPi)/dmlt + 1.0
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

              distance = (IONO_NORTH_Theta(i,j) - Loc_of_Oval(j))

              if (distance > 0.0) then
                 polarcap = .false.
              else
                 polarcap = .true.
              endif

              if (iModel.eq.4) then
                 ! Implemented Feb. 7, 2007 as modified version of iModel 5 with
                 !    new narrower fitting of the auroral oval.  DDZ

                 hall=exp(-1.0*(distance/(OvalWidthFactor*Width_of_Oval(j)))**2) * &
                      OvalStrengthFactor*( &
                      hal_a0+(hal_a1-hal_a0)*exp(-abs(iono_north_jr(i,j)*1.0e9)*hal_a2**2))
                 ped =exp(-1.0*(distance/(OvalWidthFactor*Width_of_Oval(j)))**2) * &
                      OvalStrengthFactor*( &
                      ped_a0+(ped_a1-ped_a0)*exp(-abs(iono_north_jr(i,j)*1.0e9)*ped_a2**2))
              else  ! iModel=5
                 !
                 ! We want minimal conductance lower than the oval
                 !

                 if (.not.polarcap) then
                    distance = distance/3.0
                    hal_a0 = hal_a0 * exp(-1.0*(distance/(Width_of_Oval(j)))**2)
                    ped_a0 = ped_a0 * exp(-1.0*(distance/(Width_of_Oval(j)))**2)
                 end if

                 !
                 ! A sort of correction on the fit
                 !

                 hal_a1 = hal_a0 + (hal_a1 - hal_a0)*  &
                      exp(-1.0*(distance/Width_of_Oval(j))**2)
                 ped_a1 = ped_a0 + (ped_a1 - ped_a0)*  &
                      exp(-1.0*(distance/Width_of_Oval(j))**2)

                 ! Multiply by sqrt(3) to compensate for the 3 times narrower oval
                 hall=1.7*( &
                      hal_a0-hal_a1*exp(-abs(iono_north_jr(i,j)*1.0e9)*hal_a2**2))
                 ped =1.7*( &
                      ped_a0-ped_a1*exp(-abs(iono_north_jr(i,j)*1.0e9)*ped_a2**2))
              end if

              if ((hall.gt.1.0).and.(ped.gt.0.5)) then
                 IONO_NORTH_Ave_E(i,j)  = ((hall/ped)/0.45)**(1.0/0.85)
                 IONO_NORTH_EFlux(i,j) = (ped*(16.0+IONO_NORTH_Ave_E(i,j)**2)/&
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

     case(2)
        !  Do South Next

        call Determine_Oval_Characteristics(IONO_SOUTH_JR, &
             IONO_SOUTH_Theta, &
             IONO_SOUTH_Psi, &
             Loc_of_Oval, &
             Width_of_Oval, &
             Strength_of_Oval)

        Loc_of_Oval = cPI - Loc_of_Oval

        do j = 1, IONO_nPsi
           do i = 1, IONO_nTheta


              y1 = (cPi - IONO_SOUTH_Theta(i,j))/dlat + 1.0
              if (y1 > i_cond_nlats-1) then
                 jlat = i_cond_nlats-1
                 y1   = 1.0
              else
                 jlat = y1
                 y1   = 1.0 - (y1 - jlat)
              endif
              y2 = 1.0 - y1

              x1 = mod((IONO_SOUTH_Psi(i,j) + cPi), cTwoPi)/dmlt + 1.0
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

              distance = (IONO_SOUTH_Theta(i,j)-Loc_of_Oval(j))

              if (distance < 0.0) then
                 polarcap = .false.
              else
                 polarcap = .true.
              endif

              if (iModel.eq.4) then
                 hall=exp(-1.0*(distance/(OvalWidthFactor*Width_of_Oval(j)))**2) * &
                      OvalStrengthFactor*( &
                      hal_a0+(hal_a1-hal_a0)*exp(-abs(iono_north_jr(i,j)*1.0e9)*hal_a2**2))
                 ped =exp(-1.0*(distance/(OvalWidthFactor*Width_of_Oval(j)))**2) * &
                      OvalStrengthFactor*( &
                      ped_a0+(ped_a1-ped_a0)*exp(-abs(iono_north_jr(i,j)*1.0e9)*ped_a2**2))
              else  ! iModel=5
                 if (.not.polarcap) then
                    distance = distance/3.0
                    hal_a0 = hal_a0 * exp(-1.0*(distance/(Width_of_Oval(j)))**2)
                    ped_a0 = ped_a0 * exp(-1.0*(distance/(Width_of_Oval(j)))**2)
                 endif

                 !
                 ! A sort of correction on the fit
                 !

                 hal_a1 = hal_a0 + (hal_a1 - hal_a0)*  &
                      exp(-1.0*(distance/Width_of_Oval(j))**2)
                 ped_a1 = ped_a0 + (ped_a1 - ped_a0)*  &
                      exp(-1.0*(distance/Width_of_Oval(j))**2)

                 ! Multiply by sqrt(3) to compensate for the 3 times narrower oval
                 hall=1.7*( &
                      hal_a0-hal_a1*exp(-abs(iono_south_jr(i,j)*1.0e9)*hal_a2**2))
                 ped =1.7*( &
                      ped_a0-ped_a1*exp(-abs(iono_south_jr(i,j)*1.0e9)*ped_a2**2))
              end if

              if ((hall.gt.1.0).and.(ped.gt.0.5)) then
                 IONO_SOUTH_Ave_E(i,j)  = ((hall/ped)/0.45)**(1.0/0.85)
                 IONO_SOUTH_EFlux(i,j) = (ped*(16.0+IONO_SOUTH_Ave_E(i,j)**2)/&
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
     end select
  end if

  if (iModel.eq.6) then

     MulFac = 0.75e7

     if (iBlock == 1) then 

        iono_north_eflux = 1.0e-6
        iono_north_ave_e = 1.0

        ! InvB is -1 in the pole, so we add 2 to make it not negative and not 0
        nden = iono_north_rho / (iono_north_invb+2.0) / 1.6726e-27 * 100.0**3+1
        iono_north_t = iono_north_p / (nden * 1.3807e-23)

        ! Change units from K to eV, and from ions to electrons
        iono_north_t = iono_north_t / 11000.0 * 7.0
        
        do j = 1, IONO_nPsi
           IsPeakFound = .false.
           Poleward = -1
           Equatorward = -1
           MaxP = max(maxval(iono_north_p(:,j)),1.0e-15)
           do i = 1, IONO_nTheta
              if (iono_north_p(i,j) > 0 .and. .not. IsPeakFound) then
                 iono_north_eflux(i,j) = MulFac*MaxP
                 if (iono_north_p(i,j) == MaxP) IsPeakFound = .true.
                 if (IsPeakFound) Equatorward = i
                 if (Poleward == -1) Poleward = i
              endif
           enddo
           if (.not. IsPeakFound) then
              Poleward = 10
              Equatorward = 20
           endif
           if (Poleward < 3) Poleward = 3
           Width = (Equatorward - Poleward)/2 + 1
           Center = (Equatorward + Poleward)/2
           iono_north_eflux(Poleward-2:Equatorward+2,j) = MulFac*MaxP

           ! This is the way that it should be, but it doesn't seem to work...
           MaxT = sum(iono_north_t(Poleward-2:Equatorward+2,j)) / &
                (Equatorward-Poleward+5)

           ! Divide by 1000 for eV to keV
           iono_north_ave_e(Poleward-2:Equatorward+2,j) = MaxT / 1000.0 / 300.0

           ! This is the way that seems to work
           iono_north_ave_e(:,j) = &
                iono_north_p(:,j) / (iono_north_rho(:,j)+1.0e-32) / 1.5e10
           iono_north_ave_e(Poleward-2,j) = iono_north_ave_e(Poleward,j)
           iono_north_ave_e(Poleward-1,j) = iono_north_ave_e(Poleward,j)

           do i = Poleward-2, Equatorward+2
              f = exp(-abs(float(i-Center))/float(width))
              iono_north_eflux(i,j) = iono_north_eflux(i,j) * f
           enddo
        enddo

     else

        iono_south_eflux = 1.0e-6
        iono_south_ave_e = 1.0
        
        ! T = P/nk, where n = rho/mp, but n=

        nden = iono_south_rho / (iono_south_invb+2.0) / 1.6726e-27 * 100.0**3+1
        iono_south_t = iono_south_p / nden / 1.3807e-23

        iono_south_t = iono_south_t * 11.0 * 7.0

        do j = 1, IONO_nPsi
           IsPeakFound = .false.
           Poleward = -1
           Equatorward = -1
           MaxP = max(maxval(iono_south_p(:,j)),1.0e-15)
           MaxT = max(maxval(iono_south_t(:,j)),1.0e-3)
           do i = IONO_nTheta, 1, -1
              if (iono_south_p(i,j) > 0 .and. .not. IsPeakFound) then
                 iono_south_eflux(i,j) = MulFac*MaxP
                 iono_south_ave_e(i,j) = MaxT
                 if (iono_south_p(i,j) == MaxP) IsPeakFound = .true.
                 if (Poleward == -1) Poleward = i
                 if (IsPeakFound) Equatorward = i
              endif
           enddo
           if (.not. IsPeakFound) then
              Poleward = Iono_nPsi - 10
              Equatorward = Iono_nPsi - 20
           endif
           if (Poleward > Iono_nPsi - 3) Poleward = Iono_nPsi - 3
           Width = abs(Equatorward - Poleward)/2 + 1
           Center = abs(Equatorward + Poleward)/2
           iono_south_eflux(Poleward+1,j) = MulFac*MaxP
           iono_south_eflux(Poleward+2,j) = MulFac*MaxP
           iono_south_eflux(Equatorward-1,j) = MulFac*MaxP
           iono_south_eflux(Equatorward-2,j) = MulFac*MaxP
           
           ! This is the way that seems to work
           iono_south_ave_e(:,j) = &
                iono_south_p(:,j) / (iono_south_rho(:,j)+1.0e-32) / 1.5e10
           iono_south_ave_e(Poleward+2,j) = iono_south_ave_e(Poleward,j)
           iono_south_ave_e(Poleward+1,j) = iono_south_ave_e(Poleward,j)
 
           do i = Poleward+2, Equatorward-2, -1
              f = exp(-abs(float(i-Center))/float(width))
              iono_south_eflux(i,j) = iono_south_eflux(i,j) * f
           enddo
        enddo

     endif

!!!     ! This works well for ions.  Electrons???  Hmmm....
!!!
!!!     do j = 1, IONO_nPsi
!!!        do i = 1, IONO_nTheta
!!!
!!!           if (IONO_SOUTH_rho(i,j).le.1e-21) then
!!!              IONO_SOUTH_Ave_E(i,j) = 0.0 
!!!           else
!!!              IONO_SOUTH_t(i,j) = &
!!!                   1e17*IONO_SOUTH_p(i,j)/(IONO_SOUTH_rho(i,j)*1.38)
!!!              IONO_SOUTH_Ave_E(i,j) = 0.05*IONO_SOUTH_t(i,j)*1e-27
!!!           endif
!!!
!!!           IONO_SOUTH_EFlux(i,j) = IONO_SOUTH_p(i,j)*1e5
!!!
!!!!           if (IONO_SOUTH_EFlux(i,j).le.1e-3) then
!!!!              IONO_SOUTH_EFlux(i,j) = 1e-21
!!!!           endif
!!!
!!!           if ((PolarCap_AveE > 0.0).and.(IONO_SOUTH_p(i,j).eq.0.0)) then
!!!              IONO_SOUTH_Ave_E(i,j) = max(IONO_SOUTH_Ave_E(i,j), &
!!!                   PolarCap_AveE)
!!!              IONO_SOUTH_EFlux(i,j) = max(IONO_SOUTH_EFlux(i,j),        &
!!!                   PolarCap_EFlux)
!!!           endif
!!!
!!!           if (IONO_NORTH_rho(i,j).le.1e-21) then
!!!              IONO_NORTH_Ave_E(i,j) = 0.0
!!!           else                 
!!!              IONO_NORTH_t(i,j) = &
!!!                   1e17*IONO_NORTH_p(i,j)/(IONO_NORTH_rho(i,j)*1.38)
!!!              IONO_NORTH_Ave_E(i,j) = 0.05*IONO_NORTH_t(i,j)*1e-27
!!!           endif
!!!           IONO_NORTH_EFlux(i,j) = IONO_NORTH_p(i,j)*6e5
!!!!           if (IONO_NORTH_EFlux(i,j).le.1e-3) IONO_NORTH_EFlux(i,j) = 1e-21
!!!
!!!           if ((PolarCap_AveE > 0.0).and.(IONO_NORTH_p(i,j).eq.0.0)) then
!!!              IONO_NORTH_Ave_E(i,j) = max(IONO_NORTH_Ave_E(i,j), &
!!!                   PolarCap_AveE)
!!!              IONO_NORTH_EFlux(i,j) = max(IONO_NORTH_EFlux(i,j), &
!!!                   PolarCap_EFlux)
!!!           endif
!!!
!!!        enddo
!!!     enddo

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

  use IE_ModMain
  use ModIonosphere
  use ModProcIE
  use CON_world
  use CON_comp_param
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
       oSigmaH, oSigmaP, conv_SigmaH, conv_SigmaP, cos_SZA,     &
       tmp_x, tmp_y, tmp_z

  real, dimension(1:IONO_NTheta) :: dTheta
  real, dimension(1:IONO_NPsi)   :: dPsi

  integer :: i,j
  real :: f107p53, f107p49
  real :: sn, cs, sn2, cs2, cs3, cs4, C
  real :: SigmaH_EUV, SigmaP_EUV, SigmaH_SCAT, SigmaP_SCAT, &
       SigmaH_STAR, SigmaP_STAR, SigmaH_EUV_2, SigmaP_EUV_2
  real :: SigmaH_Particles, SigmaP_Particles, tau
  logical :: old

  real    :: time_delay
  real    :: cos_limit, meeting_value_h, meeting_value_p

  !--------------------------------------------------------------------------

  if (theta(1,1) < 2.0*IONO_Theta_0) then
     north = .true.
     tmp_x = IONO_NORTH_X
     tmp_y = IONO_NORTH_Y
     tmp_z = IONO_NORTH_Z
  else
     north = .false.
     tmp_x = IONO_SOUTH_X
     tmp_y = IONO_SOUTH_Y
     tmp_z = IONO_SOUTH_Z
  endif

  if (use_comp(UA_)) then

     Sigma0 = 1000.00

     ! Use Conductances from the UAM, which are put into the following
     ! variables

     if (north) then
        SigmaH = IONO_NORTH_SigmaH
        SigmaP = IONO_NORTH_SigmaP
     else
        SigmaH = IONO_South_SigmaH
        SigmaP = IONO_South_SigmaP
     endif

  else

     if (iModel > 1) then

        cos_SZA = (tmp_x*cosTHETATilt - tmp_z*sinTHETATilt) &
             / sqrt(tmp_x**2 + tmp_y**2 + tmp_z**2)

        ! We are going to need F10.7 ^ 0.53 and F10.7 ^ 0.49 a lot,
        ! So, let's just store them straight away:

        f107p53 = f107**0.53
        f107p49 = f107**0.49
        cos_limit = cos(70.0*cDegToRad)
        meeting_value_p = f107p49*(0.34*cos_limit+0.93*sqrt(cos_limit))
        meeting_value_h = f107p53*(0.81*cos_limit+0.54*sqrt(cos_limit))

     endif

     if (iModel.eq.0) then

        do j = 1, nPsi
           do i = 1, nTheta

              Sigma0(i,j) = 1000.00
              SigmaH(i,j) = 0.00
              SigmaP(i,j) = StarLightPedConductance

           enddo
        enddo

     endif

     if (iModel.eq.1) then

        do j = 1, nPsi
           do i = 1, nTheta

              Sigma0(i,j) = 1000.00
              SigmaH(i,j) = PolarCapPedConductance
              SigmaP(i,j) = StarLightPedConductance

           enddo
        enddo

     endif

     if (iModel.eq.2) then

        do j = 1, nPsi
           do i = 1, nTheta

!!$        !\
!!$        ! Rasmussen and Schunk model B (JGR, Vol. 92, pp. 4491-4504, 1987).
!!$        !/

              Sigma0(i,j) = 1000.00

              if (cos_SZA(i,j) > 0) then
                 SigmaH_EUV=f107p53*(0.81*cos_SZA(i,j)+0.54*sqrt(cos_SZA(i,j)))
                 SigmaP_EUV=f107p49*(0.34*cos_SZA(i,j)+0.93*sqrt(cos_SZA(i,j)))
                 SigmaH_SCAT = 1.00
                 SigmaP_SCAT = 0.50
                 if (cos_SZA(i,j) < cos_limit) then
                    SigmaH_EUV_2 = meeting_value_h *   &
                         exp(-((cos_SZA(i,j)-cos_limit)**2.0)*15.0)
                    SigmaP_EUV_2 = meeting_value_p *   &
                         exp(-((cos_SZA(i,j)-cos_limit)**2.0)*15.0)
                    SigmaH_EUV = (SigmaH_EUV + SigmaH_EUV_2)/2.0
                    SigmaP_EUV = (SigmaP_EUV + SigmaP_EUV_2)/2.0
                 endif
              else
                 SigmaH_EUV = meeting_value_h *   &
                      exp(-((cos_SZA(i,j)-cos_limit)**2.0)*15.0)
                 SigmaP_EUV = meeting_value_p *   &
                      exp(-((cos_SZA(i,j)-cos_limit)**2.0)*15.0)
                 SigmaH_SCAT = 1.00*(10.00**cos_SZA(i,j))
                 SigmaP_SCAT = 0.50*(10.00**cos_SZA(i,j))
              end if

              SigmaH_STAR = StarLightPedConductance*2.0
              SigmaP_STAR = StarLightPedConductance

              SigmaH(i,j) = sqrt(SigmaH_EUV*SigmaH_EUV +                     &
                   SigmaH_SCAT*SigmaH_SCAT +                                 &
                   SigmaH_STAR*SigmaH_STAR)
              SigmaP(i,j) = sqrt(SigmaP_EUV*SigmaP_EUV +                     &
                   SigmaP_SCAT*SigmaP_SCAT +                                 &
                   SigmaP_STAR*SigmaP_STAR)

           enddo
        enddo

     endif

     if (iModel.ge.3) then

        do j = 1, nPsi
           do i = 1, nTheta

              Sigma0(i,j) = 1000.00

              if (cos_SZA(i,j) > 0) then
                 SigmaH_EUV=f107p53*(0.81*cos_SZA(i,j)+0.54*sqrt(cos_SZA(i,j)))
                 SigmaP_EUV=f107p49*(0.34*cos_SZA(i,j)+0.93*sqrt(cos_SZA(i,j)))
                 SigmaH_SCAT = 1.00
                 SigmaP_SCAT = 0.50
                 if (cos_SZA(i,j) < cos_limit) then
                    SigmaH_EUV_2 = meeting_value_h *   &
                         exp(-((cos_SZA(i,j)-cos_limit)**2.0)*15.0)
                    SigmaP_EUV_2 = meeting_value_p *   &
                         exp(-((cos_SZA(i,j)-cos_limit)**2.0)*15.0)
                    SigmaH_EUV = (SigmaH_EUV + SigmaH_EUV_2)/2.0
                    SigmaP_EUV = (SigmaP_EUV + SigmaP_EUV_2)/2.0
                 endif
              else
                 SigmaH_EUV = meeting_value_h *   &
                      exp(-((cos_SZA(i,j)-cos_limit)**2.0)*15.0)
                 SigmaP_EUV = meeting_value_p *   &
                      exp(-((cos_SZA(i,j)-cos_limit)**2.0)*15.0)
                 SigmaH_SCAT = 1.00*(10.00**cos_SZA(i,j))
                 SigmaP_SCAT = 0.50*(10.00**cos_SZA(i,j))
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

     if (north) then

        ! Subsolar is i=nTheta ; j = 0
        ! Just in case we change this in the future, add noon and midnight to
        ! test for conductance

        if (SAVE_NORTH_SigmaH(nTheta,1) +                                     &
             SAVE_NORTH_SigmaH(nTheta,nPsi/2) > 0.0) then

           ! We have "old" conductances

           old = .true.
           SigmaH = SAVE_NORTH_SigmaH
           SigmaP = SAVE_NORTH_SigmaP

        endif

     else

        ! Subsolar is i=0 ; j = 0
        ! Just in case we change this in the future, add noon and midnight to
        ! test for conductance

        if (SAVE_SOUTH_SigmaH(1,1) +                                         &
             SAVE_SOUTH_SigmaH(1,nPsi/2) > 0.0) then

           ! We have "old" conductances

           old = .true.
           SigmaH = SAVE_SOUTH_SigmaH
           SigmaP = SAVE_SOUTH_SigmaP

        endif

     endif

     if (north) then
        SAVE_NORTH_SigmaP = 0.0
        SAVE_NORTH_SigmaH = 0.0
     else
        SAVE_SOUTH_SigmaP = 0.0
        SAVE_SOUTH_SigmaH = 0.0
     endif

  end if

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

  !
  ! Reverse the Arrays for Southern Hemisphere:
  !

  if (Theta_in(1,1) < cHalfPi) then
     Current = Current_in
     Theta   = Theta_in
     Psi     = Psi_in
     north   = .true.
  else
     do i = 1, IONO_nTheta
        do j = 1, IONO_nPsi
           Current(IONO_nTheta - (i-1), j) = Current_in(i,j)
           Theta(IONO_nTheta - (i-1), j)   = cPi - Theta_in(i,j)
           Psi(IONO_nTheta - (i-1), j)     = Psi_in(i,j)
        enddo
     enddo
     north   = .false.
  endif

  !
  ! Start the Oval Determination
  !

  dJ = IONO_nPsi/9

  do n = 1, 8 

     !
     ! figure out location of auroral oval:
     ! Let's start a little ways away from the pole...
     !

     nloc = 1
     J_Start = max(1,(n-1)*dJ)
     J_End   = min(n*dJ, IONO_nPsi)
     max_fac(n) = abs(Current(3,J_Start))
     max_fac_colat(n) = Theta(3,J_Start)
     J_Save(n) = J_Start 
     do j = J_Start, J_End
        do i = 4, IONO_nTheta
           if (Current(i,j) > max_fac(n)) then
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
        if (Current(i,j) > max_fac(n)/4.0) then
           width(n) = abs(max_fac_colat(n) - Theta(i,j))
        endif
     enddo

     if (width(n).le.(theta(2,1)-theta(1,1))) width(n) = max_fac_colat(n)/5.0

  enddo

  day_colat = (max_fac_colat(1) + max_fac_colat(8))/2.0
  dusk_colat = (max_fac_colat(2) + max_fac_colat(3))/2.0
  midnight_colat = (max_fac_colat(4) + max_fac_colat(5))/2.0
  dawn_colat = (max_fac_colat(6) + max_fac_colat(7))/2.0

  day_fac = (max_fac(1) + max_fac(8))/2.0
  dusk_fac = (max_fac(2) + max_fac(3))/2.0
  midnight_fac = (max_fac(4) + max_fac(5))/2.0
  dawn_fac = (max_fac(6) + max_fac(7))/2.0

  night_width = 0.0
  sum = 0.0
  mean_colat = 0.0

  do n=1,8
     night_width = night_width + width(n) * max_fac(n)
     mean_colat = mean_colat + max_fac_colat(n) * max_fac(n)
     sum = sum + max_fac(n)
  enddo

  mean_colat = mean_colat/sum
  Night_Width = Night_Width/sum

  if (Night_Width > 6.0*cDegToRad) Night_Width=6.0*cDegToRad
  Day_Width = max(Night_Width/2.0,1.0*cDegToRad)

  if (mean_colat < 15.0*cDegToRad) then
     mean_colat = 15.0*cDegToRad
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


  !  dev_colat = 1.0*cDegToRad

  !  Day_Width = 2.5*cDegToRad
  !  Night_Width = Day_Width

  !  mean_colat = 2.0*cDegToRad

  do j=1,IONO_nPsi
     Loc_of_Oval(j)   = mean_colat - dev_colat*cos(Psi(1,j))
     Width_of_Oval(j) = Day_Width + &
          (Night_Width - Day_Width)*sin(Psi(1,j)/2.0)
     Strength_of_Oval(j) = Day_Strength + &
          (Night_Strength - Day_Strength)*sin(Psi(1,j)/2.0)
  enddo

end subroutine Determine_Oval_Characteristics


