subroutine aurora(iBlock)

  use ModGITM
  use ModSources
  use ModTime, only : tSimulation
  use ModInputs
  use ModConstants
  use ModUserGITM

  implicit none

  integer, intent(in) :: iBlock

  real :: alat, hpi, ped, hal, av_kev, eflx_ergs, a,b, maxi
  real :: Factor,temp_ED, avee, eflux, p
  integer :: i, j, k, n, iError, iED
  logical :: IsDone, IsTop

  real, dimension(nLons,nLats,nAlts) :: temp, AuroralBulkIonRate, &
       IonPrecipitationBulkIonRate, IonPrecipitationHeatingRate

  logical :: IsFirstTime(nBlocksMax) = .true.

  if (IsFirstTime(iBlock)) then
     IsFirstTime(iBlock) = .false.

     if (iBlock == 1 .and. UseIonPrecipitation) then
        call ReadIonHeat(IonIonizationFilename,  .true.) 
        call ReadIonHeat(IonHeatingRateFilename, .false.) 
     endif

  else
     if (floor((tSimulation - dT)/dTAurora) == &
          floor(tSimulation/dTAurora)) return
  endif

  AuroralBulkIonRate               = 0.0
  AuroralHeatingRate(:,:,:,iBlock) = 0.0

  call report("Aurora",1)
  call start_timing("Aurora")

  if (UseIonPrecipitation) call interpolate_ions( &
       nLons, nLats, nAlts, &
       Longitude(1:nLons,iBlock), Latitude(1:nLats,iBlock), &
       Altitude_GB(1:nLons, 1:nLats, 1:nAlts, iBlock),&
       IonPrecipitationBulkIonRate, IonPrecipitationHeatingRate)

  do i=1,nLats
     do j=1,nLons

        eflx_ergs = ElectronEnergyFlux(j,i) !/ (1.0e-7 * 100.0 * 100.0)
        av_kev    = ElectronAverageEnergy(j,i)

        UserData2d(j,i,1,2,iBlock) = av_kev
        UserData2d(j,i,1,3,iBlock) = eflx_ergs

        if (eflx_ergs > 0.02) then

           avee = av_kev * 1000.0
           eflux = eflx_ergs * 6.242e11

           a= sqrt(27.0/(2.0*3.14159)) * eflux /(avee**2.5)

           ED_flux(ED_N_energies) = 0.0 !!! But why ???
           
           do n=1,ED_N_Energies-1
              ED_flux(n) = &
                   a*sqrt(ed_energies(n))*exp(-1.5*ed_energies(n)/avee)/pi
           enddo

           call R_ELEC_EDEP (ED_Flux, 30, ED_Energies, 3, ED_Ion, 7)
           call R_ELEC_EDEP (ED_Flux, 30, ED_Energies, 3, ED_Heating, 11)

           iED = 1

           do k = 1, nAlts

              p = alog(Pressure(j,i,k,iBlock))

              IsDone = .false.
              IsTop = .false.
              do while (.not.IsDone)
                 if (ED_grid(iED) >= p .and. ED_grid(iED+1) <= p) then
                    IsDone = .true.
                    ED_Interpolation_Index(k) = iED
                    ED_Interpolation_Weight(k) = (ED_grid(iED) - p) /  &
                         (ED_grid(iED) - ED_grid(iED+1))
                 else
                    if (iED == ED_N_Alts-1) then
                       IsDone = .true.
                       IsTop = .true.
                    else
                       iED = iED + 1
                    endif
                 endif
              enddo

              if (.not.IsTop) then
                 n = ED_Interpolation_Index(k)
                 AuroralBulkIonRate(j,i,k) = ED_Ion(n) - &
                      (ED_Ion(n) - ED_Ion(n+1))*ED_Interpolation_Weight(k)
                 AuroralHeatingRate(j,i,k,iBlock) = ED_Heating(n) - &
                      (ED_Heating(n) - ED_Heating(n+1))*ED_Interpolation_Weight(k)
              else

                 ! Decrease after top of model
                 AuroralBulkIonRate(j,i,k) = ED_Ion(ED_N_Alts) * &
                      Pressure(j,i,k,iBlock) / exp(ED_grid(ED_N_Alts)) 
                 AuroralHeatingRate(j,i,k,iBlock) = ED_Heating(ED_N_Alts) * &
                      Pressure(j,i,k,iBlock) / exp(ED_grid(ED_N_Alts))

              endif

              ! write(*,*) "aurora : ",k,IsTop, AuroralBulkIonRate(j,i,k)

           enddo

        endif

     enddo
  enddo

  ! From Rees's book:

  temp = 0.92 * NDensityS(1:nLons,1:nLats,1:nAlts,iN2_,iBlock) + &
         1.00 * NDensityS(1:nLons,1:nLats,1:nAlts,iO2_,iBlock) + &
         0.56 * NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock)

  if (UseIonPrecipitation) then

     IonPrecipIonRateS(:,:,:,iO_3P_,iBlock)  = &
          0.56*IonPrecipitationBulkIonRate*&
          NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock)/temp
     IonPrecipIonRateS(:,:,:,iO2_,iBlock) = &
          1.00*IonPrecipitationBulkIonRate*&
          NDensityS(1:nLons,1:nLats,1:nAlts,iO2_,iBlock)/temp
     IonPrecipIonRateS(:,:,:,iN2_,iBlock) = &
          0.92*IonPrecipitationBulkIonRate*&
          NDensityS(1:nLons,1:nLats,1:nAlts,iN2_,iBlock)/temp

     IonPrecipHeatingRate(:,:,:,iBlock) = IonPrecipitationHeatingRate

  else

     IonPrecipIonRateS(:,:,:,:,iBlock) = 0.0
     IonPrecipHeatingRate(:,:,:,iBlock) = 0.0

  endif

  AuroralIonRateS(:,:,:,iO_3P_,iBlock)  = &
       0.56*AuroralBulkIonRate*&
       NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock)/temp
  AuroralIonRateS(:,:,:,iO2_,iBlock) = &
       1.00*AuroralBulkIonRate*&
       NDensityS(1:nLons,1:nLats,1:nAlts,iO2_,iBlock)/temp
  AuroralIonRateS(:,:,:,iN2_,iBlock) = &
       0.92*AuroralBulkIonRate*&
       NDensityS(1:nLons,1:nLats,1:nAlts,iN2_,iBlock)/temp

  call end_timing("Aurora")

end subroutine aurora




