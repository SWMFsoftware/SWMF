! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

subroutine calc_neutral_density(iBlock)

  use ModSize
  use ModGITM
  use ModPlanet
  use ModRates
  use ModEUV
  use ModSources
  use ModInputs, only: iDebugLevel

  implicit none

  integer, intent(in) :: iBlock

  real :: Sources(nSpeciesTotal)
  real :: Losses(nSpeciesTotal)
  real :: DtSub, DtOld, DtTotal, DtMin, ChemicalHeatingSub, Source

  integer :: iLon, iLat, iAlt, iNeutral, iDtReducer, nIters

  DtMin = Dt

  call report("Neutral Chemistry",2)
  call start_timing("Neutral Chemistry")

  nIters = 0

  do iLon = 1, nLons
     do iLat = 1, nLats
        do iAlt = 1, nAlts

           DtTotal = 0.0

           do while (DtTotal < Dt)

              ChemicalHeatingSub = 0.0

              DtSub = Dt - DtTotal

              Sources = 0.0
              Losses  = 0.0

              ! [O]!!!!!!!!!!

              ! Sources:
              Sources(iO_)=  &

                   RrTempDep(iLon,iLat,iAlt,iRrK1_) *   &
                   IDensityS( iLon,iLat,iAlt,iO_4SP_,iBlock)   *  &
                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)     +    &

                   RrTempInd(iRrK32_) * &
                   IDensityS(iLon,iLat,iAlt,iO2P_,iBlock) *  &
                   NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock)         +    &

                   RrTempInd(iRrK7_) * &
                   IDensityS(iLon,iLat,iAlt,iNP_,iBlock)  *  &
                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)           +    &

                   RrTempDep(iLon,iLat,iAlt,iRrA1_) * &
                   IDensityS(iLon,iLat,iAlt,iNOP_,iBlock)*  &
                   IDensityS(iLon,iLat,iAlt,ie_,iBlock)+ &

                   RrTempDep(iLon,iLat,iAlt,iRrBeta1_) * &
                   NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)  +    &

                   RrTempDep(iLon,iLat,iAlt,iRrBeta3_) * &
                   NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock)*&
                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock)  +    &

                   RrTempInd(iRrBeta6_) * &
                   NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock) *     &
                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock)    +    &

                   RrTempDep(iLon,iLat,iAlt,iRrBeta8_) * &
                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock)  +& 

                   RrTempInd(iRrBeta2_) * &
                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)       *     &
                   NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock)       +    &
 
                   1.15*RrTempDep(iLon,iLat,iAlt,iRrA2_) * &
                   IDensityS(iLon,iLat,iAlt,ie_,iBlock) * &
                   IDensityS(iLon,iLat,iAlt,iO2P_,iBlock)       + &

                   2*RrTempInd(iRrJ1_) * &
                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)


    ! Loss:
              Losses(iO_)=  &
!
                   2*RrTempDep(iLon,iLat,iAlt,iRrK12_) * &
                   NDensityS(iLon,iLat,iAlt,iO_,iBlock) * &
                   (NDensityS(iLon,iLat,iAlt,iO_,iBlock) +  &
                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)+  &
                   NDensityS(iLon,iLat,iAlt,iN2_,iBlock))+ &
!
                   RrTempDep(iLon,iLat,iAlt,iRrK13_) * &
                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock) * &
                   (NDensityS(iLon,iLat,iAlt,iO_,iBlock) +  &
                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock))+ &
!
                   RrTempDep(iLon,iLat,iAlt,iRrK3_) * &
                   IDensityS(iLon,iLat,iAlt,iN2P_,iBlock) + &
!
                   RrTempInd(iRrK8_) * &
                   IDensityS(iLon,iLat,iAlt,iNP_, iBlock) + &
!
                   EuvIonRateS(iLon,iLat,iAlt,iO_4SP_,iBlock)

!              Losses(iO_)= Losses(iO_) * &
!                   NDensityS(iLon,iLat,iAlt,iO_,iBlock)
!+&
!                   AuroralIonRateS(iLon,iLat,iAlt,iO_,iBlock)

              ! [O2]!!!!!!!!!!

              ! Sources:

              Sources(iO2_) = &
!
                   RrTempInd(iRrK5_) * &
                   IDensityS(iLon,iLat,iAlt,iO2P_, iBlock) *  &
                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock)+    &
!
                   RrTempDep(iLon,iLat,iAlt,iRrK12_) * &
                   NDensityS(iLon,iLat,iAlt,iO_,iBlock)**2 * &
                   (NDensityS(iLon,iLat,iAlt,iO_,iBlock)  +   &
                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock) +   &
                   NDensityS(iLon,iLat,iAlt,iN2_,iBlock))
              
              ! Loss:
              Losses(iO2_) = &
!
                   RrTempDep(iLon,iLat,iAlt,iRrK1_) * &
                   IDensityS(iLon,iLat,iAlt,iO_4SP_,iBlock) + &
!
                   RrTempInd(iRrK6_) * &
                   IDensityS(iLon,iLat,iAlt,iNP_, iBlock) +           &
!
                   RrTempInd(iRrK7_) * &
                   IDensityS(iLon,iLat,iAlt,iNP_, iBlock) +           &
!
                   RrTempDep(iLon,iLat,iAlt,iRrBeta1_) * &
                   NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock) +&
!
                   RrTempInd(iRrBeta2_) * &
                   NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock) +      &
!
                   RrTempDep(iLon,iLat,iAlt,iRrK13_) * &
                   NDensityS(iLon,iLat,iAlt,iO_,iBlock) *    &
                   (NDensityS(iLon,iLat,iAlt,iO_,iBlock) +    &
                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)) !+   &
!
!                   RrTempInd(iRrJ1_) 


!              Losses(iO2_) = Losses(iO2_) * &
!                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)

              ! [N2]!!!!!!!!!!
              
              ! Sources:
              Sources(iN2_) = &
                   RrTempDep(iLon,iLat,iAlt,iRrBeta3_) * &
                   NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock) + &
                   RrTempInd(iRrBeta6_) *  &
                   NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock)

              ! Loss:
              Losses(iN2_) = &
                   RrTempDep(iLon,iLat,iAlt,iRrK2_) * &
                   IDensityS(iLon,iLat,iAlt,iO_4SP_,iBlock) + & 
                   EuvIonRateS(iLon,iLat,iAlt,iN2P_,iBlock)

!              Losses(iN2_) = Losses(iN2_) * &
!                   NDensityS(iLon,iLat,iAlt,iN2_,iBlock)

              ! [N(4s)]!!!!!!!!

              ! Sources:
              Sources(iN_4S_) = &
                   RrTempDep(iLon,iLat,iAlt,iRrK2_) * &
                   IDensityS(iLon,iLat,iAlt,iO_4SP_,iBlock) *    &
                   NDensityS(iLon,iLat,iAlt,iN2_,iBlock) +       &
                   RrTempInd(iRrK6_) * &
                   IDensityS(iLon,iLat,iAlt,iNP_, iBlock) *           &
                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock) +            &
                   RrTempInd(iRrK8_) * &
                   IDensityS(iLon,iLat,iAlt,iNP_, iBlock)*            &
                   NDensityS(iLon,iLat,iAlt,iO_,iBlock) +             &
                   0.2 * RrTempDep(iLon,iLat,iAlt,iRrA1_) * &
                   IDensityS(iLon,iLat,iAlt,iNOP_,iBlock) *   &
                   IDensityS(iLon,iLat,iAlt,ie_,iBlock) +                     &
                   (2.*0.1+0.9)* RrTempDep(iLon,iLat,iAlt,iRrA3_) * &
                   IDensityS(iLon,iLat,iAlt,iN2P_,iBlock) *   &
                   IDensityS(iLon,iLat,iAlt,ie_,iBlock) +               &
                   RrTempDep(iLon,iLat,iAlt,iRrBeta5_) * &
                   NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock)* &
                   IDensityS(iLon,iLat,iAlt,ie_,iBlock) +     &
                   RrTempInd(iRrBeta7_) * &
                   NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock) +        &
                   RrTempDep(iLon,iLat,iAlt,iRrBeta8_) * &
                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock)

              Source = &
                   RrTempInd(iRrBeta4_) * &
                   NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock) *        &
                   NDensityS(iLon,iLat,iAlt,iO_,iBlock)

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + Source * 2.38

              Sources(iN_4S_) = Sources(iN_4S_) + Source

              ! Loss:
              Losses(iN_4S_) = &
                   RrTempInd(iRrK32_) * &
                   IDensityS(iLon,iLat,iAlt,iO2P_,iBlock) +                  &
                   RrTempDep(iLon,iLat,iAlt,iRrBeta1_) * &
                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock) +  &
                   RrTempDep(iLon,iLat,iAlt,iRrBeta3_) * &
                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock) + &
                   EuvIonRateS(iLon,iLat,iAlt,iNP_,iBlock)

!              Losses(iN_4S_) = Losses(iN_4S_) * &
!                   NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock)

!              Sources(iN_4S_) = 0.0
!              Losses(iN_4S_) = 0.0

              ! [NO]!!!!!!!!!!!
     

              Source = &
                   RrTempDep(iLon,iLat,iAlt,iRrBeta1_)* &
                   NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock)*&
                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + Source * 1.385

              ! Sources:
              Sources(iNO_)= Source 

              Source = &
                   RrTempInd(iRrBeta2_)* &
                   NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock)*&
                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + Source * 3.76

              Sources(iNO_)= Sources(iNO_) + Source
 

              ! Loss:
              Losses(iNO_)= &
                   RrTempInd(iRrK5_)* &
                   IDensityS(iLon,iLat,iAlt,iO2P_, iBlock)+          &
                   RrTempDep(iLon,iLat,iAlt,iRrBeta3_)*&
                   NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock)+&
                   RrTempInd(iRrBeta6_)*&
                   NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock)+      &
                   RrTempDep(iLon,iLat,iAlt,iRrBeta8_)+          &
                   RrTempDep(iLon,iLat,iAlt,iRrBeta9N_)

!              Losses(iNO_) = Losses(iNO_) * &
!                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock)

              ! [N(2D)]!!!!!!!!!
      
              ! Sources:
              Sources(iN_2D_)= &
                   RrTempDep(iLon,iLat,iAlt,iRrK3_) * &
                   IDensityS(iLon,iLat,iAlt,iN2P_,iBlock) *    &
                   NDensityS(iLon,iLat,iAlt,iO_,iBlock) +             &
                   0.8*RrTempDep(iLon,iLat,iAlt,iRrA1_) * &
                   IDensityS(iLon,iLat,iAlt,iNOP_,iBlock) *      &
                   IDensityS(iLon,iLat,iAlt,ie_,iBlock) +&
                   0.9*RrTempDep(iLon,iLat,iAlt,iRrA3_) * &
                   IDensityS(iLon,iLat,iAlt,iN2P_,iBlock) *      &
                   IDensityS(iLon,iLat,iAlt,ie_,iBlock)

              !              RrTempInd(iRrBetar2)=1.0e-20
              !      if (k > 6) RrTempInd(iRrBeta2)=5.0e-19

              ! Loss:
              Losses(iN_2D_)= &
                   RrTempInd(iRrBeta2_) * &
                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock) + &
                   RrTempInd(iRrBeta4_) * &
                   NDensityS(iLon,iLat,iAlt,iO_,iBlock) +  &
                   RrTempDep(iLon,iLat,iAlt,iRrBeta5_) * &
                   IDensityS(iLon,iLat,iAlt,ie_,iBlock) + &
                   RrTempInd(iRrBeta6_) * &
                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock) + &
                   RrTempInd(iRrBeta7_)

!              Losses(iN_2D_) = Losses(iN_2D_) * &
!                   NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock)

!              Sources(iO_) = 1.0e-10
!              Losses(iO_)  = 0.0
!              Sources(iO2_) = 1.0e-10
!              Losses(iO2_)  = 0.0
!              Sources(iN2_) = 1.0e-10
!              Losses(iN2_)  = 0.0
!
!              ! Let's trick some constituents for now....
!              ! Call then steady-state
!
              do iNeutral = 1, nSpeciesTotal
!                 if (Losses(iNeutral) * &
!                      nDensityS(iLon,iLat,iAlt,iNeutral,iBlock) > &
!                      1000.0 * Sources(iNeutral) .and. &
!                      Losses(iNeutral) > 1000.0) then
!
!                    write(*,*) "Solving steady-state : ", iNeutral
!                    call write_sources_losses(iNeutral)
!
!                    NDensityS(iLon,iLat,iAlt,iNeutral,iBlock) = &
!                         Sources(iNeutral) / Losses(iNeutral)
!                    write(*,*) NDensityS(iLon,iLat,iAlt,iNeutral,iBlock)
!                    Sources(iNeutral) = 0.0
!                    Losses(iNeutral) = 0.0
!                 else
                    Losses(iNeutral) = &
                         Losses(iNeutral) * &
                         NDensityS(iLon,iLat,iAlt,iNeutral,iBlock)
!                 endif
              enddo

              do iNeutral = 1, nSpeciesTotal
                 if (Losses(iNeutral) > Sources(iNeutral)) then
                    if (nDensityS(iLon,iLat,iAlt,iNeutral,iBlock) > 0.01) then
                       dtOld = DtSub
                       dtSub = &
                            min(0.25 * &
                            nDensityS(iLon,iLat,iAlt,iNeutral,iBlock)/ &
                            (abs(Sources(iNeutral) - &
                                 Losses(iNeutral))+0.1), DtOld)
                       if (DtSub < DtOld) iDtReducer = iNeutral

                       if (DtSub < DtOld .and. iDebugLevel > 5) &
                            write(*,*) "====> Neutral Chem, dt reduced : ", &
                            DtSub, DtOld,iAlt, iNeutral, &
                            Sources(iNeutral), Losses(iNeutral), &
                            nDensityS(iLon,iLat,iAlt,iNeutral,iBlock)
                    else
                       Sources(iNeutral) = 0.0
                       Losses(iNeutral) = 0.0
                    endif
                 endif
              enddo

              do iNeutral = 1, nSpeciesTotal
                 NDensityS(iLon,iLat,iAlt,iNeutral,iBlock) = &
                      NDensityS(iLon,iLat,iAlt,iNeutral,iBlock) + &
                      (Sources(iNeutral) - Losses(iNeutral)) * DtSub
              enddo

              ChemicalHeatingRate(iLon,iLat,iAlt) = &
                   ChemicalHeatingRate(iLon,iLat,iAlt) + &
                   ChemicalHeatingSub * DtSub

              DtTotal = DtTotal + DtSub

              if (DtSub < DtMin) DtMin = DtSub

              if (DtSub < 1.0e-7 .and. abs(DtTotal-Dt) > DtSub*1.1) then
                 write(*,*) "Neutral Chemistry is too fast!!", &
                      DtSub, abs(DtTotal-Dt), iLon, iLat, iAlt

                 call write_sources_losses(iDtReducer)

                 call stop_gitm("Neutral Chemistry is too fast!!")
              endif

              nIters = nIters + 1

           enddo

           do iNeutral = 1, nSpeciesTotal
              NDensityS(iLon,iLat,iAlt,iNeutral,iBlock) = &
                   max(0.01,NDensityS(iLon,iLat,iAlt,iNeutral,iBlock))
           enddo

           Rho(iLon,iLat,iAlt,iBlock) = &
                sum( Mass(1:nSpecies) * &
                nDensityS(iLon,iLat,iAlt,1:nSpecies,iBlock))
           nDensity(iLon,iLat,iAlt,iBlock) = &
                sum(nDensityS(iLon,iLat,iAlt,1:nSpecies,iBlock))

        enddo
     enddo
  enddo

  if (iDebugLevel > 3) then
     do iNeutral = 1, nSpeciesTotal
        write(*,*) "====> Max (Neutral) : ", iNeutral, &
             maxval(nDensityS(1:nLons,1:nLats,(nAlts*4)/5,iNeutral,iBlock))
     enddo
  endif

  if (iDebugLevel > 0) &
       write(*,*) "===> Neutral Chemistry Average Dt for this timestep : ", &
       (Dt*nLats*nLons*nAlts)/nIters

  call end_timing("Neutral Chemistry")

contains

  subroutine write_sources_losses(iSpecies)

    integer, intent(in) :: iSpecies

    write(*,*) "Constituent : ", iSpecies, &
         nDensityS(iLon,iLat,iAlt,iSpecies,iBlock), &
         sources(iSpecies), losses(iSpecies)

    if (iSpecies == iO_) then

       !  write(*,*) "k1 : ", &
            !     RrTempDep(iLon,iLat,iAlt,iRrK1_) *   &
       !     IDensityS( iLon,iLat,iAlt,iO_4SP_,iBlock)   *  &
            !     NDensityS(iLon,iLat,iAlt,iO2_,iBlock)

       !   write(*,*) "k32 : ", &
            !        RrTempInd(iRrK32_) * &
       !        IDensityS(iLon,iLat,iAlt,iO2P_,iBlock) *  &
            !        NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock)

       !   write(*,*) "k7 : ", &
            !        RrTempInd(iRrK7_) * &
       !        IDensityS(iLon,iLat,iAlt,iNP_,iBlock)  *  &
            !        NDensityS(iLon,iLat,iAlt,iO2_,iBlock)

       write(*,*) "a1 : ", &
            RrTempDep(iLon,iLat,iAlt,iRrA1_) * &
            IDensityS(iLon,iLat,iAlt,iNOP_,iBlock)*  &
            IDensityS(iLon,iLat,iAlt,ie_,iBlock)
       !
       !   write(*,*) "beta1 : ", &
            !        RrTempDep(iLon,iLat,iAlt,iRrBeta1_) * &
       !        NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock) * &
            !        NDensityS(iLon,iLat,iAlt,iO2_,iBlock)
       !
       !   write(*,*) "beta3 : ", &
            !        RrTempDep(iLon,iLat,iAlt,iRrBeta3_) * &
       !        NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock)*&
            !        NDensityS(iLon,iLat,iAlt,iNO_,iBlock)
       !
       !   write(*,*) "beta6 : ", &
            !        RrTempInd(iRrBeta6_) * &
       !        NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock) *     &
            !        NDensityS(iLon,iLat,iAlt,iNO_,iBlock)
       !
       write(*,*) "beta8 : ", &
            RrTempDep(iLon,iLat,iAlt,iRrBeta8_) * &
            NDensityS(iLon,iLat,iAlt,iNO_,iBlock)
       !
       !  write(*,*) "beta2 : ", &
            !       RrTempInd(iRrBeta2_) * &
       !       NDensityS(iLon,iLat,iAlt,iO2_,iBlock)       *     &
            !       NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock)
       !
       write(*,*) "a2 : ", &
            1.15*RrTempDep(iLon,iLat,iAlt,iRrA2_) * &
            IDensityS(iLon,iLat,iAlt,ie_,iBlock) * &
            IDensityS(iLon,iLat,iAlt,iO2P_,iBlock)

       !   write(*,*) "j1 : ", &
            !        2*RrTempInd(iRrJ1_) * &
       !        NDensityS(iLon,iLat,iAlt,iO2_,iBlock)

       write(*,*) "Loss : k12 : ", &
            2*RrTempDep(iLon,iLat,iAlt,iRrK12_) * &
            NDensityS(iLon,iLat,iAlt,iO_,iBlock) * &
            (NDensityS(iLon,iLat,iAlt,iO_,iBlock) +  &
            NDensityS(iLon,iLat,iAlt,iO2_,iBlock)+  &
            NDensityS(iLon,iLat,iAlt,iN2_,iBlock)) * &
            NDensityS(iLon,iLat,iAlt,iO_,iBlock)
    endif

    if (iSpecies == iO2_) then

       write(*,*) "o2 source, k5 : ", &
            RrTempInd(iRrK5_) * &
            IDensityS(iLon,iLat,iAlt,iO2P_, iBlock) *  &
            NDensityS(iLon,iLat,iAlt,iNO_,iBlock)
       write(*,*) "o2 source, k12 : ", &
            RrTempDep(iLon,iLat,iAlt,iRrK12_) * &
            NDensityS(iLon,iLat,iAlt,iO_,iBlock)**2 * &
            (NDensityS(iLon,iLat,iAlt,iO_,iBlock)  +   &
            NDensityS(iLon,iLat,iAlt,iO2_,iBlock) +   &
            NDensityS(iLon,iLat,iAlt,iN2_,iBlock))
       write(*,*) "o2 loss, k1 : ", &
            RrTempDep(iLon,iLat,iAlt,iRrK1_) * &
            IDensityS(iLon,iLat,iAlt,iO_4SP_,iBlock)
       write(*,*) "o2 loss, k6 : ", &
            RrTempInd(iRrK6_) * &
            IDensityS(iLon,iLat,iAlt,iNP_, iBlock)
       write(*,*) "o2 loss, k7 : ", &
            RrTempInd(iRrK7_) * &
            IDensityS(iLon,iLat,iAlt,iNP_, iBlock)
       write(*,*) "o2 loss, b1 : ", &
            RrTempDep(iLon,iLat,iAlt,iRrBeta1_) * &
            NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock)
       write(*,*) "o2 loss, b2 : ", &
            RrTempInd(iRrBeta2_) * &
            NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock)
       write(*,*) "o2 loss, k13 : ", &
            RrTempDep(iLon,iLat,iAlt,iRrK13_) * &
            NDensityS(iLon,iLat,iAlt,iO_,iBlock) *    &
            (NDensityS(iLon,iLat,iAlt,iO_,iBlock) +    &
            NDensityS(iLon,iLat,iAlt,iO2_,iBlock))
       write(*,*) "o2 loss, j1 : ", &
            RrTempInd(iRrJ1_) 

    endif

    if (iSpecies == iN2_) then

       write(*,*) "Sources(N2) : ", &
            RrTempDep(iLon,iLat,iAlt,iRrBeta3_) * &
            NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock) * &
            NDensityS(iLon,iLat,iAlt,iNO_,iBlock), &
            RrTempInd(iRrBeta6_) *  &
            NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock) * &
            NDensityS(iLon,iLat,iAlt,iNO_,iBlock)

       write(*,*) "Losses(N2) : ", &
            NDensityS(iLon,iLat,iAlt,iN2_,iBlock) * &
            RrTempDep(iLon,iLat,iAlt,iRrK2_) * &
            IDensityS(iLon,iLat,iAlt,iO_4SP_,iBlock),  & 
            NDensityS(iLon,iLat,iAlt,iN2_,iBlock) * &
            EuvIonRateS(iLon,iLat,iAlt,iN2P_,iBlock)

    endif

    if (iSpecies == iNO_) then 
       write(*,*) "NO Losses : ", &
            NDensityS(iLon,iLat,iAlt,iNO_,iBlock) * &
            RrTempInd(iRrK5_)* &
            IDensityS(iLon,iLat,iAlt,iO2P_, iBlock),          &
            NDensityS(iLon,iLat,iAlt,iNO_,iBlock) * &
            RrTempDep(iLon,iLat,iAlt,iRrBeta3_)*&
            NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock),&
            NDensityS(iLon,iLat,iAlt,iNO_,iBlock) * &
            RrTempInd(iRrBeta6_)*&
            NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock),      &
            NDensityS(iLon,iLat,iAlt,iNO_,iBlock) * &
            RrTempDep(iLon,iLat,iAlt,iRrBeta8_), &
            NDensityS(iLon,iLat,iAlt,iNO_,iBlock) * &
            RrTempDep(iLon,iLat,iAlt,iRrBeta9N_)
    endif

    if (iSpecies == iN_2D_) then

       write(*,*) "N(2D) Sources : ", &
            RrTempDep(iLon,iLat,iAlt,iRrK3_) * &
            IDensityS(iLon,iLat,iAlt,iN2P_,iBlock) *    &
            NDensityS(iLon,iLat,iAlt,iO_,iBlock),             &
            !
            0.8*RrTempDep(iLon,iLat,iAlt,iRrA1_) * &
            IDensityS(iLon,iLat,iAlt,iNOP_,iBlock) *      &
            IDensityS(iLon,iLat,iAlt,ie_,iBlock), &
            !
            0.8*RrTempDep(iLon,iLat,iAlt,iRrA1_), &
            !
       IDensityS(iLon,iLat,iAlt,iNOP_,iBlock),      &
            !
       IDensityS(iLon,iLat,iAlt,ie_,iBlock), &
            !
       0.9*RrTempDep(iLon,iLat,iAlt,iRrA3_) * &
            IDensityS(iLon,iLat,iAlt,iN2P_,iBlock) *      &
            IDensityS(iLon,iLat,iAlt,ie_,iBlock)

    endif

  end subroutine write_sources_losses

end subroutine calc_neutral_density



