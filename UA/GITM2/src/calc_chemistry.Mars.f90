
subroutine calc_chemistry(iBlock)

  use ModSizeGitm
  use ModGITM
  use ModPlanet
  use ModRates
  use ModEUV
  use ModSources
  use ModInputs, only: iDebugLevel, UseIonChemistry, UseNeutralChemistry,f107
  use ModConstants

  implicit none

  integer, intent(in) :: iBlock

  real :: IonSources(nIons), NeutralSources(nSpeciesTotal)
  real :: IonLosses(nIons), NeutralLosses(nSpeciesTotal)
  real :: DtSub, DtOld, DtTotal, DtMin, DtAve, Source, Reaction, tr, tr3, rr
  real :: te3, ti, tn, dtsubtmp, losstmp, dentmp
  real :: Ions(nIons), Neutrals(nSpeciesTotal)

  integer :: iLon, iLat, iAlt, iIon, nIters, iDtReducer, iNeutral

  real :: ChemicalHeatingSub
  real :: Emission(nEmissions), EmissionTotal(nEmissions)
  real, dimension(nLons,nLats,nAlts) :: &
       tr3d, tr33d, tr3m0443d, tr3m083d, tr3m043d, &
       te3m0393d, te3m0853d, te33d, te3m053d
  real :: te3m05, tr3m044, tr3m04, tr3m08, te3m039, te3m085, rr_opn2

  logical :: UseNeutralConstituent(nSpeciesTotal)
  logical :: UseIonConstituent(nIons)

  UseNeutralConstituent = .true.
  UseIonConstituent     = .true.

  DtMin = Dt

  if (.not.UseIonChemistry) return

  call report("Chemistry",2)
  call start_timing("calc_chemistry")

  DtAve = 0.0

  nIters=0

  tr3d = (iTemperature(1:nLons,1:nLats,1:nAlts,iBlock) &
         + Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*&
         TempUnit(1:nLons,1:nLats,1:nAlts)) / 2.0

  tr33d = tr3d/300.0
  te33d = eTemperature(1:nLons,1:nLats,1:nAlts,iBlock)/300.0

  do iLon = 1, nLons
     do iLat = 1, nLats
        do iAlt = 1, nAlts

           tr  = tr3d(iLon,iLat,iAlt)
           tr3 = tr33d(iLon,iLat,iAlt)
           te3 = te33d(iLon,iLat,iAlt)
           ti = iTemperature(iLon,iLat,iAlt,iBlock)
           tn = Temperature(iLon,iLat,iAlt,iBlock)*&
                TempUnit(iLon,iLat,iAlt)

           DtTotal = 0.0
           EmissionTotal = 0.0

           Ions = IDensityS(iLon,iLat,iAlt,:,iBlock)
           Neutrals = NDensityS(iLon,iLat,iAlt,:,iBlock)
 
           do while (DtTotal < Dt)

              ChemicalHeatingSub = 0.0
              Emission = 0.0

              DtSub = Dt - DtTotal

              IonSources = 0.0
              NeutralSources = 0.0
              IonLosses  = 0.0
              NeutralLosses = 0.0

              ! Solar EUV

              ! ----------------------------------------------------------
              ! O+
              ! ----------------------------------------------------------

              ! Solar EUV

              Reaction = EuvIonRateS(iLon,iLat,iAlt,iOP_,iBlock) * &
                   Neutrals(iO_)

              IonSources(iOP_)   = IonSources(iOP_)   + Reaction
              NeutralLosses(iO_) = NeutralLosses(iO_) + Reaction

!!!              !O2+ + N -> NO+ + O
!!!
!!!              Reaction = &
!!!                   RrTempInd(iRrR1_) * &
!!!                   Ions(iO2P_) * &
!!!                   Neutrals(iN_)
!!!
!!!              IonSources(iNOP_)   = IonSources(iNOP_)   + Reaction
!!!              NeutralSources(iO_) = NeutralSources(iO_) + Reaction
!!!              IonLosses(iO2P_)    = IonLosses(iO2P_)  + Reaction
!!!              NeutralLosses(iN_)  = NeutralLosses(iN_) + Reaction
!!!
!!!              !O2+ + NO -> NO+ + O2
!!!
!!!              Reaction = &
!!!                   RrTempInd(iRrR2_) * &
!!!                   Ions(iO2P_) * &
!!!                   Neutrals(iNO_)
!!!
!!!              IonSources(iNOP_)   = IonSources(iNOP_)   + Reaction
!!!              NeutralSources(iO2_) = NeutralSources(iO2_) + Reaction
!!!              IonLosses(iO2P_)    = IonLosses(iO2P_)  + Reaction
!!!              NeutralLosses(iNO_)  = NeutralLosses(iNO_) + Reaction
!!!
!!!!!              !N2+ + O -> NO+ + N
!!!!!
!!!!!              Reaction = &
!!!!!                   RrTempInd(iRrR3_) * &
!!!!!                   Ions(iN2P_) * &
!!!!!                   Neutrals(iO_)
!!!!!
!!!!!              IonSources(iNOP_)   = IonSources(iNOP_)   + Reaction
!!!!!              NeutralSources(iN_) = NeutralSources(iN_) + Reaction
!!!!!              IonLosses(iN2P_)    = IonLosses(iN2P_)  + Reaction
!!!!!              NeutralLosses(iO_)  = NeutralLosses(iO_) + Reaction
!!!
!!!              !CO2+ + N -> NO+ + CO
!!!
!!!              Reaction = &
!!!                   RrTempInd(iRrR4_) * &
!!!                   Ions(iCO2P_) * &
!!!                   Neutrals(iN_)
!!!
!!!              IonSources(iNOP_)   = IonSources(iNOP_)   + Reaction
!!!              NeutralSources(iCO_) = NeutralSources(iCO_) + Reaction
!!!              IonLosses(iCO2P_)    = IonLosses(iCO2P_)  + Reaction
!!!              NeutralLosses(iN_)  = NeutralLosses(iN_) + Reaction
!!!
!!!              !CO2+ + NO -> NO+ + CO2
!!!
!!!              Reaction = &
!!!                   RrTempInd(iRrR5_) * &
!!!                   Ions(iCO2P_) * &
!!!                   Neutrals(iNO_)
!!!
!!!              IonSources(iNOP_)     = IonSources(iNOP_)   + Reaction
!!!              NeutralSources(iCO2_) = NeutralSources(iCO2_) + Reaction
!!!              IonLosses(iCO2P_)     = IonLosses(iCO2P_)  + Reaction
!!!              NeutralLosses(iCO2_)  = NeutralLosses(iCO2_) + Reaction
!!!
!!!!!              !CO+ + NO -> NO+ + CO
!!!!!
!!!!!              Reaction = &
!!!!!                   RrTempInd(iRrR7_) * &
!!!!!                   Ions(iCOP_) * &
!!!!!                   Neutrals(iNO_)
!!!!!
!!!!!              IonSources(iNOP_)     = IonSources(iNOP_)   + Reaction
!!!!!              NeutralSources(iCO_) = NeutralSources(iCO_) + Reaction
!!!!!              IonLosses(iCOP_)     = IonLosses(iCOP_)  + Reaction
!!!!!              NeutralLosses(iNO_)  = NeutralLosses(iNO_) + Reaction
!!!
!!!!!              ChemicalHeatingSub = &
!!!!!                   ChemicalHeatingSub + Reaction * 1.33

              ! ----------------------------------------------------------
              ! CO2+
              ! ----------------------------------------------------------

              ! Solar EUV

              Reaction = EuvIonRateS(iLon,iLat,iAlt,iCO2P_,iBlock) * &
                   Neutrals(iCO2_)

              IonSources(iCO2P_)   = IonSources(iCO2P_)   + Reaction
              NeutralLosses(iCO2_) = NeutralLosses(iCO2_) + Reaction



!!              ! Aurora
!!
!!              Reaction = AuroralIonRateS(iLon,iLat,iAlt,iN2_, iBlock) + &
!!                   IonPrecipIonRateS(iLon,iLat,iAlt,iN2_, iBlock)
!!
!!              IonSources(iN2P_)   = IonSources(iN2P_) + Reaction
!!              NeutralLosses(iN2_) = NeutralLosses(iN2_) + Reaction
!!
              


              do iIon = 1, nIons-1

                 if (UseIonChemistry .and. &
                      UseIonConstituent(iIon)) then

                    if (IonLosses(iIon) > IonSources(iIon) .and. &
                         IonSources(iIon)>0.0) then
                       if (Ions(iIon) > 0.01) then
                          dtOld = DtSub
                          dtSub = &
                               min(0.25 * &
                               (IonSources(iIon) + &
                               Ions(iIon))/ &
                               (abs(IonLosses(iIon))+0.1), DtOld)
                          if (DtSub < DtOld) iDtReducer = iIon
                       else
                          IonSources(iIon) = 0.0
                          IonLosses(iIon) = 0.0
                       endif
                    endif
                 else
                    IonSources(iIon) = 0.0
                    IonLosses(iIon) = 0.0
                 endif
              enddo

              do iNeutral = 1, nSpeciesTotal

                 if (UseNeutralChemistry .and. &
                      UseNeutralConstituent(iNeutral)) then

                    if (NeutralLosses(iNeutral) > &
                         NeutralSources(iNeutral)) then
                       if (Neutrals(iNeutral)>0.01) then
                          dtOld = DtSub
                          dtSub = &
                               min(0.25 * &
                               Neutrals(iNeutral)/ &
                               (abs(NeutralSources(iNeutral) - &
                               NeutralLosses(iNeutral))+0.1), DtOld)
                          if (DtSub < DtOld) iDtReducer = nIons + iNeutral
                       else
                          NeutralSources(iNeutral) = 0.0
                          NeutralLosses(iNeutral) = 0.0
                       endif
                    endif
                 else
                    NeutralSources(iNeutral) = 0.0
                    NeutralLosses(iNeutral) = 0.0
                 endif
              enddo

              Ions(nIons) = 0.0
              do iIon = 1, nIons-1
                 Ions(iIon) = Ions(iIon) + &
                      (IonSources(iIon) - IonLosses(iIon)) * DtSub
                 Ions(iIon) = max(0.01,Ions(iIon))
                 
                 ! sum for e-
                 Ions(nIons) = Ions(nIons) + Ions(iIon)

                 if (Ions(iIon) < 0.0) then
                    write(*,*) "Negative Ion Density : ", &
                         iIon, iLon, iLat, iAlt, &
                         Ions(iIon), &
                         IonSources(iIon), IonLosses(iIon)
                 endif
              enddo

              do iNeutral = 1, nSpeciesTotal
                 Neutrals(iNeutral) = &
                      Neutrals(iNeutral) + &
                      (NeutralSources(iNeutral) - NeutralLosses(iNeutral)) * &
                      DtSub
              enddo

              ChemicalHeatingRate(iLon,iLat,iAlt) = &
                   ChemicalHeatingRate(iLon,iLat,iAlt) + &
                   ChemicalHeatingSub * DtSub

              EmissionTotal = EmissionTotal + Emission(:)*DtSub

              DtTotal = DtTotal + DtSub

              if (DtSub < DtMin) DtMin = DtSub

              if (DtSub < 1.0e-9 .and. abs(DtTotal-Dt) > DtSub) then
                 write(*,*) "Chemistry is too fast!!", DtSub
                 if (iDtReducer < nIons) then
                    write(*,*) "Ion Constituent : ", &
                         iDtReducer, iLon, iLat, iAlt,&
                         Ions(iDtReducer), &
                         IonSources(iDtReducer), IonLosses(iDtReducer)
                 else
                    iDtReducer = iDtReducer - nIons
                    write(*,*) "Neutral Constituent : ", &
                         iDtReducer, iLon, iLat, iAlt,&
                         Ions(iDtReducer), &
                         IonSources(iDtReducer), IonLosses(iDtReducer)
                 endif

                 call stop_gitm("Ion Chemistry is too fast!!")
              endif

              nIters = nIters + 1

           enddo

           IDensityS(iLon,iLat,iAlt,:,iBlock) = Ions
           NDensityS(iLon,iLat,iAlt,:,iBlock) = Neutrals

           Emissions(iLon, iLat, iAlt, :, iBlock) =  &
                Emissions(iLon, iLat, iAlt, :, iBlock) + EmissionTotal

        enddo
     enddo
  enddo

  if (iDebugLevel > 3) then
     do iIon = 1, nIons
        write(*,*) "====> Max Ion Density: ", iIon, &
             maxval(IDensityS(1:nLons,1:nLats,(nAlts*4)/5,iIon,iBlock))
     enddo
  endif

  if (iDebugLevel > 2) &
       write(*,*) "===> Average Dt for this timestep : ", &
       (Dt*nLats*nLons*nAlts)/nIters

  call end_timing("calc_chemistry")

end subroutine calc_chemistry
