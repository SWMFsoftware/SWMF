
subroutine calc_chemistry(iBlock)

  use ModSize
  use ModGITM
  use ModPlanet
  use ModRates
  use ModEUV
  use ModSources
  use ModInputs, only: iDebugLevel, UseIonChemistry, UseNeutralChemistry
  use ModConstants

  implicit none

  integer, intent(in) :: iBlock

  real :: IonSources(nIons), NeutralSources(nSpeciesTotal)
  real :: IonLosses(nIons), NeutralLosses(nSpeciesTotal)
  real :: DtSub, DtOld, DtTotal, DtMin, DtAve, Source, Reaction, tr, tr3, rr
  real :: te3, ti, tn, dtsubtmp, losstmp, dentmp

  integer :: iLon, iLat, iAlt, iIon, nIters, iDtReducer, iNeutral

  real :: ChemicalHeatingSub
  real :: Emission(nEmissions)
  real, dimension(nLons,nLats,nAlts) :: &
       tr3d, tr33d, tr3m0443d, tr3m083d, tr3m043d, &
       te3m0393d, te3m0853d, te33d, te3m053d
  real :: te3m05, tr3m044, tr3m04, tr3m08, te3m039, te3m085, rr_opn2

  logical :: UseNeutralConstituent(nSpeciesTotal)
  logical :: UseIonConstituent(nIons)

  UseNeutralConstituent = .true.
  UseIonConstituent     = .true.

!  UseIonConstituent(iO_4SP_) = .false.
  UseIonConstituent(iO_2PP_) = .false.
  
!  UseNeutralConstituent(iN_4S_) = .false.
!  UseNeutralConstituent(iN_2D_) = .false.
!  UseNeutralConstituent(iO2_) = .false.

  DtMin = Dt

  if (.not.UseIonChemistry) return

  call report("Chemistry",2)
  call start_timing("calc_chemistry")

  DtAve = 0.0

  nIters=0

!  AuroralIonRateS = 0.0

  tr3d = (iTemperature(1:nLons,1:nLats,1:nAlts,iBlock) &
         + Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*TempUnit) / 2.0

  tr33d = tr3d/300.0
  te33d = eTemperature(1:nLons,1:nLats,1:nAlts,iBlock)/300.0

  te3m053d   = te33d**(-0.5)
  te3m0393d  = te33d**(-0.39)
  te3m0853d  = te33d**(-0.85)

  tr3m0443d = tr33d**(-0.44)
  tr3m083d  = tr33d**(-0.8)
  tr3m043d  = tr33d**(-0.4)

  do iLon = 1, nLons
     do iLat = 1, nLats
        do iAlt = 1, nAlts

           tr  = tr3d(iLon,iLat,iAlt)
           tr3 = tr33d(iLon,iLat,iAlt)
           te3 = te33d(iLon,iLat,iAlt)
           te3m05  = te3m053d(iLon,iLat,iAlt)
           te3m085 = te3m0853d(iLon,iLat,iAlt)
           te3m039 = te3m0393d(iLon,iLat,iAlt)
           tr3m044 = tr3m0443d(iLon,iLat,iAlt)
           tr3m04  = tr3m043d(iLon,iLat,iAlt)

           tr3m08  = tr3m083d(iLon,iLat,iAlt)
           ti = iTemperature(iLon,iLat,iAlt,iBlock)
           tn = Temperature(iLon,iLat,iAlt,iBlock)*TempUnit

           rr_opn2 = min(5.0e-19,4.5e-20*tr3**2)

           DtTotal = 0.0

           do while (DtTotal < Dt)

              ChemicalHeatingSub = 0.0
              Emission = 0.0

              DtSub = Dt - DtTotal

              IonSources = 0.0
              NeutralSources = 0.0
              IonLosses  = 0.0
              NeutralLosses  = 0.0

              ! ----------------------------------------------------------
              ! N2+
              ! ----------------------------------------------------------

              ! Solar EUV

              Reaction = EuvIonRateS(iLon,iLat,iAlt,iN2P_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iN2_,iBlock)

              IonSources(iN2P_)   = IonSources(iN2P_)   + Reaction
              NeutralLosses(iN2_) = NeutralLosses(iN2_) + Reaction

              ! Aurora

              Reaction = AuroralIonRateS(iLon,iLat,iAlt,iN2_, iBlock)

              IonSources(iN2P_)   = IonSources(iN2P_) + Reaction
              NeutralLosses(iN2_) = NeutralLosses(iN2_) + Reaction

              ! O+(2D) + N2 -> N2+ + O + 1.33 eV

              Reaction = &
                   8.0e-16 * &
                   IDensityS(iLon,iLat,iAlt,iO_2DP_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iN2_,iBlock)

              IonSources(iN2P_)    = IonSources(iN2P_)   + Reaction
              NeutralSources(iO_)  = NeutralSources(iO_) + Reaction
              IonLosses(iO_2DP_)   = IonLosses(iO_2DP_)  + Reaction
              NeutralLosses(iN2_)  = NeutralLosses(iN2_) + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + Reaction * 1.33

              ! O+(2P) + N2 -> N2+ + O + 3.02 eV

              Reaction = &
                   4.8e-16 * &
                   IDensityS(iLon,iLat,iAlt,iO_2PP_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iN2_,iBlock)

              IonSources(iN2P_)    = IonSources(iN2P_)   + Reaction
              NeutralSources(iO_)  = NeutralSources(iO_) + Reaction
              IonLosses(iO_2PP_)   = IonLosses(iO_2PP_)  + Reaction
              NeutralLosses(iN2_)  = NeutralLosses(iN2_) + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + Reaction * 3.02

              ! N2+ + O2 -> O2+ + N2 + 3.02 eV

              rr = 5.0e-17 * tr3m08

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iN2P_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)

              IonSources(iO2P_)    = IonSources(iO2P_)   + Reaction
              NeutralSources(iN2_) = NeutralSources(iN2_) + Reaction
              IonLosses(iN2P_)     = IonLosses(iN2P_)  + Reaction
              NeutralLosses(iO2_)  = NeutralLosses(iO2_) + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + Reaction * 3.53

              ! N2+ + O2 -> NO+ + N(2D) + 0.70 eV
              !          -> NO+ + N(4S) + 3.08 eV

              rr = 1.4e-16 * tr3m044

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iN2P_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iO_,iBlock)

              IonSources(iNOP_)      = IonSources(iNOP_)      + Reaction
              NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + 0.5*Reaction
              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + 0.5*Reaction
              IonLosses(iN2P_)       = IonLosses(iN2P_)       + Reaction
              NeutralLosses(iO_)     = NeutralLosses(iO_)     + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   0.5 * Reaction * 0.70 + &
                   0.5 * Reaction * 3.08

              ! N2+ + e -> 2 N(2D) + 1.04 eV

              rr = 1.8e-13 * te3m039

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iN2P_,iBlock) * &
                   IDensityS(iLon,iLat,iAlt,ie_,iBlock)

              NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + 2*Reaction
              IonLosses(iN2P_)       = IonLosses(iN2P_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 1.04

              ! N2+ + O -> O+(4S) + N2 + 1.96 eV

              rr = 1.4e-16 * tr3m044

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iN2P_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iO_,iBlock)

              NeutralSources(iN2_) = NeutralSources(iN2_) + Reaction
              IonSources(iO_4SP_)  = IonSources(iO_4SP_)  + Reaction
              NeutralLosses(iO_)   = NeutralLosses(iO_) + Reaction
              IonLosses(iN2P_)     = IonLosses(iN2P_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 1.96

              ! N2+ + NO -> NO+ + N2 + 6.33 eV

              rr = 3.3e-16

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iN2P_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock)

              NeutralSources(iN2_) = NeutralSources(iN2_) + Reaction
              IonSources(iNOP_)    = IonSources(iNOP_)     + Reaction
              NeutralLosses(iNO_)  = NeutralLosses(iNO_)  + Reaction
              IonLosses(iN2P_)     = IonLosses(iN2P_)    + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 6.33

              ! ----------------------------------------------------------
              ! O2+
              ! ----------------------------------------------------------

              ! -----------
              ! Solar EUV
              ! -----------

              Reaction = EuvIonRateS(iLon,iLat,iAlt,iO2P_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)

              IonSources(iO2P_)   = IonSources(iO2P_)   + Reaction
              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction

              ! -----------
              ! Aurora
              ! -----------

              Reaction = AuroralIonRateS(iLon,iLat,iAlt,iO2_, iBlock)

              IonSources(iO2P_)   = IonSources(iO2P_)   + Reaction
              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction

              ! -----------
              ! O+(4S) + O2 -> O2+ + O + 1.55 eV
              ! -----------

              rr = 2.0e-17 * tr3m04

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iO_4SP_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)

              NeutralSources(iO_) = NeutralSources(iO_) + Reaction
              IonSources(iO2P_)   = IonSources(iO2P_)   + Reaction
              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction/1.0
              IonLosses(iO_4SP_)  = IonLosses(iO_4SP_) + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 1.55

              ! -----------
              ! O+(2D) + O2 -> O2+ + O + 4.865 eV
              ! -----------

              rr = 7.0e-16

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iO_2DP_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)

              NeutralSources(iO_) = NeutralSources(iO_) + Reaction
              IonSources(iO2P_)   = IonSources(iO2P_)   + Reaction
              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction/1.0
              IonLosses(iO_2DP_)  = IonLosses(iO_2DP_) + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 4.865

              ! -----------
              ! N+ + O2 -> O2+ + N(4S) + 2.486 eV
              ! -----------

              rr = 1.1e-16

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iNP_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)

              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              IonSources(iO2P_)      = IonSources(iO2P_)      + Reaction
              NeutralLosses(iO2_)    = NeutralLosses(iO2_)    + Reaction
              IonLosses(iNP_)        = IonLosses(iNP_)       + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 2.486

              ! -----------
              ! N+ + O2 -> O2+ + N(2D) + 0.1 eV
              ! -----------

              rr = 2.0e-16

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iNP_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)

              NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + Reaction
              IonSources(iO2P_)      = IonSources(iO2P_)      + Reaction
              NeutralLosses(iO2_)    = NeutralLosses(iO2_)    + Reaction
              IonLosses(iNP_)        = IonLosses(iNP_)       + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 0.1

              ! -----------
              ! O2+ + e -> 2 O + 5.0 eV
              ! -----------

              rr = 1.9e-13 * te3m05

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iO2P_,iBlock) * &
                   IDensityS(iLon,iLat,iAlt,ie_,iBlock)

              NeutralSources(iO_) = NeutralSources(iO_) + 2*Reaction
              IonLosses(iO2P_)    = IonLosses(iO2P_)   + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 5.0

              ! -----------
              ! O2+ + N(2D) -> N+ + O2 + 0.0 eV
              ! -----------

              rr = 2.5e-16

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iO2P_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock)

              NeutralSources(iO2_)  = NeutralSources(iO2_)  + Reaction
              IonSources(iNP_)      = IonSources(iNP_)      + Reaction
              NeutralLosses(iN_2D_) = NeutralLosses(iN_2D_) + Reaction
              IonLosses(iO2P_)      = IonLosses(iO2P_)      + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 0.0

              ! -----------
              ! O2+ + N(4S) -> NO+ + O + 4.213 eV
              ! -----------

              rr = 1.8e-16

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iO2P_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock)

              NeutralSources(iO_)   = NeutralSources(iO_)   + Reaction
              IonSources(iNOP_)     = IonSources(iNOP_)     + Reaction
              NeutralLosses(iN_4S_) = NeutralLosses(iN_4S_) + Reaction
              IonLosses(iO2P_)      = IonLosses(iO2P_)      + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 4.213

              ! -----------
              ! O2+ + NO -> NO+ + O2 + 2.813 eV
              ! -----------

              rr = 4.4e-16

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iO2P_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock)

              NeutralSources(iO2_) = NeutralSources(iO2_) + Reaction
              IonSources(iNOP_)    = IonSources(iNOP_)    + Reaction
              NeutralLosses(iNO_)  = NeutralLosses(iNO_)  + Reaction
              IonLosses(iO2P_)     = IonLosses(iO2P_)     + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 2.813

              ! -----------
              ! O2+ + N2 -> NO+ + NO + 0.9333 eV
              ! -----------

              rr = 5.0e-22

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iO2P_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iN2_,iBlock)

              NeutralSources(iNO_) = NeutralSources(iNO_) + Reaction
              IonSources(iNOP_)    = IonSources(iNOP_)    + Reaction
              NeutralLosses(iN2_)  = NeutralLosses(iN2_)  + Reaction
              IonLosses(iO2P_)     = IonLosses(iO2P_)     + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 0.9333

              ! ----------------------------------------------------------
              ! O(4S)+
              ! ----------------------------------------------------------

              ! Solar EUV

              Reaction = EuvIonRateS(iLon,iLat,iAlt,iO_4SP_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iO_,iBlock)

              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              NeutralLosses(iO_)  = NeutralLosses(iO_)  + Reaction

              ! Aurora

              Reaction = AuroralIonRateS(iLon,iLat,iAlt,iO_, iBlock)

              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              NeutralLosses(iO_)  = NeutralLosses(iO_)  + Reaction

              ! -----------
              ! O+(2D) + O -> O+(4S) + O + 3.31 eV
              ! -----------

              rr = 1.0e-17

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iO_2DP_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iO_,iBlock)

              ! We create and loose the same amount of O
              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2DP_)  = IonLosses(iO_2DP_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 3.31

              ! -----------
              ! O+(2D) + e -> O+(4S) + e + 3.31 eV
              ! -----------

              rr = 7.8e-14 * te3m05

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iO_2DP_,iBlock) * &
                   IDensityS(iLon,iLat,iAlt,ie_,iBlock)

              ! We create and loose the same amount of e
              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2DP_)  = IonLosses(iO_2DP_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 3.31

              ! -----------
              ! O+(2D) + N2 -> O+(4S) + N2 + 3.31 eV
              ! -----------

              rr = 8.0e-16

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iO_2DP_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iN2_,iBlock)

              ! We create and loose the same amount of N2
              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2DP_)  = IonLosses(iO_2DP_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 3.31

              ! -----------
              ! O+(2P) + O -> O+(4S) + O + 5.0 eV
              ! -----------

              rr = 5.2e-17

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iO_2PP_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iO_,iBlock)

              ! We create and loose the same amount of O
              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2PP_)  = IonLosses(iO_2PP_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 5.0

              ! -----------
              ! O+(2P) + e -> O+(4S) + e + 5.0 eV
              ! -----------

              rr = 4.0e-14 * te3m05

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iO_2PP_,iBlock) * &
                   IDensityS(iLon,iLat,iAlt,ie_,iBlock)

              ! We create and loose the same amount of e
              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2PP_)  = IonLosses(iO_2PP_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 5.0

              ! -----------
              ! O+(2P) -> O+(4S) + 2470A
              ! -----------

              rr = 0.047

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iO_2PP_,iBlock)

              ! We create and loose the same amount of e
              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2PP_)  = IonLosses(iO_2PP_)  + Reaction

              Emission(iE2470_) = Emission(iE2470_) + Reaction

!              ! -----------
!              ! H+ + O -> O+(4S) + H + 0.0 eV
!              ! -----------
!
!              rr = 6.0e-10/1.0e6 * (8.0/9.0) * sqrt( &
!                   (ti + tn/4) * (tn + ti/16))
!
!              Reaction = &
!                   rr * &
!                   IDensityS(iLon,iLat,iAlt,iHP_,iBlock) * &
!                   NDensityS(iLon,iLat,iAlt,iO_,iBlock)
!
!              NeutralSources(iH_) = NeutralSources(iH_) + Reaction
!              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
!              NeutralLosses(iO_)  = NeutralLosses(iO_)  + Reaction
!              IonLosses(iHP_)     = IonLosses(iHP_)     + Reaction
!
!              ChemicalHeatingSub = &
!                   ChemicalHeatingSub + &
!                   Reaction * 0.0

              ! -----------
              ! N+ + O2 -> O+(4S) + NO + 2.31 eV
              ! -----------

              rr = 3.0e-17

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iNP_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)

              NeutralSources(iNO_) = NeutralSources(iNO_) + Reaction
              IonSources(iO_4SP_)  = IonSources(iO_4SP_)  + Reaction
              NeutralLosses(iO2_)  = NeutralLosses(iO2_)  + Reaction
              IonLosses(iNP_)      = IonLosses(iNP_)      + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 2.31

              ! -----------
              ! O+(4S) + N2 -> NO+ + N(4S) + 1.10 eV
              ! -----------

              rr = rr_opn2

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iO_4SP_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iN2_,iBlock)

              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              IonSources(iNOP_)      = IonSources(iNOP_)      + Reaction
              NeutralLosses(iN2_)    = NeutralLosses(iN2_)    + Reaction
              IonLosses(iO_4SP_)     = IonLosses(iO_4SP_)     + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 1.10

              ! -----------
              ! O+(4S) + NO -> NO+ + O + 4.36 eV
              ! -----------

              rr = 8.0e-19

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iO_4SP_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock)

              NeutralSources(iO_) = NeutralSources(iO_) + Reaction
              IonSources(iNOP_)   = IonSources(iNOP_)   + Reaction
              NeutralLosses(iNO_) = NeutralLosses(iNO_) + Reaction
              IonLosses(iO_4SP_)  = IonLosses(iO_4SP_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 4.36

!              ! -----------
!              ! O+(4S) + H -> H+ + O + 0.0 eV
!              ! -----------
!
!              rr = 6.0e-10/1.0e6
!
!              Reaction = &
!                   rr * &
!                   IDensityS(iLon,iLat,iAlt,iO_4SP_,iBlock) * &
!                   NDensityS(iLon,iLat,iAlt,iH_,iBlock)
!
!              NeutralSources(iO_) = NeutralSources(iO_) + Reaction
!              IonSources(iHP_)    = IonSources(iHP_)    + Reaction
!              NeutralLosses(iH_)  = NeutralLosses(iH_)  + Reaction
!              IonLosses(iO_4SP_)  = IonLosses(iO_4SP_)  + Reaction
!
!              ChemicalHeatingSub = &
!                   ChemicalHeatingSub + &
!                   Reaction * 0.0

              ! -----------
              ! O+(4S) + N(2D) -> N+ + O + 1.45 eV
              ! -----------

              rr = 1.3e-16

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iO_4SP_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock)

              NeutralSources(iO_)   = NeutralSources(iO_)   + Reaction
              IonSources(iNP_)      = IonSources(iNP_)      + Reaction
              NeutralLosses(iN_2D_) = NeutralLosses(iN_2D_) + Reaction
              IonLosses(iO_4SP_)    = IonLosses(iO_4SP_)    + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 1.45

              ! ----------------------------------------------------------
              ! O(2D)+
              ! ----------------------------------------------------------

              ! Solar EUV

              Reaction = EuvIonRateS(iLon,iLat,iAlt,iO_2DP_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iO_,iBlock)

              IonSources(iO_2DP_) = IonSources(iO_2DP_) + Reaction
              NeutralLosses(iO_)  = NeutralLosses(iO_)  + Reaction

              ! Aurora

              !Reaction = AuroralIonRateS(iLon,iLat,iAlt,iO_, iBlock)

              !IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              !NeutralLosses(iO_)  = NeutralLosses(iO_)  + Reaction

              ! -----------
              ! O+(2P) + e -> O+(2D) + e + 1.69 eV
              ! -----------

              rr = 1.5e-13 * te3m05

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iO_2PP_,iBlock) * &
                   IDensityS(iLon,iLat,iAlt,ie_,iBlock)

              ! We create and loose the same amount of e
              IonSources(iO_2DP_) = IonSources(iO_2DP_) + Reaction
              IonLosses(iO_2PP_)  = IonLosses(iO_2PP_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 1.69

              ! -----------
              ! O+(2P) -> O+(2D) + 7320A
              ! -----------

              rr = 0.171

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iO_2PP_,iBlock)

              IonSources(iO_2DP_) = IonSources(iO_2DP_) + Reaction
              IonLosses(iO_2PP_)  = IonLosses(iO_2PP_)  + Reaction

              Emission(iE7320_) = Emission(iE7320_) + Reaction

              ! -----------
              ! O+(2D) -> O+(4S) + 3726A
              ! -----------

              rr = 7.7e-5

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iO_2DP_,iBlock)

              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2DP_)  = IonLosses(iO_2DP_)  + Reaction

              Emission(iE3726_) = Emission(iE3726_) + Reaction

              ! ----------------------------------------------------------
              ! O(2P)+
              ! ----------------------------------------------------------

              ! Solar EUV

              Reaction = EuvIonRateS(iLon,iLat,iAlt,iO_2PP_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iO_,iBlock)

              IonSources(iO_2PP_) = IonSources(iO_2PP_) + Reaction
              NeutralLosses(iO_)  = NeutralLosses(iO_)  + Reaction

              ! Aurora

              !Reaction = AuroralIonRateS(iLon,iLat,iAlt,iO_, iBlock)

              !IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              !NeutralLosses(iO_)  = NeutralLosses(iO_)  + Reaction

              ! -----------
              ! O+(2P) + N2 -> N+ + NO + 0.70 eV
              ! -----------

              rr = 1.0e-16

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iO_2PP_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iN2_,iBlock)

              NeutralSources(iNO_) = NeutralSources(iNO_) + Reaction
              IonSources(iNP_)     = IonSources(iNP_)     + Reaction
              NeutralLosses(iN2_)  = NeutralLosses(iN2_)  + Reaction
              IonLosses(iO_2PP_)   = IonLosses(iO_2PP_)   + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 0.70

              ! ----------------------------------------------------------
              ! N+
              ! ----------------------------------------------------------

              ! -----------
              ! O2+ + N(2D) -> N+ + O2 + 0.0 eV
              ! -----------

              rr = 2.5e-16

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iO2P_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock)

              NeutralSources(iO2_)  = NeutralSources(iO2_)  + Reaction
              IonSources(iNP_)      = IonSources(iNP_)      + Reaction
              NeutralLosses(iN_2D_) = NeutralLosses(iN_2D_) + Reaction
              IonLosses(iO2P_)      = IonLosses(iO2P_)      + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 0.0

!              ! -----------
!              ! He+ + N2 -> N+ + N + He + 0.28 eV
!              ! -----------
!
!              rr = 1.2e-9/1.0e6
!
!              Reaction = &
!                   rr * &
!                   IDensityS(iLon,iLat,iAlt,iHeP_,iBlock) * &
!                   NDensityS(iLon,iLat,iAlt,iN2_,iBlock)
!
!              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
!              NeutralSources(iHe_)   = NeutralSources(iHe_)   + Reaction
!              IonSources(iNP_)       = IonSources(iNP_)       + Reaction
!              NeutralLosses(iN2_)    = NeutralLosses(iN2_)    + Reaction
!              IonLosses(iHeP_)       = IonLosses(iHeP_)       + Reaction
!
!              ChemicalHeatingSub = &
!                   ChemicalHeatingSub + &
!                   Reaction * 0.28

              ! -----------
              ! O+(2P) + N -> N+ + O + 2.7 eV
              ! -----------

              rr = 1.0e-16

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iO_2PP_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock)

              NeutralSources(iO_)   = NeutralSources(iO_)   + Reaction
              IonSources(iNP_)      = IonSources(iNP_)      + Reaction
              NeutralLosses(iN_4S_) = NeutralLosses(iN_4S_) + Reaction
              IonLosses(iO_2PP_)    = IonLosses(iO_2PP_)    + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 2.7

              ! -----------
              ! O+(2D) + N -> N+ + O + 1.0 eV
              ! -----------

              rr = 7.5e-17

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iO_2DP_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock)

              NeutralSources(iO_)   = NeutralSources(iO_)   + Reaction
              IonSources(iNP_)      = IonSources(iNP_)      + Reaction
              NeutralLosses(iN_4S_) = NeutralLosses(iN_4S_) + Reaction
              IonLosses(iO_2DP_)    = IonLosses(iO_2DP_)    + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 1.0

              ! -----------
              ! N+ + O2 -> NO+ + O + 6.67 eV
              ! -----------

              rr = 2.6e-16

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iNP_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)

              NeutralSources(iO_) = NeutralSources(iO_) + Reaction
              IonSources(iNOP_)   = IonSources(iNOP_)   + Reaction
              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction
              IonLosses(iNP_)     = IonLosses(iNP_)     + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 6.67

              ! -----------
              ! N+ + O -> O+ + N + 0.93 eV
              ! -----------

              rr = 5.0e-19

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iNP_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iO_,iBlock)

              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              IonSources(iO_4SP_)    = IonSources(iO_4SP_)    + Reaction
              NeutralLosses(iO_)     = NeutralLosses(iO_)     + Reaction
              IonLosses(iNP_)        = IonLosses(iNP_)        + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 0.93

!              ! -----------
!              ! N+ + H -> H+ + N + 0.90 eV
!              ! -----------
!
!              rr = 3.6e-12/1.0e6
!
!              Reaction = &
!                   rr * &
!                   IDensityS(iLon,iLat,iAlt,iNP_,iBlock) * &
!                   NDensityS(iLon,iLat,iAlt,iH_,iBlock)
!
!              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
!              IonSources(iHP_)       = IonSources(iHP_)       + Reaction
!              NeutralLosses(iH_)     = NeutralLosses(iH_)     + Reaction
!              IonLosses(iNP_)        = IonLosses(iNP_)        + Reaction
!
!              ChemicalHeatingSub = &
!                   ChemicalHeatingSub + &
!                   Reaction * 0.9

              ! ----------------------------------------------------------
              ! NO+
              ! ----------------------------------------------------------

              ! -----------
              ! NO+ + e -> O + N(4S) + 0.38 eV  (0.78)
              ! NO+ + e -> O + N(2D) + 2.75 eV  (0.22)
              ! -----------

              rr = 4.2e-13 * te3m085

              Reaction = &
                   rr * &
                   IDensityS(iLon,iLat,iAlt,iNOP_,iBlock) * &
                   IDensityS(iLon,iLat,iAlt,ie_,iBlock)

              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + 0.78*Reaction
              NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + 0.22*Reaction
              NeutralSources(iO_)    = NeutralSources(iO_)    +      Reaction
              IonLosses(iNOP_)       = IonLosses(iNOP_)       +      Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * (0.38 * 0.78 + 2.75 * 0.22)

              ! ----------------------------------------------------------
              ! N(4S)
              ! ----------------------------------------------------------

              ! -----------
              ! N(2D) + e -> N(4S) + e + 2.38 eV
              ! -----------

              rr = 5.5e-16 * te3 ** (0.5)

              Reaction = &
                   rr * &
                   NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock) * &
                   IDensityS(iLon,iLat,iAlt,ie_,iBlock)

              ! We create and loose the same amount of e
              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              NeutralLosses(iN_2D_)  = NeutralLosses(iN_2D_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 2.38

              ! -----------
              ! N(2D) + O -> N(4S) + O + 2.38 eV
              ! -----------

              rr = 2.0e-18

              Reaction = &
                   rr * &
                   NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iO_,iBlock)

              ! We create and loose the same amount of e
              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              NeutralLosses(iN_2D_)  = NeutralLosses(iN_2D_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 2.38

              ! -----------
              ! N(2D) -> N(4S) + 5200A
              ! -----------

              rr = 1.06e-5

              Reaction = &
                   rr * &
                   NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock)

              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              NeutralLosses(iN_2D_)  = NeutralLosses(iN_2D_)  + Reaction

              Emission(iE5200_) = Emission(iE5200_) + Reaction

              ! -----------
              ! NO -> N(4S) + O
              ! -----------

              rr = 8.3e-6

              Reaction = &
                   rr * &
                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock)

              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              NeutralSources(iO_)    = NeutralSources(iO_)    + Reaction
              NeutralLosses(iNO_)    = NeutralLosses(iNO_)    + Reaction

              ! -----------
              ! N(4S) + O2 -> NO + O + 1.385 eV
              ! -----------

              rr = 4.4e-18 * exp(-3220/tn)

              Reaction = &
                   rr * &
                   NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)

              NeutralSources(iNO_)  = NeutralSources(iNO_)  + Reaction
              NeutralSources(iO_)   = NeutralSources(iO_)   + Reaction
              NeutralLosses(iN_4S_) = NeutralLosses(iN_4S_) + Reaction
              NeutralLosses(iO2_)   = NeutralLosses(iO2_)   + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 1.385

              ! -----------
              ! N(4S) + NO -> N2 + O + 3.25 eV
              ! -----------

              rr = 1.5e-18 * sqrt(tn)

              Reaction = &
                   rr * &
                   NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock)

              NeutralSources(iN2_)  = NeutralSources(iN2_)  + Reaction
              NeutralSources(iO_)   = NeutralSources(iO_)   + Reaction
              NeutralLosses(iN_4S_) = NeutralLosses(iN_4S_) + Reaction
              NeutralLosses(iNO_)   = NeutralLosses(iNO_)   + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 3.25

              ! -----------
              ! N(2P) -> N(2D) + 10400A
              ! -----------

              rr = 7.9e-2

              Reaction = &
                   rr * &
                   NDensityS(iLon,iLat,iAlt,iN_2P_,iBlock)

              NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + Reaction
              NeutralLosses(iN_2P_)  = NeutralLosses(iN_2P_)  + Reaction

              Emission(iE10400_) = Emission(iE10400_) + Reaction

              ! ----------------------------------------------------------
              ! N(2D)
              ! ----------------------------------------------------------

              ! -----------
              ! N(2D) + O2 -> NO + O + 3.76 eV
              ! -----------

              rr = 5.3e-18

              Reaction = &
                   rr * &
                   NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)

              NeutralSources(iNO_)  = NeutralSources(iNO_)  + Reaction
              NeutralSources(iO_)   = NeutralSources(iO_)   + Reaction
              NeutralLosses(iN_2D_) = NeutralLosses(iN_2D_) + Reaction
              NeutralLosses(iO2_)   = NeutralLosses(iO2_)   + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 3.76

              ! -----------
              ! N(2D) + NO -> N2 + O + 5.63 eV
              ! -----------

              rr = 7.0e-17

              Reaction = &
                   rr * &
                   NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock) * &
                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock)

              NeutralSources(iN2_)  = NeutralSources(iN2_)  + Reaction
              NeutralSources(iO_)   = NeutralSources(iO_)   + Reaction
              NeutralLosses(iN_2D_) = NeutralLosses(iN_2D_) + Reaction
              NeutralLosses(iNO_)   = NeutralLosses(iNO_)   + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 5.63

              ! ----------------------------------------------------------
              ! NO
              ! ----------------------------------------------------------
              ! -----------
              ! NO -> NO+ + e
              ! -----------

              rr = 6.0e-7

              Reaction = &
                   rr * &
                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock)

              IonSources(iNOP_)   = IonSources(iNOP_)   + Reaction
              NeutralLosses(iNO_) = NeutralLosses(iNO_) + Reaction

              do iIon = 1, nIons-1

                 if (UseIonChemistry .and. &
                      UseIonConstituent(iIon)) then

                    if (IonLosses(iIon) > IonSources(iIon) .and. &
                         IonSources(iIon)>0.0) then
                       if (IDensityS(iLon,iLat,iAlt,iIon,iBlock) > 0.01) then
                          dtOld = DtSub
                          dtSubtmp = &
                               min(0.25 * &
                               (IonSources(iIon) + &
                               IDensityS(iLon,iLat,iAlt,iIon,iBlock))/ &
                               (abs(IonLosses(iIon))+0.1), DtOld)

!                          if (DtSubtmp < 1.0e-2) then
!!                             write(*,*) "Yikes : Ion Chemistry is fast!"
!!                             write(*,*) "Loss changed from  : ",iIon, iLon, &
!!                                  Latitude(iLat,iBlock), iAlt, &
!!                                  IonLosses(iIon), IonSources(iIon), &
!!                                  IDensityS(iLon,iLat,iAlt,iIon,iBlock)
!!                             losstmp = &
!!                                  IonLosses(iIon) / &
!!                                  IDensityS(iLon,iLat,iAlt,iIon,iBlock)
!!                             dentmp = &
!!                                  (IonSources(iIon) + &
!!                                  IDensityS(iLon,iLat,iAlt,iIon,iBlock)) / &
!!                                  (1.0+losstmp)
!!                             IonSources(iIon) = 0.0
!!                             IonLosses(iIon) = &
!!                                  IDensityS(iLon,iLat,iAlt,iIon,iBlock) - &
!!                                  dentmp
!                             IonLosses(iIon) = min(IonSources(iIon)*5, &
!                                  IonLosses(iIon))
!!                             write(*,*) "to  : ",IonLosses(iIon)
!                             dtSub = &
!                                  min(0.25 * &
!                                  (IonSources(iIon) + &
!                                  IDensityS(iLon,iLat,iAlt,iIon,iBlock))/ &
!                                  (abs(IonLosses(iIon))+0.1), DtOld)
!                          else
                             DtSub = DtSubtmp
!                          endif
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
                       if (nDensityS(iLon,iLat,iAlt,iNeutral,iBlock)>0.01) then
                          dtOld = DtSub
                          dtSub = &
                               min(0.25 * &
                               nDensityS(iLon,iLat,iAlt,iNeutral,iBlock)/ &
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



              IDensityS(iLon,iLat,iAlt,nIons,iBlock) = 0.0
              do iIon = 1, nIons-1
                 IDensityS(iLon,iLat,iAlt,iIon,iBlock) = &
                      IDensityS(iLon,iLat,iAlt,iIon,iBlock) + &
                      (IonSources(iIon) - IonLosses(iIon)) * DtSub
                 IDensityS(iLon,iLat,iAlt,iIon,iBlock) = &
                      max(0.01,IDensityS(iLon,iLat,iAlt,iIon,iBlock))
                 
                 ! sum for e-
                 IDensityS(iLon,iLat,iAlt,nIons,iBlock) = &
                      IDensityS(iLon,iLat,iAlt,nIons,iBlock) + &
                      IDensityS(iLon,iLat,iAlt,iIon,iBlock)

                 if (IDensityS(iLon,iLat,iAlt,iIon,iBlock) < 0.0) then
                    write(*,*) "Negative Ion Density : ", &
                         iIon, iLon, iLat, iAlt, &
                         IDensityS(iLon,iLat,iAlt,iIon,iBlock), &
                         IonSources(iIon), IonLosses(iIon)
                 endif
              enddo

              do iNeutral = 1, nSpeciesTotal
                 NDensityS(iLon,iLat,iAlt,iNeutral,iBlock) = &
                      NDensityS(iLon,iLat,iAlt,iNeutral,iBlock) + &
                      (NeutralSources(iNeutral) - NeutralLosses(iNeutral)) * &
                      DtSub
              enddo

              ChemicalHeatingRate(iLon,iLat,iAlt) = &
                   ChemicalHeatingRate(iLon,iLat,iAlt) + &
                   ChemicalHeatingSub * DtSub

              Emissions(iLon, iLat, iAlt, :, iBlock) =  &
                   Emissions(iLon, iLat, iAlt, :, iBlock) + &
                   Emission(:)*DtSub

              DtTotal = DtTotal + DtSub

              if (DtSub < DtMin) DtMin = DtSub

              if (DtSub < 1.0e-9 .and. abs(DtTotal-Dt) > DtSub) then
                 write(*,*) "Chemistry is too fast!!", DtSub
                 if (iDtReducer < nIons) then
                    write(*,*) "Ion Constituent : ", &
                         iDtReducer, iLon, iLat, iAlt,&
                         IDensityS(iLon,iLat,iAlt,iDtReducer,iBlock), &
                         IonSources(iDtReducer), IonLosses(iDtReducer)
                 else
                    iDtReducer = iDtReducer - nIons
                    write(*,*) "Neutral Constituent : ", &
                         iDtReducer, iLon, iLat, iAlt,&
                         IDensityS(iLon,iLat,iAlt,iDtReducer,iBlock), &
                         IonSources(iDtReducer), IonLosses(iDtReducer)
                 endif

                 call stop_gitm("Ion Chemistry is too fast!!")
              endif

              nIters = nIters + 1

           enddo

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
