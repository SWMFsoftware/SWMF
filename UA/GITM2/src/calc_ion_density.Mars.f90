
subroutine calc_ion_density(iBlock)

  use ModSize
  use ModGITM
  use ModPlanet
  use ModRates
  use ModEUV
  use ModSources
  use ModInputs, only: iDebugLevel, UseIonChemistry
  use ModConstants

  implicit none

  integer, intent(in) :: iBlock

  real :: Sources(nIons)
  real :: Losses(nIons)
  real :: DtSub, DtOld, DtTotal, DtMin, DtAve

  integer :: iLon, iLat, iAlt, iIon, nIters, iDtReducer

  real :: LatDebug = 1002.5*pi/180, LonDebug = 172.5*pi/180

  DtMin = Dt

  if (.not.UseIonChemistry) return

  call report("Ion Chemistry",2)
  call start_timing("calc_ion_density")

  DtAve = 0.0

  nIters=0

  write(*,*) "Ion Chemistry does not work!!!"
  return

  do iLon = 1, nLons
     do iLat = 1, nLats
        do iAlt = 1, nAlts

           DtTotal = 0.0

           do while (DtTotal < Dt)

              DtSub = Dt - DtTotal

              Sources = 0.0
              Losses  = 0.0

!              ! N+
!
!              Sources(iNP_) = &
!                   EuvIonRateS(iLon,iLat,iAlt,iNP_,iBlock)* &
!                   NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock)
!
!              Losses(iNP_) = &
!                   rRTempInd(iRrK6_)* NDensityS(iLon,iLat,iAlt,iO2_,iBlock)+ &
!                   rRTempInd(iRrK7_)* NDensityS(iLon,iLat,iAlt,iO2_,iBlock)+ &
!                   rRTempInd(iRrK8_)* NDensityS(iLon,iLat,iAlt,iO_, iBlock)
!              Losses(iNP_) = &
!                   Losses(iNP_)*IDensityS(iLon,iLat,iAlt,iNP_,iBlock)
!
!              ! O+
!
!              Sources(iO_4SP_) = &
!                   EuvIonRateS(iLon,iLat,iAlt,iO_4SP_,iBlock) * &
!                   NDensityS(iLon,iLat,iAlt,iO_,iBlock) + &
!                   AuroralIonRateS(iLon,iLat,iAlt,iO_,iBlock)
!
!              Sources(iO_4SP_) = Sources(iO_4SP_) + &
!                   rRTempInd(iRrK8_)*IDensityS(iLon,iLat,iAlt,iNP_,iBlock) * &
!                   NDensityS(iLon,iLat,iAlt,iO_,iBlock)
!
!              Losses(iO_4SP_) =  &
!                   RrTempDep(iLon,iLat,iAlt,iRrK1_) * &
!                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock) + &
!                   RrTempDep(iLon,iLat,iAlt,iRrK2_) * &
!                   NDensityS(iLon,iLat,iAlt,iN2_,iBlock)
!
!              Losses(iO_4SP_) = Losses(iO_4SP_) * &
!                   IDensityS(iLon,iLat,iAlt,iO_4SP_,iBlock)
!
!              ! O2+
!
!              Sources(iO2P_) = &
!                   EuvIonRateS(iLon,iLat,iAlt,iO2P_,iBlock) * &
!                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock) + &
!                   AuroralIonRateS(iLon,iLat,iAlt,iO2_, iBlock)
!
!              Sources(iO2P_) = Sources(iO2P_) + &
!                   RrTempDep(iLon,iLat,iAlt,iRrK1_) * &
!                   IDensityS(iLon,iLat,iAlt,iO_4SP_,iBlock) * &
!                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)
!
!              Sources(iO2P_) = Sources(iO2P_) + &
!                   RrTempInd(iRrK6_) * &
!                   IDensityS(iLon,iLat,iAlt,iNP_,iBlock) * &
!                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)
!
!              Losses(iO2P_) = Losses(iO2P_) + &
!                   RrTempDep(iLon,iLat,iAlt,iRrA2_) * &
!                   IDensityS(iLon,iLat,iAlt,ie_,iBlock) + &
!                   RrTempInd(iRrK32_) * &
!                   NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock) + &
!                   RrTempInd(iRrK5_) * &
!                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock)
!              Losses(iO2P_) = Losses(iO2P_) * &
!                   IDensityS(iLon,iLat,iAlt,iO2P_,iBlock)
!
!              ! N2+
!
!              Sources(iN2P_) = &
!                   EuvIonRateS(iLon,iLat,iAlt,iN2P_,iBlock) * &
!                   NDensityS(iLon,iLat,iAlt,iN2_,iBlock) +  &
!                   AuroralIonRateS(iLon,iLat,iAlt,iN2_, iBlock)
!
!              Losses(iN2P_) = Losses(iN2P_) + &
!                   RrTempDep(iLon,iLat,iAlt,iRrK3_) * &
!                   NDensityS(iLon,iLat,iAlt,iO_,iBlock) + &
!                   RrTempDep(iLon,iLat,iAlt,iRrA3_) * &
!                   IDensityS(iLon,iLat,iAlt,ie_,iBlock)
!              Losses(iN2P_) = Losses(iN2P_) * &
!                   IDensityS(iLon,iLat,iAlt,iN2P_,iBlock)
!
!              ! NO+
!
!              Sources(iNOP_) = Sources(iNOP_) + &
!                   RrTempDep(iLon,iLat,iAlt,iRrK3_) * &
!                   IDensityS(iLon,iLat,iAlt,iN2P_,iBlock) * &
!                   NDensityS(iLon,iLat,iAlt,iO_,iBlock)        + &
!                   RrTempDep(iLon,iLat,iAlt,iRrK2_) * &
!                   IDensityS(iLon,iLat,iAlt,iO_4SP_,iBlock) * &
!                   NDensityS(iLon,iLat,iAlt,iN2_,iBlock)        + &
!                   RrTempDep(iLon,iLat,iAlt,iRrBeta9N_) * &
!                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock)        + &
!                   RrTempInd(iRrK32_) * &
!                   IDensityS(iLon,iLat,iAlt,iO2P_,iBlock) * &
!                   NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock)        + &
!                   RrTempInd(iRrK5_) * &
!                   IDensityS(iLon,iLat,iAlt,iO2P_,iBlock) * &
!                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock)        + &
!                   RrTempInd(iRrK7_) * &
!                   IDensityS(iLon,iLat,iAlt,iNP_,iBlock) * &
!                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)
!
!              Losses(iNOP_) = Losses(iNOP_) + &
!                   RrTempDep(iLon,iLat,iAlt,iRrA1_) * &
!                   IDensityS(iLon,iLat,iAlt,ie_,iBlock) * &
!                   IDensityS(iLon,iLat,iAlt,iNOP_,iBlock)


              do iIon = 1, nIons-1
                 if (Losses(iIon) > Sources(iIon)) then
                    if (IDensityS(iLon,iLat,iAlt,iIon,iBlock) > 0.01) then
                       dtOld = DtSub
                       dtSub = &
                            min(0.25 * &
                            IDensityS(iLon,iLat,iAlt,iIon,iBlock)/ &
                            (abs(Sources(iIon) - Losses(iIon))+0.1), DtOld)

                       if (DtSub < DtOld) iDtReducer = iIon

                       if (DtSub < DtOld .and. iDebugLevel > 5) &
                            write(*,*) "====> Ion Chem, dt reduced : ", &
                            DtSub, DtOld,iAlt, iIon, &
                            Sources(iIon), Losses(iIon), &
                            IDensityS(iLon,iLat,iAlt,iIon,iBlock)
                    else
                       Sources(iIon) = 0.0
                       Losses(iIon) = 0.0
                    endif
                 endif
              enddo

              IDensityS(iLon,iLat,iAlt,nIons,iBlock) = 0.0
              do iIon = 1, nIons-1
                 IDensityS(iLon,iLat,iAlt,iIon,iBlock) = &
                      IDensityS(iLon,iLat,iAlt,iIon,iBlock) + &
                      (Sources(iIon) - Losses(iIon)) * DtSub
                 ! sum for e-
                 IDensityS(iLon,iLat,iAlt,nIons,iBlock) = &
                      IDensityS(iLon,iLat,iAlt,nIons,iBlock) + &
                      IDensityS(iLon,iLat,iAlt,iIon,iBlock)
              enddo

              DtTotal = DtTotal + DtSub

              if (DtSub < DtMin) DtMin = DtSub

              if (DtSub < 1.0e-5 .and. abs(DtTotal-Dt) > 0.01) then
                 write(*,*) "Ion Chemistry is too fast!!", DtSub
                 write(*,*) "Constituent : ", iDtReducer, &
                      IDensityS(iLon,iLat,iAlt,iDtReducer,iBlock), &
                      sources(iDtReducer), losses(iDtReducer)
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

  call end_timing("calc_ion_density")

end subroutine calc_ion_density
