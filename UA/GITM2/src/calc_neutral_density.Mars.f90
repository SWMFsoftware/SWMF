! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

subroutine calc_neutral_density(iBlock)

  use ModSizeGitm
  use ModGITM
  use ModPlanet
  use ModRates
  use ModEUV
  use ModInputs, only: iDebugLevel

  implicit none

  integer, intent(in) :: iBlock

  real :: Sources(nSpeciesTotal)
  real :: Losses(nSpeciesTotal)
  real :: DtSub, DtOld, DtTotal, DtMin

  integer :: iLon, iLat, iAlt, iNeutral, iDtReducer, nIters

  DtMin = Dt

  call report("Neutral Chemistry",2)
  call start_timing("Neutral Chemistry")

  write(*,*) "Neutral Chemistry does not work!!!"
  return

  nIters = 0

  do iLon = 1, nLons
     do iLat = 1, nLats
        do iAlt = 1, nAlts

           DtTotal = 0.0

           do while (DtTotal < Dt)

              DtSub = Dt - DtTotal

              Sources = 0.0
              Losses  = 0.0

!              ! [O]!!!!!!!!!!
!
!              ! Sources:
!              Sources(iO_)=  &
!                   RrTempDep(iLon,iLat,iAlt,iRrK1_) *   &
!                   IDensityS( iLon,iLat,iAlt,iO_4SP_,iBlock)   *  &
!                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)     +    &
!                   RrTempInd(iRrK32_) * &
!                   IDensityS(iLon,iLat,iAlt,iO2P_,iBlock) *  &
!                   NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock)         +    &
!                   RrTempInd(iRrK7_) * &
!                   IDensityS(iLon,iLat,iAlt,iNP_,iBlock)  *  &
!                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)           +    &
!                   RrTempDep(iLon,iLat,iAlt,iRrA1_) * &
!                   IDensityS(iLon,iLat,iAlt,iNOP_,iBlock)*  &
!                   IDensityS(iLon,iLat,iAlt,ie_,iBlock)+ &
!                   RrTempDep(iLon,iLat,iAlt,iRrBeta1_) * &
!                   NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock) * &
!                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)  +    &
!                   RrTempDep(iLon,iLat,iAlt,iRrBeta3_) * &
!                   NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock)*&
!                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock)  +    &
!                   RrTempInd(iRrBeta6_) * &
!                   NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock) *     &
!                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock)    +    &
!                   RrTempDep(iLon,iLat,iAlt,iRrBeta8_) * &
!                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock)  +& 
!                   RrTempInd(iRrBeta2_) * &
!                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)       *     &
!                   NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock)       +    &
!                   2*RrTempDep(iLon,iLat,iAlt,iRrA2_) * &
!                   IDensityS(iLon,iLat,iAlt,ie_,iBlock) * &
!                   IDensityS(iLon,iLat,iAlt,iO2P_,iBlock)       + &
!                   2*RrTempInd(iRrJ1_) * &
!                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)
!
!    ! Loss:
!              Losses(iO_)=  &
!                   2*RrTempDep(iLon,iLat,iAlt,iRrK12_) * &
!                   NDensityS(iLon,iLat,iAlt,iO_,iBlock) * &
!                   (NDensityS(iLon,iLat,iAlt,iO_,iBlock) +  &
!                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)+  &
!                   NDensityS(iLon,iLat,iAlt,iN2_,iBlock))+ &
!                   RrTempDep(iLon,iLat,iAlt,iRrK13_) * &
!                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock) * &
!                   (NDensityS(iLon,iLat,iAlt,iO_,iBlock) +  &
!                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock))+ &
!                   RrTempDep(iLon,iLat,iAlt,iRrK3_) * &
!                   IDensityS(iLon,iLat,iAlt,iN2P_,iBlock) + &
!                   RrTempInd(iRrK8_) * IDensityS(iLon,iLat,iAlt,iNP_, iBlock)
!
!              Losses(iO_)= Losses(iO_) * NDensityS(iLon,iLat,iAlt,iO_,iBlock)
!
!              ! [O2]!!!!!!!!!!
!
!              ! Sources:
!
!              Sources(iO2_) = &
!                   RrTempInd(iRrK5_) * &
!                   IDensityS(iLon,iLat,iAlt,iO2P_, iBlock) *  &
!                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock)+    &
!                   RrTempDep(iLon,iLat,iAlt,iRrK12_) * &
!                   NDensityS(iLon,iLat,iAlt,iO_,iBlock)**2 * &
!                   (NDensityS(iLon,iLat,iAlt,iO_,iBlock)  +   &
!                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock) +   &
!                   NDensityS(iLon,iLat,iAlt,iN2_,iBlock))
!              
!              ! Loss:
!              Losses(iO2_) = &
!                   RrTempDep(iLon,iLat,iAlt,iRrK1_) * &
!                   IDensityS(iLon,iLat,iAlt,iO_4SP_,iBlock) + &
!                   RrTempInd(iRrK6_) * &
!                   IDensityS(iLon,iLat,iAlt,iNP_, iBlock) +           &
!                   RrTempInd(iRrK7_) * &
!                   IDensityS(iLon,iLat,iAlt,iNP_, iBlock) +           &
!                   RrTempDep(iLon,iLat,iAlt,iRrBeta1_) * &
!                   NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock) +&
!                   RrTempInd(iRrBeta2_) * &
!                   NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock) +      &
!                   RrTempDep(iLon,iLat,iAlt,iRrK13_) * &
!                   NDensityS(iLon,iLat,iAlt,iO_,iBlock) *    &
!                   (NDensityS(iLon,iLat,iAlt,iO_,iBlock) +    &
!                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock))+   &
!                   RrTempInd(iRrJ1_) 
!
!              Losses(iO2_) = Losses(iO2_) * &
!                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)
!
!              ! [N2]!!!!!!!!!!
!              
!              ! Sources:
!              Sources(iN2_) = &
!                   RrTempDep(iLon,iLat,iAlt,iRrBeta3_) * &
!                   NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock) * &
!                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock) + &
!                   RrTempInd(iRrBeta6_) *  &
!                   NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock) * &
!                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock)
!
!              ! Loss:
!              Losses(iN2_) = &
!                   RrTempDep(iLon,iLat,iAlt,iRrK2_) * &
!                   IDensityS(iLon,iLat,iAlt,iO_4SP_,iBlock)
!
!              Losses(iN2_) = Losses(iN2_) * &
!                   NDensityS(iLon,iLat,iAlt,iN2_,iBlock)
!
!              ! [N(4s)]!!!!!!!!
!
!              ! Sources:
!              Sources(iN_4S_) = &
!                   RrTempDep(iLon,iLat,iAlt,iRrK2_) * &
!                   IDensityS(iLon,iLat,iAlt,iO_4SP_,iBlock) *    &
!                   NDensityS(iLon,iLat,iAlt,iN2_,iBlock) +       &
!                   RrTempInd(iRrK6_) * &
!                   IDensityS(iLon,iLat,iAlt,iNP_, iBlock) *           &
!                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock) +            &
!                   RrTempInd(iRrK8_) * &
!                   IDensityS(iLon,iLat,iAlt,iNP_, iBlock)*            &
!                   NDensityS(iLon,iLat,iAlt,iO_,iBlock) +             &
!                   0.2 * RrTempDep(iLon,iLat,iAlt,iRrA1_) * &
!                   IDensityS(iLon,iLat,iAlt,iNOP_,iBlock) *   &
!                   IDensityS(iLon,iLat,iAlt,ie_,iBlock) +                     &
!                   (2.*0.1+0.9)* RrTempDep(iLon,iLat,iAlt,iRrA3_) * &
!                   IDensityS(iLon,iLat,iAlt,iN2P_,iBlock) *   &
!                   IDensityS(iLon,iLat,iAlt,ie_,iBlock) +               &
!                   RrTempInd(iRrBeta4_) * &
!                   NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock) *        &
!                   NDensityS(iLon,iLat,iAlt,iO_,iBlock) +           &
!                   RrTempDep(iLon,iLat,iAlt,iRrBeta5_) * &
!                   NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock)* &
!                   IDensityS(iLon,iLat,iAlt,ie_,iBlock) +     &
!                   RrTempInd(iRrBeta7_) * &
!                   NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock) +        &
!                   RrTempDep(iLon,iLat,iAlt,iRrBeta8_) * &
!                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock)
!
!              ! Loss:
!              Losses(iN_4S_) = &
!                   RrTempInd(iRrK32_) * &
!                   IDensityS(iLon,iLat,iAlt,iO2P_,iBlock) +                  &
!                   RrTempDep(iLon,iLat,iAlt,iRrBeta1_) * &
!                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock) +  &
!                   RrTempDep(iLon,iLat,iAlt,iRrBeta3_) * &
!                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock)
!
!              Losses(iN_4S_) = Losses(iN_4S_) * &
!                   NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock)
!
!              ! [NO]!!!!!!!!!!!
!     
!              ! Sources:
!              Sources(iNO_)= &
!                   RrTempDep(iLon,iLat,iAlt,iRrBeta1_)* &
!                   NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock)*&
!                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)+ &
!                   RrTempInd(iRrBeta2_)* &
!                   NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock)*&
!                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock)
!
!              ! Loss:
!              Losses(iNO_)= &
!                   RrTempInd(iRrK5_)* &
!                   IDensityS(iLon,iLat,iAlt,iO2P_, iBlock)+          &
!                   RrTempDep(iLon,iLat,iAlt,iRrBeta3_)*&
!                   NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock)+&
!                   RrTempInd(iRrBeta6_)*&
!                   NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock)+      &
!                   RrTempDep(iLon,iLat,iAlt,iRrBeta8_)+          &
!                   RrTempDep(iLon,iLat,iAlt,iRrBeta9N_)
!
!              Losses(iNO_) = Losses(iNO_) * &
!                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock)
!
!              ! [N(2D)]!!!!!!!!!
!      
!              ! Sources:
!              Sources(iN_2D_)= &
!                   RrTempDep(iLon,iLat,iAlt,iRrK3_) * &
!                   IDensityS(iLon,iLat,iAlt,iN2P_,iBlock) *    &
!                   NDensityS(iLon,iLat,iAlt,iO_,iBlock) +             &
!                   0.8*RrTempDep(iLon,iLat,iAlt,iRrA1_) * &
!                   IDensityS(iLon,iLat,iAlt,iNOP_,iBlock) *      &
!                   IDensityS(iLon,iLat,iAlt,ie_,iBlock) +&
!                   0.9*RrTempDep(iLon,iLat,iAlt,iRrA3_) * &
!                   IDensityS(iLon,iLat,iAlt,iN2P_,iBlock) *      &
!                   IDensityS(iLon,iLat,iAlt,ie_,iBlock)
!
!              ! Loss:
!              Losses(iN_2D_)= &
!                   RrTempInd(iRrBeta2_) * &
!                   NDensityS(iLon,iLat,iAlt,iO2_,iBlock) + &
!                   RrTempInd(iRrBeta4_) * &
!                   NDensityS(iLon,iLat,iAlt,iO_,iBlock) +  &
!                   RrTempDep(iLon,iLat,iAlt,iRrBeta5_) * &
!                   IDensityS(iLon,iLat,iAlt,ie_,iBlock) + &
!                   RrTempInd(iRrBeta6_) * &
!                   NDensityS(iLon,iLat,iAlt,iNO_,iBlock) + &
!                   RrTempInd(iRrBeta7_)
!
!              Losses(iN_2D_) = Losses(iN_2D_) * &
!                   NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock)

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

              DtTotal = DtTotal + DtSub

              if (DtSub < DtMin) DtMin = DtSub

              if (DtSub < 1.0e-4 .and. abs(DtTotal-Dt) > 0.01) then
                 write(*,*) "Neutral Chemistry is too fast!!", DtSub
                 write(*,*) "Constituent : ", iDtReducer, &
                      nDensityS(iLon,iLat,iAlt,iDtReducer,iBlock), &
                      sources(iDtReducer), losses(iDtReducer)
                 call stop_gitm("Neutral Chemistry is too fast!!")
              endif

              nIters = nIters + 1

           enddo

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

end subroutine calc_neutral_density



