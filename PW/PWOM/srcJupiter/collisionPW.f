
CALEX This subroutine calculates the collision frequencies and
CALEX then calculates the momentum and energy collision terms
      SUBROUTINE COLLIS(N,StateIn_GV)
      use ModCommonVariables
      use ModCommonPlanet,ONLY: HLPion1,HLPion2,HLPion3,HLPE
      use ModPhotoElectron
      integer, intent(in) :: N
      real,    intent(in) :: StateIn_GV(-1:N+2,nVar)
      
C     
C
C
C
      real :: dT_II(nIon,nSpecies),dU2_II(nIon,nSpecies)
      
      if (StateIn_GV(1,RhoH3_) < 0.0) write(*,*) 'grendel f95 is bad'

      do I=1,N
         
C**********************************************************************
C Determine the mass sources from chemistry
C**********************************************************************

C KGS Put in source
!      write(*,*) 'test5a'
      Source_CV(I,RhoH2_)=(TotalIonRate_IC(Ion1_,I) + 
     &        FFH2pc1(I)*StateIn_GV(I,RhoH2_)/Mass_I(Ion3_) +
     &        FFH2pc9(I)*StateIn_GV(I,RhoH_)/Mass_I(Ion2_))*Mass_I(Ion3_)

!      write(*,*) 'test5b'
      Source_CV(I,RhoH_)=( TotalIonRate_IC(Ion2_,I) +! s + 
     &     FFHpc2(I)*StateIn_GV(I,RhoH_)/Mass_I(Ion2_) + 
     &     FFHpc3(I)*StateIn_GV(I,RhoH_)/Mass_I(Ion2_) +
     &     FFHpc8(I)*StateIn_GV(I,RhoH_)/Mass_I(Ion2_) +
     &     FFHpc9(I)*StateIn_GV(I,RhoH_)/Mass_I(Ion2_) +
     &     FFHpr1(I)*StateIn_GV(I,RhoH_)/Mass_I(Ion2_)*
     &     StateIn_GV(I,RhoE_)/Mass_I(nIon)*
     &     (StateIn_GV(I,Te_)**(-0.7)) ) * Mass_I(Ion2_)

!      write(*,*) 'test5c'
C KGS H2+ originally folded in to H3pc1 & now put in explicitly
      Source_CV(I,RhoH3_)=(FFH3pc1(I)*StateIn_GV(I,RhoH2_)/Mass_I(Ion3_) +
     &     FFH3pc2(I)*StateIn_GV(I,RhoH_)/Mass_I(Ion2_) +
     &     FFH3pc6(I)*StateIn_GV(I,RhoH3_)/Mass_I(Ion1_) +
     &     FFH3pc7(I)*StateIn_GV(I,RhoH3_)/Mass_I(Ion1_) +
     &     FFH3pr2(I)*(StateIn_GV(I,Te_)**(-0.5))*StateIn_GV(I,RhoH3_)
     &     *StateIn_GV(I,RhoE_)/Mass_I(Ion1_)/Mass_I(nIon))*Mass_I(Ion1_)
      
!      write(*,*) 'test5d'
      Source_CV(I,RhoE_)=MassElecIon_I(Ion2_)*Source_CV(I,RhoH_)
     ;+MassElecIon_I(Ion1_)*Source_CV(I,RhoH3_)+MassElecIon_I(Ion3_)*Source_CV(I,RhoH2_)




C**********************************************************************
C Calculate collision frequencies. 
C**********************************************************************
!      write(*,*) 'test5e'
      TRHYD=0.5*(XTN(I)+StateIn_GV(I,Th_))
!      write(*,*) 'test5ea'
C These are (reduced temperatures) * (m1+m2) raised to the 1.5
C as shown on page 86 Nagy. This is for use in collision freqs
C of coulomb collisions below.       
      T1HpH3p=(StateIn_GV(I,Th3_)+3.*StateIn_GV(I,Th_))**1.5
!      write(*,*) 'test5eb',I,StateIn_GV(I,Th3_)
!      write(*,*) StateIn_GV(I,Th2_)
      T1H2pH3p=(2.*StateIn_GV(I,Th3_)+3.*StateIn_GV(I,Th2_))**1.5
!      write(*,*) 'test5ec'
      T1H2pHp=(2.*StateIn_GV(I,Th_)+StateIn_GV(I,Th2_))**1.5     
!      write(*,*) 'test5ed'

      TE32=StateIn_GV(I,Te_)**1.5
!      write(*,*) 'test5ee',I,TE32
      DTE32=StateIn_GV(I,RhoE_)/TE32
!      write(*,*) 'test5f'

C H2+, H+, and H3+
 ! KGS fix terms containing H2+ in startup
         CollisionFreq_IIC(Ion2_,Ion1_,I)=CLHpH3p(I)*StateIn_GV(I,RhoH3_)/T1HpH3p
         CollisionFreq_IIC(Ion3_,Ion1_,I)=CLH2pH3p(I)*StateIn_GV(I,RhoH3_)/T1H2pH3p
         CollisionFreq_IIC(Ion3_,Ion2_,I)=CLH2pHp(I)*StateIn_GV(I,RhoH_)/T1H2pHp
C electron H+ and electron H3+
         CollisionFreq_IIC(nIon,Ion3_,I) = CLELH2p(I)*StateIn_GV(I,RhoH2_)/TE32  ! KGS check this
         CollisionFreq_IIC(nIon,Ion2_,I) = CLELHp(I)*StateIn_GV(I,RhoH_)/TE32
         CollisionFreq_IIC(nIon,Ion1_,I)= CLELH3p(I)*StateIn_GV(I,RhoH3_)/TE32
!         write(*,*) 'test5g'
C  ion neutrals
         CollisionFreq_IIC(Ion2_,Neutral2_,I)=
     &        CLHpH(I)*SQRT(TRHYD)*(1.-.083*ALOG10(TRHYD))**2.

!         write(*,*) 'test5h'
C electron H, e H2 done in collis
         CollisionFreq_IIC(nIon,Neutral2_,I)=
     &        CLELH(I)*(1.-1.35E-4*StateIn_GV(I,Te_))*SQRT(StateIn_GV(I,Te_))
         CollisionFreq_IIC(nIon,Neutral1_,I)=
     &        getcfeh2(StateIn_GV(I,Te_),XH2(I),Mass_I(nIon),StateIn_GV(I,uE_))

!         write(*,*) 'test5i'
C Now get the inverse collision freq
         CollisionFreq_IIC(Ion2_,Ion3_,I)=
     &        StateIn_GV(I,RhoH2_)/StateIn_GV(I,RhoH_)*CollisionFreq_IIC(Ion3_,Ion2_,I)
!         write(*,*) 'test5j'
         CollisionFreq_IIC(Ion1_,Ion3_,I)=
     &        StateIn_GV(I,RhoH2_)/StateIn_GV(I,RhoH3_)*CollisionFreq_IIC(Ion3_,Ion1_,I)
!         write(*,*) 'test5k'
         CollisionFreq_IIC(Ion1_,Ion2_,I)=
     &        StateIn_GV(I,RhoH_)/StateIn_GV(I,RhoH3_)*CollisionFreq_IIC(Ion2_,Ion1_,I)
!         write(*,*) 'test5l'
         CollisionFreq_IIC(Ion3_,nIon,I)=
     &        StateIn_GV(I,RhoE_)/StateIn_GV(I,RhoH2_)*CollisionFreq_IIC(nIon,Ion3_,I)
         CollisionFreq_IIC(Ion2_,nIon,I)=
     &        StateIn_GV(I,RhoE_)/StateIn_GV(I,RhoH_)*CollisionFreq_IIC(nIon,Ion2_,I)
         CollisionFreq_IIC(Ion1_,nIon,I)=
     &        StateIn_GV(I,RhoE_)/StateIn_GV(I,RhoH3_)*CollisionFreq_IIC(nIon,Ion1_,I)

!         write(*,*) 'test5m'
C**********************************************************************
C Determin the momentum source terms
C**********************************************************************

C Velocity difference needed for source terms
      dU_32=StateIn_GV(I,uH2_)-StateIn_GV(I,uH_)
      dU_31=StateIn_GV(I,uH2_)-StateIn_GV(I,uH3_)
      dU_21=StateIn_GV(I,uH_)-StateIn_GV(I,uH3_)
      dU_3e=StateIn_GV(I,uH2_)-StateIn_GV(I,uE_)
      dU_2e=StateIn_GV(I,uH_)-StateIn_GV(I,uE_)
      dU_1e=StateIn_GV(I,uH3_)-StateIn_GV(I,uE_)

C This calculates collision source terms: 
C fclsn1=n*((u2-u1)*cf12+(u3-u1)*cf13+...)
C KGS Check H2+ here
      Source_CV(I,uH3_)=StateIn_GV(I,RhoH3_)*(dU_21*CollisionFreq_IIC(Ion1_,Ion2_,I)
     &     + dU_31*CollisionFreq_IIC(Ion1_,Ion3_,I)
     &     - dU_1e*CollisionFreq_IIC(Ion1_,nIon,I)
     &     - StateIn_GV(I,uH3_)
     &      * (CollisionFreq_IIC(Ion1_,Neutral2_,I)+CollisionFreq_IIC(Ion1_,Neutral1_,I)))
!      write(*,*) 'test5n'
      Source_CV(I,uH_)=StateIn_GV(I,RhoH_)*(-dU_21*CollisionFreq_IIC(Ion2_,Ion1_,I)
     &     + dU_32*CollisionFreq_IIC(Ion2_,Ion3_,I)
     &     - dU_2e*CollisionFreq_IIC(Ion2_,nIon,I)
     &     - StateIn_GV(I,uH_)
     &      * (CollisionFreq_IIC(Ion2_,Neutral2_,I)+CollisionFreq_IIC(Ion2_,Neutral1_,I)))
!      write(*,*) 'test5o'
      Source_CV(I,uH2_)=StateIn_GV(I,RhoH2_)*(-dU_31*CollisionFreq_IIC(Ion3_,Ion1_,I)
     &     - dU_32*CollisionFreq_IIC(Ion3_,Ion1_,I)
     &     - dU_3e*CollisionFreq_IIC(Ion3_,nIon,I)
     &     - StateIn_GV(I,uH2_)
     &      * (CollisionFreq_IIC(Ion3_,Neutral2_,I)+CollisionFreq_IIC(Ion3_,Neutral1_,I)))
!      write(*,*) 'test5p'
      Source_CV(I,uE_)=StateIn_GV(I,RhoE_)*(dU_1e*CollisionFreq_IIC(nIon,Ion1_,I)
     &     + dU_2e*CollisionFreq_IIC(nIon,Ion2_,I)
     &     + dU_3e*CollisionFreq_IIC(nIon,Ion3_,I)
     &     - StateIn_GV(I,uE_)
     &      * (CollisionFreq_IIC(nIon,Neutral2_,I)+CollisionFreq_IIC(nIon,Neutral1_,I)))
!      write(*,*) 'test5q'

C**********************************************************************
C Determine the energy source terms
C**********************************************************************
C

      dU2_II(Ion2_,Ion1_)=dU_21*dU_21
      dU2_II(Ion1_,Ion2_)=dU2_II(Ion2_,Ion1_)

!      write(*,*) 'test5r'
      dU2_II(Ion3_,Ion1_)=dU_31*dU_31
      dU2_II(Ion1_,Ion3_)=dU2_II(Ion3_,Ion1_)

!      write(*,*) 'test5s'
      dU2_II(Ion3_,Ion2_)=dU_32*dU_32
      dU2_II(Ion2_,Ion3_)=dU2_II(Ion3_,Ion2_)

!      write(*,*) 'test5t'
      dU2_II(Ion3_,nIon)=dU_3e*dU_3e
      dU2_II(nIon,Ion3_)=dU2_II(Ion3_,nIon)

!      write(*,*) 'test5u'
      dU2_II(Ion2_,nIon)=dU_2e*dU_2e
      dU2_II(nIon,Ion2_)=dU2_II(Ion2_,nIon)
      
!      write(*,*) 'test5v'
      dU2_II(Ion1_,nIon)=dU_1e*dU_1e
      dU2_II(nIon,Ion1_)=dU2_II(Ion1_,nIon)
      
!      write(*,*) 'test5w'
      dU2_II(Ion1_,Neutral1_:Neutral4_)=StateIn_GV(I,uH3_)**2
      dU2_II(Ion2_,Neutral1_:Neutral4_)=StateIn_GV(I,uH_)**2
      dU2_II(Ion3_,Neutral1_:Neutral4_)=StateIn_GV(I,uH2_)**2
      dU2_II(nIon,Neutral1_:Neutral4_)=StateIn_GV(I,uE_)**2
      
!      write(*,*) 'test5x'
CALEX these are temperature differences needed in order to calculate
CALEX the energy collision term 
      dT_II(Ion3_,Ion2_)=StateIn_GV(I,Th2_)-StateIn_GV(I,Th_)
      dT_II(Ion3_,Ion1_)=StateIn_GV(I,Th2_)-StateIn_GV(I,Th3_)
      dT_II(Ion3_,nIon)=StateIn_GV(I,Th2_)-StateIn_GV(I,Te_)
      dT_II(Ion3_,Neutral1_:Neutral4_)=StateIn_GV(I,Th2_)-XTN(I)

!      write(*,*) 'test5y'
      dT_II(Ion2_,Ion3_)=StateIn_GV(I,Th_)-StateIn_GV(I,Th2_)
      dT_II(Ion2_,Ion1_)=StateIn_GV(I,Th_)-StateIn_GV(I,Th3_)
      dT_II(Ion2_,nIon)=StateIn_GV(I,Th_)-StateIn_GV(I,Te_)
      dT_II(Ion2_,Neutral1_:Neutral4_)=StateIn_GV(I,Th_)-XTN(I)
      
!      write(*,*) 'test5z'
      dT_II(Ion1_,Ion3_)=StateIn_GV(I,Th3_)-StateIn_GV(I,Th2_)
      dT_II(Ion1_,Ion2_)=StateIn_GV(I,Th3_)-StateIn_GV(I,Th_)
      dT_II(Ion1_,nIon)=StateIn_GV(I,Th3_)-StateIn_GV(I,Te_)
      dT_II(Ion1_,Neutral1_:Neutral4_)=StateIn_GV(I,Th3_)-XTN(I)

!      write(*,*) 'test5za'
      dT_II(nIon,Ion3_)=StateIn_GV(I,Te_)-StateIn_GV(I,Th2_)
      dT_II(nIon,Ion2_)=StateIn_GV(I,Te_)-StateIn_GV(I,Th_)
      dT_II(nIon,Ion1_)=StateIn_GV(I,Te_)-StateIn_GV(I,Th3_)
      dT_II(nIon,Neutral1_:Neutral4_)=StateIn_GV(I,Te_)-XTN(I)

      Source_CV(I,pH3_) = 0.0
      Source_CV(I,pH_) = 0.0
      Source_CV(I,pH2_) = 0.0

      Source_CV(I,pE_) = 0.0

      if(UseFeedbackFromSE) then
         !add the energy deposition from SEs
         Source_CV(I,pE_) =Source_CV(I,pE_)+SeHeat_C(I)
      endif

      do jSpecies=1,nSpecies
         if(Ion1_ /= jSpecies) Source_CV(I,pH3_) = Source_CV(I,pH3_) - dT_II(Ion1_,jSpecies)
     &  *HeatFlowCoef_II(Ion1_,jSpecies)*CollisionFreq_IIC(Ion1_,jSpecies,I)
     & + dU2_II(Ion1_,jSpecies)
     &  *FricHeatCoef_II(Ion1_,jSpecies)*CollisionFreq_IIC(Ion1_,jSpecies,I)
      enddo
      Source_CV(I,pH3_) =StateIn_GV(I,RhoH3_)*Source_CV(I,pH3_)

      do jSpecies=1,nSpecies
         if(Ion2_ /= jSpecies) Source_CV(I,pH_) = Source_CV(I,pH_) - dT_II(Ion2_,jSpecies)
     &  *HeatFlowCoef_II(Ion2_,jSpecies)*CollisionFreq_IIC(Ion2_,jSpecies,I)
     & + dU2_II(Ion2_,jSpecies)
     &  *FricHeatCoef_II(Ion2_,jSpecies)*CollisionFreq_IIC(Ion2_,jSpecies,I)
      enddo
      Source_CV(I,pH_) =StateIn_GV(I,RhoH_)*Source_CV(I,pH_)

      do jSpecies=1,nSpecies
         if(Ion3_ /= jSpecies) Source_CV(I,pH2_) = Source_CV(I,pH2_) - dT_II(Ion3_,jSpecies)
     &  *HeatFlowCoef_II(Ion3_,jSpecies)*CollisionFreq_IIC(Ion3_,jSpecies,I)
     & + dU2_II(Ion3_,jSpecies)
     &  *FricHeatCoef_II(Ion3_,jSpecies)*CollisionFreq_IIC(Ion3_,jSpecies,I)
      enddo
      Source_CV(I,pH2_) =StateIn_GV(I,RhoH2_)*Source_CV(I,pH2_)

      do jSpecies=1,nSpecies
         if(nIon /= jSpecies) Source_CV(I,pE_) = Source_CV(I,pE_) - dT_II(nIon,jSpecies)
     &  *HeatFlowCoef_II(nIon,jSpecies)*CollisionFreq_IIC(nIon,jSpecies,I)
     & + dU2_II(nIon,jSpecies)
     &  *FricHeatCoef_II(nIon,jSpecies)*CollisionFreq_IIC(nIon,jSpecies,I)
      enddo
      Source_CV(I,pE_) =StateIn_GV(I,RhoE_)*Source_CV(I,pE_)


C**********************************************************************
C Calculate heat conductivities
C**********************************************************************

      HeatCon_GI(I,Ion1_)=HLPion1*(StateIn_GV(nDim,RhoH3_)/StateIn_GV(nDim,RhoE_))*StateIn_GV(I,Th3_)**2.5
      HeatCon_GI(I,nIon)=HLPE*StateIn_GV(I,Te_)**2.5
      HeatCon_GI(I,Ion2_)=HLPion2*(StateIn_GV(nDim,RhoH_)/StateIn_GV(nDim,RhoE_))*StateIn_GV(I,Th_)**2.5
      HeatCon_GI(I,Ion3_)=HLPion3*(StateIn_GV(nDim,RhoH2_)/StateIn_GV(nDim,RhoE_))*StateIn_GV(I,Th2_)**2.5
      enddo

      RETURN
      END
