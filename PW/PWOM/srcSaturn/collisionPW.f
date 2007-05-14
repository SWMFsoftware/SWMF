
CALEX This subroutine calculates the collision frequencies and
CALEX then calculates the momentum and energy collision terms
      SUBROUTINE COLLIS(N)
      use ModCommonVariables
C     
C
C
C
      real :: dT_II(nIon,nSpecies),dU2_II(nIon,nSpecies)
      
      if (State_GV(1,RhoO_) < 0.0) write(*,*) 'grendel f95 is bad'

      do I=1,N
         
C**********************************************************************
C Determine the mass sources from chemistry
C**********************************************************************

      Source_CV(I,RhoH_)=(FFHpp1(I)+FFHpp3(I)+FFHpp4(I)+FFHpc2(I)*State_GV(I,RhoH_)/Mass_I(Ion2_)
     ;+FFHpc3(I)*State_GV(I,RhoH_)/Mass_I(Ion2_)+FFHpc8(I)*State_GV(I,RhoH_)/Mass_I(Ion2_)+FFHpc9(I)*State_GV(I,RhoH_)/Mass_I(Ion2_)
     ;+FFHpr1(I)*State_GV(I,RhoH_)*State_GV(I,RhoE_)*(State_GV(I,Te_)**(-0.7))/Mass_I(Ion2_)/Mass_I(nIon))*Mass_I(Ion2_)

      Source_CV(I,RhoO_)=(FFH3pc1(I)+FFH3pc2(I)*State_GV(I,RhoH_)/Mass_I(Ion2_)+FFH3pc6(I)*State_GV(I,RhoO_)/Mass_I(Ion1_)
     ;+FFH3pc7(I)*State_GV(I,RhoO_)/Mass_I(Ion1_)
     ;+FFH3pr2(I)*(State_GV(I,Te_)**(-0.5))*State_GV(I,RhoO_)*State_GV(I,RhoE_)/Mass_I(Ion1_)/Mass_I(nIon))*Mass_I(Ion1_)
      
      Source_CV(I,RhoE_)=MassElecIon_I(Ion2_)*Source_CV(I,RhoH_)+MassElecIon_I(Ion1_)*Source_CV(I,RhoO_)




C**********************************************************************
C Calculate collision frequencies. 
C**********************************************************************

      TROX=0.5*(XTN(I)+State_GV(I,To_))
      TRHYD=0.5*(XTN(I)+State_GV(I,Th_))
      T1OX=SQRT(XTN(I)+State_GV(I,To_)/16.)

C These are (reduced temperatures) * (m1+m2) raised to the 1.5
C as shown on page 86 Nagy. This is for use in collision freqs
C of coulomb collisions below.       
      T1HpH3p=(State_GV(I,To_)+3.*State_GV(I,Th_))**1.5

      TE32=State_GV(I,Te_)**1.5
      DTE32=State_GV(I,RhoE_)/TE32


C H+ and H3+         
         CollisionFreq_IIC(Ion2_,Ion1_,I)=CLHpH3p(I)*State_GV(I,RhoO_)/T1HpH3p
C electron H+ and electron H3+
         CollisionFreq_IIC(nIon,Ion2_,I) = CLELHp(I)*State_GV(I,RhoH_)/TE32
         CollisionFreq_IIC(nIon,Ion1_,I)= CLELH3p(I)*State_GV(I,RhoO_)/TE32
                
C  ion neutrals
         CollisionFreq_IIC(Ion2_,Neutral2_,I)=
     &        CLHpH(I)*SQRT(TRHYD)*(1.-.083*ALOG10(TRHYD))**2.

C electron H, e H2 done in collis
         CollisionFreq_IIC(nIon,Neutral2_,I)=
     &        CLELH(I)*(1.-1.35E-4*State_GV(I,Te_))*SQRT(State_GV(I,Te_))
         CollisionFreq_IIC(nIon,Neutral1_,I)=
     &        getcfeh2(State_GV(I,Te_),XH2(I),Mass_I(nIon),State_GV(I,uE_))

C Now get the inverse collision freq
         CollisionFreq_IIC(Ion1_,Ion2_,I)=
     &        State_GV(I,RhoH_)/State_GV(I,RhoO_)*CollisionFreq_IIC(Ion2_,Ion1_,I)
         CollisionFreq_IIC(Ion2_,nIon,I)=
     &        State_GV(I,RhoE_)/State_GV(I,RhoH_)*CollisionFreq_IIC(nIon,Ion2_,I)
         CollisionFreq_IIC(Ion1_,nIon,I)=
     &        State_GV(I,RhoE_)/State_GV(I,RhoO_)*CollisionFreq_IIC(nIon,Ion1_,I)

C**********************************************************************
C Determin the momentum source terms
C**********************************************************************

C Velocity difference needed for source terms
      UHDOX=State_GV(I,uH_)-State_GV(I,uO_)
      UHDEL=State_GV(I,uH_)-State_GV(I,uE_)
      UOXEL=State_GV(I,uO_)-State_GV(I,uE_)

C This calculates collision source terms: 
C fclsn1=n*((u2-u1)*cf12+(u3-u1)*cf13+...)
      Source_CV(I,uO_)=State_GV(I,RhoO_)*(UHDOX*CollisionFreq_IIC(Ion1_,Ion2_,I)-
     $UOXEL*CollisionFreq_IIC(Ion1_,nIon,I)
     &-State_GV(I,uO_)
     &*(CollisionFreq_IIC(Ion1_,Neutral2_,I)+CollisionFreq_IIC(Ion1_,Neutral1_,I)))
      
      Source_CV(I,uH_)=State_GV(I,RhoH_)*(-UHDOX*CollisionFreq_IIC(Ion2_,Ion1_,I)-
     $UHDEL*CollisionFreq_IIC(Ion2_,nIon,I)-State_GV(I,uH_)
     &*(CollisionFreq_IIC(Ion2_,Neutral2_,I)+CollisionFreq_IIC(Ion2_,Neutral1_,I)))

 
      Source_CV(I,uE_)=State_GV(I,RhoE_)*(UOXEL*CollisionFreq_IIC(nIon,Ion1_,I)+
     $UHDEL*CollisionFreq_IIC(nIon,Ion2_,I)-State_GV(I,uE_)
     &*(CollisionFreq_IIC(nIon,Neutral2_,I)+CollisionFreq_IIC(nIon,Neutral1_,I)))



C**********************************************************************
C Determine the energy source terms
C**********************************************************************
C

      dU2_II(Ion2_,Ion1_)=UHDOX*UHDOX
      dU2_II(Ion1_,Ion2_)=dU2_II(Ion2_,Ion1_)

      
      dU2_II(Ion2_,nIon)=UHDEL*UHDEL
      dU2_II(nIon,Ion2_)=dU2_II(Ion2_,nIon)
      
      dU2_II(Ion1_,nIon)=UOXEL*UOXEL
      dU2_II(nIon,Ion1_)=dU2_II(Ion1_,nIon)
      
      dU2_II(Ion1_,Neutral1_:Neutral4_)=State_GV(I,uO_)**2
      dU2_II(Ion2_,Neutral1_:Neutral4_)=State_GV(I,uH_)**2
      dU2_II(nIon,Neutral1_:Neutral4_)=State_GV(I,uE_)**2
      

CALEX these are temperature differences needed in order to calculate
CALEX the energy collision term 
      dT_II(Ion2_,Ion1_)=State_GV(I,Th_)-State_GV(I,To_)
      dT_II(Ion2_,nIon)=State_GV(I,Th_)-State_GV(I,Te_)
      dT_II(Ion2_,Neutral1_:Neutral4_)=State_GV(I,Th_)-XTN(I)
      
      dT_II(Ion1_,Ion2_)=State_GV(I,To_)-State_GV(I,Th_)
      dT_II(Ion1_,nIon)=State_GV(I,To_)-State_GV(I,Te_)
      dT_II(Ion1_,Neutral1_:Neutral4_)=State_GV(I,To_)-XTN(I)

      dT_II(nIon,Ion2_)=State_GV(I,Te_)-State_GV(I,Th_)
      dT_II(nIon,Ion1_)=State_GV(I,Te_)-State_GV(I,To_)
      dT_II(nIon,Neutral1_:Neutral4_)=State_GV(I,Te_)-XTN(I)

      Source_CV(I,pO_) = 0.0
      Source_CV(I,pH_) = 0.0

      Source_CV(I,pE_) = 0.0
      do jSpecies=1,nSpecies
         if(Ion1_ /= jSpecies) Source_CV(I,pO_) = Source_CV(I,pO_) - dT_II(Ion1_,jSpecies)
     &  *HeatFlowCoef_II(Ion1_,jSpecies)*CollisionFreq_IIC(Ion1_,jSpecies,I)
     & + dU2_II(Ion1_,jSpecies)
     &  *FricHeatCoef_II(Ion1_,jSpecies)*CollisionFreq_IIC(Ion1_,jSpecies,I)
      enddo
      Source_CV(I,pO_) =State_GV(I,RhoO_)*Source_CV(I,pO_)

      do jSpecies=1,nSpecies
         if(Ion2_ /= jSpecies) Source_CV(I,pH_) = Source_CV(I,pH_) - dT_II(Ion2_,jSpecies)
     &  *HeatFlowCoef_II(Ion2_,jSpecies)*CollisionFreq_IIC(Ion2_,jSpecies,I)
     & + dU2_II(Ion2_,jSpecies)
     &  *FricHeatCoef_II(Ion2_,jSpecies)*CollisionFreq_IIC(Ion2_,jSpecies,I)
      enddo
      Source_CV(I,pH_) =State_GV(I,RhoH_)*Source_CV(I,pH_)

      do jSpecies=1,nSpecies
         if(nIon /= jSpecies) Source_CV(I,pE_) = Source_CV(I,pE_) - dT_II(nIon,jSpecies)
     &  *HeatFlowCoef_II(nIon,jSpecies)*CollisionFreq_IIC(nIon,jSpecies,I)
     & + dU2_II(nIon,jSpecies)
     &  *FricHeatCoef_II(nIon,jSpecies)*CollisionFreq_IIC(nIon,jSpecies,I)
      enddo
      Source_CV(I,pE_) =State_GV(I,RhoE_)*Source_CV(I,pE_)


C**********************************************************************
C Calculate heat conductivities
C**********************************************************************

      HeatCon_GI(I,Ion1_)=HLPO*(State_GV(nDim,RhoO_)/State_GV(nDim,RhoE_))*State_GV(I,To_)**2.5
      HeatCon_GI(I,nIon)=HLPE*State_GV(I,Te_)**2.5
      HeatCon_GI(I,Ion2_)=HLPH*(State_GV(nDim,RhoH_)/State_GV(nDim,RhoE_))*State_GV(I,Th_)**2.5
      enddo

      RETURN
      END
