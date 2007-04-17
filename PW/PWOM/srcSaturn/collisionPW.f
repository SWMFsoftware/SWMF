
CALEX This subroutine calculates the collision frequencies and
CALEX then calculates the momentum and energy collision terms
      SUBROUTINE COLLIS(N)
      use ModCommonVariables
C     
C
C
C
      real :: dT_II(nIon,nSpecies),dU2_II(nIon,nSpecies)
      
      do I=1,N
         
C**********************************************************************
C Determine the mass sources from chemistry
C**********************************************************************

      ADMSH(I)=(FFHpp1(I)+FFHpp3(I)+FFHpp4(I)+FFHpc2(I)*DHYD(I)/XMSH
     ;+FFHpc3(I)*DHYD(I)/XMSH+FFHpc8(I)*DHYD(I)/XMSH+FFHpc9(I)*DHYD(I)/XMSH
     ;+FFHpr1(I)*DHYD(I)*DELECT(I)*(Telect(I)**(-0.7))/XMSH/XMSE)*XMSH

      ADMSO(I)=(FFH3pc1(I)+FFH3pc2(I)*DHYD(I)/XMSH+FFH3pc6(I)*DOXYG(I)/XMSO
     ;+FFH3pc7(I)*DOXYG(I)/XMSO
     ;+FFH3pr2(I)*(TELECT(I)**(-0.5))*DOXYG(I)*DELECT(I)/XMSO/XMSE)*XMSO
      
      ADMSHE(I)=0.
      ADMSE(I)=RTHDEL*ADMSH(I)+RTOXEL*ADMSO(I)




C**********************************************************************
C Calculate collision frequencies. 
C**********************************************************************

      TROX=0.5*(XTN(I)+TOXYG(I))
      TRHYD=0.5*(XTN(I)+THYD(I))
      TRHEL=0.5*(XTN(I)+THEL(I))
      T1OX=SQRT(XTN(I)+TOXYG(I)/16.)

C These are (reduced temperatures) * (m1+m2) raised to the 1.5
C as shown on page 86 Nagy. This is for use in collision freqs
C of coulomb collisions below.       
      T1HpH3p=(TOXYG(I)+3.*THYD(I))**1.5
      T1OXHE=(TOXYG(I)+4.*THEL(I))**1.5
      T1HEH=(THEL(I)+4.*THYD(I))**1.5

      TE32=TELECT(I)**1.5
      DTE32=DELECT(I)/TE32


C H+ and H3+         
         CollisionFreq_IIC(Ion2_,Ion1_,I)=CLHpH3p(I)*DOXYG(I)/T1HpH3p
C electron H+ and electron H3+
         CollisionFreq_IIC(Ion4_,Ion2_,I) = CLELHp(I)*DHYD(I)/TE32
         CollisionFreq_IIC(Ion4_,Ion1_,I)= CLELH3p(I)*DOXYG(I)/TE32
                
C  ion neutrals
         CollisionFreq_IIC(Ion2_,Neutral2_,I)=
     &        CLHpH(I)*SQRT(TRHYD)*(1.-.083*ALOG10(TRHYD))**2.

C electron H, e H2 done in collis
         CollisionFreq_IIC(Ion4_,Neutral2_,I)=
     &        CLELH(I)*(1.-1.35E-4*TELECT(I))*SQRT(TELECT(I))
         CollisionFreq_IIC(Ion4_,Neutral1_,I)=
     &        getcfeh2(TELECT(I),XH2(I),XMSE,UELECT(I))

C Now get the inverse collision freq
         CollisionFreq_IIC(Ion1_,Ion2_,I)=
     &        DHYD(I)/DOXYG(I)*CollisionFreq_IIC(Ion2_,Ion1_,I)
         CollisionFreq_IIC(Ion2_,Ion4_,I)=
     &        DELECT(I)/DHYD(I)*CollisionFreq_IIC(Ion4_,Ion2_,I)
         CollisionFreq_IIC(Ion1_,Ion4_,I)=
     &        DELECT(I)/DOXYG(I)*CollisionFreq_IIC(Ion4_,Ion1_,I)

C**********************************************************************
C Determin the momentum source terms
C**********************************************************************

C Velocity difference needed for source terms
      UHDOX=UHYD(I)-UOXYG(I)
      UHDHL=UHYD(I)-UHEL(I)
      UHDEL=UHYD(I)-UELECT(I)
      UOXEL=UOXYG(I)-UELECT(I)

C This calculates collision source terms: 
C fclsn1=n*((u2-u1)*cf12+(u3-u1)*cf13+...)
      FCLSNO(I)=DOXYG(I)*(UHDOX*CollisionFreq_IIC(Ion1_,Ion2_,I)-
     $UOXEL*CollisionFreq_IIC(Ion1_,Ion4_,I)
     &-UOXYG(I)
     &*(CollisionFreq_IIC(Ion1_,Neutral2_,I)+CollisionFreq_IIC(Ion1_,Neutral1_,I)))
      
      FCLSNH(I)=DHYD(I)*(-UHDOX*CollisionFreq_IIC(Ion2_,Ion1_,I)-
     $UHDEL*CollisionFreq_IIC(Ion2_,Ion4_,I)-UHYD(I)
     &*(CollisionFreq_IIC(Ion2_,Neutral2_,I)+CollisionFreq_IIC(Ion2_,Neutral1_,I)))

      FCLSHE(I)=0.
 
      FCLSNE(I)=DELECT(I)*(UOXEL*CollisionFreq_IIC(Ion4_,Ion1_,I)+
     $UHDEL*CollisionFreq_IIC(Ion4_,Ion2_,I)-UELECT(I)
     &*(CollisionFreq_IIC(Ion4_,Neutral2_,I)+CollisionFreq_IIC(Ion4_,Neutral1_,I)))



C**********************************************************************
C Determine the energy source terms
C**********************************************************************
C

      dU2_II(Ion2_,Ion1_)=UHDOX*UHDOX
      dU2_II(Ion1_,Ion2_)=dU2_II(Ion2_,Ion1_)

      dU2_II(Ion2_,Ion3_)=UHDHL*UHDHL
      dU2_II(Ion3_,Ion2_)=dU2_II(Ion2_,Ion3_)
      
      dU2_II(Ion2_,Ion4_)=UHDEL*UHDEL
      dU2_II(Ion4_,Ion2_)=dU2_II(Ion2_,Ion4_)
      
      dU2_II(Ion3_,Ion1_)=UHEOX*UHEOX
      dU2_II(Ion1_,Ion3_)=dU2_II(Ion3_,Ion1_)

      dU2_II(Ion3_,Ion4_)=UHEEL*UHEEL
      dU2_II(Ion4_,Ion3_)=dU2_II(Ion3_,Ion4_)

      dU2_II(Ion1_,Ion4_)=UOXEL*UOXEL
      dU2_II(Ion4_,Ion1_)=dU2_II(Ion1_,Ion4_)
      
      dU2_II(Ion1_,Neutral1_:Neutral5_)=UOXYG(I)**2
      dU2_II(Ion2_,Neutral1_:Neutral5_)=UHyd(I)**2
      dU2_II(Ion3_,Neutral1_:Neutral5_)=UHel(I)**2
      dU2_II(Ion4_,Neutral1_:Neutral5_)=UElect(I)**2
      

CALEX these are temperature differences needed in order to calculate
CALEX the energy collision term 
      dT_II(Ion2_,Ion1_)=THYD(I)-TOXYG(I)
      dT_II(Ion2_,Ion3_)=THYD(I)-THEL(I)
      dT_II(Ion2_,Ion4_)=THYD(I)-TELECT(I)
      dT_II(Ion2_,Neutral1_:Neutral5_)=THYD(I)-XTN(I)
      
      dT_II(Ion1_,Ion2_)=TOxyg(I)-THyd(I)
      dT_II(Ion1_,Ion3_)=TOxyg(I)-THEL(I)
      dT_II(Ion1_,Ion4_)=TOxyg(I)-TELECT(I)
      dT_II(Ion1_,Neutral1_:Neutral5_)=TOxyg(I)-XTN(I)

      dT_II(Ion3_,Ion2_)=THel(I)-THyd(I)
      dT_II(Ion3_,Ion1_)=THel(I)-TOxyg(I)
      dT_II(Ion3_,Ion4_)=THel(I)-TELECT(I)
      dT_II(Ion3_,Neutral1_:Neutral5_)=THel(I)-XTN(I)

      dT_II(Ion4_,Ion2_)=TElect(I)-THyd(I)
      dT_II(Ion4_,Ion1_)=TElect(I)-TOxyg(I)
      dT_II(Ion4_,Ion3_)=TElect(I)-THel(I)
      dT_II(Ion4_,Neutral1_:Neutral5_)=TElect(I)-XTN(I)

      ECLSNO(I) = 0.0
      ECLSNH(I) = 0.0
      ECLSHE(I) = 0.0
      ECLSNE(I) = 0.0
      do jSpecies=1,nSpecies
         if(Ion1_ /= jSpecies) ECLSNO(I) = ECLSNO(I) - dT_II(Ion1_,jSpecies)
     &  *HeatFlowCoef_II(Ion1_,jSpecies)*CollisionFreq_IIC(Ion1_,jSpecies,I)
     & + dU2_II(Ion1_,jSpecies)
     &  *FricHeatCoef_II(Ion1_,jSpecies)*CollisionFreq_IIC(Ion1_,jSpecies,I)
      enddo
      ECLSNO(I) =dOxyg(I)*ECLSNO(I)

      do jSpecies=1,nSpecies
         if(Ion2_ /= jSpecies) ECLSNH(I) = ECLSNH(I) - dT_II(Ion2_,jSpecies)
     &  *HeatFlowCoef_II(Ion2_,jSpecies)*CollisionFreq_IIC(Ion2_,jSpecies,I)
     & + dU2_II(Ion2_,jSpecies)
     &  *FricHeatCoef_II(Ion2_,jSpecies)*CollisionFreq_IIC(Ion2_,jSpecies,I)
      enddo
      ECLSNH(I) =dHyd(I)*ECLSNH(I)

      do jSpecies=1,nSpecies
         if(Ion3_ /= jSpecies) ECLSHE(I) = ECLSHE(I) - dT_II(Ion3_,jSpecies)
     &  *HeatFlowCoef_II(Ion3_,jSpecies)*CollisionFreq_IIC(Ion3_,jSpecies,I)
     & + dU2_II(Ion3_,jSpecies)
     &  *FricHeatCoef_II(Ion3_,jSpecies)*CollisionFreq_IIC(Ion3_,jSpecies,I)
      enddo
      ECLSHE(I) =dHel(I)*ECLSHE(I)

      do jSpecies=1,nSpecies
         if(Ion4_ /= jSpecies) ECLSNE(I) = ECLSNE(I) - dT_II(Ion4_,jSpecies)
     &  *HeatFlowCoef_II(Ion4_,jSpecies)*CollisionFreq_IIC(Ion4_,jSpecies,I)
     & + dU2_II(Ion4_,jSpecies)
     &  *FricHeatCoef_II(Ion4_,jSpecies)*CollisionFreq_IIC(Ion4_,jSpecies,I)
      enddo
      ECLSNE(I) =dElect(I)*ECLSNE(I)


C**********************************************************************
C Calculate heat conductivities
C**********************************************************************

      TCONO(I)=HLPO*(DOXYG(NDIM)/DELECT(NDIM))*TOXYG(I)**2.5
      TCONE(I)=HLPE*TELECT(I)**2.5
      TCONH(I)=HLPH*(DHYD(NDIM)/DELECT(NDIM))*THYD(I)**2.5
      TCONHE(I)=0.
      enddo

      RETURN
      END
