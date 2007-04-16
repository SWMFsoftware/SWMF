
CALEX This subroutine calculates the collision frequencies and
CALEX then calculates the momentum and energy collision terms
      SUBROUTINE COLLIS(N)
      use ModCommonVariables
      
      real :: dT_II(nIon,nSpecies),dU2_II(nIon,nSpecies)
C     
C
C
C
      DO 250 I=1,N
      TOX2=TOXYG(I)*TOXYG(I)
      TOX3=TOX2*TOXYG(I)
      TOX4=TOX2*TOX2
      AAH=SQRT(THYD(I)+XTN(I)/16.+1.2E-8*UHYD(I)*UHYD(I))
      BBH=SQRT(XTN(I)+TOXYG(I)/16.+1.2E-8*UOXYG(I)*UOXYG(I))
      CCH=1.-9.149E-4*TOXYG(I)+4.228E-7*TOX2-6.790E-11*TOX3+
     ;4.225E-15*TOX4
      IF (TOXYG(I).GT.1700.) GO TO 201
      FFOXYG=FFOX4(I)*(1.-1.290E-3*TOXYG(I)+6.233E-7*TOX2)
      GO TO 205
201   FFOXYG=FFOX5(I)*(1.-1.410E-3*TOXYG(I)+6.036E-7*TOX2)
205   ADMSO(I)=FFOX1(I)+FFOX2(I)*AAH*DHYD(I)+(FFOX3(I)*BBH+
     ;FFOXYG+FFOX6(I)*CCH)*DOXYG(I)


      ADMSH(I)=FFHYD1(I)*BBH*DOXYG(I)+FFHYD2(I)*AAH*DHYD(I)
      ADMSHE(I)=FFHE1(I)+FFHE2(I)*DHEL(I)
      ADMSE(I)=RTHDEL*ADMSH(I)+RTOXEL*ADMSO(I)+RTHEEL*ADMSHE(I)
C
      TROX=0.5*(XTN(I)+TOXYG(I))
      TRHYD=0.5*(XTN(I)+THYD(I))
      TRHEL=0.5*(XTN(I)+THEL(I))
      T1OX=SQRT(XTN(I)+TOXYG(I)/16.)

CALEX These are (reduced temperatures) * (m1+m2) raised to the 1.5
CALEX as shown on page 86 Nagy. This is for use in collision freqs
CALEX of coulomb collisions below.       
      T1OXH=(TOXYG(I)+16.*THYD(I))**1.5
      T1OXHE=(TOXYG(I)+4.*THEL(I))**1.5
      T1HEH=(THEL(I)+4.*THYD(I))**1.5
      TE32=TELECT(I)**1.5
      DTE32=DELECT(I)/TE32
CALEX calculate collision frequencies. 

CALEX cf of O+ and O
!      CFOXO(I)=CLOXO(I)*SQRT(TROX)*(1.-0.064*ALOG10(TROX))**2
      CollisionFreq_IIC(Ion1_,Neutral1_,I)=
     &     CLOXO(I)*SQRT(TROX)*(1.-0.064*ALOG10(TROX))**2
CALEX this looks like cf of O+ and H but see comments for cloxh
!      CFOXH(I)=CLOXH(I)*T1OX
      CollisionFreq_IIC(Ion1_,Neutral2_,I)=CLOXH(I)*T1OX
CALEX cf of O+ and H+. coulomb collisions as eq 4.142 p.95 in Nagy
!      CFOXHD(I)=CLOXHD(I)*DHYD(I)/T1OXH
      CollisionFreq_IIC(Ion1_,Ion2_,I)=CLOXHD(I)*DHYD(I)/T1OXH
CALEX cf of O+ and He+. coulomb collisions as eq 4.142 p.95 in Nagy      
!      CFOXHL(I)=CLOXHL(I)*DHEL(I)/T1OXHE
      CollisionFreq_IIC(Ion1_,Ion3_,I)=CLOXHL(I)*DHEL(I)/T1OXHE
!      CFOXEL(I)=CLOXEL(I)*DTE32
      CollisionFreq_IIC(Ion1_,Ion4_,I)=CLOXEL(I)*DTE32
CALEX  cf of H+ and O      
!      CFHO(I)=CLHO(I)*SQRT(THYD(I))*(1.-0.047*ALOG10(THYD(I)))**2
      CollisionFreq_IIC(Ion2_,Neutral1_,I)=
     &     CLHO(I)*SQRT(THYD(I))*(1.-0.047*ALOG10(THYD(I)))**2

CALEX cf of H+ and H      
!      CFHH(I)=CLHH(I)*SQRT(TRHYD)*(1.-0.083*ALOG10(TRHYD))**2
      CollisionFreq_IIC(Ion2_,Neutral2_,I)=
     &     CLHH(I)*SQRT(TRHYD)*(1.-0.083*ALOG10(TRHYD))**2
CALEX cf of H+ and O+. coulomb collisions as eq 4.142 p.95 in Nagy      
!      CFHOX(I)=CLHOX(I)*DOXYG(I)/T1OXH
      CollisionFreq_IIC(Ion2_,Ion1_,I)=CLHOX(I)*DOXYG(I)/T1OXH
CALEX cf of H+ and He+. coulomb collisions as eq 4.142 p.95 in Nagy      
!      CFHHL(I)=CLHHL(I)*DHEL(I)/T1HEH
      CollisionFreq_IIC(Ion2_,Ion3_,I)=CLHHL(I)*DHEL(I)/T1HEH
!      CFHEL(I)=CLHEL(I)*DTE32
      CollisionFreq_IIC(Ion2_,Ion4_,I)=CLHEL(I)*DTE32
CALEX cf of He+ and HE      
!      CFHEHE(I)=CLHEHE(I)*SQRT(TRHEL)*(1.-0.093*ALOG10(TRHEL))**2
      CollisionFreq_IIC(Ion3_,Neutral5_,I)=
     &     CLHEHE(I)*SQRT(TRHEL)*(1.-0.093*ALOG10(TRHEL))**2
CALEX cf of He+ and O+      
!      CFHEOX(I)=CLHEOX(I)*DOXYG(I)/T1OXHE
      CollisionFreq_IIC(Ion3_,Ion1_,I)=CLHEOX(I)*DOXYG(I)/T1OXHE
CALEX cf of He+ and H+      
!      CFHEHD(I)=CLHEHD(I)*DHYD(I)/T1HEH
      CollisionFreq_IIC(Ion3_,Ion2_,I)=CLHEHD(I)*DHYD(I)/T1HEH
CALEX cfheel=7.43E-3/M_e * n_e / T_e^1.5, this seems like it should be
CALEX the cf between He and electrons but the formula does not match 4.144
!      CFHEEL(I)=CLHEEL(I)*DTE32
      CollisionFreq_IIC(Ion3_,Ion4_,I)=CLHEEL(I)*DTE32
CALEX cf for electrons and O, same as 4.144      
!      CFELOX(I)=CLELOX(I)*DOXYG(I)/TE32
      CollisionFreq_IIC(Ion4_,Ion1_,I)=CLELOX(I)*DOXYG(I)/TE32
CALEX cf for electrons and He      
!      CFELHL(I)=CLELHL(I)*DHEL(I)/TE32
      CollisionFreq_IIC(Ion4_,Ion3_,I)=CLELHL(I)*DHEL(I)/TE32
CALEX cf for el and H      
!      CFELHD(I)=CLELHD(I)*DHYD(I)/TE32
      CollisionFreq_IIC(Ion4_,Ion2_,I)=CLELHD(I)*DHYD(I)/TE32
C
CALEX velocity difference needed for source terms
      UHDOX=UHYD(I)-UOXYG(I)
      UHDHL=UHYD(I)-UHEL(I)
      UHDEL=UHYD(I)-UELECT(I)
      UHEOX=UHEL(I)-UOXYG(I)
      UHEEL=UHEL(I)-UELECT(I)
      UOXEL=UOXYG(I)-UELECT(I)

CALEX This calculates collision source terms: 
CALEX fclsn1=n*((u2-u1)*cf12+(u3-u1)*cf13+...)
      FCLSNO(I) = DOXYG(I)*(
     &       UHDOX*CollisionFreq_IIC(Ion1_,Ion2_,I)
     &     + UHEOX*CollisionFreq_IIC(Ion1_,Ion3_,I)
     &     - UOXEL*CollisionFreq_IIC(Ion1_,Ion4_,I)
     &     -UOXYG(I) * (
     &       CollisionFreq_IIC(Ion1_,Neutral1_,I)
     &     + CollisionFreq_IIC(Ion1_,Neutral2_,I)
     &     + CollisionFreq_IIC(Ion1_,Neutral3_,I)
     &     + CollisionFreq_IIC(Ion1_,Neutral4_,I)
     &     + CollisionFreq_IIC(Ion1_,Neutral5_,I)))
      
      FCLSNH(I) = DHYD(I)*(
     &     - UHDOX*CollisionFreq_IIC(Ion2_,Ion1_,I)
     &     - UHDHL*CollisionFreq_IIC(Ion2_,Ion3_,I)
     &     - UHDEL*CollisionFreq_IIC(Ion2_,Ion4_,I)
     &     -UHYD(I) * (
     &       CollisionFreq_IIC(Ion2_,Neutral1_,I)
     &     + CollisionFreq_IIC(Ion2_,Neutral2_,I)
     &     + CollisionFreq_IIC(Ion2_,Neutral3_,I)
     &     + CollisionFreq_IIC(Ion2_,Neutral4_,I)
     &     + CollisionFreq_IIC(Ion2_,Neutral5_,I)))

      FCLSHE(I) = DHEL(I)*(
     &       UHDHL*CollisionFreq_IIC(Ion3_,Ion2_,I)
     &     - UHEOX*CollisionFreq_IIC(Ion3_,Ion1_,I)
     &     - UHEEL*CollisionFreq_IIC(Ion3_,Ion4_,I)
     &     -UHEL(I)*(
     &       CollisionFreq_IIC(Ion3_,Neutral1_,I)
     &     + CollisionFreq_IIC(Ion3_,Neutral2_,I)
     &     + CollisionFreq_IIC(Ion3_,Neutral3_,I)
     &     + CollisionFreq_IIC(Ion3_,Neutral4_,I)
     &     + CollisionFreq_IIC(Ion3_,Neutral5_,I)))

CALEX UHLEL is not defined anywhere in this program!!! based on
CALEX UOXEL I assume it is the relative velocity between helium
CALEX and electrionsand so I define it as such below
      UHLEL=UHEEL
      

CALEX     $UHDEL,CFELHD(I),UELECT(I),CFELN2(I),CFELO2(I),CFELO(I),
CALEX     $CFELHE(I),CFELH(I)',
CALEX     $ DELECT(I),UOXEL,CFELOX(I)
     

      FCLSNE(I)=DELECT(I)*(
     &       UOXEL*CollisionFreq_IIC(Ion4_,Ion1_,I)
     &     + UHLEL*CollisionFreq_IIC(Ion4_,Ion3_,I)
     &     + UHDEL*CollisionFreq_IIC(Ion4_,Ion2_,I)
     &     -UELECT(I)*(
     &       CollisionFreq_IIC(Ion4_,Neutral1_,I)
     &     + CollisionFreq_IIC(Ion4_,Neutral2_,I)
     &     + CollisionFreq_IIC(Ion4_,Neutral3_,I)
     &     + CollisionFreq_IIC(Ion4_,Neutral4_,I)
     &     + CollisionFreq_IIC(Ion4_,Neutral5_,I)))

C
CALEX calculate square of velocity differences for use in the
CALEX energy collision term
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


CALEX These are the energy collision terms as seen in eq 4.86 in Nagy
!      ECLSNO(I)=DOXYG(I)*(-TOXN*(
!     &  HeatFlowCoef_II(Ion1_,Neutral1_)*CollisionFreq_IIC(Ion1_,Neutral1_,I)
!     & +HeatFlowCoef_II(Ion1_,Neutral2_)*CollisionFreq_IIC(Ion1_,Neutral2_,I)
!     & +HeatFlowCoef_II(Ion1_,Neutral3_)*CollisionFreq_IIC(Ion1_,Neutral3_,I)
!     & +HeatFlowCoef_II(Ion1_,Neutral4_)*CollisionFreq_IIC(Ion1_,Neutral4_,I)
!     & +HeatFlowCoef_II(Ion1_,Neutral5_)*CollisionFreq_IIC(Ion1_,Neutral5_,I))
!     & + THDOX*HeatFlowCoef_II(Ion1_,Ion2_)*CollisionFreq_IIC(Ion1_,Ion2_,I)
!     & + THEOX*HeatFlowCoef_II(Ion1_,Ion3_)*CollisionFreq_IIC(Ion1_,Ion3_,I)
!     & - TOXEL*HeatFlowCoef_II(Ion1_,Ion4_)*CollisionFreq_IIC(Ion1_,Ion4_,I)
!     & +UOXN*(
!     &  FricHeatCoef_II(Ion1_,Neutral1_)*CollisionFreq_IIC(Ion1_,Neutral1_,I)
!     & +FricHeatCoef_II(Ion1_,Neutral2_)*CollisionFreq_IIC(Ion1_,Neutral2_,I)
!     & +FricHeatCoef_II(Ion1_,Neutral3_)*CollisionFreq_IIC(Ion1_,Neutral3_,I)
!     & +FricHeatCoef_II(Ion1_,Neutral4_)*CollisionFreq_IIC(Ion1_,Neutral4_,I)
!     & +FricHeatCoef_II(Ion1_,Neutral5_)*CollisionFreq_IIC(Ion1_,Neutral5_,I))
!     $+UHDOX*FricHeatCoef_II(Ion1_,Ion2_)*CollisionFreq_IIC(Ion1_,Ion2_,I)
!     & +UHEOX*FricHeatCoef_II(Ion1_,Ion3_)*CollisionFreq_IIC(Ion1_,Ion3_,I)
!     $ +UOXEL*FricHeatCoef_II(Ion1_,Ion4_)*CollisionFreq_IIC(Ion1_,Ion4_,I))

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


C
CALEX calculate heat conductivities
      TCONO(I)=HLPO*(DOXYG(I)/DELECT(I))*TOXYG(I)**2.5
      TCONE(I)=HLPE*TELECT(I)**2.5
      TCONH(I)=HLPH*(DHYD(I)/DELECT(I))*THYD(I)**2.5
      TCONH(I)=TCONH(I)/(1.+(0.7692*
     &(CollisionFreq_IIC(Ion2_,Neutral4_,I)
     & + CollisionFreq_IIC(Ion2_,Neutral3_,I) )
     &+ 1.0962*CollisionFreq_IIC(Ion2_,Neutral1_,I) )
     & /CollisionFreq_IIC(Ion2_,Ion1_,I))

      TCONHE(I)=HLPHE*(DHEL(I)/XMSHE)*THEL(I)

      TCONHE(I)=TCONHE(I)/(0.99*CollisionFreq_IIC(Ion3_,Neutral4_,I)
     & +0.99*CollisionFreq_IIC(Ion3_,Neutral3_,I)
     & +1.02*CollisionFreq_IIC(Ion3_,Neutral1_,I)
     & +1.48*CollisionFreq_IIC(Ion3_,Neutral5_,I)
     & +2.22*CollisionFreq_IIC(Ion3_,Neutral2_,I)
     & +1.21*CollisionFreq_IIC(Ion3_,Ion1_,I)
     & +2.23*CollisionFreq_IIC(Ion3_,Ion2_,I))      


250   CONTINUE
      RETURN
      END
