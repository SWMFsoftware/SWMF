      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     ALEX(10/11/04): I THINK THIS IS A COLD START ROUTINE FOR DEALING WITH
C     A LACK OF RESTART FILE. BASICALLY THE VELOCITY IS SET TO ZERO AND 
C     THE OTHER GRID VALUES ARE SET TO VALUES THAT WON'T CRASH.
C     (11/9/04) I also see that source terms and collision terms are set
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      

      SUBROUTINE STRT1
      use ModCommonVariables
      REAL jp1,jp2,jp3,jp4,kc1,kc2,kc3,kc6,kc7,kc8,kr1,kr2
      real DensityHp,DensityH3p
CALEX define the reaction rates, label by reaction number
CALEX j is for photochemistry, k is for regular chemistry
      jp1=9.5E-11
      jp2=5.4E-10
      jp3=7.3E-10
      jp4=1.3E-10
      kc1=2.E-9
      kc2=3.2E-29
      kc3=4.15E-9
      kc6=2.4E-9
      kc7=5.3E-9
      kc8=8.2E-9
      kr1=2.E-12
      kr2=4.6E-6

C     
C     
C     
C     
C     PHIOX=7.00E-7
C     
C     DEFINE THE HE PHOTOIONIZATION RATE
C     
      PHIHE=1.30E-7
C     
C     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     C
C     DEFINE THE GAS PARAMETERS AT THE LOWER BOUNDARY                  C
C     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     C
      USURFO=0.
      USURFH=0.
      USURHE=0.
      USURFE=0.
      CALL MODATM(ALTMIN,XNH2,XNH,XNH2O,XNCH4,TEMP)

CALEX I pretend that for plasma parameters, O is H3 and HE is
CALEX chemical equilibrium value for H2+ this allow me to just
CALEX change the chemistry but leave the rest of the code the same  
      TSURFO=TEMP
      TSURFH=TEMP
      TSURHE=TEMP
      TSURFE=TEMP

      call calc_chemical_equilibrium(DensityHp,DensityH3p)
      DSURFO=XMSO*DensityH3p
      DSURFH=XMSH*DensityHp
      write(*,*) 'equilibrium',DensityHp,DensityH3p


C I have used numerically calculated chemical equilibrium
C solution for T=800k.       
!      DSURFO=XMSO*6489.69
!c      DSURHE=XMSHE*jp2/kc1
!      DSURFHE=0.
!      DSURFH=XMSH*1343.64
C I have used numerically calculated chemical equilibrium
C solution for T=1000k.       
c      DSURFO=XMSO*4725.0
c      DSURFHE=0.
c      DSURFH=XMSH*368.0

C I have used numerically calculated chemical equilibrium
C solution for T=1500k.       
c      DSURFO=XMSO*560.0
c      DSURFHE=0.
c      DSURFH=XMSH*30.0


C I have used numerically calculated chemical equilibrium
C solution for T=1500k. with enhanced water and decreased CH4      
c      DSURFO=XMSO*11435.41
c      DSURFHE=0.
c      DSURFH=XMSH*1463.48

C I have used numerically calculated chemical equilibrium
C solution for T=100k. with reduced CH4 enhanced h2o      
C      DSURFO=XMSO*5509.0
C      DSURFHE=0.
C      DSURFH=XMSH*1124.0


      DSURFE=RTHDEL*DSURFH+RTOXEL*DSURFO
      PSURFO=RGASO*TSURFO*DSURFO
      PSURFH=RGASH*TSURFH*DSURFH
      PSURHE=RGASHE*TSURHE*DSURHE
      PSURFE=RGASE*TSURFE*DSURFE
      WSURFO=SQRT(GAMMA*RGASO*TSURFO)
      WSURFH=SQRT(GAMMA*RGASH*TSURFH)
      WSURHE=SQRT(GAMMA*RGASHE*TSURHE)
      WSURFE=SQRT(GAMMA*RGASE*TSURFE)
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     

      

      DO 20 I=1,NDIM
CALEX SOURCE COEF?       
         FFHpp1(I)=jp1*XH2(I)
         FFHpp3(I)=jp3*XH(I)
         FFHpp4(I)=jp4*XH2O(I)
         FFHpc2(I)=-kc2*XH2(I)*XH2(I)
         FFHpc3(I)=-kc3*XCH4(I)
         FFHpc8(I)=-kc8*XH2O(I)
         FFHpr1(I)=-kr1
CALEX write out source coeff
CALEX         write(26,*) FFHpp1(I),FFHpp3(I),FFHpp4(I),FFHpc2(I),FFHpc3(I),FFHpc8(I),FFHpr1(I)
        

         FFH3pc1(I)=kc1*(jp2/kc1)*XH2(I)
         FFH3pc2(I)=kc2*XH2(I)*XH2(I)
         FFH3pc6(I)=-kc6*XCH4(I)
         FFH3pc7(I)=-kc7*XH2O(I)
         FFH3pr2(I)=-kr2

CALEX write out source coeff
CALEX         write(27,*) FFH3pc1(I),FFH3pc2(I),FFH3pc6(I),FFH3pc7(I),FFH3pr2(I)

CALEX CL=COLLISION COEF, CF=collision freq ?         

CALEX the coulomb collisions
CAlex H+ and H3+         
         CLHpH3p(I)=1.905*4.**1.5/XMSO
CALEX electron H+ and electron H3+
         CLELHp(I)=54.5/XMSH
         CLELH3p(I)=54.5/XMSO
       
CALEX  ion neutrals
         CLHpH(I)=2.65E-10*XH(I)
         CFHpH2(I)=2.6E-9*XH2(I)*(.82/.667)**.5
         CFH3pH(I)=2.6E-9*XH(I)*(.667/.75)**.5
         CFH3pH2(I)=2.6E-9*XH2(I)*(.82/1.2)**.5
CALEX electron H, e H2 done in collis
         CLELH(I)=4.5E-9*XH(I)

         GRAVTY(I)=-3.79E22/RAD(I)**2
         Centrifugal(I)=RAD(I)*((sin((90.-GLAT)*3.14159/180.))**2)*Omega**2
         
 20   CONTINUE
c      GRAVTY(nDim)=GRAVTY(nDim)*0.0
cAlex
      NEXP=3
      WRITE(*,*) CURR(1),RAD(1),NEXP
cendAlex      
      CRRX=CURR(1)*RAD(1)**NEXP
      DO 25 I=1,NDIM
         CURR(I)=CRRX/RAD(I)**NEXP
 25   CONTINUE
C     SGN1=1.
C     IF (CURR(1).LT.0.) SGN1=-1.
C     IF (ABS(CURR(1)).LT.1.E-4) SGN1=0.
C     DO 30 I=NDIM-250,NDIM
C     CURR(I)=SGN1*0.2998*(RAD(1)/RAD(I))**NEXP
C     30    CONTINUE
      
CALEX CT & CM = energy collision term coef?
CALEX Based on Nagy p. 83, I believe that CT is the coeff
CALEX of the term that is due to temperature difference or
CALEX heat flow between species, and CM is the term due to 
CALEX frictional heating between species moving through each other

CALEX CTOXN2 = 3*R_o*M_o/(M_o+M_{N2}) see nagy p.83
      CTHpH=3.*RGASH*XMSH/(XMSH+XMSH)
      CTHpH2=3.*RGASH*XMSH/(XMSH+2.*XMSH)
      CTHpH3p=3.*RGASH*XMSH/(XMSH+XMSO)
      CTHpEL=3.*RGASH*XMSH/(XMSH+XMSE)

      CTH3pH=3.*RGASO*XMSO/(XMSO+XMSH)
      CTH3pH2=3.*RGASO*XMSO/(XMSO+2.*XMSH)
      CTH3pHp=3.*RGASO*XMSO/(XMSO+XMSH)
      CTH3pEL=3.*RGASO*XMSO/(XMSO+XMSE)

      CTELH=3.*RGASE*XMSE/(XMSE+XMSH)
      CTELH2=3.*RGASE*XMSE/(XMSE+2.*XMSH)
      CTELHp=3.*RGASE*XMSE/(XMSE+XMSH)
      CTELH3p=3.*RGASE*XMSE/(XMSE+XMSO)
      
CALEX CMOXN2 = M_{N2}/(M_o+M_{N2}) see nagy p.83
      CMHpH=XMSH/(XMSH+XMSH)
      CMHpH2=2.*XAMU/(XMSH+2.*XAMU)
      CMHpH3p=XMSO/(XMSH+XMSO)
      CMHpEL=XMSE/(XMSH+XMSE)
       
      CMH3pH=XMSH/(XMSO+XMSH)
      CMH3pH2=2.*XAMU/(XMSO+2.*XAMU)
      CMH3pHp=XMSH/(XMSO+XMSH)
      CMH3pEL=XMSE/(XMSO+XMSE)
       
      CMELH=XMSH/(XMSE+XMSH)
      CMELH2=2.*XAMU/(XMSE+2.*XAMU)
      CMELHp=XMSH/(XMSE+XMSH)
      CMELH3p=XMSO/(XMSE+XMSO)
     
C     ALEX(10/11/04): 
C     TRY SETTING THE PLASMA PARAMETERS HERE TO THE SURFACE VALUES

      if(IsRestart) RETURN

      do K=1,NDIM
         DHYD(K)=DSURFH*exp(-(ALTD(k)-1400.E5)/5000.E5)
         UOXYG(K)=0
         POXYG(K)=PSURFO*exp(-(ALTD(k)-1400.E5)/5000.E5)
         DOXYG(K)=DSURFO*exp(-(ALTD(k)-1400.E5)/5000.E5)
         TOXYG(K)=TSURFO
         UHYD(K)=0
         PHYD(K)=PSURFH*exp(-(ALTD(k)-1400.E5)/5000.E5)
         THYD(K)=TSURFH
         UHEL(K)=0
         PHEL(K)=PSURHE
         DHEL(K)=DSURHE
         THEL(K)=TSURHE
         DELECT(K)=DSURFE*exp(-(ALTD(k)-1400.E5)/5000.E5)
         UELECT(K)=0
         PELECT(K)=PSURFE*exp(-(ALTD(k)-1400.E5)/5000.E5)
         TELECT(K)=TSURFE
         
         
      enddo
      
      
      


      RETURN
      END
