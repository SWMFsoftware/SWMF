      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     ALEX(10/11/04): I THINK THIS IS A COLD START ROUTINE FOR DEALING WITH
C     A LACK OF RESTART FILE. BASICALLY THE VELOCITY IS SET TO ZERO AND 
C     THE OTHER GRID VALUES ARE SET TO VALUES THAT WON'T CRASH.
C     (11/9/04) I also see that source terms and collision terms are set
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      

      SUBROUTINE STRT1
      use ModCommonVariables
C     
C     
C     
C     
C     PHIOX=7.00E-7
C     
C     DEFINE THE HE PHOTOIONIZATION RATE
C     
C      PHIHE=1.30E-7
      PHIHE=3.87E-8
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
      CALL MODATM(ALTMIN,XNO2,XNN2,XNO,XNH,XNHE,TEMP)
      TSURFO=TEMP
      TSURFH=TEMP
      TSURHE=TEMP
      TSURFE=TEMP
      TMP1=1.-1.290E-3*TEMP+6.233E-7*TEMP**2
      TMP2=1.-9.149E-4*TEMP+4.228E-7*TEMP**2-6.790E-11*TEMP**3+
     $     4.225E-15*TEMP**4
C     DSURFO=PHIOX*XMSO*XNO/(1.53E-12*XNN2*TMP1+2.82E-11*XNO2*TMP2)
      DSURFO=PHOTOTF(1)*XMSO*XNO/(1.53E-12*XNN2*TMP1+2.82E-11*XNO2*TMP2)
      DSURHE=PHIHE*XMSHE*XNHE/(1.10E-9*XNO2+1.60E-9*XNN2)
      DSURFH=1.136*(XNH/XNO)*(XMSH/XMSO)*DSURFO
      DSURFE=RTHDEL*DSURFH+RTOXEL*DSURFO+RTHEEL*DSURHE
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
C     FFOX1(1) BASED ON PHOTOIONIZATION FREQ FROM SUB GLOWEX;MASS OXYGEN IN GRAMS
         FFOX1(I)=PHOTOTF(I+1)*16.*1.6726E-24*XO(I)
C     FFOX1(I)=PHIOX*XMSO*XO(I)
CALEX SOURCE COEF?         
         FFOX2(I)=2.2E-11*XMSO*XO(I)/XMSH
         FFOX3(I)=-2.5E-11*XH(I)
         FFOX4(I)=-1.53E-12*XN2(I)
         FFOX5(I)=-2.73E-12*XN2(I)
         FFOX6(I)=-2.82E-11*XO2(I)
         FFHYD1(I)=2.5E-11*XMSH*XH(I)/XMSO
         FFHYD2(I)=-2.2E-11*XO(I)
         FFHE1(I)=PHIHE*XMSHE*XHE(I)
         FFHE2(I)=-(1.10E-9*XO2(I)+1.60E-9*XN2(I))


CALEX CL=COLLISION COEF, CF=collision freq ?         
CALEX cf O+ and N2         
         CFOXN2(I)=6.82E-10*XN2(I)
CALEX cf of O+ and O2         
         CFOXO2(I)=6.64E-10*XO2(I)
CALEX cf of O+ and He         
         CFOXHE(I)=1.32E-10*XHE(I)
         
CALEX cl for O+ and O
         CLOXO(I)=3.67E-11*XO(I)
CALEX seems like cl for O+ and H but I think the constant 4.63
CALEX is wrong. see nagy p.99         
         CLOXH(I)=4.63E-12*XH(I)
CALEX cl for O+ and H+         
         CLOXHD(I)=0.077*17.**1.5/XMSH
CALEX cl for O+ and He+         
         CLOXHL(I)=0.14*5.**1.5/XMSHE
CALEX cl for O and electrons??        
         CLOXEL(I)=1.86E-3/XMSE

CALEX cf of H+ and N2  
         CFHN2(I)=3.36E-9*XN2(I)
CALEX cf of H+ and O2
         CFHO2(I)=3.20E-9*XO2(I)
CALEX cf of H+ and He         
         CFHHE(I)=1.06E-9*XHE(I)

CALEX cl for ion neutral collsion       
         CLHO(I)=6.61E-11*XO(I)
         CLHH(I)=2.65E-10*XH(I)
CALEX cl for ion ion collision         
         CLHOX(I)=1.23*17.**1.5/XMSO
         CLHHL(I)=1.14*5.**1.5/XMSHE
         CLHEL(I)=2.97E-2/XMSE

CALEX cf of He+ and N2
         CFHEN2(I)=1.60E-9*XN2(I)
CALEX cf of He+ and O2         
         CFHEO2(I)=1.53E-9*XO2(I)
CALEX cl for He+ and He         
         CLHEHE(I)=8.73E-11*XHE(I)
CALEX cf of He+ and O         
         CFHEO(I)=1.01E-9*XO(I)
CALEX cf of He+ and H         
         CFHEH(I)=4.71E-10*XH(I)
         
CALEX these are cl for ion ion collisions
         CLHEOX(I)=0.57*5.**1.5/XMSO
         CLHEHD(I)=0.28*5.**1.5/XMSH
         CLHEEL(I)=7.43E-3/XMSE
        
CALEX cf of electrons and N2? this doesnt match nagy p 97 though?? 
         CFELN2(I)=3.42E-9*XN2(I)
         CFELO2(I)=3.25E-9*XO2(I)
         CFELHE(I)=1.18E-9*XHE(I)
         CFELO(I)=2.26E-9*XO(I)
         CFELH(I)=2.11E-9*XH(I)
CALEX similar to nagy p95 but I cant find the connection???         
         CLELOX(I)=54.5/XMSO
         CLELHL(I)=54.5/XMSHE
         CLELHD(I)=54.5/XMSH
         GRAVTY(I)=-3.96E20/RAD(I)**2
         Centrifugal(I)=wHorizontal**2.0*RAD(I)
 20   CONTINUE
cAlex
      NEXP=3
      WRITE(*,*) CURR(1),RAD(1),NEXP
cendAlex      
!      CRRX=CURR(1)*RAD(1)**NEXP
!      DO 25 I=1,NDIM
!         CURR(I)=CRRX/RAD(I)**NEXP
! 25   CONTINUE
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
      CTOXN2=3.*RGASO*XMSO/(XMSO+28.*XAMU)
      CTOXO2=3.*RGASO*XMSO/(XMSO+32.*XAMU)
      CTOXO=3.*RGASO*XMSO/(XMSO+XMSO)
      CTOXHE=3.*RGASO*XMSO/(XMSO+XMSHE)
      CTOXH=3.*RGASO*XMSO/(XMSO+XMSH)
      CTOXHD=CTOXH
      CTOXHL=CTOXHE
      CTOXEL=3.*RGASO*XMSO/(XMSO+XMSE)
      CTHN2=3.*RGASH*XMSH/(XMSH+28.*XAMU)
      CTHO2=3.*RGASH*XMSH/(XMSH+32.*XAMU)
      CTHHE=3.*RGASH*XMSH/(XMSH+XMSHE)
      CTHO=3.*RGASH *XMSH/(XMSH+XMSO)
      CTHH=3.*RGASH*XMSH/(XMSH+XMSH)
      CTHOX=CTHO
      CTHHL=CTHHE
      CTHEL=3.*RGASH*XMSH/(XMSH+XMSE)
      CTHEN2=3.*RGASHE*XMSHE/(XMSHE+28.*XAMU)
      CTHEO2=3.*RGASHE*XMSHE/(XMSHE+32.*XAMU)
      CTHEHE=3.*RGASHE*XMSHE/(XMSHE+XMSHE)
      CTHEO=3.*RGASHE*XMSHE/(XMSHE+XMSO)
      CTHEH=3.*RGASHE*XMSHE/(XMSHE+XMSH)
      CTHEOX=CTHEO
      CTHEHD=CTHEH
      CTHEEL=3.*RGASHE*XMSHE/(XMSHE+XMSE)
      CTELN2=3.*RGASE*XMSE/(28.*XAMU+XMSE)
      CTELO2=3.*RGASE*XMSE/(32.*XAMU+XMSE)
      CTELHE=3.*RGASE*XMSE/(XMSE+XMSHE)
      CTELO=3.*RGASE*XMSE/(XMSE+XMSO)
      CTELH=3.*RGASE*XMSE/(XMSE+XMSH)
      CTELOX=CTELO
      CTELHL=CTELHE
      CTELHD=CTELH
CALEX CMOXN2 = M_{N2}/(M_o+M_{N2}) see nagy p.83
      CMOXN2=28.*XAMU/(XMSO+28.*XAMU)
      CMOXO2=32.*XAMU/(XMSO+32.*XAMU)
      CMOXO=XMSO/(XMSO+XMSO)
      CMOXHE=XMSHE/(XMSO+XMSHE)
      CMOXH=XMSH/(XMSO+XMSH)
      CMOXHD=CMOXH
      CMOXHL=CMOXHE
      CMOXEL=XMSE/(XMSE+XMSO)
      CMHN2=28.*XAMU/(XMSH+28.*XAMU)
      CMHO2=32.*XAMU/(XMSH+32.*XAMU)
      CMHHE=XMSHE/(XMSH+XMSHE)
      CMHO=XMSO/(XMSH+XMSO)
      CMHH=XMSH/(XMSH+XMSH)
      CMHOX=CMHO
      CMHHL=CMHHE
      CMHEL=XMSE/(XMSE+XMSH)
      CMHEN2=28.*XAMU/(XMSHE+28.*XAMU)
      CMHEO2=32.*XAMU/(XMSHE+32.*XAMU)
      CMHEHE=XMSHE/(XMSHE+XMSHE)
      CMHEO=XMSO/(XMSHE+XMSO)
      CMHEH=XMSH/(XMSHE+XMSH)
      CMHEOX=CMHEO
      CMHEHD=CMHEH
      CMHEEL=XMSE/(XMSE+XMSHE)
      CMELN2=28.*XAMU/(XMSE+28.*XAMU)
      CMELO2=32.*XAMU/(XMSE+32.*XAMU)
      CMELHE=XMSHE/(XMSE+XMSHE)
      CMELO=XMSO/(XMSE+XMSO)
      CMELH=XMSH/(XMSE+XMSH)
      CMELOX=CMELO
      CMELHL=CMELHE
      CMELHD=CMELH

C     ALEX(10/11/04): 
C     TRY SETTING THE PLASMA PARAMETERS HERE TO THE SURFACE VALUES

      if(IsRestart) RETURN

      do K=1,NDIM
c         DHYD(K)=DSURFH
c         UOXYG(K)=0
c         POXYG(K)=PSURFOH
c         DOXYG(K)=DSURFOH
c         TOXYG(K)=TSURFO
c         UHYD(K)=0
c         PHYD(K)=PSURFHH
c         THYD(K)=TSURFH
c         UHEL(K)=0
c         PHEL(K)=PSURHEH
c         DHEL(K)=DSURHEH
c         THEL(K)=TSURHE
c         DELECT(K)=DSURFEH
c         UELECT(K)=0
c         PELECT(K)=PSURFEH
c         TELECT(K)=TSURFE


         DHYD(K)=1.*DSURFH*exp(-1.*real(K)/90.)
         UOXYG(K)=20.
         POXYG(K)=1.*PSURFO*exp(-1.*real(K)/90.)
         DOXYG(K)=1.*DSURFO*exp(-1.*real(K)/90.)
c         TOXYG(K)=TSURFO
         TOXYG(K)=POXYG(K)/(1.38e-16*DOXYG(K)/xmso)
         UHYD(K)=40.
         PHYD(K)=1.*PSURFH*exp(-1.*real(K)/90.)
c         THYD(K)=TSURFH*exp(-1.*real(K)/9000.)
         THYD(K)=PHYD(K)/(1.38e-16*DHYD(K)/xmsh)
         UHEL(K)=4.
         PHEL(K)=1.*PSURHE*exp(-1.*real(K)/90.)
         DHEL(K)=1.*DSURHE*exp(-1.*real(K)/90.)
c         THEL(K)=TSURHE
         THEL(K)=PHEL(K)/(1.38e-16*DHEL(K)/xmshe)
         DELECT(K)=1.*DSURFE*exp(-1.*real(K)/90.)
         UELECT(K)=1.
         PELECT(K)=1.*PSURFE*exp(-1.*real(K)/90.)
c         TELECT(K)=TSURFE
         TELECT(K)=PELECT(K)/(1.38e-16*DELECT(K)/xmse)

      enddo

      RETURN
      END
