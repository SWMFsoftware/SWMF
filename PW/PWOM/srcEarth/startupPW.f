
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     ALEX(10/11/04): I THINK THAT THIS SUBROUTINE INITIALIZES THE GRID AND
C     AND SETS CONSTANTS. IT MUST BE INITIALIZED EVERY TIME THE CODE RUNS
C     EVEN IF RESTARTING FROM A FILE.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE STRT
C
C

      use ModCommonVariables
      use ModGlow, ONLY: get_ionization
C
      use ModConst ,ONLY: cBoltzmann
      use ModPWOM  ,ONLY: UseAurora,UseIndicies, UseIE
      use ModAurora,ONLY: get_aurora,AuroralIonRateO_C
      use ModPwTime,ONLY: CurrentTime,StartTime,iStartTime,
     &                    Hour_,Minute_,Second_
      use ModIndicesInterfaces
      use ModNumConst, ONLY: cDegToRad
      use ModLatLon,   ONLY: convert_lat_lon
C     
      CurrentTime=StartTime+Time
      NPT1=14
      NPT2=16
      NPT3=30
      NPT4=35
      NPT5=60
      NPT6=70
1     FORMAT(5X,I5)
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     DEFINE THE GAS SPECIFIC HEAT AT CONSTANT PRESSURE (CP),          C
C           THE AVERAGE MOLECULAR MASS (AVMASS), THE ACTUAL            C
C           GAS CONSTANT (RGAS) AND THE SPECIFIC HEAT RATIO (GAMMA)    C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C Gas constant = k_Boltzmann/AMU
      RGAS=8.314E7
C Adiabatic index
      GAMMA=5./3.
C AMU in gramms
      XAMU=1.6606655E-24
C Mass of atomic O in gramms
      Mass_I(Ion1_)=15.994*XAMU
C Mass of atomic H in gramms
      Mass_I(Ion2_)=1.00797*XAMU
C Mass of atomic He in gramms
      Mass_I(Ion3_)=4.0026*XAMU
C Mass of electron in gramms
      Mass_I(Ion4_)=9.109534E-28
C Relative mass of atomic O to electron
      MassElecIon_I(Ion1_)=Mass_I(Ion4_)/Mass_I(Ion1_)
C Relative mass of atomic H to electron
      MassElecIon_I(Ion2_)=Mass_I(Ion4_)/Mass_I(Ion2_)
C Relative mass of atomic He to electron
      MassElecIon_I(Ion3_)=Mass_I(Ion4_)/Mass_I(Ion3_)
C kB/m_O
      RGAS_I(Ion1_)=RGAS*XAMU/Mass_I(Ion1_)
C kB/m_H
      RGAS_I(Ion2_)=RGAS*XAMU/Mass_I(Ion2_)
C kB/m_He
      RGAS_I(Ion3_)=RGAS*XAMU/Mass_I(Ion3_)
C kB/m_e
      RGAS_I(Ion4_)=RGAS*XAMU/Mass_I(Ion4_)
      GMIN1=GAMMA-1.
      GMIN2=GMIN1/2.
      GPL1=GAMMA+1.
      GPL2=GPL1/2.
      GM12=GMIN1/GAMMA/2.
      GRAR=GAMMA/GMIN2
      GREC=1./GAMMA
      CPO=GAMMA*RGAS_I(Ion1_)/GMIN1
      CPH=GAMMA*RGAS_I(Ion2_)/GMIN1
      CPHE=GAMMA*RGAS_I(Ion3_)/GMIN1
      CPE=GAMMA*RGAS_I(Ion4_)/GMIN1
      CVO=RGAS_I(Ion1_)/GMIN1
      CVH=RGAS_I(Ion2_)/GMIN1
      CVHE=RGAS_I(Ion3_)/GMIN1
      CVE=RGAS_I(Ion4_)/GMIN1

C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     READ FIELD ALIGNED CURRENT DENSITY AT LOWEST GRID POINT          C
C        (UNIT=mAMPERE/M**2)                                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C      READ (5,2) CURR(1)
C      CURR(1)=0.E-6

      call set_vertical_grid

      CURR(1)=2.998E2*CURR(1)
!      CURTIM=150.
!      CURTIM0=500.
      CURRMN=CURR(1)*(RAD(1)/(RAD(1)-DRBND))**NEXP
      CURRMX=CURR(1)*(RAD(1)/(RAD(NDIM)+DRBND))**NEXP
      
      do k=2,nDim
         CURR(k)=CURR(1)*(RAD(1)/(RAD(k)))**NEXP
      enddo


C      SGN1=1.
C      IF (CURR(1).LT.0.) SGN1=-1.
C      IF (ABS(CURR(1)).LT.1.E-4) SGN1=0.
C      CURRMX=SGN1*0.2998*(RAD(1)/(RAD(NDIM)+DRBND))**NEXP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     DEFINE THE NEUTRAL ATMOSPHERE MODEL                              C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C IF GEOMAGNETIC COORDINATES ARE SPECIFIED, SET IART TO 1 !!!!!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     DEFINE DATE                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C SOLAR MAXIMUM WINTER
C      IYD=80360
C      F107A=180.
C      F107=180.
C SOLAR MAXIMUM SUMMER
C      IYD=80172
C      F107A=180.
C      F107=180.
C SOLAR MINIMUM WINTER
C      IYD=84360
C      F107A=60.
C      F107=60.
C SOLAR MINIMUM SUMMER
c      IYD=84172
c      F107A=60.
c      F107=60.
c      SEC=43200.
c      STL=12.
c      GMLAT=80.
c      GMLONG=0.
c      IART=1
c      GLAT=0.
c      GLONG=180.
c      IART=0
C FEB 20, 1990
CALEX IYD=year_day of year

!      IYD=76183
!      STL=12.
!      F107A=60.
!      F107=60.
!      IART=1
!      GMLONG=0.
!      GMLAT=80.
C END
!      CALL GGM_PLANET(IART,GLONG,GLAT,GMLONG,GMLAT)
      CALL convert_lat_lon(Time,GMLAT,GMLONG,GLAT,GLONG)
      
      if (.not.UseStaticAtmosphere) then
         iDay  = mod(IYD,1000)+floor(Time/24.0/3600.0)
         iYear = IYD/1000
         if(iDay > 364) then
            iDay = mod(iDay,364)
            iYear= iYear+1
            IYD = iYear*1000+iDay
         else
            IYD = iYear*1000+iDay
         endif
         
         SEC=mod(iStartTime(Hour_)*3600.0 + iStartTime(Minute_)*60.0 
     &        + iStartTime(Second_) + Time, 24.0*3600.0)
         !mod((STL-GLONG/15.)*3600.,24.*3600.)
      else
         SEC=0.0
      endif
      STL=SEC/3600+GLONG/15
      DO 49 I=1,7
c      AP(I)=50.
         AP(I)=4.  
49    CONTINUE 

!     if UseIndicies = T then fill AP array appropriatly 
      if (UseIndicies) then
         AP(:) = 0.0
!     (1) DAILY AP
         do iHour =-12,12
            TempTime = real(iHour)*3600.0 + CurrentTime
            call get_ap(TempTime, TempAp, iError)
            AP(1) = AP(1) + TempAp
         enddo
         AP(1) = AP(1)/25.0
!     (2) 3 HR AP INDEX FOR CURRENT TIME
         do iHour =-1,1
            TempTime = real(iHour)*3600.0 + CurrentTime
            call get_ap(TempTime, TempAp, iError)
            AP(2) = AP(2) + TempAp
         enddo
         AP(2) = AP(2)/3.0
!     (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
          do iHour =-4,-2
            TempTime = real(iHour)*3600.0 + CurrentTime
            call get_ap(TempTime, TempAp, iError)
            AP(3) = AP(3) + TempAp
         enddo
         AP(3) = AP(3)/3.0
!     (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
         do iHour =-7,-5
            TempTime = real(iHour)*3600.0 + CurrentTime
            call get_ap(TempTime, TempAp, iError)
            AP(4) = AP(4) + TempAp
         enddo
         AP(4) = AP(4)/3.0
!     (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
         do iHour =-10,-8
            TempTime = real(iHour)*3600.0 + CurrentTime
            call get_ap(TempTime, TempAp, iError)
            AP(5) = AP(5) + TempAp
         enddo
         AP(5) = AP(5)/3.0
!             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
!                    TO CURRENT TIME
         do iHour =-35,-12
            TempTime = real(iHour)*3600.0 + CurrentTime
            call get_ap(TempTime, TempAp, iError)
            AP(6) = AP(6) + TempAp
         enddo
         AP(6) = AP(6)/24.0
!             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 59 HRS PRIOR
!                    TO CURRENT TIME
         do iHour =-59,-36
            TempTime = real(iHour)*3600.0 + CurrentTime
            call get_ap(TempTime, TempAp, iError)
            AP(7) = AP(7) + TempAp
         enddo
         AP(7) = AP(7)/24.0
!     check for errors when reading ap
         if (iError /= 0) then
            write(*,*) 'PW_ERROR: get_ap failed in startupPW_planet'
            call con_stop()
         endif              

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                     
!            Update F107 for current time 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                        
         call get_f107(CurrentTime, TempF107, iError)
         F107 = TempF107
         if(iError /=0) then
            write(*,*) 'PW_ERROR: get_f107 failed in startupPW_planet'
            call con_stop()
         endif
         
         call get_f107a(CurrentTime, TempF107a, iError)
         F107A = TempF107a
         
         if(iError /=0) then
!write(*,*) 'PW_ERROR: get_f107a failed in startupPW_planet'
            call con_stop()
         endif
         
      endif
        
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C                                                                      C
      cBoltzmannCGS = 1.0e7*cBoltzmann
      DO 50 K=1,NDIM
         CALL MODATM(ALTD(K),XO2(K),XN2(K),XO(K),XH(K),XHE(K),XTN(K))
         NDensity_CI(k,O_ ) = XO (K)
         NDensity_CI(k,O2_) = XO2(K)
         NDensity_CI(k,N2_) = XN2(K)
         NDensity_CI(k,H_ ) = XH (K)
         NDensity_CI(k,He_) = XHE(K)
         
         NeutralPressure_C(k) = cBoltzmannCGS*sum(NDensity_CI(k,:))*XTN(k)
         
50    CONTINUE


      
      if(UseAurora) then
         call get_aurora(nDim,AltD(1:nDim),NDensity_CI(1:nDim,O_:N2_),
     &        NeutralPressure_C(1:nDim)*0.1)
      endif
      	
      call  get_ionization(nDim, AltD(1:nDim), IonRateO_C(1:nDim))		

!      do K =1,nDim
!         write(*,*) AltD(K),IonRateO_C(K),AuroralIonRateO_C(K)
!      enddo
      DO 1099 J = 1,40

 9999    FORMAT(2X,1PE15.3,2X,1PE15.3)
 1099 CONTINUE
C
      CALL STRT1
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     DEFINE TOPSIDE ELECTRON HEAT FLUX AND PARAMETRIC HEAT SOURCES    C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C      READ (5,2) ETOP,ELFXIN
c      ETOP=5.0E-3

      if (UseIE .and. GmLat < 85.0 .and. UseAuroralHeatFlux) then
         ETOP = max (EfluxIE/EfluxRef * EtopAurora, EtopMin)
      else
         ETOP = EtopMin
      endif
      
      if (UsePhotoElectronHeatFlux) then
         call SOLZEN (IYD, SEC, GLAT, GLONG, SZA)
         ETOP = ETOP+EtopPhotoElectrons*max(cos(SZA*cDegToRad),0.0)
      endif

      ELFXIN=0.
C
C      ELFXIN=9.
C
2     FORMAT(6X,1PE15.4)

C      READ (5,2) HEATI1,HEATI2,ELHEAT
      HEATI1=0.
      HEATI2=0.
C
C      HEATI1=1.0E-7
C
      HEATI2=5.E-6
C
      ELHEAT=0.

      HEATA1=3.5E7
      HEATA2=1.30E7
      HEATA3=1.5E8
      HEATS1=2.*1.0E7**2
      HEATS2=2.0E6
      HEATS3=2.*1.0E7**2
      DO 53 K=1,NDIM
      HEATX1=EXP(-(ALTD(K)-HEATA1)**2/HEATS1)
      HEATX2=EXP(-(ALTD(K)-HEATA2)/HEATS2-EXP(-(ALTD(K)-HEATA2)/HEATS2))
      HEATX3=EXP(-(ALTD(K)-HEATA3)**2/HEATS3)

c      QOXYG(K)=(HEATI1*HEATX1+HEATI2*HEATX2)/
c     #         (State_GV(K,RhoO_)+State_GV(K,RhoH_)+State_GV(K,RhoHe_))
c      QHYD(K)=QOXYG(K)/16.
c      QHEL(K)=QOXYG(K)/4.

      QOXYG(K)=0.
      QHYD(K)=QOXYG(K)/16.
      QHEL(K)=QOXYG(K)/4.


C
C      QOXYG(K)=0.
C
      QELECT(K)=ELHEAT*HEATX3


53    CONTINUE
!      DO 54 K=1,NDIM,10
!      WRITE (iUnitOutput,52) ALTD(K),QOXYG(K),QHYD(K),QHEL(K),QELECT(K)
52    FORMAT(5(1PE15.4))
!54    CONTINUE
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     READ STARTING TIME, TERMINATION TIME AND BASIC TIME STEP         C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C      READ (5,4) TIME,TMAX
CCC      TIME=0.
CCC      TMAX=1.E6

            

4     FORMAT(6X,2(1PE15.4))
!!!!      READ (iUnitInput,*) Dt
C      DT=1./10.
      DTX1=DT
      DTR1=DTX1/DRBND
      DTX2=DT*NTS
      DTR2=DTX2/DRBND

      H0=0.5/DRBND
      H1E1=1./DTX1
      H1O1=1./DTX1
      H1H1=1./DTX1
      H1E2=1./DTX2
      H1O2=1./DTX2
      H1H2=1./DTX2
      H2=1./DRBND/DRBND
      H3=0.5*H2
      H4=0.5*H0
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     DEFINE THE HEAT CONDUCTION PARAMETERS                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
CALEX terms with "surf" in them refer to surface values      
      HLPE0=GMIN1/RGAS_I(Ion4_)
      HLPE=1.23E-6*GMIN1/RGAS_I(Ion4_)
      HLPO=1.24E-8*(Mass_I(Ion4_)/Mass_I(Ion1_))*GMIN1/RGAS_I(Ion1_)
      HLPHE=2.5*(RGAS_I(Ion3_)**2*Mass_I(Ion3_))*GMIN1/RGAS_I(Ion3_)
      HLPH=4.96E-8*(Mass_I(Ion4_)/Mass_I(Ion2_))*GMIN1/RGAS_I(Ion2_)
      CALL MODATM(ALTMIN,XNO2,XNN2,XNO,XNH,XNHE,TEMP)
CALEX I believe these are heat conduction coef. at lower boundary
      CXHN2=3.36E-9*XNN2
      CXHO2=3.20E-9*XNO2
      CXHO=6.61E-11*XNO*SQRT(State_GV(0,Th_))*(1.-0.047*ALOG10(State_GV(0,Th_)))**2
      CXHOX=1.23*(State_GV(0,RhoO_)/Mass_I(Ion1_))/State_GV(0,Th_)**1.5
      CXHEN2=1.60E-9*XNN2
      CXHEO2=1.53E-9*XNO2
      CXHEHE=8.73E-11*XNHE*SQRT(State_GV(0,The_))*(1.-0.093*ALOG10(State_GV(0,The_)))**2
      CXHEO=1.01E-9*XNO
      CXHEH=4.71E-10*XNH
      CXHEOX=0.57*(State_GV(0,RhoO_)/Mass_I(Ion1_))/State_GV(0,The_)**1.5
      CXHEHD=0.28*(State_GV(0,RhoH_)/Mass_I(Ion2_))/State_GV(0,The_)**1.5
CALEX I believe these are heat conductivities      
      HeatCon_GI(0,Ion1_)=HLPO*(State_GV(0,RhoO_)/State_GV(0,RhoE_))*State_GV(0,To_)**2.5
      HeatCon_GI(0,Ion4_)=HLPE*State_GV(0,Te_)**2.5
      HeatCon_GI(0,Ion2_)=HLPH*(State_GV(0,RhoH_)/State_GV(0,RhoE_))*State_GV(0,Th_)**2.5
      HeatCon_GI(0,Ion2_)=HeatCon_GI(0,Ion2_)/(1.+(0.7692*(CXHN2+CXHO2)+1.0962*CXHO)/CXHOX)
      HeatCon_GI(0,Ion3_)=HLPHE*(State_GV(0,RhoHe_)/Mass_I(Ion3_))*State_GV(0,The_)
      HeatCon_GI(0,Ion3_)=HeatCon_GI(0,Ion3_)/(0.99*CXHEN2+0.99*CXHEO2+1.02*CXHEO+1.48*CXHEHE+
     $2.22*CXHEH+1.21*CXHEOX+2.23*CXHEHD)      
      CALL MODATM(ALTMAX,XNO2,XNN2,XNO,XNH,XNHE,TEMP)
      XTNMAX=TEMP
CALEX I believe these are heat conduction coef. at upper boundary
      CZHN2=3.36E-9*XNN2
      CZHO2=3.20E-9*XNO2
      CZHO=6.61E-11*XNO
      CZHOX=1.23*17.**1.5/Mass_I(Ion1_)
      CZHEN2=1.60E-9*XNN2
      CZHEO2=1.53E-9*XNO2
      CZHEHE=8.73E-11*XNHE
      CZHEO=1.01E-9*XNO
      CZHEH=4.71E-10*XNH
      CZHEOX=0.57*5.**1.5/Mass_I(Ion1_)
      CZHEHD=0.28*5.**1.5/Mass_I(Ion2_)
      ETOP=ETOP*DRBND/1.23E-6
      CALL PW_set_upper_bc
3     FORMAT(4X,I6)
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     PRINT INITIAL AND BOUNDARY PARAMETERS                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C      READ(5,3) NCNPRT

      
      NCNPRT=0

      CALL MODATM(ALTMIN,XNO2,XNN2,XNO,XNH,XNHE,XNT)
      CALL MODATM(ALTMAX,YNO2,YNN2,YNO,YNH,YNHE,YNT)
      DO 60 I=1,NDIM
      ALTD(I)=ALTD(I)/1.E5
60    CONTINUE
      ALTMIN=ALTMIN/1.E5
      ALTMAX=ALTMAX/1.E5
      ETOP1=ETOP*1.23E-6/DRBND
      CALL COLLIS(NDIM,State_GV(-1:nDim+2,:))

      CALL PW_CALC_EFIELD(nDim,State_GV(-1:nDim+2,:))

      !write log
      if (DoLog) then
      IF (NCNPRT.NE.0) GO TO 999
      WRITE(iUnitOutput,1005) NDIM
1005  FORMAT(1H1,5X,'NUMBER OF CELLS=',I4)
      WRITE(iUnitOutput,1020) NEXP
1020  FORMAT(5X,'NEXP=',I1)
      WRITE (iUnitOutput,1008) GAMMA,RGAS_I(Ion1_),CPO,CVO,RGAS_I(Ion3_),CPHE,CVHE,
     ;RGAS_I(Ion2_),CPH,CVH,RGAS_I(Ion4_),CPE,CVE
1008  FORMAT(5X,'GAMMA=',F4.2,/5X,'RGAS(OXYGEN)=',1PE10.4,7X,
     ;'CP(OXYGEN)=',1PE10.4,7X,'CV(OXYGEN)=',1PE10.4
     ;/5X,'RGAS(HELIUM)=',1PE10.4,7X,
     ;'CP(HELIUM)=',1PE10.4,7X,'CV(HELIUM)=',1PE10.4/5X,
     ;'RGAS(HYDROGEN)=',1PE10.4,5X,'CP(HYDROGEN)=',1PE10.4,5X,
     ;'CV(HYDROGEN)=',1PE10.4,/5X,'RGAS(ELECTRON)=',1PE10.4,5X
     ;,'CP(ELECTRON)=',1PE10.4,5X,'CV(ELECTRON)=',1PE10.4)
      WRITE (iUnitOutput,1023)
1023  FORMAT(1H0,5X,'LOWER BOUNDARY PLASMA PARAMETERS:')
      WRITE(iUnitOutput,1001)
1001  FORMAT(1H ,4X,'OXYGEN:')
      WRITE (iUnitOutput,1009) State_GV(0,uO_),State_GV(0,Po_),State_GV(0,RhoO_),State_GV(0,To_),SoundSpeed_GI(0,Ion1_)
1009  FORMAT(5X,'VELOCITY=',1PE11.4,3X,'PRESSURE=',1PE10.4,3X,
     ;'MASS DENSITY=',1PE10.4,3X,'TEMPERATURE=',1PE10.4,3X,
     ;'SOUND VELOCITY=',1PE10.4)
      WRITE(iUnitOutput,10021)
10021 FORMAT(1H ,4X,'HELIUM:')
      WRITE (iUnitOutput,1009) State_GV(0,uHe_),State_GV(0,pHe_),State_GV(0,RhoHe_),State_GV(0,The_),SoundSpeed_GI(0,Ion3_)
      WRITE(iUnitOutput,1002)
1002  FORMAT(1H ,4X,'HYDROGEN:')
      WRITE (iUnitOutput,1009) State_GV(0,uH_),State_GV(0,pH_),State_GV(0,RhoH_),State_GV(0,Th_),SoundSpeed_GI(0,Ion2_)
      WRITE(iUnitOutput,1003)
1003  FORMAT(1H ,4X,'ELECTRONS:')
      WRITE (iUnitOutput,1009) State_GV(0,uE_),State_GV(0,pE_),State_GV(0,RhoE_),State_GV(0,Te_),SoundSpeed_GI(0,Ion4_)
      WRITE (iUnitOutput,1027)
1027  FORMAT(1H0,5X,'UPPER BOUNDARY INITIAL PLASMA PARAMETERS:')
      WRITE (iUnitOutput,1004)
1004  FORMAT(1H ,4X,'OXYGEN:')
      WRITE (iUnitOutput,1009) State_GV(nDim+1,uO_),State_GV(nDim+1,pO_),
     &     State_GV(nDim+1,RhoO_),State_GV(nDim+1,To_),SoundSpeed_GI(nDim+1,Ion1_)
      WRITE (iUnitOutput,1088)
1088  FORMAT(1H ,4X,'HELIUM:')
      WRITE (iUnitOutput,1009) 
     &     State_GV(nDim+1,uHe_),State_GV(nDim+1,pHe_),State_GV(nDim+1,RhoHe_),
     &     State_GV(nDim+1,The_),SoundSpeed_GI(nDim+1,Ion3_)
      WRITE (iUnitOutput,1006)
1006  FORMAT(1H ,4X,'HYDROGEN:')
      WRITE (iUnitOutput,1009) State_GV(nDim+1,uH_),State_GV(nDim+1,pH_),
     &     State_GV(nDim+1,RhoH_),State_GV(nDim+1,Th_),SoundSpeed_GI(nDim+1,Ion2_)
      WRITE(iUnitOutput,1007)
1007  FORMAT(1H ,4X,'ELECTRONS:')
      WRITE (iUnitOutput,1009) State_GV(nDim+1,uE_),State_GV(nDim+1,pE_),
     &     State_GV(nDim+1,RhoE_),State_GV(nDim+1,Te_),SoundSpeed_GI(nDim+1,Ion4_)
      WRITE (iUnitOutput,1029) ETOP1
1029  FORMAT(1H0,5X,'TOPSIDE ELECTRON HEATING RATE:',1PE10.4,
     ;' ERGS/CM**3/SEC')
      WRITE (iUnitOutput,1050) 
1050  FORMAT(1H0,5X,'ENERGY COLLISION TERM COEFFICIENTS')
      WRITE (iUnitOutput,1051) CTOXN2,CTOXO2,CTOXO,CTOXHE,CTOXH,CTOXHD,CTOXHL,
     $CTOXEL,CTHEN2,CTHEO2,CTHEO,CTHEHE,CTHEH,CTHEOX,CTHEHD,CTHEEL,
     $CTHN2,CTHO2,CTHO,CTHHE,CTHH,CTHOX,CTHHL,CTHEL,CTELN2,CTELO2,
     $CTELO,CTELHE,CTELH,CTELOX,CTELHL,CTELHD
1051  FORMAT(1H0,5X,'CTOXN2=',1PE10.4,5X,'CTOXO2=',1PE10.4,4X,
     $'CTOXO=',1PE10.4/5X,'CTOXHE=',1PE10.4,5X,'CTOXH=',1PE10.4,5X,
     $'CTOXHD=',1PE10.4,5X,'CTOXHL=',1PE10.4,5X,'CTOXEL=',1PE10.4/
     $5X,'CTHEN2=',1PE10.4,5X,'CTHEO2=',1PE10.4,4X,
     $'CTHEO=',1PE10.4/5X,'CTHEHE=',1PE10.4,5X,'CTHEH=',1PE10.4,5X,
     $'CTHEOX=',1PE10.4,5X,'CTHEHD=',1PE10.4,5X,'CTHEEL=',1PE10.4/
     $5X,'CTHN2=',1PE10.4,6X,'CTHO2=',1PE10.4,5X,'CTHO=',1PE10.4/5X,
     $'CTHHE=',1PE10.4,6X,'CTHH=',1PE10.4,6X,'CTHOX=',1PE10.4,6X,
     $'CTHHL=',1PE10.4,6X,'CTHEL=',1PE10.4/5X,
     $'CTELN2=',1PE10.4,5X,'CTELO2=',1PE10.4,4X,
     $'CTELO=',1PE10.4/5X,'CTELHE=',1PE10.4,5X,'CTELH=',1PE10.4,5X,
     $'CTELOX=',1PE10.4,5X,'CTELHL=',1PE10.4,5X,'CTELHD=',1PE10.4,5X)
      WRITE (iUnitOutput,1052) CMOXN2,CMOXO2,CMOXO,CMOXHE,CMOXH,CMOXHD,CMOXHL,
     $CMOXEL,CMHEN2,CMHEO2,CMHEO,CMHEHE,CMHEH,CMHEOX,CMHEHD,CMHEEL,
     $CMHN2,CMHO2,CMHO,CMHHE,CMHH,CMHOX,CMHHL,CMHEL,CMELN2,CMELO2,
     $CMELO,CMELHE,CMELH,CMELOX,CMELHL,CMELHD
1052  FORMAT(1H0,5X,'CMOXN2=',1PE10.4,5X,'CMOXO2=',1PE10.4,4X,
     $'CMOXO=',1PE10.4/5X,'CMOXHE=',1PE10.4,5X,'CMOXH=',1PE10.4,5X,
     $'CMOXHD=',1PE10.4,5X,'CMOXHL=',1PE10.4,5X,'CMOXEL=',1PE10.4/
     $5X,'CMHEN2=',1PE10.4,5X,'CMHEO2=',1PE10.4,4X,
     $'CMHEO=',1PE10.4/5X,'CMHEHE=',1PE10.4,5X,'CMHEH=',1PE10.4,5X,
     $'CMHEOX=',1PE10.4,5X,'CMHEHD=',1PE10.4,5X,'CMHEEL=',1PE10.4/
     $5X,'CMHN2=',1PE10.4,6X,'CMHO2=',1PE10.4,5X,'CMHO=',1PE10.4/5X,
     $'CMHHE=',1PE10.4,6X,'CMHH=',1PE10.4,6X,'CMHOX=',1PE10.4,6X,
     $'CMHHL=',1PE10.4,6X,'CMHEL=',1PE10.4/5X,
     $'CMELN2=',1PE10.4,5X,'CMELO2=',1PE10.4,4X,
     $'CMELO=',1PE10.4/5X,'CMELHE=',1PE10.4,5X,'CMELH=',1PE10.4,5X,
     $'CMELOX=',1PE10.4,5X,'CMELHL=',1PE10.4,5X,'CMELHD=',1PE10.4,5X)
      WRITE (iUnitOutput,1053)
1053  FORMAT(1H0,5X,'HEAT CONDUCTION COEFFICIENTS AT UPPER BOUNDARY')
      WRITE (iUnitOutput,1054) CZHN2,CZHO2,CZHO,CZHOX,CZHEN2,CZHEO2,CZHEHE,
     $CZHEO,CZHEH,CZHEOX,CZHEHD,XTNMAX
1054  FORMAT(1H0,5X,'CZHN2=',1PE10.4,6X,'CZHO2=',1PE10.4,6X,
     $'CZHO=',1PE10.4,7X,'CZHOX=',1PE10.4/5X,'CZHEN2=',1PE10.4,5X,
     $'CZHEO2=',1PE10.4,5X,'CZHEHE=',1PE10.4,5X,'CZHEO=',1PE10.4/
     $5X,'CZHEH=',1PE10.4,6X,'CZHEOX=',1PE10.4,5X,'CZHEHD=',
     $1PE10.4/5X,'XTNMAX=',1PE10.4)
      WRITE (iUnitOutput,1012)
1012  FORMAT(1H1,45X,'NEUTRAL ATMOSPHERE NUMBER DENSITIES')
      WRITE(iUnitOutput,1014)
1014  FORMAT(16X,'ALT',13X,'N2',13X,'O2',15X,'O',14X,'HE',
     $15X,'H',15X,'T')
      K=0
      WRITE (iUnitOutput,1022) K, ALTMIN,XNN2,XNO2,XNO,XNHE,XNH,XNT
      NDMQ=NPT1
      IF (NDIM.LT.NPT2) NDMQ=NDIM
      DO 230 K=1,NDMQ
         WRITE(iUnitOutput,1022) K,ALTD(K),XN2(K),XO2(K),XO(K),XHE(K),
     $        XH(K),XTN(K)
 230  CONTINUE
      IF (NDIM.LT.NPT2) GO TO 290
      NDMQ=NPT3
      IF (NDIM.LT.NPT4) NDMQ=NDIM
      DO 240 K=NPT2,NDMQ,2
         WRITE(iUnitOutput,1022) K,ALTD(K),XN2(K),XO2(K),XO(K),XHE(K),
     $        XH(K),XTN(K)
 240  CONTINUE
      IF (NDIM.LT.NPT4) GO TO 290
      NDMQ=NPT5
      IF (NDIM.LT.NPT6) NDMQ=NDIM
      DO 250 K=NPT4,NDMQ,5
         WRITE(iUnitOutput,1022) K,ALTD(K),XN2(K),XO2(K),XO(K),XHE(K),
     $        XH(K),XTN(K)
 250  CONTINUE
      IF (NDIM.LT.NPT6) GO TO 290
      DO 260 K=NPT6,NDIM,10
         WRITE(iUnitOutput,1022) K,ALTD(K),XN2(K),XO2(K),XO(K),XHE(K),
     $        XH(K),XTN(K)
 260  CONTINUE
 290  CONTINUE
      K=NDIM+1
      WRITE (iUnitOutput,1022) K, ALTMAX,YNN2,YNO2,YNO,YNHE,YNH,YNT
      WRITE(iUnitOutput,1015)
1015  FORMAT(1H1,55X,'SOURCE COEFFICIENTS:')
      WRITE (iUnitOutput,1016)
1016  FORMAT(12X,'ALT',6X,'GRAVTY',5X,'FFOX1',5X,'FFOX2',5X,
     ;'FFOX3',5X,'FFOX4',5X,'FFOX5',5X,'FFOX6',4X,'FFHYD1',
     ;4X,'FFHYD2',4X,'FFHE1',5X,'FFHE2')
      NDMQ=NPT1
      IF (NDIM.LT.NPT2) NDMQ=NDIM
      DO 291 K=1,NDMQ
      WRITE(iUnitOutput,1019) K,ALTD(K),GRAVTY(K),FFOX1(K),FFOX2(K),
     $FFOX3(K),FFOX4(K),FFOX5(K),FFOX6(K),FFHYD1(K),FFHYD2(K),
     $FFHE1(K),FFHE2(K)
291   CONTINUE
      IF (NDIM.LT.NPT2) GO TO 295
      NDMQ=NPT3
      IF (NDIM.LT.NPT4) NDMQ=NDIM
      DO 292 K=NPT2,NDMQ,2
      WRITE(iUnitOutput,1019) K,ALTD(K),GRAVTY(K),FFOX1(K),FFOX2(K),
     $FFOX3(K),FFOX4(K),FFOX5(K),FFOX6(K),FFHYD1(K),FFHYD2(K),
     $FFHE1(K),FFHE2(K)
292   CONTINUE
      IF (NDIM.LT.NPT4) GO TO 295
      NDMQ=NPT5
      IF (NDIM.LT.NPT6) NDMQ=NDIM
      DO 293 K=NPT4,NDMQ,5
      WRITE(iUnitOutput,1019) K,ALTD(K),GRAVTY(K),FFOX1(K),FFOX2(K),
     $FFOX3(K),FFOX4(K),FFOX5(K),FFOX6(K),FFHYD1(K),FFHYD2(K),
     $FFHE1(K),FFHE2(K)
293   CONTINUE
      IF (NDIM.LT.NPT6) GO TO 295
      DO 294 K=NPT6,NDIM,10
      WRITE(iUnitOutput,1019) K,ALTD(K),GRAVTY(K),FFOX1(K),FFOX2(K),
     $FFOX3(K),FFOX4(K),FFOX5(K),FFOX6(K),FFHYD1(K),FFHYD2(K),
     $FFHE1(K),FFHE2(K)
294   CONTINUE
295   CONTINUE
      WRITE(iUnitOutput,1017)
1017  FORMAT(1H1,50X,'COLLISION COEFFICIENTS FOR OXYGEN')
      WRITE(iUnitOutput,1018)
1018  FORMAT(12X,'ALT',7X,'CLOXN2',6X,'CLOXO2',6X,'CLOXO',7X,
     ;'CLOXHE',6X,'CLOXH',7X,'CLOXHL',6X,'CLOXHD',6X,'CLOXEL')
      NDMQ=NPT1
      IF (NDIM.LT.NPT2) NDMQ=NDIM
      DO 530 K=1,NDMQ
      WRITE (iUnitOutput,1180) K,ALTD(K),CLOXN2(K),CLOXO2(K),CLOXO(K),
     $CLOXHE(K),CLOXH(K),CLOXHL(K),CLOXHD(K),CLOXEL(K)
530   CONTINUE
      IF (NDIM.LT.NPT2) GO TO 590
      NDMQ=NPT3
      IF (NDIM.LT.NPT4) NDMQ=NDIM
      DO 540 K=NPT2,NDMQ,2
      WRITE (iUnitOutput,1180) K,ALTD(K),CLOXN2(K),CLOXO2(K),CLOXO(K),
     $CLOXHE(K),CLOXH(K),CLOXHL(K),CLOXHD(K),CLOXEL(K)
540   CONTINUE
      IF (NDIM.LT.NPT4) GO TO 590
      NDMQ=NPT5
      IF (NDIM.LT.NPT6) NDMQ=NDIM
      DO 550 K=NPT4,NDMQ,5
      WRITE (iUnitOutput,1180) K,ALTD(K),CLOXN2(K),CLOXO2(K),CLOXO(K),
     $CLOXHE(K),CLOXH(K),CLOXHL(K),CLOXHD(K),CLOXEL(K)
550   CONTINUE
      IF (NDIM.LT.NPT6) GO TO 590
      DO 560 K=NPT6,NDIM,10
      WRITE (iUnitOutput,1180) K,ALTD(K),CLOXN2(K),CLOXO2(K),CLOXO(K),
     $CLOXHE(K),CLOXH(K),CLOXHL(K),CLOXHD(K),CLOXEL(K)
560   CONTINUE
590   CONTINUE
      WRITE (iUnitOutput,1217)
1217  FORMAT(1H1,50X,'COLLISION COEFFICIENTS FOR HELIUM')
      WRITE(iUnitOutput,1218)
1218  FORMAT(12X,'ALT',7X,'CLHEN2',6X,'CLHEO2',6X,'CLHEO',7X,
     ;'CLHEHE',6X,'CLHEH',7X,'CLHEOX',6X,'CLHEHD',6X,'CLHEEL')
      NDMQ=NPT1
      IF (NDIM.LT.NPT2) NDMQ=NDIM
      DO 1230 K=1,NDMQ
      WRITE (iUnitOutput,1180) K,ALTD(K),CLHEN2(K),CLHEO2(K),CLHEO(K),
     $CLHEHE(K),CLHEH(K),CLHEOX(K),CLHEHD(K),CLHEEL(K)
1230  CONTINUE
      IF (NDIM.LT.NPT2) GO TO 1290
      NDMQ=NPT3
      IF (NDIM.LT.NPT4) NDMQ=NDIM
      DO 1240 K=NPT2,NDMQ,2
      WRITE (iUnitOutput,1180) K,ALTD(K),CLHEN2(K),CLHEO2(K),CLHEO(K),
     $CLHEHE(K),CLHEH(K),CLHEOX(K),CLHEHD(K),CLHEEL(K)
1240  CONTINUE
      IF (NDIM.LT.NPT4) GO TO 1290
      NDMQ=NPT5
      IF (NDIM.LT.NPT6) NDMQ=NDIM
      DO 1250 K=NPT4,NDMQ,5
      WRITE (iUnitOutput,1180) K,ALTD(K),CLHEN2(K),CLHEO2(K),CLHEO(K),
     $CLHEHE(K),CLHEH(K),CLHEOX(K),CLHEHD(K),CLHEEL(K)
1250  CONTINUE
      IF (NDIM.LT.NPT6) GO TO 1290
      DO 1260 K=NPT6,NDIM,10
      WRITE (iUnitOutput,1180) K,ALTD(K),CLHEN2(K),CLHEO2(K),CLHEO(K),
     $CLHEHE(K),CLHEH(K),CLHEOX(K),CLHEHD(K),CLHEEL(K)
1260  CONTINUE
1290  CONTINUE
      WRITE (iUnitOutput,2217)
2217  FORMAT(1H1,50X,'COLLISION COEFFICIENTS FOR HYDROGEN')
      WRITE(iUnitOutput,2218)
2218  FORMAT(12X,'ALT',7X,'CLHN2',7X,'CLHO2',7X,'CLHO',8X,
     ;'CLHHE',7X,'CLHH',8X,'CLHOX',7X,'CLHHL',7X,'CLHEL')
      NDMQ=NPT1
      IF (NDIM.LT.NPT2) NDMQ=NDIM
      DO 2230 K=1,NDMQ
      WRITE (iUnitOutput,1180) K,ALTD(K),CLHN2(K),CLHO2(K),CLHO(K),
     $CLHHE(K),CLHH(K),CLHOX(K),CLHHL(K),CLHEL(K)
2230  CONTINUE
      IF (NDIM.LT.NPT2) GO TO 2290
      NDMQ=NPT3
      IF (NDIM.LT.NPT4) NDMQ=NDIM
      DO 2240 K=NPT2,NDMQ,2
      WRITE (iUnitOutput,1180) K,ALTD(K),CLHN2(K),CLHO2(K),CLHO(K),
     $CLHHE(K),CLHH(K),CLHOX(K),CLHHL(K),CLHEL(K)
2240  CONTINUE
      IF (NDIM.LT.NPT4) GO TO 2290
      NDMQ=NPT5
      IF (NDIM.LT.NPT6) NDMQ=NDIM
      DO 2250 K=NPT4,NDMQ,5
      WRITE (iUnitOutput,1180) K,ALTD(K),CLHN2(K),CLHO2(K),CLHO(K),
     $CLHHE(K),CLHH(K),CLHOX(K),CLHHL(K),CLHEL(K)
2250  CONTINUE
      IF (NDIM.LT.NPT6) GO TO 2290
      DO 2260 K=NPT6,NDIM,10
      WRITE (iUnitOutput,1180) K,ALTD(K),CLHN2(K),CLHO2(K),CLHO(K),
     $CLHHE(K),CLHH(K),CLHOX(K),CLHHL(K),CLHEL(K)
2260  CONTINUE
2290  CONTINUE
      WRITE (iUnitOutput,3217)
3217  FORMAT(1H1,50X,'COLLISION COEFFICIENTS FOR ELECTRONS')
      WRITE(iUnitOutput,3218)
3218  FORMAT(12X,'ALT',6X,'CLELN2',6X,'CLELO2',6X,'CLELO',7X,
     ;'CLELHE',6X,'CLELH',7X,'CLELOX',6X,'CLELHL',6X,'CLELHD')
      NDMQ=NPT1
      IF (NDIM.LT.NPT2) NDMQ=NDIM
      DO 3230 K=1,NDMQ
      WRITE (iUnitOutput,1180) K,ALTD(K),CLELN2(K),CLELO2(K),CLELO(K),
     $CLELHE(K),CLELH(K),CLELOX(K),CLELHL(K),CLELHD(K)
3230  CONTINUE
      IF (NDIM.LT.NPT2) GO TO 3290
      NDMQ=NPT3
      IF (NDIM.LT.NPT4) NDMQ=NDIM
      DO 3240 K=NPT2,NDMQ,2
      WRITE (iUnitOutput,1180) K,ALTD(K),CLELN2(K),CLELO2(K),CLELO(K),
     $CLELHE(K),CLELH(K),CLELOX(K),CLELHL(K),CLELHD(K)
3240  CONTINUE
      IF (NDIM.LT.NPT4) GO TO 3290
      NDMQ=NPT5
      IF (NDIM.LT.NPT6) NDMQ=NDIM
      DO 3250 K=NPT4,NDMQ,5
      WRITE (iUnitOutput,1180) K,ALTD(K),CLELN2(K),CLELO2(K),CLELO(K),
     $CLELHE(K),CLELH(K),CLELOX(K),CLELHL(K),CLELHD(K)
3250  CONTINUE
      IF (NDIM.LT.NPT6) GO TO 3290
      DO 3260 K=NPT6,NDIM,10
      WRITE (iUnitOutput,1180) K,ALTD(K),CLELN2(K),CLELO2(K),CLELO(K),
     $CLELHE(K),CLELH(K),CLELOX(K),CLELHL(K),CLELHD(K)
3260  CONTINUE
3290  CONTINUE
      WRITE(iUnitOutput,1047)
1047  FORMAT(1H1,50X,'COLLISION FREQUENCIES FOR OXYGEN')
      WRITE(iUnitOutput,1048)
1048  FORMAT(12X,'ALT',7X,'CFOXN2',6X,'CFOXO2',6X,'CFOXO',7X,
     ;'CFOXHE',6X,'CFOXH',7X,'CFOXHL',6X,'CFOXHD',6X,'CFOXEL')
      NDMQ=NPT1
      IF (NDIM.LT.NPT2) NDMQ=NDIM
      DO 531 K=1,NDMQ
      WRITE (iUnitOutput,1180) K,ALTD(K),CollisionFreq_IIC(Ion1_,Neutral4_,K),
     &  CollisionFreq_IIC(Ion1_,Neutral3_,K),
     &  CollisionFreq_IIC(Ion1_,Neutral1_,K),
     &  CollisionFreq_IIC(Ion1_,Neutral5_,K),
     &  CollisionFreq_IIC(Ion1_,Neutral2_,K),
     &  CollisionFreq_IIC(Ion1_,Ion3_,K),
     &  CollisionFreq_IIC(Ion1_,Ion2_,K),
     &  CollisionFreq_IIC(Ion1_,Ion4_,K)
531   CONTINUE
      IF (NDIM.LT.NPT2) GO TO 591
      NDMQ=NPT3
      IF (NDIM.LT.NPT4) NDMQ=NDIM
      DO 541 K=NPT2,NDMQ,2
      WRITE (iUnitOutput,1180) K,ALTD(K),CollisionFreq_IIC(Ion1_,Neutral4_,K),
     &  CollisionFreq_IIC(Ion1_,Neutral3_,K),
     &  CollisionFreq_IIC(Ion1_,Neutral1_,K),
     &  CollisionFreq_IIC(Ion1_,Neutral5_,K),
     &  CollisionFreq_IIC(Ion1_,Neutral2_,K),
     &  CollisionFreq_IIC(Ion1_,Ion3_,K),
     &  CollisionFreq_IIC(Ion1_,Ion2_,K),
     &  CollisionFreq_IIC(Ion1_,Ion4_,K)
541   CONTINUE
      IF (NDIM.LT.NPT4) GO TO 591
      NDMQ=NPT5
      IF (NDIM.LT.NPT6) NDMQ=NDIM
      DO 551 K=NPT4,NDMQ,5
      WRITE (iUnitOutput,1180) K,ALTD(K),CollisionFreq_IIC(Ion1_,Neutral4_,K),
     &  CollisionFreq_IIC(Ion1_,Neutral3_,K),
     &  CollisionFreq_IIC(Ion1_,Neutral1_,K),
     &  CollisionFreq_IIC(Ion1_,Neutral5_,K),
     &  CollisionFreq_IIC(Ion1_,Neutral2_,K),
     &  CollisionFreq_IIC(Ion1_,Ion3_,K),
     &  CollisionFreq_IIC(Ion1_,Ion2_,K),
     &  CollisionFreq_IIC(Ion1_,Ion4_,K)

551   CONTINUE
      IF (NDIM.LT.NPT6) GO TO 591
      DO 561 K=NPT6,NDIM,10
      WRITE (iUnitOutput,1180) K,ALTD(K),CollisionFreq_IIC(Ion1_,Neutral4_,K),
     &  CollisionFreq_IIC(Ion1_,Neutral3_,K),
     &  CollisionFreq_IIC(Ion1_,Neutral1_,K),
     &  CollisionFreq_IIC(Ion1_,Neutral5_,K),
     &  CollisionFreq_IIC(Ion1_,Neutral2_,K),
     &  CollisionFreq_IIC(Ion1_,Ion3_,K),
     &  CollisionFreq_IIC(Ion1_,Ion2_,K),
     &  CollisionFreq_IIC(Ion1_,Ion4_,K)

561   CONTINUE
591   CONTINUE
      WRITE (iUnitOutput,1247)
1247  FORMAT(1H1,50X,'COLLISION FREQUENCIES FOR HELIUM')
      WRITE(iUnitOutput,1248)
1248  FORMAT(12X,'ALT',7X,'CFHEN2',6X,'CFHEO2',6X,'CFHEO',7X,
     ;'CFHEHE',6X,'CFHEH',7X,'CFHEOX',6X,'CFHEHD',6X,'CFHEEL')
      NDMQ=NPT1
      IF (NDIM.LT.NPT2) NDMQ=NDIM
      DO 1231 K=1,NDMQ
      WRITE (iUnitOutput,1180) K,ALTD(K),
     &  CollisionFreq_IIC(Ion3_,Neutral4_,K),
     &  CollisionFreq_IIC(Ion3_,Neutral3_,K),
     &  CollisionFreq_IIC(Ion3_,Neutral1_,K),
     &  CollisionFreq_IIC(Ion3_,Neutral5_,K),
     &  CollisionFreq_IIC(Ion3_,Neutral2_,K),
     &  CollisionFreq_IIC(Ion3_,Ion1_,K),
     &  CollisionFreq_IIC(Ion3_,Ion2_,K),
     &  CollisionFreq_IIC(Ion3_,Ion4_,K)

1231  CONTINUE
      IF (NDIM.LT.NPT2) GO TO 1291
      NDMQ=NPT3
      IF (NDIM.LT.NPT4) NDMQ=NDIM
      DO 1241 K=NPT2,NDMQ,2
      WRITE (iUnitOutput,1180) K,ALTD(K),
     &  CollisionFreq_IIC(Ion3_,Neutral4_,K),
     &  CollisionFreq_IIC(Ion3_,Neutral3_,K),
     &  CollisionFreq_IIC(Ion3_,Neutral1_,K),
     &  CollisionFreq_IIC(Ion3_,Neutral5_,K),
     &  CollisionFreq_IIC(Ion3_,Neutral2_,K),
     &  CollisionFreq_IIC(Ion3_,Ion1_,K),
     &  CollisionFreq_IIC(Ion3_,Ion2_,K),
     &  CollisionFreq_IIC(Ion3_,Ion4_,K)

1241  CONTINUE
      IF (NDIM.LT.NPT4) GO TO 1291
      NDMQ=NPT5
      IF (NDIM.LT.NPT6) NDMQ=NDIM
      DO 1251 K=NPT4,NDMQ,5
      WRITE (iUnitOutput,1180) K,ALTD(K),
     &  CollisionFreq_IIC(Ion3_,Neutral4_,K),
     &  CollisionFreq_IIC(Ion3_,Neutral3_,K),
     &  CollisionFreq_IIC(Ion3_,Neutral1_,K),
     &  CollisionFreq_IIC(Ion3_,Neutral5_,K),
     &  CollisionFreq_IIC(Ion3_,Neutral2_,K),
     &  CollisionFreq_IIC(Ion3_,Ion1_,K),
     &  CollisionFreq_IIC(Ion3_,Ion2_,K),
     &  CollisionFreq_IIC(Ion3_,Ion4_,K)

1251  CONTINUE
      IF (NDIM.LT.NPT6) GO TO 1291
      DO 1261 K=NPT6,NDIM,10
      WRITE (iUnitOutput,1180) K,ALTD(K),
     &  CollisionFreq_IIC(Ion3_,Neutral4_,K),
     &  CollisionFreq_IIC(Ion3_,Neutral3_,K),
     &  CollisionFreq_IIC(Ion3_,Neutral1_,K),
     &  CollisionFreq_IIC(Ion3_,Neutral5_,K),
     &  CollisionFreq_IIC(Ion3_,Neutral2_,K),
     &  CollisionFreq_IIC(Ion3_,Ion1_,K),
     &  CollisionFreq_IIC(Ion3_,Ion2_,K),
     &  CollisionFreq_IIC(Ion3_,Ion4_,K)

1261  CONTINUE
1291  CONTINUE
      WRITE (iUnitOutput,2247)
2247  FORMAT(1H1,50X,'COLLISION FREQUENCIES FOR HYDROGEN')
      WRITE(iUnitOutput,2248)
2248  FORMAT(12X,'ALT',7X,'CFHN2',7X,'CFHO2',7X,'CFHO',8X,
     ;'CFHHE',7X,'CFHH',8X,'CFHOX',7X,'CFHHL',7X,'CFHEL')
      NDMQ=NPT1
      IF (NDIM.LT.NPT2) NDMQ=NDIM
      DO 2231 K=1,NDMQ
      WRITE (iUnitOutput,1180) K,ALTD(K),
     &  CollisionFreq_IIC(Ion2_,Neutral4_,K),
     &  CollisionFreq_IIC(Ion2_,Neutral3_,K),
     &  CollisionFreq_IIC(Ion2_,Neutral1_,K),
     &  CollisionFreq_IIC(Ion2_,Neutral5_,K),
     &  CollisionFreq_IIC(Ion2_,Neutral2_,K),
     &  CollisionFreq_IIC(Ion2_,Ion1_,K),
     &  CollisionFreq_IIC(Ion2_,Ion3_,K),
     &  CollisionFreq_IIC(Ion2_,Ion4_,K)

2231  CONTINUE
      IF (NDIM.LT.NPT2) GO TO 2291
      NDMQ=NPT3
      IF (NDIM.LT.NPT4) NDMQ=NDIM
      DO 2241 K=NPT2,NDMQ,2
      WRITE (iUnitOutput,1180) K,ALTD(K),
     &  CollisionFreq_IIC(Ion2_,Neutral4_,K),
     &  CollisionFreq_IIC(Ion2_,Neutral3_,K),
     &  CollisionFreq_IIC(Ion2_,Neutral1_,K),
     &  CollisionFreq_IIC(Ion2_,Neutral5_,K),
     &  CollisionFreq_IIC(Ion2_,Neutral2_,K),
     &  CollisionFreq_IIC(Ion2_,Ion1_,K),
     &  CollisionFreq_IIC(Ion2_,Ion3_,K),
     &  CollisionFreq_IIC(Ion2_,Ion4_,K)

2241  CONTINUE
      IF (NDIM.LT.NPT4) GO TO 2291
      NDMQ=NPT5
      IF (NDIM.LT.NPT6) NDMQ=NDIM
      DO 2251 K=NPT4,NDMQ,5
      WRITE (iUnitOutput,1180) K,ALTD(K),
     &  CollisionFreq_IIC(Ion2_,Neutral4_,K),
     &  CollisionFreq_IIC(Ion2_,Neutral3_,K),
     &  CollisionFreq_IIC(Ion2_,Neutral1_,K),
     &  CollisionFreq_IIC(Ion2_,Neutral5_,K),
     &  CollisionFreq_IIC(Ion2_,Neutral2_,K),
     &  CollisionFreq_IIC(Ion2_,Ion1_,K),
     &  CollisionFreq_IIC(Ion2_,Ion3_,K),
     &  CollisionFreq_IIC(Ion2_,Ion4_,K)

2251  CONTINUE
      IF (NDIM.LT.NPT6) GO TO 2291
      DO 2261 K=NPT6,NDIM,10
      WRITE (iUnitOutput,1180) K,ALTD(K),
     &  CollisionFreq_IIC(Ion2_,Neutral4_,K),
     &  CollisionFreq_IIC(Ion2_,Neutral3_,K),
     &  CollisionFreq_IIC(Ion2_,Neutral1_,K),
     &  CollisionFreq_IIC(Ion2_,Neutral5_,K),
     &  CollisionFreq_IIC(Ion2_,Neutral2_,K),
     &  CollisionFreq_IIC(Ion2_,Ion1_,K),
     &  CollisionFreq_IIC(Ion2_,Ion3_,K),
     &  CollisionFreq_IIC(Ion2_,Ion4_,K)


2261  CONTINUE
2291  CONTINUE
      WRITE (iUnitOutput,3247)
3247  FORMAT(1H1,50X,'COLLISION FREQUENCIES FOR ELECTRONS')
      WRITE(iUnitOutput,3248)
3248  FORMAT(12X,'ALT',7X,'CFELN2',6X,'CFELO2',6X,'CFELO',7X,
     ;'CFELHE',6X,'CFELH',7X,'CFELOX',6X,'CFELHL',6X,'CFELHD')
      NDMQ=NPT1
      IF (NDIM.LT.NPT2) NDMQ=NDIM
      DO 3231 K=1,NDMQ
      WRITE (iUnitOutput,1180) K,ALTD(K),
     &  CollisionFreq_IIC(Ion4_,Neutral4_,K),
     &  CollisionFreq_IIC(Ion4_,Neutral3_,K),
     &  CollisionFreq_IIC(Ion4_,Neutral1_,K),
     &  CollisionFreq_IIC(Ion4_,Neutral5_,K),
     &  CollisionFreq_IIC(Ion4_,Neutral2_,K),
     &  CollisionFreq_IIC(Ion4_,Ion1_,K),
     &  CollisionFreq_IIC(Ion4_,Ion3_,K),
     &  CollisionFreq_IIC(Ion4_,Ion2_,K)

3231  CONTINUE
      IF (NDIM.LT.NPT2) GO TO 3291
      NDMQ=NPT3
      IF (NDIM.LT.NPT4) NDMQ=NDIM
      DO 3241 K=NPT2,NDMQ,2
      WRITE (iUnitOutput,1180) K,ALTD(K),
     &  CollisionFreq_IIC(Ion4_,Neutral4_,K),
     &  CollisionFreq_IIC(Ion4_,Neutral3_,K),
     &  CollisionFreq_IIC(Ion4_,Neutral1_,K),
     &  CollisionFreq_IIC(Ion4_,Neutral5_,K),
     &  CollisionFreq_IIC(Ion4_,Neutral2_,K),
     &  CollisionFreq_IIC(Ion4_,Ion1_,K),
     &  CollisionFreq_IIC(Ion4_,Ion3_,K),
     &  CollisionFreq_IIC(Ion4_,Ion2_,K)


3241  CONTINUE
      IF (NDIM.LT.NPT4) GO TO 3291
      NDMQ=NPT5
      IF (NDIM.LT.NPT6) NDMQ=NDIM
      DO 3251 K=NPT4,NDMQ,5
      WRITE (iUnitOutput,1180) K,ALTD(K),
     &  CollisionFreq_IIC(Ion4_,Neutral4_,K),
     &  CollisionFreq_IIC(Ion4_,Neutral3_,K),
     &  CollisionFreq_IIC(Ion4_,Neutral1_,K),
     &  CollisionFreq_IIC(Ion4_,Neutral5_,K),
     &  CollisionFreq_IIC(Ion4_,Neutral2_,K),
     &  CollisionFreq_IIC(Ion4_,Ion1_,K),
     &  CollisionFreq_IIC(Ion4_,Ion3_,K),
     &  CollisionFreq_IIC(Ion4_,Ion2_,K)


3251  CONTINUE
      IF (NDIM.LT.NPT6) GO TO 3291
      DO 3261 K=NPT6,NDIM,10
      WRITE (iUnitOutput,1180) K,ALTD(K),
     &  CollisionFreq_IIC(Ion4_,Neutral4_,K),
     &  CollisionFreq_IIC(Ion4_,Neutral3_,K),
     &  CollisionFreq_IIC(Ion4_,Neutral1_,K),
     &  CollisionFreq_IIC(Ion4_,Neutral5_,K),
     &  CollisionFreq_IIC(Ion4_,Neutral2_,K),
     &  CollisionFreq_IIC(Ion4_,Ion1_,K),
     &  CollisionFreq_IIC(Ion4_,Ion3_,K),
     &  CollisionFreq_IIC(Ion4_,Ion2_,K)


3261  CONTINUE
3291  CONTINUE
      WRITE (iUnitOutput,1013)
1013  FORMAT(1H1,45X,'INITIAL OXYGEN PARAMETERS')
      WRITE(iUnitOutput,1021)
1021  FORMAT(16X,'ALT',10X,'VELOCITY',8X,'MACH NO',9X,'DENSITY',9X,
     ;'PRESSURE',6X,'TEMPERATURE',/)
      K=0
      XM=State_GV(0,uO_)/SoundSpeed_GI(0,Ion1_)
      DNS1=State_GV(0,RhoO_)/Mass_I(Ion1_)
      WRITE(iUnitOutput,1022) K,ALTMIN,State_GV(0,uO_),XM,DNS1,State_GV(0,Po_),State_GV(0,To_)
      NDMQ=NPT1
      IF (NDIM.LT.NPT2) NDMQ=NDIM
      DO 630 K=1,NDMQ
      US=SQRT(GAMMA*State_GV(K,pO_)/State_GV(K,RhoO_))
      XM=State_GV(K,uO_)/US
      DNS1=State_GV(K,RhoO_)/Mass_I(Ion1_)
      WRITE(iUnitOutput,1022) K,ALTD(K),State_GV(K,uO_),XM,DNS1,State_GV(K,pO_),State_GV(K,To_)
630   CONTINUE
      IF (NDIM.LT.NPT2) GO TO 690
      NDMQ=NPT3
      IF (NDIM.LT.NPT4) NDMQ=NDIM
      DO 640 K=NPT2,NDMQ,2
      US=SQRT(GAMMA*State_GV(K,pO_)/State_GV(K,RhoO_))
      XM=State_GV(K,uO_)/US
      DNS1=State_GV(K,RhoO_)/Mass_I(Ion1_)
      WRITE(iUnitOutput,1022) K,ALTD(K),State_GV(K,uO_),XM,DNS1,State_GV(K,pO_),State_GV(K,To_)
640   CONTINUE
      IF (NDIM.LT.NPT4) GO TO 690
      NDMQ=NPT5
      IF (NDIM.LT.NPT6) NDMQ=NDIM
      DO 650 K=NPT4,NDMQ,5
      US=SQRT(GAMMA*State_GV(K,pO_)/State_GV(K,RhoO_))
      XM=State_GV(K,uO_)/US
      DNS1=State_GV(K,RhoO_)/Mass_I(Ion1_)
      WRITE(iUnitOutput,1022) K,ALTD(K),State_GV(K,uO_),XM,DNS1,State_GV(K,pO_),State_GV(K,To_)
650   CONTINUE
      IF (NDIM.LT.NPT6) GO TO 690
      DO 660 K=NPT6,NDIM,10
      US=SQRT(GAMMA*State_GV(K,pO_)/State_GV(K,RhoO_))
      XM=State_GV(K,uO_)/US
      DNS1=State_GV(K,RhoO_)/Mass_I(Ion1_)
      WRITE(iUnitOutput,1022) K,ALTD(K),State_GV(K,uO_),XM,DNS1,State_GV(K,pO_),State_GV(K,To_)
660   CONTINUE
690   CONTINUE
      K=NDIM1
      XM=State_GV(nDim+1,uO_)/SoundSpeed_GI(nDim+1,Ion1_)
      DNS1=State_GV(nDim+1,RhoO_)/Mass_I(Ion1_)
      WRITE(iUnitOutput,1022) K,ALTMAX,State_GV(nDim+1,uO_),XM,DNS1,State_GV(nDim+1,pO_),State_GV(nDim+1,To_)
      WRITE (iUnitOutput,1055)
1055  FORMAT(1H1,45X,'INITIAL HELIUM PARAMETERS')
      WRITE(iUnitOutput,1056)
1056  FORMAT(16X,'ALT',10X,'VELOCITY',8X,'MACH NO',9X,'DENSITY',9X,
     ;'PRESSURE',6X,'TEMPERATURE',/)
      K=0
      XM=State_GV(0,uHe_)/SoundSpeed_GI(0,Ion3_)
      DNS1=State_GV(0,RhoHe_)/Mass_I(Ion3_)
      WRITE(iUnitOutput,1022) K,ALTMIN,State_GV(0,uHe_),XM,DNS1,State_GV(0,pHe_),State_GV(0,The_)
      NDMQ=NPT1
      IF (NDIM.LT.NPT2) NDMQ=NDIM
      DO 639 K=1,NDMQ
      US=SQRT(GAMMA*State_GV(K,pHe_)/State_GV(K,RhoHe_))
      XM=State_GV(K,uHe_)/US
      DNS1=State_GV(K,RhoHe_)/Mass_I(Ion3_)
      WRITE(iUnitOutput,1022) K,ALTD(K),State_GV(K,uHe_),XM,DNS1,State_GV(K,pHe_),State_GV(K,The_)
639   CONTINUE
      IF (NDIM.LT.NPT2) GO TO 699
      NDMQ=NPT3
      IF (NDIM.LT.NPT4) NDMQ=NDIM
      DO 649 K=NPT2,NDMQ,2
      US=SQRT(GAMMA*State_GV(K,pHe_)/State_GV(K,RhoHe_))
      XM=State_GV(K,uHe_)/US
      DNS1=State_GV(K,RhoHe_)/Mass_I(Ion3_)
      WRITE(iUnitOutput,1022) K,ALTD(K),State_GV(K,uHe_),XM,DNS1,State_GV(K,pHe_),State_GV(K,The_)
649   CONTINUE
      IF (NDIM.LT.NPT4) GO TO 699
      NDMQ=NPT5
      IF (NDIM.LT.NPT6) NDMQ=NDIM
      DO 659 K=NPT4,NDMQ,5
      US=SQRT(GAMMA*State_GV(K,pHe_)/State_GV(K,RhoHe_))
      XM=State_GV(K,uHe_)/US
      DNS1=State_GV(K,RhoHe_)/Mass_I(Ion3_)
      WRITE(iUnitOutput,1022) K,ALTD(K),State_GV(K,uHe_),XM,DNS1,State_GV(K,pHe_),State_GV(K,The_)
659   CONTINUE
      IF (NDIM.LT.NPT6) GO TO 699
      DO 669 K=NPT6,NDIM,10
      US=SQRT(GAMMA*State_GV(K,pHe_)/State_GV(K,RhoHe_))
      XM=State_GV(K,uHe_)/US
      DNS1=State_GV(K,RhoHe_)/Mass_I(Ion3_)
      WRITE(iUnitOutput,1022) K,ALTD(K),State_GV(K,uHe_),XM,DNS1,State_GV(K,pHe_),State_GV(K,The_)
669   CONTINUE
699   CONTINUE
      K=NDIM1
      XM=State_GV(nDim+1,uHe_)/SoundSpeed_GI(nDim+1,Ion3_)
      DNS1=State_GV(nDim+1,RhoHe_)/Mass_I(Ion3_)
      WRITE(iUnitOutput,1022) K,ALTMAX,State_GV(nDim+1,uHe_),XM,DNS1,State_GV(nDim+1,pHe_),State_GV(nDim+1,The_)
      WRITE (iUnitOutput,1010)
1010  FORMAT(1H1,45X,'INITIAL HYDROGEN PARAMETERS')
      WRITE(iUnitOutput,1021)
      K=0
      XM=State_GV(0,uH_)/SoundSpeed_GI(0,Ion2_)
      DNS1=State_GV(0,RhoH_)/Mass_I(Ion2_)
      WRITE(iUnitOutput,1022) K,ALTMIN,State_GV(0,uH_),XM,DNS1,State_GV(0,pH_),State_GV(0,Th_)
      NDMQ=NPT1
      IF (NDIM.LT.NPT2) NDMQ=NDIM
      DO 730 K=1,NDMQ
      US=SQRT(GAMMA*State_GV(K,pH_)/State_GV(K,RhoH_))
      XM=State_GV(K,uH_)/US
      DNS1=State_GV(K,RhoH_)/Mass_I(Ion2_)
      WRITE(iUnitOutput,1022) K,ALTD(K),State_GV(K,uH_),XM,DNS1,State_GV(K,pH_),State_GV(K,Th_)
730   CONTINUE
      IF (NDIM.LT.NPT2) GO TO 790
      NDMQ=NPT3
      IF (NDIM.LT.NPT4) NDMQ=NDIM
      DO 740 K=NPT2,NDMQ,2
      US=SQRT(GAMMA*State_GV(K,pH_)/State_GV(K,RhoH_))
      XM=State_GV(K,uH_)/US
      DNS1=State_GV(K,RhoH_)/Mass_I(Ion2_)
      WRITE(iUnitOutput,1022) K,ALTD(K),State_GV(K,uH_),XM,DNS1,State_GV(K,pH_),State_GV(K,Th_)
740   CONTINUE
      IF (NDIM.LT.NPT4) GO TO 790
      NDMQ=NPT5
      IF (NDIM.LT.NPT6) NDMQ=NDIM
      DO 750 K=NPT4,NDMQ,5
      US=SQRT(GAMMA*State_GV(K,pH_)/State_GV(K,RhoH_))
      XM=State_GV(K,uH_)/US
      DNS1=State_GV(K,RhoH_)/Mass_I(Ion2_)
      WRITE(iUnitOutput,1022) K,ALTD(K),State_GV(K,uH_),XM,DNS1,State_GV(K,pH_),State_GV(K,Th_)
750   CONTINUE
      IF (NDIM.LT.NPT6) GO TO 790
      DO 760 K=NPT6,NDIM,10
      US=SQRT(GAMMA*State_GV(K,pH_)/State_GV(K,RhoH_))
      XM=State_GV(K,uH_)/US
      DNS1=State_GV(K,RhoH_)/Mass_I(Ion2_)
      WRITE(iUnitOutput,1022) K,ALTD(K),State_GV(K,uH_),XM,DNS1,State_GV(K,pH_),State_GV(K,Th_)
760   CONTINUE
790   CONTINUE
      K=NDIM1
      XM=State_GV(nDim+1,uH_)/SoundSpeed_GI(nDim+1,Ion2_)
      DNS1=State_GV(nDim+1,RhoH_)/Mass_I(Ion2_)
      WRITE(iUnitOutput,1022) K,ALTMAX,State_GV(nDim+1,uH_),XM,DNS1,State_GV(nDim+1,pH_),State_GV(nDim+1,Th_)
      WRITE (iUnitOutput,1011)
1011  FORMAT(1H1,45X,'INITIAL ELECTRON PARAMETERS')
      WRITE(iUnitOutput,1021)
      K=0
      XM=State_GV(0,uE_)/SoundSpeed_GI(0,Ion4_)
      DNS1=State_GV(0,RhoE_)/Mass_I(Ion4_)
      WRITE(iUnitOutput,1022) K,ALTMIN,State_GV(0,uE_),XM,DNS1,State_GV(0,pE_),State_GV(0,Te_)
      NDMQ=NPT1
      IF (NDIM.LT.NPT2) NDMQ=NDIM
      DO 830 K=1,NDMQ
      US=SQRT(GAMMA*State_GV(K,pE_)/State_GV(K,RhoE_))
      XM=State_GV(K,uE_)/US
      DNS1=State_GV(K,RhoE_)/Mass_I(Ion4_)
      WRITE(iUnitOutput,1022) K,ALTD(K),State_GV(K,uE_),XM,DNS1,State_GV(K,pE_),State_GV(K,Te_)
830   CONTINUE
      IF (NDIM.LT.NPT2) GO TO 890
      NDMQ=NPT3
      IF (NDIM.LT.NPT4) NDMQ=NDIM
      DO 840 K=NPT2,NDMQ,2
      US=SQRT(GAMMA*State_GV(K,pE_)/State_GV(K,RhoE_))
      XM=State_GV(K,uE_)/US
      DNS1=State_GV(K,RhoE_)/Mass_I(Ion4_)
      WRITE(iUnitOutput,1022) K,ALTD(K),State_GV(K,uE_),XM,DNS1,State_GV(K,pE_),State_GV(K,Te_)
840   CONTINUE
      IF (NDIM.LT.NPT4) GO TO 890
      NDMQ=NPT5
      IF (NDIM.LT.NPT6) NDMQ=NDIM
      DO 850 K=NPT4,NDMQ,5
      US=SQRT(GAMMA*State_GV(K,pE_)/State_GV(K,RhoE_))
      XM=State_GV(K,uE_)/US
      DNS1=State_GV(K,RhoE_)/Mass_I(Ion4_)
      WRITE(iUnitOutput,1022) K,ALTD(K),State_GV(K,uE_),XM,DNS1,State_GV(K,pE_),State_GV(K,Te_)
850   CONTINUE
      IF (NDIM.LT.NPT6) GO TO 890
      DO 860 K=NPT6,NDIM,10
      US=SQRT(GAMMA*State_GV(K,pE_)/State_GV(K,RhoE_))
      XM=State_GV(K,uE_)/US
      DNS1=State_GV(K,RhoE_)/Mass_I(Ion4_)
      WRITE(iUnitOutput,1022) K,ALTD(K),State_GV(K,uE_),XM,DNS1,State_GV(K,pE_),State_GV(K,Te_)
860   CONTINUE
890   CONTINUE
      K=NDIM1
      XM=State_GV(nDim+1,uE_)/SoundSpeed_GI(nDim+1,Ion4_)
      DNS1=State_GV(nDim+1,RhoE_)/Mass_I(Ion4_)
      WRITE(iUnitOutput,1022) K,ALTMAX,State_GV(nDim+1,uE_),XM,DNS1,State_GV(nDim+1,pE_),State_GV(nDim+1,Te_)
      WRITE(iUnitOutput,1024)
1024  FORMAT(1H1,40X,'INITIAL ELECTRIC FIELD AND SOURCE PARAMETERS')
      WRITE(iUnitOutput,1025)
1025  FORMAT(13X,'ALT',6X,'EFIELD',6X,'FCLSNO',6X,'ECLSNO',6X,'FCLSHE',
     ;6X,'ECLSHE',6X,'FCLSNH',6X,'ECLSNH',6X,'FCLSNE',6X,'ECLSNE')
      NDMQ=NPT1
      IF (NDIM.LT.NPT2) NDMQ=NDIM
      DO 930 K=1,NDMQ
      WRITE(iUnitOutput,1026) K,ALTD(K),EFIELD(K),Source_CV(K,uO_),Source_CV(K,pO_),Source_CV(K,uHe_),
     ;Source_CV(K,pHe_),Source_CV(K,uH_),Source_CV(K,pH_),Source_CV(K,uE_),Source_CV(K,pE_)
930   CONTINUE
      IF (NDIM.LT.NPT2) GO TO 990
      NDMQ=NPT3
      IF (NDIM.LT.NPT4) NDMQ=NDIM
      DO 940 K=NPT2,NDMQ,2
      WRITE(iUnitOutput,1026) K,ALTD(K),EFIELD(K),Source_CV(K,uO_),Source_CV(K,pO_),Source_CV(K,uHe_),
     ;Source_CV(K,pHe_),Source_CV(K,uH_),Source_CV(K,pH_),Source_CV(K,uE_),Source_CV(K,pE_)
940   CONTINUE
      IF (NDIM.LT.NPT4) GO TO 990
      NDMQ=NPT5
      IF (NDIM.LT.NPT6) NDMQ=NDIM
      DO 950 K=NPT4,NDMQ,5
      WRITE(iUnitOutput,1026) K,ALTD(K),EFIELD(K),Source_CV(K,uO_),Source_CV(K,pO_),Source_CV(K,uHe_),
     ;Source_CV(K,pHe_),Source_CV(K,uH_),Source_CV(K,pH_),Source_CV(K,uE_),Source_CV(K,pE_)
950   CONTINUE
      IF (NDIM.LT.NPT6) GO TO 990
      DO 960 K=NPT6,NDIM,10
      WRITE(iUnitOutput,1026) K,ALTD(K),EFIELD(K),Source_CV(K,uO_),Source_CV(K,pO_),Source_CV(K,uHe_),
     ;Source_CV(K,pHe_),Source_CV(K,uH_),Source_CV(K,pH_),Source_CV(K,uE_),Source_CV(K,pE_)
960   CONTINUE
990   CONTINUE
999   CONTINUE
      WRITE (iUnitOutput,2013)
2013  FORMAT(1H1,45X,'HEAT CONDUCTIVITIES')
      WRITE(iUnitOutput,2021)
2021  FORMAT(16X,'ALT',10X,'OXYGEN',10X,'HELIUM',9X,'HYDROGEN',9X,
     ;'ELECTRONS'/)
      K=0
      WRITE(iUnitOutput,1022) 
     & K,ALTMIN,HeatCon_GI(0,Ion1_),HeatCon_GI(0,Ion3_),HeatCon_GI(0,Ion2_),HeatCon_GI(0,Ion4_)
      NDMQ=NPT1
      IF (NDIM.LT.NPT2) NDMQ=NDIM
      DO 2630 K=1,NDMQ
      WRITE(iUnitOutput,1022) 
     & K,ALTD(K),HeatCon_GI(K,Ion1_),HeatCon_GI(K,Ion3_),HeatCon_GI(K,Ion2_),HeatCon_GI(K,Ion4_)
2630  CONTINUE
      IF (NDIM.LT.NPT2) GO TO 2690
      NDMQ=NPT3
      IF (NDIM.LT.NPT4) NDMQ=NDIM
      DO 2640 K=NPT2,NDMQ,2
      WRITE(iUnitOutput,1022) 
     & K,ALTD(K),HeatCon_GI(K,Ion1_),HeatCon_GI(K,Ion3_),HeatCon_GI(K,Ion2_),HeatCon_GI(K,Ion4_)
2640  CONTINUE
      IF (NDIM.LT.NPT4) GO TO 2690
      NDMQ=NPT5
      IF (NDIM.LT.NPT6) NDMQ=NDIM
      DO 2650 K=NPT4,NDMQ,5
      WRITE(iUnitOutput,1022) 
     & K,ALTD(K),HeatCon_GI(K,Ion1_),HeatCon_GI(K,Ion3_),HeatCon_GI(K,Ion2_),HeatCon_GI(K,Ion4_)
2650   CONTINUE
      IF (NDIM.LT.NPT6) GO TO 2690
      DO 2660 K=NPT6,NDIM,10
      WRITE(iUnitOutput,1022) 
     & K,ALTD(K),HeatCon_GI(K,Ion1_),HeatCon_GI(K,Ion3_),HeatCon_GI(K,Ion2_),HeatCon_GI(K,Ion4_)
2660   CONTINUE
2690   CONTINUE
      K=NDIM1
      WRITE(iUnitOutput,1022) 
     & K,ALTMAX,HeatCon_GI(nDim+1,Ion1_),HeatCon_GI(nDim+1,Ion3_),HeatCon_GI(nDim+1,Ion2_),HeatCon_GI(nDim+1,Ion4_)
1019  FORMAT(3X,I3,0PF10.2,2X,11(1PE10.2))
1022  FORMAT(3X,I3,0PF14.2,2X,6(1PE16.5))
1026  FORMAT(3X,I3,0PF11.2,2X,9(1PE12.4))
1180  FORMAT(3X,I3,0PF10.2,2X,8(1PE12.4))
1028  FORMAT(3X,I3,0PF14.2,2X,2(1PE16.5))
1031  FORMAT(3X,I3,0PF14.2,2X,4(1PE16.5))

      endif
      RETURN
      END

      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     ALEX(10/11/04): I THINK THIS IS A COLD START ROUTINE FOR DEALING WITH
C     A LACK OF RESTART FILE. BASICALLY THE VELOCITY IS SET TO ZERO AND 
C     THE OTHER GRID VALUES ARE SET TO VALUES THAT WON'T CRASH.
C     (11/9/04) I also see that source terms and collision terms are set
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      

      SUBROUTINE STRT1
      use ModCommonVariables
      use ModPWOM,ONLY: IsRestart
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
      State_GV(-1:0,uO_) =0.
      State_GV(-1:0,uh_) =0
      State_GV(-1:0,uhe_)=0.
      State_GV(-1:0,ue_) =0.
      CALL MODATM(ALTMIN-DrBnd,XNO2,XNN2,XNO,XNH,XNHE,TEMP)
      State_GV(-1,To_) =TEMP
      State_GV(-1,Th_) =TEMP
      State_GV(-1,The_)=TEMP
      State_GV(-1,Te_) =TEMP
      TMP1=1.-1.290E-3*TEMP+6.233E-7*TEMP**2
      TMP2=1.-9.149E-4*TEMP+4.228E-7*TEMP**2-6.790E-11*TEMP**3+
     $     4.225E-15*TEMP**4
      State_GV(-1,RhoO_) =IonRateO_C(1)*Mass_I(Ion1_)*XNO/(1.53E-12*XNN2*TMP1+2.82E-11*XNO2*TMP2)
      State_GV(-1,RhoHe_)=PHIHE*Mass_I(Ion3_)*XNHE/(1.10E-9*XNO2+1.60E-9*XNN2)
      State_GV(-1,RhoH_) =1.136*(XNH/XNO)*(Mass_I(Ion2_)/Mass_I(Ion1_))*State_GV(-1,RhoO_)

      CALL MODATM(ALTMIN,XNO2,XNN2,XNO,XNH,XNHE,TEMP)
      State_GV(0,To_) =TEMP
      State_GV(0,Th_) =TEMP
      State_GV(0,The_)=TEMP
      State_GV(0,Te_) =TEMP

      TMP1=1.-1.290E-3*TEMP+6.233E-7*TEMP**2
      TMP2=1.-9.149E-4*TEMP+4.228E-7*TEMP**2-6.790E-11*TEMP**3+
     $     4.225E-15*TEMP**4
C     State_GV(0,RhoO_)=PHIOX*Mass_I(Ion1_)*XNO/(1.53E-12*XNN2*TMP1+2.82E-11*XNO2*TMP2)
      State_GV(0,RhoO_) =IonRateO_C(1)*Mass_I(Ion1_)*XNO/(1.53E-12*XNN2*TMP1+2.82E-11*XNO2*TMP2)
      State_GV(0,RhoHe_)=PHIHE*Mass_I(Ion3_)*XNHE/(1.10E-9*XNO2+1.60E-9*XNN2)
      State_GV(0,RhoH_) =1.136*(XNH/XNO)*(Mass_I(Ion2_)/Mass_I(Ion1_))*State_GV(-0,RhoO_)

      State_GV(-1:0,RhoE_) =0.0 
      do iIon=1,nIon-1
         State_GV(-1:0,RhoE_) = 
     &        State_GV(-1:0,RhoE_)+MassElecIon_I(iIon)*State_GV(-1:0,iRho_I(iIon))
      enddo
      
!      State_GV(-1:0,uO_)=State_GV(1,RhoO_)*State_GV(1,uO_)/State_GV(-1:0,RhoO_)
!      State_GV(-1:0,uH_)=State_GV(1,RhoH_)*State_GV(1,uH_)/State_GV(-1:0,RhoH_)
!      State_GV(-1:0,uHe_)=State_GV(1,RhoHe_)*State_GV(1,uHe_)/State_GV(-1:0,RhoHe_)
!      State_GV(-1:0,ue_) = 0.0     

      State_GV(-1:0,pO_) =RGAS_I(Ion1_)*State_GV(-1:0,To_) *State_GV(-1:0,RhoO_)
      State_GV(-1:0,pH_) =RGAS_I(Ion2_)*State_GV(-1:0,Th_) *State_GV(-1:0,RhoH_)
      State_GV(-1:0,pHe_)=RGAS_I(Ion3_)*State_GV(-1:0,The_)*State_GV(-1:0,RhoHe_)
      State_GV(-1:0,pE_) =RGAS_I(Ion4_)*State_GV(-1:0,Te_) *State_GV(-1:0,RhoE_)
      
      do iIon=1,nIon
         SoundSpeed_GI(0,iIon) = 
     &        SQRT(
     &        GAMMA*State_GV(0,iP_I(iIon)) / State_GV(0,iRho_I(iIon))
     &        )
      enddo
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
      DO 20 I=1,NDIM
C     FFOX1(1) BASED ON PHOTOIONIZATION FREQ FROM SUB GLOWEX;MASS OXYGEN IN GRAMS
         if (I < nDim) then 
            FFOX1(I)=IonRateO_C(I+1)*16.*1.6726E-24*XO(I)
         else
            FFOX1(I)=IonRateO_C(I)*16.*1.6726E-24*XO(I)
         endif
C     FFOX1(I)=PHIOX*Mass_I(Ion1_)*XO(I)
CALEX SOURCE COEF?         
         FFOX2(I)=2.2E-11*Mass_I(Ion1_)*XO(I)/Mass_I(Ion2_)
         FFOX3(I)=-2.5E-11*XH(I)
         FFOX4(I)=-1.53E-12*XN2(I)
         FFOX5(I)=-2.73E-12*XN2(I)
         FFOX6(I)=-2.82E-11*XO2(I)
         FFHYD1(I)=2.5E-11*Mass_I(Ion2_)*XH(I)/Mass_I(Ion1_)
         FFHYD2(I)=-2.2E-11*XO(I)
         FFHE1(I)=PHIHE*Mass_I(Ion3_)*XHE(I)
         FFHE2(I)=-(1.10E-9*XO2(I)+1.60E-9*XN2(I))


CALEX CL=COLLISION COEF, CF=collision freq ?         
CALEX cf O+ and N2         
         !CFOXN2(I)=6.82E-10*XN2(I)
         CollisionFreq_IIC(Ion1_,Neutral4_,I)=6.82E-10*XN2(I)
CALEX cf of O+ and O2         
         !CFOXO2(I)=6.64E-10*XO2(I)
         CollisionFreq_IIC(Ion1_,Neutral3_,I)=6.64E-10*XO2(I)
CALEX cf of O+ and He         
         !CFOXHE(I)=1.32E-10*XHE(I)
         CollisionFreq_IIC(Ion1_,Neutral5_,I)=1.32E-10*XHE(I)
CALEX cl for O+ and O
         CLOXO(I)=3.67E-11*XO(I)
CALEX seems like cl for O+ and H but I think the constant 4.63
CALEX is wrong. see nagy p.99         
         CLOXH(I)=4.63E-12*XH(I)
CALEX cl for O+ and H+         
         CLOXHD(I)=0.077*17.**1.5/Mass_I(Ion2_)
CALEX cl for O+ and He+         
         CLOXHL(I)=0.14*5.**1.5/Mass_I(Ion3_)
CALEX cl for O and electrons??        
         CLOXEL(I)=1.86E-3/Mass_I(Ion4_)

CALEX cf of H+ and N2  
         !CFHN2(I)=3.36E-9*XN2(I)
         CollisionFreq_IIC(Ion2_,Neutral4_,I)=3.36E-9*XN2(I)
CALEX cf of H+ and O2
         !CFHO2(I)=3.20E-9*XO2(I)
         CollisionFreq_IIC(Ion2_,Neutral3_,I)=3.20E-9*XO2(I)
CALEX cf of H+ and He         
         !CFHHE(I)=1.06E-9*XHE(I)
         CollisionFreq_IIC(Ion2_,Neutral5_,I)=1.06E-9*XHE(I)
CALEX cl for ion neutral collsion       
         CLHO(I)=6.61E-11*XO(I)
         CLHH(I)=2.65E-10*XH(I)
CALEX cl for ion ion collision         
         CLHOX(I)=1.23*17.**1.5/Mass_I(Ion1_)
         CLHHL(I)=1.14*5.**1.5/Mass_I(Ion3_)
         CLHEL(I)=2.97E-2/Mass_I(Ion4_)

CALEX cf of He+ and N2
         !CFHEN2(I)=1.60E-9*XN2(I)
         CollisionFreq_IIC(Ion3_,Neutral4_,I)=1.60E-9*XN2(I)
CALEX cf of He+ and O2         
         !CFHEO2(I)=1.53E-9*XO2(I)
         CollisionFreq_IIC(Ion3_,Neutral3_,I)=1.53E-9*XO2(I)
CALEX cl for He+ and He         
         CLHEHE(I)=8.73E-11*XHE(I)
CALEX cf of He+ and O         
         !CFHEO(I)=1.01E-9*XO(I)
         CollisionFreq_IIC(Ion3_,Neutral1_,I)=1.01E-9*XO(I)
CALEX cf of He+ and H         
         !CFHEH(I)=4.71E-10*XH(I)
         CollisionFreq_IIC(Ion3_,Neutral2_,I)=4.71E-10*XH(I)
CALEX these are cl for ion ion collisions
         CLHEOX(I)=0.57*5.**1.5/Mass_I(Ion1_)
         CLHEHD(I)=0.28*5.**1.5/Mass_I(Ion2_)
         CLHEEL(I)=7.43E-3/Mass_I(Ion4_)
        
CALEX cf of electrons and N2? this doesnt match nagy p 97 though?? 
         !CFELN2(I)=3.42E-9*XN2(I)
         CollisionFreq_IIC(Ion4_,Neutral4_,I)=3.42E-9*XN2(I)
         !CFELO2(I)=3.25E-9*XO2(I)
         CollisionFreq_IIC(Ion4_,Neutral3_,I)=3.25E-9*XO2(I)
         !CFELHE(I)=1.18E-9*XHE(I)
         CollisionFreq_IIC(Ion4_,Neutral5_,I)=1.18E-9*XHE(I)
         !CFELO(I)=2.26E-9*XO(I)
         CollisionFreq_IIC(Ion4_,Neutral1_,I)=2.26E-9*XO(I)
         !CFELH(I)=2.11E-9*XH(I)
         CollisionFreq_IIC(Ion4_,Neutral2_,I)=2.11E-9*XH(I)
CALEX similar to nagy p95 but I cant find the connection???         
         CLELOX(I)=54.5/Mass_I(Ion1_)
         CLELHL(I)=54.5/Mass_I(Ion3_)
         CLELHD(I)=54.5/Mass_I(Ion2_)
         GRAVTY(I)=-3.96E20/RAD(I)**2
         Centrifugal(I)=wHorizontal**2.0*RAD(I)
 20   CONTINUE
cAlex

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
      HeatFlowCoef_II(Ion1_,Neutral4_)=3.*RGAS_I(Ion1_)*Mass_I(Ion1_)/(Mass_I(Ion1_)+28.*XAMU)
      HeatFlowCoef_II(Ion1_,Neutral3_)=3.*RGAS_I(Ion1_)*Mass_I(Ion1_)/(Mass_I(Ion1_)+32.*XAMU)
      HeatFlowCoef_II(Ion1_,Neutral1_)=3.*RGAS_I(Ion1_)*Mass_I(Ion1_)/(Mass_I(Ion1_)+Mass_I(Ion1_))
      HeatFlowCoef_II(Ion1_,Neutral5_)=3.*RGAS_I(Ion1_)*Mass_I(Ion1_)/(Mass_I(Ion1_)+Mass_I(Ion3_))
      HeatFlowCoef_II(Ion1_,Neutral2_)=3.*RGAS_I(Ion1_)*Mass_I(Ion1_)/(Mass_I(Ion1_)+Mass_I(Ion2_))
      HeatFlowCoef_II(Ion1_,Ion1_) = 0.0
 1    HeatFlowCoef_II(Ion1_,Ion2_) = 3.*RGAS_I(Ion1_)*Mass_I(Ion1_)/(Mass_I(Ion1_)+Mass_I(Ion2_))
      HeatFlowCoef_II(Ion1_,Ion3_) = 3.*RGAS_I(Ion1_)*Mass_I(Ion1_)/(Mass_I(Ion1_)+Mass_I(Ion3_))
      HeatFlowCoef_II(Ion1_,Ion4_) = 3.*RGAS_I(Ion1_)*Mass_I(Ion1_)/(Mass_I(Ion1_)+Mass_I(Ion4_))
      
      MassFracCoef_II(Ion1_,:) = HeatFlowCoef_II(Ion1_,:) / (3.0*RGAS_I(Ion1_))


      HeatFlowCoef_II(Ion2_,Neutral4_)=3.*RGAS_I(Ion2_)*Mass_I(Ion2_)/(Mass_I(Ion2_)+28.*XAMU)
      HeatFlowCoef_II(Ion2_,Neutral3_)=3.*RGAS_I(Ion2_)*Mass_I(Ion2_)/(Mass_I(Ion2_)+32.*XAMU)
      HeatFlowCoef_II(Ion2_,Neutral5_)=3.*RGAS_I(Ion2_)*Mass_I(Ion2_)/(Mass_I(Ion2_)+Mass_I(Ion3_))
      HeatFlowCoef_II(Ion2_,Neutral1_)=3.*RGAS_I(Ion2_) *Mass_I(Ion2_)/(Mass_I(Ion2_)+Mass_I(Ion1_))
      HeatFlowCoef_II(Ion2_,Neutral2_)=3.*RGAS_I(Ion2_)*Mass_I(Ion2_)/(Mass_I(Ion2_)+Mass_I(Ion2_))
      HeatFlowCoef_II(Ion2_,Ion2_) = 0.0
      HeatFlowCoef_II(Ion2_,Ion1_) = 3.*RGAS_I(Ion2_) *Mass_I(Ion2_)/(Mass_I(Ion2_)+Mass_I(Ion1_))
      HeatFlowCoef_II(Ion2_,Ion3_) = 3.*RGAS_I(Ion2_)*Mass_I(Ion2_)/(Mass_I(Ion2_)+Mass_I(Ion3_))
      HeatFlowCoef_II(Ion2_,Ion4_) =3.*RGAS_I(Ion2_)*Mass_I(Ion2_)/(Mass_I(Ion2_)+Mass_I(Ion4_))
      
      MassFracCoef_II(Ion2_,:) = HeatFlowCoef_II(Ion2_,:) / (3.0*RGAS_I(Ion2_))


      HeatFlowCoef_II(Ion3_,Neutral4_)=3.*RGAS_I(Ion3_)*Mass_I(Ion3_)/(Mass_I(Ion3_)+28.*XAMU)
      HeatFlowCoef_II(Ion3_,Neutral3_)=3.*RGAS_I(Ion3_)*Mass_I(Ion3_)/(Mass_I(Ion3_)+32.*XAMU)
      HeatFlowCoef_II(Ion3_,Neutral5_)=3.*RGAS_I(Ion3_)*Mass_I(Ion3_)/(Mass_I(Ion3_)+Mass_I(Ion3_))
      HeatFlowCoef_II(Ion3_,Neutral1_)=3.*RGAS_I(Ion3_)*Mass_I(Ion3_)/(Mass_I(Ion3_)+Mass_I(Ion1_))
      HeatFlowCoef_II(Ion3_,Neutral2_)=3.*RGAS_I(Ion3_)*Mass_I(Ion3_)/(Mass_I(Ion3_)+Mass_I(Ion2_))
      HeatFlowCoef_II(Ion3_,Ion3_) = 0.0
      HeatFlowCoef_II(Ion3_,Ion1_)=3.*RGAS_I(Ion3_)*Mass_I(Ion3_)/(Mass_I(Ion3_)+Mass_I(Ion1_))
      HeatFlowCoef_II(Ion3_,Ion2_)=3.*RGAS_I(Ion3_)*Mass_I(Ion3_)/(Mass_I(Ion3_)+Mass_I(Ion2_))
      HeatFlowCoef_II(Ion3_,Ion4_)=3.*RGAS_I(Ion3_)*Mass_I(Ion3_)/(Mass_I(Ion3_)+Mass_I(Ion4_))

      MassFracCoef_II(Ion3_,:) = HeatFlowCoef_II(Ion3_,:) / (3.0*RGAS_I(Ion3_))

      
      HeatFlowCoef_II(Ion4_,Neutral4_)=3.*RGAS_I(Ion4_)*Mass_I(Ion4_)/(28.*XAMU+Mass_I(Ion4_))
      HeatFlowCoef_II(Ion4_,Neutral3_)=3.*RGAS_I(Ion4_)*Mass_I(Ion4_)/(32.*XAMU+Mass_I(Ion4_))
      HeatFlowCoef_II(Ion4_,Neutral5_)=3.*RGAS_I(Ion4_)*Mass_I(Ion4_)/(Mass_I(Ion4_)+Mass_I(Ion3_))
      HeatFlowCoef_II(Ion4_,Neutral1_)=3.*RGAS_I(Ion4_)*Mass_I(Ion4_)/(Mass_I(Ion4_)+Mass_I(Ion1_))
      HeatFlowCoef_II(Ion4_,Neutral2_)=3.*RGAS_I(Ion4_)*Mass_I(Ion4_)/(Mass_I(Ion4_)+Mass_I(Ion2_))
      HeatFlowCoef_II(Ion4_,Ion4_) = 0.0
      HeatFlowCoef_II(Ion4_,Ion1_)=3.*RGAS_I(Ion4_)*Mass_I(Ion4_)/(Mass_I(Ion4_)+Mass_I(Ion1_))
      HeatFlowCoef_II(Ion4_,Ion3_)=3.*RGAS_I(Ion4_)*Mass_I(Ion4_)/(Mass_I(Ion4_)+Mass_I(Ion3_))
      HeatFlowCoef_II(Ion4_,Ion2_)=3.*RGAS_I(Ion4_)*Mass_I(Ion4_)/(Mass_I(Ion4_)+Mass_I(Ion2_))

      MassFracCoef_II(Ion4_,:) = HeatFlowCoef_II(Ion4_,:) / (3.0*RGAS_I(Ion4_))

CALEX CMOXN2 = M_{N2}/(M_o+M_{N2}) see nagy p.83
      FricHeatCoef_II(Ion1_,Neutral4_)=28.*XAMU/(Mass_I(Ion1_)+28.*XAMU)
      FricHeatCoef_II(Ion1_,Neutral3_)=32.*XAMU/(Mass_I(Ion1_)+32.*XAMU)
      FricHeatCoef_II(Ion1_,Neutral1_)=Mass_I(Ion1_)/(Mass_I(Ion1_)+Mass_I(Ion1_))
      FricHeatCoef_II(Ion1_,Neutral5_)=Mass_I(Ion3_)/(Mass_I(Ion1_)+Mass_I(Ion3_))
      FricHeatCoef_II(Ion1_,Neutral2_)=Mass_I(Ion2_)/(Mass_I(Ion1_)+Mass_I(Ion2_))
      FricHeatCoef_II(Ion1_,Ion2_)=Mass_I(Ion2_)/(Mass_I(Ion1_)+Mass_I(Ion2_))
      FricHeatCoef_II(Ion1_,Ion3_)=Mass_I(Ion3_)/(Mass_I(Ion1_)+Mass_I(Ion3_))
      FricHeatCoef_II(Ion1_,Ion4_)=Mass_I(Ion4_)/(Mass_I(Ion4_)+Mass_I(Ion1_))

      FricHeatCoef_II(Ion2_,Neutral4_)=28.*XAMU/(Mass_I(Ion2_)+28.*XAMU)
      FricHeatCoef_II(Ion2_,Neutral3_)=32.*XAMU/(Mass_I(Ion2_)+32.*XAMU)
      FricHeatCoef_II(Ion2_,Neutral5_)=Mass_I(Ion3_)/(Mass_I(Ion2_)+Mass_I(Ion3_))
      FricHeatCoef_II(Ion2_,Neutral1_)=Mass_I(Ion1_)/(Mass_I(Ion2_)+Mass_I(Ion1_))
      FricHeatCoef_II(Ion2_,Neutral2_)=Mass_I(Ion2_)/(Mass_I(Ion2_)+Mass_I(Ion2_))
      FricHeatCoef_II(Ion2_,Ion1_)=Mass_I(Ion1_)/(Mass_I(Ion2_)+Mass_I(Ion1_))
      FricHeatCoef_II(Ion2_,Ion3_)=Mass_I(Ion3_)/(Mass_I(Ion2_)+Mass_I(Ion3_))
      FricHeatCoef_II(Ion2_,Ion4_)=Mass_I(Ion4_)/(Mass_I(Ion4_)+Mass_I(Ion2_))

      FricHeatCoef_II(Ion3_,Neutral4_)=28.*XAMU/(Mass_I(Ion3_)+28.*XAMU)
      FricHeatCoef_II(Ion3_,Neutral3_)=32.*XAMU/(Mass_I(Ion3_)+32.*XAMU)
      FricHeatCoef_II(Ion3_,Neutral5_)=Mass_I(Ion3_)/(Mass_I(Ion3_)+Mass_I(Ion3_))
      FricHeatCoef_II(Ion3_,Neutral1_)=Mass_I(Ion1_)/(Mass_I(Ion3_)+Mass_I(Ion1_))
      FricHeatCoef_II(Ion3_,Neutral2_)=Mass_I(Ion2_)/(Mass_I(Ion3_)+Mass_I(Ion2_))
      FricHeatCoef_II(Ion3_,Ion1_)    =Mass_I(Ion1_)/(Mass_I(Ion3_)+Mass_I(Ion1_))
      FricHeatCoef_II(Ion3_,Ion2_)    =Mass_I(Ion2_)/(Mass_I(Ion3_)+Mass_I(Ion2_))
      FricHeatCoef_II(Ion3_,Ion4_)    =Mass_I(Ion4_)/(Mass_I(Ion4_)+Mass_I(Ion3_))

      FricHeatCoef_II(Ion4_,Neutral4_)=28.*XAMU/(Mass_I(Ion4_)+28.*XAMU)
      FricHeatCoef_II(Ion4_,Neutral3_)=32.*XAMU/(Mass_I(Ion4_)+32.*XAMU)
      FricHeatCoef_II(Ion4_,Neutral5_)=Mass_I(Ion3_)/(Mass_I(Ion4_)+Mass_I(Ion3_))
      FricHeatCoef_II(Ion4_,Neutral1_)=Mass_I(Ion1_)/(Mass_I(Ion4_)+Mass_I(Ion1_))
      FricHeatCoef_II(Ion4_,Neutral2_)=Mass_I(Ion2_)/(Mass_I(Ion4_)+Mass_I(Ion2_))
      FricHeatCoef_II(Ion4_,Ion1_)    =Mass_I(Ion1_)/(Mass_I(Ion4_)+Mass_I(Ion1_))
      FricHeatCoef_II(Ion4_,Ion2_)    =Mass_I(Ion2_)/(Mass_I(Ion4_)+Mass_I(Ion2_))
      FricHeatCoef_II(Ion4_,Ion3_)    =Mass_I(Ion3_)/(Mass_I(Ion4_)+Mass_I(Ion3_))



C     ALEX(10/11/04): 
C     TRY SETTING THE PLASMA PARAMETERS HERE TO THE SURFACE VALUES

      if(IsRestart) RETURN

      do K=1,NDIM
c         State_GV(K,RhoH_)=State_GV(0,RhoH)
c         State_GV(K,uO_)=0
c         State_GV(K,pO_)=State_GV(0,Po_)H
c         State_GV(K,RhoO_)=State_GV(0,RhoO_)H
c         State_GV(K,To_)=State_GV(0,To_)
c         State_GV(K,uH_)=0
c         State_GV(K,pH_)=State_GV(0,pH_)H
c         State_GV(K,Th_)=State_GV(0,Th_)
c         State_GV(K,uHe_)=0
c         State_GV(K,pHe_)=State_GV(0,pHe_)H
c         State_GV(K,RhoHe_)=State_GV(0,RhoHe_)H
c         State_GV(K,The_)=State_GV(0,The_)
c         State_GV(K,RhoE_)=State_GV(0,RhoE_)H
c         State_GV(K,uE_)=0
c         State_GV(K,pE_)=State_GV(0,pE_)H
c         State_GV(K,Te_)=State_GV(0,Te_)

         IsRestart = .true.

         State_GV(K,RhoH_)=1.*State_GV(k-1,RhoH_)*exp(-DrBnd/10000.0e5)
         State_GV(K,uO_)=1.
         State_GV(K,pO_)=1.*State_GV(k-1,Po_)*exp(-1.*DrBnd/10000.0e5)
         State_GV(K,RhoO_)=1.*State_GV(k-1,RhoO_)*exp(-1.*DrBnd/10000.0e5)
c         State_GV(K,To_)=State_GV(k-1,To_)
         State_GV(K,To_)=State_GV(K,pO_)/(1.38e-16*State_GV(K,RhoO_)/Mass_I(Ion1_))
         State_GV(K,uH_)=10.
         State_GV(K,pH_)=1.*State_GV(k-1,pH_)*exp(-1.*DrBnd/10000.0e5)
c         State_GV(K,Th_)=State_GV(k-1,Th_)*exp(-1.*DrBnd/9000.)
         State_GV(K,Th_)=State_GV(K,pH_)/(1.38e-16*State_GV(K,RhoH_)/Mass_I(Ion2_))
         State_GV(K,uHe_)=4.
         State_GV(K,pHe_)=1.*State_GV(k-1,pHe_)*exp(-1.*DrBnd/10000.0e5)
         State_GV(K,RhoHe_)=1.*State_GV(k-1,RhoHe_)*exp(-1.*DrBnd/10000.0e5)
c         State_GV(K,The_)=State_GV(k-1,The_)
         State_GV(K,The_)=State_GV(K,pHe_)/(1.38e-16*State_GV(K,RhoHe_)/Mass_I(Ion3_))
         State_GV(K,RhoE_)=1.*State_GV(k-1,RhoE_)*exp(-1.*DrBnd/10000.0e5)
         State_GV(K,uE_)=1.
         State_GV(K,pE_)=1.*State_GV(k-1,pE_)*exp(-1.*DrBnd/10000.0e5)
c         State_GV(K,Te_)=State_GV(k-1,Te_)
         if (k < nDim-20)then
            State_GV(K,Te_)=State_GV(K,pE_)/(1.38e-16*State_GV(K,RhoE_)/Mass_I(Ion4_))
         else
            State_GV(K,Te_)=(1030.-780.)/(20.)**2.*(real(k-(nDim-20)))**2.+780.
            State_GV(K,pE_)=State_GV(k,Te_)*Rgas_I(nIon)*State_GV(K,RhoE_)
         endif
         
      enddo

      RETURN
      END
