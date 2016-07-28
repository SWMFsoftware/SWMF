 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     ALEX(10/11/04): I THINK THAT THIS SUBROUTINE INITIALIZES THE GRID AND
C     AND SETS CONSTANTS. IT MUST BE INITIALIZED EVERY TIME THE CODE RUNS
C     EVEN IF RESTARTING FROM A FILE.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE STRT
C
C
      use ModCommonVariables
      use ModCommonPlanet,ONLY: HLPion1,HLPion2,HLPE,HLPE0
C
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

CALEX not I use O for H3 and he for H2
C Gas constant = k_Boltzmann/AMU
      RGAS=8.314E7
C Adiabatic index
      GAMMA=5./3.
C AMU in gramms
      XAMU=1.6606655E-24
C Mass of atomic H3 in gramms
      Mass_I(Ion1_)=3.0237*XAMU
C Mass of atomic H in gramms
      Mass_I(Ion2_)=1.00797*XAMU
C Mass of H2 in gramms
      Mass_I(Ion3_)=2.0159*XAMU
C Mass of electron in gramms
      Mass_I(nIon)=9.109534E-28
C Relative mass of H3 to electron
      MassElecIon_I(Ion1_)=Mass_I(nIon)/Mass_I(Ion1_)
C Relative mass of atomic H to electron
      MassElecIon_I(Ion2_)=Mass_I(nIon)/Mass_I(Ion2_)
C Relative mass of H2 to electron
      MassElecIon_I(Ion3_)=Mass_I(nIon)/Mass_I(Ion3_)
C kB/m_H3
      RGAS_I(Ion1_)=RGAS*XAMU/Mass_I(Ion1_)
C kB/m_H
      RGAS_I(Ion2_)=RGAS*XAMU/Mass_I(Ion2_)
C kB/m_H2
      RGAS_I(Ion3_)=RGAS*XAMU/Mass_I(Ion3_)
C kB/m_e
      RGAS_I(nIon)=RGAS*XAMU/Mass_I(nIon)
      GMIN1=GAMMA-1.
      GMIN2=GMIN1/2.
      GPL1=GAMMA+1.
      GPL2=GPL1/2.
      GM12=GMIN1/GAMMA/2.
      GRAR=GAMMA/GMIN2
      GREC=1./GAMMA
      CPion1=GAMMA*RGAS_I(Ion1_)/GMIN1
      CPion2=GAMMA*RGAS_I(Ion2_)/GMIN1
      CPion3=GAMMA*RGAS_I(Ion3_)/GMIN1
      CPE=GAMMA*RGAS_I(nIon)/GMIN1
      CVion1=RGAS_I(Ion1_)/GMIN1
      CVion2=RGAS_I(Ion2_)/GMIN1
      CVion3=RGAS_I(Ion3_)/GMIN1
      CVE=RGAS_I(nIon)/GMIN1

CALEX Set the planet radius and surface gravity, rotation freq
      RE=60268.E5
      GSURF=980.*.916
c      Omega=1./37800.

      Omega=0.
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     READ FIELD ALIGNED CURRENT DENSITY AT LOWEST GRID POINT          C
C        (UNIT=AMPERE/M**2)                                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C      READ (5,2) CURR(1)

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
C     IYD=84172
C     F107A=60.
C     F107=60.
C      SEC=43200.
C      STL=12.
C      GMLAT=80.
C      GMLONG=0.
C      IART=1
      GLAT=80.
C      GLONG=0.
C      IART=0
C FEB 20, 1990
CALEX IYD=year_day of year

c      IYD=90051
c      SEC=20.75*3600.
c      F107A=180.
c      F107=189.5
c      IART=0
c      GLONG=325.13
c      GLAT=70.47
C END
CALEX GGM determines geomagnetic lat. from geographic lat and lon
CALEX since the dipole on Saturn is aligned with the rotation axis
CALEX I have set GMLONG=GLONG
CALEX      CALL GGM(IART,GLONG,GLAT,GMLONG,GMLAT)
!      GMLONG=GLONG
!      GMLAT=GLAT

      DO 49 I=1,7
      AP(I)=50.
49    CONTINUE 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C                                                                      C
      DO 50 K=1,NDIM
         CALL MODATM(ALTD(K),XH2(K),XH(K),XH2O(K),XCH4(K),XTN(K))
         NDensity_CI(K,H2_) = XH2(K)
         NDensity_CI(K,H_)  = XH(K)
         NDensity_CI(K,H2O_)= XH2O(K)
         NDensity_CI(K,CH4_)= XCH4(K)

50    CONTINUE


CALEX I am not calling glowex now but in the future 
CALEX we might need to use this for radiative transfer etc      
CALEX      CALL GLOWEX
CALEX      DO 1099 J = 1,40
CALEX         WRITE (iUnitOutput,9999) ALTD(J),PHOTOTF(J+1) 
CALEX 9999    FORMAT(2X,1PE15.3,2X,1PE15.3)
CALEX 1099 CONTINUE
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

CALEX I don't know how to fix this heat input for Saturn.
CALEX In Dee's thesis she says you need electon heat flux to be a minimum
CALEX of 20E-3 ergs cm^2 /s.      
C      ETOP=1.0E-3
      ETOP=20.0E-3
c      ETOP=25.0E-3
      ELFXIN=0.
C
C      ELFXIN=9.
C
2     FORMAT(6X,1PE15.4)
     
C      READ (5,2) HEATI1,HEATI2,ELHEAT
      HEATI1=0.
      HEATI2=0.
C
      HEATI1=1.0E-7
calex origionally heatI1 was 0 and heat I2 was not 0, but i have a suspicion that
Calex this should not be true. I think the two refers to helium which we are not 
calex looking at.
c      HEATI2=2.5E-11
C
      ELHEAT=0.
      
      HEATA1=3.5E7
      HEATA2=2.0E8
      HEATA3=1.5E8
      HEATS1=2.*1.0E7**2
      HEATS2=2.*2.5E7**2
      HEATS3=2.*1.0E7**2
      DO 53 K=1,NDIM
      HEATX1=EXP(-(ALTD(K)-HEATA1)**2/HEATS1)
      HEATX2=EXP(-(ALTD(K)-HEATA2)**2/HEATS2)
      HEATX3=EXP(-(ALTD(K)-HEATA3)**2/HEATS3)
C KGS What's all this? delete?
      QOXYG(K)=(HEATI1*HEATX1+HEATI2*HEATX2)/
     #         (State_GV(K,RhoH3_)+State_GV(K,RhoH_))
calex origionally it was qhyd=qoxy/16 but I think 16 should be 3 since
calex I am letting oxy stand in for H3+
      QHYD(K)=QOXYG(K)/3.
      QHEL(K)=QOXYG(K)/4.
C
      QOXYG(K)=0.
C
      QOXYG(K)=0.
      QHYD(K)=0.
      QHEL(K)=0.

      QELECT(K)=ELHEAT*HEATX3
53    CONTINUE
      DO 54 K=1,NDIM,10
      
52    FORMAT(5(1PE15.4))
54    CONTINUE
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
C      READ (5,2) DT   1/20.0 is a good value
C     Read in the time step
      write(*,*) dt
c      DT=1./1.
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
      HLPE0=GMIN1/RGAS_I(nIon)
      HLPE=1.23E-6*GMIN1/RGAS_I(nIon)
      HLPion1=2.86E-8*(Mass_I(nIon)/Mass_I(Ion1_))*GMIN1/RGAS_I(Ion1_)
      HLPion2=7.37E-8*(Mass_I(nIon)/Mass_I(Ion2_))*GMIN1/RGAS_I(Ion2_)
! KGS fix this
      HLPion3=7.37E-8*(Mass_I(nIon)/Mass_I(Ion3_))*GMIN1/RGAS_I(Ion3_)

 
      
CALEX These are the heat conductivities at the lower boundary. Note:
CALEX that no allowance is made to take into account the effect of
CALEX neutrals on the heat conduction as was done at earth.     
      HeatCon_GI(0,Ion1_)=HLPion1*(State_GV(0,RhoH3_)/State_GV(0,RhoE_))*State_GV(0,Th3_)**2.5
      HeatCon_GI(0,nIon)=HLPE*State_GV(0,Te_)**2.5
      HeatCon_GI(0,Ion2_)=HLPion2*(State_GV(0,RhoH_)/State_GV(0,RhoE_))*State_GV(0,Th_)**2.5
      HeatCon_GI(0,Ion3_)=HLPion3*(State_GV(0,RhoH2_)/State_GV(0,RhoE_))*State_GV(0,Th2_)**2.5

      
C!      HeatCon_GI(0,Ion1_)=HLPion1*State_GV(0,Th3_)**2.5
C!      HeatCon_GI(0,nIon)=HLPE*State_GV(0,Te_)**2.5
C!      HeatCon_GI(0,Ion2_)=HLPion2*State_GV(0,Th_)**2.5
C!      HeatCon_GI(0,Ion2_)E=HLPHE*State_GV(0,The_)**2.5
      CALL MODATM(ALTMAX,XNH2,XNH,XNH2O,XNCH4,TEMP)
      XTNMAX=TEMP

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

      CALL MODATM(ALTMIN,XNH2,XNH,XNH2O,XNCH4,XNT)
      CALL MODATM(ALTMAX,YNH2,YNH,YNH2O,YNCH4,YNT)
      DO 60 I=1,NDIM
      ALTD(I)=ALTD(I)/1.E5
60    CONTINUE
      ALTMIN=ALTMIN/1.E5
      ALTMAX=ALTMAX/1.E5
      ETOP1=ETOP*1.23E-6/DRBND
      CALL COLLIS(NDIM,State_GV(-1:nDim+2,:))
      CALL PW_calc_efield(nDim,State_GV(-1:nDim+2,:))

      if (DoLog) then

      IF (NCNPRT.NE.0) GO TO 999

C KGS What's all this?      
      WRITE(iUnitOutput,1005) NDIM
1005  FORMAT(1H1,5X,'NUMBER OF CELLS=',I4)
      WRITE(iUnitOutput,1020) NEXP
1020  FORMAT(5X,'NEXP=',I1)
      WRITE (iUnitOutput,1008) GAMMA,RGAS_I(Ion1_),CPion1,CVion1,
     ;RGAS_I(Ion2_),CPion2,CVion2,RGAS_I(nIon),CPE,CVE
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
      WRITE (iUnitOutput,1009) State_GV(0,uH3_),State_GV(0,pH3_),State_GV(0,RhoH3_),State_GV(0,Th3_),SoundSpeed_GI(0,Ion1_)
1009  FORMAT(5X,'VELOCITY=',1PE11.4,3X,'PRESSURE=',1PE10.4,3X,
     ;'MASS DENSITY=',1PE10.4,3X,'TEMPERATURE=',1PE10.4,3X,
     ;'SOUND VELOCITY=',1PE10.4)
      WRITE(iUnitOutput,1002)
1002  FORMAT(1H ,4X,'HYDROGEN:')
      WRITE (iUnitOutput,1009) State_GV(0,uH_),State_GV(0,pH_),State_GV(0,RhoH_),State_GV(0,Th_),SoundSpeed_GI(0,Ion2_)
      WRITE(iUnitOutput,1003)
1003  FORMAT(1H ,4X,'ELECTRONS:')
      WRITE (iUnitOutput,1009) State_GV(0,uE_),State_GV(0,pE_),State_GV(0,RhoE_),State_GV(0,Te_),SoundSpeed_GI(0,nIon)
      WRITE (iUnitOutput,1027)
1027  FORMAT(1H0,5X,'UPPER BOUNDARY INITIAL PLASMA PARAMETERS:')
      WRITE (iUnitOutput,1004)
1004  FORMAT(1H ,4X,'OXYGEN:')
      WRITE (iUnitOutput,1009) 
     & State_GV(nDim+1,uH3_),State_GV(nDim+1,pH3_),State_GV(nDim+1,RhoH3_),State_GV(nDim+1,Th3_),SoundSpeed_GI(nDim+1,Ion1_)

      WRITE (iUnitOutput,1006)
1006  FORMAT(1H ,4X,'HYDROGEN:')
      WRITE (iUnitOutput,1009) 
     & State_GV(nDim+1,uH_),State_GV(nDim+1,pH_),State_GV(nDim+1,RhoH_),State_GV(nDim+1,Th_),SoundSpeed_GI(nDim+1,Ion2_)
      WRITE(iUnitOutput,1007)
1007  FORMAT(1H ,4X,'ELECTRONS:')
      WRITE (iUnitOutput,1009) 
     & State_GV(nDim+1,uE_),State_GV(nDim+1,pE_),State_GV(nDim+1,RhoE_),State_GV(nDim+1,Te_),SoundSpeed_GI(nDim+1,nIon)
      WRITE (iUnitOutput,1029) ETOP1
1029  FORMAT(1H0,5X,'TOPSIDE ELECTRON HEATING RATE:',1PE10.4,
     ;' ERGS/CM**3/SEC')

      write(iUnitOutput,*) 'CMH3pH,CMH3pH2,CMH3pHp,CMH3pEL'
      write(iUnitOutput,*) FricHeatCoef_II(Ion1_,Neutral2_),
     &     FricHeatCoef_II(Ion1_,Neutral1_),FricHeatCoef_II(Ion1_,Ion2_),
     &     FricHeatCoef_II(Ion1_,nIon)
     
      write(iUnitOutput,*) 'CMHpH,CMHpH2,CMHpH3p,CMHpEL'
      write(iUnitOutput,*) FricHeatCoef_II(Ion2_,Neutral2_),
     &     FricHeatCoef_II(Ion2_,Neutral1_),FricHeatCoef_II(Ion2_,Ion1_),
     &     FricHeatCoef_II(Ion2_,nIon)
      
      write(iUnitOutput,*) 'CMELH,CMELH2,CMELHp,CMELH3p'
      write(iUnitOutput,*) CMELH,CMELH2,CMELHp,CMELH3p
      write(iUnitOutput,*) FricHeatCoef_II(nIon,Neutral2_),
     &     FricHeatCoef_II(nIon,Neutral1_),FricHeatCoef_II(nIon,Ion2_),
     &     FricHeatCoef_II(nIon,Ion1_)

!      WRITE (iUnitOutput,1050) 
!1050  FORMAT(1H0,5X,'ENERGY COLLISION TERM COEFFICIENTS')
!      WRITE (iUnitOutput,1051)CTHpH,CTHpH2,CTHpH3p,CTHpEL,CTH3pH,CTH3pH2,CTH3pHp,
!     $ CTH3pEL,CTELH,CTELH2,CTELHp,CTELH3p
!
!1051  FORMAT(1H0,5X,'CTHpH=',1PE10.4,5X,'CTHpH2=',1PE10.4,4X,
!     $'CTHpH3p=',1PE10.4/5X,'CTHpEL=',1PE10.4,5X,'CTH3pH=',1PE10.4,5X,
!     $'CTH3pH2=',1PE10.4,5X,'CTH3pHp=',1PE10.4,5X,'CTH3pEL=',1PE10.4/
!     $5X,'CTELH=',1PE10.4,5X,'CTELH2=',1PE10.4,4X,
!     $'CTELHp=',1PE10.4/5X,'CTELH3p=',1PE10.4,5X)
!
!      WRITE (iUnitOutput,1052) CMHpH,CMHpH2,CMHpH3p,CMHpEL,CMH3pH,CMH3pH2,CMH3pHp,
!     $ CMH3pEL,CMELH,CMELH2,CMELHp,CMELH3p
!
!1052  FORMAT(1H0,5X,'CMHpH=',1PE10.4,5X,'CMHpH2=',1PE10.4,4X,
!     $'CMHpH3p=',1PE10.4/5X,'CMHpEL=',1PE10.4,5X,'CMH3pH=',1PE10.4,5X,
!     $'CMH3pH2=',1PE10.4,5X,'CMH3pHp=',1PE10.4,5X,'CMH3pEL=',1PE10.4/
!     $5X,'CMELH=',1PE10.4,5X,'CMELH2=',1PE10.4,4X,
!     $'CMELHp=',1PE10.4/5X,'CMELH3p=',1PE10.4,5X)
!
!      WRITE (iUnitOutput,1053)
!
!1053  FORMAT(1H0,5X,'HEAT CONDUCTION COEFFICIENTS AT UPPER BOUNDARY')
!
!C!      WRITE (iUnitOutput,1054) CZHN2,CZHO2,CZHO,CZHOX,CZHEN2,CZHEO2,CZHEHE,
!C!     $CZHEO,CZHEH,CZHEOX,CZHEHD,XTNMAX
!C!1054  FORMAT(1H0,5X,'CZHN2=',1PE10.4,6X,'CZHO2=',1PE10.4,6X,
!C!     $'CZHO=',1PE10.4,7X,'CZHOX=',1PE10.4/5X,'CZHEN2=',1PE10.4,5X,
!C!     $'CZHEO2=',1PE10.4,5X,'CZHEHE=',1PE10.4,5X,'CZHEO=',1PE10.4/
!C!     $5X,'CZHEH=',1PE10.4,6X,'CZHEOX=',1PE10.4,5X,'CZHEHD=',
!C!     $1PE10.4/5X,'XTNMAX=',1PE10.4)
!
      WRITE (iUnitOutput,1012)

1012  FORMAT(1H1,45X,'NEUTRAL ATMOSPHERE NUMBER DENSITIES')

      WRITE(iUnitOutput,1014)

1014  FORMAT(16X,'ALT',13X,'H2',13X,'H',15X,'H2O',14X,'CH4',
     $     15X,'T')

      K=0

      WRITE (iUnitOutput,1022) K, ALTMIN,XNH2,XNH,XNH2O,XNCH4,XNT

      NDMQ=NPT1

      IF (NDIM.LT.NPT2) NDMQ=NDIM

      do K=1,NDMQ
         WRITE(iUnitOutput,1022) K,ALTD(K),XH2(K),XH(K),XH2O(K),XCH4(K),
     $        XTN(K)
      enddo


      IF (NDIM.LT.NPT2) GO TO 290

      NDMQ=NPT3

      IF (NDIM.LT.NPT4) NDMQ=NDIM
      
      do K=NPT2,NDMQ,2
         WRITE(iUnitOutput,1022) K,ALTD(K),XH2(K),XH(K),XH2O(K),XCH4(K),
     $        XTN(K)
      enddo

      IF (NDIM.LT.NPT4) GO TO 290
      NDMQ=NPT5
      
      IF (NDIM.LT.NPT6) NDMQ=NDIM
      
      do K=NPT4,NDMQ,5
         WRITE(iUnitOutput,1022) K,ALTD(K),XH2(K),XH(K),XH2O(K),XCH4(K),
     $        XTN(K)
      enddo
      IF (NDIM.LT.NPT6) GO TO 290
      do K=NPT6,NDIM,10
         WRITE(iUnitOutput,1022) K,ALTD(K),XH2(K),XH(K),XH2O(K),XCH4(K),
     $        XTN(K)
      enddo
 290  CONTINUE
      
      K=NDIM+1
      
      WRITE (iUnitOutput,1022) K, ALTMAX,YNH2,YNH,YNH2O,YNCH4,YNT
      
      WRITE(iUnitOutput,1015)
c
calex Now output the source coef.

1015  FORMAT(1H1,55X,'SOURCE COEFFICIENTS:')
      WRITE (iUnitOutput,1016)


1016  FORMAT(12X,'ALT',6X,'GRAVTY',6X,'Centrifugal',5X,'FFHpp1',5X,
     ;'FFHpp3',5X,
     ;'FFHpp4',5X,'FFHpc2',5X,'FFHpc3',5X,'FFHpc8',4X,'FFHpr1',
     ;4X,'FFH3pc1',5X,'FFH3pc2',5X,'FFH3pc6',5X,'FFH3pc7',
     ;5X,'FFH3pr2')

      NDMQ=NPT1

      IF (NDIM.LT.NPT2) NDMQ=NDIM

      DO  K=1,NDMQ
         WRITE(iUnitOutput,1019) K,ALTD(K),GRAVTY(K),Centrifugal(K),FFHpp1(K),FFHpp3(K),
     $        FFHpp4(K),FFHpc2(K),FFHpc3(K),FFHpc8(K),FFHpr1(K),FFH3pc1(K),
     $        FFH3pc2(K), FFH3pc6(K),FFH3pc7(K),FFH3pr2(K)
      enddo


      IF (NDIM.LT.NPT2) GO TO 295
      NDMQ=NPT3
      IF (NDIM.LT.NPT4) NDMQ=NDIM

      do K=NPT2,NDMQ,2
      WRITE(iUnitOutput,1019) K,ALTD(K),GRAVTY(K),FFHpp1(K),FFHpp3(K),
     $        FFHpp4(K),FFHpc2(K),FFHpc3(K),FFHpc8(K),FFHpr1(K),FFH3pc1(K),
     $        FFH3pc2(K), FFH3pc6(K),FFH3pc7(K),FFH3pr2(K)

      enddo

      IF (NDIM.LT.NPT4) GO TO 295
      NDMQ=NPT5

      IF (NDIM.LT.NPT6) NDMQ=NDIM

      do K=NPT4,NDMQ,5
         WRITE(iUnitOutput,1019) K,ALTD(K),GRAVTY(K),FFHpp1(K),FFHpp3(K),
     $        FFHpp4(K),FFHpc2(K),FFHpc3(K),FFHpc8(K),FFHpr1(K),FFH3pc1(K),
     $        FFH3pc2(K), FFH3pc6(K),FFH3pc7(K),FFH3pr2(K)

      enddo

      IF (NDIM.LT.NPT6) GO TO 295

      do K=NPT6,NDIM,10
         WRITE(iUnitOutput,1019) K,ALTD(K),GRAVTY(K),FFHpp1(K),FFHpp3(K),
     $        FFHpp4(K),FFHpc2(K),FFHpc3(K),FFHpc8(K),FFHpr1(K),FFH3pc1(K),
     $        FFH3pc2(K), FFH3pc6(K),FFH3pc7(K),FFH3pr2(K)

      enddo
295   CONTINUE

      WRITE(iUnitOutput,1017)
1017  FORMAT(1H1,50X,'COLLISION COEFFICIENTS FOR H3p')

C!      WRITE(iUnitOutput,1018)
C!1018  FORMAT(12X,'ALT',7X,'CLOXN2',6X,'CLOXO2',6X,'CLOXO',7X,
C!     ;'CLOXHE',6X,'CLOXH',7X,'CLOXHL',6X,'CLOXHD',6X,'CLOXEL')
C!      NDMQ=NPT1
C!      IF (NDIM.LT.NPT2) NDMQ=NDIM
C!      DO 530 K=1,NDMQ
C!      WRITE (iUnitOutput,1180) K,ALTD(K),CLOXN2(K),CLOXO2(K),CLOXO(K),
C!     $CLOXHE(K),CLOXH(K),CLOXHL(K),CLOXHD(K),CLOXEL(K)
C!530   CONTINUE
C!      IF (NDIM.LT.NPT2) GO TO 590
C!      NDMQ=NPT3
C!      IF (NDIM.LT.NPT4) NDMQ=NDIM
C!      DO 540 K=NPT2,NDMQ,2
C!      WRITE (iUnitOutput,1180) K,ALTD(K),CLOXN2(K),CLOXO2(K),CLOXO(K),
C!     $CLOXHE(K),CLOXH(K),CLOXHL(K),CLOXHD(K),CLOXEL(K)
C!540   CONTINUE
C!      IF (NDIM.LT.NPT4) GO TO 590
C!      NDMQ=NPT5
C!      IF (NDIM.LT.NPT6) NDMQ=NDIM
C!      DO 550 K=NPT4,NDMQ,5
C!      WRITE (iUnitOutput,1180) K,ALTD(K),CLOXN2(K),CLOXO2(K),CLOXO(K),
C!     $CLOXHE(K),CLOXH(K),CLOXHL(K),CLOXHD(K),CLOXEL(K)
C!550   CONTINUE
C!      IF (NDIM.LT.NPT6) GO TO 590
C!      DO 560 K=NPT6,NDIM,10
C!      WRITE (iUnitOutput,1180) K,ALTD(K),CLOXN2(K),CLOXO2(K),CLOXO(K),
C!     $CLOXHE(K),CLOXH(K),CLOXHL(K),CLOXHD(K),CLOXEL(K)
C!560   CONTINUE
C!590   CONTINUE

C!1      WRITE (iUnitOutput,1217)
C!11217  FORMAT(1H1,50X,'COLLISION COEFFICIENTS FOR HELIUM')
C!1      WRITE(iUnitOutput,1218)
C!11218  FORMAT(12X,'ALT',7X,'CLHEN2',6X,'CLHEO2',6X,'CLHEO',7X,
C!1     ;'CLHEHE',6X,'CLHEH',7X,'CLHEOX',6X,'CLHEHD',6X,'CLHEEL')
C!1      NDMQ=NPT1
C!1      IF (NDIM.LT.NPT2) NDMQ=NDIM
C!1      DO 1230 K=1,NDMQ
C!1      WRITE (iUnitOutput,1180) K,ALTD(K),CLHEN2(K),CLHEO2(K),CLHEO(K),
C!1     $CLHEHE(K),CLHEH(K),CLHEOX(K),CLHEHD(K),CLHEEL(K)
C!11230  CONTINUE
C!1      IF (NDIM.LT.NPT2) GO TO 1290
C!1      NDMQ=NPT3
C!1      IF (NDIM.LT.NPT4) NDMQ=NDIM
C!1      DO 1240 K=NPT2,NDMQ,2
C!1      WRITE (iUnitOutput,1180) K,ALTD(K),CLHEN2(K),CLHEO2(K),CLHEO(K),
C!1     $CLHEHE(K),CLHEH(K),CLHEOX(K),CLHEHD(K),CLHEEL(K)
C!11240  CONTINUE
C!1      IF (NDIM.LT.NPT4) GO TO 1290
C!1      NDMQ=NPT5
C!1      IF (NDIM.LT.NPT6) NDMQ=NDIM
C!1      DO 1250 K=NPT4,NDMQ,5
C!1      WRITE (iUnitOutput,1180) K,ALTD(K),CLHEN2(K),CLHEO2(K),CLHEO(K),
C!1     $CLHEHE(K),CLHEH(K),CLHEOX(K),CLHEHD(K),CLHEEL(K)
C!11250  CONTINUE
C!1      IF (NDIM.LT.NPT6) GO TO 1290
C!1      DO 1260 K=NPT6,NDIM,10
C!1      WRITE (iUnitOutput,1180) K,ALTD(K),CLHEN2(K),CLHEO2(K),CLHEO(K),
C!1     $CLHEHE(K),CLHEH(K),CLHEOX(K),CLHEHD(K),CLHEEL(K)
C!11260  CONTINUE
C!11290  CONTINUE
C!1      WRITE (iUnitOutput,2217)

2217  FORMAT(1H1,50X,'COLLISION COEFFICIENTS FOR HYDROGEN')
      WRITE(iUnitOutput,2218)
2218  FORMAT(12X,'ALT',7X,'CLHpH3p',7X,'CLHpH')
      NDMQ=NPT1
      IF (NDIM.LT.NPT2) NDMQ=NDIM
      do K=1,NDMQ
         WRITE (iUnitOutput,1180) K,ALTD(K),CLHpH3p(K),CLHpH(K)
      enddo
      
      IF (NDIM.LT.NPT2) GO TO 2290
      NDMQ=NPT3
      
      IF (NDIM.LT.NPT4) NDMQ=NDIM
      
      do K=NPT2,NDMQ,2
         WRITE (iUnitOutput,1180) K,ALTD(K),CLHpH3p(K),CLHpH(K)
      enddo
      
      IF (NDIM.LT.NPT4) GO TO 2290
      NDMQ=NPT5
      
      IF (NDIM.LT.NPT6) NDMQ=NDIM
      
      do  K=NPT4,NDMQ,5
         WRITE (iUnitOutput,1180) K,ALTD(K),CLHpH3p(K),CLHpH(K)
      enddo
      
      IF (NDIM.LT.NPT6) GO TO 2290
      
      do K=NPT6,NDIM,10
         WRITE (iUnitOutput,1180) K,ALTD(K),CLHpH3p(K),CLHpH(K)
      enddo
2290  CONTINUE

      WRITE (iUnitOutput,3217)
3217  FORMAT(1H1,50X,'COLLISION COEFFICIENTS FOR ELECTRONS')
      WRITE(iUnitOutput,3218)
3218  FORMAT(12X,'ALT',6X,'CLELHp',6X,'CLELH3p',6X,'CLELH')
      
      NDMQ=NPT1
      
      IF (NDIM.LT.NPT2) NDMQ=NDIM
      
      do K=1,NDMQ
         WRITE (iUnitOutput,1180) K,ALTD(K),CLELHp(K),CLELH3p(K),CLELH(K)
      enddo
      
      IF (NDIM.LT.NPT2) GO TO 3290
      
      NDMQ=NPT3
      
      IF (NDIM.LT.NPT4) NDMQ=NDIM
      
      do K=NPT2,NDMQ,2
         WRITE (iUnitOutput,1180) K,ALTD(K),CLELHp(K),CLELH3p(K),CLELH(K)
      enddo
      
      IF (NDIM.LT.NPT4) GO TO 3290
      
      NDMQ=NPT5
      
      IF (NDIM.LT.NPT6) NDMQ=NDIM
      
      do K=NPT4,NDMQ,5
         WRITE (iUnitOutput,1180) K,ALTD(K),CLELHp(K),CLELH3p(K),CLELH(K)
      enddo

      IF (NDIM.LT.NPT6) GO TO 3290

      do K=NPT6,NDIM,10
         WRITE (iUnitOutput,1180) K,ALTD(K),CLELHp(K),CLELH3p(K),CLELH(K)
      enddo
3290  CONTINUE


      WRITE(iUnitOutput,1047)
1047  FORMAT(1H1,50X,'COLLISION FREQUENCIES FOR H3')
C!      WRITE(iUnitOutput,1048)
C!1048  FORMAT(12X,'ALT',7X,'CFH3pHp',6X,'CFH3pEL',6X,'CFH3pH',7X,
C!     ;'CFH3pH2')
C!
C!      NDMQ=NPT1
C!
C!      IF (NDIM.LT.NPT2) NDMQ=NDIM
C!
C!      do K=1,NDMQ
C!         WRITE (iUnitOutput,1180) K,ALTD(K),CFH3pHp(K),CFH3pEL(K),CFH3pH(K),
C!     $        CFH3pH2(K)
C!      enddo
C!      
C!      IF (NDIM.LT.NPT2) GO TO 591
C!      
C!      NDMQ=NPT3
C!      
C!      IF (NDIM.LT.NPT4) NDMQ=NDIM
C!      
C!      do K=NPT2,NDMQ,2
C!         WRITE (iUnitOutput,1180) K,ALTD(K),CFH3pHp(K),CFH3pEL(K),CFH3pH(K),
C!     $        CFH3pH2(K)
C!      enddo
C!      
C!      IF (NDIM.LT.NPT4) GO TO 591
C!      
C!      NDMQ=NPT5
C!      
C!      IF (NDIM.LT.NPT6) NDMQ=NDIM
C!      
C!      do K=NPT4,NDMQ,5
C!         WRITE (iUnitOutput,1180) K,ALTD(K),CFH3pHp(K),CFH3pEL(K),CFH3pH(K),
C!     $        CFH3pH2(K)
C!      enddo
C!      
C!      IF (NDIM.LT.NPT6) GO TO 591
C!      
C!      do K=NPT6,NDIM,10
C!         WRITE (iUnitOutput,1180) K,ALTD(K),CFH3pHp(K),CFH3pEL(K),CFH3pH(K),
C!     $        CFH3pH2(K)
C!      enddo
C!591   CONTINUE

C!1      WRITE (iUnitOutput,1247)
C!11247  FORMAT(1H1,50X,'COLLISION FREQUENCIES FOR HELIUM')
C!1      WRITE(iUnitOutput,1248)
C!11248  FORMAT(12X,'ALT',7X,'CFHEN2',6X,'CFHEO2',6X,'CFHEO',7X,
C!1     ;'CFHEHE',6X,'CFHEH',7X,'CFHEOX',6X,'CFHEHD',6X,'CFHEEL')
C!1      NDMQ=NPT1
C!1      IF (NDIM.LT.NPT2) NDMQ=NDIM
C!1      DO 1231 K=1,NDMQ
C!1      WRITE (iUnitOutput,1180) K,ALTD(K),CFHEN2(K),CFHEO2(K),CFHEO(K),
C!1     $CFHEHE(K),CFHEH(K),CFHEOX(K),CFHEHD(K),CFHEEL(K)
C!11231  CONTINUE
C!1      IF (NDIM.LT.NPT2) GO TO 1291
C!1      NDMQ=NPT3
C!1      IF (NDIM.LT.NPT4) NDMQ=NDIM
C!1      DO 1241 K=NPT2,NDMQ,2
C!1      WRITE (iUnitOutput,1180) K,ALTD(K),CFHEN2(K),CFHEO2(K),CFHEO(K),
C!1     $CFHEHE(K),CFHEH(K),CFHEOX(K),CFHEHD(K),CFHEEL(K)
C!11241  CONTINUE
C!1      IF (NDIM.LT.NPT4) GO TO 1291
C!1      NDMQ=NPT5
C!1      IF (NDIM.LT.NPT6) NDMQ=NDIM
C!1      DO 1251 K=NPT4,NDMQ,5
C!1      WRITE (iUnitOutput,1180) K,ALTD(K),CFHEN2(K),CFHEO2(K),CFHEO(K),
C!1     $CFHEHE(K),CFHEH(K),CFHEOX(K),CFHEHD(K),CFHEEL(K)
C!11251  CONTINUE
C!1      IF (NDIM.LT.NPT6) GO TO 1291
C!1      DO 1261 K=NPT6,NDIM,10
C!1      WRITE (iUnitOutput,1180) K,ALTD(K),CFHEN2(K),CFHEO2(K),CFHEO(K),
C!1     $CFHEHE(K),CFHEH(K),CFHEOX(K),CFHEHD(K),CFHEEL(K)
C!11261  CONTINUE
C!11291  CONTINUE

      WRITE (iUnitOutput,2247)
2247  FORMAT(1H1,50X,'COLLISION FREQUENCIES FOR HYDROGEN')
C!      WRITE(iUnitOutput,2248)
C!2248  FORMAT(12X,'ALT',7X,'CFHpH3p',7X,'CFHpH',7X,'CFHpEL',8X,
C!     ;'CFHpH2')
C!      
C!      NDMQ=NPT1
C!      
C!      IF (NDIM.LT.NPT2) NDMQ=NDIM
C!      
C!      do K=1,NDMQ
C!         WRITE (iUnitOutput,1180) K,ALTD(K),CFHpH3p(K),CFHpH(K),CFHpEL(K),CFHpH2(K)
C!      enddo
C!
C!      IF (NDIM.LT.NPT2) GO TO 2291
C!
C!      NDMQ=NPT3
C!
C!      IF (NDIM.LT.NPT4) NDMQ=NDIM
C!
C!      do K=NPT2,NDMQ,2
C!         WRITE (iUnitOutput,1180) K,ALTD(K),CFHpH3p(K),CFHpH(K),CFHpEL(K),CFHpH2(K)
C!      enddo
C!
C!      IF (NDIM.LT.NPT4) GO TO 2291
C!
C!      NDMQ=NPT5
C!
C!      IF (NDIM.LT.NPT6) NDMQ=NDIM
C!
C!      do K=NPT4,NDMQ,5
C!         WRITE (iUnitOutput,1180) K,ALTD(K),CFHpH3p(K),CFHpH(K),CFHpEL(K),CFHpH2(K)
C!      enddo
C!
C!      IF (NDIM.LT.NPT6) GO TO 2291
C!
C!      do K=NPT6,NDIM,10
C!         WRITE (iUnitOutput,1180) K,ALTD(K),CFHpH3p(K),CFHpH(K),CFHpEL(K),CFHpH2(K)
C!      enddo
C!2291  CONTINUE

      WRITE (iUnitOutput,3247)
3247  FORMAT(1H1,50X,'COLLISION FREQUENCIES FOR ELECTRONS')
C!      WRITE(iUnitOutput,3248)
C!3248  FORMAT(12X,'ALT',7X,'CFELHp',6X,'CFELH3p',6X,'CFELH2',7X,
C!     ;'CFELH')
C! 
C!      NDMQ=NPT1
C!      
C!      IF (NDIM.LT.NPT2) NDMQ=NDIM
C!
C!      do K=1,NDMQ
C!         WRITE (iUnitOutput,1180) K,ALTD(K),CFELHp(K),CFELH3p(K),CFELH2(K),CFELH(K)
C!      enddo
C!      IF (NDIM.LT.NPT2) GO TO 3291
C!      NDMQ=NPT3
C!      IF (NDIM.LT.NPT4) NDMQ=NDIM
C!      DO 3241 K=NPT2,NDMQ,2
C!      WRITE (iUnitOutput,1180) K,ALTD(K),CFELN2(K),CFELO2(K),CFELO(K),
C!     $CFELHE(K),CFELH(K),CFELOX(K),CFELHL(K),CFELHD(K)
C!3241  CONTINUE
C!      IF (NDIM.LT.NPT4) GO TO 3291
C!      NDMQ=NPT5
C!      IF (NDIM.LT.NPT6) NDMQ=NDIM
C!      DO 3251 K=NPT4,NDMQ,5
C!      WRITE (iUnitOutput,1180) K,ALTD(K),CFELN2(K),CFELO2(K),CFELO(K),
C!     $CFELHE(K),CFELH(K),CFELOX(K),CFELHL(K),CFELHD(K)
C!3251  CONTINUE
C!      IF (NDIM.LT.NPT6) GO TO 3291
C!      DO 3261 K=NPT6,NDIM,10
C!      WRITE (iUnitOutput,1180) K,ALTD(K),CFELN2(K),CFELO2(K),CFELO(K),
C!     $CFELHE(K),CFELH(K),CFELOX(K),CFELHL(K),CFELHD(K)
C!3261  CONTINUE
C!3291  CONTINUE

      WRITE (iUnitOutput,1013)
1013  FORMAT(1H1,45X,'INITIAL OXYGEN PARAMETERS')
      WRITE(iUnitOutput,1021)
1021  FORMAT(16X,'ALT',10X,'VELOCITY',8X,'MACH NO',9X,'DENSITY',9X,
     ;'PRESSURE',6X,'TEMPERATURE',/)
      K=0
CALEX XM stands for Mach Number      
      XM=State_GV(0,uH3_)/SoundSpeed_GI(0,Ion1_)
      DNS1=State_GV(0,RhoH3_)/Mass_I(Ion1_)
      WRITE(iUnitOutput,1022) K,ALTMIN,State_GV(0,uH3_),XM,DNS1,State_GV(0,pH3_),State_GV(0,Th3_)
      NDMQ=NPT1
      IF (NDIM.LT.NPT2) NDMQ=NDIM
      DO 630 K=1,NDMQ
      US=SQRT(GAMMA*State_GV(K,pH3_)/State_GV(K,RhoH3_))
      XM=State_GV(K,uH3_)/US
      DNS1=State_GV(K,RhoH3_)/Mass_I(Ion1_)
      WRITE(iUnitOutput,1022) K,ALTD(K),State_GV(K,uH3_),XM,DNS1,State_GV(K,pH3_),State_GV(K,Th3_)
630   CONTINUE
      IF (NDIM.LT.NPT2) GO TO 690
      NDMQ=NPT3
      IF (NDIM.LT.NPT4) NDMQ=NDIM
      DO 640 K=NPT2,NDMQ,2
      US=SQRT(GAMMA*State_GV(K,pH3_)/State_GV(K,RhoH3_))
      XM=State_GV(K,uH3_)/US
      DNS1=State_GV(K,RhoH3_)/Mass_I(Ion1_)
      WRITE(iUnitOutput,1022) K,ALTD(K),State_GV(K,uH3_),XM,DNS1,State_GV(K,pH3_),State_GV(K,Th3_)
640   CONTINUE
      IF (NDIM.LT.NPT4) GO TO 690
      NDMQ=NPT5
      IF (NDIM.LT.NPT6) NDMQ=NDIM
      DO 650 K=NPT4,NDMQ,5
      US=SQRT(GAMMA*State_GV(K,pH3_)/State_GV(K,RhoH3_))
      XM=State_GV(K,uH3_)/US
      DNS1=State_GV(K,RhoH3_)/Mass_I(Ion1_)
      WRITE(iUnitOutput,1022) K,ALTD(K),State_GV(K,uH3_),XM,DNS1,State_GV(K,pH3_),State_GV(K,Th3_)
650   CONTINUE
      IF (NDIM.LT.NPT6) GO TO 690
      DO 660 K=NPT6,NDIM,10
      US=SQRT(GAMMA*State_GV(K,pH3_)/State_GV(K,RhoH3_))
      XM=State_GV(K,uH3_)/US
      DNS1=State_GV(K,RhoH3_)/Mass_I(Ion1_)
      WRITE(iUnitOutput,1022) K,ALTD(K),State_GV(K,uH3_),XM,DNS1,State_GV(K,pH3_),State_GV(K,Th3_)
660   CONTINUE
690   CONTINUE
      K=NDIM1
      XM=State_GV(nDim+1,uH3_)/SoundSpeed_GI(nDim+1,Ion1_)
      DNS1=State_GV(nDim+1,RhoH3_)/Mass_I(Ion1_)
      WRITE(iUnitOutput,1022) K,ALTMAX,State_GV(nDim+1,uH3_),XM,DNS1,State_GV(nDim+1,pH3_),State_GV(nDim+1,Th3_)
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
      XM=State_GV(0,uE_)/SoundSpeed_GI(0,nIon)
      DNS1=State_GV(0,RhoE_)/Mass_I(nIon)
      WRITE(iUnitOutput,1022) K,ALTMIN,State_GV(0,uE_),XM,DNS1,State_GV(0,pE_),State_GV(0,Te_)
      NDMQ=NPT1
      IF (NDIM.LT.NPT2) NDMQ=NDIM
      DO 830 K=1,NDMQ
      US=SQRT(GAMMA*State_GV(K,pE_)/State_GV(K,RhoE_))
      XM=State_GV(K,uE_)/US
      DNS1=State_GV(K,RhoE_)/Mass_I(nIon)
      WRITE(iUnitOutput,1022) K,ALTD(K),State_GV(K,uE_),XM,DNS1,State_GV(K,pE_),State_GV(K,Te_)
830   CONTINUE
      IF (NDIM.LT.NPT2) GO TO 890
      NDMQ=NPT3
      IF (NDIM.LT.NPT4) NDMQ=NDIM
      DO 840 K=NPT2,NDMQ,2
      US=SQRT(GAMMA*State_GV(K,pE_)/State_GV(K,RhoE_))
      XM=State_GV(K,uE_)/US
      DNS1=State_GV(K,RhoE_)/Mass_I(nIon)
      WRITE(iUnitOutput,1022) K,ALTD(K),State_GV(K,uE_),XM,DNS1,State_GV(K,pE_),State_GV(K,Te_)
840   CONTINUE
      IF (NDIM.LT.NPT4) GO TO 890
      NDMQ=NPT5
      IF (NDIM.LT.NPT6) NDMQ=NDIM
      DO 850 K=NPT4,NDMQ,5
      US=SQRT(GAMMA*State_GV(K,pE_)/State_GV(K,RhoE_))
      XM=State_GV(K,uE_)/US
      DNS1=State_GV(K,RhoE_)/Mass_I(nIon)
      WRITE(iUnitOutput,1022) K,ALTD(K),State_GV(K,uE_),XM,DNS1,State_GV(K,pE_),State_GV(K,Te_)
850   CONTINUE
      IF (NDIM.LT.NPT6) GO TO 890
      DO 860 K=NPT6,NDIM,10
      US=SQRT(GAMMA*State_GV(K,pE_)/State_GV(K,RhoE_))
      XM=State_GV(K,uE_)/US
      DNS1=State_GV(K,RhoE_)/Mass_I(nIon)
      WRITE(iUnitOutput,1022) K,ALTD(K),State_GV(K,uE_),XM,DNS1,State_GV(K,pE_),State_GV(K,Te_)
860   CONTINUE
890   CONTINUE
      K=NDIM1
      XM=State_GV(nDim+1,uE_)/SoundSpeed_GI(nDim+1,nIon)
      DNS1=State_GV(nDim+1,RhoE_)/Mass_I(nIon)
      WRITE(iUnitOutput,1022) K,ALTMAX,State_GV(nDim+1,uE_),XM,DNS1,State_GV(nDim+1,pE_),State_GV(nDim+1,Te_)
      WRITE(iUnitOutput,1024)
1024  FORMAT(1H1,40X,'INITIAL ELECTRIC FIELD AND SOURCE PARAMETERS')
      WRITE(iUnitOutput,1025)
1025  FORMAT(13X,'ALT',6X,'EFIELD',6X,'FCLSNO',6X,'ECLSNO',
     ;6X,'ECLSHE',6X,'FCLSNH',6X,'ECLSNH',6X,'FCLSNE',6X,'ECLSNE')
      NDMQ=NPT1
      IF (NDIM.LT.NPT2) NDMQ=NDIM
      DO 930 K=1,NDMQ
      WRITE(iUnitOutput,1026) K,ALTD(K),EFIELD(K),Source_CV(K,uH3_),Source_CV(K,pH3_),
     ;Source_CV(K,uH_),Source_CV(K,pH_),Source_CV(K,uE_),Source_CV(K,pE_)
930   CONTINUE
      IF (NDIM.LT.NPT2) GO TO 990
      NDMQ=NPT3
      IF (NDIM.LT.NPT4) NDMQ=NDIM
      DO 940 K=NPT2,NDMQ,2
      WRITE(iUnitOutput,1026) K,ALTD(K),EFIELD(K),Source_CV(K,uH3_),Source_CV(K,pH3_),
     ;Source_CV(K,uH_),Source_CV(K,pH_),Source_CV(K,uE_),Source_CV(K,pE_)
940   CONTINUE
      IF (NDIM.LT.NPT4) GO TO 990
      NDMQ=NPT5
      IF (NDIM.LT.NPT6) NDMQ=NDIM
      DO 950 K=NPT4,NDMQ,5
      WRITE(iUnitOutput,1026) K,ALTD(K),EFIELD(K),Source_CV(K,uH3_),Source_CV(K,pH3_),
     ;Source_CV(K,uH_),Source_CV(K,pH_),Source_CV(K,uE_),Source_CV(K,pE_)
950   CONTINUE
      IF (NDIM.LT.NPT6) GO TO 990
      DO 960 K=NPT6,NDIM,10
      WRITE(iUnitOutput,1026) K,ALTD(K),EFIELD(K),Source_CV(K,uH3_),Source_CV(K,pH3_),
     ;Source_CV(K,uH_),Source_CV(K,pH_),Source_CV(K,uE_),Source_CV(K,pE_)
960   CONTINUE
990   CONTINUE
999   CONTINUE
      WRITE (iUnitOutput,2013)
2013  FORMAT(1H1,45X,'HEAT CONDUCTIVITIES')
      WRITE(iUnitOutput,2021)
2021  FORMAT(16X,'ALT',10X,'OXYGEN',9X,'HYDROGEN',9X,
     ;'ELECTRONS'/)
      K=0
      WRITE(iUnitOutput,1022) K,ALTMIN,HeatCon_GI(0,Ion1_),HeatCon_GI(0,Ion2_),HeatCon_GI(0,nIon)
      NDMQ=NPT1
      IF (NDIM.LT.NPT2) NDMQ=NDIM
      DO 2630 K=1,NDMQ
      WRITE(iUnitOutput,1022) K,ALTD(K),HeatCon_GI(K,Ion1_),HeatCon_GI(K,Ion2_),HeatCon_GI(K,nIon)
2630  CONTINUE
      IF (NDIM.LT.NPT2) GO TO 2690
      NDMQ=NPT3
      IF (NDIM.LT.NPT4) NDMQ=NDIM
      DO 2640 K=NPT2,NDMQ,2
      WRITE(iUnitOutput,1022) K,ALTD(K),HeatCon_GI(K,Ion1_),HeatCon_GI(K,Ion2_),HeatCon_GI(K,nIon)
2640  CONTINUE
      IF (NDIM.LT.NPT4) GO TO 2690
      NDMQ=NPT5
      IF (NDIM.LT.NPT6) NDMQ=NDIM
      DO 2650 K=NPT4,NDMQ,5
      WRITE(iUnitOutput,1022) K,ALTD(K),HeatCon_GI(K,Ion1_),HeatCon_GI(K,Ion2_),HeatCon_GI(K,nIon)
2650   CONTINUE
      IF (NDIM.LT.NPT6) GO TO 2690
      DO 2660 K=NPT6,NDIM,10
      WRITE(iUnitOutput,1022) K,ALTD(K),HeatCon_GI(K,Ion1_),HeatCon_GI(K,Ion2_),HeatCon_GI(K,nIon)
2660   CONTINUE
2690   CONTINUE
      K=NDIM1
      WRITE(iUnitOutput,1022) K,ALTMAX,HeatCon_GI(nDim+1,Ion1_),HeatCon_GI(nDim+1,Ion2_),HeatCon_GI(nDim+1,nIon)
1019  FORMAT(3X,I3,0PF10.2,2X,11(1PE11.2))
1022  FORMAT(3X,I3,0PF14.2,2X,6(1PE17.5E3))
1026  FORMAT(3X,I3,0PF11.2,2X,9(1PE13.4E3))
1180  FORMAT(3X,I3,0PF10.2,2X,8(1PE13.4E3))
1028  FORMAT(3X,I3,0PF14.2,2X,2(1PE17.5E3))
1031  FORMAT(3X,I3,0PF14.2,2X,4(1PE17.5E3))

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
      use ModPWOM, ONLY:IsRestart
!     REAL jp1,jp2,jp3,jp4,kc1,kc2,kc3,kc6,kc7,kc8,kr1,kr2
      real DensityHp,DensityH3p,DensityH2p
C     ALEX define the reaction rates, label by reaction number
C     ALEX j is for photochemistry, k is for regular chemistry
      
c     ! H2+hnu     --> H+ + H + e
      jp1=1.9E-11
!      jp1=9.5E-11
      !             --> H2+ + e
      jp2=9.9E-10
      !jp2=5.4E-10
      ! H+hnu      --> H+
      jp3=1.0E-9
      !jp3=7.3E-10
      ! H2O+hnu    --> H+ + OH +e
      jp4=4.2E-10
      !jp4=1.3E-10

      ! H2+ + H2    --> H3+ +H
      kc1=2.E-9
      ! H+ + H2 + M --> H3+ + M
      kc2=3.2E-29
      !kc2=0.0
      
      ! H+ + CH4    --> CH3+ + H2  
      !             --> CH4+ + H
      kc3=4.5E-9
      !kc3=4.15E-9
      !kc3=0.0
      ! H3+ + CH4   --> CH5+ + H2 
      kc6=2.4E-9
      !kc6=0.0
      ! H3+ + H2O   --> H3O+ + H2
      kc7=5.3E-9
      !kc7=0.0
      ! H+ + H2O    --> H2O+ + H
      kc8=8.2E-9
      !kc8=0.0
      ! H+ + e      --> H + hnu
      kr1=1.91E-10
      !kr1=2.0E-12
      !kr1=0.0
      ! H3+ + e     --> H2 + H
      !             --> 3H
      kr2=1.73E-6
      !kr2=4.6E-6
!note kr2 includes H3+ + e --> H2 + H and --> 3H

      ! H+ + H2(nu>3) --> H2+ + H
      do i=1,nDim 
         call get_rate(ALTD(i)/100000.0,kc9(i))
      enddo
      
      
      !kc9(:)=0.0*kc9(:)
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
      State_GV(-1:0,uH3_)=0.
      State_GV(-1:0,uH_)=0.
      State_GV(-1:0,uH2_)=0.
      State_GV(-1:0,uE_)=0.
      CALL MODATM(ALTMIN,XNH2,XNH,XNH2O,XNCH4,TEMP)

CALEX I pretend that for plasma parameters, O is H3 and HE is
CALEX chemical equilibrium value for H2+ this allow me to just
CALEX change the chemistry but leave the rest of the code the same  
      State_GV(-1:0,Th3_)=TEMP
      State_GV(-1:0,Th_)=TEMP
      State_GV(-1:0,Th2_)=TEMP
      State_GV(-1:0,Te_)=TEMP

C KGS this subroutine needs to be modified
      call calc_chemical_equilibrium(DensityHp,DensityH3p,DensityH2p)
      State_GV(-1:0,RhoH3_)=Mass_I(Ion1_)*DensityH3p
      State_GV(-1:0,RhoH_)=Mass_I(Ion2_)*DensityHp
      State_GV(-1:0,RhoH2_)=Mass_I(Ion3_)*DensityH2p
      write(*,*) 'H+(1400km)=',DensityHp,', H3+(1400km)=',DensityH3p,
     &     ', H2+(1400km)=',DensityH2p


C I have used numerically calculated chemical equilibrium
C solution for T=800k.       
!      State_GV(0,RhoH3_)=Mass_I(Ion1_)*6489.69
!c      State_GV(0,RhoHe_)=Mass_I(Ion3_)*jp2/kc1
!      State_GV(0,RhoH_)E=0.
!      State_GV(0,RhoH_)=Mass_I(Ion2_)*1343.64
C I have used numerically calculated chemical equilibrium
C solution for T=1000k.       
c      State_GV(0,RhoH3_)=Mass_I(Ion1_)*4725.0
c      State_GV(0,RhoH_)E=0.
c      State_GV(0,RhoH_)=Mass_I(Ion2_)*368.0

C I have used numerically calculated chemical equilibrium
C solution for T=1500k.       
c      State_GV(0,RhoH3_)=Mass_I(Ion1_)*560.0
c      State_GV(0,RhoH_)E=0.
c      State_GV(0,RhoH_)=Mass_I(Ion2_)*30.0


C I have used numerically calculated chemical equilibrium
C solution for T=1500k. with enhanced water and decreased CH4      
c      State_GV(0,RhoH3_)=Mass_I(Ion1_)*11435.41
c      State_GV(0,RhoH_)E=0.
c      State_GV(0,RhoH_)=Mass_I(Ion2_)*1463.48

C I have used numerically calculated chemical equilibrium
C solution for T=100k. with reduced CH4 enhanced h2o      
C      State_GV(0,RhoH3_)=Mass_I(Ion1_)*5509.0
C      State_GV(0,RhoH_)E=0.
C      State_GV(0,RhoH_)=Mass_I(Ion2_)*1124.0


      State_GV(-1:0,RhoE_)=MassElecIon_I(Ion3_)*State_GV(-1:0,RhoH2_) +
     &     MassElecIon_I(Ion2_)*State_GV(-1:0,RhoH_) + 
     &     MassElecIon_I(Ion1_)*State_GV(-1:0,RhoH3_)
      State_GV(-1:0,pH3_)=RGAS_I(Ion1_)*State_GV(-1:0,Th3_)*State_GV(-1:0,RhoH3_)
      State_GV(-1:0,pH_)=RGAS_I(Ion2_)*State_GV(-1:0,Th_)*State_GV(-1:0,RhoH_)
      State_GV(-1:0,pH2_)=RGAS_I(Ion3_)*State_GV(-1:0,Th2_)*State_GV(-1:0,RhoH2_)
      State_GV(-1:0,pE_)=RGAS_I(nIon)*State_GV(-1:0,Te_)*State_GV(-1:0,RhoE_)
      SoundSpeed_GI(0,Ion1_)=SQRT(GAMMA*RGAS_I(Ion1_)*State_GV(0,Th3_))
      SoundSpeed_GI(0,Ion2_)=SQRT(GAMMA*RGAS_I(Ion2_)*State_GV(0,Th_))
      SoundSpeed_GI(0,Ion3_)=SQRT(GAMMA*RGAS_I(Ion3_)*State_GV(0,Th2_))
      SoundSpeed_GI(0,nIon)=SQRT(GAMMA*RGAS_I(nIon)*State_GV(0,Te_))

C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     

! Start collision freq at zero
      CollisionFreq_IIC(:,:,:)= 0.0
      HeatFlowCoef_II(:,:)    = 0.0
      FricHeatCoef_II(:,:)    = 0.0


      DO 20 I=1,NDIM
CALEX SOURCE COEF?
C jp = photochemical reaction rate
C kc = collisional reaction rate      
         FFHpp1(I)=jp1*XH2(I)
         FFHpp3(I)=jp3*XH(I)
         FFHpp4(I)=jp4*XH2O(I)
         FFHpc2(I)=-kc2*XH2(I)*XH2(I)
         FFHpc3(I)=-kc3*XCH4(I)
         FFHpc8(I)=-kc8*XH2O(I)
         FFHpc9(I)=-kc9(I)*XH2(I)
         FFHpr1(I)=-kr1
CALEX write out source coeff
CALEX         write(26,*) FFHpp1(I),FFHpp3(I),FFHpp4(I),FFHpc2(I),FFHpc3(I),FFHpc8(I),FFHpr1(I)
        

!         FFH3pc1(I)=kc1*(jp2/kc1)*XH2(I)   ! KGS H2+ folded in here
         FFH3pc1(I)=kc1*XH2(I)              ! KGS Make sure this is multiplied by H2+ now
         FFH3pc2(I)=kc2*XH2(I)*XH2(I)
         FFH3pc6(I)=-kc6*XCH4(I)
         FFH3pc7(I)=-kc7*XH2O(I)
         FFH3pr2(I)=-kr2

CALEX write out source coeff
CALEX         write(27,*) FFH3pc1(I),FFH3pc2(I),FFH3pc6(I),FFH3pc7(I),FFH3pr2(I)

         FFH2pp2(I)=jp2*XH2(I)
         FFH2pc9(I)=kc9(I)*XH2(I)
         FFH2pc1(I)=-kc1*XH2(I)

CALEX CL=COLLISION COEF, CF=collision freq ?         
C KGS CF(1,2) = CL(1,2) * n(2)
C KGS n(2)m(2)CF(2,1) = n(1)m(1)CF(1,2)
CALEX the coulomb collisions
CAlex H+ and H3+
c KGS and H2+   from eq 4.143?
         CLHpH3p(I)=1.905*4.**1.5/Mass_I(Ion1_) ! KGS orig
         CLHpH3p(I)=1.100*4.**1.5/Mass_I(Ion1_) ! KGS correct value?
         CLH2pH3p(I)=3.17*6.**0.5/Mass_I(Ion1_) ! KGS double-check
         CLH2pHp(I)=3.81*2.**0.5/Mass_I(Ion2_)  ! KGS double-check
CALEX electron H+ and electron H3+
c KGS and electron H2+
         CLELHp(I)=54.5/Mass_I(Ion2_)
         CLELH3p(I)=54.5/Mass_I(Ion1_)
         CLELH2p(I)=54.5/Mass_I(Ion3_)

CALEX  ion neutrals
C KGS looks similar to Nagy 4.88 Maxwell Molecule Collisions
C KGS but some terms seem to be missing
         ! H+ - H2
         CLHpH(I)=2.65E-10*XH(I)
!         CFHpH2(I)=2.6E-9*XH2(I)*(.82/.667)**.5
         CollisionFreq_IIC(Ion2_,Neutral1_,I)=
     &        2.6E-9*XH2(I)*(.82/.667)**.5

! KGS 2-2 Coulomb Collision is in collisionPW.f

         ! H3+ - H
!         CFH3pH(I)=2.6E-9*XH(I)*(.667/.75)**.5
         CollisionFreq_IIC(Ion1_,Neutral2_,I)=
     &        2.6E-9*XH(I)*(.667/.75)**.5

         ! H3+ - H2
!         CFH3pH2(I)=2.6E-9*XH2(I)*(.82/1.2)**.5
         CollisionFreq_IIC(Ion1_,Neutral1_,I)=
     &        2.6E-9*XH2(I)*(.82/1.2)**.5

C KGS make sure I should be using Maxwell Molecule Collisions here
         ! H2+ - H2
         CollisionFreq_IIC(Ion3_,Neutral1_,I)=
     &        2.6E-9*XH2(I)*(.82/1.0)**.5    ! KGS check these
         ! H2+ - H
         CollisionFreq_IIC(Ion3_,Neutral2_,I)=
     &        2.6E-9*XH(I)*(.667/.667)**.5

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
      !CTHpH=3.*RGAS_I(Ion2_)*Mass_I(Ion2_)/(Mass_I(Ion2_)+Mass_I(Ion2_))
      !CTHpH2=3.*RGAS_I(Ion2_)*Mass_I(Ion2_)/(Mass_I(Ion2_)+2.*Mass_I(Ion2_))
      !CTHpH3p=3.*RGAS_I(Ion2_)*Mass_I(Ion2_)/(Mass_I(Ion2_)+Mass_I(Ion1_))
      !CTHpEL=3.*RGAS_I(Ion2_)*Mass_I(Ion2_)/(Mass_I(Ion2_)+Mass_I(nIon))
      
      HeatFlowCoef_II(Ion3_,Neutral2_)=3.*RGAS_I(Ion3_)*Mass_I(Ion3_)/(Mass_I(Ion3_)+Mass_I(Ion2_))
      HeatFlowCoef_II(Ion3_,Neutral1_)=3.*RGAS_I(Ion3_)*Mass_I(Ion3_)/(Mass_I(Ion3_)+Mass_I(Ion3_))
      HeatFlowCoef_II(Ion3_,Ion1_)=3.*RGAS_I(Ion3_)*Mass_I(Ion3_)/(Mass_I(Ion3_)+Mass_I(Ion1_))
      HeatFlowCoef_II(Ion3_,Ion2_)=3.*RGAS_I(Ion3_)*Mass_I(Ion3_)/(Mass_I(Ion3_)+Mass_I(Ion2_))
      HeatFlowCoef_II(Ion3_,nIon)=3.*RGAS_I(Ion3_)*Mass_I(Ion3_)/(Mass_I(Ion3_)+Mass_I(nIon))

      MassFracCoef_II(Ion3_,:) = HeatFlowCoef_II(Ion3_,:) / (3.0*RGAS_I(Ion3_))

      HeatFlowCoef_II(Ion2_,Neutral2_)=3.*RGAS_I(Ion2_)*Mass_I(Ion2_)/(Mass_I(Ion2_)+Mass_I(Ion2_))
      HeatFlowCoef_II(Ion2_,Neutral1_)=3.*RGAS_I(Ion2_)*Mass_I(Ion2_)/(Mass_I(Ion2_)+2.*Mass_I(Ion2_))
      HeatFlowCoef_II(Ion2_,Ion1_)=3.*RGAS_I(Ion2_)*Mass_I(Ion2_)/(Mass_I(Ion2_)+Mass_I(Ion1_))
      HeatFlowCoef_II(Ion2_,Ion3_)=3.*RGAS_I(Ion2_)*Mass_I(Ion2_)/(Mass_I(Ion2_)+Mass_I(Ion3_))
      HeatFlowCoef_II(Ion2_,nIon)=3.*RGAS_I(Ion2_)*Mass_I(Ion2_)/(Mass_I(Ion2_)+Mass_I(nIon))

      MassFracCoef_II(Ion2_,:) = HeatFlowCoef_II(Ion2_,:) / (3.0*RGAS_I(Ion2_))
!      CTH3pH=3.*RGAS_I(Ion1_)*Mass_I(Ion1_)/(Mass_I(Ion1_)+Mass_I(Ion2_))
!      CTH3pH2=3.*RGAS_I(Ion1_)*Mass_I(Ion1_)/(Mass_I(Ion1_)+2.*Mass_I(Ion2_))
!      CTH3pHp=3.*RGAS_I(Ion1_)*Mass_I(Ion1_)/(Mass_I(Ion1_)+Mass_I(Ion2_))
!      CTH3pEL=3.*RGAS_I(Ion1_)*Mass_I(Ion1_)/(Mass_I(Ion1_)+Mass_I(nIon))

      HeatFlowCoef_II(Ion1_,Neutral2_)=3.*RGAS_I(Ion1_)*Mass_I(Ion1_)/(Mass_I(Ion1_)+Mass_I(Ion2_))
      HeatFlowCoef_II(Ion1_,Neutral1_)=3.*RGAS_I(Ion1_)*Mass_I(Ion1_)/(Mass_I(Ion1_)+2.*Mass_I(Ion2_))
      HeatFlowCoef_II(Ion1_,Ion2_)=3.*RGAS_I(Ion1_)*Mass_I(Ion1_)/(Mass_I(Ion1_)+Mass_I(Ion2_))
      HeatFlowCoef_II(Ion1_,Ion3_)=3.*RGAS_I(Ion1_)*Mass_I(Ion1_)/(Mass_I(Ion1_)+Mass_I(Ion3_))
      HeatFlowCoef_II(Ion1_,nIon)=3.*RGAS_I(Ion1_)*Mass_I(Ion1_)/(Mass_I(Ion1_)+Mass_I(nIon))

      MassFracCoef_II(Ion1_,:) = HeatFlowCoef_II(Ion1_,:) / (3.0*RGAS_I(Ion1_))
!      CTELH=3.*RGAS_I(nIon)*Mass_I(nIon)/(Mass_I(nIon)+Mass_I(Ion2_))
!      CTELH2=3.*RGAS_I(nIon)*Mass_I(nIon)/(Mass_I(nIon)+2.*Mass_I(Ion2_))
!      CTELHp=3.*RGAS_I(nIon)*Mass_I(nIon)/(Mass_I(nIon)+Mass_I(Ion2_))
!      CTELH3p=3.*RGAS_I(nIon)*Mass_I(nIon)/(Mass_I(nIon)+Mass_I(Ion1_))
      
      HeatFlowCoef_II(nIon,Neutral2_)=3.*RGAS_I(nIon)*Mass_I(nIon)/(Mass_I(nIon)+Mass_I(Ion2_))
      HeatFlowCoef_II(nIon,Neutral1_)=3.*RGAS_I(nIon)*Mass_I(nIon)/(Mass_I(nIon)+2.*Mass_I(Ion2_))
      HeatFlowCoef_II(nIon,Ion3_)=3.*RGAS_I(nIon)*Mass_I(nIon)/(Mass_I(nIon)+Mass_I(Ion3_))
      HeatFlowCoef_II(nIon,Ion2_)=3.*RGAS_I(nIon)*Mass_I(nIon)/(Mass_I(nIon)+Mass_I(Ion2_))
      HeatFlowCoef_II(nIon,Ion1_)=3.*RGAS_I(nIon)*Mass_I(nIon)/(Mass_I(nIon)+Mass_I(Ion1_))

      MassFracCoef_II(nIon,:) = HeatFlowCoef_II(nIon,:) / (3.0*RGAS_I(nIon))

CALEX CMOXN2 = M_{N2}/(M_o+M_{N2}) see nagy p.83
!      CMHpH=Mass_I(Ion2_)/(Mass_I(Ion2_)+Mass_I(Ion2_))
!      CMHpH2=2.*XAMU/(Mass_I(Ion2_)+2.*XAMU)
!      CMHpH3p=Mass_I(Ion1_)/(Mass_I(Ion2_)+Mass_I(Ion1_))
!      CMHpEL=Mass_I(nIon)/(Mass_I(Ion2_)+Mass_I(nIon))

      FricHeatCoef_II(Ion3_,Neutral2_)=Mass_I(Ion3_)/(Mass_I(Ion3_)+Mass_I(Ion2_))
      FricHeatCoef_II(Ion3_,Neutral1_)=2.*XAMU/(Mass_I(Ion3_)+2.*XAMU)
      FricHeatCoef_II(Ion3_,Ion2_)=Mass_I(Ion2_)/(Mass_I(Ion3_)+Mass_I(Ion2_))
      FricHeatCoef_II(Ion3_,Ion1_)=Mass_I(Ion1_)/(Mass_I(Ion3_)+Mass_I(Ion1_))
      FricHeatCoef_II(Ion3_,nIon)=Mass_I(nIon)/(Mass_I(Ion3_)+Mass_I(nIon))

      FricHeatCoef_II(Ion2_,Neutral2_)=Mass_I(Ion2_)/(Mass_I(Ion2_)+Mass_I(Ion2_))
      FricHeatCoef_II(Ion2_,Neutral1_)=2.*XAMU/(Mass_I(Ion2_)+2.*XAMU)
      FricHeatCoef_II(Ion2_,Ion3_)=Mass_I(Ion3_)/(Mass_I(Ion2_)+Mass_I(Ion3_))
      FricHeatCoef_II(Ion2_,Ion1_)=Mass_I(Ion1_)/(Mass_I(Ion2_)+Mass_I(Ion1_))
      FricHeatCoef_II(Ion2_,nIon)=Mass_I(nIon)/(Mass_I(Ion2_)+Mass_I(nIon))
       
!      CMH3pH=Mass_I(Ion2_)/(Mass_I(Ion1_)+Mass_I(Ion2_))
!      CMH3pH2=2.*XAMU/(Mass_I(Ion1_)+2.*XAMU)
!      CMH3pHp=Mass_I(Ion2_)/(Mass_I(Ion1_)+Mass_I(Ion2_))
!      CMH3pEL=Mass_I(nIon)/(Mass_I(Ion1_)+Mass_I(nIon))
       
      FricHeatCoef_II(Ion1_,Neutral2_)=Mass_I(Ion2_)/(Mass_I(Ion1_)+Mass_I(Ion2_))
      FricHeatCoef_II(Ion1_,Neutral1_)=2.*XAMU/(Mass_I(Ion1_)+2.*XAMU)
      FricHeatCoef_II(Ion1_,Ion3_)=Mass_I(Ion3_)/(Mass_I(Ion1_)+Mass_I(Ion3_))
      FricHeatCoef_II(Ion1_,Ion2_)=Mass_I(Ion2_)/(Mass_I(Ion1_)+Mass_I(Ion2_))
      FricHeatCoef_II(Ion1_,nIon)=Mass_I(nIon)/(Mass_I(Ion1_)+Mass_I(nIon))

!      CMELH=Mass_I(Ion2_)/(Mass_I(nIon)+Mass_I(Ion2_))
!      CMELH2=2.*XAMU/(Mass_I(nIon)+2.*XAMU)
!      CMELHp=Mass_I(Ion2_)/(Mass_I(nIon)+Mass_I(Ion2_))
!      CMELH3p=Mass_I(Ion1_)/(Mass_I(nIon)+Mass_I(Ion1_))
      
      FricHeatCoef_II(nIon,Neutral2_)=Mass_I(Ion2_)/(Mass_I(nIon)+Mass_I(Ion2_))
      FricHeatCoef_II(nIon,Neutral1_)=2.*XAMU/(Mass_I(nIon)+2.*XAMU)
      FricHeatCoef_II(nIon,Ion3_)=Mass_I(Ion3_)/(Mass_I(nIon)+Mass_I(Ion3_))
      FricHeatCoef_II(nIon,Ion2_)=Mass_I(Ion2_)/(Mass_I(nIon)+Mass_I(Ion2_))
      FricHeatCoef_II(nIon,Ion1_)=Mass_I(Ion1_)/(Mass_I(nIon)+Mass_I(Ion1_))
      
C     ALEX(10/11/04): 
C     TRY SETTING THE PLASMA PARAMETERS HERE TO THE SURFACE VALUES

      print *,'Restart: ',IsRestart

      if(IsRestart) RETURN
      IsRestart = .true.
      do K=1,NDIM
         State_GV(K,uH3_)=0
         State_GV(K,pH3_)=State_GV(0,pH3_)*exp(-(ALTD(k)-1400.E5)/5000.E5)
         State_GV(K,RhoH3_)=State_GV(0,RhoH3_)*exp(-(ALTD(k)-1400.E5)/5000.E5)
         State_GV(K,Th3_)=State_GV(0,Th3_)
         State_GV(K,uH_)=0
         State_GV(K,pH_)=State_GV(0,pH_)*exp(-(ALTD(k)-1400.E5)/5000.E5)
         State_GV(K,RhoH_)=State_GV(0,RhoH_)*exp(-(ALTD(k)-1400.E5)/5000.E5)
         State_GV(K,Th_)=State_GV(0,Th_)
         State_GV(K,uH2_)=0
         State_GV(K,pH2_)=State_GV(0,pH2_)*exp(-(ALTD(k)-1400.E5)/5000.E5)
         State_GV(K,RhoH2_)=State_GV(0,RhoH2_)*exp(-(ALTD(k)-1400.E5)/5000.E5)
         State_GV(K,Th2_)=State_GV(0,Th2_)
         State_GV(K,RhoE_)=State_GV(0,RhoE_)*exp(-(ALTD(k)-1400.E5)/5000.E5)
         State_GV(K,uE_)=0
         State_GV(K,pE_)=State_GV(0,pE_)*exp(-(ALTD(k)-1400.E5)/5000.E5)
         State_GV(K,Te_)=State_GV(0,Te_)
         
         
      enddo
      
      RETURN
      END
