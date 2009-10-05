
C Subroutine GLOWEX
C
C S.C. Solomon, 6/89; Adapted for polar wind model by R. Cannata, 9/89
C
C Supplied to subroutine in labeled common /CGLOW/
C ISCALE  =0 for contrast ratio scaling of solar flux, =1 for linear interp.
C JLOCAL  =0 for electron transport calculation, =1 for local effects only
C IDATE   Date, in form yyddd
C UT      Universal Time, seconds
C GLAT    Geographic latitude, degrees
C GLONG   Geographic longitude, degrees
C F107    Solar F10.7 flux index for day being modeled (not previous day)
C F107A   Solar F10.7 flux index 81-day centered average
C EFLUX   Auroral electron fluxes: hard, soft, cusp, drizzle; erg cm-2, s-1
C EZERO   Auroral electron flux characteristic energies: h, s, c, d; keV
C ZZ      altitude array, km  (NO !!!!!! CM !!!!!)
C ZO      O number density at each altitude, cm-3
C ZN2     N2  "      "      "   "     "       "
C ZO2     O2         "
C ZNO     NO         "
C ZNS     N(4S)      "
C ZND     N(2D)      "
C ZRHO    mass density at each altitude, gm cm-3
C ZE      electron density at each alt, cm-3
C ZTN     Neutral temperature at each alt, K
C ZTI     Ion temperature at each alt, K
C ZTE     Electron temp at each alt, K (if zero, Brace/Theis model fr. IRI used)
C PHITOP  Low energy electron flux array from conjugate point, cm-2 s-1 eV-1
C RAT(JMAX)  ALTITUDE SPECIFIC FRACTION OF IONIZATION RAT WHICH COMES FROM 
C         PHOTOELECTRONS 
C
C Calculated by subroutine:
C ECALC   calculated electron density assuming photochemical equilibrium
C         below 200 km and using given O+ density at or above 200 km
C ZXDEN   array of excited and and/or ionized state densities,
C         O+(2P), O+(2D), O+(4S), N2+, N+, O2+, NO+, N2(A), N(2P), N(2D),
C         O(1S), O(1D), 8 spares, at each altitude, cm-3
C ZETA    array of volume emission rates  at 3371A, 4278A, 5200A, 5577A,
C         6300A, 7320A, 14 spares, at each altitude, cm-3 s-1
C ZCETA   array of reaction contributions to each v.e.r at each alt, cm-3 s-1
C VCB     array of vertical column brightnesses (as above), Rayleighs
C SZA     solar zenith angle, radians
C DIP     magnetic field dip angle, radians
C EFRAC   energy conservation check from electron transport routine, (out-in)/in
C IERR    error code returned from electron transport routine, should = 0
C
C Array dimensions:
C JMAX    number of altitude levels
C NBINS   number of energetic electron energy bins
C LMAX    number of wavelength intervals for solar flux
C NMAJ    number of major species
C NEX     number of ionized/excited species
C NW      number of airglow emission wavelengths
C NC      number of component production terms for each emission
C NST     number of states produced by photoionization/dissociation
C NEI     number of states produced by electron impact
C NF      number of available types of auroral fluxes
C
C
      SUBROUTINE GLOWEX
C
      use ModGlow
      use ModCommonVariables, ONLY: MaxGrid,DoLog,AltMin,GLAT,GLONG,F107,
     &     F107A,IYD,SEC,iUnitOutput, GMLAT
      use ModNumConst, ONLY: cRadToDeg,cPi

C      PARAMETER (JMAX=92)
C      PARAMETER (NBINS=84)
C      PARAMETER (LMAX=59)
C      PARAMETER (NMAJ=3)
C      PARAMETER (NEX=20)
C      PARAMETER (NW=20)
C      PARAMETER (NC=10)
C      PARAMETER (NST=6)
C      PARAMETER (NEI=10)
C      PARAMETER (NF=4)
C
      REAL ZVCD(NMAJ,JMAX)
c      DIMENSION PHOTOTF(601)
c      DIMENSION ALTD(601),RAD(601),RBOUND(601),photott(601)
      REAL photott(MaxGrid)
C      COMMON /CGLOW/
C     >  EFLUX(NF), EZERO(NF),
C     >  SZA, DIP, ISCALE, JLOCAL, EFRAC, IERR,
C     >  ZO(JMAX), ZN2(JMAX), ZO2(JMAX), ZNO(JMAX),
C     >  ZNS(JMAX), ZND(JMAX), ZRHO(JMAX), ZE(JMAX),
C     >  ZCOL(NMAJ,JMAX),ZTN(JMAX),ZMAJ(NMAJ,JMAX),ZZ(JMAX),
C     >  ZTI(JMAX), ZTE(JMAX),
C     >  WAVE1(LMAX), WAVE2(LMAX), SFLUX(LMAX),
C     >  ENER(NBINS), DEL(NBINS), PHITOP(NBINS),
C     >  PESPEC(NBINS,JMAX), SESPEC(NBINS,JMAX),
C     >  PHOTOI(NST,NMAJ,JMAX),PHOTOD(NST,NMAJ,JMAX),PHONO(NST,JMAX),
C     >  QTI(JMAX),AURI(NMAJ,JMAX),PIA(NMAJ,JMAX),SION(NMAJ,JMAX),
C     >  UFLX(NBINS,JMAX), DFLX(NBINS,JMAX), AGLW(NEI,NMAJ,JMAX),
C     >  EHEAT(JMAX), TEZ(JMAX), ECALC(JMAX),
C     >  ZXDEN(NEX,JMAX), ZETA(NW,JMAX), ZCETA(NC,NW,JMAX),VCB(NW),
C     >  NNN(NMAJ)
C
CALEX
      
CALEX origionally iyd was idate but this is inconsistant with 
CALEX common_variables.f
C      COMMON /MSIS1/ IDATE,UT,SEC,GLAT,GLONG,STL,F107A,F107,AP(7),IART
C     $,GMLAT,GMLONG
c      COMMON /SPCE/ ALTMIN,ALTMAX,ALTD,RAD,RBOUND,DRBND,DTR1,DTR2
c      COMMON /CMDN/ NDIM,NDIM1,DT,DTMX,TIME,TMAX,NSTEP,NPRINT,NSTPMX
c     ;,NDIM2,NDIMM,DTX1,DTX2
c      COMMON /PHOT/PHOTOTF(601),phototp(601)
c      COMMON /INTG/IFACTOR
      DATA IFIRST/1/
      SAVE IFIRST
      ifirst = 1

C
C UT GENERATED IN SUB MODATM IN HOURS...GLOWEX NEEDS UT IN SEC SO DEFINE "UTG"
      UTG = SEC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      I = 0
      J = 0
C  IFACTOR SCALES THE ALTITUDE INTERVAL FOR ATTENUATION CALCULATIONS
C  IN TERMS OF THE CELL THICKNESS ADOPTED IN MAIN PROGRAM
C  DMX MUST BE A POSITIVE EVEN INTEGER !
       IFACTOR = 2E5
       IDMX =INT(DrGlow/IFACTOR)
       DO K = 90,1000,IDMX
          I = I + 1
          ZKK = K*1.E05
          ZZ(I) = ZKK
       enddo 
       DO I =1,JMAX
          J = J+1     
          CALL MODATM(ZZ(I),ZO2(J),ZN2(J),ZO(J),XXX,XXX,ZTN(J))
       enddo
C
       do I = 1,JMAX
          ZMAJ(1,I) = ZO(I)
          ZMAJ(2,I) = ZO2(I)
          ZMAJ(3,I) = ZN2(I)
       enddo
       ISCALE = 0
       JLOCAL = 1
       J = 0
C     
C     First call only, set up energy grid and scale solar flux:
C     
       IF (IFIRST .EQ. 1) THEN
C     
          IFIRST = 0
          DO N=1,NBINS
             IF (N .LE. 21) THEN
                ENER(N) = 0.5 * FLOAT(N)
             ELSE
                ENER(N) = EXP (0.05 * FLOAT(N+26))
             ENDIF
          enddo
          DEL(1) = 0.5
          DO N=2,NBINS
             DEL(N) = ENER(N)-ENER(N-1)
          enddo
          DO N=1,NBINS
             ENER(N) = ENER(N) - DEL(N)/2.0
          enddo
C     
          CALL SSFLUX (ISCALE, F107, F107A, WAVE1, WAVE2, SFLUX, LMAX)
C     
       ENDIF
C     
C     
C     Find magnetic dip angle and solar zenith angle (radians):
C     
       CALL FIELDM (GLAT, GLONG, 300., XF, YF, ZF, FF, DIP, DEC, SDIP)
       DIP = ABS(DIP) * cPi/180.
C     
       CALL SOLZEN (IYD, UTG, GLAT, GLONG, SZA)
       SZA = min(SZA,85.0)
       SZA = SZA * cPi/180.
       SZAD = SZA*cRadToDeg
       if (DoLog) WRITE(iUnitOutput,9996)SZAD
 9996  FORMAT(26X,'SOLAR ZENITH ANGLE:',1F8.2,' DEGREES')
C     Calculate slant path column densities in the direction of sun of major
C     species:
C     
       
       CALL RCOLUM (SZA, ZZ, ZMAJ, ZTN, ZCOL, ZVCD, JMAX, NMAJ)
C     
C     Call Subroutine EPHOTO to calculate the photoelectron production
C     spectrum and photoionization rates as a function of altitude:
C     
       CALL EPHOTO
C     
       if (DoLog) WRITE(iUnitOutput,9997)
C     ALEX note IDATE HERE is IYD in the rest of the subroutine
 9997  FORMAT('   IDATE ','  UTG   ','   SEC   ','    GLAT    ',
     $      ' GLONG    ',' STL   ',' F107A   ','  F107   ',
     $      '   AP(3)',' IART   ',' GMLAT   ','   GMLONG   ',
     $      '   SZAD    ')
C     
!       if (DoLog) WRITE(iUnitOutput,9995)IYD,UTG,SEC,GLAT,GLONG,STL,F107A,F107,AP(3)
!     $      ,IART,GMLAT,GMLONG,SZAD
 9995  FORMAT(1X,I7,8(1X,F8.2),2X,I2,4X,F7.2,4X,F7.2,5X,F7.2)
C     
C     
!       if (GMLAT .GE. 75.0) then 
          call PRECIP(IDMX,ALTMIN,PHOTOTP)
!       else
!          PHOTOTP(:) =0.0
!       endif
C     
C     DO 109 J=291,391
C     PHOTOTP(J)=PHOTOTP(290)
       DO J=1,nCellGlow
          if (J > 14) PHOTOTP(J)=PHOTOTP(14)
       enddo
C     NOW ADD EUV AND PARTICLE IONIZATIONS!!!!!!1
       DO I = 1,nCellGlow-1
          PHOTOTT(I) = PHOTOTF(I) + PHOTOTP(I)
          
C     NOW RENAME TO MATCH MAIN PROGRAM
          PHOTOTF(I) = PHOTOTT(I)      
C     
       enddo
       PHOTOTF(nCellGlow) = PHOTOTT(nCellGlow-1)      
       RETURN
       END
      
