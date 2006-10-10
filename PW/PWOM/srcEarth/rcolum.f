C
C Subroutine RCOLUM Calculates the column density ZCOL for each species
C ZMAJ above height ZZ at solar zenith angle CHI.  Uses Chapman function
C fit to account for spherical geometry.  If CHI is less than 90 degrees,
C column densities are calculated directly; if CHI is greater than
C 90 degrees the column density at grazing height for 90 degrees is calculated
C and doubled and the column density above ZZ(J) is subtracted; if CHI is
C such that the grazing height is less than the radius of the earth the column
C densities are set to 'infinity', i.e., 1.0E30.  Densities supplied in
C array ZMAJ are used in the calculation except where grazing height is
C below the lowest level specified, in this case values are interpolated
C logarithmically from the US standard atmosphere at sea level, the tropopause,
C the stratopause, and the mesopause.
C
C
      SUBROUTINE RCOLUM (CHI, ZZ, ZMAJ, TN, ZCOL, ZVCD, JMAX, NMAJ)
C
      PARAMETER (NM=3)
      PARAMETER (NU=4)
C
      DIMENSION ZZ(JMAX),ZMAJ(NMAJ,JMAX),TN(JMAX),ZCOL(NMAJ,JMAX),
     >      ZVCD(NMAJ,JMAX), ZCG(NM), ZUS(NU), TNUS(NU), ZCUS(NM,NU)
C
      DATA PI/3.1415926535/, RE/6.37E8/
      DATA ZUS/0., 1.5E6, 5.E6, 9.E6/, TNUS/288., 217., 271., 187./
      DATA ZCUS/8.00E17, 4.54E24, 1.69E25,
     >          8.00E17, 5.46E23, 2.03E24,
     >          8.00E17, 3.63E21, 1.35E22,
     >          7.80E17, 8.48E18, 3.16E19/
C
      CALL VCD (ZZ, ZMAJ, ZVCD, JMAX, NMAJ)
C
      IF (CHI .GE. 2.) THEN 
        DO 40 I=1,NMAJ
        DO 40 J=1,JMAX
        ZCOL(I,J) = 1.0E30
   40   CONTINUE
        RETURN
      ENDIF
C
      IF (CHI .LE. PI/2.) THEN
        DO 60 I=1,NMAJ
        DO 60 J=1,JMAX
        ZCOL(I,J) = ZVCD(I,J) * CHAP(CHI,ZZ(J),TN(J),I)
   60   CONTINUE
      ELSE
        DO 220 J=1,JMAX
        GHRG=(RE+ZZ(J))*SIN(CHI) 
        GHZ=GHRG-RE 
        IF (GHZ .LE. 0.) THEN
          DO 80 I=1,NMAJ
          ZCOL(I,J) = 1.0E30
   80     CONTINUE
          GOTO 220
        ENDIF
C         PRINT*,GHZ,ZZ(1)
        IF (GHZ .GE. ZZ(1)) THEN
          DO 100 JG=1,J-1
            IF (ZZ(JG) .LE. GHZ .AND. ZZ(JG+1) .GT. GHZ) GOTO 120
  100     CONTINUE
C          IF (JG .GE. J) JG = J-1
  120     TNG = TN(JG)+(TN(JG+1)-TN(JG))*(GHZ-ZZ(JG))/(ZZ(JG+1)-ZZ(JG))
          DO 140 I=1,NMAJ
          ZCG(I) = ZVCD(I,JG) * (ZVCD(I,JG+1) / ZVCD(I,JG)) **
     >                      ((GHZ-ZZ(JG)) / (ZZ(JG+1)-ZZ(JG)))
  140     CONTINUE
        ELSE
          DO 160 JG=1,3
          IF (ZUS(JG) .LT. GHZ .AND. ZUS(JG+1) .GT. GHZ) GOTO 180
  160     CONTINUE
C          IF (JG .GE. 4) JG = 3
  180     TNG = TNUS(JG)
     >        + (TNUS(JG+1)-TNUS(JG))*(GHZ-ZUS(JG))/(ZUS(JG+1)-ZUS(JG))
C          PRINT*,JG,GHZ,ZUS(JG),ZUS(JG+1)
          DO 200 I=1,NMAJ
          ZCG(I) = ZCUS(I,JG) * (ZCUS(I,JG+1) / ZCUS(I,JG)) **
     >                      ((GHZ-ZUS(JG)) / (ZUS(JG+1)-ZUS(JG)))
  200     CONTINUE
        ENDIF
        DO 210 I=1,NMAJ
        ZCOL(I,J) = 2. * ZCG(I) * CHAP(PI/2.,GHZ,TNG,I)
     >              - ZVCD(I,J) * CHAP(CHI,ZZ(J),TN(J),I)
  210   CONTINUE
  220   CONTINUE
      ENDIF
C
      RETURN 
      END 
