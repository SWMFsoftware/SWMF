c
c##########################################################################
c#                                                                        #
c#                              GEOPACK                                   #
c#                     (MAIN SET OF FORTRAN CODES)                        #
c#                                                                        #
c##########################################################################
C

c
c This collection of subroutines is a result of several upgrades of the original
c      package by N. A. Tsyganenko. This release is dated by January 5, 2001.
c
c  Please see the additional file geopack.txt for descriptions of the individual
c    subroutines.
c
c    The following changes were made to the previous release of April 16, 1996:
c
c     (1) The latest set of IGRF expansion coefficients for the epoch 2000.0
c    was added, so that the updated version of the main field model accepts
c    dates through 2005.
c
c     (2) All explicit type statements as well as IMPLICIT NONE were taken
c    out from all subroutines, for the sake of brevity.

c     (3) To avoid potential run-time errors on some types of FORTRAN compilers,
c   the SAVE statements were added in the subroutines IGRF, DIP, RECALC, GEOMAG.
C
C     (4) In the tracing subroutine TRACE, more harmonic terms were allowed
c    to be included in the main field expansion (IGRF) for maintaining a higher
c    accuracy.
c
c     (5) An update and some correction/editing of the documentation file
c    GEOPACK.TXT was made.
c
c    Mei-Ching Fok noted (January 6, 2002):
c    subroutines dipole, locatep and trace1 was added. trace1 is a modification 
c    if TRACE. trace1 does field line tracing in sm coordinates.  Also tracing 
c    stops at magnetic equator. 
c
c-------------------------------------------------------------------------------
c
      SUBROUTINE IGRF(IY,NM,R,T,F,BR,BT,BF)
c
C     CALCULATES COMPONENTS OF THE MAIN (INTERNAL) GEOMAGNETIC FIELD IN SPHERICAL
C     GEOGRAPHICAL COORDINATE SYSTEM, USING IAGA INTERNATIONAL GEOMAGNETIC REFERENCE
C     MODEL COEFFICIENTS (e.g., http://www.ngdc.noaa.gov/IAGA/wg8/igrf2000.html)
C
C     UPDATING THE COEFFICIENTS TO A GIVEN EPOCH IS MADE AUTOMATICALLY UPON THE FIRST
C     CALL AND AFTER EVERY CHANGE OF THE PARAMETER IY.
C
C-----INPUT PARAMETERS:
C
C     IY  -  YEAR NUMBER (FOUR-DIGIT; 1965 &LE IY &LE 2005)
C     NM  -  HIGHEST ORDER OF SPHERICAL HARMONICS IN THE SCALAR POTENTIAL (NM &LE 10)
C     R,T,F -  SPHERICAL COORDINATES (RADIUS R IN UNITS RE=6371.2 KM, GEOGRAPHIC
C                COLATITUDE  T  AND LONGITUDE  F  IN RADIANS)
C
C-----OUTPUT PARAMETERS:
C
C     BR,BT,BF - SPHERICAL COMPONENTS OF THE MAIN GEOMAGNETIC FIELD IN NANOTESLA
C
C     LAST MODIFICATION:  JANUARY 5, 2001.
C     THE CODE WAS MODIFIED TO ACCEPT DATES THROUGH 2005.
C     IT HAS ALSO BEEN SLIGHTLY SIMPLIFIED BY TAKING OUT SOME REDUNDANT STATEMENTS,
C     AND A "SAVE" STATEMENT WAS ADDED, TO AVOID POTENTIAL PROBLEMS WITH SOME
C     FORTRAN COMPILERS.
c
C     WRITTEN BY: N. A. TSYGANENKO
C
      SAVE MA,IYR,IPR
C
      DIMENSION A(11),B(11),G(66),H(66),REC(66),G65(66),H65(66),G70(66),
     * H70(66),G75(66),H75(66),G80(66),H80(66),G85(66),H85(66),G90(66),
     * H90(66),G95(66),H95(66),G00(66),H00(66),DG00(45),DH00(45)
c
      DATA G65/0.,-30334.,-2119.,-1662.,2997.,1594.,1297.,-2038.,1292.,
     *856.,957.,804.,479.,-390.,252.,-219.,358.,254.,-31.,-157.,-62.,
     *45.,61.,8.,-228.,4.,1.,-111.,75.,-57.,4.,13.,-26.,-6.,13.,1.,13.,
     *5.,-4.,-14.,0.,8.,-1.,11.,4.,8.,10.,2.,-13.,10.,-1.,-1.,5.,1.,-2.,
     *-2.,-3.,2.,-5.,-2.,4.,4.,0.,2.,2.,0./
      DATA H65/0.,0.,5776.,0.,-2016.,114.,0.,-404.,240.,-165.,0.,148.,
     *-269.,13.,-269.,0.,19.,128.,-126.,-97.,81.,0.,-11.,100.,68.,-32.,
     *-8.,-7.,0.,-61.,-27.,-2.,6.,26.,-23.,-12.,0.,7.,-12.,9.,-16.,4.,
     *24.,-3.,-17.,0.,-22.,15.,7.,-4.,-5.,10.,10.,-4.,1.,0.,2.,1.,2.,
     *6.,-4.,0.,-2.,3.,0.,-6./
c
      DATA G70/0.,-30220.,-2068.,-1781.,3000.,1611.,1287.,-2091.,1278.,
     *838.,952.,800.,461.,-395.,234.,-216.,359.,262.,-42.,-160.,-56.,
     *43.,64.,15.,-212.,2.,3.,-112.,72.,-57.,1.,14.,-22.,-2.,13.,-2.,
     *14.,6.,-2.,-13.,-3.,5.,0.,11.,3.,8.,10.,2.,-12.,10.,-1.,0.,3.,
     *1.,-1.,-3.,-3.,2.,-5.,-1.,6.,4.,1.,0.,3.,-1./
      DATA H70/0.,0.,5737.,0.,-2047.,25.,0.,-366.,251.,-196.,0.,167.,
     *-266.,26.,-279.,0.,26.,139.,-139.,-91.,83.,0.,-12.,100.,72.,-37.,
     *-6.,1.,0.,-70.,-27.,-4.,8.,23.,-23.,-11.,0.,7.,-15.,6.,-17.,6.,
     *21.,-6.,-16.,0.,-21.,16.,6.,-4.,-5.,10.,11.,-2.,1.,0.,1.,1.,3.,
     *4.,-4.,0.,-1.,3.,1.,-4./
c
      DATA G75/0.,-30100.,-2013.,-1902.,3010.,1632.,1276.,-2144.,1260.,
     *830.,946.,791.,438.,-405.,216.,-218.,356.,264.,-59.,-159.,-49.,
     *45.,66.,28.,-198.,1.,6.,-111.,71.,-56.,1.,16.,-14.,0.,12.,-5.,
     *14.,6.,-1.,-12.,-8.,4.,0.,10.,1.,7.,10.,2.,-12.,10.,-1.,-1.,4.,
     *1.,-2.,-3.,-3.,2.,-5.,-2.,5.,4.,1.,0.,3.,-1./
      DATA H75/0.,0.,5675.,0.,-2067.,-68.,0.,-333.,262.,-223.,0.,191.,
     *-265.,39.,-288.,0.,31.,148.,-152.,-83.,88.,0.,-13.,99.,75.,-41.,
     *-4.,11.,0.,-77.,-26.,-5.,10.,22.,-23.,-12.,0.,6.,-16.,4.,-19.,6.,
     *18.,-10.,-17.,0.,-21.,16.,7.,-4.,-5.,10.,11.,-3.,1.,0.,1.,1.,3.,
     *4.,-4.,-1.,-1.,3.,1.,-5./
c
      DATA G80/0.,-29992.,-1956.,-1997.,3027.,1663.,1281.,-2180.,1251.,
     *833.,938.,782.,398.,-419.,199.,-218.,357.,261.,-74.,-162.,-48.,
     *48.,66.,42.,-192.,4.,14.,-108.,72.,-59.,2.,21.,-12.,1.,11.,-2.,
     *18.,6.,0.,-11.,-7.,4.,3.,6.,-1.,5.,10.,1.,-12.,9.,-3.,-1.,7.,2.,
     *-5.,-4.,-4.,2.,-5.,-2.,5.,3.,1.,2.,3.,0./
      DATA H80/0.,0.,5604.,0.,-2129.,-200.,0.,-336.,271.,-252.,0.,212.,
     *-257.,53.,-297.,0.,46.,150.,-151.,-78.,92.,0.,-15.,93.,71.,-43.,
     *-2.,17.,0.,-82.,-27.,-5.,16.,18.,-23.,-10.,0.,7.,-18.,4.,-22.,9.,
     *16.,-13.,-15.,0.,-21.,16.,9.,-5.,-6.,9.,10.,-6.,2.,0.,1.,0.,3.,
     *6.,-4.,0.,-1.,4.,0.,-6./
c
      DATA G85/0.,-29873.,-1905.,-2072.,3044.,1687.,1296.,-2208.,1247.,
     *829.,936.,780.,361.,-424.,170.,-214.,355.,253.,-93.,-164.,-46.,
     *53.,65.,51.,-185.,4.,16.,-102.,74.,-62.,3.,24.,-6.,4.,10.,0.,21.,
     *6.,0.,-11.,-9.,4.,4.,4.,-4.,5.,10.,1.,-12.,9.,-3.,-1.,7.,1.,-5.,
     *-4.,-4.,3.,-5.,-2.,5.,3.,1.,2.,3.,0./
      DATA H85/0.,0.,5500.,0.,-2197.,-306.,0.,-310.,284.,-297.,0.,232.,
     *-249.,69.,-297.,0.,47.,150.,-154.,-75.,95.,0.,-16.,88.,69.,-48.,
     *-1.,21.,0.,-83.,-27.,-2.,20.,17.,-23.,-7.,0.,8.,-19.,5.,-23.,11.,
     *14.,-15.,-11.,0.,-21.,15.,9.,-6.,-6.,9.,9.,-7.,2.,0.,1.,0.,3.,
     *6.,-4.,0.,-1.,4.,0.,-6./
c
      DATA G90/0., -29775.,  -1848.,  -2131.,   3059.,   1686.,   1314.,
     *     -2239.,   1248.,    802.,    939.,    780.,    325.,   -423.,
     *       141.,   -214.,    353.,    245.,   -109.,   -165.,    -36.,
     *        61.,     65.,     59.,   -178.,      3.,     18.,    -96.,
     *        77.,    -64.,      2.,     26.,     -1.,      5.,      9.,
     *         0.,     23.,      5.,     -1.,    -10.,    -12.,      3.,
     *         4.,      2.,     -6.,      4.,      9.,      1.,    -12.,
     *         9.,     -4.,     -2.,      7.,      1.,     -6.,     -3.,
     *        -4.,      2.,     -5.,     -2.,      4.,      3.,      1.,
     *         3.,      3.,      0./

      DATA H90/0.,      0.,   5406.,      0.,  -2279.,   -373.,      0.,
     *      -284.,    293.,   -352.,      0.,    247.,   -240.,     84.,
     *      -299.,      0.,     46.,    154.,   -153.,    -69.,     97.,
     *         0.,    -16.,     82.,     69.,    -52.,      1.,     24.,
     *         0.,    -80.,    -26.,      0.,     21.,     17.,    -23.,
     *        -4.,      0.,     10.,    -19.,      6.,    -22.,     12.,
     *        12.,    -16.,    -10.,      0.,    -20.,     15.,     11.,
     *        -7.,     -7.,      9.,      8.,     -7.,      2.,      0.,
     *         2.,      1.,      3.,      6.,     -4.,      0.,     -2.,
     *         3.,     -1.,     -6./

      DATA G95/0., -29682.,  -1789.,  -2197.,   3074.,   1685.,   1329.,
     *     -2268.,   1249.,    769.,    941.,    782.,    291.,   -421.,
     *       116.,   -210.,    352.,    237.,   -122.,   -167.,    -26.,
     *        66.,     64.,     65.,   -172.,      2.,     17.,    -94.,
     *        78.,    -67.,      1.,     29.,      4.,      8.,     10.,
     *        -2.,     24.,      4.,     -1.,     -9.,    -14.,      4.,
     *         5.,      0.,     -7.,      4.,      9.,      1.,    -12.,
     *         9.,     -4.,     -2.,      7.,      0.,     -6.,     -3.,
     *        -4.,      2.,     -5.,     -2.,      4.,      3.,      1.,
     *         3.,      3.,      0./

      DATA H95/0.,      0.,   5318.,      0.,  -2356.,   -425.,      0.,
     *      -263.,    302.,   -406.,      0.,    262.,   -232.,     98.,
     *      -301.,      0.,     44.,    157.,   -152.,    -64.,     99.,
     *         0.,    -16.,     77.,     67.,    -57.,      4.,     28.,
     *         0.,    -77.,    -25.,      3.,     22.,     16.,    -23.,
     *        -3.,      0.,     12.,    -20.,      7.,    -21.,     12.,
     *        10.,    -17.,    -10.,      0.,    -19.,     15.,     11.,
     *        -7.,     -7.,      9.,      7.,     -8.,      1.,      0.,
     *         2.,      1.,      3.,      6.,     -4.,      0.,     -2.,
     *         3.,     -1.,     -6./

      DATA G00/0., -29615.,  -1728.,  -2267.,   3072.,   1672.,   1341.,
     *     -2290.,   1253.,    715.,    935.,    787.,    251.,   -405.,
     *       110.,   -217.,    351.,    222.,   -131.,   -169.,    -12.,
     *        72.,     68.,     74.,   -161.,     -5.,     17.,    -91.,
     *        79.,    -74.,      0.,     33.,      9.,      7.,      8.,
     *        -2.,     25.,      6.,     -9.,     -8.,    -17.,      9.,
     *         7.,     -8.,     -7.,      5.,      9.,      3.,    - 8.,
     *         6.,     -9.,     -2.,      9.,     -4.,     -8.,     -2.,
     *        -6.,      2.,     -3.,      0.,      4.,      1.,      2.,
     *         4.,      0.,     -1./

      DATA H00/0.,      0.,   5186.,      0.,  -2478.,   -458.,      0.,
     *      -227.,    296.,   -492.,      0.,    272.,   -232.,    119.,
     *      -304.,      0.,     44.,    172.,   -134.,    -40.,    107.,
     *         0.,    -17.,     64.,     65.,    -61.,      1.,     44.,
     *         0.,    -65.,    -24.,      6.,     24.,     15.,    -25.,
     *        -6.,      0.,     12.,    -22.,      8.,    -21.,     15.,
     *         9.,    -16.,     -3.,      0.,    -20.,     13.,     12.,
     *        -6.,     -8.,      9.,      4.,     -8.,      5.,      0.,
     *         1.,      0.,      4.,      5.,     -6.,     -1.,     -3.,
     *         0.,     -2.,     -8./

      DATA DG00/0.0,  14.6,    10.7,   -12.4,     1.1,    -1.1,     0.7,
     *         -5.4,   0.9,    -7.7,    -1.3,     1.6,    -7.3,     2.9,
     *         -3.2,   0.0,    -0.7,    -2.1,    -2.8,    -0.8,     2.5,
     *          1.0,  -0.4,     0.9,     2.0,    -0.6,    -0.3,     1.2,
     *         -0.4,  -0.4,    -0.3,     1.1,     1.1,    -0.2,     0.6,
     *         -0.9,  -0.3,     0.2,    -0.3,     0.4,    -1.0,     0.3,
     *         -0.5,  -0.7,    -0.4/

      DATA DH00/0.0,   0.0,   -22.5,     0.0,   -20.6,    -9.6,     0.0,
     *          6.0,  -0.1,   -14.2,     0.0,     2.1,     1.3,     5.0,
     *          0.3,   0.0,    -0.1,     0.6,     1.7,     1.9,     0.1,
     *          0.0,  -0.2,    -1.4,     0.0,    -0.8,     0.0,     0.9,
     *          0.0,   1.1,     0.0,     0.3,    -0.1,    -0.6,    -0.7,
     *          0.2,   0.0,     0.1,     0.0,     0.0,     0.3,     0.6,
     *         -0.4,   0.3,     0.7/
c
c
      DATA MA,IYR,IPR/0,0,0/

      IF(MA.NE.1) GOTO 10
      IF(IY.NE.IYR) GOTO 30
      GOTO 130
10    MA=1
C
      DO 20 N=1,11
         N2=2*N-1
         N2=N2*(N2-2)
         DO 20 M=1,N
            MN=N*(N-1)/2+M
20    REC(MN)=FLOAT((N-M)*(N+M-2))/FLOAT(N2)
C
30    IYR=IY
      IF (IYR.LT.1965) IYR=1965
      IF (IYR.GT.2005) IYR=2005
      IF (IY.NE.IYR.AND.IPR.EQ.0) WRITE (*,999)IY,IYR
      IF (IYR.NE.IY) IPR=1
      IF (IYR.LT.1970) GOTO 50          !INTERPOLATE BETWEEN 1965 - 1970
      IF (IYR.LT.1975) GOTO 60          !INTERPOLATE BETWEEN 1970 - 1975
      IF (IYR.LT.1980) GOTO 70          !INTERPOLATE BETWEEN 1975 - 1980
      IF (IYR.LT.1985) GOTO 80          !INTERPOLATE BETWEEN 1980 - 1985
      IF (IYR.LT.1990) GOTO 90          !INTERPOLATE BETWEEN 1985 - 1990
      IF (IYR.LT.1995) GOTO 100         !INTERPOLATE BETWEEN 1990 - 1995
      IF (IYR.LT.2000) GOTO 110         !INTERPOLATE BETWEEN 1995 - 2000
C
C       EXTRAPOLATE BEYOND 2000:
C
      DT=FLOAT(IYR)-2000.
      DO 40 N=1,66
         G(N)=G00(N)
         H(N)=H00(N)
         IF (N.GT.45) GOTO 40
         G(N)=G(N)+DG00(N)*DT
         H(N)=H(N)+DH00(N)*DT
40    CONTINUE
      GOTO 300
C
C       INTERPOLATE BETWEEEN 1965 - 1970:
C
50    F2=(IYR-1965)/5.
      F1=1.-F2
      DO 55 N=1,66
         G(N)=G65(N)*F1+G70(N)*F2
55       H(N)=H65(N)*F1+H70(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 1970 - 1975:
C
60    F2=(IYR-1970)/5.
      F1=1.-F2
      DO 65 N=1,66
         G(N)=G70(N)*F1+G75(N)*F2
65       H(N)=H70(N)*F1+H75(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 1975 - 1980:
C
70    F2=(IYR-1975)/5.
      F1=1.-F2
      DO 75 N=1,66
         G(N)=G75(N)*F1+G80(N)*F2
75       H(N)=H75(N)*F1+H80(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 1980 - 1985:
C
80    F2=(IYR-1980)/5.
      F1=1.-F2
      DO 85 N=1,66
         G(N)=G80(N)*F1+G85(N)*F2
85       H(N)=H80(N)*F1+H85(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 1985 - 1990:
C
90    F2=(IYR-1985)/5.
      F1=1.-F2
      DO 95 N=1,66
         G(N)=G85(N)*F1+G90(N)*F2
95       H(N)=H85(N)*F1+H90(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 1990 - 1995:
C
100   F2=(IYR-1990)/5.
      F1=1.-F2
      DO 105 N=1,66
         G(N)=G90(N)*F1+G95(N)*F2
105      H(N)=H90(N)*F1+H95(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 1995 - 2000:
C
110   F2=(IYR-1995)/5.
      F1=1.-F2
      DO 115 N=1,66
         G(N)=G95(N)*F1+G00(N)*F2
115      H(N)=H95(N)*F1+H00(N)*F2
      GOTO 300
C
C   COEFFICIENTS FOR A GIVEN YEAR HAVE BEEN CALCULATED; NOW MULTIPLY
C   THEM BY SCHMIDT NORMALIZATION FACTORS:
C
300   S=1.
      DO 120 N=2,11
         MN=N*(N-1)/2+1
         S=S*FLOAT(2*N-3)/FLOAT(N-1)
         G(MN)=G(MN)*S
         H(MN)=H(MN)*S
         P=S
         DO 120 M=2,N
            AA=1.
            IF (M.EQ.2) AA=2.
            P=P*SQRT(AA*FLOAT(N-M+1)/FLOAT(N+M-2))
            MNN=MN+M-1
            G(MNN)=G(MNN)*P
120         H(MNN)=H(MNN)*P
C
C     NOW CALCULATE THE FIELD COMPONENTS
C     (IN CASE OF MULTIPLE INVOCATIONS WITH THE SAME VALUES OF IY AND NM,
C      CALCULATIONS START RIGHT HERE):
C
130   PP=1./R
      P=PP

      K=NM+1
      DO 150 N=1,K
         P=P*PP
         A(N)=P
150      B(N)=P*N
      P=1.
      D=0.
      BBR=0.
      BBT=0.
      BBF=0.
      U=T
      CF=COS(F)
      SF=SIN(F)
      C=COS(U)
      S=SIN(U)
      DO 200 M=1,K
         IF(M.EQ.1) GOTO 160
         MM=M-1
         W=X
         X=W*CF+Y*SF
         Y=Y*CF-W*SF
         GOTO 170
160      X=0.
         Y=1.
170      Q=P
         Z=D
         BI=0.
         P2=0.
         D2=0.
         DO 190 N=M,K
            AN=A(N)
            MN=N*(N-1)/2+M
            E=G(MN)
            HH=H(MN)
            W=E*Y+HH*X
            BBR=BBR+B(N)*W*Q
            BBT=BBT-AN*W*Z
            IF(M.EQ.1) GOTO 180
            QQ=Q
            IF(S.LT.1.E-5) QQ=Z
            BI=BI+AN*(E*X-HH*Y)*QQ
180         XK=REC(MN)
            DP=C*Z-S*Q-XK*D2
            PM=C*Q-XK*P2
            D2=Z
            P2=Q
            Z=DP
190        Q=PM
         D=S*D+C*P
         P=S*P
         IF(M.EQ.1) GOTO 200
         BI=BI*MM
         BBF=BBF+BI
200   CONTINUE
C
      BR=BBR
      BT=BBT
      IF(S.LT.1.E-5) GOTO 210
      BF=BBF/S
      RETURN
210   IF(C.LT.0.) BBF=-BBF
      BF=BBF

      RETURN
C
999   FORMAT(//1X,
     * 'IGRF WARNS:**** YEAR IS OUT OF INTERVAL 1965-2005: IY =',I5,/,
     *',        CALCULATIONS WILL BE DONE FOR IYR =',I5,' ****'//)
      END
C

c
       SUBROUTINE DIP(PS,X,Y,Z,BX,BY,BZ)
C
C
C  CALCULATES GSM COMPONENTS OF A GEODIPOLE FIELD WITH THE DIPOLE MOMENT
C  CORRESPONDING TO THE EPOCH OF 2000.
C
C----INPUT PARAMETERS:
C     PS - GEODIPOLE TILT ANGLE IN RADIANS,
C     X,Y,Z - GSM COORDINATES IN RE (1 RE = 6371.2 km)
C
C----OUTPUT PARAMETERS:
C     BX,BY,BZ - FIELD COMPONENTS IN GSM SYSTEM, IN NANOTESLA.
C
C  LAST MODIFICATION: JAN. 5, 2001. THE VALUE OF THE DIPOLE MOMENT WAS UPDATED TO 2000.
C    AND A "SAVE" STATEMENT HAS BEEN ADDED, TO AVOID POTENTIAL PROBLEMS WITH SOME
C    FORTRAN COMPILERS
C
C  WRITTEN BY: N. A. TSYGANENKO
C
      SAVE M,PSI
      DATA M,PSI/0,5./
      IF(M.EQ.1.AND.ABS(PS-PSI).LT.1.E-5) GOTO 1
      SPS=SIN(PS)
      CPS=COS(PS)
      PSI=PS
      M=1
  1   P=X**2
      U=Z**2
      V=3.*Z*X
      T=Y**2
      Q=30115./SQRT(P+T+U)**5
      BX=Q*((T+U-2.*P)*SPS-V*CPS)
      BY=-3.*Y*Q*(X*SPS+Z*CPS)
      BZ=Q*((P+T-2.*U)*CPS-V*SPS)
      RETURN
      END

c
      SUBROUTINE SUN(IYR,IDAY,IHOUR,MIN,ISEC,GST,SLONG,SRASN,SDEC)
C
C  CALCULATES FOUR QUANTITIES NECESSARY FOR COORDINATE TRANSFORMATIONS
C  WHICH DEPEND ON SUN POSITION (AND, HENCE, ON UNIVERSAL TIME AND SEASON)
C
C-------  INPUT PARAMETERS:
C  IYR,IDAY,IHOUR,MIN,ISEC -  YEAR, DAY, AND UNIVERSAL TIME IN HOURS, MINUTES,
C    AND SECONDS  (IDAY=1 CORRESPONDS TO JANUARY 1).
C
C-------  OUTPUT PARAMETERS:
C  GST - GREENWICH MEAN SIDEREAL TIME, SLONG - LONGITUDE ALONG ECLIPTIC
C  SRASN - RIGHT ASCENSION,  SDEC - DECLINATION  OF THE SUN (RADIANS)
C  ORIGINAL VERSION OF THIS SUBROUTINE HAS BEEN COMPILED FROM:
C  RUSSELL, C.T., COSMIC ELECTRODYNAMICS, 1971, V.2, PP.184-196.
C
C     LAST MODIFICATION:  JAN 5, 2001 (NO ESSENTIAL CHANGES, BUT
C     SOME REDUNDANT STATEMENTS TAKEN OUT FROM THE PREVIOUS VERSION)
C
C     ORIGINAL VERSION WRITTEN BY:    Gilbert D. Mead
C
      DOUBLE PRECISION DJ,FDAY
      DATA RAD/57.295779513/
C
      IF(IYR.LT.1901.OR.IYR.GT.2099) RETURN
      FDAY=REAL(IHOUR*3600+MIN*60+ISEC)/86400.D0
      DJ=365*(IYR-1900)+(IYR-1901)/4+IDAY-0.5D0+FDAY
      T=DJ/36525.
      VL=DMOD(279.696678+0.9856473354*DJ,360.D0)
      GST=DMOD(279.690983+.9856473354*DJ+360.*FDAY+180.,360.D0)/RAD
      G=DMOD(358.475845+0.985600267*DJ,360.D0)/RAD
      SLONG=(VL+(1.91946-0.004789*T)*SIN(G)+0.020094*SIN(2.*G))/RAD
      IF(SLONG.GT.6.2831853) SLONG=SLONG-6.2831853
      IF (SLONG.LT.0.) SLONG=SLONG+6.2831853
      OBLIQ=(23.45229-0.0130125*T)/RAD
      SOB=SIN(OBLIQ)
      SLP=SLONG-9.924E-5
C
C   THE LAST CONSTANT IS A CORRECTION FOR THE ANGULAR ABERRATION  DUE TO
C   THE ORBITAL MOTION OF THE EARTH
C
      SIND=SOB*SIN(SLP)
      COSD=SQRT(1.-SIND**2)
      SC=SIND/COSD
      SDEC=ATAN(SC)
      SRASN=3.141592654-ATAN2(COS(OBLIQ)/SOB*SC,-COS(SLP)/COSD)
      RETURN
      END

c
      SUBROUTINE SPHCAR(R,TETA,PHI,X,Y,Z,J)
C
C   CONVERTS SPHERICAL COORDS INTO CARTESIAN ONES AND VICA VERSA
C    (TETA AND PHI IN RADIANS).
C
C                  J &Gt 0            J &Lt 0
C-----INPUT:   J,R,TETA,PHI     J,X,Y,Z
C----OUTPUT:      X,Y,Z        R,TETA,PHI
C
C   LAST MOFIFICATION:  JAN 5, 2001 (NO ESSENTIAL CHANGES, BUT
C                        SOME REDUNDANT STATEMENTS TAKEN OUT)
C
C   WRITTEN BY:  N. A. TSYGANENKO
C
      IF(J.GT.0) GOTO 3
      SQ=X**2+Y**2
      R=SQRT(SQ+Z**2)
      IF (SQ.NE.0.) GOTO 2
      PHI=0.
      IF (Z.LT.0.) GOTO 1
      TETA=0.
      RETURN
  1   TETA=3.141592654
      RETURN
  2   SQ=SQRT(SQ)
      PHI=ATAN2(Y,X)
      TETA=ATAN2(SQ,Z)
      IF (PHI.LT.0.) PHI=PHI+6.28318531
      RETURN
  3   SQ=R*SIN(TETA)
      X=SQ*COS(PHI)
      Y=SQ*SIN(PHI)
      Z=R*COS(TETA)
      RETURN
      END

c
      SUBROUTINE BSPCAR(TETA,PHI,BR,BTET,BPHI,BX,BY,BZ)
C
C   CALCULATES CARTESIAN FIELD COMPONENTS FROM SPHERICAL ONES
C-----INPUT:   TETA,PHI - SPHERICAL ANGLES OF THE POINT IN RADIANS
C              BR,BTET,BPHI -  SPHERICAL COMPONENTS OF THE FIELD
C-----OUTPUT:  BX,BY,BZ - CARTESIAN COMPONENTS OF THE FIELD
C
C   LAST MOFIFICATION:  JAN 5 2001 (NO ESSENTIAL CHANGES, BUT
C                        SOME REDUNDANT STATEMENTS TAKEN OUT)
C
C   WRITTEN BY:  N. A. TSYGANENKO
C
      S=SIN(TETA)
      C=COS(TETA)
      SF=SIN(PHI)
      CF=COS(PHI)
      BE=BR*S+BTET*C
      BX=BE*CF-BPHI*SF
      BY=BE*SF+BPHI*CF
      BZ=BR*C-BTET*S
      RETURN
      END
c

C
      SUBROUTINE RECALC(IYR,IDAY,IHOUR,MIN,ISEC)
C
C  PREPARES ELEMENTS OF ROTATION MATRICES FOR TRANSFORMATIONS OF VECTORS BETWEEN
C  SEVERAL COORDINATE SYSTEMS, MOST FREQUENTLY USED IN SPACE PHYSICS.
C
C  THIS SUBROUTINE SHOULD BE INVOKED BEFORE USING THE FOLLOWING SUBROUTINES:
C         GEOGSM, MAGSM, SMGSM, GSMGSE, GEIGEO.
C
C  THERE IS NO NEED TO REPEATEDLY INVOKE RECALC, IF MULTIPLE CALCULATIONS ARE MADE
C    FOR THE SAME DATE AND TIME.
C
C-----INPUT PARAMETERS:
C
C     IYR   -  YEAR NUMBER (FOUR DIGITS)
C     IDAY  -  DAY OF YEAR (DAY 1 = JAN 1)
C     IHOUR -  HOUR OF DAY (00 TO 23)
C     MIN   -  MINUTE OF HOUR (00 TO 59)
C     ISEC  -  SECONDS OF MINUTE (00 TO 59)
C
C-----OUTPUT PARAMETERS:   NONE (ALL OUTPUT QUANTITIES ARE PLACED
C                                 INTO THE COMMON BLOCK /GEOPACK/)
C
C    OTHER SUBROUTINES CALLED BY THIS ONE: SUN
C
C   ################################################
C   #  WRITTEN BY  N.A. TSYGANENKO ON DEC.1, 1991  #
C   ################################################
c
c    Last modification:  January 5, 2001.
c    The code has been modified to accept dates through 2005.
c    Also, a "save" statement was added, to avoid potential problems
c    with some Fortran compilers.
C
      COMMON /GEOPACK/ ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI,SPS,
     * CPS,SHI,CHI,HI,PSI,XMUT,A11,A21,A31,A12,A22,A32,A13,A23,A33,DS3,
     * K,IY,CGST,SGST,BA(6)
c
      SAVE IYE,IDE,IPR
      DATA IYE,IDE,IPR/3*0/
C
      IF (IYR.EQ.IYE.AND.IDAY.EQ.IDE) GOTO 5
C
C   IYE AND IDE ARE THE CURRENT VALUES OF YEAR AND DAY NUMBER
C
      IY=IYR
      IDE=IDAY
      IF(IY.LT.1965) IY=1965
      IF(IY.GT.2005) IY=2005
C
C  WE ARE RESTRICTED BY THE INTERVAL 1965-2005,
C  FOR WHICH THE IGRF COEFFICIENTS ARE KNOWN; IF IYR IS OUTSIDE THIS INTERVAL
C  THE SUBROUTINE PRINTS A WARNING (BUT DOES NOT REPEAT IT AT NEXT INVOCATIONS)
C
      IF(IY.NE.IYR.AND.IPR.EQ.0) WRITE (*,10) IYR,IY
      IF(IY.NE.IYR) IPR=1
      IYE=IY
C
C  LINEAR INTERPOLATION OF THE GEODIPOLE MOMENT COMPONENTS BETWEEN THE
C  VALUES FOR THE NEAREST EPOCHS:
C
        IF (IY.LT.1970) THEN                            !1965-1970
           F2=(FLOAT(IY)+FLOAT(IDAY)/365.-1965.)/5.
           F1=1.D0-F2
           G10=30334.*F1+30220.*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
           G11=-2119.*F1-2068.*F2
           H11=5776.*F1+5737.*F2
        ELSEIF (IY.LT.1975) THEN                        !1970-1975
           F2=(FLOAT(IY)+FLOAT(IDAY)/365.-1970.)/5.
           F1=1.D0-F2
           G10=30220.*F1+30100.*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
           G11=-2068.*F1-2013.*F2
           H11=5737.*F1+5675.*F2
        ELSEIF (IY.LT.1980) THEN                        !1975-1980
           F2=(REAL(IY)+REAL(IDAY)/365.-1975.)/5.
           F1=1.D0-F2
           G10=30100.*F1+29992.*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
           G11=-2013.*F1-1956.*F2
           H11=5675.*F1+5604.*F2
        ELSEIF (IY.LT.1985) THEN                        !1980-1985
           F2=(FLOAT(IY)+FLOAT(IDAY)/365.-1980.)/5.
           F1=1.D0-F2
           G10=29992.*F1+29873.*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
           G11=-1956.*F1-1905.*F2
           H11=5604.*F1+5500.*F2
        ELSEIF (IY.LT.1990) THEN                        !1985-1990
           F2=(FLOAT(IY)+FLOAT(IDAY)/365.-1985.)/5.
           F1=1.D0-F2
           G10=29873.*F1+29775.*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
           G11=-1905.*F1-1848.*F2
           H11=5500.*F1+5406.*F2
        ELSEIF (IY.LT.1995) THEN                        !1990-1995
           F2=(FLOAT(IY)+FLOAT(IDAY)/365.-1990.)/5.
           F1=1.D0-F2
           G10=29775.*F1+29682.*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
           G11=-1848.*F1-1789.*F2
           H11=5406.*F1+5318.*F2
        ELSEIF (IY.LT.2000) THEN                        !1995-2000
c          F2=(FLOAT(IY)+FLOAT(IDAY)/365.-1990.)/5.
           F2=(FLOAT(IY)+FLOAT(IDAY)/365.-1995.)/5.     ! Kolya April 3, 2001
           F1=1.D0-F2
           G10=29682.*F1+29615.*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
           G11=-1789.*F1-1728.*F2
           H11=5318.*F1+5186.*F2
        ELSE                                            !2000-2005
C
C   LINEAR EXTRAPOLATION BEYOND 2000 BY USING SECULAR VELOCITY COEFFICIENTS:
C
           DT=FLOAT(IY)+FLOAT(IDAY)/366.-2000.
           G10=29615.-14.6*DT      ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
           G11=-1728.+10.7*DT
           H11=5186.-22.5*DT
        ENDIF
C
C  NOW CALCULATE THE COMPONENTS OF THE UNIT VECTOR EzMAG IN GEO COORD.SYSTEM:
C   SIN(TETA0)*COS(LAMBDA0), SIN(TETA0)*SIN(LAMBDA0), AND COS(TETA0)
C         ST0 * CL0                ST0 * SL0                CT0
C
      SQ=G11**2+H11**2
      SQQ=SQRT(SQ)
      SQR=SQRT(G10**2+SQ)
      SL0=-H11/SQQ
      CL0=-G11/SQQ
      ST0=SQQ/SQR
      CT0=G10/SQR
      STCL=ST0*CL0
      STSL=ST0*SL0
      CTSL=CT0*SL0
      CTCL=CT0*CL0
C
C      THE CALCULATIONS ARE TERMINATED IF ONLY GEO-MAG TRANSFORMATION
C       IS TO BE DONE  (IHOUR>24 IS THE AGREED INDICATOR FOR THIS CASE):
C
   5   IF (IHOUR.GT.24) RETURN
C
      CALL SUN(IY,IDAY,IHOUR,MIN,ISEC,GST,SLONG,SRASN,SDEC)
C
C  S1,S2, AND S3 ARE THE COMPONENTS OF THE UNIT VECTOR EXGSM=EXGSE IN THE
C   SYSTEM GEI POINTING FROM THE EARTH'S CENTER TO THE SUN:
C
      S1=COS(SRASN)*COS(SDEC)
      S2=SIN(SRASN)*COS(SDEC)
      S3=SIN(SDEC)
      CGST=COS(GST)
      SGST=SIN(GST)
C
C  DIP1, DIP2, AND DIP3 ARE THE COMPONENTS OF THE UNIT VECTOR EZSM=EZMAG
C   IN THE SYSTEM GEI:
C
      DIP1=STCL*CGST-STSL*SGST
      DIP2=STCL*SGST+STSL*CGST
      DIP3=CT0
C
C  NOW CALCULATE THE COMPONENTS OF THE UNIT VECTOR EYGSM IN THE SYSTEM GEI
C   BY TAKING THE VECTOR PRODUCT D x S AND NORMALIZING IT TO UNIT LENGTH:
C
      Y1=DIP2*S3-DIP3*S2
      Y2=DIP3*S1-DIP1*S3
      Y3=DIP1*S2-DIP2*S1
      Y=SQRT(Y1*Y1+Y2*Y2+Y3*Y3)
      Y1=Y1/Y
      Y2=Y2/Y
      Y3=Y3/Y
C
C   THEN IN THE GEI SYSTEM THE UNIT VECTOR Z = EZGSM = EXGSM x EYGSM = S x Y
C    HAS THE COMPONENTS:
C
      Z1=S2*Y3-S3*Y2
      Z2=S3*Y1-S1*Y3
      Z3=S1*Y2-S2*Y1
C
C    THE VECTOR EZGSE (HERE DZ) IN GEI HAS THE COMPONENTS (0,-SIN(DELTA),
C     COS(DELTA)) = (0.,-0.397823,0.917462); HERE DELTA = 23.44214 DEG FOR
C   THE EPOCH 1978 (SEE THE BOOK BY GUREVICH OR OTHER ASTRONOMICAL HANDBOOKS).
C    HERE THE MOST ACCURATE TIME-DEPENDENT FORMULA IS USED:
C
      DJ=FLOAT(365*(IY-1900)+(IY-1901)/4 +IDAY)
     * -0.5+FLOAT(IHOUR*3600+MIN*60+ISEC)/86400.
      T=DJ/36525.
      OBLIQ=(23.45229-0.0130125*T)/57.2957795
      DZ1=0.
      DZ2=-SIN(OBLIQ)
      DZ3=COS(OBLIQ)
C
C  THEN THE UNIT VECTOR EYGSE IN GEI SYSTEM IS THE VECTOR PRODUCT DZ x S :
C
      DY1=DZ2*S3-DZ3*S2
      DY2=DZ3*S1-DZ1*S3
      DY3=DZ1*S2-DZ2*S1
C
C   THE ELEMENTS OF THE MATRIX GSE TO GSM ARE THE SCALAR PRODUCTS:
C   CHI=EM22=(EYGSM,EYGSE), SHI=EM23=(EYGSM,EZGSE), EM32=(EZGSM,EYGSE)=-EM23,
C     AND EM33=(EZGSM,EZGSE)=EM22
C
      CHI=Y1*DY1+Y2*DY2+Y3*DY3
      SHI=Y1*DZ1+Y2*DZ2+Y3*DZ3
      HI=ASIN(SHI)
C
C    TILT ANGLE: PSI=ARCSIN(DIP,EXGSM)
C
      SPS=DIP1*S1+DIP2*S2+DIP3*S3
      CPS=SQRT(1.-SPS**2)
      PSI=ASIN(SPS)
C
C    THE ELEMENTS OF THE MATRIX MAG TO SM ARE THE SCALAR PRODUCTS:
C CFI=GM22=(EYSM,EYMAG), SFI=GM23=(EYSM,EXMAG); THEY CAN BE DERIVED AS FOLLOWS:
C
C IN GEO THE VECTORS EXMAG AND EYMAG HAVE THE COMPONENTS (CT0*CL0,CT0*SL0,-ST0)
C  AND (-SL0,CL0,0), RESPECTIVELY.    HENCE, IN GEI THE COMPONENTS ARE:
C  EXMAG:    CT0*CL0*COS(GST)-CT0*SL0*SIN(GST)
C            CT0*CL0*SIN(GST)+CT0*SL0*COS(GST)
C            -ST0
C  EYMAG:    -SL0*COS(GST)-CL0*SIN(GST)
C            -SL0*SIN(GST)+CL0*COS(GST)
C             0
C  THE COMPONENTS OF EYSM IN GEI WERE FOUND ABOVE AS Y1, Y2, AND Y3;
C  NOW WE ONLY HAVE TO COMBINE THE QUANTITIES INTO SCALAR PRODUCTS:
C
      EXMAGX=CT0*(CL0*CGST-SL0*SGST)
      EXMAGY=CT0*(CL0*SGST+SL0*CGST)
      EXMAGZ=-ST0
      EYMAGX=-(SL0*CGST+CL0*SGST)
      EYMAGY=-(SL0*SGST-CL0*CGST)
      CFI=Y1*EYMAGX+Y2*EYMAGY
      SFI=Y1*EXMAGX+Y2*EXMAGY+Y3*EXMAGZ
C
      XMUT=(ATAN2(SFI,CFI)+3.1415926536)*3.8197186342
C
C  THE ELEMENTS OF THE MATRIX GEO TO GSM ARE THE SCALAR PRODUCTS:
C
C   A11=(EXGEO,EXGSM), A12=(EYGEO,EXGSM), A13=(EZGEO,EXGSM),
C   A21=(EXGEO,EYGSM), A22=(EYGEO,EYGSM), A23=(EZGEO,EYGSM),
C   A31=(EXGEO,EZGSM), A32=(EYGEO,EZGSM), A33=(EZGEO,EZGSM),
C
C   ALL THE UNIT VECTORS IN BRACKETS ARE ALREADY DEFINED IN GEI:
C
C  EXGEO=(CGST,SGST,0), EYGEO=(-SGST,CGST,0), EZGEO=(0,0,1)
C  EXGSM=(S1,S2,S3),  EYGSM=(Y1,Y2,Y3),   EZGSM=(Z1,Z2,Z3)
C                                                           AND  THEREFORE:
C
      A11=S1*CGST+S2*SGST
      A12=-S1*SGST+S2*CGST
      A13=S3
      A21=Y1*CGST+Y2*SGST
      A22=-Y1*SGST+Y2*CGST
      A23=Y3
      A31=Z1*CGST+Z2*SGST
      A32=-Z1*SGST+Z2*CGST
      A33=Z3
C
 10   FORMAT(//1X,
     * '****RECALC WARNS:  YEAR IS OUT OF INTERVAL 1965-2005: IYR=',I4,
     * /,6X,'CALCULATIONS WILL BE DONE FOR IYR=',I4,/)
      RETURN
      END

C
      SUBROUTINE GEOMAG(XGEO,YGEO,ZGEO,XMAG,YMAG,ZMAG,J,IYR)
C
C    CONVERTS GEOGRAPHIC (GEO) TO DIPOLE (MAG) COORDINATES OR VICA VERSA.
C    IYR IS YEAR NUMBER (FOUR DIGITS).
C
C                           J &Gt  0                J &Lt 0
C-----INPUT:  J,XGEO,YGEO,ZGEO,IYR   J,XMAG,YMAG,ZMAG,IYR
C-----OUTPUT:    XMAG,YMAG,ZMAG        XGEO,YGEO,ZGEO
C
C   OTHER SUBROUTINES CALLED BY THIS ONE:  RECALC
C
C   LAST MOFIFICATION:  JAN 5, 2001 (NO ESSENTIAL CHANGES, BUT
C                         SOME REDUNDANT STATEMENTS TAKEN OUT)
C
C   WRITTEN BY:  N. A. TSYGANENKO
C
      COMMON /GEOPACK/ ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,AB(19),
     * K,IY,BB(8)
      SAVE II
      DATA II/1/
C
      IF(IYR.EQ.II) GOTO 1
      II=IYR
      CALL RECALC(II,0,25,0,0)
  1   CONTINUE
      IF(J.LT.0) GOTO 2
      XMAG=XGEO*CTCL+YGEO*CTSL-ZGEO*ST0
      YMAG=YGEO*CL0-XGEO*SL0
      ZMAG=XGEO*STCL+YGEO*STSL+ZGEO*CT0
      RETURN
  2   XGEO=XMAG*CTCL-YMAG*SL0+ZMAG*STCL
      YGEO=XMAG*CTSL+YMAG*CL0+ZMAG*STSL
      ZGEO=ZMAG*CT0-XMAG*ST0
      RETURN
      END
c

c
      SUBROUTINE GEIGEO(XGEI,YGEI,ZGEI,XGEO,YGEO,ZGEO,J)
C
C   CONVERTS EQUATORIAL INERTIAL (GEI) TO GEOGRAPHICAL (GEO) COORDS
C   OR VICA VERSA.
C                    J &Gt 0                J &Lt 0
C----INPUT:  J,XGEI,YGEI,ZGEI    J,XGEO,YGEO,ZGEO
C----OUTPUT:   XGEO,YGEO,ZGEO      XGEI,YGEI,ZGEI
C
C  ATTENTION:  SUBROUTINE  RECALC  MUST BE INVOKED BEFORE GEIGEO IN TWO CASES:
C     /A/  BEFORE THE FIRST INVOCATION OF GEIGEO
C     /B/  IF THE CURRENT VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC HAVE BEEN CHANGED
C
C     LAST MODIFICATION:  JAN 5, 2001 (NO ESSENTIAL CHANGES, BUT
C                         SOME REDUNDANT STATEMENTS TAKEN OUT)

C     WRITTEN BY:  N. A. TSYGANENKO
C
      COMMON /GEOPACK/ A(27),KKKK(2),CGST,SGST,B(6)
C
      IF(J.LT.0) GOTO 1
      XGEO=XGEI*CGST+YGEI*SGST
      YGEO=YGEI*CGST-XGEI*SGST
      ZGEO=ZGEI
      RETURN
  1   XGEI=XGEO*CGST-YGEO*SGST
      YGEI=YGEO*CGST+XGEO*SGST
      ZGEI=ZGEO
      RETURN
      END
C

C
      SUBROUTINE MAGSM(XMAG,YMAG,ZMAG,XSM,YSM,ZSM,J)
C
C  CONVERTS DIPOLE (MAG) TO SOLAR MAGNETIC (SM) COORDINATES OR VICA VERSA
C
C                    J &Gt 0              J &Lt 0
C----INPUT: J,XMAG,YMAG,ZMAG     J,XSM,YSM,ZSM
C----OUTPUT:    XSM,YSM,ZSM       XMAG,YMAG,ZMAG
C
C  ATTENTION:  SUBROUTINE  RECALC  MUST BE INVOKED BEFORE MAGSM IN TWO CASES:
C     /A/  BEFORE THE FIRST INVOCATION OF MAGSM
C     /B/  IF THE VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC HAVE BEEN CHANGED
C
C     LAST MODIFICATION:  JAN 5, 2001 (NO ESSENTIAL CHANGES, BUT
C                            SOME REDUNDANT STATEMENTS TAKEN OUT)
C
C     WRITTEN BY:  N. A. TSYGANENKO
C
      COMMON /GEOPACK/ A(8),SFI,CFI,B(7),AB(10),K,IY,BA(8)
C
      IF (J.LT.0) GOTO 1
      XSM=XMAG*CFI-YMAG*SFI
      YSM=XMAG*SFI+YMAG*CFI
      ZSM=ZMAG
      RETURN
  1   XMAG=XSM*CFI+YSM*SFI
      YMAG=YSM*CFI-XSM*SFI
      ZMAG=ZSM
      RETURN
      END

C
       SUBROUTINE GSMGSE(XGSM,YGSM,ZGSM,XGSE,YGSE,ZGSE,J)
C
C CONVERTS GEOCENTRIC SOLAR MAGNETOSPHERIC (GSM) COORDS TO SOLAR ECLIPTIC (GSE) ONES
C   OR VICA VERSA.
C                    J &Gt 0                J &Lt 0
C-----INPUT: J,XGSM,YGSM,ZGSM    J,XGSE,YGSE,ZGSE
C----OUTPUT:   XGSE,YGSE,ZGSE      XGSM,YGSM,ZGSM
C
C  ATTENTION:  SUBROUTINE  RECALC  MUST BE INVOKED BEFORE GSMGSE IN TWO CASES:
C     /A/  BEFORE THE FIRST INVOCATION OF GSMGSE
C     /B/  IF THE VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC HAVE BEEN CHANGED
C
C     LAST MODIFICATION:  JAN 5, 2001  (NO ESSENTIAL CHANGES, BUT
C                             SOME REDUNDANT STATEMENTS TAKEN OUT)
C
C     WRITTEN BY:  N. A. TSYGANENKO
C
      COMMON /GEOPACK/ A(12),SHI,CHI,AB(13),K,IY,BA(8)
C
      IF (J.LT.0) GOTO 1
      XGSE=XGSM
      YGSE=YGSM*CHI-ZGSM*SHI
      ZGSE=YGSM*SHI+ZGSM*CHI
      RETURN
1     XGSM=XGSE
      YGSM=YGSE*CHI+ZGSE*SHI
      ZGSM=ZGSE*CHI-YGSE*SHI
      RETURN
      END

C
       SUBROUTINE SMGSM(XSM,YSM,ZSM,XGSM,YGSM,ZGSM,J)
C
C CONVERTS SOLAR MAGNETIC (SM) TO GEOCENTRIC SOLAR MAGNETOSPHERIC
C   (GSM) COORDINATES OR VICA VERSA.
C                  J &Gt 0                 J &Lt 0
C-----INPUT: J,XSM,YSM,ZSM        J,XGSM,YGSM,ZGSM
C----OUTPUT:  XGSM,YGSM,ZGSM       XSM,YSM,ZSM
C
C  ATTENTION:  SUBROUTINE  RECALC  MUST BE INVOKED BEFORE SMGSM IN TWO CASES:
C     /A/  BEFORE THE FIRST INVOCATION OF SMGSM
C     /B/  IF THE VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC HAVE BEEN CHANGED
C
C     LAST MODIFICATION:  JAN 5, 2001  (NO ESSENTIAL CHANGES, BUT
C                             SOME REDUNDANT STATEMENTS TAKEN OUT)
C
C     WRITTEN BY:  N. A. TSYGANENKO
C
      COMMON /GEOPACK/ A(10),SPS,CPS,B(15),K,IY,AB(8)
      IF (J.LT.0) GOTO 1
      XGSM=XSM*CPS+ZSM*SPS
      YGSM=YSM
      ZGSM=ZSM*CPS-XSM*SPS
      RETURN
  1   XSM=XGSM*CPS-ZGSM*SPS
      YSM=YGSM
      ZSM=XGSM*SPS+ZGSM*CPS
      RETURN
      END

C
      SUBROUTINE GEOGSM(XGEO,YGEO,ZGEO,XGSM,YGSM,ZGSM,J)
C
C CONVERTS GEOGRAPHIC (GEO) TO SOLAR MAGNETOSPHERIC (SM) COORDINATES OR VICA VERSA.
C
C                   J &Gt 0                   J &Lt 0
C----- INPUT:  J,XGEO,YGEO,ZGEO    J,XGSM,YGSM,ZGSM
C---- OUTPUT:    XGSM,YGSM,ZGSM      XGEO,YGEO,ZGEO
C
C  ATTENTION:  SUBROUTINE  RECALC  MUST BE INVOKED BEFORE GEOGSM IN TWO CASES:
C     /A/  BEFORE THE FIRST INVOCATION OF GEOGSM
C     /B/  IF THE VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC  HAVE BEEN CHANGED
C
C     LAST MODIFICATION:  JAN 5, 2001 (NO ESSENTIAL CHANGES, BUT
C                            SOME REDUNDANT STATEMENTS TAKEN OUT)
C
C     WRITTEN BY:  N. A. TSYGANENKO
C
      COMMON /GEOPACK/ AA(17),A11,A21,A31,A12,A22,A32,A13,A23,A33,
     *  D,K,IY,B(8)
C
      IF (J.LT.0) GOTO 1
      XGSM=A11*XGEO+A12*YGEO+A13*ZGEO
      YGSM=A21*XGEO+A22*YGEO+A23*ZGEO
      ZGSM=A31*XGEO+A32*YGEO+A33*ZGEO
      RETURN
  1   XGEO=A11*XGSM+A21*YGSM+A31*ZGSM
      YGEO=A12*XGSM+A22*YGSM+A32*ZGSM
      ZGEO=A13*XGSM+A23*YGSM+A33*ZGSM
      RETURN
      END

C
      SUBROUTINE RHAND(X,Y,Z,R1,R2,R3,IOPT,PARMOD,EXNAME)
C
C  COMPUTES RIGHT HAND EXPRESSIONS IN THE GEOMAGNETIC FIELD LINE EQUATION
C      (a subsidiary subroutine for the subr. STEP)
C
C     LAST MODIFICATION:  JAN 5, 2001 (NO ESSENTIAL CHANGES, BUT
C                            SOME REDUNDANT STATEMENTS TAKEN OUT)
C
C     WRITTEN BY:  N. A. TSYGANENKO
C
      DIMENSION PARMOD(10)
      EXTERNAL EXNAME
C
      COMMON /GEOPACK/ A(15),PSI,AA(10),DS3,K,IY,BB(8)
      CALL EXNAME(IOPT,PARMOD,PSI,X,Y,Z,BX,BY,BZ)
      IF (K.EQ.0) GOTO 1
      CALL GEOGSM(XG,YG,ZG,X,Y,Z,-1)
      CALL SPHCAR(R,T,F,XG,YG,ZG,-1)
      CALL IGRF(IY,K,R,T,F,BR,BT,BF)
      CALL BSPCAR(T,F,BR,BT,BF,FX,FY,FZ)
      CALL GEOGSM(FX,FY,FZ,HX,HY,HZ,1)
      GOTO 2
  1   CALL DIP(PSI,X,Y,Z,HX,HY,HZ)
  2   BX=BX+HX
      BY=BY+HY
      BZ=BZ+HZ
      B=DS3/SQRT(BX**2+BY**2+BZ**2)
      R1=BX*B
      R2=BY*B
      R3=BZ*B
      RETURN
      END

C
      SUBROUTINE STEP(N,X,Y,Z,DS,ERRIN,IOPT,PARMOD,EXNAME)
C
C RE-CALCULATES COORDS X,Y,Z FOR ONE STEP ALONG FIELD LINE. N IS MAXIMUM
C ORDER OF HARMONICS IN MAIN FIELD EXPANSION, DS IS STEP SIZE, ERRIN IS
C PERMISSIBLE ERROR VALUE, IOPT AND EXNAME - SEE COMMENTS TO SUBROUTINE TRACE
C  ALL THE PARAMETERS ARE INPUT ONES; OUTPUT IS THE RENEWED TRIPLET X,Y,Z
C
C     LAST MODIFICATION:  JAN 5, 2001 (NO ESSENTIAL CHANGES, BUT
C                            SOME REDUNDANT STATEMENTS TAKEN OUT)
C
C     WRITTEN BY:  N. A. TSYGANENKO
C
C
      DIMENSION PARMOD(10)
      COMMON /GEOPACK/ A(26),DS3,K,IY,B(8)
      EXTERNAL EXNAME
      K=N
  1   DS3=-DS/3.
      CALL RHAND(X,Y,Z,R11,R12,R13,IOPT,PARMOD,EXNAME)
      CALL RHAND(X+R11,Y+R12,Z+R13,R21,R22,R23,IOPT,PARMOD,EXNAME)
      CALL RHAND(X+.5*(R11+R21),Y+.5*(R12+R22),Z+.5*
     *(R13+R23),R31,R32,R33,IOPT,PARMOD,EXNAME)
      CALL RHAND(X+.375*(R11+3.*R31),Y+.375*(R12+3.*R32
     *),Z+.375*(R13+3.*R33),R41,R42,R43,IOPT,PARMOD,EXNAME)
      CALL RHAND(X+1.5*(R11-3.*R31+4.*R41),Y+1.5*(R12-
     *3.*R32+4.*R42),Z+1.5*(R13-3.*R33+4.*R43),
     *R51,R52,R53,IOPT,PARMOD,EXNAME)
      ERRCUR=ABS(R11-4.5*R31+4.*R41-.5*R51)+ABS(R12-4.5*R32+4.*R42-.5*
     *R52)+ABS(R13-4.5*R33+4.*R43-.5*R53)
      IF (ERRCUR.LT.ERRIN) GOTO 2
      DS=DS*.5
      GOTO 1
  2   X=X+.5*(R11+4.*R41+R51)
      Y=Y+.5*(R12+4.*R42+R52)
      Z=Z+.5*(R13+4.*R43+R53)
      IF(ERRCUR.LT.ERRIN*.04.AND.ABS(DS).LT.1.33) DS=DS*1.5
      RETURN
      END

C
      SUBROUTINE TRACE(XI,YI,ZI,DIR,RLIM,R0,IHARM,NP,IOPT,PARMOD,EXNAME,
     *XF,YF,ZF,XX,YY,ZZ,L)
C
C   TRACES A FIELD LINE FROM AN ARBITRARY POINT OF SPACE TO THE EARTH'S
C   SURFACE OR TO A MODEL LIMITING BOUNDARY.
C
C------------- INPUT PARAMETERS:
C
C   XI,YI,ZI - GSM COORDS OF INITIAL POINT (IN EARTH RADII),
C
C   DIR - SIGN OF TRACING DIRECTION: IF DIR=1. THEN ANTIPARALLEL TO
C     B VECTOR (E.G. FROM NORTHERN TO SOUTHERN CONJUGATE POINT),
C     AND IF DIR=-1. THEN PARALLEL TO B.
C
C   R0 -  RADIUS OF A SPHERE (IN RE) FOR WHICH FOOTPOINT POINT COORDINATES
C     XF,YF,ZF  ARE TO BE CALCULATED.
C
C   RLIM - UPPER LIMIT OF THE GEOCENTRIC DISTANCE, WHERE THE TRACING IS TERMINATED.
C
C   IHARM - THE HIGHEST ORDER OF SPHERICAL HARMONICS IN THE MAIN FIELD EXPANSION,
C     THIS PARAMETER DEPENDS ON A VERSION OF THE MAIN GEOMAGNETIC FIELD MODEL.
C     THE SUBROUTINE IGRF INCLUDED IN THIS PACKAGE HAS IHARM=10.
C     IF THE MAIN FIELD CAN BE ASSUMED PURELY DIPOLAR, THEN TAKE IHARM=0
C     AND SPECIFY A VALUE OF THE GEODIPOLE TILT ANGLE PSI (IN RADIANS)
C     IN THE 16-TH ELEMENT OF THE COMMON BLOCK /GEOPACK/ BEFORE INVOKING TRACE.
C     OTHERWISE INVOKE  RECALC  BEFORE USING  TRACE.
C
C   NP - UPPER ESTIMATE OF THE NUMBER OF STEPS ALONG A FIELD LINE
C     (IN MOST CASES, NP=500 IS SUFFICIENT).
C
C   IOPT - A MODEL INDEX; CAN BE USED FOR SPECIFYING AN OPTION OF THE EXTERNAL FIELD
C       MODEL (E.G., INTERVAL OF THE KP-INDEX). ALTERNATIVELY, ONE CAN USE THE ARRAY
C       PARMOD FOR THAT PURPOSE (SEE BELOW); IN THAT CASE IOPT CAN BE JUST IGNORED.
C
C    PARMOD -  A 10-ELEMENT ARRAY CONTAINING MODEL PARAMETERS, NEEDED FOR A UNIQUE
C      SPECIFICATION OF THE EXTERNAL FIELD. THE CONCRETE MEANING OF THE COMPONENTS
C      OF PARMOD DEPENDS ON A SPECIFIC VERSION OF THE EXTERNAL FIELD MODEL.
C
C    EXNAME - NAME OF A SUBROUTINE PROVIDING COMPONENTS OF THE EXTERNAL MAGNETIC FIELD.
C
C-------------- OUTPUT PARAMETERS:
C
C   XF,YF,ZF - GSM COORDS OF THE LAST CALCULATED POINT OF A FIELD LINE
C   XX,YY,ZZ - ARRAYS, CONTAINING COORDS OF FIELD LINE POINTS. THEIR SIZE SHOULD BE
C      NOT LESS THAN  NP.
C   L - ACTUAL NUMBER OF THE CALCULATED FIELD LINE POINTS. IF L EXCEEDS NP, THE TRACING
C     TERMINATES, AND A WARNING IS DISPLAYED.
C
C
C     LAST MODIFICATION:  JAN 5, 2001 (NO ESSENTIAL CHANGES)
C
C     WRITTEN BY:  N. A. TSYGANENKO
C
      DIMENSION XX(NP),YY(NP),ZZ(NP), PARMOD(10)
      COMMON /GEOPACK/ AA(26),DD,K1,K2,BB(8)
      EXTERNAL EXNAME

      J=IHARM
      ERR=0.0005
      L=0
      DS=0.5*DIR
      X=XI
      Y=YI
      Z=ZI
      DD=DIR
      K1=IHARM
      AL=0.
      CALL RHAND(X,Y,Z,R1,R2,R3,IOPT,PARMOD,EXNAME)
      AD=SIGN(0.01,X*R1+Y*R2+Z*R3)
      RR=SQRT(X**2+Y**2+Z**2)+AD
  1   L=L+1
      IF(L.GT.NP) GOTO 7
      XX(L)=X
      YY(L)=Y
      ZZ(L)=Z
      RYZ=Y**2+Z**2
      R2=X**2+RYZ
      R=SQRT(R2)
      IF(R.GT.RLIM.OR.RYZ.GT.1600..OR.X.GT.20.) GOTO 8
      IF(R.LT.R0.AND.RR.GT.R) GOTO 6
      IF(R.GE.RR) GOTO 5
      IF(R.GT.5.) GOTO 5
      IF(R.GE.3.) GOTO 3
      FC=0.2
      IF(R-R0.LT.0.05) FC=0.05
      AL=FC*(R-R0+0.2)
      DS=DIR*AL
      GOTO 4
  3   DS=DIR
      AL=1.
  4   XR=X
      YR=Y
      ZR=Z
  5   RR=R
      IF(IHARM.NE.0) J=1+20./(R-AL)
      IF(J.GT.IHARM) J=IHARM
      CALL STEP(J,X,Y,Z,DS,ERR,IOPT,PARMOD,EXNAME)
      GOTO 1
  6   R1=(R0-R)/(RR-R)
      X=X-(X-XR)*R1
      Y=Y-(Y-YR)*R1
      Z=Z-(Z-ZR)*R1
      GOTO 8
c MCF
c 7   WRITE (*,10)
  7   continue     
c MCF end
      L=NP
      RETURN
  8   XF=X
      YF=Y
      ZF=Z
      DO 9 I=L,NP
      XX(I)=XF
      YY(I)=YF
  9   ZZ(I)=ZF
      RETURN
 10   FORMAT(//,1X,'**** COMPUTATIONS IN THE SUBROUTINE TRACE',
     *' ARE TERMINATED: NP IS TOO SMALL ****'//)
      END


      subroutine trace1(xism,yism,zism,DIR,RLIM,R0,IHARM,NP,IOPT,
     *  PARMOD,EXNAME,XFsm,YFsm,ZFsm,XXsm,YYsm,ZZsm,L,iout)
c
c   A modification of TRACE. This routine input sm xi, yi, zi, output iout, 
c   the tracing stops at the magnetic equator, and output fieldline points in 
c   sm coordinates. All the changes are in lower case.
c   (Mei-Ching Fok Jan. 11, 2001. Last modification: July 30, 2001)
C
C   TRACES A FIELD LINE FROM AN ARBITRARY POINT OF SPACE TO THE EARTH'S
C   SURFACE OR TO A MODEL LIMITING BOUNDARY.
C
C------------- INPUT PARAMETERS:
C
C   XI,YI,ZI - GSM COORDS OF INITIAL POINT (IN EARTH RADII),
C
C   DIR - SIGN OF TRACING DIRECTION: IF DIR=1. THEN ANTIPARALLEL TO
C     B VECTOR (E.G. FROM NORTHERN TO SOUTHERN CONJUGATE POINT),
C     AND IF DIR=-1. THEN PARALLEL TO B.
C
C   R0 -  RADIUS OF A SPHERE (IN RE) FOR WHICH FOOTPOINT POINT COORDINATES
C     XF,YF,ZF  ARE TO BE CALCULATED.
C
C   RLIM - UPPER LIMIT OF THE GEOCENTRIC DISTANCE, WHERE THE TRACING IS TERMINATED.
C
C   IHARM - THE HIGHEST ORDER OF SPHERICAL HARMONICS IN THE MAIN FIELD EXPANSION,
C     THIS PARAMETER DEPENDS ON A VERSION OF THE MAIN GEOMAGNETIC FIELD MODEL.
C     THE SUBROUTINE IGRF INCLUDED IN THIS PACKAGE HAS IHARM=10.
C     IF THE MAIN FIELD CAN BE ASSUMED PURELY DIPOLAR, THEN TAKE IHARM=0
C     AND SPECIFY A VALUE OF THE GEODIPOLE TILT ANGLE PSI (IN RADIANS)
C     IN THE 16-TH ELEMENT OF THE COMMON BLOCK /GEOPACK/ BEFORE INVOKING TRACE.
C     OTHERWISE INVOKE  RECALC  BEFORE USING  TRACE.
C
C   NP - UPPER ESTIMATE OF THE NUMBER OF STEPS ALONG A FIELD LINE
C     (IN MOST CASES, NP=500 IS SUFFICIENT).
C
C   IOPT - A MODEL INDEX; CAN BE USED FOR SPECIFYING AN OPTION OF THE EXTERNAL FIELD
C       MODEL (E.G., INTERVAL OF THE KP-INDEX). ALTERNATIVELY, ONE CAN USE THE ARRAY
C       PARMOD FOR THAT PURPOSE (SEE BELOW); IN THAT CASE IOPT CAN BE JUST IGNORED.
C
C    PARMOD -  A 10-ELEMENT ARRAY CONTAINING MODEL PARAMETERS, NEEDED FOR A UNIQUE
C      SPECIFICATION OF THE EXTERNAL FIELD. THE CONCRETE MEANING OF THE COMPONENTS
C      OF PARMOD DEPENDS ON A SPECIFIC VERSION OF THE EXTERNAL FIELD MODEL.
C
C    EXNAME - NAME OF A SUBROUTINE PROVIDING COMPONENTS OF THE EXTERNAL MAGNETIC FIELD.
C
C-------------- OUTPUT PARAMETERS:
C
C   XF,YF,ZF - GSM COORDS OF THE LAST CALCULATED POINT OF A FIELD LINE
C   XX,YY,ZZ - ARRAYS, CONTAINING COORDS OF FIELD LINE POINTS. THEIR SIZE SHOULD BE
C      NOT LESS THAN  NP.
C   L - ACTUAL NUMBER OF THE CALCULATED FIELD LINE POINTS. IF L EXCEEDS NP, THE TRACING
C     TERMINATES, AND A WARNING IS DISPLAYED.
C
C
C     LAST MODIFICATION:  JAN 5, 2001 (NO ESSENTIAL CHANGES)
C
C     WRITTEN BY:  N. A. TSYGANENKO
C
      REAL XX(NP),YY(NP),ZZ(NP), PARMOD(10),XXsm(NP),YYsm(NP),ZZsm(NP)
      COMMON /GEOPACK/ AA(26),DD,K1,K2,BB(8)
      EXTERNAL EXNAME

      iout=0
      xn_pd=parmod(1)    ! solar wind ram pressure in nPa
      vel=-1.
      J=IHARM
      ERR=0.0005
      L=0
      DS=0.3*DIR
      xrsm=xism
      yrsm=yism
      zrsm=zism
      j1=1
      call smgsm(xism,yism,zism,xi,yi,zi,j1)
      X=XI
      Y=YI
      Z=ZI
      DD=DIR
      K1=IHARM
      AL=0.
      CALL RHAND(X,Y,Z,R1,R2,R3,IOPT,PARMOD,EXNAME)
      AD=SIGN(0.01,X*R1+Y*R2+Z*R3)
      RR=SQRT(X**2+Y**2+Z**2)+AD
  1   L=L+1
      IF(L.GT.NP) GOTO 7
      XX(L)=X
      YY(L)=Y
      ZZ(L)=Z
      j1=-1
      call smgsm(xsm,ysm,zsm,x,y,z,j1)
      XXsm(L)=Xsm
      YYsm(L)=Ysm
      ZZsm(L)=Zsm
      RYZ=Y**2+Z**2
      R2=X**2+RYZ
      R=SQRT(R2)
      call locatep(xn_pd,vel,x,y,z,xmgnp,ymgnp,zmgnp,dist,id)
c     IF(R.GT.RLIM.OR.RYZ.GT.1600..OR.X.GT.20.) GOTO 8
      if(r.gt.rlim.or.ryz.gt.1600..or.x.gt.20..or.id.eq.-1) then
         iout=1         ! hit open field lines or beyond simulation box
         goto 8
      endif
      zsign=zsm*zrsm
      if (zsign.le.0.) goto 99            ! cross magnetic equator
      IF(R.LT.R0.AND.RR.GT.R) GOTO 6
      IF(R.GE.RR) GOTO 5
      IF(R.GT.5.) GOTO 5
      IF(R.GE.3.) GOTO 3
      FC=0.2
      IF(R-R0.LT.0.05) FC=0.05
      AL=FC*(R-R0+0.2)
      DS=DIR*AL
      GOTO 4
  3   DS=DIR
      AL=1.
  4   XR=X
      YR=Y
      ZR=Z
  5   RR=R
      xrsm=xsm
      yrsm=ysm
      zrsm=zsm
      IF(IHARM.NE.0) J=1+20./(R-AL)
      IF(J.GT.IHARM) J=IHARM
      CALL STEP(J,X,Y,Z,DS,ERR,IOPT,PARMOD,EXNAME)
      GOTO 1
 99   r1=-zsm/(zrsm-zsm)
      Xsm=Xsm-(Xsm-xrsm)*R1
      Ysm=Ysm-(Ysm-yrsm)*R1
      Zsm=Zsm-(Zsm-zrsm)*R1
      j1=1
      call smgsm(xsm,ysm,zsm,x,y,z,j1)
      goto 8
  6   R1=(R0-R)/(RR-R)
      X=X-(X-XR)*R1
      Y=Y-(Y-YR)*R1
      Z=Z-(Z-ZR)*R1
      GOTO 8
  7   WRITE (*,10)
      L=NP
      RETURN
  8   XF=X
      YF=Y
      ZF=Z
      XFsm=Xsm
      YFsm=Ysm
      ZFsm=Zsm
      DO 9 I=L,NP
      XXsm(I)=XFsm
      YYsm(I)=YFsm
      ZZsm(I)=ZFsm
      XX(I)=XF
      YY(I)=YF
  9   ZZ(I)=ZF
      RETURN
 10   FORMAT(//,1X,'**** COMPUTATIONS IN THE SUBROUTINE TRACE',
     *' ARE TERMINATED: NP IS TOO SMALL ****'//)
      END


C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C                                                                          &
C          SUPPLEMENTARY CODES:   LOCATE  and   CROSSING                   &
C                                                                          &
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

c  The name of this routine has changed to LOCATEp, Mei-Ching Fok, Jan. 11, 2001

      SUBROUTINE LOCATEp(XN_PD,VEL,XGSM,YGSM,ZGSM,XMGNP,
     &                   YMGNP,ZMGNP,DIST,ID)
C
C     THIS SUBROUTINE DEFINES THE POSITION OF A POINT (XMGNP,YMGNP,ZMGNP)
C            AT THE MODEL MAGNETOPAUSE, CLOSEST TO A GIVEN POINT OF SPACE
C               (XGSM,YGSM,ZGSM),   AND THE DISTANCE BETWEEN THEM (DIST)
C
C INPUT:  XN_PD - EITHER SOLAR WIND PROTON NUMBER DENSITY (PER C.C.) (IF VEL>0)
C                    OR THE SOLAR WIND RAM PRESSURE IN NANOPASCALS   (IF VEL<0)
C         VEL - EITHER SOLAR WIND VELOCITY (KM/SEC)
C                  OR ANY NEGATIVE NUMBER, WHICH INDICATES THAT XN_PD STANDS
C                     FOR THE SOLAR WIND PRESSURE, RATHER THAN FOR THE DENSITY
C
C         XGSM,YGSM,ZGSM - COORDINATES OF THE OBSERVATION POINT IN EARTH RADII
C
C  OUTPUT: XMGNP,YMGNP,ZMGNP - COORDINATES OF A POINT AT THE MAGNETOPAUSE,
C                                    CLOSEST TO THE POINT  XGSM,YGSM,ZGSM
C          DIST -  THE DISTANCE BETWEEN THE ABOVE TWO POINTS, IN RE,
C          ID -    INDICATOR; ID=+1 AND ID=-1 MEAN THAT THE POINT
C                       (XGSM,YGSM,ZGSM)  LIES INSIDE OR OUTSIDE
C                        THE MODEL MAGNETOPAUSE, RESPECTIVELY
C
C   THE PRESSURE-DEPENDENT MAGNETOPAUSE IS THAT USED IN THE T96_01 MODEL
C
c            CODED BY:  N.A. TSYGANENKO, AUG.1, 1995;  REVISED  JUNE 22, 1996
C
      IF (VEL.LT.0.) THEN
       PD=XN_PD
                  ELSE
       PD=1.94E-6*XN_PD*VEL**2  ! PD IS THE SOLAR WIND DYNAMIC PRESSURE
C                                      (IN NANOPASCALS)
      ENDIF
C
      RAT=PD/2.0  ! RATIO OF PD TO THE AVERAGE PRESSURE, ASSUMED AS 2 nPa
C
      RAT16=RAT**0.14 ! THE POWER IN THE SCALING FACTOR IS THE BEST-FIT VALUE
C                         OBTAINED FROM DATA IN THE T96_01 VERSION OF THE MODEL
C
      A0=70.
      S00=1.08
      X00=5.48    !  VALUES OF THE MAGNETOPAUSE PARAMETERS FOR  PD = 2 nPa
C
      A=A0/RAT16
      S0=S00
      X0=X00/RAT16   !  VALUES OF THE MAGNETOPAUSE PARAMETERS, SCALED TO THE
C                         ACTUAL PRESSURE
C
       XM=X0-A    !  THIS IS THE X-COORDINATE OF THE "SEAM" BETWEEN THE
C                             ELLIPSOID AND THE CYLINDER
C
C        (FOR DETAILS ON THE ELLIPSOIDAL COORDINATES, SEE THE PAPER:
C            N.A.TSYGANENKO, SOLUTION OF CHAPMAN-FERRARO PROBLEM FOR AN
C             ELLIPSOIDAL MAGNETOPAUSE, PLANET.SPACE SCI., V.37, P.1037, 1989).
C
          IF (YGSM.NE.0..OR.ZGSM.NE.0.) THEN
             PHI=ATAN2(YGSM,ZGSM)
              ELSE
             PHI=0.
          ENDIF
C
          RHO=SQRT(YGSM**2+ZGSM**2)
C
         IF (XGSM.LT.XM) THEN
           XMGNP=XGSM
           RHOMGNP=A*SQRT(S0**2-1)
           YMGNP=RHOMGNP*SIN(PHI)
           ZMGNP=RHOMGNP*COS(PHI)
           DIST=SQRT((XGSM-XMGNP)**2+(YGSM-YMGNP)**2+(ZGSM-ZMGNP)**2)
           IF (RHOMGNP.GT.RHO) ID=+1
           IF (RHOMGNP.LT.RHO) ID=-1
           RETURN
         ENDIF
C
          XKSI=(XGSM-X0)/A+1.
          XDZT=RHO/A
          SQ1=SQRT((1.+XKSI)**2+XDZT**2)
          SQ2=SQRT((1.-XKSI)**2+XDZT**2)
          SIGMA=0.5*(SQ1+SQ2)
          TAU=0.5*(SQ1-SQ2)
C
C  NOW CALCULATE (X,Y,Z) FOR THE CLOSEST POINT AT THE MAGNETOPAUSE
C
          XMGNP=X0-A*(1.-S0*TAU)
          RHOMGNP=A*SQRT((S0**2-1.)*(1.-TAU**2))
          YMGNP=RHOMGNP*SIN(PHI)
          ZMGNP=RHOMGNP*COS(PHI)
C
C  NOW CALCULATE THE SHORTEST DISTANCE BETWEEN THE POINT XGSM,YGSM,ZGSM AND THE
C            MAGNETOPAUSE
C
      DIST=SQRT((XGSM-XMGNP)**2+(YGSM-YMGNP)**2+(ZGSM-ZMGNP)**2)
C
      IF (SIGMA.GT.S0) ID=-1   !  ID=-1 MEANS THAT THE POINT LIES OUTSIDE
      IF (SIGMA.LT.S0) ID=+1   !  ID=+1 MEANS THAT THE POINT LIES INSIDE
C                                           THE MAGNETOSPHERE
      RETURN
      END
