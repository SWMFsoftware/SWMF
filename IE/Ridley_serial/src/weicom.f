C          Routines common to both 1996 and 2001 versions of Weimer's
C          Ionospheric Electrostatic Potential Model

************************ Copyright 1996,2001 Dan Weimer/MRC ***********************
	FUNCTION FSVal(omega,MaxN,FSC)
* Fourier Series Value
* Return value of Sine/Cosine Fourier series for N terms up to MaxN
* at angle omega, given the F.S. coeficients in array FSC
	REAL omega,FSC(0:1,0:*)
	INTEGER MaxN,n
	REAL Y,theta
	Y=0.
	DO n=0,MaxN
	  theta=omega*n
	  Y=Y + FSC(0,n)*COS(theta) + FSC(1,n)*SIN(theta)
	ENDDO
	FSVal=Y
	RETURN
	END
************************ Copyright 1996,2001 Dan Weimer/MRC ***********************
* COORDINATE TRANSFORMATION UTILITIES

CNCAR      Feb 01:  Changed TRANS to GET_TILT s.t. the dipole tilt angle is
C          returned.

	FUNCTION GET_TILT (YEAR,MONTH,DAY,HOUR)
C       SUBROUTINE TRANS(YEAR,MONTH,DAY,HOUR,IDBUG)
CNCAR

	INTEGER YEAR,MONTH,DAY,IDBUG
	REAL HOUR
C         
C      THIS SUBROUTINE DERIVES THE ROTATION MATRICES AM(I,J,K) FOR 11
C      TRANSFORMATIONS, IDENTIFIED BY K.
C          K=1 TRANSFORMS GSE to GEO
C          K=2     "      GEO to MAG
C          K=3     "      GSE to MAG
C          K=4     "      GSE to GSM
C          K=5     "      GEO to GSM
C          K=6     "      GSM to MAG
C          K=7     "      GSE to GEI
C          K=8     "      GEI to GEO
C          K=9     "      GSM to SM 
C	   K=10    "      GEO to SM 
C	   K=11    "      MAG to SM 
C
C      IF IDBUG IS NOT 0, THEN OUTPUTS DIAGNOSTIC INFORMATION TO
C      FILE UNIT=IDBUG
C	
	INTEGER GSEGEO,GEOGSE,GEOMAG,MAGGEO
	INTEGER GSEMAG,MAGGSE,GSEGSM,GSMGSE
	INTEGER GEOGSM,GSMGEO,GSMMAG,MAGGSM
	INTEGER GSEGEI,GEIGSE,GEIGEO,GEOGEI
	INTEGER GSMSM,SMGSM,GEOSM,SMGEO,MAGSM,SMMAG

	PARAMETER (GSEGEO= 1,GEOGSE=-1,GEOMAG= 2,MAGGEO=-2)
	PARAMETER (GSEMAG= 3,MAGGSE=-3,GSEGSM= 4,GSMGSE=-4)
	PARAMETER (GEOGSM= 5,GSMGEO=-5,GSMMAG= 6,MAGGSM=-6)
	PARAMETER (GSEGEI= 7,GEIGSE=-7,GEIGEO= 8,GEOGEI=-8)
	PARAMETER (GSMSM = 9,SMGSM =-9,GEOSM =10,SMGEO=-10)
	PARAMETER (MAGSM =11,SMMAG =-11)
C
C      The formal names of the coordinate systems are:
C	GSE - Geocentric Solar Ecliptic
C	GEO - Geographic
C	MAG - Geomagnetic
C	GSM - Geocentric Solar Magnetospheric
C	SM  - Solar Magnetic
C	
C      THE ARRAY CX(I) ENCODES VARIOUS ANGLES, STORED IN DEGREES
C      ST(I) AND CT(I) ARE SINES & COSINES.       
C
C      Program author:  D. R. Weimer
C
C      Some of this code has been copied from subroutines which had been
C      obtained from D. Stern, NASA/GSFC.  Other formulas are from "Space 
C      Physics Coordinate Transformations: A User Guide" by M. Hapgood (1991).
C
C      The formulas for the calculation of Greenwich mean sidereal time (GMST)
C      and the sun's location are from "Almanac for Computers 1990",
C      U.S. Naval Observatory.
C
	REAL UT,T0,GMSTD,GMSTH,ECLIP,MA,LAMD,SUNLON

	COMMON/MFIELD/EPOCH,TH0,PH0,DIPOLE
	COMMON/TRANSDAT/CX(9),ST(6),CT(6),AM(3,3,11)
C         
C UM - Changed Data statement to the following:

	EPOCH = 1980. 
	TH0 = 11.19
	PH0 = -70.76
	DIPOLE = .30574

C UM

CNCAR      Feb 01:  Prevent debug printing to a file
	IDBUG = 0
CNCAR

	IF(YEAR.LT.1900)THEN
	  IYR=1900+YEAR
	ELSE
	  IYR=YEAR
	ENDIF
	UT=HOUR
	JD=JULDAY(MONTH,DAY,IYR)
	MJD=JD-2400001
	T0=(FLOAT(MJD)-51544.5)/36525.0
	GMSTD=100.4606184 + 36000.770*T0 + 3.87933E-4*T0*T0 +
     $        15.0410686*UT
	CALL ADJUST(GMSTD)
	GMSTH=GMSTD*24./360.
	ECLIP=23.439 - 0.013*T0
	MA=357.528 + 35999.050*T0 + 0.041066678*UT
	CALL ADJUST(MA)
	LAMD=280.460 + 36000.772*T0 + 0.041068642*UT
	CALL ADJUST(LAMD)
	SUNLON=LAMD + (1.915-0.0048*T0)*SIND(MA) + 0.020*SIND(2.*MA)
	CALL ADJUST(SUNLON)
	IF(IDBUG.NE.0)THEN
	  WRITE(IDBUG,*) YEAR,MONTH,DAY,HOUR
	  WRITE(IDBUG,*) 'MJD=',MJD
	  WRITE(IDBUG,*) 'T0=',T0
	  WRITE(IDBUG,*) 'GMSTH=',GMSTH
	  WRITE(IDBUG,*) 'ECLIPTIC OBLIQUITY=',ECLIP
	  WRITE(IDBUG,*) 'MEAN ANOMALY=',MA
	  WRITE(IDBUG,*) 'MEAN LONGITUDE=',LAMD
	  WRITE(IDBUG,*) 'TRUE LONGITUDE=',SUNLON
	ENDIF

	CX(1)= GMSTD
	CX(2) = ECLIP
	CX(3) = SUNLON
	CX(4) = TH0
	CX(5) = PH0
c Derived later:
c       CX(6) = Dipole tilt angle  
c       CX(7) = Angle between sun and magnetic pole
c       CX(8) = Subsolar point latitude
c       CX(9) = Subsolar point longitude

	DO I=1,5
	  ST(I) = SIND(CX(I))
	  CT(I) = COSD(CX(I))
	ENDDO
C         
      AM(1,1,GSEGEI) = CT(3)
      AM(1,2,GSEGEI) = -ST(3)
      AM(1,3,GSEGEI) = 0.         
      AM(2,1,GSEGEI) = ST(3)*CT(2)
      AM(2,2,GSEGEI) = CT(3)*CT(2)
      AM(2,3,GSEGEI) = -ST(2)
      AM(3,1,GSEGEI) = ST(3)*ST(2)
      AM(3,2,GSEGEI) = CT(3)*ST(2)
      AM(3,3,GSEGEI) = CT(2)      
C         
      AM(1,1,GEIGEO) = CT(1)      
      AM(1,2,GEIGEO) = ST(1)      
      AM(1,3,GEIGEO) = 0.         
      AM(2,1,GEIGEO) = -ST(1)     
      AM(2,2,GEIGEO) = CT(1)      
      AM(2,3,GEIGEO) = 0.         
      AM(3,1,GEIGEO) = 0.         
      AM(3,2,GEIGEO) = 0.         
      AM(3,3,GEIGEO) = 1.         
C         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GSEGEO) = AM(I,1,GEIGEO)*AM(1,J,GSEGEI) +
     $    AM(I,2,GEIGEO)*AM(2,J,GSEGEI) +  AM(I,3,GEIGEO)*AM(3,J,GSEGEI)
      ENDDO
      ENDDO
C         
      AM(1,1,GEOMAG) = CT(4)*CT(5) 
      AM(1,2,GEOMAG) = CT(4)*ST(5) 
      AM(1,3,GEOMAG) =-ST(4)       
      AM(2,1,GEOMAG) =-ST(5)       
      AM(2,2,GEOMAG) = CT(5)       
      AM(2,3,GEOMAG) = 0.
      AM(3,1,GEOMAG) = ST(4)*CT(5) 
      AM(3,2,GEOMAG) = ST(4)*ST(5) 
      AM(3,3,GEOMAG) = CT(4)       
C         
      DO I=1,3   
      DO J=1,3   
       AM(I,J,GSEMAG) = AM(I,1,GEOMAG)*AM(1,J,GSEGEO) +
     $   AM(I,2,GEOMAG)*AM(2,J,GSEGEO) +  AM(I,3,GEOMAG)*AM(3,J,GSEGEO)
      ENDDO
      ENDDO
C         
      B32 = AM(3,2,GSEMAG)         
      B33 = AM(3,3,GSEMAG)         
      B3  = SQRT(B32*B32+B33*B33)       
      IF (B33.LE.0.) B3 = -B3    
C         
      AM(2,2,GSEGSM) = B33/B3      
      AM(3,3,GSEGSM) = AM(2,2,GSEGSM)   
      AM(3,2,GSEGSM) = B32/B3      
      AM(2,3,GSEGSM) =-AM(3,2,GSEGSM)   
      AM(1,1,GSEGSM) = 1.
      AM(1,2,GSEGSM) = 0.
      AM(1,3,GSEGSM) = 0.
      AM(2,1,GSEGSM) = 0.
      AM(3,1,GSEGSM) = 0.
C         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GEOGSM) = AM(I,1,GSEGSM)*AM(J,1,GSEGEO) +
     $    AM(I,2,GSEGSM)*AM(J,2,GSEGEO) + AM(I,3,GSEGSM)*AM(J,3,GSEGEO)
      ENDDO
      ENDDO
C         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GSMMAG) = AM(I,1,GEOMAG)*AM(J,1,GEOGSM) +
     $   AM(I,2,GEOMAG)*AM(J,2,GEOGSM) + AM(I,3,GEOMAG)*AM(J,3,GEOGSM)
      ENDDO
      ENDDO
C
	ST(6) = AM(3,1,GSEMAG)       
	CT(6) = SQRT(1.-ST(6)*ST(6))      
	CX(6) = ASIND(ST(6))     

        AM(1,1,GSMSM) = CT(6)
        AM(1,2,GSMSM) = 0.
        AM(1,3,GSMSM) = -ST(6)
        AM(2,1,GSMSM) = 0.
        AM(2,2,GSMSM) = 1.
        AM(2,3,GSMSM) = 0.
        AM(3,1,GSMSM) = ST(6)
        AM(3,2,GSMSM) = 0.
        AM(3,3,GSMSM) = CT(6)
C         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GEOSM) = AM(I,1,GSMSM)*AM(1,J,GEOGSM) +
     $    AM(I,2,GSMSM)*AM(2,J,GEOGSM) +  AM(I,3,GSMSM)*AM(3,J,GEOGSM)
      ENDDO
      ENDDO
C         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,MAGSM) = AM(I,1,GSMSM)*AM(J,1,GSMMAG) +
     $   AM(I,2,GSMSM)*AM(J,2,GSMMAG) + AM(I,3,GSMSM)*AM(J,3,GSMMAG)
      ENDDO
      ENDDO
C
      CX(7)=ATAN2D( AM(2,1,11) , AM(1,1,11) )
      CX(8)=ASIND( AM(3,1,1) )
      CX(9)=ATAN2D( AM(2,1,1) , AM(1,1,1) )

      IF(IDBUG.NE.0)THEN
	  WRITE(IDBUG,*) 'Dipole tilt angle=',CX(6)
	  WRITE(IDBUG,*) 'Angle between sun and magnetic pole=',CX(7)
	  WRITE(IDBUG,*) 'Subsolar point latitude=',CX(8)
	  WRITE(IDBUG,*) 'Subsolar point longitude=',CX(9)

        DO K=1,11
         WRITE(IDBUG,1001) K
         DO I=1,3
           WRITE(IDBUG,1002) (AM(I,J,K),J=1,3)
         ENDDO
        ENDDO
 1001   FORMAT(' ROTATION MATRIX ',I2)
 1002   FORMAT(3F9.5)
      ENDIF

CNCAR      Mar 96: return the dipole tilt from this function call.
      GET_TILT = CX(6)
CNCAR

      RETURN
      END 
******************************************************************************
CNCAR      Feb 01:  Eliminate unused routines from translib.for: ROTATE,
C          ROTATEV, FROMCART, TOCART, MLT, MAGLONG, SUNLOC.  Remaining
C          are ADJUST and JULDAY
CNCAR
	SUBROUTINE ADJUST(ANGLE)
C	ADJUST AN ANGLE IN DEGREES TO BE IN RANGE OF 0 TO 360.
 10	CONTINUE
	IF(ANGLE.LT.0.)THEN
	  ANGLE=ANGLE+360.
	  GOTO 10
	ENDIF
 20	CONTINUE
 	IF(ANGLE.GE.360.)THEN
	  ANGLE=ANGLE-360.
	  GOTO 20
	ENDIF
	RETURN
	END
******************************************************************************
      FUNCTION JULDAY(MM,ID,IYYY)
      PARAMETER (IGREG=15+31*(10+12*1582))
      IF (IYYY.EQ.0) PAUSE 'There is no Year Zero.'
      IF (IYYY.LT.0) IYYY=IYYY+1
      IF (MM.GT.2) THEN
        JY=IYYY
        JM=MM+1
      ELSE
        JY=IYYY-1
        JM=MM+13
      ENDIF
      JULDAY=INT(365.25*JY)+INT(30.6001*JM)+ID+1720995
      IF (ID+31*(MM+12*IYYY).GE.IGREG) THEN
        JA=INT(0.01*JY)
        JULDAY=JULDAY+2-JA+INT(0.25*JA)
      ENDIF
      RETURN
      END

CNCAR      Routines added to work around non-ANSI trig functions which
C          input degrees instead of radians:  SIND, COSD, ASIND, ATAN2D

      FUNCTION SIND (DEG)
      PARAMETER ( D2R =  0.0174532925199432957692369076847 ,
     +            R2D = 57.2957795130823208767981548147)
      SIND = SIN (DEG * D2R)
      RETURN
      END

      FUNCTION COSD (DEG)
      PARAMETER ( D2R =  0.0174532925199432957692369076847 ,
     +            R2D = 57.2957795130823208767981548147)

      COSD = COS (DEG * D2R)
      RETURN
      END

      FUNCTION ASIND (RNUM)
      PARAMETER ( D2R =  0.0174532925199432957692369076847 ,
     +            R2D = 57.2957795130823208767981548147)
      ASIND = R2D * ASIN (RNUM)
      RETURN
      END

      FUNCTION ATAN2D (RNUM1,RNUM2)
      PARAMETER ( D2R =  0.0174532925199432957692369076847 ,
     +            R2D = 57.2957795130823208767981548147)
      ATAN2D = R2D * ATAN2 (RNUM1,RNUM2)
      RETURN
      END
CNCAR

