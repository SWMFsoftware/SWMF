      SUBROUTINE APEX(DATE,DLAT,DLON,ALT
     1 ,A,ALAT,ALON,BMAG,XMAG,YMAG,ZMAG,V)
C Calculate apex radius, latitude, longitude; and magnetic field and potential
C 940822 A. D. Richmond
C
C  Input:
C
C DATE = Year and fraction (1990.0 = 1990 January 1, 0 UT)
C DLAT = Geodetic latitude in degrees
C DLON = Geodetic longitude in degrees
C ALT = Altitude in km
C
C  Output:
C
C A = Apex height in RE (similar to L)
C ALAT = Apex latitude in degrees (negative in S. magnetic hemisphere)
C ALON = Apex longitude (Geomagnetic longitude of apex)
C BMAG = Magnitude of geomagnetic field in nT
C XMAG,YMAG,ZMAG = North, east, and downward geomagnetic field components in nT
C V    = Magnetic potential, in T.m
C
C     COMMON/DIPOLE/COLAT,ELON,VP,CTP,STP                           
C COLAT = Geocentric colatitude of geomagnetic dipole north pole (deg)
C ELON = East longitude of geomagnetic dipole north pole (deg)
C VP = Magnitude, in T.m, of dipole component of magnetic potential at 
C        geomagnetic pole and geocentric radius of 6371.2 km
C CTP,STP = cosine, sine of COLAT
C
      PARAMETER (RTOD=5.72957795130823E1,DTOR=1.745329251994330E-2)
      PARAMETER (RE=6371.2,REQ=6378.160)
      COMMON/DIPOLE/COLAT,ELON,VP,CTP,STP                           
      CALL COFRM(DATE)
      CALL DYPOL(CLATP,POLON,VPOL)
      COLAT = CLATP
      CTP = COS(CLATP*DTOR)
      STP = SQRT(1. - CTP*CTP)
      ELON = POLON
      VP = VPOL
      CALL LINAPX(DLAT,DLON,ALT
     2 ,A,ALAT,ALON,XMAG,YMAG,ZMAG,BMAG)
      XMAG = XMAG*1.E5
      YMAG = YMAG*1.E5
      ZMAG = ZMAG*1.E5
      BMAG = BMAG*1.E5
      CALL GD2CART(DLAT,DLON,ALT,X,Y,Z)
      CALL FELDG(3,X/RE,Y/RE,Z/RE,BX,BY,BZ,V)
      RETURN
      END
C
       SUBROUTINE LINAPX(GDLAT,GLON,ALT
     2 ,A,ALAT,ALON,XMAG,YMAG,ZMAG,F)
C***BEGIN PROLOGUE  LINAPX                                                      
C***DATE WRITTEN   731029   (YYMMDD)                                            
C***AUTHOR  CLARK, W., N.O.A.A. ERL LAB.                                        
C***REVISION DATE  880201   (YYMMDD) Harsh Anand Passi, NCAR
C***LATEST REVISION 940803 A. D. Richmond
C***PURPOSE  Transforms the geographic coordinates to apex coordinates.         
C***DESCRIPTION                                                                 
C     The method used is as follow:                                             
C       1. Calculates step size as a function of the geomagnetic                
C          dipole coordinates of the starting point.                            
C       2. Determine direction of trace                                         
C       3. Convert the geodetic coordinates of the starting point               
C          to the cartesian coordinates for tracing.                            
C       Loop:                                                                   
C       i)   Increment step count, if count > 200,                              
C            assume it is dipole field, call DIPAPX to                          
C            determine Apex coordinates else continue.                          
C       ii)  Get field components in X, Y and Z direction                       
C       iii) Call ITRACE to trace field line.                                   
C       iv)  Test if Apex passed call FNDAPX to determine Apex coordinates      
C            else loop:                                                         
C                                                                               
C   INPUT                                                                       
C     GDLAT  Latitude of starting point (Deg)                                   
C     GLON   Longitude (East=+) of starting point (Deg)                         
C     ALT    Ht. of starting point (Km)                                         
C                                                                               
C  OUTPUT
C          A  Apex radius (km)                                                  
C       ALAT  Apex Lat. (deg)                                                   
C       ALON  Apex Lon. (deg)                                                   
C       XMAG  North component of magnetic field at starting point
C       YMAG  East component of magnetic field at starting point
C       ZMAG  Down component of magnetic field at starting point
C          F  Magnetic field magnitude at starting point
C
C     COMMON Blocks Used                                                        
C     /APXIN/ YAPX(3,3)                                                         
C       YAPX    Matrix of cartesian coordinates (loaded columnwise)             
C               of the 3 points about APEX. Set in subprogram ITRACE.           
C                                                                               
C   /DIPOLE/COLAT,ELON,VP,CTP,STP                           
C     COLAT   Geographic colatitude of the north pole of the                    
C             earth-centered dipole (Deg).                                      
C     ELON    Geographic longitude of the north pole of the                     
C             earth-centered dipole (Deg).                                      
C     VP      Magnetic potential magnitude at geomagnetic pole (T.m)
C     CTP     cos(COLAT*DTOR)
C     STP     sin(COLAT*DTOR)
C                                                                               
C     /ITRA/ NSTP, Y(3), YOLD(3), SGN, DS                                       
C     NSTP      Step count. Incremented in subprogram LINAPX.                   
C     Y         Array containing current tracing point cartesian coordinates.   
C     YOLD      Array containing previous tracing point cartesian coordinates.
C     SGN       Determines direction of trace. Set in subprogram LINAPX         
C     DS        Step size (Km) Computed in subprogram LINAPX.                   
C                                                                               
C     /FLDCOMD/ BX, BY, BZ, BB                                                  
C     BX        X comp. of field vector at the current tracing point (Gauss)    
C     BY        Y comp. of field vector at the current tracing point (Gauss)    
C     BZ        Z comp. of field vector at the current tracing point (Gauss)    
C     BB        Magnitude of field vector at the current tracing point (Gauss)  
C                                                                               
C***REFERENCES  Stassinopoulos E. G. , Mead Gilbert D., X-841-72-17             
C                 (1971) GSFC, Greenbelt, Maryland                              
C                                                                               
C***ROUTINES CALLED  GD2CART,CONVRT,FELDG,
C                    ITRACE,DIPAPX,FNDAPX                                       
C                                                                               
C***END PROLOGUE  LINAPX                                                        
       PARAMETER ( MAXS = 200 )
      PARAMETER (RTOD=5.72957795130823E1,DTOR=1.745329251994330E-2)
      PARAMETER (RE=6371.2,REQ=6378.160)
       COMMON /FLDCOMD/ BX, BY, BZ, BB                                          
       COMMON /APXIN/ YAPX(3,3)                                                 
      COMMON/DIPOLE/COLAT,ELON,VP,CTP,STP                           
       COMMON /ITRA/ NSTP, Y(3), YP(3),  SGN, DS                                
C                                                                               
C ** DETERMINE STEP SIZE AS A FUNCTION OF GEOMAGNETIC DIPOLE COORDINATES
C    OF THE STARTING POINT                                                    
C **                                                                            
       CALL CONVRT(2,GDLAT,ALT,GCLAT,R)
C       SINGML = .98*SIN(GCLAT*DTOR) + .199*COS(GCLAT*DTOR)*              
       SINGML = CTP*SIN(GCLAT*DTOR) + STP*COS(GCLAT*DTOR)*                     
     +  COS((GLON-ELON)*DTOR)
       DS = .06*R/(1.-SINGML*SINGML) - 370.                                     
C  Limit DS to its value at 60 deg gm latitude.
       IF (DS .GT. 1186.) DS = 1186.                                            

C     Initialize YAPX array                                                     

       DO 4 J = 1,3                                                             
       DO 4 I =1,3                                                              
        YAPX(I,J) = 0.                                                          
4      CONTINUE                                                                 

C     Convert from geodetic to earth centered cartesian coordinates             

      CALL GD2CART(GDLAT,GLON,ALT,Y(1),Y(2),Y(3))
       NSTP = 0                                                                 
C                                                                               
C **GET FIELD COMPONENTS TO DETERMINE THE DIRECTION OF TRACING THE FIELD LINE
C                                                                               
       CALL FELDG(1,GDLAT,GLON,ALT,XMAG,YMAG,ZMAG,F)                   
       SGN = SIGN(1.,-ZMAG)                                   
C                                                                               
   10  CONTINUE                                                                 
C
C **   GET FIELD COMPONENTS AND TRACE ALONG FIELD LINE                          
C
       IENTRY=2                                                                 
       CALL FELDG(IENTRY,Y(1)/RE,Y(2)/RE,Y(3)/RE,BX,BY,BZ,BB)                   
C                                                                               
       NSTP = NSTP + 1                                                          
C                                                                               
C **   QUIT IF TOO MANY STEPS                                                   
C **                                                                            
      IF (NSTP .GE. MAXS) THEN                                                 
        RHO = SQRT(Y(1)*Y(1) + Y(2)*Y(2))
        CALL CONVRT(3,XLAT,HT,RHO,Y(3))
        XLON = RTOD*ATAN2(Y(2),Y(1))
        CALL FELDG(1,XLAT,XLON,HT,BNRTH,BEAST,BDOWN,BABS)                   
        CALL DIPAPX(XLAT,XLON,HT,BNRTH,BEAST,BDOWN,A,ALON)
        ALAT = -SGN*RTOD*ACOS(SQRT(1./A))
        RETURN                                                                
      END IF                                                                   
C                                                                               
C **   FIND NEXT POINT USING ADAMS ALGORITHM AFTER 7 POINTS                     
C **                                                                            
       CALL ITRACE( IAPX)                                                       
C                                                                               
C                                                                               
       GO TO (10, 20) IAPX                                                      
   20  CONTINUE                                                                 
C                                                                               
C **    MAXIMUM RADIUS JUST PASSED.  FIND APEX COORDS                           
C **                                                                            
       CALL FNDAPX(ALT,ZMAG,A,ALAT,ALON)
       RETURN                                                                   
       END

       SUBROUTINE ITRACE( IAPX)                                                 
      SAVE
C***BEGIN PROLOGUE  ITRACE                                                      
C***DATE WRITTEN   731029   (YYMMDD)                                            
C***AUTHOR  CLARK, W. N.O.A.A. ERL LAB.                                         
C***REVISION DATE  880201, H. Passi 
C***PURPOSE  Field line integration routine.                                    
C***DESCRIPTION                                                                 
C     It uses 4-point ADAMS formula after initialization.                       
C     First 7 iterations advance point by 3 steps.                              
C                                                                               
C   INPUT                                                                       
C            Passed in through the common blocks ITRA, FLDCOMD.                 
C            See the description below.                                         
C   OUTPUT                                                                      
C    IAPX    Flag set to 2 when APEX passed, otherwise set to 1.                
C                                                                               
C            Passed out through the common block APXIN.                         
C            See the description below.                                         
C                                                                               
C     COMMON Blocks Used                                                        
C     /APXIN/ YAPX(3,3)                                                         
C     YAPX      Matrix of cartesian coordinates (loaded columnwise)             
C               of the 3 points about APEX. Set in subprogram ITRACE.           
C                                                                               
C     /FLDCOMD/ BX, BY, BZ, BB                                                  
C     BX        X comp. of field vector at the current tracing point (Gauss)    
C     BY        Y comp. of field vector at the current tracing point (Gauss)    
C     BZ        Z comp. of field vector at the current tracing point (Gauss)    
C     BB        Magnitude of field vector at the current tracing point (Gauss)  
C                                                                               
C     /ITRA/ NSTP, Y(3), YOLD(3), SGN, DS                                       
C     NSTP      Step count for line integration.                                
C               Incremented in subprogram LINAPX.                               
C     Y         Array containing current tracing point cartesian coordinates.   
C     YOLD      Array containing previous tracing point cartesian coordinates.  
C     SGN       Determines direction of trace. Set in subprogram LINAPX         
C     DS        Integration step size (arc length Km)                           
C               Computed in subprogram LINAPX.                                  
C                                                                               
C***REFERENCES  reference 1                                                     
C                                                                               
C***ROUTINES CALLED  None                                                       
C***COMMON BLOCKS    APXIN,FLDCOMD,ITRA                                         
C***END PROLOGUE  ITRACE                                                        
C **                                                                            
       COMMON /ITRA/ NSTP, Y(3), YOLD(3), SGN, DS                               
       COMMON /FLDCOMD/ BX, BY, BZ, BB                                          
       COMMON /APXIN/ YAPX(3,3)                                                 
       DIMENSION  YP(3, 4)                                                      
C      Statement function                                                       
       RDUS(D,E,F) = SQRT( D**2 + E**2 + F**2 )                                 
C                                                                               
       IAPX = 1                                                                 
C      Field line is defined by the following differential equations            
C      in cartesian coordinates:                                                
       YP(1, 4) = SGN*BX/BB                                                     
       YP(2, 4) = SGN*BY/BB                                                     
       YP(3, 4) = SGN*BZ/BB                                                     
       IF (NSTP .GT. 7) GO TO 90                                                
C **  FIRST SEVEN STEPS USE THIS BLOCK                                          
       DO 80 I = 1, 3                                                           
         GO TO (10, 20, 30, 40, 50, 60, 70) NSTP                                
C                                                                               
   10   D2 = DS/2.                                                              
         D6 = DS/6.                                                             
         D12 = DS/12.                                                           
         D24 = DS/24.                                                           
         YP(I, 1) = YP(I, 4)                                                    
         YOLD(I) = Y(I)                                                         
         YAPX(I, 1) = Y(I)                                                      
         Y(I) = YOLD(I) + DS*YP(I, 1)                                           
         GO TO 80                                                               
C                                                                               
   20   YP(I, 2) = YP(I, 4)                                                     
         Y(I) = YOLD(I) + D2*(YP(I,2)+YP(I,1))                                  
         GO TO 80                                                               
C                                                                               
   30   Y(I) = YOLD(I) + D6*(2.*YP(I,4)+YP(I,2)+3.*YP(I,1))                     
         GO TO 80                                                               
C                                                                               
   40   YP(I, 2) = YP(I, 4)                                                     
         YAPX(I, 2) = Y(I)                                                      
         YOLD(I) = Y(I)                                                         
         Y(I) = YOLD(I) + D2*(3.*YP(I,2)-YP(I,1))                               
         GO TO 80                                                               
C                                                                               
  50   Y(I) = YOLD(I) + D12*(5.*YP(I,4)+8.*YP(I,2)-YP(I,1))                     
         GO TO 80                                                               
C                                                                               
   60   YP(I, 3) = YP(I, 4)                                                     
         YOLD(I) = Y(I)                                                         
         YAPX(I, 3) = Y(I)                                                      
         Y(I) = YOLD(I) + D12*(23.*YP(I,3)-16.*YP(I,2)+5.*YP(I,1))              
         GO TO 80                                                               
C                                                                               
   70  YAPX(I, 1) = YAPX(I, 2)                                                  
        YAPX(I, 2) = YAPX(I, 3)                                                 
        Y(I) = YOLD(I) + D24*(9.*YP(I,4)+19.*YP(I,3)-5.*YP(I,2)+YP(I,1))        
         YAPX(I, 3) = Y(I)                                                      
   80 CONTINUE                                                                  
C   **   SIGNAL IF APEX PASSED                                                  
       IF ( NSTP .EQ. 6 .OR. NSTP .EQ. 7) THEN                                  
        RC = RDUS( YAPX(1,3), YAPX(2,3), YAPX(3,3))                             
        RP = RDUS( YAPX(1,2), YAPX(2,2), YAPX(3,2))                             
        IF ( RC .LT. RP) IAPX=2                                                 
       ENDIF                                                                    
       RETURN                                                                   
C                                                                               
C **   STEPPING BLOCK FOR NSTEP .GT. 7                                          
   90 DO 110 I = 1, 3                                                           
         YAPX(I, 1) = YAPX(I, 2)                                                
         YAPX(I, 2) = Y(I)                                                      
         YOLD(I) = Y(I)                                                         
         TERM = 55.*YP(I, 4) - 59.*YP(I, 3) + 37.*YP(I, 2) - 9.*YP(I, 1)        
         Y(I) = YOLD(I) + D24*TERM                                              
         YAPX(I, 3) = Y(I)                                                      
         DO 100 J = 1, 3                                                        
           YP(I, J) = YP(I, J+1)                                                
  100   CONTINUE                                                                
  110 CONTINUE                                                                  
        RC =RDUS( Y(1), Y(2), Y(3))                                             
        RP = RDUS( YOLD(1), YOLD(2), YOLD(3))                                   
        IF ( RC .LT. RP) IAPX=2                                                 
       RETURN                                                                   
C                                                                               
       END                                                                      

       SUBROUTINE FNDAPX(ALT,ZMAG,A,ALAT,ALON)
C***BEGIN PROLOGUE  FNDAPX                                                      
C***DATE WRITTEN   731023   (YYMMDD)                                            
C***AUTHOR  CLARK, W., NOAA BOULDER                                             
C***REVISION DATE  940803, A. D. Richmond, NCAR 
C***PURPOSE  Finds apex coords once tracing has signalled that the apex         
C            has been passed.  
C***DESCRIPTION                                                                 
C                                                                               
C     It uses second degree interpolation routine, FINT, to find                
C     apex latitude and apex longtitude.                                        
C   INPUT                                                                       
C     ALT    Altitude of starting point
C     ZMAG   Downward component of geomagnetic field at starting point
C     NMAX   Order of IGRF coefficients being used
C     G      Array of coefficients from COFRM
C   OUTPUT                                                                      
C          A  Apex radius (km)                                                  
C       ALAT  Apex Lat. (deg)                                                   
C       ALON  Apex Lon. (deg)                                                   
C                                                                               
C***LONG DESCRIPTION                                                            
C                                                                               
C     COMMON Blocks Used                                                        
C     /APXIN/ YAPX(3,3)                                                         
C     YAPX      Matrix of cartesian coordinates (loaded columnwise)             
C               of the 3 points about APEX. Set in subprogram ITRACE.           
C                                                                               
C   /DIPOLE/COLAT,ELON,VP,CTP,STP                           
C     COLAT   Geocentric colatitude of the north pole of the                    
C             earth-centered dipole (Deg).                                      
C     ELON    Geographic longitude of the north pole of the                     
C             earth-centered dipole (Deg).                                      
C     CTP     cos(COLAT*DTOR)
C     STP     sin(COLAT*DTOR)
C                                                                               
C***ROUTINES CALLED  FINT                                                       
C***COMMON BLOCKS    APXIN,DIPOLE                                         
C***END PROLOGUE  FNDAPX                                                        
C                                                                               
      PARAMETER (RTOD=5.72957795130823E1,DTOR=1.745329251994330E-2)
      PARAMETER (RE=6371.2,REQ=6378.160)
       COMMON /APXIN/ YAPX(3,3)                                                 
      COMMON/DIPOLE/COLAT,ELON,VP,CTP,STP
       DIMENSION Z(3), HT(3), Y(3)

C       ****
C       ****     GET GEODETIC FIELD COMPONENTS
C       ****
      IENTRY = 1
      DO 2 I = 1,3
	RHO = SQRT(YAPX(1,I)**2+YAPX(2,I)**2)
	CALL CONVRT(3,GDLT,HT(I),RHO,YAPX(3,I))
        GDLN = RTOD*ATAN2(YAPX(2,I),YAPX(1,I))
	CALL FELDG(IENTRY,GDLT,GDLN,HT(I),X,YDUM,Z(I),F)
    2 CONTINUE
C     ****
C     ****     FIND CARTESIAN COORDINATES AT DIP EQUATOR BY INTERPOLATION
C     ****
      DO 3 I = 1,3
	CALL FINT(Z(1),Z(2),Z(3),YAPX(I,1),YAPX(I,2),YAPX(I,3),0.,Y(I))
    3 CONTINUE
C     ****
C     ****     FIND APEX HEIGHT BY INTERPOLATION
C     ****
      CALL FINT(Z(1),Z(2),Z(3),HT(1),HT(2),HT(3),0.,XINTER)
C Ensure that XINTER is not less than original starting altitude (ALT)
      XINTER = AMAX1(ALT,XINTER)
      A = (REQ+XINTER)/(REQ)
C     ****
C     ****     FIND APEX COORDINATES , GIVING ALAT SIGN OF DIP AT
C     ****       STARTING POINT.  ALON IS THE VALUE OF THE GEOMAGNETIC
C     ****       LONGITUDE AT THE APEX.
C     ****
      IF (A.LT.1.) THEN
	WRITE(6,20)
  20    FORMAT (' BOMBED! THIS MAKES A LESS THAN ONE')
	call stop_gitm("died in apex.f")
      ENDIF
      RASQ = RTOD*ACOS(SQRT(1./A))
      ALAT = SIGN(RASQ,ZMAG)
C Algorithm for ALON:
C   Use spherical coordinates.
C   Let GP be geographic pole.
C   Let GM be geomagnetic pole (colatitude COLAT, east longitude ELON).
C   Let XLON be longitude of apex.
C   Let TE be colatitude of apex.
C   Let ANG be longitude angle from GM to apex.
C   Let TP be colatitude of GM.
C   Let TF be arc length between GM and apex.
C   Let PA = ALON be geomagnetic longitude, i.e., Pi minus angle measured 
C     counterclockwise from arc GM-apex to arc GM-GP.
C   Then, using notation C=cos, S=sin, spherical-trigonometry formulas 
C     for the functions of the angles are as shown below.  Note: STFCPA,
C     STFSPA are sin(TF) times cos(PA), sin(PA), respectively.

      XLON = ATAN2(Y(2),Y(1))
      ANG = XLON-ELON*DTOR
      CANG = COS(ANG)
      SANG = SIN(ANG)
      R = SQRT(Y(1)**2+Y(2)**2+Y(3)**2)
      CTE = Y(3)/R
      STE = SQRT(1.-CTE*CTE)
      STFCPA = STE*CTP*CANG - CTE*STP
      STFSPA = SANG*STE
      ALON = ATAN2(STFSPA,STFCPA)*RTOD
      RETURN
      END                                                                      

      SUBROUTINE DIPAPX(GDLAT,GDLON,ALT,BNORTH,BEAST,BDOWN,A,ALON)
C Compute a, alon from local magnetic field using dipole and spherical approx.
C 940501 A. D. Richmond
C Input:
C   GDLAT  = geodetic latitude, degrees
C   GDLON  = geodetic longitude, degrees
C   ALT    = altitude, km
C   BNORTH = geodetic northward magnetic field component (any units)
C   BEAST  = eastward magnetic field component
C   BDOWN  = geodetic downward magnetic field component
C Output:
C   A      = apex radius, 1 + h_A/R_eq
C   ALON   = apex longitude, degrees
C
C Algorithm:
C   Use spherical coordinates.
C   Let GP be geographic pole.
C   Let GM be geomagnetic pole (colatitude COLAT, east longitude ELON).
C   Let G be point at GDLAT,GDLON.
C   Let E be point on sphere below apex of dipolar field line passing through G.
C   Let TD be dipole colatitude of point G, found by applying dipole formula
C     for dip angle to actual dip angle.
C   Let B be Pi plus local declination angle.  B is in the direction 
C     from G to E.
C   Let TG be colatitude of G.
C   Let ANG be longitude angle from GM to G.
C   Let TE be colatitude of E.
C   Let TP be colatitude of GM.
C   Let A be longitude angle from G to E.
C   Let APANG = A + ANG
C   Let PA be geomagnetic longitude, i.e., Pi minus angle measured 
C     counterclockwise from arc GM-E to arc GM-GP.
C   Let TF be arc length between GM and E.
C   Then, using notation C=cos, S=sin, COT=cot, spherical-trigonometry formulas 
C     for the functions of the angles are as shown below.  Note: STFCPA,
C     STFSPA are sin(TF) times cos(PA), sin(PA), respectively.

      PARAMETER (RTOD=5.72957795130823E1,DTOR=1.745329251994330E-2)
      PARAMETER (RE=6371.2,REQ=6378.160)
      COMMON/DIPOLE/COLAT,ELON,VP,CTP,STP
      BHOR = SQRT(BNORTH*BNORTH + BEAST*BEAST)
      IF (BHOR.EQ.0.) THEN
	ALON = 0.
	A = 1.E34
	RETURN
      ENDIF
      COTTD = BDOWN*.5/BHOR
      STD = 1./SQRT(1.+COTTD*COTTD)
      CTD = COTTD*STD
      SB = -BEAST/BHOR
      CB = -BNORTH/BHOR
      CTG = SIN(GDLAT*DTOR) 
      STG = COS(GDLAT*DTOR)
      ANG = (GDLON-ELON)*DTOR 
      SANG = SIN(ANG)
      CANG = COS(ANG)
      CTE = CTG*STD + STG*CTD*CB
      STE = SQRT(1. - CTE*CTE)
      SA = SB*CTD/STE
      CA = (STD*STG - CTD*CTG*CB)/STE
      CAPANG = CA*CANG - SA*SANG
      SAPANG = CA*SANG + SA*CANG
      STFCPA = STE*CTP*CAPANG - CTE*STP
      STFSPA = SAPANG*STE
      ALON = ATAN2(STFSPA,STFCPA)*RTOD
      R = ALT + RE
      HA = ALT + R*COTTD*COTTD
      A = 1. + HA/REQ
      RETURN
      END

       SUBROUTINE FINT(A1, A2, A3, A4, A5, A6, A7, RESULT)                      
C***PURPOSE  Second degree interpolation routine                                
C***REFER TO  FNDAPX                                                            
       RESULT = ((A2-A3)*(A7-A2)*(A7-A3)*A4-(A1-A3)*(A7-A1)*(A7-A3)*A5+         
     +  (A1-A2)*(A7-A1)*(A7-A2)*A6)/((A1-A2)*(A1-A3)*(A2-A3))
       RETURN                                                                   
       END                                                                      
      SUBROUTINE COFRM (DATE)
C          Assign DGRF/IGRF spherical harmonic coefficients, to degree and
C          order NMAX, for DATE, yyyy.fraction, into array G.  Coefficients
C          are interpolated from the DGRF dates through the current IGRF year.
C          Coefficients for a later DATE are extrapolated using the IGRF
C          initial value and the secular change coefficients.  A warning
C          message is issued if DATE is later than the last recommended
C          (5 yrs later than the IGRF).  An DATE input earlier than the
C          first DGRF (EPOCH(1)), results in a diagnostic and a STOP.
C
C          Output in COMMON /MAGCOF/ NMAX,GB(144),GV(144),ICHG
C             NMAX = Maximum order of spherical harmonic coefficients used
C             GB   = Coefficients for magnetic field calculation
C             GV   = Coefficients for magnetic potential calculation
C             ICHG = Flag indicating when GB,GV have been changed in COFRM
C
C          HISTORY (blame):
C          COFRM and FELDG originated 15 Apr 83 by Vincent B. Wickwar
C          (formerly at SRI. Int., currently at Utah State).  Although set
C          up to accomodate second order time derivitives, the IGRF
C          (GTT, HTT) have been zero.  The spherical harmonic coefficients
C          degree and order is defined by NMAX (currently 10).
C
C          Jun 86:  Updated coefficients adding DGRF 1980 & IGRF 1985, which
C          were obtained from Eos Vol. 7, No. 24.  Common block MAG was
C          replaced by MAGCOF, thus removing variables not used in subroutine
C          FELDG.  (Roy Barnes)
C
C          Apr 1992 (Barnes):  Added DGRF 1985 and IGRF 1990 as described
C          in EOS Vol 73 Number 16 Apr 21 1992.  Other changes were made so
C          future updates should:
C            (1) Increment NDGY;
C            (2) Append to EPOCH the next IGRF year;
C            (3) Append the next DGRF coefficients to G1DIM and H1DIM; and
C            (4) Replace the IGRF initial values (G0, GT) and rates of
C                change indices (H0, HT).
C
C          Apr 94 (Art Richmond): Computation of GV added, for finding
C          magnetic potential.
C
C          Aug 95 (Barnes):  Added DGRF for 1990 and IGRF for 1995, which were
C          obtained by anonymous ftp geomag.gsfc.nasa.gov (cd pub, mget table*)
C          as per instructions from Bob Langel (langel@geomag.gsfc.nasa.gov),
C          but, problems are to be reported to baldwin@geomag.gsfc.nasa.gov
 
C          Oct 95 (Barnes):  Correct error in IGRF-95 G 7 6 and H 8 7 (see
C          email in folder).  Also found bug whereby coefficients were not being
C          updated in FELDG when IENTY did not change. ICHG was added to flag
C          date changes.  Also, a vestigial switch (IS) was removed from COFRM:
C          It was always 0 and involved 3 branch if statements in the main
C          polynomial construction loop (now numbered 200).
 
      PARAMETER (RTOD=5.72957795130823E1,DTOR=1.745329251994330E-2)
      DOUBLE PRECISION F,F0
      COMMON /MAGCOF/ NMAX,G(144),GV(144),ICHG
 
      PARAMETER (NDGY=6 , NYT = NDGY+1 , NGH = 144*NDGY)
C          NDGY = Number of DGRF years of sets of coefficients
C          NYT  = Add one for the IGRF set (and point to it).
C          NGH  = Dimension of the equivalenced arrays
      DIMENSION GYR(12,12,NYT) , HYR(12,12,NYT), EPOCH(NYT) ,
     +          G1DIM(NGH)     , H1DIM(NGH) ,
     +          G0(12,12) , GT(12,12) , GTT(12,12) ,
     +          H0(12,12) , HT(12,12) ,  HTT(12,12)
      EQUIVALENCE (GYR(1,1,1),G1DIM(1))  , (HYR(1,1,1),H1DIM(1)) ,
     +            (GYR(1,1,NYT),G0(1,1)) , (HYR(1,1,NYT),H0(1,1))
 
      SAVE DATEL,GYR,HYR,EPOCH,G1DIM,H1DIM,G0,H0,GT,HT,GTT,HTT
C     SAVE DATEL,GYR,HYR,EPOCH,            G0,H0,GT,HT,GTT,HTT
      DATA DATEL /-999./
      DATA EPOCH /1965. , 1970. , 1975. , 1980. , 1985. , 1990. , 1995./
C          D_/Dtime2 coefficients are 0
      DATA GTT/144*0./,HTT/144*0./
C          DGRF g(n,m) for 1965:
C          The "column" corresponds to "n" and
C          the "line" corresponds to "m" as indicated in column 6;
C          e.g., for 1965 g(0,3) = 1297. or g(6,6) = -111.
      DATA (G1DIM(I),I=1,144) /0.,
     O  -30334.,-1662., 1297., 957.,-219., 45., 75., 13.,  8.,-2., 2*0.,
     1   -2119., 2997.,-2038., 804., 358., 61.,-57.,  5., 10.,-3., 3*0.,
     2           1594., 1292., 479., 254.,  8.,  4., -4.,  2., 2., 4*0.,
     3                   856.,-390., -31.,-228.,13.,-14.,-13.,-5., 5*0.,
     4                         252.,-157.,  4.,-26.,  0., 10.,-2., 6*0.,
     5                               -62.,  1., -6.,  8., -1., 4., 7*0.,
     6                                   -111., 13., -1., -1., 4., 8*0.,
     7                                           1., 11.,  5., 0., 9*0.,
     8                                                4.,  1., 2.,10*0.,
     9                                                    -2., 2.,11*0.,
     A                                                         0.,13*0./
C          DGRF g(n,m) for 1970:
       DATA (G1DIM(I),I=145,288)/ 0.,
     O  -30220.,-1781., 1287., 952.,-216., 43., 72., 14.,  8.,-3., 2*0.,
     1   -2068., 3000.,-2091., 800., 359., 64.,-57.,  6., 10.,-3., 3*0.,
     2           1611., 1278., 461., 262., 15.,  1., -2.,  2., 2., 4*0.,
     3                   838.,-395., -42.,-212.,14.,-13.,-12.,-5., 5*0.,
     4                         234.,-160.,  2.,-22., -3., 10.,-1., 6*0.,
     5                               -56.,  3., -2.,  5., -1., 6., 7*0.,
     6                                   -112., 13.,  0.,  0., 4., 8*0.,
     7                                          -2., 11.,  3., 1., 9*0.,
     8                                                3.,  1., 0.,10*0.,
     9                                                    -1., 3.,11*0.,
     A                                                        -1.,13*0./
C          DGRF g(n,m) for 1975:
      DATA (G1DIM(I),I=289,432)/ 0.,
     O  -30100.,-1902., 1276., 946.,-218., 45., 71., 14.,  7.,-3., 2*0.,
     1   -2013., 3010.,-2144., 791., 356., 66.,-56.,  6., 10.,-3., 3*0.,
     2           1632., 1260., 438., 264., 28.,  1., -1.,  2., 2., 4*0.,
     3                   830.,-405., -59.,-198.,16.,-12.,-12.,-5., 5*0.,
     4                         216.,-159.,  1.,-14., -8., 10.,-2., 6*0.,
     5                               -49.,  6.,  0.,  4., -1., 5., 7*0.,
     6                                   -111., 12.,  0., -1., 4., 8*0.,
     7                                          -5., 10.,  4., 1., 9*0.,
     8                                                1.,  1., 0.,10*0.,
     9                                                    -2., 3.,11*0.,
     A                                                        -1.,13*0./
C          DGRF g(n,m) for 1980:
      DATA (G1DIM(I),I=433,576)/ 0.,
     O  -29992.,-1997., 1281., 938.,-218., 48., 72., 18.,  5.,-4., 2*0.,
     1   -1956., 3027.,-2180., 782., 357., 66.,-59.,  6., 10.,-4., 3*0.,
     2           1663., 1251., 398., 261., 42.,  2.,  0.,  1., 2., 4*0.,
     3                   833.,-419., -74.,-192.,21.,-11.,-12.,-5., 5*0.,
     4                         199.,-162.,  4.,-12., -7.,  9.,-2., 6*0.,
     5                               -48., 14.,  1.,  4., -3., 5., 7*0.,
     6                                   -108., 11.,  3., -1., 3., 8*0.,
     7                                          -2.,  6.,  7., 1., 9*0.,
     8                                               -1.,  2., 2.,10*0.,
     9                                                    -5., 3.,11*0.,
     A                                                         0.,13*0./
C          DGRF g(n,m) for 1985:
      DATA (G1DIM(I),I=577,720)/ 0.,
     O  -29873.,-2072., 1296., 936.,-214., 53., 74., 21.,  5.,-4., 2*0.,
     1   -1905., 3044.,-2208., 780., 355., 65.,-62.,  6., 10.,-4., 3*0.,
     2           1687., 1247., 361., 253., 51.,  3.,  0.,  1., 3., 4*0.,
     3                   829.,-424., -93.,-185.,24.,-11.,-12.,-5., 5*0.,
     4                         170.,-164.,  4., -6., -9.,  9.,-2., 6*0.,
     5                               -46., 16.,  4.,  4., -3., 5., 7*0.,
     6                                   -102., 10.,  4., -1., 3., 8*0.,
     7                                           0.,  4.,  7., 1., 9*0.,
     8                                               -4.,  1., 2.,10*0.,
     9                                                    -5., 3.,11*0.,
     A                                                         0.,13*0./
C          DGRF g(n,m) for 1990:
      DATA (G1DIM(I),I=721,864)/ 0.,
     O -29775.,-2131., 1314., 939., -214., 61., 77., 23., 4., -3., 2*0.,
     1  -1848., 3059.,-2239., 780.,  353., 65.,-64.,  5., 9., -4., 3*0.,
     2          1686., 1248., 325.,  245., 59.,  2., -1., 1.,  2., 4*0.,
     3                  802.,-423., -109.,-178.,26.,-10.,-12.,-5., 5*0.,
     4                        141., -165.,  3., -1.,-12., 9., -2., 6*0.,
     5                               -36., 18.,  5.,  3.,-4.,  4., 7*0.,
     6                                    -96.,  9.,  4.,-2.,  3., 8*0.,
     7                                           0.,  2., 7.,  1., 9*0.,
     8                                                -6., 1., 3.,10*0.,
     9                                                    -6., 3.,11*0.,
     A                                                         0.,13*0./
C          DGRF h(n,m) for 1965:
      DATA (H1DIM(I),I=1,144)/13*0.,
     1    5776,-2016., -404., 148.,  19.,-11.,-61.,  7.,-22., 2., 3*0.,
     2           114.,  240.,-269., 128.,100.,-27.,-12., 15., 1., 4*0.,
     3                  -165., 13.,-126., 68., -2.,  9.,  7., 2., 5*0.,
     4                        -269.,-97.,-32.,  6.,-16., -4., 6., 6*0.,
     5                               81., -8., 26.,  4., -5.,-4., 7*0.,
     6                                    -7.,-23., 24., 10., 0., 8*0.,
     7                                        -12., -3., 10.,-2., 9*0.,
     8                                             -17., -4., 3.,10*0.,
     9                                                    1., 0.,11*0.,
     A                                                       -6.,13*0./
C          DGRF h(n,m) for 1970:
      DATA (H1DIM(I),I=145,288)/ 13*0.,
     1    5737.,-2047., -366., 167.,  26.,-12.,-70.,  7.,-21., 1., 3*0.,
     2             25.,  251.,-266., 139.,100.,-27.,-15., 16., 1., 4*0.,
     3                  -196.,  26.,-139., 72., -4.,  6.,  6., 3., 5*0.,
     4                        -279., -91.,-37.,  8.,-17., -4., 4., 6*0.,
     5                                83., -6., 23.,  6., -5.,-4., 7*0.,
     6                                      1.,-23., 21., 10., 0., 8*0.,
     7                                         -11., -6., 11.,-1., 9*0.,
     8                                              -16., -2., 3.,10*0.,
     9                                                     1., 1.,11*0.,
     A                                                        -4.,13*0./
C          DGRF h(n,m) for 1975:
      DATA (H1DIM(I),I=289,432)/ 13*0.,
     1    5675.,-2067., -333., 191.,  31.,-13.,-77.,  6.,-21., 1., 3*0.,
     2            -68.,  262.,-265., 148., 99.,-26.,-16., 16., 1., 4*0.,
     3                  -223.,  39.,-152., 75., -5.,  4.,  7., 3., 5*0.,
     4                        -288., -83.,-41., 10.,-19., -4., 4., 6*0.,
     5                                88., -4., 22.,  6., -5.,-4., 7*0.,
     6                                     11.,-23., 18., 10.,-1., 8*0.,
     7                                         -12.,-10., 11.,-1., 9*0.,
     8                                              -17., -3., 3.,10*0.,
     9                                                     1., 1.,11*0.,
     A                                                        -5.,13*0./
C          DGRF h(n,m) for 1980:
      DATA (H1DIM(I),I=433,576)/ 13*0.,
     1    5604.,-2129., -336., 212.,  46.,-15.,-82.,  7.,-21., 1., 3*0.,
     2           -200.,  271.,-257., 150., 93.,-27.,-18., 16., 0., 4*0.,
     3                  -252.,  53.,-151., 71., -5.,  4.,  9., 3., 5*0.,
     4                        -297., -78.,-43., 16.,-22., -5., 6., 6*0.,
     5                                92., -2., 18.,  9., -6.,-4., 7*0.,
     6                                     17.,-23., 16.,  9., 0., 8*0.,
     7                                         -10.,-13., 10.,-1., 9*0.,
     8                                              -15., -6., 4.,10*0.,
     9                                                     2., 0.,11*0.,
     A                                                        -6.,13*0./
C          DGRF h(n,m) for 1985:
      DATA (H1DIM(I),I=577,720)/ 13*0.,
     1    5500.,-2197., -310., 232.,  47.,-16.,-83.,  8.,-21., 1., 3*0.,
     2           -306.,  284.,-249., 150., 88.,-27.,-19., 15., 0., 4*0.,
     3                  -297.,  69.,-154., 69., -2.,  5.,  9., 3., 5*0.,
     4                        -297., -75.,-48., 20.,-23., -6., 6., 6*0.,
     5                                95., -1., 17., 11., -6.,-4., 7*0.,
     6                                     21.,-23., 14.,  9., 0., 8*0.,
     7                                          -7.,-15.,  9.,-1., 9*0.,
     8                                              -11., -7., 4.,10*0.,
     9                                                     2., 0.,11*0.,
     A                                                        -6.,13*0./
C          DGRF h(n,m) for 1990:
      DATA (H1DIM(I),I=721,864)/ 13*0.,
     1   5406., -2279., -284., 247.,  46.,-16.,-80., 10.,-20., 2., 3*0.,
     2           -373.,  293.,-240., 154., 82.,-26.,-19., 15., 1., 4*0.,
     3                  -352.,  84.,-153., 69.,  0.,  6., 11., 3., 5*0.,
     4                        -299., -69.,-52., 21.,-22., -7., 6., 6*0.,
     5                                97.,  1., 17., 12., -7.,-4., 7*0.,
     6                                     24.,-23., 12.,  9., 0., 8*0.,
     7                                          -4.,-16.,  8.,-2., 9*0.,
     8                                              -10., -7., 3.,10*0.,
     9                                                     2.,-1.,11*0.,
     A                                                        -6.,13*0./
 
 
C          Initial coefficients g0 (IGRF for 1995):
      DATA G0/0.,
     O -29682., -2197., 1329., 941.,-210., 66., 78., 24.,  4.,-3., 2*0.,
     1  -1789.,  3074.,-2268., 782., 352., 64.,-67.,  4.,  9.,-4., 3*0.,
     2           1685., 1249., 291., 237., 65.,  1., -1.,  1., 2., 4*0.,
     3                   769.,-421.,-122.,-172.,29., -9.,-12.,-5., 5*0.,
     4                         116.,-167.,  2.,  4.,-14.,  9.,-2., 6*0.,
     5                               -26., 17.,  8.,  4., -4., 4., 7*0.,
     6                                    -94., 10.,  5., -2., 3., 8*0.,
     7                                          -2.,  0.,  7., 1., 9*0.,
     8                                               -7.,  0., 3.,10*0.,
     9                                                    -6., 3.,11*0.,
     A                                                         0.,13*0./
C          D_/Dtime coefficients gt (IGRF for 1995-2000):
      DATA GT/0.,
     O      17.6, -13.2,  1.5,   .8,   .8,  .5, -.2,  .3, 0.0,0.0, 2*0.,
     1      13.0,   3.7, -6.4,   .9,   .1, -.4, -.8, -.2, 0.0,0.0, 3*0.,
     2              -.8,  -.2, -6.9, -1.5,  .6, -.6,  .1, 0.0,0.0, 4*0.,
     3                   -8.1,   .5, -2.0, 1.9,  .6,  .4, 0.0,0.0, 5*0.,
     4                         -4.6,  -.1, -.2, 1.2,-1.1, 0.0,0.0, 6*0.,
     5                                2.3, -.2,  .1,  .3, 0.0,0.0, 7*0.,
     6                                      .0,  .2,  .2, 0.0,0.0, 8*0.,
     7                                          -.6, -.9, 0.0,0.0, 9*0.,
     8                                               -.3, 0.0,0.0,10*0.,
     9                                                    0.0,0.0,11*0.,
     A                                                        0.0,13*0./
C          Initial coefficients h0 (IGRF for 1995-2000):
      DATA H0/13*0.,
     1    5318., -2356., -263., 262., 44.,-16.,-77., 12.,-19., 2., 3*0.,
     2            -425.,  302.,-232.,157., 77.,-25.,-20., 15., 1., 4*0.,
     3                   -406.,  98.,-152., 67., 3.,  7., 11., 3., 5*0.,
     4                         -301.,-64.,-57., 22.,-21., -7., 6., 6*0.,
     5                                99.,  4., 16., 12., -7.,-4., 7*0.,
     6                                     28.,-23., 10.,  9., 0., 8*0.,
     7                                          -3.,-17.,  7.,-2., 9*0.,
     8                                              -10., -8., 3.,10*0.,
     9                                                     1.,-1.,11*0.,
     A                                                        -6.,13*0./
C          D_/Dtime coefficients ht (IGRF for 1995-2000):
      DATA HT/13*0.,
     1        -18.3, -15.0,  4.1, 1.8, .2,  .3,  .8,  .4, 0.0,0.0, 3*0.,
     2                -8.8,  2.2, 1.2,1.2,-1.6,  .2, -.2, 0.0,0.0, 4*0.,
     3                     -12.1, 2.7, .3, -.2,  .6,  .2, 0.0,0.0, 5*0.,
     4                           -1.0,1.8, -.9, -.4,  .7, 0.0,0.0, 6*0.,
     5                                 .9, 1.0,  .0,  .0, 0.0,0.0, 7*0.,
     6                                     2.2, -.3,-1.2, 0.0,0.0, 8*0.,
     7                                           .0, -.7, 0.0,0.0, 9*0.,
     8                                               -.6, 0.0,0.0,10*0.,
     9                                                    0.0,0.0,11*0.,
     A                                                        0.0,13*0./
 
      nmax = 10
      ichg = -99999
C      DATA NMAX,ICHG /10,-99999/


C          Do not need to load new coefficients if date has not changed
      ICHG = 0
      IF (DATE .EQ. DATEL) GO TO 300
      DATEL = DATE
      ICHG = 1
 
C          Trap out of range date:
      IF (DATE .LT. EPOCH(1)) GO TO 9100
      IF (DATE .GT. EPOCH(NYT)+5.) WRITE(6,9200) DATE
 
      DO 100 I=1,NYT
      IF (DATE .LT. EPOCH(I)) GO TO 110
      IY = I
  100 CONTINUE
  110 CONTINUE
 
      TIME = DATE
      T = TIME-EPOCH(IY)
      G(1) = 0.0
      I = 2
      F0 = -1.0D-5
      DO 200 N=1,NMAX
      RN = REAL (N)
      F0 = F0*RN*(2.*RN-1.) / (4.*RN-2.)
      F  = 0.5*SQRT(2.0)*F0
      NN = N+1
      MM = 1
      IF (IY .EQ. NYT) THEN
C          Extrapolate coefficients
        G(I) = ((GTT(NN,MM)*T + GT(NN,MM))*T + G0(NN,MM)) * F0
      ELSE
C          Interpolate coefficients
        G(I) = (GYR(NN,MM,IY) +
     +         T/5.0 * (GYR(NN,MM,IY+1)-GYR(NN,MM,IY))) * F0
      ENDIF
      GV(I) = G(I)/FLOAT(NN)
      I = I+1
      DO 200 M=1,N
      TMP1 = REAL (N+M)
      TMP2 = REAL (N-M+1)
      F = F*SQRT(TMP2/TMP1)*TMP1/TMP2
      NN = N+1
      MM = M+1
      I1 = I+1
      IF (IY .EQ. NYT) THEN
C          Extrapolate coefficients
        G(I)  = ((GTT(NN,MM)*T + GT(NN,MM))*T + G0(NN,MM)) * F
        G(I1) = ((HTT(NN,MM)*T + HT(NN,MM))*T + H0(NN,MM)) * F
      ELSE
C          Interpolate coefficients
        G(I)  = (GYR(NN,MM,IY) +
     +          T/5.0 * (GYR(NN,MM,IY+1)-GYR(NN,MM,IY))) * F
        G(I1) = (HYR(NN,MM,IY) +
     +          T/5.0 * (HYR(NN,MM,IY+1)-HYR(NN,MM,IY))) * F
      ENDIF
      GV(I)  = G(I)  / REAL(NN)
      GV(I1) = G(I1) / REAL(NN)
  200 I = I+2
 
  300 CONTINUE

      RETURN
 
C          Error trap diagnostics:
 9100 WRITE(6,'('' '',/,
     +'' COFRM:  DATE'',F9.3,'' preceeds DGRF coefficients'',
     +                          '' presently coded.'')') DATE
      call stop_gitm("died in apex.f")
 9200 FORMAT(' ',/,
     +' COFRM:  DATE',F9.3,' is after the maximum',
     +                          ' recommended for extrapolation.')
      END
 
      SUBROUTINE DYPOL (COLAT,ELON,VP)
C          Computes parameters for dipole component of geomagnetic field.
C          COFRM must be called before calling DYPOL!
C          940504 A. D. Richmond
C
C          INPUT from COFRM through COMMON /MAGCOF/ NMAX,GB(144),GV(144),ICHG
C            NMAX = Maximum order of spherical harmonic coefficients used
C            GB   = Coefficients for magnetic field calculation
C            GV   = Coefficients for magnetic potential calculation
C            ICHG = Flag indicating when GB,GV have been changed
C
C          RETURNS:
C            COLAT = Geocentric colatitude of geomagnetic dipole north pole
C                    (deg)
C            ELON  = East longitude of geomagnetic dipole north pole (deg)
C            VP    = Magnitude, in T.m, of dipole component of magnetic
C                    potential at geomagnetic pole and geocentric radius
C                    of 6371.2 km
 
      PARAMETER (RTOD=5.72957795130823E1 , DTOR=1.745329251994330E-2,
     +           RE=6371.2, REQ=6378.160)
      COMMON /MAGCOF/ NMAX,G(144),GV(144),ICHG
 
C          Compute geographic colatitude and longitude of the north pole of
C          earth centered dipole
      GPL = SQRT( G(2  )**2+ G(3  )**2+ G(4  )**2)
      CTP = G(2  )/GPL
      STP = SQRT(1. - CTP*CTP)
      COLAT = (ACOS(CTP))*RTOD
      ELON = ATAN2( G(4  ), G(3  ))*RTOD
 
C          Compute magnitude of magnetic potential at pole, radius Re.
      VP = .2*GPL*RE
C          .2 = 2*(10**-4 T/gauss)*(1000 m/km) (2 comes through F0 in COFRM).
 
      RETURN
      END
 
      SUBROUTINE FELDG (IENTY,GLAT,GLON,ALT, BNRTH,BEAST,BDOWN,BABS)
C          Compute the DGRF/IGRF field components at the point GLAT,GLON,ALT.
C          COFRM must be called to establish coefficients for correct date
C          prior to calling FELDG.
C
C          IENTY is an input flag controlling the meaning and direction of the
C                remaining formal arguments:
C          IENTY = 1
C            INPUTS:
C              GLAT = Latitude of point (deg)
C              GLON = Longitude (east=+) of point (deg)
C              ALT  = Ht of point (km)
C            RETURNS:
C              BNRTH  north component of field vector (Gauss)
C              BEAST  east component of field vector  (Gauss)
C              BDOWN  downward component of field vector (Gauss)
C              BABS   magnitude of field vector (Gauss)
C
C          IENTY = 2
C            INPUTS:
C              GLAT = X coordinate (in units of earth radii 6371.2 km )
C              GLON = Y coordinate (in units of earth radii 6371.2 km )
C              ALT  = Z coordinate (in units of earth radii 6371.2 km )
C            RETURNS:
C              BNRTH = X component of field vector (Gauss)
C              BEAST = Y component of field vector (Gauss)
C              BDOWN = Z component of field vector (Gauss)
C              BABS  = Magnitude of field vector (Gauss)
C          IENTY = 3
C            INPUTS:
C              GLAT = X coordinate (in units of earth radii 6371.2 km )
C              GLON = Y coordinate (in units of earth radii 6371.2 km )
C              ALT  = Z coordinate (in units of earth radii 6371.2 km )
C            RETURNS:
C              BNRTH = Dummy variable
C              BEAST = Dummy variable
C              BDOWN = Dummy variable
C              BABS  = Magnetic potential (T.m)
C
C          INPUT from COFRM through COMMON /MAGCOF/ NMAX,GB(144),GV(144),ICHG
C            NMAX = Maximum order of spherical harmonic coefficients used
C            GB   = Coefficients for magnetic field calculation
C            GV   = Coefficients for magnetic potential calculation
C            ICHG = Flag indicating when GB,GV have been changed
C
C          HISTORY:
C          COFRM and FELDG originated 15 Apr 83 by Vincent B. Wickwar
C          (formerly at SRI. Int., currently at Utah State).
C
C          May 94 (A.D. Richmond): Added magnetic potential calculation
C
C          Oct 95 (Barnes): Added ICHG
 
      PARAMETER (RTOD=57.2957795130823, DTOR=0.01745329251994330,
     +           RE=6371.2, REQ=6378.160)
      COMMON /MAGCOF/ NMAX,GB(144),GV(144),ICHG
      DIMENSION XI(3),H(144),G(144)
      SAVE IENTYP, G
      DATA IENTYP/-10000/
 
      IF (IENTY .EQ. 1) THEN
        IS=1
        RLAT = GLAT*DTOR
        CT   = SIN(RLAT)
        ST   = COS(RLAT)
        RLON = GLON*DTOR
        CP   = COS(RLON)
        SP   = SIN(RLON)
        CALL GD2CART (GLAT,GLON,ALT,XXX,YYY,ZZZ)
        XXX = XXX/RE
        YYY = YYY/RE
        ZZZ = ZZZ/RE
      ELSE
        IS   = 2
        XXX  = GLAT
        YYY  = GLON
        ZZZ  = ALT
      ENDIF
      RQ    = 1./(XXX**2+YYY**2+ZZZ**2)
      XI(1) = XXX*RQ
      XI(2) = YYY*RQ
      XI(3) = ZZZ*RQ
      IHMAX = NMAX*NMAX+1
      LAST  = IHMAX+NMAX+NMAX
      IMAX  = NMAX+NMAX-1
 
      IF (IENTY .NE. IENTYP .OR. ICHG .EQ. 1) THEN
        IENTYP = IENTY
	ICHG = 0
        IF (IENTY .NE. 3) THEN
          DO 12 I=1,LAST
   12     G(I) = GB(I)
        ELSE
          DO 18 I=1,LAST
   18     G(I) = GV(I)
        ENDIF
      ENDIF
 
      DO 8 I=IHMAX,LAST
    8 H(I) = G(I)
      MK = 3
      IF ( IMAX .EQ. 1) MK=1
      DO 6 K=1,MK,2
      I  = IMAX
      IH = IHMAX
    1 IL = IH-I
      F = 2./FLOAT(I-K+2)
      X = XI(1)*F
      Y = XI(2)*F
      Z = XI(3)*(F+F)
      I = I-2
      IF(I-1)5,4,2
    2 DO 3 M=3,I,2
      IHM = IH+M
      ILM = IL+M
      H(ILM+1) = G(ILM+1)+ Z*H(IHM+1) + X*(H(IHM+3)-H(IHM-1))
     +                                        -Y*(H(IHM+2)+H(IHM-2))
    3 H(ILM)   = G(ILM)  + Z*H(IHM)   + X*(H(IHM+2)-H(IHM-2))
     +                                        +Y*(H(IHM+3)+H(IHM-1))
    4 H(IL+2) = G(IL+2) + Z*H(IH+2) + X*H(IH+4) - Y*(H(IH+3)+H(IH))
      H(IL+1) = G(IL+1) + Z*H(IH+1) + Y*H(IH+4) + X*(H(IH+3)-H(IH))
    5 H(IL)   = G(IL)   + Z*H(IH)   + 2.*(X*H(IH+1)+Y*H(IH+2))
      IH = IL
      IF (I-K) 6,1,1
    6 CONTINUE
 
      S = .5*H(1)+2.*(H(2)*XI(3)+H(3)*XI(1)+H(4)*XI(2))
      T = (RQ+RQ)*SQRT(RQ)
      BXXX = T*(H(3)-S*XXX)
      BYYY = T*(H(4)-S*YYY)
      BZZZ = T*(H(2)-S*ZZZ)
      BABS = SQRT(BXXX**2+BYYY**2+BZZZ**2)
      IF (IS .EQ. 1) THEN            ! (convert back to geodetic)
        BEAST = BYYY*CP-BXXX*SP
        BRHO  = BYYY*SP+BXXX*CP
        BNRTH =  BZZZ*ST-BRHO*CT
        BDOWN = -BZZZ*CT-BRHO*ST
      ELSEIF (IS .EQ. 2) THEN        ! (leave in earth centered cartesian)
        BNRTH = BXXX
        BEAST = BYYY
        BDOWN = BZZZ
      ENDIF
 
C          Magnetic potential computation makes use of the fact that the
C          calculation of V is identical to that for r*Br, if coefficients
C          in the latter calculation have been divided by (n+1) (coefficients
C          GV).  Factor .1 converts km to m and gauss to tesla.
      IF (IENTY.EQ.3) BABS = (BXXX*XXX + BYYY*YYY + BZZZ*ZZZ)*RE*.1
 
      RETURN
      END
 
      SUBROUTINE GD2CART (GDLAT,GLON,ALT,X,Y,Z)
C          Convert geodetic to cartesian coordinates by calling CONVRT
C          940503 A. D. Richmond
      PARAMETER (RTOD=57.2957795130823E1, DTOR=0.01745329251994330)
      CALL CONVRT (1,GDLAT,ALT,RHO,Z)
      ANG = GLON*DTOR
      X = RHO*COS(ANG)
      Y = RHO*SIN(ANG)
      RETURN
      END
 
      SUBROUTINE CONVRT (I,GDLAT,ALT,X1,X2)
C          Convert space point from geodetic to geocentric or vice versa.
C
C          I is an input flag controlling the meaning and direction of the
C            remaining formal arguments:
C
C          I = 1  (convert from geodetic to cylindrical)
C            INPUTS:
C              GDLAT = Geodetic latitude (deg)
C              ALT   = Altitude above reference ellipsoid (km)
C            RETURNS:
C              X1    = Distance from Earth's rotation axis (km)
C              X2    = Distance above (north of) Earth's equatorial plane (km)
C
C          I = 2  (convert from geodetic to geocentric spherical)
C            INPUTS:
C              GDLAT = Geodetic latitude (deg)
C              ALT   = Altitude above reference ellipsoid (km)
C            RETURNS:
C              X1    = Geocentric latitude (deg)
C              X2    = Geocentric distance (km)
C
C          I = 3  (convert from cylindrical to geodetic)
C            INPUTS:
C              X1    = Distance from Earth's rotation axis (km)
C              X2    = Distance from Earth's equatorial plane (km)
C            RETURNS:
C              GDLAT = Geodetic latitude (deg)
C              ALT   = Altitude above reference ellipsoid (km)
C
C          I = 4  (convert from geocentric spherical to geodetic)
C            INPUTS:
C              X1    = Geocentric latitude (deg)
C              X2    = Geocentric distance (km)
C            RETURNS:
C              GDLAT = Geodetic latitude (deg)
C              ALT   = Altitude above reference ellipsoid (km)
C
C
C          HISTORY:
C          940503 (A. D. Richmond):  Based on a routine originally written
C          by V. B. Wickwar.
C
C          REFERENCE:  ASTRON. J. VOL. 66, p. 15-16, 1961.
 
      PARAMETER (RTOD=57.2957795130823, DTOR=0.01745329251994330 ,
     +  RE=6371.2 , REQ=6378.160 , FLTNVRS=298.25 ,
     +  E2=(2.-1./FLTNVRS)/FLTNVRS , E4=E2*E2 , E6=E4*E2 , E8=E4*E4 ,
     +  OME2REQ = (1.-E2)*REQ ,
     +     A21 =     (512.*E2 + 128.*E4 + 60.*E6 + 35.*E8)/1024. ,
     +     A22 =     (                        E6 +     E8)/  32. ,
     +     A23 = -3.*(                     4.*E6 +  3.*E8)/ 256. ,
     +     A41 =    -(           64.*E4 + 48.*E6 + 35.*E8)/1024. ,
     +     A42 =     (            4.*E4 +  2.*E6 +     E8)/  16. ,
     +     A43 =                                   15.*E8 / 256. ,
     +     A44 =                                      -E8 /  16. ,
     +     A61 =  3.*(                     4.*E6 +  5.*E8)/1024. ,
     +     A62 = -3.*(                        E6 +     E8)/  32. ,
     +     A63 = 35.*(                     4.*E6 +  3.*E8)/ 768. ,
     +     A81 =                                   -5.*E8 /2048. ,
     +     A82 =                                   64.*E8 /2048. ,
     +     A83 =                                 -252.*E8 /2048. ,
     +     A84 =                                  320.*E8 /2048. )
C          E2 = Square of eccentricity
 
      IF (I .GE. 3) GO TO 300
 
C          Geodetic to geocentric
 
C          Compute RHO,Z
      SINLAT = SIN(GDLAT*DTOR)
      COSLAT = SQRT(1.-SINLAT*SINLAT)
      D      = SQRT(1.-E2*SINLAT*SINLAT)
      Z      = (ALT+OME2REQ/D)*SINLAT
      RHO    = (ALT+REQ/D)*COSLAT
      X1 = RHO
      X2 = Z
      IF (I .EQ. 1) RETURN
 
C          Compute GCLAT,RKM
      RKM   = SQRT(Z*Z + RHO*RHO)
      GCLAT = RTOD*ATAN2(Z,RHO)
      X1 = GCLAT
      X2 = RKM
      RETURN
 
C          Geocentric to geodetic
  300 IF (I .EQ. 3) THEN
         RHO = X1
         Z = X2
         RKM = SQRT(Z*Z+RHO*RHO)
         SCL = Z/RKM
         GCLAT = ASIN(SCL)*RTOD
      ELSEIF (I .EQ. 4) THEN
         GCLAT = X1
         RKM = X2
         SCL = SIN(GCLAT*DTOR)
      ELSE
         RETURN
      ENDIF
 
      RI = REQ/RKM
      A2 = RI*(A21+RI*(A22+RI* A23))
      A4 = RI*(A41+RI*(A42+RI*(A43+RI*A44)))
      A6 = RI*(A61+RI*(A62+RI* A63))
      A8 = RI*(A81+RI*(A82+RI*(A83+RI*A84)))
      CCL = SQRT(1.-SCL*SCL)
      S2CL = 2.*SCL*CCL
      C2CL = 2.*CCL*CCL-1.
      S4CL = 2.*S2CL*C2CL
      C4CL = 2.*C2CL*C2CL-1.
      S8CL = 2.*S4CL*C4CL
      S6CL = S2CL*C4CL+C2CL*S4CL
      DLTCL = S2CL*A2+S4CL*A4+S6CL*A6+S8CL*A8
      GDLAT = DLTCL*RTOD+GCLAT
      SGL = SIN(GDLAT*DTOR)
      ALT = RKM*COS(DLTCL)-REQ*SQRT(1.-E2*SGL*SGL)
      RETURN
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine magloctm(alon,sbsllat,sbsllon,clatp,polon,mlt)
C  Computes magnetic local time from magnetic longitude, subsolar coordinates,
C   and geomagnetic pole coordinates.
C  950302 A. D. Richmond, NCAR
C  Algorithm:  MLT is calculated from the difference of the apex longitude,
C   alon, and the geomagnetic dipole longitude of the subsolar point.
C
C   Inputs:
C    alon = apex magnetic longitude of the point (deg)
C    sbsllat = geographic latitude of subsolar point (degrees)
C    sbsllon = geographic longitude of subsolar point (degrees)
C    clatp = Geocentric colatitude of geomagnetic dipole north pole (deg)
C    polon = East longitude of geomagnetic dipole north pole (deg)
C
C   Output:
C    mlt (real) = magnetic local time for the apex longitude alon (hours)
C
C To go from mlt to alon (see comments following Entry mlt2alon for definition 
C  of variables), use:
C
C     call mlt2alon(mlt,sbsllat,sbsllon,clatp,polon,alon)
C
C  NOTE: If the calling routine uses subroutine magloctm in conjunction with 
C   file magfld.f (which is used by subroutine APEX), then clatp and polon can 
C   be found by invoking
C
C     call DYPOL(clatp,polon,vp)
C
C   where vp is an unneeded variable.  (Note that subroutine COFRM must have
C   been called before DYPOL, in order to set up the coefficients for the
C   desired epoch.)  Alternatively, if subroutine apxntrp is
C   used to get alon from previously computed arrays, then
C   clatp and polon can be obtained for use in magloctm by adding
C
C     common/dipol/clatp,polon,dum1,dum2,dum3
C
C   to the calling routine (where dum1,dum2,dum3 are unneeded dummy variables). 
C
	real mlt
	call SOLGMLON(sbsllat,sbsllon,clatp,polon,smlon)
        mlt = (alon - smlon)/15.0 + 12.0
	if (mlt.ge.24.0) mlt = mlt - 24.0
	if (mlt.lt.0.) mlt = mlt + 24.0
	return
C
      Entry mlt2alon(xmlt,sbsllat,sbsllon,clatp,polon,alonx) 
C
C   Inputs:
C    xmlt (real) = magnetic local time for the apex longitude alon (hours, 
C                 0. to 24.)
C    sbsllat = geographic latitude of subsolar point (degrees)
C    sbsllon = geographic longitude of subsolar point (degrees)
C    clatp = Geocentric colatitude of geomagnetic dipole north pole (deg)
C    polon = East longitude of geomagnetic dipole north pole (deg)
C
C   Output:
C    alonx = apex magnetic longitude of the point (deg, -180. to 180.)
C
	call SOLGMLON(sbsllat,sbsllon,clatp,polon,smlon)
	alonx = 15.*(xmlt - 12.0) + smlon
	if (alonx.gt.180.) alonx = alonx - 360.0
	if (alonx.le.-180.) alonx = alonx + 360.0
	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE SOLGMLON(XLAT,XLON,COLAT,ELON,MLON)
C Computes geomagnetic longitude of the point with geocentric spherical
C  latitude and longitude of XLAT and XLON, respectively.
C 940719 A. D. Richmond, NCAR
C Inputs:
C   XLAT = geocentric spherical latitude (deg)
C   XLON = geocentric spherical longitude (deg)
C   COLAT = Geocentric colatitude of geomagnetic dipole north pole (deg)
C   ELON = East longitude of geomagnetic dipole north pole (deg)
C Output:
C   MLON = Geomagnetic dipole longitude of the point (deg, -180. to 180.)

      REAL MLON
      PARAMETER (RTOD=5.72957795130823E1,DTOR=1.745329251994330E-2)

C Algorithm:
C   Use spherical coordinates.
C   Let GP be geographic pole.
C   Let GM be geomagnetic pole (colatitude COLAT, east longitude ELON).
C   Let XLON be longitude of point P.
C   Let TE be colatitude of point P.
C   Let ANG be longitude angle from GM to P.
C   Let TP be colatitude of GM.
C   Let TF be arc length between GM and P.
C   Let PA = MLON be geomagnetic longitude, i.e., Pi minus angle measured
C     counterclockwise from arc GM-P to arc GM-GP.
C   Then, using notation C=cos, S=sin, spherical-trigonometry formulas
C     for the functions of the angles are as shown below.  Note: STFCPA,
C     STFSPA are sin(TF) times cos(PA), sin(PA), respectively.

      CTP = COS(COLAT*DTOR)
      STP = SQRT(1. - CTP*CTP)
      ANG = (XLON-ELON)*DTOR
      CANG = COS(ANG)
      SANG = SIN(ANG)
      CTE = SIN(XLAT*DTOR)
      STE = SQRT(1.-CTE*CTE)
      STFCPA = STE*CTP*CANG - CTE*STP
      STFSPA = SANG*STE
      MLON = ATAN2(STFSPA,STFCPA)*RTOD
      RETURN
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE SUBSOLR (IYR,IDAY,IHR,IMN,SEC,SBSLLAT,SBSLLON)
C          Find subsolar geographic latitude and longitude from date and time.
C	   Based on formulas in Astronomical Almanac for the year 1996, p. C24.
C	    (U.S. Government Printing Office, 1994).
C	   Usable for years 1601-2100, inclusive.  According to the Almanac, 
C	    results are good to at least 0.01 degree latitude and 0.025 degree 
C	    longitude between years 1950 and 2050.  Accuracy for other years 
C	    has not been tested.  Every day is assumed to have exactly
C	    86400 seconds; thus leap seconds that sometimes occur on December
C	    31 are ignored:  their effect is below the accuracy threshold of
C	    the algorithm.
C          961026 A. D. Richmond, NCAR
C
C          INPUTS:
C            IYR  = year (e.g., 1994).  1600 < IYR < 2101
C            IDAY = day of the year (1 = January 1; 365 [366 in leap year] = 
C		December 31).  No constraint on permitted values, i.e.,
C		may be negative or greater than 366.
C            IHR  = hour (e.g., 13 for 13:49). 
C		No constraint on permitted values.
C            IMN  = minute (e.g., 49 for 13:49).
C		No constraint on permitted values.
C            SEC  = second and fraction after the hour/minute.
C		No constraint on permitted values.
C          OUTPUT:
C            SBSLLAT = geographic latitude of subsolar point (degrees)
C            SBSLLON = geographic longitude of subsolar point (degrees,
C                      between -180 and +180)
      PARAMETER (D2R=0.0174532925199432957692369076847 ,
     +           R2D=57.2957795130823208767981548147)
      PARAMETER (MSGUN=6)
      INTEGER IYR,YR,IDAY,IHR,IMN,NLEAP,NCENT,NROT
      REAL SEC,SBSLLAT,SBSLLON,L0,G0,DF,LF,GF,L,G,LAMBDA,EPSILON,N
     1   ,GRAD,LAMRAD,SINLAM,EPSRAD,DELTA,UT,ETDEG
C
C Number of years from 2000 to IYR (negative if IYR < 2000):
      YR = IYR - 2000
C
C NLEAP (final) = number of leap days from (2000 January 1) to (IYR January 1)
C                 (negative if IYR is before 1997)
      NLEAP = (IYR-1601)/4
      NLEAP = NLEAP - 99
      IF (IYR.LE.1900) THEN
	IF (IYR.LE.1600) THEN
	 WRITE(MSGUN,*) 'SUBSOLR INVALID BEFORE 1601: INPUT YEAR = ',IYR       
	 call stop_gitm("died in apex.f")
	ENDIF
	NCENT = (IYR-1601)/100
	NCENT = 3 - NCENT 
	NLEAP = NLEAP + NCENT
      ENDIF
      IF (IYR.GE.2101) THEN
	WRITE(MSGUN,*) 'SUBSOLR INVALID AFTER 2100:  INPUT YEAR = ',IYR
	call stop_gitm("died in apex.f")
      ENDIF
C
C L0 = Mean longitude of Sun at 12 UT on January 1 of IYR:
C     L0 = 280.461 + .9856474*(365*(YR-NLEAP) + 366*NLEAP) 
C	   - (ARBITRARY INTEGER)*360.
C        = 280.461 + .9856474*(365*(YR-4*NLEAP) + (366+365*3)*NLEAP) 
C	   - (ARBITRARY INTEGER)*360.
C        = (280.461 - 360.) + (.9856474*365 - 360.)*(YR-4*NLEAP) 
C	   + (.9856474*(366+365*3) - 4*360.)*NLEAP,
C  where ARBITRARY INTEGER = YR+1.  This gives:
      L0 = -79.549 + (-.238699*(YR-4*NLEAP) + 3.08514E-2*NLEAP)
C
C G0 = Mean anomaly at 12 UT on January 1 of IYR:
C     G0 = 357.528 + .9856003*(365*(YR-NLEAP) + 366*NLEAP) 
C	   - (ARBITRARY INTEGER)*360.
C        = 357.528 + .9856003*(365*(YR-4*NLEAP) + (366+365*3)*NLEAP) 
C	   - (ARBITRARY INTEGER)*360.
C        = (357.528 - 360.) + (.9856003*365 - 360.)*(YR-4*NLEAP) 
C	   + (.9856003*(366+365*3) - 4*360.)*NLEAP,
C  where ARBITRARY INTEGER = YR+1.  This gives:
      G0 = -2.472 + (-.2558905*(YR-4*NLEAP) - 3.79617E-2*NLEAP)
C
C Universal time in seconds:
      UT = FLOAT(IHR*3600 + IMN*60) + SEC
C
C Days (including fraction) since 12 UT on January 1 of IYR:
      DF = (UT/86400. - 1.5) + IDAY
C
C Addition to Mean longitude of Sun since January 1 of IYR:
      LF = .9856474*DF
C
C Addition to Mean anomaly since January 1 of IYR:
      GF = .9856003*DF
C
C Mean longitude of Sun:
      L = L0 + LF
C
C Mean anomaly:
      G = G0 + GF
      GRAD = G*D2R
C
C Ecliptic longitude:
      LAMBDA = L + 1.915*SIN(GRAD) + .020*SIN(2.*GRAD)
      LAMRAD = LAMBDA*D2R
      SINLAM = SIN(LAMRAD)
C
C Days (including fraction) since 12 UT on January 1 of 2000:
      N = DF + FLOAT(365*YR + NLEAP)
C
C Obliquity of ecliptic:
      EPSILON = 23.439 - 4.E-7*N
      EPSRAD = EPSILON*D2R
C
C Right ascension:
      ALPHA = ATAN2(COS(EPSRAD)*SINLAM,COS(LAMRAD))*R2D
C
C Declination:
      DELTA = ASIN(SIN(EPSRAD)*SINLAM)*R2D
C
C Subsolar latitude:
      SBSLLAT = DELTA
C
C Equation of time (degrees):
      ETDEG = L - ALPHA
      NROT = NINT(ETDEG/360.)
      ETDEG = ETDEG - FLOAT(360*NROT)
C
C Apparent time (degrees):
      APTIME = UT/240. + ETDEG
C          Earth rotates one degree every 240 s.
C
C Subsolar longitude:
      SBSLLON = 180. - APTIME
      NROT = NINT(SBSLLON/360.)
      SBSLLON = SBSLLON - FLOAT(360*NROT)
C
      RETURN
      END
