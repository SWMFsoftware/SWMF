************************ Copyright 1996,2001 Dan Weimer/MRC ***********************
*
* Subroutines to calculate the electric potentials from the "Weimer 2K" model of
* the polar cap ionospheric electric potentials described in the publication: 
* Weimer, D. R., An improved model of ionospheric electric potentials including
* substorm perturbations and application to the Geospace Environment Modeling
* November 24, 1996 event, Journal of Geophysical Research, Vol. 106, p. 407, 2001.
*
* To use, first call procedure SETMODEL01 with the specified input parameters:
*   angle: IMF Y-Z clock angle in degrees, 0=northward, 180=southward
*   Bt: Magnitude of IMF in Y-Z plane in nT
*   Tilt: dipole tilt angle in degrees.
*   SWVel: solar wind velocity in km/sec
*   SWDen: solar wind density in #/cc
*   ALindex: (optional) AL index in nT
*
* The function EPOTVAL01(gLAT,gMLT) can then be used repeatively to get the
* electric potential in kV at the desired location.
* Input coordinates assume use of 'altitude adjusted' corrected geomagnetic
* coordinates for R=1, also refered to as AACGM0.
*
* The function BOUNDARYLAT(gMLT) can be used to get the latitude of the boundary
*   where the potential goes to zero.  This boundary is a function of MLT, and
*   varies with the SETMODEL01 parameters.  The potential is zero everywhere below
*   this boundary.
*
* Two data files are provided:
*	'w2klittle.dat' for LITTLE_ENDIAN machines.
*	'w2kbig.dat'    for    BIG_ENDIAN machines.
* You must copy or rename the correct one to the file 'w2k.dat'
*
* This code is protected by copyright and is distributed
* for research or educational use only.
* Commerical use without written permission from Dan Weimer/MRC is prohibited.

CNCAR      Revisions for use at NCAR:
C            (1) Change behavior at minimum magnetic latitude.  When approaching
C                the model equatorial edge (which varies with MLT) the electric
C                potential returned used to go to zero discontinuously; although
C                intended as a flag, it created artificial gradients in the
C                electric field calculation.  Now the potential returned is
C                constant (that of the minimum latitude) for any latitude at
C                or equatorward of the minimum.
C            (2) Accomodate running simultaneously 1996 and 2001 versions.  To
C                avoid name collisions this required: (i) revising names (e.g.,
C                adding '01') for differing subprograms, and (ii) relocating
C                common routines into another file (weicom.f).
C            (3) Pass the coefficients file name and unit number into READCOEF01
C                rather than using hard coded values.
C            (4) Add wrapper subroutines for non-ANSI trig functions which
C                input angles in degrees.
C            (5) Add electric field routine (GECMP01) to deterine the electric
C                potential gradient.
C            (6) Add wrapper routine (WEIEPOT01) for use with AMIE; this is a
C                substitute for calling SETMODEL01 and EPOTVAL01.
C            (7) Remove blanks in some statements to fit in 72 columns to
C                adhere to ANSI Fortran 77.
CNCAR      NCAR changes are delimited by "CNCAR"
*
************************ Copyright 1996,2001 Dan Weimer/MRC ***********************
	SUBROUTINE DLEGENDRE(x,lmax,mmax,Plm,dPlm,derivative)
* compute Double Precision Associate Legendre Function P_l^m(x)
* for all l up to lmax and all m up to mmax.
* Returns results in array Plm.
* If the LOGICAL flag "derivative" is set, then the first derivatives are also
* computed, and put into array dPlm.
* The recursion formulas keep a count of the exponents of the factor SQRT(1-x^2)
*  in both the numerator and denominator, which may cancel each other in the
*  final evaluation, particularly with the derivatives.  This prevents infinities
*  at x=-1 or +1. 
* If X is out of range ( abs(x)>1 ) then value is returns as if x=1.
	DOUBLE PRECISION x,xx,Plm(0:10,0:10),P(0:10,0:10,0:1),fact,sfact
	DOUBLE PRECISION dPlm(0:10,0:10),dP(0:10,0:10,0:2),anum,term
	LOGICAL derivative

	DO l=0,lmax
	    DO m=0,mmax
		  Plm(l,m)=0.D0
		  P(l,m,0)=0.D0
		  P(l,m,1)=0.D0
	    ENDDO
	ENDDO
	IF(lmax .LT. 0 .OR. mmax .LT. 0 .OR. mmax .GT. lmax )THEN
	  Print *,'Bad arguments to DLegendre'
	  RETURN
	ENDIF

* Copy x to xx, and make sure it is in range of -1. to +1.
	xx=MIN(x,1.D0)
	xx=MAX(xx,-1.D0)

	P(0,0,1)=1.D0
	IF(lmax.GT.0) P(1,0,1)=xx
	IF(lmax.GT.1)THEN
	   DO L=2,lmax
	    P(L,0,1)=( (2.D0*L-1)*xx*P(L-1,0,1) - (L-1)*P(L-2,0,1) ) / L
	   ENDDO
	ENDIF

	fact=1.D0-xx**2
	sfact=DSQRT(fact)

	IF(mmax .GT. 0)THEN
		DO M=1,mmax
		  DO L=M,lmax
			L2=MAX( L-2 ,  0 )
			P(L,M,1)= P(L2,M,1) -(2*L-1)*P(L-1,M-1,0)*fact
			P(L,M,0)= P(L2,M,0) -(2*L-1)*P(L-1,M-1,1)
		  ENDDO
	    ENDDO
	ENDIF

	IF(Derivative)Then
* First zero arrays
		DO l=0,lmax
			DO m=0,mmax
				dPlm(l,m)=0.D0
				dP(l,m,0)=0.D0
				dP(l,m,1)=0.D0
				dP(l,m,2)=0.D0
			ENDDO
		ENDDO

		IF(lmax .GT. 0) dP(1,0,1)=1.D0

		IF(lmax .GT. 1)THEN
			DO L=2,lmax  
				dP(L,0,1)=( (2*L-1)*P(L-1,0,1) + 
     $                  (2*L-1)*xx*dP(L-1,0,1) - 
     $                  (L-1)*dP(L-2,0,1) ) / L
			ENDDO
		ENDIF

		IF(mmax .GT. 0)THEN
		  DO M=1,mmax  
		    DO L=M,lmax
		      L2=MAX( L-2 ,  0 )
		      dP(L,M,1)= dP(L2,M,1) - (2*L-1)*fact*dP(L-1,M-1,0)
     $                 - (2*L-1)*dP(L-1,M-1,2) + (2*L-1)*xx*P(L-1,M-1,0)
		      dP(L,M,0)= dP(L2,M,0) - (2*L-1)*dP(L-1,M-1,1)
		      dP(L,M,2)=dP(L2,M,2) +(2*L-1)*xx*P(L-1,M-1,1)
		    ENDDO
		  ENDDO
		ENDIF

		DO L=0,lmax  
	      mlimit=MIN(mmax,L)
		  DO M=0,mlimit
* Prevent a divide by zero
		    anum=dP(L,M,2) !numerator
			IF(sfact.NE.0.)Then !denominator is OK
			  term=anum/sfact 
			ELSE !denominator is zero
			  IF(DABS(anum).LT.1.D-7)THEN
				term=0.D0 !return 0 in cases where numerator is near zero
			  ELSE !return nearly infinity with same sign as numerator
				term=DSIGN(1.D36,anum) 
			  ENDIF
			ENDIF
			dPlm(L,M)=dP(L,M,1) + dP(L,M,0)*sfact + term
		  ENDDO
		ENDDO

	ENDIF !End doing derivative

	DO L=0,lmax
	    mlimit=MIN(mmax,L)
	    DO M=0,mlimit
		  Plm(L,M)=P(L,M,1) + P(L,M,0)*sfact
	    ENDDO
	ENDDO	

	RETURN
	END
************************ Copyright 1996,2001 Dan Weimer/MRC ***********************

	SUBROUTINE READCOEF01 (udat)
CNCAR      Feb 01: Make a subroutine (a.k.a. ReadCoef) from FIRST=TRUE block
C          of SETMODEL01 (a.k.a. SetModel) so that one can load the coefficients
C          independently of defining geophysical conditions.  This mimics the
C          1996 model.
C          Sep 01: Change argument list to pass both unit number and file name
C          rather than hard coding
CNCAR

	INTEGER udat
	CHARACTER*30 Copyright
	PARAMETER (MJ=3,ML=4,MM=3,MN=2,MO=2)
	REAL  CS( 0:MJ, 0:1 , 0:MO, 0:1 , 0:ML, 0:MM)
	REAL BCS( 0:MJ, 0:1 , 0:MO, 0:1 , 0:MN)
	REAL  SS( 0:1 , 0:1 , 0:MO, 0:1 , 0:ML, 0:MM)
	REAL BSS( 0:1 , 0:1 , 0:MO, 0:1 , 0:MN)
	REAL Coef(0:1,0:5,0:5),BoundFit(0:1,0:5),pi
	DOUBLE PRECISION dpi
	INTEGER*4 maxj,MaxL,MaxM,MaxN,MaxO
	COMMON /AllW2kCoefs/MaxJ,MaxO,CS,BCS,SS,BSS
	COMMON /SetW2kCoef/MaxL,MaxM,MaxN,Coef,BoundFit,pi,dpi

CNCAR      Feb 01:  Initialize constants used in GECMP01
C          Sep 01:  Omit unneeded min lat variables because of switch to
C                   constant potential at and below min lat (hardcoded
C                   in EPOTVAL01).
      COMMON /CECMP/ ALAMX,STPD,STP2,CSTP,SSTP
C            ALAMX = Absolute max latitude (deg) for normal gradient calc.
C            STPD  = Angular dist (deg) of step @ 300km above earth (r=6371km)
C            STP2  = Denominator in gradient calc
      PARAMETER ( D2R =  0.0174532925199432957692369076847 ,
     +            R2D = 57.2957795130823208767981548147)
      STEP = 10.
      STPR = STEP/6671.
      STPD = STPR*R2D
      STP2 = 2.*STEP
      CSTP = COS (STPR)
      SSTP = SQRT (1. - CSTP*CSTP)
      ALAMX = 90. - STPD
CNCAR

CNCAR      Sep 01:  udat,cfile are now a formal arguments
C  Assume file has been opened at unit elsewhere
c     cfile = 'w2k.dat' !make sure correct ENDIAN type file is used.
c     udat  = 99
CNCAR
c     OPEN (UNIT=udat,FILE=cfile,STATUS='OLD',form='UNFORMATTED')
      READ (udat) Copyright
      PRINT *,Copyright
      READ (udat) Maxj,MaxL,MaxM,MaxN,MaxO
      If (maxj.NE.MJ .OR. MaxL.NE.ML .OR. MaxM.NE.MM .OR.
     $    MaxN.NE.MN .OR. MaxO.NE.MO) Then
	      PRINT *,'Data File Error'
	      STOP !Data file did not match size expected for arrays
      Endif
      READ (udat) CS
      READ (udat) BCS
      READ (udat) SS
      READ (udat) BSS
      CLOSE (udat)
      pi=2.*ASIN(1.)
      dpi=2.D0*DASIN(1.D0)

      RETURN
      END
************************ Copyright 1996,2001 Dan Weimer/MRC ***********************
CNCAR      Sep 01: Change name from SetModel to be distinct

	SUBROUTINE SetModel01 (angle,Bt,Tilt,SWVel,SWDen,ALindex,UseAL)
C       SUBROUTINE SetModel   (angle,Bt,Tilt,SWVel,SWDen,ALindex,UseAL)
CNCAR
*
* Calculate the complete set of spherical harmonic coeficients,
* given an aribitrary IMF angle (degrees from northward toward +Y),
* magnitude Bt (nT), dipole tilt angle (degrees), 
* solar wind velocity (km/sec), SWDen (#/cc),
* ALindex (nT), and Logical flag to use optional AL index.
*
* Sets the value of Coef and Boundfit in the common block SetW2kCoef.
*
	REAL angle,Bt,Tilt,SWVel,SWDen,ALindex
C       LOGICAL First,UseAL
	LOGICAL UseAL

CNCAR      Feb 01: Move logic block to subroutine RdCoef01
C       DATA First/.TRUE./
C       SAVE First
C       INTEGER unit
C       CHARACTER*15 cfile
C       CHARACTER*30 Copyright
	PARAMETER (MJ=3,ML=4,MM=3,MN=2,MO=2)
	REAL  CS( 0:MJ, 0:1 , 0:MO, 0:1 , 0:ML, 0:MM)
	REAL BCS( 0:MJ, 0:1 , 0:MO, 0:1 , 0:MN)
	REAL  SS( 0:1 , 0:1 , 0:MO, 0:1 , 0:ML, 0:MM)
	REAL BSS( 0:1 , 0:1 , 0:MO, 0:1 , 0:MN)
	REAL XA(0:MJ),XB(0:MJ),FSC(0:1,0:4),PSS(0:1)
	REAL Coef(0:1,0:5,0:5),BoundFit(0:1,0:5),pi
	DOUBLE PRECISION dpi
	INTEGER*4 i,j,k,l,m,n,o
	INTEGER*4 maxj,MaxL,MaxM,MaxN,MaxO
	COMMON /AllW2kCoefs/MaxJ,MaxO,CS,BCS,SS,BSS
	COMMON /SetW2kCoef/MaxL,MaxM,MaxN,Coef,BoundFit,pi,dpi

CNCAR      Feb 01: Move logic block to subroutine RdCoef01
C  All First=TRUE in new subroutine ReadCoef
C       If(First)Then
C       cfile='w2k.dat' !make sure correct ENDIAN type file is used.
C       unit=99
C       OPEN(UNIT=unit,FILE=cfile,STATUS='OLD',form='UNFORMATTED')
C       READ(unit) Copyright
C       PRINT *,Copyright
C       READ(unit) Maxj,MaxL,MaxM,MaxN,MaxO
C       If(maxj.NE.MJ .OR. MaxL.NE.ML .OR. MaxM.NE.MM .OR.
C    $   MaxN.NE.MN .OR. MaxO.NE.MO)Then
C               PRINT *,'Data File Error'
C               STOP !Data file did not match sixe expected for arrays
C       Endif
C       READ(unit) CS
C       READ(unit) BCS
C       READ(unit) SS
C       READ(unit) BSS
C       CLOSE(unit)
C       pi=2.*ASIN(1.)
C       dpi=2.D0*DASIN(1.D0)
C       First=.FALSE.
C       Endif
CNCAR

	SinTilt=SIND(Tilt)
	omega=angle*pi/180.
	XA(0)=1.
	XA(1)=Bt**(2./3.) *SWvel
	XA(2)=SinTilt
	XA(3)=SWvel**2 *SWDen
	XB(0)=1.
	XB(1)=Bt
	XB(2)=SinTilt
	XB(3)=SWvel**2 *SWDen

	DO l=0,MaxL
	    mlimit=MIN(l,MaxM)
	    DO m=0,mlimit
		  klimit=MIN(m,1)
		  DO k=0,klimit
			acoef=0. !rezero for summation
			DO j=0,MaxJ
			    DO o=0,MaxO
				  DO i=0,1
					FSC(i,o)=CS(j,i,o,k,l,m)
				  ENDDO
			    ENDDO
     			    acoef=acoef+ XA(j)*FSVAL(omega,MaxO,fsc)
			ENDDO
			IF(UseAL)THEN
			    DO j=0,1
				  DO o=0,MaxO
					DO i=0,1
					    FSC(i,o)=SS(j,i,o,k,l,m)
					ENDDO
				  ENDDO
				  PSS(j)=FSVAL(omega,MaxO,fsc)
			    ENDDO
			    acoef=acoef + PSS(0) + PSS(1)*ALindex
			ENDIF
			Coef(k,l,m)=acoef
		  ENDDO
	    ENDDO
	ENDDO

	DO n=0,MaxN
	    klimit=MIN(n,1)
	    DO k=0,klimit
		  acoef=0. !rezero for summation
		  DO j=0,MaxJ
			DO o=0,MaxO
			    DO i=0,1
				  FSC(i,o)=BCS(j,i,o,k,n)
			    ENDDO
			ENDDO
     			acoef=acoef+ XB(j)*FSVAL(omega,MaxO,fsc)
		  ENDDO
		  IF(UseAL)THEN
			DO j=0,1
			    DO o=0,MaxO
				  DO i=0,1
					FSC(i,o)=BSS(j,i,o,k,n)
				  ENDDO
			    ENDDO
			    PSS(j)=FSVAL(omega,MaxO,fsc)
			ENDDO
			acoef=acoef + PSS(0) + PSS(1)*ALindex
		  ENDIF
		  BoundFit(k,n)=acoef
	    ENDDO
	ENDDO
	RETURN
	END
****************** Copyright 1996, 2001, Dan Weimer/MRC ***********************
	FUNCTION BoundaryLat(gmlt)
	REAL gmlt
	REAL Coef(0:1,0:5,0:5),BoundFit(0:1,0:5),pi
	DOUBLE PRECISION dpi
	INTEGER MaxL,MaxM,MaxN
	COMMON /SetW2kCoef/MaxL,MaxM,MaxN,Coef,BoundFit,pi,dpi
	BoundaryLat=FSVal(gmlt*pi/12.,MaxN,BoundFit)
	RETURN
	END
****************** Copyright 1996, 2001, Dan Weimer/MRC ***********************
CNCAR      Sep 01: Change name from EpotVal to be distinct

	FUNCTION EPOTVAL01 (gLAT,gMLT)
C       FUNCTION EpotVal   (gLAT,gMLT)
CNCAR

* Return the value of the electric potential in kV at
* corrected geomagnetic coordinates gLAT (degrees) and gMLT (hours).
*
* Must first call ReadCoef and SetModel to set up the model coeficients for
* the desired values of Bt, IMF clock angle, Dipole tilt angle, and SW Vel.
*
      REAL gLAT,gMLT
      DOUBLE PRECISION Phi,Z,O,x,ct,Phim
      DOUBLE PRECISION Plm(0:10,0:10), OPlm(0:10,0:10), DPLM(0:10,0:10)

      REAL Coef(0:1,0:5,0:5), BoundFit(0:1,0:5),pi
      DOUBLE PRECISION dpi
      COMMON /SetW2kCoef/ MaxL,MaxM,MaxN,Coef,BoundFit,pi,dpi

      blat = BoundaryLat (gMLT)
CNCAR      Feb 01:  For latitudes at and equatorward of the model minimum
C          (a function of MLT), limit the latitude used to never be less
C          than the model minimum, thus returning a constant potential for
C          points at and equatorward of BoundaryLat (gmlt).
C     IF (glat .GT. blat) THEN
      glatlim = amax1 (glat,blat)
CNCAR

      Phi = DBLE (gMLT)*dpi/12.D0

CNCAR
      colat = 90.-glatlim
C     colat = 90.-glat
CNCAR

      bcolat  = 90.-blat
      x  = DBLE (colat)*dpi/DBLE(bcolat)
      DC = DBLE (Coef(0,0,0))
      Z  = DC
      O  = DC
      ct = DCOS(x)
      DPLM(0:10,0:10)=0.D0
      CALL DLegendre (ct,   MaxL,MaxM,Plm ,DPLM,.FALSE.)
C          Also find value at outer boundary at angle Pi, or cos(pi)=-1.
      CALL DLegendre (-1.D0,MaxL,MaxM,OPlm,DPLM,.FALSE.)
      DO l=1,MaxL
	Z = Z +  Plm(l,0)*DBLE(Coef(0,l,0))
	O = O + OPlm(l,0)*DBLE(Coef(0,l,0))
	mlimit = MIN(l,MaxM)
	DO m=1,mlimit
	  phim = phi*m
	  Z = Z + Plm(l,m)*(DBLE(Coef(0,l,m))*DCOS(Phim) +
     $                      DBLE(Coef(1,l,m))*DSIN(Phim) )
	  O = O +OPlm(l,m)*DBLE(Coef(0,l,m))
	ENDDO
      ENDDO
CNCAR
      EPOTVAL01 = SNGL(Z-O)
C     EpotVal   = SNGL(Z-O)
C     ELSE
C       EpotVal=0.
C     ENDIF
CNCAR
      RETURN
      END

CNCAR
      SUBROUTINE GECMP01 (AMLA,RMLT,ET,EP)
C          Get Electric field components for the 2001 Weimer electrostatic
C          potential model.  Before use, first load coefficients (CALL
C          READCOEF01) and initialize model conditions (CALL SETMODEL01).
C          The electric field, the gradient of the electic potential,
C          is determined using finite differences over a distance of
C          STEP km, defined in READCOEF01 and accessed here via common
C          block CECMP.
C
C          It has been useful to modify STEP.  One application discovered
C          the Weimer electrostatic potential can have minor unrealistic
C          perturbations over short distances.  To avoid these ripples
C          STEP was increased to 5 degrees arc (582 km at 300 km altitude,
C          R=6671 km).
C
C          INPUTS:
C            AMLA = Absolute value of magnetic latitude (deg)
C            RMLT = Magnetic local time (hours).
C          RETURNS:
C            ET = Etheta (magnetic equatorward*) E field component (V/m)
C            EP = Ephi   (magnetic eastward)     E field component (V/m)
C
C          * ET direction is along the magnetic meridian away from the
C            current hemisphere; i.e., when ET > 0, the direction is
C              southward when in northern magnetic hemisphere
C              northward when in southern magnetic hemisphere
C          Since the Weimer model uses only NH (positive) latitudes, the
C          sign of ET for SH should be changed outside of this routine.
C
C          HISTORY:
C          Jan 97:  Initial implementation at NCAR for the 1996 version.
C          Feb 01:  Error corrected in determining DPHI.  Old version
C          used wrong spherical right triangle formula (inadvertently
C          computing along a latitude line); more importantly, it
C          neglected to convert mlt from degrees to hour angle input
C          to EPOTVAL01.
C          Sep 01:  Revised equatorward boundary logic to the scheme Barbara
C          uses for AMIE:  As absolute magnetic latitude decreases to the
C          model minimum, the electric potential from EPOTVAL01 now goes to
C          a non-zero constant (rather than zero), thus obviating the need
C          for special handling here.  Therefore, special logic for determining
C          gradients at lower limit has been removed.

      PARAMETER ( D2R =  0.0174532925199432957692369076847 ,
     +            R2D = 57.2957795130823208767981548147)

C          CECMP contains constants initialized in READCOEF01
      COMMON /CECMP/ ALAMX,STPD,STP2,CSTP,SSTP

      ET = -99999.
      EP = -99999.
      IF (AMLA .LT. 0.) GO TO 100

C          Calculate -(latitude gradient) by stepping along the magnetic
C          latitude line in each direction (flipping coordinates when
C          going over pole to keep lat <= 90).
      KPOL  = 0
      XMLT  = RMLT
   10 XMLT1 = XMLT
      AMLA1 = AMLA + STPD
      IF (AMLA1 .GT. 90.) THEN
	AMLA1 = 180. - AMLA1
	XMLT1 = XMLT1 + 12.
      ENDIF
      P1 = EPOTVAL01 (AMLA1    ,XMLT1)
      P2 = EPOTVAL01 (AMLA-STPD,XMLT )
      IF (KPOL .EQ. 1) GO TO 20
      ET = (P1 - P2) / STP2

C          Calculate -(lon gradient).  Step along the magnetic meridion
C          in both directions to obtain the electric potential
      IF (AMLA .LT. ALAMX) THEN
	AMLA1 = ASIN(SIN(AMLA*D2R)*CSTP)
	DPHI  = ASIN (SSTP/COS(AMLA1))*R2D/15.      ! 15 converts from degrees to hours
	AMLA1 = AMLA1*R2D
	P1 = EPOTVAL01 (AMLA1,XMLT+DPHI)
	P2 = EPOTVAL01 (AMLA1,XMLT-DPHI)
      ELSE
	AMLA = 90.
	XMLT = XMLT + 6.
	KPOL = 1
	GO TO 10
      ENDIF
   20 EP = (P2 - P1) / STP2
      IF (KPOL .EQ. 1) EP = -EP

  100 RETURN
      END   

      SUBROUTINE WEIEPOT01 (IYR,IMO,IDA,IHR,IMN,SWS,SWD,BY,BZ,AL,USEAL,
     +                     IHEM,ISETM,RMLA,RMLT, ET,EP,EPOT)
C          Interface to Weimer-01 (a.k.a w2k) for AMIE's combined electrostatic
C          potential models.  This replaces two calls (SETMODEL01 and
C          EPOTVAL01) and automates sign changes for southern hemisphere.
C          INPUTS:
C            IYR   = UT year ((4-digit)
C            IMO   = month of IYR
C            IDA   = day of IMO
C            IHR   = hour of day
C            IMN   = min of hour
C            SWS   = Solar wind speed (km/s)
C            SWD   = Solar wind density (#/cm3)
C            BY    = IMF By component in GSM coordinates (nt)
C            BZ    = IMF Bz component in GSM coordinates (nt)
C            AL    = AL index (nT)
C            USEAL = logical (T or F) to use AL index or not
C            IHEM  = Hemisphere flag: (-1) southern, (1) northern
C            ISETM = Model conditions change flag: (0) no-change
C                                                  (1) time, IMF or SW changed
C            RMLA  = Magnetic latitude  (deg) of point to determine Potential
C                    (RMLA should be NH positive latitudes only, since SH
C                    values found by changing sign of By (or ANGL) and tilt.)
C            RMLT  = Magnetic longitude (hrs) of point to determine Potential
C          RETURNS:
C            ET   = Etheta (magnetic equatorward*) E field component (V/m)
C            EP   = Ephi   (magnetic eastward)     E field component (V/m)
C            EPOT = Electrostatic potential (kV)
C
C          * ET direction is along the magnetic meridian away from the
C            current hemisphere; i.e., when ET > 0, the direction is
C              southward when in northern magnetic hemisphere
C              northward when in southern magnetic hemisphere
C
C          Since the AMIE model assumes a NH solution although the
C          latitudes are negative for SH data, the sign of ET
C          for the SH should be changed outside of this routine.

C          HISTORY:
C          Jan 97:  Initial implementation. B. Emery.
C          Feb 01:  Remove bug (HR was not defined)
C          Sep 01:  Add common block s.t. computed values are available to
C          the calling routine without changing the argument list
      COMMON /WEIAT/ ANGL, TILT
      PARAMETER (R2D=57.2957795130823208767981548147)
      LOGICAL USEAL

      H = REAL (IHEM)

      IF (ISETM .EQ. 1) THEN
	ANGL = ATAN2 (BY,BZ)*R2D
	BT   = SQRT (BY*BY + BZ*BZ)
	HR = FLOAT(IHR) + FLOAT(IMN)/60.
	TILT = GET_TILT (IYR,IMO,IDA,HR)

	HANGL = H * ANGL
	HTILT = H * TILT
	CALL SETMODEL01 (HANGL,BT,HTILT,SWS,SWD,AL,USEAL)
C     WRITE (6,'('' WEIMER-2001: '',I4,5I3,8F7.2,L)') IYR,IMO,IDA,IHR,
C    +                    IMN,IHEM,SWS,SWD,AL,BY,BZ,BT,HANGL,HTILT,USEAL
      ENDIF

C          NH assumed latitudes only
      AMLA  = ABS (RMLA)
      EPOT = EPOTVAL01 (AMLA,RMLT)
      CALL GECMP01 (AMLA,RMLT,ET,EP)

      RETURN
      END
CNCAR
