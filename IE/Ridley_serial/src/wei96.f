************************ Copyright 1996, Dan Weimer/MRC ***********************
*
* Subroutines to calculate the electric potentials from the Weimer '96 model of
* the polar cap ionospheric electric potentials.
*
* To use, first call subroutine ReadCoef once.
* Next, call SetModel with the specified input parameters.
* The function EpotVal(gLAT,gMLT) can then be used repeatively to get the
* electric potential at the desired location in geomagnetic coordinates.
*
* This code is protected by copyright and is
* distributed for research or educational use only.
* Commerical use without written permission from Dan Weimer/MRC is prohibited.
*

CNCAR      Revisions for use at NCAR:
C            (1) Change behavior at minimum magnetic latitude.  When approaching
C                the model minimum absolute latitude (45 degrees) the electric
C                potential returned used to go to zero discontinuously; although
C                intended as a flag, it created artificial gradients in the
C                electric field calculation.  Now the potential returned is
C                constant (that of the minimum latitude) for any latitude at
C                or equatorward of the minimum.
C            (2) Accommodate running sumultaneously 1996 and 2001 versions.  To
C                avoid name collisions this required: (i) revising names (e.g.,
C                adding '96') for differing subprograms, and (ii) relocating
C                common routines into another file (weicom.f).
C            (3) Pass the coefficients file name and unit number into RDCOEF96
C                rather than using hard coded values.
C            (4) Add wrapper subroutines for non-ANSI trig functions which
C                input angles in degrees.
C            (5) Add electric field routine (GECMP96) to deterine the electric
C                potential gradient.
C            (6) Add wrapper routine (WEIEPOT96) for use with AMIE; this is a
C                substitute for calling SETMODEL96 and EPOTVAL96.
CNCAR      NCAR changes are delimited by "CNCAR"

************************ Copyright 1996, Dan Weimer/MRC ***********************

CNCAR      Oct 01:  Distinguish from the 2001 version

	FUNCTION EPOTVAL96 (gLAT,gMLT)
C       FUNCTION EPOTVAL   (gLAT,gMLT)
CNCAR
* Return the value of the electric potential in kV at
* corrected geomagnetic coordinates gLAT (degrees) and gMLT (hours).
*
* Must first call ReadCoef and SetModel to set up the model coeficients for
* the desired values of Bt, IMF clock angle, Dipole tilt angle, and SW Vel.
*
	REAL gLAT,gMLT
	Real Theta,Phi,Z,ct,Phim
	REAL Plm(0:20,0:20)

	REAL Coef(0:1,0:8,0:3),pi
	INTEGER ML,MM
	COMMON/SetCoef/ML,MM,Coef,pi

	r = 90.-gLAT

CNCAR      Sep 01: Limit the colatitude (r) such that an input absolute
C          magnetic latitude of 45 degrees or less always calculates the
C          potential at 45.
	R = AMIN1 (R,45.)
C       IF(r .LT. 45.)THEN
CNCAR
	  Theta=r*pi/45.
          Phi=gMLT*pi/12.
	  Z=Coef(0,0,0)
	  ct=COS(Theta)
	  CALL Legendre(ct,ML,MM,Plm)
	  DO l=1,ML
	    Z=Z + Coef(0,l,0)*Plm(l,0)
	    IF(l.LT.MM)THEN
	      limit=l
	    ELSE
	      limit=MM
	    ENDIF
	    DO m=1,limit
	      phim=phi*m
	      Z=Z + Coef(0,l,m)*Plm(l,m)*COS(phim) +
     $		   Coef(1,l,m)*Plm(l,m)*SIN(phim)
	    ENDDO
	  ENDDO
CNCAR
C       ELSE
C         Z=0.
C       ENDIF
	EPOTVAL96 = Z
C       EpotVal   = Z
CNCAR
	RETURN
	END

************************ Copyright 1996, Dan Weimer/MRC ***********************

CNCAR      Sep 01:  Pass in the unit no. and file name rather than hard coding
C          inside.  Also change name to distinguish 2001 version.

	SUBROUTINE READCOEF96 (udat)
C       SUBROUTINE ReadCoef
CNCAR

*
* Read in the data file with the model coefficients
*
	INTEGER udat
CNCAR
	CHARACTER    skip*15
C       CHARACTER*15 cfile,skip
CNCAR
	REAL C(0:3)
	REAL Cn( 0:3 , 0:1 , 0:4 , 0:1 , 0:8 , 0:3 )
	INTEGER MaxL,MaxM,MaxN
	COMMON /AllCoefs/MaxL,MaxM,MaxN,Cn


CNCAR      Jan 97:  Initialize constants used in GECMP96
C          Sep 01:  Omit unneeded min lat variables because of switch to
C                   constant potential at and below min lat (hardcoded
C                   in EPOTVAL96).
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
C       cfile='wei96.cofcnts'
C       udat=99
CNCAR
C	OPEN(udat,FILE=cfile,STATUS='OLD')
  900   FORMAT(A15)
 1000	FORMAT(3I8)
 2000	FORMAT(3I2)
 3000	FORMAT(2I2,4E15.6)

	READ(udat,900) skip
	READ(udat,1000) MaxL,MaxM,MaxN
	DO l=0,MaxL
	  IF(l.LT.MaxM)THEN
	    mlimit=l
	  ELSE
	    mlimit=MaxM
	  ENDIF
	  DO m=0,mlimit
	    IF(m.LT.1)THEN
	      klimit=0
	    ELSE
	      klimit=1
	    ENDIF
	    DO k=0,klimit
	      READ(udat,2000) ll,mm,kk
	      IF(ll.NE.l .OR. mm.NE.m .OR. kk.NE.k)THEN
		PRINT *,'Data File Format Error'
		STOP
	      ENDIF
	      DO n=0,MaxN
	        IF(n.LT.1)THEN
	          ilimit=0
	        ELSE
	          ilimit=1
	        ENDIF
		DO i=0,ilimit
		  READ(udat,3000) nn,ii,C
	          IF(nn.NE.n .OR. ii.NE.i)THEN
		    PRINT *,'Data File Format Error'
		    STOP
		  ENDIF
		  Cn(0,i,n,k,l,m)=C(0)
		  Cn(1,i,n,k,l,m)=C(1)
		  Cn(2,i,n,k,l,m)=C(2)
		  Cn(3,i,n,k,l,m)=C(3)
		ENDDO
	      ENDDO
	    ENDDO
	  ENDDO
	ENDDO

	CLOSE(udat)
	RETURN
	END
************************ Copyright 1996, Dan Weimer/MRC ***********************
CNCAR      Oct 01:  Change name from SetModel to be distinct

	SUBROUTINE SETMODEL96 (angle,Bt,Tilt,SWVel)
C       SUBROUTINE SetModel   (angle,Bt,Tilt,SWVel)
CNCAR
*
* Calculate the complete set of spherical harmonic coefficients,
* given an arbitrary IMF angle (degrees from northward toward +Y),
* magnitude Bt (nT), dipole tilt angle (degrees),
* and solar wind velocity (km/sec).
* Returns the Coef in the common block SetCoef.
*
	REAL angle,Bt,Tilt,SWVel
	REAL FSC(0:1,0:4)
	REAL Cn( 0:3 , 0:1 , 0:4 , 0:1 , 0:8 , 0:3 )
	INTEGER MaxL,MaxM,MaxN
	COMMON /AllCoefs/MaxL,MaxM,MaxN,Cn

	REAL Coef(0:1,0:8,0:3),pi
	INTEGER ML,MM
	COMMON/SetCoef/ML,MM,Coef,pi

	pi=2.*ASIN(1.)
	ML=MaxL
	MM=MaxM
	SinTilt=SIND(Tilt)

	omega=angle*pi/180.
	DO l=0,MaxL
	  IF(l.LT.MaxM)THEN
	    mlimit=l
	  ELSE
	    mlimit=MaxM
	  ENDIF
	  DO m=0,mlimit
	    IF(m.LT.1)THEN
	      klimit=0
	    ELSE
	      klimit=1
	    ENDIF
	    DO k=0,klimit
* Retrieve the regression coefficients and evaluate the function
* as a function of Bt,Tilt,and SWVel to get each Fourier coefficient.
	      DO n=0,MaxN
	        IF(n.LT.1)THEN
	          ilimit=0
	        ELSE
	          ilimit=1
	        ENDIF
		DO i=0,ilimit
		  FSC(i,n)=Cn(0,i,n,k,l,m) + Bt*Cn(1,i,n,k,l,m) +
     $		   SinTilt*Cn(2,i,n,k,l,m) + SWVel*Cn(3,i,n,k,l,m)
		ENDDO
	      ENDDO
* Next evaluate the Fourier series as a function of angle.
      	      Coef(k,l,m)=FSVal(omega,MaxN,FSC)
	    ENDDO
	  ENDDO
	ENDDO
	RETURN
	END
************************ Copyright 1996, Dan Weimer/MRC ***********************
	SUBROUTINE LEGENDRE(x,lmax,mmax,Plm)
* compute Associate Legendre Function P_l^m(x)
* for all l up to lmax and all m up to mmax.
* returns results in array Plm
* if X is out of range ( abs(x)>1 ) then value is returned as if x=1.
	DIMENSION Plm(0:20,0:20)
	  DO l=0,20
	    DO m=0,20
		Plm(l,m)=0.
	    ENDDO
	  ENDDO
	xx=MIN(x,1.)
	xx=MAX(xx,-1.)
	IF(lmax .LT. 0 .OR. mmax .LT. 0 .OR. mmax .GT. lmax )THEN
	  Print *,'Bad arguments to Legendre'
	  RETURN
	ENDIF
* First calculate all Pl0 for l=0 to l
	Plm(0,0)=1.
	IF(lmax.GT.0)Plm(1,0)=xx
	IF (lmax .GT. 1 )THEN
	  DO L=2,lmax
	    Plm(L,0)=( (2.*L-1)*xx*Plm(L-1,0) - (L-1)*Plm(L-2,0) )/L
	  ENDDO
	ENDIF
	IF (mmax .EQ. 0 )RETURN
	fact=SQRT( (1.-xx)*(1.+xx) )
	DO M=1,mmax
	  DO L=m,lmax
	    lm2=MAX(L-2,0)
	    Plm(L,M)=Plm(lm2,M) - ( 2*L-1)*fact*Plm(L-1,M-1)
	  ENDDO
	ENDDO
	RETURN
	END

CNCAR
      SUBROUTINE GECMP96 (AMLA,RMLT,ET,EP)
C          Get Electric field components for the 1996 Weimer electrostatic
C          potential model.  Before use, first load coefficients (CALL
C          READCOEF96) and initialize model conditions (CALL SETMODEL96).
C          The electric field, the gradient of the electic potential,
C          is determined using finite differences.
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
C
C          HISTORY:
C          Jan 97:  Initial implementation at NCAR.  R.Barnes
C          Feb 01:  Error corrected in determining DPHI.  Old version
C          used wrong spherical right triangle formula (inadvertently
C          computing along a latitude line); more importantly, it
C          neglected to convert mlt from degrees to hour angle input
C          to EPOTVAL96.
C          Sep 01:  Revised equatorward boundary logic to the scheme Barbara
C          uses for AMIE:  As absolute magnetic latitude decreases to the
C          model minimum, the electric potential from EPOTVAL96 now goes to
C          a non-zero constant (rather than zero), thus obviating the need
C          for special handling here.  Therefore, special logic for determining
C          gradients at lower limit has been removed.

      PARAMETER ( D2R =  0.0174532925199432957692369076847 ,
     +            R2D = 57.2957795130823208767981548147)

C          CECMP contains constants initialized in READCOEF
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
      P1 = EPOTVAL96 (AMLA1    ,XMLT1)
      P2 = EPOTVAL96 (AMLA-STPD,XMLT )
      IF (KPOL .EQ. 1) GO TO 20
      ET = (P1 - P2) / STP2

C          Calculate -(lon gradient).  Step along the magnetic meridion
C          in both directions to obtain the electric potential
      IF (AMLA .LT. ALAMX) THEN
	AMLA1 = ASIN(SIN(AMLA*D2R)*CSTP)

C       DPHI  = ASIN (SSTP/SIN(AMLA1))*R2D          ! old and wrong! (changed Feb 01)
	DPHI  = ASIN (SSTP/COS(AMLA1))*R2D/15.      ! 15 converts from degrees to hours

	AMLA1 = AMLA1*R2D
	P1 = EPOTVAL96 (AMLA1,XMLT+DPHI)
	P2 = EPOTVAL96 (AMLA1,XMLT-DPHI)
      ELSE
C          At the pole.  Avoid a divide by zero by using Art's trick
C          where Ephi(90,lon) = Etheta(90,lon+90).
	AMLA = 90.
	XMLT = XMLT + 6.
	KPOL = 1
	GO TO 10
      ENDIF
   20 EP = (P2 - P1) / STP2
      IF (KPOL .EQ. 1) EP = -EP

  100 RETURN
      END   

      SUBROUTINE WEIEPOT96 (IYR,IMO,IDA,IHR,IMN,SWS,BY,BZ,IHEM,ISETM,
     +                     RMLA,RMLT, ET,EP,EPOT)
C          Interface to Weimer-96 for AMIE's combined electrostatic
C          potential models.  This replaces two calls (SETMODEL96 and
C          EPOTVAL96) and automates sign changes for southern hemisphere.
C          INPUTS:
C            IYR   = UT year ((4-digit)
C            IMO   = month of IYR
C            IDA   = day of IMO
C            IHR   = hour of day
C            IMN   = min of hour
C            SWS   = Solar wind speed (km/s)
C            BY    = IMF By component in GSM coordinates (nt)
C            BZ    = IMF Bz component in GSM coordinates (nt)
C            IHEM  = Hemisphere flag: (-1) southern, (1) northern
C            ISETM = Model conditions change flag: (0) no-change
C                                                  (1) time, IMF or SW changed
C            RMLA  = Magnetic latitude (deg) of point to determine Potential.
C                    RMLA should be positive (NH) only, since SH values are
C                    obtained by changing sign of By (or ANGL) and TILT.
C            RMLT  = Magnetic longitude (hrs) of point to determine Potential
C          RETURNS:
C            ET    = Etheta (magnetic equatorward*) E field component (V/m)
C            EP    = Ephi   (magnetic eastward)     E field component (V/m)
C            EPOT  = Electric potential (kV)
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
C          Oct 01:  Revise routine names s.t. this can run with the 2001 version
      COMMON /WEIAT/ ANGL, TILT
      PARAMETER (R2D=57.2957795130823208767981548147)

      H = REAL (IHEM)

      IF (ISETM .EQ. 1) THEN
	ANGL = ATAN2 (BY,BZ)*R2D
	BT   = SQRT (BY*BY + BZ*BZ)
	HR   = FLOAT(IHR) + FLOAT(IMN)/60.       ! define HR (bug fix Feb 01)
	TILT = GET_TILT (IYR,IMO,IDA,HR)

	HANGL = H * ANGL
	HTILT = H * TILT
	CALL SETMODEL96 (HANGL,BT,HTILT,SWS)
C       WRITE (6,'('' WEIMER-1996: sw by bz bt angle tilt='',6f7.2)')
C    +                            SWS,BY,BZ,BT,ANGL,TILT
      ENDIF

C          NH assumed latitudes only
      AMLA = ABS (RMLA)
      EPOT = EPOTVAL96 (AMLA, RMLT)
      CALL GECMP96 (AMLA,RMLT,ET,EP)

      RETURN
      END
CNCAR
