c
c     The code has been made to implicit real*8 by Mei-Ching Fok on
c     Jan. 30, 2002

************************ Copyright 1996,2001 Dan Weimer/MRC ***********************

*

* Subroutines to calculate the electric potentials from the "Weimer 2K" model of

* the polar cap ionospheric electric potentials described in the publication: 

* Weimer, D. R., An improved model of ionospheric electric potentials including

* substorm perturbations and application to the Geospace Environment Modeling

* November 24, 1996 event, Journal of Geophysical Research, Vol. 106, p. 407, 2001.

*

* To use, first call procedure SETMODEL00 with the specified input parameters:

*   angle: IMF Y-Z clock angle in degrees, 0=northward, 180=southward

*   Bt: Magnitude of IMF in Y-Z plane in nT

*   Tilt: dipole tilt angle in degrees.

*   SWVel: solar wind velocity in km/sec

*   SWDen: solar wind density in #/cc

*   ALindex: (optional) AL index in nT

*

* The function EPOTVAL00(gLAT,gMLT) can then be used repeatively to get the

* electric potential in kV at the desired location.

* Input coordinates assume use of 'altitude adjusted' corrected geomagnetic

* coordinates for R=1, also refered to as AACGM0.

*

* The function BOUNDARYLAT00(gMLT) can be used to get the latitude of the boundary

*   where the potential goes to zero.  This boundary is a function of MLT, and

*   varies with the SETMODEL00 parameters.  The potential is zero everywhere below

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

*

************************ Copyright 1996,2001 Dan Weimer/MRC ***********************

	SUBROUTINE DLEGENDRE00(x,lmax,mmax,Plm,dPlm)

* compute Double Precision Associate Legendre Function P_l^m(x), as well as optional

* derivatives, for all L up to lmax and all M up to mmax.

* Returns results in array Plm, which should be dimensioned as double precision

* with size (0:10,0:10).

* If the first element of dPlm is not zero, then the first derivatives are also

* computed, and put into array dPlm.  To skip the derivatives it is only necessary

* to put a scalar 0.D0 in the dPlm parameter. Otherwise, dPlm should also be a

* double precision array of size (0:10,0:10).

* The recursion formulas keep a count of the exponents of the factor SQRT(1-x^2)

*  in both the numerator and denominator, which may cancel each other in the

*  final evaluation, particularly with the derivatives.  This prevents infinities

*  at x=-1 or +1. 

* If X is out of range ( abs(x)>1 ) then value is returns as if x=1.
       
        implicit real*8 (a-h,o-z)
	real*8           x,xx,Plm(0:10,0:10),P(0:10,0:10,0:1),fact,sfact

	real*8           dPlm(0:10,0:10),dP(0:10,0:10,0:2),anum,term

	LOGICAL DodPlm



	DodPlm=dPlm(0,0) .NE. 0.d0



	DO l=0,lmax

	    DO m=0,mmax

		  Plm(l,m)=0.D0

		  P(l,m,0)=0.D0

		  P(l,m,1)=0.D0

	    ENDDO

	ENDDO

	IF(lmax .LT. 0 .OR. mmax .LT. 0 .OR. mmax .GT. lmax )THEN

	  Print *,'Bad arguments to DLegendre00'

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



	IF(DodPlm)Then !do derivatives

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

			  dP(L,M,1)= dP(L2,M,1) - (2*L-1)*fact*dP(L-1,M-1,0) - 

     $                  (2*L-1)*dP(L-1,M-1,2) + (2*L-1)*xx*P(L-1,M-1,0)

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

*	FUNCTION FSVal00(omega,MaxN,FSC)
*
** Fourier Series Value
*
** Return value of Sine/Cosine Fourier series for N terms up to MaxN
*
** at angle omega, given the F.S. coeficients in array FSC
*
*        implicit real*8 (a-h,o-z)
*	REAL*8 omega,FSC(0:1,0:*)
*
*	INTEGER MaxN,n
*
*	REAL*8 Y,theta
*
*	Y=0.
*
*	DO n=0,MaxN
*
*	  theta=omega*n
*
*	  Y=Y + FSC(0,n)*COS(theta) + FSC(1,n)*SIN(theta)
*
*	ENDDO
*
*	FSVal=Y
*
*	RETURN
*
*	END

************************ Copyright 1996,2001 Dan Weimer/MRC ***********************

	SUBROUTINE SetModel00(angle,Bt,Tilt,SWVel,SWDen,ALindex,UseAL)

*

* Calculate the complete set of spherical harmonic coeficients,

* given an aribitrary IMF angle (degrees from northward toward +Y),

* magnitude Bt (nT), dipole tilt angle (degrees), 

* solar wind velocity (km/sec), SWDen (#/cc),

* ALindex (nT), and Logical flag to use optional AL index.

*

* Sets the value of Coef and Boundfit in the common block SetW00Coef.

*
        implicit real*8 (a-h,o-z)
	REAL*8 angle,Bt,Tilt,SWVel,SWDen,ALindex

	LOGICAL First,UseAL



	DATA First/.TRUE./

	SAVE First

	INTEGER unit

	CHARACTER*15 cfile

	CHARACTER*30 Copyright

	PARAMETER (MJ=3,ML=4,MM=3,MN=2,MO=2)

	integer :: iMJ,iMO,iML,iMM,iMN

	REAL*4  CS( 0:MJ, 0:1 , 0:MO, 0:1 , 0:ML, 0:MM)

	REAL*4 BCS( 0:MJ, 0:1 , 0:MO, 0:1 , 0:MN)

	REAL*4  SS( 0:1 , 0:1 , 0:MO, 0:1 , 0:ML, 0:MM)

	REAL*4 BSS( 0:1 , 0:1 , 0:MO, 0:1 , 0:MN)

	REAL*8 XA(0:MJ),XB(0:MJ),FSC(0:1,0:4),PSS(0:1)

	REAL*8 Coef(0:1,0:5,0:5),BoundFit(0:1,0:5),pi

	real*8           dpi

	INTEGER*4 i,j,k,l,m,n,o

	INTEGER*4 maxj,MaxL,MaxM,MaxN,MaxO

	COMMON /AllW00Coefs/MaxJ,MaxO,CS,BCS,SS,BSS

	COMMON /SetW00Coef/Coef,BoundFit,pi,dpi,MaxL,MaxM,MaxN



	If(First)Then

	cfile='w2k.dat' !make sure correct ENDIAN type file is used.

	unit=9

c	OPEN(UNIT=unit,FILE=cfile,STATUS='OLD',form='UNFORMATTED')
	OPEN(UNIT=unit,FILE=cfile,STATUS='old',form='FORMATTED')


	READ(unit,*) Copyright
	
	
	
c       PRINT *,Copyright

	
	read(unit,*) Maxj,MaxL,MaxM,MaxN,MaxO


	If(maxj.NE.MJ .OR. MaxL.NE.ML .OR. MaxM.NE.MM .OR.

     $   MaxN.NE.MN .OR. MaxO.NE.MO)Then

		PRINT *,'Data File Error'

		STOP !Data file did not match sixe expected for arrays

	Endif   

	
	do iMJ=0,MJ
	   do i=0,1
	      do iMO=0,MO
		 do j=0,1
		    do iML=0,ML
		       do iMM=0,MM
			  read(unit,*) CS( iMJ, i , iMO, j , iML, iMM)
		       enddo
		    enddo
		 enddo
	      enddo
	   enddo
	enddo

	
	do iMJ=0,MJ
	   do i=0,1
	      do iMO=0,MO
		 do j=0,1
		    do iMN=0,MN
		       read(unit,*) BCS( iMJ, i , iMO, j , iMN)
		    enddo
		 enddo
	      enddo
	   enddo
	enddo

	
	do i=0,1
	   do j=0,1
	      do iMO=0,MO
		 do k=0,1
		    do iML=0,ML
		       do iMM=0,MM
			  read(unit,*) SS( i , j , iMO, k , iML, iMM)
		       enddo
		    enddo
		 enddo
	      enddo
	   enddo
	enddo
	

	
	do i=0,1
	   do j=0,1
	      do iMO = 0,MO
		 do k=0,1
		    do iMN=0,MN
		       read(unit,*) BSS( i , j , iMO, k , iMN)
		    enddo
		 enddo
	      enddo
	   enddo
	enddo

	CLOSE(unit)
	

	pi=2.*ASIN(1.)

	dpi=2.D0*DASIN(1.D0)

	First=.FALSE.

	Endif



	SinTilt=SIN(Tilt)

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

	FUNCTION BoundaryLat00(gmlt)

        implicit real*8 (a-h,o-z)
	REAL*8 gmlt

	REAL*8 Coef(0:1,0:5,0:5),BoundFit(0:1,0:5),pi

	real*8           dpi

        integer*4 MaxL,MaxM,MaxN
	COMMON /SetW00Coef/Coef,BoundFit,pi,dpi,MaxL,MaxM,MaxN

	BoundaryLat00 = FSVal(gmlt*pi/12.,MaxN,BoundFit)

	RETURN

	END

****************** Copyright 1996, 2001, Dan Weimer/MRC ***********************

	FUNCTION EpotVal00(gLAT,gMLT)

* Return the value of the electric potential in kV at

* corrected geomagnetic coordinates gLAT (degrees) and gMLT (hours).

*

* Must first call SetModel00 to set up the model coeficients for

* the desired values of Bt, IMF clock angle, Dipole tilt angle, SW Vel,

* number density, and (optional) AL index.

*
        implicit real*8 (a-h,o-z)
	REAL*8 gLAT,gMLT

	real*8           Phi,Z,O,x,ct,Phim

	real*8            Plm(0:10,0:10),OPlm(0:10,0:10)
	real*8           dPlm(0:10,0:10)


	REAL*8 Coef(0:1,0:5,0:5),BoundFit(0:1,0:5),pi

	real*8           dpi

        integer*4 MaxL,MaxM,MaxN
	COMMON /SetW00Coef/Coef,BoundFit,pi,dpi,MaxL,MaxM,MaxN

	blat = BoundaryLat00(gmlt)

	IF(glat .GT. blat)THEN

        Phi=DBLE(gMLT)*dpi/12.D0

	  colat=90.-glat

	  bcolat=90.-blat

	  x=DBLE(colat)*dpi/DBLE(bcolat)

	  DC=DBLE(Coef(0,0,0))

	  Z=DC

	  O=DC

	  ct=DCOS(x)
	  dPlm(:,:)=0.D0
	  CALL DLegendre00(ct,MaxL,MaxM,Plm,dPlm)

!Also find value at outer boundary at angle Pi, or cos(pi)=-1.

	  CALL DLegendre00(-1.D0,MaxL,MaxM,OPlm,dPlm)

	  DO l=1,MaxL

	    Z=Z +  Plm(l,0)*DBLE(Coef(0,l,0))

	    O=O + OPlm(l,0)*DBLE(Coef(0,l,0))

	    mlimit=MIN(l,MaxM)

	    DO m=1,mlimit

	      phim=phi*m

	      Z=Z + Plm(l,m)*(DBLE(Coef(0,l,m))*DCOS(phim) +

     $			    DBLE(Coef(1,l,m))*DSIN(phim) )

	      O=O +OPlm(l,m)*DBLE(Coef(0,l,m))

	    ENDDO

	  ENDDO

	  EpotVal00 = SNGL(Z-O)

	ELSE

	  EpotVal00 = 0.

	ENDIF

	RETURN

	END

************************ Copyright 1996,2001 Dan Weimer/MRC ***********************
