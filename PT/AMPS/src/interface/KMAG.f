C*********************************************************************************************
C*********************************************************************************************
	SUBROUTINE 	KMAG(TIME,EPOCH,R,THETA,PHI,RLT,BY_IMF,BZ_IMF,
     ~	BRM,BTM,BPM,Dp)
c	R,	THETA and PHI are radial distance, colatitude (in radians) and longitude (in radians)
c	of the spacecraft. The program returns, local time (LT, in hours) and the model
c	magnetic field in system 3 shperical coordinates.



	IMPLICIT REAL*8(A-H,O-Z)
	Dimension pvecin(3),pvecou(3)
	Dimension bvecin(3),bvecou(3)
	Real*8 time
	character*5 Epoch
	DIMENSION C(4)
	DATA NL/4/ !number of current sheet modes
	DATA dql /1.41/

C	Strength of tail modes


	DATA C/
     ~  -0.553754,
     ~   0.714655,
     ~   0.119748,
     ~   1.237026/




	DATA Dp0/0.017/ !nominal dynamic pressure


	DATA a1/9.696/
	DATA a2/0.2351/
	DATA dK/0.71516/ !for selfsimilar mp

	DATA dhD/2.0/ !nominal half thickness of the cs


	r0=a1*Dp0**(-a2)


	r0or=r0/(a1*Dp**(-a2))


	D=dhD/r0or !Half thickness of the CS


C
c	FIRST CALCULATE THE INTERNAL FIELD HIGHER THAN DIPOLE

      CALL  KRONIAN_HIGHER(3,R,THETA,PHI,BRI,BTI,BPI)
C Jon V. - To remove non-axisymmetric terms, zero these out
C and don't call KRONIAN_HIGHER
C      BRI = 0.d0
C      BTI = 0.d0
C      BPI = 0.d0

c debug
c      write(6,*) 'BI(r,t,p): ',BRI,BTI,BPI
c end debug

C
C	NOW CALCULATE THE AXISYMMETRIC DIPOLE AND ITS SHIELD FIELD
C
      XS3=R*DSIN(THETA)*DCOS(PHI)
      YS3=R*DSIN(THETA)*DSIN(PHI)
      ZS3=R*DCOS(THETA)
C
C	ROTATE TRAJECTORY INTO DIS COORDINATES
C
	pvecin(1)=xs3
	pvecin(2)=ys3
	pvecin(3)=zs3
	CALL KROT('S3C','DIS',pvecin,pvecou,time,Epoch)
	XDIS=pvecou(1)
	YDIS=pvecou(2)
	ZDIS=pvecou(3)
c	Calculate the mapped location
	call mapit(time,XDIS,YDIS,ZDIS,XMAP,YMAP,ZMAP,r0or,Epoch)

c debug
c      write(6,*) 'X,Y,ZMAP: ', XMAP, YMAP, ZMAP
c end debug
	RMAP=dsqrt(XMAP**2+YMAP**2+ZMAP**2)

	CALL checkIfInsideMpSATURN(RMAP,XMAP,is_inside_mp,Dp,a1,a2,dK)


	CALL dipole_shielded(time,XMAP,YMAP,ZMAP,BXD,BYD,BZD,r0or)
c debug
c      write(6,*) 'BXD,BYD,BZD: ', BXD,BYD,BZD
c end debug
C
c
c	NOW MAP IT AND ROTATE INTO SYSTEM III
C
	CALL Mapped_field(time,XMAP,YMAP,ZMAP,BXD,BYD,BZD,BXDMAP,BYDMAP,
     ~BZDMAP,r0or,Epoch)
c debug
c      write(6,*) 'BXDMAP,BYDMAP,BZDMAP: ',BXDMAP,BYDMAP,BZDMAP
c end debug
	bvecin(1)=BXDMAP
	bvecin(2)=BYDMAP
	bvecin(3)=BZDMAP
	CALL KROT('DIS','S3C',bvecin,bvecou,time,Epoch)
	BXDS3=bvecou(1)
	BYDS3=bvecou(2)
	BZDS3=bvecou(3)
      BRD=BXDS3*DSIN(THETA)*DCOS(PHI)+BYDS3*DSIN(THETA)*
     .DSIN(PHI)+BZDS3*DCOS(THETA)
      BTD=BXDS3*DCOS(THETA)*DCOS(PHI)+BYDS3*DCOS(THETA)*DSIN(PHI)-
     .BZDS3*DSIN(THETA)
      BPD=-BXDS3*DSIN(PHI)+BYDS3*DCOS(PHI)

C	NOW CALCULATE THE FIELD FROM CURRENT SHEET AND ITS SHIELDING
	call shielded_csheetfield(time,R,THETA,PHI,RLT,Brscs,Btscs,Bpscs,
     ~	NL,C,D,r0or,dql,Epoch)

c debug
c      write(6,*) 'Brscs,Btscs,Bpscs: ',Brscs,Btscs,Bpscs
c end debug

C     NOW get the influence of the IMF

	CALL getIMF_penetration(BY_IMF,BZ_IMF,BY_p,BZ_p)

c debug
c      write(6,*) 'BY_p,BZ_p: ',BY_p,BZ_p
c end debug
	bvecin(1)=0.0
	bvecin(2)=BY_p
	bvecin(3)=BZ_p
	CALL KROT('KSM','S3C',bvecin,bvecou,time,Epoch)
	BX_ps3=bvecou(1)
	BY_ps3=bvecou(2)
	BZ_ps3=bvecou(3)

	CALL CAR2SPH_MAG(BX_ps3,BY_ps3,BZ_ps3,BR_p,BT_p,BP_p,THETA,PHI)

c debug
c      write(6,*)
c      write(6,*) 'R,T,P: ', R, THETA, PHI
c      write(6,*) 'BI(r,t,p): ',BRI,BTI,BPI
c      write(6,*) 'BD(r,t,p): ',BRD,BTD,BPD
c      write(6,*) 'Bscs(r,t,p): ',Brscs,Btscs,Bpscs
c      write(6,*) 'B_p(r,t,p): ',Br_p,Bt_p,Bp_p
c end debug

c	Now add all the fields

	IF (is_inside_mp .eq. 1) THEN


		Brm=BRI+BRD + Brscs + Br_p
		Btm=BTI+BTD + Btscs + Bt_p
		Bpm=BPI+BPD + Bpscs + Bp_p


	!Brm=BRD + Brscs + Br_p
	!Btm=BTD + Btscs + Bt_p
	!Bpm=BPD + Bpscs + Bp_p


	!Brm=BRI
	!Btm=BTI
	!Bpm=BPI


!	Brm=Brscs
!	Btm=Btscs
!	Bpm=Bpscs
	ELSE

		Brm=Br_p
		Btm=Bt_p
		Bpm=Bp_p

	ENDIF

c debug
c      write(6,*) 'inside KMAG: BR,BT,BP: ', Brm,Btm,Bpm
c debug
	Return
	End
C*********************************************************************************
C*********************************************************************************
C
      SUBROUTINE KRONIAN_HIGHER (NM,R,T,F,BR,BT,BF)
C
C     BASED ON THE SUBROUTINE IGRF WRITTEN BY N.A. TSYGANENKO (1979)
C     MODIFIED BY KRISHAN KHURANA, JULY, 1996. AND   NOV. 2004.
C
C     CALCULATES COMPONENTS OF MAIN KRONIAN FIELD IN RIGHT HANDED SPHERICAL
C     COORD SYSTEM. BASED ON THE  SPHERICAL HARMONIC COEFFICIENTS GIVEN BY
C     ACUNA ET AL. [1983] (Z3 MODEL,  V1+V2 DATA)
C     MAXIMUM ORDER OF HARMONICS TAKEN INTO ACCOUNT (NOT MORE THAN 0RDER 3)
C
C
C     IT IS ASSUMED THAT THE TRAJECTORY IS IS IN RIGHT HANDED S III COORDINATES.
C     THE OUTPUT IS ALSO IN RTP (RH) COORDINATES.
C
C            INPUT:  NM (INTEGER)- MAXIMUM ORDER OF HARMONICS TAKEN
C                                  INTO ACCOUNT (NM.LE.12)
C
C                    R,T,F (REAL)- POSITION OF DESIRED FIELD VALUE IN
C                                  SPHERICAL JOVIGRAPHIC COORDINATE SYSTEM
C                                  (R IN PLANET RADII, COLATITUDE T AND
C                                   RH LONGITUDE F IN RADIANS)
C
C           OUTPUT: BR,BT,BF (REAL)- COMPONENTS OF THE INTERNAL PORTION
C                                    OF THE MAIN MAGNETIC FIELD IN
C                                    SPHERICAL S III COORD SYSTEM
C                                    (VALUES GIVEN IN GAMMA)
C
C***************************************************************************
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL BK,BM
      DIMENSION A(13),B(13),G(91),H(91),REC(91)
	DATA REC/ 91 * 1.0/


C      DATA G/0.,0.,0.,1449.,0,0.,2012.,
C     ~       53.,92.,-40.,
C     ~       81*0.0/

       DATA G/0.,0.,24.,1449.,-179.,-35.,2012.,
     ~       53.,92.,-40.,
     ~       81*0.0/

C      DATA H/0.,0.,0.,0.,43.5680,0.,0.,
C     ~       140.9273,-92.8876,68.2759,
C     ~       81*0.0/

       DATA H/0.,0.,-0.9110,0.,43.5680,-44.7762,0.,
     ~       140.9273,-92.8876,68.2759,
     ~       81*0.0/


      DATA FIRSTI/0.0/
C     WRITE(1,'(5F15.3)')(G(I),I=1,10)
C     WRITE(1,'(5F14.3)')(H(I),I=1,10)
      IF(FIRSTI.EQ.0.0) GO TO 1
      GO TO 6
 1    FIRSTI=1.0
      G(1)=0.
      H(1)=0.
      KNM=15
      DO 2 N=1,13
      N2=2*N-1
      N2=N2*(N2-2)
      DO 2 M=1,N
      MN=N*(N-1)/2+M
  2   REC(MN)=DBLE((N-M)*(N+M-2))/DBLE(N2)
      S=1.
      DO 5 N=2,13
      MN=N*(N-1)/2+1
      S=S*DBLE(2*N-3)/DBLE(N-1)
      G(MN)=G(MN)*S
      H(MN)=H(MN)*S
      P=S
      DO 5 M=2,N
      AA=1.
      IF (M.EQ.2) AA=2.
      P=P*DSQRT(AA*DBLE(N-M+1)/DBLE(N+M-2))
      MNN=MN+M-1
      G(MNN)=G(MNN)*P
  5   H(MNN)=H(MNN)*P
  6   IF(KNM.EQ.NM) GOTO 61
      KNM=NM
      K=KNM+1
 61   PP=1./R
      P=PP
      DO 7 N=1,K
      P=P*PP
      A(N)=P
  7   B(N)=P*N
      P=1.
      D=0.
      BBR=0.
      BBT=0.
      BBF=0.
      U=T
      CF=DCOS(F)
      SF=DSIN(F)
      C=DCOS(U)
      S=DSIN(U)
      BK=(S.LT.1.D-5)
      DO 12 M=1,K
      BM=(M.EQ.1)
      IF(BM) GOTO 8
      MM=M-1
      W=X
      X=W*CF+Y*SF
      Y=Y*CF-W*SF
      GOTO 9
  8   X=0.
      Y=1.
  9   Q=P
      Z=D
      BI=0.
      P2=0.
      D2=0.
      DO 11 N=M,K
      AN=A(N)
      MN=N*(N-1)/2+M
      E=G(MN)
      HH=H(MN)
      W=E*Y+HH*X
      IF (DABS(P2).LT.1.D-38) P2=0.0
      IF (DABS(Q).LT.1.D-38) Q=0.0
      BBR=BBR+B(N)*W*Q
      BBT=BBT-AN*W*Z
      IF(BM) GOTO 10
      QQ=Q
      IF(BK) QQ=Z
      BI=BI+AN*(E*X-HH*Y)*QQ
  10  XK=REC(MN)
      DP=C*Z-S*Q-XK*D2
      PM=C*Q-XK*P2
      D2=Z
      P2=Q
      Z=DP
      Q=PM
   11 CONTINUE
      D=S*D+C*P
      P=S*P
      IF(BM) GOTO 12
      BI=BI*MM
      BBF=BBF+BI
  12  CONTINUE
      BR=BBR
      BT=BBT
      IF(BK) GOTO 13
      BF=BBF/S
      GOTO 14
  13  IF(C.LT.0.) BBF=-BBF
      BF=BBF
  14  CONTINUE
      RETURN
      END
C
C

C*************************************************************************

	Subroutine dipole_shielded(time,x,y,z,Bx,By,Bz,r0or)
C           WRITTEN BY K. K. KHURANA     11/2004

C           x,y,z input position
C		  OUTPUT: mag .filed Bx,By,Bz at x,y,z
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION pvecin(3),pvecou(3)

	DATA B0x/0.0/
	DATA B0y/0.0/
	DATA B0z/21160.0/


c	We will  first calculate the field of the dipole



	call dipole(B0x,B0y,B0z,x,y,z,Bxd,Byd,Bzd)
c
c	Now calculate the perpendicular dipole shielding field
c
	if( (z .eq. 0) .and. (y .eq. 0)) Then
	phi=0.0
	go to 1
	end if
	phi=datan2(z,y)
    1	rho=y*dcos(phi)+z*dsin(phi)



	Call B_mp_perp(rho,phi,x,Brho1,Bphi1,Bx1,8,r0or)


	By1=Brho1*dcos(phi)-Bphi1*dsin(phi)
	Bz1=Brho1*dsin(phi)+Bphi1*dcos(phi)



	Bx=Bxd -Bx1
	By=Byd -By1
	Bz=Bzd -Bz1



	Return
	End
c******************************************************************
c
	Subroutine dipole(B0x,B0y,B0z,x,y,z,Bx,By,Bz)
	IMPLICIT REAL*8(A-H,O-Z)

c	Calculates the field of Jupiter's dipole for shielding in the magnetopause
c	(B0x, B0y, B0z) is the dipole moment
c	x, y, z is the position vector
c     Bx, By, Bz is the output field vector
c
	r=dsqrt(x**2+y**2+z**2)
	a11=(3*x**2-r**2)/r**5
	a12=(3*x*y)/r**5
	a13=(3*x*z)/r**5
	a21=(3*x*y)/r**5
	a22=(3*y**2-r**2)/r**5
	a23=(3*y*z)/r**5
	a31=(3*x*z)/r**5
	a32=(3*y*z)/r**5
	a33=(3*z**2-r**2)/r**5
	Bx=a11*B0x+a12*B0y+a13*B0z
	By=a21*B0x+a22*B0y+a23*B0z
	Bz=a31*B0x+a32*B0y+a33*B0z
	Return
	End
C*************************************************************************
	Subroutine B_mp_perp(rho,phi,x,Bperpr,Bperpf,Bperpx,M,r0or)
	IMPLICIT REAL*8(A-H,O-Z)
	Dimension a(8),b(8)
C

	Data a/3.32125816385087D-004, -0.784083155434928,
     ~	-0.416208789540178, -3.60127324820496D-002, -0.273188407067209,
     ~   0.736497845500708, -0.939713175408542,0.433616725727916/
      Data b/9.24095494266823, 30.4656008995161,
     ~	   18.1008047927882, 60.6014301282221,  125.347578447821,
     ~	   253.707136094045, 509.877831023220, 1021.96079781869/


	Bperpr=0.
	Bperpf=0.
	Bperpx=0.
	Do 4 K=1,M
	Term1=Dsin(phi)*Dexp(x/(b(K)/r0or))
	Term2=bessj1(rho/(b(K)/r0or))/(rho/(b(K)/r0or))
     ~	-bessj0(rho/(b(K)/r0or))
	Term3=-Dcos(phi)*Dexp(x/(b(K)/r0or))
	Term4=bessj1(rho/(b(K)/r0or))/(rho/(b(K)/r0or))
	Term5=-Dsin(phi)*Dexp(x/(b(K)/r0or))
	Term6=bessj1(rho/(b(K)/r0or))
    	Bperpr=Bperpr+a(K)*r0or**3*Term1*Term2
	Bperpf=Bperpf+a(K)*r0or**3*Term3*Term4
      Bperpx=Bperpx+a(K)*r0or**3*Term5*Term6
    4 continue

      Return
	End
C***************************************************************
      Function DBSJ2(x)
	IMPLICIT REAL*8(A-H,O-Z)
	DBSJ2=(2.0/x)*bessj1(x)-bessj0(x)
	Return
	End
C***********************************************************************

      Function DBSJ3(x)
	IMPLICIT REAL*8(A-H,O-Z)
	DBSJ3=(4.0/x)*DBSJ2(x)-bessj1(x)
	Return
	End




C***********************************************************************

C     Hannes needs his own Bessl functions

C***********************************************************************


C     Taken from: http://audiolab.uwaterloo.ca/~jeffb/thesis/node50.html
	function bessj1(x)
      IMPLICIT REAL*8(A-H,O-Z)
      data r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,
     #242396853.1d0,
     #-2972611.439d0,15704.48260d0,-30.160366606d0/,
     #s1,s2,s3,s4,s5,s6/144725228442.d0,2300535178.d0,
     #18583304.74d0,99447.43394d0,376.9991397d0,1.d0/
      data p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,
     #.2457520174d-5,
     #-.240337019d-6/, q1,q2,q3,q4,q5
     #/.04687499995d0,-.2002690873d-3,
     #.8449199096d-5,-.88228987d-6,.105787412d-6/
      if(dabs(x).lt.8) then
	    y=x**2
	    bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))
     #      /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
      else
	    ax=dabs(x)
	    z=8./ax
	    y=z**2
	    xx=ax-2.356194491
	    bessj1=dsqrt(.636619772/ax)*(dcos(xx)*(p1+y*(p2+y*(p3+y*
     #      (p4+y*p5))))-z*dsin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
     #      *dsign(dble(1.),x)
      endif
      return
      end

C***********************************************************************

C     Taken from: http://audiolab.uwaterloo.ca/~jeffb/thesis/node50.html
      function bessj0(x)
      IMPLICIT REAL*8(A-H,O-Z)
      data p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,
     #-.2073370639d-5,.2093887211d-6/
     #,q1,q2,q3,q4,q5/-.1562499995d-1,
     #.1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
      data r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,
     #651619640.7d0,
     #-11214424.18d0,77392.33017d0,-184.9052456d0/,
     #s1,s2,s3,s4,s5,s6/57568490411.d0,1029532985.d0,
     #9494680.718d0,59272.64853d0,267.8532712d0,1.d0/
      if(abs(x).lt.8) then
	    y=x**2
	    bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))
     #      /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
      else
	    ax=dabs(x)
	    z=8./ax
	    y=z**2
	    xx=ax-.785398164
	    bessj0=dsqrt(.636619772/ax)*(dcos(xx)*(p1+y*(p2+y*(p3+
     #      y*(p4+y*p5))))-z*dsin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      return
      end

C**********************************************************************
c
	Subroutine KROT(From,To,Vecin,Vecout,time,epoch)
c
c	INPUTS: From,	a chracter*3 variable denoting the incoming coordinate system
c			To:		a chracter*3 variable denoting the outgoing coordinate system
c			Vecin	a variable of dimension 3 containing the incoming vector
c			Vecout:	a variable of dimension 3 containing the outgoing vector
c			time	A double precision variable denoting number of seconds from an epoch
c			epoch   A five letter character which can be either 'ctime', 'etime' or 'J2000'/upper or lower case
c
c
C	Rotates the vector vecin from the coordinate system "From" to the
c	system "To". The rotated vector is placed in the vector Vecout.
c
c	We always first rotate the Vecin in "from" system to system III. Then we rotate
c	it to the "To"coordinate system
c
c	The supported coordinate systems are:
c			S3C System III Cartesian (right-handed)
c			KSO Kronian  -Sun-Orbital
c			KSM Kronian-Sun-Magnetic (So far same as KSO)
c	        DIP Dipole (cartesian)
	Implicit Real*8(A-H,O-Z)
	Dimension Vecin(3),Vecout(3),vector(3)
	Character*3 From, To
	Character*5 epoch
	Real*8 dipole(3,3),first(3,3),second(3,3),dummy(3,3)
c	Real*4 stheta,sphi
	Real*8 time,jtime,eetime
	Parameter (PI=3.1415927,twopi=2.*PI,radian=PI/180.,degree=180./PI)
c      Parameter (TH_DIP=0.0,PH_DIP=0.0) !The dipole parameters
	Data dipole(3,1),dipole(3,2)/0.,0./
	Data dipole(3,3)/1.0/
	Data dipole(2,1),dipole(2,2),dipole(2,3)/0.,1.0,0/
	data dipole(1,1),dipole(1,2)/1.0,0./
	data dipole(1,3)/0.0/
!        write(*,*) time
!	write(6,*)"dipole: ",dipole

c	First ascertain the time epoch

	if((epoch(1:5).eq."j2000").or.(epoch(1:5).eq."J2000")) then
	jtime=time
	Go to 101
	end if


	if((epoch(1:5).eq."i2000").or.(epoch(1:5).eq."I2000")) then
	jtime=time+32.184D0
	Go to 101
	end if


	  if(((epoch(1:1) .eq. "c") .or. epoch(1:1) .eq. "C") .and.
     ~     ((epoch(2:2) .eq. "t") .or. epoch(2:2) .eq. "T") .and.
     ~     ((epoch(3:3) .eq. "i") .or. epoch(3:3) .eq. "I") .and.
     ~     ((epoch(4:4) .eq. "m") .or. epoch(4:4) .eq. "M") .and.
     ~     ((epoch(5:5) .eq. "e") .or. epoch(5:5) .eq. "E")) then
           jtime=eetime(time)
	      go to 101
	  end if
c
	  if(((epoch(1:1) .eq. "e") .or. epoch(1:1) .eq. "E") .and.
     ~     ((epoch(2:2) .eq. "t") .or. epoch(2:2) .eq. "T") .and.
     ~     ((epoch(3:3) .eq. "i") .or. epoch(3:3) .eq. "I") .and.
     ~     ((epoch(4:4) .eq. "m") .or. epoch(4:4) .eq. "M") .and.
     ~     ((epoch(5:5) .eq. "e") .or. epoch(5:5) .eq. "E")) then
		 jtime=time
	     go to 101
	  end if
	Write(*,*)" I do not understand the time epoch."
	Go to 100

c
c	Next ascertain the incoming coordinate system
c	Is it system 3?
  101	If(((From(1:1) .eq. "s") .or. from(1:1) .eq. "S") .and.
     ~   ((From(2:2) .eq. "3") .or. from(2:2) .eq. "3") .and.
     ~   ((From(3:3) .eq. "c") .or. from(3:3) .eq. "C")) Go to 1
c	Is it JSO?
	If(((From(1:1) .eq. "k") .or. from(1:1) .eq. "K") .and.
     ~   ((From(2:2) .eq. "s") .or. from(2:2) .eq. "S") .and.
     ~   ((From(3:3) .eq. "o") .or. from(3:3) .eq. "O")) Go to 2
c	Is it JSM?
	If(((From(1:1) .eq. "k") .or. from(1:1) .eq. "K") .and.
     ~   ((From(2:2) .eq. "s") .or. from(2:2) .eq. "S") .and.
     ~   ((From(3:3) .eq. "m") .or. from(3:3) .eq. "M")) Go to 3
c	Is it Dipole?
	If(((From(1:1) .eq. "d") .or. from(1:1) .eq. "D") .and.
     ~   ((From(2:2) .eq. "i") .or. from(2:2) .eq. "I") .and.
     ~   ((From(3:3) .eq. "p") .or. from(3:3) .eq. "P")) Go to 4
c	Is it Dipole Sun?
	If(((From(1:1) .eq. "d") .or. from(1:1) .eq. "D") .and.
     ~   ((From(2:2) .eq. "i") .or. from(2:2) .eq. "I") .and.
     ~   ((From(3:3) .eq. "s") .or. from(3:3) .eq. "S")) Go to 5

	Write(*,*)" I do not understand the incoming coordinate system"
	Go to 100
    1 continue
c	Write(*,*) " Incoming coordinate system is S3 Cartesian"
	First(1,1)=1.0
	First(1,2)=0.0
	First(1,3)=0.0
	First(2,1)=0.0
	First(2,2)=1.0
	First(2,3)=0.0
	First(3,1)=0.0
	First(3,2)=0.0
	First(3,3)=1.0
	Go to 60
    2 continue
c	Write(*,*) " Incoming coordinate system is KSO"
c	Calculate the rotation matrix to go from KSO to S3R
	call KSun(jtime,stheta,sphi,Ztheta,Zphi)
c	Define X component in system III of XKSO unit vector etc.
	dummy(1,1)=cos(stheta)*cos(sphi)
	dummy(1,2)=cos(stheta)*sin(sphi)
	dummy(1,3)=sin(stheta)
C	Calculate the Z axis of the KSO coordinate system.
c
c	Define X component in system III of ZKSO unit vector etc.
	dummy(3,1)=cos(Ztheta)*cos(Zphi)
	dummy(3,2)=cos(Ztheta)*sin(Zphi)
	dummy(3,3)=sin(Ztheta)
c
c	Define X component in system III of YKSO unit vector etc.

	dummy(2,1)=dummy(3,2)*dummy(1,3)-dummy(3,3)*dummy(1,2)
	dummy(2,2)=dummy(3,3)*dummy(1,1)-dummy(3,1)*dummy(1,3)
	dummy(2,3)=dummy(3,1)*dummy(1,2)-dummy(3,2)*dummy(1,1)
c	The transpose of the dummy matrix takes KSO to System 3.
	Do 21 i=1,3
	Do 21 j=1,3
   21	first(i,j)=dummy(j,i)
	go to 60

    3 continue
c	Write(*,*) " Incoming coordinate system is KSM"

	call KSun(jtime,stheta,sphi,Ztheta,Zphi)

c	Now define the KSM transpose vector
c	Define X component in system III of XJSM unit vector etc.

	dummy(1,1)=cos(stheta)*cos(sphi)
	dummy(1,2)=cos(stheta)*sin(sphi)
	dummy(1,3)=sin(stheta)
c	Now define the Y vector so that it is perpendicular to the dipole vector and X
	dummy(2,1)=dummy(1,3)*dipole(3,2)-dummy(1,2)*dipole(3,3)
	dummy(2,2)=dummy(1,1)*dipole(3,3)-dummy(1,3)*dipole(3,1)
	dummy(2,3)=dummy(1,2)*dipole(3,1)-dummy(1,1)*dipole(3,2)
	denom=sqrt(dummy(2,1)**2+dummy(2,2)**2+dummy(2,3)**2)
	dummy(2,1)=dummy(2,1)/denom
	dummy(2,2)=dummy(2,2)/denom
	dummy(2,3)=dummy(2,3)/denom
c	Now define the z vector
	dummy(3,1)=dummy(1,2)*dummy(2,3)-dummy(1,3)*dummy(2,2)
	dummy(3,2)=dummy(1,3)*dummy(2,1)-dummy(1,1)*dummy(2,3)
	dummy(3,3)=dummy(1,1)*dummy(2,2)-dummy(1,2)*dummy(2,1)
	Do 30 i=1,3
	Do 30 j=1,3
   30 first(i,j)=dummy(j,i)
	Go to 60
c

    4 continue
c	Write(*,*) " Incoming coordinate system is DIP"
	Do 40 i=1,3
	Do 40 j=1,3
   40 first(i,j)=dipole(j,i)
	go to 60
c
    5 continue
c	Write(*,*) " Incoming coordinate system is DIS"

c	Define  ZDIS unit vector
	dummy(3,1)=	dipole(3,1)
	dummy(3,2)=	dipole(3,2)
	dummy(3,3)=	dipole(3,3)

	call KSun(jtime,stheta,sphi,Ztheta,Zphi)


c	Define  XKSO unit vector
	dummy(1,1)=cos(stheta)*cos(sphi)
	dummy(1,2)=cos(stheta)*sin(sphi)
	dummy(1,3)=sin(stheta)
c
C	Define YDIS unit vector	(as normal to ZDIS and XKSO)
	dummy(2,1)=dummy(1,3)*dummy(3,2)-dummy(1,2)*dummy(3,3)
	dummy(2,2)=dummy(1,1)*dummy(3,3)-dummy(1,3)*dummy(3,1)
	dummy(2,3)=dummy(1,2)*dummy(3,1)-dummy(1,1)*dummy(3,2)
	denom=sqrt(dummy(2,1)**2+dummy(2,2)**2+dummy(2,3)**2)
	dummy(2,1)=dummy(2,1)/denom
	dummy(2,2)=dummy(2,2)/denom
	dummy(2,3)=dummy(2,3)/denom

c	Define XDIS unit vector as normal to YDIS and ZDIS
c
	dummy(1,1)=dummy(2,2)*dummy(3,3)-dummy(2,3)*dummy(3,2)
	dummy(1,2)=dummy(2,3)*dummy(3,1)-dummy(2,1)*dummy(3,3)
	dummy(1,3)=dummy(2,1)*dummy(3,2)-dummy(2,2)*dummy(3,1)


c	The transpose of the dummy matrix takes DIS to System 3.
	Do 22 i=1,3
	Do 22 j=1,3
   22	first(i,j)=dummy(j,i)


c
c
   60 continue
c
c	Now Figure out the outgoing coordinate system

	If(((To(1:1) .eq. "s") .or. To(1:1) .eq. "S") .and.
     ~   ((To(2:2) .eq. "3") .or. To(2:2) .eq. "3") .and.
     ~   ((To(3:3) .eq. "c") .or. To(3:3) .eq. "C")) Go to 6
c	Is it JSO?
	If(((To(1:1) .eq. "k") .or. To(1:1) .eq. "K") .and.
     ~   ((To(2:2) .eq. "s") .or. To(2:2) .eq. "S") .and.
     ~   ((To(3:3) .eq. "o") .or. To(3:3) .eq. "O")) Go to 7
c	Is it JSM?
	If(((To(1:1) .eq. "k") .or. To(1:1) .eq. "K") .and.
     ~   ((To(2:2) .eq. "s") .or. To(2:2) .eq. "S") .and.
     ~   ((To(3:3) .eq. "m") .or. To(3:3) .eq. "M")) Go to 8
c	Is it Dipole?
	If(((To(1:1) .eq. "d") .or. To(1:1) .eq. "D") .and.
     ~   ((To(2:2) .eq. "i") .or. To(2:2) .eq. "I") .and.
     ~   ((To(3:3) .eq. "p") .or. To(3:3) .eq. "P")) Go to 9
c	Is it DIS?
	If(((To(1:1) .eq. "d") .or. To(1:1) .eq. "D") .and.
     ~   ((To(2:2) .eq. "i") .or. To(2:2) .eq. "I") .and.
     ~   ((To(3:3) .eq. "s") .or. To(3:3) .eq. "S")) Go to 10
c
	Write(*,*)" I do not understand the outgoing coordinate system"
	Go to 100

    6	Continue
c	Write(*,*) " The outgoing coordinate system is S3 cartesian"
	Second(1,1)=1.0
	Second(1,2)=0.0
	Second(1,3)=0.0
	Second(2,1)=0.0
	Second(2,2)=1.0
	Second(2,3)=0.0
	Second(3,1)=0.0
	Second(3,2)=0.0
	Second(3,3)=1.0
	Go to 70
c
    7 continue
c	Write(*,*) " The outgoing coordinate system is KSO"
c	Get the matrix that Rotates from System III cartesian into KSO.
	call KSun(jtime,stheta,sphi,Ztheta,Zphi)
c	Define X component in system III of XKSO unit vector etc.
	Second(1,1)=cos(stheta)*cos(sphi)
	Second(1,2)=cos(stheta)*sin(sphi)
	Second(1,3)=sin(stheta)
c
c	Define X component in system III of ZKSO unit vector etc.
	Second(3,1)=cos(Ztheta)*cos(Zphi)
	Second(3,2)=cos(Ztheta)*sin(Zphi)
	Second(3,3)=sin(Ztheta)
c
c	Check orthogonality of the new coordinate system
	dotprod=second(3,1)*second(1,1)+second(3,2)*second(1,2)+
     .second(3,3)*second(1,3)
C	Write(*,*) dotprod
c	Define X component in system III of YKSO unit vector etc.

	Second(2,1)=Second(3,2)*Second(1,3)-Second(3,3)*Second(1,2)
	Second(2,2)=Second(3,3)*Second(1,1)-Second(3,1)*Second(1,3)
	Second(2,3)=Second(3,1)*Second(1,2)-Second(3,2)*Second(1,1)
c	seconm=	sqrt(second(2,1)**2+second(2,2)**2+second(2,3)**2)
c	Write(*,*)seconm
	Go to 70
    8 continue
c	Write(*,*) " The outgoing coordinate system is KSM"

c	Get the matrix that Rotates from System III cartesian into KSM.
	call KSun(jtime,stheta,sphi,Ztheta,Zphi)
c	Define X component in system III of XKSM unit vector etc.
	Second(1,1)=cos(stheta)*cos(sphi)
	Second(1,2)=cos(stheta)*sin(sphi)
	Second(1,3)=sin(stheta)
c	Now define the Y vector so that it is perpendicular to  the dipole vector and X
	Second(2,1)=second(1,3)*dipole(3,2)-second(1,2)*dipole(3,3)
	Second(2,2)=second(1,1)*dipole(3,3)-second(1,3)*dipole(3,1)
	Second(2,3)=second(1,2)*dipole(3,1)-second(1,1)*dipole(3,2)
	denom=sqrt(second(2,1)**2+second(2,2)**2+second(2,3)**2)
	Second(2,1)=second(2,1)/denom
	Second(2,2)=second(2,2)/denom
	Second(2,3)=second(2,3)/denom
c	Now define the z vector
	Second(3,1)=second(1,2)*second(2,3)-second(1,3)*second(2,2)
	Second(3,2)=second(1,3)*second(2,1)-second(1,1)*second(2,3)
	Second(3,3)=second(1,1)*second(2,2)-second(1,2)*second(2,1)
	Go to 70
c
    9 continue
c	Write(*,*) " The outgoing coordinate system is Dipole"
	Do 80 i=1,3
	Do 80 j=1,3
   80	Second(i,j)=dipole(i,j)
	go to 70
   10 Continue
c	Write(*,*) " Outgoing coordinate system is DIS"

c	Define  ZDIS unit vector
	dummy(3,1)=	dipole(3,1)
	dummy(3,2)=	dipole(3,2)
	dummy(3,3)=	dipole(3,3)


	!write(6,*)"DUMMY1: ",dummy

	call KSun(jtime,stheta,sphi,Ztheta,Zphi)
c	Define  XKSO unit vector
	dummy(1,1)=cos(stheta)*cos(sphi)
	dummy(1,2)=cos(stheta)*sin(sphi)
	dummy(1,3)=sin(stheta)
c

	!write(6,*)"stheta: ",stheta
	!write(6,*)"sphi: ",sphi

	!write(6,*)"DUMMY2: ",dummy

C	Define YDIS unit vector	(as normal to ZDIS and XKSO)
	dummy(2,1)=dummy(1,3)*dummy(3,2)-dummy(1,2)*dummy(3,3)
	dummy(2,2)=dummy(1,1)*dummy(3,3)-dummy(1,3)*dummy(3,1)
	dummy(2,3)=dummy(1,2)*dummy(3,1)-dummy(1,1)*dummy(3,2)
	denom=sqrt(dummy(2,1)**2+dummy(2,2)**2+dummy(2,3)**2)
	dummy(2,1)=dummy(2,1)/denom
	dummy(2,2)=dummy(2,2)/denom
	dummy(2,3)=dummy(2,3)/denom

      !write(6,*)"DUMMY3: ",dummy


c	Define XDIS unit vector as normal to YDIS and ZDIS
c
	dummy(1,1)=dummy(2,2)*dummy(3,3)-dummy(2,3)*dummy(3,2)
	dummy(1,2)=dummy(2,3)*dummy(3,1)-dummy(2,1)*dummy(3,3)
	dummy(1,3)=dummy(2,1)*dummy(3,2)-dummy(2,2)*dummy(3,1)


c	The  dummy matrix takes  System 3 to DIS.
	Do 23 i=1,3
	Do 23 j=1,3
   23	second(i,j)=dummy(i,j)
c
   70	continue
c	Now multimply vecin with first and second matrices to get the vecout
c	Write(*,*)" calling matmult 1"
	call matmult(first,vecin,vector,3,3,1)
c	Write(*,*)" calling matmult 2"
	call matmult(second,vector,vecout,3,3,1)

	!write(6,*)"dipole: ",dipole
	!write(6,*)"DUMMY: ",dummy
	!write(6,*)"FIRST: ",first
	!write(6,*)"sECOND: ",second



c
  100 Continue
	Return
	End
c***********************************************************************

c***********************************************************************
	Subroutine KSun(time,stheta,sphi,Ztheta,Zphi)
c	INPUT:  J2000 time of the data point
C	OUTPUTS: stheta, sphi, latitude and longitude  (in radians) of the Sun in system III (RH).
C	OUTPUTS: Ztheta, Zphi, latitude and longitude  (in radians) of the Zkso in system III (RH).
c	The equations are written in J2000 convention followed by PDS.
c

c	Last updated August 26, 2004.
c
c	theta=a1*cos(omegay*t)+a2*sin(omegay*t)+
c	a3*cos(2.*omegay*t)+a4*sin(2.*omegay*t)+
c	a5*cos(3.*omegay*t)+a6*sin(3.*omegay*t)+
c	a7*cos(4.*omegay*t)+a8*sin(4.*omegay*t)+
c     a9*t**2+a10*t+a11
c	fphi=b1*cos(omegay*t)+b2*sin(omegay*t)+
c     b3*cos(2.*omegay*t)+b4*sin(2.*omegay*t)+
c     b5*cos(3.*omegay*t)+b6*sin(3.*omegay*t)+
c     b7*cos(4.*omegay*t)+b8*sin(4.*omegay*t)+
c	b9*t**2+b10*t+b11	(I assume I have despun Saturn first.)

c	Then we rotate into the System III coordinates

c	time is a double precision variable and is assumed to be in J2000. theta and phi are single precision
c	variables.
c

c
c	omega is Saturn's rotation rate
c	omegaz is saturn's rotation rate for use of Zkso lat and lon
c	omegay is Saturn's yearly orbital rate

c	Initialize variables
	Implicit Real*8(A-H,O-Z)
	Real*8 aa(11),bb(11), cc(11),x(11)
	Real*8 time,jtime1,fphi,yrsat,omega,D360,omegay,t,year,omegaz
c	Real*4 stheta,sphi,ztheta,Zphi
	Parameter (PI=3.141592654d0)
	PARAMETER (twopi=2.*PI,radian=PI/180.,degree=180./PI)
	Parameter (yrsat=10759.22d0*86400.0d0)
	Parameter (omega=360.D0/(10.65622222D0*3600.D0))
	PARAMETER (omegay=2.D0*PI/yrsat)
	Parameter (omegaz=360.D0/(10.65622222D0*3600.D0)+3.880688963D-7)
	Parameter (year=86400D0*.36525D3,D360=360.D0)
	Parameter (jtime1=-.825767955816D+09)
	Parameter (Zthetd=63.271d0)
c	Parameter (Zthetd=-26.729)
	Data aa  /-.26237934D+02,.30525104D+01,-.13686733D+01,
     +.10182726D+00,-.30897805D+00,.89277949D-01,-.39638771D-01,
     +.97499653D-02,.40974159D-04,.14075684D-02,.13505775D+01/
c
	Data bb /-.32108124D+00,.58569866D+01,.71266272D+00,.33244539D+01,
     +.47529267D-01,.32770362D+00,.46935622D-01,.10720536D+00,
     + -.50594764D-03,.29439612D-01,.26423581D+03/
c
	Data cc/-.13200649D-02,-.18358429D-02,.64658927D-03,.12800487D-02,
     +.17618936D-04,-.58790898D-03,.49804081D-04,.42372137D-03,
     +.14744891D-04,.25369744D-01,.77943328D+02/

c	First calculate the latitude and longitude in non-rotating Saturn coordinates.

c	Calculate the best fit theta and fphi
	t=time-jtime1
	x(1)=cos(omegay*t)
	x(2)=sin(omegay*t)
	x(3)=cos(2.*omegay*t)
	x(4)=sin(2.*omegay*t)
	x(5)=cos(3.*omegay*t)
	x(6)=sin(3.*omegay*t)
	x(7)=cos(4.*omegay*t)
	x(8)=sin(4.*omegay*t)
	x(9)=(t/year)**2
	x(10)=t/year
	x(11)=1.0
	stheta=0.
c	fphi is phi in Saturn fixed (non-rotating) coordinate
	fphi=0.0
	Zfphi=0.0
	Do 1 j=1,11
	fphi=fphi+bb(j)*x(j)
	Zfphi=Zfphi+cc(j)*x(j)
    1 stheta=stheta+aa(j)*x(j)
C	Now rotate the longitude to Saturn System III
c	First Add the rotation of Saturn around the Sun.
c	fphi is the phi of the Sun as unspinning Saturn goes around the sun
	fphi=DMod(fphi+t/yrsat*360.d0, D360)
c	Next add the rotation of Saturn around its axis.
	sphi=DMod(fphi-t*omega, D360)
	if (sphi .lt. 0.0) sphi = sphi+360.0
	sphi=sphi*radian
	stheta=stheta*radian
c
C	Next rotate the longitude of Zfphi to Saturn System III
c	First Add the rotation of Saturn around the Sun.
c	Zfphi is the phi of Zkso as unspinning Saturn goes around the sun
	 Zfphi=DMod(Zfphi+t/yrsat*360.d0+180.D0, D360)
c	Next add the rotation of Saturn around its axis.
	Zphi=DMod(Zfphi-t*omegaz, D360)
	if (Zphi .lt. 0.0) Zphi = Zphi+360.0
	Zphi=Zphi*radian
	Ztheta=Zthetd*radian
	Return
      end
c*******************************************************************
	Function eetime(ctime)
	IMPLICIT real*8(A-H,O-z)
	Real*8 ctime, eetime
	Tcor=0.
	if(ctime .ge. 189302400.000) Tcor = Tcor+10
	if(ctime .ge. 205027200.000) Tcor = Tcor+1
	if(ctime .ge. 220924800.000) Tcor = Tcor+1
	if(ctime .ge. 252460800.000) Tcor = Tcor+1
	if(ctime .ge. 283996800.000) Tcor = Tcor+1
	if(ctime .ge. 315532800.000) Tcor = Tcor+1
	if(ctime .ge. 347155200.000) Tcor = Tcor+1
	if(ctime .ge. 378691200.000) Tcor = Tcor+1
	if(ctime .ge. 410227200.000) Tcor = Tcor+1
	if(ctime .ge. 441763200.000) Tcor = Tcor+1
	if(ctime .ge. 489024000.000) Tcor = Tcor+1
	if(ctime .ge. 520560000.000) Tcor = Tcor+1
	if(ctime .ge. 552096000.000) Tcor = Tcor+1
	if(ctime .ge. 615254400.000) Tcor = Tcor+1
	if(ctime .ge. 694224000.000) Tcor = Tcor+1
	if(ctime .ge. 757382400.000) Tcor = Tcor+1
	if(ctime .ge. 788918400.000) Tcor = Tcor+1
	if(ctime .ge. 836179200.000) Tcor = Tcor+1
	if(ctime .ge. 867715200.000) Tcor = Tcor+1
	if(ctime .ge. 899251200.000) Tcor = Tcor+1
	if(ctime .ge. 946684800.000) Tcor = Tcor+1
	if(ctime .ge. 993945600.000) Tcor = Tcor+1
	if(ctime .ge. 1041379200.000) Tcor = Tcor+1
	eetime=ctime+dble(Tcor)-.1072958367816D10
	Return
	End
C*************************************************************************


      SUBROUTINE matmult(xx,yy,zz,nxrow,nxcol,nycol)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION xx(nxrow,nxcol),yy(nxcol,nycol),zz(nxrow,nycol)
	Do 1 i=1,nxrow
	Do 1 j=1,nycol
	zz(i,j)=0.
	Do 1 k=1,nxcol
	zz(i,j)=zz(i,j)+xx(i,k)*yy(k,j)
    1 continue
      RETURN
      END
C******************************************************************
	Function ctimer(iyr,imon,iday,ihr,imin,sec)
	integer iyr,imon,iday,ihr,imin
	integer*8 ndays,doy
	Real*8 sec,ctimer

c	First calculate the number of days from Jan 1, 1966
		ndays=0
		If	(iyr .ge. 1966) then
		Do 1 i=1966,iyr-1
		ndays=ndays+365
		if(((mod(i,4) .eq. 0) .and. (mod(i,100) .ne. 0)) .or.
     +	(mod(i,400)) .eq. 0) ndays=ndays+1
   1		continue
c	Now add the number of days of the current year
    	ndays=ndays+doy(iyr,imon,iday)-1
    	ctimer = dble(ndays)*86400D0+dble(ihr)*3600D0+dble(imin)*60D0+
     +dble(sec)
	go to 4
	end if

c	Calculate the seconds for the negative years
	If	(iyr .lt. 1966) then
	Do 2 i=iyr, 1965
	ndays=ndays-365
		if(((mod(i,4) .eq. 0) .and. (mod(i,100) .ne. 0)) .or.
     +	(mod(i,400)) .eq. 0) ndays=ndays-1
    2 continue
c	Now subtract the number of days of the current year
	ndays=ndays+doy(iyr,imon,iday)
	end if
    	ctimer = dble(ndays-1)*86400D0+dble(ihr)*3600D0+dble(imin)*60D0+
     +dble(sec)
    4 continue
	idyr=doy(iyr,imon,iday)
	Write(*,*) " doy is ",idyr
	Write(*,*)" ndays is ",ndays
	Return
	End
c*******************************************************************
	Function doy(iyr,imon,iday)
	Integer*8 mon(12),doy
	data mon /31,28,31,30,31,30,31,31,30,31,30,31/
	doy=0
	Do 1 i = 2,imon
    	doy=doy+mon(i-1)
c	Add an extra day for February
	if(i .eq. 3) then
		if(((mod(iyr,4) .eq. 0) .and. (mod(iyr,100) .ne. 0)) .or.
     +	(mod(iyr,400)) .eq. 0) doy=doy+1
	end if

    1 continue
	doy=doy+iday
	Return
	End
C********************************************************************************
C********************************************************************************
	Subroutine mapit(time,XDIS,YDIS,ZDIS,XMAP,YMAP,ZMAP,r0or,Epoch)
	IMPLICIT REAL*8(A-H,O-Z)
	data RH0/23.0/
	Dimension pvecin(3),pvecou(3)
	character*5 Epoch
C	r=sqrt(XKSM**2+YKSM**2+ZKSM**2)

C	r0or=1.0

	RH=RH0/r0or !hinging as function of pressure
c	First calculate the dipole tilt (si) in the KSM coordinates
	pvecin(1)=0.
	pvecin(2)=0.
	pvecin(3)=1.
	CALL KROT('S3C','KSM',pvecin,pvecou,time,epoch)
c Jon V - comment from Krishan: si = Saturn inclination
	si=	dacos(pvecou(3))
	if(pvecou(1) .lt. 0) si=-si
c	Now calculate the height of the current sheet
	rhoDIS=dsqrt(XDIS**2+YDIS**2)
	phiLT=datan2(YDIS,XDIS)
	theta=datan2(ZDIS,rhoDIS)
	Zcs=(rhoDIS-RH*dtanh(rhoDIS/RH))*dtan(-si) !sign of si?
	Thcs=datan2(Zcs,rhoDIS) !sign of Thcs?
	rhomap= rhoDIS*dcos(thcs)+ZDIS*dsin(thcs)
	Zmap=  -rhoDIS*dsin(thcs)+ZDIS*dcos(thcs)
	Xmap=rhomap*dcos(phiLT)
	Ymap=rhomap*dsin(phiLT)
	Return
	End
C********************************************************************************
C********************************************************************************
	SUBROUTINE Mapped_field(
     ~	time,X,Y,Z,BX,BY,BZ,BXMAP,BYMAP,BZMAP,r0or,Epoch)
	IMPLICIT REAL*8(A-H,O-Z)
	character*5 Epoch
C	DIMENSION pvecin(3),pvecou(3)
c
C	r0or=1.0
c
c
c	Now define the nine derivatives
c	These are calculated at the original location
c
	dx=0.01
	dy=0.01
	dz=0.01
	xp=x+dx
	xm=x-dx
	yp=y+dy
	ym=y-dy
	zp=z+dz
	zm=z-dz
c
c	We begin with the x derivatives

	call mapit(time,xp,y,z,xpp,ypp,zpp,r0or,Epoch)
	call mapit(time,xm,y,z,xpm,ypm,zpm,r0or,Epoch)

	dxpdx=(xpp-xpm)/(2.*dx)
	dypdx=(ypp-ypm)/(2.*dx)
	dzpdx=(zpp-zpm)/(2.*dx)

c	Next the y derivatives
c

C	dxpdy=0.
C	dypdy=1.
C	dzpdy=0.


	call mapit(time,x,yp,z,xpp,ypp,zpp,r0or,Epoch)
	call mapit(time,x,ym,z,xpm,ypm,zpm,r0or,Epoch)


	dxpdy=(xpp-xpm)/(2.*dy)
	dypdy=(ypp-ypm)/(2.*dy)
	dzpdy=(zpp-zpm)/(2.*dy)


c	Next the z gradients

	call mapit(time,x,y,zp,xpp,ypp,zpp,r0or,Epoch)
	call mapit(time,x,y,zm,xpm,ypm,zpm,r0or,Epoch)

	dxpdz=(xpp-xpm)/(2.*dz)
	dypdz=(ypp-ypm)/(2.*dz)
	dzpdz=(zpp-zpm)/(2.*dz)


c
C	Now calculate the T matrix
c  	Calculate the mapped location
c
	Txx=dypdy*dzpdz-dypdz*dzpdy
	Txy=dxpdz*dzpdy-dxpdy*dzpdz
	Txz=dxpdy*dypdz-dxpdz*dypdy
	Tyx=dypdz*dzpdx-dypdx*dzpdz
	Tyy=dxpdx*dzpdz-dxpdz*dzpdx
	Tyz=dxpdz*dypdx-dxpdx*dypdz
	Tzx=dypdx*dzpdy-dypdy*dzpdx
	Tzy=dxpdy*dzpdx-dxpdx*dzpdy
	Tzz=dxpdx*dypdy-dxpdy*dypdx

C	Now calculate the field at the mapped location
c
      Bxmap=Txx*Bx+Txy*By+Txz*Bz
      Bymap=Tyx*Bx+Tyy*By+Tyz*Bz
      Bzmap=Tzx*Bx+Tzy*By+Tzz*Bz
	return
	End

C******************************************************************
	Subroutine shielded_csheetfield(time,R1,T1,P1,RLT,Brm1,Btm1,
     ~Bpm1,NL,C,D,r0or,dql,epoch)
      IMPLICIT REAL*8(A-H,O-Z)
c	This subroutine calculates the best fit current sheet by modeling only
c	Brho and Bz components. Once the coefficients for the 7 modes are known,
c	the other subroutines MAPIT and Mapped_field introduce the stretch transformations
c	and calculate the new values of Brho,Bphi and Bz.
c	Real*4 stheta,sphi,ztheta,Zphi

	Dimension pvecin(3),pvecou(3)
	Dimension bvecin(3),bvecou(3)
	Real*8 time
	character*5 Epoch
	Dimension f(6,6),beta(6,6),C(NL)
	Dimension Xrho(NL),Xz(NL)


C     RING CURRENT
C     ============



c     PEAK AT 6

C      Data (f (1, j), j = 1, 6) /
C     ~   -96480.9970060387,
C     ~	 -733.2997300981,
C     ~   -17474.3949850542,
C     ~	13878.3552501633,
C     ~	99937.8587320000,
C     ~	  988.4579332559/
C      Data (beta (1, j), j = 1, 6) /
C     ~	4.5216256036,
C     ~	2.5132995869,
C     ~	3.2135379771,
C     ~	3.0509806322,
C     ~	4.4629794036,
C     ~	6.6822375235/
C L=4
c     PEAK AT 4
	Data (f (1, j), j = 1, 6) /
     ~ 	   69605.0534869644,
     ~ 	   149.0772338679,
     ~ 	   78909.3792557924,
     ~ 	-54703.2221035385,
     ~ 	-93857.8599750019,
     ~ 	   -59.8779082335/
C bi
      Data (beta (1, j), j = 1, 6) /
     ~     3.7051342090,
     ~     2.3672348371,
     ~     3.2236662867,
     ~     3.1391528825,
     ~     3.6275112423,
     ~     7.8003148590/

      Data (f (2, j), j = 1, 6) /
     ~        -60.6568564967,
     ~     -73849.9534670401,
     ~      97442.8437444251,
     ~     -23178.2043425704,
     ~         66.3747390965,
     ~         -1.4848752731/
      Data (beta (2, j), j = 1, 6) /
     ~    48.1885954914,
     ~     5.7107159217,
     ~     5.6345236716,
     ~     5.3879744819,
     ~    89.0583961619,
     ~    22.6415533854/
C L=8
      Data (f (3, j), j = 1, 6) /
     ~      310.2640136545,
     ~    24690.5295574541,
     ~      786.0977568405,
     ~   -24122.7960738661,
     ~     -259.2663294173,
     ~     -598.5252855461/
      Data (beta (3, j), j = 1, 6) /
     ~  48.9404306370,
     ~   5.2832855659,
     ~  16.8568900868,
     ~   5.2434883472,
     ~  89.2258188514,
     ~  25.7902108924/
C L=16
C      Data (f (4, j), j = 1, 6) /
C     ~    -359.4857268465,
C     ~ -718318.9196872473,
C     ~     952.9729049448,
C     ~  718848.3348118553,
C     ~      29.9449422503,
C     ~    1422.1605650097/
C      Data (beta (4, j), j = 1, 6) /
C     ~    48.9404373087,
C     ~     5.4911324770,
C     ~    13.9191935499,
C     ~     5.4924655928,
C     ~    89.2258188514,
C     ~    25.7560973761/


C L=32
	Data (f (4, j), j = 1, 6) /
     ~     4285.5052331319,
     ~  -600859.2054973751,
     ~     1283.3234583067,
     ~   601269.3489620079,
     ~     -421.4303718400,
     ~     3031.9399845112/

	Data (beta (4, j), j = 1, 6) /
     ~   48.9404306370,
     ~    5.1544511385,
     ~   13.9220169737,
     ~    5.1557936291,
     ~   89.2258188514,
     ~   25.7556569433/

C L=32
	Data (f (5, j), j = 1, 6) /
     ~      4285.5052331319,
     ~   -600859.2054973751,
     ~      1283.3234583067,
     ~    601269.3489620079,
     ~      -421.4303718400,
     ~      3031.9399845112/
	Data (beta (5, j), j = 1, 6) /
     ~   48.9404306370,
     ~    5.1544511385,
     ~   13.9220169737,
     ~    5.1557936291,
     ~   89.2258188514,
     ~   25.7556569433/


C L=64
	Data (f (6, j), j = 1, 6) /
     ~    9938.2744323159,
     ~ -500764.5439833679,
     ~    1287.6877590235,
     ~  501147.4432669706,
     ~   14912.2836873664,
     ~    3963.5336376826/
	Data (beta (6, j), j = 1, 6) /
     ~    48.9368813776,
     ~     5.1508684107,
     ~    13.9644917853,
     ~     5.1524290713,
     ~    89.2282691400,
     ~    25.7556037045/

C L=64
c      f(6,1)=      9938.2744323159
c      f(6,2)=   -500764.5439833679
c      f(6,3)=      1287.6877590235
c      f(6,4)=    501147.4432669706
c      f(6,5)=     14912.2836873664
c      f(6,6)=      3963.5336376826

c      beta(6,1)=   48.9368813776
c      beta(6,2)=    5.1508684107
c      beta(6,3)=   13.9644917853
c      beta(6,4)=    5.1524290713
c      beta(6,5)=   89.2282691400
C      beta(6,6)=   25.7556037045



	DATA NI/6/
C     NL number of modes
C     NI mumber of coefficients per mode


c	We work in the KSM coordinate system

c	In this subroutine, index I=1-6 corresponds to each component of a mode
c	Index L = 1,7 denotes the seven L modes
c	Index j sums over observations (Br and Btheta) and is twice the number of
c	observations
c	index M is used to obtain least squares equations (same dimension as L)
      PI=4.0*datan(1.D0)
	twopi=2.*PI
	PI2=PI/2.
      DEGREE=180./PI
	RADIAN=PI/180.
	drho=0.05
	dz=0.05
c	D is the half thickness of current sheet

	M=5
	M1=25
	M2=2*M
C
C     END DEFINITIONS
C..........................................................................
c......................................................................
c debug
c      write(6,*) 'inputs:time,R1,T1,P1,RLT,NL,C,D,r0or,dql,epoch: ',
c     +time,R1,T1,P1,RLT,NL,C,D,r0or,dql,epoch
c end debug
	R_S3=R1
	THETA=T1
	PHI=P1
	rho=R_S3*DSIN(THETA)
	time=time
c	RLT=(posLT(j)-12.0)*15.0
c	If (RLT .LT. 0.0) RLT=RLT+360.0
c	RLT=RLT*Radian
      XS3=R_S3*DSIN(THETA)*DCOS(PHI)
      YS3=R_S3*DSIN(THETA)*DSIN(PHI)
      ZS3=R_S3*DCOS(THETA)
C
C   NOW ROTATE the trajectory INTO KSM COORDINATES
C

	!write(6,*)"Xs3: ",Xs3
	!write(6,*)"Ys3: ",Ys3
	!write(6,*)"Zs3: ",Zs3

      !write(6,*)"time: ",time,epoch


	pvecin(1)=xs3
	pvecin(2)=ys3
	pvecin(3)=zs3
	CALL KROT('S3C','DIS',pvecin,pvecou,time,epoch)
	XDIS=pvecou(1)
	YDIS=pvecou(2)
	ZDIS=pvecou(3)

	RLT=datan2(YDIS,XDIS)


c	Start calculating the field for each unit mode	at the mapped location
c	First map the trajectory

	!write(6,*)"XDIS: ",XDIS
	!write(6,*)"YDIS: ",YDIS
	!write(6,*)"ZDIS: ",ZDIS



	call mapit(time,XDIS,YDIS,ZDIS,XMAP,YMAP,ZMAP,r0or,Epoch)


	Z=ZMAP
	RHOMAG=sqrt(XMAP**2+YMAP**2)
c debug
c      write(6,*) 'XMAP,YMAP,ZMAP,RHOMAG: ',XMAP,YMAP,ZMAP,RHOMAG
c end debug

	Do 2 L = 1, NL

	ZM=dabs(Z-dz)
	If (ZM < D) ZM= 0.5*(ZM**2/D+D)
	ZP=dabs(Z+dz)
	If (ZP < D) ZP= 0.5*(ZP**2/D+D)
	xlpp=0.
	xlpm=0.
	Do 3 i=1,NI
	S1p=dsqrt((beta(L,i)/r0or+ZP)**2+(RHOMAG+Beta(L,i)/r0or)**2)
	S2p=dsqrt((beta(L,i)/r0or+ZP)**2+(RHOMAG-Beta(L,i)/r0or)**2)
	S1m=dsqrt((beta(L,i)/r0or+ZM)**2+(RHOMAG+Beta(L,i)/r0or)**2)
	S2m=dsqrt((beta(L,i)/r0or+ZM)**2+(RHOMAG-Beta(L,i)/r0or)**2)
	tp=2*(Beta(L,i)/r0or)/(S1p+S2p)
	tm=2*(Beta(L,i)/r0or)/(S1m+S2m)
	AAp=tp*dsqrt(1.-tp**2)/(S1p*S2p)
	AAm=tm*dsqrt(1.-tm**2)/(S1m*S2m)
	xlpp=xlpp+f(L,i)*AAp*rhomag
	xlpm=xlpm+f(L,i)*AAm*rhomag
c debug
c      write(6,*) 'S1p,S2p,S1m,S2m: ',
c     +  S1p,S2p,S1m,S2m
c      write(6,*) 'tp,tm,AAp,AAm,xlpp,xlpm: ',
c     +  tp,tm,AAp,AAm,xlpp,xlpm
c end debug
    3 continue
	dxpldz=(xlpp-xlpm)/(2.*dz)
	!write(6,*)"dxpldz: ",dxpldz
c debug
c	write(6,*)'L,i,dxpldz: ',L,i,dxpldz
c end debug
	Xrho(L)=-dxpldz
    2	continue

c	Now CALCULATE THE BZ COMPONENT
	Do 20 L = 1, NL
	rhom=RHOMAG-drho
	rhop=RHOMAG+drho

	xi=dabs(Z)
	If (dabs(Z) .le. D) xi= 0.5*(Z**2/D+D)

c
	xlpp=0.
	xlpm=0.
	Do 30 i=1,NI
	S1p=dsqrt((beta(L,i)/r0or+xi)**2+(rhop+Beta(L,i)/r0or)**2)
	S2p=dsqrt((beta(L,i)/r0or+xi)**2+(rhop-Beta(L,i)/r0or)**2)
	S1m=dsqrt((beta(L,i)/r0or+xi)**2+(rhom+Beta(L,i)/r0or)**2)
	S2m=dsqrt((beta(L,i)/r0or+xi)**2+(rhom-Beta(L,i)/r0or)**2)
	tp=2*(Beta(L,i)/r0or)/(S1p+S2p)
	tm=2*(Beta(L,i)/r0or)/(S1m+S2m)
	AAp=tp*dsqrt(1.-tp**2)/(S1p*S2p)
	AAm=tm*dsqrt(1.-tm**2)/(S1m*S2m)
	xlpp=xlpp+f(L,i)*AAp*rhop
	xlpm=xlpm+f(L,i)*AAm*rhom
   30 continue
	dxpldr=(rhop*xlpp-rhom*xlpm)/(2.*drho)
	!write(6,*)"dxpldr: ",dxpldr
c debug
c	write(6,*)"L,i,dxpldr: ",L,i,dxpldr
c end debug
	Xz(L)=dxpldr/RHOMAG
	!write(6,*)"RHOMAG: ",RHOMAG
c debug
c	write(6,*)"L,i, RHOMAG: ",L,i,RHOMAG
c end debug
   20	continue


c
c	Now calculate the new mapped field for each  UNIT mode
	Bx=0.0
	By=0.0
	Bz=0.0

C	Do 60 L=5,5
	Do 60 L=1,NL
	if(( ymap .eq. 0.) .and. (xmap .eq. 0.)) then
	phimap=0.
	go to 69
	end if
	phimap=atan2(YMAP,XMAP)
   69	Bx1=Xrho(L)*Dcos(phimap)
	By1=Xrho(L)*Dsin(phimap)
	Bz1=Xz(L)

	!write(6,*)"Bx1: ",Bx1
	!write(6,*)"By1: ",By1
	!write(6,*)"Bz1: ",Bz1



c	Now Add the ACTUAL shield field to the current sheet field
	XMAP1=XMAP
	YMAP1=YMAP
	ZMAP1=ZMAP
	call shieldfield_csheet(XMAP1,YMAP1,ZMAP1,L,M,M1,M2,
     ~BX2,BY2,BZ2,r0or)


!      BX2=0.0
!	BY2=0.0
!	BZ2=0.0




	Bx=Bx+C(L)*r0or**dql*(Bx1*(1/(r0or))**2.0-Bx2/r0or)
	By=By+C(L)*r0or**dql*(By1*(1/(r0or))**2.0-By2/r0or)
	Bz=Bz+C(L)*r0or**dql*(Bz1*(1/(r0or))**2.0-Bz2/r0or)


	!write(6,*)"Bx: ",Bx
	!write(6,*)"By: ",By
	!write(6,*)"Bz: ",Bz


   60 Continue
c
	!write(6,*)"XMAP: ",XMAP
	!write(6,*)"YMAP: ",YMAP
	!write(6,*)"ZMAP: ",ZMAP


	call Mapped_field(
     ~	time,XMAP,YMAP,ZMAP,Bx,By,Bz,Bxdis,Bydis,Bzdis,r0or,Epoch)
c
	bvecin(1)=Bxdis
	bvecin(2)=Bydis
	bvecin(3)=Bzdis

	!write(6,*)"Bxdis: ",Bxdis
	!write(6,*)"Bydis: ",Bydis
	!write(6,*)"Bzdis: ",Bzdis

c	Now rotate back into spherical coordinate system


	CALL KROT('DIS','S3C',bvecin,bvecou,time,epoch)

	Bxs3=bvecou(1)
	Bys3=bvecou(2)
	Bzs3=bvecou(3)

	!write(6,*)"Bxs3: ",Bxs3
	!write(6,*)"Bys3: ",Bys3
	!write(6,*)"Bzs3: ",Bzs3

	Brm1=Bxs3*sin(theta)*cos(phi)+Bys3*sin(theta)*sin(phi)+
     ~Bzs3*cos(theta)
	Btm1=Bxs3*cos(theta)*cos(phi)+Bys3*cos(theta)*sin(phi)-
     ~Bzs3*sin(theta)
	Bpm1=-Bxs3*sin(phi)+Bys3*cos(phi)

!	!write(6,*)"Brm1: ",Brm1
!	!write(6,*)"Btm1: ",Btm1
!	!write(6,*)"Bpm1: ",Bpm1

	Return
	End

C********************************************************************
	Subroutine shieldfield_csheet(X,Y,Z,L,M,M1,M2,
     ~Bxm,Bym,Bzm,r0or)
	IMPLICIT REAL*8(A-H,O-Z)
c

c
	Dimension ABCD(6,25),PQ(6,10)
	Dimension a(M,M),p(M2)

	Dimension b(M,M),c(M,M),d(M,M)
	Dimension q(M2),r(M2),s(M2)
	Dimension XXP(M1),XXM(M1)


C	MODE 1 Peak at 6
C	Data (abcd (1, j), j = 1, 25) /
C     ~    0.57900389912988492469D-05,
C     ~   -0.33867095260193265104D-02,
C     ~   -0.17885402929184579079D+00,
C     ~   -0.76105861802226257850D-01,
C     ~   -0.21707320837598329532D-01,
C     ~   -0.22806696336541141256D-04,
C     ~    0.13605855599644953723D-01,
C     ~    0.72655745785709626716D+00,
C     ~    0.33979077925435658968D+00,
C     ~    0.84015518189067677212D-01,
C     ~    0.38251019250246102387D-04,
C     ~   -0.23331745620323913747D-01,
C     ~   -0.12632616230609083896D+01,
C     ~   -0.66501059469106635901D+00,
C     ~   -0.11557383129450049530D+00,
C     ~   -0.46761660834837321942D-04,
C     ~    0.29188863483653651798D-01,
C     ~    0.16062192067320260058D+01,
C     ~    0.10008527972927729887D+01,
C     ~    0.96874076440719001368D-01,
C     ~    0.25527293306380216542D-04,
C     ~   -0.16076374426793411664D-01,
C     ~   -0.89077972230583615242D+00,
C     ~   -0.61419050551211071820D+00,
C     ~   -0.48351231840305022302D+00/

C     MODE 1 Peak at 4
	Data (abcd (1, j), j = 1, 25) /
     ~   -0.14004329253036342350D-07,
     ~    0.37688393599679170797D-05,
     ~    0.49241640240943258532D-04,
     ~    0.27365751661153279172D-03,
     ~   -0.29637572764194883845D-03,
     ~    0.23608493437412536586D-06,
     ~   -0.99590524882565070186D-04,
     ~   -0.13716499714981711388D-02,
     ~   -0.68651650806999411358D-02,
     ~    0.13086363094019834996D-01,
     ~   -0.17741658850262338020D-05,
     ~    0.95936069671529971003D-03,
     ~    0.14635178751353254966D-01,
     ~    0.13851244022798143618D-01,
     ~   -0.50732363273520704183D-01,
     ~    0.59728889093050563374D-05,
     ~   -0.35847147114331114892D-02,
     ~   -0.60417634422173911445D-01,
     ~   -0.24014471471244385192D-02,
     ~   -0.87454103956755879778D-02,
     ~   -0.44206961724389151058D-05,
     ~    0.27227339013033198256D-02,
     ~    0.47812648728699302935D-01,
     ~   -0.49208189730297080544D-02,
     ~    0.84257347502468000755D-02/


C	MODE 2
	Data (abcd (2, j), j = 1, 25) /
     ~    0.77668182761180606377D-08,
     ~   -0.20196291280988356575D-05,
     ~   -0.49273322968357806672D-04,
     ~    0.21539440444312103473D-03,
     ~   -0.83705446882959897436D-04,
     ~    0.33536606594698534777D-04,
     ~   -0.12866462057380196881D-01,
     ~   -0.31621223228494135248D+00,
     ~    0.61757501135689674143D+00,
     ~   -0.11476830385717948779D+00,
     ~   -0.17431745682254278229D-04,
     ~    0.65663659425452749474D-02,
     ~    0.15888921816088288352D+00,
     ~   -0.33062421999133104755D+00,
     ~    0.13245199895207815377D+00,
     ~   -0.26812881293819263994D-04,
     ~    0.10672183150391700845D-01,
     ~    0.27287667133660016283D+00,
     ~   -0.42942834594719494489D+00,
     ~   -0.32044936014889375819D+00,
     ~    0.10700204255233096706D-04,
     ~   -0.43707712411099670646D-02,
     ~   -0.11613927454426622443D+00,
     ~    0.71015669683247804044D-01,
     ~   -0.13462690359270563789D+01/



C	MODE 3
	Data (abcd (3, j), j = 1, 25) /
     ~    0.92353191216470111868D-02,
     ~   -0.47065512098254735917D+01,
     ~   -0.22554133829293983026D+03,
     ~   -0.23457137432873151894D+02,
     ~    0.46111664280590520803D+02,
     ~   -0.12239552063973810902D-01,
     ~    0.62483152611269074938D+01,
     ~    0.29997204397613383974D+03,
     ~    0.31604341484015141539D+02,
     ~   -0.57823703323961463951D+02,
     ~    0.41170674967848803760D-02,
     ~   -0.21219401709444536408D+01,
     ~   -0.10293144491195860279D+03,
     ~   -0.11947843043877295343D+02,
     ~    0.14488934852695538602D+02,
     ~   -0.18939743624590780868D-02,
     ~    0.99324763584517850034D+00,
     ~    0.49141758650519626883D+02,
     ~    0.77481048808290280405D+01,
     ~   -0.39887477569330695104D+01,
     ~    0.78113976972219205663D-03,
     ~   -0.41307234546100435323D+00,
     ~   -0.20642069925114081563D+02,
     ~   -0.40978583394053940125D+01,
     ~   -0.35575024161355539575D+01/



C	MODE 4
C	Data (abcd (4, j), j = 1, 25) /
C     ~    0.21992617916451395743D+02,
C     ~   -0.11399307544820429516D+05,
C     ~   -0.76667036088618019107D+06,
C     ~   -0.39185539354281475610D+06,
C     ~   -0.59239538019716295735D+03,
C     ~    0.47524457523937337910D+00,
C     ~   -0.24323272755420672908D+03,
C     ~   -0.16297746944208286734D+05,
C     ~   -0.78974933554931041612D+04,
C     ~   -0.34576773059694034806D+02,
C     ~   -0.67509882718360936237D+00,
C     ~    0.34607166925607515395D+03,
C     ~    0.23198896555695553978D+05,
C     ~    0.11312410520290614446D+05,
C     ~    0.45910310222559074233D+02,
C     ~   -0.21888994250814901576D+02,
C     ~    0.11346764342961157545D+05,
C     ~    0.76316098499373987706D+06,
C     ~    0.39024560870819433588D+06,
C     ~    0.58817978798588868016D+03,
C     ~    0.96230585996413289251D-01,
C     ~   -0.50295745279677515071D+02,
c     ~   -0.33917785686196171290D+04,
C     ~   -0.18056934409611232084D+04,
C     ~   -0.22829873613259561437D+02/


C     L=32
	Data (abcd (4, j), j = 1, 25)  /
     ~       0.38370706169389450224D+02,
     ~      -0.32479962500340833209D+05,
     ~      -0.94660043127081774372D+06,
     ~      -0.56864021527375161468D+06,
     ~      -0.27814993022411260703D+04,
     ~       0.20600370116089345984D+00,
     ~      -0.17025940328754780139D+03,
     ~      -0.49395843669987780089D+04,
     ~      -0.30000800342043403290D+04,
     ~      -0.97079695288155960497D+01,
     ~      -0.14579821620379142643D+01,
     ~       0.12204100420071475330D+04,
     ~       0.35488135357577208495D+05,
     ~       0.21388387530527799285D+05,
     ~       0.76407552757476420168D+02,
     ~      -0.37682332024326319341D+02,
     ~       0.31912372589853306159D+05,
     ~       0.93015110639191078689D+06,
     ~       0.55873046694301766734D+06,
     ~       0.28055544715537447331D+04,
     ~       0.56360430778218297831D+00,
     ~      -0.48256080155038478807D+03,
     ~      -0.14099258801109653127D+05,
     ~      -0.84807776710395010866D+04,
     ~      -0.13800109969955343203D+03/

C     L=32
C	Data (abcd (5, j), j = 1, 25)  /
C     ~       0.38370706169389450224D+02,
C     ~      -0.32479962500340833209D+05,
C     ~      -0.94660043127081774372D+06,
C     ~      -0.56864021527375161468D+06,
C     ~      -0.27814993022411260703D+04,
C     ~       0.20600370116089345984D+00,
C     ~      -0.17025940328754780139D+03,
C     ~      -0.49395843669987780089D+04,
C     ~      -0.30000800342043403290D+04,
C     ~      -0.97079695288155960497D+01,
C     ~      -0.14579821620379142643D+01,
C     ~       0.12204100420071475330D+04,
C     ~       0.35488135357577208495D+05,
C     ~       0.21388387530527799285D+05,
C     ~       0.76407552757476420168D+02,
C     ~      -0.37682332024326319341D+02,
C     ~       0.31912372589853306159D+05,
C     ~       0.93015110639191078689D+06,
C     ~       0.55873046694301766734D+06,
C     ~       0.28055544715537447331D+04,
C     ~       0.56360430778218297831D+00,
C     ~      -0.48256080155038478807D+03,
c     ~      -0.14099258801109653127D+05,
C     ~      -0.84807776710395010866D+04,
C     ~      -0.13800109969955343203D+03/


C     L=64
	Data (abcd (5, j), j = 1, 25)  /
     ~   -0.88135632848555722773D+01,
     ~    0.11514098078690224724D+05,
     ~    0.21560988124946391941D+06,
     ~    0.45832673059045454522D+05,
     ~    0.21572565347949170799D+03,
     ~    0.54760001008211061090D+00,
     ~   -0.70240562027308603987D+03,
     ~   -0.13144093980819422373D+05,
     ~   -0.29206319879254337479D+04,
     ~    0.16213317181975170910D+02,
     ~   -0.19639794571317899851D+03,
     ~    0.25945417987723824460D+06,
     ~    0.48620565446258918740D+07,
     ~    0.10126003837642234550D+07,
     ~    0.84424842477639039373D+04,
     ~    0.18558502462981778080D+03,
     ~   -0.24493795838444620827D+06,
     ~   -0.45896954969516352562D+07,
     ~   -0.95728945603260022778D+06,
     ~   -0.75585692609229706562D+04,
     ~    0.19078884327836625800D+02,
     ~   -0.25327914385527452445D+05,
     ~   -0.47482700183867816434D+06,
     ~   -0.98233311870422994616D+05,
     ~   -0.12493968682442764261D+04/


C     L=64
	Data (abcd (6, j), j = 1, 25)  /
     ~   -0.88135632848555722773D+01,
     ~    0.11514098078690224724D+05,
     ~    0.21560988124946391941D+06,
     ~    0.45832673059045454522D+05,
     ~    0.21572565347949170799D+03,
     ~    0.54760001008211061090D+00,
     ~   -0.70240562027308603987D+03,
     ~   -0.13144093980819422373D+05,
     ~   -0.29206319879254337479D+04,
     ~    0.16213317181975170910D+02,
     ~   -0.19639794571317899851D+03,
     ~    0.25945417987723824460D+06,
     ~    0.48620565446258918740D+07,
     ~    0.10126003837642234550D+07,
     ~    0.84424842477639039373D+04,
     ~    0.18558502462981778080D+03,
     ~   -0.24493795838444620827D+06,
     ~   -0.45896954969516352562D+07,
     ~   -0.95728945603260022778D+06,
     ~   -0.75585692609229706562D+04,
     ~    0.19078884327836625800D+02,
     ~   -0.25327914385527452445D+05,
     ~   -0.47482700183867816434D+06,
     ~   -0.98233311870422994616D+05,
     ~   -0.12493968682442764261D+04/


C	Peak at 6
C	Data (pq (1, j), j = 1, 10) /
C     ~       0.14671022415161132812D+02,
C     ~       0.17176361083984375000D+02,
C     ~       0.22174074172973634588D+02,
C     ~       0.38154026031494140625D+02,
C     ~       0.73457305908203123223D+02,
C     ~       0.11755837202072143554D+01,
C     ~       0.19926652908325195312D+01,
C     ~       0.38097207546234130859D+01,
C     ~       0.89784660339355468750D+01,
C     ~       0.32120651245117186611D+02/

C     MODE1 PEAK at 4
	Data (pq (1, j), j = 1, 10) /
     ~    0.67008705139160156250D+01,
     ~    0.10257417678833007368D+02,
     ~    0.17875274658203124111D+02,
     ~    0.36698616027832029473D+02,
     ~    0.73508399963378909802D+02,
     ~    0.14384461641311645507D+01,
     ~    0.27737910747528076171D+01,
     ~    0.63922395706176757812D+01,
     ~    0.35446662902832031250D+02,
     ~    0.73347290039062498223D+02/


C	MODE 2
	Data (pq (2, j), j = 1, 10) /
     ~      0.62806210517883300781D+01,
     ~      0.19055501937866210937D+02,
     ~      0.16761466979980468750D+02,
     ~      0.28816335678100584161D+02,
     ~      0.77478836059570310723D+02,
     ~      0.12525197267532348632D+01,
     ~      0.21577339172363281250D+01,
     ~      0.42458128929138183593D+01,
     ~      0.10150924682617188388D+02,
     ~      0.32478603363037108486D+02/
C	MODE 3
	Data (pq (3, j), j = 1, 10) /
     ~     -0.20120578765869141513D+02,
     ~      0.20586240768432615411D+02,
     ~      0.23876194000244139736D+02,
     ~      0.38100570678710937500D+02,
     ~      0.71865959167480468750D+02,
     ~      0.11369531154632568359D+01,
     ~      0.18996990919113159179D+01,
     ~      0.35767178535461425781D+01,
     ~      0.82503862380981445312D+01,
     ~      0.30094644546508790838D+02/
C	MODE 4
C	Data (pq (4, j), j = 1, 10) /
C     ~      0.38784885406494140625D+02,
C     ~      0.26744506835937498223D+02,
C     ~      0.27683023452758788174D+02,
C     ~      0.38958595275878904473D+02,
C     ~      0.72836189270019531250D+02,
C     ~      0.11952087879180908203D+01,
C     ~      0.20218055248260498046D+01,
C     ~      0.38115730285644531250D+01,
C     ~      0.88158550262451171875D+01,
C     ~      0.33938060760498047763D+02/

C     L=32
	Data (pq (4, j), j = 1, 10) /
     ~      0.38836956024169921875D+02,
     ~      0.24311889648437499111D+02,
     ~      0.29589153289794921875D+02,
     ~      0.39441062927246095526D+02,
     ~      0.74735847473144527697D+02,
     ~      0.13389585018157958984D+01,
     ~      0.23527803421020507812D+01,
     ~      0.45696606636047363281D+01,
     ~      0.10911170959472655894D+02,
     ~      0.44013973236083980822D+02/


C     L=32
C	Data (pq (5, j), j = 1, 10) /
C     ~      0.38836956024169921875D+02,
C     ~      0.24311889648437499111D+02,
C     ~      0.29589153289794921875D+02,
C     ~      0.39441062927246095526D+02,
C     ~      0.74735847473144527697D+02,
C     ~      0.13389585018157958984D+01,
C     ~      0.23527803421020507812D+01,
C     ~      0.45696606636047363281D+01,
C     ~      0.10911170959472655894D+02,
C     ~      0.44013973236083980822D+02/
C
C     L=64
	Data (pq (5, j), j = 1, 10) /
     ~     0.36148849487304688388D+02,
     ~     0.26422071456909179687D+02,
     ~     0.52445343017578123223D+02,
     ~     0.50159130096435546875D+02,
     ~     0.72792282104492187500D+02,
     ~     0.14043723344802856445D+01,
     ~     0.26591603755950927734D+01,
     ~     0.55725831985473632812D+01,
     ~     0.15859733581542969638D+02,
     ~     0.70801498413085939276D+02/



C     L=64
	Data (pq (6, j), j = 1, 10) /
     ~     0.36148849487304688388D+02,
     ~     0.26422071456909179687D+02,
     ~     0.52445343017578123223D+02,
     ~     0.50159130096435546875D+02,
     ~     0.72792282104492187500D+02,
     ~     0.14043723344802856445D+01,
     ~     0.26591603755950927734D+01,
     ~     0.55725831985473632812D+01,
     ~     0.15859733581542969638D+02,
     ~     0.70801498413085939276D+02/



c	Separate out the a. b, c, d and p, q.
	KA=1
	Do 10 I=1,M
	Do 10 J=1,M
	a(i,j)=ABCD(L,KA)
	KA=KA+1
   10 continue
c	Now the PQs
	Do 14 I=1,M
	J=I+M
	p(i)=pq(L,I)/r0or
	p(j)=pq(L,J)/r0or
   14 continue

c	Now compute the magnetic field from differentiating the potential U.

	call gradU(Bxm,Bym,Bzm,x,y,z,a,p,M,M1,M2)

c
    3 continue
      Return
	End
C*************************************************************************
	Subroutine gradU(Bxm1,Bym1,Bzm1,x1,y1,z1,a,p,M,M1,M2)
	IMPLICIT REAL*8(A-H,O-Z)
c
c	Calculates the gradient of the potential.
c
c
	Dimension a(M,M),p(M2)
	M1=M1
	Bxm1=0.
	Bym1=0.
	Bzm1=0.
	Do 111 i=1,M
	Do 222 k=M+1,M2
	J=K-M
c	Write(*,*)M,M1,M2,p(i),p(k)
	Term1=dsqrt(1.D0/p(i)**2+1.D0/p(k)**2)
	Term2=dexp(Term1*x1)
	Bxm1=Bxm1+Term1*a(i,j)*term2*dcos(y1/p(i))*dsin(z1/p(k))
	Bym1=Bym1- a(i,j)/p(i)*term2*dsin(y1/p(i))*dsin(z1/p(k))
	Bzm1=Bzm1+ a(i,j)/p(k)*term2*dcos(y1/p(i))*dcos(z1/p(k))
  222 continue
  111 continue
	Bxm1=-Bxm1
	Bym1=-Bym1
	Bzm1=-Bzm1
	Return
	End
C*********************************************************************
C******************************************************************
	SUBROUTINE CAR2SPH_MAG(BX,BY,BZ,BR,BTH,BPHI,TH,PHI)
	IMPLICIT REAL*8(A-H,O-Z)
	! ARFKEN
	BR  = BX*dsin(TH)*dcos(PHI)+BY*dsin(TH)*dsin(PHI)+BZ*dcos(TH)
	BTH = BX*dcos(TH)*dcos(PHI)+BY*dcos(TH)*dsin(PHI)-BZ*dsin(TH)
	BPHI=-BX*dsin(PHI)+BY*dcos(PHI)
	RETURN
	END
C******************************************************************
	SUBROUTINE getIMF_penetration(BY_IMF,BZ_IMF,BY_p,BZ_p)
	IMPLICIT REAL*8(A-H,O-Z)
	DATA e1/0.06816/
	DATA e2/0.55417/

	pi=3.14159265358979


	IF ((BY_IMF .eq. 0.0) .AND. (BZ_IMF .eq. 0.0)) THEN
		THETA=0.0
	ELSE
		THETA=DATAN2(BY_IMF,BZ_IMF)
		IF (THETA .lt. 0.0 ) THEN
			THETA = THETA + 2*pi
		ENDIF
	ENDIF

	BY_p=(e1+e2*(dcos(THETA/2.0))**2)*BY_IMF
	BZ_p=(e1+e2*(dcos(THETA/2.0))**2)*BZ_IMF

	RETURN
	END

C******************************************************************************
	SUBROUTINE checkIfInsideMpSATURN(RR,x,NOUT,Dp,a1,a2,dK)
	IMPLICIT REAL*8(A-H,O-Z)



	r0=a1*Dp**(-a2)

	cosTheta=x/RR
	Rmp=r0*(2/(1+cosTheta))**dK

	IF (Rmp .lt. RR) THEN
		NOUT = 0
	ELSE
		NOUT = 1
	ENDIF


	RETURN
	END
c*********************************************************************

C*********************************************************************************
C
      SUBROUTINE KRONIAN(NM,R,T,F,BR,BT,BF)
C
C     BASED ON THE SUBROUTINE IGRF WRITTEN BY N.A. TSYGANENKO (1979)
C     MODIFIED BY KRISHAN KHURANA, JULY, 1996. AND   NOV. 2004.
C
C     CALCULATES COMPONENTS OF MAIN KRONIAN FIELD IN RIGHT HANDED SPHERICAL
C     COORD SYSTEM. BASED ON THE  SPHERICAL HARMONIC COEFFICIENTS GIVEN BY
C     ACUNA ET AL. [1983] (Z3 MODEL,  V1+V2 DATA)
C     MAXIMUM ORDER OF HARMONICS TAKEN INTO ACCOUNT (NOT MORE THAN 0RDER 3)
C
C
C     IT IS ASSUMED THAT THE TRAJECTORY IS IS IN RIGHT HANDED S III COORDINATES.
C     THE OUTPUT IS ALSO IN RTP (RH) COORDINATES.
C
C            INPUT:  NM (INTEGER)- MAXIMUM ORDER OF HARMONICS TAKEN
C                                  INTO ACCOUNT (NM.LE.12)
C
C                    R,T,F (REAL)- POSITION OF DESIRED FIELD VALUE IN
C                                  SPHERICAL JOVIGRAPHIC COORDINATE SYSTEM
C                                  (R IN PLANET RADII, COLATITUDE T AND
C                                   RH LONGITUDE F IN RADIANS)
C
C           OUTPUT: BR,BT,BF (REAL)- COMPONENTS OF THE INTERNAL PORTION
C                                    OF THE MAIN MAGNETIC FIELD IN
C                                    SPHERICAL S III COORD SYSTEM
C                                    (VALUES GIVEN IN GAMMA)
C
C***************************************************************************
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL BK,BM
      DIMENSION A(13),B(13),G(91),H(91),REC(91)



	DATA G/0.,21144.,24.,1449.,-179.,-35.,2012.,
     ~       53.,92.,-40.,
     ~       81*0.0/

	DATA H/0.,0.,-0.9110,0.,43.5680,-44.7762,0.,
     ~       140.9273,-92.8876,68.2759,
     ~       81*0.0/



      DATA FIRSTI/0.0/
C     WRITE(1,'(5F15.3)')(G(I),I=1,10)
C     WRITE(1,'(5F14.3)')(H(I),I=1,10)
      IF(FIRSTI.EQ.0.0) GO TO 1
      GO TO 6
 1    FIRSTI=1.0
      G(1)=0.
      H(1)=0.
      KNM=15
      DO 2 N=1,13
      N2=2*N-1
      N2=N2*(N2-2)
      DO 2 M=1,N
      MN=N*(N-1)/2+M
  2   REC(MN)=DBLE((N-M)*(N+M-2))/DBLE(N2)
      S=1.
      DO 5 N=2,13
      MN=N*(N-1)/2+1
      S=S*DBLE(2*N-3)/DBLE(N-1)
      G(MN)=G(MN)*S
      H(MN)=H(MN)*S
      P=S
      DO 5 M=2,N
      AA=1.
      IF (M.EQ.2) AA=2.
      P=P*DSQRT(AA*DBLE(N-M+1)/DBLE(N+M-2))
      MNN=MN+M-1
      G(MNN)=G(MNN)*P
  5   H(MNN)=H(MNN)*P
  6   IF(KNM.EQ.NM) GOTO 61
      KNM=NM
      K=KNM+1
 61   PP=1./R
      P=PP
      DO 7 N=1,K
      P=P*PP
      A(N)=P
  7   B(N)=P*N
      P=1.
      D=0.
      BBR=0.
      BBT=0.
      BBF=0.
      U=T
      CF=DCOS(F)
      SF=DSIN(F)
      C=DCOS(U)
      S=DSIN(U)
      BK=(S.LT.1.D-5)
      DO 12 M=1,K
      BM=(M.EQ.1)
      IF(BM) GOTO 8
      MM=M-1
      W=X
      X=W*CF+Y*SF
      Y=Y*CF-W*SF
      GOTO 9
  8   X=0.
      Y=1.
  9   Q=P
      Z=D
      BI=0.
      P2=0.
      D2=0.
      DO 11 N=M,K
      AN=A(N)
      MN=N*(N-1)/2+M
      E=G(MN)
      HH=H(MN)
      W=E*Y+HH*X
      IF (DABS(P2).LT.1.D-38) P2=0.0
      IF (DABS(Q).LT.1.D-38) Q=0.0
      BBR=BBR+B(N)*W*Q
      BBT=BBT-AN*W*Z
      IF(BM) GOTO 10
      QQ=Q
      IF(BK) QQ=Z
      BI=BI+AN*(E*X-HH*Y)*QQ
  10  XK=REC(MN)
      DP=C*Z-S*Q-XK*D2
      PM=C*Q-XK*P2
      D2=Z
      P2=Q
      Z=DP
      Q=PM
   11 CONTINUE
      D=S*D+C*P
      P=S*P
      IF(BM) GOTO 12
      BI=BI*MM
      BBF=BBF+BI
  12  CONTINUE
      BR=BBR
      BT=BBT
      IF(BK) GOTO 13
      BF=BBF/S
      GOTO 14
  13  IF(C.LT.0.) BBF=-BBF
      BF=BBF
  14  CONTINUE
      RETURN
      END
C
C

C*************************************************************************


