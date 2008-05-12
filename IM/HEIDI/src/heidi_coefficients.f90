! File name: heidi_coefficients_010.f90
!
! Contains: drift and diffusion coefficient definition routines for HEIDI
!	CEPARA
!	OTHERPARA
!	MAGCONV
!
! Last Modified: March 2006, Mike Liemohn
!
!***********************************************************************
!			       CEPARA
! Calculates the cross-section of charge exchange of ring current
! ion species with the neutral hydrogen and the charge exchange rate
! Also calculates the atmospheric loss rate (out the loss cone)
!***********************************************************************
	SUBROUTINE CEPARA

	use ModHeidiSize
	use ModHeidiMain
	use ModHeidiDrifts

	implicit none

	integer I,K,L,I1,LUP,ier
	real x,y		!,taub
	real funt,FACI,fac,cosd,acosd
	REAL RLAMBDA(71),PA(71),HDNS(NR,NPA),LH(20),HDNSIN(20,71)
	CHARACTER*80 TITLE
	external cosd,acosd

!.......Open a file with the bounce-averaged H dens [m-3] in geocorona
!	obtained through bounce-averaging the Chamberlain model
!	(Mei-Ching's fit by a poly of order 4 good for L from 1 to 6.5)
	S=0
	DO I=2,NS
	  S=S+SCALC(I)   ! S=0, no ions in calc, S>0, ions in calc
	END DO
	IF (S.GE.1) THEN	! Only needed for ion charge exchange
! Read in from file instead
!	  DO I=1,20
!	    LH(I)=1.5+.25*REAL(I)
!	  END DO
!cc	  IF(LO.EQ.31) THEN
!cc	    OPEN(UNIT=8,FILE='hgeo.in',STATUS='OLD')
!cc	    LUP=LO-1
!cc	  ELSE
!cc	    IF(LO.EQ.50) OPEN(UNIT=8,FILE='hgeo50.in',STATUS='OLD')
!cc	    IF(LO.EQ.60) OPEN(UNIT=8,file='hgeo64.in',STATUS='OLD')
!cc	    IF(LO.EQ.71) OPEN(UNIT=8,file='hgeo71.in',STATUS='OLD')
	    OPEN(UNIT=8,file='hgeo71.in',STATUS='OLD')
!cc	    IF(LO.EQ.91) OPEN(UNIT=8,file='hgeo91.in',STATUS='OLD')
	    LUP=71
!cc	  ENDIF
	  do I=1,20
	    READ(8,3) TITLE,LH(i)
	    do L=1,LUP
	      READ(8,4) PA(L),RLAMBDA(L),HDNSIN(I,L)
	    end do	! L loop
!cc	    if (LO.EQ.31) HDNSIN(I,LUP+1)=HDNSIN(I,LUP)
	  end do	! I loop
	  CLOSE(8)
!cc	  If (LO.EQ.31) LUP=LUP+1
	  do L=1,LUP
	    PA(L)=COSD(PA(L))
	  end do
          DO I=1,IO
	   do L=1,LO
	     CALL LINTP2(LH,PA,HDNSIN,20,71,LZ(i),MU(L),FAC,IER)
	     HDNS(i,l)=FAC
           enddo
          END DO
	    
	END IF
3	FORMAT(A61,F5.2)
4	FORMAT(2F7.3,2X,1PE12.5)

!.......Calculate charge exchange cross-section of species S with H
!         and then the charge exchange decay rate ACHAR

!.......H+ charge exchange
	IF (SCALC(2).EQ.1) THEN
	 do K=2,KO
	  X=ALOG10(EKEV(K))
	  if (X.LT.-2.) X=-2.
	  Y=-18.767-0.11017*X-3.8173e-2*X**2-0.1232*X**3-5.0488e-2*X**4
!............10**Y is cross section of H+ in m2
	  do L=2,Lo
	    achar(2:io,k,l,2)=exp(-(10.**Y*V(K,2)*HDNS(2:io,L)*DT))
	  end do	! L loop	
	 end do	! K loop
	END IF

!.......He+ charge exchange
	IF (SCALC(3).EQ.1) THEN
	 do K=2,KO
	  X=ALOG10(EKEV(K))
	  if (X.LT.-2.) X=-2.
	  Y=-20.789+0.92316*X-0.68017*X**2+0.66153*X**3-0.20998*X**4
!..........10**Y is cross sect of He+ in m2
	  do L=2,Lo
	    ACHAR(2:io,K,L,3)=EXP(-(10.**Y*V(K,3)*HDNS(2:io,L)*DT))
	  end do	! L loop
	 end do	! K loop
	END IF

!.......O+ charge exchange
	IF (SCALC(4).EQ.1) THEN
	 do K=2,KO
	  X=ALOG10(EKEV(K))
	  if (X.LT.-2.) X=-2.
	  Y=-18.987-0.10613*X-5.4841E-3*X**2-1.6262E-2*X**3   &
            -7.0554E-3*X**4
!...........10**Y is cross sect of O+ in m2
	  do L=2,Lo
	    achar(2:io,k,l,4)=exp(-(10.**Y*V(K,4)*HDNS(2:io,L)*DT))
	  end do	! L loop
	 end do	! K loop
	END IF

!.......Calculate the losses due to the collis with atmosphere
	do s=1,NS
	 do I=2,IO
	  do K=2,KO
	   do L=UPA(I),LO
	     ATLOS(I,K,L,S)=EXP(-DT/(29.1347*LZ(I)*FUNT(MU(L))*SQRT(M1(S)/EKEV(K))))
	   end do	! L loop
	  end do	! K loop
	 end do		! I loop
	end do		! S loop

	RETURN
	END
!
! End of subroutine CEPARA
!

!***********************************************************************
!				OTHERPARA
!	  Calculate parameters used in drifts and Coulomb drags
!***********************************************************************
	SUBROUTINE OTHERPARA

	use ModHeidiSize
	use ModHeidiIO
	use ModHeidiMain
	use ModHeidiDrifts

	implicit none

	integer i,j,k,l,is,iss,ier
	real eps,dln,sqm,qe,gama,cco,ccd,edrco,x,xd,g,cci,cdi,edri,ccde,ccdi
        real bane,badif,c,gpa,cce,cde,edre
	real funi,funt,erf
	REAL COULDE(NE,NPA),COULDI(NE,NPA),AFIR,ASEC
	REAL VF(NSTH),RA(NSTH),MUBOUN,MULC,TMAS(NSTH),TM1(NSTH)
	external funi,funt,erf

	DATA TM1/5.4462E-4,1.0,4.0,16.0/  ! Thermal mass in AMU
	DATA RA/1.,.77,.2,.03/     ! Thermal proportions in plasmasphere

!.......Parameters used in calculating drifts at boundaries
	C=1.44E-2*RE**2		! Constant of corotation, [C]=V*m
	ISS=-1			! sign of specie's charge
	IF (S.GE.2) ISS=1

!.......Determine the dayside
	J6=0
	J=0
	do WHILE (J6.EQ.0)
	  J=J+1
	  IF (MLT(J).GE.5.99) J6=J
	end do	! While J6 loop
	J18=0
	J=0
	do WHILE (J18.EQ.0)
	  J=J+1
	  IF (MLT(J).GE.17.99) J18=J
	end do	! While J18 loop

	do K=1,KO
	 do i=1,io
	   do L=1,UPA(I)    ! Kp independent part of azimuthal drift
	    GPA=1.-FUNI(MU(L))/6./FUNT(MU(L))
	    P2(I,K,L)=(C-ISS*3.*EKEV(K)*1000.*Z(I)*GPA)*DT/DPHI/ME
!	    P2(I,K,L)=C*DT/DPHI/ME
	   end do	! 1st L loop
	   do L=UPA(I)+1,LO
	    P2(I,K,L)=(C-ISS*3.*EKEV(K)*1000.*Z(I)*GPA)*DT/DPHI/ME
!	    P2(I,K,L)=C*DT/DPHI/ME
	   end do	! 2nd L loop
	 end do
	end do	! K loop
	do I=1,IO
	  VR(I,J6)=0.
	  VR(I,J18)=0.
	end do	! I loop processing VR

!.......Pitch angle and energy time derivatives at boundaries of grids
	do I=1,IO
	 do J=1,JO
	  do L=1,UPA(I)
	   MUBOUN=MU(L)+0.5*WMU(L)           ! MU at boundary of grid
!	   write(*,*) "muboun = ", muboun
	   MUDOT(I,J,L)=(1.-MUBOUN**2)*FUNI(MUBOUN)/LZ(I)  &
     	     /MUBOUN/4./FUNT(MUBOUN)*DL1/DMU(L)
	   GPA=1.-FUNI(MU(L))/6./FUNT(MU(L))
	   do K=1,KO
	    EDOT(I,J,K,L)=-3.*EBND(K)/LZ(I)*GPA*DL1/DE(K)
	   end do	! K loop
	  end do 	! L loop
	  MULC=MU(UPA(I))+0.5*WMU(UPA(I))
	  do L=UPA(I)+1,LO-1
	   MUBOUN=MU(L)+0.5*WMU(L)
	   MUDOT(I,J,L)=(1.-MUBOUN**2)*FUNI(MULC)/LZ(I)  &
     	     /MUBOUN/4./FUNT(MULC)*DL1/DMU(L)
	   do K=1,KO
	    EDOT(I,J,K,L)=-3.*EBND(K)/LZ(I)*GPA*DL1/DE(K)
	   end do	! K loop
	  end do	! L loop
	  MUDOT(I,J,LO)=0.
	  do K=1,KO
            EDOT(I,J,K,LO)=-3.*EBND(K)/LZ(I)*GPA*DL1/DE(K)
	  end do	! K loop
	 end do 	! J loop
	end do	! I loop

!.......We assume Te=Ti=1eV (kT=1eV)
	do I=1,NSTH
	  TMAS(I)=TM1(I)*MP
	  VF(I)=SQRT(2.*Q/TMAS(I))	! [m/s]
	end do	! I loop
	EPS=8.854E-12                ! [F/m]
	DLN=21.5                     ! Coulomb logarithm
	QE=(Q**2/EPS)
	DO S=1,NS
	  IF (SCALC(S).EQ.1) THEN
	    SQM=SQRT(MAS(S))
!!! GAMA and EDRCO divided by an extra 2 to match Khazanov defn of Q
	    GAMA=DLN/8./PI*QE*1E6/SQM*QE ! *1E6 to transform ne into m-3
	    CCO=GAMA/Q*SQRT(2.)*DT
	    CCD=GAMA*2*DT/16./SQRT(2.)   ! 2*DT due to double the time step
	    EDRCO=2.*DLN*QE*RE/MAS(S)*QE/MAS(S)*(1.E-10) ! to invert in eV/cm2/s
! There was an autoscope directive here, but I removed it
	    do K=1,KO
	      X=VBND(K,S)/VF(1)
	      XD=V(K,S)/VF(1)
	      CCE=RA(1)*G(X)		     ! electrons
	      CDE=RA(1)*(ERF(XD,IER)-G(XD))
	      EDRE=RA(1)*G(XD)
	      CCI=0.
	      CDI=0.
	      EDRI=0.
	      do IS=2,NSTH
	       X=VBND(K,S)/VF(IS)
	       XD=V(K,S)/VF(IS)
	       CCI=CCI+RA(IS)*G(X)	     ! ions
	       EDRI=EDRI+RA(IS)*G(XD)
	       CDI=CDI+RA(IS)*(ERF(XD,IER)-G(XD))
	      end do	! IS loop
!.......Collisions with plasmaspheric e- (minus: dF/dt-d(cF)/dt=0)
	      COULE(1,K,1,S)=-CCE*SQRT(EBND(K)/Q/1000.)/DE(K)*CCO
!.......Collisions with plasmaspheric ions
	      COULI(1,K,1,S)=-CCI*SQRT(EBND(K)/Q/1000.)/DE(K)*CCO
!.......Energy deposition rate [eV/cm2/s]
	      CEDR(K,1,S)=EDRCO*EKEV(K)*EDRE
	      CIDR(K,1,S)=EDRCO*EKEV(K)*EDRI
!.......Coulomb scattering
	      CCDE=CCD*CDE/(EKEV(K)*Q*1000)**(1.5)
	      CCDI=CCD*CDI/(EKEV(K)*Q*1000)**(1.5)
	      COULDE(K,1)=0.
	      GTAE(K,2,S)=0.
	      COULDI(K,1)=0.
	      GTAI(K,2,S)=0.
	      do L=2,LO-1
	       BANE=(1.-FUNI(MU(L))/2./FUNT(MU(L)))/(1.-MU(L)*MU(L))
	       COULE(1,K,L,S)=COULE(1,K,1,S)*BANE
	       COULI(1,K,L,S)=COULI(1,K,1,S)*BANE
	       CEDR(K,L,S)=CEDR(K,1,S)*BANE*FACMU(L)*WMU(L)
	       CIDR(K,L,S)=CIDR(K,1,S)*BANE*FACMU(L)*WMU(L)
	       MUBOUN=MU(L)+0.5*WMU(L)
	       BADIF=(1.-MUBOUN*MUBOUN)*FUNI(MUBOUN)/MUBOUN
!.......Collisions with plasmaspheric electrons
	       COULDE(K,L)=CCDE*BADIF
	       ATAE(K,L,S)=COULDE(K,L)/FACMU(L+1)/DMU(L)/WMU(L)
	       AFIR=COULDE(K,L)/FACMU(L)/DMU(L)/WMU(L)
	       ASEC=COULDE(K,L-1)/FACMU(L)/DMU(L-1)/WMU(L)
	       BTAE(K,L,S)=AFIR+ASEC
	       GTAE(K,L+1,S)=COULDE(K,L)/FACMU(L)/WMU(L+1)/DMU(L)
!.......Collisions with plasmaspheric ions
	       COULDI(K,L)=CCDI*BADIF
	       ATAI(K,L,S)=COULDI(K,L)/FACMU(L+1)/DMU(L)/WMU(L)
	       AFIR=COULDI(K,L)/FACMU(L)/DMU(L)/WMU(L)
	       ASEC=COULDI(K,L-1)/FACMU(L)/DMU(L-1)/WMU(L)
	       BTAI(K,L,S)=AFIR+ASEC
	       GTAI(K,L+1,S)=COULDI(K,L)/FACMU(L)/WMU(L+1)/DMU(L)
	      end do	! L loop
	      COULDE(K,LO)=0.
	      ATAE(K,LO,S)=0.
	      COULDI(K,LO)=0.
	      ATAI(K,LO,S)=0.
	    end do	! K loop
	    do K=1,KO
 	     do I=IO,1,-1		! Set to edge of loss cone when
	      do L=1,UPA(I)-1	! in the loss cone
	        COULE(I,K,L,S)=COULE(1,K,L,S)
	        COULI(I,K,L,S)=COULI(1,K,L,S)
	      end do	! L loop	
	      do L=UPA(I),LO
	        COULE(I,K,L,S)=COULE(I,K,L-1,S)
	        COULI(I,K,L,S)=COULI(I,K,L-1,S)
	      end do	! 2nd L loop
	     end do	! I loop
	    end do	! K loop
	  END IF	! SCALC check
	END DO		! S loop

!...Set continuous source/loss integrals to zero
	RNS=0.	! Radial drift particle source
	RNL=0.	! Radial drift particle loss
	RES=0.	! Radial drift energy source
	REL=0.	! Radial drift energy loss
	ESN=0.	! Energy drift particle gain
	ELN=0.	! Energy drift particle loss
	ESE=0.	! Energy drift energy gain
	ELE=0.	! Energy drift energy loss
	ECN=0.	! Coulomb collision particle gain/loss
	ECE=0.	! Coulomb collision energy gain/loss
	ALN=0.	! Atmospheric loss particle gain/loss
	ALE=0.	! Atmospheric loss energy gain/loss
	CEN=0.	! Charge exchange particle gain/loss
	CEE=0.	! Charge exchange energy gain/loss
!...Set "magnetopause" at outer boundary
	ILMP(1:JO)=IO
	LMP(1:JO)=LZ(IO)+DL1


	RETURN
	END
!
! End of subroutine OTHERPARA
!

! **********************************************************************
!                             MAGCONV
!	Calculates magnetospheric convection strength and puts it
!	into the drift terms 
!	(added to already-calculated corotation and gradient-curvature)
!  Calculation options:
!	IA=0  No convection field
!	IA=1  Kp driven V-S with Maynard and Chen activity dependence
!	IA=2  MBI driven V-S, acitivity from force balance at dusk
!	IA=3  McIlwain field
!	IA=4  McIlwain field, Eo driven by Cross Polar Cap Potential
!	IA=5  Kp V-S field plus Burke-Wygant penetration field
!	IA=6  McIlwain field plus Burke-Wygant penetration field
!	IA=7  McIlwain(CPCP) field plus B-W penetration field
!	IA=8  Kp V-S field plus self-consistent penetration field
!	IA=9  McIlwain field plus self-consistent penetration field
!	IA=10 McIlwain(CPCP) field plus S-C penetration field
!	IA=11 Unshielded V-S(CPCP) field plus penetration field
!	IA=12 Unshielded VS(CPCP) field (no Epen)
!       IA=13 W96 field plus S-C penetration field
!       IA=14 W96 field everywhere (done in S-C field block)
!       IA=15 AMIE field plus S-C penetration field
!       IA=16 AMIE field everywhere (done in S-C field block)
!       IA=20 Read in electric field from a file
!       IA=21 AMIE potentials everywhere
!       IA=22 Weimer-96 potentials everywhere
!       IA=23 Foster potentials everywhere
!       IA=?? More field models in Aaron's new subroutine? Put them in...
!  Outputs (both are unitless advection coefficients):
!	VR => Radial drift
!	P1 => Activity dependent azimuthal drift
! **********************************************************************
	SUBROUTINE MAGCONV(I3,NST)

	use ModHeidiSize
	use ModHeidiMain
	use ModHeidiIO
	use ModHeidiDrifts
	use ModHeidiCurrents
	use ModHeidiDGCPM

	implicit none

	integer ABASE(24),AEPEN(24),IER,j,i,ii,I4,I3,NST,jj,bc_choice
	integer L,IOpot,JOpot,iproc
	real KR,SPJ,CPJ,PJ,DLMAG,Ppc0,Eom,EOJ,EOJ1,NY(4),FAC,CRo
	real KGAM,KS(4),KALP,KBETA,KPHI,KEX,KEY,LP,DLP,PCO,PHIPOFF
	real Ro,RR,RoRg,KEPS,KSD,SJ,CJ,PHIP,EPR,EPP,ERF,DP1,FACP
        real DPP,DP2,sLP,dLPdphi,Jreducer
	real sind,cosd
	real jfac_temp(NR+3,NT)
	real FPOT1(3*NR,NT),FPOT2(3*NR,NT),FPOT12(3*NR,NT)
	real TP1,TP2,LZpot(3*NR),Data(3*NR),Fmin,Fmax,MLTpot(NT)
	real t_sawtooth(9),tdiff,dt_saw
	character*80 header
	DATA ABASE/1,1,1,2,3,1,2,3,1,2,3,4,4,5,5,6,6,1,1,1,10,11,12,13/
	DATA AEPEN/0,0,0,0,0,1,1,1,2,2,2,2,0,2,3,2,3,0,0,0, 0, 0, 0, 0/
	DATA KGAM/8./, KS/9.8,-1.4,-0.9,-0.3/
	DATA KALP/0.3/, KBETA/0.8/, KPHI/3000./
	DATA KEX/3.E-5/, KEY/1.2E-4/, PCO/9.17E4/
	DATA t_sawtooth/9.6e4,9.84e4,1.062e5,1.152e5,1.278e5,   &
                       1.374e5,1.458e5,1.554e5,1.626e5/
!cc Use the next t_sawtooth definition for 4.5 h of 2*Vsw early on 4/18
!cc	DATA t_sawtooth/8.64e4,8.82e4,9.0e4,9.18e4,9.36e4,
!cc     &                  9.54e4,9.72e4,9.9e4,1.008e5/
	DATA dt_saw/600./
	EXTERNAL ERF,sind,cosd

! Variables to send to Aaron's subroutine
        real eyear, eday, ehour, eminute, esecond, eby, ebz, evsw
        integer nemlts1, nelats1, nemlts2, nelats2, nechoice
        real emlts1(NT,NR+3), elats1(NT,NR+3), epots1(NT,NR+3)
        real emlts2(nphicells,nthetacells), elats2(nphicells,nthetacells),  &
         epots2(nphicells,nthetacells)

	SAVE I4,TP1,TP2,FPOT1,FPOT2,IOpot,LZpot,JOpot,MLTpot

	integer edayplus
	real univ_time

	PHIPOFF=0.
	DP1=.4*PI
	DP2=2.*PI/3.
	FACP=AMIN1((.25*KP)**2,1.5)
	Kr=KP/(1.+.1*KP)
!**  The next line is to cut down the size of the field-aligned currents...
	Jreducer=1.0

!** Find polar cap potential value (if needed)
        print *, 'Is ABASE 3 or 4?'
	IF (ABASE(IA+1).EQ.3 .or. ABASE(IA+1).EQ.4) THEN  
	  IF (I3.EQ.NST) THEN
	    I4=0
	    ii=2
	    DO WHILE (I4.eq.0)
		IF (TPPC(ii).GE.T) THEN
		  I4=ii-1
		ELSE IF (ii.EQ.IPPC-1) THEN
		  I4=IPPC-1
		ELSE
		  ii=ii+1
		END IF
	    END DO
	  ELSE IF (T+2.*DT.GE.TPPC(IPPC)) THEN
	    I4=IPPC-1
	  ELSE IF (T.GE.TPPC(I4+1)) THEN
	    I4=I4+1
	  END IF
!	  print *, 'PPC:',I4,IPPC,TPPC(I4),PPC(I4)
	  IF (T.LE.TPPC(1)) THEN
	    PPC0=PPC(1)
	  ELSE IF (T.GE.TPPC(IPPC)) THEN
	    Ppc0=PPC(IPPC)
	  ELSE
	    Ppc0=PPC(I4)+(PPC(I4+1)-PPC(I4))*(T-TPPC(I4))   &
             /AMAX1(.1*DT,TPPC(I4+1)-TPPC(I4))
	  END IF
!	print *, 'Ppc0:',I4,TPPC(I4),TPPC(I4+1),PPC(I4),PPC(I4+1)
	END IF

	DLMAG=LMP(7)+LMP(19)  ! Width across magnetosphere at X=0
	Eom=sqrt(KEX**2+KEY**2)  ! Eo factor for McIlwain field

	IF (IA.EQ.0) A=0.     ! No convection field

!  Calculate base convection electric field

!CC Shielded Volland-Stern
        print *, 'Is ABASE 1-4?'

	IF (ABASE(IA+1).EQ.1) THEN

	 do J=1,JO   ! Fill in drift values
	  do i=1,io
	   VR(I,J)=-A*COS(PHI(J))*(LZ(I)+0.5*DL1)**(LAMGAM+2.)*   &
     		DT/DL1*(RE*RE/ME)
	   P1(I,J)=A*LAMGAM*LZ(I)**(LAMGAM+1.)*SIN(PHI(J)+0.5*DPHI)*   &
     		DT/DPHI*(RE*RE/ME)
	   BASEPOT(i+1,j)=A*RE*LZ(i)**(LAMGAM)*SIN(phi(j))
	  end do
	   BASEPOT(1,j)=A*RE*(LZ(i)-DL1)**(LAMGAM)*SIN(phi(j))
	   BASEPOT(i+1,j)=A*RE*(LZ(io)+dl1)**(LAMGAM)*SIN(phi(j))
	   BASEPOT(i+1,j)=A*RE*(LZ(io)+2*dl1)**(LAMGAM)*SIN(phi(j))	  
	 end do
	 do J=1,nphicells   ! Fill in DGCPM potentials
	  do I=1,nthetacells
	   potdgcpm(i,j)=A*RE*vlzcells(I)**(LAMGAM)*SIND(vphicells(J))
	  enddo
	 enddo

!CC McIlwain or modified McIlwain field
	ELSE IF (ABASE(IA+1).EQ.2 .OR. ABASE(IA+1).EQ.3) THEN

	 KEPS=1.+KALP*Kr   ! McIlwain E5D potential structure
	 IF (ABASE(IA+1).EQ.3) KEPS=Ppc0/(DLMAG*Re*Eom)
	 do J=1,JO    ! Fill in radial drift values
	  SJ=SIN(PHI(J))
	  CJ=COS(PHI(J))
	  Ro=KBETA*(KS(1)+KS(2)*CJ+Kr*(KS(3)+KS(4)*CJ))
	  If (LAMGAM.EQ.0.) Ro=0.
	  EOJ=RE*(KEY*SJ+KEX*CJ)
	  EOJ1=RE*(KEY*CJ-KEX*SJ)
	  do i=1,io
	   RR=LZ(i)+0.5*DL1
	   RoRg=(Ro/RR)**(KGAM-1.)
	   KSD=1.+RoRg*(Ro/RR)
	   VR(I,J)=-DT/DL1*RE/ME*RR**3*KEPS/KSD*(EOJ1+KGAM*SJ/KSD   &
             *(EOJ+KPHI/RR)*RoRg*KBETA/RR*(KS(2)+Kr*KS(4)))
	  end do
	 end do
	 do J=1,JO   ! Fill in azimuthal drift values
	  SJ=SIN(PHI(J)+0.5*DPHI)
	  CJ=COS(PHI(J)+0.5*DPHI)
	  Ro=KBETA*(KS(1)+KS(2)*CJ+Kr*(KS(3)+KS(4)*CJ))
	  If (LAMGAM.EQ.0.) Ro=0.
	  EOJ=RE*(KEY*SJ+KEX*CJ)
	  do i=1,io
	   RR=LZ(i)
	   RoRg=(Ro/RR)**KGAM
	   KSD=1.+RoRg
	   P1(I,J)=DT/DPHI*RE/ME*RR**2*KEPS*EOJ/KSD*(1.+KGAM/KSD   &
             *(1.+KPHI/EOJ/RR)*RoRg)
	  end do
	 end do
	 do j=1,jo          ! Fill in potentials for output
	   SJ=sin(phi(j))
	   CJ=cos(phi(j))
	   EOJ=RE*(KEY*sj+KEX*cj)
	   Ro=KBETA*(KS(1)+KS(2)*CJ+Kr*(KS(3)+KS(4)*CJ))
	   do i=1,io+3
	     RoRg=(Ro/Lsh(i))**KGAM
	     KSD=1.+RoRg
	     BASEPOT(i,j)=KEPS*(Lsh(i)*EOJ+KPHI)/KSD
	  enddo
	 enddo	     
 	 do j=1,nphicells   ! Fill in DGCPM potentials
	  SJ=sind(vphicells(j))
	  CJ=cosd(vphicells(j))
	  EOJ=RE*(KEY*sj+KEX*cj)
	  Ro=KBETA*(KS(1)+KS(2)*CJ+Kr*(KS(3)+KS(4)*CJ))
	  If (LAMGAM.EQ.0) Ro=0.
	  do i=1,nthetacells
	   RoRg=(Ro/vlzcells(i))**KGAM
	   KSD=1.+RoRg
	   potdgcpm(i,j)=KEPS*(vlzcells(i)*EOJ+KPHI)/KSD
	  enddo
	 enddo

!CC Unshielded Volland-Stern field
	ELSE IF (ABASE(IA+1).EQ.4) THEN

	  KEPS=Ppc0/DLMAG
	  do J=1,JO   !  Fill in drift values
	   do i=1,io
	     VR(I,J)=-KEPS*COS(PHI(J))*(LZ(I)+0.5*DL1)**3   &
     		*DT/DL1*(RE/ME)
	     P1(I,J)=KEPS*LZ(I)**2*SIN(PHI(J)+0.5*DPHI)*   &
     		DT/DPHI*(RE/ME)
	   end do
	   do i=1,io+3
	      BASEPOT(i,j)=KEPS*Lsh(i)*SIN(phi(j))
	   end do
	  end do
	 do J=1,nphicells   ! Fill in DGCPM potentials
	  do I=1,nthetacells
	   potdgcpm(i,j)=KEPS*vlzcells(I)*SIND(vphicells(J))
	  enddo
	 enddo
!	print *, 'Conv:',KEPS,VR(IO,1),P1(IO,7)

!CC Base field is specified during the self-consistent calculation
	ELSE IF (ABASE(IA+1).GE.5) THEN
           print *, 'ABASE is 5.'

	  do J=1,JO   !  Zero out drift values
	   do i=1,io
	     VR(I,J) = 0.0
	     P1(I,J) = 0.0
	   end do
	  end do
	 do J=1,nphicells   ! Zero out DGCPM potentials
	  do I=1,nthetacells
	   potdgcpm(i,j) = 0.0
	  enddo
	 enddo
	   
!CC Read in potentials from a file
         print *, 'Is ABASE >=10?'
        ELSE IF (ABASE(IA+1).EQ.10) THEN

	  if (T.EQ.TIME) then
	    TP2=TIME-1.
	    TP1=TP2
	    OPEN(UNIT=18,FILE=name//'_pot.in',status='old')
	    do i=1,5                   ! lines of header material
	      read (18,*) HEADER
	    end do
	    read (18,*) JOpot
	    read (18,*) (MLTpot(j),j=1,JOpot)
	    read (18,*) header
	    read (18,*) IOpot
	    read (18,*) (LZpot(I),i=1,IOpot)
	    do j=JOpot+1,NT
	       MLTpot(j)=MLTpot(j-1)+(MLTpot(JOpot)-MLTpot(JOpot-1))
	    enddo
	    do i=IOpot+1,3*NR
	       LZpot(i)=LZpot(i-1)+(LZpot(IOpot)-LZpot(IOpot-1))
	    enddo
	    FPOT2(1:3*NR,1:NT)=0.
	  end if
	  if (TP2.LT.T) then			! Read in potentials
	    do while (TP2.LE.T)		! Best if final TP2 > final T
	      FPOT1(1:IOpot,1:JOpot)=FPOT2(1:IOpot,1:JOpot)
	      READ (18,*,IOSTAT=L) TP2
	      do j=1,JOpot
		READ (18,*,IOSTAT=L) (DATA(I),I=1,IOpot)
		FPOT2(1:IOpot,J)=DATA(1:IOpot)  ! Potential in Volts
	      end do
	      FPOT2(1:IOpot,JOpot+1)=FPOT2(1:IOpot,1)
	      IF (L.LT.0) TP2=TIME+2*DT*(NSTEP+1)
	      IF (T.EQ.TIME) THEN		! In case T2>T already
		TP1=TIME
		FPOT1(1:IOpot,1:JOpot+1)=FPOT2(1:IOpot,1:JOpot+1)
	      END IF
	    END DO
	  END IF
	  FAC=(T-TP1)/(TP2-TP1)			! Linearly interpolate
	  Fmax=0.
 	  Fmin=0.
	  do j=1,JOpot+1
	    FPOT12(1:IOpot,J)=FAC*FPOT2(1:IOpot,J)+(1.-FAC)*FPOT1(1:IOpot,J)
	    BASEPOT(1:IO+3,J)=FPOT12(1:IO+3,J)
	    Do i=1,IOpot
	       If (FPOT12(i,j).gt.Fmax) Fmax=FPOT12(i,j)
	       If (FPOT12(i,j).lt.Fmin) Fmin=FPOT12(i,j)
	    Enddo
	  enddo
	  print *, 'FPOT12 Max, Min:',Fmin,Fmax,Fmax-Fmin
!         transform FPOT to equatorial plane drifts
	  do j=1,jo
	    SJ=SIN(PHI(J))
	    CJ=COS(PHI(J))
	    jj=j+1
!	    IF (j.eq.jo) jj=1
	    do i=1,io   ! Note this is shifted from FPOT12's I grid by 1
	      RR=LZ(I)+0.5*DL1
	      EPP=-(FPOT12(i+1,jj)-FPOT12(i+1,j))/RR/DPHI
	      VR(i,j)=VR(i,j)+EPP*DT/DL1*RR**3*RE/ME
	    end do
	    SJ=SIN(PHI(J)+.5*DPHI)
	    CJ=COS(PHI(J)+.5*DPHI)
	    do i=1,io
	      RR=LZ(I)
	      EPR=-(FPOT12(i+2,j)-FPOT12(i+1,j))/DL1
	      P1(i,j)=P1(i,j)-EPR*DT/DPHI*RR**2*RE/ME
	    end do
	  end do
!         potdgcpm filled in from FPOT12
	  Fmax=0.
 	  Fmin=0.
	  do j=1,nphicells   !   Fill in DGCPM potentials
	   do i=1,nthetacells
	    CALL LINTP2(LZpot,MLT,FPOT12,3*NR,NT,   &
     	      vlzcells(i),vmltcells(J),FAC,IER)
	    potdgcpm(i,j)=FAC
	   enddo
	   Do i=1,nthetacells
	       If (potdgcpm(i,j).gt.Fmax) Fmax=potdgcpm(i,j)
	       If (potdgcpm(i,j).lt.Fmin) Fmin=potdgcpm(i,j)
	   Enddo
	  enddo
	  print *, 'potdgcpm Max, Min:',Fmin,Fmax,Fmax-Fmin

!CC Calculate potentials using an ionospheric model (AMIE, Weimer, etc.)
        ELSE IF (ABASE(IA+1).GE.11) THEN

!         Transfer date and time to new variables
          eyear=year
          eday=day
          edayplus=AINT(t/86400.)
          eday=eday+edayplus
          univ_time=t-edayplus*86400.
          ehour=AINT(univ_time/3600.)
          univ_time=univ_time-ehour*3600.
          eminute=AINT(univ_time/60.)
          univ_time=univ_time-eminute*60.
          esecond=univ_time

!         Transfer solar wind values to new variables
	  eby=BYSW
	  ebz=BZSW
	  evsw=USW
!         Transfer grids to new variables
!         Note that we have 2 sets of grids, for RAM and for DGCPM
          nemlts1=NT
          nelats1=NR+3
          do i=1,nelats1
            emlts1(1:nemlts1,i)=MLT(1:nemlts1)
          end do
          do j=1,nemlts1
            elats1(j,1:nelats1)=lats(1:nelats1)
          end do
          nemlts2=nphicells
          nelats2=nthetacells
          do i=1,nelats2
             emlts2(1:nemlts2,i)=vmltcells(1:nemlts2)
          end do
          do j=1,nemlts2
             elats2(j,1:nelats2)=90.-vthetacells(1:nelats2)
          end do
!         Transfer the potential model choice to a new variable
!         nechoice = 1 means AMIE potentials
!                  = 2 means Weimer-96 potentials
!                  = 3 means Foster potentials
!                  = ?? we have to add in more to the DATA defs above
	  nechoice=ABASE(IA+1) - 10

!       Aaron: insert your call here, using the "e" variables

        call get_potential(nechoice,eby,ebz,evsw,eyear,eday,ehour,   &
         eminute,esecond,elats1,emlts1,epots1,elats2,emlts2,epots2)

!         Use epots1 to update the advection coefficients:
	  do j=1,jo  ! transform pots1 to equatorial plane drifts
	    jj=j+1
	    IF (j.eq.jo) jj=1
	    BASEPOT(1:io+3,j)=epots1(j,1:io+3)
	    do i=1,io   ! Note this is shifted from epots1's I grid
	      RR=LZ(I)+0.5*DL1
	      EPP=-(epots1(jj,i+1)-epots1(j,i+1))/RR/DPHI
	      VR(i,j)=EPP*DT/DL1*RR**3*RE/ME
	    end do
	    do i=1,io
	      RR=LZ(I)
	      EPR=-(epots1(j,i+2)-epots1(j,i+1))/DL1
	      P1(i,j)=-EPR*DT/DPHI*RR**2*RE/ME
	    end do
	  end do
!         write pots2 into potdgcpm (for Dan Ober's code)
          do j=1,nphicells
             potdgcpm(1:nthetacells,j)=epots2(j,1:nthetacells)
          end do

!CC Done with base field calculations
	END IF

!  Calculate Burke-Wygant penetration electric field
        print *, 'Burke-Wygant pen field?'

	IF (AEPEN(IA+1).EQ.1) THEN
	IF (ABASE(IA+1).EQ.1) THEN
	  PHIP=SQRT(PCO*RE*A*(1.5*DLMAG)**(LAMGAM-1))
	ELSE IF (ABASE(IA+1).GE.2) THEN
	  PHIP=SQRT(PCO*RE*Eom*KEPS)
	END IF
	PHIP=AMAX1(0.,1.4*PHIP-3200.)  ! From Burke's analysis, fig12
!	print *, PHIP,FACP,PCO,RE*Eom,KEPS
	PHIP=-PHIP*FACP    ! To have the Burke value at MLT=18
	do J=1,JO   ! Fill in radial drift values
	 CJ=COS(PHI(J))
	 SJ=SIN(PHI(J))
	 PJ=PHI(J)-PHIPOFF
	 CPJ=ABS(PJ)
	 SPJ=1.
	 IF (exp(-(CPJ/DP1)**3) .lt. exp(-((2.*pi-CPJ)/DP2)**3)) THEN
	   CPJ=2.*PI-CPJ
           IF (PJ.GT.0.) SPJ=-1.
	   DPP=DP2
	 ELSE
	   IF (PJ.LT.0.) SPJ=-1.
	   DPP=DP1
	 ENDIF
	 CPJ=CPJ/DPP
	 LP=0.5+KBETA*(KS(1)+KS(2)*CJ+Kr*(KS(3)+KS(4)*CJ))
	 dLPdphi=-KBETA*SJ*(KS(2)+Kr*KS(4))
	 sLP=sqrt(LP)
	 do i=1,io
	   RR=LZ(I)+0.5*DL1
	   DLP=10.
	   IF (RR.LT.LP) DLP=sLP-1.
	   CRo=(RR-LP)/DLP
	   RoRg=EXP(-CRo**2-CPJ**3)
	   KSD=3.*SPJ*CPJ**2-2.*CRo*dLPdphi/DLP*(1.+CRo/sLP)
	   EPP=PHIP*RoRg/RR*KSD
	   FAC=EPP*DT/DL1*RR**3*RE/ME
	   VR(I,J)=VR(I,J)+FAC
!	if (i.eq.io) print 51,RR,MLT(J),CPJ,CRo,SPJ,RoRg,KSD,EPP,FAC,VR(I,j) 
!	if (i.eq.2 .and. mod(j-1,6).eq.0) print 51,RR,MLT(J),EPP,FAC,VR(I,j) 
	 end do
	end do
	do J=1,JO   ! Fill in azimuthal drift values
	 CJ=COS(PHI(J)+.5*DPHI)
	 PJ=PHI(J)-PHIPOFF+.5*DPHI
	 CPJ=ABS(PJ)
	 IF (exp(-(CPJ/DP1)**3) .lt. exp(-((2.*pi-CPJ)/DP2)**3)) THEN
	   CPJ=2.*PI-CPJ
	   DPP=DP2
	 ELSE
	   DPP=DP1
	 ENDIF
	 CPJ=CPJ/DPP
	 LP=0.5+KBETA*(KS(1)+KS(2)*CJ+Kr*(KS(3)+KS(4)*CJ))
	 sLP=sqrt(LP)
	 do i=1,io
	   RR=LZ(I)
	   DLP=10.
	   IF (RR.LT.LP) DLP=sLP-1.
	   CRo=(RR-LP)/DLP
	   RoRg=EXP(-CRo**2-CPJ**3)
	   KSD=2.*CRo/DLP
	   EPR=PHIP*RoRg*KSD
	   FAC=-EPR*DT/DPHI*RR**2*RE/ME
	   P1(I,J)=P1(I,J)+FAC
!	if (i.eq.io) print 51,RR,MLT(J),CPJ,CRo,PHIP,RoRg,KSD,EPR,FAC,P1(I,j) 
!	if (i.eq.2 .and. mod(j-1,6).eq.0) print 51,RR,MLT(J),EPR,FAC,P1(I,j) 
	 end do
	end do
	do j=1,nphicells   ! Fill in DGCPM potentials
	 CJ=cosd(vphicells(j))
	 PJ=vphicells(j)*pi/180.+phipoff
	 CPJ=ABS(PJ)
	 IF (exp(-(CPJ/DP1)**3) .lt. exp(-((2.*pi-CPJ)/DP2)**3)) THEN
	   CPJ=2.*PI-CPJ
	   DPP=DP2
	 ELSE
	   DPP=DP1
	 ENDIF
	 CPJ=CPJ/DPP
	 LP=0.5+KBETA*(KS(1)+KS(2)*CJ+Kr*(KS(3)+KS(4)*CJ))
	 sLP=sqrt(LP)
	 do i=1,nthetacells
	  RR=vlzcells(i)
	  DLP=10.
	  If (RR.LT.LP) DLP=sLP-1.
	  CRo=(RR-LP)/DLP
	  RoRg=EXP(-CRo**2-CPJ**3)
	  potdgcpm(i,j)=potdgcpm(i,j)+PHIP*RoRg
	 enddo
	enddo
	END IF
51	format(f6.3,f6.2,1p,10e11.3)

!  Calculate self-consistent penetration electric field
        print *, 'Self-consistent pen field?'
	IF (AEPEN(IA+1).GE.2) THEN
           print *, '...Calling PRESSURES'
	  CALL PRESSURES    ! Updates the bulk ring current parameters
          print *, '...Calling CURRENTCALC'
	  CALL CURRENTCALC  ! Finds the perpendicular "ring current"
!			    ! and the field-aligned closure currents
          print *, '...More setup'
	  Irfac=Ir
	  Latfac(1:Ir)=Lats(1:Ir)

	  Ilfac=JO
	  Lonfac(1:JO)=MLT(1:JO)

	  bc_choice=0       ! High-lat boundary has pot=0.
	  if (ABASE(IA+1).eq.5) bc_choice=1  ! W-96 high-lat BC
	  if (ABASE(IA+1).eq.6) bc_choice=2  ! W-96 high-lat BC
	  if (bc_choice.ne.0 .and. AEPEN(IA+1).EQ.3) bc_choice=-bc_choice
!             !! W-96/AMIE everywhere, done in S-C potential call

!CC The Jfac adjustment that's commented out is if you only do 1
!CC species and want to proportionally increase the Jfac value
!CC to reflect the current from all species
!	  NY(2)=0.34*EXP(0.054*KP)	! From Young et al 1982
!	  NY(3)=5.1E-3*EXP(6.6E-3*F107)
!	  NY(4)=0.011*EXP(0.24*KP+0.011*F107)
!	  FAC=(NY(2)+NY(3)+NY(4))/NY(S)
!	  DO j=1,JO
!	    Jfac(1:Ir,j)=Jion1(1:Ir,j)*FAC
!	  end do
!CC Here the Jfac from the various species are summed together

	  if (numprocs.gt.1) then

	     call MPI_BARRIER(MPI_COMM_WORLD,ira)

	     Jfac = 0.0

	     do iproc=1,numprocs

		if (iproc-1.eq.me_world) then
		   jfac_temp(:,:) = Jion1(:,:,parallel_species(iproc))
		endif
		   
		call MPI_Bcast(jfac_temp,(NR+3)*NT,MPI_Real,   &
		     iproc-1,MPI_COMM_WORLD,ira)

		Jion1(:,:,parallel_species(iproc)) = jfac_temp(:,:)

		JFAC(:,:)=Jfac(:,:)+Jreducer*Jion1(:,:,parallel_species(iproc))

	     enddo

	  else

	     DO J=1,JO
		Jfac(1:Ir,J)=0.
		DO S=1,NS
		   IF (SCALC(S).EQ.1) JFAC(1:Ir,J)=Jfac(1:Ir,j)+Jreducer*Jion1(1:Ir,j,s)
		END DO
	     END DO

	  endif

	  if (me_world.eq.0) then

	print *, 'FACs being sent to potential solver'
!	do j=1,JO 
!	print 50, MLT(J),(JFAC(i,j),i=1,Ir)
!	end do
50	format(F5.1,1P,25E12.4)

!cc The next few lines are for doubling convection during sawtooth "snaps"
	  evsw=USW
	  do ii=1,9
 	    tdiff = t - t_sawtooth(ii)
	    if (tdiff.ge.0. .and. tdiff.le.dt_saw) evsw=2.*evsw
	  end do
!cc End the sawtooth snap block
          print *, '...Calling EPENCALC'

	     CALL EPENCALC(t,f107,bc_choice,BYSW,BZSW,evsw) 
!CC                             ! Aaron Ridley's solver 
!CC                             ! for the potential
!CC				! from the field-aligned currents
!CC         ! Note: By and Bz in nT, Usw in m/s

	print *, 'Potentials returned from the solver'
!	do j=1,JO 
!	print 50, MLT(J),(FPOT(i,j),i=1,Ir)
!	end do

	  endif

	  call MPI_BARRIER(MPI_COMM_WORLD,ira)

	  call MPI_Bcast(fpot,(NR+3)*NT,MPI_Real,0,MPI_COMM_WORLD,ira)

! While EPENCALC is commented out, use these lines to zero out FPOT
!	do j=1,jo
!	 do i=1,ir
!	  FPOT(i,j)=0.
!	 enddo
!	end do  ! End FPOT zeroing

	If (AEPEN(IA+1).GT.1) THEN  !! >2 Prohibits inclusion of Epen

	  do j=1,jo  ! transform FPOT to equatorial plane drifts
	    SJ=SIN(PHI(J))
	    CJ=COS(PHI(J))
	    jj=j+1
	    IF (j.eq.jo) jj=1
	    do i=1,io   ! Note this is shifted from FPOT's I grid
	      RR=LZ(I)+0.5*DL1
	      EPP=-(FPOT(i+1,jj)-FPOT(i+1,j))/RR/DPHI
	      VR(i,j)=VR(i,j)+EPP*DT/DL1*RR**3*RE/ME
	    end do
	    SJ=SIN(PHI(J)+.5*DPHI)
	    CJ=COS(PHI(J)+.5*DPHI)
	    do i=1,io
	      RR=LZ(I)
	      EPR=-(FPOT(i+2,j)-FPOT(i+1,j))/DL1
	      P1(i,j)=P1(i,j)-EPR*DT/DPHI*RR**2*RE/ME
	    end do
	  end do
!CC potdgcpm filled in write_ring_current from epencalc
!	  do j=1,nphicells   !   Fill in DGCPM potentials
!	   do i=1,nthetacells
!	    CALL LINTP2(LZ,MLT,FPOT,NR,NT,
!     &	      vlzcells(I),vmltcells(J),FAC,IER)
!	    potdgcpm(i,j)=potdgcpm(i,j)+FAC
!	   enddo
!	  enddo

	END IF  !! Prohibits inclusion of Epen

	END IF  
	print *, 'Done with MAGCONV'

	RETURN
	END
!
! End of subroutine MAGCONV
!


