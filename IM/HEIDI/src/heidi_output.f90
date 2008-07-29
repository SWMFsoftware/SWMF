! File name: heidi_output_010.f90
!
! Contains: output routines for HEIDI
!	ECFL
!	WRESULT
!	PSRCLOSS
!
! Last Modified: March 2006, Mike Liemohn
!
! **********************************************************************
!				ECFL
!	Checks the energy advection CFL numbers
! **********************************************************************

	SUBROUTINE ECFL

	use ModHeidiSize
	use ModHeidiIO
	use ModHeidiMain
	use ModHeidiDrifts

	implicit none

	integer :: i,j,k,l,Im,Km,Jm,Lm,ibad
	real :: cfl,cmax
	real :: acosd
	INTEGER :: NTC
	CHARACTER*5 ST1,ST3
	CHARACTER*3 SUF
	CHARACTER*2 ST2,SUF2
	CHARACTER*1 SUF1
	SAVE NTC

        integer :: iUnitOut = 46

!.......Define parts of the output file names
	IF (T.EQ.0) THEN
	  NTC=0
	  SUF='000'
	ELSE
	  NTC=NTC+1
	  WRITE (SUF,13) NTC
	  IF (NTC.LT.10) THEN
	    WRITE (SUF1,11) NTC
	    SUF='00'//SUF1
	  ELSE IF (NTC.LT.100) THEN
	    WRITE (SUF2,12) NTC
	    SUF='0'//SUF2
	  END IF
	END IF
11	FORMAT (I1)
12	FORMAT (I2)
13	FORMAT (I3)
	ST1=NAME
	ST3='_edr.'

	DO S=1,NS
	IF (SCALC(S).EQ.1) THEN

	IF (S.EQ.1) ST2='_e'
	IF (S.EQ.2) ST2='_h'
	IF (S.EQ.3) ST2='he'
	IF (S.EQ.4) ST2='_o'
	OPEN(iUnitOut,FILE=cOutputDir//ST1//ST2//ST3//SUF,STATUS='UNKNOWN')
	WRITE (iUnitOut,*) 'Filename: '//ST1//ST2//ST3//SUF
	WRITE (iUnitOut,*) 'T,KP,A:',T,KP,A
	WRITE (iUnitOut,38) 'L-shell','MLT(hr)','E(keV)','PA(deg)','n(cm-3)',   &
          'ABS CFL'
	WRITE (iUnitOut,*) 'Energy advection |CFL| > 1 (entire grid)'
	Ibad=0
	CMAX=0.
	Im=0
	Jm=0
	Km=0
	Lm=0
	do I=1,IO
	 do J=1,JO
	  do K=1,KO
	   do L=1,LO	     ! UPA(I)-1 , changed to include the l.c.
	    CFL=ABS(EDOT(I,J,K,L)*VR(I,J)+(COULI(I,K,L,S)+COULE(I,K,L,S))*XNE(I,J))
	    IF (CFL.GT.1.) THEN
	      WRITE (iUnitOut,39) LZ(I),MLT(J),EKEV(K),ACOSD(MU(L)),XNE(I,J),CFL
	     Ibad=Ibad +1
	    END IF
	    IF (CFL.GT.CMAX) THEN
		CMAX=CFL
		Im=I
		Jm=J
		Km=K
		Lm=L
	    END IF
	    IF (Ibad.GT.1000) GOTO 50
	   end do	! L loop
	  end do	! K loop
	 end do	! J loop
	end do	! I loop
50	continue
	WRITE (iUnitOut,*) 'Bad CFLs: ',Ibad
	WRITE (iUnitOut,*) 'Max CFL: ',CMAX,Im,Jm,Km,Lm

	WRITE (iUnitOut,*) 'Radial advection |CFL| > 1 (entire grid)'
	Ibad=0
	CMAX=0.
	do I=1,IO
	 do J=1,JO
	    CFL=ABS(VR(I,J))
	    IF (CFL.GT.1.) THEN
	      WRITE (iUnitOut,39) LZ(I),MLT(J),A,VR(I,J),CFL
	     Ibad=Ibad +1
	    END IF
	    IF (CFL.GT.CMAX) THEN
		CMAX=CFL
		Im=I
		Jm=J
	    END IF
	    IF (Ibad.GT.1000) GOTO 51
	 end do	! J loop
	end do	! I loop
51	continue
	WRITE (iUnitOut,*) 'Bad CFLs: ',Ibad
	WRITE (iUnitOut,*) 'Max CFL: ',CMAX,Im,Jm

	Ibad=0
	CMAX=0.
	WRITE (iUnitOut,*) 'Azimuthal advection |CFL| > 1 (entire grid)'
	do I=1,IO
	 do J=1,JO
	  do K=1,KO
	   do L=1,LO	     
	    CFL=ABS(P1(I,J)+P2(I,K,L))
	    IF (CFL.GT.1.) THEN
	      WRITE (iUnitOut,39) LZ(I),MLT(J),EKEV(K),ACOSD(MU(L)),A,P1(I,J),P2(I,K,L),CFL
	     Ibad=Ibad +1
	    END IF
	    IF (CFL.GT.CMAX) THEN
		CMAX=CFL
		Im=I
		Jm=J
		Km=K
		Lm=L
	    END IF
	    IF (Ibad.GT.1000) GOTO 52
	   end do	! L loop
	  end do	! K loop
	 end do	! J loop
	end do	! I loop
52	continue
	WRITE (iUnitOut,*) 'Bad CFLs: ',Ibad
	WRITE (iUnitOut,*) 'Max CFL: ',CMAX,Im,Jm,Km,Lm

	WRITE (iUnitOut,*) 'Mu advection |CFL| > 1 (entire grid)'
	Ibad=0
	CMAX=0.
	do I=1,IO
	 do J=1,JO
	   do L=1,LO	     
	    CFL=ABS(MUDOT(I,J,L)*VR(I,J))
	    IF (CFL.GT.1.) THEN
	      WRITE (iUnitOut,39) LZ(I),MLT(J),ACOSD(MU(L)),A,CFL
	     Ibad=Ibad +1
	    END IF
	    IF (CFL.GT.CMAX) THEN
		CMAX=CFL
		Im=I
		Jm=J
		Lm=L
	    END IF
	    IF (Ibad.GT.1000) GOTO 53
	   end do	! L loop
	 end do	! J loop
	end do	! I loop
53	continue
	WRITE (iUnitOut,*) 'Bad CFLs: ',Ibad
	WRITE (iUnitOut,*) 'Max CFL: ',CMAX,Im,Jm,Lm

	CLOSE(iUnitOut)

	END IF		! SCALC Check
	END DO		! S LOOP

37	FORMAT(4(2X,A3,2X),7(2X,A7,2X))
38	FORMAT(20(2X,A7,2X))
39	FORMAT(10(1PE11.3))

	RETURN
	END
!
! End of subroutine ECFL
!

! **********************************************************************
!				WRESULT
!       Routine prints all the results at time T after injection
!	IRES(1)  'psd'	F throughout magnetosphere
!	IRES(2)  'etf'	Equ. trapped F throughout magnetosphere
!	IRES(3)  'dep'	Flux tube energy deposition
!	IRES(4)  '*pf'	Total precipitation flux (3 E ranges)
!	IRES(5)  'flx'	Differential precipitation flux
!	IRES(6)  'los'	Particle and energy losses
!	IRES(7)  'pla'	Thermal plasma densities
!	IRES(8)  'cfl'	CFL numbers for advection operators
!	IRES(9)  'drf'	Drift velocities for advection
!	IRES(10) 'evl'	E vs. L distributions at given MLT and PA
!	IRES(11) 'lft'	Particle lifetimes
!	IRES(12) 'prs'	Pressures, densities, and Dst
!	IRES(13) 'fun'  Unformatted output of all F2
!	IRES(14) 'sal'	Continuous sources and losses of number/energy
!	IRES(15) 'fbc'	Nightside boundary condition distribution
!***********************************************************************
	SUBROUTINE WRESULT(LNC,XN,IFIR)

	use ModHeidiSize
	use ModHeidiIO
	use ModHeidiMain
	use ModHeidiDrifts
	use ModHeidiCurrents
	use ModHeidiWaves
	use ModHeidiDGCPM

	implicit none

	integer :: i,j,k,l,nec,nec1,nec2,il,ibad,ii,nlc,ifn,ir2,IFIR
	real :: flux,esum,csum,psum,erate,xlec,xlnc,cfl,weight,xr2,   &
     		edr,xe,xn1,FUNI,FUNT
	EXTERNAL :: FUNI,FUNT
	real :: acosd
	REAL :: LNC(NR,NS),XN(NR,NS),XNO(NR)
     	REAL :: NSUM,TAUD,TAUBO,TAUCHE,EO(NR),AVEFL(NR,NT,NE),RFAC,   &
     	     FZERO(NR,NT,NE),DEP(NR,NT),NBC(NT)
	INTEGER :: NTC,NIC(3),IFM(38),PAV(3),EV(3)
	CHARACTER*5 ST1,ST3,IPF(3)
	CHARACTER*3 SUF
	CHARACTER*2 ST2,SUF2
	CHARACTER*1 SUF1
 	CHARACTER*80 filename
	SAVE NTC

        integer :: iUnitOut = 46

	DATA IPF/'_lpf.','_mpf.','_hpf.'/
	DATA IFM/2,7,13,20,28,35,42,47,50,52,54,56,58,60,62,64,66,68,70,   &
                2,11,21,31,41,51,61,70,75,79,82,83,84,85,86,87,88,89,90/
	DATA PAV,EV/2,20,42,5,19,31/

!.......Define parts of the output file names
	IF (IFIR.EQ.1) THEN
	  NTC=INT(NINT(TIME/TINT))
	ELSE
	  NTC=NTC+1
	END IF
	WRITE (SUF,13) NTC
	IF (NTC.LT.10) THEN
	  WRITE (SUF1,11) NTC
	  SUF='00'//SUF1
	ELSE IF (NTC.LT.100) THEN
	  WRITE (SUF2,12) NTC
	  SUF='0'//SUF2
	END IF
11	FORMAT (I1)
12	FORMAT (I2)
13	FORMAT (I3)
	ST1=NAME

!.......Find output counters
	NLC=NINT(REAL(LO-1)/10.)
	IF (NLC.LT.1) NLC=1
	NEC=NINT(REAL(KO-1)/25.)
	IF (NEC.LT.1) NEC=1

!.......Calculate bulk quantities
        print *, 'Calling PRESSURES'
	CALL PRESSURES
        print *, 'Calling CURRENTCALC'
!	IF (IA.GE.8) CALL CURRENTCALC
	CALL CURRENTCALC  ! Do it all the time
	print *, 'WRESULT DST VALUES: ',(DST(S),S=1,NS)

!.......L counter offset in PAD outputs
	IFN=0
	IF (IPA.EQ.0) IFN=19

!.......Write the plasmaspheric thermal densities (IRES(7), 'pla' & 'dgcpm')
	IF (IRES(7).EQ.1 .and. me_world.eq.0) THEN
	print *, 'Printing plasmasphere'
!	  First create Dan's output file for his plotting software
	  filename=cOutputDir//ST1//'_dgcpm_'//SUF//'.dat'
	call getdensity(vthetacells,nthetacells,vphicells,nphicells,   &
              dendgcpm)
!	call saveit(vthetacells,nthetacells,vphicells,nphicells,   &
!              dendgcpm,gridx,gridy,gridoc,filename)
	call saveplasmasphere(filename)
!	  Next create an output like we have made before
	  ST3='_pla.'
        OPEN(UNIT=iUnitOut,file=cOutputDir//ST1//ST3//SUF,STATUS='UNKNOWN')
	  WRITE (iUnitOut,*) 'Filename: '//ST1//ST3//SUF
	  WRITE (iUnitOut,*) 'Plasmaspheric thermal densities from the DGCPM'
        WRITE (iUnitOut,16) T,KP
	  WRITE (iUnitOut,31) (MLT(J),J=1,JO)
	  do I=2,IO
	    WRITE(iUnitOut,29) LZ(I),(XNE(I,J),J=1,JO)
	  end do
	  CLOSE (iUnitOut)
	END IF

!CC Output from Aaron's ionosphere code
!CC Didn't work...we need Aaron to diagnose the problem
!	print *, 'Saving Aaron''s ionosphere model output'
	if ((IA.ge.8 .and. IA.le.11) .or. IA.ge.13) then
	 if (me_world.eq.0) call IonoHeidiWriteOutput(1,t,ST1,SUF)
	end if

!CC Start the main loop over species we're calculating
	DO S=1,NS
!	print *, 'WRESULT: ',S,SCALC(S)
	IF (SCALC(S).EQ.1) THEN

	IF (S.EQ.1) ST2='_e'
	IF (S.EQ.2) ST2='_h'
	IF (S.EQ.3) ST2='he'
	IF (S.EQ.4) ST2='_o'

!.......Find the energy and particle losses
	DO I=2,IO
	  XNO(I)=XN(I,S)
	  XN(I,S)=0
	  EO(I)=ENER(I,S)
	  ENER(I,S)=0.
	  do K=2,KO
	    do L=2,LO-1
	      DO J=1,JO
		  IF(L.LT.UPA(I)) THEN
		   WEIGHT=F2(I,J,K,L,S)*WE(K)*WMU(L)
	 	   XN(I,S)=XN(I,S)+WEIGHT  		      ! N in LZ
		   ENER(I,S)=ENER(I,S)+EKEV(K)*WEIGHT     ! E in LZ
		  ENDIF
	      END DO
	    end do	! K loop
	  end do	! L loop
	  LNC(I,S)=XNO(I)-XN(I,S)
	  LEC(I,S)=EO(I)-ENER(I,S)
	end do

!.......Start the output routines

!.......Write the phase space distribution, F (IRES(1), 'psd')
!	IF (MOD(T,21600.).LT.2*DT) THEN	! Only every 6 hours
	IF (IRES(1).EQ.1) THEN
	  ST3='_psd.'
        OPEN(UNIT=iUnitOut,file=cOutputDir//ST1//ST2//ST3//SUF,STATUS='UNKNOWN')
	  WRITE (iUnitOut,*) 'Filename: '//ST1//ST2//ST3//SUF
	  IF (IFAC.EQ.1) THEN
	    WRITE (iUnitOut,*) 'Phase space flux function, PHI = 2EF/m^2'
	  ELSE
	    WRITE (iUnitOut,*) 'Phase space distribution function, F'
	  END IF
        do I=2,IO,2  ! ,4
	    do J=1,JO  ! ,3
		WRITE(iUnitOut,45) T,LZ(I),MLT(J),KP,XNE(I,J)
		WRITE(iUnitOut,44) (ACOSD(MU(IFM(L))),L=1+IFN,19+IFN)
		do K=2,KO,NEC
		  WRITE(iUnitOut,43) EKEV(K),(F2(I,J,K,IFM(L),S)/   &
     		    FFACTOR(I,K,IFM(L)),L=1+IFN,19+IFN)
		end do	! K loop
	    end do		! J loop
	  end do		! I loop
	  CLOSE(iUnitOut)
	END IF
!	END IF					! 6 hour check

!.......Write equatorially trapped distribution (IRES(2), 'etf')
	IF (IRES(2).EQ.1) THEN
        ST3='_etf.'
        OPEN(UNIT=iUnitOut,FILE=cOutputDir//ST1//ST2//ST3//SUF,STATUS='UNKNOWN')
	  WRITE (iUnitOut,*) 'Filename: '//ST1//ST2//ST3//SUF
	  WRITE (iUnitOut,*) 'Equatorially trapped distribution function, F'
	  NIC(1)=2
	  NIC(2)=IO/2
	  NIC(3)=IO
	  do II=1,3
	    I=NIC(II)
	    WRITE(iUnitOut,33) T,LZ(I),KP,ACOSD(MU(2))
	    WRITE(iUnitOut,31) (MLT(J),J=1,JO,3)
	    do K=1,KO,NEC
		WRITE(iUnitOut,29) EKEV(K),(F2(I,J,K,2,S)/FFACTOR(I,K,L),   &
     		  J=1,JO,3)
	    end do
	  end do
	  CLOSE(iUnitOut)
	END IF

!.......Write the plasmaspheric heating (IRES(3), 'dep')
	IF (IRES(3).EQ.1) THEN
	 XR2=IO/2.
	 IR2=IO/2
	 ST3='_dep.'
	 OPEN(UNIT=iUnitOut,file=cOutputDir//ST1//ST2//ST3//SUF,STATUS='UNKNOWN')
	 WRITE (iUnitOut,*) 'Filename: '//ST1//ST2//ST3//SUF
	 do I=2,IO			! Electron heating
	  do J=1,JO
	    EDR=0.
	    do K=2,KO
	      FZERO(I,J,K)=0.
	      do L=2,UPA(I)-1
		  FZERO(I,J,K)=FZERO(I,J,K)+F2(I,J,K,L,S)*CEDR(K,L,S)   &
      		    /FFACTOR(I,K,L)
		end do
	      do L=UPA(I),LO
		  FZERO(I,J,K)=FZERO(I,J,K)+F2(I,J,K,L,S)   &
     		    *CEDR(K,UPA(I)-1,S)/FFACTOR(I,K,L)
	      end do
	      IF (IFAC.EQ.1) FZERO(I,J,K)=FZERO(I,J,K)/EKEV(K)/FLUXFACT(S)
	      EDR=EDR+FZERO(I,J,K)*WE(K)
	    end do
	    EDR=EDR+EKEV(2)*FZERO(I,J,2)-EKEV(KO)*FZERO(I,J,KO)
	    DEP(I,J)=EDR*XNE(I,J)*LZ(I)
	  end do
	 end do
	 WRITE (iUnitOut,*) 'Energy deposition rates of thermal electrons'
	 WRITE (iUnitOut,*) 'through Coulomb collisions [eV/cm2/s]'
         WRITE (iUnitOut,15) T,KP
	 WRITE (iUnitOut,31) (LZ(I),I=2,IO)
	 do J=1,JO
	   WRITE (iUnitOut,29) MLT(J),(DEP(I,J),I=2,IO)
	 end do
	 do I=2,IO			! Ion heating
	  do J=1,JO
	    EDR=0.
	    do K=2,KO
	      FZERO(I,J,K)=0.
	      do L=2,UPA(I)-1
		  FZERO(I,J,K)=FZERO(I,J,K)+F2(I,J,K,L,S)*CIDR(K,L,S)   &
     		    /FFACTOR(I,K,L)
		end do
	      do L=UPA(I),LO
		  FZERO(I,J,K)=FZERO(I,J,K)+F2(I,J,K,L,S)   &
     		    *CIDR(K,UPA(I)-1,S)/FFACTOR(I,K,L)
	      end do
	      IF (IFAC.EQ.1) FZERO(I,J,K)=FZERO(I,J,K)/EKEV(K)/FLUXFACT(S)
	      EDR=EDR+FZERO(I,J,K)*WE(K)
	    end do
	    EDR=EDR+EKEV(2)*FZERO(I,J,2)-EKEV(KO)*FZERO(I,J,KO)
	    DEP(I,J)=EDR*XNE(I,J)*LZ(I)
	  end do
	 end do
	 WRITE (iUnitOut,*) 'Heating rates of thermal ions'
	 WRITE (iUnitOut,*) 'through Coulomb collisions [eV/cm2/s]'
         WRITE (iUnitOut,15) T,KP
	 WRITE (iUnitOut,31) (LZ(I),I=2,IO)
	 do J=1,JO
	   WRITE (iUnitOut,30) MLT(J),(DEP(I,J),I=2,IO)
	 end do
	 CLOSE(iUnitOut)
	END IF

!.......Write the precipitation flux (IRES(4), '*pf')
	if (IRES(4).EQ.1) then	
	 do II=1,3
	  ST3=IPF(II)
	  NEC1=(II-1)*KO/3+1
	  NEC2=II*KO/3
          OPEN(UNIT=iUnitOut,file=cOutputDir//ST1//ST2//ST3//SUF,STATUS='UNKNOWN')
	  WRITE (iUnitOut,*) 'Filename: '//ST1//ST2//ST3//SUF
	  WRITE (iUnitOut,*) 'Total precipitation flux [1/cm2/s]'
	  WRITE (iUnitOut,*) 'Integrated over the energy range:'
	  WRITE (iUnitOut,*) 'Lower edge (keV): ',EKEV(NEC1)
	  WRITE (iUnitOut,*) 'Upper edge (keV): ',EKEV(NEC2)
	  WRITE(iUnitOut,71) T,KP
	  WRITE(iUnitOut,72)
	  do I=2,IO
	    do J=1,JO
	      FLUX=0.
	      do K=NEC1,NEC2
		AVEFL(I,J,K)=0.
		do L=UPA(I),LO-1
		  AVEFL(I,J,K)=AVEFL(I,J,K)+F2(I,J,K,L,S)*WMU(L)   &
     		    /FFACTOR(I,K,L)
		end do ! L loop
  	        AVEFL(I,J,K)=AVEFL(I,J,K)/(MU(LO)-MU(UPA(I)))
                IF (IFAC.EQ.2) AVEFL(I,J,K)=AVEFL(I,J,K)*FLUXFACT(S)*EKEV(K)
		FLUX=FLUX+AVEFL(I,J,K)*PI*WE(K)
	      end do ! K loop
	      WRITE(iUnitOut,70) LZ(I),PHI(J),FLUX
	    end do	! J loop
	  end do	! I loop
	  CLOSE (iUnitOut)
	 end do
	endif

!.......Write the differential precip flux (IRES(5), 'flx')
	IF (IRES(5).EQ.1) THEN	
	  ST3='_flx.'
          OPEN(UNIT=iUnitOut,file=cOutputDir//ST1//ST2//ST3//SUF,STATUS='UNKNOWN')
	  WRITE (iUnitOut,*) 'Filename: '//ST1//ST2//ST3//SUF
	  WRITE (iUnitOut,*) 'Differential precipitation fluxes'
          DO I=4,IO,2
!	   IF(I.EQ.4.OR.I.EQ.6.OR.I.EQ.10.OR.I.EQ.16) THEN
	    IL=2
!	    IF(I.EQ.16) IL=1
	    DO J=1,JO,3
!	     IF(J.EQ.1.OR.J.EQ.10.OR.J.EQ.16) THEN
	      WRITE(iUnitOut,34) T,LZ(I),MLT(J)
	      WRITE(iUnitOut,40) (ACOSD(MU(L)),L=UPA(I),LO-1,IL)
	      WRITE(iUnitOut,*)'     EKEV  \ FLUX[1/cm2/s/ster/keV]  AVEFL'
	      IF (IFAC.EQ.1) THEN
 	        DO K=2,KO,NEC
		  WRITE(iUnitOut,30) EKEV(K),(F2(I,J,K,L,S)/FFACTOR(I,K,L),   &
                     L=UPA(I),LO-1,IL),AVEFL(I,J,K)
	        END DO
	      ELSE
 	        DO K=2,KO,NEC
		  WRITE(iUnitOut,30) EKEV(K),(F2(I,J,K,L,S)*FLUXFACT(S)*EKEV(K)   &
     		    /FFACTOR(I,K,L),L=UPA(I),LO-1,IL),AVEFL(I,J,K)
	        END DO
	      END IF
!	     END IF
            END DO
!	   END IF
          END DO
	  CLOSE (iUnitOut)
        END IF

!..Write the particle & energy losses (IRES(6), 'los')
	IF (IRES(6).EQ.1) THEN
	  ST3='_los.'
	  OPEN(UNIT=iUnitOut,file=cOutputDir//ST1//ST2//ST3//SUF,STATUS='UNKNOWN')
	  WRITE (iUnitOut,*) 'Filename: '//ST1//ST2//ST3//SUF
	  WRITE (iUnitOut,*) 'Particle and energy losses since last output'
	  WRITE(iUnitOut,71) T,KP
	  ESUM=0
	  NSUM=0
	  CSUM=0
	  PSUM=0
	  ERATE=0
        WRITE(iUnitOut,35)
        do I=2,IO
          XE=ENER(I,S)*FACTOR(S)
          XN1=XN(I,S)*FACTOR(S)
          XLEC=LEC(I,S)*FACTOR(S)
	    XLNC=LNC(I,S)*FACTOR(S)
          ESUM=ESUM+XE
          NSUM=NSUM+XN1
	    CSUM=CSUM+XLEC
	    PSUM=PSUM+XLNC
	    WRITE(iUnitOut,30) LZ(I),XE,XN1,XLEC,XLNC
	    do J=1,JO
	      do K=2,KO
	        ERATE=ERATE+EKEV(K)*AVEFL(I,J,K)*WE(K)
		end do
	    end do
	    ERATE=ERATE*2*PI*Q*ECOF(I)*1e7	!to convert in cm and J
	  end do
	  WRITE(iUnitOut,60) ESUM,NSUM,CSUM,PSUM
	  WRITE(iUnitOut,61) ERATE
	  CLOSE (iUnitOut)
	END IF

!.......Print out the CFLs for the advection operators (IRES(8), 'cfl')
	IF (IRES(8).EQ.1) THEN
	 ST3='_cfl.'
	 OPEN (iUnitOut,FILE=cOutputDir//ST1//ST2//ST3//SUF,STATUS='UNKNOWN')
	 WRITE (iUnitOut,*) 'Filename: '//ST1//ST2//ST3//SUF
	 WRITE (iUnitOut,*) 'CFL numbers for the advection operators'
	 WRITE (iUnitOut,71) T,KP
	 WRITE (iUnitOut,*)
	 WRITE (iUnitOut,*) 'Coulomb collsion Courant numbers (with e-)'
	 do I=2,IO,IO-2
	  do J=1,JO,3
	    WRITE (iUnitOut,37) LZ(I),MLT(J),XNE(I,J)
	    WRITE (iUnitOut,40) (ACOSD(MU(L)),L=2,LO,NLC),ACOSD(MU(LO-1))
	    do K=2,KO,KO-2
	      WRITE (iUnitOut,30) EKEV(K),(COULE(I,K,L,S)*XNE(I,J),L=2,LO,NLC),   &
     		COULE(I,K,LO-1,S)*XNE(I,J)
	    end do
	  end do
	 end do
	 WRITE (iUnitOut,*)
	 WRITE (iUnitOut,*) 'Coulomb collsion Courant numbers (with ions)'
	 do I=2,IO,IO-2
	  do J=1,JO,3
	    WRITE (iUnitOut,37) LZ(I),MLT(J),XNE(I,J)
	    WRITE (iUnitOut,40) (ACOSD(MU(L)),L=2,LO,NLC),ACOSD(MU(LO-1))
	    do K=2,KO,KO-2
	      WRITE (iUnitOut,30) EKEV(K),(COULI(I,K,L,S)*XNE(I,J),L=2,LO,NLC),   &
     		COULI(I,K,LO-1,S)*XNE(I,J)
	    end do
	  end do
	 end do
	 WRITE (iUnitOut,*)
	 WRITE (iUnitOut,*) 'Radial drift Courant numbers'
	 WRITE (iUnitOut,41) ' L \ MLT =',(MLT(J),J=1,JO,3)
	 do I=2,IO,4
	  WRITE (iUnitOut,42) LZ(I),(VR(I,J),J=1,JO,3)
	 end do
	 WRITE (iUnitOut,*)
	 WRITE (iUnitOut,*) 'Azimuthal drift Courant numbers'
	 do I=2,IO,IO-2
	  do J=1,JO,3
	    WRITE (iUnitOut,36) LZ(I),MLT(J)
	    WRITE (iUnitOut,40) (ACOSD(MU(L)),L=2,LO,NLC),ACOSD(MU(LO-1))
	    do K=2,KO,KO-2
	      WRITE (iUnitOut,30) EKEV(K),(P1(I,J)+P2(I,K,L),L=2,LO,NLC),   &
               P1(I,J)+P2(I,K,LO-1)
	    end do
	  end do
	 end do
	 WRITE (iUnitOut,*)
	 WRITE (iUnitOut,*) 'Energy drift Courant numbers'
	 do I=2,IO,IO-2
	  do J=1,JO,3
	    WRITE (iUnitOut,36) LZ(I),MLT(J)
	    WRITE (iUnitOut,40) (ACOSD(MU(L)),L=2,LO,NLC),ACOSD(MU(LO-1))
	    do K=2,KO,KO-2
	      WRITE (iUnitOut,30) EKEV(K),(VR(I,J)*EDOT(I,J,K,L),L=2,LO,NLC),   &
             VR(I,J)*EDOT(I,J,K,LO-1)
	    end do
	  end do
	 end do
	 WRITE (iUnitOut,*)
	 WRITE (iUnitOut,*) 'Mu drift Courant numbers'
	 do I=2,IO,IO-2
	  WRITE (iUnitOut,*) ' L =',LZ(I)
	  WRITE (iUnitOut,40) (ACOSD(MU(L)),L=2,LO,NLC),ACOSD(MU(LO-1))
	  do J=1,JO,3
	   WRITE (iUnitOut,30) MLT(J),(VR(I,J)*MUDOT(I,J,L),L=2,LO,NLC),   &
           VR(I,J)*MUDOT(I,J,LO-1)
	  end do
	 end do
	 WRITE (iUnitOut,*) 'Energy advection |CFL| > 1 (entire grid)'
	 WRITE (iUnitOut,38) 'L-shell','MLT(hr)','E(keV)','PA(deg)','n(cm-3)',   &
         'ABS CFL'
	 Ibad=0
	 do I=1,IO
	  do J=1,JO
	   do K=1,KO
	    do L=1,LO	     ! UPA(I)-1 , changed to include the l.c.
	     CFL=ABS(EDOT(I,J,K,L)*VR(I,J)+   &
     		(COULI(I,K,L,S)+COULE(I,K,L,S))*XNE(I,J))
	     if (CFL.GT.1.) then
	      WRITE (iUnitOut,39)LZ(I),MLT(J),EKEV(K),ACOSD(MU(L)),XNE(I,J),CFL
	      Ibad=Ibad +1
	     end if
	    end do	! L loop
	   end do	! K loop
	  end do	! J loop
	 end do	! I loop
	 WRITE (iUnitOut,*) 'Bad CFLs: ',Ibad
	 CLOSE(iUnitOut)
	END IF

!.......Print out the drift velocities (IRES(9), 'drf')
	IF (IRES(9).EQ.1) THEN
	 ST3='_drf.'
	 OPEN (iUnitOut,FILE=cOutputDir//ST1//ST2//ST3//SUF,STATUS='UNKNOWN')
	 WRITE (iUnitOut,*) 'Filename: '//ST1//ST2//ST3//SUF
	 WRITE (iUnitOut,*) 'Spatial coordinate drift velocities'
	 WRITE (iUnitOut,71) T,KP
	 WRITE (iUnitOut,*)
	 WRITE (iUnitOut,*) 'Radial drift component'
	 WRITE (iUnitOut,41) ' L \ MLT =',(MLT(J),J=1,JO,2)
	 do I=2,IO,2
	  WRITE (iUnitOut,42) LZ(I),(VR(I,J),J=1,JO,2)
	 end do
	 WRITE (iUnitOut,*)
	 WRITE (iUnitOut,*) 'Azimuthal drift component at various E,mu'
	 do L=1,3
	  do K=1,3
	    WRITE (iUnitOut,36) EKEV(EV(K)),ACOSD(MU(PAV(L)))
	    WRITE (iUnitOut,41) ' L \ MLT =',(MLT(J),J=1,JO,2)
	    do I=2,IO,2
	     WRITE (iUnitOut,42) LZ(I),(P1(I,J)+P2(I,EV(K),PAV(L)),J=1,JO,2)
	    end do
	  end do
	 end do
	 CLOSE(iUnitOut)
	END IF

!.......Print out midnight energy vs. L distributions (IRES(10), 'evl')
!.......Also noon distributions, at PA=90 deg and at UPA(I) boundary
	IF (IRES(10).EQ.1) THEN
	 ST3='_evl.'
	 OPEN (iUnitOut,FILE=cOutputDir//ST1//ST2//ST3//SUF,STATUS='UNKNOWN')
	 WRITE (iUnitOut,*) 'Filename: '//ST1//ST2//ST3//SUF
	 WRITE (iUnitOut,*) 'Energy-Lshell spectra at MLT=0 and PA=90'
	 WRITE (iUnitOut,71) T,KP
	 WRITE(iUnitOut,41) '   Ne =',(XNE(I,1),I=2,IO)
	 WRITE (iUnitOut,41) ' E \ L =',(LZ(I),I=2,IO)
	 do K=2,KO
	   WRITE(iUnitOut,43) EKEV(K),(F2(I,1,K,2,S)/FFACTOR(I,K,2),I=2,IO)
	 end do		! K loop
	 WRITE (iUnitOut,*) 'Energy-Lshell spectra at MLT=12 and PA=90'
	 WRITE (iUnitOut,71) T,KP
	 WRITE(iUnitOut,41) '   Ne =',(XNE(I,13),I=2,IO)
	 WRITE (iUnitOut,41) ' E \ L =',(LZ(I),I=2,IO)
	 do K=2,KO
	   WRITE(iUnitOut,43) EKEV(K),(F2(I,13,K,2,S)/FFACTOR(I,K,2),I=2,IO)
	 end do		! K loop
	 WRITE (iUnitOut,*) 'Energy-Lshell spectra at MLT=0 and PA=UPA'
	 WRITE (iUnitOut,71) T,KP
	 WRITE(iUnitOut,41) '   Ne =',(XNE(I,1),I=2,IO)
	 WRITE (iUnitOut,41) ' E \ L =',(LZ(I),I=2,IO)
	 do K=2,KO
	   WRITE(iUnitOut,43) EKEV(K),(F2(I,1,K,UPA(I)-1,S)/   &
     	     FFACTOR(I,K,UPA(I)-1),I=2,IO)
	 end do		! K loop
	 WRITE (iUnitOut,*) 'Energy-Lshell spectra at MLT=12 and PA=UPA'
	 WRITE (iUnitOut,71) T,KP
	 WRITE(iUnitOut,41) '   Ne =',(XNE(I,13),I=2,IO)
	 WRITE (iUnitOut,41) ' E \ L =',(LZ(I),I=2,IO)
	 do K=2,KO
	   WRITE(iUnitOut,43) EKEV(K),(F2(I,13,K,UPA(I)-1,S)/   &
             FFACTOR(I,K,UPA(I)-1),I=2,IO)
	 end do		! K loop
	 CLOSE(iUnitOut)
	END IF

!.......Print lifetimes [hr] and Coulomb diff coeff (IRES(11), 'lft')
	IF (IRES(11).EQ.1) THEN
	 ST3='_lft.'
	 OPEN(UNIT=iUnitOut,FILE=cOutputDir//ST1//ST2//ST3//SUF,STATUS='UNKNOWN')
	 WRITE (iUnitOut,*) 'Filename: '//ST1//ST2//ST3//SUF
	 WRITE (iUnitOut,71) T,KP
	DO  I=4,IO,6
         L=UPA(I)-1
	 DO  J=1,1
          WRITE(iUnitOut,*)' Lifetimes for ',ST2,' rc ion; PA=',   &
            ACOSD(MU(L))
          WRITE(iUnitOut,*)' L   E[KEV]  TAUBO[HR] TAUCHE[HR] TAUD[HR]'
	  DO K=2,KO
	   TAUD=2*pi*ME/ABS(1.44E-2*RE-3*EKEV(K)*1000*LZ(I))/RE/   &
      		(1-FUNI(MU(L))/6/FUNT(MU(L)))/3600.
           TAUBO=4*LZ(I)*RE/V(K,S)*FUNT(MU(L))/3600.
	   TAUCHE=-DT/ALOG(ACHAR(I,K,L,S))/3600.
	   WRITE(iUnitOut,80) LZ(I),EBND(K),TAUBO,TAUCHE,TAUD
	  ENDDO
	 ENDDO
	ENDDO
	CLOSE (iUnitOut)
	END IF

!.......Print out pressures, densities, and Dst (IRES(12), 'prs')
	IF (IRES(12).EQ.1) THEN
	 ST3='_prs.'
	 OPEN (iUnitOut,FILE=cOutputDir//ST1//ST2//ST3//SUF,STATUS='UNKNOWN')
	 WRITE (iUnitOut,*) 'Filename: '//ST1//ST2//ST3//SUF
	 WRITE (iUnitOut,*) 'Pressures, densities, etc. for RC species ',ST2
	 WRITE (iUnitOut,71) T,KP
	 WRITE (iUnitOut,*) 'Total Energy [keV]   Total Particles   Dst [nT]'
	 WRITE (iUnitOut,73) ETOT(S),NTOT(S),Dst(S)
	print *, 'WRESULT Dst: ',S,Dst(S)
	 WRITE (iUnitOut,*) 'Equatorial density [cm-3]'
	 WRITE (iUnitOut,31) (LZ(I),I=2,IO)
	 do J=1,JO
	   WRITE (iUnitOut,29) MLT(J),(RNHT(I,J,S),I=2,IO)
	 end do
	 WRITE (iUnitOut,*) 'Equatorial perpendicular pressure [keV cm-3]'
	 WRITE (iUnitOut,31) (LZ(I),I=2,IO)
	 do J=1,JO
	   WRITE (iUnitOut,29) MLT(J),(PPER(I,J,S),I=2,IO)
	 end do
	 WRITE (iUnitOut,*) 'Equatorial parallel pressure [keV cm-3]'
	 WRITE (iUnitOut,31) (LZ(I),I=2,IO)
	 do J=1,JO
	   WRITE (iUnitOut,29) MLT(J),(PPAR(I,J,S),I=2,IO)
	 end do
	 WRITE (iUnitOut,*) 'Equatorial anisotropy [Tper/Tpar - 1]'
	 WRITE (iUnitOut,31) (LZ(I),I=2,IO)
	 do J=1,JO
	   WRITE (iUnitOut,29) MLT(J),(ANIS(I,J,S),I=2,IO)
	 end do
	 WRITE (iUnitOut,*) 'Equatorial energy density [keV cm-3]'
	 WRITE (iUnitOut,31) (LZ(I),I=2,IO)
	 do J=1,JO
	   WRITE (iUnitOut,29) MLT(J),(EDEN(I,J,S),I=2,IO)
	 end do
	 WRITE (iUnitOut,*) 'Equatorial azimuthal current [A m-2]'
	 WRITE (iUnitOut,31) (LZ(I),I=2,IO)
	 do J=1,JO
	   WRITE (iUnitOut,29) MLT(J),(JPER(I,J,S),I=2,IO)
	 end do
	 WRITE (iUnitOut,*) 'Total particle count in the spatial volume [ions]'
	 WRITE (iUnitOut,31) (LZ(I),I=2,IO)
	 do J=1,JO
	   WRITE (iUnitOut,29) MLT(J),(Nspace(I,J,S),I=2,IO)
	 end do
	 WRITE (iUnitOut,*) 'Total energy in the spatial volume [keV]'
	 WRITE (iUnitOut,31) (LZ(I),I=2,IO)
	 do J=1,JO
	   WRITE (iUnitOut,29) MLT(J),(Espace(I,J,S),I=2,IO)
	 end do
	 WRITE (iUnitOut,*) 'Base convection potentials [kV]'
	 WRITE (iUnitOut,31) (Lsh(I),I=1,Ir)
	 do J=1,JO
	   WRITE (iUnitOut,29) MLT(J),(BASEPOT(I,J)*1.E-3,I=1,Ir)
	 end do
!	 IF (IA.GE.8) THEN
	   WRITE (iUnitOut,*) 'Total azimuthal current [A]'
	   WRITE (iUnitOut,31) (Lsh(I),I=1,Ir)
	   do J=1,JO
	     WRITE (iUnitOut,29) MLT(J),(Iphi(I,J,S),I=1,Ir)
	   end do
	   WRITE (iUnitOut,*) 'Total radial current [A]'
	   WRITE (iUnitOut,31) (Lsh(I),I=1,Ir)
	   do J=1,JO
	     WRITE (iUnitOut,29) MLT(J),(Irad(I,J,S),I=1,Ir)
	   end do
	   WRITE (iUnitOut,*) 'Field-aligned current density into one ',   &
             'hemisphere [A m-2]'
	   WRITE (iUnitOut,31) (Lsh(I),I=1,Ir)
	   do J=1,JO
	     WRITE (iUnitOut,29) MLT(J),(Jion1(I,J,S),I=1,Ir)
	   end do
	   WRITE (iUnitOut,*) 'Potentials from RC FACs [kV]'
	   WRITE (iUnitOut,31) (Lsh(I),I=1,Ir)
	   do J=1,JO
	     WRITE (iUnitOut,29) MLT(J),(FPOT(I,J)*1.E-3,I=1,Ir)
	   end do
	   WRITE (iUnitOut,*) 'FAC density from all species into one ',   &
             'hemisphere [A m-2]'
	   WRITE (iUnitOut,31) (Lsh(I),I=1,Ir)
	   do J=1,JO
	     WRITE (iUnitOut,29) MLT(J),(Jfac(I,J),I=1,Ir)
	   end do
!	 END IF
	 CLOSE (iUnitOut)
	END IF

!.......Write out F2 to an unformatted file (IRES(13), 'unff')
	IF (IRES(13).EQ.1) THEN
	 ST3='.unff'
	 OPEN(UNIT=iUnitOut,FILE=cOutputDir//ST1//ST2//ST3,status='unknown',   &
           form='unformatted')
	 DO L=1,NPA
	  DO K=1,NE
	   DO J=1,NT
	    write(iUnitOut) (f2(I,J,K,L,S),I=1,NR)
	   END DO
	  END DO
	 END DO
	 close(iUnitOut)
	END IF

!.......Open file for source/loss continual output (IRES(14), 'sal')
	 IF (IRES(14).EQ.1) THEN
!	print *, 'Continuous loss file:',15+S,ST1//ST2//ST3//SUF
	  CLOSE (15+S)
	  ST3='_sal.'
	  OPEN (15+S,FILE=cOutputDir//ST1//ST2//ST3//SUF,STATUS='UNKNOWN')
	  WRITE (15+S,*) 'Filename: '//ST1//ST2//ST3//SUF
	  WRITE (15+S,*) 'Sources and losses: continuous output'
	  WRITE (15+S,71) T,KP
	  WRITE (15+S,74) 'T','LMP6','LMP12','RNS','RNL','RES','REL',   &
            'ESN','ELN','ESE','ELE','ECN','ECE','ALN','ALE','CEN',   &
            'CEE','DNT','DET'
	END IF

!.......Print out nightside boundary condition (IRES(15), 'fbc')
	IF (IRES(15).EQ.1) THEN
	 NBC(1:JO)=0.
	 RFAC=4.*PI*SQRT(2.)*1.E-24*(Q*1.E3/MP/M1(S))**1.5
	 DO L=1,LO
	  DO K=2,KO
	   DO J=1,JO
	     NBC(J)=NBC(J)+FGEOS(J,K,L,S)*ERNM(IO,K,L)*ERNH(K,S)*RFAC
	   END DO
	  END DO
	 END DO
	 ST3='_fbc.'
	 OPEN (iUnitOut,FILE=cOutputDir//ST1//ST2//ST3//SUF,STATUS='UNKNOWN')
	 WRITE (iUnitOut,*) 'Filename: '//ST1//ST2//ST3//SUF
	 WRITE (iUnitOut,*) 'Nightside boundary conditions'
	 WRITE (iUnitOut,38) 'Ninj','Einj','Kinj','NSWB','USWB'
	 WRITE (iUnitOut,39) Ninj,Einj,Kinj,NSWB,USWB
	 WRITE (iUnitOut,*) 'Boundary densities'
	 WRITE (iUnitOut,38) 'MLT','N,cm-3'
	 DO J=1,JO
	  WRITE (iUnitOut,30) MLT(J),NBC(J)
	 END DO
	 IF (IFAC.EQ.1) THEN
	   WRITE (iUnitOut,*) 'Phase space flux function, PHI = 2EF/m^2'
	 ELSE
	   WRITE (iUnitOut,*) 'Phase space distribution function, F'
	 END IF
	 do J=-3,5,2
	   II=J
	   IF (II.LT.0) II=II+JO
	   WRITE(iUnitOut,45) T,LZ(IO)+DL1,MLT(II),KP,0.
	   WRITE(iUnitOut,44) (ACOSD(MU(IFM(L))),L=1+IFN,19+IFN)
	   do K=2,KO,NEC
	     WRITE(iUnitOut,43) EKEV(K),(FGEOS(II,K,IFM(L),S)/   &
       	       FFACTOR(IO,K,IFM(L)),L=1+IFN,19+IFN)
	   end do	! K loop
	 end do		! J loop
	 CLOSE (iUnitOut)
	END IF

	END IF		! SCALC check
	END DO		! S loop

15      FORMAT(' EKEV \ T =',F8.0,2X,'Kp =',F6.2,10X)
16	FORMAT('LSHELL\ T =',F8.0,2X,'Kp =',F6.2,10X)
29	FORMAT(F7.3,1P,40E11.3)
30	FORMAT(F7.3,1P,20E10.3)
31	FORMAT(' MLT\L ',40(3X,F8.2))
32	FORMAT(' EKEV \ T=',F8.0,' L=',F6.2,' Kp=',F6.2,' MLT=',F4.1)
33	FORMAT(' EKEV \ T=',F8.0,' L=',F6.2,' Kp=',F6.2,' PA=',F6.2)
34	FORMAT(' T=',F8.0,' L=',F6.2,' MLT=',F6.2)
35	FORMAT(2X,'L',4X,'Energy[keV]',3X,'# Part',3X,'Ener Loss',   &
          2X,'Part Loss')
36      FORMAT(' EKEV \ L=',F6.2,' MLT=',F4.1)
37      FORMAT(' EKEV \ L=',F6.2,' MLT=',F4.1,' Ne=',1PE11.3)
38	FORMAT(20(2X,A7,1X))
39	FORMAT(1P,20E10.3)
40	FORMAT(' PA =',20(2X,F8.2))
41	FORMAT(A10,1P,20E10.3)
42	FORMAT(F6.2,4X,1P,20E10.3)
43	FORMAT(F7.3,1P,20E9.2)
44	FORMAT(' PA =',20(1X,F8.2))
45	FORMAT(' EKEV \ T=',F8.0,' L=',F6.2,' MLT=',F4.1,' Kp=',F6.2,   &
          ' Ne=',1PE11.3)
46	FORMAT('EKEV = ',F7.3,' PA = ',F6.2)
60       FORMAT(/,4X,'Total',1P,5(2X,E10.3))
61	FORMAT(/,4X,'DEP ENERGY AT 1000 KM =',1PE11.4,' J/s')
70	FORMAT(F5.2,F10.6,1P,E13.3)
71      FORMAT(2X,'T =',F8.0,2X,'Kp =',F6.2)
72      FORMAT(2X,'L',4X,'AZIMUTH',4X,'I P FLUX')
73	FORMAT(1P,3E16.4)
74	FORMAT(A8,2A6,16A10)
80	FORMAT(F6.2,F9.3,1P,8E11.3)

      RETURN
      END
!
! End of subroutine WRESULT
!

!***********************************************************************
!			     PSRCLOSS
!     Continuous printing of energy and particle sources and losses
!	Also resets values for the next time step
!***********************************************************************
	SUBROUTINE PSRCLOSS(T)

	use ModHeidiSize
	use ModHeidiIO

	implicit none

	REAL :: T,DN,DE

	DN=RNS-RNL+ESN-ELN+ECN-ALN-CEN
	DE=RES-REL+ESE-ELE+ECE-ALE-CEE
	WRITE (15+S,50) T,LMP(7),LMP(13),RNS,RNL,RES,REL,ESN,ELN,ESE,   &
          ELE,ECN,ECE,ALN,ALE,CEN,CEE,DN,DE
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
50	FORMAT(F8.1,2F6.2,1P,16E10.2)
	RETURN
	END
!
! End of subroutine PSRCLOSS
!
