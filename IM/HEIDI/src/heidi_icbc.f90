! File name: heidi_initial.f90
!
! Contains: initial and boundary condition definition routines for HEIDI
!	INITIAL
!	LMPLOSS
!	GEOSB
!	FBC
!	FINJ
!
! Last Modified: March 2006, Mike Liemohn
!
! **********************************************************************
!				INITIAL
!   	Initial set up of distribution functions (F2), energy (ENER)
!  	                and number of particle (N)
!***********************************************************************
SUBROUTINE INITIAL(LNC,XN,J6,J18)

  use ModHeidiSize
  use ModHeidiIO
  use ModHeidiMain

  implicit none

  integer :: j6,j18,i,j,k,l,ifn,jj,ii,kk,ig1,ig2,ig3,ig4,  &
       kg,kg1,kg2,IER,Iin,Kin
  real :: weight,elat1,elat2,etotal,efractn,pg,sg,y,yz,y10,y12,x,xlt
  real :: fbc,Cst1,Cst2,GAMMLN,esum
  REAL :: XN(NR,NS),LNC(NR,NS),N,FAC
  integer ::I1,J1,K1,L1
  PARAMETER (I1=5, J1=9, K1=25, L1=19)	! From restart.bcf sizes
  !	PARAMETER (I1=5, J1=9, K1=22, L1=19)	! From restart.bcf sizes
  REAL :: FI(6,11),LI(6),EI(11),  &
       LIN(3),EIN(5),NI(3,5),E1(NE),E2(NE),F0(K1,L1),DUMMY,  &
       MU0(L1),E0(K1),MLT0(J1),R0(I1),G0(I1,J1)
  INTEGER ::IFM(38)
  CHARACTER*80 HEADER
  CHARACTER*5 ST1
  CHARACTER*2 ST2
  CHARACTER*3 ST3

  ! Common block for tests
  REAL :: XL1(4),XL2(2),XL3
  COMMON /CLIN2/ XL1,XL2,XL3

  DATA LI/2.,3.,3.5,4.5,5.5,6.5/  ! for moder storm input
  DATA EI/0,0.25,0.5,0.75,1.,1.25,1.5,1.75,2.,2.25,2.5/
  !		EI => LOG(E[keV])
  DATA LIN/2.,4.25,6.5/           ! LIN and EIN => model PA shape
  DATA EIN/2.05,21.6,50.15,153.,350./	! EIN => E[keV]
  DATA IFM/2,7,13,20,28,35,42,47,50,52,54,56,58,60,62,64,66,68,70,  &
       2,11,21,31,41,51,61,70,75,79,82,83,84,85,86,87,88,89,90/
  external :: GAMMLN

  !.......Start the loop over RC species
  DO S=1,NS

     !.......Zero out F2
     do I=1,IO
        do J=1,JO
           do K=1,KO
              do L=1,LO
                 F2(I,J,K,L,S)=1.E-25*FFACTOR(I,K,L)
              end do
           end do
        end do
     end do

     !.......Do the rest only if we're calculating this species
     IF (SCALC(S).EQ.1) THEN

        !.......Define the input file names
        IF (ISTORM.EQ.1) ST1='major'
        IF (ISTORM.EQ.2) ST1='moder'
        IF (ISTORM.EQ.3) ST1='tests'
        IF (S.EQ.1) ST2='_e'
        IF (S.EQ.2) ST2='_h'
        IF (S.EQ.3) ST2='he'
        IF (S.EQ.4) ST2='_o'

        !.......Loss cone distribution (INI=1)
        IF (INI(S).EQ.1) THEN
           IBC=1
           IF (Ab.GT.0.) THEN	! Maxwellian distribution
              do I=2,IO
                 do L=UPA(I),LO
                    do K=2,KO
                       do J=J6,J18
                          F2(I,J,K,L,S)=FBC(EKEV(K),FFACTOR(I,K,L),FINI(K)*CHI(I,J))
                       end do	! J loop
                    end do	! K loop
                 end do	! L loop
              end do	! I loop
           ELSE			! Read in from input file: 'cone.bcf'
              OPEN(unit=20,FILE='cone.bcf',STATUS='OLD')
              READ (20,101) HEADER
              do K=2,KO
                 READ (20,*) DUMMY,FINI(K)
              end do	! K loop
              CLOSE(20)
              !	  CHI0=1./CHI0	! Only do 1 divide, then multiply
              !	  CALL ZENITH(J6,J18)
              do J=J6,J18
                 do I=2,IO
110                 CHI(I,J)=1.		! Same BC flux over entire dayside
                    !110	      CHI(I,J)=CHI(I,J)*CHI0
                 end do	! I loop
              end do	! J loop
              do I=2,IO
                 do L=UPA(I),LO
                    do K=2,KO
                       do J=J6,J18
                          F2(I,J,K,L,S)=FINI(K)*CHI(I,J)*FFACTOR(I,K,L)
                       end do	! J loop
                    end do	! K loop
                 end do		! L loop
              end do		! I loop
           END IF	! end of INI=1 case

           !.......Gaussian in R and PHI about some location (INI=2)
           !.......Distribution from input files (INI=3)
        ELSE IF ((INI(S).EQ.2).OR.(INI(S).EQ.3)) THEN
           !.......Read in FI and NI from files
           ST3='.in '
           OPEN(UNIT=20,file=ST1//ST2//ST3,STATUS='OLD')
           READ(20,*) YEAR,DAY,R,AP,KP,ETOTAL,EFRACTN  !Replaces input.glo
           READ(20,101) HEADER
           READ(20,101) HEADER
           do  Kin=1,11
              READ(20,*) (FI(Iin,Kin),Iin=1,6)
           end do
           CLOSE(20)
           ST1='testn'
           OPEN(UNIT=22,file=ST1//ST2//ST3,STATUS='OLD')
           READ(22,101) HEADER
           READ(22,101) HEADER
           do KK=1,5	
              READ(22,*) (NI(II,KK),II=1,3)
           end do
           CLOSE(22)
           !.......Determine parameters for specific initial condition
           IF (INI(S).EQ.2) THEN
              ig1=4			! Gausian centered at L=LZ(ig1)
              ig2=ig1
              ig3=ig1-2
              ig4=ig1+2
              kg=24			! Gaussian centered at E=EKEV(kg)
              kg1=kg-7
              kg2=kg+5
              do k=kg1,kg2
                 E1(k)=ALOG10(EKEV(kg))
                 E2(k)=10*EXP(-(EKEV(k)-EKEV(kg))**2/15.)
              end do
              pg=3.0		! Gaussian extends to |PHI|<pg
              sg=0.05		! Variance of PHI Gaussian
           ELSE
              ig1=2
              ig2=IO
              ig3=2
              ig4=IO
              kg1=2
              kg2=KO
              do k=kg1,kg2
                 E1(k)=ALOG10(EKEV(k))
                 E2(k)=1.
              end do
              pg=5.5
              sg=18.0
           END IF
           !.......Find distribution at local midnight
           do I=ig1,ig2
              do K=kg2,kg1,-1
                 IF (EKEV(K).GT.1.) THEN
                    CALL LINTP2(LI,EI,FI,6,11,LZ(I),E1(K),Y,IER)
                    CALL LINTP2(LIN,EIN,NI,3,5,LZ(I),10**E1(K),YZ,IER)
                    Y10=(10**Y)*E2(k)
                    do L=2,UPA(I)-1
                       F2(I,1,K,L,S)=Y10*(1.-MU(L)**2)**(YZ/2.)*FFACTOR(I,K,L)
                    end do	! L loop
                    KK=K
                 ELSE			! Maxwellian below 1 keV
                    X=EKEV(k)/EKEV(KK)
                    Y12=Y10*X*EXP(1.-X)
                    do L=2,UPA(I)-1
                       F2(I,1,K,L,S)=Y12*(1.-MU(L)**2)**(YZ/2.)*FFACTOR(I,K,L)
                    end do ! L loop
                 END IF
              end do	! K loop
           end do	! I loop

           !.......Gaussian in R for INI=2
           IF (INI(S).EQ.2) THEN
              do l=2,UPA(i)
                 do k=kg1,kg2
                    do i=ig3,ig4
                       F2(i,1,k,l,s)=F2(ig1,1,k,l,s)*EXP(-(LZ(i)-LZ(ig1))**2/0.005)
                    end do
                 end do
              end do
           END IF
           !.......Gaussian in PHI
           do J=1,JO
              XLT=MLT(J)
              IF (XLT.GT.12) XLT=XLT-24.
              XLT=ABS(XLT)
              do i=ig3,ig4
                 do l=2,UPA(i)-1
                    do k=kg1,kg2
                       if (XLT.LT.pg) then
                          F2(i,j,k,l,s)=F2(i,1,k,l,s)*EXP(-XLT**2/sg)
                       end if
                    end do	! K loop
                 end do		! L loop
              end do		! I loop
           end do 		! J loop

           !.......MICS quiet RC, constant with R, PHI, and MU (INI=4)
        ELSE IF (INI(S).EQ.4) THEN
           JJ=0
           K=1
           DO WHILE (JJ.EQ.0) 
              IF (EKEV(K).GT.Eob) THEN
                 JJ=K
              ELSE IF (K.EQ.KO) THEN
                 JJ=K+1
              END IF
              K=K+1
           END DO
           II=0
           K=1
           DO WHILE (II.EQ.0) 
              IF (EKEV(K).GT.40.) THEN
                 II=K
              ELSE IF (K.EQ.KO) THEN
                 II=K+1
              END IF
              K=K+1
           END DO
           FAC=1.
           do L=1,LO
              do K=2,JJ
                 IF (IFAC.NE.1) FAC=1./FLUXFACT(S)/EKEV(K)
                 DUMMY=(ALOG(EKEV(K))-ALOG(1.))/(ALOG(Eob)-ALOG(1.))
                 WEIGHT=EXP(DUMMY*ALOG(Ab)+(1.-DUMMY)*ALOG(0.1*Ab))
                 do J=1,JO
                    do I=2,6
                       F2(I,J,K,L,S)=Ab*FAC*FFACTOR(I,K,L)*(.5)**(6-I)
                    end do		! I loop
                    do I=11,IO
                       F2(I,J,K,L,S)=Ab*FAC*FFACTOR(I,K,L)
                    end do		! I loop
                 end do		! J loop
              end do		! K Loop
              do K=JJ+1,II
                 IF (IFAC.NE.1) FAC=1./FLUXFACT(S)/EKEV(K)
                 DUMMY=(ALOG(EKEV(K))-ALOG(Eob))/(ALOG(40.)-ALOG(Eob))
                 WEIGHT=EXP(DUMMY*ALOG(0.1*Ab)+(1.-DUMMY)*ALOG(Ab))
                 do J=1,JO
                    do I=2,6
                       F2(I,J,K,L,S)=Ab*FAC*FFACTOR(I,K,L)*(.5)**(6-I)
                    end do		! I loop
                    do I=11,IO
                       F2(I,J,K,L,S)=Ab*FAC*FFACTOR(I,K,L)
                    end do		! I loop
                 end do		! J loop
              end do		! K Loop
              do K=II+1,KO
                 IF (IFAC.NE.1) FAC=1./FLUXFACT(S)/EKEV(K)
                 DUMMY=(ALOG(EKEV(K))-ALOG(40.))/(ALOG(300.)-ALOG(40.))
                 WEIGHT=EXP(DUMMY*ALOG(Ab)+(1.-DUMMY)*ALOG(0.1*Ab))
                 do J=1,JO
                    do I=2,6
                       F2(I,J,K,L,S)=WEIGHT*FAC*FFACTOR(I,K,L)*(.5)**(6-I)
                    end do		! I loop
                    do I=7,10
                       F2(I,J,K,L,S)=WEIGHT*FAC*FFACTOR(I,K,L)
                    end do		! I loop
                 end do		! J loop
                 WEIGHT=EXP(DUMMY*ALOG(0.1*Ab)+(1.-DUMMY)*ALOG(Ab))
                 do J=1,JO
                    do I=11,IO
                       F2(I,J,K,L,S)=WEIGHT*FAC*FFACTOR(I,K,L)
                    end do		! I loop
                 end do		! J loop
              end do		! K Loop
           end do		! L loop

           !.......Read in F from a file: 'restart.bcf' (INI=5)
        ELSE IF (INI(S).EQ.5) THEN
           OPEN(unit=20,FILE='restart.bcf',STATUS='OLD')
           IFN=0
           IF (IPA.EQ.0) IFN=L1
           do L=1,L1
              MU0(L)=MU(IFM(L+IFN))
           end do
           READ (20,101) HEADER
           READ (20,101) HEADER
           do I=1,I1		! Read in and perform first 2D interp.
              II=I*4-2
              do J=1,J1-1
                 JJ=J*3-2
                 READ (20,102) R0(I),MLT0(J)
                 READ (20,101) HEADER
                 do K=1,K1
                    READ (20,*) E0(K),(F0(K,L),L=1,L1)
                 end do
                 !	  IF (I+J.EQ.2) THEN
                 !	   PRINT 50,'Inputs :',1,2,2,18,E0(1),E0(2),MU0(2),MU0(18)
                 !	   PRINT 50,'My grid:',2,5,7,68,EKEV(2),EKEV(5),MU(7),MU(68)
                 !	  END IF
                 !	  PRINT 50, 'F0:',I,II,J,JJ,F0(1,2),F0(1,18),F0(2,2),F0(2,18)
50               FORMAT(A,4I4,1P,4E12.4)
                 do L=1,L1		! Convert to F2 and log scale
                    do K=1,K1
                       if (F0(K,L) .LT. 1.E-30) then
                          F0(K,L)=1.E-30
                       end if
                       F0(K,L)=ALOG10(F0(K,L))
                    end do	! K loop
                 end do	! L loop
                 do K=2,KO
                    do L=1,LO-1
                       CALL LINTP2(E0,MU0,F0,K1,L1,EKEV(K),MU(L),F2(II,JJ,K,L,S),IER)
                       IF (EKEV(K).GT.E0(K1)) &
                            F2(II,JJ,K,L,S)=AMIN1(F2(II,JJ,K,L,S),F2(II,JJ,K-1,L,S))
                    end do	! L loop
                 end do	! K loop
                 !	  PRINT 50, 'F2:',I,II,J,JJ,10**F2(II,JJ,2,7,S),10**F2(II,JJ,2,68,S),  &
                 !     		10**F2(II,JJ,5,7,S),10**F2(II,JJ,5,68,S)
                 !	  PRINT 51, 'LIN:',(XL1(K),K=1,4),(XL2(L),L=1,2),XL3
51               FORMAT (A,1P,7E10.2)
                 !	  IF (I+J.EQ.2) THEN
                 !	   PRINT 50,'Inputs :',1,2,2,18,E0(1),E0(2),MU0(2),MU0(18)
                 !	   PRINT 50,'My grid:',2,5,7,68,EKEV(2),EKEV(5),MU(7),MU(68)
                 !	   PRINT 50, 'F0#2:',I,II,J,JJ,10**F0(1,2),10**F0(1,18),10**F0(2,2),  &
                 !     		10**F0(2,18)
                 !	   PRINT 50,'F2in:',2,3,4,5,(10**F2(II,JJ,K,7,S),K=2,5)
                 !	  END IF
              end do	! J loop
           end do		! I Loop
           CLOSE(20)
           !	STOP			!TEST RUNS ONLY
           MLT0(J1)=MLT(1)
           do K=2,KO		! Perform second 2D interpolation
              do L=1,LO-1
                 do J=1,J1
                    JJ=J*3-2
                    if (J.EQ.J1) then
                       JJ=1
                    end if
                    do I=1,I1
                       II=I*4-2
                       G0(I,J)=F2(II,JJ,K,L,S)
                    end do	! I loop
                 end do	! J loop
                 do J=1,JO
                    do I=2,IO
                       CALL LINTP2(R0,MLT0,G0,I1,J1,LZ(I),MLT(J),F2(I,J,K,L,S),IER)
                       IF (LZ(I).GT.R0(I1)) &
                            F2(I,J,K,L,S)=AMIN1(F2(I,J,K,L,S),F2(I-1,J,K,L,S))
                    end do 	! I loop
                 end do 	! J loop
              end do	! L loop
           end do		! K loop
           do L=1,LO-1		! Remove negative and near zero values
              do K=2,KO		! Convert from log scale
                 do J=1,JO
                    do I=2,IO
                       F2(I,J,K,L,S)=10**(F2(I,J,K,L,S))*FFACTOR(I,K,L)
                       if (F2(I,J,K,L,S).LE.1.E-30*FFACTOR(I,K,L)) then
                          F2(I,J,K,L,S)=1.E-30*FFACTOR(I,K,L)
                       end if
                    end do	! I loop
                 end do	! J loop
              end do	! K loop
           end do	! L loop
           do K=2,KO		! Pitch angle boundary condition
              do J=1,JO
                 do I=2,IO
                    F2(I,J,K,LO,S)=F2(I,J,K,LO-1,S)*CONMU2
                 end do	! I loop
              end do	! J loop
           end do		! K loop

           !  Nightside plasmasheet injection (INI=6)
        ELSE IF (INI(S).EQ.6) THEN
           IBC=6
           IF (S.EQ.1) THEN
              Einj=.2			! Characteristic E, keV
              Kinj=6.			! Kappa value
              Ninj=.1			! Density, cm-3
           ELSE 
              Einj=1.4			! Characteristic E, keV
              Kinj=5.5			! Kappa value
              Ninj=.4			! Density, cm-3
           END IF
           Cst1=SQRT(Q*1.E4/(2.*MAS(S)*(PI*Kinj*Einj)**3))
           Cst2=EXP(GAMMLN(Kinj+1.,IER)-GAMMLN(Kinj-0.5,IER))
           do L=2,UPA(IO)-1
              do K=1,KO
                 do J=1,JO
                    !	    F2(IO,J,K,L,S)=Ninj*EKEV(K)*Cst1*Cst2*FFACTOR(IO,K,L)*   &
                    !             (1.+EKEV(K)/(Kinj*Einj))**(-Kinj-1.)
                    F2(IO,J,K,L,S)=1.E-30	! Set in GEOSB subroutine
                 end do
              end do
           end do
           do K=1,KO
              !	 PRINT 5,EKEV(K),F2(IO,1,K,2,S)/FFACTOR(IO,K,2),   &
              !           F2(IO,1,K,UPA(IO)-1,S)/FFACTOR(IO,K,UPA(IO)-1)
              do j=1,JO
                 F2(IO,J,K,1,S)=F2(IO,J,K,2,S)
              end do
           end do
5          FORMAT (1P,5E12.4)

           !  Read in from a unformatted file (INI=7)
        ELSE IF (INI(S).EQ.7) THEN
           OPEN(UNIT=1,FILE=NAME//ST2//'.unff',status='old',   &
                form='unformatted')
           DO L=1,NPA
              !	  DO K=8,NE  ! Changed the Egrid for runs "e" and "f" !1,NE
              DO K=1,NE  ! Change back to this for restarts
                 DO J=1,NT 
                    read(1) (f2(I,J,K,L,S),I=1,NR)
                    !	    f2(1:NR,J,K,L,S)=0.5*f2(1:NR,J,K,L,S)  ! Special restart line
                 END DO
              END DO
              !	  DO K=8,1,-1   ! New loop to fill in low-E grid
              !	   DO J=1,NT    ! Comment out for restarts 
              !	    F2(1:NR,J,K,L,S)=F2(1:NR,J,K+1,L,S)
              !	   END DO
              !	  END DO        ! to here
           END DO
           close(1)

        END IF
        !.......Done with initial particle distribution set up

        IF (IBC(S).EQ.1) THEN	! This is to redo INI=1 above
           Ib=1
           OPEN(unit=20,FILE='cone.bcf',STATUS='OLD')
           READ (20,101) HEADER
           do K=2,KO
              READ (20,*) DUMMY,FINI(K)
           end do
           CLOSE(20)
           do J=J6,J18
              do I=2,IO
                 CHI(I,J)=1.		! Same BC flux over entire dayside
              end do
           end do
        END IF

        !.......Calculate the total # of particles and energy of this species
        N=0				! total # dens of RC specie "s"
        ESUM=0				! total E of RC for specie "s"
        do I=1,IO
           ENER(I,S)=0			! E of RC specie for some LZ
        end do
        do I=2,IO
           do J=2,JO
              do K=2,KO
                 do L=2,UPA(I)-1
                    WEIGHT=F2(I,J,K,L,S)*WE(K)*WMU(L) ! F2 - average in cell
                    XN(I,S)=XN(I,S)+WEIGHT		 !      (i,j,k,l)
                    ENER(I,S)=ENER(I,S)+EKEV(K)*WEIGHT  ! for some LZ
                 end do	! L loop
              end do		! K loop
           end do 		! J loop
           ESUM=ESUM+ENER(I,S)
           N=N+XN(I,S)	
        end do

        !.......FACTOR, a scaling for ESUM and N, depends on IFAC:
        IF (IFAC.EQ.1) FACTOR(S)=8.6474E13/M1(S)**1.5*DR*DPHI
        IF (IFAC.EQ.2) FACTOR(S)=4.3237E13/M1(S)**1.5*DR*DPHI

        !.......Calculate the characteristics (mean energy) at T=0
        ST1='rnsc1'
        IF(T/3600.EQ.48) ST1='rnsc3'
        IF(T/3600.EQ.96) ST1='wpa10'
        !	OPEN(1,FILE=ST1//ST2//'.l')
        !	WRITE(1,*)' Losses due to some processes'
        !	WRITE(1,15) KP
15      FORMAT(2X,5HT = 0,2X,4HKp =,F6.2,/,6X,1HL,10X,11HENIGHT[keV],2X,   &
             9HEDAY[keV],2X,6HNNIGHT,2X,4HNDAY,4X,6HMean E)
        !	DO I=2,IO
        !	 EMEAN=ENER(I,S)/XN(I,S)
        !	 WRITE(1,17) LZ(I),ENER(I,S)*FACTOR(S),XN(I,S)*FACTOR(S),EMEAN
        !	END DO
17      FORMAT(2X,F7.2,5(2X,1PE12.5))
        !	AMEAN=ESUM/N
        !	WRITE(1,19)ESUM*FACTOR(S),N*FACTOR(S),AMEAN
19      FORMAT(/,4X,5HTotal,3(2X,1PE12.5))
        !	CLOSE(1)

        !.......Initial loss is zero
        LNC(1:io,S)=0       ! Loss of particles due to all processes
        LEC(1:io,S)=0       ! Loss of energy due to     "       "

     END IF		! SCALC Check
  END DO		! S loop

  !.......Calculate the energy coeff
  X=RE+HMIN
  do I=1,IO
     ELAT1=SQRT(X/(LZ(I)-DL1/2.)/RE)
     ELAT2=SQRT(X/(LZ(I)+DL1/2.)/RE)
     ECOF(I)=DPHI*X**2*(ACOS(ELAT2)-ACOS(ELAT1))*(ELAT1+ELAT2)/2.
  end do	! I loop

101 FORMAT(A80)
102 FORMAT(21X,F6.2,5X,F4.1)
103 FORMAT(F7.3,20(1PE9.2))

  RETURN
END SUBROUTINE INITIAL
!
! End of subroutine INITIAL
!

! **********************************************************************
!				 LMPLOSS
!   		If LMP moved inward, then lose the particles 
!  	                
!***********************************************************************
SUBROUTINE LMPLOSS

  use ModHeidiSize
  use ModHeidiIO
  use ModHeidiMain

  IMPLICIT NONE

  REAL :: FLO
  INTEGER :: I,K,L,J

  DO J=1,JO
     IF (ILMP(J).LT.ILold(J)) THEN
        DO L=1,LO		! Lose everything beyond magnetopause
	   DO K=2,KO
              DO I=ILMP(J)+1,ILold(J)
                 FLO=1.E-30*FFACTOR(I,K,L)
                 RNL=RNL+(F2(I,J,K,L,S)-FLO)*CONSL(K,S)*WE(K)*WMU(L)*DR*DPHI
                 REL=REL+(F2(I,J,K,L,S)-FLO)*CONSL(K,S)*EKEV(K)   &
                      *WE(K)*WMU(L)*DR*DPHI
                 F2(I,J,K,L,S)=FLO
              END DO
	   END DO
        END DO
     END IF
  END DO

  RETURN
END SUBROUTINE LMPLOSS
!
! End of subroutine LMPLOSS
!

! **********************************************************************
!				 GEOSB
!   			Boundary conditions set up 
!  	                
!***********************************************************************
SUBROUTINE GEOSB

  use ModHeidiSize
  use ModHeidiIO
  use ModHeidiMain
  use ModHeidiDrifts

  IMPLICIT NONE

  CHARACTER*80 HEADER
  REAL :: Fkapb(NE),Foutb,X,DATA(9),fac
  REAL :: TM1,TM2,NM1,NM2,TFM1,TFM2,TCM1,TCM2,NM,TFM,TCM,TS1,TS2,  &
       Flanl(NE,NPA,NS),FS1(7),FS2(7),ES(7),FS(7),NY(NS),NEL,  &
       NE1,NE2,TEF1,TEF2,TEC1,TEC2,NMFAC(NS),FSFAC(NS),TEF,TEC,  &
       Ekap,Kappa,GAMMLN
  integer :: I,J,K,L,IER,GetBCI,BCI,I2,KES(NS),I6,I7,I9,IG7
  integer :: SKAPPA(4)
  external :: GetBCI,GAMMLN
  DATA ES/45.,60.,94.,141.,210.,325.,535./, Kappa/5.00/
  DATA SKAPPA/1,1,1,1/
  SAVE KES,TM1,TM2,NM1,NM2,TFM1,TFM2,TCM1,TCM2,TS1,TS2,FS1,FS2,  &
       I2,I6,I7,I9,IG7,NE1,NE2,TEF1,TEF2,TEC1,TEC2

  print *, 'Resetting the outer boundary condition'
  !CCC Create a few flags and open a few files
  IF (T.EQ.TIME) THEN
     I6=0
     I7=0
     IG7=0
     I9=0
     DO S=1,NS
        IF (IBC(S).EQ.6) I6=1
        IF (IBC(S).EQ.7) I7=1
        IF (IBC(S).EQ.9) I9=1
        IF (IBC(S).GT.7) IG7=1
     END DO
     IF (I7.EQ.1 .OR. IG7.EQ.1) THEN
        DO S=1,NS
           IF (S.EQ.1 .OR. IBC(S).EQ.7) THEN    ! no SOPA data
              KES(S)=KO
           ELSE
              KES(S)=0
              K=0
              DO WHILE (KES(S).EQ.0)
                 K=K+1
                 IF (EKEV(K).GE.ES(1)) THEN
                    KES(S)=K-1
                 ELSE
                    IF (K.EQ.KO) KES(S)=KO
                 END IF
              END DO
           END IF
        END DO
        IF (IG7.EQ.1) THEN
           TS2=TIME-1.		! Prepare SOPA input file
           TS1=TS2
           FS2(1:7)=0.
           OPEN(UNIT=14,FILE=NAME//'_sopa.in',status='old')
           DO I=1,3
              READ(14,*) HEADER
           END DO
        END IF
        TM2=TIME-1.		! Prepare MPA input file
        TM1=TM2
        NM2=0.
        TFM2=0.
        TCM2=0.
        NE2=0.
        TEC2=0.
        TEF2=0.
        OPEN(UNIT=16,FILE=NAME//'_mpa.in',status='old')
        DO I=1,3			! 3 lines of header material
           READ(16,*) HEADER
        END DO
        I2=0
        IF (S.EQ.1) I2=3
     END IF
  END IF

  DO S=1,NS

     DO L=1,LO
        DO K=1,KO
           DO J=1,JO
              FGEOS(J,K,L,S)=0. 
           enddo
        enddo
     enddo
     IF (SCALC(S).EQ.1) THEN

        IF (TINJ.GT.TIME+2.*DT*NSTEP) THEN ! No injection, use IC for BC
           DO L=1,LO
              DO K=1,KO
                 DO J=1,JO
                    FGEOS(J,K,L,S)=F2(IO,J,K,L,S)
                 END DO
              END DO
           END DO
        END IF

     END IF   ! SCALC check
  END DO	 ! S loop

  IF (I7.EQ.1 .OR. IG7.EQ.1) THEN	! LANL data injection
     IF (TM2.LT.T) THEN			! MPA DATA
        DO WHILE (TM2.LE.T)		! Best if final TM2 > final T
           TM1=TM2
           NM1=NM2
           TFM1=TFM2
           TCM1=TCM2
           READ (16,*,IOSTAT=L) (DATA(I),I=1,9)
           TM2=DATA(2)
           NM2=DATA(4)
           TFM2=DATA(6)
           TCM2=DATA(5)
           NE2=DATA(7)
           TEF2=DATA(9)
           TEC2=DATA(8)
           IF (L.LT.0) TM2=TIME+2*DT*(NSTEP+1)
           IF (T.EQ.TIME) THEN		! In case T2>T already
              TM1=TIME
              NM1=NM2
              TFM1=TFM2
              TCM1=TCM2
              NE1=NE2
              TEC1=TEC2
              TEF1=TEF2
           END IF
        END DO
     END IF
     FAC=(T-TM1)/(TM2-TM1)			! Linearly interpolate
     NM=FAC*NM2+(1.-FAC)*NM1		! in cm-3
     TFM=(FAC*TFM2+(1.-FAC)*TFM1)*1.E-3	! in keV
     TCM=(FAC*TCM2+(1.-FAC)*TCM1)*1.E-3	! in keV
     NEL=FAC*NE2+(1.-FAC)*NE1		! in cm-3
     TEF=(FAC*TEF2+(1.-FAC)*TEF1)*1.E-3	! in keV
     TEC=(FAC*TEC2+(1.-FAC)*TEC1)*1.E-3	! in keV
     NY(2)=0.34*EXP(0.054*KP)	! From Young et al 1982
     NY(3)=5.1E-3*EXP(6.6E-3*F107)
     NY(4)=0.011*EXP(0.24*KP+0.011*F107)
     FAC=0.
     DO I=2,4 			! Weights corrects MPA moments
        FAC=FAC+NY(I)/SQRT(M1(I))
     END DO
     NMFAC(1)=1.
     DO S=2,NS
        NMFAC(S)=NY(S)/FAC
        !CCC The next line is for the idealized test runs only, no BC variation:
        !CCC	    NMFAC(S)=1.
     END DO
     IF (IG7.EQ.1) THEN		! SOPA DATA
        IF (TS2.LT.T) THEN
           DO WHILE (TS2.LE.T)	! Best if final TS2 > final T
              TS1=TS2
              FS1(2:7)=FS2(2:7)
              READ (14,*,IOSTAT=I) (DATA(I),I=1,8)
              TS2=DATA(1)
              FS2(2:7)=DATA(3:8)
              IF (I.LT.0) TS2=TIME+2*DT*(NSTEP+1)
              IF (T.EQ.TIME) THEN		! In case T2>T already
                 TS1=TIME
                 FS1(2:7)=FS2(2:7)
              END IF
           END DO
        END IF
        FAC=(T-TS1)/(TS2-TS1)		! Linearly interpolate
        FS(2:7)=FAC*FS2(2:7)+(1.-FAC)*FS1(2:7)	! flux units
     END IF
     FAC=(NY(2)+NY(3)+NY(4))  
     FSFAC(1)=1.
     DO S=2,NS
        FSFAC(S)=NY(S)/FAC
        !CCC The next line sometimes gives bad results...use with care.
        IF (NEL*FSFAC(S).GT.NM*NMFAC(S)) NMFAC(S)=NEL*FSFAC(S)/NM
     END DO
     DO S=1,NS
        IF (SCALC(S).EQ.1) THEN
           !CCC _hsopa: H+ has a different function than all other species
           IF (SKAPPA(S).EQ.1) THEN  ! bi-kappa=5,  SOPA (if designated)
              Ekap=TFM*(Kappa-1.5)/Kappa
              FAC=(TFM/TCM)*(MAS(S)*1.E13/(2.*PI*Kappa*Ekap*Q))**1.5
              FAC=FAC*EXP(GAMMLN(Kappa+1.,IER)-GAMMLN(Kappa-0.5,IER))
              DO L=1,LO
                 do K=1,KES(S)          ! MPA moments
                    Flanl(K,L,S)=NM*NMFAC(S)*FAC*(1.+EKEV(K)*(MU(L)**2+   &
                         (TFM/TCM)*(1.-MU(L)**2))/(Kappa*Ekap))**(-Kappa-1.)
                 end do
              END DO
              FS(1)=AMAX1(.5*FS(2),Flanl(KES(S),10,S)   &
                   *FLUXFACT(S)*EKEV(KES(S)))
              DO K=KES(S)+1,KO                        ! SOPA fluxes
                 CALL LINTP(ES,FS,7,EKEV(K),FAC,IER)
                 Flanl(K,1:LO,S)=FAC/FLUXFACT(S)/EKEV(K)
              END DO
           ELSE  ! Maxwellian everywhere, SOPA (if designated)
              DO L=1,LO
                 DO K=2,KES(S)			! MPA moments
                    Flanl(K,L,S)=NM*NMFAC(S)*(MAS(S)*1.E13/Q/2./PI)**1.5/  &
                         SQRT(TFM)/TCM*EXP(-EKEV(K)*((1.-MU(L)*MU(L))/TCM    &
                         +MU(L)*MU(L)/TFM))
                 END DO  ! K loop for MPA
              END DO   ! L loop for MPA
              FS(1)=AMAX1(.5*FS(2)*FSFAC(S),Flanl(KES(S),10,S)   &
                   *FLUXFACT(S)*EKEV(KES(S)))/FSFAC(S)
              DO K=KES(S)+1,KO			! SOPA fluxes
                 CALL LINTP(ES,FS,7,EKEV(K),FAC,IER)
                 Flanl(K,1:LO,S)=FSFAC(S)*FAC/FLUXFACT(S)/EKEV(K)
              END DO ! K loop for SOPA
           END IF  ! Block for S=2 and all others
        END IF  ! SCALC check
     END DO   ! S loop
  ELSE 
     RETURN
  END IF

  !...Find injection boundary and fill in (Note: IBC=9, inject everywhere)
  DO S=1,NS
     IF (SCALC(S).EQ.1) THEN
        IF (IBC(S).EQ.6) THEN
           CALL FINJ(Fkapb)
           DO L=1,LO
              DO K=2,KO
                 DO J=1,JO
                    FGEOS(J,K,L,S)=Fkapb(K)*FFACTOR(IO,K,L)
                 END DO	! J LOOP
              END DO	! K LOOP
           END DO	! L LOOP
        ELSE IF (IBC(s).GE.7) THEN
           DO L=1,LO
              DO K=2,KO
                 DO J=1,JO
                    FGEOS(J,K,L,S)=Flanl(K,L,S)*FFACTOR(IO,K,L)
                 END DO	! J LOOP
              END DO	! K LOOP
           END DO	! L LOOP
        END IF
     END IF		! SCALC check
  END DO		! S loop

  RETURN
END SUBROUTINE GEOSB

!
! End of subroutine GEOSB
!

!***********************************************************************
!				FBC
!	This function calculates the dayside loss cone boundary flux
!	It uses either a Maxwellian or a distribution from a file
!!***********************************************************************
REAL FUNCTION FBC(E,FAC,F1)

  use ModHeidiSize
  use ModHeidiIO

  implicit none

  real :: e,fac,f1

  IF (Ab.GT.0.) THEN
     FBC=Ab*EXP(-E/Eob)*FAC
  ELSE
     FBC=F1*FAC
  END IF
  RETURN
END FUNCTION FBC
!
! End of subroutine FBC
!

!***********************************************************************
!			     FINJ
!     Calculates the plasma sheet distribution for the boundary
!***********************************************************************
SUBROUTINE FINJ(F)

  use ModHeidiSize
  use ModHeidiIO
  use ModHeidiMain

  IMPLICIT NONE

  INTEGER :: K,IER,I
  REAL :: F(NE),Cst1,Cst2,GAMMLN,CONV
  real :: T1,T2,BZ1,BZ2,MD1,MD2,U1,U2,FAC,ERF,NY(NS),TLAG
  character header*80
  SAVE T1,T2,BZ1,BZ2,MD1,MD2,U1,U2
  external :: GAMMLN,ERF

  TLAG=4.*3600.			! From Borovsky et al, Aug 98
  IF (ISWB.EQ.1) THEN
     IF (T.EQ.TIME) THEN
        T2=TIME-1.
        T1=T2
        OPEN(UNIT=15,FILE=NAME//'_sw2.in',status='old')
        DO I=1,6                      ! 6 lines of header material
           READ(15,*) HEADER
        END DO
     END IF
     IF (T2.LT.T) THEN
        DO WHILE (T2.LT.T)    ! Best if final T2 > final T
           T1=T2
           BZ1=BZ2
           MD1=MD2
           U1=U2
           READ (15,*,IOSTAT=I) T2,BZ2,MD2,U2
           T2=T2+TLAG
           IF (I.LT.0) T2=TIME+2*DT*(NSTEP+1)+TLAG
           IF (T.EQ.TIME) THEN                 ! In case T2>T already
              T1=TIME
              BZ1=BZ2
              MD1=MD2
              U1=U2
           END IF
        END DO
     END IF
     FAC=(T-T1)/(T2-T1)                      ! Linearly interpolate
     NSWB=(FAC*MD2+(1.-FAC)*MD1)	        ! in cm-3
     USWB=(FAC*U2+(1.-FAC)*U1)                ! in km/s
     Einj=2.17+0.0223*USWB		! From Borovsky et al 1998
     Ninj=0.292*NSWB**0.49
  ELSE
     NSWB=0.
     USWB=0.
     Ninj=.6			! Average ion values
     Einj=10.
  END IF
  IF (S.EQ.1) THEN
     Einj=Einj/7.			! From Huang and Frank 1986
     Kinj=6.
     Ninj=Ninj*.5			! From various sources
  ELSE
     !	  Kinj=5.*SQRT(Einj)		! From Huang and Frank 1986
     Kinj=8.-5.*ERF(Einj/23.,IER)	! Gets harder with Einj
     NY(2)=0.34*EXP(0.054*KP)	! From Young et al 1982
     NY(3)=5.1E-3*EXP(6.6E-3*F107)
     NY(4)=0.011*EXP(0.24*KP+0.011*F107)
     NY(1)=0.
     DO I=2,4 			! Weights corrects LANL data
        NY(1)=NY(1)+NY(I)/SQRT(M1(I))
     END DO
     Ninj=Ninj*NY(S)/NY(1)
  END IF
  Cst1=(MAS(S)*1.E13/(2.*PI*Kinj*Einj*Q))**1.5
  Cst2=EXP(GAMMLN(Kinj+1.,IER)-GAMMLN(Kinj-0.5,IER))
  CONV=1.
  do K=1,KO
     IF (IFAC.EQ.1) CONV=FLUXFACT(S)*EKEV(K)
     F(K)=Ninj*CONV*Cst1*Cst2*(1.+EKEV(K)/(Kinj*Einj))**(-Kinj-1.)
  end do
  RETURN
END SUBROUTINE FINJ
!
! End of subroutine FINJ
!


