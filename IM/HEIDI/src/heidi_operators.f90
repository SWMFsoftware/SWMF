! File name: heidi_operators_010.f90
!
! Contains: phase space density calculation routines for HEIDI
!	LLIMIT
!	FCHECK
!	DRIFTR
!	DRIFTP
!	DRECOUL
!	DRIFTMU
!	CHAREXCHANGE
!	COULMU
!
! Last Modified: March 2006, Mike Liemohn
!
!***********************************************************************
!				LLIMIT
!	Sets very small  values to zero to avoid computer zero
!***********************************************************************
SUBROUTINE LLIMIT

  use ModHeidiSize
  use ModHeidiMain
  
  implicit none
  
  integer :: i,j,k,l
  
  !	do s=1,ns
  !	 if (SCALC(S).EQ.1) THEN
  do L=1,LO
     do K=1,KO
        do J=1,JO
           do I=1,IO
              if (ABS(F2(I,J,K,L,S)).LE.1.E-29*FFACTOR(I,K,L))   &
                   F2(I,J,K,L,S)=1.E-30*FFACTOR(I,K,L)
           end do
        end do
     end do
  end do
  !	 end if
  !	end do

  RETURN
END SUBROUTINE LLIMIT

!***********************************************************************
!			   	FCHECK
!	Checks the distribution for negative or infinite results
!***********************************************************************
SUBROUTINE FCHECK(MARK)
  
  use ModHeidiSize
  use ModHeidiMain
  use ModHeidiIo, ONLY: write_prefix, iUnitStdOut
  implicit none
  
  integer :: mark
  
  integer :: Ibad,i,j,k,l,i1,j1,k1,l1,Ibadt,m
  
  Ibadt=0
  DO m=1,NS
     IF (SCALC(m).EQ.1) THEN
        
	Ibad=0
	do L=1,Lo
           do k=2,ko
              do j=1,jo
                 do i=2,io
                    IF ((F2(I,J,K,L,m).LT.-1.E-12).OR.   &   ! -1E-29
                         (F2(I,J,K,L,m)-F2(I,J,K,L,m).NE.0.)) THEN
                       call write_prefix; write(iUnitStdOut,*) 'Bad F,T,I,J,K,L,MARK:',F2(I,J,K,L,m),   &
                            T,m,I,J,K,L,MARK,Ibad
                       call write_prefix; write(iUnitStdOut,*) 'Radius:',(F2(I1,J,K,L,m),I1=2,IO,IO/10)
                       call write_prefix; write(iUnitStdOut,*) 'Azimuth:',(F2(I,J1,K,L,m),J1=1,JO,JO/10)
                       call write_prefix; write(iUnitStdOut,*) 'Energy:',(F2(I,J,K1,L,m),K1=2,KO,KO/10)
                       call write_prefix; write(iUnitStdOut,*) 'Pitch angle:',(F2(I,J,K,L1,m),L1=1,LO,LO/10)
                       Ibad=Ibad+1
                       IF (Ibad.EQ.20) GOTO 151 
                    END IF
                 end do
              end do
           end do
	end do
151 	CONTINUE
	Ibadt=Ibadt+Ibad
     end if
  end do
  
  IF (Ibadt.GT.0) THEN
     call write_prefix; write(iUnitStdOut,*) 'Ibad=',Ibadt
     call write_prefix; write(iUnitStdOut,*) 'A=',A
     CALL ECFL
     STOP
  END IF
  
99 FORMAT(A,1PE11.3,0PF8.0,7(2X,I4))
98 FORMAT(A,200(1PE11.3))
  RETURN
END SUBROUTINE FCHECK

!***********************************************************************
!			     DRIFTR
!     Routine calculate the change of distribution function due to
! 			    radial drift
!***********************************************************************
SUBROUTINE DRIFTR

  use ModHeidiSize
  use ModHeidiIO
  use ModHeidiMain
  use ModHeidiDrifts

  implicit none

  REAL :: f01(nt,ne,npa,ns),f02(nt,ne,npa,ns)
  REAL ::F(NR+2),FBND(NR),C(NR),LIMITER
  INTEGER :: UR,isign
  integer ::i,j,k,l,n
  real :: fup,RR,corr,x
  real :: fbc
  save f01,f02

  !...Set up injection boundary fluxes
  IF (AMOD(T,TINJ).LT.2.*DT .OR. T.EQ.TIME) THEN
     DO L=1,LO
        DO K=1,KO
           DO J=1,JO
              f01(j,k,l,S)=FGEOS(J,K,L,S)*CONF1
              f02(j,k,l,S)=FGEOS(J,K,L,S)*CONF2
           end do
        end do
     end do
  END IF				! End BC update

  do J=1,JO
     do K=2,KO		!	< FBND(I) >
        do L=2,LO		!  |____.____|____.____|____.____|
	   do I=1,IO		!  <   F(I)  >
              F(I) = F2(I,J,K,L,S)  ! F - average in cell(i,j,k,l)
	   end do
	   IF (VR(IO,J).GE.0.) THEN
              FBND(1)=0.	          ! outflow b. c.
              FBND(IO)=F(IO)       ! upwind for the side with no b.c.
              C(IO)=VR(IO,J)
              UR=IO-1
	   ELSE
              FBND(1)=F(2)         ! upwind for the side with no b.c
              C(1)=VR(1,J)
              UR=IO
              F(IO+1)=f01(j,k,l,s)
              F(IO+2)=f02(j,k,l,s)
	   END IF
           !	   C(1)=AMIN1(0.99,AMAX1(-0.99,C(1)))
	   do I=2,UR
              C(I)=VR(I,J)
              !	     C(I)=AMIN1(0.99,AMAX1(-0.99,C(I)))
              X=F(I+1)-F(I)
              ISIGN=1
              IF(C(I).NE.ABS(C(I))) ISIGN=-1
              FUP=0.5*(F(I)+F(I+1)-ISIGN*X) 	! upwind
              IF (ABS(X).LE.1.E-27) FBND(I)=FUP
              IF (ABS(X).GT.1.E-27) THEN
                 N=I+1-ISIGN
                 RR=(F(N)-F(N-1))/X
                 IF (RR.LE.0) FBND(I)=FUP
                 IF (RR.GT.0) THEN
                    LIMITER=AMAX1(AMIN1(2.*RR,1.),AMIN1(RR,2.))
                    CORR=-0.5*(C(I)-ISIGN)*X
                    FBND(I)=FUP+LIMITER*CORR	! at boundary of cell
                 END IF
              END IF
              !	 IF (I.EQ.9 .AND. J.EQ.3 .AND. K.EQ.2 .AND. L.EQ.2)    &
              !         PRINT 50,I,J,ISIGN,N,FBND(I),FUP,RR,CORR,LIMITER,X,F(I),F(I+1),F(N),F(N-1)
	   end do ! I loop

           !........update the solution for next time step
	   !C(1) = -1.5951970E-05
           do I=2,ILMP(J)
              F2(I,J,K,L,S)=F2(I,J,K,L,S)-C(I)*FBND(I)+C(I-1)*FBND(I-1)
	   end do ! I loop 
	   do I=ILMP(J)+1,IO
              F2(I,J,K,L,S)=1.E-30*FFACTOR(I,K,L)
	   end do
50         FORMAT(4I3,1P,10E10.2)

           !.......Count sources and losses out the spatial boundaries
	   I=ILMP(J)			! Choose I_magnetopause
	   IF (C(I).LT.0.) THEN	! Gain at outer boundary
              RNS=RNS-C(I)*FBND(I)*CONSL(K,S)*WE(K)*WMU(L)*DR*DPHI
              RES=RES-C(I)*FBND(I)*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
	   ELSE				! Loss at outer boundary
              RNL=RNL+C(I)*FBND(I)*CONSL(K,S)*WE(K)*WMU(L)*DR*DPHI
              REL=REL+C(I)*FBND(I)*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
	   END IF
	   IF (C(1).GT.0.) THEN		! Gain at inner boundary
              RNS=RNS+C(1)*FBND(1)*CONSL(K,S)*WE(K)*WMU(L)*DR*DPHI
              RES=RES+C(1)*FBND(1)*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
	   ELSE				! Loss at inner boundary
              RNL=RNL-C(1)*FBND(1)*CONSL(K,S)*WE(K)*WMU(L)*DR*DPHI
              REL=REL-C(1)*FBND(1)*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
	   END IF

        end do	! End L loop

! dayside BC for IBC=1
        IF (Ib.EQ.1 .AND. S.EQ.1 .AND. J.GE.J6 .AND. J.LE.J18) THEN
           do L=UPA(I),LO
              do I=2,IO
                 F2(I,J,K,L,S)=FBC(EKEV(K),FFACTOR(I,K,L),FINI(K)*CHI(I,J))
              end do
           end do
        END IF
     end do	! K loop
  end do	! J loop
  RETURN
END SUBROUTINE DRIFTR

!***********************************************************************
!			     DRIFTP
!     Routine calculate the change of distribution function due to
! 			 azimuthal drift
!***********************************************************************
	SUBROUTINE DRIFTP

	use ModHeidiSize
	use ModHeidiIO
	use ModHeidiMain
	use ModHeidiDrifts

	implicit none

	REAL :: FBND(0:NT),F(0:NT+2),C(0:NT),LIMITER
	integer :: i,j,k,l,n,isign,imag
	real :: fup,RR,corr,x
	real :: fbc

  	do I=2,IO
	 do K=2,KO
	  do L=2,LO
	   do J=1,JO
	     F(J)=F2(I,J,K,L,S)
	   end do
	   do J=1,2
	     F(JO+J)=F2(I,J,K,L,S)
	   end do
	   F(0)=F2(I,JO,K,L,S)

	   do J=1,JO
	    C(J)=P1(I,J)+P2(I,K,L)
!	    C(J)=AMIN1(0.99,AMAX1(-0.99,C(J)))
	    ISIGN=1
	    IF(C(J).NE.ABS(C(J))) ISIGN=-1
	    X=F(J+1)-F(J)
	    FUP=0.5*(F(J)+F(J+1)-ISIGN*X)
	    IF (ABS(X).LE.1.E-27) FBND(J)=FUP
	    IF (ABS(X).GT.1.E-27) THEN
	      N=J+1-ISIGN
	      RR=(F(N)-F(N-1))/X
	      IF (RR.LE.0) FBND(J)=FUP
	      IF (RR.GT.0) THEN
	        LIMITER=AMAX1(AMIN1(2.*RR,1.),AMIN1(RR,2.))
	        CORR=-0.5*(C(J)-ISIGN)*X
	        FBND(J)=FUP+LIMITER*CORR
              END IF
	    END IF
	   end do ! End J loop
	   C(0)=C(JO)
	   FBND(0)=FBND(JO)

	   do J=1,JO
	     F2(I,J,K,L,S)=F2(I,J,K,L,S)-C(J)*FBND(J)+C(J-1)*FBND(J-1)
	   end do
	  end do	! End L loop

	  IF (ISW.GT.0) THEN		! SW can compress magnetopause
	  imag=1
	  do J=1,JO
	    IF (ILMP(J+1).LT.I .AND. imag.EQ.1) THEN ! Leaving m'sphere
	     imag=0
	     IF (C(J).LT.0.) THEN	! Gain at outer boundary
	      RNS=RNS-C(J)*FBND(J)*CONSL(K,S)*WE(K)*WMU(L)*DR*DPHI
	      RES=RES-C(J)*FBND(J)*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
	     ELSE			! Loss at outer boundary
	      RNL=RNL+C(J)*FBND(J)*CONSL(K,S)*WE(K)*WMU(L)*DR*DPHI
	      REL=REL+C(J)*FBND(J)*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
	     END IF
	    ELSE IF (ILMP(J+1).GE.I .AND. imag.EQ.0) THEN ! Reentering
	     imag=1
	     IF (C(J).GT.0.) THEN	! Gain at outer boundary
	      RNS=RNS+C(J)*FBND(J)*CONSL(K,S)*WE(K)*WMU(L)*DR*DPHI
	      RES=RES+C(J)*FBND(J)*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
	     ELSE			! Loss at outer boundary
	      RNL=RNL-C(J)*FBND(J)*CONSL(K,S)*WE(K)*WMU(L)*DR*DPHI
	      REL=REL-C(J)*FBND(J)*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
	     END IF
	    END IF
	  end do
	  END IF

	  IF (Ib.EQ.1 .AND. S.EQ.1) THEN
	  do J=J6,J18
	    do L=UPA(I),LO
	     F2(I,J,K,L,S)=FBC(EKEV(K),FFACTOR(I,K,L),FINI(K)*CHI(I,J))
	    end do
	  end do	! J loop
	  END IF

	 end do	! big K loop
	end do	! big I loop

	RETURN
	END
! 
! End of subroutine DRIFTP
!

!***********************************************************************
!			DRECOUL
!	Routine calculates the change of distribution function due to
!	energization along the drift path and Coulomb energy decay
!***********************************************************************
	SUBROUTINE DRECOUL

	use ModHeidiSize
	use ModHeidiIO
	use ModHeidiMain
	use ModHeidiDrifts

	implicit none

	REAL :: FBND(NE),F(0:NE+2),C(NE),LIMITER,corr,fup,RR,x,CD
	integer :: i,j,k,l,n,isign
	real :: fbc

	do J=1,JO
	 do I=2,ILMP(J)
	  do L=2,LO
	   do k=0,KO+2
	    F(K)=0.
	   end do
	   do k=1,KO
	    C(k)=0.
	   end do 
!	   do K=1,KO
!		F(K)=F2(I,J,K,L,S)
		f(1:ko)=f2(i,j,1:ko,l,s)
!	   end do
	   do K=1,KO
	    C(K)=EDOT(I,J,K,L)*VR(I,J)+(COULE(I,K,L,S)+COULI(I,K,L,S))*XNE(I,J)
!	    C(K)=AMIN1(0.99,AMAX1(-0.99,C(K)))
	    ISIGN=1
	    IF (C(K).NE.ABS(C(K))) ISIGN=-1
	    X=F(K+1)-F(K)
	    FUP=0.5*(F(K)+F(K+1)-ISIGN*X)
	    IF (ABS(X).LE.1.E-27) FBND(K)=FUP
	    IF (ABS(X).GT.1.E-27) THEN
		N=K+1-ISIGN
 		RR=(F(N)-F(N-1))/X
		IF (RR.LE.0) FBND(K)=FUP
		IF (RR.GT.0) THEN
		  LIMITER=AMAX1(AMIN1(2.*RR,1.),AMIN1(RR,2.))
		  CORR=-0.5*(C(K)-ISIGN)*X
		  FBND(K)=FUP+LIMITER*CORR
		END IF
	    END IF
	   end do	! K loop

	   do K=2,KO
		F2(I,J,K,L,S)=F2(I,J,K,L,S)-C(K)*FBND(K)*DE(K)/WE(K)   &
     		  +C(K-1)*FBND(K-1)*DE(K-1)/WE(K)
	   end do

!.......Keep track of sources and losses, separately for drift and CC
!	Actually a net values for each, adding the sources and
!	subtracting the losses
!	Note: particle changes only at the Erange boundaries, while
!	energy changes throughout the range; E endpoints done in loop
	   CD=EDOT(I,J,1,L)*VR(I,J)*DE(1)/WE(2)	! Drift at K=1,2 bnd
	   IF (CD.GT.0) THEN
	     ESN=ESN+CD*FBND(1)*CONSL(2,S)*WE(2)*WMU(L)*DR*DPHI
	   ELSE
	     ELN=ELN-CD*FBND(1)*CONSL(2,S)*WE(2)*WMU(L)*DR*DPHI
	   END IF
	   CD=C(1)*DE(1)/WE(2)-CD		! CC at K=1,2 bnd
	   ECN=ECN+CD*FBND(1)*CONSL(2,S)*WE(2)*WMU(L)*DR*DPHI
	   CD=EDOT(I,J,KO,L)*VR(I,J)*DE(KO)/WE(KO)	! Drift at K=KO
	   IF (CD.LE.0) THEN
	     ESN=ESN-CD*FBND(KO)*CONSL(KO,S)*WE(KO)*WMU(L)*DR*DPHI
	   ELSE 
	     ELN=ELN+CD*FBND(KO)*CONSL(KO,S)*WE(KO)*WMU(L)*DR*DPHI
	   END IF
	   CD=C(KO)*DE(KO)/WE(KO)-CD		! CC at K=KO
	   ECN=ECN-CD*FBND(KO)*CONSL(KO,S)*WE(KO)*WMU(L)*DR*DPHI
	   DO K=2,KO				! Now do energy changes
	     CD=EDOT(I,J,K-1,L)*VR(I,J)*DE(K-1)/WE(K)	! Drift at lower bnd
	     IF (CD.GT.0) THEN
	       ESE=ESE+CD*FBND(K-1)*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
	     ELSE
	       ELE=ELE-CD*FBND(K-1)*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
	     END IF
	     CD=C(K-1)*DE(K-1)/WE(K)-CD		! CC at lower bnd
	     ECE=ECE+CD*FBND(K-1)*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
	     CD=EDOT(I,J,K,L)*VR(I,J)*DE(K)/WE(K)	! Drift at upper bnd
	     IF (CD.LE.0) THEN
	       ESE=ESE-CD*FBND(K)*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
	     ELSE
	       ELE=ELE+CD*FBND(K)*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
	     END IF
	     CD=C(K)*DE(K)/WE(K)-CD		! CC at upper bnd 
	     ECE=ECE-CD*FBND(K)*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
	   END DO

	  end do	! End L loop

! dayside BC for IBC=1
	  IF (Ib.EQ.1 .AND. S.EQ.1 .AND. J.GE.J6 .AND. J.LE.J18) THEN
	    do L=UPA(I),LO
	     do K=2,KO
	      F2(I,J,K,L,S)=FBC(EKEV(K),FFACTOR(I,K,L),FINI(K)*CHI(I,J))
	     end do
	    end do
	  END IF

	 end do	! I loop
	end do	! J loop

	RETURN
	END
!
! End of subroutine DRECOUL
!

!***********************************************************************
!			DRIFTMU
!	Routine calculates the change of distribution function due to
!		pitch angle changes along the drift path
!***********************************************************************
	SUBROUTINE DRIFTMU

	use ModHeidiSize
	use ModHeidiIO
	use ModHeidiMain
	use ModHeidiDrifts

	implicit none

	REAL :: FBND(NPA),F(NPA),C(NPA),LIMITER,corr,fup,RR,x
	INTEGER :: UL,ULL
	integer :: i,j,k,l,n,isign
	real :: fbc

	do J=1,JO
	 do I=2,ILMP(J)
	  do K=2,KO
	   do L=2,LO
		F(L)=F2(I,J,K,L,S)
	   end do	! L loop
	   C(1)=0.
	   FBND(1)=0.
	   C(LO-1)=0.
	   FBND(LO-1)=0.
	   F(1)=0.
	   UL=LO-2
	   ULL=LO-1
!  dayside BC for IBC=1
	   IF (Ib.EQ.1 .AND. S.EQ.1 .AND. J.GE.J6 .AND. J.LE.J18) THEN
	      F(UPA(I))=FBC(EKEV(K),FFACTOR(I,K,L),FINI(K)*CHI(I,J))
	      UL=UPA(I)-1
	      ULL=UPA(I)-1
	    END IF

	   do L=2,UL
	    C(L)=MUDOT(I,J,L)*VR(I,J)
!	    C(L)=AMIN1(0.99,AMAX1(-0.99,C(L)))
	    X=F(L+1)-F(L)
	    ISIGN=1
	    IF(C(L).NE.ABS(C(L))) ISIGN=-1
	    FUP=0.5*(F(L)+F(L+1)-ISIGN*X)
	    IF (ABS(X).LE.1.E-27) FBND(L)=FUP
	    IF (ABS(X).GT.1.E-27) THEN
	      N=L+1-ISIGN
 	      RR=(F(N)-F(N-1))/X
	      IF (RR.LE.0) FBND(L)=FUP
	      IF (RR.GT.0) THEN
	        LIMITER=AMAX1(AMIN1(2.*RR,1.),AMIN1(RR,2.))
	        CORR=-0.5*(C(L)-ISIGN)*X
	        FBND(L)=FUP+LIMITER*CORR
	      END IF
	    END IF
	   end do	! second L loop

	   do L=2,ULL		! f(i,j,k,1)=f(i,j,k,2)
	    F2(I,J,K,L,S)=F2(I,J,K,L,S)-C(L)*FBND(L)*DMU(L)/WMU(L)   &
     	      +C(L-1)*FBND(L-1)*DMU(L-1)/WMU(L)
	   end do	! third L loop
! dayside BC for IBC=1
	   IF (Ib.EQ.1 .AND. S.EQ.1 .AND. J.GE.J6 .AND. J.LE.J18) THEN 
	     do L=UPA(I),LO
	      F2(I,J,K,L,S)=FBC(EKEV(K),FFACTOR(I,K,L),FINI(K)*CHI(I,J))
	     end do	! yet another L loop
	   ENDIF
!	   F2(I,J,K,LO,S)=F2(I,J,K,LO-1,S)*CONMU2
	  end do	! K loop
	 end do	! I loop
	end do	! J loop

	RETURN
	END
!
! End of subroutine DRIFTMU
!

!***********************************************************************
!		CHAREXCHANGE & ATMOSPHERIC LOSSES
!  Routine calculates the decay of distributions due to charge exchange
!***********************************************************************
	SUBROUTINE CHAREXCHANGE

	use ModHeidiSize
	use ModHeidiIO
	use ModHeidiMain
	use ModHeidiDrifts

	implicit none

	integer :: i,j,k,l
	real :: FL,fbc,FN

	IF (S.GE.2) THEN
	 do L=2,LO
	  do K=2,KO
	   do J=1,JO
	    do I=2,ILMP(J)
	      FN=AMAX1(F2(I,J,K,L,S)*ACHAR(I,K,L,S),1.E-30*FFACTOR(I,K,L))
	      FL=F2(I,J,K,L,S)-FN
	      CEN=CEN+FL*CONSL(K,S)*WE(K)*WMU(L)*DR*DPHI
	      CEE=CEE+FL*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
	      F2(I,J,K,L,S)=FN
	    end do
	   end do
	  end do
	 end do
	END IF
	do K=2,KO
	 do J=1,JO
	  do I=2,ILMP(J)
	   do L=UPA(I),LO
		FN=AMAX1(F2(I,J,K,L,S)*ATLOS(I,K,L,S),1.E-30*FFACTOR(I,K,L))
		FL=F2(I,J,K,L,S)-FN
		ALN=ALN+FL*CONSL(K,S)*WE(K)*WMU(L)*DR*DPHI
		ALE=ALE+FL*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
		F2(I,J,K,L,S)=FN
	   end do
	  end do
	 end do
	end do
	IF (Ib.EQ.1 .AND. S.EQ.1) THEN  ! dayside BC for IBC=1
	  do K=2,KO
	   do J=J6,J18
	    do I=2,IO
	     do L=UPA(I),LO-1
	      F2(I,J,K,L,S)=FBC(EKEV(K),FFACTOR(I,K,L),FINI(K)*CHI(I,J))
	     end do
	    end do
	   end do	! J loop
	  end do	! K loop
	ENDIF

	RETURN
	END
! 
! End of subroutine CHAREXCHANGE
!

!***********************************************************************
!			 	COULMU
!     Routine calculates the Coulomb decay of the distribution function
!            due to pitch angle diffusion
!***********************************************************************
	SUBROUTINE COULMU

	use ModHeidiSize
	use ModHeidiIO
	use ModHeidiMain
	use ModHeidiDrifts
	use ModHeidiWaves

	implicit none

	REAL :: F(NPA),AN,BN,GN,RP,DENOM,RK(NPA),RL(NPA),FBC,BTAW
	INTEGER :: IL,i,j,k,l,ll,UL

	il=1				! =0 for fully implicit
	do J=1,JO
	 do I=2,ILMP(J)
	  do K=2,KO
	   do L=1,LO
	      F(L)=F2(I,J,K,L,S)
	   end do	! 1st L loop

	   RK(2)=0.			! lower b.c.
	   RL(2)=-CONMU1		!   "     "
	   UL=LO-1
	   IF (Ib.EQ.1 .AND. S.EQ.1 .AND. J.GE.J6 .AND. J.LE.J18)    &
     		UL=UPA(I)-1 
	   do L=3,UL
	    BTAW=ATAW(I,J,K,L)+GTAW(I,J,K,L)
	    AN=(ATAI(K,L,S)+ATAE(K,L,S))*XNE(I,J)+IWPI*ATAW(I,J,K,L)
	    BN=(BTAI(K,L,S)+BTAE(K,L,S))*XNE(I,J)+IWPI*BTAW
	    GN=(GTAI(K,L,S)+GTAE(K,L,S))*XNE(I,J)+IWPI*GTAW(I,J,K,L)
	    IF (L.GE.UPA(I)) THEN
	     ll=UPA(I)-1
	     BTAW=ATAW(I,J,K,LL)+GTAW(I,J,K,LL)
	     AN=(ATAI(K,LL,S)+ATAE(K,LL,S))*XNE(I,J)+IWPI*ATAW(I,J,K,LL)
	     BN=(BTAI(K,LL,S)+BTAE(K,LL,S))*XNE(I,J)+IWPI*BTAW
	     GN=(GTAI(K,LL,S)+GTAE(K,LL,S))*XNE(I,J)+IWPI*GTAW(I,J,K,LL)
	    END IF
	    IF (il.EQ.0) THEN
	      AN=2.*AN
	      BN=2.*BN
	      GN=2.*GN
	      RP=F(L)
	    ELSE
	      RP=AN*F(L+1)+(1.-BN)*F(L)+GN*F(L-1)
	    END IF
	    IF (RP.LT.0.) il=0
	    DENOM=BN+GN*RL(L-1)+1.
	    RK(L)=(RP+GN*RK(L-1))/DENOM
	    RL(L)=-AN/DENOM
!	    IF (RK(L).LT.-1.E-29 .OR. ABS(RL(L)+.5).GT.0.50001) THEN
!	       LL=L
!	       IF (L.GE.UPA(I)) LL=UPA(I)-1
!	       print 51,I,J,K,L,LL,IWPI,RK(L),RL(L),AN,BN,GN
!	       print 52,ATAI(K,LL),ATAE(K,LL),ATAW(I,J,K,LL),BTAI(K,LL),   &
!                 BTAE(K,LL),GTAI(K,LL),GTAE(K,LL),GTAW(I,J,K,LL)
!	       print 52,F(L-1),F(L),F(L+1),DENOM,RP,RL(L-1),RK(L-1),i,   &
!      		 XNE(I,J)
!	       STOP
!51		FORMAT (6I4,1P,5E11.3)
!52		FORMAT (1P,8E10.2)
!	    END IF
	   end do	! 2nd L loop

! dayside BC for IBC=1
	   IF (Ib.EQ.1 .AND. S.EQ.1 .AND. J.GE.J6 .AND. J.LE.J18) THEN
	     do L=UPA(I),LO-1
	      F2(I,J,K,L,S)=FBC(EKEV(K),FFACTOR(I,K,L),FINI(K)*CHI(I,J))
	     end do	! 3rd L loop
	   ELSE 				  ! day or night, IBC>1
	     F2(I,J,K,LO-1,S)=RK(LO-1)/(1+RL(LO-1)*CONMU2)
	     do L=LO-2,UPA(I),-1
		F2(I,J,K,L,S)=RK(L)-RL(L)*F2(I,J,K,L+1,S)
	     end do	! 4th L loop
!	     IF (F2(I,J,K,UPA(I),S).LT.0.2*F2(I,J,K,UPA(I)+1,S)) THEN
!	      PRINT 10, 'Bad F2:',I,J,K,UPA(I),T,F2(I,J,K,UPA(I),S),   &
!     		0.2*F2(I,J,K,UPA(I)+1,S),RK(UPA(I)),RL(UPA(I))
!    	      F2(I,J,K,UPA(I),S)=0.2*F2(I,J,K,UPA(I)+1,S)
!	     END IF
	   ENDIF
	   do L=UPA(I)-1,2,-1
		F2(I,J,K,L,S)=RK(L)-RL(L)*F2(I,J,K,L+1,S)
	   end do	! 4th L loop
	   F2(I,J,K,LO,S)=F2(I,J,K,LO-1,S)*CONMU2
	  end do	! K loop
	 end do		! I loop
	end do		! J loop

10	FORMAT(A10,4I4,F10.1,1P,5E11.3)
	RETURN
	END
! 
! End of subroutine COULMU
!

