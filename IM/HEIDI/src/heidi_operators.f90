! File name: heidi_operators_010.f90
! Contains: phase space density calculation routines for HEIDI
!	LLIMIT
!	FCHECK
!	DRIFTR
!	DRIFTP
!	DRECOUL
!	DRIFTMU
!	CHAREXCHANGE
!	COULMU
!======================================================================
!				LLIMIT
!	Sets very small  values to zero to avoid computer zero
!======================================================================
subroutine LLIMIT

  use ModHeidiSize
  use ModHeidiMain

  implicit none

  integer :: i,j,k,l
  !----------------------------------------------------------------------
  do L=1,LO
     do K=1,KO
        do J=1,JO
           do I=1,IO
              if (abs(F2(I,J,K,L,S)).le.1.E-29*FFACTOR(I,j,K,L))   &
                   F2(I,J,K,L,S)=1.E-30*FFACTOR(I,j,K,L)
           end do
        end do
     end do
  end do

end subroutine LLIMIT
!======================================================================
!			   	FCHECK
!	Checks the distribution for negative or infinite results
!======================================================================
subroutine FCHECK(MARK)

  use ModHeidiSize
  use ModHeidiMain
  use ModHeidiIo, ONLY: write_prefix, iUnitStdOut

  implicit none

  integer :: mark
  integer :: Ibad,i,j,k,l,i1,j1,k1,l1,Ibadt,m
  !----------------------------------------------------------------------  
  Ibadt=0

  do m=1,NS
     if (SCALC(m).eq.1) then
     	Ibad=0
	do L=1,Lo
           do k=2,ko
              do j=1,jo
                 do i=2,io
                    if ((F2(I,J,K,L,m).lt.-1.E-29).or.   &   ! -1E-29
                         (F2(I,J,K,L,m)-F2(I,J,K,L,m).ne.0.)) then
!                         (F2(I,J,K,L,m)-F2(I,J,K,L,m).ge.1e-12)) then
                      call write_prefix; write(iUnitStdOut,*) 'Bad F,T,I,J,K,L,MARK:',F2(I,J,K,L,m),   &
                            T,m,I,J,K,L,MARK,Ibad
!                       call write_prefix; write(iUnitStdOut,*) 'Radius:',(F2(I1,J,K,L,m),I1=2,IO,IO/10)
!                       call write_prefix; write(iUnitStdOut,*) 'Azimuth:',(F2(I,J1,K,L,m),J1=1,JO,JO/10)
!                       call write_prefix; write(iUnitStdOut,*) 'Energy:',(F2(I,J,K1,L,m),K1=2,KO,KO/10)
!                       call write_prefix; write(iUnitStdOut,*) 'Pitch angle:',(F2(I,J,K,L1,m),L1=1,LO,LO/10)
                       Ibad=Ibad+1
                       if (Ibad.eq.20) goto 151 
                    end if
                 end do
              end do
           end do
	end do
151 	continue
	Ibadt=Ibadt+Ibad
     end if
  end do

  if (Ibadt.gt.0) then
     call write_prefix; write(iUnitStdOut,*) 'Ibad=',Ibadt
     call write_prefix; write(iUnitStdOut,*) 'A=',A
     call ECFL
     call CON_stop('ERROR in heidi_operators.f90')
  end if

99 format(A,1PE11.3,0PF8.0,7(2X,I4))
98 format(A,200(1PE11.3))

end subroutine FCHECK
!======================================================================
!			     DRIFTR
!     Routine calculate the change of distribution function due to
! 			    radial drift
!======================================================================
subroutine DRIFTR

  use ModHeidiSize
  use ModHeidiIO
  use ModHeidiMain
  use ModHeidiDrifts

  implicit none

  real    :: f01(nt,ne,npa,ns),f02(nt,ne,npa,ns)
  real    :: F(NR+2),FBND(NR),C(NR),LIMITER
  integer :: UR,isign
  integer :: i,j,k,l,n
  real    :: fup,RR,corr,x
  real    :: fbc
  save f01,f02
  !----------------------------------------------------------------------
  !\
  ! Set up injection boundary fluxes
  !/
  if (AMOD(T,TINJ).lt.2.*DT .or. T.eq.TIME) then
     do L=1,LO
        do K=1,KO
           do J=1,JO
              f01(j,k,l,S)=FGEOS(J,K,L,S)*CONF1
              f02(j,k,l,S)=FGEOS(J,K,L,S)*CONF2
           end do
        end do
     end do
  end if				! End BC update

  do J=1,JO
     do K=2,KO		!	< FBND(I) >
        do L=2,LO		!  |____.____|____.____|____.____|
	   do I=1,IO		!  <   F(I)  >
              F(I) = F2(I,J,K,L,S)  ! F - average in cell(i,j,k,l)
	   end do
	   if (VR(IO,J,k,l).ge.0.) then
              FBND(1)=0.	          ! outflow b. c.
              FBND(IO)=F(IO)       ! upwind for the side with no b.c.
              C(IO)=VR(IO,J,k,l)
              UR=IO-1
           else
              FBND(1)=F(2)         ! upwind for the side with no b.c
              C(1)=VR(1,J,k,l)
              UR=IO
              f01(j,k,l,S)=FGEOS(J,K,L,S)*CONF1
              f02(j,k,l,S)=FGEOS(J,K,L,S)*CONF2
              F(IO+1)=f01(j,k,l,s)
              F(IO+2)=f02(j,k,l,s)
           end if
           !	   C(1)=AMIN1(0.99,AMAX1(-0.99,C(1)))
	   do I=2,UR
              C(I)=VR(I,J,k,l)
              !	     C(I)=AMIN1(0.99,AMAX1(-0.99,C(I)))
              X=F(I+1)-F(I)
              ISIGN=1
              if(C(I).ne.abs(C(I))) ISIGN=-1
              FUP=0.5*(F(I)+F(I+1)-ISIGN*X) 	! upwind
              if (abs(X).le.1.E-27) FBND(I)=FUP
              if (abs(X).gt.1.E-27) then
                 N=I+1-ISIGN
                 RR=(F(N)-F(N-1))/X
                 if (RR.le.0) FBND(I)=FUP
                 if (RR.gt.0) then
                    LIMITER=AMAX1(AMIN1(2.*RR,1.),AMIN1(RR,2.))
                    CORR=-0.5*(C(I)-ISIGN)*X
                    FBND(I)=FUP+LIMITER*CORR	! at boundary of cell
                 end if
              end if
           end do ! I loop

           !\
           ! Update the solution for next time step
           !/

           do I=2,ILMP(J)
              F2(I,J,K,L,S)=F2(I,J,K,L,S)-C(I)*FBND(I)+C(I-1)*FBND(I-1)
	   end do ! I loop 
	   do I=ILMP(J)+1,IO
              F2(I,J,K,L,S)=1.E-30*FFACTOR(I,j,K,L)
	   end do
50         format(4I3,1P,10E10.2)

           !\
           ! Count sources and losses out the spatial boundaries
	   !/

           I=ILMP(J)			! Choose I_magnetopause
	   if (C(I).lt.0.) then	! Gain at outer boundary
              RNS=RNS-C(I)*FBND(I)*CONSL(K,S)*WE(K)*WMU(L)*DR*DPHI
              RES=RES-C(I)*FBND(I)*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
	   else				! Loss at outer boundary
              RNL=RNL+C(I)*FBND(I)*CONSL(K,S)*WE(K)*WMU(L)*DR*DPHI
              REL=REL+C(I)*FBND(I)*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
	   end if
	   if (C(1).gt.0.) then		! Gain at inner boundary
              RNS=RNS+C(1)*FBND(1)*CONSL(K,S)*WE(K)*WMU(L)*DR*DPHI
              RES=RES+C(1)*FBND(1)*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
	   else				! Loss at inner boundary
              RNL=RNL-C(1)*FBND(1)*CONSL(K,S)*WE(K)*WMU(L)*DR*DPHI
              REL=REL-C(1)*FBND(1)*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
	   end if

        end do	! End L loop
        !\
        ! Dayside BC for IBC=1
        !/
        if (Ib.eq.1 .and. S.eq.1 .and. J.ge.J6 .and. J.le.J18) then
           do L=UPA(I),LO
              do I=2,IO
                 F2(I,J,K,L,S)=FBC(EKEV(K),FFACTOR(I,j,K,L),FINI(K)*CHI(I,J))
              end do
           end do
        end if
     end do	! K loop
  end do	! J loop
end subroutine DRIFTR
!======================================================================
!			     DRIFTP
!     Routine calculate the change of distribution function due to
! 			 azimuthal drift
!======================================================================
subroutine DRIFTP

  use ModHeidiSize
  use ModHeidiIO
  use ModHeidiMain
  use ModHeidiDrifts

  implicit none

  real    :: FBND(0:NT),F(0:NT+2),C(0:NT),LIMITER
  integer :: i,j,k,l,n,isign,imag
  real    :: fup,RR,corr,x
  real    :: fbc
  !----------------------------------------------------------------------
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
              C(J)=P1(I,J)+P2(I,j,K,L)
              !	    C(J)=AMIN1(0.99,AMAX1(-0.99,C(J)))
              ISIGN=1
              if(C(J).ne.abs(C(J))) ISIGN=-1
              X=F(J+1)-F(J)
              FUP=0.5*(F(J)+F(J+1)-ISIGN*X)
              if (abs(X).le.1.E-27) FBND(J)=FUP
              if (abs(X).gt.1.E-27) then
                 N=J+1-ISIGN
                 RR=(F(N)-F(N-1))/X
                 if (RR.le.0) FBND(J)=FUP
                 if (RR.gt.0) then
                    LIMITER=AMAX1(AMIN1(2.*RR,1.),AMIN1(RR,2.))
                    CORR=-0.5*(C(J)-ISIGN)*X
                    FBND(J)=FUP+LIMITER*CORR
                 end if
              end if
	   end do ! End J loop
	   C(0)=C(JO)
	   FBND(0)=FBND(JO)

	   do J=1,JO
              F2(I,J,K,L,S)=F2(I,J,K,L,S)-C(J)*FBND(J)+C(J-1)*FBND(J-1)
	   end do
        end do	! End L loop

        if (ISW.gt.0) then		! SW can compress magnetopause
           imag=1
           do J=1,JO
              if (ILMP(J+1).lt.I .and. imag.eq.1) then ! Leaving m'sphere
                 imag=0
                 if (C(J).lt.0.) then	! Gain at outer boundary
                    RNS=RNS-C(J)*FBND(J)*CONSL(K,S)*WE(K)*WMU(L)*DR*DPHI
                    RES=RES-C(J)*FBND(J)*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
                 else			! Loss at outer boundary
                    RNL=RNL+C(J)*FBND(J)*CONSL(K,S)*WE(K)*WMU(L)*DR*DPHI
                    REL=REL+C(J)*FBND(J)*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
                 end if
              else if (ILMP(J+1).ge.I .and. imag.eq.0) then ! Reentering
                 imag=1
                 if (C(J).gt.0.) then	! Gain at outer boundary
                    RNS=RNS+C(J)*FBND(J)*CONSL(K,S)*WE(K)*WMU(L)*DR*DPHI
                    RES=RES+C(J)*FBND(J)*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
                 else			! Loss at outer boundary
                    RNL=RNL-C(J)*FBND(J)*CONSL(K,S)*WE(K)*WMU(L)*DR*DPHI
                    REL=REL-C(J)*FBND(J)*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
                 end if
              end if
           end do
        end if

        if (Ib.eq.1 .and. S.eq.1) then
           do J=J6,J18
              do L=UPA(I),LO
                 F2(I,J,K,L,S)=FBC(EKEV(K),FFACTOR(I,j,K,L),FINI(K)*CHI(I,J))
              end do
           end do	! J loop
        end if
     end do	! big K loop
  end do	! big I loop

end subroutine DRIFTP
!======================================================================
!			DRECOUL
!	Routine calculates the change of distribution function due to
!	energization along the drift path and Coulomb energy decay
!======================================================================
subroutine DRECOUL

  use ModHeidiSize
  use ModHeidiIO
  use ModHeidiMain
  use ModHeidiDrifts
  use ModIoUnit, ONLY : UNITTMP_

  implicit none

  real    :: FBND(NE),F(0:NE+2),C(NE),LIMITER,corr,fup,RR,x,CD
  real    :: fbc
  integer :: i,j,k,l,n,isign
  !---------------------------------------------------------------------- 

  do J=1,JO
     do I=2,ILMP(J)
        do L=2,LO
           do k=0,KO+2
              F(K)=0.
           end do
           do k=1,KO
              C(k)=0.
           end do
           f(1:ko)=f2(i,j,1:ko,l,s)
           do K=1,KO
!              C(K)=EDOT(I,J,K,L)*VR(I,J)+(COULE(I,j,K,L,S)+COULI(I,j,K,L,S))*XNE(I,J)
              C(K)=EDOT(I,J,K,L)+(COULE(I,j,K,L,S)+COULI(I,j,K,L,S))*XNE(I,J)
              !	    C(K)=AMIN1(0.99,AMAX1(-0.99,C(K)))
              ISIGN=1
              if (C(K).ne.abs(C(K))) ISIGN=-1
              X=F(K+1)-F(K)
              FUP=0.5*(F(K)+F(K+1)-ISIGN*X)
              if (abs(X).le.1.E-27) FBND(K)=FUP
              if (abs(X).gt.1.E-27) then
                 N=K+1-ISIGN
                 RR=(F(N)-F(N-1))/X
                 if (RR.le.0) FBND(K)=FUP
                 if (RR.gt.0) then
                    LIMITER=AMAX1(AMIN1(2.*RR,1.),AMIN1(RR,2.))
                    CORR=-0.5*(C(K)-ISIGN)*X
                    FBND(K)=FUP+LIMITER*CORR
                 end if
              end if
           end do	! K loop


           do K=2,KO
              F2(I,J,K,L,S)=F2(I,J,K,L,S)-C(K)*FBND(K)*DE(K)/WE(K)   &
                   +C(K-1)*FBND(K-1)*DE(K-1)/WE(K)

           end do

           !\                 
           ! Keep track of sources and losses, separately for drift and CC
           ! Actually a net values for each, adding the sources and
           ! subtracting the losses
           ! Note: particle changes only at the Erange boundaries, while
           ! energy changes throughout the range; E endpoints done in loop
           !/  
!           CD=EDOT(I,J,1,L)*VR(I,J)*DE(1)/WE(2)	! Drift at K=1,2 bnd
           CD=EDOT(I,J,1,L)*DE(1)/WE(2)	! Drift at K=1,2 bnd

           if (CD.gt.0) then
              ESN=ESN+CD*FBND(1)*CONSL(2,S)*WE(2)*WMU(L)*DR*DPHI
           else
              ELN=ELN-CD*FBND(1)*CONSL(2,S)*WE(2)*WMU(L)*DR*DPHI
           end if
           CD=C(1)*DE(1)/WE(2)-CD		! CC at K=1,2 bnd
           ECN=ECN+CD*FBND(1)*CONSL(2,S)*WE(2)*WMU(L)*DR*DPHI
!           CD=EDOT(I,J,KO,L)*VR(I,J)*DE(KO)/WE(KO)	! Drift at K=KO
           CD=EDOT(I,J,KO,L)*DE(KO)/WE(KO)	! Drift at K=KO
           
          if (CD.le.0) then
              ESN=ESN-CD*FBND(KO)*CONSL(KO,S)*WE(KO)*WMU(L)*DR*DPHI
           else 
              ELN=ELN+CD*FBND(KO)*CONSL(KO,S)*WE(KO)*WMU(L)*DR*DPHI
           end if
           CD=C(KO)*DE(KO)/WE(KO)-CD		! CC at K=KO
           ECN=ECN-CD*FBND(KO)*CONSL(KO,S)*WE(KO)*WMU(L)*DR*DPHI
           do K=2,KO				! Now do energy changes
!              CD=EDOT(I,J,K-1,L)*VR(I,J)*DE(K-1)/WE(K)	! Drift at lower bnd
              CD=EDOT(I,J,K-1,L)*DE(K-1)/WE(K)	! Drift at lower bnd 

             if (CD.gt.0) then
                 ESE=ESE+CD*FBND(K-1)*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
              else
                 ELE=ELE-CD*FBND(K-1)*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
              end if
              CD=C(K-1)*DE(K-1)/WE(K)-CD		! CC at lower bnd
              ECE=ECE+CD*FBND(K-1)*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
!              CD=EDOT(I,J,K,L)*VR(I,J)*DE(K)/WE(K)	! Drift at upper bnd
              CD=EDOT(I,J,K,L)*DE(K)/WE(K)	! Drift at upper bnd

             if (CD.le.0) then
                 ESE=ESE-CD*FBND(K)*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
              else
                 ELE=ELE+CD*FBND(K)*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
              end if
              CD=C(K)*DE(K)/WE(K)-CD		! CC at upper bnd 
              ECE=ECE-CD*FBND(K)*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
           end do

        end do	! End L loop
        !\
        ! Dayside BC for IBC=1
        !/
        if (Ib.eq.1 .and. S.eq.1 .and. J.ge.J6 .and. J.le.J18) then
           do L=UPA(I),LO
              do K=2,KO
                 F2(I,J,K,L,S)=FBC(EKEV(K),FFACTOR(I,j,K,L),FINI(K)*CHI(I,J))
              end do
           end do
        end if
     end do	! I loop
  end do	! J loop

end subroutine DRECOUL

!======================================================================
!			DRIFTMU
!	Routine calculates the change of distribution function due to
!		pitch angle changes along the drift path
!======================================================================
subroutine DRIFTMU

  use ModHeidiSize
  use ModHeidiIO
  use ModHeidiMain
  use ModHeidiDrifts

  implicit none

  real :: FBND(NPA),F(NPA),C(NPA),LIMITER,corr,fup,RR,x
  integer :: UL,ULL
  integer :: i,j,k,l,n,isign
  real :: fbc

  !---------------------------------------------------------------------- 
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
           !\
           ! Dayside BC for IBC=1
	   !/
           if (Ib.eq.1 .and. S.eq.1 .and. J.ge.J6 .and. J.le.J18) then
	      F(UPA(I))=FBC(EKEV(K),FFACTOR(I,j,K,L),FINI(K)*CHI(I,J))
	      UL=UPA(I)-1
	      ULL=UPA(I)-1
           end if

	   do L=2,UL
!              C(L)=MUDOT(I,J,L)*VR(I,J)
              C(L)=MUDOT(I,J,K,L) 
             !	    C(L)=AMIN1(0.99,AMAX1(-0.99,C(L)))
              X=F(L+1)-F(L)
              ISIGN=1
              if(C(L).ne.abs(C(L))) ISIGN=-1
              FUP=0.5*(F(L)+F(L+1)-ISIGN*X)
              if (abs(X).le.1.E-27) FBND(L)=FUP
              if (abs(X).gt.1.E-27) then
                 N=L+1-ISIGN
                 RR=(F(N)-F(N-1))/X
                 if (RR.le.0) FBND(L)=FUP
                 if (RR.gt.0) then
                    LIMITER=AMAX1(AMIN1(2.*RR,1.),AMIN1(RR,2.))
                    CORR=-0.5*(C(L)-ISIGN)*X
                    FBND(L)=FUP+LIMITER*CORR
                 end if
              end if
	   end do	! second L loop

	   do L=2,ULL		! f(i,j,k,1)=f(i,j,k,2)
              F2(I,J,K,L,S)=F2(I,J,K,L,S)-C(L)*FBND(L)*DMU(L)/WMU(L)   &
                   +C(L-1)*FBND(L-1)*DMU(L-1)/WMU(L)
	   end do	! third L loop
           !\
           ! Dayside BC for IBC=1
	   !/
           if (Ib.eq.1 .and. S.eq.1 .and. J.ge.J6 .and. J.le.J18) then 
              do L=UPA(I),LO
                 F2(I,J,K,L,S)=FBC(EKEV(K),FFACTOR(I,j,K,L),FINI(K)*CHI(I,J))
              end do	! yet another L loop
	   endif
        end do	! K loop
     end do	! I loop
  end do	! J loop

end subroutine DRIFTMU
!======================================================================
!		CHAREXCHANGE & ATMOSPHERIC LOSSES
!  Routine calculates the decay of distributions due to charge exchange
!======================================================================
subroutine CHAREXCHANGE

  use ModHeidiSize
  use ModHeidiIO
  use ModHeidiMain
  use ModHeidiDrifts

  implicit none

  integer :: i,j,k,l
  real    :: FL,fbc,FN
  !---------------------------------------------------------------------- 
  if (S.ge.2) then
     do L=2,LO
        do K=2,KO
	   do J=1,JO
              do I=2,ILMP(J)
                 FN=AMAX1(F2(I,J,K,L,S)*achar(I,j,K,L,S),1.E-30*FFACTOR(I,j,K,L))
                 FL=F2(I,J,K,L,S)-FN
                 CEN=CEN+FL*CONSL(K,S)*WE(K)*WMU(L)*DR*DPHI
                 CEE=CEE+FL*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
                 F2(I,J,K,L,S)=FN
              end do
	   end do
        end do
     end do
  end if
  do K=2,KO
     do J=1,JO
        do I=2,ILMP(J)
	   do L=UPA(I),LO
              FN=AMAX1(F2(I,J,K,L,S)*ATLOS(I,j,K,L,S),1.E-30*FFACTOR(I,j,K,L))
              FL=F2(I,J,K,L,S)-FN
              ALN=ALN+FL*CONSL(K,S)*WE(K)*WMU(L)*DR*DPHI
              ALE=ALE+FL*CONSL(K,S)*EKEV(K)*WE(K)*WMU(L)*DR*DPHI
              F2(I,J,K,L,S)=FN
	   end do
        end do
     end do
  end do
  if (Ib.eq.1 .and. S.eq.1) then  ! dayside BC for IBC=1
     do K=2,KO
        do J=J6,J18
           do I=2,IO
              do L=UPA(I),LO-1
                 F2(I,J,K,L,S)=FBC(EKEV(K),FFACTOR(I,j,K,L),FINI(K)*CHI(I,J))
              end do
           end do
        end do	! J loop
     end do	! K loop
  endif

end subroutine CHAREXCHANGE

!======================================================================
!			 	COULMU
!     Routine calculates the Coulomb decay of the distribution function
!            due to pitch angle diffusion
!======================================================================
subroutine COULMU

  use ModHeidiSize
  use ModHeidiIO
  use ModHeidiMain
  use ModHeidiDrifts
  use ModHeidiWaves

  implicit none

  real    :: F(NPA),AN,BN,GN,RP,DENOM,RK(NPA),RL(NPA),FBC,BTAW
  integer :: IL,i,j,k,l,ll,UL
  !---------------------------------------------------------------------- 

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
           if (Ib.eq.1 .and. S.eq.1 .and. J.ge.J6 .and. J.le.J18)    &
                UL=UPA(I)-1 
           do L=3,UL
              BTAW=ATAW(I,J,K,L)+GTAW(I,J,K,L)
              AN=(ATAI(i,j,K,L,S)+ATAE(i,j,K,L,S))*XNE(I,J)+IWPI*ATAW(I,J,K,L)
              BN=(BTAI(i,j,K,L,S)+BTAE(i,j,K,L,S))*XNE(I,J)+IWPI*BTAW
              GN=(GTAI(i,j,K,L,S)+GTAE(i,j,K,L,S))*XNE(I,J)+IWPI*GTAW(I,J,K,L)

              if (L.ge.UPA(I)) then
                 ll=UPA(I)-1
                 BTAW=ATAW(I,J,K,LL)+GTAW(I,J,K,LL)
                 AN=(ATAI(i,j,K,LL,S)+ATAE(i,j,K,LL,S))*XNE(I,J)+IWPI*ATAW(I,J,K,LL)
                 BN=(BTAI(i,j,K,LL,S)+BTAE(i,j,K,LL,S))*XNE(I,J)+IWPI*BTAW
                 GN=(GTAI(i,j,K,LL,S)+GTAE(i,j,K,LL,S))*XNE(I,J)+IWPI*GTAW(I,J,K,LL)

              end if
              if (il.eq.0) then
                 AN=2.*AN
                 BN=2.*BN
                 GN=2.*GN
                 RP=F(L)
              else
                 RP=AN*F(L+1)+(1.-BN)*F(L)+GN*F(L-1)
              end if
              if (RP.lt.0.) il=0
              DENOM=BN+GN*RL(L-1)+1.
              RK(L)=(RP+GN*RK(L-1))/DENOM
              RL(L)=-AN/DENOM
           end do	! 2nd L loop
           !\
           ! Dayside BC for IBC=1
           !/
           if (Ib.eq.1 .and. S.eq.1 .and. J.ge.J6 .and. J.le.J18) then
              do L=UPA(I),LO-1
                 F2(I,J,K,L,S)=FBC(EKEV(K),FFACTOR(I,j,K,L),FINI(K)*CHI(I,J))
              end do	! 3rd L loop
           else 				  ! day or night, IBC>1
              F2(I,J,K,LO-1,S)=RK(LO-1)/(1+RL(LO-1)*CONMU2)
              do L=LO-2,UPA(I),-1
                 DENOM=BN+GN*RL(L-1)+1.
                 RK(L)=(RP+GN*RK(L-1))/DENOM
                 RL(L)=-AN/DENOM
                 F2(I,J,K,L,S)=RK(L)-RL(L)*F2(I,J,K,L+1,S)
              end do	! 4th L loop

           endif

           do L=UPA(I)-1,2,-1
              F2(I,J,K,L,S)=RK(L)-RL(L)*F2(I,J,K,L+1,S)
           end do	! 4th L loop

           F2(I,J,K,LO,S)=F2(I,J,K,LO-1,S)*CONMU2
        end do	! K loop
     end do       ! I loop
  end do		! J loop


10 format(A10,4I4,F10.1,1P,5E11.3)
end subroutine COULMU
!======================================================================
