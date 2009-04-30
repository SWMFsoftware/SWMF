! File name: heidi_setup.f90
!
! Contains: input and array setup routines for HEIDI
!	HEIDI_READ
!	CONSTANT
!       BFIELD_SETUP
!	ARRAYS
!	GETKPA
!	GETSWIND
!	THERMAL
!	G
!	FUNT
!	FUNI
!	APPX
!	ACOSD
!	ASIND (unused)
!	COSD (unused)
!	SIND (unused)
!
! Change from 010 to 011: Dipole field definition now a subroutine
!	instead of embedded in ARRAYS
!
! Last Modified: March 2006, Mike Liemohn
!
!======================================================================
!                             HEIDI_READ
!	Read parameters: DT, TMAX, Species, Storm, .......
!======================================================================

subroutine heidi_read

  use ModHeidiSize
  use ModHeidiIO
  use ModHeidiMain
  use ModIoUnit, only : UNITTMP_
  use ModHeidiInput, only: set_parameters
  use ModProcIM, only: iComm

  implicit none

  integer          :: k,i
  character(len=1) :: header

  !------------------------------------------------------------------------

  call set_parameters

  Dt = DTmax
  call write_prefix; write(iUnitStdOut,*) ' year,month,day,UT',year,month,day,UT

  write(*,*) 'DT,TMAX,TINT,TIME',DT,TMAX,TINT,TIME
  write(*,*) 'IO,JO,KO,LO,ISO',IO,JO,KO,LO,ISO
  write(*,*)' ELB,SWE,RW,HMIN',ELB,SWE,RW,HMIN
  write(*,*) 'ISTORM,IKP,IPA,IFAC,IST,IWPI,ISW,IA,ITHERMINIT',&
       ISTORM,IKP,IPA,IFAC,IST,IWPI,ISW,IA,ITHERMINIT
  write(*,*) '(SCALC(k),k=1,NS)',(SCALC(k),k=1,NS)
  write(*,*) 'YEAR,month,day,UT,R,AP,KP',YEAR,month,day,UT,R,AP,KP
  write(*,*) '(INI(k),k=1,NS)',(INI(k),k=1,NS)
  write(*,*) '(IBC(k),k=1,NS)',(IBC(k),k=1,NS)
  write(*,*) 'TINJ,Ab,Eob',TINJ,Ab,Eob
  write(*,*) '(IRES(k),k=1,15)',(IRES(k),k=1,15)
  write(*,*) 'NAME',NameRun


  TimeArray(1) = year
  TimeArray(2) = month
  TimeArray(3) = day
  TimeArray(4) = UT
  TimeArray(5) = 0.0
  TimeArray(6) = 0.0
  TimeArray(7) = 0.0

  ISWB=ISW

  NSTEP=nint(TMAX/DT/2.)        ! time splitting
  Ib=0
  if (IBC(1).eq.1) then
     Ib=1			! Loss cone BC
  end if

  ithermfirst=1		! So we do the setup routines in THERMAL

  if (IKP.eq.4 .or. IA.eq.2) then  ! Read in MBI file
     open(UNITTMP_,FILE=NameInputDirectory//trim(NameRun)//'_Le.dat',status='old')
     do I=1,3
        read (UNITTMP_,*) header
     end do
     read (UNITTMP_,*) ILAMBE,LAMGAM
     read (UNITTMP_,*) header
     read (UNITTMP_,*) header
     do I=1,ilambe
        read (UNITTMP_,*) TLAME(I),LAMBE(I)
     end do
     close(UNITTMP_)
     TLAME(1:ILAMBE)=TLAME(1:ILAMBE)*86400.
     TLAME(1:ILAMBE)=TLAME(1:ILAMBE)-3600.  ! Forward shift 1 h
  else
     LAMGAM=2.
  end if

  if (IA.eq.4 .or. IA.eq.7 .or. IA.ge.10) then ! Read in PC Potential File

     open(UNITTMP_,FILE=NameInputDirectory//trim(NameRun)//'_ppc.dat',status='old')
     do I=1,3
        read (UNITTMP_,*) header
     end do
     read (UNITTMP_,*) IPPC
     read (UNITTMP_,*) header
     read (UNITTMP_,*) header
     do I=1,ippc
        read (UNITTMP_,*) TPPC(I),PPC(I)
     end do
     close(UNITTMP_)
     TPPC(1:IPPC)=TPPC(1:IPPC)*86400. ! Convert to seconds
     PPC(1:IPPC)=PPC(1:IPPC)*1.E3       ! Convert to Volts
     call write_prefix; write(iUnitStdOut,*) 'PPC:',IPPC,TPPC(1),TPPC(IPPC),PPC(1),PPC(IPPC)
  end if

end subroutine heidi_read
!======================================================================
!				CONSTANT
!			    Define constants
!======================================================================
subroutine CONSTANT(NKP)

  use ModHeidiSize
  use ModHeidiIO
  use ModHeidiMain
  use ModIoUnit, only : UNITTMP_

  implicit none

  real              :: kpt,DUT
  integer           :: I,NKP
  character(len=80) :: header
!------------------------------------------------------------------------
  !.......Read Kp history of the modeled storm
  if (IKP.ge.3) then
     open(UNITTMP_,FILE=NameInputDirectory//trim(NameRun)//'_kp.in',STATUS='OLD') 
     read(UNITTMP_,10) HEADER
10   format(A80)
     do I=1,NSTEP/NKP+2
        read(UNITTMP_,*) DAYR(I),DUT,RKPH(I),F107R(I),APR(I),RSUNR(I)
     enddo
     close(UNITTMP_)
  end if

  KPT=-1./3./3600.   ! model rate of decay of Kp in s-1
  DKP=KPT*DT*2.	   ! discrete step size of Kp
  ME=7.9E15	   ! Magnetic moment of the earth
  RE=6.371E6	   ! Earth's radius (m)
  MP=1.673E-27       ! Mass of H+ in kg
  Q=1.6E-19	   ! Electron charge, eV -> J conversion
  PI=3.141592654

  !	Conversion factor from phase space distribution F[s^3/km^6] into
  !	the equatorial directional flux: FLUX[1/cm^2/s/keV/ster]
  !	FLUXFACT = FLUX/F/E[keV]

  FLUXFACT(1)=6.1847E6
  FLUXFACT(2)=1.833847
  FLUXFACT(3)=0.1146154
  FLUXFACT(4)=7.16346E-3

end subroutine CONSTANT

!======================================================================
!				BFIELD_SETUP
!			Set up all the arrays
!======================================================================

subroutine BFIELD_SETUP(cone)

  use ModHeidiSize
  use ModHeidiMain
  use ModHeidiWaves
  use ModHeidiIO

  implicit none

  real     :: CONE(NR+4),degrad,camlra,asind
  integer  :: i,iml
  external :: asind
  !------------------------------------------------------------------------
  degrad=pi/180.
  do i=1,IO
     do IML=1,ISO
        camlra=amla(iml)*degrad
        BE(I,IML)=0.32/LZ(I)**3*sqrt(1.+3.*sin(camlra)**2)  &
             /cos(camlra)**6				! in gauss
     enddo
     BFC(I)=ME/Z(I)**3/40./sqrt(PI*Q)		! in SI units
  end do

  !.......CONE - pitch angle loss cone in degree
  do i=1,io	
     CONE(I)=ASIND(sqrt(((RE+HMIN)/Z(I))**3   &
          /sqrt(4.-3.*((RE+HMIN)/Z(I)))))
  end do
  CONE(IO+1)=2.5    ! to calcul PA grid near 0 deg for IPA=1
  CONE(IO+2)=1.5
  CONE(IO+3)=1.
  CONE(IO+4)=0.

end subroutine BFIELD_SETUP

!======================================================================
!				ARRAYS
!			Set up all the arrays
!======================================================================
subroutine ARRAYS

  use ModHeidiSize
  use ModHeidiIO
  use ModHeidiMain
  use ModHeidiWaves

  implicit none

  integer :: I,J,K,L,IDL1,IC,IML,LL
  real :: DPA,SEG1,SEG2,RWU,CONV,degrad,camlra,muboun,spa
  real :: ASIND,COSD,ACOSD,FUNT,CON1
  real :: sind
  external :: sind, ASIND,COSD,ACOSD, funt
  real :: amla0(Slen)

  real :: CONE(NR+4),PA(NPA),MUB(NPA),sumd,sumw,LZMAX,LZMIN

  !	DATA M1/5.4462E-4,1.,4.,16./	! Mi/MAS(H+) for e-, H+, He+, O+
  data amla0/0.0,0.2,0.4,0.6,0.8,1.0,1.5,2.0,2.5,3.0,3.5, &
       4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0/
  !------------------------------------------------------------------------
  
  M1(1) = 5.4462E-4
  M1(2) = 1.0
  M1(3) = 4.0
  M1(4) = 16.0

  do i=1,Slen
     AMLA(i) = AMLA0(i)
  enddo

  LZMIN=2.
  LZMAX=6.5
  DL1=((LZMAX-LZMIN)/(IO-2))	! Grid size of L shell
  !        new code for Cray  ===(length 1 line)====ab==
  IDL1=nint(DL1*1000.)
  if(DL1.ge.0.25 .and. mod(IDL1,250).ne.0) then
     call write_prefix; write(iUnitStdOut,*)  'Error : 0.25 is not a factor of DL1 '
     !WRITE(6,*) ' Error : 0.25 is not a factor of DL1 '
     call CON_stop('ERROR in heidi_setup.f90')
  else if (DL1.lt.0.25 .and. mod(250,IDL1).ne.0) then
     call write_prefix; write(iUnitStdOut,*)  'Error : DL1 is not a factor of 0.25'
     !WRITE (6,*) ' Error : DL1 is not a factor of 0.25'
     call CON_stop('ERROR in heidi_setup.f90')
  end if

  DR=DL1*RE               !Grid size for Z=RO
  do I=1,IO
     LZ(I)=LZMIN+(I-2)*DL1
  end do
  Z(1:IO)=RE*LZ(1:IO)
  degrad=pi/180.

  call BFIELD_SETUP(cone)

  DPHI=2.*PI/JO		      ! Grid size for local time [rad]
  if(mod(NLT,JO).ne.0) then
     call write_prefix; write(iUnitStdOut,*)  'Error : JO is not a factor of NLT '
     !WRITE(6,*) ' Error : JO is not a factor of NLT '
     call CON_stop('ERROR in heidi_setup.f90')
  end if

  do J=1,JO+1
     PHI(J)=(J-1)*DPHI	! Magnetic local time in radian
  end do
  MLT(1:JO+1)=PHI(1:JO+1)*12./PI    ! Magnetic local time in hour 

  do I=1,NS
     MAS(I)=MP*M1(I)		! Mass of each species (kg)
  end do
  !.......Calculate Kinetic Energy EKEV [keV] at cent
  WE(1)=abs(SWE)	!  |_._|___.___|____.____|______.______|
  !    .     <   DE   >    <      WE     >
  !   EKEV                EBND
  EBND(1)=ELB+WE(1)
  !	IF (SWE.GT.0.) THEN	! Growing WE
  do K=1,KO-1
     WE(K+1)=WE(K)*RW	! WE(K) [keV] is a power series
  end do
  if (IST.eq.1) then	! Wide dE at top for PE runs
     RW=1.35			 
     do K=KO-10,KO-1		
        WE(K+1)=WE(K)*RW
     end do
  else if (IST.eq.2) then ! Wide dE at bottom for PSE runs
     RW=1.057			
     do K=KO,16,-1	
        WE(K)=WE(K-15)
     end do
     do K=15,1,-1		
        WE(K)=WE(K+1)/RW
        EBND(1)=EBND(1)-WE(K+1)
     end do
     ELB=EBND(1)-WE(1)
  end if
  !	ELSE			! Set WE from PE source
  !	   IF (SWE.EQ.-.001) THEN
  !	      WE(2:39)=0.001
  !	      WE(40:45)=0.01
  !	      WE(46:70)=0.02
  !	      do k=71,KO
  !		 WE(K)=WE(K-1)*RW
  !	      end do
  !	   ELSE
  !	      WE(2:20)=0.002
  !	      WE(21:26)=0.01
  !	      WE(27:51)=0.02
  !	      do k=52,KO
  !		 WE(K)=WE(K-1)*RW
  !	      end do
  !	   END IF
  !	END IF

  do k=1,KO-1
     EBND(K+1)=EBND(K)+WE(K+1)	! E[keV] at bound of grid
  end do
  DE(1:KO-1)=0.5*(WE(1:KO-1)+WE(2:KO))
  DE(KO)=WE(KO)			! 0.5*WE(KO)*(1.+RW)
  EKEV(1)=ELB+0.5*WE(1)
  do K=1,KO-1
     EKEV(K+1)=EKEV(K)+DE(K)	! E[keV] at cent of grid
  end do
  do S=1,NS
     V(1:ko,s)=sqrt(2.*EKEV(1:ko)*1000.*Q/MAS(S))    ! Vel [m/s] at cent
     VBND(1:ko,s)=sqrt(2.*EBND(1:ko)*1000.*Q/MAS(S)) ! Vel [m/s] at bound
  end do

  !.......PA is equatorial pitch angle in deg - PA(1)=90, PA(LO)=0.
  !.......MU is cosine of equatorial PA
  if (IPA.ne.1) then	! Constant DPA based on LO
     DPA=-90./(LO-1)	! Grid size for pitch angle in degree
     PA(1)=90.
     MU(1)=0.
     do L=1,LO-2
        PA(L+1)=PA(L)+DPA
     end do
     do l=1,lo-1
        MUB(l)=COSD(PA(l)+DPA/2.)
     end do
     WMU(1)=2*MUB(1)
     WMU(2:lo-1)=MUB(2:lo-1)-MUB(2:lo-1)
     DMU(1:lo-2)=0.5*(WMU(2:lo-1)+WMU(1:lo-2))
     do l=1,lo-2
        MU(L+1)=MU(L)+DMU(L)
     end do
     do l=2,LO-1
        PA(l)=ACOSD(MU(l))
     end do
     PA(LO)=0.
     MU(LO)=1.
     WMU(LO)=2*(1.-MUB(LO-1))
     DMU(LO-1)=MU(LO)-MU(LO-1)
     DMU(LO)=WMU(LO)
     !.......Find the pitch angle nearest to the edge of the loss cone
     do I=1,IO
        UPA(I)=0
        L=0
        do while (UPA(I).eq.0)
           L=L+1
           if(PA(L).le.CONE(I)) then
              UPA(I)=L
              SEG1=CONE(I)-PA(L)
              SEG2=PA(L-1)-CONE(I)
              if (SEG2.lt.SEG1) UPA(I)=L-1
           end if
        end do	! While loop
     end do	! I loop

  else		! IPA(1): PA grid w/ pts on the loss cone edges
     LO=71		! For this grid, LO must be 71
     PA(1)=90.	! |__.__|___.___|____.____|_____._____|
     MU(1)=0.	!    MU     <  DMU   >    <    WMU    >
     PA(LO)=0.
     MU(LO)=1.
     RWU=0.98
     WMU(1)=(MU(LO)-MU(1))/32
     do L=1,46
        WMU(L+1)=WMU(L)*RWU
     end do
     DMU(1:46)=0.5*(WMU(1:46)+WMU(2:47))
     do l=1,46
        MU(L+1)=MU(L)+DMU(L)
     end do
     do l=2,47
        PA(l)=ACOSD(MU(l))
     end do
     PA(48)=18.65
     MU(48)=COSD(PA(48))
     DMU(47)=(MU(48)-MU(47))
     UPA(1)=44	! Loss cone indices determined along with grid
     UPA(2)=47
     IC=3 
     do L=48,LO-1
        PA(L+1)=CONE(IC)
        if(L.eq.49) then
           PA(50)=16.
        else 
           if (IC.le.IO) UPA(IC)=L+1
           IC=IC+1
        endif
        MU(l+1)=COSD(PA(l+1))
     end do
     DMU(48:lo-1)=(MU(49:lo)-MU(48:lo-1))	! Grid size in cos(PA) 
     do l=48,55
        WMU(L)=2.*(DMU(L-1)-0.5*WMU(L-1))
     end do
     WMU(56:lo-1)=0.5*(DMU(56:lo-1)+DMU(55:lo-2))
     DMU(LO)=0.				! DMU(LO-1)
     WMU(LO)=MU(LO)-MU(LO-1)-0.5*WMU(LO-1)	! DMU(LO-1)
     WMU(1)=0.5*WMU(1)
  end if

  do L=1,LO-1
     MUBOUN=MU(L)+0.5*DMU(L)
     PAbn(L)=ACOSD(MUBOUN)
  end do
  PAbn(LO)=0.
  SUMW=0.
  SUMD=0.
  do L=1,LO
     SUMW=SUMW+WMU(L)
     SUMD=SUMD+DMU(L)
  end do


  call write_prefix; write(iUnitStdOut,*) 'SUMS:',SUMW,SUMD
  call write_prefix; write(iUnitStdOut,*) 'LO:',MU(LO-1),MU(LO),DMU(LO),WMU(LO),MU(LO-1)+.5*DMU(LO-1)

  !print *, '1:',MU(1),MU(2),DMU(1),WMU(1),MU(2)-.5*DMU(2)
  !**  calculate pitch angles for mlat
  do I=1,IO
     do IML=1,ISO
        do L=1,LO
           spa=sqrt(SIND(PAbn(L))**2*BE(i,iml)/BE(i,1))
           if (spa.gt.1.0) spa=1.0
           ZRpabn(i,L,iml)=asin(spa)
           if (spa.eq.1.0) ZRpabn(i,L,iml)=-1.0
        end do
     end do
  end do

  !.......Define conversion factors
  !.......FFACTOR is ratio of phase space F to F2 in conservative space
  !.......IFAC indicates conversion to flux or distribution function
  if (IFAC.eq.1) then	! Flux function
     do K=1,KO
        FFACTOR(1:io,K,1)=0.
        do L=2,LO
           do i=1,io
              FFACTOR(I,K,L)=LZ(i)*LZ(i)/sqrt(EKEV(K))*MU(L)*FUNT(MU(L))
           end do
        end do
     end do
  else 			! Distribution function
     do K=1,KO
        FFACTOR(1:io,K,1)=0.
        do L=2,LO
           do i=1,io
              FFACTOR(I,K,L)=LZ(i)*LZ(i)*sqrt(EKEV(K))*MU(L)*FUNT(MU(L))
           end do
        end do
     end do
  end if
  do K=1,KO
     do i=1,io
        FFACTOR(I,K,1)=FFACTOR(I,K,2)
     end do
  end do

  !.......Variables for pressure and anisotropy calculations
  do L=1,LO
     do K=1,KO
        do i=1,io
           ERNM(I,K,L)=WMU(L)/FFACTOR(I,K,L)
           EPMA(I,K,L)=ERNM(I,K,L)*MU(L)*MU(L)
           EPME(I,K,L)=ERNM(I,K,L)-EPMA(I,K,L)
        end do
     end do
  end do
  CONV=1.
  do s=1,NS
     do k=1,KO
        if (IFAC.eq.1) CONV=1./(EKEV(K)*FLUXFACT(S))
        ERNH(K,S)=WE(K)*sqrt(EKEV(K))*CONV
        EPP(K,S)=ERNH(K,S)*EKEV(K)
     end do
  end do

  !.......zero out wave diff coeff
  !  DO 20 I=1,NR
  !     DO 20 J=1,NT
  !	DO 20 K=1,NE
  !           DO 20 L=1,NPA
  !              ATAW(I,J,K,L)=0.
  !              GTAW(I,J,K,L)=0.
  !20            CONTINUE
  do I=1,NR
     do J=1,NT
	do K=1,NE
           do L=1,NPA
              ATAW(I,J,K,L)=0.
              GTAW(I,J,K,L)=0.
              !.......to keep F constant at boundary
              CONF1=((LZ(IO)+DL1)/LZ(IO))**2
              CONF2=((LZ(IO)+2.*DL1)/LZ(IO))**2

              !.......FACMU is the PA dependent factor due to conversion to F2
              do LL=1,lo
                 FACMU(LL)=FUNT(MU(LL))*MU(LL)
              end do

              !.......to keep F const at 90 & at 0 deg
              CONMU1=FACMU(2)/FACMU(3)
              CONMU2=FACMU(LO)/FACMU(LO-1)

              !	CONSL is for source/loss calculations, integrating F to give
              !	particles or keV
              !	CONSL is 8*PI*SQRT(2)*(kg->amu * keV->J)^1.5 *(m->km)^6 * Re^2 
              !	  which should be ~4.3E13/M(amu)^1.5
              do S=1,NS
                 CON1=8.*sqrt(2.)*PI*(Q*1.E3/MP/M1(S))**1.5*1.E-18*RE**2
                 CONSL(1:KO,S)=CON1
                 if (IFAC.eq.1) CONSL(1:KO,S)=CONSL(1:KO,S)/FLUXFACT(S)/EKEV(1:KO)
              end do

              return
           end do
        end do
     end do
  end do
end subroutine ARRAYS

!======================================================================
!                             GETKPA
!	Finds KP and other parameters as well as convection strength A
!	Note [A]= V/m * 1/Re^(lamgam-1), last term =1 (in units of L)
!  IKP choices:
!	IKP=0  Static Kp
!	IKP=1  Decrease by DKP until Kp=1
!	IKP=2  Read from table
!	IKP=3  Read from file
!	IKP=4  Read from MBI file (also read kp file for other inputs)
!======================================================================

subroutine GETKPA(i3,nst,i2,nkp)

  use ModHeidiSize
  use ModHeidiIO
  use ModHeidiMain

  implicit none

  integer  :: i3,nst,i4,i2,nkp,JKP,ii
  real     :: KPN,KPO,lambdae,KPtab(48),tol,RSUN,KPP,XKP
  real     :: cosd
  external :: cosd
  save tol,i4,XKP,KPP,KPN,KPO

  !	Data for simulation of the CRRES observations in 1/1991
  data KPtab/1.667,2.000,2.000,2.333,1.777,3.333,4.667,3.000, &
       2.667,3.667,2.333,2.000,2.000,1.667,1.000,3.333, &
       0.333,1.333,1.000,1.667,2.000,1.333,2.333,3.000, &
       1.000,1.667,1.333,1.333,0.667,1.667,1.000,2.000, &
       2.000,0.667,1.333,1.667,0.667,0.333,0.333,1.667, &
       0.333,0.333,0.333,0.667,1.333,2.000,2.000,1.000/
  !	DATA KPtab/1.667,2.000,2.000,2.333,1.777,3.333,4.667,3.000, &
  !                  2.667,3.667,2.333,2.000,2.000,1.667,1.000,3.333, &
  !     		   0.333,1.000,1.000,1.000,1.000,1.000,1.000,1.000, &
  !     		   1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000, &
  !     		   1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000, &
  !     		   1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000/

  !  IKP=0, keep Kp constant, otherwise...
  !------------------------------------------------------------------------
  
  if (IKP.eq.1) then	! IKP=1, Degrades until Kp<=1
     if(KP.gt.1.) KP=KP+DKP
  else if (IKP.eq.2) then	! IKP=2, read from table
     JKP=min(48,int(T/10800.)+1)
     KP=KPtab(JKP) 
  else if (IKP.ge.3) then	! IKP=3 Read from file
     if (mod(I3-1,NKP).eq.0 .or. I3.eq.NST) then
        KPP=RKPH(MAX0(1,I2-1))
        KPO=RKPH(I2)
        KPN=RKPH(MIN0(NSTEP/NKP+2,I2+1))
        XKP=.125*(10.*KPO-KPP-KPN)
        TOL=T
        DAY=DAYR(I2)
        AP=APR(I2)
        F107=F107R(I2)
        RSUN=RSUNR(I2)
        I2=I2+1
        KPN=RKPH(I2)	
     end if
     if (T-TOL.lt.3600.) then
        KP=.5*(KPP+KPO)+(XKP-.5*(KPP+KPO))*(T-TOL)/3600.
     else if (T-TOL.lt.7200.) then
        KP=XKP
     else
        KP=XKP+(.5*(KPO+KPN)-XKP)*(T-TOL-7200.)/3600.
     end if

     !            KP=KPO+(KPN-KPO)*(T-TOL)/10800.
  end if
  if (IKP.eq.4 .or. IA.eq.2) then   ! Midnight Boundary Index
     if (I3.eq.NST) then
        I4=0
        ii=2
        do while (I4.eq.0)
           if (TLAME(ii).ge.T) then
              I4=ii-1
           else if (ii.eq.ILAMBE-1) then
              I4=ILAMBE-1
           else
              ii=ii+1
           end if
        end do
     else if (T+2.*DT.ge.TLAME(ILAMBE)) then
        I4=ILAMBE-1
     else if (T.ge.TLAME(I4+1)) then
        I4=I4+1
     end if
     if (T.le.TLAME(1)) then
        LAMBDAE=LAMBE(1)
     else if (T.ge.TLAME(ILAMBE)) then
        LAMBDAE=LAMBE(ILAMBE)
     else
        LAMBDAE=LAMBE(I4)+(LAMBE(I4+1)-LAMBE(I4))*(T-TLAME(I4)) &
             /AMAX1(.1*DT,TLAME(I4+1)-TLAME(I4))
     end if
     if (IKP.eq.4) KP=LAMBDAE
  end if
  A=7.05E-6/(1.-0.159*KP+0.0093*KP**2)**3
  if (IA.eq.2) then  ! A from MBI transformation
     A=1.44E-2/LAMGAM*COSD(LAMBDAE)**(2.*(LAMGAM+1.))
  end if
end subroutine GETKPA

!======================================================================
!                             GETSWIND
!	Reads swind.dat each time step for solar wind parameters
!======================================================================
subroutine GETSWIND

  use ModHeidiSize
  use ModHeidiIO
  use ModHeidiMain
  use ModIoUnit, only : io_unit_new 


  implicit none

  real             :: T1,T2,BZ1,BZ2,MD1,MD2
  real             :: U1,U2,FAC,FLO,LMPO,AMPO
  real             :: BY1,BY2,BZ,BT
  integer          :: i, j, k, l
  character(len=8) :: header
  save T1,T2,BZ1,BZ2,MD1,MD2,U1,U2,BY1,BY2

  real :: bx
  integer :: iUnit
  !------------------------------------------------------------------------
  iUnit = io_unit_new()
  ILold(1:JO)=ILMP(1:JO)
  if (T.eq.TIME) then
     T2=TIME-1.
     T1=T2
     open(UNIT=iUnit,FILE=NameInputDirectory//trim(NameRun)//'_sw1.in',status='old')
     do I=1,6			! 6 lines of header material
        read (iUnit,*) HEADER
     end do
  end if
  if (T2.lt.T) then
     do while (T2.le.T)	! Best if final T2 > final T
        T1=T2
        BY1=BY2
        BZ1=BZ2
        MD1=MD2
        U1=U2
        read (iUnit,*,IOSTAT=I) T2,BT,BX,BY2,BZ2,MD2,U2
        if (I.lt.0) T2=TIME+2*DT*(NSTEP+1)
        if (T.eq.TIME) then			! In case T2>T already
           T1=TIME
           BY1=BY2
           BZ1=BZ2
           MD1=MD2
           U1=U2
        end if
     end do
  end if
  FAC=(T-T1)/(T2-T1)			! Linearly interpolate
  BYSW=FAC*BY2+(1.-FAC)*BY1               ! in nT
  BZSW=FAC*BZ2+(1.-FAC)*BZ1		! in nT
  MDSW=(FAC*MD2+(1.-FAC)*MD1)*MP*1.E6	! in kg/m3
  USW=(FAC*U2+(1.-FAC)*U1)*1.E3		! in m/s
  DPSW=(FAC*MD2*U2*U2+(1.-FAC)*MD1*U1*U1)*(MP*1.E21)  ! in nPa

  !...Set up for source/loss calculation, magnetopause location vs MLT
  LMPO=(10.22+1.29*tanh(0.184*(BZSW+8.14)))*DPSW**(-0.15152)
  AMPO=(0.58-0.007*BZSW)*(1.+0.024*ALOG(DPSW))
  LMP(2:JO)=LMPO*(2./(1.-cos(PHI(2:JO))))**AMPO
  LMP(1)=LMP(2)
  do J=1,JO
     if (LMP(J).lt.LZ(IO)) then
        ILMP(J)=0		! ILMP=last closed I at this J
        I=1
        do while (ILMP(J).eq.0)
           I=I+1
           if (LMP(J).lt.LZ(I)) ILMP(J)=I-1
        end do
     else
        ILMP(J)=IO
     end if
  end do

end subroutine GETSWIND

!======================================================================
!				THERMAL
!	Converts thermal plasma densities from Craig's grid to ours
!	Rewritten for mram04.f to use Dan Ober's plasmasphere code:
!	  the Dynamic Global Core Plasma Model (DGCPM)
!	  including passage of our electric potentials to the DGCPM
!======================================================================
subroutine THERMAL

  use ModHeidiSize
  use ModHeidiMain
  use ModHeidiDGCPM
  use ModHeidiIO
  use ModProcIM

  implicit none

  integer           :: i,j,i1,j1,l,jj,jjj
  integer           :: j2,i2,ier,ntimestep,n
  real              :: EO,FAC,FACI,sind
  real              :: delt, par(2)
  real              ::chi1, kpinit, dtimestep
  character(len=80) :: filename
  external          :: sind
  !------------------------------------------------------------------------

  par(1)=-1.0
  par(2)=-1.0
  ntimestep=10
  dtimestep=1./ntimestep

  !  If ITHERMFIRST=1, then do initial setup for the plasmasphere code and return
  if (ithermfirst.eq.1) then
     if (me_world.eq.0) then
        call initmain()
        call getgrid(vthetacells,nthetacells,vphicells,nphicells)
        call getxydipole(nthetacells,nphicells,vthetacells,vphicells, &
             gridx,gridy,gridoc)
        do i=1,nthetacells
           vlzcells(i)=1./sind(vthetacells(i))**2
        enddo
        do j=1,nphicells
           vmltcells(j)=vphicells(j)/15.
        enddo
     endif

     call MPI_BARRIER(iComm,iError)
     call MPI_Bcast(vthetacells,nthetacells,MPI_Real,0,iComm,iError)
     call MPI_Bcast(vlzcells,nthetacells,MPI_Real,0,iComm,iError)
     call MPI_Bcast(vphicells,nphicells,MPI_Real,0,iComm,iError)
     call MPI_Bcast(vmltcells,nphicells,MPI_Real,0,iComm,iError)
     call MPI_Bcast(gridx,nthetacells*nphicells,MPI_Real,0,iComm,iError)
     call MPI_Bcast(gridy,nthetacells*nphicells,MPI_Real,0,iComm,iError)
     call MPI_Bcast(gridoc,nthetacells*nphicells,MPI_Real,0,iComm,iError)
     ithermfirst=2
     return

  end if

  !  Call the plasmasphere code and interpolate onto our spatial grid
  if (me_world.eq.0) then
     delt=2*DT*dtimestep
     if (ithermfirst.eq.2) then
        if (itherminit.eq.1) then  ! Read initial plasmasphere from file
           filename=trim(NameRun)//'_dgcpm.in'

           call write_prefix; write(iUnitStdOut,*) 'Reading in plasmasphere file: ',filename
           call loadplasmasphere(filename)
           call getdensity(vthetacells,nthetacells,vphicells,nphicells, &
                dendgcpm)
           delt=0.
        else	 ! initialize plasmasphere with 48 h of low activity
           kpinit = 1.0 
           chi1 = 7350.0 / (9.0 - kpinit)
           delt=48.0*60.0*60.0
           do i = 1, nthetacells
              do j = 1, nphicells
                 fac = vlzcells(i) * sind(vphicells(j))
                 potdgcpm(i,j) = chi1 * fac
              enddo
           enddo
        endif  ! itherminit
     endif    ! ithermfirst
     if (delt.gt.0.) then
        call setpot(vthetacells,nthetacells,vphicells,nphicells,potdgcpm)
        call write_prefix; write(iUnitStdOut,*)  'Calling plasmasphere:',potdgcpm(nthetacells,nphicells/4), &
             potdgcpm(nthetacells,3*nphicells/4)
        do n=1,ntimestep
           call plasmasphere(delt,par)
        end do
        call getdensity(vthetacells,nthetacells,vphicells,nphicells,dendgcpm)
     endif
     do j=1,JO
        do i=1,IO
           call LINTP2(vlzcells,vmltcells,dendgcpm,nthetacells,nphicells, &
                LZ(I),MLT(J),XNE(I,J),IER)
           XNE(I,J)=1.E-6*XNE(i,J)
        end do	! I loop
     end do	! J loop
  endif		! me_world=0

  call MPI_BARRIER(iComm,iError)
  call MPI_Bcast(XNE,NR*NT,MPI_Real,0,iComm,iError)

  !.......Set F2(K=1) to the thermal plasma flux level (PE runs)
  if (IST.eq.1 .and. SCALC(1).eq.1) then
     EO=0.001	! Te=1 eV
     do l=1,LO
        do j=1,JO
           do i=1,io
              F2(I,J,1,L,1)=sqrt(5.*Q/MAS(1)*(PI*EO)**3)*XNE(I,J)*EKEV(1)* &
                   exp(-EKEV(1)/EO)*FFACTOR(I,1,L)
           end do
        end do
     end do
  end if

  ithermfirst = 0

end subroutine THERMAL

!======================================================================
!                        Function G(X) 
!======================================================================
real function G(X)

  use ModHeidiSize
  use ModHeidiMain
  use ModConst, only : cPi

  implicit none

  real    :: x,G1,ERF
  integer :: IER
  external :: ERF
 !------------------------------------------------------------------------ 

  G1=ERF(X,IER)-2.*X/sqrt(cPi)*exp(-X*X)
  G=G1/2./X/X

end function G
!======================================================================
!                       Function T(X) 
!	 function f(y) taken from Ejiri, JGR,1978
!======================================================================
real function FUNT(X)		! X is cos of equat pa = MU

  use ModHeidiSize
  use ModHeidiMain

  implicit none

  real :: Y,X,ALPHA,BETA,A1,A2,A3,A4
  !------------------------------------------------------------------------ 

  Y=sqrt(1-X*X)
  ALPHA=1.+ALOG(2.+sqrt(3.))/2./sqrt(3.)
  BETA=ALPHA/2.-PI*sqrt(2.)/12.
  A1=0.055
  A2=-0.037
  A3=-0.074
  A4=0.056
  FUNT=ALPHA-BETA*(Y+sqrt(Y))+A1*Y**(1./3.)+A2*Y**(2./3.)+  &
       A3*Y+A4*Y**(4./3.)

end function FUNT
!======================================================================
!                        Function I(X) 
!	 function I(y) taken from Ejiri, JGR,1978
!======================================================================
real function FUNI(X)		! X is cos of equat pa = MU

  use ModHeidiSize
  use ModHeidiMain

  implicit none

  real :: Y,X,ALPHA,BETA,A1,A2,A3,A4
  !------------------------------------------------------------------------
  Y=sqrt(1.-X*X)
  ALPHA=1.+ALOG(2.+sqrt(3.))/2./sqrt(3.)
  BETA=ALPHA/2.-PI*sqrt(2.)/12.
  A1=0.055
  A2=-0.037
  A3=-0.074
  A4=0.056
  FUNI=2.*ALPHA*(1.-Y)+2.*BETA*Y*ALOG(Y)+4.*BETA*(Y-sqrt(Y))+  &
       3.*A1*(Y**(1./3.)-Y)+6.*A2*(Y**(2./3.)-Y)+6.*A4*(Y-Y**(4./3.))  &
       -2.*A3*Y*ALOG(Y)

end function FUNI
!============================================================================
!       new code to adjust the Cray's floating point operations====
!	Couldn't we just ROUND() it?
!============================================================================
real function APPX(T)

  implicit none
  real :: EPSLN, T
  !------------------------------------------------------------------------
  
  EPSLN=0.00000001
  if (T.lt.0.) then
     APPX=T-EPSLN
  else
     APPX=T+EPSLN
  end if
  if(int(APPX).ne.int(T)) return
  APPX=T

  return
end function APPX
!============================================================================
!                    end of new code     adjust floating
!        new code for missing functions on cray         begin    ab
!============================================================================           
real function ACOSD(X)

  use ModHeidiSize
  use ModHeidiMain
  use ModConst, only : cPi
  
  implicit none
  real :: x
  !------------------------------------------------------------------------
  
  acosd=180.0/cPi*acos(X)

end function ACOSD
!============================================================================          
real function ASIND(X)

  use ModHeidiSize
  use ModHeidiMain
  use ModConst, only : cPi
  
  implicit none
  real :: x
  !------------------------------------------------------------------------

  ASIND=180.0/cPi*asin(X)

end function ASIND

!============================================================================           
real function COSD(X)

  use ModHeidiSize
  use ModHeidiMain
  use ModConst, only : cPi
  
  implicit none
  real :: x
  !------------------------------------------------------------------------
  
  COSD=cos(cPi/180.0*X)

end function COSD

!============================================================================           
real function SIND(X)

  use ModHeidiSize
  use ModHeidiMain
  use ModConst, only : cPi
  
  implicit none
  real :: x
  !------------------------------------------------------------------------
  SIND=sin(cPi/180.0*X)

end function SIND


