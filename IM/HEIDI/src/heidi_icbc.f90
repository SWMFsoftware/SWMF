! File name: heidi_icbc.f90
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
!=======================================================================
!				INITIAL
!   	Initial set up of distribution functions (F2), energy (ENER)
!  	                and number of particle (N)
!=======================================================================

subroutine INITIAL(LNC,XN,J6,J18)

  use ModHeidiSize
  use ModHeidiIO
  use ModHeidiMain
  use ModIoUnit, only : io_unit_new, UNITTMP_
  use ModHeidiIO, only: NameInputDirectory
  use ModPlotFile, only: read_plot_file
  use ModHeidiInput, only:TypeFile

  implicit none

  integer :: j6,j18,i,j,k,l,ifn,jj,ii,kk,ig1,ig2,ig3,ig4,  &
       kg,kg1,kg2,IER,Iin,Kin
  real :: weight,elat1,elat2,etotal,efractn,pg,sg,y,yz,y10,y12,x,xlt
  real :: fbc,Cst1,Cst2,GAMMLN,esum
  real :: XN(NR,NS),LNC(NR,NS),N,FAC
  integer ::I1,J1,K1,L1
  parameter (I1=5, J1=9, K1=25, L1=19)	! From restart.bcf sizes
  !	PARAMETER (I1=5, J1=9, K1=22, L1=19)	! From restart.bcf sizes
  real :: FI(6,11),LI(6),EI(11),  &
       LIN(3),EIN(5),NI(3,5),E1(NE),E2(NE),F0(K1,L1),DUMMY,  &
       MU0(L1),E0(K1),MLT0(J1),R0(I1),G0(I1,J1)
  integer ::IFM(38)
  character(len=80) :: HEADER

  ! Common block for tests
  real :: XL1(4),XL2(2),XL3
  common /CLIN2/ XL1,XL2,XL3

  data LI/2.,3.,3.5,4.5,5.5,6.5/  ! for moder storm input
  data EI/0,0.25,0.5,0.75,1.,1.25,1.5,1.75,2.,2.25,2.5/
  !		EI => LOG(E[keV])
  data LIN/2.,4.25,6.5/           ! LIN and EIN => model PA shape
  data EIN/2.05,21.6,50.15,153.,350./	! EIN => E[keV]
  data IFM/2,7,13,20,28,35,42,47,50,52,54,56,58,60,62,64,66,68,70,  &
       2,11,21,31,41,51,61,70,75,79,82,83,84,85,86,87,88,89,90/
  external :: GAMMLN

  integer :: iUnit
  integer :: iUnit_tst
  character(len=30) :: NameOutputSpecies
  character(len=5)  :: NameSpecies
  character(len=5)  :: NameStormType
  character(len=3)  :: NamePrefix
  character(LEN=500):: StringVarName, StringHeader, NameFile
  real              :: acosd
  integer           :: ntc, IFIR

  real, allocatable :: Var_IIV(:,:,:) 
  integer           :: nDim, nParam,nVar       
  integer           :: n1, n2
  character(len=500):: NameVar
  real              :: Param_I(4)
  real              :: reshape
  character(len=20) :: TypePosition
  save ntc
  !-----------------------------------------------------------------------
  if (IFIR.eq.1) then
     ntc=int(nint(TIME/TINT))
  else
     ntc=ntc+1
  end if
  
  !.......Start the loop over RC species
  do S=1,NS

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
     if (SCALC(S).eq.1) then

        !.......Define the input file names
        if (ISTORM.eq.1) NameStormType='major'
        if (ISTORM.eq.2) NameStormType='moder'
        if (ISTORM.eq.3) NameStormType='tests'

        if (S.eq.1) NameSpecies='_e' 
	if (S.eq.2) NameSpecies='_h' 
	if (S.eq.3) NameSpecies='_he' 
	if (S.eq.4) NameSpecies='_o' 



        !.......Loss cone distribution (INI=1)
        if (INI(S).eq.1) then
           IBC=1
           if (Ab.gt.0.) then	! Maxwellian distribution
              do I=2,IO
                 do L=UPA(I),LO
                    do K=2,KO
                       do J=J6,J18
                          F2(I,J,K,L,S)=FBC(EKEV(K),FFACTOR(I,K,L),FINI(K)*CHI(I,J))
                       end do	! J loop
                    end do	! K loop
                 end do	! L loop
              end do	! I loop
           else			! Read in from input file: 'cone.bcf'

              open(UNITTMP_,FILE=NameInputDirectory//'cone.bcf',STATUS='OLD')
              read (UNITTMP_,101) HEADER
              do K=2,KO
                 read (UNITTMP_,*) DUMMY,FINI(K)
              end do	! K loop
              close(UNITTMP_)
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
           end if	! end of INI=1 case

           !.......Gaussian in R and PHI about some location (INI=2)
           !.......Distribution from input files (INI=3)
        else if ((INI(S).eq.2).or.(INI(S).eq.3)) then
           !.......Read in FI and NI from files
           NamePrefix='.in '
           open(UNITTMP_,file=NameStormType//trim(NameSpecies)//NamePrefix,STATUS='OLD')
           read(UNITTMP_,*) YEAR,DAY,R,AP,KP,ETOTAL,EFRACTN  !Replaces input.glo
           read(UNITTMP_,101) HEADER
           read(UNITTMP_,101) HEADER
           do  Kin=1,11
              read(UNITTMP_,*) (FI(Iin,Kin),Iin=1,6)
           end do
           close(UNITTMP_)
           NameStormType='testn'
           open(UNITTMP_,file=NameStormType//trim(NameSpecies)//NamePrefix,STATUS='OLD')
           read(UNITTMP_,101) HEADER
           read(UNITTMP_,101) HEADER
           do KK=1,5	
              read(UNITTMP_,*) (NI(II,KK),II=1,3)
           end do
           close(UNITTMP_)
           !.......Determine parameters for specific initial condition
           if (INI(S).eq.2) then
              ig1=4			! Gausian centered at L=LZ(ig1)
              ig2=ig1
              ig3=ig1-2
              ig4=ig1+2
              kg=24			! Gaussian centered at E=EKEV(kg)
              kg1=kg-7
              kg2=kg+5
              do k=kg1,kg2
                 E1(k)=ALOG10(EKEV(kg))
                 E2(k)=10*exp(-(EKEV(k)-EKEV(kg))**2/15.)
              end do
              pg=3.0		! Gaussian extends to |PHI|<pg
              sg=0.05		! Variance of PHI Gaussian
           else
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
           end if
           !.......Find distribution at local midnight
           do I=ig1,ig2
              do K=kg2,kg1,-1
                 if (EKEV(K).gt.1.) then
                    call LINTP2(LI,EI,FI,6,11,LZ(I),E1(K),Y,IER)
                    call LINTP2(LIN,EIN,NI,3,5,LZ(I),10**E1(K),YZ,IER)
                    Y10=(10**Y)*E2(k)
                    do L=2,UPA(I)-1
                       F2(I,1,K,L,S)=Y10*(1.-MU(L)**2)**(YZ/2.)*FFACTOR(I,K,L)
                    end do	! L loop
                    KK=K
                 else			! Maxwellian below 1 keV
                    X=EKEV(k)/EKEV(KK)
                    Y12=Y10*X*exp(1.-X)
                    do L=2,UPA(I)-1
                       F2(I,1,K,L,S)=Y12*(1.-MU(L)**2)**(YZ/2.)*FFACTOR(I,K,L)
                    end do ! L loop
                 end if
              end do	! K loop
           end do	! I loop

           !.......Gaussian in R for INI=2
           if (INI(S).eq.2) then
              do l=2,UPA(i)
                 do k=kg1,kg2
                    do i=ig3,ig4
                       F2(i,1,k,l,s)=F2(ig1,1,k,l,s)*exp(-(LZ(i)-LZ(ig1))**2/0.005)
                    end do
                 end do
              end do
           end if
           !.......Gaussian in PHI
           do J=1,JO
              XLT=MLT(J)
              if (XLT.gt.12) XLT=XLT-24.
              XLT=abs(XLT)
              do i=ig3,ig4
                 do l=2,UPA(i)-1
                    do k=kg1,kg2
                       if (XLT.lt.pg) then
                          F2(i,j,k,l,s)=F2(i,1,k,l,s)*exp(-XLT**2/sg)
                       end if
                    end do	! K loop
                 end do		! L loop
              end do		! I loop
           end do 		! J loop

           !.......MICS quiet RC, constant with R, PHI, and MU (INI=4)
        else if (INI(S).eq.4) then
           JJ=0
           K=1
           do while (JJ.eq.0) 
              if (EKEV(K).gt.Eob) then
                 JJ=K
              else if (K.eq.KO) then
                 JJ=K+1
              end if
              K=K+1
           end do
           II=0
           K=1
           do while (II.eq.0) 
              if (EKEV(K).gt.40.) then
                 II=K
              else if (K.eq.KO) then
                 II=K+1
              end if
              K=K+1
           end do
           FAC=1.
           do L=1,LO
              do K=2,JJ
                 if (IFAC.ne.1) FAC=1./FLUXFACT(S)/EKEV(K)
                 DUMMY=(ALOG(EKEV(K))-ALOG(1.))/(ALOG(Eob)-ALOG(1.))
                 WEIGHT=exp(DUMMY*ALOG(Ab)+(1.-DUMMY)*ALOG(0.1*Ab))
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
                 if (IFAC.ne.1) FAC=1./FLUXFACT(S)/EKEV(K)
                 DUMMY=(ALOG(EKEV(K))-ALOG(Eob))/(ALOG(40.)-ALOG(Eob))
                 WEIGHT=exp(DUMMY*ALOG(0.1*Ab)+(1.-DUMMY)*ALOG(Ab))
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
                 if (IFAC.ne.1) FAC=1./FLUXFACT(S)/EKEV(K)
                 DUMMY=(ALOG(EKEV(K))-ALOG(40.))/(ALOG(300.)-ALOG(40.))
                 WEIGHT=exp(DUMMY*ALOG(Ab)+(1.-DUMMY)*ALOG(0.1*Ab))
                 do J=1,JO
                    do I=2,6
                       F2(I,J,K,L,S)=WEIGHT*FAC*FFACTOR(I,K,L)*(.5)**(6-I)
                    end do		! I loop
                    do I=7,10
                       F2(I,J,K,L,S)=WEIGHT*FAC*FFACTOR(I,K,L)
                    end do		! I loop
                 end do		! J loop
                 WEIGHT=exp(DUMMY*ALOG(0.1*Ab)+(1.-DUMMY)*ALOG(Ab))
                 do J=1,JO
                    do I=11,IO
                       F2(I,J,K,L,S)=WEIGHT*FAC*FFACTOR(I,K,L)
                    end do		! I loop
                 end do		! J loop
              end do		! K Loop
           end do		! L loop

           !.......Read in F from a file: 'restart.bcf' (INI=5)
        else if (INI(S).eq.5) then
           iUnit = io_unit_new()
           open(unit=iUnit,FILE='restart.bcf',STATUS='OLD')
           IFN=0
           if (IPA.eq.0) IFN=L1
           do L=1,L1
              MU0(L)=MU(IFM(L+IFN))
           end do
           read (iUnit,101) HEADER
           read (iUnit,101) HEADER
           do I=1,I1		! Read in and perform first 2D interp.
              II=I*4-2
              do J=1,J1-1
                 JJ=J*3-2
                 read (iUnit,102) R0(I),MLT0(J)
                 read (iUnit,101) HEADER
                 do K=1,K1
                    read (iUnit,*) E0(K),(F0(K,L),L=1,L1)
                 end do
                 !	  IF (I+J.EQ.2) THEN
                 !	   PRINT 50,'Inputs :',1,2,2,18,E0(1),E0(2),MU0(2),MU0(18)
                 !	   PRINT 50,'My grid:',2,5,7,68,EKEV(2),EKEV(5),MU(7),MU(68)
                 !	  END IF
                 !	  PRINT 50, 'F0:',I,II,J,JJ,F0(1,2),F0(1,18),F0(2,2),F0(2,18)
50               format(A,4I4,1P,4E12.4)
                 do L=1,L1		! Convert to F2 and log scale
                    do K=1,K1
                       if (F0(K,L) .lt. 1.E-30) then
                          F0(K,L)=1.E-30
                       end if
                       F0(K,L)=ALOG10(F0(K,L))
                    end do	! K loop
                 end do	! L loop
                 do K=2,KO
                    do L=1,LO-1
                       call LINTP2(E0,MU0,F0,K1,L1,EKEV(K),MU(L),F2(II,JJ,K,L,S),IER)
                       if (EKEV(K).gt.E0(K1)) &
                            F2(II,JJ,K,L,S)=AMIN1(F2(II,JJ,K,L,S),F2(II,JJ,K-1,L,S))
                    end do	! L loop
                 end do	! K loop
                 !	  PRINT 50, 'F2:',I,II,J,JJ,10**F2(II,JJ,2,7,S),10**F2(II,JJ,2,68,S),  &
                 !     		10**F2(II,JJ,5,7,S),10**F2(II,JJ,5,68,S)
                 !	  PRINT 51, 'LIN:',(XL1(K),K=1,4),(XL2(L),L=1,2),XL3
51               format (A,1P,7E10.2)
                 !	  IF (I+J.EQ.2) THEN
                 !	   PRINT 50,'Inputs :',1,2,2,18,E0(1),E0(2),MU0(2),MU0(18)
                 !	   PRINT 50,'My grid:',2,5,7,68,EKEV(2),EKEV(5),MU(7),MU(68)
                 !	   PRINT 50, 'F0#2:',I,II,J,JJ,10**F0(1,2),10**F0(1,18),10**F0(2,2),  &
                 !     		10**F0(2,18)
                 !	   PRINT 50,'F2in:',2,3,4,5,(10**F2(II,JJ,K,7,S),K=2,5)
                 !	  END IF
              end do	! J loop
           end do		! I Loop
           close(iUnit)
           !	STOP			!TEST RUNS ONLY
           MLT0(J1)=MLT(1)
           do K=2,KO		! Perform second 2D interpolation
              do L=1,LO-1
                 do J=1,J1
                    JJ=J*3-2
                    if (J.eq.J1) then
                       JJ=1
                    end if
                    do I=1,I1
                       II=I*4-2
                       G0(I,J)=F2(II,JJ,K,L,S)
                    end do	! I loop
                 end do	! J loop
                 do J=1,JO
                    do I=2,IO
                       call LINTP2(R0,MLT0,G0,I1,J1,LZ(I),MLT(J),F2(I,J,K,L,S),IER)
                       if (LZ(I).gt.R0(I1)) &
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
                       if (F2(I,J,K,L,S).le.1.E-30*FFACTOR(I,K,L)) then
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
        else if (INI(S).eq.6) then
           IBC=6
           if (S.eq.1) then
              Einj=.2			! Characteristic E, keV
              Kinj=6.			! Kappa value
              Ninj=.1			! Density, cm-3
           else 
              Einj=1.4			! Characteristic E, keV
              Kinj=5.5			! Kappa value
              Ninj=.4			! Density, cm-3
           end if
           Cst1=sqrt(Q*1.E4/(2.*MAS(S)*(PI*Kinj*Einj)**3))
           Cst2=exp(GAMMLN(Kinj+1.,IER)-GAMMLN(Kinj-0.5,IER))
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
5          format (1P,5E12.4)

           !  Read in from a unformatted file (INI=7)
!!$        else if (INI(S).eq.7) then
!!$           open(UNIT=UnitTmp_,FILE=NameRestartInDir//trim(NameRun)//trim(NameSpecies)//'.unff',status='old',   &
!!$               form='unformatted')
!!$           do L=1,NPA
!!$              do K=1,NE  ! Change back to this for restarts
!!$                 do J=1,NT 
!!$                    read(UnitTmp_) (f2(I,J,K,L,S),I=1,NR)
!!$                 end do
!!$              end do
!!$           end do
!!$           close(UnitTmp_)
!!$
!!$        end if

           !  Read in from a restart file (INI=7) 
        else if (INI(S).eq.7) then
           NameFile       = trim(NameRestartInDir)//'restart'//trim(NameSpecies)//'.out'
           StringHeader   = &
                'Phase space distribution function for all pitch angles, energies and locations.'
           
           StringVarName ='R   MLT   F   E   PA'
           iUnit = io_unit_new()
           n1 = NR
           n2 = NT

           do L=1,NPA 
              do K=1,NE
                 ! Read header info 
                 call read_plot_file(NameFile,&
                      iUnitIn           = iUnit,&
                      TypeFileIn        = TypeFile, &
                      StringHeaderOut   = StringHeader,&
                      nStepOut          = ntc, &
                      TimeOut           = T,&
                      nDimOut           = nDim, &
                      nParamOut         = nParam, &
                      nVarOut           = nVar,&
                      n1Out             = n1,&
                      n2Out             = n2, &
                      ParamOut_I        = Param_I)
                                
                 ! Determine the shape of arrays from the header
                 !allocate(Coord_DII(2,NR,NT),Var_IIV(NR,NT,1))
                 
                 allocate(Var_IIV(NR,NT,1))
                 ! Read the coord and var arrays
                 call read_plot_file(NameFile,&
                      iUnitIn       = iUnit,&
                      TypeFileIn    = TypeFile, &
                      !CoordOut_DII  = Coord_DII,&
                      VarOut_IIV    = Var_IIV)
                 f2(:,:,K,L,S:S) = Var_IIV

           deallocate(Var_IIV)  
        end do
     end do


end if
        !.......Done with initial particle distribution set up

        if (IBC(S).eq.1) then	! This is to redo INI=1 above
           Ib=1
           open(unit=iUnit,FILE='cone.bcf',STATUS='OLD')
           read (iUnit,101) HEADER
           do K=2,KO
              read (iUnit,*) DUMMY,FINI(K)
           end do
           close(iUnit)
           do J=J6,J18
              do I=2,IO
                 CHI(I,J)=1.		! Same BC flux over entire dayside
              end do
           end do
        end if

        !.......Calculate the total # of particles and energy of this species
        N=0				! total # dens of RC specie "s"
        ESUM=0				! total E of RC for specie "s"
        !comment out when crashes from here

!!$        do I=1,IO
!!$           ENER(I,S)=0			! E of RC specie for some LZ
!!$        end do
!!$        do I=2,IO
!!$           do J=2,JO
!!$              do K=2,KO
!!$                 do L=2,UPA(I)-1
!!$                    WEIGHT=F2(I,J,K,L,S)*WE(K)*WMU(L) ! F2 - average in cell
!!$                    XN(I,S)=XN(I,S)+WEIGHT		 !      (i,j,k,l)
!!$                    ENER(I,S)=ENER(I,S)+EKEV(K)*WEIGHT  ! for some LZ
!!$                 end do	        ! L loop
!!$              end do		! K loop
!!$           end do 		! J loop
!!$           ESUM=ESUM+ENER(I,S)
!!$           N=N+XN(I,S)	
!!$        end do

        ! to here
        !.......FACTOR, a scaling for ESUM and N, depends on IFAC:
        if (IFAC.eq.1) FACTOR(S)=8.6474E13/M1(S)**1.5*DR*DPHI
        if (IFAC.eq.2) FACTOR(S)=4.3237E13/M1(S)**1.5*DR*DPHI

        !.......Calculate the characteristics (mean energy) at T=0
        NameStormType='rnsc1'
        if(T/3600.eq.48) NameStormType='rnsc3'
        if(T/3600.eq.96) NameStormType='wpa10'
        !	OPEN(1,FILE=NameStormType//NameSpecies//'.l')
        !	WRITE(1,*)' Losses due to some processes'
        !	WRITE(1,15) KP
15      format(2X,5HT = 0,2X,4HKp =,F6.2,/,6X,1HL,10X,11HENIGHT[keV],2X,   &
             9HEDAY[keV],2X,6HNNIGHT,2X,4HNDAY,4X,6HMean E)
        !	DO I=2,IO
        !	 EMEAN=ENER(I,S)/XN(I,S)
        !	 WRITE(1,17) LZ(I),ENER(I,S)*FACTOR(S),XN(I,S)*FACTOR(S),EMEAN
        !	END DO
17      format(2X,F7.2,5(2X,1PE12.5))
        !	AMEAN=ESUM/N
        !	WRITE(1,19)ESUM*FACTOR(S),N*FACTOR(S),AMEAN
19      format(/,4X,5HTotal,3(2X,1PE12.5))
        !	CLOSE(1)

        !.......Initial loss is zero
        LNC(1:io,S)=0       ! Loss of particles due to all processes
        LEC(1:io,S)=0       ! Loss of energy due to     "       "

     end if		! SCALC Check
  end do		! S loop

  !.......Calculate the energy coeff
  X=RE+HMIN
  do I=1,IO
     ELAT1=sqrt(X/(LZ(I)-DL1/2.)/RE)
     ELAT2=sqrt(X/(LZ(I)+DL1/2.)/RE)
     ECOF(I)=DPHI*X**2*(acos(ELAT2)-acos(ELAT1))*(ELAT1+ELAT2)/2.
  end do	! I loop

101 format(A80)
102 format(21X,F6.2,5X,F4.1)
103 format(F7.3,20(1PE9.2))

end subroutine INITIAL
!=======================================================================
!				 LMPLOSS
!   		If LMP moved inward, then lose the particles 
!=======================================================================
subroutine LMPLOSS

  use ModHeidiSize
  use ModHeidiIO
  use ModHeidiMain
  use ModIoUnit, only : io_unit_new

  implicit none

  real :: FLO
  integer :: I,K,L,J
  !---------------------------------------------------------------------
  do J=1,JO
     if (ILMP(J).lt.ILold(J)) then
        do L=1,LO		! Lose everything beyond magnetopause
	   do K=2,KO
              do I=ILMP(J)+1,ILold(J)
                 FLO=1.E-30*FFACTOR(I,K,L)
                 RNL=RNL+(F2(I,J,K,L,S)-FLO)*CONSL(K,S)*WE(K)*WMU(L)*DR*DPHI
                 REL=REL+(F2(I,J,K,L,S)-FLO)*CONSL(K,S)*EKEV(K)   &
                      *WE(K)*WMU(L)*DR*DPHI
                 F2(I,J,K,L,S)=FLO
              end do
	   end do
        end do
     end if
  end do

end subroutine LMPLOSS

!=======================================================================
!				 GEOSB
!   			Boundary conditions set up 
!=======================================================================

subroutine GEOSB

  use ModHeidiSize
  use ModHeidiIO
  use ModHeidiMain
  use ModHeidiDrifts
  use ModIonoHeidi
  use ModNumConst, only: cDegToRad
  use ModIoUnit, only : io_unit_new

  implicit none

  character(len=80) :: HEADER
  real :: Fkapb(NE),Foutb,X,data(9),fac
  real :: TM1,TM2,NM1,NM2,TFM1,TFM2,TCM1,TCM2,NM,TFM,TCM,TS1,TS2,  &
       Flanl(NE,NPA,NS),FS1(7),FS2(7),ES(7),FS(7),NY(NS),NEL,  &
       NE1,NE2,TEF1,TEF2,TEC1,TEC2,NMFAC(NS),FSFAC(NS),TEF,TEC,  &
       Ekap,Kappa,GAMMLN
  integer :: I,J,K,L,IER,GetBCI,BCI,I2,KES(NS),I6,I7,I9,IG7
  integer :: SKAPPA(4)
  external :: GetBCI,GAMMLN
  data ES/45.,60.,94.,141.,210.,325.,535./, Kappa/5.00/
  data SKAPPA/1,1,1,1/
  save KES,TM1,TM2,NM1,NM2,TFM1,TFM2,TCM1,TCM2,TS1,TS2,FS1,FS2,  &
       I2,I6,I7,I9,IG7,NE1,NE2,TEF1,TEF2,TEC1,TEC2

  integer :: iLatBoundary=-1, iLonBoundary=-1
  !---------------------------------------------------------------------

  call write_prefix; write(iUnitStdOut,*) 'Resetting the outer boundary condition'

  !CCC Create a few flags and open a few files
  if (T.eq.TIME) then
     I6=0
     I7=0
     IG7=0
     I9=0
     do S=1,NS
        if (IBC(S).eq.6) I6=1
        if (IBC(S).eq.7) I7=1
        if (IBC(S).eq.9) I9=1
        if (IBC(S).gt.7) IG7=1
     end do
     if (I7.eq.1 .or. IG7.eq.1) then
        do S=1,NS
           if (S.eq.1 .or. IBC(S).eq.7) then    ! no SOPA data
              KES(S)=KO
           else
              KES(S)=0
              K=0
              do while (KES(S).eq.0)
                 K=K+1
                 if (EKEV(K).ge.ES(1)) then
                    KES(S)=K-1
                 else
                    if (K.eq.KO) KES(S)=KO
                 end if
              end do
           end if
        end do
        if (IG7.eq.1) then
           TS2=TIME-1.		! Prepare SOPA input file
           TS1=TS2
           FS2(1:7)=0.
           iUnitSopa = io_unit_new()
           open(UNIT=iUnitSopa,FILE=NameInputDirectory//trim(NameRun)//'_sopa.in',status='old')
           do I=1,3
              read(iUnitSopa,*) HEADER
           end do
        end if
        TM2=TIME-1.		! Prepare MPA input file
        TM1=TM2
        NM2=0.
        TFM2=0.
        TCM2=0.
        NE2=0.
        TEC2=0.
        TEF2=0.
        iUnitMpa = io_unit_new()
        open(UNIT=iUnitMpa,FILE=NameInputDirectory//trim(NameRun)//'_mpa.in',status='old')
        do I=1,3			! 3 lines of header material
           read(iUnitMpa,*) HEADER
        end do
        I2=0
        if (S.eq.1) I2=3
     end if
  end if

  do S=1,NS

     do L=1,LO
        do K=1,KO
           do J=1,JO
              FGEOS(J,K,L,S)=0. 
           enddo
        enddo
     enddo

     if (SCALC(S).eq.1) then

        if (TINJ.gt.TIME+2.*DT*NSTEP) then ! No injection, use IC for BC
           do L=1,LO
              do K=1,KO
                 do J=1,JO
                    FGEOS(J,K,L,S)=F2(IO,J,K,L,S)
                 end do
              end do
           end do
        end if

     end if   ! SCALC check
  end do	 ! S loop

  if (I7.eq.1 .or. IG7.eq.1) then	! LANL data injection
     if (TM2.lt.T) then			! MPA DATA
        do while (TM2.le.T)		! Best if final TM2 > final T
           TM1=TM2
           NM1=NM2
           TFM1=TFM2
           TCM1=TCM2
           read (iUnitMpa,*,IOSTAT=L) (data(I),I=1,9)
           TM2=data(2)
           NM2=data(4)
           TFM2=data(6)
           TCM2=data(5)
           NE2=data(7)
           TEF2=data(9)
           TEC2=data(8)
           if (L.lt.0) TM2=TIME+2*DT*(NSTEP+1)
           if (T.eq.TIME) then		! In case T2>T already
              TM1=TIME
              NM1=NM2
              TFM1=TFM2
              TCM1=TCM2
              NE1=NE2
              TEC1=TEC2
              TEF1=TEF2
           end if
        end do
     end if

     ! inputs are in /cc and eV

     FAC=(T-TM1)/(TM2-TM1)		! Linearly interpolate
     NM=FAC*NM2+(1.-FAC)*NM1		! in cm-3
     TFM=(FAC*TFM2+(1.-FAC)*TFM1)*1.E-3	! in keV
     TCM=(FAC*TCM2+(1.-FAC)*TCM1)*1.E-3	! in keV
     NEL=FAC*NE2+(1.-FAC)*NE1		! in cm-3
     TEF=(FAC*TEF2+(1.-FAC)*TEF1)*1.E-3	! in keV
     TEC=(FAC*TEC2+(1.-FAC)*TEC1)*1.E-3	! in keV

     !-------------------
     ! DONE reading in data

     ! This is a total hack - we need to figure out how to NOT read in the
     ! data above.

     if (maxval(IonoGmDensity) > 0.0) then

        !write(*,*) "---------------------------------------------------"
        !write(*,*) "Ignoring LANL data, and overwriting with GM Data!!!"

        ! Find location to take boundary condition from : 67 degrees
        iLonBoundary = IONO_nPsi/2.0
        !        if (iLatBoundary < 0) then
        iLatBoundary = 1
        do while (IONO_NORTH_Theta(iLatBoundary,1) < (90.0-67.0)*cDegToRad .and. &
             IonoGmDensity(iLatBoundary, iLonBoundary) == 0.0)
           iLatBoundary = iLatBoundary + 1
        enddo
        !        endif

	!write(*,*) "Taking boundary condition from location : ", &
        !     iLatBoundary, iLonBoundary

        NM = IonoGmDensity(iLatBoundary, iLonBoundary)
        TCM = IonoGmTemperature(iLatBoundary, iLonBoundary)
        TFM = TCM

        !write(*,*) "Density : ", NM, " /cc"
        !write(*,*) "Temperature : ", TCM, " eV"

        !write(*,*) "---------------------------------------------------"

     endif

     NY(2)=0.34*exp(0.054*KP)	! From Young et al 1982
     NY(3)=5.1E-3*exp(6.6E-3*F107)
     NY(4)=0.011*exp(0.24*KP+0.011*F107)
     FAC=0.
     do I=2,4 			! Weights corrects MPA moments
        FAC=FAC+NY(I)/sqrt(M1(I))
     end do
     NMFAC(1)=1.
     do S=2,NS
        NMFAC(S)=NY(S)/FAC
        !CCC The next line is for the idealized test runs only, no BC variation:
        !CCC	    NMFAC(S)=1.
     end do
     if (IG7.eq.1) then		! SOPA DATA
        if (TS2.lt.T) then
           do while (TS2.le.T)	! Best if final TS2 > final T
              TS1=TS2
              FS1(2:7)=FS2(2:7)
              read (iUnitSopa,*,IOSTAT=I) (data(I),I=1,8)
              TS2=data(1)
              FS2(2:7)=data(3:8)
              if (I.lt.0) TS2=TIME+2*DT*(NSTEP+1)
              if (T.eq.TIME) then		! In case T2>T already
                 TS1=TIME
                 FS1(2:7)=FS2(2:7)
              end if
           end do
        end if
        FAC=(T-TS1)/(TS2-TS1)		! Linearly interpolate
        FS(2:7)=FAC*FS2(2:7)+(1.-FAC)*FS1(2:7)	! flux units
     end if
     FAC=(NY(2)+NY(3)+NY(4))  
     FSFAC(1)=1.
     do S=2,NS
        FSFAC(S)=NY(S)/FAC
        !CCC The next line sometimes gives bad results...use with care.
        if (NEL*FSFAC(S).gt.NM*NMFAC(S)) NMFAC(S)=NEL*FSFAC(S)/NM
     end do
     do S=1,NS
        if (SCALC(S).eq.1) then
           !CCC _hsopa: H+ has a different function than all other species
           if (SKAPPA(S).eq.1) then  ! bi-kappa=5,  SOPA (if designated)
              Ekap=TFM*(Kappa-1.5)/Kappa
              FAC=(TFM/TCM)*(MAS(S)*1.E13/(2.*PI*Kappa*Ekap*Q))**1.5
              FAC=FAC*exp(GAMMLN(Kappa+1.,IER)-GAMMLN(Kappa-0.5,IER))
              do L=1,LO
                 do K=1,KES(S)          ! MPA moments
                    Flanl(K,L,S)=NM*NMFAC(S)*FAC*(1.+EKEV(K)*(MU(L)**2+   &
                         (TFM/TCM)*(1.-MU(L)**2))/(Kappa*Ekap))**(-Kappa-1.)
                 end do
              end do
              if (IG7.eq.1) then	
                 FS(1)=AMAX1(.5*FS(2),Flanl(KES(S),10,S)   &
                      *FLUXFACT(S)*EKEV(KES(S)))
                 do K=KES(S)+1,KO                        ! SOPA fluxes
                    call LINTP(ES,FS,7,EKEV(K),FAC,IER)
                    Flanl(K,1:LO,S)=FAC/FLUXFACT(S)/EKEV(K)
                 end do
              endif

           else  ! Maxwellian everywhere, SOPA (if designated)
              do L=1,LO
                 do K=2,KES(S)			! MPA moments
                    ! Converts Density & Temperature Moments to Fluxes
                    ! TCM = Temperature (par or perp)
                    ! TFM = Temperature (perp or par) - can have same
                    ! NM  = density
                    Flanl(K,L,S)=NM*NMFAC(S)*(MAS(S)*1.E13/Q/2./PI)**1.5/  &
                         sqrt(TFM)/TCM*exp(-EKEV(K)*((1.-MU(L)*MU(L))/TCM    &
                         +MU(L)*MU(L)/TFM))
                 end do  ! K loop for MPA
              end do   ! L loop for MPA

              if (IG7.eq.1) then
                 FS(1)=AMAX1(.5*FS(2)*FSFAC(S),Flanl(KES(S),10,S)   &
                      *FLUXFACT(S)*EKEV(KES(S)))/FSFAC(S)
                 do K=KES(S)+1,KO			! SOPA fluxes
                    call LINTP(ES,FS,7,EKEV(K),FAC,IER)
                    Flanl(K,1:LO,S)=FSFAC(S)*FAC/FLUXFACT(S)/EKEV(K)
                 end do ! K loop for SOPA
              endif
           end if  ! Block for S=2 and all others
        end if  ! SCALC check
     end do   ! S loop
  else 
     return
  end if

  !...Find injection boundary and fill in (Note: IBC=9, inject everywhere)

  ! FGEOS is outer boundary condition on plasma
  ! j - azimuth (mlt)
  ! k - energy
  ! l - pitch angle
  ! s - species (1=electrons, 2=H+, 3=He+, 4=O+)

  do S=1,NS
     if (SCALC(S).eq.1) then
        if (IBC(S).eq.6) then
           call FINJ(Fkapb)
           do L=1,LO
              do K=2,KO
                 do J=1,JO
                    FGEOS(J,K,L,S)=Fkapb(K)*FFACTOR(IO,K,L)
                 end do	! J LOOP
              end do	! K LOOP
           end do	! L LOOP
        else if (IBC(s).ge.7) then
           do L=1,LO
              do K=2,KO
                 do J=1,JO
                    FGEOS(J,K,L,S)=Flanl(K,L,S)*FFACTOR(IO,K,L)
                 end do	! J LOOP
              end do	! K LOOP
           end do	! L LOOP
        end if
     end if		! SCALC check
  end do		! S loop

end subroutine GEOSB
!=======================================================================
!				FBC
!	This function calculates the dayside loss cone boundary flux
!	It uses either a Maxwellian or a distribution from a file
!=======================================================================
real function FBC(E,FAC,F1)

  use ModHeidiSize
  use ModHeidiIO

  implicit none

  real :: e,fac,f1
  !-----------------------------------------------------------------------
  if (Ab.gt.0.) then
     FBC=Ab*exp(-E/Eob)*FAC
  else
     FBC=F1*FAC
  end if
  return
end function FBC
!=======================================================================
!			     FINJ
!     Calculates the plasma sheet distribution for the boundary
!=======================================================================

subroutine FINJ(F)

  use ModHeidiSize
  use ModHeidiIO
  use ModHeidiMain
  use ModIoUnit, only : io_unit_new

  implicit none

  integer :: K,IER,I
  real    :: F(NE),Cst1,Cst2,GAMMLN,CONV
  real    :: T1,T2,BZ1,BZ2,MD1,MD2,U1,U2,FAC,ERF,NY(NS),TLAG
  character header*80
  save T1,T2,BZ1,BZ2,MD1,MD2,U1,U2
  external :: GAMMLN,ERF

  !-----------------------------------------------------------------------

  TLAG=4.*3600.			! From Borovsky et al, Aug 98
  if (ISWB.eq.1) then
     if (T.eq.TIME) then
        T2=TIME-1.
        T1=T2
        iUnitSw2 = io_unit_new()
        open(UNIT=iUnitSw2,FILE=NameInputDirectory//trim(NameRun)//'_sw2.in',status='old')
        do I=1,6                      ! 6 lines of header material
           read(iUnitSw2,*) HEADER
        end do
     end if
     if (T2.lt.T) then
        do while (T2.lt.T)    ! Best if final T2 > final T
           T1=T2
           BZ1=BZ2
           MD1=MD2
           U1=U2
           read (iUnitSw2,*,IOSTAT=I) T2,BZ2,MD2,U2
           T2=T2+TLAG
           if (I.lt.0) T2=TIME+2*DT*(NSTEP+1)+TLAG
           if (T.eq.TIME) then                 ! In case T2>T already
              T1=TIME
              BZ1=BZ2
              MD1=MD2
              U1=U2
           end if
        end do
     end if
     FAC=(T-T1)/(T2-T1)                      ! Linearly interpolate
     NSWB=(FAC*MD2+(1.-FAC)*MD1)	        ! in cm-3
     USWB=(FAC*U2+(1.-FAC)*U1)                ! in km/s
     Einj=2.17+0.0223*USWB		! From Borovsky et al 1998
     Ninj=0.292*NSWB**0.49
  else
     NSWB=0.
     USWB=0.
     Ninj=.6			! Average ion values
     Einj=10.
  end if
  if (S.eq.1) then
     Einj=Einj/7.			! From Huang and Frank 1986
     Kinj=6.
     Ninj=Ninj*.5			! From various sources
  else
     !	  Kinj=5.*SQRT(Einj)		! From Huang and Frank 1986
     Kinj=8.-5.*ERF(Einj/23.,IER)	! Gets harder with Einj
     NY(2)=0.34*exp(0.054*KP)	! From Young et al 1982
     NY(3)=5.1E-3*exp(6.6E-3*F107)
     NY(4)=0.011*exp(0.24*KP+0.011*F107)
     NY(1)=0.
     do I=2,4 			! Weights corrects LANL data
        NY(1)=NY(1)+NY(I)/sqrt(M1(I))
     end do
     Ninj=Ninj*NY(S)/NY(1)
  end if
  Cst1=(MAS(S)*1.E13/(2.*PI*Kinj*Einj*Q))**1.5
  Cst2=exp(GAMMLN(Kinj+1.,IER)-GAMMLN(Kinj-0.5,IER))
  CONV=1.
  do K=1,KO
     if (IFAC.eq.1) CONV=FLUXFACT(S)*EKEV(K)
     F(K)=Ninj*CONV*Cst1*Cst2*(1.+EKEV(K)/(Kinj*Einj))**(-Kinj-1.)
  end do
end subroutine FINJ



