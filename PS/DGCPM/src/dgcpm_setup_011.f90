! File name: dgcpm_setup_011.f90
!
! Contains: input and array setup routines for DGCPM
!	READPARA
!	CONSTANT
!	ARRAYS
!	GETKPA
!	THERMAL
!
! Change from 010 to 011: Dipole field definition now a subroutine
!	instead of embedded in ARRAYS
!
! Last Modified: December 2006, Mike Liemohn
!
! **********************************************************************
!                             READPARA
!	Read parameters: DT, TMAX, Species, Storm, .......
! **********************************************************************
subroutine readpara

  use ModSizeDGCPM
  use ModIoDGCPM
  use ModMainDGCPM

  implicit none

  real tmax
  integer k,i
  character*1 header

!  OPEN (UNIT=1,FILE='input.dgcpm',STATUS='OLD')
!  READ (1,*) DT,TMAX,TINT,TIME
!  READ (1,*) IO,JO,KO,LO,ISO
!  READ (1,*) ELB,SWE,RW,HMIN
!  READ (1,*) ISTORM,IKP,IPA,IFAC,IST,IWPI,ISW,IA,ITHERMINIT
!  READ (1,*) (SCALC(k),k=1,NS)
!  READ (1,*) YEAR,DAY,UT,R,AP,KP
!  READ (1,*) (INI(k),k=1,NS)
!  READ (1,*) (IBC(k),k=1,NS)
!  READ (1,*) TINJ,Ab,Eob
!  READ (1,*) (IRES(k),k=1,15)
!  READ (1,*) NAME

  i = 1
     
  write(*,*) "Setting Parameters"

  READ (cInputText(i),*) DT,TMAX,TINT,TIME
  i=i+1
  READ (cInputText(i),*) IO,JO,KO,LO,ISO
  i=i+1
  READ (cInputText(i),*) ELB,SWE,RW,HMIN
  i=i+1
  READ (cInputText(i),*) ISTORM,IKP,IPA,IFAC,IST,IWPI,ISW,IA,ITHERMINIT
  i=i+1
  READ (cInputText(i),*) (SCALC(k),k=1,NS)
  i=i+1
  READ (cInputText(i),*) YEAR,DAY,UT,R,AP,KP
  i=i+1
  READ (cInputText(i),*) (INI(k),k=1,NS)
  i=i+1
  READ (cInputText(i),*) (IBC(k),k=1,NS)
  i=i+1
  READ (cInputText(i),*) TINJ,Ab,Eob
  i=i+1
  READ (cInputText(i),*) (IRES(k),k=1,15)
  i=i+1
  READ (cInputText(i),*) NAME
  i=i+1

!  CLOSE(1)

  NSTEP=NINT(TMAX/DT/2.)                 ! time splitting

  ithermfirst=1		! So we do the setup routines in THERMAL

  LAMGAM=2.       ! Empirical E-field shielding parameter

  RETURN

end subroutine readpara
!
! End of subroutine READPARA
!

!***********************************************************************
!				CONSTANT
!			    Define constants
!***********************************************************************
SUBROUTINE CONSTANT

  use ModSizeDGCPM
  use ModIoDGCPM
  use ModMainDGCPM

  implicit none

  real kpt,DUT
  integer I
  character*80 header

  !.......Read Kp history of the modeled storm
  IF (IKP.GE.3) THEN
     OPEN(1,FILE=cInputDir//NAME//'_kp.in',STATUS='OLD') 
     READ(1,10) HEADER
10   FORMAT(A80)
     DO I=1,NSTEP/NKP+2
        READ(1,*) DAYR(I),DUT,RKPH(I),F107R(I),APR(I),RSUNR(I)
     ENDDO
     CLOSE(1)
  END IF

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

  RETURN
END SUBROUTINE CONSTANT
!
! End of subroutine CONSTANT
!

!***********************************************************************
!				BFIELD_SETUP
!			Set up all the arrays
!***********************************************************************
SUBROUTINE BFIELD_SETUP(cone)

  use ModSizeDGCPM
  use ModMainDGCPM
  use ModDgcpmWaves
  use ModIoDGCPM

  Implicit None

  REAL CONE(NR+4),degrad,camlra
  integer i,iml

  degrad=pi/180.
  do i=1,IO
     DO IML=1,ISO
        camlra=amla(iml)*degrad
        BE(I,IML)=0.32/LZ(I)**3*SQRT(1.+3.*SIN(camlra)**2)  &
             /COS(camlra)**6				! in gauss
     ENDDO
     BFC(I)=ME/Z(I)**3/40./SQRT(PI*Q)		! in SI units
  END DO

  !.......CONE - pitch angle loss cone in degree
  do i=1,io	
     CONE(I)=ASIN(SQRT(((RE+HMIN)/Z(I))**3   &
          /SQRT(4.-3.*((RE+HMIN)/Z(I)))))/degrad
  end do
  CONE(IO+1)=2.5    ! to calcul PA grid near 0 deg for IPA=1
  CONE(IO+2)=1.5
  CONE(IO+3)=1.
  CONE(IO+4)=0.

  RETURN
END SUBROUTINE BFIELD_SETUP

!
! End of subroutine BFIELD_SETUP
!

!***********************************************************************
!				ARRAYS
!			Set up all the arrays
!***********************************************************************
SUBROUTINE ARRAYS

  use ModSizeDGCPM
  use ModIoDGCPM
  use ModMainDGCPM

  implicit none

  integer I,J,K,L,IDL1,IC,IML
  real DPA,SEG1,SEG2,RWU,CONV,degrad,camlra,muboun,spa
  real FUNT,CON1
  external funt
  real amla0(Slen)

  REAL CONE(NR+4),PA(NPA),MUB(NPA),sumd,sumw,LZMAX,LZMIN

  DATA amla0/0.0,0.2,0.4,0.6,0.8,1.0,1.5,2.0,2.5,3.0,3.5, &
       4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0/

  LZMIN=2.
  LZMAX=6.5
  DL1=((LZMAX-LZMIN)/(IO-2))	! Grid size of L shell
  IDL1=NINT(DL1*1000.)
  IF(DL1.GE.0.25 .AND. MOD(IDL1,250).NE.0) THEN
     WRITE(6,*) ' Error : 0.25 is not a factor of DL1 '
     STOP
  ELSE IF (DL1.LT.0.25 .AND. MOD(250,IDL1).NE.0) THEN
     WRITE (6,*) ' Error : DL1 is not a factor of 0.25'
     STOP
  END IF

  DR=DL1*RE               !Grid size for Z=RO
  do I=1,IO
     LZ(I)=LZMIN+(I-2)*DL1
  end do
  Z(1:IO)=RE*LZ(1:IO)
  degrad=pi/180.

  CALL BFIELD_SETUP(cone)

  DPHI=2.*PI/JO		      ! Grid size for local time [rad]
  IF(MOD(NLT,JO).NE.0) THEN
     WRITE(6,*) ' Error : JO is not a factor of NLT '
     STOP
  END IF

  do J=1,JO+1
     PHI(J)=(J-1)*DPHI	! Magnetic local time in radian
  end do
  MLT(1:JO+1)=PHI(1:JO+1)*12./PI    ! Magnetic local time in hour 

  RETURN
END SUBROUTINE ARRAYS
!
! End of subroutine ARRAYS
!

! **********************************************************************
!                             GETKPA
!	Finds KP and other parameters as well as convection strength A
!	Note [A]= V/m * 1/Re^(lamgam-1), last term =1 (in units of L)
!  IKP choices:
!	IKP=0  Static Kp
!	IKP=1  Decrease by DKP until Kp=1
!	IKP=2  Read from table
!	IKP=3  Read from file
!	IKP=4  Read from MBI file (also read kp file for other inputs)
! **********************************************************************
!! SUBROUTINE GETKPA(i3,nst,i2,nkp)
SUBROUTINE GETKPA(i3)

  use ModSizeDGCPM
  use ModIoDGCPM
  use ModMainDGCPM

  implicit none

  integer i3,i4,JKP,ii
  real KPN,KPO,lambdae,KPtab(48),tol,RSUN,KPP,XKP
  save tol,i4,XKP,KPP,KPN,KPO

  !	Data for simulation of the CRRES observations in 1/1991
  DATA KPtab/1.667,2.000,2.000,2.333,1.777,3.333,4.667,3.000, &
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
  IF (IKP.EQ.1) THEN	! IKP=1, Degrades until Kp<=1
     IF(KP.GT.1.) KP=KP+DKP
  ELSE IF (IKP.EQ.2) THEN	! IKP=2, read from table
     JKP=MIN(48,INT(T/10800.)+1)
     KP=KPtab(JKP) 
  ELSE IF (IKP.GE.3) THEN	! IKP=3 Read from file
     IF (MOD(I3-1,NKP).EQ.0 .OR. I3.EQ.NST) THEN
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
     END IF
     IF (T-TOL.LT.3600.) THEN
        KP=.5*(KPP+KPO)+(XKP-.5*(KPP+KPO))*(T-TOL)/3600.
     ELSE IF (T-TOL.LT.7200.) THEN
        KP=XKP
     ELSE
        KP=XKP+(.5*(KPO+KPN)-XKP)*(T-TOL-7200.)/3600.
     END IF
     !            KP=KPO+(KPN-KPO)*(T-TOL)/10800.
  END IF
  IF (IKP.EQ.4 .OR. IA.EQ.2) THEN   ! Midnight Boundary Index
     IF (I3.EQ.NST) THEN
        I4=0
        ii=2
        DO WHILE (I4.eq.0)
           IF (TLAME(ii).GE.T) THEN
              I4=ii-1
           ELSE IF (ii.EQ.ILAMBE-1) THEN
              I4=ILAMBE-1
           ELSE
              ii=ii+1
           END IF
        END DO
     ELSE IF (T+2.*DT.GE.TLAME(ILAMBE)) THEN
        I4=ILAMBE-1
     ELSE IF (T.GE.TLAME(I4+1)) THEN
        I4=I4+1
     END IF
     IF (T.LE.TLAME(1)) THEN
        LAMBDAE=LAMBE(1)
     ELSE IF (T.GE.TLAME(ILAMBE)) THEN
        LAMBDAE=LAMBE(ILAMBE)
     ELSE
        LAMBDAE=LAMBE(I4)+(LAMBE(I4+1)-LAMBE(I4))*(T-TLAME(I4)) &
             /AMAX1(.1*DT,TLAME(I4+1)-TLAME(I4))
     END IF
     IF (IKP.EQ.4) KP=LAMBDAE
  END IF
  A=7.05E-6/(1.-0.159*KP+0.0093*KP**2)**3
  IF (IA.EQ.2) THEN  ! A from MBI transformation
     A=1.44E-2/LAMGAM*COS(LAMBDAE)**(2.*(LAMGAM+1.))
  END IF
  RETURN
END SUBROUTINE GETKPA
!
! End of subroutine GETKPA
!

! **********************************************************************
!				THERMAL
!	Converts thermal plasma densities from Craig's grid to ours
!	Rewritten for mram04.f to use Dan Ober's plasmasphere code:
!	  the Dynamic Global Core Plasma Model (DGCPM)
!	  including passage of our electric potentials to the DGCPM
! **********************************************************************
SUBROUTINE THERMAL

  use ModSizeDGCPM
  use ModMainDGCPM
  use ModHeidiDGCPM
  use ModIoDGCPM
  use ModConstants

  implicit none

  integer i,j,i1,j1,l,jj,jjj,j2,ier,ntimestep,n
  real EO,FAC,FACI
  real delt, par(2)
  real chi1, kpinit, dtimestep
  character*80 filename

  par(1)=-1.0
  par(2)=-1.0
  ntimestep=10
  dtimestep=1./ntimestep

  !  If ITHERMFIRST=1, then do initial setup for the plasmasphere code and return
  IF (ithermfirst.eq.1) then
     call initmain()
     call getgrid(vthetacells,nthetacells,vphicells,nphicells)
     call getxydipole(nthetacells,nphicells,vthetacells,vphicells, &
          gridx,gridy,gridoc)
     do i=1,nthetacells
        vlzcells(i)=1./sin(dtor*vthetacells(i))**2
     enddo
     do j=1,nphicells
        vmltcells(j)=vphicells(j)/15.
     enddo

     ithermfirst=2
     RETURN

  END IF

  !  Call the plasmasphere code and interpolate onto our spatial grid
  delt=2*DT*dtimestep
  If (ithermfirst.eq.2) then
     If (itherminit.eq.1) then  ! Read initial plasmasphere from file
        filename=cInputDir//name//'_dgcpm.in'
        print *, 'Reading in plasmasphere file: ',filename
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
              fac = vlzcells(i) * sin(dtor*vphicells(j))
              potdgcpm(i,j) = chi1 * fac
           enddo
        enddo
     endif  ! itherminit
  endif    ! ithermfirst
  If (delt.gt.0.) then
     call setpot(vthetacells,nthetacells,vphicells,nphicells,potdgcpm)
     print *, 'Calling plasmasphere:',potdgcpm(nthetacells,nphicells/4), &
          potdgcpm(nthetacells,3*nphicells/4)
     do n=1,ntimestep
        call plasmasphere(delt,par)
     end do
     call getdensity(vthetacells,nthetacells,vphicells,nphicells,dendgcpm)
  endif

  ithermfirst = 0

  RETURN
END SUBROUTINE THERMAL
!
! End of subroutine THERMAL
!

