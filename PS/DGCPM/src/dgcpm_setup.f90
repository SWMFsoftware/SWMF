!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
! File name: dgcpm_setup.f90
!
! Contains: Routines for DGCPM
!	GETKPA
!       SHUE
!	THERMAL
!
! Change from 011 to 012: Supports SWMF Style Input, removes legacy
!        inoperable HEIDI code.
!
! Last Modified: January 2012, Aron Dodger
!
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
!   IKP=5  Read from file (Direct file values, no interpolation)
! **********************************************************************
!! SUBROUTINE GETKPA(i3,nst,i2,nkp)
SUBROUTINE GETKPA(i3)

  use ModIoUnit,  ONLY: UnitTMP_
  use ModSizeDGCPM
  use ModIoDGCPM
  use ModMainDGCPM
  use ModTimeDGCPM
  use ModIndicesInterfaces

  implicit none

  integer, intent(in) :: i3

  integer, save :: i4
  real,    save :: tol=0.0, XKP, KPP, KPN, KPO

  integer :: JKP,ii, i, ierror
  real :: lambdae, KPtab(48), RSUN
  real :: kpt, dut, dkp
  character(len=80) ::  header

  ! Setup Variables
  KPT=-1./3./3600. ! model rate of decay of Kp in s-1
  DKP=KPT*DT*2.    ! discrete step size of Kp

  !	Data for simulation of the CRRES observations in 1/1991
  DATA KPtab/1.667,2.000,2.000,2.333,1.777,3.333,4.667,3.000, &
       2.667,3.667,2.333,2.000,2.000,1.667,1.000,3.333, &
       0.333,1.333,1.000,1.667,2.000,1.333,2.333,3.000, &
       1.000,1.667,1.333,1.333,0.667,1.667,1.000,2.000, &
       2.000,0.667,1.333,1.667,0.667,0.333,0.333,1.667, &
       0.333,0.333,0.333,0.667,1.333,2.000,2.000,1.000/ 

  ! Read KP file if Uninitialized
  IF (IsUninitialized) THEN
    If ((IKP.EQ.3).OR.(IKP.EQ.4).OR.(IKP.EQ.5)) THEN
     OPEN(UnitTMP_,FILE=cInputDir//NAME//'_kp.in',STATUS='OLD') 
     read(UnitTMP_,*) Header
     DO I=1,NSTEP/NKP+2
        READ(UnitTMP_,*) DAYR(I),DUT,RKPH(I),F107R(I),APR(I),RSUNR(I)
     ENDDO
     CLOSE(UnitTMP_)
    ENDIF
  ELSE 
  
  ! IKP=0, keep Kp constant, otherwise...
  IF (IKP.EQ.1) THEN	! IKP=1, Degrades until Kp<=1
     IF(KP.GT.1.) KP=KP+DKP
  
  ELSE IF (IKP.EQ.2) THEN	! IKP=2, read from table
     JKP=MIN(48,INT(MaxTime/10800.)+1)
     KP=KPtab(JKP) 
  
  ELSE IF (IKP.EQ.3 .OR. IKP.EQ.4) THEN	! IKP=3 Read from file
     !IF (MOD(I3-1,NKP).EQ.0 .OR. I3.EQ.NST) THEN
     KPP=RKPH(MAX0(1,I2-1))
     KPO=RKPH(I2)
     KPN=RKPH(MIN0(NSTEP/NKP+2,I2+1))
     XKP=.125*(10.*KPO-KPP-KPN)
     AP=APR(I2)
     F107=F107R(I2)
     RSUN=RSUNR(I2)
     I2=I2+1
     KPN=RKPH(I2)	
     !END IF
     IF (T-TOL.LT.3600.) THEN
        KP=.5*(KPP+KPO)+(XKP-.5*(KPP+KPO))*(T-TOL)/3600.
     ELSE IF (T-TOL.LT.7200.) THEN
        KP=XKP
     ELSE
        KP=XKP+(.5*(KPO+KPN)-XKP)*(T-TOL-7200.)/3600.
     END IF
     !            KP=KPO+(KPN-KPO)*(T-TOL)/10800.
  
  ELSE IF (IKP.EQ.5) THEN ! IKP=5 USE EXACT FILE VALES
    KP=RKPH(MIN0(I2, NSTEP/NKP+2))
    IF (KP.GT.10.) KP=10.
    IF (KP.LT.0.)  KP=0.
  
  ELSE IF (IKP.EQ.6) THEN ! IKP=6 Use NGDC Index File
    call get_kp(CurrentTime, KP, iError)
    write(*,*) 'DGCPM: Kp = ', KP
    if (iError /= 0) then
           write(*,*) "Can not find IMF Bz."
            kp=0.
    ENDIF
  ENDIF

  if (ilambe == 0) ilambe = 3
  IF (IKP.EQ.4) THEN   ! Midnight Boundary Index
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
  
  END IF

    RETURN
END SUBROUTINE GETKPA

!
! End of subroutine GETKPA
!

! **********************************************************************
!                               SHUE
!   Calculates the magnetopause location using the equation described
!   in Shue et al. 1998. Equation is a function of the solar wind 
!   dynamic pressure (in nPa) and Bz field strength (in nT)
! **********************************************************************

    subroutine Shue()

    use ModSizeDGCPM 
    use ModMainDGCPM
    use ModTimeDGCPM
    use ModIndicesInterfaces
    
        implicit none

! r0 is the subsolar standoff distance in Re
        real r0, r, theta, temp_r
! alpha is the power factor for the angular dependance
        real alpha
! Dp is the solar wind dynamic pressure, in nPa
        real Dp

        integer i, j, iError

        call get_IMF_Bz(CurrentTime, IMF_Bz, iError)
        if (iError /= 0) then
           write(*,*) "Can not find IMF Bz."
        endif

        call get_SW_V(CurrentTime, SW_V, iError)
        if (iError /= 0) then
           write(*,*) "Can not find SW_V."
        endif

        call get_SW_N(CurrentTime, SW_N, iError)
        if (iError /= 0) then
           write(*,*) "Can not find IMF Bz."
        endif
        
        SW_V = SW_V * 1000.
        SW_N = SW_N * (1./100.) ** 3.

        Dp = ((SW_V ** 2.0) * SW_N * 1.67262158E-27) * (1E21)

        r0 = (10.22 + 1.29 * tanh( 0.184 * ( IMF_Bz + 8.14 ) ) ) * &
               (Dp ** (- (1.0 / 6.6 ) ) )

        alpha = (0.58 - 0.007 * IMF_Bz) * (1.0 + 0.024 * alog(Dp)) 

        do i=1, nthetacells
            do j=1, nphicells
                temp_r = vrcells(i)
                theta  = atan(mgridx(i,j) / mgridy(i,j))
                r = r0 * ( 2.0 / (1.0 + cos(theta))) ** alpha

                if (r.lt.temp_r) then
                    mgridoc(i,j) = 1. 
                else 
                    mgridoc(i,j) = 0.
                endif
            enddo
        enddo

        return
      END SUBROUTINE SHUE
!
! End of subroutine SHUE
!


! **********************************************************************
!				THERMAL
!	  Dynamic Global Core Plasma Model (DGCPM) is called
!         through Thermal via the plasmasphere command. This routine
!         also performs most of the necessary setup and initialization
!         commands necessary for the model itself. 
! **********************************************************************
SUBROUTINE THERMAL

  use ModSizeDGCPM
  use ModMainDGCPM
  use ModIoDGCPM
  use ModConstants
  use ModCoupleDGCPM
  use ModTimeDGCPM

  implicit none

  integer i,j,i1,j1,l,jj,jjj,j2,ier,ntimestep,n,testi,testj
  real EO,FAC,FACI
  real delt, par(2)
  real chi1, kpinit, dtimestep
  character*80 filename 
  integer ierror

  par(1)=-1.0
  par(2)=-1.0
  ntimestep=10
  dtimestep=1./ntimestep

  !  If ITHERMFIRST=1, then do initial setup for the plasmasphere code
  !  and RETURN
  IF (ithermfirst.eq.1) then
    
     call initmain()
     call getgrid(vthetacells,nthetacells,vphicells,nphicells)
     call getxydipole()
     do i=1,nthetacells
        vlzcells(i)=1./sin(dtor*vthetacells(i))**2.
     enddo
     do j=1,nphicells
        vmltcells(j)=vphicells(j)/15.
     enddo

     ithermfirst=2
     RETURN

  END IF

  ! Call the plasmasphere code and interpolate onto our spatial grid
  delt=2*DT*dtimestep
  If (ithermfirst.eq.2) then
     If (itherminit.eq.1) then  ! Read initial plasmasphere from file
        filename=cInputDir//name//'_dgcpm.in'
        if (debug .gt. 0) write(*,*) 'Reading in plasmasphere file: ',filename
        call loadplasmasphere(filename)
        call getdensity(vthetacells,nthetacells,vphicells,nphicells, &
             dendgcpm)
        delt=0.        
     else	 ! initialize plasmasphere with 48 h of low activity
        if (debug .gt. 0) write(*,*) 'Initializing Plasmasphere with 48h of Low Activity'
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
     if (isCoupled) then
         mgridpot = coupled_potential
     else
         call magconv
     endif

     call setpot(vthetacells,nthetacells,vphicells,nphicells,mgridpot)
    
    ! Shue Magnetopause - Still in testing
    If (UseShue) then
        call shue()
    endif
     
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

