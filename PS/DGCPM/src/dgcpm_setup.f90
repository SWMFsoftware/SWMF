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
!
!       Previously, five Kp choices were available via "IKP".  These have
!       been reduced to two options, set via the #KP param command.
!
!  NameSourceKp Options:
!       'const'   Constant Kp.
!       'file'    Read from NGDC-formatted file; interpolate.
!
!  OBSOLETE IKP choices [removed]:
!	IKP=0  Static Kp
!	IKP=1  Decrease by DKP until Kp=1
!	IKP=2  Read from table
!	IKP=3  Read from file
!	IKP=4  Read from MBI file (also read kp file for other inputs)
!       IKP=5  Read from file (Direct file values, no interpolation)
!       IKP=6  Read from NGDC formatted file with interpolation.
! **********************************************************************
!! SUBROUTINE GETKPA(i3,nst,i2,nkp)
SUBROUTINE GETKPA

  use ModIoDGCPM,   ONLY: Kp, NameSourceKp, KpConst
  use ModMainDGCPM, ONLY: A
  use ModTimeDGCPM, ONLY: CurrentTime
  use ModIndicesInterfaces

  implicit none

  integer :: ierror

  logical :: DoTest, DoTestMe
  character(len=*), parameter :: NameSub = 'getkpa'  
  !--------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  ! Get Kp values either constant or from file:
  if (NameSourceKp .eq. 'const') then
     Kp=kpConst
  else if (index(NameSourceKp, 'file')>0) then
     ! Interpolate to current time:
     call get_kp(CurrentTime, Kp, iError)
     if (iError /= 0) call CON_stop(NameSub//' DGCPM error getting Kp')
  else
     call CON_stop(NameSub//' Unrecognized NameSourceKp='//NameSourceKp)
  end if

  A=7.05E-6/(1.-0.159*KP+0.0093*KP**2)**3

  ! Write test info to screen as necessary:
  if (DoTestMe) write(*,'(a,f15.1,a, f4.2, E12.3)') &
       'PS/DGCPM: T=', CurrentTime, ' Kp, A = ', Kp, A
  
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

  integer :: i,j,i1,j1,l,jj,jjj,j2,ier,ntimestep,n,testi,testj
  real :: EO,FAC,FACI
  real :: delt, par(2)
  real :: chi1, kpinit, dtimestep
  character(len=100) ::  filename =''
  integer :: ierror

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
        filename=cRestartIn//'dgcpm_restart.dat'
        if (debug .gt. 0) write(*,*) &
             'Reading in plasmasphere file: ', trim(filename)
        call load_restart_file(trim(filename))
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

