Module ModElecTrans
  PRIVATE
  
  real,allocatable :: ALPHA(JMAX), BETA(JMAX), GAMA(JMAX), PSI(JMAX), &
       DELZ(JMAX), DEL2(JMAX), DELA(JMAX), DELP(JMAX), &
       DELM(JMAX), DELS(JMAX), DEN(JMAX)
  real :: FAC
  
  logical :: IsDebug=.true.
  
  public :: ETRANS
  ! Subroutine ETRANS
  !
  ! M.H.: this appears to calculate electron precipitation impact outcomes on various species
  !
  ! This software is part of the GLOW model.  Use is governed by the Open Source
  ! Academic Research License Agreement contained in the file glowlicense.txt.
  ! For more information see the file glow.txt.
  !
  ! Banks & Nagy 2-stream electron transport code
  ! Adapted by Stan Solomon, 1986, 1988
  ! Uses variable altitude and energy grids
  ! LMAX obtained from glow.h, SMB, 1994
  ! Updated comments and removed artifacts, SCS, 2005
  !
  ! Subroutine EXSECT called first time only to calculate electron impact
  ! cross sections.
  !
  ! Definitions:
  ! COMMON /CGLOW/:  see subroutine GLOW
  ! COMMON /CXSECT/: see subroutine EXSECT
  ! COMMON /CXPARS/: see subroutine EXSECT
  ! PSI    first term of parabolic d.e., = 1
  ! ALPHA  second term "; cm-1
  ! BETA   third term  "; cm-2
  ! GAMA  forth term  "; cm-4 s-1 eV-1
  ! DELZ   altitude increments; cm
  ! DEL2   sum of altitude increment and next higher increment; cm
  ! DELA   average of "
  ! DELP   product of DELA and next higher DELZ
  ! DELM   product of DELA and DELZ
  ! DELS   product of DELZ and next higer DELZ
  ! DEN    dummy array for transfer of calculated downward flux
  ! FAC    factor for extrapolating production rate, = 0
  ! PROD   sum of photo- and secondary electron production; cm-3 s-1 eV-1
  ! EPROD  energy of "; eV cm-3
  ! T1     elastic collision term; cm-1
  ! T2     elastic + inelastic collision term; cm-1
  ! TSA    total energy loss cross section for each species; cm2
  ! PRODUP upward   cascade + secondary production; cm-3 s-1 eV-1
  ! PRODWN downward cascade + secondary production; cm-3 s-1 eV-1
  ! PHIUP  upward   flux; cm-2 s-1 eV-1
  ! PHIDWN downward flux; cm-2 s-1 eV-1
  ! TSIGNE thermal electron collision term; cm-1
  ! SECION total ionization rate; cm-3 s-1
  ! SECP   secondary electron production; cm-3 s-1 eV-1
  ! R1     ratio term for calculating upward flux; cm-2 s-1 eV-1
  ! EXPT2  exponential term for calculating upward flux
  ! PRODUA collection array for calculating PRODUP; cm-3 s-1 eV-1
  ! PRODDA  "                               PRODWN
  ! PHIINF downward flux at top of atmos., divided by AVMU; cm-2 s-1 eV-1
  ! POTION ionizaition potential for each species; eV
  ! AVMU   cosine of the average pitch angle
  !
  ! Array dimensions:
  ! JMAX    number of altitude levels
  ! NBINS   number of energetic electron energy bins
  ! LMAX    number of wavelength intervals for solar flux
  ! NMAJ    number of major species
  ! NEX     number of ionized/excited species
  ! NW      number of airglow emission wavelengths
  ! NC      number of component production terms for each emission
  ! NST     number of states produced by photoionization/dissociation
  ! NEI     number of states produced by electron impact
  ! NF      number of available types of auroral fluxes
  !
  !
  SUBROUTINE ETRANS
    use ModSeGrid, only: NBINS=>nEnergy, JMAX=>nAlt,Alt_C,&
         ENER=>EnergyGrid_I, DEL=>DeltaE_I
    use ModSeProduction, only: PESPEC=>ePhotoProdSpec_IC
    use ModSeCross,only: IIMAXX,SIGI,SIGA,SIGS,SIGEX,NEI,WW,PE,PIN
    use ModBackground,only: dip, NMAJ=>nNeutralSepcies, ZE=>eThermalDensity_C,&
         ZTE=>eThermalTemp_C
    !      use cglow, only: nmaj,jmax,nw,nst,nbins,nei,nf,lmax,nc,nex
    ! *************************
    ! TODO using "implicit none" causes this file to give error with f2py/f2py3
    !constructing wrapper function "aurora"...
    !         pyion,pyecalc,pypi,pysi,pyisr = aurora(z,pyidate,pyut,pyglat,pyglong,pyf107a,pyf107,pyf107p,pyap,pyphitop)
    !{}
    !analyzevars: charselector={'len': '4'} unhandled.analyzevars: charselector={'len': '4'} unhandled.analyzevars: charselector={'len': '4'} unhandled.getctype: No C-type found in "{}", assuming void.
    !
    ! issue with common block? Leave it alone till Stan has new F90 code available.
    ! ************************
    
!    use, intrinsic :: iso_fortran_env, only : stdout=>output_unit, &
    !         stderr=>error_unit
    !      implicit none
!    INCLUDE 'cglow.h'
    
    logical isfinite
    logical :: IsLocal = .false.
    integer  IERR
    
    !phitop should be set
    real PHITOP(NBINS),EFRAC, SION(NMAJ,JMAX), &
       UFLX(NBINS,JMAX), DFLX(NBINS,JMAX), AGLW(NEI,NMAJ,JMAX), &
       EHEAT(JMAX), TEZ(JMAX)
  
  real    PROD(JMAX), EPROD(JMAX), T1(JMAX), T2(JMAX), TSA(NMAJ), &
       PRODUP(JMAX,NBINS), PRODWN(JMAX,NBINS), &
       PHIUP(JMAX), PHIDWN(JMAX), TSIGNE(JMAX), TAUE(JMAX), &
       SECION(JMAX), SECP(NMAJ,JMAX), R1(JMAX), EXPT2(JMAX), &
       PRODUA(JMAX), PRODDA(JMAX), PHIINF(NBINS)
  
  real  APROD,DAG,EDEP,EET,ein,eout,epe,ephi,et,fluxj,phiout, &
       rmusin, sindip
  integer  i,ib,ibb,ii,im,iq,iv,j,jj,jjj4,k,kk,ll,n
  
  integer :: IFIRST=1
  real,parameter :: AVMU=0.5
  
  ! this should be adjusted...planet specific
  real potion(nmaj)
  DATA POTION/16.,16.,18./
  !------------------------------------------------------------------------
  ! allocate solver arrays
  if (.not.allocated(ALPHA)) allocate(ALPHA(JMAX), BETA(JMAX), GAMA(JMAX), &
       PSI(JMAX), DELZ(JMAX), DEL2(JMAX), DELA(JMAX), DELP(JMAX), &
       DELM(JMAX), DELS(JMAX), DEN(JMAX))


  IERR = 0
  FAC = 0.
  SINDIP = SIN(DIP)
  RMUSIN = 1. / SINDIP / AVMU
  
  ! First call only:  calculate cross-sectons:
  !
  !      IF (IFIRST .EQ. 1) THEN
  !        CALL EXSECT (ENER, DEL)
  !        IFIRST = 0
  !      ENDIF
  !
  !
  ! Zero variables:
  !
  !      DO 100 II=1,JMAX
  !        DO 100 IB=1,NMAJ
  !          DO 100 IBB=1,NEI
  !            AGLW(IBB,IB,II) = 0.0
  !  100 CONTINUE
  !
  PSI(1)   = 1.
  ALPHA(1) = 0.
  BETA(1) = 0.
  GAMA(1) = 0.
  PHIOUT = 0.0
  !
  DO I = 1, JMAX
     EHEAT(I) = 0.0
     EPROD(I) = 0.0
     SECION(I) = 0.0
     DO N = 1, NMAJ
        SION(N,I) = 0.0
     End Do
  End Do
  !
  DO JJ = 1, NBINS
     DO I = 1, JMAX
        PRODUP(I,JJ) = 1.0E-20
        PRODWN(I,JJ) = 1.0E-20
     enddo
  enddo
  !
  !
  ! Divide downward flux at top of atmos. by average pitch angle cosine:
  !
  DO  J=1,NBINS
     PHIINF(J) = PHITOP(J) / AVMU
     if(IsDebug .and. .not.isfinite(phiinf(j))) &
          error stop 'etrans: nonfinite PhiInf'
  End Do
  !
  !
  ! Calcualte delta z's:
  !
  DELZ(1) = Alt_C(2)-Alt_C(1)
  DO I=2,JMAX
     DELZ(I) = Alt_C(I)-Alt_C(I-1)
  End Do
  
  DO I=1,JMAX-1
     DEL2(I) = DELZ(I)+DELZ(I+1)
     DELA(I) = DEL2(I)/2.
     DELP(I) = DELA(I)*DELZ(I+1)
     DELM(I) = DELA(I)*DELZ(I)
     DELS(I) = DELZ(I)*DELZ(I+1)
  End Do
  
  DEL2(JMAX) = DEL2(JMAX-1)
  DELA(JMAX) = DELA(JMAX-1)
  DELP(JMAX) = DELP(JMAX-1)
  DELM(JMAX) = DELP(JMAX-1)
  DELS(JMAX) = DELS(JMAX-1)
  !
  !
  !
  ! Top of Energy loop:
  !
  ENERGY_LOOP: do J=NBINS,1,-1
     !
     !
     ! Calculate production:
     !
     DO I = 1, JMAX
        !SESPEC is obsolete
        ! PROD(I) = (PESPEC(J,I)+SESPEC(J,I)) * RMUSIN / DEL(J)
        PROD(I) = (PESPEC(J,I)) * RMUSIN / DEL(J)
        EPROD(I) = EPROD(I) + PROD(I) * ENER(J) * DEL(J) / RMUSIN
     End Do
     !
     !
     ! Total energy loss cross section for each species:
     !
     do I = 1, NMAJ
        TSA(I) = 0.0
     enddo
     IF (J .GT. 1) THEN
        DO K = 1, J-1
           DO I = 1, NMAJ
              TSA(I) = TSA(I) + SIGA(I,K,J) * (DEL(J-K)/DEL(J))
           enddo
        enddo
     ELSE
        DO I=1,NMAJ
           TSA(I) = TSA(I) + SIGA(I,1,J) + 1.E-18
        End Do
     ENDIF
     !
     !
     ! Thermal electron energy loss:
     !
     JJJ4 = J - 1
     IF (J .EQ. 1) JJJ4 = 1
     DAG = ENER(J) - ENER(JJJ4)
     IF (DAG .LE. 0.0) DAG = DEL(1)
     !
     LOOP: DO I = 1, JMAX
        ET = 8.618E-5 * ZTE(I)
        EET = ENER(J) - ET
        IF (EET .LE. 0.0) then
           TSIGNE(I) = 0.0
           cycle LOOP
        endif
        TSIGNE(I) = ((3.37E-12*ZE(I)**0.97)/(ENER(J)**0.94)) &
             * ((EET)/(ENER(J) - (0.53*ET))) ** 2.36
        TSIGNE(I) = TSIGNE(I) * RMUSIN / DAG
     enddo LOOP
     !
     !
     ! Collision terms:
     !
     DO I = 1, JMAX
        T1(I) = 0.0
        T2(I) = 0.0
        DO IV = 1, NMAJ
           T1(I) = T1(I) + ZMAJ(IV,I) * SIGS(IV,J) * PE(IV,J)
           T2(I) = T2(I) + ZMAJ(IV,I) * (SIGS(IV,J)*PE(IV,J) + TSA(IV))
        End Do
        T1(I) = T1(I) * RMUSIN
        T2(I) = T2(I) * RMUSIN + TSIGNE(I)
     enddo
     !
     !
     ! Bypass next section if local calculation was specified:
     !
     IF (.not.IsLocal) then
        !
        !
        ! Solve parabolic d.e. by Crank-Nicholson method to find downward flux:
        !
        DO I = 2, JMAX-1
           PSI(I) = 1.
           ALPHA(I) = (T1(I-1) - T1(I+1)) / (DEL2(I) * T1(I))
           BETA(I) = T2(I) * (T1(I+1) - T1(I-1)) / (T1(I) * DEL2(I)) &
                - (T2(I+1) - T2(I-1)) / DEL2(I) &
                - T2(I)**2 + T1(I)**2
           IF (PROD(I) .LT. 1.E-30) PROD(I) = 1.E-30
           IF (PRODWN(I,J) .LT. 1.E-30) PRODWN(I,J) = 1.E-30
           GAMA(I) = (PROD(I)/2.0) &
                * (-T1(I) - T2(I) - ALPHA(I) &
                - (PROD(I+1) - PROD(I-1))/PROD(I)/DEL2(I)) &
                + PRODWN(I,J) &
                * (-ALPHA(I) - T2(I) &
                - (PRODWN(I+1,J) - PRODWN(I-1,J)) &
                /PRODWN(I,J)/DEL2(I)) &
                - PRODUP(I,J) * T1(I)
           !if (.not.isfinite(GAMA(I))) GAMA(I) = 0.
           if (IsDebug .and. .not.isfinite(GAMA(I))) then
              write (*,*),'etrans.f: GAMA PRODWNn1 PRODWN PRODWNp1', &
                   ' PRODUP PRODn1 PROD PRODp1 T1 T2 ALPHA DEL2'
              write(stderr,*), GAMA(I),PRODWN(I-1,J),PRODWN(I,J), &
                   PRODWN(I+1,J), PRODUP(I,J),PROD(I-1),PROD(I),PROD(I+1), &
                   T1(I),T2(I),ALPHA(I), DEL2(I)
              call con_stop('etran: non-finite GAMA')
           end if
        End DO
        
        IF (ABS(BETA(2)) .LT. 1.E-20) THEN
           BETA(2) = 1.E-20
           IERR = 2
        ENDIF
        PHIDWN(2) = GAMA(2) / BETA(2)
        DEN(1) = PHIDWN(2)
        FLUXJ = PHIINF(J)
        CALL IMPIT(FLUXJ) !computes DEN via module
        DO I = 1, JMAX
           PHIDWN(I) = DEN(I)
           if(IsDebug .and. phidwn(i).gt.1e30) then
              write(*,*) 'GAMA(i), beta(i)',gama(i),beta(i)
              call con_stop('etrans: very large PHIDWN (jlocal != 1)')
           endif
        End Do
        !
        !
        ! Apply lower boundary condition: PHIUP=PHIDWN.  Should be nearly zero.
        ! Then integrate back upward to calculate upward flux:
        !
        PHIUP(1) = PHIDWN(1)
        DO I = 2, JMAX
           R1(I) = (T1(I)*PHIDWN(I) + (PROD(I)+2.*PRODUP(I,J))/2.) / T2(I)
           TAUE(I) = T2(I)*DELZ(I)
           IF (TAUE(I) .GT. 60.) TAUE(I)=60.
           EXPT2(I) = EXP(-TAUE(I))
        enddo
        DO I=2,JMAX
           PHIUP(I) = R1(I) + (PHIUP(I-1)-R1(I)) * EXPT2(I)
           if (IsDebug .and. .not.isfinite(phiup(i))) &
                call con_stop('etrans: nonfinite PHIUP  (.not.IsLocal)')
        End DO
        
        
     else
        !local calculation
        DO I = 1, JMAX
           IF (T2(I) .LE. T1(I)) THEN
              IERR = 1
              T2(I) = T1(I) * 1.0001
           ENDIF
           PHIUP(I)  = (PROD(I)/2.0 + PRODUP(I,J)) / (T2(I) - T1(I))
           PHIDWN(I) = (PROD(I)/2.0 + PRODWN(I,J)) / (T2(I) - T1(I))
           if (IsDebug .and. .not. isfinite(phiup(i))) &
                call con_stop('etrans: nonfinite PHIUP  (jlocal=1)')
           if (IsDebug .and. phidwn(i).gt.1e30) &
                call con_stop('etrans: very large PHIDWN (jlocal=1)')
        End Do
        !
     endif
     
     !
     !
     ! Multiply fluxes by average pitch angle cosine and put in arrays,
     ! and calculate outgoing electron energy flux for conservation check:
     !
     DO I=1,JMAX
        UFLX(J,I) = PHIUP(I) * AVMU
        DFLX(J,I) = PHIDWN(I) * AVMU
     End do
     !
     PHIOUT = PHIOUT + PHIUP(JMAX) * DEL(J) * ENER(J)
     !
     !
     ! Cascade production:
     !
     do K = 1, J-1
        LL = J - K
        DO I=1,JMAX
           PRODUA(I) = 0.0
           PRODDA(I) = 0.0
        enddo
        do N = 1, NMAJ
           do I=1,JMAX
              PRODUA(I) = PRODUA(I) &
                   + ZMAJ(N,I) * (SIGA(N,K,J)*PIN(N,J)*PHIDWN(I) &
                   + (1. - PIN(N,J))*SIGA(N,K,J)*PHIUP(I))
              PRODDA(I) = PRODDA(I) &
                   + ZMAJ(N,I) * (SIGA(N,K,J)*PIN(N,J)*PHIUP(I) &
                   + (1. - PIN(N,J))*SIGA(N,K,J)*PHIDWN(I))
           enddo
        enddo
        do I=1,JMAX
           PRODUP(I,LL) = PRODUP(I,LL) + PRODUA(I) * RMUSIN
           PRODWN(I,LL) = PRODWN(I,LL) + PRODDA(I) * RMUSIN
        enddo
     enddo
     !
     KK = J - 1
     IF (KK>0) then
        do I = 1, JMAX
           PRODUP(I,KK) = PRODUP(I,KK) + TSIGNE(I) * PHIUP(I) * (DEL(J) &
                / DEL(KK))
           PRODWN(I,KK) = PRODWN(I,KK) + TSIGNE(I) * PHIDWN(I) * (DEL(J) &
                / DEL(KK))
        enddo
     endif
     !
     !
     ! Electron heating rate:
     !
     DAG = DEL(J)
     DO I = 1, JMAX
        EHEAT(I) = EHEAT(I) + TSIGNE(I) * (PHIUP(I)+PHIDWN(I)) * DAG**2
     End Do
     !
     !
     ! Electron impact excitation rates:
     !
     DO II = 1, JMAX
        DO I = 1, NMAJ
           DO IBB = 1, NEI
              AGLW(IBB,I,II) = AGLW(IBB,I,II) + (PHIUP(II) + PHIDWN(II)) &
                   * SIGEX(IBB,I,J) * DEL(J) * ZMAJ(I,II)
           enddo
        enddo
     enddo
     !
     !
     ! Calculate production of secondaries into K bin for energy J bin and
     ! add to production:
     !
     DO  K = 1, IIMAXX(J) ! iimaxx set near exsect.f:424
        DO  N = 1, NMAJ
           DO  I = 1, JMAX
              !            SECP(N,I) = SEC(N,K,J) * ZMAJ(N,I) * (PHIUP(I) + PHIDWN(I))
              SECP(N,I) = SIGI(N,K,J) * ZMAJ(N,I) * (PHIUP(I) + PHIDWN(I))
              !             call get_secprod(J,N,ZMAJ(N,I),PHIUP(I) + PHIDWN(I),SECP(N,I)
              !            if (isnan(sec(n,k,j))) stop 'etrans: NaN in SEC'
              !            if (isnan(zmaj(n,i))) stop 'etrans: NaN in ZMAJ'
              !            if (isnan(phiup(i))) stop 'etrans: NaN in PHIUP'
              if (IsDebug .and. phidwn(i).gt.1e30) &
                   call con_stop('etrans.f: very large PHIDWN')
              if (IsDebug .and. .not.isfinite(secp(n,i))) &
                   call con_stop('etrans: nonfinite SECP')
              SION(N,I) = SION(N,I) + SECP(N,I) * DEL(K)
              
              if (IsDebug .and. .not.isfinite(sion(n,i)))  then
                 write(*,*)'etrans: nonfinite impact ioniz SION.', &
                      '  n,i=',n,i
                 write(*,*) 'del(k)=',del(k)
                 write(*,*) 'secp(n,i)=',secp(n,i)
                 write(*,*) 'sion(n,i)=',sion(n,i)
                 call con_stop('stoping in secondary production')
              endif
              
              SECION(I) = SECION(I) + SECP(N,I) * DEL(K)
              PRODUP(I,K) = PRODUP(I,K) + (SECP(N,I)*.5*RMUSIN)
              PRODWN(I,K) = PRODWN(I,K) + (SECP(N,I)*.5*RMUSIN)
           enddo
        enddo
     enddo
     
     !
  enddo ENERGY_LOOP  ! Bottom of Energy loop
  
  DO I = 1, JMAX
     EHEAT(I) = EHEAT(I) / RMUSIN
  end DO
  !
  !
  ! Calculate energy deposited as a function of altitude
  ! and total energy deposition:
  !
  EDEP = 0.
  DO IM=1,JMAX
     TEZ(IM) = EHEAT(IM)
     DO II=1,NMAJ
        TEZ(IM) = TEZ(IM) + SION(II,IM)*POTION(II)
        DO IQ=1,NEI
           TEZ(IM) = TEZ(IM) + AGLW(IQ,II,IM)*WW(IQ,II)
        enddo
     enddo
     EDEP = EDEP + TEZ(IM) * DELA(IM)
  enddo
  !
  !
  ! Calculate energy input, output, and fractional conservation:
  !
  EPE = 0.0
  EPHI = 0.0
  DO I = 2, JMAX
     APROD = SQRT(EPROD(I)*EPROD(I - 1))
     EPE = EPE + APROD * DELZ(I)
  enddo
  DO JJ = 1, NBINS
     EPHI = EPHI + PHIINF(JJ) * ENER(JJ) * DEL(JJ) / RMUSIN
  End Do
  EIN = EPHI + EPE
  PHIOUT = PHIOUT / RMUSIN
  EOUT = EDEP + PHIOUT
  EFRAC = (EOUT - EIN) / EIN

END SUBROUTINE ETRANS
!
!
!
!
! Subroutine IMPIT solves parabolic differential equation by implicit
! Crank-Nicholson method
!
      SUBROUTINE IMPIT(FLUXJ)
  use ModSeGrid, only: NBINS=>nEnergy, JMAX=>nAlt,Alt_C,&
       ENER=>EnergyGrid_I, DEL=>DeltaE_I
!      use, intrinsic :: iso_fortran_env, only : stdout=>output_unit, &
!                                                stderr=>error_unit
!      use cglow,only: jmax
      Implicit None
!      include 'cglow.h'
!Args:
      Real, Intent(In) :: FLUXJ
!Local:
      logical isfinite
      Real dem
      real,dimension(jmax) ::K, L, A, B, C, D
      Integer i,i1,jk,kk

!      COMMON /CIMPIT/ ALPHA, BETA, GAMA, PSI, DELZ, DEL2, DELA, DELP, &
!                      DELM, DELS, DEN, FAC
!
      I1 = JMAX - 1

      DO I = 1, I1
        A(I) = PSI(I) / DELP(I) + ALPHA(I) / DEL2(I)
        B(I) = -2. * PSI(I) / DELS(I) + BETA(I)
        C(I) = PSI(I) / DELM(I) - ALPHA(I) / DEL2(I)
        D(I) = GAMA(I)
       if(IsDebug .and. .not.isfinite(c(i))) &
                    error stop 'etrans:impit nonfinite C(I)'
       if(IsDebug .and. .not.isfinite(d(i))) &
                    error stop 'etrans:impit nonfinite D(I)'
      End Do

      K(2) = (D(2) - C(2)*DEN(1)) / B(2)
      L(2) = A(2) / B(2)
      if (IsDebug .and. .not.isfinite(k(2))) then
        write(stderr,*) '***********    **********'
        write(stderr,*) 'K(2) D(2) C(2) DEN(1) B(2)',K(2),D(2),C(2), &
         DEN(1),B(2)
        error stop 'etrans:impit non-finite K(2)'
      end if
!      if (isnan(k(2))) stop 'etrans:impit NaN in K(2)'

      DO I = 3, I1
        DEM = B(I) - C(I) * L(I-1)
        if (IsDebug .and. .not.isfinite(dem)) &
              error stop 'etrans:impit nonfinite DEM'
        K(I) = (D(I) - C(I)*K(I-1)) / DEM
        L(I) = A(I) / DEM

        if (IsDebug .and. .not. isfinite(K(i))) then
         write(stderr,*) k
         error stop 'etrans:impit NaN in K(i)'
        end if
!        if (isnan(L(i))) stop'etrans:impit NaN in L(i)'
      End DO

      DEN(I1) = (K(I1) - L(I1)*FLUXJ) / (1. + L(I1)*FAC)
!      if (isnan(K(i1))) stop'etrans:impit NaN in K(i1)'
!      if (isnan(L(i1))) stop'etrans:impit NaN in L(i1)'
      if(IsDebug .and. .not.isfinite(den(i1))) &
                error stop 'impit nonfinite DEN(I1)'
      DEN(JMAX) = DEN(I1)

      Do KK = 1, JMAX-3
        JK = I1 - KK
        DEN(JK) = K(JK) - L(JK) * DEN(JK + 1)
        if (IsDebug .and. .not.isfinite(den(jk))) &
          error stop 'etrans:impit: non-finite DEN'

      End Do

      END Subroutine IMPIT


      logical Function isfinite(x)
      implicit none

      real,intent(in) :: x

      if (abs(x) <= huge(x)) then
        isfinite = .true.
      else
        isfinite = .false.
      end if

      end function isfinite

!      !========================================================================
!      SUBROUTINE get_secprod(iEnergyIn,iNeutral,NeutralDens,NetFlux,SecProdTotal)
!        use ModSeGrid,      ONLY: nEnergy, DeltaE_I,EnergyGrid_I,BINNUM
!        use ModMath,        ONLY: midpnt_int
!        use ModNumConst,    ONLY: cPi
!        use ModSeCross,     ONLY: SIGS,SIGI,SIGA
!
!        IMPLICIT NONE
!        integer, intent(in) :: iEnergyIn,iNeutral
!        real   , intent(in) :: NeutralDens,NetFlux
!        !incomming crossections
!!        real   , intent(in) :: SIGS(nNeutral,nEnergy)
!!        real   , intent(in) :: SIGI(nNeutral,nEnergy,nEnergy)
!!        real   , intent(in) :: SIGA(nNeutral,nEnergy,nEnergy)
!        !outgoing electron production rate
!        real   , intent(out):: SecProdTotal
!        !electron production rate from secondary production for each neutral
!        real   :: NetFlux, SecProd(nEnergy)
!        
!        ! Minimum ionization threshold
!        real,parameter :: Eplus = 12.0 
!        INTEGER m,n,LL,jj
!        !--------------------------------------------------------------------------
!        
!        SecProd(:)=0.0
!        SecProdTotal=0.0
!        ! this calculates the electron total production
!        IF (2*EnergyGrid_I(iEnergyIn)+Eplus &
!             < EnergyGrid_I(nEnergy)+.5*DeltaE_I(nEnergy)) THEN
!           LL=BINNUM(2*EnergyGrid_I(iEnergyIn)+Eplus)
!           
!           DO jj=LL,nEnergy
!              SecProd(jj)=&
!                   SecProd(jj)+NeutralDens_I(n)*SIGI(n,iEnergyIn,jj)*NetFlux
!           enddo
!        END IF
!        ! energy-integrated secondary production per species into Qe
!        CALL midpnt_int(SecProdTotal,SecProd(:),DeltaE_I,1,nEnergy,nEnergy,2)
!        
!        RETURN
!      end SUBROUTINE get_secprod
!
!      SUBROUTINE midpnt_int(sum,f,x,a,b,N,ii)
!        IMPLICIT NONE
!        INTEGER a,b,N,j,ii
!        REAL f(N),x(N),sum
!        
!        sum=0.
!        if ((b-a.LT.0).OR.(b.GT.N).OR.(a.LE.0)) RETURN
!        if (ii.EQ.1) then
!           do  j=a,b-1
!              sum=sum+(x(j+1)-x(j))*(f(j)+f(j+1))*0.5
!           enddo
!        else      ! ii=2
!           do j=a,b
!              sum=sum+x(j)*f(j)
!           enddo
!        END IF
!        RETURN
!      END SUBROUTINE midpnt_int
