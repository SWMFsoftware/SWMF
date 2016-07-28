Module ModElecTrans
  implicit none

  PRIVATE

  real,allocatable :: ALPHA(:), BETA(:), GAMA(:), PSI(:), &
       DELZ(:), DEL2(:), DELA(:), DELP(:), &
       DELM(:), DELS(:), DEN(:)
  real :: FAC
  
  logical :: IsDebug=.true.

  real,public,allocatable :: HeatingRate_C(:)! volume heating rate [eV/cm3/s]
  real,public,allocatable :: NumberDens_C(:) ! number density of SE [/cm3]
  real,public,allocatable :: NumberFlux_C(:) ! number flux of SE [/cm2/s]
  real,public,allocatable :: SecondaryIonRate_IC(:,:) 
  real,public,allocatable :: TotalIonizationRate_C(:) 

  !Precipitation info
  real,allocatable :: PrecipCombinedPhi_I(:)
  logical, public :: UsePrecipitation = .false.
  real   , public :: PrecipEmin, PrecipEmax,PrecipEmean,PrecipEflux
  !Polar Rain info
  logical, public :: UsePolarRain = .false.
  real   , public :: PolarRainEmin=80.0,  PolarRainEmax=300.0, &
       PolarRainEmean=100.0, PolarRainEflux=1.0e-2
  !Ovation Precipitation info
  logical, public :: UseOvation
  real   , public :: OvationEmin,OvationEmax
  real   , public :: EMeanDiff,EFluxDiff,EMeanWave,EFluxWave,EMeanMono,EFluxMono

  !time, used for output only
  real,public ::Time=0.0
  
  public :: ETRANS
  
contains
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
         ENER=>EnergyGrid_I, DEL=>DeltaE_I,nAltExtended,DeltaPot_C,IsVerbose
    use ModSeCross,only: IIMAXX,SIGIX,SIGA,SIGS,SIGEX,SEC,NEI,WW,PE,PIN
!    use ModSeCross,only: IIMAXX,SIGIX=>SIGI,SIGA,SIGS,SIGEX,NEI,WW,PE,PIN
    use ModSeBackground,only:dip, NMAJ=>nNeutralSpecies, ZE=>eThermalDensity_C,&
         ZTE=>eThermalTemp_C,PESPEC=>ePhotoProdSpec_IC,ZMAJ=>NeutralDens_IC
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
    
!    logical,external :: isfinite
    logical :: IsLocal = .false.
    integer  IERR
    
    !phitop should be set
    real PHITOP(NBINS),EFRAC, SION(NMAJ,JMAX), &
         UFLX(NBINS,JMAX), DFLX(NBINS,JMAX), AGLW(NEI,NMAJ,JMAX), &
         EHEAT(JMAX), TEZ(JMAX)
    real :: PolarRainPhi_I(NBINS),PrecipPhi_I(NBINS)
    real :: OvationDiffPhi_I(NBINS),OvationWavePhi_I(NBINS),&
         OvationMonoPhi_I(NBINS)

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
    
    ! For convergence: 
    ! DiffMax is the maximum difference between successive solutions
    ! needed for convergence.
    real :: DiffMax, LastPHIUP(NBINS,JMAX), LastPHIDWN(NBINS,JMAX),DiffMax_I(NBINS)
    integer :: iEnergy
    
    !logical to define if we should iterate flux calculation to reflect portion 
    ! of flux below high altitude potential drop
    logical :: DoIterateFlux=.true.
    ! this should be adjusted...planet specific
    real,allocatable:: potion(:)

    !  DATA potion/16.,16.,18./
    !------------------------------------------------------------------------
    !DeltaPot_C(:)=0.0
    DiffMax=1.0
    
    !set the precipition arrays
    if (.not.allocated(PrecipCombinedPhi_I))&
         allocate(PrecipCombinedPhi_I(NBINS))
    if (UsePrecipitation) then
       call maxt(PrecipEflux, PrecipEmean, ENER, DEL,0, 0.0, 0.0, PrecipPhi_I)
    else
       PrecipPhi_I(:)=0.0
    endif
    
    if (UsePolarRain) then
       call maxt(PolarRainEflux, PolarRainEmean, ENER, DEL,0, 0.0, 0.0, &
            PolarRainPhi_I)
    else
       PolarRainPhi_I(:)=0.0
    endif
    
    if (UseOvation) then
       call maxt(EfluxDiff, EMeanDiff, ENER, DEL,0, 0.0, 0.0, &
            OvationDiffPhi_I)
       call maxt(EfluxWave, EMeanWave, ENER, DEL,0, 0.0, 0.0, &
            OvationWavePhi_I)
       call maxt(0.0, 1.0, ENER, DEL,0, EfluxMono, EmeanMono, &
            OvationMonoPhi_I)
    else
       OvationDiffPhi_I(:)=0.0
       OvationWavePhi_I(:)=0.0
       OvationMonoPhi_I(:)=0.0
    endif
    
    !now combine the different types of precipitation
    PrecipCombinedPhi_I(:)=0.0
    do iEnergy=1,nBINS
       PrecipCombinedPhi_I(iEnergy)=0.0
       if (ENER(iEnergy)>PrecipEmin .and. ENER(iEnergy)<PrecipEmax) then
          PrecipCombinedPhi_I(iEnergy)=&
               PrecipCombinedPhi_I(iEnergy)+PrecipPhi_I(iEnergy)
       endif
       if (ENER(iEnergy)>PolarRainEmin .and. ENER(iEnergy)<PolarRainEmax) then
          PrecipCombinedPhi_I(iEnergy)=&
               PrecipCombinedPhi_I(iEnergy)+PolarRainPhi_I(iEnergy)
       endif
       if (ENER(iEnergy)>OvationEmin .and. ENER(iEnergy)<OvationEmax) then
          PrecipCombinedPhi_I(iEnergy)=&
               PrecipCombinedPhi_I(iEnergy)+OvationDiffPhi_I(iEnergy)&
               +OvationWavePhi_I(iEnergy)+OvationMonoPhi_I(iEnergy)
       endif
    enddo
    if (UseOvation .or. UsePrecipitation .or. UsePolarRain) then
       PHITOP(:)=PrecipCombinedPhi_I(:)
    else
       PHITOP(:)=0.0
    endif
    
    write(*,*)2.0*3.14*sum(PHITOP(:)*DEL(:))*AVMU
    
    do while (DiffMax >0.1 .and. DoIterateFlux)
       
       !reflect solution below max potential drop
       do iEnergy=1,nBINS
          if (ENER(iEnergy)<DeltaPot_C(nAltExtended)) then 
             PHITOP(iEnergy)=LastPhiUp(iEnergy,JMAX)
             !write(*,*) 'E,Phi',ENER(iEnergy),PhiUp(iEnergy)
          else
             if (UseOvation .or. UsePrecipitation .or. UsePolarRain) then
                PHITOP(iEnergy)=PrecipCombinedPhi_I(iEnergy)
             else
                PHITOP(iEnergy)=0.0
             endif
          end if
       enddo
       
       !allocate return vars and initiate to zero
       if (.not.allocated(HeatingRate_C))allocate(HeatingRate_C(nAltExtended))
       if (.not.allocated(SecondaryIonRate_IC))&
            allocate(SecondaryIonRate_IC(NMAJ,nAltExtended))
       if (.not.allocated(TotalIonizationRate_C))&
            allocate(TotalIonizationRate_C(nAltExtended))
       
       HeatingRate_C(:)=0.0
       SecondaryIonRate_IC(:,:)=0.0
       TotalIonizationRate_C(:)=0.0
       
       if(.not.allocated(potion)) then
          allocate(potion(nmaj))
          potion(1)=16.
          potion(2)=16.
          potion(3)=18.
       end if
       
       ! allocate solver arrays
       if (.not.allocated(ALPHA)) allocate(ALPHA(JMAX), BETA(JMAX), GAMA(JMAX), &
            PSI(JMAX), DELZ(JMAX), DEL2(JMAX), DELA(JMAX), DELP(JMAX), &
            DELM(JMAX), DELS(JMAX), DEN(JMAX))
       
       
       IERR = 0
       FAC = 0.
       SINDIP = SIN(DIP)
       RMUSIN = 1. / SINDIP / AVMU
       
       ! set phitop
       !    PHITOP(:)=0.0
       
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
          !       if(IsDebug .and. .not.isfinite(phiinf(j))) &
          !            error stop 'etrans: nonfinite PhiInf'
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
          
          !write(*,*) 'Working on energy bin',J
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
                !  if (IsDebug .and. .not.isfinite(GAMA(I))) then
                !     write (*,*),'etrans.f: GAMA PRODWNn1 PRODWN PRODWNp1', &
                !          ' PRODUP PRODn1 PROD PRODp1 T1 T2 ALPHA DEL2'
                !     write(*,*), GAMA(I),PRODWN(I-1,J),PRODWN(I,J), &
                !     PRODWN(I+1,J), PRODUP(I,J),PROD(I-1),PROD(I),PROD(I+1), &
                !          T1(I),T2(I),ALPHA(I), DEL2(I)
                !     call con_stop('etran: non-finite GAMA')
                !  end if
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
             ! Apply lower boundary condition: PHIUP=PHIDWN.  
             ! Should be nearly zero.
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
                !             if (IsDebug .and. .not.isfinite(phiup(i))) &
                !      call con_stop('etrans: nonfinite PHIUP  (.not.IsLocal)')
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
                !             if (IsDebug .and. .not. isfinite(phiup(i))) &
                !           call con_stop('etrans: nonfinite PHIUP  (jlocal=1)')
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
                   SECP(N,I) = SEC(N,K,J) * ZMAJ(N,I) * (PHIUP(I) + PHIDWN(I))
                   !SECP(N,I) = SIGIX(N,K,J) * ZMAJ(N,I) * (PHIUP(I) + PHIDWN(I))
                   !call get_secprod(J,N,ZMAJ(N,I),PHIUP(I) + PHIDWN(I),SECP(N,I)
                   !    if (isnan(sec(n,k,j))) stop 'etrans: NaN in SEC'
                   !    if (isnan(zmaj(n,i))) stop 'etrans: NaN in ZMAJ'
                   !     if (isnan(phiup(i))) stop 'etrans: NaN in PHIUP'
                   if (IsDebug .and. phidwn(i).gt.1e30) &
                        call con_stop('etrans.f: very large PHIDWN')
                   !                if (IsDebug .and. .not.isfinite(secp(n,i))) &
                   !                     call con_stop('etrans: nonfinite SECP')
                   SION(N,I) = SION(N,I) + SECP(N,I) * DEL(K)
                   
                   !     if (IsDebug .and. .not.isfinite(sion(n,i)))  then
                   !       write(*,*)'etrans: nonfinite impact ioniz SION.', &
                   !            '  n,i=',n,i
                   !       write(*,*) 'del(k)=',del(k)
                   !       write(*,*) 'secp(n,i)=',secp(n,i)
                   !       write(*,*) 'sion(n,i)=',sion(n,i)
                   !       call con_stop('stoping in secondary production')
                   !      endif
                   
                   SECION(I) = SECION(I) + SECP(N,I) * DEL(K)
                   SecondaryIonRate_IC(N,I)=SecondaryIonRate_IC(N,I)&
                        +SECP(N,I) * DEL(K)
                   PRODUP(I,K) = PRODUP(I,K) + (SECP(N,I)*.5*RMUSIN)
                   PRODWN(I,K) = PRODWN(I,K) + (SECP(N,I)*.5*RMUSIN)
                enddo
             enddo
          enddo
          
          !check convergence
          DiffMax_I(J)= &
               max(maxval(PhiUp(:)-LastPhiUp(J,:)),&
               maxval(PhiDwn(:)-LastPhiDwn(J,:)))
          LastPhiUp(J,:)=PhiUp
          LastPhiDwn(J,:)=PhiDwn
       enddo ENERGY_LOOP  ! Bottom of Energy loop
       
       !add secondary production to totalionization rate
       !add photoproduction to total ionization rate
       do I=1,JMAX
          TotalIonizationRate_C(I)=sum(PESPEC(:,I)*DEL(:))+SECION(I)
          !       write(*,*) 'Alt_C(I)*1e-5,ionrate',&
          !            Alt_C(I)*1e-5,TotalIonizationRate_C(I)
          
       enddo
       
       DO I = 1, JMAX
          EHEAT(I) = EHEAT(I) / RMUSIN
          !       write(*,*) 'Alt_C(I)*1e-5,EHEAT(I)',Alt_C(I)*1e-5,EHEAT(I)
          HeatingRate_C(I)=EHEAT(I)
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
          !       write(*,*)'total edep=', EDEP
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
       
       !check convergence
       DiffMax= &
            maxval(DiffMax_I)
       if (IsVerbose)write(*,*) 'DiffMax',DiffMax
       
    end do !end while
    !call plot_omni_iono(time,uFlux_IC,dFlux_IC)
    call plot_omni_iono(time,uFlx,dFlx)
    
    call calc_integrated_values(AVMU,uFlx,dFlx)
    !call map_flux(AVMU,uFlx)
    call plot_integrated(time)
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
!    logical,external:: isfinite
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
!       if(IsDebug .and. .not.isfinite(c(i))) &
!            error stop 'etrans:impit nonfinite C(I)'
!       if(IsDebug .and. .not.isfinite(d(i))) &
!            error stop 'etrans:impit nonfinite D(I)'
    End Do
    
    K(2) = (D(2) - C(2)*DEN(1)) / B(2)
    L(2) = A(2) / B(2)
!    if (IsDebug .and. .not.isfinite(k(2))) then
!       write(*,*) '***********    **********'
!       write(*,*) 'K(2) D(2) C(2) DEN(1) B(2)',K(2),D(2),C(2), &
!            DEN(1),B(2)
!       error stop 'etrans:impit non-finite K(2)'
!    end if
    !      if (isnan(k(2))) stop 'etrans:impit NaN in K(2)'
    
    DO I = 3, I1
       DEM = B(I) - C(I) * L(I-1)
!       if (IsDebug .and. .not.isfinite(dem)) &
!            error stop 'etrans:impit nonfinite DEM'
       K(I) = (D(I) - C(I)*K(I-1)) / DEM
       L(I) = A(I) / DEM
       
!       if (IsDebug .and. .not. isfinite(K(i))) then
!          write(*,*) k
!          error stop 'etrans:impit NaN in K(i)'
!       end if
       !        if (isnan(L(i))) stop'etrans:impit NaN in L(i)'
    End DO
    
    DEN(I1) = (K(I1) - L(I1)*FLUXJ) / (1. + L(I1)*FAC)
    !      if (isnan(K(i1))) stop'etrans:impit NaN in K(i1)'
    !      if (isnan(L(i1))) stop'etrans:impit NaN in L(i1)'
!    if(IsDebug .and. .not.isfinite(den(i1))) &
!         error stop 'impit nonfinite DEN(I1)'
    DEN(JMAX) = DEN(I1)
    
    Do KK = 1, JMAX-3
       JK = I1 - KK
       DEN(JK) = K(JK) - L(JK) * DEN(JK + 1)
!       if (IsDebug .and. .not.isfinite(den(jk))) &
!            error stop 'etrans:impit: non-finite DEN'
       
    End Do
    
  END Subroutine IMPIT
  
  
!  function isfinite(x)
!    implicit none
!    logical :: isfinite
!    real,intent(in) :: x
!    
!    if (abs(x) <= huge(x)) then
!       isfinite = .true.
!    else
!       isfinite = .false.
!    end if
!    
!  end function isfinite
  
  !      !========================================================================
  !      SUBROUTINE get_secprod(iEnergyIn,iNeutral,NeutralDens,NetFlux,SecProdTotal)
  !        use ModSeGrid,      ONLY: nEnergy, DeltaE_I,EnergyGrid_I,BINNUM
  !        use ModMath,        ONLY: midpnt_int
  !        use ModNumConst,    ONLY: cPi
  !        use ModSeCross,     ONLY: SIGS,SIGIX,SIGA
  !
  !        IMPLICIT NONE
  !        integer, intent(in) :: iEnergyIn,iNeutral
  !        real   , intent(in) :: NeutralDens,NetFlux
  !        !incomming crossections
  !!        real   , intent(in) :: SIGS(nNeutral,nEnergy)
  !!        real   , intent(in) :: SIGIX(nNeutral,nEnergy,nEnergy)
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
  !                   SecProd(jj)+NeutralDens_I(n)*SIGIX(n,iEnergyIn,jj)*NetFlux
  !           enddo
  !        END IF
  !        ! energy-integrated secondary production per species into Qe
  !        CALL midpnt_int(SecProdTotal,SecProd(:),DeltaE_I,1,nEnergy,nEnergy,2)
  !        
  !        RETURN
  !      end SUBROUTINE get_secprod
  !
  SUBROUTINE midpnt_int(sum,f,x,a,b,N,ii)
    IMPLICIT NONE
    INTEGER a,b,N,j,ii
    REAL f(N),x(N),sum
    
    sum=0.
    if ((b-a.LT.0).OR.(b.GT.N).OR.(a.LE.0)) RETURN
    if (ii.EQ.1) then
       do  j=a,b-1
          sum=sum+(x(j+1)-x(j))*(f(j)+f(j+1))*0.5
       enddo
    else      ! ii=2
       do j=a,b
          sum=sum+x(j)*f(j)
       enddo
    END IF
    RETURN
  END SUBROUTINE midpnt_int

    ! plot omnidirectional flux in the ionosphere
  !============================================================================
  subroutine plot_omni_iono(time,uFlux_IC,dFlux_IC)
    use ModSeGrid,     ONLY: Alt_C, nAlt, nEnergy,&
         EnergyGrid_I,iLineGlobal
    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModNumConst,   ONLY: cRadToDeg,cPi
    
    real,    intent(in) :: time
    real,    intent(in) :: uFlux_IC(nEnergy,nAlt),&
         dFlux_IC(nEnergy,nAlt)
    
    real, allocatable   :: Coord_DII(:,:,:), PlotState_IIV(:,:,:)

    !grid parameters
    integer, parameter :: nDim =2, nVar=3, S_=2, E_=1
    !variable parameters
    integer, parameter :: Flux_=1, FluxUp_=2, FluxDn_=3
    
    character(len=100),parameter :: NamePlotVar='E[eV] Alt[km] NetFlux[cm-2eV-1s-1] FluxUp[cm-2eV-1s-1] FluxDn[cm-2eV-1s-1] g r'
    character(len=100) :: NamePlot='OmniIono.out'
    
    character(len=*),parameter :: NameHeader='SE output iono'
    character(len=5) :: TypePlot='ascii'
    integer :: iAlt,iEnergy

    logical :: IsFirstCall1=.true.,IsFirstCall2=.true.
    !--------------------------------------------------------------------------
    
    allocate(Coord_DII(nDim,nEnergy,nAlt),PlotState_IIV(nEnergy,nAlt,nVar))
    
    !    do iLine=1,nLine
    PlotState_IIV = 0.0
    Coord_DII     = 0.0
    
    !Set values
    do iEnergy=1,nEnergy
       do iAlt=1,nAlt
             Coord_DII(E_,iEnergy,iAlt) = EnergyGrid_I(iEnergy)             
             Coord_DII(S_,iEnergy,iAlt) = Alt_C(iAlt)/1e5
             ! set plot state
             PlotState_IIV(iEnergy,iAlt,Flux_)  = &
                  uFlux_IC(iEnergy,iAlt)-dFlux_IC(iEnergy,iAlt)
             PlotState_IIV(iEnergy,iAlt,FluxUp_)  = &
                  uFlux_IC(iEnergy,iAlt)
             PlotState_IIV(iEnergy,iAlt,FluxDn_)  = &
                  dFlux_IC(iEnergy,iAlt)
          enddo
       enddo
       
       ! set name for plotfile
       write(NamePlot,"(a,i4.4,a)") 'OmniIono_',iLineGlobal,'.out'
  
       !Plot grid for given line
       if(IsFirstCall1) then
          call save_plot_file(NamePlot, TypePositionIn='rewind', &
               TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
               NameVarIn = NamePlotVar, nStepIn=0,TimeIn=time,     &
               nDimIn=nDim,CoordIn_DII=Coord_DII,                &
               VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
          IsFirstCall1 = .false.
       else
          call save_plot_file(NamePlot, TypePositionIn='append', &
               TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
               NameVarIn = NamePlotVar, nStepIn=0,TimeIn=time,     &
               nDimIn=nDim,CoordIn_DII=Coord_DII,                &
               VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
       endif
    
    deallocate(Coord_DII, PlotState_IIV)
  end subroutine plot_omni_iono

    ! plot omnidirectional flux in the ionosphere
  !============================================================================
  subroutine plot_integrated(time)
    use ModSeGrid,     ONLY: AltExtended_C, nAltExtended, iLineGlobal,DeltaPot_C
    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModNumConst,   ONLY: cRadToDeg,cPi
    
    real,    intent(in) :: time
    
    real, allocatable   :: Coord_I(:), PlotState_IV(:,:)

    !grid parameters
    integer, parameter :: nDim =1, nVar=4, S_=1
    !variable parameters
    integer, parameter :: n_=1, nflux_=2, qe_=3, pot_=4
    
    character(len=100),parameter :: NamePlotVar='Alt[km] n_se[cm-3] nflux[cm-2s-1] HeatingRate[eV/cm3/s] Pot[eV] g r'
    character(len=100) :: NamePlot='integrated.out'
    
    character(len=*),parameter :: NameHeader='SE output iono'
    character(len=5) :: TypePlot='ascii'
    integer :: iAlt,iEnergy

    logical :: IsFirstCall1=.true.,IsFirstCall2=.true.
    !--------------------------------------------------------------------------
    
    allocate(Coord_I(nAltExtended),PlotState_IV(nAltExtended,nVar))
    
    !    do iLine=1,nLine
    PlotState_IV = 0.0
    Coord_I     = 0.0
    
    !Set values
    do iAlt=1,nAltExtended
       Coord_I(iAlt) = AltExtended_C(iAlt)/1e5
       ! set plot state
       PlotState_IV(iAlt,n_)  = NumberDens_C(iAlt)
       PlotState_IV(iAlt,nflux_)  = NumberFlux_C(iAlt)
       PlotState_IV(iAlt,qe_)  = HeatingRate_C(iAlt)
       PlotState_IV(iAlt,Pot_)  = DeltaPot_C(iAlt)
       
    enddo
    
    ! set name for plotfile
    write(NamePlot,"(a,i4.4,a)") 'se_integrated_',iLineGlobal,'.out'
  
    !Plot grid for given line
    if(IsFirstCall1) then
       call save_plot_file(NamePlot, TypePositionIn='rewind', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn=0,TimeIn=time,     &
            nDimIn=nDim,CoordIn_I=Coord_I,                &
            VarIn_IV = PlotState_IV, ParamIn_I = (/1.6, 1.0/))
       IsFirstCall1 = .false.
    else
       call save_plot_file(NamePlot, TypePositionIn='append', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn=0,TimeIn=time,     &
            nDimIn=nDim,CoordIn_I=Coord_I,                &
            VarIn_IV = PlotState_IV, ParamIn_I = (/1.6, 1.0/))
    endif
    
    deallocate(Coord_I, PlotState_IV)
  end subroutine plot_integrated

  !=============================================================================
  
  subroutine calc_integrated_values(AVMU,uFlux_IC,dFlux_IC)
    use ModSeGrid,     ONLY: AltExtended_C,Alt_C, nAlt,nAltExtended, nEnergy,&
         EnergyGrid_I,DeltaE_I    
    use ModNumConst,    ONLY: cPi
    real,    intent(in) :: AVMU, uFlux_IC(nEnergy,nAlt),&
         dFlux_IC(nEnergy,nAlt)
    real,allocatable :: NumDensIntegrand_I(:) !integrand of number density
    integer :: iEnergy,iAlt


    !--------------------------------------------------------------------------
    !\
    ! calculate number density
    !/
    if (.not.allocated(NumDensIntegrand_I))allocate(NumDensIntegrand_I(nEnergy))
    if (.not.allocated(NumberDens_C)) allocate(NumberDens_C(nAltExtended))
    if (.not.allocated(NumberFlux_C))allocate(NumberFlux_C(nAltExtended))
    
    ALT: do iAlt=1,nAlt
       ENERGY: do iEnergy=1,nEnergy
          NumDensIntegrand_I(iEnergy) = &
               (uFlux_IC(iEnergy,iAlt)/AVMU+dFlux_IC(iEnergy,iAlt)/AVMU)&
               /sqrt(EnergyGrid_I(iEnergy))
       enddo ENERGY
       CALL midpnt_int(NumberDens_C(iAlt),NumDensIntegrand_I,&
            EnergyGrid_I,1,nEnergy,nEnergy,1)
       NumberDens_C(iAlt)=4.*cPi*1.7E-8*NumberDens_C(iAlt)
       
       ! calculate the number flux
!       CALL midpnt_int(NumberFlux_C(iAlt),&
!            uFlux_IC(:,iAlt)-dFlux_IC(:,iAlt),&
!            DeltaE_I,1,nEnergy,nEnergy,2)
       NumberFlux_C(iAlt)=2.0*cPi*sum((uFlux_IC(:,iAlt)-dFlux_IC(:,iAlt))*DeltaE_I(:))
       
    enddo ALT

    ! Map flux above calculation using Liouville's theorem
!    call map_flux(AVMU,uFlux_IC-dFlux_IC)
    call map_flux(AVMU,uFlux_IC)
!    do iAlt=1,nAltExtended
!       write(*,*) 'Alt_C(iAlt)*1.0e-5,NumberDens_C(iAlt),NumberFlux_C(iAlt)',&
!            AltExtended_C(iAlt)*1.0e-5,NumberDens_C(iAlt),NumberFlux_C(iAlt)
!    enddo

  end subroutine calc_integrated_values


  !=============================================================================
  
  subroutine map_flux(AVMU,uFlux_IC)
    use ModSeGrid,     ONLY: nAlt, nAltExtended,nEnergy,&
         EnergyGrid_I,DeltaE_I,AltExtended_C,DeltaPot_C    
    use ModNumConst,    ONLY: cPi
    
    real,    intent(in) :: AVMU, uFlux_IC(nEnergy,nAlt)
    integer :: iAlt,iEnergy,iPA, nPAforInt
    real :: B0,Biono
    integer, parameter :: nPA=90 ! number of PA
    real :: PA_I(nPa) !array for PA
    real :: mu0_I(nPA), mu_I(nPA),dmu0_I(nPA)
    real :: dPA,alpha,IntMu0,IntdMu0,UpFlux_I(nEnergy),DnFlux_I(nEnergy)
    real :: UpOmni_I(nEnergy),DnOmni_I(nEnergy)
    real :: KE_I(nEnergy),dKE_I(nEnergy)
    real,allocatable :: NumDensIntegrand_I(:) !integrand of number density
    real,allocatable :: NumFluxIntegrand_I(:)!integrand of number flux
    real, parameter :: cCmToM=1e-2
    logical :: DoReflect
    integer :: nReflect
    !-----------------------------------------------------------------------

    if (.not.allocated(NumDensIntegrand_I))allocate(NumDensIntegrand_I(nEnergy))
    if (.not.allocated(NumFluxIntegrand_I))allocate(NumFluxIntegrand_I(nEnergy))
    ! fill PA_I
    dPA = 0.5*cPi/real(nPA)
    PA_I(1) = 0.0
    mu_I(1) = cos(PA_I(1))
    do iPA =2,nPA
       PA_I(iPA) = PA_I(iPA-1)+dPA
       mu_I(iPA) = cos(PA_I(iPA))
    enddo
    
    DeltaPot_C(:)=0.0
!    DeltaPot_C(:)=0.5*EnergyGrid_I(1)

    !mu = sqrt((mu0^2-1) * (B(s)[E-e(dPhi(s)-dPhi0)])/(B0*E)+1)
    
    !set B at top of iono
    call get_b(AltExtended_C(nAlt)*cCmToM,Biono)
    ALT_LOOP: do iAlt=nAlt+1,nAltExtended
       !for a given pitchangle at iAlt find coresponding PA at nAlt 
       call get_b(AltExtended_C(iAlt)*cCmToM,B0)
       !B0=Biono
       ENERGY_LOOP: do iEnergy=1,nEnergy
          nPAforInt=0
          dmu0_I(:)=0.0
          mu0_I(:)=0.0
          nReflect=0
          PA_LOOP: do iPA =1,nPA
             !mu0=sqrt(1-(B0*E/(B(s)[E-e(dPfhi(s)-dPhi0)])) (1-mu^2)))
             !mu0=sqrt(1-(B0*(B(s))) (1-mu^2)))
             !mu0=sqrt(1-alpha)
             if (EnergyGrid_I(iEnergy)-DeltaPot_C(iAlt)>0)then
                alpha = (1.0-mu_I(iPA)**2.0)*(B0*EnergyGrid_I(iEnergy))&
                     /(Biono*(EnergyGrid_I(iEnergy)-DeltaPot_C(iAlt)))
             else
                alpha=2.0
             endif
             if ((1.0-alpha)>0.0) then
                mu0_I(iPA) = sqrt(1.0-alpha)
                nPAforInt=nPAforInt+1
                !save dmu0 spacing
                if (iPA>1) dmu0_I(iPA-1)=mu0_I(iPA-1)-mu0_I(iPA)
             else
                mu0_I(iPA) = -99.0 
             endif
             !write(*,*) 'mu_I(iPA),mu0_I(iPA)',mu_I(iPA),mu0_I(iPA)
             !check if this pitch angle and energy will reflect
             call check_will_reflect(iEnergy,iAlt,mu_I(iPA),DoReflect)
 !            write(*,*) 'EnergyGrid_I(iEnergy),mu_I(iPA),DoReflect',&
 !                 EnergyGrid_I(iEnergy),mu_I(iPA),DoReflect
             if (DoReflect) then
                nReflect=nReflect+1
             endif
          enddo PA_LOOP
         
          !fill in last delta for mu
          if (nPAforInt>1) dmu0_I(nPAforInt)=0.0!dmu0_I(nPAforInt-1)
          !set the Kinetic Energy array for this altitude
          KE_I(iEnergy) = EnergyGrid_I(iEnergy)-DeltaPot_C(iAlt)
          if (KE_I(iEnergy) <0.0 .or. nPAforInt<=1) then
             !when KE is less than 0 then set KE to 0 and there is no 
             ! contribution to the flux
             KE_I(iEnergy)=0.0
             dKE_I(iEnergy)=0.0
             UpFlux_I(iEnergy) = 0.0
             DnFlux_I(iEnergy) = 0.0
             NumDensIntegrand_I(iEnergy) = 0.0
             NumFluxIntegrand_I(iEnergy) = 0.0
             cycle ENERGY_LOOP

          else
             !get mu0 integral
!             CALL midpnt_int(IntMu0,&
!                  mu0_I(1:nPAforInt),dmu0_I(1:nPAforInt),1,nPAforInt,&
!                  nPAforInt,2)
             IntMu0=sum(mu0_I(1:nPAforInt)*dmu0_I(1:nPAforInt))
             IntdMu0=sum(dmu0_I(1:nPAforInt))
             UpOmni_I(iEnergy) = IntdMu0*uFlux_IC(iEnergy,nAlt)/AVMU
             UpFlux_I(iEnergy) = IntMu0*uFlux_IC(iEnergy,nAlt)/AVMU

             !set downward flux due to reflection and PA change
             if(nReflect==nPA) then
                 DnFlux_I(iEnergy) = 0.0
                 DnOmni_I(iEnergy) = 0.0
              else
                 !nPAforInt=nPA-nReflect
                 if(nPA-nReflect>=nPAforInt) then
                    DnFlux_I(iEnergy) = 0.0
                    DnOmni_I(iEnergy) = 0.0
                 else
                    IntMu0=sum(mu0_I(nPA-nReflect:nPAforInt)&
                         *dmu0_I(nPA-nReflect:nPAforInt))
                    IntdMu0=sum(dmu0_I(nPA-nReflect:nPAforInt))
                    DnOmni_I(iEnergy) = &
                         IntdMu0*(uFlux_IC(iEnergy,nAlt)/AVMU)
                    DnFlux_I(iEnergy) = &
                         IntMu0*(uFlux_IC(iEnergy,nAlt)/AVMU)
                 endif
              endif
          endif
 !         write(*,*) 'IntMu0',IntMu0
 !         write(*,*) 'nPAforInt',nPAforInt
          NumDensIntegrand_I(iEnergy) = &
               (UpOmni_I(iEnergy)+DnOmni_I(iEnergy))&
               /sqrt(KE_I(iEnergy))
          NumFluxIntegrand_I(iEnergy) = &
               (UpFlux_I(iEnergy)-DnFlux_I(iEnergy))

          !set dKE_I
          if (iEnergy>1.and.iEnergy<nEnergy) &
               dKE_I(iEnergy-1)=KE_I(iEnergy)-KE_I(iEnergy-1)
       enddo ENERGY_LOOP
       dKE_I(nEnergy)=dKE_I(nEnergy-1)

       !now integrate quantities for density and number flux
       NumberDens_C(iAlt)=4.0*cPi*1.7E-8*sum(NumDensIntegrand_I(:)*dKE_I(:))
       NumberFlux_C(iAlt)=2.0*cPi*sum(NumFluxIntegrand_I(:)*dKE_I(:))

    enddo ALT_LOOP
    
  end subroutine map_flux

  !=============================================================================
  ! check if a given energy and PA will reflect above iAltIn
  subroutine check_will_reflect(iEnergyIn,iAltIn,MuIn,DoReflect)
    use ModSeGrid,     ONLY: nAlt, nAltExtended,nEnergy,&
         EnergyGrid_I,DeltaE_I,AltExtended_C,DeltaPot_C    
    use ModNumConst,    ONLY: cPi
    
    integer,    intent(in) :: iEnergyIn,iAltIn
    real,    intent(in) :: MuIn
    logical, intent(out):: DoReflect
    integer :: iAlt,iEnergy,iPA, nPAforInt
    real :: B0,Biono
    integer, parameter :: nPA=90 ! number of PA
    real :: PA_I(nPa) !array for PA
    real :: mu0_I(nPA), mu_I(nPA),dmu0_I(nPA),alpha
    real, parameter :: cCmToM=1e-2
    !-----------------------------------------------------------------------
    DoReflect=.false.

    !mu = sqrt((mu0^2-1) * (B(s)[E-e(dPhi(s)-dPhi0)])/(B0*E)+1)
    !set B at iono point
    call get_b(AltExtended_C(iAltIn)*cCmToM,Biono)
    ALT_LOOP: do iAlt=iAltIn+1,nAltExtended
       !for a given pitchangle at iAlt find coresponding PA at nAlt 
       call get_b(AltExtended_C(iAlt)*cCmToM,B0)
       !B0=Biono
       !mu0=sqrt(1-(B0*E/(B(s)[E-e(dPfhi(s)-dPhi0)])) (1-mu^2)))
       !mu0=sqrt(1-(B0*(B(s))) (1-mu^2)))
       !mu0=sqrt(1-alpha)
       if (EnergyGrid_I(iEnergyIn)-DeltaPot_C(iAlt)>0)then
          alpha = (1.0-MuIn**2.0)*(B0*EnergyGrid_I(iEnergyIn))&
               /(Biono*(EnergyGrid_I(iEnergyIn)-DeltaPot_C(iAlt)))
       else
          alpha=2.0
          DoReflect=.true.
          return
       endif
       
       if ((1.0-alpha)> 0.0) then
          DoReflect=.false.
       else
          DoReflect=.true.
          return
       endif
    enddo ALT_LOOP
     
    
  end subroutine check_will_reflect

  !============================================================================
  subroutine get_b(AltRef, B0ref)
    
    use ModPlanetConst,     ONLY: Planet_,DipoleStrengthPlanet_I,rPlanet_I
    use ModNumConst,        ONLY: cDegToRad
    use ModSeBackground,ONLY:mLat
    real, intent(in):: AltRef !incomming reference alt [m]
    real, intent(out)   :: B0ref

    real, parameter :: cCmToM=1.0e-2
    real    :: Lshell, rPlanet, dipmom, Lat
    real    :: rRef ! reference radius
    !--------------------------------------------------------------------------

    rPlanet = rPlanet_I(Planet_)                          ! planet's radius (m)
    dipmom  = abs(DipoleStrengthPlanet_I(Planet_)*rPlanet**3)! planet's dipole 
    
    !set the reference radius
    rRef = rPlanet+AltRef
    ! find corresponding l-shell
    Lshell = 1.0/(cos(mLat*cDegToRad))**2.0
    
    ! find corresponding latitude for location on l-shell
    Lat = acos(sqrt(rRef/(Lshell*rPlanet)))
    
    ! get the magnetic field of the reference altitude
    B0ref = &
         dipmom*sqrt(1+3.0*(sin(Lat))**2.0)/(rRef)**3.0

  end subroutine get_b

  ! Subroutine MAXT
  !
  ! This software is part of the GLOW model.  Use is governed by the Open Source
  ! Academic Research License Agreement contained in the file glowlicense.txt.
  ! For more information see the file glow.txt.
  !
  ! Stan Solomon, 11/89, 9/91, 1/94, 3/05
  !
  ! Generates Maxwellian electron spectra with, optionally, a low energy tail
  ! of the form used by Meier et al., JGR 94, 13541, 1989.
  !
  ! Supplied by calling routine:
  !     EFLUX  total energy flux in erg cm-2 s-1
  !     EZER   characteristic energy in eV
  !     ENER   energy grid in eV
  !     dE   energy bin width in eV
  !     NBINS  number of energy bins (dimension of ENER, dE, and PHI)
  !     ITAIL  1 = Maxwellian with low-energy tail, 0 = regular Maxwellian
  !     FMONO  additional monoenergetic energy flux in erg cm-2 s-1
  !     EMONO  characteristic enerngy of FMONO in eV
  !
  ! Returned by subroutine:
  !     MAXT   Hemispherical flux in cm-2 s-1 eV-1
  !
  Subroutine maxt(EFLUX, EZER, ENER, dE, ITAIL, FMONO, EMONO, phi)
    !      use cglow, only: nbins
    use ModSeGrid, only: NBINS=>nEnergy
    implicit none

    !Args:
    Real,Intent(Out) :: phi(NBINS)
    Real, Intent(In)  :: EFLUX, EZER, ENER(NBINS), dE(NBINS), &
         FMONO,EMONO
    Integer,Intent(In):: ITAIL
    !Local:
    Real B,TE, PHIMAX, ERAT
    Integer K
    
    TE = 0.
    !
    IF (EZER < 500.) THEN
       B = 0.8*EZER
    ELSE
       B = 0.1*EZER + 350.
    ENDIF
    !
    PHIMAX = EXP(-1.)
    !
    DO K=1,NBINS
       ERAT = ENER(K) / EZER
       IF (ERAT > 60.) ERAT = 60.
       PHI(K) = ERAT * EXP(-ERAT)
       IF (ITAIL > 0) &
            PHI(K) = PHI(K) + 0.4*PHIMAX*(EZER/ENER(K))*EXP(-ENER(K)/B)
       TE = TE + PHI(K) * dE(K) * ENER(K) * 1.6022E-12
    enddo
    !
    DO K=1,NBINS
       PHI(K) = PHI(K) * EFLUX / TE
    enddo
    !
    !
    !
    if (fmono > 0.) then
       do  k=1,nbins
          if (emono > ener(k)-dE(k)/2..and.emono < ener(k)+dE(k)/2.) &
               phi(k)=phi(k)+fmono/(1.6022E-12*dE(k)*ener(k))
       enddo
    endif
  END Subroutine maxt
  
  
end Module ModElecTrans
