Module ModSeProduction
  ! This module contains all of the converted awfulness 
  save
  private !except
  public :: RCOLUM
  public :: RCOLUM_ABOVE
  public :: ephoto
  public :: SOLZEN
  public :: SSFLUX
  public :: init_production
  public :: unit_test_plas_resonant_scattering

  !number of wavelengthincrements
  integer, parameter :: LMAX=123
  !number of excited states
!  integer, parameter :: NEI=10
  
  !number of major neutral species
  integer:: nNeutral

  ! WAVE1   Array of minimum boundaries for radiation intervals; Ang
  ! WAVE2   Array of maximum boundaries for radiation intervals: Ang
  ! SFLUX   Array of solar flux intensities in the intervals defined 
  real  :: WAVE1(LMAX),WAVE2(LMAX),SFLUX(LMAX)

  ! store the slant column density
  real, allocatable ::ZCOL(:,:)

  ! store the column density above for starlight input calculation
  real, allocatable ::ZCOLabove(:,:)

  ! Store the crossections
  real, allocatable :: SIGION(:,:),SIGABS(:,:)

  ! Branching ratios
  real, allocatable :: PROB(:,:,:)
  
  !atomic mass of neutrals
  real,allocatable :: NeutralMassAMU_I(:)

  !Auger e- energies and thresholds
  real,    allocatable :: AugEnergy_I(:)
  real, allocatable :: AugWave_I(:)

  !surface gravity in cm/2^2
  real :: gSurface = 0.0
  
  integer, allocatable :: nStatesPerSpecies_I(:)
  integer :: nStatesMax 

  ! Store the ionization potentials
  real,public,allocatable :: TPOT(:,:)

  ! named constants for ion indicies in photoionrate array
  !Jupiter
  integer, parameter :: H2plus_=1,Heplus_=2,Hplus_=3,CH4plus_=4,&
                        CH3plus_=5,CH2plus_=6,CHplus_=7
  
  !named constants for neutrals
  !Earth
  integer, parameter :: O_=1, O2_=2, N2_=3
  !Jupiter
  integer, parameter :: H2_=1, He_=2, H_=3, CH4_=4
  
contains

  ! Subroutine RCOLUM
  !
  ! This software is part of the GLOW model.  Use is governed by the Open Source
  ! Academic Research License Agreement contained in the file glowlicense.txt.
  ! For more information see the file glow.txt.
  !
  ! Stan Solomon, 1988, 1991
  !
  ! Calculates the column density ZCOL for each species ZMAJ above height
  ! ZZ at zenith angle CHI.  Uses a Chapman function fit [Smith and Smith,
  ! JGR 77, 3592, 1972].  If CHI is less than 90 degrees, column
  ! densities are calculated directly; if CHI is greater than 90 degrees
  ! the column density at grazing height for 90 degrees is calculated and
  ! doubled and the column density above ZZ(J) is subtracted; if CHI is
  ! such that the grazing height is less than the radius of the earth the
  ! column densities are set to 'infinity', i.e., 1.0E30.  Densities
  ! supplied in array ZMAJ are used in the calculation except where
  ! grazing height is below the lowest level specified, in this case
  ! values are interpolated logarithmically from a US Standard Atmosphere
  ! at sea level, the tropopause, the stratopause, and the mesopause.
  SUBROUTINE RCOLUM (CHI, ZZ, ZMAJ, TN, IONO)

    use ModSeGrid,      ONLY: rPlanetCM
    use ModNumConst,    ONLY: cPi
    PARAMETER (NM=3) 
    PARAMETER (NU=4) !heights in us standard model
    integer, intent(in) :: IONO
    real,    intent(in) :: CHI, ZMAJ(nNeutral,IONO)
    !
    DIMENSION ZZ(IONO), TN(IONO)
    real, allocatable :: ZVCD(:,:),ZCG(:)
    
    ! *** planet specific
    DIMENSION  ZUS(NU), TNUS(NU), ZCUS(NM,NU)

    !

    ! These are the heights, temperatures and densities for US standard model
    DATA ZUS/0., 1.5E6, 5.E6, 9.E6/, TNUS/288., 217., 271., 187./
    DATA ZCUS/8.00E17, 4.54E24, 1.69E25, &
         8.00E17, 5.46E23, 2.03E24, &
         8.00E17, 3.63E21, 1.35E22, &
         7.80E17, 8.48E18, 3.16E19/
    !
    !
    if (.not.allocated(ZVCD)) allocate(ZVCD(nNeutral,IONO))
    if (.not.allocated(ZCG)) allocate(ZCG(nNeutral))

    if (.not.allocated(ZCOL)) allocate(ZCOL(nNeutral,IONO))
    CALL VCD (ZZ, ZMAJ, ZVCD, IONO, nNeutral)
    !
    IF (CHI .GE. 2.) THEN
       do I=1,nNeutral
          do J=1,IONO
             ZCOL(I,J) = 1.0E30
          enddo
       enddo
       RETURN
    ENDIF
    !
    IF (CHI .LE. cPi/2.) THEN
       DO I=1,nNeutral
          DO J=1,IONO
             ZCOL(I,J) = ZVCD(I,J) * CHAP(CHI,ZZ(J),TN(J),I)
          enddo
       enddo
    ELSE
       do J=1,IONO
          GHRG=(rPlanetCM+ZZ(J))*SIN(CHI)
          GHZ=GHRG-rPlanetCM
          IF (GHZ .LE. 0.) THEN
             do I=1,nNeutral
                ZCOL(I,J) = 1.0E30
             end do
             cycle
          ENDIF
          IF (GHZ .GE. ZZ(1)) THEN
             DO JG=1,J-1
                IF (ZZ(JG) .LE. GHZ .AND. ZZ(JG+1) .GT. GHZ) GOTO 120
             enddo
120          TNG = TN(JG)+(TN(JG+1)-TN(JG))*(GHZ-ZZ(JG))/(ZZ(JG+1)-ZZ(JG))
             do I=1,nNeutral
                ZCG(I) = ZVCD(I,JG) * (ZVCD(I,JG+1) / ZVCD(I,JG)) ** &
                     ((GHZ-ZZ(JG)) / (ZZ(JG+1)-ZZ(JG)))
             enddo
          ELSE
             !Here the grazing altitude is less than the bottom of the model 
             !therefore the values are interpolated from US standard atmosphere 
             !at sea level. Only good for Earth, Not suitable for Jupiter!!!
             ! for Jupiter we have densities down to 0 altitude so just use 
             ! GITM values. For jupiter zz should go to 0 alt.
             do JG=1,3
                IF (ZUS(JG) .LT. GHZ .AND. ZUS(JG+1) .GT. GHZ) GOTO 180
             enddo
180          TNG = TNUS(JG) &
                  + (TNUS(JG+1)-TNUS(JG))*(GHZ-ZUS(JG))/(ZUS(JG+1)-ZUS(JG))
             do I=1,nNeutral
                ZCG(I) = ZCUS(I,JG) * (ZCUS(I,JG+1) / ZCUS(I,JG)) ** &
                     ((GHZ-ZUS(JG)) / (ZUS(JG+1)-ZUS(JG)))
             end do
          ENDIF
          do I=1,nNeutral
             ZCOL(I,J) = 2. * ZCG(I) * CHAP(CPI/2.,GHZ,TNG,I) &
                  - ZVCD(I,J) * CHAP(CHI,ZZ(J),TN(J),I)
          end do
       end do
    ENDIF
    !
    RETURN
  END SUBROUTINE RCOLUM

  !============================================================================
  !variant of rcolum to get column density above the location of interest
  SUBROUTINE RCOLUM_ABOVE (ZZ, ZMAJ, TN, IONO)

    use ModSeGrid,      ONLY: rPlanetCM
    use ModNumConst,    ONLY: cPi
    PARAMETER (NM=3) 
    PARAMETER (NU=4) !heights in us standard model
    integer, intent(in) :: IONO
    real,    intent(in) :: ZMAJ(nNeutral,IONO)
    !
    DIMENSION ZZ(IONO), TN(IONO)
    real, allocatable :: ZVCD(:,:),ZCG(:)
    real,parameter ::CHI=0.0
    ! *** planet specific
    DIMENSION  ZUS(NU), TNUS(NU), ZCUS(NM,NU)

    !

    ! These are the heights, temperatures and densities for US standard model
    DATA ZUS/0., 1.5E6, 5.E6, 9.E6/, TNUS/288., 217., 271., 187./
    DATA ZCUS/8.00E17, 4.54E24, 1.69E25, &
         8.00E17, 5.46E23, 2.03E24, &
         8.00E17, 3.63E21, 1.35E22, &
         7.80E17, 8.48E18, 3.16E19/
    !
    !
    if (.not.allocated(ZVCD)) allocate(ZVCD(nNeutral,IONO))
    if (.not.allocated(ZCG)) allocate(ZCG(nNeutral))

    if (.not.allocated(ZCOLabove)) allocate(ZCOLabove(nNeutral,IONO))
    CALL VCD (ZZ, ZMAJ, ZVCD, IONO, nNeutral)
    !
    IF (CHI .GE. 2.) THEN
       do I=1,nNeutral
          do J=1,IONO
             ZCOLabove(I,J) = 1.0E30
          enddo
       enddo
       RETURN
    ENDIF
    !
    IF (CHI .LE. cPi/2.) THEN
       DO I=1,nNeutral
          DO J=1,IONO
             ZCOLabove(I,J) = ZVCD(I,J) * CHAP(CHI,ZZ(J),TN(J),I)
          enddo
       enddo
    ELSE
       do J=1,IONO
          GHRG=(rPlanetCM+ZZ(J))*SIN(CHI)
          GHZ=GHRG-rPlanetCM
          IF (GHZ .LE. 0.) THEN
             do I=1,nNeutral
                ZCOLabove(I,J) = 1.0E30
             end do
             cycle
          ENDIF
          IF (GHZ .GE. ZZ(1)) THEN
             DO JG=1,J-1
                IF (ZZ(JG) .LE. GHZ .AND. ZZ(JG+1) .GT. GHZ) GOTO 120
             enddo
120          TNG = TN(JG)+(TN(JG+1)-TN(JG))*(GHZ-ZZ(JG))/(ZZ(JG+1)-ZZ(JG))
             do I=1,nNeutral
                ZCG(I) = ZVCD(I,JG) * (ZVCD(I,JG+1) / ZVCD(I,JG)) ** &
                     ((GHZ-ZZ(JG)) / (ZZ(JG+1)-ZZ(JG)))
             enddo
          ELSE
             !Here the grazing altitude is less than the bottom of the model 
             !therefore the values are interpolated from US standard atmosphere 
             !at sea level. Only good for Earth, Not suitable for Jupiter!!!
             ! for Jupiter we have densities down to 0 altitude so just use 
             ! GITM values. For jupiter zz should go to 0 alt.
             do JG=1,3
                IF (ZUS(JG) .LT. GHZ .AND. ZUS(JG+1) .GT. GHZ) GOTO 180
             enddo
180          TNG = TNUS(JG) &
                  + (TNUS(JG+1)-TNUS(JG))*(GHZ-ZUS(JG))/(ZUS(JG+1)-ZUS(JG))
             do I=1,nNeutral
                ZCG(I) = ZCUS(I,JG) * (ZCUS(I,JG+1) / ZCUS(I,JG)) ** &
                     ((GHZ-ZUS(JG)) / (ZUS(JG+1)-ZUS(JG)))
             end do
          ENDIF
          do I=1,nNeutral
             ZCOLabove(I,J) = 2. * ZCG(I) * CHAP(CPI/2.,GHZ,TNG,I) &
                  - ZVCD(I,J) * CHAP(CHI,ZZ(J),TN(J),I)
          end do
       end do
    ENDIF
    !
    RETURN
  END SUBROUTINE RCOLUM_ABOVE
  !============================================================================
  !
  !
  !
  !
  FUNCTION CHAP (CHI, Z, T, I)
    use ModSeGrid,      ONLY: rPlanetCM
    use ModNumConst,    ONLY: cPi
    GR=gSurface*(rPlanetCM/(rPlanetCM+Z))**2
    HN=1.38E-16*T/(NeutralMassAMU_I(I)*1.662E-24*GR)
    HG=(rPlanetCM+Z)/HN
    HF=0.5*HG*(COS(CHI)**2)
    SQHF=SQRT(HF)
    CHAP=SQRT(0.5*cPi*HG)*SPERFC(SQHF)
    RETURN
  END FUNCTION CHAP
  !
  !
  !
  !
  FUNCTION SPERFC(DUMMY)
! error function erfc(y) = 1 - erf(y) 
! from Smith & Smith [1972], eq 12
    IF (DUMMY .LE. 8.) THEN
       SPERFC = (1.0606963+0.55643831*DUMMY) / &
            (1.0619896+1.7245609*DUMMY+DUMMY*DUMMY)
    ELSE
       SPERFC=0.56498823/(0.06651874+DUMMY)
    ENDIF
    RETURN
  END FUNCTION SPERFC
  !
  !
  !
  !
  SUBROUTINE VCD(ZZ,ZMAJ,ZVCD,IONO,nNeutralIn)
    DIMENSION ZZ(IONO), ZMAJ(nNeutralIn,IONO), ZVCD(nNeutralIn,IONO)
    !
    DO I=1,nNeutralIn
       ZVCD(I,IONO) =   ZMAJ(I,IONO) &
            * (ZZ(IONO)-ZZ(IONO-1)) &
            / ALOG(ZMAJ(I,IONO-1)/ZMAJ(I,IONO))

       DO J=IONO-1,1,-1
          RAT = ZMAJ(I,J+1) / ZMAJ(I,J)
          ZVCD(I,J) =   ZVCD(I,J+1) &
               + ZMAJ(I,J) * (ZZ(J)-ZZ(J+1)) / ALOG(RAT) * (1.-RAT)
       enddo
    end DO
    RETURN
  END SUBROUTINE VCD


!
!
! Subroutine EPHOTO
!
! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.
!
! Adapted from Banks & Nagy 2-stream input code by Stan Solomon, 6/88
! Modified to handle Auger electrons, Stan Solomon, 7/90
! Reads cross sectons from files (for 1-nm bins), Scott Bailey, ~1994
! Modified bin structure, fixed CIII problem, Stan Solomon, 12/2000
! Corrected additional Auger problem, Liying Qian, 11/2002
! Converged above three branches, Stan Solomon, 3/2005
! Removed LIMIN, wavelength loop now runs from 1 to LMAX, SCS, 3/2005
!
! This subroutine calculates photoionization, rates, certain
! photodissociative excitation rates, and the photoelectron production
! spectrum as a function of altitude.  Uses continuously variable energy
! grid.  3 major species: O, O2, N2; NO is treated as a minor (non-
! absorbing) specie.
!
! Modified for multiple planets and to separate file reading - Alex Glocer2016
!  For Earth - 3 major neutral species : O, O2, N2
! NO is no longer tracked
! For Jupiter/Saturn - 4 major neutral species : H2, H, CH4, He

!
! Supplied by calling routine:
! WAVE1   wavelength array, upper bound; Angstroms
! WAVE2   wavelength array, lower bound; Angstroms
! SFLUX   solar flux array; photons cm-2 sec-1
! ZZ      altitude array; cm above earth
! ZMAJ    density array for species O, O2, N2, altitude; cm-3
! ZNO     density of NO at each altitude; cm-3
! ZCOL    slant column density for species O, O2, N2, altitude; cm-2
! ENER    energy grid for photoelectrons; eV
! DEL     array of energy grid increments; eV
!
! Calculated by subroutine:
! PESPEC  photoelectron production spectrum for each altitude; cm-3 s-1
! PHOTOI  photoionization rates for state, species, altitude; cm-3 s-1
! PHOTOD  photodissoc./exc. rates for state, species, alt.; cm-3 s-1
! PHONO   photoionization/dissoc./exc. rates for NO; cm-3 s-1
!
! Other definitions:
! DSPECT  ionization rate in particular wavelength bin; cm-3 s-1
! TAU     optical depth, dimensionless
! FLUX    solar flux at altitude; cm-2 s-1
! SIGABS  photoabsorption cross sections, O, O2, N2; cm2
! SIGION  photoionization cross sections, O, O2, N2; cm2
! SIGAO, SIGAO2, SIGAN2, SIGIO, SIGIO2, SIGIN2; cross sect. data arrays
! NNN     number of states for each species
! TPOT    ionization potentials for each species, state; eV
! PROB    branching ratios for each state, species, and wavelength bin:
!         O+ states: 4S, 2Do, 2Po, 4Pe, 2Pe
!         O2+ states: X, a+A, b, dissoc.
!         N2+ states: X, A, B, C, F, dissoc.
! PROBO, PROBO2, PROBN2; branching ratio data arrays
! BSO2    yield of O(1S) from dissociation of O2
! EPSIL1  energy loss lower bound for state, species, wavelength; eV
! EPSIL2  energy loss upper bound for state, species, wavelength; eV
! SIGNO   NO photoionization xsect at Ly-alpha
! AUGE    Mean energy of Auger electrons for each species; eV
! AUGL    Wavelength threshold for Auger electrons; Angstroms
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
  SUBROUTINE EPHOTO(ZMAJ,PESPEC,nAlt,Ioff,SZA,AltKm_C,nIons,PhotoIonRate_IC)

    !
    !
    
    use ModSeGrid,   only:nEnergy,del=>DeltaE_I,ener=>EnergyGrid_I,Emin=>EnergyMin
!    use ModMath,     only:midpnt_int
    use ModPlanetConst, ONLY: Planet_, NamePlanet_I
    use ModNumConst, only:cPi
!    PARAMETER (NEX=20)
!    PARAMETER (NW=20)
!    PARAMETER (NC=10)
!    PARAMETER (NF=4)

    integer, intent(in) :: nAlt
    real, intent(in) :: SZA, AltKm_C(nAlt)
    real :: FluxRes
    !Set named constants for particular wavelength bins
    integer,parameter :: LyAlpha_= 12, LyBeta_=18 , HeI_=62, HeII_=90
    !for starlight bin1 is 1000-1050A, bin2 is 950-1000A and bin3 is 900-950A
    integer,parameter :: Starlight1a_= 16,Starlight1b_= 20,&
         Starlight2a_= 21,Starlight2b_= 25,Starlight3a_=26,Starlight3b_=30
    
    !incident starlight intensity for starlight bins (Thitheridge 2000)
    real, parameter :: StarLightIntensity = 5.0e6 !photons/cm2/s

!    COMMON /CGLOW/ &
!         ZCOL(nNeutral,nAlt),WAVE1(LMAX),WAVE2(LMAX),SFLUX(LMAX)
    !
    !      COMMON /CENERGY/ ener(Elen),del(Elen),Emin,Jo
    !
    DIMENSION FLUX(LMAX,nAlt), &
         PESPEC(NEnergy,nAlt), PeSpectrumSpecies_IIIC(nStatesMax,nNeutral,NEnergy,nAlt), &
         ZMAJ(nNeutral,nAlt), PhotoIonRate_IC(nIons,nAlt),&
         PHOTOI(nStatesMax,nNeutral,nAlt), PHOTOD(nStatesMax,nNeutral,nAlt), &
         BSO2(LMAX),PHONO(nStatesMax,nAlt), RION(LMAX,nNeutral,nAlt)

    real, allocatable :: EPSIL1(:,:,:), EPSIL2(:,:,:)
    real :: DSPECT(nAlt)
    real    :: TAU
    integer :: iAlt

    !
    SAVE EPSIL1, EPSIL2, EPB1, EPB2
    !
    DATA  IFIRST/1/, LIMIN/16/ ! No PROBs below L=16
    !
    !

    !
    !DATA BSO2/12*0.,.01,.03,7*.10,8*.07,5*.03,5*.01,84*0./
    !
    !
    !
    ! First time only:  pack photoabsorption, photoioniation cross sections
    ! and convert to cm2; pack branching ratios; calculate energy losses:
    !
    !for some reason this is set to zero later so just do it now
    BSO2(:)=0.0


    IF (IFIRST .EQ. 1) THEN
       IFIRST = 0


       if (.not.allocated(EPSIL1)) allocate(EPSIL1(nStatesMax,nNeutral,LMAX))
       if (.not.allocated(EPSIL2)) allocate(EPSIL2(nStatesMax,nNeutral,LMAX))

       DO  L=1,LMAX
          DO  I=1,nNeutral
             DO  K=1,nStatesPerSpecies_I(I)
                !wavelength to energy - tpot
                EPSIL1(K,I,L)=12397.7/WAVE1(L)-TPOT(K,I)
                EPSIL2(K,I,L)=12397.7/WAVE2(L)-TPOT(K,I)
                IF (WAVE1(L) .LE. AugWave_I(I)) THEN
                   EPSIL1(K,I,L) = EPSIL1(K,I,L) - AugEnergy_I(I)
                   EPSIL2(K,I,L) = EPSIL2(K,I,L) - AugEnergy_I(I)
                endif
             enddo
          end do
       end do

    ENDIF

    !
    !
    ! Zero arrays:
    !
    DO J=1,nAlt
       DO I=1,nNeutral
          DO K=1,nStatesMax
             PHONO(K,J) = 0.
             PHOTOI(K,I,J) = 0.
             PHOTOD(K,I,J) = 0.
          end do
       end do
       DO  M=1,nEnergy
          PESPEC(M,J) = 0.
          PeSpectrumSpecies_IIIC(:,:,M,J)=0.
       end do
    end do

!    write(*,*) '!!! SZA',SZA



    !
    !
    ! Calculate attenuated solar flux at all altitudes and wavelengths:
    !
    DO  L=1,LMAX
       DO  J=1,nAlt
          TAU=0.
          DO  I=1,nNeutral
             TAU=TAU+SIGABS(I,L)*ZCOL(I,J)
          end do
          IF (TAU .LT. 20.) THEN
             FLUX(L,J)=SFLUX(L)*EXP(-TAU)
          ELSE
             FLUX(L,J) = 0.0
          ENDIF

          if(NamePlanet_I(Planet_)=='EARTH')then
             ! add in the resonant scattering from the plasmasphere from 
             ! strobel et al 1974, only for earth now
             if (L==LyAlpha_)then
                call get_plas_resonant_scattering(&
                     (/SZA,AltKm_C(J)/),'LyAlpha',FluxRes)
                FLUX(L,J) = FLUX(L,J) + FluxRes
             endif
             if (L==LyBeta_)then
                call get_plas_resonant_scattering(&
                     (/SZA,AltKm_C(J)/),'LyBeta',FluxRes)
                FLUX(L,J) = FLUX(L,J) + FluxRes
             endif
             if (L==HeI_)then
                call get_plas_resonant_scattering(&
                     (/SZA,AltKm_C(J)/),'HeI',FluxRes)
                FLUX(L,J) = FLUX(L,J) + FluxRes
             endif
             if (L==HeII_)then
                call get_plas_resonant_scattering(&
                     (/SZA,AltKm_C(J)/),'HeII',FluxRes)
                FLUX(L,J) = FLUX(L,J) + FluxRes
             endif
             
          endif

          ! add in the starlight source
          if ((Starlight1a_>=L .and. Starlight1b_<=L)&
               .or. (Starlight2a_>=L .and. Starlight2b_<=L)&
               .or. (Starlight3a_>=L .and. Starlight3b_<=L)) then

             ! recalculate tau for starlight which arrives over multiple 
             ! directions using approach of Titheridge 2000
             TAU=0.
             DO  I=1,nNeutral
                TAU=TAU+SIGABS(I,L)*ZCOLabove(I,J)
             end do
             FLUX(L,J) = FLUX(L,J) &
                  + StarLightIntensity * 0.25* (exp(-TAU/0.9) &
                  + exp(-TAU/0.7)+exp(-TAU/0.5)+exp(-TAU/0.3))
          endif


          !
          !
          ! Calculate SRC photodissociation of O2, dissociative excitation of
          ! O(1S), photodissociation of N2, and photoionization of NO by solar
          ! Ly-alpha:
          ! *** planet specific, but doesn't seem to be used currently
          if(NamePlanet_I(Planet_)=='EARTH')then
             IF (WAVE1(L) .LT. 1751. .AND. WAVE2(L) .GT. 1349.) &
                  PHOTOD(1,2,J) = PHOTOD(1,2,J)+ZMAJ(2,J)*SIGABS(2,L)*FLUX(L,J)
             PHOTOD(2,2,J) = PHOTOD(2,2,J) + ZMAJ(2,J)*SIGABS(2,L)*FLUX(L,J) &
                  * BSO2(L)
             PHOTOD(1,3,J) = PHOTOD(1,3,J) + &
                  ZMAJ(3,J)*(SIGABS(3,L)-SIGION(3,L))*FLUX(L,J)
             !for now leave out the NO
             !IF (WAVE1(L) .LT. 1221. .AND. WAVE2(L) .GT. 1209.) &
             !     PHONO(1,J) = PHONO(1,J) + ZNO(J)*SIGNO*FLUX(L,J)
          endif
!          IF (WAVE1(L) .LT. 1751 .AND. WAVE2(L) .GT. 1349.) &
!               PHOTOD(1,2,J) = PHOTOD(1,2,J)+ZMAJ(2,J)*SIGABS(2,L)*FLUX(L,J)
!          PHOTOD(2,2,J) = PHOTOD(2,2,J) + ZMAJ(2,J)*SIGABS(2,L)*FLUX(L,J) &
!               * BSO2(L)
!          PHOTOD(1,3,J) = PHOTOD(1,3,J) + &
!               ZMAJ(3,J)*(SIGABS(3,L)-SIGION(3,L))*FLUX(L,J)
          
       end do
    end do
    

    !
    !
!    Emax=ENER(nEnergy)+DEL(nEnergy)/2.
    !
    ! Calculate ionization rates and photoelectron production:
    !
    !
    ! Loop over ionizing wavelengths:
    !
    WAVE_LENGTH_LOOP: DO  L=LIMIN,LMAX
       !
       !
       ! Loop over species:
       !
       SPECIES_LOOP: DO  I=1,nNeutral
          !
          ! Calculate total ionization rates for all species and altitudes:
          !
          DO J=1,nAlt
             RION(L,I,J)=ZMAJ(I,J)*SIGION(I,L)*FLUX(L,J)
          enddo
          !
          ! Loop over states:
          !
          DO  K=1,nStatesPerSpecies_I(I)
             !
             E1= EPSIL1(K,I,L)
             E2= EPSIL2(K,I,L)
             IF (E2 .LT. 0.0) cycle
             IF (E1 .LT. 0.0) E1=0.0
             
             ! Calculate state-specific ionization rates at all altitudes:
             !
             DO J=1,nAlt
                DSPECT(J) = RION(L,I,J)*PROB(K,I,L)
                PHOTOI(K,I,J) = PHOTOI(K,I,J) + DSPECT(J)
             enddo

             !
             !
             ! Find box numbers M1, M2 corresponding to energies E1, E2:
             !
             CALL BOXNUM (E1, E2, M1, M2, R1, R2, DEL, ENER)
             IF (M1 .GT. nEnergy) cycle 
             
             !
             !
             ! Fill the boxes from M1 to M2 at all altitudes:
             !
             Y = E2 - E1
             DO  N=M1,M2
                IF (M1 .EQ. M2) THEN
                   FAC = 1.
                ELSE
                   IF (N .EQ. M1) THEN
                      FAC = (R1-E1) / Y
                   ELSE
                      IF (N .EQ. M2) THEN
                         FAC = (E2-R2) / Y
                      ELSE
                         FAC = DEL(N) / Y
                      ENDIF
                   ENDIF
                ENDIF
                DO  J=1,nAlt
                   PESPEC(N,J) = PESPEC(N,J) + DSPECT(J) * FAC
                   PeSpectrumSpecies_IIIC(K,I,N,J)=DSPECT(J) * FAC
                enddo
             enddo

             !
             !
             ! Generate Auger electrons if energy is sufficient:
             !
             IF (WAVE1(L) .LE. AugWave_I(I)) THEN
                E1 = AugEnergy_I(I)
                E2 = AugEnergy_I(I)
                CALL BOXNUM (E1, E2, M1, M2, R1, R2, DEL, ENER)
                IF (M1.GT.nEnergy .OR. M2.GT.nEnergy) CYCLE
                DO J=1,nAlt
                   PESPEC(M1,J) = PESPEC(M1,J) + RION(L,I,J)
                   PeSpectrumSpecies_IIIC(K,I,N,J)=DSPECT(J) * FAC
                enddo
             ENDIF
          end DO
       enddo SPECIES_LOOP
    enddo WAVE_LENGTH_LOOP
    
    !integrate now to get photoionization rates for each species
    select case(NamePlanet_I(Planet_))
    case('EARTH')
       !for now just set rate for Earth to zero
       PhotoIonRate_IC(1,:) = 0.0
    case('JUPITER')
       !uncomment for debugging
       !CALL plot_pespecspecies(PeSpectrumSpecies_IIIC)
       !CALL plot_flux(FLUX)
       !CALL plot_crossec(SIGION,PROB)
       
       do iAlt=1,nAlt
          !
          !integrate to get production rate for each species as a function of 
          ! altitude
          !H2+
          PhotoIonRate_IC(H2plus_,iAlt) = 4.0*cPi*&
               sum(PeSpectrumSpecies_IIIC(1,H2_,:,iAlt)*del(:))
          
          !He+
          PhotoIonRate_IC(Heplus_,iAlt) = 4.0*cPi*&
               sum(PeSpectrumSpecies_IIIC(1,He_,:,iAlt)*del(:))

          !H+
          PhotoIonRate_IC(Hplus_,iAlt) = 4.0*cPi*&
               (sum(PeSpectrumSpecies_IIIC(2,H2_,:,iAlt)*del(:)) &
               + sum(PeSpectrumSpecies_IIIC(1,H_,:,iAlt)*del(:)) &
               + sum(PeSpectrumSpecies_IIIC(5,CH4_,:,iAlt)*del(:)))

          !CH4+
          PhotoIonRate_IC(CH4plus_,iAlt) = 4.0*cPi*&
               sum(PeSpectrumSpecies_IIIC(1,CH4_,:,iAlt)*del(:))
          
          !CH3+
          PhotoIonRate_IC(CH3plus_,iAlt) = 4.0*cPi*&
               sum(PeSpectrumSpecies_IIIC(2,CH4_,:,iAlt)*del(:))
          
          !CH2+
          PhotoIonRate_IC(CH2plus_,iAlt) = 4.0*cPi*&
               sum(PeSpectrumSpecies_IIIC(3,CH4_,:,iAlt)*del(:))
          
          !CH+
          PhotoIonRate_IC(CHplus_,iAlt) = 4.0*cPi*&
               sum(PeSpectrumSpecies_IIIC(4,CH4_,:,iAlt)*del(:))

       enddo
    end select

  end SUBROUTINE EPHOTO

  !=============================================================================
  ! save ion production plot for verification
  ! doesn't work for Earth yet
  ! debugging subroutines to plot epec params. only for one fieldline
  subroutine plot_pespecspecies(PeSpectrumSpecies_IIIC)
    use ModSeGrid,     ONLY: Alt_C,nAlt,nEnergy, &
         DeltaE_I,EnergyGrid_I
    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModNumConst,   ONLY: cRadToDeg,cPi
    use ModPlanetConst, only: Planet_, NamePlanet_I

    real, intent(in) :: PeSpectrumSpecies_IIIC(nStatesMax,nNeutral,nEnergy,nAlt)

    real, allocatable   :: Coord_DII(:,:,:), PlotState_IIV(:,:,:)
    !grid parameters
    integer, parameter :: nDim =2,E_=1, S_=2
    integer :: nVar
    integer, parameter ::iLine=1

    !Jupiter
    integer, parameter :: H2plus_=1,Heplus_=2,Hplus_=3,CH4plus_=4,&
         CH3plus_=5,CH2plus_=6,CHplus_=7
    
    !Earth
    integer, parameter :: Oplus_=1
    
    character(len=100),parameter :: NamePlotVarEarth=&
         'Alt[km] O+[cm-3s-1] g r'
    character(len=100),parameter :: NamePlotVarJupiter=&
         'E[eV] Alt[km] He+[cm-3s-1] g r'

    character(len=100) :: NamePlotVar
    character(len=*),parameter :: NameHeader='Photoionization Rates'
    character(len=5) :: TypePlot='ascii'
    integer :: iIon,iAlt
    character(len=100) :: NamePlot
    logical,save :: IsFirstCall =.true.
    !--------------------------------------------------------------------------

!    nVar=nIons!+2
    nVar=1
    allocate(Coord_DII(nDim,nEnergy,nAlt),PlotState_IIV(nEnergy,nAlt,nVar))

    PlotState_IIV = 0.0
    Coord_DII     = 0.0
    
    select case(NamePlanet_I(Planet_))
    case('EARTH')
       NamePlotVar=NamePlotVarEarth
    case('JUPITER')
              NamePlotVar=NamePlotVarJupiter
    end select
       

       !Set Coordinates along field line and PA
       do iEnergy=1,nEnergy
          do iAlt=1,nAlt
             Coord_DII(E_,iEnergy,iAlt) = EnergyGrid_I(iEnergy)             
             Coord_DII(S_,iEnergy,iAlt) = Alt_C(iAlt)/1e5
             PlotState_IIV(iEnergy,iAlt,1)  = &
                  PeSpectrumSpecies_IIIC(1,He_,iEnergy,iAlt)
          enddo
       enddo       

    ! set name for plotfile
    write(NamePlot,"(a)") 'PeSpecHe.out'
    
    !Plot grid for given line
    if(IsFirstCall) then
       call save_plot_file(NamePlot, TypePositionIn='rewind', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
            nDimIn=nDim,CoordIn_DII=Coord_DII,                &
            VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/)) !***
       IsFirstCall = .false.
    else
       call save_plot_file(NamePlot, TypePositionIn='append', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
            nDimIn=nDim,CoordIn_DII=Coord_DII,                &
            VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/)) !***
    endif
    
    deallocate(Coord_DII, PlotState_IIV)
  end subroutine plot_pespecspecies

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  ! save solar flux plot for verification
  !     currently plots only the lower wavelength of each bin
  !     doesn't treat lines appropriately
  !     only for one line and debugging
  !--------------------------------------------------------------------------
  subroutine plot_flux(FLUX)
    use ModSeGrid,     ONLY: Alt_C,nAlt,nEnergy,  &
         DeltaE_I,EnergyGrid_I
    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModNumConst,   ONLY: cRadToDeg,cPi
    use ModPlanetConst, only: Planet_, NamePlanet_I

    real, intent(in) :: FLUX(LMAX,nAlt)

    real, parameter :: C1=12397.7       ! Converting wavelengths to energy

    real, allocatable   :: Coord_DII(:,:,:), PlotFlux_IIV(:,:,:)
    !grid parameters
    integer, parameter :: nDim =2,E_=1, S_=2
    integer :: nVar
    integer, parameter ::iLine=1
    !Jupiter
    integer, parameter :: H2plus_=1,Heplus_=2,Hplus_=3,CH4plus_=4,&
         CH3plus_=5,CH2plus_=6,CHplus_=7
    
    !Earth
    integer, parameter :: Oplus_=1
    
    character(len=100),parameter :: NamePlotVarEarth=&
         'Alt[km] O+[cm-3s-1] g r'
    character(len=100),parameter :: NamePlotVarJupiter=&
         'Energy[eV] Alt[km] PhotonFlux g r'

    character(len=100) :: NamePlotVar
    character(len=*),parameter :: NameHeader='Solar Flux'
    character(len=5) :: TypePlot='ascii'
    integer :: iIon,iAlt
    character(len=100) :: NamePlot
    logical,save :: IsFirstCall =.true.
    !--------------------------------------------------------------------------

    nVar=1
    allocate(Coord_DII(nDim,LMAX,nAlt),PlotFlux_IIV(LMAX,nAlt,nVar))

    PlotFlux_IIV = 0.0
    Coord_DII     = 0.0
    
    select case(NamePlanet_I(Planet_))
    case('EARTH')
       NamePlotVar=NamePlotVarEarth
    case('JUPITER')
       NamePlotVar=NamePlotVarJupiter
    end select
       
       !Set Flux Coordinates along field line and wavelength converted to eV
       do iLambda=1,LMAX
          do iAlt=1,nAlt
             Coord_DII(E_,iLambda,iAlt) = C1/Wave1(iLambda)
             Coord_DII(S_,iLambda,iAlt) = Alt_C(iAlt)/1e5
             PlotFlux_IIV(iLambda,iAlt,1)  = FLUX(iLambda,iAlt)
          enddo
       enddo       
       
    ! set name for plotfile
    write(NamePlot,"(a)") 'Flux.out'
    
    !Plot grid for given line
    if(IsFirstCall) then
       call save_plot_file(NamePlot, TypePositionIn='rewind', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
            nDimIn=nDim,CoordIn_DII=Coord_DII,                &
            VarIn_IIV = PlotFlux_IIV, ParamIn_I = (/1.6, 1.0/)) !***
       IsFirstCall = .false.
    else
       call save_plot_file(NamePlot, TypePositionIn='append', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
            nDimIn=nDim,CoordIn_DII=Coord_DII,                &
            VarIn_IIV = PlotFlux_IIV, ParamIn_I = (/1.6, 1.0/)) !***
    endif
    
    deallocate(Coord_DII, PlotFlux_IIV)
  end subroutine plot_flux

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  ! save crossection plot for verification
  ! doesn't work for Earth yet
  ! for debugging only one line
  !--------------------------------------------------------------------------
  subroutine plot_crossec(SIGION,PROB)
    use ModSeGrid,     ONLY: Alt_C,nAlt,nEnergy,  &
         DeltaE_I,EnergyGrid_I
    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModNumConst,   ONLY: cRadToDeg,cPi
    use ModPlanetConst, only: Planet_, NamePlanet_I

    real, intent(in) :: SIGION(nNeutral,LMAX),PROB(nStatesMax,nNeutral,LMAX)

    real, parameter :: C1=12397.7       ! Converting wavelengths to energy

    real, allocatable   :: Coord_I(:), PlotState_IV(:,:)
    !grid parameters
    integer, parameter :: nDim =1
    integer :: nVar
    integer, parameter ::iLine=1
    !Jupiter
    integer, parameter :: H2plus_=1,Heplus_=2,Hplus_=3,CH4plus_=4,&
         CH3plus_=5,CH2plus_=6,CHplus_=7
    
    !Earth
    integer, parameter :: Oplus_=1
    
    character(len=100),parameter :: NamePlotVarEarth=&
         'E[eV] Var g r'
    character(len=100),parameter :: NamePlotVarJupiter=&
         'E[eV] SigmaIon Prob g r'

    character(len=100) :: NamePlotVar
    character(len=*),parameter :: NameHeader='He+ Cross Section'
    character(len=5) :: TypePlot='ascii'
    integer :: iIon,iAlt
    character(len=100) :: NamePlot
    logical,save :: IsFirstCall =.true.
    !--------------------------------------------------------------------------

!    nVar=nIons!+2
    nVar=2
    allocate(Coord_I(LMAX),PlotState_IV(LMAX,nVar))

    PlotState_IV = 0.0
    Coord_I     = 0.0
    
    select case(NamePlanet_I(Planet_))
    case('EARTH')
       NamePlotVar=NamePlotVarEarth
    case('JUPITER')
       NamePlotVar=NamePlotVarJupiter
    end select
       
       !Set Flux Coordinates along field line and wavelength converted to eV
       do iLambda=1,LMAX
             Coord_I(iLambda) = C1/Wave1(iLambda)
             PlotState_IV(iLambda,1) = SIGION(He_,iLambda)
             PlotState_IV(iLambda,2) = PROB(1,He_,iLambda)
       enddo       
       
    ! set name for plotfile
    write(NamePlot,"(a,i4.4,a)") 'CrossSec_iLine',iLine,'.out'
    
    !Plot grid for given line
    if(IsFirstCall) then
       call save_plot_file(NamePlot, TypePositionIn='rewind', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
            nDimIn=nDim,CoordIn_I=Coord_I,                &
            VarIn_IV = PlotState_IV, ParamIn_I = (/1.6, 1.0/)) !***
       IsFirstCall = .false.
    else
       call save_plot_file(NamePlot, TypePositionIn='append', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
            nDimIn=nDim,CoordIn_I=Coord_I,                &
            VarIn_IV = PlotState_IV, ParamIn_I = (/1.6, 1.0/)) !***
    endif
    
    deallocate(Coord_I, PlotState_IV)
  end subroutine plot_crossec


  !
  !
  !
  !
  
!-----------------------------------------------------------------------
  SUBROUTINE BOXNUM (E1, E2, M1, M2, R1, R2, DEL, ENER)
    !      use cglow,only: nbins
    use ModSeGrid,only:NBINS=>nEnergy
    implicit none
    
    ! This subroutine finds the box numbers corresponding to
    ! energies E1 and E2, and calls them M1 and M2
    !
    ! R1 is the upper edge of the lower box, R2 is the lower edge of the
    ! upper box.
    ! Args:
    Real,Intent(in)    :: DEL(NBINS), ENER(NBINS)
    Real,Intent(in)    :: E1,E2
    Real,Intent(out)   :: R1,R2
    Integer,Intent(out):: M1,M2
    !Local:
    Integer I
    
    IF (E1 > ENER(NBINS)+DEL(NBINS)/2.) then
       M1 = NBINS+1
       RETURN
    endif
    
    DO I=1,NBINS
       IF (E1 .LT. ENER(I)+DEL(I)/2.) then
          M1 = I
          R1 = ENER(I) + DEL(I)/2.
          exit
       endif
   enddo
!

   !
   IF (E2 > ENER(NBINS)+DEL(I)/2.) then
      M2 = NBINS
      R2 = E2 - DEL(NBINS)
      RETURN
   endif

   DO I=1,NBINS
      IF (E2 .LT. ENER(I)+DEL(I)/2.) then
         M2 = I
         R2 = ENER(I) - DEL(I)/2.
         exit
      endif
   end DO

 END Subroutine BOXNUM

 
  !=============================================================================
  ! Subroutine SOLZEN
  !
  ! Stan Solomon, 1988
  !
  ! Returns Solar Zenith Angle SZA in degrees for specified date in form
  ! yyddd, universal time in seconds, geographic latitude and longitude
  ! in degrees.
  !
  ! note for jupiter or another planet will need planet specific calculation
  SUBROUTINE SOLZEN (IDATE, UT, GLAT, GLONG, SZA)
    !
    DATA PI/3.1415926536/
    !
    RLAT = GLAT * PI/180.
    RLONG = GLONG * PI/180.

    ! Convert IDATE to 4-digit year and 3-digit day of year
    ! for dates between 1950 and 2049
    IYR=IDATE/1000
    IDAY=IDATE-IYR*1000
    if (IYR.LT.50) IYR = IYR+100
    IYR = IYR+1900
    
    CALL SunCoordsGEI (IYR, IDAY, UT, SDEC, SRASN, GST)
    write(*,*) IDATE, SDEC,SRASN,GST
    RH = SRASN - (GST+RLONG)
    COSSZA = SIN(SDEC)*SIN(RLAT) + COS(SDEC)*COS(RLAT)*COS(RH)
    SZA = ACOS(COSSZA) * 180./PI
    RETURN
  END SUBROUTINE SOLZEN
  !
  !
  !
  ! Subroutine SunCoordsGEI returns the declination SDEC and right ascension
  ! SRASN of the sun in GEI coordinates, radians, for a given date IDATE
  ! in yyddd format and universal time UT in seconds.  Greenwich Sidereal
  ! Time GST in radians is also returned. For years 1901-2099. Accuracy 0.006
  ! degree. Reference:  C.T. Russell, Geophysical Coordinate Transforms.
  !
  Subroutine SunCoordsGEI (iYear, iDOY, UTseconds, Dec, RA, GST)
    use ModNumConst, ONLY: cPi,cDegToRad
    integer, intent(in) :: iYear, iDOY
    real, intent(in)    :: UTseconds

    real, intent(out)   :: Dec, RA, GST

    real :: SecPerDay,DaysPerYr,DaysPer100Yr,FractionOfDay,JD1900
    real :: CenturiesElapsed,MeanLongitude,MeanAnomaly,EclipticLongitude
    real :: ObliquityOfEcliptic,SinDec,CosDec
    !
    SecPerDay = 86400.
    DaysPerYr = 365.
    DaysPer100Yr=36525.
    FractionOfDay=UTseconds/SecPerDay
    JD1900=DaysPerYr*(iYear-1900)+(iYear-1901)/4+iDOY+FractionOfDay-0.5
    CenturiesElapsed=JD1900/DaysPer100Yr
    MeanLongitude=AMOD(279.696678+.9856473354*JD1900,360.)
    GST=AMOD(279.696678+.9856473354*JD1900+360.*FractionOfDay+180.,360.) &
         * cDegToRad
    MeanAnomaly=AMOD(358.475845+.985600267*JD1900,360.) * cDegtoRad
    EclipticLongitude = MeanLongitude + &
         (1.91946-.004789*CenturiesElapsed)*SIN(MeanAnomaly) + &
         .020094*SIN(2.*MeanAnomaly)
    ObliquityOfEcliptic=(23.45229-0.0130125*CenturiesElapsed) * cDegtoRad
    EclipticLongitude=(EclipticLongitude-.005686) * cDegtoRad
    SinDec=SIN(ObliquityOfEcliptic)*SIN(EclipticLongitude)
    CosDec=SQRT(1.-SinDec**2)
    Dec=ATAN(SinDec/CosDec)
    RA=cPi-ATAN2(1./TAN(ObliquityOfEcliptic)*SinDec/CosDec, &
         -COS(EclipticLongitude)/CosDec)
    RETURN
  END SUBROUTINE SunCoordsGEI


  ! Subroutine SSFLUX
  !
  ! This software is part of the GLOW model.  Use is governed by the Open Source
  ! Academic Research License Agreement contained in the file glowlicense.txt.
  ! For more information see the file glow.txt.
  !
  ! Subroutine SSFLUX calculates the solar EUV and FUV flux in the range
  ! 0.5 to 1750 Angstroms for a specified level of solar activity.
  !
  ! The calling routine supplies a scaling switch ISCALE, the daily 10.7
  ! cm flux F107, its 81-day centered average F107A, and, optionally,
  ! the H Lyman-Beta 1026A ratio to its solar minimum value HLYBR, the
  ! Fe XVI 335A ratio to its solar minimum value FEXVIR, the H Lyman-alpha
  ! 1216A flux HLYA, the He I 10830A equivalent width HEIEW, and an XUV
  ! enhancement factor XUVFAC.  Any optional calling parameters not used
  ! should be set to zero.
  !
  ! XUVFAC is applied from 18-250 A for the Hinteregger model (ISCALE=0),
  ! from 18-50 A for the EUVAC model (ISCALE=1), and not at all for
  ! user-supplied data (ISCALE=2)
  !
  ! The subroutine returns the longwave boundary WAVE1 and shortwave
  ! boundary WAVE2 of the wavelenth bins, and the solar flux in each bin
  ! SFLUX.  Bins are arranged in energy multiples of 2 from 0.5 to 8 A,
  ! aligned with k-shell boundaries at 23, 32, and 44 A from 18 to 44 nm,
  ! 10 A in width from 60 to 1050 A, and 50 A in width from 1050 to
  ! 1750 A with the exception of Lyman-alpha which has its own bin from
  ! 1210 to 1220 A.
  !
  ! Methods used:
  !   If ISCALE=0 the flux is scaled using parameterization methods based
  ! on F107 and F107A.  For ionizing EUV, Hinteregger's contrast ratio
  ! method (Hinteregger et al., GRL, 8, 1147, 1981) is used, based on the
  ! reference spectrum SC#21REFW at 1 nm bin resolution.  If the
  ! H Lyman-Beta (1026A) or Fe XVI (335A) enhancement ratios are provided
  ! (>0) as calling arguments, they are used to scale the EUV spectrum.
  ! Otherwise, enhancement ratios for H Ly-B and Fe XVI are calculated
  ! from F107 and F107A using Hinteregger's formula, employing
  ! coefficients which reduce to the reference values at F107=67.6,
  ! F107A=71.5.  The 'best fit' coefficients are not used as they produce
  ! some negative values at low solar activity, but remain in a
  ! 'commented out' data statement for reference.  The EUV spectrum is
  ! then scaled from these modeled ratios.  Scaling factors were
  ! calculated from contrast ratios in the SC#21REFW data file.
  !   If ISCALE=1, the EUV flux (50-1050A) is scaled using the EUVAC model
  ! (Richards et al., JGR 99, 8981, 1994) re-binned onto ~1 nm intervals.
  ! The Hinteregger spectrum, scaled using the EUVAC algorithm, is used
  ! from 18 to 50A.
  !   Neither of these models extends shortward of 18A, so from 1-18 A
  ! an amalgam of sources are used to derive an estimated flux, e.g.,
  ! DeJager, in Astronomical Observations from Space Vehicles, Steinberg,
  ! ed., 1964; Smith & Gottlieb, SSR 16, 771, 1974; Manson, in The Solar
  ! Output and its Variation, White, ed., 1977; Kreplin et al, ibid;
  ! Horan & Kreplin, Solar Physics 74, 265, 1981; Wagner, Adv. Space Res.
  ! 8, (7)67, 1988.
  !    For FUV from 1050A-1750A, 50A interval bins from the Woods and
  ! Rottman [2002] reference spectrum and scale factors based on
  ! UARS SOLSTICE data are used.  The scaling method follows the
  ! Hinteregger or EUVAC algorithm, whichever is selected, so as to
  ! linearly scale the spectrum between the reference value and maximum
  ! value calculated with F10.7=F10.7A=200.  If a value for Lyman-alpha
  ! (HLYA>0) is provided by the calling program, it is subsituted into
  ! the spectrum.
  !   If ISCALE=2, the solar flux (0-1750A) is read from a file named
  ! ssflux_user.dat in the current working directory.  The file must
  ! contain three columns:  WAVES, WAVEL, SFLUX (Angstroms and cm-2 s-1)
  ! in order of increasing wavelength.  The number of lines in the file
  ! must match the value of LMAX in glow.h.
  !
  ! Modification history:
  !   Stan Solomon, 12/88  Basic Hinteregger EUV, approx. SME FUV
  !   Chris Gaskill, 7/89  Added early Tobiska model
  !   Stan Solomon,  8/89  Corrections to above
  !   Stan Solomon,  1/90  Tobiska SERF2; added W & R spectra
  !   Stan Solomon,  6/91  Tobiska EUV 91; Hntggr Ly-B, Fe XVI scaling
  !   Stan Solomon,  2/92  Updated Tobiska EUV91; corrected SME FUV
  !   Scott Bailey, 12/93  Initial one-nm bins version
  !   Stan Solomon,  6/04  Added EUVAC option, cleaned up artifacts
  !   Stan Solomon,  9/04  Added ability to specify input data file
  !   Stan Solomon,  3/05  Changed all to photon units
  !   Stan Solomon,  1/15  Updated for f90; only read file on first call
  !
  ! Calling parameters:
  ! ISCALE   =0 for Hinteregger contrast ratio method
  !          =1 for EUVAC
  !          =2 for user-supplied data
  ! F107     daily 10.7 cm flux (1.E-22 W m-2 Hz-1)
  ! F107A    81-day centered average 10.7 cm flux
  ! HLYBR    ratio of H Ly-b 1026A flux to solar minimum value (optional)
  ! FEXVIR   ratio of Fe XVI 335A flux to solar minimum value (optional)
  ! HLYA     H Lyman-alpha flux (photons cm-2 s-1) (optional)
  ! HEIEW    He I 10830A equivalent width (mAngstroms) (obsolete)
  ! XUVFAC   factor for scaling flux 18-250A or 18-50A (optional)
  
  ! Returned parameters:
  ! WAVE1    longwave bound of spectral intervals (Angstroms)
  ! WAVE2    shortwave bound of intervals
  ! SFLUX    scaled solar flux returned by subroutine (photons cm-2 s-1)
  !
  ! Other definitions:
  ! LMAX     dimension of flux and scaling arrays, currently = 123
  ! WAVEL    = WAVE1
  ! WAVES    = WAVE2
  ! RFLUX    low solar activity flux
  ! XFLUX    high solar activity flux
  ! SCALE1   scaling factors for H LyB-keyed chromospheric emissions
  ! SCALE2   scaling factors for FeXVI-keyed coronal emissions
  ! B1       fit coefficients for H LyB
  ! B2       fit coefficients for FeXVI
  ! R1       enhancement ratio for H LyB
  ! R2       enhancement ratio for FeXVI
  ! P107     average of F107 and F107A
  ! A        scaling factor for EUVAC model
  !
  !
  SUBROUTINE SSFLUX (ISCALE, F107, F107A, HLYBR, FEXVIR, HLYA, &
       HEIEW, XUVFAC)
    use ModPlanetConst, ONLY: Planet_, NamePlanet_I
    use ModIoUnit,     ONLY: UnitTmp_
    !      use cglow,only: lmax
    implicit none
    !      include 'cglow.h'
    
    ! Args:
    integer,intent(in) :: iscale
    real,intent(in)    :: f107, f107a,HLYBR, FEXVIR,HLYA,XUVFAC,heiew

    ! Local:
    real WAVEL(LMAX), WAVES(LMAX), RFLUX(LMAX), UFLUX(LMAX), &
         SCALE1(LMAX), SCALE2(LMAX), A(LMAX),p107,r1,r2
    integer l
    real,parameter :: epsil=1.0E-6
    !      integer :: islast=-1  !this didn't work right in f2py
    
    
    ! regression coefficients which reduce to solar min. spectrum:
    real :: B1(3),B2(3)
    DATA B1/1.0, 0.0138, 0.005/
    DATA B2/1.0, 0.59425, 0.3811/
    
    ! 'best fit' regression coefficients, commented out, for reference:
    !     DATA B1/1.31, 0.01106, 0.00492/, B2/-6.618, 0.66159, 0.38319/
    
    
    ! Hinteregger contrast ratio method:
    
    IF (iscale .eq. 0) then
       !        if (islast .ne. iscale) then
       open(unit=UnitTmp_,file='PW/ssflux_hint.dat',status='old')
       read(UnitTmp_,*)
       do l=lmax,1,-1
          read(UnitTmp_,*) waves(l),wavel(l),rflux(l),scale1(l),scale2(l)
       enddo
       close(unit=UnitTmp_)
       !         endif
       
       IF (HLYBR .GT. EPSIL) THEN
          R1 = HLYBR
       ELSE
          R1 =  B1(1) + B1(2)*(F107A-71.5) + B1(3)*(F107-F107A+3.9)
       ENDIF
       IF (FEXVIR .GT. EPSIL) THEN
          R2 = FEXVIR
       ELSE
          R2 =  B2(1) + B2(2)*(F107A-71.5) + B2(3)*(F107-F107A+3.9)
       ENDIF
       
       do l=1,lmax
          SFLUX(L) = RFLUX(L) + (R1-1.)*SCALE1(L) + (R2-1.)*SCALE2(L)
          IF (SFLUX(L) .LT. 0.0) SFLUX(L) = 0.0
          IF (XUVFAC .GT. EPSIL .AND. &
               WAVEL(L).LT.251.0 .AND. WAVES(L).GT.17.0) &
               SFLUX(L)=SFLUX(L)*XUVFAC
          !write(*,*) 'L,WAVEL(L),SFLUX(L)',L,WAVEL(L),SFLUX(L)
       enddo
    endif
    
    ! EUVAC Method:
    
    IF (iscale .eq. 1) then
       !        if (islast .ne. iscale) then
       open(unit=UnitTmp_,file='PW/ssflux_euvac.dat',status='old')
       read(UnitTmp_,*)
       do l=lmax,1,-1
          read(UnitTmp_,*) waves(l),wavel(l),rflux(l),a(l)
       enddo
       close(unit=UnitTmp_)
       !        endif
       
       P107 = (F107+F107A)/2.
       do l=1,lmax
          SFLUX(L) = RFLUX(L) * (1. + A(L)*(P107-80.))
          IF (SFLUX(L) .LT. 0.8*RFLUX(L)) SFLUX(L) = 0.8*RFLUX(L)
          IF (XUVFAC .GT. EPSIL .AND. &
               WAVEL(L).LT.51.0 .AND. WAVES(L).GT.17.0) &
               SFLUX(L)=SFLUX(L)*XUVFAC
       enddo
    ENDIF
    
    ! User-supplied data:
    
    if (iscale .eq. 2) then
       !        if (islast .ne. iscale) then
       open(unit=UnitTmp_,file='PW/ssflux_user.dat',status='old')
       read(UnitTmp_,*)
       do l=lmax,1,-1
          read(UnitTmp_,*) waves(l),wavel(l),uflux(l)
       enddo
       close(unit=UnitTmp_)
       !        endif
       do l=1,lmax
          sflux(l)=uflux(l)
       enddo
    endif
    
    ! Fill wavelength arrays, substitute in H Lyman-alpha if provided:
    
    do l=1,lmax
       WAVE1(L) = WAVEL(L)
       WAVE2(L) = WAVES(L)
       IF (HLYA .GT. EPSIL .AND. &
            WAVEL(L).LT.1221. .AND. WAVES(L).GT.1209.) &
            SFLUX(L) = HLYA
    enddo

    ! Fluxes assumed for Earth, but should be scaled for other planets
    select case(NamePlanet_I(Planet_))
    case('JUPITER')
       SFLUX(:) = SFLUX(:) * (1.0/5.2)**2
    case('SATURN')
       SFLUX(:) = SFLUX(:) * (1.0/9.5)**2
    end select
    
    !      islast=iscale
    
  End SUBROUTINE SSFLUX




  !============================================================================
  subroutine init_production
    use ModPlanetConst, ONLY: Planet_, NamePlanet_I, rPlanet_I
    use ModIoUnit,     ONLY: UnitTmp_
    !\
    ! Earth
    !/
    integer, parameter :: nStatesMaxEarth=6 !only for setting data arrays
    ! branching ratios (O,O2,N2)
    real :: PROBO(nStatesMaxEarth,LMAX), PROBO2(nStatesMaxEarth,LMAX) 
    real :: PROBN2(nStatesMaxEarth,LMAX)
    ! absorption crossections (O,O2,N2)
    real :: SIGAO(LMAX), SIGAO2(LMAX), SIGAN2(LMAX)
    ! ionization crossections (O,O2,N2)
    real :: SIGIO(LMAX), SIGIO2(LMAX), SIGIN2(LMAX)
    
    !\
    ! Jupiter/Saturn
    !/
    ! branching ratios
    real,allocatable :: ProbSpecies(:,:)
    ! ionization and absorption crossections
    real :: SigIonSpecies(LMAX),SigAbsSpecies(LMAX)


    select case(NamePlanet_I(Planet_))
    case('EARTH')
       ! Earth
       ! set the surface gravity
       gSurface = 978.1

       ! set the neutral parameters
       nNeutral=3

       ! now that nNeutral is set allocate arrays       
       call allocate_neutral_arrays

       ! set the neutral mass per species
       NeutralMassAMU_I = (/16.0,32.0,28.0/)

       ! set auger e- mean energies and wavelength thresholds
       AugEnergy_I = (/500., 500., 360./)
       
       ! set wave length thresholds for Auger production
       AugWave_I = (/24.,24.,33./)

       ! set number of states per species
       nStatesPerSpecies_I=(/5,4,6/)

       nStatesMax = maxval(nStatesPerSpecies_I)

       call allocate_state_arrays(ProbSpecies)

       ! set ionization potentials per state, per species
       TPOT(:,O_) = (/13.61, 16.93, 18.63, 28.50, 40.00,  0.00/)
       TPOT(:,O2_)= (/12.07, 16.10, 18.20, 20.00,  0.00,  0.00/) 
       TPOT(:,N2_)= (/15.60, 16.70, 18.80, 30.00, 34.80, 25.00/)
       
       !read ionization and absorption states
       open(unit=UnitTmp_,file='PW/ephoto_xn2.dat',status='old')
       read(UnitTmp_,*)
       read(UnitTmp_,*)
       read(UnitTmp_,*)
       read(UnitTmp_,*)
       do l=lmax,1,-1
          read(UnitTmp_,*) aa,bb,(probn2(n,l),n=1,nStatesMax),sigin2(l),sigan2(l)
          !bso2(l)=0.0
       enddo
      close(UnitTmp_)
!
      open(unit=UnitTmp_,file='PW/ephoto_xo2.dat',status='old')
      read(UnitTmp_,*)
      read(UnitTmp_,*)
      read(UnitTmp_,*)
      read(UnitTmp_,*)
      do l=lmax,1,-1
         read(UnitTmp_,*) aa,bb,(probo2(n,l),n=1,nStatesMax),sigio2(l),sigao2(l)
      enddo
      close(UnitTmp_)
!
      open(unit=UnitTmp_,file='PW/ephoto_xo.dat',status='old')
      read(UnitTmp_,*)
      read(UnitTmp_,*)
      read(UnitTmp_,*)
      read(UnitTmp_,*)
      do l=lmax,1,-1
         read(UnitTmp_,*) aa,bb,(probo(n,l),n=1,nStatesMax),sigio(l),sigao(l)
      enddo
      close(UnitTmp_)

       !set branching ratios and crossections
       DO  L=1,LMAX
          DO  K=1,nStatesMax
             PROB(K,1,L) = PROBO(K,L)
             PROB(K,2,L) = PROBO2(K,L)
             PROB(K,3,L) = PROBN2(K,L)
          end do
       end DO
       
       DO  L=1,LMAX
          SIGABS(1,L) = SIGAO(L)  * 1.E-18
          SIGABS(2,L) = SIGAO2(L) * 1.E-18
          SIGABS(3,L) = SIGAN2(L) * 1.E-18
          SIGION(1,L) = SIGIO(L)  * 1.E-18
          SIGION(2,L) = SIGIO2(L) * 1.E-18
          SIGION(3,L) = SIGIN2(L) * 1.E-18
       end DO
       
       
    case('JUPITER')
       ! Jupiter
       ! set the surface gravity
       gSurface = 2479.0

       ! set the neutral parameters
       nNeutral=4

       ! now that nNeutral is set allocate arrays
       call allocate_neutral_arrays

       ! set the neutral mass per species
       NeutralMassAMU_I = (/2.0,4.0,1.0,16.0/)

       ! set auger e- mean energies and wavelength thresholds
       ! both set to zero for no Auger production
       AugEnergy_I = (/0.0, 0.0, 0.0, 0.0/)
       
       ! set iAugWaveBin_I to zero since we have no Auger production for now
       AugWave_I = (/0., 0., 0., 0./)

       ! set number of states per species
!       nStatesPerSpecies_I=(/1,1,1,3/)
       nStatesPerSpecies_I=(/2,1,1,5/)

       nStatesMax = maxval(nStatesPerSpecies_I)

       call allocate_state_arrays(ProbSpecies)

       ! set ionization potentials per state, per species
       !H2+      H+ H  
       TPOT(:,H2_) = (/15.427, 18.08, 0.00, 0.00, 0.00/) 
       !He+
       TPOT(:,He_) = (/24.6, 0.00, 0.00, 0.00, 0.00/) 
       !H
       TPOT(:,H_)  = (/13.5, 0.00, 0.00, 0.00, 0.00/) 
       !CH4+     CH3+H    CH2+H2   CH+H2/H  H+CH3
       TPOT(:,CH4_)= (/13.12, 14.23, 14.99, 19.70, 18.08/) 

!       TPOT(:,H2_) = (/15.42589, 0.00, 0.00/) 
!       TPOT(:,He_) = (/24.6, 0.00, 0.00/) 
!       TPOT(:,H_)  = (/13.5, 0.00, 0.00/) 
!       TPOT(:,CH4_)= (/12.98, 24.0, 27.55/)

!       call read_data_array('PW/PhotoH2.dat',nStatesPerSpecies_I(1),2, &
!            ProbSpecies,SigAbsSpecies,SigIonSpecies)
!       PROB(:,1,:) = ProbSpecies(:,:)
!       SIGABS(1,:) = SigAbsSpecies(:) * 1.e-18
!       SIGION(1,:) = SigIonSpecies(:) * 1.e-18

       call read_phidrates('H2',nStatesPerSpecies_I(1), &
            ProbSpecies,SigAbsSpecies,SigIonSpecies)
       PROB(:,1,:) = ProbSpecies(:,:)
       SIGABS(1,:) = SigAbsSpecies(:)
       SIGION(1,:) = SigIonSpecies(:)

       call read_data_array('PW/PhotoHe.dat',nStatesPerSpecies_I(2),1, &
            ProbSpecies,SigAbsSpecies,SigIonSpecies)
       PROB(:,2,:) = ProbSpecies(:,:)
       SIGABS(2,:) = SigAbsSpecies(:) * 1.e-18
       SIGION(2,:) = SigIonSpecies(:) * 1.e-18
       call read_data_array('PW/PhotoH.dat',nStatesPerSpecies_I(3),1, &
            ProbSpecies,SigAbsSpecies,SigIonSpecies)
       PROB(:,3,:) = ProbSpecies(:,:)
       SIGABS(3,:) = SigAbsSpecies(:) * 1.e-18
       SIGION(3,:) = SigIonSpecies(:) * 1.e-18
!       call read_data_array('PW/PhotoCH4.dat',nStatesPerSpecies_I(4),2, &
!            ProbSpecies,SigAbsSpecies,SigIonSpecies)
!       PROB(:,4,:) = ProbSpecies(:,:)
!       SIGABS(4,:) = SigAbsSpecies(:) * 1.e-18
!       SIGION(4,:) = SigIonSpecies(:) * 1.e-18

       call read_phidrates('CH4',nStatesPerSpecies_I(4), &
            ProbSpecies,SigAbsSpecies,SigIonSpecies)
       PROB(:,4,:) = ProbSpecies(:,:)
       SIGABS(4,:) = SigAbsSpecies(:)
       SIGION(4,:) = SigIonSpecies(:)


       !set branching ratios and crossections
!       DO  L=1,LMAX
!          DO  K=1,nStatesMax
!             PROB(K,1,L) = PROBH2(K,L)
!             PROB(K,2,L) = PROBH(K,L)
!             PROB(K,3,L) = PROBCH4(K,L)
!             PROB(K,4,L) = PROBHe(K,L)
!          end do
!       end DO
       
!       DO  L=1,LMAX
!          SIGABS(1,L) = SIGAH2(L)  * 1.E-18
!          SIGABS(2,L) = SIGAH(L) * 1.E-18
!          SIGABS(3,L) = SIGACH4(L) * 1.E-18
!          SIGABS(4,L) = SIGAHe(L) * 1.E-18
!          SIGION(1,L) = SIGIH2(L)  * 1.E-18
!          SIGION(2,L) = SIGIH(L) * 1.E-18
!          SIGION(3,L) = SIGICH4(L) * 1.E-18
!          SIGION(4,L) = SIGIHe(L) * 1.E-18
!       end DO

    end select


  end subroutine init_production


  subroutine allocate_neutral_arrays
    ! array to hold neutral mass
    if (.not.allocated(NeutralMassAMU_I)) &
         allocate(NeutralMassAMU_I(nNeutral))
    ! states per species
    if (.not.allocated(nStatesPerSpecies_I)) &
         allocate(nStatesPerSpecies_I(nNeutral))
    ! auger e- mean energy per species
    if (.not.allocated(AugEnergy_I)) &
         allocate(AugEnergy_I(nNeutral))
    ! auger e- wavelength threshold per species
    if (.not.allocated(AugWave_I)) &
         allocate(AugWave_I(nNeutral))
    ! ionization and absorption crossections
    if (.not.allocated(SigIon)) &
         allocate(SigIon(nNeutral,LMAX))
    if (.not.allocated(SigAbs)) &
         allocate(SigAbs(nNeutral,LMAX))
    
  end subroutine allocate_neutral_arrays
!=======

  subroutine allocate_state_arrays(ProbSpecies)

    real, allocatable :: ProbSpecies(:,:)

    ! allocate ionization potentials per state, per species
    if (.not.allocated(TPOT)) &
         allocate(TPOT(nStatesMax,nNeutral))
    ! allocate generic species branching probability array
    if (.not.allocated(ProbSpecies)) &
         allocate(ProbSpecies(nStatesMax,LMAX))
    ! allocate global branching probability array
    if (.not.allocated(Prob)) &
         allocate(Prob(nStatesMax,nNeutral,LMAX))

  end subroutine allocate_state_arrays
  !=============================================================================
  subroutine read_data_array(DatafileName,nStates,GridType, &
       ProbSpecies,SigAbsSpecies,SigIonSpecies)
    use ModInterpolate, ONLY: linear
    use ModIoUnit,      ONLY: UnitTmp_
    character (len=*), intent(in) :: DatafileName
    integer, intent(in) :: nStates
    integer, intent(in) :: GridType !  (1 for points, 2 for bins)

    real, intent(out) :: ProbSpecies(nStatesMax,LMAX), &
         SigAbsSpecies(LMAX),SigIonSpecies(LMAX)

    character (len=8) :: inputfmt ! e.g. (8E10.6)
    integer :: nLines
    real :: StandardWaveGrid(LMAX)
    real, allocatable :: SigAbsIn(:),SigIonIn(:),ProbIn(:,:),WaveGrid(:,:)
    real, allocatable :: WaveGridCenters(:)

    open(UnitTmp_,FILE=DatafileName,STATUS='OLD')

    read(UnitTmp_,*) nLines
    
    allocate(SigAbsIn(nLines))
    allocate(SigIonIn(nLines))
    allocate(ProbIn(nStates,nLines))
    allocate(WaveGrid(GridType,nLines))
    allocate(WaveGridCenters(nLines))

    print *,'GridType: ',GridType,'nLines: ',nLines

    do i=1,nLines
       write(inputfmt,'("(", I0, "E10.6)")') 2+nStates+GridType
       read(UnitTmp_,*) WaveGrid(:,i),SigAbsIn(i),ProbIn(:,i), &
            SigIonIn(i)
!5001   format(<2+nStates+GridType>E10.6)
    enddo
    close(UnitTmp_)

    print *,WaveGrid(:,nLines),SigAbsIn(nLines),ProbIn(:,nLines), &
         SigIonIn(nLines)

    if (GridType == 2) then
       WaveGridCenters(:) = (WaveGrid(1,:) + WaveGrid(2,:))/2.
    else
       WaveGridCenters(:) = WaveGrid(1,:)
    endif

    StandardWaveGrid = (WAVE1 + WAVE2)/2.
    
    ! interpolate all variables to new wavelength grid
    ! extrapolate set to true
    do i=1,LMAX
       !when outside of range set to zero otherwise interpolate
       if (StandardWaveGrid(i) > maxval(WaveGridCenters) &
            .or. StandardWaveGrid(i) < minval(WaveGridCenters)) then
          SigAbsSpecies(i) = 0.0
          SigIonSpecies(i) = 0.0
          ProbSpecies(:,i) = 0.0
       else
          SigAbsSpecies(i) = linear(SigAbsIn,1,nLines, &
               StandardWaveGrid(i),WaveGridCenters)
          SigIonSpecies(i) = linear(SigIonIn,1,nLines, &
               StandardWaveGrid(i),WaveGridCenters)
          do l=1,nStates
             ProbSpecies(l,i) = linear(ProbIn(l,:),1,nLines, &
                  StandardWaveGrid(i),WaveGridCenters)/100.0
          end do
       endif
    enddo

    write(*,*) 'For DatafileName=',DataFileName
    write(*,*) 'maxval(SigIonSpecies),minval(SigIonSpecies)',&
         maxval(SigIonSpecies),minval(SigIonSpecies)
    
    deallocate(SigAbsIn)
    deallocate(SigIonIn)
    deallocate(ProbIn)
    deallocate(WaveGrid)
    deallocate(WaveGridCenters)

  end subroutine read_data_array

  !=============================================================================
  subroutine read_phidrates(NameNeutral,nStates, &
       ProbSpecies,SigAbsSpecies,SigIonSpecies)
    use ModInterpolate, ONLY: linear
    use ModIoUnit,      ONLY: UnitTmp_
    character (len=*), intent(in) :: NameNeutral
    integer, intent(in) :: nStates

    real, intent(out) :: ProbSpecies(nStatesMax,LMAX), &
         SigAbsSpecies(LMAX),SigIonSpecies(LMAX)

    integer :: nLines
    real :: StandardWaveGrid(LMAX)
    real, allocatable :: SigAbsIn(:),SigIonIn(:),ProbIn(:,:),WaveGrid(:,:)
    real, allocatable :: WaveGridCenters(:)
    character(len=100) :: junk
    real :: junk1,junk2,junk3,junk4,cross1,cross2,cross3,cross4,cross5
    
    if(NameNeutral == 'H2') then
       open(UnitTmp_,FILE='PW/phidratesH2.dat',STATUS='OLD')
       !read and discard header
       do iLine=1,2
          read(UnitTmp_,*) junk
       enddo
       nLines=102
    elseif(NameNeutral == 'CH4') then
       open(UnitTmp_,FILE='PW/phidratesCH4.dat',STATUS='OLD')
        do iLine=1,2
          read(UnitTmp_,*) junk
       enddo
       nLines=137
    else
       call con_stop(NameNeutral//' not available')
    endif

    allocate(SigAbsIn(nLines))
    allocate(SigIonIn(nLines))
    allocate(ProbIn(nStates,nLines))
    allocate(WaveGridCenters(nLines))

    do i=1,nLines
       if(NameNeutral == 'H2') then
          read(UnitTmp_,*) WaveGridCenters(i),SigAbsIn(i),junk1,junk2,&
               cross1,cross2
          SigIonIn(i) = cross1+cross2
          if (SigIonIn(i) == 0.0 ) then
             ProbIn(1,i) = 0.0
             ProbIn(2,i) = 0.0
          else
             ProbIn(1,i) = cross1/SigIonIn(i)
             ProbIn(2,i) = cross2/SigIonIn(i)
          endif
       elseif(NameNeutral == 'CH4') then
          read(UnitTmp_,*) WaveGridCenters(i),SigAbsIn(i),junk1,junk2,junk3,&
               cross1,cross2,cross3,cross4,cross5,junk4
          SigIonIn(i) = cross1+cross2+cross3+cross4+cross5
          if (SigIonIn(i) == 0.0 ) then
             ProbIn(1,i) = 0.0
             ProbIn(2,i) = 0.0
             ProbIn(3,i) = 0.0
             ProbIn(4,i) = 0.0
             ProbIn(5,i) = 0.0
          else
             ProbIn(1,i) = cross1/SigIonIn(i)
             ProbIn(2,i) = cross2/SigIonIn(i)
             ProbIn(3,i) = cross3/SigIonIn(i)
             ProbIn(4,i) = cross4/SigIonIn(i)
             ProbIn(5,i) = cross5/SigIonIn(i)
          endif
          
       endif
              
!5001   format(<2+nStates+GridType>E10.6)
    enddo
    close(UnitTmp_)

    StandardWaveGrid = (WAVE1 + WAVE2)/2.
    
    ! interpolate all variables to new wavelength grid
    ! extrapolate set to true
    do i=1,LMAX
       !when outside of range set to zero otherwise interpolate
       if (StandardWaveGrid(i) > maxval(WaveGridCenters) &
            .or. StandardWaveGrid(i) < minval(WaveGridCenters)) then
          SigAbsSpecies(i) = 0.0
          SigIonSpecies(i) = 0.0
          ProbSpecies(:,i) = 0.0
       else
          SigAbsSpecies(i) = linear(SigAbsIn,1,nLines, &
               StandardWaveGrid(i),WaveGridCenters)
          SigIonSpecies(i) = linear(SigIonIn,1,nLines, &
               StandardWaveGrid(i),WaveGridCenters)
          do l=1,nStates
             ProbSpecies(l,i) = linear(ProbIn(l,:),1,nLines, &
                  StandardWaveGrid(i),WaveGridCenters)
          end do
       endif
    enddo

    write(*,*) 'For NameNeutral=',NameNeutral
    write(*,*) 'maxval(SigIonSpecies),minval(SigIonSpecies)',&
         maxval(SigIonSpecies),minval(SigIonSpecies)
    
    deallocate(SigAbsIn)
    deallocate(SigIonIn)
    deallocate(ProbIn)
    deallocate(WaveGridCenters)

  end subroutine read_phidrates

  !============================================================================
  ! subroutine to get the flux contribution in four wavelengths for 
  ! resonant scattering in plasmasphere
  ! Xy_D is input: Xy_D(1)=SZA and Xy_D(2) = Alt of required point
  ! Flux is output
  subroutine get_plas_resonant_scattering(Xy_D,NameLamda,Flux)
    use ModTriangulate,ONLY:calc_triangulation, find_triangle
    use ModIoUnit,     ONLY: UnitTmp_
    implicit none
    real, intent(in)  :: Xy_D(2)
    character(len=*), intent(in) :: NameLamda
    real, intent(out) :: Flux
    
    integer, parameter   :: nPointLyAlpha=56
    real,save,    allocatable :: CoordLyAlphaXy_DI(:,:),FluxLyAlpha_I(:)
    integer,save, allocatable :: iNodeTriangleLyAlpha_II(:,:)
    integer,save :: nTriangleLyAlpha
    integer, parameter   :: nPointLyBeta=43
    real,save,    allocatable :: CoordLyBetaXy_DI(:,:),FluxLyBeta_I(:)
    integer,save, allocatable :: iNodeTriangleLyBeta_II(:,:)
    integer,save :: nTriangleLyBeta

    integer, parameter   :: nPointHeI=41
    real,save,    allocatable :: CoordHeIXy_DI(:,:),FluxHeI_I(:)
    integer,save, allocatable :: iNodeTriangleHeI_II(:,:)
    integer,save :: nTriangleHeI

    integer, parameter   :: nPointHeII=46
    real,save,    allocatable :: CoordHeIIXy_DI(:,:),FluxHeII_I(:)
    integer, save,allocatable :: iNodeTriangleHeII_II(:,:)
    integer,save :: nTriangleHeII


    integer, parameter :: nCoord=2
    integer    :: iPoint
    logical    :: IsTriangleFound
    integer    :: iNode1, iNode2, iNode3
    real       :: Area1, Area2, Area3

    logical, save ::IsFirstCall = .true.
    !---------------------------------------------------------------------------
    ! on first call allocate arrays and read all data files and calc the 
    ! triangulations
    
    if (IsFirstCall) then
       !\
       ! LyAlpha
       !/
       ! allocate arrays
       allocate(iNodeTriangleLyAlpha_II(3, 2*nPointLyAlpha),&
            CoordLyAlphaXy_DI(nCoord, nPointLyAlpha),&
            FluxLyAlpha_I(nPointLyAlpha))
       
       ! read data
       open(UNIT=UnitTmp_, FILE='PW/strobel_LyAlpha.dat', STATUS='OLD')
       do iPoint=1,nPointLyAlpha
          read(UnitTmp_,*) CoordLyAlphaXy_DI(1,iPoint),&
               CoordLyAlphaXy_DI(2,iPoint),FluxLyAlpha_I(iPoint)
       end do
       close(UnitTmp_)
       ! find a triangulation on points
       call calc_triangulation(nPointLyAlpha, CoordLyAlphaXy_DI(:,1:nPointLyAlpha), &
            iNodeTriangleLyAlpha_II, nTriangleLyAlpha)

       !\
       ! LyBeta
       !/
       ! allocate arrays
       allocate(iNodeTriangleLyBeta_II(3, 2*nPointLyBeta),&
            CoordLyBetaXy_DI(nCoord, nPointLyBeta),&
            FluxLyBeta_I(nPointLyBeta))
       
       ! read data
       open(UNIT=UnitTmp_, FILE='PW/strobel_LyBeta.dat', STATUS='OLD')
       do iPoint=1,nPointLyBeta
          read(UnitTmp_,*) CoordLyBetaXy_DI(1,iPoint),&
               CoordLyBetaXy_DI(2,iPoint),FluxLyBeta_I(iPoint)
       end do
       close(UnitTmp_)
       ! find a triangulation on points
       call calc_triangulation(nPointLyBeta, CoordLyBetaXy_DI(:,1:nPointLyBeta), &
            iNodeTriangleLyBeta_II, nTriangleLyBeta)

       !\
       ! HeI
       !/
       ! allocate arrays
       allocate(iNodeTriangleHeI_II(3, 2*nPointHeI),&
            CoordHeIXy_DI(nCoord, nPointHeI),&
            FluxHeI_I(nPointHeI))
       
       ! read data
       open(UNIT=UnitTmp_, FILE='PW/strobel_HeI.dat', STATUS='OLD')
       do iPoint=1,nPointHeI
          read(UnitTmp_,*) CoordHeIXy_DI(1,iPoint),&
               CoordHeIXy_DI(2,iPoint),FluxHeI_I(iPoint)
       end do
       close(UnitTmp_)
       ! find a triangulation on points
       call calc_triangulation(nPointHeI, CoordHeIXy_DI(:,1:nPointHeI), &
            iNodeTriangleHeI_II, nTriangleHeI)

       !\
       ! HeII
       !/
       ! allocate arrays
       allocate(iNodeTriangleHeII_II(3, 2*nPointHeII),&
            CoordHeIIXy_DI(nCoord, nPointHeII),&
            FluxHeII_I(nPointHeII))
       
       ! read data
       open(UNIT=UnitTmp_, FILE='PW/strobel_HeII.dat', STATUS='OLD')
       do iPoint=1,nPointHeII
          read(UnitTmp_,*) CoordHeIIXy_DI(1,iPoint),&
               CoordHeIIXy_DI(2,iPoint),FluxHeII_I(iPoint)
       end do
       close(UnitTmp_)
       ! find a triangulation on points
       call calc_triangulation(nPointHeII, CoordHeIIXy_DI(:,1:nPointHeII), &
            iNodeTriangleHeII_II, nTriangleHeII)
       
       IsFirstCall=.false.
    end if

    ! select case based on name of wavelength considered
    select case(NameLamda)
    case('LyAlpha')
       call find_triangle(&
            nPointLyAlpha, nTriangleLyAlpha, Xy_D, &
            CoordLyAlphaXy_DI(:,1:nPointLyAlpha), &
            iNodeTriangleLyAlpha_II(:,1:nTriangleLyAlpha), &
            iNode1, iNode2, iNode3, Area1, Area2, Area3, IsTriangleFound)
       if(IsTriangleFound) then
          Flux = Area1*FluxLyAlpha_I(iNode1) + Area2*FluxLyAlpha_I(iNode2) &
               + Area3*FluxLyAlpha_I(iNode3)
       else
          Flux=0.0
       endif
    case('LyBeta')
       call find_triangle(&
            nPointLyBeta, nTriangleLyBeta, Xy_D, &
            CoordLyBetaXy_DI(:,1:nPointLyBeta), &
            iNodeTriangleLyBeta_II(:,1:nTriangleLyBeta), &
            iNode1, iNode2, iNode3, Area1, Area2, Area3, IsTriangleFound)
       if(IsTriangleFound) then
          Flux = Area1*FluxLyBeta_I(iNode1) + Area2*FluxLyBeta_I(iNode2) &
               + Area3*FluxLyBeta_I(iNode3)
          
!          write(*,*) 'details of found triangle for Xy_D=',Xy_D
!          write(*,*) IsTriangleFound,nTriangleLyBeta
!          write(*,*)CoordLyBetaXy_DI(:,iNode1)
!          write(*,*)CoordLyBetaXy_DI(:,iNode2)
!          write(*,*)CoordLyBetaXy_DI(:,iNode3)
          !stop
       else
          Flux=0.0
!          write(*,*) IsTriangleFound,Xy_D
       endif
       

    case('HeI')
       call find_triangle(&
            nPointHeI, nTriangleHeI, Xy_D, &
            CoordHeIXy_DI(:,1:nPointHeI), &
            iNodeTriangleHeI_II(:,1:nTriangleHeI), &
            iNode1, iNode2, iNode3, Area1, Area2, Area3, IsTriangleFound)
       if(IsTriangleFound) then
          Flux = Area1*FluxHeI_I(iNode1) + Area2*FluxHeI_I(iNode2) &
               + Area3*FluxHeI_I(iNode3)
       else
          Flux=0.0
       endif
    case('HeII')
       call find_triangle(&
            nPointHeII, nTriangleHeII, Xy_D, &
            CoordHeIIXy_DI(:,1:nPointHeII), &
            iNodeTriangleHeII_II(:,1:nTriangleHeII), &
            iNode1, iNode2, iNode3, Area1, Area2, Area3, IsTriangleFound)
       if(IsTriangleFound) then
          Flux = Area1*FluxHeII_I(iNode1) + Area2*FluxHeII_I(iNode2) &
               + Area3*FluxHeII_I(iNode3)
       else
          Flux=0.0
       endif
    end select
    !locate triangle containing input point
 


    !interpolate flux to input point


    
  end subroutine get_plas_resonant_scattering

  !============================================================================
  subroutine unit_test_plas_resonant_scattering
    use ModIoUnit,     ONLY: UnitTmp_
    implicit none
    real :: Xy_D(2),Flux
    integer :: iAlt,iSZA
    integer :: nAlt=100,nSZA=18
    !--------------------------------------------------------------------------
    
    ! test one point
    Xy_D(1)=100.0
    Xy_D(2)=400.0
    
    call get_plas_resonant_scattering(Xy_D,'LyBeta',Flux)
    
    write(*,*) 'For SZA=',Xy_D(1),' and Alt=',Xy_D(2),' Flux=',Flux

    call get_plas_resonant_scattering(Xy_D,'LyBeta',Flux)
    
    write(*,*) 'AGAIN!For SZA=',Xy_D(1),' and Alt=',Xy_D(2),' Flux=',Flux

    ! now test file
    open(UnitTmp_,FILE='StrobelLyBeta.dat')
    write(UnitTmp_,'(a)') &
         'VARIABLES = "Sza", "Alt", "Flux [photons/cm2/s]"'
    write(UnitTmp_,'(a,i3,a,i3,a)') 'Zone I=', nSZA, &
         ', J=', nAlt,', DATAPACKING=POINT'
    do iAlt=1,nAlt
       do iSZA=1,nSZA
          Xy_D(2) = real(iAlt-1)*10.0+90.0
          Xy_D(1) = real(iSZA-1)*5.0+90.0
          call get_plas_resonant_scattering(Xy_D,'LyBeta',Flux)
 !         write(*,"(100es18.10)") Xy_D(1),&
 !                Xy_D(2),Flux
          write(UnitTmp_,"(100es18.10)") Xy_D(1),&
                 Xy_D(2),Flux
       end do
    end do
    close(UnitTmp_)


    ! now test file
    open(UnitTmp_,FILE='StrobelLyAlpha.dat')
    write(UnitTmp_,'(a)') &
         'VARIABLES = "Sza", "Alt", "Flux [photons/cm2/s]"'
    write(UnitTmp_,'(a,i3,a,i3,a)') 'Zone I=', nSZA, &
         ', J=', nAlt,', DATAPACKING=POINT'
    do iAlt=1,nAlt
       do iSZA=1,nSZA
          Xy_D(2) = real(iAlt-1)*10.0+90.0
          Xy_D(1) = real(iSZA-1)*5.0+90.0
          call get_plas_resonant_scattering(Xy_D,'LyAlpha',Flux)
 !         write(*,"(100es18.10)") Xy_D(1),&
 !                Xy_D(2),Flux
          write(UnitTmp_,"(100es18.10)") Xy_D(1),&
                 Xy_D(2),Flux
       end do
    end do
    close(UnitTmp_)

    ! now test file
    open(UnitTmp_,FILE='StrobelHeI.dat')
    write(UnitTmp_,'(a)') &
         'VARIABLES = "Sza", "Alt", "Flux [photons/cm2/s]"'
    write(UnitTmp_,'(a,i3,a,i3,a)') 'Zone I=', nSZA, &
         ', J=', nAlt,', DATAPACKING=POINT'
    do iAlt=1,nAlt
       do iSZA=1,nSZA
          Xy_D(2) = real(iAlt-1)*10.0+90.0
          Xy_D(1) = real(iSZA-1)*5.0+90.0
          call get_plas_resonant_scattering(Xy_D,'HeI',Flux)
 !         write(*,"(100es18.10)") Xy_D(1),&
 !                Xy_D(2),Flux
          write(UnitTmp_,"(100es18.10)") Xy_D(1),&
                 Xy_D(2),Flux
       end do
    end do
    close(UnitTmp_)

    ! now test file
    open(UnitTmp_,FILE='StrobelHeII.dat')
    write(UnitTmp_,'(a)') &
         'VARIABLES = "Sza", "Alt", "Flux [photons/cm2/s]"'
    write(UnitTmp_,'(a,i3,a,i3,a)') 'Zone I=', nSZA, &
         ', J=', nAlt,', DATAPACKING=POINT'
    do iAlt=1,nAlt
       do iSZA=1,nSZA
          Xy_D(2) = real(iAlt-1)*10.0+90.0
          Xy_D(1) = real(iSZA-1)*5.0+90.0
          call get_plas_resonant_scattering(Xy_D,'HeII',Flux)
 !         write(*,"(100es18.10)") Xy_D(1),&
 !                Xy_D(2),Flux
          write(UnitTmp_,"(100es18.10)") Xy_D(1),&
                 Xy_D(2),Flux
       end do
    end do
    close(UnitTmp_)

  end subroutine unit_test_plas_resonant_scattering

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
end Module ModSeProduction
