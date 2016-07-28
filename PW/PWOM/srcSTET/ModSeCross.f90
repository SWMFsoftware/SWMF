! A module that contains all the routines for setting the cross sections 
Module ModSeCross
  private !except
  
  public :: cross
  public :: cross_jupiter
  !number of excited states
  integer, parameter :: NEI=10
  
  !number of major neutral species
  integer, parameter :: NMAJ = 3


  ! crossesction that are calculated here for other parts of the code
  real, allocatable,public :: SIGS(:,:), &
       SIGI(:,:,:), SIGA(:,:,:)

  ! WW-energy threshold for each excited state,species in eV
  ! WW,AO,OMEG,ANU,BB - revised excitation crossesction parameters from
  ! Green & Stolarski (1972) formula)
  ! AUTO autoionization coefs
  ! AK, AJ, TS, TA, TB, GAMS, GAMB:  Jackman et al (1977) ioniz. params
  real :: WW(NEI,NMAJ), AO(NEI,NMAJ), OMEG(NEI,NMAJ), &
       ANU(NEI,NMAJ), BB(NEI,NMAJ), AUTO(NEI,NMAJ), &
       THI(NEI,NMAJ),  AK(NEI,NMAJ),   AJ(NEI,NMAJ), &
       TS(NEI,NMAJ),   TA(NEI,NMAJ),   TB(NEI,NMAJ), &
       GAMS(NEI,NMAJ), GAMB(NEI,NMAJ)


    !
    !     O states:   1D,   1S, 3s5S, 3s3S, 3p5P, 3p3P, 3d3D, 3s'3D
    !     O2 states:   a,    b, AA'c,    B,  9.9, Ryds,  vib
    !     N2 states: ABW,   B',    C, aa'w,  1Pu,   b', Ryds,  vib
    !
    DATA WW  /1.96, 4.17, 9.29, 9.53,10.76,10.97,12.07,12.54, 0.,0., &
         0.98, 1.64, 4.50, 8.44, 9.90,13.50, 0.25, 0.00, 0.,0., &
         6.17, 8.16,11.03, 8.40,12.85,14.00,13.75, 1.85, 0.,0./
    DATA AO /.0100,.0042,.1793,.3565,.0817,.0245,.0293,.1221, 0.,0., &
         .0797,.0211,.0215,.3400,.0657,1.110,3.480, 0.00, 0.,0., &
         2.770,.1140,.1790,.0999,.8760,.6010,1.890,1.350, 0.,0./
    DATA OMEG/1.00, 1.00, 3.00, 0.75, 3.00, 0.85, 0.75, 0.75, 0.,0., &
         2.00, 2.00, 1.15, 0.75, 0.75, 0.75, 7.00, 0.00, 0.,0., &
         3.00, 3.00, 3.00, 1.00, 0.75, 0.75, 0.75, 8.00, 0.,0./
    DATA ANU /2.00, 1.04, 2.53, 0.54, 2.43, 2.87, 0.93, 0.72, 0.,0., &
         6.18, 4.14, 1.00, 1.05, 1.60, 3.00,10.87, 0.00, 0.,0., &
         4.53, 4.78, 4.32, 4.05, 1.47, 1.27, 3.00, 1.58, 0.,0./
    DATA BB / 1.00, 0.50, 1.02, 0.01, 4.19, 4.88, 0.66, 0.17, 0.,0., &
         0.53, 0.51, 0.98, 0.99, 1.86, 1.00, 1.00, 0.00, 0.,0., &
         1.42, 3.54,12.70, 5.20, 0.86, 0.45, 1.00, 1.00, 0.,0./
    DATA AUTO/30*0.0/
    DATA THI/13.60,16.90,18.50, 0.00, 0.00, 0.00, 0.00, 0.,0.,0., &
         12.10,16.10,16.90,18.20,20.00,23.00,37.00, 0.,0.,0., &
         15.58,16.73,18.75,22.00,23.60,40.00, 0.00, 0.,0.,0./
    DATA AK / 1.13, 1.25, 0.67, 0.00, 0.00, 0.00, 0.00, 0.,0.,0., &
         0.47, 1.13, 1.13, 1.01, 0.65, 0.95, 0.59, 0.,0.,0., &
         2.42, 1.06, 0.55, 0.37, 0.37, 0.53, 0.00, 0.,0.,0./
    DATA AJ / 1.81, 1.79, 1.78, 0.00, 0.00, 0.00, 0.00, 0.,0.,0., &
         3.76, 3.76, 3.76, 3.76, 3.76, 3.76, 3.76, 0.,0.,0., &
         1.74, 1.74, 1.74, 1.74, 1.74, 1.74, 0.00, 0.,0.,0./
    DATA TS / 6.41, 6.41, 6.41, 0.00, 0.00, 0.00, 0.00, 0.,0.,0., &
         1.86, 1.86, 1.86, 1.86, 1.86, 1.86, 1.86, 0.,0.,0., &
         4.71, 4.71, 4.71, 4.71, 4.71, 4.71, 0.00, 0.,0.,0./
    DATA TA /3450.,3450.,3450.,   0.,   0.,   0.,   0., 0.,0.,0., &
         1000.,1000.,1000.,1000.,1000.,1000.,1000., 0.,0.,0., &
         1000.,1000.,1000.,1000.,1000.,1000.,   0., 0.,0.,0./
    DATA TB/162.00,162.0,162.0, 0.00, 0.00, 0.00, 0.00, 0.,0.,0., &
         24.20,32.20,33.80,36.40,40.60,46.00,74.00, 0.,0.,0., &
         31.16,33.46,37.50,44.00,47.20,80.00, 0.00, 0.,0.,0./
    DATA GAMS/13.00,13.0,13.00, 0.00, 0.00, 0.00, 0.00, 0.,0.,0., &
         18.50,18.5,18.50,18.50,18.50,18.50,18.50, 0.,0.,0., &
         13.80,13.8,13.80,13.80,13.80,13.80, 0.00, 0.,0.,0./
    DATA GAMB/-.815,-.815,-.815,0.00, 0.00, 0.00, 0.00, 0.,0.,0., &
         12.10,16.10,16.90,18.2,20.30,23.00,37.00, 0.,0.,0., &
         15.58,16.73,18.75,22.0,23.60,40.00, 0.00, 0.,0.,0./
  
contains
  
  ! Subroutine CROSS, an adaptation from EXSECT to fit our "black
  !*  box" approach to the problem.
  !*  This subroutine can now handle a nonzero Emin (which it handled
  !*  incorrectly before): 3/29/95      (M W Liemohn)
  !
  ! Adapted from Banks & Nagy 2-stream code by Stan Solomon, 1988
  !
  ! Calculates electron impact cross sections
  !
  ! Definitions:
  ! SIGS   elastic cross sections for each species, energy; cm2
  ! PE     elastic backscatter probabilities for each species, energy
  ! PI     inelastic  "
  ! SIGA   energy loss cross section for each species, loss, energy; cm2
  ! WW     energy threshold for each excited state, species; eV
  ! WW, AO, OMEG, ANU, BB: revised excitation cross section parameters,
  !        from Green & Stolarski (1972) formula (W, A, omega, nu, gamma)
  ! AUTO   autoionization coefs (= 0 as autoion. included in ion xsects)
  ! THI    energy threshold for each ionized state, species; eV
  ! AK, AJ, TS, TA, TB, GAMS, GAMB:  Jackman et al (1977) ioniz. params
  ! ENER   energy grid; eV
  ! DEL    energy grid spacing; eV
  ! NNN    number of excited states for each species
  ! NINN   number of ionized states for each species
  ! NUM    number of points on elastic data trid for each species
  ! EC     data energy grid of elastic xsects and backscatter ratios
  !        for each species; eV
  ! CC     elastic xsects on data grid for each species, cm2
  ! CE     elastic backscat. probs on data grid for each species; cm2
  ! CI     inelastic "
  !
  ! Array dimensions:
  ! Elen   number of energy levels
  ! NMAJ   number of major species
  ! NEI    number of slots for excited and ionized states
  !
  !
  SUBROUTINE CROSS
    !
!    INCLUDE 'numbers.h'
    !
!    COMMON /CXSECT/ SIGS(NMAJ,Elen), &
!         SIGI(NMAJ,Elen,Elen), SIGA(NMAJ,Elen,Elen)
    !
!    COMMON /CXPARS/ WW(NEI,NMAJ), AO(NEI,NMAJ), OMEG(NEI,NMAJ), &
!         ANU(NEI,NMAJ), BB(NEI,NMAJ), AUTO(NEI,NMAJ), &
!         THI(NEI,NMAJ),  AK(NEI,NMAJ),   AJ(NEI,NMAJ), &
!         TS(NEI,NMAJ),   TA(NEI,NMAJ),   TB(NEI,NMAJ), &
!         GAMS(NEI,NMAJ), GAMB(NEI,NMAJ)
    !
    use ModSeGrid,only:nEnergy,del=>DeltaE_I,ener=>EnergyGrid_I, &
         Emin=>EnergyMin, BINNUM
!    COMMON /CENERGY/ENER(Elen), DEL(Elen), Emin, Jo
    !

    DIMENSION NNN(NMAJ), NINN(NMAJ), NUM(NMAJ), &
         EC(31,NMAJ), CC(31,NMAJ)
    !
    !      DATA QQN/6.51E-14/, NNN/8,7,7/, NINN/3,7,6/, NUM/31,28,28/
    !  Changed to remove the N2 vibrational state.
    DATA QQN/6.51E-14/, NNN/8,7,8/, NINN/3,7,6/, NUM/31,28,28/

    DATA EC /      1.00,     2.00,     4.00,     6.00,     8.00, &
         10.00,    12.00,    14.00,    16.00,    18.00, &
         20.00,    30.00,    40.00,    50.00,    60.00, &
         70.00,    80.00,    90.00,   100.00,   150.00, &
         200.00,   300.00,   500.00,  1000.00,  2000.00, &
         3000.00,  5000.00, 10000.00, 20000.00, 40000.00, &
         50000.00, &
         1.00,     2.00,     3.00,     5.00,     7.00, &
         10.00,    15.00,    20.00,    30.00,    40.00, &
         50.00,    70.00,   100.00,   150.00,   200.00, &
         300.00,   400.00,   500.00,   600.00,   700.00, &
         1000.00,  2000.00,  3000.00,  5000.00, 10000.00, &
         20000.00, 40000.00, 50000.00,     0.00,     0.00, &
         0.00, &
         1.00,     2.00,     2.50,     3.00,     4.00, &
         5.00,     6.00,     8.00,    10.00,    15.00, &
         20.00,    30.00,    40.00,    50.00,    70.00, &
         100.00,   200.00,   300.00,   500.00,   700.00, &
         1000.00,  2000.00,  3000.00,  5000.00, 10000.00, &
         20000.00, 40000.00, 50000.00,     0.00,     0.00, 0.0/
    DATA CC /  5.00E-16, 6.00E-16, 7.50E-16, 7.60E-16, 7.70E-16, &
         7.80E-16, 7.50E-16, 7.20E-16, 6.90E-16, 6.70E-16, &
         6.50E-16, 5.60E-16, 4.60E-16, 4.00E-16, 3.50E-16, &
         3.20E-16, 2.90E-16, 2.70E-16, 2.50E-16, 1.90E-16, &
         1.50E-16, 1.20E-16, 8.00E-17, 5.00E-17, 3.02E-17, &
         1.99E-17, 1.20E-17, 6.08E-18, 3.06E-18, 1.55E-18, &
         1.24E-18, &
         5.50E-16, 6.90E-16, 7.50E-16, 8.50E-16, 9.60E-16, &
         1.00E-15, 1.00E-15, 9.00E-16, 8.30E-16, 7.70E-16, &
         6.90E-16, 5.70E-16, 4.40E-16, 3.30E-16, 2.70E-16, &
         2.10E-16, 1.80E-16, 1.60E-16, 1.40E-16, 1.30E-16, &
         1.10E-16, 7.00E-17, 5.00E-17, 3.00E-17, 1.53E-17, &
         7.72E-18, 3.90E-18, 3.13E-18, 0.00E+00, 0.00E+00, &
         0.00E+00, &
         9.00E-16, 2.27E-15, 2.52E-15, 1.93E-15, 1.32E-15, &
         1.15E-15, 1.16E-15, 1.17E-15, 1.18E-15, 1.14E-15, &
         1.13E-15, 9.50E-16, 8.60E-16, 7.30E-16, 5.90E-16, &
         4.70E-16, 3.30E-16, 2.50E-16, 1.60E-16, 1.30E-16, &
         1.10E-16, 6.35E-17, 4.18E-17, 2.54E-17, 1.28E-17, &
         6.44E-18, 3.27E-18, 2.62E-18, 0.00E+00, 0.00E+00, 0.0/

    !allocate sig arrays if not already done
    if (.not.allocated(SIGS)) allocate(SIGS(NMAJ,nEnergy)) 
    if (.not.allocated(SIGI)) allocate(SIGI(NMAJ,nEnergy,nEnergy)) 
    if (.not.allocated(SIGA)) allocate(SIGA(NMAJ,nEnergy,nEnergy)) 
    

    !
    !
    ! Zero energy loss xsect and secondary production xsect arrays:
    !
    DO I3=1,nEnergy
       DO I2=1,nEnergy
          DO I1=1,NMAJ
             SIGA(I1,I2,I3)=0.0
             SIGI(I1,I2,I3)=0.0
          enddo
       enddo
    enddo
    !
    !
    ! Interpolate elastic cross sections and backscatter ratios:
    !
    DO  IJ=1,NMAJ
       DO   IV=1,nEnergy
          EX=ENER(IV)
          DO  II=1,NUM(IJ)
             IF (EC(II,IJ) .GT. EX) GOTO 60
          end DO
          SIGS(IJ,IV)=CC(NUM(IJ),IJ)*(EC(NUM(IJ),IJ)/EX)**0.8
          IF(IJ.EQ.1) &
               SIGS(IJ,IV)=CC(NUM(IJ),IJ)*(EC(NUM(IJ),IJ)/EX)**2
          cycle
60        I=II-1
          IF (I .LE. 0) THEN
             SIGS(IJ,IV)=CC(II,IJ)
          ELSE
             FAC = ALOG (EX/EC(I,IJ)) / ALOG (EC(II,IJ)/EC(I,IJ))
             SIGS(IJ,IV) = EXP (ALOG (CC(I,IJ)) &
                  + ALOG (CC(II,IJ)/CC(I,IJ)) * FAC)
          ENDIF
       enddo
    enddo
    !
    !
    ! Calculate electron impact excitation cross sections (put into SIGA):
    !
    DO  JY=1,nEnergy
       ETJ=ENER(JY)
       DO  I=1,NMAJ
          DO  J=1,NNN(I)
             ! difference between energy and threshold
             ETA = ETJ - WW(J,I)
             IF (ETA .GT. 0.) THEN
                !when energy exceeds the threshold 
                WE = WW(J,I) / ETJ
                SIGG = QQN * AO(J,I) * (WE**OMEG(J,I) / WW(J,I)**2) &
                     * (1.0 - WE**BB(J,I)) ** ANU(J,I)
                IF (SIGG .LT. 1.E-30) SIGG = 0.0
                !find closest energy bins bracketing (Energy-threshold)
                IE = INV (nEnergy,ETA,JY,ENER,Emin)
                IEE = IE - 1
                ! K is shifted index. Note that second index in SIGA is measured relative 
                ! to the index of the energy bin cooresponding to E-threshold
                K = JY - IE
                KK = JY - IEE
                IF (IE .EQ. JY) THEN
                   !when shifted index is at 0 (when energybin cooresponds to index)
                   IF (JY .EQ. 1) THEN
                      !special case when energy index is at bottom of energy grid
                      SIGA(I,1,JY)=SIGA(I,1,JY)+SIGG*WW(J,I)/(.5*DEL(JY))
                   ELSE
                      !exactly on an energy bin but not at the bottom
                      SIGA(I,1,JY)=SIGA(I,1,JY) &
                           +SIGG*WW(J,I)/(ETJ-ENER(JY-1))
                   END IF
                ELSE
                   IF (IE .LE. 1 .OR. ETA.LE.Emin) THEN
                      !special case when index is less than bottom or energy grid
                      ! when difference between threshold and energy is very small
                      SIGA(I,K,JY) = SIGA(I,K,JY) + SIGG
                   ELSE
                      ! usual case where you are interpolating between two bins
                      FF = (ENER(IE)-ETA) / (ENER(IE)-ENER(IEE))
                      FF = 1.0 - ABS(FF)
                      SIGA(I,K,JY) = SIGA(I,K,JY) + SIGG * FF
                      SIGA(I,KK,JY) = SIGA(I,KK,JY) + SIGG * (1.0-FF)
                   END IF
                END IF
             END IF
          enddo
       enddo
    enddo
    !
    !
    ! Loop over energy:
    !
    DO  JY=1,nEnergy
       !
       ETJ=ENER(JY)
       DETJ=DEL(JY)
       !
       !
       ! Loop over species:
       !
       DO  I=1,NMAJ
          !
          !
          ! Loop over ion states:
          !
          DO  ML=1,NINN(I)
             !
             !
             ! Calculate cross-section for production of secondaries into each
             ! bin from 1 to ITMAX and store in SIGI(II). Also store the average
             ! energy of the secondaries in T120:
             !
             WAG = THI(ML,I)
             TMAX = (ETJ-WAG) / 2.
             IF (TMAX .LE. 0.) cycle
             ITMAX = INV (nEnergy,TMAX,JY,ENER,Emin)
             IF (.not.(ITMAX.EQ.0 .OR. TMAX.LE.Emin))then
                TMT=0.
                I1=0
                DO  WHILE (TMT.EQ.0.)
                   I1=I1+1
                   TMT = ENER(I1)+DEL(I1)/2.0
                   IF (TMT.LE.Emin) TMT=0.
                   IF (I1.EQ.nEnergy .AND. TMT.EQ.0.) TMT=Emin
                enddo
                IF (TMAX .LT. TMT) TMT = TMAX
                SIGI(I,I1,JY) = SIGI(I,I1,JY) &
                     + SIGION(I,ML,ETJ,Emin,TMT,T120)
                TMT = ENER(I1) + DEL(I1) / 2.
                IF (TMAX .GT. TMT) THEN
                   IF (TMAX .LE. ENER(I1+1)) ITMAX = I1+1
                   DO  II=I1+1,ITMAX
                      E1 = ENER(II) - DEL(II) / 2.
                      E2 = E1 + DEL(II)
                      IF (E2 .GT. TMAX) E2 = TMAX
                      IF (E1 .LE. E2) SIGI(I,II,JY) = SIGI(I,II,JY) &
                           + SIGION(I,ML,ETJ,E1,E2,T120)
                   enddo
                ENDIF
             endif
             !
             ! Put degradation cross sections into SIGA
             !
             ETA=TMAX               ! Lowest primary energy after collision
             ETB=ETJ-WAG            ! Highest primary energy after collision
             ITMAX=BINNUM(ETA)
             ITMIN=BINNUM(ETB)
             IF (ITMAX.EQ.0 .OR. ETA.LE.Emin) THEN
                E1=ETB-Emin
                IF (E1.LT.0.) E1=0.
                E2=TMAX
                SIGA(I,JY,JY)=SIGA(I,JY,JY)+SIGION(I,ML,ETJ,E1,E2,T120)
                ITMAX=1
                TMAX=ETB-Emin
             END IF
             IF (ITMIN.EQ.0 .OR. ETB.LE.Emin) cycle
             DO  II=ITMIN,ITMAX,-1
                E1=ETB-(ENER(II)+.5*DEL(II))
                E2=E1+DEL(II)
                IF (E1.LT.0.) E1=0.
                IF (II.EQ.ITMAX) E2=TMAX
                IF (E2.LE.1.E-10) cycle
                SIGG=SIGION(I,ML,ETJ,E1,E2,T120)
                IF (II.EQ.JY) THEN
                   SIGA(I,1,JY)=SIGA(I,1,JY) &
                        +SIGG*(WAG+T120)/(ETJ-ENER(JY-1))
                ELSE
                   SIGA(I,JY-II,JY)=SIGA(I,JY-II,JY)+SIGG
                END IF
             enddo
             !
             !
          enddo           ! ion states
          !
       enddo           ! species
       !
    enddo            ! energy
    !

!    call plot_diffion_cross
    RETURN
  END SUBROUTINE CROSS
  !
  !
  !
  !
  ! Function SIGION calculates ionization cross section for species I,
  ! state ML, primary energy E, secondary energy from E1 to E2
  !
  FUNCTION SIGION(I,ML,E,E1,E2,T12)
    !
    PARAMETER (NMAJ=3)
    PARAMETER (NEI=10)
    !
!    COMMON /CXPARS/ WW(NEI,NMAJ), AO(NEI,NMAJ), OMEG(NEI,NMAJ), &
!         ANU(NEI,NMAJ), BB(NEI,NMAJ), AUTO(NEI,NMAJ), &
!         THI(NEI,NMAJ),  AK(NEI,NMAJ),   AJ(NEI,NMAJ), &
!         TS(NEI,NMAJ),   TA(NEI,NMAJ),   TB(NEI,NMAJ), &
!         GAMS(NEI,NMAJ), GAMB(NEI,NMAJ)
    !
    DATA QQ/1.E-16/
    !
    !
    IF (.not.(E .LE. THI(ML,I))) then
       IF (.not.(E2.LE.E1)) then
          !
          AK1=AK(ML,I)
          AJ1=AJ(ML,I)
          TS1=TS(ML,I)
          TA1=TA(ML,I)
          TB1=TB(ML,I)
          GAMS1=GAMS(ML,I)
          GAMB1=GAMB(ML,I)
          S=QQ*AK1*ALOG(E/AJ1)
          A=S/E
          TZ=TS1-TA1/(E+TB1)
          GG=(GAMS1*E)/(E+GAMB1)
          TTL=(E-THI(ML,I))/2.0
          TTL1=TTL-0.01
          IF(.not.(E1.GE.TTL1)) then
             IF(E2.GE.TTL)E2=TTL
             ABB=(E2-TZ)/GG
             ABC=(E1-TZ)/GG
             AL2=GG*GG*(ABB*ABB+1.0)
             AL1=GG*GG*(ABC*ABC+1.0)
             ABB=ATAN(ABB)-ATAN(ABC)
             T12=TZ+0.5*GG*(ALOG(AL2)-ALOG(AL1))/ABB
             SIGION=A*GG*ABB
             RETURN
          endif
       endif
    endif
    !
    SIGION=0.0
    T12=0.0
    RETURN
    !
  END FUNCTION SIGION
  !
  !
  !
  !
  ! Function INV finds the bin number closest to energy ETA on grid ENER.
  ! Bin INV or INV-1 will contain ETA.
  !
  FUNCTION INV (nEnergy,ETA,JY,ENER,Emin)
    !
!    INCLUDE 'numbers.h'
    DIMENSION ENER(nEnergy)
    !
    IF (ETA .LE. Emin) THEN
       INV = 0
    ELSE
       DO  IV=1,JY
          IF (ETA .LE. ENER(IV)) then 
             INV=IV
             return
          else
             INV=JY
          endif
       enddo
    ENDIF
    !
    RETURN
  END FUNCTION INV

  !=============================================================================
  subroutine cross_jupiter(nNeutralSpecies)
    
    use ModSeGrid,only:nEnergy,DeltaE_I,EnergyGrid_I, &
         EnergyMin, BINNUM
    !jupiter
    integer,parameter  :: H2_=1, He_=2, H_=3, CH4_=4
    
    !allocate sig arrays if not already done
    if (.not.allocated(SIGS)) allocate(SIGS(nNeutralSpecies,nEnergy)) 
    if (.not.allocated(SIGI)) allocate(SIGI(nNeutralSpecies,nEnergy,nEnergy)) 
    if (.not.allocated(SIGA)) allocate(SIGA(nNeutralSpecies,nEnergy,nEnergy)) 
    
    ! currently just set crossections to zero. 
    SIGS(:,:)=0.0
    SIGI(:,:,:)=0.0
    SIGA(:,:,:)=0.0

    write(*,*) 'Neutrals: ',nNeutralSpecies
    
    do iNeutral = 1, nNeutralSpecies

       write(*,*) 'getting H2 cross'
       call get_excitation_crossection('H2', SIGA(H2_,:,:))
       call get_crossection_diffion('H2', SIGI(H2_,:,:), SIGA(H2_,:,:))
       call get_scattering_crossection('H2',SIGS(H2_,:))
       write(*,*) 'getting He cross'
       call get_excitation_crossection('He', SIGA(He_,:,:))
       call get_crossection_diffion('He', SIGI(He_,:,:), SIGA(He_,:,:))
       call get_scattering_crossection('He',SIGS(He_,:))
       write(*,*) 'getting H cross'
       call get_excitation_crossection('H', SIGA(H_,:,:))
       call get_crossection_diffion('H',  SIGI(H_,:,:), SIGA(H_,:,:))
       call get_scattering_crossection('H',SIGS(H_,:))
       write(*,*) 'getting CH4 cross'
       call get_excitation_crossection('CH4', SIGA(CH4_,:,:))
       call get_crossection_diffion('CH4',SIGI(CH4_,:,:), SIGA(CH4_,:,:))
       call get_scattering_crossection('CH4',SIGS(CH4_,:))
       
    end do
    
    ! plot the differential ion crossection (for testing)
    call plot_diffion_cross
  end subroutine cross_jupiter

  !=============================================================================
  subroutine get_excitation_crossection(NameNeutralSpecies,SigA)
    use ModSeGrid,only:nEnergy,DeltaE_I,EnergyGrid_I,EnergyMin,BINNUM
    
    character(len=*), intent(in) :: NameNeutralSpecies
    real, intent(out) :: SigA(nEnergy,nEnergy)
    integer :: nStates

    real, allocatable :: SigExcitation(:,:),Threshold(:)

    integer :: iEnergy,iState, iBinHigh,iBinLow,nBinShiftHigh,nBinShiftLow
    real :: Energy, DeltaE, Sigma, UpperBinFrac
    !
    !
    ! Calculate electron impact excitation cross sections (put into SIGA):
    !
    write(*,*) 'Getting excitation crossections for ',NameNeutralSpecies
    
    select case(NameNeutralSpecies)
    case('O')
       write(*,*) 'WARNING: Species ',NameNeutralSpecies&
            ,' not yet supported. Using 0 for crossection.' 
       nStates = 0
    case('N2')
       write(*,*) 'WARNING: Species ',NameNeutralSpecies&
            ,' not yet supported. Using 0 for crossection.' 
       nStates = 0
    case('O2')
       write(*,*) 'WARNING: Species ',NameNeutralSpecies&
            ,' not yet supported. Using 0 for crossection.' 
       nStates = 0
    case('H2')
       nStates = 12
       if(.not.allocated(SigExcitation)) allocate(SigExcitation(nEnergy,nStates))
       if(.not.allocated(Threshold)) allocate(Threshold(nStates))
       call read_excitation_crossection(NameNeutralSpecies,nStates,Threshold,SigExcitation)
    case('H')
       nStates = 6
       if(.not.allocated(SigExcitation)) allocate(SigExcitation(nEnergy,nStates))
       if(.not.allocated(Threshold)) allocate(Threshold(nStates))
       call read_excitation_crossection(NameNeutralSpecies,nStates,Threshold,SigExcitation)
   case('CH4')
       nStates = 5
       if(.not.allocated(SigExcitation)) allocate(SigExcitation(nEnergy,nStates))
       if(.not.allocated(Threshold)) allocate(Threshold(nStates))
       call read_excitation_crossection(NameNeutralSpecies,nStates,Threshold,SigExcitation)
    case('He')
       nStates = 18
       if(.not.allocated(SigExcitation)) allocate(SigExcitation(nEnergy,nStates))
       if(.not.allocated(Threshold)) allocate(Threshold(nStates))
       call read_excitation_crossection(NameNeutralSpecies,nStates,Threshold,SigExcitation)
    case default
       write(*,*) 'WARNING: Species ',NameNeutralSpecies&
            ,' not yet supported. Using 0 for crossection.' 
       nStates = 0
       
    end select

    DO  iEnergy=1,nEnergy
       Energy=EnergyGrid_I(iEnergy)
       DO  iState=1,nStates ! was J
          ! difference between energy and threshold
          DeltaE = Energy - Threshold(iState) ! was ETA = ETJ-WW
          IF (DeltaE .GT. 0.) THEN
             !when energy exceeds the threshold
             ! *************************
             !Calculate crossection here:
             ! *************************
             ! GSinv_eps = Threshold(iState) / Energy
             ! from Green and Stolarski, 1972, eqs 1 & 7
             ! Sigma = GSq0 * GSA(iState) * (GSinv_eps**GSomega(iState,iSpecies) / &
             ! Threshold(iState)**2) * (1.0 - GSinv_eps**GSgamma(iState,iSpecies)) &
             ! ** GSnu(iState,iSpecies)
             Sigma = SigExcitation(iEnergy,iState)
             IF (Sigma .LT. 1.E-30) Sigma = 0.0
             !find closest energy bins bracketing (Energy-threshold)
             !could prob use iBin = BINNUM(DeltaE) but would make a small difference
             !in interpolation
             iBinHigh = INV(nEnergy,DeltaE,iEnergy,EnergyGrid_I,Emin)
             iBinLow  = iBinHigh - 1
             ! Find shifted indices. Note that second index in SIGA is measured relative 
             ! to the index of the energy bin cooresponding to E-threshold
             nBinShiftHigh = iEnergy - iBinHigh
             nBinShiftLow  = iEnergy - iBinLow
             IF (iBinHigh .EQ. iEnergy) THEN
                !when shifted index is at 0 (when energybin corresponds to index)
                IF (iEnergy .EQ. 1) THEN
                   !special case when energy index is at bottom of energy grid
                   SigA(1,iEnergy)=SigA(1,iEnergy)+ &
                        Sigma*Threshold(iState)/(.5*DeltaE_I(iEnergy))
                ELSE
                   !exactly on an energy bin but not at the bottom
                   SigA(1,iEnergy)=SigA(1,iEnergy) &
                        +Sigma*Threshold(iState)/(Energy-EnergyGrid_I(iEnergy-1))
                END IF
             ELSE
                IF (iBinHigh .LE. 1 .OR. DeltaE.LE.Emin) THEN
                   !special case when index is less than bottom or energy grid
                   ! when difference between threshold and energy is very small
                   SigA(nBinShiftHigh,iEnergy) = SigA(nBinShiftHigh,iEnergy) &
                        + Sigma
                ELSE
                   ! usual case where you are interpolating between two bins
                   UpperBinFrac = (EnergyGrid_I(iBinHigh)-DeltaE) / &
                        (EnergyGrid_I(iBinHigh)-EnergyGrid_I(iBinLow))
                   UpperBinFrac = 1.0 - ABS(UpperBinFrac)
                   SigA(nBinShiftHigh,iEnergy) = SigA(nBinShiftHigh,iEnergy) &
                        + Sigma * UpperBinFrac
                   SigA(nBinShiftLow,iEnergy)  = SigA(nBinShiftLow,iEnergy) &
                        + Sigma * (1.0-UpperBinFrac)
                END IF
             END IF
          END IF
       enddo
    enddo
  end subroutine get_excitation_crossection
  
  !=============================================================================
  subroutine get_scattering_crossection(NameNeutralSpecies,SigS)
    use ModSeGrid,only:nEnergy,DeltaE_I,EnergyGrid_I,EnergyMin,BINNUM

    character(len=*), intent(in) :: NameNeutralSpecies
    real, intent(out) :: SigS(nEnergy)

    real, allocatable :: SigTotalI(:)
    integer :: iEnergy,iEnergySec

    real :: NewPrimaryELow,NewPrimaryEHigh,EnergySecondary,SIGG
    integer ::iEnergyLow,iEnergyHigh,iNewEnergy

    select case(NameNeutralSpecies)
    case('O')
       write(*,*) 'WARNING: Species ',NameNeutralSpecies&
            ,' not yet supported. Using 0 for crossection.' 
       !placeholder set to 0
       SigS(:) = 0.0
    case('N2')
       write(*,*) 'WARNING: Species ',NameNeutralSpecies&
            ,' not yet supported. Using 0 for crossection.' 
       !placeholder set to 0
       SigS(:) = 0.0
    case('O2')
       write(*,*) 'WARNING: Species ',NameNeutralSpecies&
            ,' not yet supported. Using 0 for crossection.' 
       !placeholder set to 0
       SigS(:) = 0.0
    case('H2')
       call read_scattering_crossection(NameNeutralSpecies,SigS)
    case('H')
       call read_scattering_crossection(NameNeutralSpecies,SigS)       
    case('CH4')
       ! crossections from Shirai et al, 2002
       ! Analytic Cross Sections for Electron Collisions with Hydrocarbons
       Sig0 = 1.0e-16  ! cm**2
       ER = 1.361e-2 ! keV
       
       a1 = 2.93e-2
       a2 = -1.03
       a3 = 3.28e2
       a4 = 2.19
       a5 = 5.4e-3
       a6 = 7.92e-1

       do iEnergy = 1,nEnergy
          E1 = 1.0e-3*EnergyGrid_I(iEnergy) ! convert from eV to keV
          
          f1 = Sig0*a1*(E1/ER)**a2
          f12 = Sig0*a3*(E1/ER)**a4
          f2 = f12/(1.0+(E1/a5)**(a4+a6))
          
          SigS(iEnergy) = f1 + f2
       end do
       
    case('He')
       call read_scattering_crossection(NameNeutralSpecies,SigS)
       
    case default
       write(*,*) 'WARNING: Species ',NameNeutralSpecies&
            ,' not yet supported. Using 0 for crossection.' 
       SigS(:) = 0.0
    end select

    
  end subroutine get_scattering_crossection
  
 !=============================================================================
  subroutine get_crossection_diffion(NameNeutralSpecies,SigDiffI,SigDiffA)
    use ModSeGrid,only:nEnergy,DeltaE_I,EnergyGrid_I,EnergyMin,BINNUM

    character(len=*), intent(in) :: NameNeutralSpecies
    real, intent(out) :: SigDiffI(nEnergy,nEnergy), SigDiffA(nEnergy,nEnergy)

    real :: Ethreshold,Ebar,OpalCoef
    real :: S_kr,t_kr,coef_kr,coef2_kr,w_kr,dfdw_kr,Term1_kr,Term2_kr,Term3_kr
    real, allocatable :: SigTotalI(:)
    integer :: iEnergy,iEnergySec

    real :: NewPrimaryELow,NewPrimaryEHigh,EnergySecondary,SIGG
    integer ::iEnergyLow,iEnergyHigh,iNewEnergy
   
    !allocate arrays to hold the total ionization crossection
    if(.not.allocated(SigTotalI)) allocate(SigTotalI(nEnergy))

    !for testing
    SigTotalI(:) = 0.0

    write(*,*) 'Getting ',NameNeutralSpecies
    
    select case(NameNeutralSpecies)
    case('O')
       write(*,*) 'WARNING: Species ',NameNeutralSpecies&
            ,' not yet supported. Using 0 for crossection.' 
       !placeholder set to 0 
       SigTotalI(:)=0.0
       Ebar=1
       Ethreshold=0.0
    case('N2')
       call read_total_ionization_crossection(NameNeutralSpecies,SigTotalI)
       Ethreshold = 15.6
       Ebar = 13.0
    case('O2')
       call read_total_ionization_crossection(NameNeutralSpecies,SigTotalI)
       Ethreshold = 12.2
       Ebar = 17.4
    case('H2')
       call read_total_ionization_crossection(NameNeutralSpecies,SigTotalI)
       Ethreshold = 15.4
       Ebar = 8.3
    case('H')
       call read_diff_ionization_crossection(NameNeutralSpecies,SigDiffI)
 !      SigDiffI(:,:) = 0.0

   case('CH4')
       call read_total_ionization_crossection(NameNeutralSpecies,SigTotalI)
       Ethreshold = 13.0
       Ebar = 7.3
    case('He')
       call read_total_ionization_crossection(NameNeutralSpecies,SigTotalI)
       Ethreshold = 24.6
       Ebar = 15.8
    case default
       write(*,*) 'WARNING: Species ',NameNeutralSpecies&
            ,' not yet supported. Using 0 for crossection.' 
       SigTotalI(:)=0.0
       Ebar=1
       Ethreshold=0.0
    end select
    
    !loop over primary and secondary energies to fill diff crossection
    do iEnergy = 1,nEnergy

       if (NameNeutralSpecies.NE.'H') then
          !use the paper by Opal et al 1971 to set
          ! the differential ionization crossection
          if(EnergyGrid_I(iEnergy)>=Ethreshold) then
             !set the Opal coef
             OpalCoef = SigTotalI(iEnergy)&
                  /(Ebar * atan((EnergyGrid_I(iEnergy)-Ethrehold)/(2.0*Ebar)))
          else
             ! when energy of primary is below threshold make sure coef is zero
             OpalCoef=0.0
          end if
          do iEnergySec=1,nEnergy
             SigDiffI(iEnergySec,iEnergy) = &
                  OpalCoef / (1.0+(EnergyGrid_I(iEnergySec)/Ebar)**2.0)
             !SigDiffI(iEnergySec,iEnergy) = 0.0
          end do
       endif

       ! Put degradation cross sections into SIGA
       ! JY is energy loop iEnergy
       ! ETJ is EnergyGrid_I(iEnergy)

       ! Lowest primary energy after collision
       NewPrimaryELow = (EnergyGrid_I(iEnergy)-Ethreshold) / 2.
       ! Highest primary energy after collision
       NewPrimaryEHigh = EnergyGrid_I(iEnergy)-Ethreshold
       
       IF (NewPrimaryELow .LE. 0.) cycle
       
       iEnergyLow=BINNUM(NewPrimaryELow)
       iEnergyHigh=BINNUM(NewPrimaryEHigh)
       
       ! if NewPrimaryELow is below the EnergyMin then say that the primary
       ! loses all of its energy for all New Primary Energies from
       ! NewPrimaryELow to EnergyMin, which means Secondary Energies from
       ! NewPrimaryEHigh-EnergyMin to NewPrimaryEHigh-NewPrimaryELow,
       ! or NewPrimaryELow (because NewPrimaryELow*2 = NewPrimaryEHigh)
       
       IF (iEnergyLow.EQ.0 .OR. NewPrimaryELow.LE.EnergyMin) THEN
          !          E1Secondary=NewPrimaryEHigh-EnergyMin
          !          IF (E1Secondary.LT.0.) E1Secondary=0.
          !          E2Secondary=NewPrimaryELow
          ! Not sure about the right way to convert this to a single secondary energy...
          EnergySecondary=NewPrimaryELow
          iEnergySec=BINNUM(EnergySecondary)
          SigDiffA(iEnergy,iEnergy) = SigDiffA(iEnergy,iEnergy) + &
               SigDiffI(iEnergySec,iEnergy)
          iEnergyLow=1
          NewPrimaryELow=NewPrimaryEHigh-EnergyMin
       END IF
       IF (iEnergyHigh.EQ.0 .OR. NewPrimaryEHigh.LE.EnergyMin) cycle
       DO  iNewEnergy=iEnergyHigh,iEnergyLow,-1
          ! Secondary energy = new total energy - new primary energy
          EnergySecondary = NewPrimaryEHigh - EnergyGrid_I(iNewEnergy)
          IF (iNewEnergy.EQ.iEnergyLow) EnergySecondary=NewPrimaryELow
          IF (EnergySecondary.LE.1.E-10) cycle
          iEnergySec = BINNUM(EnergySecondary)
          SIGG=SigDiffI(iEnergySec,iEnergy)
          IF (iNewEnergy.EQ.iEnergy) THEN
             SigDiffA(1,iEnergy)=SigDiffA(1,iEnergy) &
                  + SIGG*(Ethreshold+EnergySecondary) / &
                  (EnergyGrid_I(iEnergy)-EnergyGrid_I(iEnergy-1))
          ELSE
             ! Put secondary energy into energy
             ! SigDiffA(PrimaryEnergyLoss,OriginalPrimaryEnergy)
             SigDiffA(iEnergy-iNewEnergy,iEnergy) = &
                  SigDiffA(iEnergy-iNewEnergy,iEnergy)+SIGG
          END IF
       enddo
       !
       !
    end do
    
    !deallocate to save memory
    deallocate(SigTotalI)

  end subroutine get_crossection_diffion

  !=============================================================================
  ! crossection for H taken from Shyn 1992
  subroutine read_diff_ionization_crossection(NameSpecies,SigDiffI)
    use ModSeGrid,only:nEnergy,DeltaE_I,EnergyGrid_I
    use ModIoUnit, ONLY : UnitTmp_
    use ModInterpolate, ONLY: bilinear
    
    character(len=*), intent(in) :: NameSpecies
    real :: SigDiffI(nEnergy,nEnergy)
    character(len=100) :: DatafileName
    character(len=31) :: TmpStr
    integer, parameter :: DataLen = 126, nE1 = 6, nE2 = 21
    integer :: iEnergy1,iEnergy2
    real :: DataArray(3,DataLen)
    real :: Energy1Array(nE1,nE2),Energy2Array(nE1,nE2),CrossSecArray(nE1,nE2)
    
    write(*,*) 'getting total ionization crossection for species ', &
         NameSpecies

    write(DatafileName,"(3a)") 'PW/DiffIon',NameSpecies,'.dat'

    open(UnitTmp_,FILE=DatafileName,STATUS='OLD')
    
    read(UnitTmp_,*) TmpStr

    ! Input energy arrays in eV and differential crossections in cm2
    read(UnitTmp_,*) DataArray
    
    close(UnitTmp_)

    SigDiffI(:,:) = 0.0

    do i=1,nE1
       Energy1Array(i,:)  = DataArray(1,(i-1)*nE2+1:i*nE2)
       Energy2Array(i,:)  = DataArray(2,(i-1)*nE2+1:i*nE2)
       if (i<nE1-1) then
          CrossSecArray(i,:) = DataArray(3,(i-1)*nE2+1:i*nE2)*1.e-17
       else
          CrossSecArray(i,:) = DataArray(3,(i-1)*nE2+1:i*nE2)*1.e-18
       endif
    enddo

    do iEnergy1 = 1,nEnergy
       do iEnergy2 = 1,nEnergy
          if (EnergyGrid_I(iEnergy1) < Energy1Array(1,1) &
               .or. EnergyGrid_I(iEnergy1) > Energy1Array(nE1,1)  &
               .or. EnergyGrid_I(iEnergy2) < Energy2Array(1,1)  &
               .or. EnergyGrid_I(iEnergy2) > Energy2Array(1,nE2)) then
             ! if outside of data range set crossection to 0
             SigDiffI(iEnergy2,iEnergy1)=0.0
          else
             ! when inside of data range interpolate
             SigDiffI(iEnergy2,iEnergy1) = &
                  bilinear(CrossSecArray(:,:),1,nE1,1,nE2, &
                  [EnergyGrid_I(iEnergy1),EnergyGrid_I(iEnergy2)], &
                  Energy1Array(:,1),Energy2Array(1,:))
          end if
       end do
    end do
    
  end subroutine read_diff_ionization_crossection
  !=============================================================================
  subroutine read_total_ionization_crossection(NameSpecies,SigTotalI)
    use ModSeGrid,only:nEnergy,DeltaE_I,EnergyGrid_I
    use ModIoUnit, ONLY : UnitTmp_
    use ModInterpolate, ONLY: linear
    
    character(len=*), intent(in) :: NameSpecies
    real :: SigTotalI(:)  ! should be allocated and passed in
    character(len=100) :: DatafileName
    integer :: DataLen
    real, allocatable :: EnergyArray(:),CrossSecArray(:)
    
    write(*,*) 'getting total ionization crossection for species ', &
         NameSpecies

    write(DatafileName,"(3a)") 'PW/IonCross',NameSpecies,'.dat'

    open(UnitTmp_,FILE=DatafileName,STATUS='OLD')
    
    read(UnitTmp_,*) DataLen

    if(.not.allocated(EnergyArray)) allocate(EnergyArray(DataLen))
    if(.not.allocated(CrossSecArray)) allocate(CrossSecArray(DataLen))
    
    ! Input energy array in eV and total crossections in cm2
    read(UnitTmp_,*) EnergyArray
    read(UnitTmp_,*) CrossSecArray
    
    close(UnitTmp_)

    SigTotalI(:) = 0.0

    do iEnergy = 1,nEnergy
       if (EnergyGrid_I(iEnergy) < EnergyArray(1) &
            .or. EnergyGrid_I(iEnergy) > EnergyArray(DataLen)) then
          ! if outside of data range set crossection to 0
          SigTotalI(iEnergy)=0.0
       else
          ! when inside of data range interpolate
          SigTotalI(iEnergy) = linear(CrossSecArray(:),1,DataLen, &
               EnergyGrid_I(iEnergy),EnergyArray(:))
       end if
    end do
    
    deallocate(EnergyArray)
    deallocate(CrossSecArray)

  end subroutine read_total_ionization_crossection
  !============================================================================
  subroutine read_excitation_crossection(NameSpecies,nStates,Threshold,SigEx)
    use ModSeGrid,only:nEnergy,DeltaE_I,EnergyGrid_I
    use ModIoUnit, ONLY : UnitTmp_
    use ModInterpolate, ONLY: linear
    
    character(len=*), intent(in) :: NameSpecies
    integer, intent(in) :: nStates
    !real :: Threshold(:),SigEx(:,:)  ! should be allocated and passed in
    real :: Threshold(nStates),SigEx(nEnergy,nStates)
    character(len=100) :: DatafileName
    integer :: DataLen,iState,iEnergy
    real, allocatable :: EnergyArray(:),CrossSecArray(:)
    
    write(*,*) 'getting excitation/absorption crossections for species ', &
         NameSpecies

    write(DatafileName,"(3a)") 'PW/',NameSpecies,'Across.dat'

    SigEx(:,:) = 0.0

    open(UnitTmp_,FILE=DatafileName,STATUS='OLD')

    do iState=1,nStates
       read(UnitTmp_,*) Threshold(iState)
       read(UnitTmp_,*) DataLen
       
       if(.not.allocated(EnergyArray)) allocate(EnergyArray(DataLen))
       if(.not.allocated(CrossSecArray)) allocate(CrossSecArray(DataLen))
    
       ! Input energy array in eV and excitation crossections in cm2
       read(UnitTmp_,*) EnergyArray
       read(UnitTmp_,*) CrossSecArray

       do iEnergy = 1,nEnergy
          if (EnergyGrid_I(iEnergy) < EnergyArray(1) &
               .or. EnergyGrid_I(iEnergy) > EnergyArray(DataLen)) then
             ! if outside of data range set crossection to 0
             SigEx(iEnergy,iState)=0.0
          else
             ! when inside of data range interpolate
             SigEx(iEnergy,iState) = linear(CrossSecArray(:),1,DataLen, &
                  EnergyGrid_I(iEnergy),EnergyArray(:))
          end if
       end do
    
       deallocate(EnergyArray)
       deallocate(CrossSecArray)

    end do
    
    close(UnitTmp_)
  end subroutine read_excitation_crossection
  !============================================================================
  subroutine read_scattering_crossection(NameSpecies,SigS)
    use ModSeGrid,only:nEnergy,DeltaE_I,EnergyGrid_I
    use ModIoUnit, ONLY : UnitTmp_
    use ModInterpolate, ONLY: linear
    
    character(len=*), intent(in) :: NameSpecies
    real :: SigS(:)  ! should be allocated and passed in
    character(len=100) :: DatafileName
    integer :: DataLen
    real, allocatable :: EnergyArray(:),CrossSecArray(:)
    
    write(*,*) 'getting elastic scattering crossection for species ', &
         NameSpecies

    write(DatafileName,"(3a)") 'PW/',NameSpecies,'Scross.dat'

    open(UnitTmp_,FILE=DatafileName,STATUS='OLD')
    
    read(UnitTmp_,*) DataLen

    if(.not.allocated(EnergyArray)) allocate(EnergyArray(DataLen))
    if(.not.allocated(CrossSecArray)) allocate(CrossSecArray(DataLen))
    
    ! Input energy array in eV and elastic scattering crossections in cm2
    read(UnitTmp_,*) EnergyArray
    read(UnitTmp_,*) CrossSecArray
    
    close(UnitTmp_)

    SigS(:) = 0.0

    do iEnergy = 1,nEnergy
       if (EnergyGrid_I(iEnergy) < EnergyArray(1) &
            .or. EnergyGrid_I(iEnergy) > EnergyArray(DataLen)) then
          ! if outside of data range set crossection to 0
          SigS(iEnergy)=0.0
       else
          ! when inside of data range interpolate
          SigS(iEnergy) = linear(CrossSecArray(:),1,DataLen, &
               EnergyGrid_I(iEnergy),EnergyArray(:))
       end if
    end do
    
    deallocate(EnergyArray)
    deallocate(CrossSecArray)

  end subroutine read_scattering_crossection
  !============================================================================
  ! plot differential ionization crossection
  subroutine plot_diffion_cross
    use ModSeGrid,     ONLY: nEnergy, EnergyGrid_I
    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModPlanetConst, ONLY: Planet_, NamePlanet_I
    
    real, allocatable   :: Coord_DII(:,:,:), PlotState_IIV(:,:,:)

    real    :: time=0
    integer :: nStep =0
    
    !grid parameters
    integer, parameter :: nDim =2, E1_=1, E2_=2
    
    ! planet-specific values
    integer :: nVar, nNeutral
    character(len=100) :: NamePlotVar

    integer :: iNeutral
    
    character(len=100) :: NamePlot = 'DiffIonCross.out'
    
    character(len=*),parameter :: NameHeader='SE output iono'
    character(len=5) :: TypePlot='ascii'
    integer :: iEnergyPrimary,iEnergySecondary

    logical :: IsFirstCall=.true.
    
    !--------------------------------------------------------------------------

    select case(NamePlanet_I(Planet_))
    case('EARTH')
       nVar=3
       nNeutral=3
       NamePlotVar= &
            'Es[eV] Ep[eV]  sigmaO[/cc/eV] sigmaO2[/cc/eV] sigmaN2[/cc/eV] g r'

    case('JUPITER')
       nVar=4
       nNeutral=4
       NamePlotVar= 'Es[eV] Ep[eV]  sigmaH2[/cc/eV] sigmaHe[/cc/eV] sigmaH[/cc/eV] sigmaCH4[/cc/eV] g r'

    end select
    
    allocate(Coord_DII(nDim,nEnergy,nEnergy),PlotState_IIV(nEnergy,nEnergy,nVar))
    
    !    do iLine=1,nLine
    PlotState_IIV = 0.0
    Coord_DII     = 0.0
    
    !Set values
    do iEnergyPrimary=1,nEnergy
       do iEnergySecondary=1,nEnergy
          Coord_DII(E1_,iEnergySecondary,iEnergyPrimary) = EnergyGrid_I(iEnergySecondary)             
          Coord_DII(E2_,iEnergySecondary,iEnergyPrimary) = EnergyGrid_I(iEnergyPrimary)             
          ! set plot state
          do iNeutral=1,nNeutral
             PlotState_IIV(iEnergySecondary,iEnergyPrimary,iNeutral) = &
                  SIGI(iNeutral,iEnergySecondary,iEnergyPrimary)
          enddo
       enddo
    enddo
    
    !Plot crossection
    if(IsFirstCall) then
       call save_plot_file(NamePlot, TypePositionIn='rewind', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
            nDimIn=nDim,CoordIn_DII=Coord_DII,                &
            VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
       IsFirstCall = .false.
    else
       call save_plot_file(NamePlot, TypePositionIn='append', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
            nDimIn=nDim,CoordIn_DII=Coord_DII,                &
            VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
    endif
    
    deallocate(Coord_DII, PlotState_IIV)
  end subroutine plot_diffion_cross
    

end Module ModSeCross
