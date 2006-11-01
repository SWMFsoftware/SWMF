C
C
C
C
C Subroutine EPHOTO calculates the photoelectron production spectrum
C
C Adapted from Banks & Nagy 2-stream input code by Stan Solomon, 6/88
C
C Uses continuously variable energy grid
C 3 species:  O, O2, N2
C
C Supplied by calling routine:
C WAVE1   wavelength array, upper bound, Angstroms
C WAVE2   wavelength array, lower bound, Angstroms
C SFLUX   solar flux array, photons cm-2 sec-1
C ZZ      altitude array, cm above earth
C ZMAJ    density array for species O, O2, N2, altitude, cm-3
C ZNO     density of NO at each altitude, cm-3
C ZCOL    slant column density array for species O, O2, N2, altitude, cm-2
C ENER    energy grid for photoelectrons, eV
C DEL     array of energy grid increments, eV
C
C Calculated by subroutine:
C PESPEC  photoelectron production spectrum for each altitude, cm-3 s-1
C PHOTOI  photoionization rates for each state, species, altitude, cm-3 s-1
C PHOTOD  photodissociation/exc. rates for state, species, alt., cm-3 s-1
C PHONO   photoionization/dissociation/excitation rates for NO, cm-3 s-1
C
C Other definitions:
C SPECT   photoelectron production spectrum at current z, cm-3 s-1
C FLUX    solar flux at altitude
C TPOT    ionization potentials for each species, state, eV
C SIGABS  photoabsorption cross sections, O, O2, N2; cm2
C SIGION  photoionization cross sections, O, O2, N2; cm2
C SIGAO, SIGAO2, SIGAN2, SIGIO, SIGIO2, SIGIN2; cross section data arrays
C NNN     number of states for each species
C PROB    branching ratios for each state, species, and wavelength interval
C PROBO, PROBO2, PROBN2; branching ratio data arrays
C EPSIL1  energy loss for each state, species, wavelength, lower bound; eV
C EPSIL2  energy loss for each state, species, wavelength, upper bound; eV
C SIGNO   NO photoionization xsect at Ly-alpha
C
C Array dimensions:
C JMAX    number of altitude levels for GLOW's internal calculations
C NBINS   number of energetic electron energy bins
C LMAX    number of wavelength intervals for solar flux
C NMAJ    number of major species
C NEX     number of ionized/excited species
C NW      number of airglow emission wavelengths
C NC      number of component production terms for each emission
C NST     number of states produced by photoionization/dissociation
C NEI     number of states produced by electron impact
C NF      number of available types of auroral fluxes
C
C
      SUBROUTINE EPHOTO
C
      use ModCommonVariables

      DIMENSION SPECT(NBINS), FLUX(LMAX),
     >          SIGION(NMAJ,LMAX), SIGABS(NMAJ,LMAX),
     >          TPOT(NST,NMAJ), PROB(NST,NMAJ,LMAX),
     >          EPSIL1(NST,NMAJ,LMAX), EPSIL2(NST,NMAJ,LMAX),
     >          SIGAO(LMAX), SIGAO2(LMAX), SIGAN2(LMAX),
     >          SIGIO(LMAX), SIGIO2(LMAX), SIGIN2(LMAX),
     >          PROBO(NST,LMAX), PROBO2(NST,LMAX), PROBN2(NST,LMAX),
     >          PHOTOTE(NMAJ,JMAX),PHOTOTR(JMAX),TAU1(JMAX),
     >          RAT(JMAX), PHOTOT(JMAX)

      DATA (RAT(L),L=1,JMAX)/JMAX*0./

C      SAVE IFIRST
      DATA SIGNO/2.0 E-18/, IFIRST/1/
      ifirst =1
CCC   NNN/5,4,6/
C
      DATA TPOT/13.60, 16.90, 18.60, 28.50, 40.00,  0.00,
     >          12.10, 16.10, 18.20, 20.00,  0.00,  0.00,
     >          15.60, 16.70, 18.80, 30.00, 34.80, 25.00/
C NB - absorption and ionization cross sections are multiplied by 1.E-18
C on first call.
      DATA SIGAO /  18 * 0.00,
     >             0.00, 0.00, 2.12, 4.18, 4.38, 4.23,
     >             4.28, 4.18, 4.18, 8.00,11.35,10.04,
     >            12.21,12.22,12.23,11.90,12.17,12.13,
     >            11.91,11.64,11.25,11.21, 9.64, 9.95,
     >             8.67, 7.70, 7.68, 6.61, 7.13, 6.05,
     >             5.30, 2.90, 1.60, 0.59, 0.16, 0.05,
     >             0.51, 0.07, .012, .002, .0002/
C
      DATA SIGAO2/ 0.50, 1.50, 3.40, 6.00,10.00,13.00,
     >            15.00,12.00, 2.20, 0.30, 3.00, 0.01,
     >             0.30, 0.10, 1.00, 1.10, 1.00, 1.60,
     >            16.53, 4.00,15.54, 9.85,20.87,27.09,
     >            26.66,25.18,21.96,29.05,25.00,26.27,
     >            26.02,25.80,26.10,25.04,22.00,25.59,
     >            24.06,21.59,20.40,19.39,18.17,18.40,
     >            17.19,16.80,16.80,15.10,15.70,13.20,
     >            10.60, 7.10, 4.00, 1.18, 0.32, 0.10,
     >             1.02, 0.14, .024, .004, .0004/
C
      DATA SIGAN2/  18 * 0.00,
     >            36.16, 0.70,16.99,46.63,15.05,30.71,
     >            19.26,26.88,35.46,30.94,26.30,29.75,
     >            23.22,23.20,23.10,22.38,23.20,24.69,
     >            24.53,21.85,21.80,21.07,17.51,18.00,
     >            13.00,11.60,11.60,10.30,10.60, 9.70,
     >             8.00, 4.40, 1.90, 0.60, 0.24, 1.16,
     >             0.48, 0.09, .015, .003, .0003/
C
      DATA SIGIO /  18 * 0.00,
     >             0.00, 0.00, 2.12, 4.18, 4.38, 4.23,
     >             4.28, 4.18, 4.18, 8.00,11.35,10.04,
     >            12.21,12.22,12.23,11.90,12.17,12.13,
     >            11.91,11.64,11.25,11.21, 9.64, 9.95,
     >             8.67, 7.70, 7.68, 6.61, 7.13, 6.05,
     >             5.30, 2.90, 1.60, 0.59, 0.16, 0.05,
     >             0.51, 0.07, .012, .002, .0002/
C
      DATA SIGIO2/ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     >             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     >             0.00, 0.00, 0.00, 0.27, 0.00, 1.00,
     >            12.22, 2.50, 9.34, 4.69, 6.12, 9.39,
     >            11.05, 9.69, 8.59,23.81,23.00,22.05,
     >            25.94,25.80,26.10,25.04,22.00,25.59,
     >            24.06,21.59,20.40,19.39,18.17,18.40,
     >            17.19,16.80,16.80,15.10,15.70,13.20,
     >            10.60, 7.10, 4.00, 1.18, 0.32, 0.10,
     >             1.02, 0.14, .024, .004, .0004/
C
      DATA SIGIN2/  18 * 0.00,
     >             0.00, 0.00, 0.00, 0.00, 0.00,16.75,
     >            10.18,18.39,23.77,23.20,23.00,25.06,
     >            23.22,23.20,23.10,22.38,23.20,24.69,
     >            24.53,21.85,21.80,21.07,17.51,18.00,
     >            13.00,11.60,11.60,10.30,10.60, 9.70,
     >             8.00, 4.40, 1.90, 0.60, 0.24, 1.16,
     >             0.48, 0.09, .015, .003, .0003/
C
      DATA ((PROBO(K,L),K=1,6),L=1,38)
     >           / 120 * 0.00,
     >             1.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     >             1.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     >             1.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     >             1.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     >             1.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     >             1.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     >             1.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     >             0.52, 0.48, 0.00, 0.00, 0.00, 0.00,
     >             0.43, 0.57, 0.00, 0.00, 0.00, 0.00,
     >             0.40, 0.55, 0.05, 0.00, 0.00, 0.00,
     >             0.30, 0.45, 0.25, 0.00, 0.00, 0.00,
     >             0.41, 0.45, 0.24, 0.00, 0.00, 0.00,
     >             0.30, 0.45, 0.25, 0.00, 0.00, 0.00,
     >             0.29, 0.45, 0.26, 0.00, 0.00, 0.00,
     >             0.29, 0.45, 0.25, 0.00, 0.00, 0.00,
     >             0.29, 0.45, 0.26, 0.00, 0.00, 0.00,
     >             0.28, 0.45, 0.26, 0.00, 0.00, 0.00,
     >             0.28, 0.45, 0.27, 0.00, 0.00, 0.00/
      DATA ((PROBO(K,L),K=1,6),L=39,52)
     >           / 0.28, 0.45, 0.27, 0.00, 0.00, 0.00,
     >             0.27, 0.42, 0.26, 0.05, 0.00, 0.00,
     >             0.26, 0.40, 0.25, 0.08, 0.00, 0.00,
     >             0.26, 0.40, 0.25, 0.08, 0.00, 0.00,
     >             0.26, 0.40, 0.25, 0.09, 0.00, 0.00,
     >             0.25, 0.37, 0.24, 0.09, 0.04, 0.00,
     >             0.25, 0.37, 0.24, 0.09, 0.04, 0.00,
     >             0.25, 0.36, 0.23, 0.10, 0.06, 0.00,
     >             0.25, 0.37, 0.23, 0.10, 0.05, 0.00,
     >             0.25, 0.36, 0.23, 0.10, 0.06, 0.00,
     >             0.30, 0.31, 0.20, 0.11, 0.07, 0.00,
     >             0.37, 0.26, 0.17, 0.13, 0.06, 0.00,
     >             0.29, 0.32, 0.21, 0.10, 0.08, 0.00,
     >             0.30, 0.32, 0.21, 0.09, 0.08, 0.00/
      DATA ((PROBO(K,L),K=1,6),L=53,59)
     >           / 0.30, 0.32, 0.21, 0.09, 0.08, 0.00,
     >             0.30, 0.32, 0.21, 0.09, 0.08, 0.00,
     >             0.30, 0.32, 0.21, 0.09, 0.08, 0.00,
     >             0.30, 0.32, 0.21, 0.09, 0.08, 0.00,
     >             0.30, 0.32, 0.21, 0.09, 0.08, 0.00,
     >             0.30, 0.32, 0.21, 0.09, 0.08, 0.00,
     >             0.30, 0.32, 0.21, 0.09, 0.08, 0.00/
C
      DATA ((PROBO2(K,L),K=1,6),L=1,33)
     >           /  90 * 0.00,
     >             1.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     >             1.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     >             1.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     >             1.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     >             1.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     >             1.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     >             1.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     >             1.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     >             1.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     >             1.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     >             0.96, 0.04, 0.00, 0.00, 0.00, 0.00,
     >             0.90, 0.10, 0.00, 0.00, 0.00, 0.00,
     >             0.56, 0.44, 0.00, 0.00, 0.00, 0.00,
     >             0.56, 0.40, 0.00, 0.00, 0.00, 0.00,
     >             0.40, 0.46, 0.14, 0.00, 0.00, 0.00,
     >             0.24, 0.36, 0.35, 0.05, 0.00, 0.00,
     >             0.24, 0.36, 0.35, 0.05, 0.00, 0.00,
     >             0.23, 0.38, 0.31, 0.07, 0.00, 0.00/
      DATA ((PROBO2(K,L),K=1,6),L=34,52)
     >           / 0.35, 0.28, 0.21, 0.16, 0.00, 0.00,
     >             0.30, 0.33, 0.21, 0.16, 0.00, 0.00,
     >             0.36, 0.24, 0.22, 0.18, 0.00, 0.00,
     >             0.36, 0.28, 0.13, 0.23, 0.00, 0.00,
     >             0.42, 0.25, 0.12, 0.21, 0.00, 0.00,
     >             0.42, 0.25, 0.12, 0.21, 0.00, 0.00,
     >             0.42, 0.24, 0.12, 0.22, 0.00, 0.00,
     >             0.40, 0.22, 0.12, 0.26, 0.00, 0.00,
     >             0.37, 0.21, 0.12, 0.30, 0.00, 0.00,
     >             0.36, 0.20, 0.12, 0.32, 0.00, 0.00,
     >             0.35, 0.19, 0.11, 0.35, 0.00, 0.00,
     >             0.35, 0.19, 0.11, 0.35, 0.00, 0.00,
     >             0.34, 0.18, 0.11, 0.37, 0.00, 0.00,
     >             0.34, 0.18, 0.11, 0.37, 0.00, 0.00,
     >             0.34, 0.18, 0.11, 0.37, 0.00, 0.00,
     >             0.33, 0.18, 0.10, 0.39, 0.00, 0.00,
     >             0.30, 0.16, 0.09, 0.45, 0.00, 0.00,
     >             0.20, 0.11, 0.07, 0.62, 0.00, 0.00,
     >             0.10, 0.06, 0.04, 0.80, 0.00, 0.00/
      DATA ((PROBO2(K,L),K=1,6),L=53,59)
     >           / 0.10, 0.06, 0.04, 0.80, 0.00, 0.00,
     >             0.10, 0.06, 0.04, 0.80, 0.00, 0.00,
     >             0.10, 0.06, 0.04, 0.80, 0.00, 0.00,
     >             0.10, 0.06, 0.04, 0.80, 0.00, 0.00,
     >             0.10, 0.06, 0.04, 0.80, 0.00, 0.00,
     >             0.10, 0.06, 0.04, 0.80, 0.00, 0.00,
     >             0.10, 0.06, 0.04, 0.80, 0.00, 0.00/
C
      DATA ((PROBN2(K,L),K=1,6),L=1,38)
     >           / 120 * 0.00,
     >             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     >             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     >             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     >             1.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     >             1.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     >             1.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     >             1.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     >             0.53, 0.47, 0.00, 0.00, 0.00, 0.00,
     >             0.64, 0.36, 0.00, 0.00, 0.00, 0.00,
     >             0.34, 0.66, 0.00, 0.00, 0.00, 0.00,
     >             0.31, 0.59, 0.10, 0.00, 0.00, 0.00,
     >             0.31, 0.59, 0.10, 0.00, 0.00, 0.00,
     >             0.31, 0.59, 0.10, 0.00, 0.00, 0.00,
     >             0.33, 0.57, 0.10, 0.00, 0.00, 0.00,
     >             0.32, 0.58, 0.10, 0.00, 0.00, 0.00,
     >             0.35, 0.55, 0.10, 0.00, 0.00, 0.00,
     >             0.38, 0.53, 0.09, 0.00, 0.00, 0.00,
     >             0.40, 0.47, 0.09, 0.00, 0.00, 0.04/
      DATA ((PROBN2(K,L),K=1,6),L=39,52)
     >           / 0.42, 0.46, 0.08, 0.00, 0.00, 0.04,
     >             0.44, 0.44, 0.08, 0.00, 0.00, 0.04,
     >             0.35, 0.46, 0.10, 0.03, 0.00, 0.06,
     >             0.33, 0.46, 0.10, 0.03, 0.00, 0.07,
     >             0.25, 0.43, 0.10, 0.05, 0.01, 0.16,
     >             0.22, 0.38, 0.09, 0.06, 0.05, 0.20,
     >             0.22, 0.38, 0.09, 0.06, 0.05, 0.20,
     >             0.20, 0.33, 0.07, 0.03, 0.10, 0.26,
     >             0.20, 0.36, 0.07, 0.04, 0.09, 0.24,
     >             0.19, 0.27, 0.07, 0.04, 0.12, 0.31,
     >             0.18, 0.21, 0.07, 0.04, 0.15, 0.35,
     >             0.17, 0.18, 0.07, 0.04, 0.17, 0.36,
     >             0.17, 0.18, 0.07, 0.04, 0.17, 0.36,
     >             0.17, 0.18, 0.07, 0.04, 0.17, 0.36/
      DATA ((PROBN2(K,L),K=1,6),L=53,59)
     >           / 0.17, 0.18, 0.07, 0.04, 0.17, 0.36,
     >             0.17, 0.18, 0.07, 0.04, 0.17, 0.36,
     >             0.17, 0.18, 0.07, 0.04, 0.17, 0.36,
     >             0.17, 0.18, 0.07, 0.04, 0.17, 0.36,
     >             0.17, 0.18, 0.07, 0.04, 0.17, 0.36,
     >             0.17, 0.18, 0.07, 0.04, 0.17, 0.36,
     >             0.17, 0.18, 0.07, 0.04, 0.17, 0.36/
C
C
C First time only:  pack photoabsorption, photoioniation cross sections and
C convert to cm2; pack branching ratios; calculate energy losses:
C
      NNN(1)=5
      NNN(2)=4
      NNN(3)=6


      IF (IFIRST .EQ. 1) THEN
      IFIRST = 0
        DO 10 L=1,LMAX
        SIGABS(1,L) = SIGAO(L)  * 1.E-18
        SIGABS(2,L) = SIGAO2(L) * 1.E-18
        SIGABS(3,L) = SIGAN2(L) * 1.E-18
        SIGION(1,L) = SIGIO(L)  * 1.E-18
        SIGION(2,L) = SIGIO2(L) * 1.E-18
        SIGION(3,L) = SIGIN2(L) * 1.E-18
  
   10   CONTINUE
C
        DO 20 L=1,LMAX
        DO 20 K=1,NST
        PROB(K,1,L) = PROBO(K,L)
        PROB(K,2,L) = PROBO2(K,L)
        PROB(K,3,L) = PROBN2(K,L)
   20   CONTINUE
C
        DO 40 L=1,LMAX 
        DO 40 I=1,NMAJ 
        DO 40 K=1,NNN(I) 
        EPSIL1(K,I,L)=12397.7/WAVE1(L)-TPOT(K,I) 
        EPSIL2(K,I,L)=12397.7/WAVE2(L)-TPOT(K,I) 
   40   CONTINUE 
C
      ENDIF
C
C
C Main loop over altitude:
C
      DO 700 J=1,JMAX
C
      DO 53 I=1,NMAJ
      PHONO(I,J) = 0.
      DO 53 K=1,NST
      PHOTOI(K,I,J) = 0.
      PHOTOD(K,I,J) = 0.
   53 CONTINUE
      DO 55 M=1,NBINS
      SPECT(M) = 0.0 
   55 CONTINUE
C
C loop over wavelengths:
C
      DO 400 L=1,LMAX 
C
C
C Calculate attenuated solar flux:
C
      TAU=0. 
      DO 57 I=1,NMAJ 
         TAU=TAU+SIGABS(I,L)*ZCOL(I,J) 
C      TAU=TAU+SIGABS(I,L)*ZCOL(I,J)*1.E5 
C  NOW CALCULATE FRACTION OF IONIZATION WHICH ARISES DUE TO PHOTOELECTRONS
C  THIS SCHEME DESCRIBED BY RICHARDS AND TORR, JGR, MAY 1988, pp 4063-4064
         IF(I.EQ.1 .AND. WAVE1(L) .LE. 250. .AND. WAVE2(L) .GE. 200.
     >        .AND. ZZ(J) .GT. 0.90E7 .AND. TAU .LT. 9.0)THEN
            RATD=EXP(-TAU)+2*(EXP(-1.3*TAU)+
     >           EXP(-2.0*TAU)+EXP(-2.5*TAU))
            RAT(J)=2.4*(EXP(-TAU))/RATD
         ENDIF
C     PRINT*,TAU
 57   CONTINUE
      IF (TAU .LT. 20) THEN
         FLUX(L)=SFLUX(L)*EXP(-TAU) 
      ELSE
         FLUX(L) = 0.0
      ENDIF
C
C
C Calculate SRC photodissociation of O2, photodissociation of N2, and
C photoionization of NO by solar Ly-alpha:
C
      IF (WAVE1(L) .LT. 1751. AND. WAVE2(L) .GT. 1349.)
     >   PHOTOD(1,2,J) = PHOTOD(1,2,J) + ZMAJ(2,J)*SIGABS(2,L)*FLUX(L)
      PHOTOD(1,3,J) = PHOTOD(1,3,J) +
     >                ZMAJ(3,J)*(SIGABS(3,L)-SIGION(3,L))*FLUX(L)
      IF (WAVE1(L) .LT. 1216. .AND. WAVE2(L) .GT. 1215.)
     >   PHONO(1,J) = PHONO(1,J) + ZNO(J)*SIGNO*FLUX(L)
C
C Loop over species:
C
      DO 300 I=1,NMAJ 
C
C
C Loop over states:
C
      DO 200 K=1,NNN(I) 
C
      E1= EPSIL1(K,I,L) 
      E2= EPSIL2(K,I,L) 
      IF(E2.LT.0.) GO TO 200 
      IF(E1.LT.0.) THEN
        E1=0. 
        DSPECT=(ZMAJ(I,J)*SIGION(I,L)*FLUX(L)*PROB(K,I,L))* 
     >         (E2/(EPSIL2(K,I,L)-EPSIL1(K,I,L))) 
      ELSE
        DSPECT= ZMAJ(I,J)*SIGION(I,L)*FLUX(L)*PROB(K,I,L) 

      ENDIF

      Y=E2-E1 
C SUMMATION OVER ENTIRE WAVELENGTH INTERVAL 
      PHOTOI(K,I,J) = PHOTOI(K,I,J) + DSPECT

C
  200 CONTINUE 
C
  300 CONTINUE
C
  400 CONTINUE 
C
C Put photoelectron spectrum into return array:
C
      DO 500 M=1,NBINS 
      PESPEC(M,J)=SPECT(M)
  500 CONTINUE
C
  700 CONTINUE
C
C    SUM UP PHOTOIONIZATION RATE FOR OXYGEN; ALL STATES FOR EACH ALTITUDE
      DO 402 J = 1,JMAX
C    PHOTOT IS PHOTOIONIZATION RATE (CM-3 SEC-1)
      PHOTOT(J) = 0.0
      DO 401 K  =  1,NNN(1)
      PHOTOT(J) = PHOTOT(J) + PHOTOI(K,1,J)

  401 CONTINUE 
      
C   SAVE THE PRIMARY PHOTOIONIZATION RATE 'PHOTOT'
          PHOTOTR(J) = PHOTOT(J)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C VERTICAL TRANSPORT (ESCAPE) OF PHOTOELECTRONS ABOVE 250 KM REDUCES
C IMPACT IONIZATION EFFECTIVENESS AND FURTHER LOWERS THE RATIO OF P.E.
C IONIZATION TO EUV IONIZATION.   THEREFORE, APPLY LINEAR CORRECTION
C TO MATCH RICHARDS AND TORR PROFILE, JGR,MAY,88,PG 4063 FIG 5.
C ALTITUDE OF INTERPOLATION 250 TO 750 KM
C
      IF (ZZ(J) .GT. 0.250E8 .AND. ZZ(J) .LT. 0.760E8)THEN
          RAT(J)=RAT(J)*(ZZ(17)/ZZ(J))
         ENDIF
         IF (ZZ(J) .GT. 0.750E8) RAT(J)= RAT(67) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   FOLD IN THE CONTRIBUTION FROM PHOTOELECTRONS
C         write(44,*) 'ok'
C         write(44,*) 'PHOTOTR(J)+PHOTOTR(J)*RAT(J)',
C     ;    PHOTOTR(J),PHOTOTR(J),RAT(J)
         
         PHOTOTR(J)= PHOTOTR(J)+PHOTOTR(J)*RAT(J)
C   COMPUTE THE PRIMARY PHOTOIONIZATION FREQUENCY (SEC-1) FOR OPTICAL DEPTH CALC
          PHOTOTF(J) =  PHOTOT(J)/ZMAJ(1,J)
         
       ZZ(J) = ZZ(J)/1.E5
  402 CONTINUE
C
       WRITE(iUnitOutput,9991)
 9991 FORMAT(5X,'      ALT   ','      TAU  ','     PHOT RATE',
     $ '      PHOTOE ','     PHOTO FREQ  ')  
C
C CALCULATE TAU1 AS A FUNCTION OF ALTITIUDE 
      DO 403 J = 1,JMAX
         IF (PHOTOTF(J).GE.1.E-20.AND.PHOTOTF(JMAX).GE.1.E-20)THEN
            TAU1(J) =ALOG((PHOTOTF(JMAX)/PHOTOTF(J)))
         ELSE
            TAU1(J)= 99.999
         ENDIF
C  NOW  COMPUTE THE TOTAL PHOTOIONIZATION FREQUENCY (SEC-1)
         PHOTOTF(J) =  PHOTOTR(J)/ZMAJ(1,J)
         IF(PHOTOTF(J) .LT. 1.E-12) PHOTOTF(J)=9.99E-13   
C     
 403  CONTINUE
C       J=0
C REDUCE TOTAL PHOTIONIZATION ARRAY AND OTHER ARRAYS TO "PW SIZE"
       DO 399 K=1,JMAX
          IF (ZZ(K) .LE. (ALTMIN/1.E5))THEN
            J = 1 
            PHOTOTF(J)=PHOTOTF(K)
            ZZ(J)=ZZ(K)
            TAU1(J)=TAU1(K)
            RAT(J)=RAT(K)
            PHOTOTR(J)=PHOTOTR(K)
            KKK =INT( K + IFACTOR/1E5)
           ENDIF
 399   CONTINUE
C
             J = 1
            DO 499 L= KKK,(JMAX-1),IFACTOR/1E5
            PHOTOTF(J+1)=PHOTOTF(L)
            ZZ(J+1)=ZZ(L)
            TAU1(J+1)=TAU1(L)
            RAT(J+1)=RAT(L)
            PHOTOTR(J+1)=PHOTOTR(L)
            ISUM=J
            J= J + 1
 499  CONTINUE
C NOW FILL PHOTOIONIZATION ARRAY TO MATCH PW CODE SIZE
       DO 404 J=ISUM,NDIM1
       PHOTOTF(J)=PHOTOTF(ISUM+1)
  404  CONTINUE
C
       DO 405 J = 1,40
C PRINT LOW ALT PHOTONIZATION VALUES
      WRITE(iUnitOutput,9994)ZZ(J),TAU1(J),PHOTOTR(J),RAT(J),PHOTOTF(J)
 9994 FORMAT(6X,F9.1,5X,F7.3,5X,F9.3,6X,F7.5,1X,1PE15.3)
  405  CONTINUE
C
      RETURN
      END
