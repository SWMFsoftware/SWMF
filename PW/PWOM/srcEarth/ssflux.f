C
C    Subroutine SSFLUX scales the solar flux according to the 10.7 cm flux
C F107 and the 81-day centered average 10.7 cm flux F107A.  The longwave
C boundary WAVE1 and shortwave boundary WAVE2 of the wavelenth bins are
C returned (Angstroms), and the solar flux in photons cm-2 s-1 returned
C in SFLUX.
C   If ISCALE=0 the flux is scaled using parameterization methods based on
C F107 and F107A.  For ionizing EUV, Hinteregger's contrast ratio method
C is used, based on the Torr and Torr (JGR 90, 6675, 1985) bin structure
C for reference spectrum SC#21REFW.  The 1026A (H LyB) and 335A (FeXVI)
C enhancement ratios are calculated from Hinteregger's formula, using the
C coefficients which reduce to the reference values at F107=71.5,
C F107A=75.4.  The 'best fit' coefficients are not used as they produce some
C negative values at low solar activity, but remain in a 'commented out'
C data statement for reference.  The rest of the spectrum is then calculated
C from these key emissions using Hinteregger's method.  Scaling factors
C were calculated by the author from contrast ratios in the original
C spectrum data file.  For FUV in the 1050A-1350A region, 50A interval
C averaging was done by the author from SC#21REFW, and scaling factors
C also calculated.  This is a mere place holder since 50A bins are not
C adequate for actual band calculations in this region.  For Lyman alpha,
C which is treated seperately as an individual line, Rottman's
C paraterization cited by Bossy (PSS 31, 977, 1983) is used.  For the
C SR continuum, the Torr et al. (JGR 80, 6063, 1980) parameterization
C is used with the coefficients adjusted to reflect the measurements of Rottman
C (JGR 86, 6697, 1981), Mount et al. (JGR 85, 4271, 1980), and Mount and
C Rottman (JGR 86, 9193, 1981; JGR 88, 5403, 1983; JGR 90, 13031, 1985).
C (N.b. - the evidence from SME indicates that the change with solar activity
C of the SR continuum flux is much smaller than what is used here - personal
C communication from Gary Rottman.) 
C   If ISCALE=1 linear interpolation between high and low activity spectra
C is used, based on F107 alone, and assuming that the low activity spectrum
C corresponds to F107=68 and the high activitiy spectrum to F107=243.
C The Hinteregger SC#21REFW and F79050 spectra as binned by Torr and Torr
C are used for ionizing EUV.  For the 1050A-1350A region, the SC#21REFW
C spectrum  averaged into 50A intervals is used for low solar activity;
C the high activity spectrum was obtained by scaling this spectrum
C using the contrast ratios.  For Lyman alpha and the SR continuum,
C linear interpolation amounts to the same thing as the aformentioned
C parameterization.
C    In either case, the EUV fluxes between 250A and 50A are normalized
C upwards (Richards and Torr, 1984).  The normalization coefficient is 2
C at F107=68 and reduces linearly to 1 at F107=243.
C    X-ray fluxes shortwards of 50A are included, 10/88.  The Hinteregger
C fluxes were used to 18A; shortwards of there approximations taken from 495
C notes are employed on a temporary basis.
C
C Definitions:
C ISCALE   =0 for contrast ratio method, =1 for linear interpolation
C F107     daily 10.7 cm flux
C F107A    81-day centered average 10.7 cm flux
C WAVE1    longwave bound of spectral intervals
C WAVE2    shortwave bound of spectral intervals (= WAVE1 for indiv. lines)
C SFLUX    scaled solar flux returned by subroutine
C LMAX     dimension of WAVE1, WAVE2, and SFLUX arrays, must be <= LM
C WAVEL    = WAVE1
C WAVES    = WAVE2
C RFLUX    low solar activity reference flux
C XFLUX    high solar activity flux
C SCALE1   scaling factors for H LyB-keyed chromospheric emissions
C SCALE2   scaling factors for FeXVI-keyed coronal emissions
C LM       dimension of above arrays, currently = 59
C SRA      'A' value for S-R continuum scaling formula
C SRB      'B' value for S-R continuum scaling formula
C B1       fit coefficients for H LyB
C B2       fit coefficients for FeXVI
C R1       enhancement ratio for H LyB
C R2       enhancement ratio for FeXVI
C SFNORM   normalization factor for scaling flux shortwards of 250A
C
C
      SUBROUTINE SSFLUX(ISCALE, F107, F107A, WAVE1, WAVE2, SFLUX, LMAX)
C
      PARAMETER (LM=59)
C
      DIMENSION WAVE1(LMAX), WAVE2(LMAX), SFLUX(LMAX),
     >          WAVEL(LM), WAVES(LM), RFLUX(LM), XFLUX(LM),
     >          SCALE1(LM), SCALE2(LM), SRA(8), SRB(8), B1(3), B2(3)
C
C new B's:
      DATA B1/1.0, 0.0138, 0.005/, B2/1.0, 0.59425, 0.3811/
C
C old B's, commented out:
C     DATA B1/1.31, 0.01106, 0.00492/, B2/-6.618, 0.66159, 0.38319/
C
      DATA WAVEL/ 1750.00, 1700.00, 1650.00, 1600.00, 1550.00, 1500.00,
     >            1450.00, 1400.00, 1350.00, 1300.00, 1250.00, 1215.67,
     >            1200.00, 1150.00, 1100.00, 1050.00, 1031.91, 1025.72,
     >            1000.00,  977.02,  950.00,  900.00,  850.00,  800.00,
     >             789.36,  770.41,  765.15,  750.00,  703.31,  700.00,
     >             650.00,  629.73,  609.76,  600.00,  584.33,  554.37,
     >             550.00,  500.00,  465.22,  450.00,  400.00,  368.07,
     >             350.00,  303.78,  303.31,  300.00,  284.15,  256.30,
     >             250.00,  200.00,  150.00,  100.00,   50.00,   32.00,
     >              23.00,   16.00,    8.00,    4.00,    2.00/
      DATA WAVES/ 1700.00, 1650.00, 1600.00, 1550.00, 1500.00, 1450.00,
     >            1400.00, 1350.00, 1300.00, 1250.00, 1200.00, 1215.67,
     >            1150.00, 1100.00, 1050.00, 1000.00, 1031.91, 1025.72,
     >             950.00,  977.02,  900.00,  850.00,  800.00,  750.00,
     >             789.36,  770.41,  765.15,  700.00,  703.31,  650.00,
     >             600.00,  629.73,  609.76,  550.00,  584.33,  554.37,
     >             500.00,  450.00,  465.22,  400.00,  350.00,  368.07,
     >             300.00,  303.78,  303.31,  250.00,  284.15,  256.30,
     >             200.00,  150.00,  100.00,   50.00,   32.00,   23.00,
     >              16.00,    8.00,    4.00,    2.00,    1.00/
      DATA RFLUX/  370.45,  203.69,   96.00,   69.71,   50.70,   26.67,
     >              17.21,    8.26,   12.86,    4.10,    5.20,  333.80,
     >               2.78,    0.70,    3.07,    3.64,    3.18,    4.38,
     >               1.78,    5.96,    4.22,    4.43,    1.93,    0.87,
     >               0.79,    0.24,    0.20,    0.17,    0.39,    0.22,
     >               0.17,    1.50,    0.45,    0.48,    1.58,    0.80,
     >               0.51,    0.31,    0.18,    0.39,    0.21,    0.74,
     >               0.87,    6.00,    0.24,    0.84,    0.10,    0.27,
     >               0.92,    1.84,    0.13,    0.38,  0.0215,  0.0067,
     >             0.0009,  0.0003,   1.E-6,   3.E-9,   1.E-11/
      DATA XFLUX/  464.20,  241.50,  131.50,  101.90,   81.32,   48.71,
     >              37.16,   21.14,   30.70,   11.20,   12.00,  438.80,
     >               6.50,    1.60,    6.40,    8.66,    9.04,   13.12,
     >               4.42,   13.18,   12.03,   13.29,    5.01,    2.18,
     >               1.59,    0.67,    0.43,    0.43,    0.72,    0.46,
     >               0.48,    3.02,    1.46,    1.02,    4.86,    1.59,
     >               1.57,    1.67,    0.36,    0.99,    2.20,    1.39,
     >               5.63,   11.28,    2.50,    4.14,    3.16,    0.59,
     >               3.70,    4.85,    0.34,    1.15,    0.18,    0.08,
     >              0.025,    0.03,   1.E-3,   3.E-5,   1.E-6/
      DATA SCALE1/35347.5, 33095.6, 18040.6, 13733.0, 12564.2, 7121.38,
     >            6608.74, 5779.89, 8009.80, 3186.34, 3033.78,  47555.,
     >            1692.09,  405.95, 1516.20, 2731.70, 3314.57, 4375.00,
     >            1316.91, 3621.91, 3908.56, 4432.54, 1541.21,  531.73,
     >             364.83,    0.00,  116.00,  129.41,  162.48,   94.07,
     >              41.29,  709.50,    0.00,  268.47, 1561.05,  367.64,
     >             290.06,  184.36,    0.00,   86.15,    7.50,    0.00,
     >               0.00, 2220.00,    0.00,   61.00,    0.00,   86.95,
     >             206.00,  135.89,   60.35,  157.12,    7.06,    0.75,
     >               0.00,    0.00,    0.00,    0.00,    0.00/
      DATA SCALE2/   0.00,    0.00,    0.00,    0.00,    0.00,    0.00,
     >               0.00,    0.00,    0.00,    0.00,    0.00,    0.00,
     >               0.00,    0.00,    0.00,    0.00,    0.00,    0.00,
     >               0.00,    0.00,    0.00,    0.00,    0.00,    0.00,
     >               0.00,    5.34,    0.00,    0.00,    0.00,    0.54,
     >               3.30,    0.00,   12.60,    0.00,    0.00,    0.00,
     >               5.34,   11.63,    2.28,    5.56,   24.93,    8.16,
     >              60.69,    0.00,   28.20,   45.90,   40.80,    1.27,
     >              35.47,   42.80,    1.12,    6.19,    1.26,    0.69,
     >               0.23,    0.30,    0.01,   3.E-4,   1.E-5/
      DATA SRA/     0.536,   0.216,   0.203,   0.184,   0.175,   0.126,
     >              0.114,   0.073/
      DATA SRB/     334.0,   189.0,    82.2,    57.2,    38.8,    18.1,
     >               9.46,    3.30/
C
C
      IF (ISCALE .EQ. 0) THEN
        R1 =  B1(1) + B1(2)*(F107-71.5) + B1(3)*(F107-F107A+3.9)
        R2 =  B2(1) + B2(2)*(F107-71.5) + B2(3)*(F107-F107A+3.9)
        DO 100 L=1,LMAX
        IF (L .LT. 9) THEN
          SFLUX(L) = SRA(L) * F107 + SRB(L)
        ELSE
          IF (L .EQ. 12) THEN
            SFLUX(L) = 332. + 0.6 * (F107-65.)
          ELSE
            SFLUX(L) = (RFLUX(L) + ((R1-1.)*SCALE1(L)
     >                            + (R2-1.)*SCALE2(L)) / 1000.)
            IF (SFLUX(L) .LT. 0.01) SFLUX(L) = 0.01
          ENDIF
        ENDIF
  100   CONTINUE
      ELSE
        FRAT = (F107-68.) / (243.-68.)
        DO 200 L=1,LMAX
        SFLUX(L) = RFLUX(L) + (XFLUX(L)-RFLUX(L)) * FRAT
  200   CONTINUE
      ENDIF
C
      SFNORM = 2. - (F107-68.) / (243.-68.)
      IF (SFNORM .LT. 1.0) SFNORM = 1.0
C
      DO 300 L=1,LMAX
      WAVE1(L) = WAVEL(L)
      WAVE2(L) = WAVES(L)
      IF (WAVE1(L) .LT. 251. .AND. WAVE2(L) .GT. 0.1)
     >   SFLUX(L) = SFLUX(L) * SFNORM
      SFLUX(L) = SFLUX(L) * 1.E9
  300 CONTINUE
      RETURN
      END
