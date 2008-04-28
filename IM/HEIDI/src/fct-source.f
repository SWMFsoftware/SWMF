      SUBROUTINE  ETBFCT(RHOO, RHON, N, RBC, LBC)
C
      REAL   RHOO(*),  RHON(*)
CD
CD    ***************************************************************
CD
CD    ETBFCT  ( RHOO, RHON, N, RBC, LBC )              D.3
CD    ORIGINATOR - J.P. BORIS   CODE  7706,  NRL       OCT. 1975
CD
CD    DESCRIPTION:  THIS  ROUTINE  SOLVES  GENERALISED  CONTINUITY  
CD    EQUATIONS  OF  THE  FORM  DRHO/DT = -DIV (RHO*V) + SOURCES IN
CD    EITHER  CARTESIAN,  CYLINDRICAL, OR  SPHERICAL  GEOMETRY.  THE
CD    FINITE-DIFFERENCE  GRID  CAN  BE  EULERIAN, SLIDING  REZONE
CD    OR  LAGRANGIAN   AND  CAN  BE  ARBITRARILY  SPACED.  THE
CD    ALGORITHM  USED  IS  A  LOW-PHASE  -ERROR  FCT  ALGORITHM,
CD    VECTORISED  AND  OPTIMISED  FOR  SPEED.
CD
CD    ARGUMENTS : IN  THIS  ROUTINE  THE  RIGHT  BOUNDARY  AT  RADR  IS
CD    HALF  A   CELL  BEYOND  THE  LAST  GRID  POINT  N  AT  RADN(N) AND
CD    THE  LEFT  BOUNDARY  AT  RADL  IS  HALF  A  CELL  BEFORE  THE
CD    FIRST  GRID  POINT  AT  RADN(1).
CD    RHON      REAL  ARRAY(N)  GRID  POINT  DENSITIES  AT  END OF  STEP.
CD    RHON      REAL  ARRAY(N)  GRID  POINT  DENSITIES  AT START OF STEP.
CD    N         INTEGER         NUMBER  OF  INTERIOR  GRID  POINTS.
CD    RBC       REAL            RIGHT  BOUNDARY  CONDITION  FACTOR.
CD    LBC       REAL            LEFT   BOUNDARY  CONDITION  FACTOR.
CD    
CD    LANGUAGE  AND  LIMITATIONS: THE  SUBROUTINE  ETBFCT  IS  A
CD    MULTIPLE - ENTRY  FORTRAN  ROUTINE  IN  SINGLE  PRECISION ( 32
CD    BITS  ASC). THE  ASC  PARAMETER  STATEMENT  IS  USED  TO  SET
CD    SYMBOLICALLY  THE  INTERNAL  ARRAY  DIMENSIONS.  UNDERFLOWS  ARE
CD    POSSIBLE  WHEN  THE  FUNCTION  BEING  TRANSFORMED  HAS  MANY  
CD    ZEROES.  THE  CALCULATIONS  GENERALLY  MISCONSERVE  BY  ONE  OR
CD    TWO  BITS  PER  CYCLE.  THE  RELATIVE  PHASE  AND  AMPLITUDE  
CD    ERRORS ( FOR  SMOOTH  FUNCTIONS )  ARE  TYPICALLY  SEVERAL  PERCENT
CD    FOR  CHARACTERISTIC  LENGTHS  OF  1-2  CELLS ( WAVELENGTHS  OF
CD    ORDER  10  CELLS ).  SHOCKS  ARE  GENERALLY  ACCURATE  TO  BETTER
CD    THEN  1  PERCENT.  THIS  SUBROUTINE  MUST  BE  COMPILED  WITH  THE
CD    Y  OPTION  TO  FORCE  STORAGE  AND  RETENTION  OF  INTERNAL  
CD    VARIABLES.  ALTERNATIVELY  A  COMMON  BLOCK  CAN  BE  ADDED  TO
CD    ACCOMPLISH  THE  SAME  END.
CD
CD    ENTRY  POINTS : OGRIDE, NGRIDE, VELOCE, SOURCE, CONSRE, SEE #11
CD    OF  THE  DETAILED  DOCUMENTATION  ( OR  THE  LISTING  BELOW) FOR
CD    THE  EXPLANATIONS  AND  USE  OF  THE  ARGUMENTS  TO  THESE  ENTRIES.
CD
CD    NO  AUXILIARY  OR  LIBRARY  ROUTINES  ARE  CALLED  BY  ETBFCT.
CD
CD
CD    *****************************************************************
C
      LOGICAL       LSOURC, first			! initialize variables (cer)
      REAL          RBC, LBC
      REAL          SOURCE(350), SCRH(350), RHOT(350), DIFF(350)
      REAL          ADUGTH(350), FLXH(350), NULH(350), MULH(350)
      REAL          LNRHOT(350), FSGN(350), EPSH(350)
      REAL          LORHOT(350), ADUDTH(350)
      REAL          LO(350), LN(350), LH(350), RLO(350), RLN(350)
      REAL          RNH(350), ROH(350), RLH(350), AH(350)
      REAL          U(*), RADN(*), C(*), D(*), RHO(*)
      REAL          TERP(350), TERM(350), FABS(350)
      INTEGER       ALPHA, N
      EQUIVALENCE   (EPSH(1), SCRH(1)),     (LNRHOT(1), LORHOT(1))
      EQUIVALENCE   (FLXH(1), SCRH(1)),     (FSGN(1),  RHOT(1))
      EQUIVALENCE   (FABS(1), SCRH(1)),     (TERP(1),  SOURCE(1))
      EQUIVALENCE   (TERM(1), SOURCE(1))
      DATA          PI, FTPI/3.1415927, 4.1887902 /
      DATA          LSOURC/.FALSE./, first/.true./
      save					! needed for MPW FORTRAN (cer)
C
C     CALCULATE  THE  DIFFUSIVE  AND  CONVECTIVE  FLUXES.
C
      NP  =  N  +  1
      DO  11  I  = 2, N
      FLXH(I) = 0.5*ADUDTH(I)*(RHOO(I) + RHOO(I-1))
 11   DIFF(I) = NULH(I)*(RHOO(I) - RHOO(I-1))
      RHOL    = RHOO(1)*LBC
      RHOR    = RHOO(N)*RBC
      DIFF(1) = NULH(1)*(RHOO(1) - RHOL )
      DIFF(NP)= NULH(NP)*(RHOR - RHOO(N) )
      FLXH(1) = 0.5*ADUDTH(1)*(RHOO(1) + RHOL )
      FLXH(NP)= 0.5*ADUDTH(NP)*(RHOR + RHOO(N) )
C
C     CALCULATE  LAMBDAO*RHOT,  THE  TRANSPORTED  MASS  ELEMENTS.
C
      DO  12  I  = 1, N
 12   LORHOT(I) = LO(I)*RHOO(I) - FLXH(I+1) + FLXH(I)
C
C     ADD  IN  THE  SOURCE  TERMS  AS  APPROPRIATE.
C
      IF  (.NOT. LSOURC ) GOTO  14
      DO  13  I = 1, N
 13   LORHOT(I) = LORHOT(I) + SOURCE(I)
C
C     CALCULATE  THE  PHOENICAL  ANTIDIFFUSIVE  FLUXES  HERE.
C
 14   DO  16  I = 1, N
 16   RHOT (I) = LORHOT(I)*RLO(I)
      DO  17  I = 2, N
 17   FLXH(I) = MULH(I)*(RHOT(I) - RHOT(I-1))
      FLXH(1) = MULH(1)*(RHOT(1) - LBC*RHOT(1))
      FLXH(NP)= MULH(NP)*(RBC*RHOT(N) - RHOT(N))
C
C     CALCULATE  THE  TRANSPORTED/ DIFFUSED  DENSITY  AND  GRID  
C     DIFFERENCES.
C
      DO  18  I = 1, N
 18   LNRHOT(I) = LORHOT(I) + DIFF(I+1) - DIFF(I)
      DO  19  I = 1, N
 19   RHOT(I) = LNRHOT(I)*RLN(I)
      DO  20  I = 2, N
 20   DIFF(I) = RHOT(I) - RHOT(I-1)
      DIFF(1) = RHOT(1) - LBC*RHOT(1)
      DIFF(NP)= RBC*RHOT(N) - RHOT(N)
C
C     CALCULATE  THE  SIGN  AND  MAGNITUDE  OF  THE  ANTIDIFFUSIVE  FLUX.
C
      DO  21  I  = 1, NP
      FSGN(I)  =  +1.0
      IF ( DIFF(I).LT.0 ) FSGN(I) = -1.0
      FABS(I)  =  FLXH(I)
      IF ( FLXH(I) .LT. 0 ) FABS(I) = -FLXH(I)
 21   CONTINUE
C
C     CALCULATE THE FLUX-LIMITING CHANGES ON THE RIGHT AND THE LEFT.
C
      DO  23  I = 1, N
 23   TERP(I) = FSGN(I)*LN(I)*DIFF(I+1)
      TERP(NP)= 1.0E+35
      DO 24 I = 1, NP
 24   FABS(I) = AMIN1(TERP(I), FABS(I) )
      DO  25  I = 2, NP
 25   TERM(I) = FSGN(I)*LN(I-1)*DIFF(I-1)
      TERM(1) = 1.0E+35
C
C     CORRECT  THE  FLUXES  COMPLETELY  NOW.
C
      DO  35  I = 1, NP
 35   DIFF(I) = AMIN1(FABS(I), TERM(I))
      DO  36  I = 1, NP
 36   FLXH(I) = AMAX1( 0.0, DIFF(I))
      DO  37  I = 1, NP
 37   FLXH(I) = FSGN(I)*FLXH(I)
C
C     CALCULATE  THE  NEW  FLUX-CORRECTED  DENSITIES.
C
      DO  40  I = 1,N
 40   LNRHOT(I) = LNRHOT(I) - FLXH(I+1) + FLXH(I)
      DO  41  I = 1, N
      SOURCE(I) = 0.0
 41   RHON(I)   =  LNRHOT(I)*RLN(I)
      LSOURC    =  .FALSE.
      GO TO 600
C
C
C    -------------------------------------------------------------------
C
      ENTRY   VELOCE  ( U, N, UR, UL, DT )
C
C
CD    ******************************************************************
CD
CD    VELOCE  (U, N, UR, UL, DT )
CD    DESCRIPTION : THIS  ENTRY  CALCULATES  ALL  VELOCITY-DEPENDANT
CD                  COEFFS.
CD
CD    ARGUMENTS :
CD    U    REAL  ARRAY(N)   FLOW  VELOCITY  AT  THE  GRID  POINTS
CD    N    INTEGER          NUMBER  OF  INTERIOR  GRID  POINTS
CD    UR   REAL             VELOCITY  OF  FLOW  AT  RIGHT  BOUNDARY
CD    UL   REAL             VELOCITY  OF  FLOW  AT  LEFT   BOUNDARY
CD    DT   REAL             STEPSIZE  FOR THE  TIME  INTEGRATION
CD
CD    ******************************************************************
C
C     CALCULATES  THE  INTERFACE  AREA  X  VELOCITY  DIFFERENTIAL  X  DT.
C
      NP = N + 1
      DTH= 0.5*DT
      DO  101  I = 2, N
 101  ADUDTH(I)  = AH(I)*DTH*(U(I) + U(I-1)) - ADUGTH(I)
      ADUDTH(1)  = AH(1)*DT*UL - ADUGTH(1)
      ADUDTH(NP) = AH(NP)*DT*UR- ADUGTH(NP)
C
C     CALCULATE  THE  HALF-CELL  EPSILON  (V*DT/DX)
C
      DO  102  I = 1, NP
 102  EPSH(I) = ADUDTH(I)*RLH(I)
C
C     NEXT  CALCULATE  THE  DIFFUSION  AND  ANTIDIFFUSION  COEFFICIENTS.
C     VARIATION  WITH  EPSILON  MEANS  FOURTH-ORDER  ACCURATE  PHASES.
C
      DO  103  I = 1, NP
      NULH(I) = 0.16666667 + 0.3333333*EPSH(I)*EPSH(I)
 103  MULH(I) = 0.25 - 0.5*NULH(I)
      DO 104 I = 1, NP
      NULH(I) = LH(I)*NULH(I)
 104  MULH(I) = LH(I)*MULH(I)
      RETURN
C
C
C     ------------------------------------------------------------------
C
      ENTRY  NGRIDE  ( RADN, N, RADR, RADL, ALPHA )
C
C
CD    ******************************************************************
CD
CD    NGRIDE  (RADN, N, RADR, RADL, ALPHA)
CD    DESCRIPTION :  THIS  ENTRY  SETS  NEW  GEOMETRY  VARIABLES  AND
CD                   COEFFS.
CD    ARGUMENTS :
CD    RADN      REAL  ARRAY(N)  NEW  GRID  POSITIONS
CD    N         INTEGER         NUMBER  OF  INTERIOR  GRID  POSITIONS
CD    RADR      REAL            POSITION  OF  THE  RIGHT  BOUNDARY
CD    RADL      REAL            POSITION  OF  THE  LEFT   BOUNDARY
CD    ALPHA     INTEGER         = 1 FOR  CARTESIAN  GEOMETRY
CD                              = 2 FOR  CYLINDRICAL GEOMETRY
CD                              = 3 FOR  SPHERICAL   GEOMETRY
CD
CD    *******************************************************************
C
C     CALCULATES  THE  NEW  HALF-CELL  POSITIONS  AND  GRID  CHANGES.
C
      NP = N + 1
      DO  202  I = 2, N
      RNH(I) = 0.5*(RADN(I) + RADN(I-1) )
C     WRITE ( 3,* ) ( RADN(I), RADN(I-1),RNH(I) )
 202  CONTINUE
      RNH(1) = RADL
      RNH(NP)= RADR
C
C     CALCULATE  THE  THREE  COORDINATE  SYSTEMS.
C
      GO TO ( 203, 206, 209 ), ALPHA
C
C     CARTESIAN  COORDINATES.
C
 203  DO  204  I = 1, NP
 204  AH(I) = 1.0
      DO  205  I = 1, N
      LN(I) = RNH(I+1) - RNH(I)
C     WRITE ( 3,* )( RNH(I),RNH(I+1),LN(I) )
 205  CONTINUE
      GOTO  213
C
C     CYLINDRICAL  COORDINATES.
C
 206  DO  207  I = 1, NP
      DIFF(I) = RNH(I)*RNH(I)
 207  AH(I) = PI*(ROH(I) + RNH(I))
      DO 208 I = 1, N
 208  LN(I) = PI*(DIFF(I+1) - DIFF(I) )
      GOTO  213
C
C     SPHERICAL   COORDINATES.
C
 209  DO  210  I = 1, NP
      DIFF(I) = RNH(I)*RNH(I)*RNH(I)
 210  SCRH(I) = (ROH(I) + RNH(I))*ROH(I)
      DO  211  I = 1, NP
 211  AH(I)  =  FTPI*(SCRH(I) + RNH(I)*RNH(I))
      DO  212  I = 1, N
 212  LN(I)  =  FTPI*(DIFF(I+1) - DIFF(I) )
C
C     NOW  THE  GEOMETRIC  VARIABLES  WHICH  ARE  SYSTEM  INDEPENDENT.
C
 213  DO  214  I = 2, N
 214  LH(I) = 0.5*(LN(I) + LN(I-1))
      LH(1) = LN(1)
      LH(NP)= LN(N)
      DO  215  I = 1, N
 215  RLN(I) = 1.0/LN(I)
      DO  216  I = 1, NP
 216  ADUGTH(I) = AH(I)*(RNH(I) - ROH(I))
      DO  217  I = 2, N
 217  RLH(I)  = 0.5*(RLN(I) + RLN(I-1))
      RLH(1)  = RLN(1)
      RLH(NP) = RLN(N)
      RETURN
C
C
C     -----------------------------------------------------------------
C
      ENTRY  SOURCQ ( N, DT, MODES, C, D, DR, DL )
C
C
CD    ******************************************************************
CD
CD    SOURCE ( N, DT, MODES, C, D, DR, DL )
CD    DESCRIPTION: THIS  ENTRY  ACCUMULATES  DIFFERENT  SOURCE  TERMS.
CD
CD    ARGUMENTS:
CD    N   INTEGER        NUMBER  OF  INTERIOR  GRID  POINTS
CD    DT  REAL           STEPSIZE  FOR  THE  TIME  INTEGRATION
CD    MODES  INTEGER     = 1  COMPUTES  +  DIV (D)
CD                       = 2  COMPUTES  +  C*GRAD(D)
CD                       = 3  ADDS  +  D TO THE  SOURCES
CD    C   REAL  ARRAY(N) ARRAY  OF  SOURCE  VARIABLES  AT  GRID  POINTS
CD    D   REAL  ARRAY(N) ARRAY  OF  SOURCE  VARIABLES  AT  GRID  POINTS
CD    DR  REAL           RIGHT  BOUNDARY  VALUE  OF  D
CD    DL  REAL           LEFT   BOUNDARY  VALUE  OF  D
CD
CD    *****************************************************************
C
      if (first) then					! initialize to 0.0 (CER)
          do i = 1, N
	      source(i) = 0.0
	  end do
	  first = .false.
      end if
      NP = N + 1
      DTH = 0.5*DT
      DTQ = 0.25*DT
      GOTO  ( 310, 320, 330 ), MODES
C
C     DIV (D) IS  COMPUTED  CONSERVATIVELY  AND  ADDED  TO  THE  SOURCES.
C
 310  DO  311  I = 2, N
 311  SCRH(I) = DTH*AH(I)*(D(I) + D(I-1))
      SCRH(1) = DT*AH(1)*DL
      SCRH(NP)= DT*AH(NP)*DR
      DO  312  I = 1, N
 312  SOURCE(I) = SOURCE(I) + SCRH(I+1) - SCRH(I)
      LSOURC  = .TRUE.
      RETURN
C
C     C*GRAD(D) IS  COMPUTED  EFFICIENTLY  AND  ADDED  TO  THE  SOURCES.
C
 320  DO  321  I = 2, N
 321  SCRH(I) = DTQ*(D(I) + D(I-1))
      SCRH(1) = DTH*DL
      SCRH(NP)= DTH*DR
      DO  322  I = 1, N
 322  DIFF(I) = SCRH(I+1) - SCRH(I)
      DO  323  I = 1, N
 323  SOURCE(I) = SOURCE(I) + C(I)*(AH(I+1) +AH(I))*DIFF(I)
      LSOURC = .TRUE.
      RETURN
C
C     D  IS  ADDED  TO  THE  SOURCES  IN  AN  EXPLICIT  FORMULATION.
C
 330  DO  331  I = 1, N
 331  SOURCE(I) = SOURCE(I) + DT*LO(I)*D(I)
      LSOURC  =  .TRUE.
      RETURN
C
C
C     -----------------------------------------------------------------
C
      ENTRY   OGRIDE (N)
C
CD    *****************************************************************
CD
CD    OGRIDE   (N)
CD    DESCRIPTION:  THIS  ENTRY  COPIES  OLD  GRID  AND  GEOMETRY  VARIABLES.
CD
CD    ARGUMENTS : 
CD    N   INTEGER     NUMBER  OF  INTERIOR  GRID  POINTS
CD
CD    ******************************************************************
C
C      COPY  THE  PREVIOUSLY  NEW  GRID  VALUES  TO  BE  USED  AS  THE
C            OLD  GRID.
C
       NP = N + 1
       DO  401  I = 1, N
       LO(I) = LN(I)
 401   RLO(I) = RLN(I)
       DO  402  I = 1, NP
 402   ROH(I) = RNH(I)
       RETURN
C
C
C      ------------------------------------------------------------------
C
       ENTRY  CONSRE  ( RHO, N, CSUM )
C
C
CD     *****************************************************************
CD
CD     CONSRE  ( RHO, N, CSUM )
CD     DESCRIPTION : THIS  ENTRY  COMPUTES  THE  OSTENSIBLY  CONSERVED
CD                   SUM.
CD
CD     ARGUMENTS:
CD     RHO     REAL  ARRAY(N)  GRID  PT. VALUES  FOR  CONSERVATION  SUM.
CD     N       INTEGER         NUMBER  OF  INTERIOR  GRID  POINTS.
CD     CSUM    REAL            VALUE   OF  THE  CONSERVATION  SUM OF RHO
CD
CD     *****************************************************************
CD
C      COMPUTE  THE  OSTENSIBLY  CONSERVED  TOTAL  MASS ( BEWARE YOUR B.C. )
C
       CSUM  =  0.0
       DO  501  I = 1, N
 501   CSUM  =  CSUM  +  LN(I)*RHO(I)
 600   RETURN
        END

      SUBROUTINE  ETBFCT1(RHOO, RHON, N, RBC, LBC)
C
      REAL   RHOO(*),  RHON(*)
CD
CD    ***************************************************************
CD
CD    ETBFCT  ( RHOO, RHON, N, RBC, LBC )              D.3
CD    ORIGINATOR - J.P. BORIS   CODE  7706,  NRL       OCT. 1975
CD
CD    DESCRIPTION:  THIS  ROUTINE  SOLVES  GENERALISED  CONTINUITY  
CD    EQUATIONS  OF  THE  FORM  DRHO/DT = -DIV (RHO*V) + SOURCES IN
CD    EITHER  CARTESIAN,  CYLINDRICAL, OR  SPHERICAL  GEOMETRY.  THE
CD    FINITE-DIFFERENCE  GRID  CAN  BE  EULERIAN, SLIDING  REZONE
CD    OR  LAGRANGIAN   AND  CAN  BE  ARBITRARILY  SPACED.  THE
CD    ALGORITHM  USED  IS  A  LOW-PHASE  -ERROR  FCT  ALGORITHM,
CD    VECTORISED  AND  OPTIMISED  FOR  SPEED.
CD
CD    ARGUMENTS : IN  THIS  ROUTINE  THE  RIGHT  BOUNDARY  AT  RADR  IS
CD    HALF  A   CELL  BEYOND  THE  LAST  GRID  POINT  N  AT  RADN(N) AND
CD    THE  LEFT  BOUNDARY  AT  RADL  IS  HALF  A  CELL  BEFORE  THE
CD    FIRST  GRID  POINT  AT  RADN(1).
CD    RHON      REAL  ARRAY(N)  GRID  POINT  DENSITIES  AT  END OF  STEP.
CD    RHON      REAL  ARRAY(N)  GRID  POINT  DENSITIES  AT START OF STEP.
CD    N         INTEGER         NUMBER  OF  INTERIOR  GRID  POINTS.
CD    RBC       REAL            RIGHT  BOUNDARY  CONDITION  FACTOR.
CD    LBC       REAL            LEFT   BOUNDARY  CONDITION  FACTOR.
CD    
CD    LANGUAGE  AND  LIMITATIONS: THE  SUBROUTINE  ETBFCT  IS  A
CD    MULTIPLE - ENTRY  FORTRAN  ROUTINE  IN  SINGLE  PRECISION ( 32
CD    BITS  ASC). THE  ASC  PARAMETER  STATEMENT  IS  USED  TO  SET
CD    SYMBOLICALLY  THE  INTERNAL  ARRAY  DIMENSIONS.  UNDERFLOWS  ARE
CD    POSSIBLE  WHEN  THE  FUNCTION  BEING  TRANSFORMED  HAS  MANY  
CD    ZEROES.  THE  CALCULATIONS  GENERALLY  MISCONSERVE  BY  ONE  OR
CD    TWO  BITS  PER  CYCLE.  THE  RELATIVE  PHASE  AND  AMPLITUDE  
CD    ERRORS ( FOR  SMOOTH  FUNCTIONS )  ARE  TYPICALLY  SEVERAL  PERCENT
CD    FOR  CHARACTERISTIC  LENGTHS  OF  1-2  CELLS ( WAVELENGTHS  OF
CD    ORDER  10  CELLS ).  SHOCKS  ARE  GENERALLY  ACCURATE  TO  BETTER
CD    THEN  1  PERCENT.  THIS  SUBROUTINE  MUST  BE  COMPILED  WITH  THE
CD    Y  OPTION  TO  FORCE  STORAGE  AND  RETENTION  OF  INTERNAL  
CD    VARIABLES.  ALTERNATIVELY  A  COMMON  BLOCK  CAN  BE  ADDED  TO
CD    ACCOMPLISH  THE  SAME  END.
CD
CD    ENTRY  POINTS : OGRIDE1, NGRIDE1, VELOCE1, SOURCE1, CONSRE1, SEE #11
CD    OF  THE  DETAILED  DOCUMENTATION  ( OR  THE  LISTING  BELOW) FOR
CD    THE  EXPLANATIONS  AND  USE  OF  THE  ARGUMENTS  TO  THESE  ENTRIES.
CD
CD    NO  AUXILIARY  OR  LIBRARY  ROUTINES  ARE  CALLED  BY  ETBFCT.
CD
CD
CD    *****************************************************************
C
      LOGICAL       LSOURC, first			! initialize variables (cer)
      REAL          RBC, LBC
      REAL          SOURCE(350), SCRH(350), RHOT(350), DIFF(350)
      REAL          ADUGTH(350), FLXH(350), NULH(350), MULH(350)
      REAL          LNRHOT(350), FSGN(350), EPSH(350)
      REAL          LORHOT(350), ADUDTH(350)
      REAL          LO(350), LN(350), LH(350), RLO(350), RLN(350)
      REAL          RNH(350), ROH(350), RLH(350), AH(350)
      REAL          U(*), RADN(*), C(*), D(*), RHO(*)
      REAL          TERP(350), TERM(350), FABS(350)
      INTEGER       ALPHA, N
      EQUIVALENCE   (EPSH(1), SCRH(1)),     (LNRHOT(1), LORHOT(1))
      EQUIVALENCE   (FLXH(1), SCRH(1)),     (FSGN(1),  RHOT(1))
      EQUIVALENCE   (FABS(1), SCRH(1)),     (TERP(1),  SOURCE(1))
      EQUIVALENCE   (TERM(1), SOURCE(1))
      DATA          PI, FTPI/3.1415927, 4.1887902 /
      DATA          LSOURC/.FALSE./, first/.true./
      save					! needed for MPW FORTRAN (cer)
C
C     CALCULATE  THE  DIFFUSIVE  AND  CONVECTIVE  FLUXES.
C
      NP  =  N  +  1
      DO  11  I  = 2, N
      FLXH(I) = 0.5*ADUDTH(I)*(RHOO(I) + RHOO(I-1))
 11   DIFF(I) = NULH(I)*(RHOO(I) - RHOO(I-1))
      RHOL    = RHOO(1)*LBC
      RHOR    = RHOO(N)*RBC
      DIFF(1) = NULH(1)*(RHOO(1) - RHOL )
      DIFF(NP)= NULH(NP)*(RHOR - RHOO(N) )
      FLXH(1) = 0.5*ADUDTH(1)*(RHOO(1) + RHOL )
      FLXH(NP)= 0.5*ADUDTH(NP)*(RHOR + RHOO(N) )
C
C     CALCULATE  LAMBDAO*RHOT,  THE  TRANSPORTED  MASS  ELEMENTS.
C
      DO  12  I  = 1, N
 12   LORHOT(I) = LO(I)*RHOO(I) - FLXH(I+1) + FLXH(I)
C
C     ADD  IN  THE  SOURCE  TERMS  AS  APPROPRIATE.
C
      IF  (.NOT. LSOURC ) GOTO  14
      DO  13  I = 1, N
 13   LORHOT(I) = LORHOT(I) + SOURCE(I)
C
C     CALCULATE  THE  PHOENICAL  ANTIDIFFUSIVE  FLUXES  HERE.
C
 14   DO  16  I = 1, N
 16   RHOT (I) = LORHOT(I)*RLO(I)
      DO  17  I = 2, N
 17   FLXH(I) = MULH(I)*(RHOT(I) - RHOT(I-1))
      FLXH(1) = MULH(1)*(RHOT(1) - LBC*RHOT(1))
      FLXH(NP)= MULH(NP)*(RBC*RHOT(N) - RHOT(N))
C
C     CALCULATE  THE  TRANSPORTED/ DIFFUSED  DENSITY  AND  GRID  
C     DIFFERENCES.
C
      DO  18  I = 1, N
 18   LNRHOT(I) = LORHOT(I) + DIFF(I+1) - DIFF(I)
      DO  19  I = 1, N
 19   RHOT(I) = LNRHOT(I)*RLN(I)
      DO  20  I = 2, N
 20   DIFF(I) = RHOT(I) - RHOT(I-1)
      DIFF(1) = RHOT(1) - LBC*RHOT(1)
      DIFF(NP)= RBC*RHOT(N) - RHOT(N)
C
C     CALCULATE  THE  SIGN  AND  MAGNITUDE  OF  THE  ANTIDIFFUSIVE  FLUX.
C
      DO  21  I  = 1, NP
      FSGN(I)  =  +1.0
      IF ( DIFF(I).LT.0 ) FSGN(I) = -1.0
      FABS(I)  =  FLXH(I)
      IF ( FLXH(I) .LT. 0 ) FABS(I) = -FLXH(I)
 21   CONTINUE
C
C     CALCULATE THE FLUX-LIMITING CHANGES ON THE RIGHT AND THE LEFT.
C
      DO  23  I = 1, N
 23   TERP(I) = FSGN(I)*LN(I)*DIFF(I+1)
      TERP(NP)= 1.0E+35
      DO 24 I = 1, NP
 24   FABS(I) = AMIN1(TERP(I), FABS(I) )
      DO  25  I = 2, NP
 25   TERM(I) = FSGN(I)*LN(I-1)*DIFF(I-1)
      TERM(1) = 1.0E+35
C
C     CORRECT  THE  FLUXES  COMPLETELY  NOW.
C
      DO  35  I = 1, NP
 35   DIFF(I) = AMIN1(FABS(I), TERM(I))
      DO  36  I = 1, NP
 36   FLXH(I) = AMAX1( 0.0, DIFF(I))
      DO  37  I = 1, NP
 37   FLXH(I) = FSGN(I)*FLXH(I)
C
C     CALCULATE  THE  NEW  FLUX-CORRECTED  DENSITIES.
C
      DO  40  I = 1,N
 40   LNRHOT(I) = LNRHOT(I) - FLXH(I+1) + FLXH(I)
      DO  41  I = 1, N
      SOURCE(I) = 0.0
 41   RHON(I)   =  LNRHOT(I)*RLN(I)
      LSOURC    =  .FALSE.
      GO TO 600
C
C
C    -------------------------------------------------------------------
C
      ENTRY   VELOCE1  ( U, N, UR, UL, DT )
C
C
CD    ******************************************************************
CD
CD    VELOCE  (U, N, UR, UL, DT )
CD    DESCRIPTION : THIS  ENTRY  CALCULATES  ALL  VELOCITY-DEPENDANT
CD                  COEFFS.
CD
CD    ARGUMENTS :
CD    U    REAL  ARRAY(N)   FLOW  VELOCITY  AT  THE  GRID  POINTS
CD    N    INTEGER          NUMBER  OF  INTERIOR  GRID  POINTS
CD    UR   REAL             VELOCITY  OF  FLOW  AT  RIGHT  BOUNDARY
CD    UL   REAL             VELOCITY  OF  FLOW  AT  LEFT   BOUNDARY
CD    DT   REAL             STEPSIZE  FOR THE  TIME  INTEGRATION
CD
CD    ******************************************************************
C
C     CALCULATES  THE  INTERFACE  AREA  X  VELOCITY  DIFFERENTIAL  X  DT.
C
      NP = N + 1
      DTH= 0.5*DT
      DO  101  I = 2, N
 101  ADUDTH(I)  = AH(I)*DTH*(U(I) + U(I-1)) - ADUGTH(I)
      ADUDTH(1)  = AH(1)*DT*UL - ADUGTH(1)
      ADUDTH(NP) = AH(NP)*DT*UR- ADUGTH(NP)
C
C     CALCULATE  THE  HALF-CELL  EPSILON  (V*DT/DX)
C
      DO  102  I = 1, NP
 102  EPSH(I) = ADUDTH(I)*RLH(I)
C
C     NEXT  CALCULATE  THE  DIFFUSION  AND  ANTIDIFFUSION  COEFFICIENTS.
C     VARIATION  WITH  EPSILON  MEANS  FOURTH-ORDER  ACCURATE  PHASES.
C
      DO  103  I = 1, NP
      NULH(I) = 0.16666667 + 0.3333333*EPSH(I)*EPSH(I)
 103  MULH(I) = 0.25 - 0.5*NULH(I)
      DO 104 I = 1, NP
      NULH(I) = LH(I)*NULH(I)
 104  MULH(I) = LH(I)*MULH(I)
      RETURN
C
C
C     ------------------------------------------------------------------
C
      ENTRY  NGRIDE1  ( RADN, N, RADR, RADL, ALPHA )
C
C
CD    ******************************************************************
CD
CD    NGRIDE  (RADN, N, RADR, RADL, ALPHA)
CD    DESCRIPTION :  THIS  ENTRY  SETS  NEW  GEOMETRY  VARIABLES  AND
CD                   COEFFS.
CD    ARGUMENTS :
CD    RADN      REAL  ARRAY(N)  NEW  GRID  POSITIONS
CD    N         INTEGER         NUMBER  OF  INTERIOR  GRID  POSITIONS
CD    RADR      REAL            POSITION  OF  THE  RIGHT  BOUNDARY
CD    RADL      REAL            POSITION  OF  THE  LEFT   BOUNDARY
CD    ALPHA     INTEGER         = 1 FOR  CARTESIAN  GEOMETRY
CD                              = 2 FOR  CYLINDRICAL GEOMETRY
CD                              = 3 FOR  SPHERICAL   GEOMETRY
CD
CD    *******************************************************************
C
C     CALCULATES  THE  NEW  HALF-CELL  POSITIONS  AND  GRID  CHANGES.
C
      NP = N + 1
      DO  202  I = 2, N
      RNH(I) = 0.5*(RADN(I) + RADN(I-1) )
C     WRITE ( 3,* ) ( RADN(I), RADN(I-1),RNH(I) )
 202  CONTINUE
      RNH(1) = RADL
      RNH(NP)= RADR
C
C     CALCULATE  THE  THREE  COORDINATE  SYSTEMS.
C
      GO TO ( 203, 206, 209 ), ALPHA
C
C     CARTESIAN  COORDINATES.
C
 203  DO  204  I = 1, NP
 204  AH(I) = 1.0
      DO  205  I = 1, N
      LN(I) = RNH(I+1) - RNH(I)
C     WRITE ( 3,* )( RNH(I),RNH(I+1),LN(I) )
 205  CONTINUE
      GOTO  213
C
C     CYLINDRICAL  COORDINATES.
C
 206  DO  207  I = 1, NP
      DIFF(I) = RNH(I)*RNH(I)
 207  AH(I) = PI*(ROH(I) + RNH(I))
      DO 208 I = 1, N
 208  LN(I) = PI*(DIFF(I+1) - DIFF(I) )
      GOTO  213
C
C     SPHERICAL   COORDINATES.
C
 209  DO  210  I = 1, NP
      DIFF(I) = RNH(I)*RNH(I)*RNH(I)
 210  SCRH(I) = (ROH(I) + RNH(I))*ROH(I)
      DO  211  I = 1, NP
 211  AH(I)  =  FTPI*(SCRH(I) + RNH(I)*RNH(I))
      DO  212  I = 1, N
 212  LN(I)  =  FTPI*(DIFF(I+1) - DIFF(I) )
C
C     NOW  THE  GEOMETRIC  VARIABLES  WHICH  ARE  SYSTEM  INDEPENDENT.
C
 213  DO  214  I = 2, N
 214  LH(I) = 0.5*(LN(I) + LN(I-1))
      LH(1) = LN(1)
      LH(NP)= LN(N)
      DO  215  I = 1, N
 215  RLN(I) = 1.0/LN(I)
      DO  216  I = 1, NP
 216  ADUGTH(I) = AH(I)*(RNH(I) - ROH(I))
      DO  217  I = 2, N
 217  RLH(I)  = 0.5*(RLN(I) + RLN(I-1))
      RLH(1)  = RLN(1)
      RLH(NP) = RLN(N)
      RETURN
C
C
C     -----------------------------------------------------------------
C
      ENTRY  SOURCQ1 ( N, DT, MODES, C, D, DR, DL )
C
C
CD    ******************************************************************
CD
CD    SOURCE ( N, DT, MODES, C, D, DR, DL )
CD    DESCRIPTION: THIS  ENTRY  ACCUMULATES  DIFFERENT  SOURCE  TERMS.
CD
CD    ARGUMENTS:
CD    N   INTEGER        NUMBER  OF  INTERIOR  GRID  POINTS
CD    DT  REAL           STEPSIZE  FOR  THE  TIME  INTEGRATION
CD    MODES  INTEGER     = 1  COMPUTES  +  DIV (D)
CD                       = 2  COMPUTES  +  C*GRAD(D)
CD                       = 3  ADDS  +  D TO THE  SOURCES
CD    C   REAL  ARRAY(N) ARRAY  OF  SOURCE  VARIABLES  AT  GRID  POINTS
CD    D   REAL  ARRAY(N) ARRAY  OF  SOURCE  VARIABLES  AT  GRID  POINTS
CD    DR  REAL           RIGHT  BOUNDARY  VALUE  OF  D
CD    DL  REAL           LEFT   BOUNDARY  VALUE  OF  D
CD
CD    *****************************************************************
C
      if (first) then					! initialize to 0.0 (CER)
          do i = 1, N
	      source(i) = 0.0
	  end do
	  first = .false.
      end if
      NP = N + 1
      DTH = 0.5*DT
      DTQ = 0.25*DT
      GOTO  ( 310, 320, 330 ), MODES
C
C     DIV (D) IS  COMPUTED  CONSERVATIVELY  AND  ADDED  TO  THE  SOURCES.
C
 310  DO  311  I = 2, N
 311  SCRH(I) = DTH*AH(I)*(D(I) + D(I-1))
      SCRH(1) = DT*AH(1)*DL
      SCRH(NP)= DT*AH(NP)*DR
      DO  312  I = 1, N
 312  SOURCE(I) = SOURCE(I) + SCRH(I+1) - SCRH(I)
      LSOURC  = .TRUE.
      RETURN
C
C     C*GRAD(D) IS  COMPUTED  EFFICIENTLY  AND  ADDED  TO  THE  SOURCES.
C
 320  DO  321  I = 2, N
 321  SCRH(I) = DTQ*(D(I) + D(I-1))
      SCRH(1) = DTH*DL
      SCRH(NP)= DTH*DR
      DO  322  I = 1, N
 322  DIFF(I) = SCRH(I+1) - SCRH(I)
      DO  323  I = 1, N
 323  SOURCE(I) = SOURCE(I) + C(I)*(AH(I+1) +AH(I))*DIFF(I)
      LSOURC = .TRUE.
      RETURN
C
C     D  IS  ADDED  TO  THE  SOURCES  IN  AN  EXPLICIT  FORMULATION.
C
 330  DO  331  I = 1, N
 331  SOURCE(I) = SOURCE(I) + DT*LO(I)*D(I)
      LSOURC  =  .TRUE.
      RETURN
C
C
C     -----------------------------------------------------------------
C
      ENTRY   OGRIDE1 (N)
C
CD    *****************************************************************
CD
CD    OGRIDE   (N)
CD    DESCRIPTION:  THIS  ENTRY  COPIES  OLD  GRID  AND  GEOMETRY  VARIABLES.
CD
CD    ARGUMENTS : 
CD    N   INTEGER     NUMBER  OF  INTERIOR  GRID  POINTS
CD
CD    ******************************************************************
C
C      COPY  THE  PREVIOUSLY  NEW  GRID  VALUES  TO  BE  USED  AS  THE
C            OLD  GRID.
C
       NP = N + 1
       DO  401  I = 1, N
       LO(I) = LN(I)
 401   RLO(I) = RLN(I)
       DO  402  I = 1, NP
 402   ROH(I) = RNH(I)
       RETURN
C
C
C      ------------------------------------------------------------------
C
       ENTRY  CONSRE1  ( RHO, N, CSUM )
C
C
CD     *****************************************************************
CD
CD     CONSRE  ( RHO, N, CSUM )
CD     DESCRIPTION : THIS  ENTRY  COMPUTES  THE  OSTENSIBLY  CONSERVED
CD                   SUM.
CD
CD     ARGUMENTS:
CD     RHO     REAL  ARRAY(N)  GRID  PT. VALUES  FOR  CONSERVATION  SUM.
CD     N       INTEGER         NUMBER  OF  INTERIOR  GRID  POINTS.
CD     CSUM    REAL            VALUE   OF  THE  CONSERVATION  SUM OF RHO
CD
CD     *****************************************************************
CD
C      COMPUTE  THE  OSTENSIBLY  CONSERVED  TOTAL  MASS ( BEWARE YOUR B.C. )
C
       CSUM  =  0.0
       DO  501  I = 1, N
 501   CSUM  =  CSUM  +  LN(I)*RHO(I)
 600   RETURN
        END
C#####################################################################
C     THIS SUBROUTINE CALCULATES AND THEN PRINTS THE BALANCE OF TERMS.
C#####################################################################
      SUBROUTINE SUMTRM
C
      DIMENSION ROLD(252),RNEW(252),VOLD(252),VNEW(252),RADO(252)
      COMMON/CDAT/ROLD,RNEW,VOLD,VNEW,RADO,TRATIO,GAMA,REOL2,JSTEP,
     >   DELTAX,DELTAT,NX,MAXSTP,RIN
C
      DELX2=2.*DELTAX
C
      WRITE(8,10)JSTEP,DELTAX,DELTAT,GAMA,TRATIO,NX,MAXSTP,RIN
   10 FORMAT(////,1X,'FIRST EQUATION - BALANCE OF TERMS - TIME STEP='
     >  ,I5,/,1X,'OTHER PARAMETERS: DELTAX,DELTAT,GAMA,TRATIO,NX,
     >MAXSTP,RIN=',/,1P,4E10.2,0P,5X,2I6,5X,1P,E10.2,0P,//)
      WRITE(8,15)
   15 FORMAT(5X,'    I',9X,'RADIUS',10X,'DN/DT',8X,'N(DV/DR)',8X,
     >   'V(DN/DR)',7X,'3NV/R',10X,'SUM')
C
      DO 20 I=2,(NX-1)
      IM1=I-1
      IP1=I+1
      TERM1=(RNEW(I)-ROLD(I))/DELTAT
      TERM2=RNEW(I)*(VNEW(IP1)-VNEW(IM1))/DELX2
      TERM3=VNEW(I)*(RNEW(IP1)-RNEW(IM1))/DELX2
      TERM4=3.*RNEW(I)*VNEW(I)/(RIN+RADO(I)-RADO(1))
      SUM=TERM1+TERM2+TERM3+TERM4
C
      WRITE(8,30)I,RADO(I),TERM1,TERM2,TERM3,TERM4,SUM
   30 FORMAT(5X,I5,1P,6E15.2,0P)
C
   20 CONTINUE
C
C
      WRITE(8,39)
   39 FORMAT(///,1X,'SECOND EQUATION - BALANCE OF TERMS',/)
      WRITE(8,40)
   40 FORMAT(5X,'    I',9X,'R',14X,'DV/DT',7X,'V(DV/DR)',7X,
     >   'F/N(DN/DR)',5X,'REOL2/(3+R)**2',4X,'SUM')
C
      DO 50 I=2,(NX-1)
      IM1=I-1
      IP1=I+1
      TERM1=(VNEW(I)-VOLD(I))/DELTAT
      TERM2=VNEW(I)*(VNEW(IP1)-VNEW(IM1))/DELX2
      TERM3=(1.+TRATIO*GAMA*RNEW(I)**GAMA/RNEW(I))/RNEW(I)
     >   *(RNEW(IP1)-RNEW(IM1))/DELX2
      FR=RIN+RADO(I)-RADO(1)
      TERM4=REOL2/FR**2
      SUM=TERM1+TERM2+TERM3+TERM4
C
      WRITE(8,60) I,RADO(I),TERM1,TERM2,TERM3,TERM4,SUM
   60 FORMAT(5X,I5,1P,6E15.2,0P)
   50 CONTINUE
C
      RETURN
      END


c@@@@@@@@@@@@@@@@@@  Periodic routines @@@@@@@@@@@@@@@@@@@@@

      SUBROUTINE  PRBFCT(RHOO, RHON, NP)
C
      REAL   RHOO(*),  RHON(*)
CD
CD    ***************************************************************
CD
CD    PRBFCT(RHOO, RHON, NP)                           D.3
CD    ORIGINATOR - J.P. BORIS   CODE  7706,  NRL       OCT. 1975
CD
CD    DESCRIPTION:  THIS  ROUTINE  SOLVES  GENERALISED  CONTINUITY  
CD    EQUATIONS  OF  THE  FORM  DRHO/DT = -DIV (RHO*V) + SOURCES WITH
CD    PERIODIC BOUNDARY CONDITIONS IN CARTESIAN GEOMETRY.  THE
CD    FINITE-DIFFERENCE  GRID  CAN  BE  EULERIAN, SLIDING  REZONE
CD    OR  LAGRANGIAN   AND  CAN  BE  ARBITRARILY  SPACED.  THE
CD    ALGORITHM  USED  IS  A  LOW-PHASE  -ERROR  FCT  ALGORITHM,
CD    VECTORISED  AND  OPTIMISED  FOR  SPEED.
CD
CD    ARGUMENTS: IN THIS ROUTINE ALL INPUT VARIABLES MUST BE GIVEN AT
CD    NP = N + 1 POINTS WHERE RHO(NP) = RHO(1) AND THE LENGTH OF THE
CD    SYSTEM IS RADN(NP) - RADN(1)
CD    RHOO      REAL  ARRAY(NP)  GRID POINT DENSITIES AT START OF STEP
CD    RHON      REAL  ARRAY(NP)  GRID POINT DENSITIES AT END OF STEP
CD    NP        INTEGER          NUMBER OF GRID POINTS PLUS ONE
CD    
CD    LANGUAGE  AND  LIMITATIONS: THE  SUBROUTINE  PRBFCT  IS  A
CD    MULTIPLE - ENTRY  FORTRAN  ROUTINE  IN  SINGLE  PRECISION ( 32
CD    BITS  ASC). THE  ASC  PARAMETER  STATEMENT  IS  USED  TO  SET
CD    SYMBOLICALLY  THE  INTERNAL  ARRAY  DIMENSIONS.  UNDERFLOWS  ARE
CD    POSSIBLE  WHEN  THE  FUNCTION  BEING  TRANSFORMED  HAS  MANY  
CD    ZEROES.  THE  CALCULATIONS  GENERALLY  MISCONSERVE  BY  ONE  OR
CD    TWO  BITS  PER  CYCLE.  THE  RELATIVE  PHASE  AND  AMPLITUDE  
CD    ERRORS ( FOR  SMOOTH  FUNCTIONS )  ARE  TYPICALLY  SEVERAL  PERCENT
CD    FOR  CHARACTERISTIC  LENGTHS  OF  1-2  CELLS ( WAVELENGTHS  OF
CD    ORDER  10  CELLS ).  SHOCKS  ARE  GENERALLY  ACCURATE  TO  BETTER
CD    THEN  1  PERCENT.  THIS  SUBROUTINE  MUST  BE  COMPILED  WITH  THE
CD    Y  OPTION  TO  FORCE  STORAGE  AND  RETENTION  OF  INTERNAL  
CD    VARIABLES.  ALTERNATIVELY  A  COMMON  BLOCK  CAN  BE  ADDED  TO
CD    ACCOMPLISH  THE  SAME  END.
CD
CD    ENTRY  POINTS : OGRIDP, NGRIDP, VELOCP, SOURCP, CONSRP, SEE #11
CD    OF  THE  DETAILED  DOCUMENTATION  ( OR  THE  LISTING  BELOW) FOR
CD    THE  EXPLANATIONS  AND  USE  OF  THE  ARGUMENTS  TO  THESE  ENTRIES.
CD
CD    NO  AUXILIARY  OR  LIBRARY  ROUTINES  ARE  CALLED  BY  PRBFCT.
CD
CD
CD    *****************************************************************
C
      PARAMETER (NPT = 350)
      LOGICAL       LSOURC, first		! initialize variables (cer)
C      REAL	    MASK1, MASK2, MASK3
      REAL          SOURCE(NPT), SCRH(NPT), RHOT(NPT), DIFF(NPT)
      REAL          ADUGTH(NPT), FLXH(NPT), NULH(NPT), MULH(NPT)
      REAL          LNRHOT(NPT), FSGN(NPT), FABS(NPT), EPSH(NPT)
      REAL          LORHOT(NPT), TERP(NPT), TERM(NPT), ADUDTH(NPT)
      REAL          LO(NPT), LN(NPT), LH(NPT), RLO(NPT), RLN(NPT)
      REAL          RNH(NPT), ROH(NPT), RLH(NPT)
      REAL	    U(*)
      REAL	    RADN(*)
      REAL	    C(*), D(*)
      REAL          RHO(*)
C      DATA	    MASK1, MASK2, MASK3 /$80000000, 1.0, $7FFFFFFF/
      DATA	    ROH /NPT*1.0/, SOURCE/NPT*0.0/, LSOURC/.false./
      EQUIVALENCE   (EPSH(1), SCRH(1)),     (LNRHOT(1), LORHOT(1))
      EQUIVALENCE   (FLXH(1), SCRH(1)),     (FSGN(1),  RHOT(1))
      EQUIVALENCE   (FABS(1), SCRH(1)),     (TERP(1),  SOURCE(1))
      EQUIVALENCE   (TERM(1), SOURCE(1))
      save					! needed for MPW FORTRAN (cer)
C
C     CALCULATE  THE  DIFFUSIVE  AND  CONVECTIVE  FLUXES.
C
      N = NP - 1
      DO  11  I  = 2, NP
      FLXH(I) = 0.5*ADUDTH(I)*(RHOO(I) + RHOO(I-1))
 11   DIFF(I) = NULH(I)*(RHOO(I) - RHOO(I-1))
      FLXH(1) = FLXH(NP)
      DIFF(1) = DIFF(NP)
C
C     CALCULATE  LAMBDAO*RHOT,  THE  TRANSPORTED  MASS  ELEMENTS.
C
      DO  12  I  = 1, N
 12   LORHOT(I) = LO(I)*RHOO(I) - FLXH(I+1) + FLXH(I)
C
C     ADD  IN  THE  SOURCE  TERMS  AS  APPROPRIATE.
C
      IF  (.NOT. LSOURC ) GO TO 14
      DO  13  I = 1, N
 13   LORHOT(I) = LORHOT(I) + SOURCE(I)
C
C     CALCULATE  THE  PHOENICAL  ANTIDIFFUSIVE  FLUXES  HERE.
C
 14   DO  16  I = 1, N
 16   RHOT(I) = LORHOT(I)*RLO(I)
      RHOT(NP) = RHOT(1)
      DO  17  I = 2, NP
 17   FLXH(I) = MULH(I)*(RHOT(I) - RHOT(I-1))
      FLXH(1) = FLXH(NP)
C
C  DIFFUSE THE SOLUTION RHOT USING OLD FLUXES.   
C
      DO  18  I = 1, N
 18   LNRHOT(I) = LORHOT(I) + DIFF(I+1) - DIFF(I)
C
C  CALCULATE THE TRANSPORTED/DIFFUSED DENSITY AND GRID DIFFERENCES.
C
      DO  19  I = 1, N
 19   RHOT(I) = LNRHOT(I)*RLN(I)
      RHOT(NP) = RHOT(1)
      DO  20  I = 2, NP
 20   DIFF(I) = RHOT(I) - RHOT(I-1)
      DIFF(1) = DIFF(NP)
C
C     CALCULATE THE SIGN AND MAGNITUDE OF THE ANTIDIFFUSIVE FLUX.
C
c...      DO 21 I = 1, NP
c... 21   FSGN(I) = AND(MASK1, DIFF(I))
c...      DO 22 I = 1, NP
c...      FABS(I) = AND(MASK3, FLXH(I))
c... 22   FSGN(I) = OR(MASK2, FSGN(I))      
C
C  (CER) replace the above two loops
C
      do 21 i = 1, np
          if (diff(i) .lt. 0) then
	      fsgn(i) = -1.0
	  else
	      fsgn(i) =  1.0
	  end if
	  if (flxh(i) .lt. 0) then
	      fabs(i) = -flxh(i)
	  else
	      fabs(i) =  flxh(i)
	  end if
 21   continue
C
C     CALCULATE  THE  FLUX-LIMITING  CHANGES  ON  THE RIGHT AND  THE  LEFT.
C
      DO 23 I = 1, N
 23   TERP(I) = FSGN(I)*LN(I)*DIFF(I+1)
      TERP(NP)= TERP(1)
      DO 24 I = 1, NP
 24   FABS(I) = AMIN1(TERP(I), FABS(I) )
      DO 25 I = 2, NP
 25   TERM(I) = FSGN(I)*LN(I-1)*DIFF(I-1)
      TERM(1) = TERM(NP)
C
C     CORRECT  THE  FLUXES  COMPLETELY  NOW.
C
      DO 35 I = 1, NP
 35   DIFF(I) = AMIN1(FABS(I), TERM(I))
      DO 36 I = 1, NP
 36   FLXH(I) = AMAX1( 0.0, DIFF(I))
      DO 37 I = 1, NP
 37   FLXH(I) = FSGN(I)*FLXH(I)
C
C     CALCULATE  THE  NEW  FLUX-CORRECTED  DENSITIES.
C
      DO 40 I = 1, N
 40   LNRHOT(I) = LNRHOT(I) - FLXH(I+1) + FLXH(I)
      LNRHOT(NP) = LNRHOT(1)
      DO 41 I = 1, NP
      SOURCE(I) = 0.0
 41   RHON(I) =  LNRHOT(I)*RLN(I)
      LSOURC = .FALSE.
      RETURN
C
C
C    -------------------------------------------------------------------
C
      ENTRY VELOCP(U, NP, DT)
C
CD    ******************************************************************
CD
CD    VELOCP (U, NP, DT )
CD    DESCRIPTION : THIS ENTRY CALCULATES ALL VELOCITY-DEPENDANT COEFFS.
CD
CD    ARGUMENTS :
CD    U    REAL  ARRAY(NP)   FLOW VELOCITY AT THE GRID POINTS
CD    NP   INTEGER           NUMBER OF GRID POINTS PLUS ONE
CD    DT   REAL              STEPSIZE FOR THE TIME INTEGRATION
CD
CD    ******************************************************************
C
C     CALCULATES  THE  INTERFACE  AREA  X  VELOCITY  DIFFERENTIAL  X  DT.
C
      N = NP - 1
      DTH= 0.5*DT
      DO 101 I = 2, NP
 101  ADUDTH(I) = DTH*(U(I) + U(I-1)) - ADUGTH(I)
      ADUDTH(1) = ADUDTH(NP)
C
C     CALCULATE  THE  HALF-CELL  EPSILON  (V*DT/DX)
C
      DO 102 I = 1, NP
 102  EPSH(I) = ADUDTH(I)*RLH(I)
C
C     NEXT  CALCULATE  THE  DIFFUSION  AND  ANTIDIFFUSION  COEFFICIENTS.
C     VARIATION  WITH  EPSILON  MEANS  FOURTH-ORDER  ACCURATE  PHASES.
C
      DO 103 I = 1, NP
      NULH(I) = 0.16666667 + 0.3333333*EPSH(I)*EPSH(I)
 103  MULH(I) = 0.25 - 0.5*NULH(I)
      DO 104 I = 1, NP
      NULH(I) = LH(I)*NULH(I)
 104  MULH(I) = LH(I)*MULH(I)
      RETURN
C
C
C     ------------------------------------------------------------------
C
      ENTRY  NGRIDP(RADN, NP)
C
C
CD    ******************************************************************
CD
CD    NGRIDP (RADN, NP)
CD    DESCRIPTION :  THIS  ENTRY  SETS  NEW  GEOMETRY  VARIABLES  AND
CD                   COEFFS.
CD    ARGUMENTS :
CD    RADN      REAL  ARRAY(NP) NEW GRID POINT POSITIONS
CD    NP        INTEGER         NUMBER OF GRID POINTS PLUS ONE
CD
CD    *******************************************************************
C
C     CALCULATES  THE  NEW  HALF-CELL  POSITIONS  AND  GRID  CHANGES.
C
      N = NP - 1
      DO 202 I = 2, NP
 202  RNH(I) = 0.5*(RADN(I) + RADN(I-1) )
      RNH(1) = RADN(1) - 0.5*(RADN(NP) - RADN(N))
C
C  CALCULATE THE CELL VOLUMES AND OTHER GEOMETRIC VARIABLES FOR
C  RETENTION AND GLOBAL USE.
C
      DO 205 I = 1, N
 205  LN(I) = RNH(I+1) - RNH(I)
      LN(NP) = LN(1)
      DO 213 I = 2, NP
 213  LH(I) = 0.5*(LN(I) + LN(I-1))
      LH(1) = LN(NP)
      DO 214 I = 1, NP
 214  RLN(I) = 1.0/LN(I)
      DO  215  I = 1, NP
 215  ADUGTH(I) = (RNH(I) - ROH(I))
      DO 216 I = 2, NP
 216  RLH(I)  = 0.5*(RLN(I) + RLN(I-1))
      RLH(1)  = RLH(NP)
      RETURN
C
C
C     -----------------------------------------------------------------
C
      ENTRY  SOURCP(NP, DT, MODES, C, D)
C
CD    ******************************************************************
CD
CD    SOURCP (NP, DT, MODES, C, D)
CD    DESCRIPTION: THIS  ENTRY  ACCUMULATES  DIFFERENT  SOURCE  TERMS.
CD
CD    ARGUMENTS:
CD    NP  INTEGER        NUMBER OF GRID POINTS PLUS ONE.
CD    DT  REAL           STEPSIZE FOR THE TIME INTEGRATION
CD    MODES  INTEGER     = 1  COMPUTES  +  DIV (D)
CD                       = 2  COMPUTES  +  C*GRAD(D)
CD                       = 3  ADDS  +  D TO THE  SOURCES
CD    C   REAL  ARRAY(N) ARRAY  OF  SOURCE  VARIABLES  AT  GRID  POINTS
CD    D   REAL  ARRAY(N) ARRAY  OF  SOURCE  VARIABLES  AT  GRID  POINTS
CD
CD    *****************************************************************
C
      if (first) then					! initialize to 0.0 (CER)
          do i = 1, NP
	      source(i) = 0.0
	  end do
	  first = .false.
      end if
      N = NP - 1
      DTH = 0.5*DT
      GOTO  ( 310, 320, 330 ), MODES
C
C     DIV(D) IS  COMPUTED  CONSERVATIVELY  AND  ADDED  TO  THE  SOURCES.
C
 310  DO 311 I = 2, NP
 311  SCRH(I) = DTH*(D(I) + D(I-1))
      SCRH(1) = SCRH(NP)
      DO  312  I = 1, N
 312  SOURCE(I) = SOURCE(I) + SCRH(I+1) - SCRH(I)
      SOURCE(NP) = SOURCE(1)
      LSOURC  = .TRUE.
      RETURN
C
C     C*GRAD(D) IS  COMPUTED  EFFICIENTLY  AND  ADDED  TO  THE  SOURCES.
C
 320  DO  321  I = 2, NP
 321  SCRH(I) = DTH*(D(I) + D(I-1))
      SCRH(1) = SCRH(NP)
      DO  322  I = 1, N
 322  SOURCE(I) = SOURCE(I) + C(I)*(SCRH(I+1) - SCRH(I))
      SOURCE(NP) = SOURCE(1)
      LSOURC = .TRUE.
      RETURN
C
C     D  IS  ADDED  TO  THE  SOURCES  IN  AN  EXPLICIT  FORMULATION.
C
 330  DO  331  I = 1, NP
 331  SOURCE(I) = SOURCE(I) + DT*LO(I)*D(I)
      LSOURC  =  .TRUE.
      RETURN
C
C
C     -----------------------------------------------------------------
C
      ENTRY OGRIDP(NP)
C
CD    *****************************************************************
CD
CD    OGRIDP (NP)
CD    DESCRIPTION:  THIS  ENTRY  COPIES  OLD  GRID  AND  GEOMETRY  VARIABLES.
CD
CD    ARGUMENTS : 
CD    NP   INTEGER     NUMBER OF GRID POINTS PLUS ONE
CD
CD    ******************************************************************
C
C      COPY  THE  PREVIOUSLY  NEW  GRID  VALUES  TO  BE  USED  AS  THE
C            OLD  GRID.
C
       N = NP - 1
       DO  401  I = 1, NP
       ROH(I) = RNH(I)
       LO(I) = LN(I)
 401   RLO(I) = RLN(I)
       RETURN
C
C
C      ------------------------------------------------------------------
C
       ENTRY CONSRP(RHO, NP, CSUM)
C
CD     *****************************************************************
CD
CD     CONSRP(RHO, NP, CSUM)
CD     DESCRIPTION : THIS  ENTRY  COMPUTES  THE  OSTENSIBLY  CONSERVED
CD                   SUM.
CD
CD     ARGUMENTS:
CD     RHO     REAL  ARRAY(N)  GRID  PT. VALUES  FOR  CONSERVATION  SUM.
CD     NP      INTEGER         NUMBER OF GRID POINTS PLUS ONE.
CD     CSUM    REAL            VALUE   OF  THE  CONSERVATION  SUM OF RHO
CD
CD     *****************************************************************
CD
C      COMPUTE  THE  OSTENSIBLY  CONSERVED  TOTAL  MASS ( BEWARE YOUR B.C. )
C
       N = NP - 1
       CSUM  =  0.0
       DO  501  I = 1, N
 501   CSUM  =  CSUM  +  LN(I)*RHO(I)
       RETURN
       END


c@@@@@@@@@@@@@@@@@@  Periodic routines @@@@@@@@@@@@@@@@@@@@@

      SUBROUTINE  PRBFCT1(RHOO, RHON, NP)
C
      REAL   RHOO(*),  RHON(*)
CD
CD    ***************************************************************
CD
CD    PRBFCT(RHOO, RHON, NP)                           D.3
CD    ORIGINATOR - J.P. BORIS   CODE  7706,  NRL       OCT. 1975
CD
CD    DESCRIPTION:  THIS  ROUTINE  SOLVES  GENERALISED  CONTINUITY  
CD    EQUATIONS  OF  THE  FORM  DRHO/DT = -DIV (RHO*V) + SOURCES WITH
CD    PERIODIC BOUNDARY CONDITIONS IN CARTESIAN GEOMETRY.  THE
CD    FINITE-DIFFERENCE  GRID  CAN  BE  EULERIAN, SLIDING  REZONE
CD    OR  LAGRANGIAN   AND  CAN  BE  ARBITRARILY  SPACED.  THE
CD    ALGORITHM  USED  IS  A  LOW-PHASE  -ERROR  FCT  ALGORITHM,
CD    VECTORISED  AND  OPTIMISED  FOR  SPEED.
CD
CD    ARGUMENTS: IN THIS ROUTINE ALL INPUT VARIABLES MUST BE GIVEN AT
CD    NP = N + 1 POINTS WHERE RHO(NP) = RHO(1) AND THE LENGTH OF THE
CD    SYSTEM IS RADN(NP) - RADN(1)
CD    RHOO      REAL  ARRAY(NP)  GRID POINT DENSITIES AT START OF STEP
CD    RHON      REAL  ARRAY(NP)  GRID POINT DENSITIES AT END OF STEP
CD    NP        INTEGER          NUMBER OF GRID POINTS PLUS ONE
CD    
CD    LANGUAGE  AND  LIMITATIONS: THE  SUBROUTINE  PRBFCT  IS  A
CD    MULTIPLE - ENTRY  FORTRAN  ROUTINE  IN  SINGLE  PRECISION ( 32
CD    BITS  ASC). THE  ASC  PARAMETER  STATEMENT  IS  USED  TO  SET
CD    SYMBOLICALLY  THE  INTERNAL  ARRAY  DIMENSIONS.  UNDERFLOWS  ARE
CD    POSSIBLE  WHEN  THE  FUNCTION  BEING  TRANSFORMED  HAS  MANY  
CD    ZEROES.  THE  CALCULATIONS  GENERALLY  MISCONSERVE  BY  ONE  OR
CD    TWO  BITS  PER  CYCLE.  THE  RELATIVE  PHASE  AND  AMPLITUDE  
CD    ERRORS ( FOR  SMOOTH  FUNCTIONS )  ARE  TYPICALLY  SEVERAL  PERCENT
CD    FOR  CHARACTERISTIC  LENGTHS  OF  1-2  CELLS ( WAVELENGTHS  OF
CD    ORDER  10  CELLS ).  SHOCKS  ARE  GENERALLY  ACCURATE  TO  BETTER
CD    THEN  1  PERCENT.  THIS  SUBROUTINE  MUST  BE  COMPILED  WITH  THE
CD    Y  OPTION  TO  FORCE  STORAGE  AND  RETENTION  OF  INTERNAL  
CD    VARIABLES.  ALTERNATIVELY  A  COMMON  BLOCK  CAN  BE  ADDED  TO
CD    ACCOMPLISH  THE  SAME  END.
CD
CD    ENTRY  POINTS : OGRIDP1, NGRIDP1, VELOCP1, SOURCP1, CONSRP1, SEE #11
CD    OF  THE  DETAILED  DOCUMENTATION  ( OR  THE  LISTING  BELOW) FOR
CD    THE  EXPLANATIONS  AND  USE  OF  THE  ARGUMENTS  TO  THESE  ENTRIES.
CD
CD    NO  AUXILIARY  OR  LIBRARY  ROUTINES  ARE  CALLED  BY  PRBFCT.
CD
CD
CD    *****************************************************************
C
      PARAMETER (NPT = 350)
      LOGICAL       LSOURC, first		! initialize variables (cer)
C      REAL	    MASK1, MASK2, MASK3
      REAL          SOURCE(NPT), SCRH(NPT), RHOT(NPT), DIFF(NPT)
      REAL          ADUGTH(NPT), FLXH(NPT), NULH(NPT), MULH(NPT)
      REAL          LNRHOT(NPT), FSGN(NPT), FABS(NPT), EPSH(NPT)
      REAL          LORHOT(NPT), TERP(NPT), TERM(NPT), ADUDTH(NPT)
      REAL          LO(NPT), LN(NPT), LH(NPT), RLO(NPT), RLN(NPT)
      REAL          RNH(NPT), ROH(NPT), RLH(NPT)
      REAL	    U(*)
      REAL	    RADN(*)
      REAL	    C(*), D(*)
      REAL          RHO(*)
C      DATA	    MASK1, MASK2, MASK3 /$80000000, 1.0, $7FFFFFFF/
      DATA	    ROH /NPT*1.0/, SOURCE/NPT*0.0/, LSOURC/.false./
      EQUIVALENCE   (EPSH(1), SCRH(1)),     (LNRHOT(1), LORHOT(1))
      EQUIVALENCE   (FLXH(1), SCRH(1)),     (FSGN(1),  RHOT(1))
      EQUIVALENCE   (FABS(1), SCRH(1)),     (TERP(1),  SOURCE(1))
      EQUIVALENCE   (TERM(1), SOURCE(1))
      save					! needed for MPW FORTRAN (cer)
C
C     CALCULATE  THE  DIFFUSIVE  AND  CONVECTIVE  FLUXES.
C
      N = NP - 1
      DO  11  I  = 2, NP
      FLXH(I) = 0.5*ADUDTH(I)*(RHOO(I) + RHOO(I-1))
 11   DIFF(I) = NULH(I)*(RHOO(I) - RHOO(I-1))
      FLXH(1) = FLXH(NP)
      DIFF(1) = DIFF(NP)
C
C     CALCULATE  LAMBDAO*RHOT,  THE  TRANSPORTED  MASS  ELEMENTS.
C
      DO  12  I  = 1, N
 12   LORHOT(I) = LO(I)*RHOO(I) - FLXH(I+1) + FLXH(I)
C
C     ADD  IN  THE  SOURCE  TERMS  AS  APPROPRIATE.
C
      IF  (.NOT. LSOURC ) GO TO 14
      DO  13  I = 1, N
 13   LORHOT(I) = LORHOT(I) + SOURCE(I)
C
C     CALCULATE  THE  PHOENICAL  ANTIDIFFUSIVE  FLUXES  HERE.
C
 14   DO  16  I = 1, N
 16   RHOT(I) = LORHOT(I)*RLO(I)
      RHOT(NP) = RHOT(1)
      DO  17  I = 2, NP
 17   FLXH(I) = MULH(I)*(RHOT(I) - RHOT(I-1))
      FLXH(1) = FLXH(NP)
C
C  DIFFUSE THE SOLUTION RHOT USING OLD FLUXES.   
C
      DO  18  I = 1, N
 18   LNRHOT(I) = LORHOT(I) + DIFF(I+1) - DIFF(I)
C
C  CALCULATE THE TRANSPORTED/DIFFUSED DENSITY AND GRID DIFFERENCES.
C
      DO  19  I = 1, N
 19   RHOT(I) = LNRHOT(I)*RLN(I)
      RHOT(NP) = RHOT(1)
      DO  20  I = 2, NP
 20   DIFF(I) = RHOT(I) - RHOT(I-1)
      DIFF(1) = DIFF(NP)
C
C     CALCULATE THE SIGN AND MAGNITUDE OF THE ANTIDIFFUSIVE FLUX.
C
c...      DO 21 I = 1, NP
c... 21   FSGN(I) = AND(MASK1, DIFF(I))
c...      DO 22 I = 1, NP
c...      FABS(I) = AND(MASK3, FLXH(I))
c... 22   FSGN(I) = OR(MASK2, FSGN(I))      
C
C  (CER) replace the above two loops
C
      do 21 i = 1, np
          if (diff(i) .lt. 0) then
	      fsgn(i) = -1.0
	  else
	      fsgn(i) =  1.0
	  end if
	  if (flxh(i) .lt. 0) then
	      fabs(i) = -flxh(i)
	  else
	      fabs(i) =  flxh(i)
	  end if
 21   continue
C
C     CALCULATE  THE  FLUX-LIMITING  CHANGES  ON  THE RIGHT AND  THE  LEFT.
C
      DO 23 I = 1, N
 23   TERP(I) = FSGN(I)*LN(I)*DIFF(I+1)
      TERP(NP)= TERP(1)
      DO 24 I = 1, NP
 24   FABS(I) = AMIN1(TERP(I), FABS(I) )
      DO 25 I = 2, NP
 25   TERM(I) = FSGN(I)*LN(I-1)*DIFF(I-1)
      TERM(1) = TERM(NP)
C
C     CORRECT  THE  FLUXES  COMPLETELY  NOW.
C
      DO 35 I = 1, NP
 35   DIFF(I) = AMIN1(FABS(I), TERM(I))
      DO 36 I = 1, NP
 36   FLXH(I) = AMAX1( 0.0, DIFF(I))
      DO 37 I = 1, NP
 37   FLXH(I) = FSGN(I)*FLXH(I)
C
C     CALCULATE  THE  NEW  FLUX-CORRECTED  DENSITIES.
C
      DO 40 I = 1, N
 40   LNRHOT(I) = LNRHOT(I) - FLXH(I+1) + FLXH(I)
      LNRHOT(NP) = LNRHOT(1)
      DO 41 I = 1, NP
      SOURCE(I) = 0.0
 41   RHON(I) =  LNRHOT(I)*RLN(I)
      LSOURC = .FALSE.
      RETURN
C
C
C    -------------------------------------------------------------------
C
      ENTRY VELOCP1(U, NP, DT)
C
CD    ******************************************************************
CD
CD    VELOCP (U, NP, DT )
CD    DESCRIPTION : THIS ENTRY CALCULATES ALL VELOCITY-DEPENDANT COEFFS.
CD
CD    ARGUMENTS :
CD    U    REAL  ARRAY(NP)   FLOW VELOCITY AT THE GRID POINTS
CD    NP   INTEGER           NUMBER OF GRID POINTS PLUS ONE
CD    DT   REAL              STEPSIZE FOR THE TIME INTEGRATION
CD
CD    ******************************************************************
C
C     CALCULATES  THE  INTERFACE  AREA  X  VELOCITY  DIFFERENTIAL  X  DT.
C
      N = NP - 1
      DTH= 0.5*DT
      DO 101 I = 2, NP
 101  ADUDTH(I) = DTH*(U(I) + U(I-1)) - ADUGTH(I)
      ADUDTH(1) = ADUDTH(NP)
C
C     CALCULATE  THE  HALF-CELL  EPSILON  (V*DT/DX)
C
      DO 102 I = 1, NP
 102  EPSH(I) = ADUDTH(I)*RLH(I)
C
C     NEXT  CALCULATE  THE  DIFFUSION  AND  ANTIDIFFUSION  COEFFICIENTS.
C     VARIATION  WITH  EPSILON  MEANS  FOURTH-ORDER  ACCURATE  PHASES.
C
      DO 103 I = 1, NP
      NULH(I) = 0.16666667 + 0.3333333*EPSH(I)*EPSH(I)
 103  MULH(I) = 0.25 - 0.5*NULH(I)
      DO 104 I = 1, NP
      NULH(I) = LH(I)*NULH(I)
 104  MULH(I) = LH(I)*MULH(I)
      RETURN
C
C
C     ------------------------------------------------------------------
C
      ENTRY  NGRIDP1(RADN, NP)
C
C
CD    ******************************************************************
CD
CD    NGRIDP (RADN, NP)
CD    DESCRIPTION :  THIS  ENTRY  SETS  NEW  GEOMETRY  VARIABLES  AND
CD                   COEFFS.
CD    ARGUMENTS :
CD    RADN      REAL  ARRAY(NP) NEW GRID POINT POSITIONS
CD    NP        INTEGER         NUMBER OF GRID POINTS PLUS ONE
CD
CD    *******************************************************************
C
C     CALCULATES  THE  NEW  HALF-CELL  POSITIONS  AND  GRID  CHANGES.
C
      N = NP - 1
      DO 202 I = 2, NP
 202  RNH(I) = 0.5*(RADN(I) + RADN(I-1) )
      RNH(1) = RADN(1) - 0.5*(RADN(NP) - RADN(N))
C
C  CALCULATE THE CELL VOLUMES AND OTHER GEOMETRIC VARIABLES FOR
C  RETENTION AND GLOBAL USE.
C
      DO 205 I = 1, N
 205  LN(I) = RNH(I+1) - RNH(I)
      LN(NP) = LN(1)
      DO 213 I = 2, NP
 213  LH(I) = 0.5*(LN(I) + LN(I-1))
      LH(1) = LN(NP)
      DO 214 I = 1, NP
 214  RLN(I) = 1.0/LN(I)
      DO  215  I = 1, NP
 215  ADUGTH(I) = (RNH(I) - ROH(I))
      DO 216 I = 2, NP
 216  RLH(I)  = 0.5*(RLN(I) + RLN(I-1))
      RLH(1)  = RLH(NP)
      RETURN
C
C
C     -----------------------------------------------------------------
C
      ENTRY  SOURCP1(NP, DT, MODES, C, D)
C
CD    ******************************************************************
CD
CD    SOURCP (NP, DT, MODES, C, D)
CD    DESCRIPTION: THIS  ENTRY  ACCUMULATES  DIFFERENT  SOURCE  TERMS.
CD
CD    ARGUMENTS:
CD    NP  INTEGER        NUMBER OF GRID POINTS PLUS ONE.
CD    DT  REAL           STEPSIZE FOR THE TIME INTEGRATION
CD    MODES  INTEGER     = 1  COMPUTES  +  DIV (D)
CD                       = 2  COMPUTES  +  C*GRAD(D)
CD                       = 3  ADDS  +  D TO THE  SOURCES
CD    C   REAL  ARRAY(N) ARRAY  OF  SOURCE  VARIABLES  AT  GRID  POINTS
CD    D   REAL  ARRAY(N) ARRAY  OF  SOURCE  VARIABLES  AT  GRID  POINTS
CD
CD    *****************************************************************
C
      if (first) then					! initialize to 0.0 (CER)
          do i = 1, NP
	      source(i) = 0.0
	  end do
	  first = .false.
      end if
      N = NP - 1
      DTH = 0.5*DT
      GOTO  ( 310, 320, 330 ), MODES
C
C     DIV(D) IS  COMPUTED  CONSERVATIVELY  AND  ADDED  TO  THE  SOURCES.
C
 310  DO 311 I = 2, NP
 311  SCRH(I) = DTH*(D(I) + D(I-1))
      SCRH(1) = SCRH(NP)
      DO  312  I = 1, N
 312  SOURCE(I) = SOURCE(I) + SCRH(I+1) - SCRH(I)
      SOURCE(NP) = SOURCE(1)
      LSOURC  = .TRUE.
      RETURN
C
C     C*GRAD(D) IS  COMPUTED  EFFICIENTLY  AND  ADDED  TO  THE  SOURCES.
C
 320  DO  321  I = 2, NP
 321  SCRH(I) = DTH*(D(I) + D(I-1))
      SCRH(1) = SCRH(NP)
      DO  322  I = 1, N
 322  SOURCE(I) = SOURCE(I) + C(I)*(SCRH(I+1) - SCRH(I))
      SOURCE(NP) = SOURCE(1)
      LSOURC = .TRUE.
      RETURN
C
C     D  IS  ADDED  TO  THE  SOURCES  IN  AN  EXPLICIT  FORMULATION.
C
 330  DO  331  I = 1, NP
 331  SOURCE(I) = SOURCE(I) + DT*LO(I)*D(I)
      LSOURC  =  .TRUE.
      RETURN
C
C
C     -----------------------------------------------------------------
C
      ENTRY OGRIDP1(NP)
C
CD    *****************************************************************
CD
CD    OGRIDP (NP)
CD    DESCRIPTION:  THIS  ENTRY  COPIES  OLD  GRID  AND  GEOMETRY  VARIABLES.
CD
CD    ARGUMENTS : 
CD    NP   INTEGER     NUMBER OF GRID POINTS PLUS ONE
CD
CD    ******************************************************************
C
C      COPY  THE  PREVIOUSLY  NEW  GRID  VALUES  TO  BE  USED  AS  THE
C            OLD  GRID.
C
       N = NP - 1
       DO  401  I = 1, NP
       ROH(I) = RNH(I)
       LO(I) = LN(I)
 401   RLO(I) = RLN(I)
       RETURN
C
C
C      ------------------------------------------------------------------
C
       ENTRY CONSRP1(RHO, NP, CSUM)
C
CD     *****************************************************************
CD
CD     CONSRP(RHO, NP, CSUM)
CD     DESCRIPTION : THIS  ENTRY  COMPUTES  THE  OSTENSIBLY  CONSERVED
CD                   SUM.
CD
CD     ARGUMENTS:
CD     RHO     REAL  ARRAY(N)  GRID  PT. VALUES  FOR  CONSERVATION  SUM.
CD     NP      INTEGER         NUMBER OF GRID POINTS PLUS ONE.
CD     CSUM    REAL            VALUE   OF  THE  CONSERVATION  SUM OF RHO
CD
CD     *****************************************************************
CD
C      COMPUTE  THE  OSTENSIBLY  CONSERVED  TOTAL  MASS ( BEWARE YOUR B.C. )
C
       N = NP - 1
       CSUM  =  0.0
       DO  501  I = 1, N
 501   CSUM  =  CSUM  +  LN(I)*RHO(I)
       RETURN
       END


