CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE BOXNUM (E1, E2, M1, M2, R1, R2, NBINS, DEL, ENER)
C
C This subroutine finds the box numbers corresponding to
C energies E1 and E2, and calls them M1 and M2
C 
C R1 is the upper edge of the lower box, R2 is the lower edge of the
C upper box.
C
C
      DIMENSION DEL(NBINS), ENER(NBINS)
C
      DO 100 I=1,NBINS
      IF (E1.GE.ENER(I)-DEL(I)/2. .AND. E1.LT.ENER(I)+DEL(I)/2.)
     >  GOTO 200
  100 CONTINUE
C
      IF (E1 .LT. 0.25) THEN
        M1 = 1
        R1 = 0.75
      ELSE
        M1 = NBINS+1
        R1 = ENER(I) + DEL(I) * 2.
      ENDIF
      GOTO 300
C
  200 M1 = I
      R1 = ENER(I) + DEL(I) / 2.
C
  300 DO 400 I=1,NBINS
      IF (E2.GE.ENER(I)-DEL(I)/2. .AND. E2.LT.ENER(I)+DEL(I)/2.)
     >  GOTO 500
  400 CONTINUE
C
      IF (E2 .LT. 0.25) THEN
        M2 = 1
        R2 = 0.25
      ELSE
        M2 = NBINS+1
        R2 = ENER(I) + DEL(I)
      ENDIF
      GOTO 600
C
  500 M2 = I
      R2 = ENER(I) - DEL(I) / 2.
C
  600 RETURN 
C
      END
