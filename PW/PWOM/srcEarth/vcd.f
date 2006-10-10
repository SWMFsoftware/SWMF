C
C
C
C
      SUBROUTINE VCD(ZZ,ZMAJ,ZVCD,JMAX,NMAJ)
      DIMENSION ZZ(JMAX), ZMAJ(NMAJ,JMAX), ZVCD(NMAJ,JMAX)
C
      DO 200 I=1,NMAJ
      ZVCD(I,JMAX) =   ZMAJ(I,JMAX)
     >               * (ZZ(JMAX)-ZZ(JMAX-1))
     >               / ALOG(ZMAJ(I,JMAX-1)/ZMAJ(I,JMAX))
      DO 200 J=JMAX-1,1,-1
      RAT = ZMAJ(I,J+1) / ZMAJ(I,J)
      ZVCD(I,J) =   ZVCD(I,J+1)
     >            + ZMAJ(I,J) * (ZZ(J)-ZZ(J+1)) / ALOG(RAT) * (1.-RAT)
C       PRINT*,RAT,ZMAJ(I,J+1),ZVCD(I,J)
 
  200 CONTINUE
      RETURN
      END
