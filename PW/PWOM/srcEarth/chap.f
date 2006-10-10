 
C
C
C
C
      FUNCTION CHAP (CHI, Z, T, I)
      DIMENSION AM(3)
      DATA AM/16., 32., 28./, PI/3.1415926535/, RE/6.37E8/, G/978.1/
      GR=G*(RE/(RE+Z))**2 
      HN=1.38E-16*T/(AM(I)*1.662E-24*GR)
      HG=(RE+Z)/HN 
      HF=0.5*HG*(COS(CHI)**2) 
      SQHF=SQRT(HF) 
      CHAP=SQRT(0.5*PI*HG)*SPERFC(SQHF) 
      RETURN
      END
