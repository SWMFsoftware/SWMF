      SUBROUTINE CLFME2W
      use ModCommonVariables
      REAL C1(MaxGrid),C2(MaxGrid),C3(MaxGrid),C4(MaxGrid),
     $     D1(MaxGrid),D2(MaxGrid),
     $     D3(MaxGrid),D4(MaxGrid)
      REAL YL(MaxGrid),YK(MaxGrid)
C     
C
      C2(1)=H0*(TCONE(2)-TCSFE)
      C2(NDIM)=H0*(TCBGE-TCONE(NDIM2))
      DO 8715 K=2,NDIM2
 8715    C2(K)=H0*(TCONE(K+1)-TCONE(K-1))
      DO 8720 K=1,NDIM
         C1(K)=TCONE(K)/DELECT(K)
         C2(K)=C2(K)/DELECT(K)
 8720 CONTINUE
      C3(1)=H0*(UELECT(2)-USURFE)
      C3(NDIM)=H0*(UBGNDE-UELECT(NDIM2))
      DO 8725 K=2,NDIM2
 8725    C3(K)=H0*(UELECT(K+1)-UELECT(K-1))
      DO 8730 K=1,NDIM
      C2(K)=C2(K)-UELECT(K)
      C2(K)=C2(K)+DAREA(K)*C1(K)
 8730 C3(K)=-GMIN1*(DAREA(K)*UELECT(K)+C3(K))
C
C
      XHLP=EXP(-(TIME-300.)**2/2./150./150.)
      DO 8735 K=1,NDIM
         C3(K)=C3(K)-ADMSE(K)/DELECT(K)
         C4(K)=HLPE0*(ECLSNE(K)+XHLP*QELECT(K))/DELECT(K)
 8735 CONTINUE
C
C
C      DO 8738 K=1,NCL
C      XX1=H3*C1(K)
C      XX2=H4*C2(K)
C      XX3=H2*C1(K)
C      XX4=0.5*C3(K)-XX3
C      D1(K)=-XX1-XX2
C      D2(K)=H1E1-XX4
C      D3(K)=-XX1+XX2
C8738  D4(K)=H1E1+XX4
C      DO 8740 K=NCL+1,NDIM
      DO 8740 K=1,NDIM
      XX1=H3*C1(K)
      XX2=H4*C2(K)
      XX3=H2*C1(K)
      XX4=0.5*C3(K)-XX3
      D1(K)=-XX1-XX2
      D2(K)=H1E2-XX4
      D3(K)=-XX1+XX2
8740  D4(K)=H1E2+XX4
      D4(1)=C4(1)-D3(1)*TSURFE+D4(1)*TELECT(1)-D1(1)*TELECT(2)
      D4(NDIM)=C4(NDIM)-D3(NDIM)*TELECT(NDIM2)+
     $D4(NDIM)*TELECT(NDIM)-D1(NDIM)*TBGNDE
      DO 8750 K=2,NDIM2
         D4(K)=C4(K)-D3(K)*TELECT(K-1)+D4(K)*TELECT(K)-
     $        D1(K)*TELECT(K+1)
 8750 CONTINUE
      YL(1)=D1(1)/D2(1)
      YK(1)=(D4(1)-D3(1)*TSURFE)/D2(1)
CALEX is the syntax for these do loops correct?      
      DO 8760 K=2,NDIM
      XX1=D2(K)-D3(K)*YL(K-1)
      YK(K)=(D4(K)-D3(K)*YK(K-1))/XX1
8760  YL(K)=D1(K)/XX1
      TELECT(NDIM)=YK(NDIM)-YL(NDIM)*TBGNDE
      PELECT(NDIM)=RGASE*DELECT(NDIM)*TELECT(NDIM)
      DO 8770 K=NDIM2,1,-1
      TELECT(K)=YK(K)-YL(K)*TELECT(K+1)
8770  PELECT(K)=RGASE*DELECT(K)*TELECT(K)
      RETURN
      END
