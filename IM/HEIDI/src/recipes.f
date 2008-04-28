      SUBROUTINE LINTP(XX,YY,N,X,Y,IER)
CD
      DIMENSION XX(N),YY(N)
      IER = 0
c
c    initialize lower and upper values
c
      JL=1
      JU=N
c
c    if not done compute a midpoint
c
10    IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
c
c    now replace lower or upper limit
c
        IF((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM)))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
c
c    try again
c
      GO TO 10
      ENDIF
c
c    this is J
c
      J=JL	! If X.LE.XX(1) then J=1
c		  If X.GT.X(J).AND.X.LE.X(J+1) then J=J
c		  If X.GT.X(N) then J=N-1	
      D=XX(J+1)-XX(J)
      Y=(YY(J)*(XX(J+1)-X)+YY(J+1)*(X-XX(J)))/D
      RETURN
      END
C
C
      SUBROUTINE LINTP2(X1A,X2A,YA,M,N,X1,X2,Y,IER)
c
c    NN...maximum expected dimension of N or M
c
      PARAMETER (NN=500)
      DIMENSION X1A(M),X2A(N),YA(M,N),YTMP(NN),YYTMP(NN)
      IER = 0
c
c    do M evaluations of row constructed using 1-dimensional evaluator LINTP
c
      DO 12 J=1,M
        DO 11 K=1,N
          YTMP(K)=YA(J,K)
11      CONTINUE
        CALL LINTP(X2A,YTMP,N,X2,YYTMP(J),IER)
        IER = 10 * IER
        IF (IER.EQ.10) RETURN
12    CONTINUE
c
c    evaluate it
c
      CALL LINTP(X1A,YYTMP,M,X1,Y,IER)
      IER = IER * 10
      RETURN
      END
C
C
C
      FUNCTION GAMMLN(XX,IER)
c
c    double precision for recurrences
c
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      IER = 0
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END
C
C
C
      SUBROUTINE GSER(GAMSER,A,X,GLN,IER)
c
c    ITMAX...max iterations
c    EPS...small number
c
      PARAMETER (ITMAX=100,EPS=3.E-7)
      IER = 0
      GLN=GAMMLN(A,IER)
      IF(X.LE.0.)THEN
        IF(X.LT.0.)PAUSE
        GAMSER=0.
        RETURN
      ENDIF
      AP=A
      SUM=1./A
      DEL=SUM
      DO 11 N=1,ITMAX
        AP=AP+1.
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(ABS(DEL).LT.ABS(SUM)*EPS)GO TO 1
11    CONTINUE
c
c    too many iterations
c
      IER = 1
      RETURN
1     GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
      RETURN
      END
C
C
C
      SUBROUTINE GCF(GAMMCF,A,X,GLN,IER)
c
c    ITMAX...max iterations to preform
c    EPS...a small number
c
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A,IER)
      IER = 0
c
c    previous value to check for convergence
c
      GOLD=0.
c
c    setting up to evaluate continuous fraction
c
      A0=1.
      A1=X
      B0=0.
      B1=1.
c
c    renormalized factor preventing overflow
c
      FAC=1.
      DO 11 N=1,ITMAX
        AN=FLOAT(N)
        ANA=AN-A
c
c    one step of the recurrence
c
        A0=(A1+A0*ANA)*FAC
        B0=(B1+B0*ANA)*FAC
c
c    next step
c
        ANF=AN*FAC
        A1=X*A0+ANF*A1
        B1=X*B0+ANF*B1
c
c    time to renormalize ?
c
        IF(A1.NE.0.)THEN
          FAC=1./A1
          G=B1*FAC
c
c    converged ?
c
          IF(ABS((G-GOLD)/G).LT.EPS)GO TO 1
          GOLD=G
        ENDIF
11    CONTINUE
c
c    error
      IER = 1
      RETURN
1     GAMMCF=EXP(-X+A*ALOG(X)-GLN)*G
      RETURN
      END
C
C
C
      FUNCTION GAMMP(A,X,IER)
c
      IER = 0
      IF(X.LT.0..OR.A.LE.0.)PAUSE
c
c    use series representation
c
      IF(X.LT.A+1.)THEN
        CALL GSER(GAMMP,A,X,GLN,IER)
        IER = IER * 20
        IF (IER.EQ.20) RETURN
c
c    use continued fraction representation
c
      ELSE
        CALL GCF(GAMMCF,A,X,GLN,IER)
        IER = 10 * IER
        IF (IER.EQ.10) RETURN
c
c    take its components
c
        GAMMP=1.-GAMMCF
      ENDIF
      RETURN
      END
C
C
C
      FUNCTION ERF(X,IER)
c
      IER = 0
      IF(X.LT.0.)THEN
        ERF=-GAMMP(.5,X**2,IER)
        IF (IER.NE.0) RETURN
      ELSE
        ERF=GAMMP(.5,X**2,IER)
        IF (IER.NE.0) RETURN
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE TRIDAG(A,B,C,R,U,N,IER)
CD
      PARAMETER (NMAX=100)
      DIMENSION GAM(NMAX),A(N),B(N),C(N),R(N),U(N)
c
c    problem can be simplified to N-1
c
      IF(B(1).EQ.0.)THEN
        IER = 1
        RETURN
      ENDIF
      IER = 0
      BET=B(1)
      U(1)=R(1)/BET
c
c    decomposition and forward substitution
c
      DO 11 J=2,N
        GAM(J)=C(J-1)/BET
        BET=B(J)-A(J)*GAM(J)
c
c    algotithm fails
c
        IF(BET.EQ.0.)THEN
          IER = 2
          RETURN
        ENDIF
        U(J)=(R(J)-A(J)*U(J-1))/BET
11    CONTINUE
c
c    back substitution
c
      DO 12 J=N-1,1,-1
        U(J)=U(J)-GAM(J+1)*U(J+1)
12    CONTINUE
      RETURN
      END
C
C
C
      subroutine bessel2(jn,arg1,bs)
      external bessj,bessj0,bessj1
cc this version of bessel2 returns  bs(1),n=jn-1;bs(2),n=jn;bs(3),n=jn+1.
      real arg,arg1,bs(3),bessj
      integer jj(3)
      arg=arg1
      if(arg.lt.0)arg=-arg
      jj(1)=jn-1
      jj(2)=jn
      jj(3)=jn+1
      do 5 in=1,3
      ij=jj(in)
      if(ij.lt.0)ij=-ij
      if(ij.ne.0)go to 6
      bs(in)=bessj0(arg)
      go to 5
    6 if(ij.ne.1)go to 7
      bs(in)=bessj1(arg)
      go to 5
    7 bs(in)=bessj(ij,arg)
    5  continue
      if(arg1.ge.0.)go to 16
      do 15 in=1,3
      ij=jj(in)
      if(ij.lt.0)ij=-ij
      ie=ij
      if(ie.eq.0)go to 15
      efac=1.0
      do 17 ii=1,ie
      efac=efac*(-1)
   17 continue
      bs(in)=efac*bs(in)
   15 continue
   16 continue
cc  check jj for negative values and calculate negative harmonics
      do 18 in=1,3
      ij=jj(in)
      if(ij.ge.0)go to 18
      ie=-ij
      efac=1.0
      do 19 ii=1,ie
      efac=efac*(-1)
   19 continue
      bs(in)=efac*bs(in)
   18 continue
      return
      end
      real function bessj(n,x)
      parameter (iacc=40,bigno=1.e10,bigni=1.e-10)
      if(n.lt.2) then
      write(59,777)
 777  format(1x,'bad argument n in bessj-exiting routine ')
      end if
      tox=2./x
      if(x.gt.float(n))then
        bjm=bessj0(x)
        bj=bessj1(x)
        do 11 j=1,n-1
          bjp=j*tox*bj-bjm
          bjm=bj
          bj=bjp
11      continue
        bessj=bj
      else
        m=2*((n+int(sqrt(float(iacc*n))))/2)
        bessj=0.
        jsum=0
        sum=0.
        bjp=0.
        bj=1.
        do 12 j=m,1,-1
          bjm=j*tox*bj-bjp
          bjp=bj
          bj=bjm
          if(abs(bj).gt.bigno)then
            bj=bj*bigni
            bjp=bjp*bigni
            bessj=bessj*bigni
            sum=sum*bigni
          endif
          if(jsum.ne.0)sum=sum+bj
          jsum=1-jsum
          if(j.eq.n)bessj=bjp
12      continue
        sum=2.*sum-bj
        bessj=bessj/sum
      endif
      return
      end
C
C
C
      real function bessj0(x)
      real y,p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     *    s1,s2,s3,s4,s5,s6
      data p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,
     *    -.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-
     *1,
     *    .1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
      data r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,651619640.7d
     *0,
     *    -11214424.18d0,77392.33017d0,-184.9052456d0/,
     *    s1,s2,s3,s4,s5,s6/57568490411.d0,1029532985.d0,
     *    9494680.718d0,59272.64853d0,267.8532712d0,1.d0/
      if(abs(x).lt.8.)then
        y=x**2
        bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))
     *      /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8./ax
        y=z**2
        xx=ax-.785398164
        bessj0=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y
     *      *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      return
      end
C
C
C
      real function bessj1(x)
      real y,p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     *    s1,s2,s3,s4,s5,s6
      data r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,242396853.1d0
     *,
     *    -2972611.439d0,15704.48260d0,-30.16036606d0/,
     *    s1,s2,s3,s4,s5,s6/144725228442.d0,2300535178.d0,
     *    18583304.74d0,99447.43394d0,376.9991397d0,1.d0/
      data p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,.2457520174d-5
     *,
     *    -.240337019d-6/, q1,q2,q3,q4,q5/.04687499995d0,-.2002690873d-3
     *,
     *    .8449199096d-5,-.88228987d-6,.105787412d-6/
      if(abs(x).lt.8.)then
        y=x**2
        bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))
     *      /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8./ax
        y=z**2
        xx=ax-2.356194491
        bessj1=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y
     *      *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
     *      *sign(1.,x)
      endif
      return
      end
C
C
C
