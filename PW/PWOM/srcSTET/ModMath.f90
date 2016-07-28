!*  This is a collection of math subroutines based on versions origionally 
!*  coded by Michael Liemohn (1994). Updated and put in this module by 
!*  Alex Glocer (2014)

Module ModMath
  private !except
  
  public :: midpnt_int
  public :: erf
  public :: G

contains

  !----------------------------------------------------------------------*
  !*  SUBROUTINE midpnt_int
  !*  This subroutine integrates a function array over an abscissa array
  !*  using the midpoint integration technique from the Second Euler-
  !*  MacLaurin Summation Formula:
  !*        INT[f(x),x=a..b] = SUM[dx_j*f_j+1/2,j=0,n-1]
  !*        for n points along [a,b].
  !*  Essentially, compute the function at the midpoint of the interval
  !*  of each step along x.
  !*  The index ii determines what the x array is: abscissa values (1) or
  !*  abscissa intervals (2). For ii=1, f is assumed to be on the values,
  !*  for ii=2, f is assumed to be on the midpoints (already averaged).
  !*  This distinction is for the PA/s (1) and the energy (2) integrals.
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
       
!==============================================================================
  REAL FUNCTION erf(y)
    !*  specifications for arguments
    REAL y
    !*  specifications for local variables
    INTEGER isw,i
    DIMENSION p(3),q(2),p1(5),q1(4),p2(3),q2(2)
    REAL p,q,p1,q1,p2,q2,xmin,xlarge,sqrpi,x,res,xsq,xnum,xden,xi
    !*  coefficients for 0.0 .le. y .lt. 0.477
    DATA p(1)/.3166529/, p(2)/1.722276/, p(3)/21.38533/
    DATA q(1)/7.843746/, q(2)/18.95226/
    !*  coefficients for .477 .le. y .le. 4.0
    DATA p1(1)/.5631696/, p1(2)/3.031799/, p1(3)/6.865018/, &
         p1(4)/7.373888/, p1(5)/4.318779e-5/
    DATA q1(1)/5.354217/, q1(2)/12.79553/, q1(3)/15.18491/, &
         q1(4)/7.373961/
    !*  coefficients for 4.0 .lt. y
    DATA p2(1)/-5.168823e-2/, p2(2)/-.1960690/, p2(3)/-4.257996e-2/
    DATA q2(1)/.9214524/, q2(2)/.1509421/
    !*  constants
    DATA xmin/1.0e-5/,xlarge/4.1875e0/
    DATA sqrpi/.5641896/
    !*  first executable statement
    x = y
    isw = 1
    IF (x.LT.0.0e0) THEN
       isw = -1
       x = -x
    ELSE IF (x.LT.xmin) THEN
       res = x*p(3)/q(2)
    ELSE IF (x.LT..477e0) THEN
       xsq = x*x
       xnum = (p(1)*xsq+p(2))*xsq+p(3)
       xden = (xsq+q(1))*xsq+q(2)
       res = x*xnum/xden
    ELSE IF (x.GE.xlarge) THEN
       res = 1.0e0
    ELSE
       IF (x.LE.4.0e0) THEN
          xsq = x*x
          xnum = p1(5)*x+p1(1)
          xden = x+q1(1)
          DO  i=2,4
             xnum = xnum*x+p1(i)
             xden = xden*x+q1(i)
          end do
          res = xnum/xden
       ELSE
          xsq = x*x
          xi = 1.0e0/xsq
          xnum = (p2(1)*xi+p2(2))*xi+p2(3)
          xden = (xi+q2(1))*xi+q2(2)
          res = (sqrpi+xi*xnum/xden)/x
       END IF
       res = res*EXP(-xsq)
       res = 1.0e0-res
    END IF
    IF (isw.EQ.-1) res = -res
    erf = res
    RETURN
  END FUNCTION erf

  !* ------------------------------------------------------------------ **
  !*  The G function
  FUNCTION G(X)
    REAL G,X,pi
    DATA pi/3.1415926/
    X2=X*X
    !write(*,*) 'X, X2',X, X2
    G=.5*erf(X)/X2 - EXP(-X2)/SQRT(pi)/X
    RETURN
  END FUNCTION G


end Module ModMath
