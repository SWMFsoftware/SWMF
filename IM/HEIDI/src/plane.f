C  plane()  --  calculate flux-tube content in the equatorial plane
C
      subroutine plane(yearIn, dayIn, secIn, Kp, ap, R, dt, n)
	implicit none
	include 'plane.h'
	integer yearIn, dayIn
        real secIn, Kp, ap, R, dt, n(NL,0:*)
     &     , dtH , time, L(0:NL+1), mlt(0:NLT)
     &	   , tau0(NL,0:NLT), n0(0:NL+1), nsw(0:NLT,NL)
     &	   , uL(0:NL+1,0:NLT), ut(0:NLT,0:NL+1), vol(NL)
     &	   , gLong, gLat, mLong0, mLat, update
	integer  IOcount, Tcount   ! k, nk, IOskip, Tskip
C       character*1 ans
        logical first
	data first/.true./, Tcount/0/, IOcount/0/
	save
cv        save first
        integer year, month, day
        real sec
	common /DateTime/ year, month, day, sec

        if (first) then
           year = yearIn
           day = dayIn
           sec = secIn
           time = sec/3600.
           update = time
	   month = .5 + day/30.
	   gLat = 45.                       ! get mLong at 0 mlt and
	   gLong = 270.                     ! at initial time (American sector)
	   call ggmxx(0, gLong, gLat, mLong0, mLat)
*
*  initialize grid
*
	   call initGrid(L, mlt)

	   call ngride1(L(1), NL, L(NL+1), L(0), CARTESIAN)
	   call ogride1(NL)
	   call ngride1(L(1), NL, L(NL+1), L(0), CARTESIAN)
	   call ngridp1(mlt, NLT+1)
	   call ogridp1(NLT+1)
	   call ngridp1(mlt, NLT+1)
*
*  initial data
*
           call initContent(R, day, L, n0)
	   call volume(L, vol)
           call TauFill(R, ap, mLong0, mlt, L, vol, n0, tau0) ! get time scales
           update = time + 3.                  ! change fluxes every () hours
           first = .false.
        end if

c.        if (update .le. time) then
c.           update = time + 3.                  ! change fluxes every () hours
c.           call TauFill(R, ap, mLong0, mlt, L, vol, n0, tau0) ! get time scales
c.        end if

        dtH = dt/2./3600.                             ! also convert to hr.
        call velocity(Kp, L, mlt, uL, ut)

	call advLtime(time, dtH, L, tau0, n0, n, uL)
	call nTOnsw(n, nsw)
	call advTtime(time, dtH, mlt, nsw, ut)

	time = time + dtH

	call advTtime(time, dtH, mlt, nsw, ut)
	call nswTOn(nsw, n)
	call advLtime(time, dtH, L, tau0, n0, n, uL)

	time = time + dtH

	return
      end


*
*  advLtime()  --  advance one time step in L direction
*

      subroutine advLtime(time, dt, L, tau0, n0, n, uL)
	implicit none
	include 'plane.h'
	real time, dt, L(0:*), tau0(NL,0:*), n0(0:*), n(NL,0:*)
     &	   , uL(0:NL+1,0:*), S(NL,0:NLT), null(1), sum(NL)
     &	   , bcB, bcE
	integer i, j

	call IonSource(n0, n, tau0, S)
	bcB = n0(0)/n0(1)
	do j = 0, NLT
	    call veloce1(uL(1,j), NL, uL(NL+1,j), uL(0,j), dt)
	    call sourcq1(NL, dt, 3, null, S(1,j), 0.0, 0.0)
	    do i = 1, NL
	        sum(i) = -4.*n(i,j)*uL(i,j)/L(i)
	    end do
	    call sourcq1(NL, dt, 3, null, sum, 0.0, 0.0)
	    call sourcq1(NL,dt,2,n(1,j),uL(1,j),uL(NL+1,j),uL(0,j))
	    if (uL(NL+1,j) .ge. 0.0) then
	        bcE = n0(NL+1)/n0(NL)
		call etbFCT1(n(1,j), n(1,j), NL, bcE, bcB)
	    else
	        bcE = 0.
		call etbFCT1(n(1,j), n(1,j), NL, bcE, bcB)
		n(NL,j) = 0.
	    end if
	end do
	
	return
      end


*
*  advTtime()  --  advance one time step in azimuthal direction
*

      subroutine advTtime(time, dt, mlt, nsw, ut)
	implicit none
	include 'plane.h'
	real time, dt, mlt(0:NLT), nsw(0:NLT,NL), ut(0:NLT,0:NL+1)
	integer i

	do i = 1, NL
	    call velocp1(ut(0,i), NLT+1, dt)
	    call sourcp1(NLT+1, dt, 2, nsw(0,i), ut(0,i))
	    call prbfct1(nsw(0,i), nsw(0,i), NLT+1)
	end do

	return
      end


*
*  initContent()  --  initialize tube-content
*
      subroutine initContent(R, day, L, n0)
	implicit none
	include 'plane.h'
	real R, L(0:*), n0(0:*), x, fac, dy
     &	   , a, Lc, nLc
	integer day, i   ! ,j

*
*  calculate saturation levels [from, Carpenter and Anderson, 1992]
*
	dy = 2.*PI*(day + 9.)/365.
	fac = 0.15*(cos(dy) - 0.5*cos(2.*dy)) + 0.00127*R - 0.0635
	do i = 0, NL+1
	    x = exp(-(L(i) - 2.)/1.5)
	    n0(i) = 10**(3.9043 - 0.3145*L(i) + fac*x)
	end do

	a = 3.				! change to L^-a dependence
	Lc = 5.2			! at L = 5.2
	x = exp(-(Lc - 2.)/1.5)
	nLc = 10**(3.9043 - 0.3145*Lc + fac*x)
	fac = nLc * Lc**a
	do i = 0, NL+1
	    if (L(i) .gt. Lc) n0(i) = fac * L(i)**(-a)
	end do
*
	return
      end


*
*  volume()  --  calculate normalized flux-tube volume (X B0/r0)
*
      subroutine volume(L, vol)
	implicit none
	include 'plane.h'
	real L(0:*), vol(*), a, x
	integer i

	a = (1. + 250./6371.)
	do i = 1, NL
	    x = a/L(i)
	    vol(i) = (32./35.) * L(i)**4 * sqrt(1. - x)
     &		* (1. + (1./2.)*x + (3./8.)*x**2 + (5./16.)*x**3)
	end do
	
	return
      end

*
*  TauFill()  --  calculate refilling time constants (from ionospheric fluxes)
*
      subroutine TauFill(R, ap, mLong0, mlt, L, vol, n0, tau0)
	
        use ModIoUnit, ONLY : io_unit_new

        implicit none
	include 'plane.h'
	real R, ap, mLong0, mlt(0:*), L(0:*), vol(*), n0(0:*),RE
     &     ,tau0(NL,0:NLT), n(NL,0:*), S(NL,0:*) 
C    &     ,a,F10p7,sqroot,flux10p7,flux,fluxN,fluxS,gLong,gLat,mLat,mLong
	integer i, j
        parameter(RE = 6371.e5)
	real sec
	integer year, month, day
	common /DateTime/ year, month, day, sec

        integer :: iUnitIn! = 32
        
        iUnitIn = io_unit_new()
        
        open(unit=iUnitIn, file='tau0.dat', status='old', err=10
     &		, form='unformatted')
	read(iUnitIn) tau0
        close(iUnitIn)
        return

 10     stop ' error on open to tau0.dat'

 99	return

      entry IonSource(n0, n, tau0, S)
	
	do i = 1, NL
	    do j = 0, NLT
	        S(i,j) = -(n(i,j) - n0(i)) / tau0(i,j)
c...	        if (S(i,j) .ne. 0.)        ! speed up refilling near saturation
c...     &		    S(i,j) = S(i,j) / sqrt( abs(1. - n(i,j)/n0(i)) )
	    end do
	end do

	return
      end

*
*  initGrid()  --  initialize grid (L in units of RE)
*
      subroutine initGrid(L, mlt)
	implicit none
	include 'plane.h'
	real  L(0:*), mlt(0:*)
	integer i

	L(1) = 1.5
	L(0) = L(1) - 0.5*DL
	do i = 2, NL
	    L(i) = L(i-1) + DL
	end do
	L(NL+1) = L(NL) + 0.5*DL
	
	mlt(0) = 0.0
	do i = 1, NLT
	    mlt(i) = mlt(i-1) + DLT
	end do

	return
      end

*
*  fluxInf()  --  calculate maximum ionospheric flux at 1000 km
*
      subroutine fluxInf(year,day,month,sec,gLat,gLong,f10p7,ap,flux)
      implicit none
      include 'plane.h'
      integer year, day, month
      real sec, gLat, gLong, f10p7, ap, flux
     &     , altr, alt0, glt, glng, sLT, Tn, Ti, Te, nO, nH, nOp
     &     , rf, r0, g, HOp, Ho, Hc, Hp, d(8), exoT(2), Rs, magbr, Rsun
      real RE, KB, MHP, MOP
      parameter(RE=6371.0e5, KB=1.38e-16, MHP=1.67e-24, MOP=2.67e-23)
      integer iyd, iMass, jMag
      data jMag/0/

      iyd = 1000*(year - 1900) + day
*                                        altr 50 km less than Richards and Torr
      altr = 450. + (750. - 450.)*(f10p7 - 70.)/(230. - 70.)
      glt = gLat					! degrees
      glng = gLong					! degrees
      sLT = amod(sec/3600. + glng/15., 24.)
      Rs = Rsun(f10p7)
      
      iMass = 0
      call gts5(iyd, sec, altr, glt, glng, sLT
     &     , f10p7, f10p7, ap, iMass, d, exoT)
      Tn = exoT(2)
      iMass = 16
      call gts5(iyd, sec, altr, glt, glng, sLT
     &     , f10p7, f10p7, ap, iMass, d, exoT)
      nO = d(2)                                                     ! O
      call iri(altr, jMag, glt, glng, Rs, month, sLT
     &       , nOp, Ti, Te, magbr)                                  ! Ti, Te
      
      rf = RE + altr*1.e5
      g = 980.*(RE/rf)**2
      HOp = KB*(Ti + Te)/(MOP*g)
      Ho  = KB*Tn/(MOP*g)
      Hc = HOp*Ho / (HOp - Ho)
      Hp = ( HOp*Ho / (HOp + Ho) ) / 1.e5

      alt0 = altr - Hp*log( (2.e19/(nOp*nO)) * (Ti/Hc)**2 )
      iMass = 0
      call gts5(iyd, sec, alt0, glt, glng, sLT
     &     , f10p7, f10p7, ap, iMass, d, exoT)
      Tn = exoT(2)
      iMass = 1
      call gts5(iyd, sec, alt0, glt, glng, sLT
     &     , f10p7, f10p7, ap, iMass, d, exoT)
      nH = d(7)                                                    ! H
      call iri(alt0, jMag, glt, glng, Rs, month, sLT
     &       , nOp, Ti, Te, magbr)                                 ! O^+
      r0 = RE + alt0*1.e5
      g = 980.*(RE/r0)**2
      HOp = KB*(Ti+Te)/(MOP*g)

      flux = 2.5e-11*sqrt(Tn)*nH*nOp*HOp                           ! at alt0
      flux = flux * ( (6371.+alt0)/(6371.+1000.) )**3              ! at 1000 km

      return
      end


*
*  Rsun()  --  calculate sunspot number Rs from F_10.7
*
      real function Rsun(f10p7)
      implicit none
      real f10p7, a, b
      data a/.9261497/, b/57.36373/

      Rsun = (f10p7 - b)/a

      return
      end

*
*  flux10p7()  --  calculate F_10.7 from sunspot number Rs
*
      real function flux10p7(Rs)
      implicit none
      real Rs, a, b
      data a/.9261497/, b/57.36373/

      flux10p7 = a*Rs + b

      return
      end

C
C                     
C************************************************************
C*************** EARTH MAGNETIC FIELD ***********************
C************************************************************
C
C
      SUBROUTINE GGMXX(ART,LONG,LATI,MLONG,MLAT)
C CALCULATES GEOMAGNETIC LONGITUDE (MLONG) AND LATITUDE (MLAT)
C FROM GEOGRAFIC LONGITUDE (LONG) AND LATITUDE (LATI) FOR ART=0
C AND REVERSE FOR ART=1. ALL ANGLES IN DEGREE.
C LATITUDE:-90 TO 90. LONGITUDE:0 TO 360 EAST.
      INTEGER ART
      REAL MLONG,MLAT,LONG,LATI

      faktor = ATAN(1.0)*4./180.
      ZPI=FAKTOR*360.                    
      CBG=11.4*FAKTOR                              
      CI=COS(CBG)     
      SI=SIN(CBG)
      IF(ART.EQ.0) GOTO 10                         
      CBM=COS(MLAT*FAKTOR)                           
      SBM=SIN(MLAT*FAKTOR)                           
      CLM=COS(MLONG*FAKTOR)                          
      SLM=SIN(MLONG*FAKTOR)
      SBG=SBM*CI-CBM*CLM*SI                     
      LATI=ASIN(SBG)   
      CBG=COS(LATI)     
      SLG=(CBM*SLM)/CBG  
      CLG=(SBM*SI+CBM*CLM*CI)/CBG
        IF(ABS(CLG).GT.1.) CLG=SIGN(1.,CLG)                  
      LONG=ACOS(CLG)  
      IF(SLG.LT.0.0) LONG=ZPI-ACOS(CLG)            
      LATI=LATI/FAKTOR    
      LONG=LONG/FAKTOR  
      LONG=LONG-69.8    
      IF(LONG.LT.0.0) LONG=LONG+360.0                 
      RETURN          
10    YLG=LONG+69.8    
      CBG=COS(LATI*FAKTOR)                           
      SBG=SIN(LATI*FAKTOR)                           
      CLG=COS(YLG*FAKTOR)                          
      SLG=SIN(YLG*FAKTOR)                          
      SBM=SBG*CI+CBG*CLG*SI                        
      MLAT=ASIN(SBM)   
      CBM=COS(MLAT)     
      SLM=(CBG*SLG)/CBM                            
      CLM=(-SBG*SI+CBG*CLG*CI)/CBM
        IF(CLM.GT.1.) CLM=1.                 
      MLONG=ACOS(CLM)  
      IF(SLM.LT..0) MLONG=ZPI-ACOS(CLM)             
      MLAT=MLAT/FAKTOR    
      MLONG=MLONG/FAKTOR  
      RETURN          
      END             

*
*  nTOnsw()  --  copy n(i,j) to nsw(j,i)
*
      subroutine nTOnsw(n, nsw)
	implicit none
	include 'plane.h'
	real n(NL,0:*), nsw(0:NLT,*)
	integer i, j
	
	do i = 1, NL
	    do j = 0, NLT
		nsw(j,i) = n(i,j)
	    end do
	end do

	return
      end
      
*
*  nswTOn()  --  copy nsw(j,i) to n(i,j)
*
      subroutine nswTOn(nsw, n)
	implicit none
	include 'plane.h'
	real nsw(0:NLT,*), n(NL,0:*)
	integer i, j
	
	do i = 1, NL
	    do j = 0, NLT
		n(i,j) = nsw(j,i)
	    end do
	end do

	return
      end

*
*  velocity()  --  get velocity from Volland-Stern model (assume GAMMA=2)
*
      subroutine velocity(Kp, L, mlt, uL, ut)
	implicit none
	include 'plane.h'
	real Kp, L(0:*), mlt(0:*), uL(0:NL+1,0:*), ut(0:NLT,0:*)
     &	   , A, B, phi, vr, vp, facL, facT, RE
	integer i, j
	parameter(RE = 6371.e3)
        logical corotation
        data corotation/.true./

	A = 7.05e-6 / (1. - 0.159*Kp + 0.0093*Kp*Kp)**3
	facL = 3600. / RE	                        ! per hour
	facT = 3600. * (24.0/(2.*PI)) / RE	        ! per hour

	do i = 0, NL+1
            B = 0.3e-4 / L(i)**3                         ! weber/m^2
            vr = -(A*L(i)/B) * facL
            vp = (2.0*A/B) * facT
	    do j = 0, NLT
                phi = mlt(j)*15.*DTR
		uL(i,j) = vr*cos(phi)		         ! dL/dt
		ut(j,i) = vp*sin(phi)		         ! d(mlt)/dt
                if (corotation) ut(j,i) = ut(j,i) + 1.
	    end do
	end do

	return
      end
C
C
C Routine CARPENTER is a plasmaspheric model of electron densities established
C by Carpenter and Anderson [JGR,1992]
C
        SUBROUTINE CARPENTER(NE,KPMAX,DAY,R)
        INCLUDE 'plane.h'
        EXTERNAL FUNC
        REAL NE(NL,0:NLT),L,KPMAX,LPPI,LPPO
	INTEGER DAY
        COMMON /PARAS1/LPPI,X,C1,C2,T1,IDAY,XR
	IDAY=DAY
	XR=R
        LPPI=5.6-0.46*KPMAX
        DO 3 J=0,NLT
          T1=24./NLT*J
          IF(T1.LT.6.) THEN
            X=0.1
            C1=5800.
            C2=300.
          ENDIF
          IF(T1.GE.6.0.AND.T1.LE.19.) THEN
            X=0.1+0.011*(T1-6.)
            C1=-800.
            C2=1400.
          ENDIF
          IF(T1.GT.19..AND.T1.LT.20.) THEN    ! Day-night transition
            X=0.1-0.143*(T1-20.)
            C1=428600.
            C2=-21200.
          ENDIF
          IF(T1.GE.20.) THEN
            X=0.1
            C1=-1400.
            C2=300.
          ENDIF
C
C Find Lppo
        LPPO=RTBIS(FUNC,LPPI,10.,0.025,IER)
        DO 3 I=1,NL
          L=1.5+(I-1)*0.25
          IF(L.LE.LPPI) NE(I,J)=SPS(L,DAY,R)
          IF(L.GT.LPPI.AND.L.LE.LPPO) NE(I,J)=PS(L,LPPI,X,DAY,R)
          IF(L.GT.LPPO) THEN
            XN=PS(LPPO,LPPI,X,DAY,R)*(L/LPPO)**(-4.5)
            NE(I,J)=XN+1.-EXP((2.-L)/10.)
          ENDIF
3       CONTINUE
C
        RETURN
        END

C
C ************************** Functions ***********************************
C
C Function SPS calculate the electron density in the saturated plasmsphere
C segment 
        FUNCTION SPS(L,DAY,R)
        INCLUDE 'plane.h'
        REAL L
        INTEGER DAY
        X1=-0.3145*L+3.9043
        X2=0.15*(COS(2.*PI*(DAY+9)/365.)-0.5*COS(4.*PI*(DAY+9)/365.))
        X3=0.00127*R-0.0635
        XT=X1+(X2+X3)*EXP(-(L-2.)/1.5)
        SPS=10**XT
        RETURN
        END
C
C Function PS calculate the electron density at the plasmapause segment
        FUNCTION PS(L,LPPI,X,DAY,R)
        REAL L,LPPI
	INTEGER DAY
        PS=SPS(LPPI,DAY,R)*10.**((LPPI-L)/X)
        RETURN
        END
C
C Function FUNC is used to locate Lppo
        FUNCTION FUNC(L)
        REAL L,LPPI
	INTEGER DAY
        COMMON /PARAS1/LPPI,X,C1,C2,T1,DAY,R
        Y=(C1+C2*T1)*L**(-4.5)+1.-EXP((2.-L)/10.)
        FUNC=Y-PS(L,LPPI,X,DAY,R)
        RETURN
        END
C
      FUNCTION RTBIS(FUNC,X1,X2,XACC,IER)
	external func	! vania
c
c    JMAX...max iterations
c
      PARAMETER (JMAX=40)

      COMMON /PARAS1/LPPI,X,C1,C2,T1,DAY,R
      COMMON /CFUNY/Z1,B1,SI1
      COMMON /Cyy/yy

      IER = 0
      FMID=FUNC(X2)
      F=FUNC(X1)
      IF(FMID.EQ.0) THEN
        RTBIS=X2
	RETURN
      ENDIF	
      IF(F.EQ.0) THEN
        RTBIS=X1
	RETURN
      ENDIF	
c
c    not bracketed
c
      IF(F*FMID.GT.0.) THEN
        IER = 1
        RETURN
      ENDIF
c
c    orient F>0 for X+DX
c
      IF(F.LT.0.)THEN
        RTBIS=X1
        DX=X2-X1
      ELSE
        RTBIS=X2
        DX=X1-X2
      ENDIF
c
c    bisection loop
c
      DO 11 J=1,JMAX
        DX=DX*.5
        XMID=RTBIS+DX
        FMID=FUNC(XMID)
        IF(FMID.LE.0.)RTBIS=XMID
        IF(ABS(DX).LT.XACC .OR. FMID.EQ.0.) RETURN
11    CONTINUE
      IER = 2
      RETURN
      END
