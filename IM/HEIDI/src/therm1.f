c      INCLUDE 'plane.f'		! Craig's thermal density model

      PROGRAM Global_thermal_densities
**  Driver program for Craig's plane.f subroutine
      
      use ModIoUnit, ONLY : io_unit_new
      use ModHeidiMain, ONLY: Re

      INCLUDE 'numv.h'
      REAL T,KP,AP,R,DT,NECR(NL,0:NLT),KPN,KPO
        DIMENSION DAYR(48),RKPH(48),F107R(48),APR(48),RSUNR(48)
      INTEGER YEAR,DAY,I,NSTEP,NPR
      CHARACTER*5 NAME
      CHARACTER*80 HEADER

      integer :: iUnitKp

	print *, 'Starting therm1.f'
      DT=20.
      DKP=-DT/3./3600.	! Kp decreases by 1 every 3 hours
      !RE=6.378E6
	IKP=1 ! IKP=0 means no KP update from IC, else update from file
      NAME='d0301'
c      YEAR=1999	!1998	!1991  !1998	!1997	!1998	!1997	!1999
c      DAY=294	!267	!155   !291	!134	!29	!282	!294
c      R=97.	!117.	!179.  !105.	!13.	!47.	!20.	!97.
c      AP=18.	!18.	!7.    !5.	!3.	!5.	!15.	!18.
c      KP=3.3	!3.3	!2.    !1.3	!0.7	!1.3	!3.0	!3.3
      YEAR=2001	!2000	!2000	!2001
      DAY=89	!97	!196	!89
      R=231.	!108.	!164.	!231.
      AP=9.	!6.	!15.	!9.
      KP=2.3	!1.7	!3.0	!2.3
c      TMAX=172800.	! 172800 = 2 days of run time
      TMAX=25.*3600.
      TST=0.
      NST=NINT(TST/DT/2.)+1
      NSTEP=TMAX/DT/2.
      TPR=14400.	! 14400 = every 4 hours
      NPR=TPR/DT/2.
      T=0.
      NKP=NINT(10800./DT/2.)
      I2=(NST-1)/NKP +1
	print *, 'Done with setup'
      IF (IKP.NE.0) THEN
c.......Read Kp history of the modeled storm
         
         iUnitKp = io_unit_new()
         OPEN(iUnitKp,FILE=NAME//'_kp.in',STATUS='OLD')
         READ(iUnitKp,10) HEADER
10       FORMAT(A80)
	PRINT *, NSTEP/NKP+2,INT(NSTEP/NKP)+2
         DO I=1,INT(NSTEP/NKP)+2
          READ(iUnitKp,*) DAYR(I),DUT,RKPH(I),F107R(I),APR(I),RSUNR(I)
         ENDDO
         CLOSE(iUnitKp)
	KPN=KP
	PRINT *, I-1,NKP,KPN,NSTEP
      END IF
      CALL GETDENS(NECR)
c      CALL WRESULT(T,DT,KP,NECR)
      DO 100 I=NST,NSTEP
	IF (IKP.NE.0) THEN
        IF (MOD(I-1,NKP).EQ.0 .OR. I.EQ.NST) THEN
          KPO=KPN
          TOL=T
          DAY=DAYR(I2)
          AP=APR(I2)
          F107=F107R(I2)
          R=RSUNR(I2)
          I2=I2+1
          KPN=RKPH(I2)
	PRINT *, I,I2,T,KPO,KPN,DAY,AP,F107,R
        END IF
        KP=KPO+(KPN-KPO)*(T-TOL)/10800.
	ELSE
	IF (MOD(T+DT,86400.).LE.2.*DT) PRINT *, I,T,KP
	END IF
        A=7.05E-6/(1.-0.159*KP+0.0093*KP**2)**3/RE
	CALL PLANE(YEAR,DAY,T,KP,AP,R,2.*DT,NECR)
c	IF (MOD(I,NPR).EQ.0) CALL WRESULT(T,DT,KP,NECR)
C  Comented out the print because we only need the final output
	T=T+DT*2.
100   CONTINUE
	PRINT *, T,KP,KPO,KPN,TOL
      CALL WRESULT(T,0.,KP,NECR)

      call CON_stop('ERROR in therm1.f')
      END

**--------------------------------------------------------------------**
**  Subroutine GETDENS
**  Reads in the initial densities
      SUBROUTINE GETDENS(NECR)
      
      use ModIoUnit, ONLY : UNITTMP_

      INCLUDE 'numv.h'
      REAL NECR(NL,0:NLT)
      CHARACTER HEADER*80
      
      OPEN (UNITTMP_,FILE='ne.dat',STATUS='OLD')
      READ (UNITTMP_,101) HEADER
      READ (UNITTMP_,*) ((NECR(I,J),I=1,NL),J=0,NLT)
      CLOSE (UNITTMP_)
101   FORMAT (A80)
      RETURN
      END

**--------------------------------------------------------------------**
**  Subroutine WRESULT
**  Prints out the thermal densities to a file
      SUBROUTINE WRESULT(T,DT,KP,NECR)
      
      use ModIoUnit, ONLY : io_unit_new

      INCLUDE 'numv.h'
      REAL T,NECR(NL,0:NLT),DT,KP
      INTEGER NTC
      CHARACTER name*7,SUF*2,SUF1*1

      integer ::iUnitOut != 32

      SAVE NTC

c.....Define the output name
      name='thermal'
      IF (T.EQ.0) THEN
	NTC=0
	SUF='00'
      ELSE
	NTC=NTC+1
	WRITE (SUF,11) NTC
	IF (NTC.LT.10) THEN
	  WRITE (SUF1,12) NTC
	  SUF='0'//SUF1
	END IF
      END IF
11    FORMAT (I2)
12    FORMAT (I1)

c.....Write out results
      IF (DT.GT.0.) THEN
         iUnitOut = io_unit_new()
       OPEN (iUnitOut,FILE=name//SUF,STATUS='UNKNOWN')
       WRITE (iUnitOut,*) 'Filename: '//name//SUF
       WRITE (iUnitOut,*) 'Thermal densities in the plasmasphere from plane.f'
       WRITE (iUnitOut,50) T,DT,KP
       WRITE (iUnitOut,51) (REAL(I)*.25+1.25,I=1,NL-1,2)
       DO 100 J=0,NLT
100 	 WRITE (iUnitOut,52) REAL(J)*.5,(NECR(I,J),I=1,NL-1,2)
       CLOSE (iUnitOut)
      ELSE
       OPEN (iUnitOut,FILE='ne_new.dat',STATUS='UNKNOWN')
       WRITE (iUnitOut,53) KP
       DO 110 J=0,NLT
	 WRITE (iUnitOut,54) (NECR(I,J),I=1,5)
	 WRITE (iUnitOut,54) (NECR(I,J),I=6,10)
	 WRITE (iUnitOut,54) (NECR(I,J),I=11,15)
	 WRITE (iUnitOut,54) (NECR(I,J),I=16,20)
	 WRITE (iUnitOut,54) (NECR(I,J),I=21,25)
	 WRITE (iUnitOut,54) (NECR(I,J),I=26,30)
	 WRITE (iUnitOut,54) (NECR(I,J),I=31,35)
110    CONTINUE
       CLOSE (iUnitOut)
      END IF

50    FORMAT ('T =',F10.1,'  DT =',F5.1,'  Kp =',F6.3)
51    FORMAT ('MLT\L',50(2X,F6.2,2X))
52    FORMAT (F5.1,50(1PG10.3))
53    FORMAT ('Ne in cm-3, Kp =',F12.7)
54    FORMAT (5(1PG15.7))

      RETURN
      END
