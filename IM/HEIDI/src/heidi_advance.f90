! File name: heidi_advance.f90
!
! Contains: the main time loop program for HEIDI
!
! Last Modified: April 2008, Raluca Ilie
!
!************************************************************************

	DO I3=NST,NSTEP			! Begin time loop
	  CALL GETKPA(i3,nst,i2,nkp)
	  IF (ISW.GT.0) CALL GETSWIND
          print *, 'Calling MAGCONV'
	  CALL MAGCONV(i3,nst)
          print *, 'Calling THERMAL'
	  CALL THERMAL
          print *, 'Calling WRESULT'
	  IF (I3.EQ.NST) CALL WRESULT(LNC,XN,1)
	  DO S=1,NS
	    IF (SCALC(S).EQ.1) THEN
	print *, '   S: ',S
	      CALL LMPLOSS
!.......Call routines to calculate the changes of distribution function
!               considering drifts, charge exchange and Coulomb drag
!	CALL FCHECK(1)
	      CALL DRIFTR
	      CALL DRIFTP
  	      CALL DRECOUL
	      CALL DRIFTMU
	      CALL CHAREXCHANGE
	      CALL COULMU
!	CALL FCHECK(2)
	      CALL CHAREXCHANGE
	      CALL DRIFTMU
	      CALL DRECOUL
	      CALL DRIFTP
	      CALL DRIFTR
	      CALL LLIMIT	  ! Truncates results below 1E-30
!	CALL FCHECK(3)
!.......Print sources and losses (continuous output stream)
	      IF (IRES(14).EQ.1) CALL PSRCLOSS(T)
	    END IF  ! SCALC check
	  END DO  ! S loop
!.......Increment time
	  T=T+2.*DT
!.......Print desired result files at every TINT sec 
	CALL FCHECK(11)
          IF (MOD(I3,NPR).EQ.0 .OR. I3.EQ.NSTEP) THEN
	    CALL FCHECK(10)	! Checks for negative results
	    CALL WRESULT(LNC,XN,0)
	    IF (iwpi.GT.0) CALL ANISCH
	  END IF
!	Update boundary condition
	  IF (MOD(I3,NIBC).EQ.0) CALL GEOSB
	end do			! End time loop

!************************  END OF MAIN  *******************************
