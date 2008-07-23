subroutine heidi_run
  use ModInit
  use ModHeidiSize
  use ModHeidiIO
  use ModHeidiMain
  
  implicit none 

  CALL GETKPA(iStep,iStepStart,nCounterForCurrents,nFrequencyKp)

  print *, 'I3:',iStep,T,KP
  IF (ISW.GT.0) CALL GETSWIND

  print *, 'Calling MAGCONV'
  CALL MAGCONV(iStep,iStepStart)

  print *, 'Calling THERMAL'
  CALL THERMAL

  print *, 'Calling WRESULT'
  IF (iStep.EQ.iStepStart ) CALL WRESULT(LNC,XN,1)

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
  IF (MOD(iStep,nOutputFrequency).EQ.0 .OR. iStep.EQ.NSTEP) THEN
     CALL FCHECK(10)	! Checks for negative results
     CALL WRESULT(LNC,XN,0)
     IF (iwpi.GT.0) CALL ANISCH
  END IF
  !	Update boundary condition

  write(*,*) "kp(heidi_run) : ",kp
  IF (MOD(iStep,nFrequencyBC).EQ.0) CALL GEOSB
  
  iStep = iStep + 1

end subroutine heidi_run
