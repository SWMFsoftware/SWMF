subroutine heidi_init

  use ModHeidiSize
  use ModHeidiMain
  use ModHeidiDrifts
  use ModHeidiWaves
  use ModHeidiIO
  use ModProcIM
  use ModInit
  
  implicit none

  !integer:: nst,npr,nkp,NIBC,i2

  if (me_world.eq.0) then
     call IonoHeidiInit(year,day,ut)
  endif
  
  T=TIME
  write(*,*) 'TIME =',TIME
  NST=NINT(TIME/DT/2.) + 1
  NKP=NINT(10800./DT/2.)
  NIBC=NINT(TINJ/DT/2.)
  
  write (*,*) 'nSteps, nSteps in KP, nSteps IBC:',NST,NKP,NIBC

  CALL CONSTANT(NKP)
  I2=(NST-1)/NKP + 1
  write (*,*) 'I2', I2

  IF (IKP.GE.3) F107=F107R(I2)
  write(*,*) 'F107', F107
  
  write(*,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  write(*,*) 'Calling ARRAYS'
  CALL ARRAYS
  write(*,*) 'Calling CEPARA'
  CALL CEPARA
  write(*,*) 'Calling OTHERPARA'
  CALL OTHERPARA
  write(*,*) 'Calling CURRENTSETUP'
  CALL CURRENTSETUP
  write(*,*) 'Calling INITIAL'
  CALL INITIAL(LNC,XN,J6,J18)
  write(*,*) 'LNC,XN,J6,J18', LNC,XN,J6,J18
  write(*,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  
  write(*,*) 'Calling THERMAL'
  CALL THERMAL   ! Setup only
  write(*,*) 'Calling GEOSB'
  CALL GEOSB
  write(*,*) 'IWPI', iwpi
  IF (iwpi.GT.0) THEN
     CALL WAPARA
     CALL ANISCH
  END IF
  !	IF (ISW.GT.0) CALL GETSWIND  ! Do it inside loop
  !.......Start the calculation
  NPR=NINT(TINT/DT/2.)
  if(NSTEP.LT.45) write(*,*) NSTEP
  write(*,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  print *, 'Times:',NST,NSTEP,NPR,NKP,NIBC,I2,DT,NSTEP*2.*DT

  write(*,*) 'Initialization done!'

!!$  T=TIME
!!$
!!$  iStepStart=NINT(TIME/DT/2.) + 1
!!$  iStep = iStepStart-1
!!$  nFrequencyKp=NINT(10800./DT/2.)
!!$  nFrequencyBC = NINT(TINJ/DT/2.)
!!$
!!$  write (*,*) 'iStepStart,nFrequencyKp,nFrequencyBC:',&
!!$       iStepStart,nFrequencyKp,nFrequencyBC
!!$
!!$  CALL CONSTANT(nFrequencyKp)
!!$  nCounterForCurrents = (iStepStart-1)/nFrequencyKp + 1
!!$  IF (IKP.GE.3) F107=F107R(nCounterForCurrents)
!!$
!!$  CALL ARRAYS
!!$  CALL CEPARA
!!$  CALL OTHERPARA
!!$  CALL CURRENTSETUP
!!$  CALL INITIAL(LNC,XN,J6,J18)
!!$  CALL THERMAL   ! Setup only
!!$  CALL GEOSB
!!$
!!$  IF (iwpi.GT.0) THEN
!!$     CALL WAPARA
!!$     CALL ANISCH
!!$  END IF
!!$  !	IF (ISW.GT.0) CALL GETSWIND  ! Do it inside loop
!!$  !.......Start the calculation
!!$
!!$  nOutputFrequency=NINT(TINT/DT/2.)
!!$  if(NSTEP.LT.45) write(*,*) nStep
!!$
!!$  print *, 'Times:',iStepStart,&
!!$                    nStep,     &
!!$                    nOutputFrequency,&
!!$                    nFrequencyKp, &
!!$                    nOutputFrequency,nCounterForCurrents,DT,nStep*2.*DT

end subroutine heidi_init
