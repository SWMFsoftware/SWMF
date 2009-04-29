subroutine heidi_init

  use ModHeidiSize
  use ModHeidiMain
  use ModHeidiDrifts
  use ModHeidiWaves
  use ModHeidiIO
  use ModProcIM
  use ModInit

  implicit none

  !--------------------------------------------------------------------------

  if (me_world.eq.0) then
     call IonoHeidiInit(year,day,ut)
  endif

  T=TIME
  call write_prefix; write(iUnitStdOut,*) 'TIME =',TIME
  NST=nint(TIME/DT/2.) + 1
  NKP=nint(10800./DT/2.)
  NIBC=nint(TINJ/DT/2.)

  call write_prefix; write(iUnitStdOut,*)'nSteps, nSteps in KP, nSteps IBC:',NST,NKP,NIBC
  call CONSTANT(NKP)
  I2=(NST-1)/NKP + 1
  if (IKP.ge.3) F107=F107R(I2)
  call ARRAYS
  call CEPARA
  call OTHERPARA
  call CURRENTSETUP
  call INITIAL(LNC,XN,J6,J18)
  call THERMAL   ! Setup only
  call GEOSB
  if (iwpi.gt.0) then
     call WAPARA
     call ANISCH
  end if
  !	IF (ISW.GT.0) CALL GETSWIND  ! Do it inside loop
  !.......Start the calculation
  NPR=nint(TINT/DT/2.)
  if(NSTEP.lt.45) call write_prefix; write(iUnitStdOut,*) 'NSTEP',NSTEP

  call write_prefix; write(iUnitStdOut,*) 'Times:',NST,NSTEP,NPR,NKP,NIBC,I2,DT,NSTEP*2.*DT

end subroutine heidi_init
