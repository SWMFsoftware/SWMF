subroutine heidi_init

  use ModHeidiMain, only: T, f107R, Re, DipoleFactor
  use ModHeidiDrifts, only: j18,j6
  use ModProcIM, only: iProc
  use ModHeidiIO
  use ModInit

  use CON_planet, ONLY: get_planet
  implicit none

  !--------------------------------------------------------------------------

  if (iProc == 0) then
     call IonoHeidiInit(year,day,ut)
  endif

  if (.not. IsFramework) T = TIME

  call get_planet(RadiusPlanetOut = Re, DipoleStrengthOut = DipoleFactor)
  DipoleFactor = DipoleFactor/Re**3

  if (iProc == 0) then
     call write_prefix; write(iUnitStdOut,*) 'TIME =',TIME
  end if
  NST=nint(TIME/DT/2.) + 1
  NKP=nint(10800./DT/2.)
  NIBC=nint(TINJ/DT/2.)

  if (iProc==0) then
     call write_prefix; write(iUnitStdOut,*)'nSteps, nSteps in KP, nSteps IBC:',NST,NKP,NIBC
  end if

  call CONSTANT(NKP)
  I2=(NST-1)/NKP + 1
  if (IKP.ge.3) F107=F107R(I2)
  call ARRAYS
  call heidi_cepara
  call OTHERPARA
  call CURRENTSETUP
  call heidi_initial(LNC,XN,J6,J18)
  call THERMAL   ! Setup only
  call GEOSB
  if (iwpi.gt.0) then
     call WAPARA
     call ANISCH
  end if
  !	IF (ISW.GT.0) CALL GETSWIND  ! Do it inside loop
  !.......Start the calculation
  NPR=nint(TINT/DT/2.)
  if (NSTEP.lt.45) then
     if (iProc==0) then
        !call write_prefix; write(iUnitStdOut,*) 'NSTEP',NSTEP
        call write_prefix; write(iUnitStdOut,*) 'Times:',NST,NSTEP,NPR,NKP,NIBC,I2,DT,NSTEP*2.*DT
     end if
  end if

end subroutine heidi_init
