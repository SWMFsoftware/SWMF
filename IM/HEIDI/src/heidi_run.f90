subroutine heidi_run
  
  use ModInit
  use ModHeidiIO,   ONLY: isw, nStep,iwpi,ires,&
       write_prefix,iUnitStdOut
  use ModHeidiMain, ONLY: s,t,dt,scalc
  use ModProcIM, ONLY: iProc
  
  implicit none 
  !------------------------------------------------------------------------------

!  call OTHERPARA
  call GETKPA(i3,nst,i2,nkp)
  if (ISW.gt.0) call GETSWIND
  call MAGCONV(i3,nst)
  call get_E_mu_dot
  call THERMAL
  if (I3.eq.NST) call WRESULT(LNC,XN,1)
  do S = 1,NS
     if (SCALC(S).eq.1) then
        call LMPLOSS
        !\
        ! Call routines to calculate the changes of distribution function
        ! considering drifts, charge exchange and Coulomb drag
        !/
        call DRIFTR
        call DRIFTP
        call DRECOUL
        call DRIFTMU
        call CHAREXCHANGE
        call COULMU
        call CHAREXCHANGE
        call DRIFTMU
        call DRECOUL
        call DRIFTP
        call DRIFTR
        call LLIMIT	  ! Truncates results below 1E-30
        !\
        ! Print sources and losses (continuous output stream)
        !/
        if (IRES(14).eq.1) call PSRCLOSS(T)
     end if  ! SCALC check
  end do  ! S loop
  
  !\
  ! Increment time
  !/
  T = T + 2.*DT
  if (iProc==0) then
     call write_prefix; write(iUnitStdOut,*)&
          'SimulationTime=', T
  end if

  !\
  ! Print desired result files at every TINT sec 
  !/

  call FCHECK(11)
  if (mod(I3,NPR).eq.0 .or. I3.eq.NSTEP) then
     call FCHECK(10)	! Checks for negative results
     call WRESULT(LNC,XN,0)
     if (iwpi.gt.0) call ANISCH
  end if
  
  !\
  ! Update boundary condition
  !/
  if (mod(I3,NIBC).eq.0) call GEOSB
 
end subroutine heidi_run
