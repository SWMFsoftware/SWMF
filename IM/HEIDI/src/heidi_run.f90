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
  call FCHECK(1)
  call MAGCONV(i3,nst)
  call FCHECK(2)
  call get_E_mu_dot
  call FCHECK(3)
  call THERMAL
  call FCHECK(4)
  if (I3.eq.NST) call WRESULT(LNC,XN,1)
  call FCHECK(5)
  do S = 1,NS
     if (SCALC(S).eq.1) then
        call LMPLOSS
        !\
        ! Call routines to calculate the changes of distribution function
        ! considering drifts, charge exchange and Coulomb drag
        !/
        
        call FCHECK(5)

        call DRIFTR
        call FCHECK(6)
        call heidi_driftp
        call FCHECK(7)
        call DRECOUL
        call FCHECK(8)
        call DRIFTMU
        call FCHECK(9)
        call heidi_charexchange
        call FCHECK(10)
        call COULMU
        call FCHECK(11)
        call heidi_charexchange
        call FCHECK(12)
        call DRIFTMU
        call FCHECK(13)
        call DRECOUL
        call FCHECK(14)
        call heidi_driftp
        call FCHECK(15)
        call DRIFTR
        call FCHECK(16)
        call LLIMIT	  ! Truncates results below 1E-30
        call FCHECK(17)
        !\
        ! Print sources and losses (continuous output stream)
        !/
        if (IRES(14).eq.1) call PSRCLOSS(T)
          call FCHECK(18)
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

  call FCHECK(19)
  if (mod(I3,NPR).eq.0 .or. I3.eq.NSTEP) then
     call FCHECK(20)	! Checks for negative results
     call WRESULT(LNC,XN,0)
     call FCHECK(21)
     if (iwpi.gt.0) call ANISCH
     call FCHECK(22)
  end if
  
  !\
  ! Update boundary condition
  !/
  if (mod(I3,NIBC).eq.0) call GEOSB
    call FCHECK(23)

end subroutine heidi_run
