!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
! File name: dgcpm_main_010.f90
!
! Contains: the main driver program for DGCPM
!	DGCPM
!
! Last Modified: December 2006, Mike Liemohn
!
!************************************************************************
! PROGRAM DGCPM - 
!   ***  Hot Electron and Ion Drift Integrator  ***
! A program to calculate the growth and decay of a given population of 
! ring current ions and electrons solving the bounce-averaged 
! Boltzmann equation considering drifts, charge exchange, atmospheric 
! loss, wave-particle interactions, electric and magnetic field effects
! and feedback, and Coulomb collisions.
!
! Numerical scheme for cons terms: Lax-Wendroff + superbee flux limiter
! Converted from mram05.f into heidi010.f90, March 2006, Mike Liemohn
!***********************************************************************

!.......NR=no. grids in radial direction, NT=no. grids in azimuth,
!	NE=no. of energy grids, NS=no. of species (e-, H+, He+, O+),
!	NPA=no. of grids in equatorial pitch angle

program dgcpm

  use ModMainDGCPM
  use ModIoDGCPM
  use ModTimeDGCPM

  implicit none

  integer npr,i3

  !.......Preparation

  call DGCPC_read_inputs(cInputFile)

  call readpara

  t=time
  nst=nint(time/dt/2.) + 1
  nkp=nint(10800./dt/2.)
  nibc=nint(tinj/dt/2.)

  write (*,*) 'nst,nkp,nibc:',nst,nkp,nibc

  call constant(nkp)

  i2=(nst-1)/nkp + 1

  if (ikp >= 3) f107=f107r(i2)

  call arrays
!  call otherpara
  call thermal   ! setup only

!.......start the calculation

  npr=nint(tint/dt/2.)
  write(*,*) 'times:',nst,nstep,npr,nkp,nibc,i2,dt,nstep*2.*dt

  do i3=nst,nstep			! begin time loop
    
     open(unit=17, file='test.txt', format='formatted', access='APPEND') 
     write(17,*) 'i3:', i3, ' T:', t, ' kp:', kp
     close(17)

     call getkpa(i3,nst,i2,nkp)

     call magconv()
     call thermal

     if (i3.eq.nst) call wresult(1)

!.......increment time

     t=t+2.*dt
     CurrentTime = CurrentTime + 2.*dt

!.......print desired result files at every tint sec 

     if (mod(i3,npr) == 0 .or. i3 == nstep) then
        call wresult(0)
     end if

  enddo

! ************************  END OF MAIN  *******************************

end program dgcpm
