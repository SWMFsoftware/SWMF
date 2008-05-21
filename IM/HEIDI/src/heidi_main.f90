! File name: heidi_main.f90
!
! Contains: the main driver program for HEIDI
!	HEIDI
!
! Last Modified: March 2006, Mike Liemohn
!
!************************************************************************
! PROGRAM HEIDI - 
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

program heidi_main
  use ModHeidiSize
  use ModHeidiMain
  use ModHeidiDrifts
  use ModHeidiWaves
  use ModHeidiIO
  use ModProcIM

  implicit none

  integer :: nst,npr,i3,nkp,NIBC,i2
  real :: XN(NR,NS),LNC(NR,NS)

  !.......Preparation

  call MPI_INIT(iError)

  iComm= MPI_COMM_WORLD

  call MPI_COMM_RANK(iComm, iProc, iError)
  call MPI_COMM_SIZE(iComm, nProc, iError)

  CALL READPARA

  !\
  ! Check to make sure that we are running on the "correct" number of processors
  !/

  if (numprocs.gt.1) then

     total_species = 0

     DO S=1,NS
        IF (SCALC(S).EQ.1) total_species = total_species + 1
     enddo

     if (numprocs.ne.total_species) then

        if (me_world.eq.0) then
           write(*,*) "In this version, numprocs must = total_species"
           write(*,*) "numprocs      : ", numprocs
           write(*,*) "total_species : ", total_species
           write(*,*) "scalc         : ", scalc
        endif

        call MPI_abort(iComm, erno, iError)

     else

        !\
        ! Modify scalc so each processor knows which species to do
        !/

        total_species = 0
        do s=1,ns
           if (scalc(s).eq.1) then
              if (total_species.ne.me_world) then
                 scalc(s) = 0
              else
                 ispecies = s
              endif
              parallel_species(total_species+1) = s
              total_species = total_species + 1
           endif
        enddo

     endif

  else

     total_species = 0
     do s=1,ns
        if (scalc(s).eq.1) then
           parallel_species(total_species+1) = s
           total_species = total_species + 1
        endif
     enddo

  endif

  write(*,*) "me_world : ", me_world
  write(*,*) "scalc : ", scalc
  write(*,*) "total_species : ", total_species
  write(*,*) "parallel_species : ", parallel_species

  if (me_world.eq.0) then
     call ionosphere_fine_grid
     call ionosphere_init(year,day,ut)
  endif

  T=TIME
  NST=NINT(TIME/DT/2.) + 1
  NKP=NINT(10800./DT/2.)
  NIBC=NINT(TINJ/DT/2.)
  write (*,*) 'NST,NKP,NIBC:',NST,NKP,NIBC
  CALL CONSTANT(NKP)
  I2=(NST-1)/NKP + 1
  IF (IKP.GE.3) F107=F107R(I2)
  CALL ARRAYS
  CALL CEPARA
  CALL OTHERPARA
  CALL CURRENTSETUP
  CALL INITIAL(LNC,XN,J6,J18)
  CALL THERMAL   ! Setup only
  CALL GEOSB
  IF (iwpi.GT.0) THEN
     CALL WAPARA
     CALL ANISCH
  END IF
  !	IF (ISW.GT.0) CALL GETSWIND  ! Do it inside loop
  !.......Start the calculation
  NPR=NINT(TINT/DT/2.)
  if(NSTEP.LT.45) write(*,*) NSTEP
  print *, 'Times:',NST,NSTEP,NPR,NKP,NIBC,I2,DT,NSTEP*2.*DT

  DO I3=NST,NSTEP			! Begin time loop
     CALL GETKPA(i3,nst,i2,nkp)
     print *, 'I3:',I3,T,KP
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
  DO S=1,NS
     CLOSE(15+S)		! Closes continuous output file
  END DO
  CLOSE(13)		! Closes sw1 input file
  CLOSE(15)		! Closes sw2 input file
  CLOSE(14)               ! Closes MPA input file
  CLOSE(16)               ! Closes SOPA input file
  CLOSE(18)               ! Closes FPOT input file

  call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------  
  call MPI_finalize(iError)

end program heidi_main
! ************************  END OF MAIN  *******************************


