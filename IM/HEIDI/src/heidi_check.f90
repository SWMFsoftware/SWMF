!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine heidi_check

  use ModHeidiSize, ONLY: ns, s, scalc
  use ModHeidiMain, ONLY: nSpecies, nParallelSpecies
  use ModProcIM

  !---------------------------------------------------------------------------
  !\
  ! Check to make sure that we are running on the "correct" number of processors
  !/

  if (nProc.gt.1) then

     nSpecies = 0

     do s = 1, ns
        if (scalc(s).eq.1) nSpecies = nSpecies + 1
     enddo

     !if (nProc > sum(scalc)) then
     if (nProc.ne.nSpecies) then

        if (iProc.eq.0) then
           write(*,*) "In this version, nProc must = nSpecies"
           write(*,*) "nProc         : ", nProc
           write(*,*) "nSpecies      : ", nSpecies
           write(*,*) "scalc         : ", scalc
        endif

        call MPI_abort(iComm,nError,iError)

     else

        !\
        ! Modify scalc so each processor knows which species to do
        !/

        nSpecies = 0
        do s = 1, ns
           if (scalc(s).eq.1) then
              if (nSpecies.ne.iProc) then
                 scalc(s) = 0
              else
                 iSpecies = s
              endif
              nParallelSpecies(nSpecies+1) = s
              nSpecies = nSpecies + 1
           endif
        enddo

     endif

  else

     nSpecies = 0
     do s=1,ns
        if (scalc(s).eq.1) then
           nParallelSpecies(nSpecies+1) = s
           nSpecies = nSpecies + 1
        endif
     enddo

  endif

  write(*,*) "iProc : ", iProc
  write(*,*) "scalc (which species): ", scalc
  write(*,*) "total species : ", nSpecies
  write(*,*) "parallel species : ", nParallelSpecies

end subroutine heidi_check
