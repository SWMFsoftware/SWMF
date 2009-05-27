subroutine heidi_check

  use ModHeidiSize, ONLY: ns, s, scalc
  use ModHeidiMain, ONLY: TotalSpecies, ParallelSpecies
  use ModProcIM

  !---------------------------------------------------------------------------
  !\
  ! Check to make sure that we are running on the "correct" number of processors
  !/

  if (nProc.gt.1) then

     TotalSpecies = 0

     do s = 1, ns
        if (scalc(s).eq.1) TotalSpecies = TotalSpecies + 1
     enddo

     if (nProc.ne.TotalSpecies) then

        if (iProc.eq.0) then
           write(*,*) "In this version, nProc must = TotalSpecies"
           write(*,*) "nProc         : ", nProc
           write(*,*) "TotalSpecies  : ", TotalSpecies
           write(*,*) "scalc         : ", scalc
        endif

        call MPI_abort(iComm,nError,iError)

     else

        !\
        ! Modify scalc so each processor knows which species to do
        !/

        TotalSpecies = 0
        do s = 1, ns
           if (scalc(s).eq.1) then
              if (TotalSpecies.ne.iProc) then
                 scalc(s) = 0
              else
                 iSpecies = s
              endif
              ParallelSpecies(TotalSpecies+1) = s
              TotalSpecies = TotalSpecies + 1
           endif
        enddo

     endif

  else

     TotalSpecies = 0
     do s=1,ns
        if (scalc(s).eq.1) then
           ParallelSpecies(TotalSpecies+1) = s
           TotalSpecies = TotalSpecies + 1
        endif
     enddo

  endif

  write(*,*) "iProc : ", iProc
  write(*,*) "scalc (which species): ", scalc
  write(*,*) "total species : ", TotalSpecies
  write(*,*) "parallel species : ", ParallelSpecies

end subroutine heidi_check
