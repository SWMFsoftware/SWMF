subroutine heidi_check
  use ModHeidiSize
  use ModHeidiMain
  use ModHeidiDrifts
  use ModHeidiWaves
  use ModHeidiIO
  use ModProcIM


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

end subroutine heidi_check
