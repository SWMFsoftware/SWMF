
subroutine init_energy_deposition

  use ModInputs, only: iDebugLevel
  use ModGITM, only: Altitude
  use ModSources

  implicit none

  integer :: ierr, iAlt, i
  real :: a
  logical :: IsDone

  !--------------------------------------
  ! Start doing energy deposition stuff
  !--------------------------------------

  call report("init_energy_deposition",2)

  if (iDebugLevel > 2) write(*,*) "===> ED_Init"
  call ED_Init(ierr)

  if (ierr /= 0) then
     call stop_gitm("Error in initilizing the energy deposition tables")
  endif

  if (iDebugLevel > 2) write(*,*) "===> ED_Get_Grid_Size"
  call ED_Get_Grid_Size(ED_N_Alts)
  allocate(ED_grid(ED_N_Alts), stat=ierr)
  if (ierr /= 0) then
     call stop_gitm("Error allocating array ED_grid")
  endif

  allocate(ED_Ion(ED_N_Alts), stat=ierr)
  if (ierr /= 0) then
     call stop_gitm("Error allocating array ED_Ion")
  endif

  allocate(ED_Heating(ED_N_Alts), stat=ierr)
  if (ierr /= 0) then
     call stop_gitm("Error allocating array ED_Heating")
  endif

  if (iDebugLevel > 2) write(*,*) "===> ED_Get_Number_of_Energies"
  call ED_Get_Number_of_Energies(ED_N_Energies)
  allocate(ED_Energies(ED_N_Energies), stat=ierr)
  if (ierr /= 0) then
     call stop_gitm("Error allocating array ED_Energies")
  endif

  allocate(ED_Flux(ED_N_Energies), stat=ierr)
  if (ierr /= 0) then
     call stop_gitm("Error allocating array ED_Flux")
  endif

  if (iDebugLevel > 2) write(*,*) "===> ED_Get_Grid"
  call ED_Get_Grid(ED_grid, .false., ierr)

  do iAlt = 1, nAlts

     a = Altitude(iAlt)

     i = 1
     IsDone = .false.
     do while (.not.IsDone)
        if (ED_grid(i) <= a .and. ED_grid(i+1) >= a) then
           IsDone = .true.
           ED_Interpolation_Index(iAlt) = i
           ED_Interpolation_Weight(iAlt) = (ED_grid(i) - a) /  &
                (ED_grid(i) - ED_grid(i+1))
        else
           if (i == ED_N_Alts-1) then
              IsDone = .true.
              ED_Interpolation_Weight(iAlt) = 1.0
              ED_Interpolation_Index(iAlt) = ED_N_Alts - 1
           else
              i = i + 1
           endif
        endif
     enddo

  enddo

  call ED_Get_Energies(ED_Energies)

end subroutine init_energy_deposition
