
program GITM

  use ModInputs
  use ModTime
  use ModGITM

  implicit none

  integer :: iBlock

  ! ------------------------------------------------------------------------
  ! initialize stuff
  ! ------------------------------------------------------------------------

  call init_mpi
  call start_timing("GITM")
  call delete_stop

  call init_planet
  call set_defaults

  call read_inputs(cInputFile)
  call set_inputs

  call initialize_gitm

  call write_output

  call report("Starting Main Time Loop",0)

  ! ------------------------------------------------------------------------
  ! Run for a few iterations
  ! ------------------------------------------------------------------------

  do while (CurrentTime < EndTime)

     call calc_pressure

     !!! We may have to split cMax and Dt calculation!!!
     Dt = 1.e32

     call calc_timestep_vertical
     if (.not. Is1D) call calc_timestep_horizontal

     call advance

     if (.not.IsFramework) call check_stop

     iStep = iStep + 1

     call write_output

  enddo

  ! ------------------------------------------------------------------------
  ! Finish run
  ! ------------------------------------------------------------------------

  call finalize_gitm

end program GITM

!============================================================================
! The following subroutines are here so that we can use SWMF library routines
! Also some features available in SWMF mode only require empty subroutines
! for compilation of the stand alone code.
!============================================================================

subroutine CON_stop(StringError)

  implicit none
  character (len=*), intent(in) :: StringError
  call stop_gitm(StringError)

end subroutine CON_stop


subroutine CON_io_unit_new(iUnit)

  implicit none
  integer, intent(in) :: iUnit

  return

end subroutine CON_io_unit_new
