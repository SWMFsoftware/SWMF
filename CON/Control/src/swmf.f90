!^CMP COPYRIGHT UM
!===========================================================!
!           SWMF: Space Weather Modeling Framework          |
!                   University of Michigan                  |
!============================================================
!INTERFACE:
program SWMF

  !USES:
  use CON_main
  use CON_variables, ONLY: iErrorSwmf
  use ModMpi

  !REVISION HISTORY:
  ! This main program is for the stand alone SWMF.
  ! It uses methods from CON_main which can also be used when 
  ! the SWMF is a library, e.g. when it is run by an ESMF application.
  !EOP

  implicit none
  integer :: iError
  !---------------------------------------------------------------------------

  call MPI_init(iError)
  if(iError==MPI_SUCCESS) call initialize
  if(iErrorSwmf==0)       call run
  if(iErrorSwmf==0)       call finalize
  call MPI_finalize(iError)

end program SWMF

