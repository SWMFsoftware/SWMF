!^CMP COPYRIGHT UM
!===========================================================!
!           SWMF: Space Weather Modeling Framework          |
!                   University of Michigan                  |
!============================================================
!INTERFACE:
program SWMF

  !USES:
  use CON_main

  !REVISION HISTORY:
  ! This main program is for the stand alone SWMF.
  ! It uses methods from CON_main which can also be used when 
  ! the SWMF is a library, e.g. when it is run by an ESMF application.
  !EOP

  implicit none
  !---------------------------------------------------------------------------

  call initialize
  call run
  call finalize

end program SWMF

