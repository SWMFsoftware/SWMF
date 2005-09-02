!^CMP COPYRIGHT UM
!===========================================================!
!           SWMF: Space Weather Modeling Framework          |
!                   University of Michigan                  |
!============================================================
!BOP 
!
!QUOTE: \chapter{Control Module}
!QUOTE: \section{CON/Control: Main Executable and Control}
!MODULE: swmf - the main program for the stand-alone SWMF executable
!INTERFACE:
program SWMF

  !USES:
  use CON_main,      ONLY: initialize, run, finalize
  use CON_variables, ONLY: iErrorSwmf
  use ModMpi,        ONLY: MPI_SUCCESS

  implicit none

  !LOCAL VARIABLES:
  integer :: iErrorMpi ! MPI error code

  !DESCRIPTION:
  ! This main program is for the stand alone SWMF.
  ! It uses methods from CON\_main which can also be used when 
  ! the SWMF is a library, e.g. when it is run by an ESMF application.

  !REVISION HISTORY:
  ! 09/01/05 G.Toth - initial version
  !EOP
  !---------------------------------------------------------------------------
  !BOC
  call MPI_init(iErrorMpi)
  if(iErrorMpi  == MPI_SUCCESS) call initialize
  if(iErrorSwmf == 0          ) call run
  if(iErrorSwmf == 0          ) call finalize
  call MPI_finalize(iErrorMpi)
  !EOC
end program SWMF
