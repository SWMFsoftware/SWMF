!BOP ======================================================================
!ROUTINE: example_routine - put short description here

!INTERFACE:
subroutine example_routine(InputVar, OutputVar, InputOutputVar)

  !USES:
  use ModExample

  !INPUT ARGUMENTS: 
  real, intent(in) :: InputVar          ! short description of InputVar

  !OUTPUT ARGUMENTS:
  logical, intent(out) :: OutputVar     ! short description of OutputVar

  !INPUT/OUTPUT ARGUMENTS: 
  real, intent(inout) :: InputOutputVar ! short description of InputOutputVar

  !DESCRIPTION: 
  ! This is a template for a subroutine that is not part of another
  ! module or procedure.
  !
  ! Provide here long description of subroutine example\_routine 
  ! in Latex format.
  ! Do not repeat information already provided by the above tags.

  !LOCAL VARIABLES:
  real :: AnImportantLocalVariable

  !REVISION HISTORY: 
  ! 04/27/2004 G. Toth <myemail@umich.edu> - initial version
  ! 04/28/2004 G. Toth <myemail@umich.edu> - fixed some typos
  !EOP  

  ! local variables not worth of documenting come here

  character(len=*), parameter:: NameSub='example_routine'

  !------------------------------------------------------------------------
  !BOC
  write(*,*)'Executable statement worth documenting come here'
  !EOC

  write(*,*)'This part should not appear in the documentation'

end subroutine example_routine

