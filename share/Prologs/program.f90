!BOP
!
!MODULE: example_program - the main executable
!
!DESCRIPTION:
!
! The main executable is documented as if it was a module, because
! Protex does not provide a !PROGRAM: tag. Use LaTex syntax to 
! describe the purpose of this program.
!
!INTERFACE:
program example_program

  !USES:
  use ModExample, ONLY: ExampleParameter
  use ModExample2

  implicit none

  !LOCAL VARIABLES:
  real :: VeryImpotantVariable

  !REVISION HISTORY:
  ! 04/27/2004 G. Toth <myemail@umich.edu> - initial version
  !EOP
  ! other local variables not worth documenting
  !...
  !---------------------------------------------------------------------------
  !BOC
  write(*,*)'First executable statement of example_program worth documenting'
  !EOC
  write(*,*)'Executable statement of example_program not worth documenting'

end program example_program
