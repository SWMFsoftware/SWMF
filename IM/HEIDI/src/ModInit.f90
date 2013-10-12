!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModInit
  use ModHeidiSize,ONLY:nR,nS
  implicit none
  SAVE

  integer :: i3
  integer :: nst
  integer :: npr

  !To count the currents              
  integer :: i2

  !How often the boundary conditions should b called
  integer :: NIBC
  real,dimension(nR,nS)::XN,LNC
  integer:: nkp

end module ModInit
