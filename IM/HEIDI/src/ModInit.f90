module ModInit
  use ModHeidiSize,ONLY:nR,nS
  implicit none
  SAVE

  integer :: i3
  
  integer :: nst

! To be verified
  integer :: npr

  !To count the currents              
  integer :: i2

  !How often the boundary conditions should b called
  integer :: NIBC

  real,dimension(nR,nS)::XN,LNC
  
  integer:: nkp
end module ModInit
