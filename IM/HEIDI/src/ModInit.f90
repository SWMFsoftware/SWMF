module ModInit
  use ModHeidiSize,ONLY:nR,nS
  implicit none
  SAVE

  integer :: iStep           !i3
  
  integer :: iStepStart      !nst

! To be verified
  integer :: nOutputFrequency !npr

  !To count the currents              
  integer :: nCounterForCurrents !i2

  !How often the boundary conditions should b called
  integer :: nFrequencyBC        !NIBC

  real,dimension(nR,nS)::XN,LNC
  
  integer:: nFrequencyKp  !Freqiency of output of the Kp index
end module ModInit
