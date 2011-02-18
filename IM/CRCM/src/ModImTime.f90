Module ModImTime
  use ModKind,only:real8_

  real (kind=real8_)   ::  StartTime,CurrentTime !needed for SWMF stuff
  integer              ::  iStartTime_I(7)  =(/1976,6,28,0,0,0,0/)
  integer              ::  iCurrentTime_I(7)=(/1976,6,28,0,0,0,0/)
!  integer              ::  iDOY,IYD
end Module ModImTime
