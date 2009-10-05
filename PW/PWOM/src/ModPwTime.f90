Module ModPwTime
  use ModKind,only:real8_
  
  real (kind=real8_)   ::  StartTime,CurrentTime !needed for SWMF stuff
  integer              ::  iStartTime(7)=(/1976,6,28,0,0,0,0/)
  integer, parameter   ::  Year_=1, Month_=2, Day_=3, Hour_=4, Minute_=5, &
                           Second_=6, mSecond_=7

end Module ModPwTime
