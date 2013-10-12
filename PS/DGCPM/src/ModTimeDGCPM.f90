!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModTimeDGCPM

! Time related variables only.

use ModKind, only: Real8_

real (Real8_) :: StartTime, CurrentTime, OldTime, MaxTime
integer :: YEAR, MONTH, DAY, HOUR, MINUTE, SECOND

end Module ModTimeDGCPM

