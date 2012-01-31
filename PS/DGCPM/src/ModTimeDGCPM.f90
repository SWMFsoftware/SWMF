module ModTimeDGCPM

! Time related variables only.

use ModKind, only: Real8_

real (Real8_) :: StartTime, CurrentTime, OldTime, MaxTime
integer :: YEAR, MONTH, DAY, HOUR, MINUTE, SECOND

end Module ModTimeDGCPM

