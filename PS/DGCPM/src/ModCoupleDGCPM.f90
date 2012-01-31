module ModCoupleDGCPM

! Couple related variables only.

use ModKind, only: Real8_
use ModSizeDGCPM, only: nthetacells, nphicells

real (Real8_) :: coupled_potential(nthetacells, nphicells)
logical :: isCoupled = .false.

end Module ModCoupleDGCPM

