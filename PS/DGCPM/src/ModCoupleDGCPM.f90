!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModCoupleDGCPM

! Couple related variables only.

use ModKind, only: Real8_
use ModSizeDGCPM, only: nthetacells, nphicells

real (Real8_) :: coupled_potential(nthetacells, nphicells)
logical :: isCoupled = .false.

end Module ModCoupleDGCPM

