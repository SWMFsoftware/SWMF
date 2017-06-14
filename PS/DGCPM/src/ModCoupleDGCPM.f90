!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModCoupleDGCPM

! Couple related variables only.

use ModKind, only: Real8_
use ModSizeDGCPM, only: nthetacells, nphicells

logical :: isCoupled = .false.

integer :: nThetaIe=0, nPhiIe=0
real (Real8_), allocatable :: IePot_II(:,:), IeTheta_I(:), IePhi_I(:)

end Module ModCoupleDGCPM

