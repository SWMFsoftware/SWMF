!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module BATL_size

  implicit none

  SAVE

  ! Maximum dimensionality of grid is 3 (cannot be modified)
  integer, parameter :: MaxDim = 3

  ! Number of cells per block in each direction. 
  ! These values are set by the Config.pl script.
  ! Set 1 for ignored directions! 
  integer, parameter :: nI = 8, nJ = 8, nK = 8

  ! 0 or 1 if a given dimension is ignnored or used
  integer, parameter:: iDim_ = min(nI - 1, 1)
  integer, parameter:: jDim_ = min(nJ - 1, 1)
  integer, parameter:: kDim_ = min(nK - 1, 1)

  ! Number of nodes per block in each direction
  integer, parameter:: nINode = nI + iDim_
  integer, parameter:: nJNode = nJ + jDim_
  integer, parameter:: nKNode = nK + kDim_

  ! Index values for ghost cells that may or may not exist
  integer, parameter:: j0_   =  1 - jDim_ ! j=0
  integer, parameter:: j2_   =  1 + jDim_ ! j=2
  integer, parameter:: nJp1_ = nJ + jDim_ ! j=nJ+1
  integer, parameter:: nJm1_ = nJ - jDim_ ! j=nJ-1
  integer, parameter:: k0_   =  1 - kDim_ ! k=0
  integer, parameter:: k2_   =  1 + kDim_ ! k=2
  integer, parameter:: nKp1_ = nK + kDim_ ! k=nK+1
  integer, parameter:: nKm1_ = nK - kDim_ ! k=nK-1

  ! Number of not ignored dimensions
  integer, parameter:: nDim = iDim_ + jDim_ + kDim_

  ! Refinement ratios in the 3 dimensions. Either 1 or 2.
  ! The values are set by the Config.pl script.
  integer, parameter:: &
       iRatio = min(2,nI), jRatio = min(2,nJ), kRatio = min(2,nK)

  ! Number of dimensions in which grid adaptation is done
  integer, parameter :: nDimAmr = iRatio + jRatio + kRatio - 3

  ! Indexes of AMR dimensions. 
  ! The magic formulas should be correct from 1 to nDimAmr. 
  integer, parameter, private :: iDimAmrTmp_D(MaxDim) = &
       (/ 1 + (2-iRatio)*(3-jRatio), 6-iRatio-jRatio, 3 /)

  integer, parameter :: iDimAmr_D(nDimAmr) = iDimAmrTmp_D(1:nDimAmr)

  ! Maximum number of ghost cells
  integer, parameter :: nG = 2

  ! Number of ghost cells in each direction
  integer, parameter:: nGI = nG*iDim_
  integer, parameter:: nGJ = nG*jDim_
  integer, parameter:: nGK = nG*kDim_

  ! Cell index ranges including ghost cells
  integer, parameter :: &
       MinI = 1 - nGI, MaxI = nI + nGI, &
       MinJ = 1 - nGJ, MaxJ = nJ + nGJ, &
       MinK = 1 - nGK, MaxK = nK + nGK

  ! Number of cells per block
  integer, parameter:: nIJK = nI*nJ*nK

  ! Number of nodes per block
  integer, parameter:: nIJKNode = nINode*nJNode*nKNode

  ! Array for block size
  integer, parameter:: nIJK_D(MaxDim) = (/ nI, nJ, nK /)

  ! Array for block size
  integer, parameter:: nIJKNode_D(MaxDim) = (/ nINode, nJNode, nKNode /)

  ! Array of refinement ratios
  integer, parameter:: iRatio_D(MaxDim) = (/ iRatio, jRatio, kRatio /)

  ! Ratio of cells in coarsening and prolongation
  integer, parameter:: IjkRatio = iRatio*jRatio*kRatio
  
  ! Inverse volume ratio for Cartesian case
  real, parameter:: InvIjkRatio = 1.0/IjkRatio

  ! Maximum number of blocks per processor is set dynamically
  integer :: MaxBlock = 0

  ! Largest used block index on a processor at any given time
  integer :: nBlock = 0

end module BATL_size
