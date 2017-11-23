!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module BATL_size

  implicit none

  SAVE

  ! Number of cells per block in each direction. 
  ! These values are set by the Config.pl script.
  ! Set 1 for ignored directions! 
  integer, parameter :: nI = 8, nJ = 8, nK = 8

  ! Maximum number of ghost cells set by Config.pl script.
  ! Valid values are 0,1,2,3,4,5
  integer, parameter :: nG = 2

  ! Refinement ratios in the 3 dimensions. Either 1 or 2.
  ! The values are set by the Config.pl script.
  integer, parameter:: &
       iRatio = min(2, nI), jRatio = min(2, nJ), kRatio = min(2, nK)

  ! Maximum dimensionality of grid is 3 (cannot be modified)
  integer, parameter :: MaxDim = 3

  ! 0 or 1 if a given dimension is ignored or used
  integer, parameter:: iDim_ = min(nI - 1, 1)
  integer, parameter:: jDim_ = min(nJ - 1, 1)
  integer, parameter:: kDim_ = min(nK - 1, 1)

  ! Number of not ignored dimensions
  integer, parameter:: nDim = iDim_ + jDim_ + kDim_

  ! Index names limited by nDim
  integer, parameter:: Dim1_=1, Dim2_=min(nDim,2), Dim3_=min(nDim,3)

  ! Number of dimensions in which grid adaptation is done
  integer, parameter :: nDimAmr = iRatio + jRatio + kRatio - 3

  ! Number of nodes per block in each direction
  integer, parameter:: nINode = nI + iDim_
  integer, parameter:: nJNode = nJ + jDim_
  integer, parameter:: nKNode = nK + kDim_

  ! Number of ghost cells in each direction
  integer, parameter:: nGI = nG*iDim_
  integer, parameter:: nGJ = nG*jDim_
  integer, parameter:: nGK = nG*kDim_

  ! Cell index ranges including ghost cells
  integer, parameter :: &
       MinI = 1 - nGI, MaxI = nI + nGI, &
       MinJ = 1 - nGJ, MaxJ = nJ + nGJ, &
       MinK = 1 - nGK, MaxK = nK + nGK

  ! Index values for cells limited to the existing range
  integer, parameter:: im2_  = max(    -2, MinI) ! i=-2
  integer, parameter:: im1_  = max(    -1, MinI) ! i=-1
  integer, parameter:: i0_   = max(     0, MinI) ! i=0
  integer, parameter:: i2_   = min(     2, nI)   ! i=2
  integer, parameter:: i3_   = min(     3, nI)   ! i=3
  integer, parameter:: nIm2_ = max(nI - 2, 1)    ! i=nI-2
  integer, parameter:: nIm1_ = max(nI - 1, 1)    ! i=nI-1
  integer, parameter:: nIp1_ = min(nI + 1, MaxI) ! i=nI+1
  integer, parameter:: nIp2_ = min(nI + 2, MaxI) ! i=nI+2
  integer, parameter:: nIp3_ = min(nI + 3, MaxI) ! i=nI+3

  integer, parameter:: jm2_  = max(    -2, MinJ) ! j=-2
  integer, parameter:: jm1_  = max(    -1, MinJ) ! j=-1
  integer, parameter:: j0_   = max(     0, MinJ) ! j=0
  integer, parameter:: j2_   = min(     2, nJ)   ! j=2
  integer, parameter:: j3_   = min(     3, nJ)   ! j=3
  integer, parameter:: nJm2_ = max(nJ - 2, 1)    ! j=nJ-2
  integer, parameter:: nJm1_ = max(nJ - 1, 1)    ! j=nJ-1
  integer, parameter:: nJp1_ = min(nJ + 1, MaxJ) ! j=nJ+1
  integer, parameter:: nJp2_ = min(nJ + 2, MaxJ) ! j=nJ+2
  integer, parameter:: nJp3_ = min(nJ + 3, MaxJ) ! j=nJ+3

  integer, parameter:: km2_  = max(    -2, MinK) ! k=-2
  integer, parameter:: km1_  = max(    -1, MinK) ! k=-1
  integer, parameter:: k0_   = max(     0, MinK) ! k=0
  integer, parameter:: k2_   = min(     2, nK)   ! k=2
  integer, parameter:: k3_   = min(     3, nK)   ! k=3
  integer, parameter:: nKm2_ = max(nK - 2, 1)    ! k=nK-2
  integer, parameter:: nKm1_ = max(nK - 1, 1)    ! k=nK-1
  integer, parameter:: nKp1_ = min(nK + 1, MaxK) ! k=nK+1
  integer, parameter:: nKp2_ = min(nK + 2, MaxK) ! k=nK+2
  integer, parameter:: nKp3_ = min(nK + 3, MaxK) ! k=nK+3

  ! Indexes of AMR dimensions. 
  ! The magic formulas should be correct from 1 to nDimAmr. 
  integer, parameter, private :: iDimAmrTmp_D(MaxDim) = &
       (/ 1 + (2-iRatio)*(3-jRatio), 6-iRatio-jRatio, 3 /)

  integer, parameter :: iDimAmr_D(nDimAmr) = iDimAmrTmp_D(1:nDimAmr)

  ! Number of cells per block
  integer, parameter:: nIJK = nI*nJ*nK

  ! Number of cells per block including ghost cells
  integer, parameter:: MaxIJK = (MaxI-MinI+1)*(MaxJ-MinJ+1)*(MaxK-MinK+1)

  ! Number of nodes per block
  integer, parameter:: nIJKNode = nINode*nJNode*nKNode

  ! Arrays for block size
  integer, parameter:: nIJK_D(MaxDim) = (/ nI, nJ, nK /)
  integer, parameter:: MinIJK_D(MaxDim) = (/MinI, MinJ, MinK/)
  integer, parameter:: MaxIJK_D(MaxDim) = (/MaxI, MaxJ, MaxK/)

  ! Array for block size
  integer, parameter:: nIJKNode_D(MaxDim) = (/ nINode, nJNode, nKNode /)

  ! Array of refinement ratios
  integer, parameter:: iRatio_D(MaxDim) = (/ iRatio, jRatio, kRatio /)

  ! Ratio of cells in coarsening and prolongation
  integer, parameter:: IjkRatio = iRatio*jRatio*kRatio
  
  ! Inverse volume ratio for Cartesian case
  real, parameter:: InvIjkRatio = 1.0/IjkRatio

  ! Maximum number of blocks per processor is set dynamically
  integer :: MaxBlock = 400

  ! Largest used block index on a processor at any given time
  integer :: nBlock = 0

  ! Number of different kinds of particles
  integer, parameter:: nKindParticle = 8

end module BATL_size
