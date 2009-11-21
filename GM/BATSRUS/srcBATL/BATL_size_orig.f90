module BATL_size

  implicit none

  SAVE

  ! Maximum dimensionality of grid is 3 (cannot be modified)
  integer, parameter :: MaxDim = 3

  ! Number of cells per block in each direction. 
  ! These values are set by the Config.pl script.
  ! Set 1 for ignored directions! 
  integer, parameter :: nI = 8, nJ = 8, nK = 8

  ! Number of not ignored dimensions
  integer, parameter :: nDim = min(1,nI-1) + min(1,nJ-1) + min(1,nK-1)

  ! Refinement ratios in the 3 dimensions. Either 1 or 2.
  ! The values are set by the Config.pl script.
  integer, parameter:: &
       iRatio = min(2,nI), jRatio = min(2,nJ), kRatio = min(2,nK)

  ! Number of dimensions in which grid adaptation is done
  integer, parameter :: nDimAmr = iRatio + jRatio + kRatio - 3

  ! Maximum number of ghost cells
  integer, parameter :: nG = 2

  ! Number of ghost cells in each direction
  integer, parameter:: nGI = nG*min(1,nI-1)
  integer, parameter:: nGJ = nG*min(1,nJ-1)
  integer, parameter:: nGK = nG*min(1,nK-1)

  ! Cell index ranges including ghost cells
  integer, parameter :: &
       MinI = 1 - nGI, MaxI = nI + nGI, &
       MinJ = 1 - nGJ, MaxJ = nJ + nGJ, &
       MinK = 1 - nGK, MaxK = nK + nGK

  ! Number of cells per block
  integer, parameter:: nIJK = nI*nJ*nK

  ! Array for block size
  integer, parameter:: nIJK_D(MaxDim) = (/ nI, nJ, nK /)

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

