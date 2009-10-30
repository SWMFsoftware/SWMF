module BATL_size

  implicit none

  SAVE

  ! Dimensionality of grid and AMR
  integer, parameter :: MaxDim = 3    ! This has to be 3 all the time
  integer, parameter :: nDim   = 3    ! Number of not ignored dimensions
  integer, parameter :: nDimAmr = 3  ! Number of refined dimensions

  ! Note: 3 = MaxDim >= nDim >= nDimAmr >= 1

  ! Maximum number of blocks per processor
  integer :: MaxBlock = 0

  ! Largest block index
  integer :: nBlock = 0

  ! Number of cells per block in each direction
  integer, parameter :: nI = 8, nJ = 4, nK = 2

  ! Number of cells per block
  integer, parameter :: nIJK = nI*nJ*nK

  ! Array for block size
  integer, parameter:: &
       nIJK_D(MaxDim) = (/ nI, nJ, nK /)

  ! number of ghost cells
  integer, parameter :: nG = 2

  integer, parameter:: nGI = nG
  integer, parameter:: nGJ = nG*min(1,nDim-1)
  integer, parameter:: nGK = nG*max(0,nDim-2)
  
  integer, parameter :: &
       MinI = 1 - nGI, MaxI = nI + nGI, &
       MinJ = 1 - nGJ, MaxJ = nJ + nGJ, &
       MinK = 1 - nGK, MaxK = nK + nGK

  ! Refinement ratios in the 3 dimensions (depends on nDimAmr)
  integer, parameter:: &
       iRatio = 2, jRatio = min(2,nDimAmr), kRatio = max(1,nDimAmr-1)

  ! Array of refinement ratios
  integer, parameter:: iRatio_D(MaxDim) = (/ iRatio, jRatio, kRatio /)

  ! Ratio of cells in coarsening and prolongation
  integer, parameter:: IjkRatio = iRatio*jRatio*kRatio
  
  ! Inverse volume ratio for Cartesian case
  real, parameter:: InvIjkRatio = 1.0/IjkRatio

end module BATL_size

