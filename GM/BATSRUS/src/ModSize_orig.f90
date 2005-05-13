!^CFG COPYRIGHT UM
module ModSize
  !\
  ! Block parameters.
  !/
  ! Maximum number of blocks per processor
  integer, parameter :: nBLK=400

  integer, parameter :: MaxBlock = nBLK

  ! Maximum number of implicit blocks                   !^CFG IF IMPLICIT
  integer, parameter :: MaxImplBLK = min(MaxBlock, 100) !^CFG IF IMPLICIT

  !\
  ! Block size in terms of number of cells in X, Y and Z directions.
  !   Minimum is 2x2x2 for fixed grid, 4x4x4 with AMR.
  !   Use even numbers only.
  !/
  integer, parameter :: nI=8, nJ=8, nK=8

  integer, parameter :: nIJK=nI*nJ*nK
  integer, parameter, dimension(3) :: nCells=(/nI,nJ,nK/)

  !\
  ! number of dimensions and ghostcells. KEEP FIXED AT 3 and 2!
  !/
  integer, parameter :: nDim=3, gcn=2

  ! Parameters to hold direction names
  integer, parameter:: east_=1, west_=2, south_=3, north_=4, bot_=5, top_=6

end module ModSize
