!^CFG COPYRIGHT UM
module ModSize
  
  use BATL_size, MaxBlockBatl => MaxBlock
  use BATL_geometry, ONLY: r_, Phi_, Theta_, Lon_, Lat_

  implicit none

  ! number of ghostcells. KEEP FIXED AT 2 for NOW!
  integer, parameter :: gcn=2

  ! Maximum number of blocks per processor (set by Config.pl)
  integer, parameter :: MaxBlock = 400

  integer, parameter :: nBLK = MaxBlock ! alias for convenience

  ! Maximum number of implicit blocks (set by Config.pl) !^CFG IF IMPLICIT
  integer, parameter :: MaxImplBLK = min(MaxBlock, 100)  !^CFG IF IMPLICIT

  ! Named indexes for Cartesian directions (these are limited by nDim in BATL)
  integer, parameter :: x_ = 1, y_ = 2, z_ = 3

end module ModSize
