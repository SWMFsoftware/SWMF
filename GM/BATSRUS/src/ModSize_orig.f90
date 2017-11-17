module ModSize

  use BATL_size, MaxBlockBatl => MaxBlock

  ! Named indexes for coordinate directions
  use BATL_geometry, ONLY: x_, y_, z_, r_, Phi_, Theta_, Lon_, Lat_

  implicit none

  ! Maximum number of blocks per processor (set by Config.pl)
  integer, parameter :: MaxBlock = 400

  ! Maximum number of implicit blocks (set by Config.pl)
  integer, parameter :: MaxImplBLK = min(MaxBlock, 100)

end module ModSize
!==============================================================================
