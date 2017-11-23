module ModSize

  use BATL_size

  ! Named indexes for coordinate directions
  use BATL_geometry, ONLY: x_, y_, z_, r_, Phi_, Theta_, Lon_, Lat_

  implicit none

  ! Maximum number of implicit blocks (set by Config.pl)
  integer, parameter :: MaxImplBLK = 100

end module ModSize
!==============================================================================
