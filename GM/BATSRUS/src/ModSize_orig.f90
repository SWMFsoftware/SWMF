!^CFG COPYRIGHT UM
module ModSize
  
  use BATL_size, nDimBatl => nDim, MaxBlockBatl => MaxBlock, &
       nBlockBatl => nBlock

  implicit none

  ! number of dimensions and ghostcells. KEEP FIXED AT 3 and 2 for NOW!
  integer, parameter :: nDim=3, gcn=2

  !\
  ! Block parameters.
  !/
  ! Maximum number of blocks per processor (set by Config.pl)
  integer, parameter :: MaxBlock = 400

  integer, parameter :: nBLK = MaxBlock ! alias for convenience

  ! Maximum number of implicit blocks (set by Config.pl) !^CFG IF IMPLICIT
  integer, parameter :: MaxImplBLK = min(MaxBlock, 100)  !^CFG IF IMPLICIT

  ! Named indexes for directions and sides
  integer, parameter :: x_ = 1, y_ = 2, z_ = 3
  integer, parameter :: R_ = x_, Phi_ = y_, Theta_ = z_

  integer, parameter:: East_=1, West_=2, South_=3, North_=4, Bot_=5, Top_=6

end module ModSize
