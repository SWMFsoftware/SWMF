
module ModUserGITM

  use ModGITM

  implicit none

  integer, parameter :: nUserOutputs = 40
  real :: UserData3D(-1:nLons+2,-1:nLats+2,-1:nAlts+2,nUserOutputs,nBlocksMax)=0.0
  real :: UserData2D(-1:nLons+2,-1:nLats+2,         1,nUserOutputs,nBlocksMax)=0.0

  integer :: nVarsUser3d=0
  integer :: nVarsUser2d=0

end module ModUserGITM
