!^CFG COPYRIGHT UM
!^CFG FILE CONSTRAINB
module ModCT

  use ModSize
  implicit none
  SAVE

  !\
  ! Variables for Constrained Transport
  !/

  ! Face centered magnetic field components, Bxface is centered on face X etc.

  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK) :: &
       Bxface_BLK,Byface_BLK,Bzface_BLK

  ! Face centered normal magnetic field from the 4 finer neighbors
  ! on the two sides of the block (1=east/south/bot, 2=west/north/top)

  real, dimension(nJ,nK,4,2,nBLK) :: BxFaceFine_XQSB
  real, dimension(nI,nK,4,2,nBLK) :: ByFaceFine_YQSB
  real, dimension(nI,nJ,4,2,nBLK) :: BzFaceFine_ZQSB

  ! VxB stored at edges

  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK):: &
       VxB_x,VxB_y,VxB_z

  logical :: DoInitConstrainB = .true.

end module ModCT
