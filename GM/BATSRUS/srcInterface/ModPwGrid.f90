module ModPwGrid

  use ModNumConst, ONLY: cDegToRad

  implicit none
  save

  integer, parameter :: nCoord=2

  real, allocatable :: CoordXyPw_DI(:,:), StatePw_VI(:,:)
  integer, allocatable :: iNodeTriangle_II(:, :)
  integer, parameter :: Ub_=1, RhoPw_=2
  integer, parameter :: nOuterPoint=12
  real,    parameter :: dThetaOuter=5.0*cDegToRad

  integer :: nVarPw=-1, nSpeciesPw=-1, nLinePw=-1, nPoint=-1, nTriangle=-1
  integer :: iRhoPwFirst=-1, iRhoPwLast=-1, iUPwFirst=-1, iUPwLast=-1


end module ModPwGrid
