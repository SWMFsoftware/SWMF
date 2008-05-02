module ModPwGrid

  use ModNumConst, ONLY: cDegToRad

  implicit none
  save

  character (len=3) :: NamePwCoord = '???'
  integer, parameter :: nCoord=2

  ! StateGm_VI stores the PW variables that occur in GM too
  ! These are interpolated then to the GM grid
  real, allocatable :: CoordXyzPw_DI(:,:), CoordXyPw_DI(:,:), StateGm_VI(:,:)
  integer, allocatable :: iNodeTriangle_II(:, :)

  integer, parameter :: nOuterPoint=12
  real,    parameter :: dThetaOuter=5.0*cDegToRad

  integer :: nVarPw=-1, nSpeciesPw=-1, nLinePw=-1, nPoint=-1, nTriangle=-1
  integer :: iRhoPwFirst=-1, iRhoPwLast=-1, iUPwFirst=-1, iUPwLast=-1
  integer :: iRhoGmFirst=-1, iRhoGmLast=-1, iUGmFirst=-1, iUGmLast=-1


end module ModPwGrid
