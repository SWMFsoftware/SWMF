!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModPwGrid

  use ModNumConst,   ONLY: cDegToRad
  use ModMultiFluid, ONLY: nIonFluid
  
  implicit none
  save

  character (len=3) :: NamePwCoord = '???'
  integer, parameter :: nCoord=3

  ! StateGm_VI stores the PW variables that occur in GM too
  ! These are interpolated then to the GM grid
  real, allocatable :: CoordXyzPw_DI(:,:), CoordXyPw_DI(:,:), StateGm_VI(:,:),&
       CoordXyzPw1_DI(:,:), CoordXyzPw2_DI(:,:), &
       StateGm1_VI(:,:),StateGm2_VI(:,:)
  
  integer, allocatable :: list1_I(:),lptr1_I(:), lend1_I(:),&
       list2_I(:),lptr2_I(:), lend2_I(:)  

  integer, parameter :: nOuterPoint=12
  real,    parameter :: dThetaOuter=5.0*cDegToRad

  integer :: iCoupleFluid(nIonFluid) = -1

  integer :: nVarPw=-1, nSpeciesPw=-1, nLinePw=-1, nPoint=-1, nTriangle=-1
  integer :: iRhoPwFirst=-1, iRhoPwLast=-1, iUPwFirst=-1, iUPwLast=-1
  integer :: iRhoGmFirst=-1, iRhoGmLast=-1, iUGmFirst=-1, iUGmLast=-1
  integer :: nPoint1=-1, nPoint2=-1,nTriangle1=-1,nTriangle2=-1
  integer :: nLinePw1 = -1, nLinePw2 = -1

end module ModPwGrid
