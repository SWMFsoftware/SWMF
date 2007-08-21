module ModPWOM

  use ModParameters
  use ModCommonPlanet, ONLY: nVar
  implicit none

  logical :: IsStandAlone = .false.

  integer :: iUnitOut
  character (len=7) :: StringPrefix=''

  integer :: iProc, nProc, iComm,errcode

  integer :: iProcTest = 0, iLineTest = 1
  character (len=100):: StringTest = ''

  integer, parameter :: MaxLine = 500
  
  integer :: nTotalLine=1
  integer :: iTheta, iPhi, iUnitSouth,iUnitNorth,i,iLine=0
  integer, dimension(MaxLine)::iLineGlobal

  real    ::   Bcoef,MagMoment,rPlanet,OmegaPlanet,Dtheta,Dphi,rLowerBoundary

  integer :: nLine,nLog=0

  ! The number of lines on each processor and on the processors with lower rank
  integer, allocatable :: nLine_P(:), nLineBefore_P(:)

  ! ionosphere variables
  integer :: nTheta = -1, nPhi = -1
  real, dimension(:,:), allocatable :: &
       Theta_G,Phi_G,SigmaH_G,SigmaP_G, &
       Jr_G,Potential_G,Br_G,Btheta_G,  &
       BmagnitudeSquared_G, &
       Ephi_C,Etheta_C,Er_C, uExBtheta_C,&
       uExBphi_C,uExBr_C


  real, dimension(MaxLine) ::            &
       GeoMagLat_I,GeoMagLon_I,          &
       ThetaLine_I, PhiLine_I,           &
       xLine_I,yLine_I,zLine_I,          &
       xLineOld_I,yLineOld_I,zLineOld_I, &
       UthetaLine_I,UphiLine_I,          &
       UxLine_I,UyLine_I,UzLine_I,       &
       OmegaLine_I,                      &
       JrLine_I
  
  integer :: nStep=0
  
  integer, dimension(maxLine)        ::  iThetaLine_I,iPhiLine_I
  real   ::  DtHorizontalOrig = 50.0, DtHorizontal=50.0, Time, TimeMax

  logical::  DoMoveLine=.true., UseJr=.true., UseCentrifugal=.true.
  logical::  UseIE=.false.
  logical::  DoPlotElectrodynamics=.false.
  logical::  DoSavePlot=.true.
  character(len=100) :: NamePhiNorth, NamePhiSouth

  character(len=100) :: NameInput,  &
                   NameCollision,NameSourceGraphics

  character(len=100),dimension(MaxLine):: &
       NameRestartIn, NameRestart, NameGraphics,NameOutput

  integer       :: iUnitInput,iUnitSourceGraphics,&
                   iUnitCollision,nAlt=390
                 
  integer, dimension(maxLine) :: &
       iUnitRestart,iUnitRestartIn,iUnitGraphics,iUnitOutput
  
  real :: r_C(1:maxGrid)
  real :: State_CVI(maxGrid,nVar,maxLine)
  real, dimension(1:maxGrid,0:maxLine):: &
       dOxyg_CI, uOxyg_CI, pOxyg_CI, TOxyg,     &
       dHel_CI, uHel_CI, pHel_CI, THel,         &
       dHyd_CI, uHyd_CI, pHyd_CI, THyd,         &
       dElect_CI, uElect_CI, pElect_CI, TElect
  real :: DToutput=50.0, DtVertical=0.05, Tmax=100.0,DtPlotElectrodynamics=10.0

  logical:: &
       IsFullyImplicit    = .false.,  &
       IsPointImplicit    = .false.,  & 
       IsPointImplicitAll = .false.,  & ! Including ion-ion friction
       IsRestart          = .true.,   &
       IsVariableDt       = .true.

  character(7) :: TypeSolver='Godunov'
  
  real ::  Beta = 1.0              ! 1 <= Beta <= 2

contains
  !==========================================================================
  subroutine allocate_ie_variables(iSize, jSize)
    integer, intent(in):: iSize, jSize

    nPhi = iSize
    nTheta = jSize

    allocate( &
         Theta_G(0:iSize+1, 0:jSize+1), &
         Phi_G(0:iSize+1, 0:jSize+1), &
         Potential_G(0:iSize+1, 0:jSize+1), &
         Jr_G(0:iSize+1, 0:jSize+1), &
         SigmaH_G(0:iSize+1, 0:jSize+1), &
         SigmaP_G(0:iSize+1, 0:jSize+1), &
         Br_G(0:iSize+1, 0:jSize+1), &
         Btheta_G(0:iSize+1, 0:jSize+1),  &
         BmagnitudeSquared_G(0:iSize+1, 0:jSize+1), &
         Ephi_C(iSize, jSize), &
         Etheta_C(iSize, jSize), &
         Er_C(iSize, jSize), & 
         uExBtheta_C(iSize, jSize), &
         uExBphi_C(iSize, jSize), &
         uExBr_C(iSize, jSize))

  end subroutine allocate_ie_variables

end module ModPWOM
