module ModPWOM

  use ModParameters
  implicit none

  logical :: IsStandAlone = .false.

  integer :: iUnitOut
  character (len=7) :: StringPrefix=''

  integer :: iProc, nProc, iComm,errcode

  integer, parameter :: MaxLine = 500
  
  integer :: nTotalLine=1
  integer ::   iTheta, iPhi, iUnitSouth,iUnitNorth,i,iLine
  integer, dimension(MaxLine)::iLineGlobal

  real    ::   Bcoef,MagMoment,rPlanet,Dtheta,Dphi,rLowerBoundary

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
  
  
  integer, dimension(maxLine)        ::  iThetaLine_I,iPhiLine_I
  real                               ::  DtHorizontal=50.0, Time, TimeMax

  logical::  DoMoveLine=.true., UseJr=.true., UseCentrifugal=.true.
  logical::  UseIE=.false.
  logical::  DoPlotElectrodynamics=.false.
  character(len=100) :: NamePhiNorth, NamePhiSouth

  character(len=100) :: NameInput,  &
                   NameCollision,NameSourceGraphics

  character(len=100),dimension(MaxLine):: &
       NameRestartIn, NameRestart, NameGraphics,NameOutput

  integer       :: iUnitInput,iUnitSourceGraphics,&
                   iUnitCollision,nAlt=390
                 
  integer, dimension(maxLine) :: iUnitRestart,iUnitRestartIn,iUnitGraphics,iUnitOutput
  
  real, dimension(1:maxGrid,0:maxLine):: dOxyg_CI, uOxyg_CI, pOxyg_CI, TOxyg,     &
                                          dHel_CI, uHel_CI, pHel_CI, THel,         &
                                          dHyd_CI, uHyd_CI, pHyd_CI, THyd,         &
                                          dElect_CI, uElect_CI, pElect_CI, TElect
  real :: DToutput=10.0, DtVertical=0.05, Tmax=100.0,DtPlotElectrodynamics=10.0

  logical:: IsImplicit=.false., IsRestart=.true., IsVariableDt=.true.

  character(7) :: TypeSolver='Godunov'

contains

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
