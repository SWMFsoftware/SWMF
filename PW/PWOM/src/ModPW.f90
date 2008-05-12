module ModPWOM

  use ModParameters
  use ModCommonPlanet, ONLY: nVar
  implicit none

  logical :: IsStandAlone = .false.
  logical :: UseIonHeat=.true.,UseEleHeat=.true.,UseExplicitHeat=.false.
  logical :: UseIndicies
  integer :: iUnitOut
  character (len=7) :: StringPrefix=''

  integer :: iProc, nProc, iComm,errcode

  integer :: iProcTest = 0, iLineTest = 1
  character (len=100):: StringTest = ''

  integer, parameter :: MaxLine = 500
  
  integer :: nTotalLine=1
  integer :: iTheta, iPhi, iUnitSouth,iUnitNorth,i,iLine=0
  integer, dimension(:), allocatable::iLineGlobal

  real    ::   Bcoef,MagMoment,rPlanet,OmegaPlanet,Dtheta,Dphi,rLowerBoundary

  integer :: nLine,nLog=0

  ! The number of lines on each processor and on the processors with lower rank
  integer, allocatable :: nLine_P(:), nLineBefore_P(:)

  !Joule Heating parameters
  real    :: uJoule2
  logical :: UseJouleHeating=.false.

  ! ionosphere variables
  integer :: nTheta = -1, nPhi = -1
  real, dimension(:,:), allocatable :: &
       Theta_G,Phi_G,SigmaH_G,SigmaP_G, &
       Jr_G,Potential_G,Br_G,Btheta_G,  &
       BmagnitudeSquared_G, &
       Ephi_C,Etheta_C,Er_C, uExBtheta_C,&
       uExBphi_C,uExBr_C,ElectronEnergyFlux_C,&
       ElectronAverageEnergy_C


  real, dimension(:), allocatable ::            &
       GeoMagLat_I,GeoMagLon_I,          &
       ThetaLine_I, PhiLine_I,           &
       xLine_I,yLine_I,zLine_I,          &
       xLineOld_I,yLineOld_I,zLineOld_I, &
       UthetaLine_I,UphiLine_I,          &
       UxLine_I,UyLine_I,UzLine_I,       &
       OmegaLine_I,                      &
       JrLine_I
  
  integer :: nStep=0
  
  integer, dimension(:), allocatable        ::  iThetaLine_I,iPhiLine_I
  real   ::  DtHorizontalOrig = 50.0, DtHorizontal=50.0, Time, TimeMax
  real, allocatable   ::  Dt_I(:)
  logical::  DoMoveLine=.true., UseJr=.true., UseCentrifugal=.true.
  logical::  UseIE=.false.,UseAurora=.false.
  logical::  DoPlotElectrodynamics=.false.
  logical::  DoSavePlot=.true.
  character(len=100) :: NamePhiNorth, NamePhiSouth

  character(len=100) :: NameInput,  &
                   NameCollision,NameSourceGraphics

  character(len=100),dimension(:), allocatable:: &
       NameRestartIn, NameRestart, NameGraphics,NameOutput

  integer       :: iUnitInput,iUnitSourceGraphics,&
                   iUnitCollision,nAlt=390
                 
  integer, dimension(:), allocatable :: &
       iUnitRestart,iUnitRestartIn,iUnitGraphics,iUnitOutput
  
  real,allocatable :: r_C(:)
  real,allocatable :: State_CVI(:,:,:)

  real :: DToutput=50.0, DtVertical=0.05, Tmax=100.0,DtPlotElectrodynamics=10.0
  integer :: MaxStep = -1, DnOutput=-1
  logical :: DoTimeAccurate = .true.

  logical:: &
       IsFullyImplicit    = .false.,  &
       IsPointImplicit    = .false.,  & 
       IsPointImplicitAll = .false.,  & ! Including ion-ion friction
       IsRestart          = .true.,   &
       IsVariableDt       = .false.

  character(7)   :: TypeSolver='Godunov'
  character(13)  :: TypeDiffusion='LaxFriedrichs'
  real ::  BetaIn = 1.0  ! limiter beta: 1 <= Beta <= 2, 0 for first order
  real ::  Beta = 1.0    ! actual beta used (changes in implicit scheme) 

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
