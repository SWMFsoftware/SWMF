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

  integer :: nLine

  ! The number of lines on each processor and on the processors with lower rank
  integer, allocatable :: nLine_P(:), nLineBefore_P(:)

  ! ionosphere variables
  integer :: nTheta = -1, nPhi = -1
  real, dimension(:,:), allocatable :: &
       Theta_G,Phi_G,SigmaH_G,SigmaP_G, &
       Jr_G,Potential_G,Br_G,Btheta_G,  &
       BmagnitudeSquared_G, &
       Ephi,Etheta,Er, VelocityExBtheta,&
       VelocityExBphi,VelocityExBr,     &
       Ux,Uy,Uz,x,y,z,Ex,Ey,Ez


  real,    dimension(maxLine)        ::  FieldLineTheta, FieldLinePhi,     &
                                          FieldLineX,FieldLineY,FieldLineZ, &
                                          OldVelocityTheta,OldVelocityPhi,  &
                                          AccelerationTheta,AccelerationPhi,&
                                          FieldLineVelTheta,FieldLineVelPhi,&
                                          FieldLineVelx,FieldLineVely,      &
                                          FieldLineVelz,                    &
                                          OmegaHorFieldline,                &
                                          OldFieldLineTheta,OldFieldLinePhi,&
                                          OldFieldLineX,OldFieldLineY,      &
                                          OldFieldLineZ,GeoMagLat,GeoMagLon,&
                                          FieldLineJr
  
  
  integer, dimension(maxLine)        ::  iFieldLineTheta,iFieldLinePhi
  real                               ::  DtMax=50.0, Dt, Time, maxTime

  logical::  IsMoveFluxTube=.true., IsUseJr=.true., IsCentrifugal=.true.
  logical::  UseIE=.false.
  character(len=100) :: NamePhiNorth, NamePhiSouth

  character(len=100) :: NameInput, NameOutput, &
                   NameCollision,NameSourceGraphics

  character(len=100),dimension(MaxLine):: &
       NameRestartIn, NameRestart, NameGraphics

  integer       :: iUnitInput,iUnitOutput,iUnitSourceGraphics,&
                   iUnitCollision,nDim=390
                 
  integer, dimension(maxLine) :: iUnitRestart,iUnitRestartIn,iUnitGraphics
  
  real, dimension(1:maxGrid,0:maxLine):: dOxyg, uOxyg, pOxyg, TOxyg,     &
                                          dHel, uHel, pHel, THel,         &
                                          dHyd, uHyd, pHyd, THyd,         &
                                          dElect, uElect, pElect, TElect
  real   :: DToutput=10.0, DTpolarwind=0.05, Tmax=100.0

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
         Ephi(iSize, jSize), &
         Etheta(iSize, jSize), &
         Er(iSize, jSize), & 
         VelocityExBtheta(iSize, jSize), &
         VelocityExBphi(iSize, jSize), &
         VelocityExBr(iSize, jSize), &
         Ux(iSize, jSize), &
         Uy(iSize, jSize), &
         Uz(iSize, jSize), &
         x(iSize, jSize), &
         y(iSize, jSize), &
         z(iSize, jSize), &
         Ex(iSize, jSize), &
         Ey(iSize, jSize), &
         Ez(iSize, jSize))
  end subroutine allocate_ie_variables

end module ModPWOM
