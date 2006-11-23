module ModPWOM

  use ModParameters
  implicit none

  logical :: IsStandAlone = .false.

  integer :: iUnitOut
  character (len=7) :: StringPrefix=''

  integer :: iProc, nProc, iComm,errcode

  integer ::   nTheta, nPhi, maxLine
  parameter (nTheta     = 65 )
  parameter (nPhi       = 257)
  parameter (maxLine    = 500)
  
  integer :: nTotalLine=1
  integer ::   iTheta, iPhi, iUnitSouth,iUnitNorth,i,iLine
  integer, dimension(MaxLine)::iLineGlobal

  real    ::   Bcoef,MagMoment,rPlanet,Dtheta,Dphi,rLowerBoundary

  integer :: nLine

  real, dimension(0:nPhi+1,0:nTheta+1)::  Theta_G,Phi_G,SigmaH_G, SigmaP_G,&
                                          Jr_G,Potential_G,Br_G,Btheta_G,  &
                                          BmagnitudeSquared_G
 
  real, dimension(nPhi,nTheta)        ::  Ephi,Etheta,Er, VelocityExBtheta,&
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

end module ModPWOM
