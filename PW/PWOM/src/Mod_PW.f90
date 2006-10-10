Module Mod_PW
  use ModParameters
  integer ::   nTheta, nPhi,maxLine
  Parameter (nTheta     = 65 )
  Parameter (nPhi       = 257)
  Parameter (maxLine    = 1  )
  
  integer ::   iTheta, iPhi, iUnitSouth,iUnitNorth,i,iLine
  integer, dimension(maxLine)::iLineGlobal

  real    ::   Bcoef,MagMoment,rPlanet,Dtheta,Dphi,rLowerBoundary,nLine

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
  real                               ::  Dt, Time, maxTime

  Logical                            ::  IsMoveFluxTube,IsUseJr,IsCentrifugal
  
  Character*100 :: NamePhiNorth, NamePhiSouth

  CHARACTER*100 :: NameInput, NameOutput, &
                   NameCollision,NameSourceGraphics

  CHARACTER*100,dimension(maxline) ::NameRestartIn, NameRestart,NameGraphics

  integer       :: iUnitInput,iUnitOutput,iUnitSourceGraphics,&
                   iUnitCollision,nDim
                 
  integer, dimension(maxLine) :: iUnitRestart,iUnitRestartIn,iUnitGraphics
  
  real, dimension(1:maxGrid,0:maxLine):: dOxyg, uOxyg, pOxyg, TOxyg,     &
                                          dHel, uHel, pHel, THel,         &
                                          dHyd, uHyd, pHyd, THyd,         &
                                          dElect, uElect, pElect, TElect
  real   :: DToutput,DTpolarwind, Tmax

  logical:: IsImplicit, IsRestart, IsVariableDt

  CHARACTER(7) :: TypeSolver
end Module Mod_PW
