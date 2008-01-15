
subroutine conductance_gradients

  use ModSizeRIM
  use ModRIM
  use ModNumConst, only : cPi

  implicit none

  real, dimension(0:nLons+1,nLats) :: &
       sn, cs, sn2, cs2, cs3, cs4, Theta, C

  real, dimension(0:nLons+1,nLats) :: &
       TermTheta2, TermTheta1, TermPsi2, TermPsi1

  real, dimension(0:nLons+1,nLats) :: &
       dLatitude2, dLongitude2

  Theta = cPi - Latitude

  sn = sin(Theta)
  cs = cos(Theta)
  sn2= sn*sn
  cs2 = cs*cs
  cs3 = 1.00 + 3.00*cs2
  cs4 = sqrt(cs3)
  C = 4.00*Sigma0*cs2 + SigmaP*sn2

  SigmaThTh = Sigma0*SigmaP*cs3/C
  SigmaThPs = 2.00*Sigma0*SigmaH*cs*cs4/C
  SigmaPsPs = SigmaP+SigmaH*SigmaH*sn2/C

  ! Edge 1

  dSigmaThTh_dLatitude(0:nLons+1,1) = &
       (SigmaThTh(0:nLons+1,2)-SigmaThTh(0:nLons+1,1))/ &
       (dLatitude(0:nLons+1,1))

  dSigmaThTh_dLongitude(0,1:nLats) = &
       (SigmaThTh(1,1:nLats)-SigmaThTh(0,1:nLats))/ &
       (dLongitude(1,1:nLats))

  dSigmaThPs_dLatitude(0:nLons+1,1) = &
       (SigmaThPs(0:nLons+1,2)-SigmaThPs(0:nLons+1,1))/ &
       (dLatitude(0:nLons+1,1))
  dSigmaThPs_dLongitude(0,1:nLats) = &
       (SigmaThPs(1,1:nLats)-SigmaThPs(0,1:nLats))/ &
       (dLongitude(1,1:nLats))

  dSigmaPsPs_dLatitude(0:nLons+1,1) = &
       (SigmaPsPs(0:nLons+1,2)-SigmaPsPs(0:nLons+1,1))/ &
       (dLatitude(0:nLons+1,1))
  dSigmaPsPs_dLongitude(0,1:nLats) = &
       (SigmaPsPs(1,1:nLats)-SigmaPsPs(0,1:nLats))/ &
       (dLongitude(1,1:nLats))

  ! Middle

  dSigmaThTh_dLatitude(0:nLons+1,2:nLats-1) = &
       (SigmaThTh(0:nLons+1,3:nLats)-SigmaThTh(0:nLons+1,1:nLats-2))/ &
       (dLatitude(0:nLons+1,2:nLats-1))

  dSigmaThTh_dLongitude(1:nLons,1:nLats) = &
       (SigmaThTh(2:nLons+1,1:nLats)-SigmaThTh(0:nLons-1,1:nLats))/ &
       (dLongitude(1:nLons,1:nLats))

  dSigmaThPs_dLatitude(0:nLons+1,2:nLats-1) = &
       (SigmaThPs(0:nLons+1,3:nLats)-SigmaThPs(0:nLons+1,1:nLats-2))/ &
       (dLatitude(0:nLons+1,2:nLats-1))
  dSigmaThPs_dLongitude(1:nLons,1:nLats) = &
       (SigmaThPs(2:nLons+1,1:nLats)-SigmaThPs(0:nLons-1,1:nLats))/ &
       (dLongitude(1:nLons,1:nLats))

  dSigmaPsPs_dLatitude(0:nLons+1,2:nLats-1) = &
       (SigmaPsPs(0:nLons+1,3:nLats)-SigmaPsPs(0:nLons+1,1:nLats-2))/ &
       (dLatitude(0:nLons+1,2:nLats-1))
  dSigmaPsPs_dLongitude(1:nLons,1:nLats) = &
       (SigmaPsPs(2:nLons+1,1:nLats)-SigmaPsPs(0:nLons-1,1:nLats))/ &
       (dLongitude(1:nLons,1:nLats))

  ! Edge 2

  dSigmaThTh_dLatitude(0:nLons+1,nLats) = &
       (SigmaThTh(0:nLons+1,nLats)-SigmaThTh(0:nLons+1,nLats-1))/ &
       (dLatitude(0:nLons+1,nLats))
  dSigmaThTh_dLongitude(nLons+1,1:nLats) = &
       (SigmaThTh(nLons+1,1:nLats)-SigmaThTh(nLons,1:nLats))/ &
       (dLongitude(nLons+1,1:nLats))

  dSigmaThPs_dLatitude(0:nLons+1,nLats) = &
       (SigmaThPs(0:nLons+1,nLats)-SigmaThPs(0:nLons+1,nLats-1))/ &
       (dLatitude(0:nLons+1,nLats))
  dSigmaThPs_dLongitude(nLons+1,1:nLats) = &
       (SigmaThPs(nLons+1,1:nLats)-SigmaThPs(nLons,1:nLats))/ &
       (dLongitude(nLons+1,1:nLats))

  dSigmaPsPs_dLatitude(0:nLons+1,nLats) = &
       (SigmaPsPs(0:nLons+1,nLats)-SigmaPsPs(0:nLons+1,nLats-1))/ &
       (dLatitude(0:nLons+1,nLats))
  dSigmaPsPs_dLongitude(nLons+1,1:nLats) = &
       (SigmaPsPs(nLons+1,1:nLats)-SigmaPsPs(nLons,1:nLats))/ &
       (dLongitude(nLons+1,1:nLats))

  dLatitude2          = (dLatitude/2       )*(dLatitude/2       )
  dLatitude2(:,    1) = (dLatitude(:,    1))*(dLatitude(:,    1))
  dLatitude2(:,nLats) = (dLatitude(:,nLats))*(dLatitude(:,nLats))
  
  dLongitude2            = (dLongitude/2         )*(dLongitude/2         )
  dLongitude2(      0,:) = (dLongitude(      0,:))*(dLongitude(      0,:))
  dLongitude2(nLons+1,:) = (dLongitude(nLons+1,:))*(dLongitude(nLons+1,:))
  
  TermTheta2 = SigmaThTh*sn2/dLatitude2
  TermTheta1 = ( SigmaThTh*sn*cs   &
       + dSigmaThTh_dLatitude*sn2 &
       - dSigmaThPs_dLongitude*sn) / dLatitude
  TermPsi2 = SigmaPsPs / dLongitude2
  TermPsi1 = (dSigmaThPs_dLatitude*sn + dSigmaPsPs_dLongitude) / dLongitude

  ! Form the complete matrix
  SolverA =  -2.0 * (TermTheta2 + TermPsi2)
  SolverB =          TermTheta2 - TermTheta1
  SolverC =          TermTheta2 + TermTheta1
  SolverD =          TermPsi2   - TermPsi1
  SolverE =          TermPsi2   + TermPsi1

!  SolverA = 4.0
!  SolverB = 1.0
!  SolverC = 1.0
!  SolverD = 1.0
!  SolverE = 1.0

end subroutine conductance_gradients
