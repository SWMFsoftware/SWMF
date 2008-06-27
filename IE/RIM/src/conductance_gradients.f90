
subroutine conductance_gradients

  use ModSizeRIM
  use ModRIM
  use ModNumConst, only : cPi
  use ModParamRIM

  implicit none

  real, dimension(0:nLons+1,nLats) :: &
       sn, cs, sn2, cs2, cs3, cs4, Theta, C

  real, dimension(0:nLons+1,nLats) :: &
       TermTheta2, TermTheta1, TermPsi2, TermPsi1

  real, dimension(0:nLons+1,nLats) :: &
       dLongitude2, dTheta, dTheta2

  real :: r
  integer :: iLon, iLat, iLh

  Theta = cPi/2 - Latitude

  sn = sin(Theta)
  cs = cos(Theta)

  sn2 = sn*sn
  cs2 = cs*cs
  cs3 = 1.00 + 3.00*cs2
  cs4 = sqrt(cs3)
  C   = 4.00*Sigma0*cs2 + SigmaP*sn2

  if (DoFold) then 

     write(*,*) "Folding conductances!"

     do iLon = 1, nLons
        do iLat = 1, nLats
           if (Latitude(iLon,iLat) < -OCFLB(iLon)-OCFLBBuffer .or. &
                Latitude(iLon,iLat) > OCFLB(iLon)+OCFLBBuffer) then
              SigmaThTh(iLon,iLat) = &
                   Sigma0(iLon,iLat)*SigmaP(iLon,iLat)*cs3(iLon,iLat)/&
                   C(iLon,iLat)
              SigmaThPs(iLon,iLat) = &
                   2.00*Sigma0(iLon,iLat)*SigmaH(iLon,iLat)*&
                   cs(iLon,iLat)*cs4(iLon,iLat)/C(iLon,iLat)
              SigmaPsPs(iLon,iLat) = SigmaP(iLon,iLat)+ &
                   SigmaH(iLon,iLat)*SigmaH(iLon,iLat)*sn2(iLon,iLat)/&
                   C(iLon,iLat)
           else

              ! linearly change from the OCFLB through the buffer region
              if (Latitude(iLon,iLat) > -OCFLB(iLon) .and. &
                   Latitude(iLon,iLat) < OCFLB(iLon)) then
                 r = 1.0
              else
                 r = 1-(abs(Latitude(iLon,iLat))-OCFLB(iLon))/OCFLBBuffer
              endif

              ! iLh = latitude of other hemisphere
              iLh = nLats - iLat + 1

              ! We add the conductance from one hemisphere to the other
              ! hemisphere.  Here we simply add them together, while you
              ! could imagine having the magnetic field geometry incorporated
              ! into this somehow....

              SigmaThTh(iLon,iLat) = &
                   (Sigma0(iLon,iLat)+r*Sigma0(iLon,iLh)) * &
                   (SigmaP(iLon,iLat)+r*SigmaP(iLon,iLh))*cs3(iLon,iLat)/&
                   (C(iLon,iLat)+r*C(iLon,iLh))

              SigmaThPs(iLon,iLat) = &
                   2.0*(Sigma0(iLon,iLat)+r*Sigma0(iLon,iLh))* &
                   (SigmaH(iLon,iLat)+r*SigmaH(iLon,iLh))*&
                   cs(iLon,iLat)*cs4(iLon,iLat)/&
                   (C(iLon,iLat)+r*C(iLon,iLh))

              SigmaPsPs(iLon,iLat) = &
                   (SigmaP(iLon,iLat) + r*SigmaP(iLon,iLh)) + &
                   (SigmaH(iLon,iLat) + r*SigmaH(iLon,iLh))**2 * &
                   sn2(iLon,iLat)/&
                   (C(iLon,iLat)+r*C(iLon,iLh))

           endif
        enddo
     enddo

  else

     SigmaThTh = Sigma0*SigmaP*cs3/C
     SigmaThPs = 2.00*Sigma0*SigmaH*cs*cs4/C
     SigmaPsPs = SigmaP+SigmaH*SigmaH*sn2/C

  endif

  ! Edge 1

  dTheta = -dLatitude

  dSigmaThTh_dLatitude(0:nLons+1,1) = &
       (SigmaThTh(0:nLons+1,2)-SigmaThTh(0:nLons+1,1))/ &
       (dTheta(0:nLons+1,1))

  dSigmaThTh_dLongitude(0,1:nLats) = &
       (SigmaThTh(1,1:nLats)-SigmaThTh(0,1:nLats))/ &
       (dLongitude(1,1:nLats))

  dSigmaThPs_dLatitude(0:nLons+1,1) = &
       (SigmaThPs(0:nLons+1,2)-SigmaThPs(0:nLons+1,1))/ &
       (dTheta(0:nLons+1,1))
  dSigmaThPs_dLongitude(0,1:nLats) = &
       (SigmaThPs(1,1:nLats)-SigmaThPs(0,1:nLats))/ &
       (dLongitude(1,1:nLats))

  dSigmaPsPs_dLatitude(0:nLons+1,1) = &
       (SigmaPsPs(0:nLons+1,2)-SigmaPsPs(0:nLons+1,1))/ &
       (dTheta(0:nLons+1,1))
  dSigmaPsPs_dLongitude(0,1:nLats) = &
       (SigmaPsPs(1,1:nLats)-SigmaPsPs(0,1:nLats))/ &
       (dLongitude(1,1:nLats))

  ! Middle

  dSigmaThTh_dLatitude(0:nLons+1,2:nLats-1) = &
       (SigmaThTh(0:nLons+1,3:nLats)-SigmaThTh(0:nLons+1,1:nLats-2))/ &
       (dTheta(0:nLons+1,2:nLats-1))

  dSigmaThTh_dLongitude(1:nLons,1:nLats) = &
       (SigmaThTh(2:nLons+1,1:nLats)-SigmaThTh(0:nLons-1,1:nLats))/ &
       (dLongitude(1:nLons,1:nLats))

  dSigmaThPs_dLatitude(0:nLons+1,2:nLats-1) = &
       (SigmaThPs(0:nLons+1,3:nLats)-SigmaThPs(0:nLons+1,1:nLats-2))/ &
       (dTheta(0:nLons+1,2:nLats-1))
  dSigmaThPs_dLongitude(1:nLons,1:nLats) = &
       (SigmaThPs(2:nLons+1,1:nLats)-SigmaThPs(0:nLons-1,1:nLats))/ &
       (dLongitude(1:nLons,1:nLats))

  dSigmaPsPs_dLatitude(0:nLons+1,2:nLats-1) = &
       (SigmaPsPs(0:nLons+1,3:nLats)-SigmaPsPs(0:nLons+1,1:nLats-2))/ &
       (dTheta(0:nLons+1,2:nLats-1))
  dSigmaPsPs_dLongitude(1:nLons,1:nLats) = &
       (SigmaPsPs(2:nLons+1,1:nLats)-SigmaPsPs(0:nLons-1,1:nLats))/ &
       (dLongitude(1:nLons,1:nLats))

  ! Edge 2

  dSigmaThTh_dLatitude(0:nLons+1,nLats) = &
       (SigmaThTh(0:nLons+1,nLats)-SigmaThTh(0:nLons+1,nLats-1))/ &
       (dTheta(0:nLons+1,nLats))
  dSigmaThTh_dLongitude(nLons+1,1:nLats) = &
       (SigmaThTh(nLons+1,1:nLats)-SigmaThTh(nLons,1:nLats))/ &
       (dLongitude(nLons+1,1:nLats))

  dSigmaThPs_dLatitude(0:nLons+1,nLats) = &
       (SigmaThPs(0:nLons+1,nLats)-SigmaThPs(0:nLons+1,nLats-1))/ &
       (dTheta(0:nLons+1,nLats))
  dSigmaThPs_dLongitude(nLons+1,1:nLats) = &
       (SigmaThPs(nLons+1,1:nLats)-SigmaThPs(nLons,1:nLats))/ &
       (dLongitude(nLons+1,1:nLats))

  dSigmaPsPs_dLatitude(0:nLons+1,nLats) = &
       (SigmaPsPs(0:nLons+1,nLats)-SigmaPsPs(0:nLons+1,nLats-1))/ &
       (dTheta(0:nLons+1,nLats))
  dSigmaPsPs_dLongitude(nLons+1,1:nLats) = &
       (SigmaPsPs(nLons+1,1:nLats)-SigmaPsPs(nLons,1:nLats))/ &
       (dLongitude(nLons+1,1:nLats))

  dTheta2          = (dTheta/2       )*(dTheta/2       )
  dTheta2(:,    1) = (dTheta(:,    1))*(dTheta(:,    1))
  dTheta2(:,nLats) = (dTheta(:,nLats))*(dTheta(:,nLats))
  
  dLongitude2            = (dLongitude/2         )*(dLongitude/2         )
  dLongitude2(      0,:) = (dLongitude(      0,:))*(dLongitude(      0,:))
  dLongitude2(nLons+1,:) = (dLongitude(nLons+1,:))*(dLongitude(nLons+1,:))
  

  TermTheta2 = SigmaThTh*sn2/dTheta2
  TermTheta1 = ( SigmaThTh*sn*cs   &
       + dSigmaThTh_dLatitude*sn2 &
       - dSigmaThPs_dLongitude*sn) / dTheta
  TermPsi2 = SigmaPsPs / dLongitude2
  TermPsi1 = (dSigmaThPs_dLatitude*sn + dSigmaPsPs_dLongitude) / dLongitude

  ! Form the complete matrix
  SolverA =  -2.0 * (TermTheta2 + TermPsi2)
  SolverB =          TermTheta2 - TermTheta1
  SolverC =          TermTheta2 + TermTheta1
  SolverD =          TermPsi2   - TermPsi1
  SolverE =          TermPsi2   + TermPsi1

end subroutine conductance_gradients
