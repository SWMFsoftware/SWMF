

subroutine get_jouleheating

  use ModRIM
  use ModParamRIM
  use ModNumConst, only : cPi

  implicit none


  integer :: iLon, iLat
  real, dimension(0:nLons+1, nLats) :: Eth, Eps, sinTheta,dTheta,dPsi

  dTheta = dLatitude
  dPsi = dLongitude
  sinTheta = sin(cPi/2 - Latitude)
  JouleHeating(:,:) = 0.0
  Eth = 0.0
  Eps = 0.0

  do iLon = 0, nLons+1
     if(iLon > 0 .and. iLon < nLons+1)then
        do iLat=2, nLats-1
           Eth(iLon,iLat) = -(Potential(iLon,iLat-1)-Potential(iLon,iLat+1))/ &
                (dTheta(iLon,iLat)*Radius)
           Eps(iLon,iLat) = -(Potential(iLon+1,iLat)-Potential(iLon-1,iLat))/ &
                (dPsi(iLon,iLat)*Radius*sinTheta(iLon,iLat))
        end do
        Eth(iLon,1) = -(Potential(iLon,1)-Potential(iLon,2))/ &
             (dTheta(iLon,1)*Radius)
        Eps(iLon,1) = Eps(iLon,2)
        Eth(iLon,nLats) = -(Potential(iLon,nLats-1)-Potential(iLon,nLats))/&
             (dTheta(iLon,nLats)*Radius)
        Eps(iLon,nLats) = Eps(iLon,nLats-1)
     else if (iLon ==0)then
        do iLat=2, nLats-1
           Eth(iLon,iLat) = -(Potential(iLon,iLat-1)-Potential(iLon,iLat+1))/ &
                (dTheta(iLon,iLat)*Radius)
           Eps(iLon,iLat) = -(Potential(iLon+1,iLat)-Potential(nLons,iLat))/ &
                (dPsi(iLon,iLat)*Radius*sinTheta(iLon,iLat))
        end do
        Eth(iLon,1) = -(Potential(iLon,1)-Potential(iLon,2))/ &
             (dTheta(iLon,1)*Radius)
        Eps(iLon,1) = Eps(iLon,2)
        Eth(iLon,nLats) = -(Potential(iLon,nLats-1)-Potential(iLon,nLats))/&
             (dTheta(iLon,nLats)*Radius)
        Eps(iLon,nLats) = Eps(iLon,nLats-1)
     else
        do iLat=2, nLats-1
           Eth(iLon,iLat) = -(Potential(iLon,iLat-1)-Potential(iLon,iLat+1))/ &
                (dTheta(iLon,iLat)*Radius)
           Eps(iLon,iLat) = -(Potential(1,iLat)-Potential(iLon,iLat))/&
                (dPsi(iLon,iLat)*Radius*sinTheta(iLon,iLat))
        end do
        Eth(iLon,1) = -(Potential(iLon,1)-Potential(iLon,2))/ &
             (dTheta(iLon,1)*Radius)
        Eps(iLon,1) = Eps(iLon,2)
        Eth(iLon,nLats) = -(Potential(iLon,nLats-1)-Potential(iLon,nLats))/&
             (dTheta(iLon,nLats)*Radius)
        Eps(iLon,nLats) = Eps(iLon,nLats-1)
     end if
  end do

  JouleHeating = SigmaP * (Eth**2+Eps**2)

end subroutine get_jouleheating
