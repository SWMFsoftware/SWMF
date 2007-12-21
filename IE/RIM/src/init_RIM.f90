subroutine init_RIM()
  use ModProcIE
  use ModRIM
  use ModNumConst

  implicit none
  integer :: i


  ! Set actual grid values
  do i=1,nLats
     Latitude(:,i)=-cHalfPi+(cPi/nLats)*(real(i)-cHalf)
  end do
  do i=0,nLons+1
     Longitude(i,:)=((cTwoPi/nLons)/real(nProc))*(iProc*nLons+real(i)-cHalf) 
  end do
  allocate(AllLons(nLons*nProc))
  do i=1,nLons*nProc
     AllLons(i)=((cTwoPi/nLons)/real(nProc))*(real(i)-cHalf) 
  end do

  dLatitude(:,1)=(Latitude(:,2)-Latitude(:,1))
  do i=2,nLats-1
     dLatitude(:,i)=(Latitude(:,i+1)-Latitude(:,i-1))*cHalf
  end do
  dLatitude(:,nLats)=(Latitude(:,nLats)-Latitude(:,nLats-1))

  dLongitude(0,:)=(Longitude(1,:)-Longitude(0,:))
  do i=1,nLons
     dLongitude(i,:)=(Longitude(i+1,:)-Longitude(i-1,:))*cHalf
  end do
  dLongitude(i,:)=(Longitude(nLons+1,:)-Longitude(nLons,:))

  Potential=cZero

end subroutine init_RIM
