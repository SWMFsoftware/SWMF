program dipole11

  ! Create a magnetogram with spherical harmonics (1,1)
  ! proportional to sin(Theta)*cos(Phi)

  use ModNumConst, ONLY: cTwoPi

  implicit none

  real,parameter::Amplitude = 5.0
  integer, parameter:: nPhi = 360, nTheta = 180
  real, parameter:: dPhi= cTwoPi/nPhi, dSinTheta = 2.0/nTheta

  real :: Theta
  integer:: iPhi, iTheta
  !--------------------------------------------------------------

  open(11,file='fitsfile.dat',status='replace')
  write(11,'(a)')'#CR'
  write(11,'(a)')'1200'
  write(11,'(a)')'#nMax'
  write(11,'(a)')'90'
  write(11,'(a)')'#ARRAYSIZE'
  write(11,*)nPhi
  write(11,*)nTheta
  write(11,'(a)')'#START'

  do iTheta = 0, nTheta-1
     Theta = acos( (iTheta+0.5)*dSinTheta - 1.0 )
     do iPhi = 0, nPhi-1
        write(11,'(e15.6)')Amplitude * sin(Theta) *cos(iPhi*dPhi)
     end do
  end do
  close(11)
  stop

end program dipole11
!=================================================================
subroutine CON_stop(String)
  character(LEN=*),intent(in):: String
  write(*,*) String
  stop
end subroutine CON_stop
