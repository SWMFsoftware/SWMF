program dipole
  use ModMagHarmonics
  implicit none

 
  real,parameter::Amplitude = 5.0
  nPhi   = 360
  nTheta = 180
  dPhi   = 2.0*cPi/nPhi
  dTheta =     cPi/nTheta
  dSinTheta =  2.0/nTheta
  open(11,file='fitsfile.dat',status='replace')
  write(11,'(a)')'#CR'
  write(11,'(a)')'1200'
  write(11,'(a)')'#nMax'
  write(11,'(a)')'90'
  write(11,'(a)')'#ARRAYSIZE'
  write(11,*)nPhi
  write(11,*)nTheta
  write(11,'(a)')'#START'
  do iTheta=0,nTheta-1
     Theta = colatitude(iTheta)
     write(*,*) Theta*180/cPi
     do iPhi = 0,nPhi-1
       ! write(11,'(e15.6)')Amplitude * cos(Theta)
        write(11,'(e15.6)')Amplitude * sin(Theta) *cos(iPhi*dPhi)
     end do
  end do
  close(11)
  stop
end program dipole
subroutine CON_stop(TypeMessage)
  character(LEN=*),intent(in)::TypeMessage
  stop
end subroutine CON_stop
