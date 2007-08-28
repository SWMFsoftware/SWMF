subroutine get_weimer_potential
  use ModIndicesInterfaces
  use ModInputs
  use ModPwTime,ONLY:CurrentTime,StartTime
  use ModNumConst,ONLY:cPi,cTwoPi,cHalfPi
  character (len=100), dimension(100):: Lines
  !----------------------------------------------------------------------------

  call allocate_ie_variables(257, 65)
  
  CurrentTime=StartTime+Time
  
  !Setup Theta and Phi Grids
  dTheta=cHalfPi/nTheta
  dPhi  =cTwoPi /nPhi
  
  Theta_G(:,1) = 0.0
  Phi_G  (1,:) = 0.0
  do iPhi=1,nPhi
     do iTheta=1,nTheta
        Theta_G(iPhi,iTheta)=    Theta_G(iPhi,iTheta-1)+dTheta
        Phi_G  (iPhi,iTheta)=mod(Phi_G(iPhi-1,iTheta)  +dPhi,cTwoPi)
     enddo
  enddo

  MLT_C(1:nPhi,1:nTheta)=Phi_G(1:nPhi,1:nTheta)*12.0/cPi
  MLatitude_C(1:nPhi,1:nTheta)=cHalfPi-Theta_G(1:nPhi,1:nTheta)
  
  Lines(1) = "#BACKGROUND"
  Lines(2) = "PW/"
  
  UseHPI = .true.
  call get_IMF_Bz(CurrentTime, temp, iError)
  call IO_SetIMFBz(temp)
  write(*,*) "Setting potential to Weimer [1996]."
  Lines(3) = "weimer96"    ! Change to "zero" if you want
  UseIMF = .true.
  Lines(4) = "ihp"
  Lines(5) = "idontknow"
  Lines(6) = ""
  
  call EIE_set_inputs(Lines)
  
  call EIE_Initialize(iError)

  call get_IMF_Bz(CurrentTime, temp, iError)
  call IO_SetIMFBz(temp)

  call get_IMF_By(CurrentTime, temp, iError)
  call IO_SetIMFBy(temp)
  
  call get_SW_V(CurrentTime, temp, iError)
  call IO_SetSWV(temp)

  call IO_SetGrid(                    &
       MLT(1:nPhi,1:nTheta), &
       MLatitude(1:nPhi,1:nTheta), iError)
  
  call IO_GetPotential(Potential_G(1:nPhi,1:nTheta), iError)
  
end subroutine get_weimer_potential
