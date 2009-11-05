subroutine get_weimer_potential
  use ModIndicesInterfaces
  use ModPwTime,  ONLY:CurrentTime,StartTime
  use ModNumConst,ONLY:cPi,cTwoPi,cHalfPi,cRadToDeg
  use ModPWOM,    ONLY:Theta_G,Phi_G,nTheta,nPhi,&
       SigmaH_G,SigmaP_G,Jr_G,Potential_G,&
       Time,allocate_ie_variables
  implicit none
  character (len=100), dimension(100):: Lines_I
  integer :: iError,iPhi,iTheta
  real    :: temp, dTheta, dPhi
  logical,save :: UseIMF, IsFirst=.true.
  !real :: MLT_C(257,65), MLatitude_C(257,65),TempPotential_C(257,65)
  real,save :: MLT_C(360,90), MLatitude_C(360,90)
  !----------------------------------------------------------------------------

  call allocate_ie_variables(360, 90)

  CurrentTime=StartTime+Time

  if (IsFirst) then 
     !Setup Theta and Phi Grids
     dTheta=cHalfPi/(real(nTheta)-1.0)
     dPhi  =cTwoPi /(real(nPhi)-1.0)

     Theta_G(:,1) = 0.0
     Phi_G  (1,:) = 0.0
     do iPhi=2,nPhi
        Phi_G(iPhi,:)=mod(Phi_G(iPhi-1,:) + dPhi,cTwoPi)
     enddo

     do iTheta=2,nTheta
        Theta_G(:,iTheta)=Theta_G(:,iTheta-1)+dTheta
     enddo

     MLT_C(1:nPhi,1:nTheta)=mod(Phi_G(1:nPhi,1:nTheta)*12.0/cPi+12.0,24.0)
     MLatitude_C(1:nPhi,1:nTheta)=(cHalfPi-Theta_G(1:nPhi,1:nTheta))*cRadToDeg

     Lines_I(1) = "#BACKGROUND"
     Lines_I(2) = "EIE/"


     call get_IMF_Bz(CurrentTime, temp, iError)
     
     call IO_SetIMFBz(temp)
     if (iError /= 0) then
        call con_stop("PW_ERROR:Can not find IMF Bz.")
     else
        write(*,*) "Setting potential to Weimer [1996]."
        Lines_I(3) = "weimer96"    ! Change to "zero" if you want
     endif
     Lines_I(4) = "ihp"
     Lines_I(5) = "idontknow"
     Lines_I(6) = ""
     Lines_I(7) = "#DEBUG"
     Lines_I(8) = "0"
     Lines_I(9) = "0"
     Lines_I(10) = ""
     Lines_I(11) = "#END"


     call EIE_set_inputs(Lines_I)

     call EIE_Initialize(iError)
     if (iError /= 0) call &
          con_stop('PW_ERROR: EIE_Initialize failed at get_weimer_potential')

     IsFirst = .false.
  endif

  call IO_SetnMLTs(nPhi)
  call IO_SetnLats(nTheta)
  call IO_SetTime(CurrentTime)

  call IO_SetNorth


  call get_IMF_Bz(CurrentTime, temp, iError)
  call IO_SetIMFBz(temp)

  call get_IMF_By(CurrentTime, temp, iError)
  call IO_SetIMFBy(temp)

  call get_SW_V(CurrentTime, temp, iError)
  call IO_SetSWV(temp)
  if (iError /= 0) call &
       con_stop('PW_ERROR: get_SW_V failed in get_weimer_potential')

  !  call get_kp(CurrentTime, temp, iError)
  !  call IO_Setkp(temp)
  !  call IO_Setkp(1.0)
  !  if (iError /= 0) &
  !      call con_stop('PW_ERROR: get_kp failed in get_weimer_potential')


  call IO_SetGrid(                    &
       MLT_C(1:nPhi,1:nTheta), &
       MLatitude_C(1:nPhi,1:nTheta), iError)
  if (iError /= 0) call &
       con_stop('PW_ERROR: IO_SetGrid failed in get_weimer_potential')

  call IO_GetPotential(Potential_G(1:nPhi,1:nTheta),iError)

  if (iError /= 0) then
     write(*,*) "Error in get_potential (IO_GetPotential):"
     write(*,*) iError
     call con_stop("Stopping in get_weimer_potential")
  endif


  !quantities not set by weimer are set to zero
  SigmaH_G(:,:)=0.0
  SigmaP_G(:,:)=0.0
  Jr_G    (:,:)=0.0

  ! Set Values for Aurora


end subroutine get_weimer_potential
