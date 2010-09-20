!******************************************************************************
!  Determine the interpolated velocity at point theta3, and phi3. 
!  Uses Shepard Interpolation given on page 430 of 
!  Numerical Analysis by Kincaid
!******************************************************************************


subroutine interpolate_velocity(  Theta1,Theta2,Theta3, &
                                  Phi1,  Phi2,  Phi3,   &
                                  VelocityTheta11, VelocityTheta12, &
                                  VelocityTheta21, VelocityTheta22, &
                                  VelocityPhi11,   VelocityPhi12,   &
                                  VelocityPhi21,   VelocityPhi22,   &
                                  VelocityThetaOut,VelocityPhiOut)


  real, intent(in)  ::           Theta1,Theta2,Theta3, &
                                Phi1,  Phi2,  Phi3,   &
                                VelocityTheta11, VelocityTheta12, &
                                VelocityTheta21, VelocityTheta22, &
                                VelocityPhi11,   VelocityPhi12,   &
                                VelocityPhi21,   VelocityPhi22   
  
  real, intent(out)  ::          VelocityThetaOut,VelocityPhiOut
  

  real,dimension(4) ::          Theta, Phi,VelocityTheta,VelocityPhi,numerator
  
  real              ::          multiplyer
  real, dimension(4,4) ::       d

  integer           ::          i,j

  Theta(1)         = Theta1
  Phi(1)           = Phi1
  VelocityTheta(1) = VelocityTheta11
  VelocityPhi(1)   = VelocityPhi11

  Theta(2)         = Theta1
  Phi(2)           = Phi2
  VelocityTheta(2) = VelocityTheta21
  VelocityPhi(2)   = VelocityPhi21

  Theta(3)         = Theta2
  Phi(3)           = Phi1
  VelocityTheta(3) = VelocityTheta12
  VelocityPhi(3)   = VelocityPhi12

  Theta(4)         = Theta2
  Phi(4)           = Phi2
  VelocityTheta(4) = VelocityTheta22
  VelocityPhi(4)   = VelocityPhi22


  do i=1,4
     
     Numerator(i) = (Theta3-Theta(i))**2.0 + (Phi3-Phi(i))**2.0

     do j=1,4
        d(i,j)= (Theta(i)-Theta(j))**2.0 + (Phi(i)-Phi(j))**2.0
     enddo
  enddo
     
  VelocityThetaOut = 0.0
  VelocityPhiOut   = 0.0
  
  do i=1,4
     multiplyer = 1.0
     do j=1,4
        if (j .ne. i) then
           multiplyer=multiplyer*numerator(j)/(d(i,j))
        endif
     enddo
     VelocityThetaOut = VelocityThetaOut+VelocityTheta(i)*multiplyer
     VelocityPhiOut   = VelocityPhiOut  +VelocityPhi(i)  *multiplyer
     
  enddo
     
end subroutine interpolate_velocity


subroutine New_Interpolate_Velocity(Theta1,Theta2,Theta3,uTheta1,&
                                    uTheta2,uPhi1,uPhi2,uTheta3,uPhi3)

real, intent(in)            ::   Theta1,Theta2,Theta3,uTheta1,&
                                 uTheta2,uPhi1,uPhi2

real, intent(out)           ::   uTheta3,uPhi3
  
  
uTheta3 = (uTheta2 - uTheta1) / (Theta2-Theta1) * (Theta3-Theta1) + uTheta1

uPhi3   = (uPhi2   - uPhi1)   / (Theta2-Theta1) * (Theta3-Theta1) + uPhi1
  
  
end subroutine New_Interpolate_Velocity


subroutine New_Interpolate_Jr(Theta1,Theta2,Theta3,Jr1,&
                                    Jr2,Jr3)

real, intent(in)            ::   Theta1,Theta2,Theta3,Jr1,&
                                 Jr2

real, intent(out)           ::   Jr3
  
  
Jr3 = (Jr2 - Jr1) / (Theta2-Theta1) * (Theta3-Theta1) + Jr1


  
  
end subroutine New_Interpolate_Jr

!******************************************************************************

subroutine initial_line_location
  
  use ModNumConst, ONLY: cTwoPi
  use ModPWOM
  implicit none

  ! set global variables Dtheta and Dphi
  dTheta = maxval( Theta_G(1:nPhi,1:nTheta)) / (nTheta - 1)
  dPhi   = cTwoPi / (nPhi - 1)

  do iLine=1,nLine
     iThetaLine_I(iLine) = floor(ThetaLine_I(iLine)/ dTheta) + 1
     
     iPhiLine_I(iLine)   = mod(floor(PhiLine_I(iLine)/ dPhi),nPhi-1) + 1
     
     xLine_I(iLine)      = &
          rLowerBoundary*sin(ThetaLine_I(iLine))*cos(PhiLine_I(iLine))
     
     yLine_I(iLine)      = &
          rLowerBoundary*sin(ThetaLine_I(iLine))*sin(PhiLine_I(iLine))
     
     zLine_I(iLine)      = &
          rLowerBoundary*cos(ThetaLine_I(iLine))
     
     xLineOld_I(iLine)   = xLine_I(iLine)
     yLineOld_I(iLine)   = yLine_I(iLine)
     zLineOld_I(iLine)   = zLine_I(iLine)
  enddo
  
end subroutine initial_line_location

!=============================================================================

subroutine move_line
  use ModNumConst, ONLY: cTwoPi, cRadToDeg
  use ModPWOM
  use ModIoUnit, ONLY: UnitTmp_, io_unit_new
  use ModInterpolate, ONLY: bilinear
  implicit none
  real :: a
  character(len=*), parameter :: NameSub = 'PW_move_line'
  logical :: DoTest, DoTestMe
  !---------------------------------------------------------------------------
  call PW_set_do_test(NameSub, DoTest, DoTestMe)

  ! Get the velocity of field line advection from a
  ! bilinear interpolation.
 
  UthetaLine_I(iLine) = bilinear(uExBtheta_C,1,nPhi,1,nTheta, &
       (/ PhiLine_I(iLine)/Dphi+1.0,ThetaLine_I(iLine)/Dtheta+1.0 /) )
  
  UphiLine_I(iLine)   = bilinear(uExBphi_C  ,1,nPhi,1,nTheta, &
       (/ PhiLine_I(iLine)/Dphi+1.0,ThetaLine_I(iLine)/Dtheta+1.0 /) )

  JrLine_I(iLine)     = bilinear(Jr_G, 0,nPhi+1,0,nTheta+1, &
       (/ PhiLine_I(iLine)/Dphi+1.0,ThetaLine_I(iLine)/Dtheta+1.0 /) )
  AvELine_I(iLine)     = bilinear(AvE_G, 0,nPhi+1,0,nTheta+1, &
       (/ PhiLine_I(iLine)/Dphi+1.0,ThetaLine_I(iLine)/Dtheta+1.0 /) )
  EfluxLine_I(iLine)     = bilinear(Eflux_G, 0,nPhi+1,0,nTheta+1, &
       (/ PhiLine_I(iLine)/Dphi+1.0,ThetaLine_I(iLine)/Dtheta+1.0 /) )

  ! save ExB velocity to get joule heating
  if (UseJouleHeating) then
     uJoule2 = (&
       UthetaLine_I(iLine)**2.0&
       +(UphiLine_I(iLine) - OmegaPlanet*rPlanet*sin(ThetaLine_I(iLine)))**2.0)
  endif
  
  xLineOld_I(iLine) = xLine_I(iLine)
  yLineOld_I(iLine) = yLine_I(iLine)
  zLineOld_I(iLine) = zLine_I(iLine)
  
  UxLine_I(iLine)   = &
       UthetaLine_I(iLine)*cos(ThetaLine_I(iLine)) &
       * cos(PhiLine_I(iLine)) &
       - UphiLine_I(iLine)*sin(PhiLine_I(iLine)) 
  
  UyLine_I(iLine)   = &
       UthetaLine_I(iLine)*cos(ThetaLine_I(iLine)) &
       * sin(PhiLine_I(iLine)) & 
       + UphiLine_I(iLine)*cos(PhiLine_I(iLine))

  UzLine_I(iLine)   = &
       -UthetaLine_I(iLine)*sin(ThetaLine_I(iLine))

  if (DoMoveLine) then
     xLine_I(iLine) = UxLine_I(iLine)*DtHorizontal + xLineOld_I(iLine)
     yLine_I(iLine) = UyLine_I(iLine)*DtHorizontal + yLineOld_I(iLine)
     zLine_I(iLine) = UzLine_I(iLine)*DtHorizontal + zLineOld_I(iLine) 
  else
     xLine_I(iLine) =  xLineOld_I(iLine)
     yLine_I(iLine) =  yLineOld_I(iLine)
     zLine_I(iLine) =  zLineOld_I(iLine) 
  endif
  
  if(DoTestMe)then
     write(*,*) NameSub, ': iLine=',iLine
     write(*,*) NameSub, ': Theta, Phi=',ThetaLine_I(iLine)*cRadToDeg, &
          PhiLine_I(iLine)*cRadToDeg
     write(*,*) NameSub, ': xyzLineOld=',xLineOld_I(iLine), &
          yLineOld_I(iLine), zLineOld_I(iLine)
     write(*,*) NameSub, ': UxyzLine  =',UxLine_I(iLine), &
          UyLine_I(iLine), UzLine_I(iLine)
     write(*,*) NameSub, ': DtHorizontal=',DtHorizontal
     write(*,*) NameSub, ': xyzLine   =',xLine_I(iLine), &
          yLine_I(iLine), zLine_I(iLine)
  end if
  
  ! Get new angles from new XYZ positions.Use Theta=arccos(Z/R) and
  ! use Phi = arctan(y/x). make sure to use case statements for phi
  ! depending on what quadrant the position lies in. 

  a = rLowerBoundary / &
       sqrt(xLine_I(iLine)**2 + yLine_I(iLine)**2 + zLine_I(iLine)**2)

  xLine_I(iLine) = xLine_I(iLine)*a
  yLine_I(iLine) = yLine_I(iLine)*a
  zLine_I(iLine) = zLine_I(iLine)*a
  
  ThetaLine_I(iLine) = acos(max(-1.0,min(1.0, zLine_I(iLine)/rLowerBoundary)))
  PhiLine_I(iLine)   = modulo(atan2(yLine_I(iLine), xLine_I(iLine)), cTwoPi)
  
  ! Deal with posibility that phi is negative
  if (PhiLine_I(iLine) .lt. 0.0) then
     write(*,*) 'TTT',PhiLine_I(iLine)
     PhiLine_I(iLine) = PhiLine_I(iLine) + 6.283185
     write(*,*) 'TTTT',PhiLine_I(iLine)
     write(*,*) xLine_I(iLine),yLine_I(iLine),zLine_I(iLine)
     call con_stop('Error: Phi is negative')
  endif
  
  ! Extract the nearest grid point greater then or equal to the 
  ! actual location. The result is used to get the velocity. 
  iThetaLine_I(iLine) = floor( ThetaLine_I(iLine) / Dtheta) + 1
  iPhiLine_I(iLine)   = mod(floor(PhiLine_I(iLine) / Dphi),nPhi-1)+1

  if(DoTestMe)then
     write(*,*) NameSub, ': Factor a=',a
     write(*,*) NameSub, ': xyzLine   =',xLine_I(iLine), &
          yLine_I(iLine), zLine_I(iLine)
     write(*,*) NameSub, ': Theta, Phi=',ThetaLine_I(iLine)*cRadToDeg, &
          PhiLine_I(iLine)*cRadToDeg
  end if

end subroutine move_line

!==============================================================================

subroutine PW_advance_line

  !  This routine advances a single field line

  use ModNumConst, ONLY: cRadToDeg
  use ModPWOM
  use ModFieldLine
  use ModAurora,   ONLY: set_aurora, set_auroral_rates
  implicit none

  real XXX,MaxLineTime
  logical :: DoLog
  !---------------------------------------------------------------------------
  ! Get the GeoMagnetic latitude and longitude 
  GeoMagLat_I(iLine) = 90.0 - ThetaLine_I(iLine)*cRadToDeg
  GeoMagLon_I(iLine) = PhiLine_I(iLine)*cRadToDeg
  
  !Calculate the horizontal fieldline velocity in cm/s to be used for 
  !centrifugal force in model
  if (UseCentrifugal) then
     OmegaLine_I(iLine)= &
          sqrt(UthetaLine_I(iLine)**2.0+UphiLine_I(iLine)**2.0)/rLowerBoundary
  else
     OmegaLine_I(iLine) = 0.0
  endif
  
  MaxLineTime=Time+DtHorizontal
  if (MaxLineTime > Tmax) MaxLineTime = Tmax
  
  if (nLog ==-1 .or. iLineGlobal(iLine) == nLog) then
     DoLog=.true.
  else
     DoLog=.False.
  endif

!  if (UseAurora) then
!     call set_aurora
!     call set_auroral_rates
!  endif

  call put_field_line(nAlt,&
       State_CVI(:,:,iLine),&
       GeoMagLat_I(iLine),GeoMagLon_I(iLine),JrLine_I(iLine),              &
       OmegaLine_I(iLine), uJoule2=uJoule2,iUnitOutput=iUnitOutput(iLine), &
       NameRestart=NameRestart(iLine),  &
       iLine=iLine,Time=Time,MaxLineTime=MaxLineTime,TypeSolver=TypeSolver,&
       DToutput=DToutput,    &
       DoLog=DoLog,nStep=nStep,Dt=Dt_I(iLine),&
       AvE=AvELine_I(iLine),Eflux=EfluxLine_I(iLine))
    
  call polar_wind
  
  if (iLine == nLine) then
     call get_field_line( nAlt,&
       State_CVI(:,:,iLine),&
       GeoMagLat_I(iLine),GeoMagLon_I(iLine),JrLine_I(iLine),               &
       OmegaLine_I(iLine),       &
       iLine=iLine,Time=Time,MaxLineTime=MaxLineTime,nStep=nStep,r_C=r_C, &
       Dt=Dt_I(iLine))
     DoSavePlot=.true.
  else
     call get_field_line( nAlt,&
       State_CVI(:,:,iLine),&
       GeoMagLat_I(iLine),GeoMagLon_I(iLine),JrLine_I(iLine),               &
       OmegaLine_I(iLine),       &
       iLine=iLine,MaxLineTime=MaxLineTime,Dt=Dt_I(iLine))
  endif
  
end subroutine PW_advance_line

!============================================================================
subroutine PW_set_do_test(String, DoTest, DoTestMe)
  use ModPWOM, ONLY : iProc, iProcTest, iLine, iLineTest, StringTest
  implicit none
  
  character (len=*), intent(in) :: String
  logical, intent(out) :: DoTest, DoTestMe

  if(iLineTest == 0 .or. iLine == iLineTest)then
     DoTest = index(StringTest, String) > 0
     DoTestMe = DoTest .and. iProc == iProcTest
  else
     DoTest = .false.; DoTestMe = .false.
  end if
end subroutine PW_set_do_test
