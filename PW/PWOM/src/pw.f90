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
  
  use ModPWOM
  implicit none

  ! set global variables Dtheta and Dphi
  DTheta = maxval( Theta_G(:,:)) / nTheta
  DPhi   = maxval( Phi_G(:,:)  ) / nPhi

  do iLine=1,nLine
     if (.not. IsRestart) then
        FieldLineTheta (iLine) = 10.0 * 3.141592653589/180.0
        FieldLinePhi   (iLine) = 0.0  * 3.141592653589/180.0
     else
        FieldLineTheta (iLine) = 1.57079632679-GeoMagLat(iLine)&
             * 3.141592653589/180.0
        FieldLinePhi   (iLine) = GeoMagLon(iLine)              &
             * 3.141592653589/180.0
     endif
     iFieldLineTheta(iLine) = &
          ceiling(  FieldLineTheta(iLine)  / Dtheta)
     
     iFieldLinePhi(iLine)   = &
          mod(ceiling(FieldLinePhi(iLine)  / Dphi),nPhi-1)
     
     if (FieldLinePhi(iLine) .eq. 0.0)  iFieldLinePhi(iLine) = 1
     
     
     FieldLineX(iLine)      = &
          rLowerBoundary*sin(FieldLineTheta(iLine))*cos(FieldLinePhi(iLine))
     
     FieldLineY(iLine)      = &
          rLowerBoundary*sin(FieldLineTheta(iLine))*sin(FieldLinePhi(iLine))
     
     FieldLineZ(iLine)      = &
          rLowerBoundary*cos(FieldLineTheta(iLine))
     
     OldFieldLineX(iLine)   = FieldLineX(iLine)
     OldFieldLineY(iLine)   = FieldLineY(iLine)
     OldFieldLineZ(iLine)   = FieldLineZ(iLine)
  enddo
  
end subroutine initial_line_location


!******************************************************************************

subroutine MoveFluxTube
  use ModPWOM
  use ModIoUnit, ONLY: UnitTmp_, io_unit_new
  implicit none
  

  ! if the field line is in the polar circle then set velocity equal
  ! to that of the nearest point. Otherwise, get the velocity from an
  !interpolation.
  
  if (iFieldLineTheta(iLine) .le. 1) then
     write(*,*) 'TT',FieldLineTheta(iLine),Dtheta,FieldLineZ(iLine),OldFieldLineZ(iLine)
     
     
     
     !stop
     FieldLineVelTheta(iLine)= &
          VelocityExBTheta(iFieldLinePhi(iLine),iFieldLineTheta(iLine))
     FieldLineVelPhi(iLine)= &
          VelocityExBPhi(iFieldLinePhi(iLine),iFieldLineTheta(iLine))
     
     FieldLineJr(iLine)    = &
          Jr_G(iFieldLinePhi(iLine),iFieldLineTheta(iLine))
     write(*,*) iFieldLineTheta(iline),time,FieldLineVelX(iLine),FieldLineX(iLine)
  else
     Call New_Interpolate_Velocity( &
          Theta_G(iFieldLinePhi(iLine),iFieldLineTheta(iLine)),           &
          Theta_G(iFieldLinePhi(iLine),iFieldLineTheta(iLine)-1),         &
          FieldLineTheta(iLine),                                          &
          VelocityExBTheta(iFieldLinePhi(iLine),iFieldLineTheta(iLine)),  &
          VelocityExBTheta(iFieldLinePhi(iLine),iFieldLineTheta(iLine)-1),&
          VelocityExBPhi(iFieldLinePhi(iLine),iFieldLineTheta(iLine)),    &
          VelocityExBPhi(iFieldLinePhi(iLine),iFieldLineTheta(iLine)-1),  &
          FieldLineVelTheta(iLine),FieldLineVelPhi(iLine))
     
     Call New_Interpolate_Jr( &
          Theta_G(iFieldLinePhi(iLine),iFieldLineTheta(iLine)),           &
          Theta_G(iFieldLinePhi(iLine),iFieldLineTheta(iLine)-1),         &
          FieldLineTheta(iLine),                                          &
          Jr_G(iFieldLinePhi(iLine),iFieldLineTheta(iLine)),  &
          Jr_G(iFieldLinePhi(iLine),iFieldLineTheta(iLine)-1),&
          FieldLineJr(iLine))
  endif
  
  
  
  OldFieldLineTheta(iLine) = FieldLineTheta(iLine)
  OldFieldLinePhi(iLine)   = FieldLinePhi(iLine)
  
  
  OldFieldLineX(iLine) = FieldLineX(iLine)
  OldFieldLineY(iLine) = FieldLineY(iLine)
  OldFieldLineZ(iLine) = FieldLineZ(iLine)
  
  
  FieldLineVelx(iLine)   = &
       FieldLineVelTheta(iLine)*cos(OldFieldLineTheta(iLine)) &
       * cos(OldFieldLinePhi(iLine)) &
       - FieldLineVelPhi(iLine)*sin(OldFieldLinePhi(iLine)) 
  
  FieldLineVely(iLine)   = &
       FieldLineVelTheta(iLine)*cos(OldFieldLineTheta(iLine)) &
       * sin(OldFieldLinePhi(iLine)) & 
       + FieldLineVelPhi(iLine)*cos(OldFieldLinePhi(iLine))
  
  FieldLineVelz(iLine)   = &
       -FieldLineVelTheta(iLine)*sin(OldFieldLineTheta(iLine))
  
  FieldLineX(iLine) =  &
       (FieldLineVelx(iLine)) *Dt     &
       + OldFieldLineX(iLine)                                   
  
  FieldLineY(iLine) =  &
       (FieldLineVely(iLine)) *Dt     &
       + OldFieldLineY(iLine)                                   
  
  FieldLineZ(iLine) =  &
       (FieldLineVelz(iLine)) *Dt     &
       + OldFieldLineZ(iLine) 
  
  if (iFieldLineTheta(iLine) .le. 1) write(*,*) FieldLineZ(iLine)/rLowerBoundary,Time
  
  ! Get new angles from new XYZ positions.Use Theta=arccos(Z/R) and
  ! use Phi = arctan(y/x). make sure to use case statements for phi
  ! depending on what quadrant the position lies in. 
  
  if (FieldLineZ(iLine)/rLowerBoundary .ge.1) then
     FieldLineTheta(iLine) = 1.57079632679 &
          - acos(sqrt(FieldLineX(iLine)**2+FieldLineY(iLine)**2) & 
          / rLowerBoundary)
  else
     FieldLineTheta(iLine) = &
          abs(acos(FieldLineZ(iLine)/rLowerBoundary))
  endif
  
  
  
  
  if (FieldLineX(iLine) .gt. 0.0 .and. FieldLineY(iLine) .ge. 0.0) then
     
     FieldLinePhi(iLine)   = &
          (atan(FieldLineY(iLine)/FieldLineX(iLine)))
  else if(FieldLineX(iLine) .lt. 0.0 .and. FieldLineY(iLine).ge.0.0) then
     FieldLinePhi(iLine)   = &
          3.14159-(atan(abs(FieldLineY(iLine)/FieldLineX(iLine))))
  else if(FieldLineX(iLine) .lt.0.0.and. FieldLineY(iLine) .le. 0.0) then
     FieldLinePhi(iLine)   = &
          3.14159+(atan(abs(FieldLineY(iLine)/FieldLineX(iLine))))
  else if(FieldLineX(iLine) .gt.0.0.and. FieldLineY(iLine) .le. 0.0) then
     FieldLinePhi(iLine)   = &
          6.283185-(atan(abs(FieldLineY(iLine)/FieldLineX(iLine))))
  else if(FieldLineX(iLine) .eq.0.0.and. FieldLineY(iLine) .ge.0.0) then
     FieldLinePhi(iLine)   = &
          3.14159/2.0
  else if(FieldLineX(iLine) .eq.0.0.and. FieldLineY(iLine) .le.0.0) then
     FieldLinePhi(iLine)   = &
          3.0*3.14159/2.0
  endif
  
  
  ! Deal with posibility that phi is negative
  if (FieldLinePhi(iLine) .lt. 0.0) then
     write(*,*) 'TTT',FieldLinePhi(iLine)
     FieldLinePhi(iLine) = FieldLinePhi(iLine) + 6.283185
     write(*,*) 'TTTT',FieldLinePhi(iLine)
     write(*,*) FieldLineX(iLine),FieldLineY(iLine),FieldLineZ(iLine)
     call con_stop('Error: Phi is negative')
  endif
  
  
  
  ! Extract the nearest grid point greater then or equal to the 
  ! actual location. The result is used to get the velocity. 
  iFieldLineTheta(iLine) = &
       ceiling(  FieldLineTheta(iLine)  / Dtheta)
  
  iFieldLinePhi(iLine)   = &
       mod(ceiling(FieldLinePhi(iLine)  / Dphi),nPhi-1)
  if (iFieldLinePhi(iLine) .eq. 0)  iFieldLinePhi(iLine) = 1
  if (iFieldLineTheta(iLine) .eq. 0)  iFieldLineTheta(iLine) = 1
  
  ! Put field line locations back on the sphere
  
  if (FieldLineTheta(iLine) .ne.0.0 ) then
     FieldLineX(iLine)      = &
          rLowerBoundary*sin(FieldLineTheta(iLine))     &
          * cos(FieldLinePhi(iLine))
     
     FieldLineY(iLine)      = &
          rLowerBoundary*sin(FieldLineTheta(iLine))     &
          * sin(FieldLinePhi(iLine))
     
     FieldLineZ(iLine)      = &
          rLowerBoundary*cos(FieldLineTheta(iLine))
     
  endif
  

!******************************************************************************
!  Write output, use cartesian coords for output
!******************************************************************************
! 
!  do iPhi=1,nPhi
!     do iTheta=1,nTheta
!        ux(iPhi,iTheta) =  & 
!             VelocityExBtheta(iPhi,iTheta)*cos(Theta_G(iPhi,iTheta)) &
!             * cos(Phi_G(iPhi,iTheta)) &
!             - VelocityExBphi(iPhi,iTheta)*sin(Phi_G(iPhi,iTheta))   &
!             + VelocityExBr(iPhi,iTheta)  *sin(Theta_G(iPhi,iTheta)) &
!             * cos(Phi_G(iPhi,iTheta))
!        
!        uy(iPhi,iTheta) =  &
!             VelocityExBtheta(iPhi,iTheta)*cos(Theta_G(iPhi,iTheta)) &
!             * sin(Phi_G(iPhi,iTheta)) &
!             + VelocityExBphi(iPhi,iTheta)*cos(Phi_G(iPhi,iTheta))   &
!             + VelocityExBr(iPhi,iTheta)  *sin(Theta_G(iPhi,iTheta)) &
!             * sin(Phi_G(iPhi,iTheta))
!        
!        uz(iPhi,iTheta) =  &
!             -VelocityExBtheta(iPhi,iTheta)*sin(Theta_G(iPhi,iTheta)) &
!             + VelocityExBr(iPhi,iTheta)  *cos(Theta_G(iPhi,iTheta)) 
!        
!        x(iPhi,iTheta)  =  &
!             rLowerBoundary*sin(Theta_G(iPhi,iTheta))*cos(Phi_G(iPhi,iTheta))
!        
!        y(iPhi,iTheta)  =  &
!             rLowerBoundary*sin(Theta_G(iPhi,iTheta))*sin(Phi_G(iPhi,iTheta))
!        
!        z(iPhi,iTheta)  =  &
!             rLowerBoundary*cos(Theta_G(iPhi,iTheta))
!        
!        Ex(iPhi,iTheta) =  & 
!             Etheta(iPhi,iTheta)*cos(Theta_G(iPhi,iTheta)) &
!             * cos(Phi_G(iPhi,iTheta))                     &
!             - Ephi(iPhi,iTheta)*sin(Phi_G(iPhi,iTheta))   &
!             + Er(iPhi,iTheta)  *sin(Theta_G(iPhi,iTheta)) &
!             * cos(Phi_G(iPhi,iTheta))
!        
!        Ey(iPhi,iTheta) =  &
!             Etheta(iPhi,iTheta)*cos(Theta_G(iPhi,iTheta)) &
!             * sin(Phi_G(iPhi,iTheta))                     &
!             + Ephi(iPhi,iTheta)*cos(Phi_G(iPhi,iTheta))   &
!             + Er(iPhi,iTheta)  *sin(Theta_G(iPhi,iTheta)) &
!             * sin(Phi_G(iPhi,iTheta))
!        
!        Ez(iPhi,iTheta) =  &
!             -Etheta(iPhi,iTheta)*sin(Theta_G(iPhi,iTheta))&
!             + Er(iPhi,iTheta)  *cos(Theta_G(iPhi,iTheta)) 
!     enddo
!  enddo
! 
!
! open(UnitTmp_,FILE='PW/Test.out')
!  write(UnitTmp_,*) &
!     'VARIABLES = "X", "Y", "Z", "Ux", "Uy", "Uz", "V", "Ex", "Ey", "Ez", "Jr"'
!  
!  write(UnitTmp_,*) 'Zone I=', nPhi, ', J=', nTheta,', DATAPACKING=POINT'
!  
!  do iTheta=1,nTheta
!     do iPhi=1,nPhi
!        write(UnitTmp_,*) &
!             x(iPhi,iTheta),y(iPhi,iTheta),z(iPhi,iTheta),     &
!             ux(iPhi,iTheta),uy(iPhi,iTheta),uz(iPhi,iTheta),  &
!             Potential_G(iPhi,iTheta),                         &
!             Ex(iPhi,iTheta), Ey(iPhi,iTheta), Ez(iPhi,iTheta),&
!             Jr_G(iPhi,iTheta)
!        
!     enddo
!  enddo
!  
!  close(UnitTmp_)
end subroutine MoveFluxTube


!******************************************************************************


!******************************************************************************
!  This routine advances a single field line
!******************************************************************************

subroutine PW_advance_line
  use ModPWOM
  use ModFieldLine
  implicit none

  real XXX,MaxLineTime

  ! Get the GeoMagnetic latitude and longitude 
  GeoMagLat(iLine) = (1.57079632679 - OldFieldLineTheta(iLine)) &
       * 180.0 / 3.141592653589
  
  GeoMagLon(iLine) = (OldFieldLinePhi(iLine)) &
       * 180.0 / 3.141592653589
  
  !Calculate the horizontal fieldline velocity in cm/s to be used for 
  !centrifugal force in model
  if (IsCentrifugal) then
     OmegaHorFieldLine(iLine)= sqrt(( &
          FieldLineVelTheta(iLine)**2.0+FieldLineVelPhi(iLine)**2.0) &
          /rLowerBoundary**2.0)
  else
     OmegaHorFieldLine(iLine) = 0.0
  endif
  
  MaxLineTime=Time+DT
  if (MaxLineTime > Tmax) MaxLineTime = Tmax
  

  call put_field_line(&
       dOxyg(:,iLine), uOxyg(:,iLine), pOxyg(:,iLine), TOxyg(:,iLine),     &
       dHel(:,iLine), uHel(:,iLine), pHel(:,iLine), THel(:,iLine),         &
       dHyd(:,iLine), uHyd(:,iLine), pHyd(:,iLine), THyd(:,iLine),         &
       dElect(:,iLine), uElect(:,iLine), pElect(:,iLine), TElect(:,iLine), &
       GeoMagLat(iLine),GeoMagLon(iLine),FieldLineJr(iLine),               &
       OmegaHorFieldLine(iLine), iUnitOutput=iUnitOutput,                  &
       iUnitGraphics=iUnitGraphics(iLine),NameRestart=NameRestart(iLine),  &
       iLine=iLine,Time=Time,MaxLineTime=MaxLineTime,TypeSolver=TypeSolver,&
       IsVariableDt=IsVariableDt,IsRestart=IsRestart,DToutput=DToutput,    &
       nDim=nDim)
    
  call polar_wind
  
  if (iLine == nLine) then
     call get_field_line( &
       dOxyg(:,iLine), uOxyg(:,iLine), pOxyg(:,iLine), TOxyg(:,iLine),     &
       dHel(:,iLine), uHel(:,iLine), pHel(:,iLine), THel(:,iLine),         &
       dHyd(:,iLine), uHyd(:,iLine), pHyd(:,iLine), THyd(:,iLine),         &
       dElect(:,iLine), uElect(:,iLine), pElect(:,iLine), TElect(:,iLine), &
       GeoMagLat(iLine),GeoMagLon(iLine),FieldLineJr(iLine),               &
       OmegaHorFieldLine(iLine), iUnitGraphics=iUnitGraphics(iLine),       &
       iUnitOutput=iUnitOutput,iLine=iLine,Time=Time,MaxLineTime=MaxLineTime)
  else
     call get_field_line( &
       dOxyg(:,iLine), uOxyg(:,iLine), pOxyg(:,iLine), TOxyg(:,iLine),     &
       dHel(:,iLine), uHel(:,iLine), pHel(:,iLine), THel(:,iLine),         &
       dHyd(:,iLine), uHyd(:,iLine), pHyd(:,iLine), THyd(:,iLine),         &
       dElect(:,iLine), uElect(:,iLine), pElect(:,iLine), TElect(:,iLine), &
       GeoMagLat(iLine),GeoMagLon(iLine),FieldLineJr(iLine),               &
       OmegaHorFieldLine(iLine), iUnitGraphics=iUnitGraphics(iLine),       &
       iUnitOutput=iUnitOutput,iLine=iLine,MaxLineTime=MaxLineTime)
  endif
  
end subroutine PW_advance_line

