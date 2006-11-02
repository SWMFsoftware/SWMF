program pw

  use ModPwom
  use ModFieldLine
  use ModMpi
  use ModReadParam
  implicit none

  !****************************************************************************
  ! Initiallize MPI and get number of processors and rank of given processor
  !****************************************************************************

  !---------------------------------------------------------------------------
  call MPI_INIT(errcode)
  iComm = MPI_COMM_WORLD

  call MPI_COMM_RANK(iComm,iProc,errcode)
  call MPI_COMM_SIZE(iComm,nProc,errcode)

  !****************************************************************************
  ! Read the input file
  !****************************************************************************
  IsStandAlone      = .true.
  NameThisComp = 'PW'
  NameInput          = 'pw.input'

  call read_file('pw.input',iComm)
  call read_init('  ',iSessionIn=1,iLineIn=0)

  call PW_set_parameters('READ')
  ! call PW_set_parameters('CHECK')

  call PW_initialize

  !****************************************************************************
  ! Use Get_GITM to bring in neutral atmosphere from GITM
  !****************************************************************************
  !call GetNeutralData

  !****************************************************************************
  !  Set parameters for reading in potential and time of simulation
  !****************************************************************************


  Dtheta  = 0.0242
  Dphi    = 0.0245

  Dt      =    50.0
  !  Dt      =    0.2
  !  Dt      =    Tmax
  !maxTime = 10000.0
  maxTime = Tmax
  Time    =     0.0


  !****************************************************************************
  ! Read information from IE file, and get the velocities
  !****************************************************************************

  call Get_ElectrodynamicPW  



  !****************************************************************************
  !  Move flux tube around
  !****************************************************************************

  !initialize field line locations

  call Initial_Line_Location

  !****************************************************************************
  ! Move the flux tube, solve each fieldline, and advance the time
  !****************************************************************************

  TIMELOOP:do
     if (Time >= Tmax) exit TIMELOOP
     do iLine=1,nLine

        ! MoveFluxTube moves the flux tube, then we can use the angular
        !position to get the lat and lon

        call MoveFluxTube

        !  Call the flux tube to be solved

        call AdvancePWline
     enddo
  enddo TIMELOOP


  !****************************************************************************
  !  Write output, use cartesian coords for output
  !****************************************************************************


  do iLine=1,nLine
     CLOSE(UNIT=iUnitGraphics(iLine))
  enddo
  close(UNIT=iUnitOutput)

  call MPI_FINALIZE(errcode)

end program pw



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
!  This subroutine gets the electrodynamic parameters, and gets the convection
!  velocity. 
!******************************************************************************
subroutine Get_ElectrodynamicPW

  use ModIoUnit, ONLY: UnitTmp_
  use ModPWOM
  implicit none
  !---------------------------------------------------------------------------
  open(UnitTmp_, FILE=NamePhiNorth)  

  do iPhi=1,nPhi
     do iTheta=1,nTheta
        read(unit=UnitTmp_,fmt='(6(1PE13.5))') &
             Theta_G(iPhi,iTheta),Phi_G(iPhi,iTheta),SigmaH_G(iPhi,iTheta),&
             SigmaP_G(iPhi,iTheta),Jr_G(iPhi,iTheta),Potential_G(iPhi,iTheta)

     enddo
  enddo

  close(UnitTmp_)


!******************************************************************************
!  Change angles to radians
!******************************************************************************
  Theta_G(:,:) = Theta_G(:,:)*3.1416/180.0
  Phi_G  (:,:) = Phi_G  (:,:)*3.1416/180.0

!******************************************************************************
!  Convert potential from kilovolts to Volts
!******************************************************************************
  Potential_G(:,:) = Potential_G(:,:)*1.0e3
!******************************************************************************
!  Calc Bfield components
!******************************************************************************
  rPlanet        = 6378000.0
  rLowerBoundary = rPlanet+110.0e3
  MagMoment      = 7.84e15
  Bcoef          = MagMoment/(rLowerBoundary)**3 
  do iPhi=1,nPhi
     do iTheta=1,nTheta
        Br_G(iPhi,iTheta)     = -Bcoef*2.0*sin(Theta_G(iPhi,iTheta))
        Btheta_G(iPhi,iTheta) = Bcoef*cos(Theta_G(iPhi,iTheta))
     enddo
  enddo
  
  BmagnitudeSquared_G(:,:) = Br_G(:,:)**2 + Btheta_G(:,:)**2
!******************************************************************************
!  Fill Ghost cells
!******************************************************************************
  Theta_G    (nPhi+1,:) = Theta_G    (2,:) 
  Phi_G      (nPhi+1,:) = Phi_G      (2,:)
  SigmaH_G   (nPhi+1,:) = SigmaH_G   (2,:)
  SigmaP_G   (nPhi+1,:) = SigmaP_G   (2,:) 
  Jr_G       (nPhi+1,:) = Jr_G       (2,:)
  Potential_G(nPhi+1,:) = Potential_G(2,:)

  Theta_G    (0,:) = Theta_G    (nPhi-1,:) 
  Phi_G      (0,:) = Phi_G      (nPhi-1,:)
  SigmaH_G   (0,:) = SigmaH_G   (nPhi-1,:)
  SigmaP_G   (0,:) = SigmaP_G   (nPhi-1,:) 
  Jr_G       (0,:) = Jr_G       (nPhi-1,:)
  Potential_G(0,:) = Potential_G(nPhi-1,:)

  Theta_G    (:,nTheta+1) = Theta_G    (:,nTheta) 
  Phi_G      (:,nTheta+1) = Phi_G      (:,nTheta)
  SigmaH_G   (:,nTheta+1) = SigmaH_G   (:,nTheta)
  SigmaP_G   (:,nTheta+1) = SigmaP_G   (:,nTheta) 
  Jr_G       (:,nTheta+1) = Jr_G       (:,nTheta)
  Potential_G(:,nTheta+1) = Potential_G(:,nTheta)

  do iPhi=1,nPhi
     Theta_G    (iPhi,0) = Theta_G    (mod(iPhi+128,nPhi-1),2) 
     Phi_G      (iPhi,0) = Phi_G      (mod(iPhi+128,nPhi-1),2)
     SigmaH_G   (iPhi,0) = SigmaH_G   (mod(iPhi+128,nPhi-1),2)
     SigmaP_G   (iPhi,0) = SigmaP_G   (mod(iPhi+128,nPhi-1),2) 
     Jr_G       (iPhi,0) = Jr_G       (mod(iPhi+128,nPhi-1),2)
     Potential_G(iPhi,0) = Potential_G(mod(iPhi+128,nPhi-1),2)
  enddo



!******************************************************************************
!  Calc electric field from E=-grad Potential
!******************************************************************************
 
  
  do iPhi=1,nPhi
     do iTheta=1,nTheta
        Etheta(iPhi,iTheta) = &
             (Potential_G(iPhi,iTheta)-Potential_G(iPhi,iTheta+1))&
             /((rLowerBoundary)*Dtheta)
        if (iTheta == 1) then
           Ephi(iPhi,iTheta) = 0.0
        else
           Ephi(iPhi,iTheta)   = &
                (Potential_G(iPhi-1,iTheta)-Potential_G(iPhi+1,iTheta))&
                / ((rLowerBoundary)*2.0*Dphi  &
                * sin(Theta_G(iPhi,iTheta)))
        endif
    enddo
  enddo
  Er(:,:) = 0.0
!******************************************************************************
!  Calc VelocityExB drift from E and B
!******************************************************************************
  

  do iPhi=1,nPhi
     do iTheta=1,nTheta
        VelocityExBtheta(iPhi,iTheta) = &
             (Ephi(iPhi,iTheta)*Br_G(iPhi,iTheta)) &
             / BmagnitudeSquared_G(iPhi,iTheta)
        
        VelocityExBphi(iPhi,iTheta)   = &
             (-Etheta(iPhi,iTheta)*Br_G(iPhi,iTheta)) &
             / BmagnitudeSquared_G(iPhi,iTheta)
        
        VelocityExBr(iPhi,iTheta)     = &
             (-Ephi(iPhi,iTheta)*Btheta_G(iPhi,iTheta)) &
             / BmagnitudeSquared_G(iPhi,iTheta)


     enddo
  enddo
  VelocityExBtheta(:,1) = VelocityExBtheta(:,2)
  VelocityExBphi  (:,1) = VelocityExBphi  (:,2)
  
  if (.not.IsMoveFluxTube) Then 
     VelocityExBtheta(:,:) = 0.0
     VelocityExBPhi  (:,:) = 0.0
     VelocityExBr    (:,:) = 0.0
  endif

  if (.not.IsUseJr) Then 
     Jr_G(:,:) = 0.0
  endif


end subroutine Get_ElectrodynamicPW



!******************************************************************************

subroutine Initial_Line_Location
  
  use ModPWOM
  implicit none
  
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
  
end subroutine Initial_Line_Location


!******************************************************************************

subroutine MoveFluxTube
  use ModPWOM
  implicit none
  

  ! if the field line is in the polar circle then set velocity equal
  ! to that of the nearest point. Otherwise, get the velocity from an
  !interpolation.
  
  if (iFieldLineTheta(iLine) .le. 1) then
     write(*,*) 'TT',FieldLineTheta(iLine)/Dtheta,FieldLineZ(iLine),OldFieldLineZ(iLine)
     
     
     
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
  
end subroutine MoveFluxTube


!******************************************************************************


!******************************************************************************
!  This routine advances a single field line
!******************************************************************************

subroutine AdvancePWline
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
  
  call put_field_line(&
       dOxyg(:,iLine), uOxyg(:,iLine), pOxyg(:,iLine), TOxyg(:,iLine),     &
       dHel(:,iLine), uHel(:,iLine), pHel(:,iLine), THel(:,iLine),         &
       dHyd(:,iLine), uHyd(:,iLine), pHyd(:,iLine), THyd(:,iLine),         &
       dElect(:,iLine), uElect(:,iLine), pElect(:,iLine), TElect(:,iLine), &
       GeoMagLat(iLine),GeoMagLon(iLine),FieldLineJr(iLine),               &
       OmegaHorFieldLine(iLine), iUnitOutput=iUnitOutput,                  &
       iUnitGraphics=iUnitGraphics(iLine),NameRestart=NameRestart(iLine),  &
       iLine=iLine,Time=Time,MaxLineTime=MaxLineTime,TypeSolver=TypeSolver,&
       IsVariableDt=IsVariableDt,IsRestart=IsRestart,DToutput=DToutput)
    
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
  
end subroutine AdvancePWline

!============================================================================
! The following subroutines are here so that we can use SWMF library routines
! Also some features available in SWMF mode only require empty subroutines
! for compilation of the stand alone code.
!============================================================================
subroutine CON_stop(StringError)
  implicit none
  character (len=*), intent(in) :: StringError
  call stop_mpi(StringError)
end subroutine CON_stop

subroutine stop_mpi(str)

  use ModPWOM, ONLY : NameThisComp,IsStandAlone,iProc,iComm
  use ModMpi
  implicit none

  character (len=*), intent(in) :: str

  ! Local variables:
  integer :: iError,nError

  !----------------------------------------------------------------------------

  if(IsStandAlone)then
     write(*,*)'Stopping execution! me=',iProc,' at iteration=',&
          ' with msg:'
     write(*,*)str
     call MPI_abort(iComm, nError, iError)
     stop
  else
     write(*,*)NameThisComp,': stopping execution! me=',iProc
     call CON_stop(NameThisComp//':'//str)
  end if

end subroutine stop_mpi

