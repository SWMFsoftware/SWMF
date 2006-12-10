subroutine PW_print_electrodynamics

  use ModIoUnit, ONLY: UnitTmp_
  use ModPWOM
  use ModNumConst, ONLY:cDegToRad
  implicit none

  integer :: TimeOut
  Character(len=100) :: NameElectrodynamics
!-----------------------------------------------------------------------------


!******************************************************************************
!  Write output, use cartesian coords for output
!******************************************************************************
 
  do iPhi=1,nPhi
     do iTheta=1,nTheta
        ux(iPhi,iTheta) =  & 
             VelocityExBtheta(iPhi,iTheta)*cos(Theta_G(iPhi,iTheta)) &
             * cos(Phi_G(iPhi,iTheta)) &
             - VelocityExBphi(iPhi,iTheta)*sin(Phi_G(iPhi,iTheta))   &
             + VelocityExBr(iPhi,iTheta)  *sin(Theta_G(iPhi,iTheta)) &
             * cos(Phi_G(iPhi,iTheta))
        
        uy(iPhi,iTheta) =  &
             VelocityExBtheta(iPhi,iTheta)*cos(Theta_G(iPhi,iTheta)) &
             * sin(Phi_G(iPhi,iTheta)) &
             + VelocityExBphi(iPhi,iTheta)*cos(Phi_G(iPhi,iTheta))   &
             + VelocityExBr(iPhi,iTheta)  *sin(Theta_G(iPhi,iTheta)) &
             * sin(Phi_G(iPhi,iTheta))
        
        uz(iPhi,iTheta) =  &
             -VelocityExBtheta(iPhi,iTheta)*sin(Theta_G(iPhi,iTheta)) &
             + VelocityExBr(iPhi,iTheta)  *cos(Theta_G(iPhi,iTheta)) 
        
        x(iPhi,iTheta)  =  &
             rLowerBoundary*sin(Theta_G(iPhi,iTheta))*cos(Phi_G(iPhi,iTheta))
        
        y(iPhi,iTheta)  =  &
             rLowerBoundary*sin(Theta_G(iPhi,iTheta))*sin(Phi_G(iPhi,iTheta))
        
        z(iPhi,iTheta)  =  &
             rLowerBoundary*cos(Theta_G(iPhi,iTheta))
        
        Ex(iPhi,iTheta) =  & 
             Etheta(iPhi,iTheta)*cos(Theta_G(iPhi,iTheta)) &
             * cos(Phi_G(iPhi,iTheta))                     &
             - Ephi(iPhi,iTheta)*sin(Phi_G(iPhi,iTheta))   &
             + Er(iPhi,iTheta)  *sin(Theta_G(iPhi,iTheta)) &
             * cos(Phi_G(iPhi,iTheta))
        
        Ey(iPhi,iTheta) =  &
             Etheta(iPhi,iTheta)*cos(Theta_G(iPhi,iTheta)) &
             * sin(Phi_G(iPhi,iTheta))                     &
             + Ephi(iPhi,iTheta)*cos(Phi_G(iPhi,iTheta))   &
             + Er(iPhi,iTheta)  *sin(Theta_G(iPhi,iTheta)) &
             * sin(Phi_G(iPhi,iTheta))
        
        Ez(iPhi,iTheta) =  &
             -Etheta(iPhi,iTheta)*sin(Theta_G(iPhi,iTheta))&
             + Er(iPhi,iTheta)  *cos(Theta_G(iPhi,iTheta)) 
     enddo
  enddo
 
  TimeOut=int(Time)
  write(NameElectrodynamics,"(a,i8.8,a)") &
       'PW/Electrodynamics_Time',TimeOut,'.dat'
  open(UnitTmp_,FILE=NameElectrodynamics)
  write(UnitTmp_,*) &
     'VARIABLES = "X", "Y", "Z", "Ux", "Uy", "Uz", "V", "Ex", "Ey", "Ez", "Jr"'
  
  write(UnitTmp_,*) 'Zone I=', nPhi, ', J=', nTheta,', DATAPACKING=POINT'
  
  do iTheta=1,nTheta
     do iPhi=1,nPhi
        write(UnitTmp_,*) &
             x(iPhi,iTheta),y(iPhi,iTheta),z(iPhi,iTheta),     &
             ux(iPhi,iTheta),uy(iPhi,iTheta),uz(iPhi,iTheta),  &
             Potential_G(iPhi,iTheta),                         &
             Ex(iPhi,iTheta), Ey(iPhi,iTheta), Ez(iPhi,iTheta),&
             Jr_G(iPhi,iTheta)
        
     enddo
  enddo
  
  close(UnitTmp_)


end subroutine PW_print_electrodynamics
