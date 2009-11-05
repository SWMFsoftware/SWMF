subroutine PW_print_electrodynamics

  use ModIoUnit, ONLY: UnitTmp_
  use ModPWOM
  use ModNumConst, ONLY:cDegToRad,cRadToDeg,cPi
  implicit none

  real,dimension(:,:),allocatable  :: Ux,Uy,Uz,x,y,z,Ex,Ey,Ez
  real :: Lat,Lon
  integer :: TimeOut
  Character(len=100) :: NameElectrodynamics
!-----------------------------------------------------------------------------


  !Alocate output arrays
  allocate(&
       Ux(nPhi, nTheta), &
       Uy(nPhi, nTheta), &
       Uz(nPhi, nTheta), &
       x (nPhi, nTheta), &
       y (nPhi, nTheta), &
       z (nPhi, nTheta), &
       Ex(nPhi, nTheta), &
       Ey(nPhi, nTheta), &
       Ez(nPhi, nTheta))


!******************************************************************************
!  Write output, use cartesian coords for output
!******************************************************************************
 
  do iPhi=1,nPhi
     do iTheta=1,nTheta
        ux(iPhi,iTheta) =  & 
             uExBtheta_C(iPhi,iTheta)*cos(Theta_G(iPhi,iTheta)) &
             * cos(Phi_G(iPhi,iTheta)) &
             - uExBphi_C(iPhi,iTheta)*sin(Phi_G(iPhi,iTheta))   &
             + uExBr_C(iPhi,iTheta)  *sin(Theta_G(iPhi,iTheta)) &
             * cos(Phi_G(iPhi,iTheta))
        
        uy(iPhi,iTheta) =  &
             uExBtheta_C(iPhi,iTheta)*cos(Theta_G(iPhi,iTheta)) &
             * sin(Phi_G(iPhi,iTheta)) &
             + uExBphi_C(iPhi,iTheta)*cos(Phi_G(iPhi,iTheta))   &
             + uExBr_C(iPhi,iTheta)  *sin(Theta_G(iPhi,iTheta)) &
             * sin(Phi_G(iPhi,iTheta))
        
        uz(iPhi,iTheta) =  &
             -uExBtheta_C(iPhi,iTheta)*sin(Theta_G(iPhi,iTheta)) &
             + uExBr_C(iPhi,iTheta)  *cos(Theta_G(iPhi,iTheta)) 
        
        x(iPhi,iTheta)  =  &
             1.0*sin(Theta_G(iPhi,iTheta))*cos(Phi_G(iPhi,iTheta))
        
        y(iPhi,iTheta)  =  &
             1.0*sin(Theta_G(iPhi,iTheta))*sin(Phi_G(iPhi,iTheta))
        
        z(iPhi,iTheta)  =  &
             1.0*cos(Theta_G(iPhi,iTheta))
        Lat = (0.5*cPi-Theta_G(iPhi,iTheta))*cRadToDeg
        Lon = Phi_G(iPhi,iTheta)*cRadToDeg
        Ex(iPhi,iTheta) =  & 
             Etheta_C(iPhi,iTheta)*cos(Theta_G(iPhi,iTheta)) &
             * cos(Phi_G(iPhi,iTheta))                     &
             - Ephi_C(iPhi,iTheta)*sin(Phi_G(iPhi,iTheta))   &
             + Er_C(iPhi,iTheta)  *sin(Theta_G(iPhi,iTheta)) &
             * cos(Phi_G(iPhi,iTheta))
        
        Ey(iPhi,iTheta) =  &
             Etheta_C(iPhi,iTheta)*cos(Theta_G(iPhi,iTheta)) &
             * sin(Phi_G(iPhi,iTheta))                     &
             + Ephi_C(iPhi,iTheta)*cos(Phi_G(iPhi,iTheta))   &
             + Er_C(iPhi,iTheta)  *sin(Theta_G(iPhi,iTheta)) &
             * sin(Phi_G(iPhi,iTheta))
        
        Ez(iPhi,iTheta) =  &
             -Etheta_C(iPhi,iTheta)*sin(Theta_G(iPhi,iTheta))&
             + Er_C(iPhi,iTheta)  *cos(Theta_G(iPhi,iTheta)) 
     enddo
  enddo
 
  TimeOut=int(Time)
  write(NameElectrodynamics,"(a,i8.8,a)") &
       'PW/Electrodynamics_Time',TimeOut,'.dat'
  open(UnitTmp_,FILE=NameElectrodynamics)
  write(UnitTmp_,*) &
     'VARIABLES = "X", "Y", "Z", "Ux [m/s]", "Uy [m/s]", "Uz [m/s]", '&
     // '"Pot [volts]", "Ex", "Ey", "Ez", "Jr", "Eflux [Ergs/cm2/s]", "AvE [keV]"'
  write(UnitTmp_,*) 'Zone I=', nPhi, ', J=', nTheta,', DATAPACKING=POINT'
  
  do iTheta=1,nTheta
     do iPhi=1,nPhi
        write(UnitTmp_,*) &
             x(iPhi,iTheta),y(iPhi,iTheta),z(iPhi,iTheta),     &
             ux(iPhi,iTheta),uy(iPhi,iTheta),uz(iPhi,iTheta),  &
             Potential_G(iPhi,iTheta),                         &
             Ex(iPhi,iTheta), Ey(iPhi,iTheta), Ez(iPhi,iTheta),&
             Jr_G(iPhi,iTheta),Eflux_G(iPhi,iTheta),AvE_G(iPhi,iTheta)
        
     enddo
  enddo
  
  close(UnitTmp_)

    deallocate(&
       Ux, &
       Uy, &
       Uz, &
       x , &
       y , &
       z , &
       Ex, &
       Ey, &
       Ez)

end subroutine PW_print_electrodynamics
