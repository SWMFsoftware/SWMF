!******************************************************************************
!  This subroutine gets the electrodynamic parameters, and gets the convection
!  velocity. 
!******************************************************************************
subroutine PW_get_electrodynamics
  use ModNumConst, ONLY: cTwoPi
  use ModIoUnit, ONLY: UnitTmp_
  use ModPWOM
  use ModNumConst, ONLY:cDegToRad
  use ModAurora , ONLY: set_aurora
  use ModCommonVariables,ONLY:Ap
  implicit none

  character (len=100), dimension(100):: Lines_I
  integer :: iError
  logical,save :: IsFirst = .true.
  real :: dTheta1, dPhi1
  !---------------------------------------------------------------------------
  
  if ((IsStandAlone .or. .not. UseIE) .and. .not.UseWeimer) then
     open(UnitTmp_, FILE=NamePhiNorth)  
     if(IsFirst)then
        call allocate_ie_variables(257, 65)
        do iPhi=1,nPhi
           do iTheta=1,nTheta
              read(unit=UnitTmp_,fmt='(6(1PE13.5))') &
                   Theta_G(iPhi,iTheta),Phi_G(iPhi,iTheta),SigmaH_G(iPhi,iTheta),&
                   SigmaP_G(iPhi,iTheta),Jr_G(iPhi,iTheta),Potential_G(iPhi,iTheta)
              
           enddo
        enddo
        !  Change angles to radians
        Theta_G(1:nPhi,1:nTheta) = Theta_G(1:nPhi,1:nTheta)*cDegToRad
        Phi_G  (1:nPhi,1:nTheta) = Phi_G  (1:nPhi,1:nTheta)*cDegToRad
        !  Convert potential from kilovolts to Volts
        Potential_G(1:nPhi,1:nTheta) = Potential_G(1:nPhi,1:nTheta)*1.0e3
        !  Convert microA/m^2 --> A/m^2
        Jr_G(1:nPhi,1:nTheta) = Jr_G(1:nPhi,1:nTheta) * 1.0e-6 
        close(UnitTmp_)   
     endif
  elseif (UseWeimer) then
     call get_weimer_potential
  endif
  
  ! Set EIE for Aurora when weimer not used
  if (.not.UseWeimer .and. UseAurora .and. IsFirst) then
          Lines_I(1) = "#BACKGROUND"
          Lines_I(2) = "EIE/"
          Lines_I(3) = "zero"
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
     if (iError /= 0) then
        write(*,*) 'PW_ERROR: EIE_Initialize failed at get_electrodynamic'
        call con_stop()
     endif

  end if
  
  ! Set AvE_G and Eflux_G when using empirical Aurora
  if (UseAurora) then
     call set_aurora
     AvE_G  (1:nPhi,1:nTheta) = ElectronAverageEnergy_C(1:nPhi,1:nTheta)
     Eflux_G(1:nPhi,1:nTheta) = ElectronEnergyFlux_C   (1:nPhi,1:nTheta)
  end if

!******************************************************************************
!  Calc Bfield components
!******************************************************************************
  OmegaPlanet    = cTwoPi/24.0/3600.
  rPlanet        = 6378000.0
  rLowerBoundary = rPlanet+110.0e3
  MagMoment      = 7.84e15
  Bcoef          = MagMoment/(rLowerBoundary)**3 
  do iPhi=1,nPhi
     do iTheta=1,nTheta
        Br_G(iPhi,iTheta)     = -Bcoef*2.0*cos(Theta_G(iPhi,iTheta))
        Btheta_G(iPhi,iTheta) = -Bcoef*sin(Theta_G(iPhi,iTheta))
    enddo
  enddo
  
  BmagnitudeSquared_G(1:nPhi,1:nTheta) = &
       Br_G(1:nPhi,1:nTheta)**2 + Btheta_G(1:nPhi,1:nTheta)**2
!******************************************************************************
!  Fill Ghost cells
!******************************************************************************
  Theta_G    (nPhi+1,1:nTheta) = Theta_G    (2,1:nTheta) 
  Phi_G      (nPhi+1,1:nTheta) = Phi_G      (2,1:nTheta)
  SigmaH_G   (nPhi+1,1:nTheta) = SigmaH_G   (2,1:nTheta)
  SigmaP_G   (nPhi+1,1:nTheta) = SigmaP_G   (2,1:nTheta) 
  Jr_G       (nPhi+1,1:nTheta) = Jr_G       (2,1:nTheta)
  Potential_G(nPhi+1,1:nTheta) = Potential_G(2,1:nTheta)

  Theta_G    (0,1:nTheta) = Theta_G    (nPhi-1,1:nTheta) 
  Phi_G      (0,1:nTheta) = Phi_G      (nPhi-1,1:nTheta)
  SigmaH_G   (0,1:nTheta) = SigmaH_G   (nPhi-1,1:nTheta)
  SigmaP_G   (0,1:nTheta) = SigmaP_G   (nPhi-1,1:nTheta) 
  Jr_G       (0,1:nTheta) = Jr_G       (nPhi-1,1:nTheta)
  Potential_G(0,1:nTheta) = Potential_G(nPhi-1,1:nTheta)

  Theta_G    (0:nPhi+1,nTheta+1) = Theta_G    (0:nPhi+1,nTheta) 
  Phi_G      (0:nPhi+1,nTheta+1) = Phi_G      (0:nPhi+1,nTheta)
  SigmaH_G   (0:nPhi+1,nTheta+1) = SigmaH_G   (0:nPhi+1,nTheta)
  SigmaP_G   (0:nPhi+1,nTheta+1) = SigmaP_G   (0:nPhi+1,nTheta) 
  Jr_G       (0:nPhi+1,nTheta+1) = Jr_G       (0:nPhi+1,nTheta)
  Potential_G(0:nPhi+1,nTheta+1) = Potential_G(0:nPhi+1,nTheta)

  do iPhi=1,nPhi
     Theta_G    (iPhi,0) = Theta_G    (mod(iPhi+floor(nPhi/2.0),nPhi-1),2) 
     Phi_G      (iPhi,0) = Phi_G      (mod(iPhi+floor(nPhi/2.0),nPhi-1),2)
     SigmaH_G   (iPhi,0) = SigmaH_G   (mod(iPhi+floor(nPhi/2.0),nPhi-1),2)
     SigmaP_G   (iPhi,0) = SigmaP_G   (mod(iPhi+floor(nPhi/2.0),nPhi-1),2) 
     Jr_G       (iPhi,0) = Jr_G       (mod(iPhi+floor(nPhi/2.0),nPhi-1),2)
     Potential_G(iPhi,0) = Potential_G(mod(iPhi+floor(nPhi/2.0),nPhi-1),2)
  enddo

  ! Fill ghost cells for AvE and Eflux from IE when UseIE=T, else set to 0
  if (UseIE .or. UseAurora) then
     AvE_G  (nPhi+1,1:nTheta)  = AvE_G  (2,1:nTheta)     
     Eflux_G(nPhi+1,1:nTheta)  = Eflux_G(2,1:nTheta)     
     AvE_G  (0,1:nTheta)       = AvE_G  (nPhi-1,1:nTheta)
     Eflux_G(0,1:nTheta)       = Eflux_G(nPhi-1,1:nTheta)
     AvE_G  (0:nPhi+1,nTheta+1)= AvE_G  (0:nPhi+1,nTheta)
     Eflux_G(0:nPhi+1,nTheta+1)= Eflux_G(0:nPhi+1,nTheta)
     do iPhi=1,nPhi
        AvE_G  (iPhi,0) = AvE_G  (mod(iPhi+floor(nPhi/2.0),nPhi-1),2)
        Eflux_G(iPhi,0) = Eflux_G(mod(iPhi+floor(nPhi/2.0),nPhi-1),2)
     enddo
  else
     AvE_G (:,:) = 0.0
     Eflux_G(:,:) = 0.0
  endif

!******************************************************************************
!  Calc electric field from E=-grad Potential
!******************************************************************************
  do iPhi=1,nPhi
     do iTheta=1,nTheta
        DTheta1=abs(Theta_G(iPhi,iTheta)-Theta_G(iPhi,iTheta-1))
        DPhi1=abs(Phi_G(iPhi,iTheta)    -Phi_G(iPhi-1,iTheta))
        Etheta_C(iPhi,iTheta) = &
             (Potential_G(iPhi,iTheta)-Potential_G(iPhi,iTheta+1))&
             /((rLowerBoundary)*DTheta1)
        if (iTheta == 1) then
           Ephi_C(iPhi,iTheta) = 0.0
        else
           Ephi_C(iPhi,iTheta)   = &
                (Potential_G(iPhi-1,iTheta)-Potential_G(iPhi+1,iTheta))&
                / ((rLowerBoundary)*2.0*DPhi1  &
                * sin(Theta_G(iPhi,iTheta)))
        endif
    enddo
  enddo
  Er_C(:,:) = 0.0

!******************************************************************************
!  Calc VelocityExB drift from E and B, add corotation velocity to uExBphi_C
!******************************************************************************
  

  do iPhi=1,nPhi
     do iTheta=1,nTheta
        uExBtheta_C(iPhi,iTheta) = &
             (Ephi_C(iPhi,iTheta)*Br_G(iPhi,iTheta)) &
             / BmagnitudeSquared_G(iPhi,iTheta)
        
        uExBphi_C(iPhi,iTheta)   = &
             (-Etheta_C(iPhi,iTheta)*Br_G(iPhi,iTheta)) &
             / BmagnitudeSquared_G(iPhi,iTheta) &
             + OmegaPlanet*rPlanet*sin(Theta_G(iPhi,iTheta)) !corotation
        
        uExBr_C(iPhi,iTheta)     = &
             (-Ephi_C(iPhi,iTheta)*Btheta_G(iPhi,iTheta)) &
             / BmagnitudeSquared_G(iPhi,iTheta)


     enddo
  enddo
  uExBtheta_C(:,1) = uExBtheta_C(:,2)
  uExBphi_C  (:,1) = uExBphi_C  (:,2)
  
!  if (.not.DoMoveLine) Then 
!     uExBtheta_C(:,:) = 0.0
!     uExBphi_C  (:,:) = 0.0
!     uExBr_C    (:,:) = 0.0
!  endif

  if (.not.UseJr) Then 
     Jr_G(:,:) = 0.0
  endif

  IsFirst = .false.

  !set auroral heating 
!  call set_theta0(nPhi,nTheta,uExBphi_C,uExBtheta_C,Theta_G)
!  call set_Emax(Ap(1))
end subroutine PW_get_electrodynamics
