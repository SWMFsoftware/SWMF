Module ModAurora
  implicit none
  save
  
  real :: Theta0,dtheta=0.174532925199 !10 degrees
  real,parameter :: E0=30.0e-3 !ergs/cm**2
  real :: Emax


contains
  !============================================================================
  
  real function get_etop(GmLat,GmLon)
    !returns topside electron energy input
    use ModNumConst, ONLY:cPi,cDegToRad
    
    real,intent(in) :: GmLat,GmLon
    real :: Theta,Phi,DeltaTheta,DeltaPhi
    real :: eTopMax
    !--------------------------------------------------------------------------
    Theta=cPi/2-GmLat*cDegToRad
    Phi = GmLon*cDegToRad
    eTopMax = 10.0*E0
    ! If point is inside aurora return auroral etop value, else return default
    
    if(Theta > Theta0-dtheta/2 .and. Theta < Theta0+dtheta/2 &
         .and. phi > 0.5*cPi .and. phi < 1.5*cPi) then
       DeltaTheta = abs(Theta-Theta0)
       DeltaPhi = abs(Phi - cPi)
       get_etop = &
            (1 - 2.0*DeltaPhi/cPi)*(eTopMax-(eTopMax-E0)*DeltaTheta/dtheta)
       
    else
       get_etop = E0
    endif
    
  end function get_etop
  
  !============================================================================
  !============================================================================
  
  real function get_eflux(GmLat,GmLon)
    !returns topside electron energy input
    use ModNumConst, ONLY:cPi,cDegToRad
    
    real,intent(in) :: GmLat,GmLon
    real :: Theta,Phi,DeltaTheta,DeltaPhi
    !--------------------------------------------------------------------------
    Theta=cPi/2-GmLat*cDegToRad
    Phi = GmLon*cDegToRad
    
    ! If point is inside aurora return auroral etop value, else return default
    
    if(Theta > Theta0-dtheta/2 .and. Theta < Theta0+dtheta/2 &
         .and. phi > 0.5*cPi .and. phi < 1.5*cPi) then
       DeltaTheta = abs(Theta-Theta0)
       DeltaPhi = abs(Phi - cPi)
       get_eflux = (1 - 2.0*DeltaPhi/cPi)*(Emax-(Emax-E0)*DeltaTheta/dtheta)
       
    else
       get_eflux = E0
    endif
    
  end function get_eflux
  
  !============================================================================
  !============================================================================
  
  real function get_eAverageE(GmLat,GmLon)
    !returns average precipitating energy
    use ModNumConst, ONLY:cPi,cDegToRad
    
    real,intent(in) :: GmLat,GmLon
    real :: Theta,Phi,DeltaTheta,DeltaPhi
    real,parameter:: AveE0=.2, AveEmax=.75
    !--------------------------------------------------------------------------
    Theta=cPi/2-GmLat*cDegToRad
    Phi = GmLon*cDegToRad
    
    ! If point is inside aurora return auroral etop value, else return default
    
    if(Theta > Theta0-dtheta/2 .and. Theta < Theta0+dtheta/2 &
         .and. phi > 0.5*cPi .and. phi < 1.5*cPi) then
       DeltaTheta = abs(Theta-Theta0)
       DeltaPhi = abs(Phi - cPi)
       get_eAverageE = &
            (1 - 2.0*DeltaPhi/cPi)*(AveEmax-(AveEmax-AveE0)*DeltaTheta/dtheta)
       
    else
       get_eAverageE = AveE0
    endif
    
  end function get_eAverageE
  
  !============================================================================
  
  subroutine set_Emax(Ap)
    real,intent(in) :: Ap
    !--------------------------------------------------------------------------
    
    Emax=1.0+7.0*Ap/207.0
  end subroutine set_Emax
  
  !============================================================================
  
  subroutine set_theta0(nPhi,nTheta,uExBphi_C,uExBtheta_C,Theta_G)
    use ModPWOM, ONLY: DoMoveLine
    integer,intent(in) :: nPhi,nTheta
    real   ,intent(in) :: uExBphi_C(nPhi,nTheta),uExBtheta_C(nPhi,nTheta)
    real   ,intent(in) :: Theta_G(0:nPhi+1,0:nTheta+1)
    real,parameter :: AlphaMax = 0.5            !30 degrees  
    real,parameter :: ThetaMin = 0.174532925199 !10 degrees  
    integer:: iPhi,iTheta
    real :: Alpha
    
    !--------------------------------------------------------------------------
    iPhi = nPhi/2
    do iTheta=nTheta,1,-1
       Theta0=Theta_G(iPhi,iTheta)
       if (.not. DoMoveLine) return
       Alpha = atan(uExBtheta_C(iPhi,iTheta)/uExBphi_C(iPhi,iTheta))
       
       if (Alpha >= AlphaMax .or. Theta0 <= ThetaMin ) then
          return
       endif
    enddo
  end subroutine set_theta0
  
  !============================================================================
end Module ModAurora
