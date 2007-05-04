
SUBROUTINE PW_set_upper_bc
  use ModCommonVariables
  use ModGmPressure
  use ModPWOM, ONLY: iLine
  use ModConst,ONLY: cBoltzmann

  real :: ScaleHeightO,ScaleHeightH,ScaleHeightHe
  !----------------------------------------------------------------------------

  if (TypeSolver == 'Godunov') then
     
     UBGNDO=UOXYG(NDIM)
     UBGNHE=UHEL(NDIM)
     UBGNDH=UHYD(NDIM)
     
     TBGNDO=TOXYG(NDIM)
     TBGNHE=THEL(NDIM)
     TBGNDH=THYD(NDIM)
     
     PBGNDO = (POXYG(nDim-1)+ POXYG(nDim)- POXYG(nDim-2))*.8
     PBGNDH = (PHYD(nDim-1) + PHYD(nDim) - PHYD(nDim-2)) *.8
     PBGNHE = (PHEL(nDim-1) + PHEL(nDim) - PHEL(nDim-2)) *.8
     
     DBGNDO=PBGNDO/RGASO/TBGNDO
     DBGNHE=PBGNHE/RGASHE/TBGNHE
     DBGNDH=PBGNDH/RGASH/TBGNDH
  ENDIF
  
  if (TypeSolver == 'Rusanov') then
     
     TBGNDO=TOXYG(NDIM)
     TBGNHE=THEL(NDIM)
     TBGNDH=THYD(NDIM)
     
     TBGNDO2=TBGNDO
     TBGNHE2=TBGNHE
     TBGNDH2=TBGNDH
     
     !alternative pressure condition using free expansion
     !      PBGNDO=POXYG(NDIM)*(CellVolume_C(nDim)/CellVolumeTop(1)  
     !     &     *                   (Rad(nDim)/(Rad(nDim)+DrBnd))**3.0)
     !      PBGNHE=PHEL(NDIM)*(CellVolume_C(nDim)/CellVolumeTop(1)  
     !     &     *                   (Rad(nDim)/(Rad(nDim)+DrBnd))**3.0)
     !      PBGNDH=PHYD(NDIM)*(CellVolume_C(nDim)/CellVolumeTop(1)  
     !     &     *                   (Rad(nDim)/(Rad(nDim)+DrBnd))**3.0)
     
     !      PBGNDO2=PBGNDO
     !      PBGNHE2=PBGNHE
     !      PBGNDH2=PBGNDH
     
     !      PBGNDO2=PBGNDO*(CellVolumeTop(1)/CellVolumeTop(2)  
     !     &     *        ((Rad(nDim)+DrBnd)/(Rad(nDim)+2.0*DrBnd))**3.0)
     !      PBGNHE2=PBGNHE*(CellVolumeTop(1)/CellVolumeTop(2)  
     !     &     *        ((Rad(nDim)+DrBnd)/(Rad(nDim)+2.0*DrBnd))**3.0)
     !      PBGNDH2=PBGNDH*(CellVolumeTop(1)/CellVolumeTop(2)  
     !     &     *        ((Rad(nDim)+DrBnd)/(Rad(nDim)+2.0*DrBnd))**3.0)
     
     ! cBoltzmann [cgs] = 1.0e7 * cBoltzmann [SI] 
     ScaleHeightO =1.0e7*cBoltzmann*(TOXYG(nDim)+Telect(nDim))&
          /(abs(Gravty(nDim))*xmso)
     
     ScaleHeightH =1.0e7*cBoltzmann*(THYD(nDim)+Telect(nDim))&
          /(abs(Gravty(nDim))*xmsh)
     
     if(NamePlanet =='Earth ')&
          ScaleHeightHe=1.0e7*cBoltzmann*(THEL(nDim)+Telect(nDim))&
          /(abs(Gravty(nDim))*xmshe)
     
     PBGNDO=POXYG(nDim)*exp(-DrBnd/ScaleHeightO)
     PBGNDH=PHYD(nDim) *exp(-DrBnd/ScaleHeightH)
     if(NamePlanet =='Earth ')&
          PBGNHE=PHEL(nDim) *exp(-DrBnd/ScaleHeightHe)
     
     ! PBGNDO2=PBGNDO*exp(-DrBnd/ScaleHeightO)
     ! PBGNHE2=PBGNHE*exp(-DrBnd/ScaleHeightH)
     ! PBGNDH2=PBGNDH*exp(-DrBnd/ScaleHeightHe)
     
     PBGNDO2=PBGNDO
     if(NamePlanet =='Earth ')PBGNHE2=PBGNHE
     PBGNDH2=PBGNDH
     
     DBGNDO=PBGNDO/RGASO/TBGNDO
     if(NamePlanet =='Earth ')DBGNHE=PBGNHE/RGASHE/TBGNHE
     DBGNDH=PBGNDH/RGASH/TBGNDH
     
     DBGNDO2=PBGNDO2/RGASO/TBGNDO2
     if(NamePlanet =='Earth ')DBGNHE2=PBGNHE2/RGASHE/TBGNHE2
     DBGNDH2=PBGNDH2/RGASH/TBGNDH2         
     
     UBGNDO=UOXYG(nDim)
     if(NamePlanet =='Earth ')UBGNHE=UHEL(nDim)
     UBGNDH=UHYD(nDim)
     
     UBGNDO2=UBGNDO
     if(NamePlanet =='Earth ')UBGNHE2=UBGNHE
     UBGNDH2=UBGNDH

!      UBGNDO=DOXYG(nDim)/DBGNDO*AR12TOP(1)/AR23(nDim)*UOXYG(nDim)
!      UBGNHE=DHEL(nDim) /DBGNHE*AR12TOP(1)/AR23(nDim)*UHEL(nDim)
!      UBGNDH=DHYD(nDim) /DBGNDH*AR12TOP(1)/AR23(nDim)*UHYD(nDim)
!      
!      UBGNDO2=DBGNDO/DBGNDO2*AR12Top(2)/AR23Top(1)*UBGNDO
!      UBGNHE2=DBGNHE/DBGNHE2*AR12Top(2)/AR23Top(1)*UBGNHE
!      UBGNDH2=DBGNDH/DBGNDH2*AR12Top(2)/AR23Top(1)*UBGNDH
     
  endif
  
  
      
  
  XHTM=1.+ELFXIN*EXP(-(TIME-300.)**2/2./150./150.)
  
  ! Top boundary electron temperature is calculated from
  ! Te(nDim) + Tgradient*Dx where Tgradient = Etop/kappa
  ! and Etop is the electron heat flux or proportional to that
  TBGNDE=TELECT(nDim)+XHTM*ETOP/TELECT(nDim)**2.5
  
  ! If using GM --> PW coupling, overwrite ghostcell pressure
  ! with fraction of MHD pressure using 
  ! P_total = P_oxyg+P_e + smaller terms.
  ! n_e ~= n_oxyg hence P_oxyg ~= T_oxyg/(T_oxyg+T_e) * P_total
  if (UseGmToPw) then
     PBGNDO = (TBGNDO / (TBGNDO+TBGNDE))*p_I(iLine)
     DBGNDO = PBGNDO/RGASO/TBGNDO
  endif
  
  ! n_e   = n_H+  + n_O+ + n_He+
  ! rho_e = rho_H*(m_e/m_H) + rho_O*(m_e/m_O) + rho_He*(m_e/m_He)
  
  
  
  DBGNDE=RTHDEL*DBGNDH+RTOXEL*DBGNDO
  if(NamePlanet =='Earth ') DBGNDE=DBGNDE+RTHEEL*DBGNHE
  PBGNDE=RGASE*TBGNDE*DBGNDE
  
  UBGNDE=(RTHDEL*DBGNDH*UBGNDH+RTOXEL*DBGNDO*UBGNDO &
       +RTHEEL*DBGNHE*UBGNHE-1.8965E-18*CURRMX)/DBGNDE

  !these are sound speeds      
  WBGNDO=SQRT(GAMMA*PBGNDO/DBGNDO)
  if(NamePlanet =='Earth ') WBGNHE=SQRT(GAMMA*PBGNHE/DBGNHE)
  WBGNDH=SQRT(GAMMA*PBGNDH/DBGNDH)
  WBGNDE=SQRT(GAMMA*PBGNDE/DBGNDE)

  if(NamePlanet =='Earth ') then
     TCBGO=HLPO*(DBGNDO/DBGNDE)*TBGNDO**2.5
     TCBGE=HLPE*TBGNDE**2.5
     TCBGH=HLPH*(DBGNDH/DBGNDE)*TBGNDH**2.5
     THLP=SQRT(TBGNDH)*(1.-0.047*ALOG10(TBGNDH))**2
     TKLP=DBGNDO/(TBGNDO+16.*TBGNDH)**1.5
     TCBGH=TCBGH/(1.+(0.7692*(CZHN2+CZHO2)+1.0962*CZHO*THLP)/CZHOX &
          /TKLP)
     THLP=0.5*(XTNMAX+TBGNHE)
     THLP=SQRT(THLP)*(1.-0.093*ALOG10(THLP))**2
     TKLP=DBGNDO/(TBGNDO+4.*TBGNHE)**1.5
     TLLP=DBGNDH/(TBGNHE+4.*TBGNDH)**1.5
     TCBGHE=HLPHE*(DBGNHE/XMSHE)*TBGNHE
     TCBGHE=TCBGHE/(0.99*CZHEN2+0.99*CZHEO2+1.02*CZHEO+1.48*CZHEHE*&
          THLP+2.22*CZHEH+1.21*CZHEOX*TKLP+2.23*CZHEHD*TLLP)      
  else
     TCBGO=HLPO*(DBGNDO/DBGNDE)*TBGNDO**2.5
     TCBGE=HLPE*TBGNDE**2.5
     TCBGH=HLPH*(DBGNDH/DBGNDE)*TBGNDH**2.5
     TCBGHE=HLPHE*(DBGNHE/XMSHE)*TBGNHE**2.5
  endif
  RETURN
END SUBROUTINE PW_set_upper_bc
