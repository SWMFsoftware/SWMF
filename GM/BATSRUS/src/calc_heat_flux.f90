!^CFG FILE DISSFLUX
!^CFG COPYRIGHT UM
subroutine add_heat_flux(DoResChangeOnly)
  use ModSize,     ONLY:nI,nJ,nK,gcn
  use ModMain,     ONLY:nIFace,nJFace,nKFace,&
       iMinFaceY,iMaxFaceY,iMinFaceZ,iMaxFaceZ, &
       jMinFaceX,jMaxFaceX,jMinFaceZ,jMaxFaceZ, &
       kMinFaceX,kMaxFaceX,kMinFaceY,kMaxFaceY, &
       globalBLK,UseSpitzerForm,x_,y_,z_,UseB0
  use ModVarIndexes,ONLY:rho_,P_,Energy_,&
       Bx_,By_,Bz_
  use ModGeometry, ONLY: dx_BLK, &
       dy_BLK,dz_BLK,fAx_BLK,fAy_BLK,fAz_BLK
  use ModParallel, ONLY:neiLeast,neiLwest,neiLsouth, &
       neiLnorth,neiLtop,neiLbot
  use ModAdvance,  ONLY:State_VGB, &
       B0_DGB, & 
       VdtFace_x,VdtFace_y,VdtFace_z, &
       Flux_VX,Flux_VY,Flux_VZ
  use ModConst,    ONLY:cBoltzmann,cProtonMass, &
       cElectronMass,cElectronCharge,cEps
  use ModNumConst, ONLY:cOne,cHalf,cZero,cTiny,cTwoPi,cPi
  use ModPhysics,  ONLY:Kappa0Heat,ExponentHeat,g,gm1, &
       No2Si_V, UnitX_, UnitT_, UnitEnergyDens_, &
       UnitTemperature_, UnitB_, UnitRho_, &
       UseDefaultUnits
  use ModMpi
  implicit none
  
  logical,intent(in):: DoResChangeOnly
  integer:: i,j,k
  real:: CU_x,CU_t,CU_temperature,CU_energydens
  real:: CristophCoefficient,HeatFluxCoefficient
  real:: PressureUp,PressureDown,DensityUp,DensityDown
  real:: KappaParHeat,KappaParConstHeat
  real:: KappaPerpHeat,KappaPerpConstHeat
  real:: LogLambdaHeat,Omega_iTau_ii2Heat
  real,dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn):: &
       BX2,BY2,BZ2,B2Inverted
  !\
  ! Introduce some units for dimensionalization:: 
  ! This is not necessary, but makes things more clear!
  !/
  CU_x = No2Si_V(UnitX_)                          ! in [m]
  CU_t = No2Si_V(UnitT_)                          ! in [s]
  if (UseDefaultUnits) then
     CU_temperature = No2Si_V(UnitTemperature_)   ! in [K]
  else
     CU_temperature = No2Si_V(UnitTemperature_)/g ! in [K]
  end if
  CU_energydens     = No2Si_V(UnitEnergydens_)    ! in [J/m^3]
  !
  if(UseB0)then
     BX2(:,:,:) = (State_VGB(Bx_,:,:,:,globalBLK)+&
          B0_DGB(x_,:,:,:,globalBLK))**2
     BY2(:,:,:) = (State_VGB(By_,:,:,:,globalBLK)+&
          B0_DGB(y_,:,:,:,globalBLK))**2
     BZ2(:,:,:) = (State_VGB(Bz_,:,:,:,globalBLK)+&
          B0_DGB(z_,:,:,:,globalBLK))**2
  else
     BX2(:,:,:) = State_VGB(Bx_,:,:,:,globalBLK)**2
     BY2(:,:,:) = State_VGB(By_,:,:,:,globalBLK)**2
     BZ2(:,:,:) = State_VGB(Bz_,:,:,:,globalBLK)**2
  end if
  B2Inverted = cOne/max(BX2+BY2+BZ2,cTiny**2)
  !\
  ! Compute and add the x_heat_flux to the x-face fluxes 
  !/
  CristophCoefficient = fAx_BLK(globalBLK)/dx_BLK(globalBLK)
  if (.not.DoResChangeOnly) then
     do k=kMinFaceX,kMaxFaceX
        do j=jMinFaceX,jMaxFaceX
           do i=1,nIFace
              call add_heat_flux_x
           end do
        end do
     end do
  else if (neiLeast(globalBLK)==+1) then
     i=1
     do k=1,nK
        do j=1,nJ
           call add_heat_flux_x
        end do
     end do
  else if ( neiLwest(globalBLK)==+1) then
     i=nIFace
     do k=1,nK
        do j=1,nJ
           call add_heat_flux_x
        end do
     end do
  end if
  !\  
  ! Compute and add the y_heat_flux to the y-face fluxes 
  !/
  CristophCoefficient = fAy_BLK(globalBLK)/dy_BLK(globalBLK)
  if (.not.DoResChangeOnly) then
     do k=kMinFaceY,kMaxFaceY
        do j=1,nJFace
           do i=iMinFaceY,iMaxFaceY
              call add_heat_flux_y
           end do
        end do
     end do
  else if(neiLsouth(globalBLK)==+1)then
     j=1
     do k=1,nK
        do i=1,nI
           call add_heat_flux_y
        end do
     end do
  else if (neiLnorth(globalBLK)==+1) then
     j=nJFace 
     do k=1,nK
        do i=1,nI
           call add_heat_flux_y
        end do
     end do
  end if
  !\
  ! Compute and add the z_heat_flux to the z-face fluxes  
  !/
  CristophCoefficient = fAz_BLK(globalBLK)/dz_BLK(globalBLK)
  if (.not.DoResChangeOnly) then
     do k=1,nKFace
        do j=jMinFaceZ,jMaxFaceZ
           do i=iMinFaceZ,iMaxFaceZ
              call add_heat_flux_z
           end do
        end do
     end do
  else if (neiLbot(globalBLK)==+1) then
     k=1
     do j=1,nJ
        do i=1,nI
           call add_heat_flux_z
        end do
     end do
  else if ( neiLtop(globalBLK)==+1) then
     k=nKFace            
     do j=1,nJ
        do i=1,nI
           call add_heat_flux_z
        end do
     end do
  end if
Contains
  
  subroutine add_heat_flux_x
    PressureUp   = State_VGB(P_,i  ,j,k,globalBLK)
    PressureDown = State_VGB(P_,i-1,j,k,globalBLK)
    DensityUp    = State_VGB(rho_,i  ,j,k,globalBLK)
    DensityDown  = State_VGB(rho_,i-1,j,k,globalBLK)
    ! Compute the parallel heat conduction coefficient
    call compute_kappa_coefficient
    ! Compute the heat flux coefficient
    HeatFluxCoefficient = CristophCoefficient*KappaParHeat * &
         min(max(BX2(i  ,j,k),cTiny**2)*B2Inverted(i  ,j,k), &
             max(BX2(i-1,j,k),cTiny**2)*B2Inverted(i-1,j,k))
    ! Subtract the heat flux from the x-face fluxes

    Flux_VX(Energy_,i,j,k) = &
         Flux_VX(Energy_,i,j,k)-HeatFluxCoefficient  * &
         (PressureUp/DensityUp-PressureDown/DensityDown)
    ! Update the CFL condition
    VdtFace_x(i,j,k) = VdtFace_x(i,j,k)+ &
         2*HeatFluxCoefficient*gm1  / &
         min(DensityUp,DensityDown)
    !\
    ! If Spitzer form is used, then add the cross-field conductive flux   
    !/
    if (UseSpitzerForm) then
       HeatFluxCoefficient = CristophCoefficient*KappaPerpHeat* &
            min(max(BY2(i  ,j,k)+BZ2(i  ,j,k),cTiny**2)       * &
                B2Inverted(i  ,j,k), &
                max(BY2(i-1,j,k)+BZ2(i-1,j,k),cTiny**2)       * &
                B2Inverted(i-1,j,k))
       ! Subtract the cross-field heat flux from the x-face fluxes 
       Flux_VX(Energy_,i,j,k) = &
            Flux_VX(Energy_,i,j,k)-HeatFluxCoefficient  * &
            (PressureUp/DensityUp-PressureDown/DensityDown)
       ! Update the CFL condition
       VdtFace_x(i,j,k) = VdtFace_x(i,j,k)+ &
            2*HeatFluxCoefficient*gm1  / &
            min(DensityUp,DensityDown)  
    end if
  end subroutine add_heat_flux_x

  subroutine add_heat_flux_y
    PressureUp   = State_VGB(P_,i,j  ,k,globalBLK)
    PressureDown = State_VGB(P_,i,j-1,k,globalBLK)
    DensityUp    = State_VGB(rho_,i,j  ,k,globalBLK)
    DensityDown  = State_VGB(rho_,i,j-1,k,globalBLK)
    ! Compute the parallel heat conduction coefficient
    call compute_kappa_coefficient
    ! Compute the heat flux coefficient
    HeatFluxCoefficient = CristophCoefficient*KappaParHeat * &
         min(max(BY2(i,j  ,k),cTiny**2)*B2Inverted(i,j  ,k), &
             max(BY2(i,j-1,k),cTiny**2)*B2Inverted(i,j-1,k))
    ! Subtract the heat flux from the y-face fluxes

    Flux_VY(Energy_,i,j,k) = &
         Flux_VY(Energy_,i,j,k)-HeatFluxCoefficient  * &
         (PressureUp/DensityUp-PressureDown/DensityDown)
    ! Update the CFL condition
    VdtFace_y(i,j,k) = VdtFace_y(i,j,k)+ &
         2*HeatFluxCoefficient*gm1  / &
         min(DensityUp,DensityDown)   
    !\
    ! If Spitzer form is used, then add the cross-field conductive flux   
    !/
    if (UseSpitzerForm) then
       HeatFluxCoefficient = CristophCoefficient*KappaPerpHeat* &
            min(max(BX2(i,j  ,k)+BZ2(i,j  ,k),cTiny**2)       * &
                B2Inverted(i,j  ,k), &
                max(BX2(i,j-1,k)+BZ2(i,j-1,k),cTiny**2)       * &
                B2Inverted(i,j-1,k))
       ! Subtract the cross-field heat flux from the y-face fluxes 
       Flux_VY(Energy_,i,j,k) = &
            Flux_VY(Energy_,i,j,k)-HeatFluxCoefficient  * &
            (PressureUp/DensityUp-PressureDown/DensityDown)
       ! Update the CFL condition
       VdtFace_y(i,j,k) = VdtFace_y(i,j,k)+ &
            2*HeatFluxCoefficient*gm1  / &
            min(DensityUp,DensityDown)  
    end if
  end subroutine add_heat_flux_y
  
  subroutine add_heat_flux_z
    PressureUp   = State_VGB(P_,i,j,k  ,globalBLK)
    PressureDown = State_VGB(P_,i,j,k-1,globalBLK)
    DensityUp    = State_VGB(rho_,i,j,k  ,globalBLK)
    DensityDown  = State_VGB(rho_,i,j,k-1,globalBLK)
    ! Compute the parallel heat conduction coefficient
    call compute_kappa_coefficient
    ! Compute the heat flux coefficient
    HeatFluxCoefficient = CristophCoefficient*KappaParHeat * &
         min(max(BZ2(i,j,k  ),cTiny**2)*B2Inverted(i,j,k  ), &
             max(BZ2(i,j,k-1),cTiny**2)*B2Inverted(i,j,k-1))
    ! Subtract the heat flux from the z-face fluxes 
    Flux_VZ(Energy_,i,j,k) = &
         Flux_VZ(Energy_,i,j,k)-HeatFluxCoefficient  * &
         (PressureUp/DensityUp-PressureDown/DensityDown)
    ! Update the CFL condition
    VdtFace_z(i,j,k) = VdtFace_z(i,j,k)+ &
         2*HeatFluxCoefficient*gm1  / &
         min(DensityUp,DensityDown)
    !\
    ! If Spitzer form is used, then add the cross-field conductive flux   
    !/
    if (UseSpitzerForm) then
       HeatFluxCoefficient = CristophCoefficient*KappaPerpHeat* &
            min(max(BX2(i,j,k  )+BY2(i,j,k  ),cTiny**2)       * &
                B2Inverted(i,j,k  ), &
                max(BX2(i,j,k-1)+BY2(i,j,k-1),cTiny**2)       * &
                B2Inverted(i,j,k-1))
       ! Subtract the cross-field heat flux from the z-face fluxes 
       Flux_VZ(Energy_,i,j,k) = &
            Flux_VZ(Energy_,i,j,k)-HeatFluxCoefficient  * &
            (PressureUp/DensityUp-PressureDown/DensityDown)
       ! Update the CFL condition
       VdtFace_z(i,j,k) = VdtFace_z(i,j,k)+ &
            2*HeatFluxCoefficient*gm1  / &
            min(DensityUp,DensityDown)  
    end if
  end subroutine add_heat_flux_z
  
  subroutine compute_kappa_coefficient
    if (UseSpitzerForm) then
       !\
       ! Assume ln(Lambda) = 20 (solar corona) -- coded
       !                  or 30 (magnetosphere)
       ! _________________or 10 (lab. plasma)
       !/
       LogLambdaHeat = 20.0 ! Coronal value
       !\
       ! KappaParConstHeat is in units of [W/m/K]=[m^2*(J/m^3)/(s*K)]
       ! KappaParConstHeat = 1.820173178845863049E-10
       !/
       KappaParConstHeat = 3.16*3*(sqrt(cElectronMass)   / &
            cElectronCharge**2)*(cTwoPi*cBoltzmann)**(cOne+cHalf)/ &
            cElectronCharge**2*cEps**2*cBoltzmann                  / &
            cElectronMass*cBoltzmann
       !\
       ! Take into account the Coulomb logarithm::
       ! KappaParConstHeat  = KappaParConstHeat/LogLambdaHeat
       !/
       KappaParConstHeat  = KappaParConstHeat/LogLambdaHeat
       !\
       ! Finally, compute KappaParHeat::
       ! KappaParHeat = KappaParConstHeat*Te^2.5
       !/
       KappaParHeat = KappaParConstHeat*(No2Si_V(UnitTemperature_)* &
            (PressureUp+PressureDown)                      / &
            (DensityUp+DensityDown))**2.5
       ! Dimensionalize KappaParHeat::
       KappaParHeat = KappaParHeat*CU_t*CU_temperature     / &
            (CU_energydens*CU_x**2)
       !\
       ! Take into account the cross-field heat conduction.
       ! This is mostly due to i-i collisions::
       ! KappaPerpHeat = KappaPerpConstHeat*(Ti^2.5/LogLambdaHeat)/
       ! (1+Omega_iTau_ii2Heat), where
       ! Omega_iTau_ii2Heat = 2.552789521363004029E+30*B^2*T^3/
       ! (Ni*LogLambdaHeat)
       !/
       !\
       ! KappaPerpConstHeat is in units of [W/m/K]=[m^2*(J/m^3)/(s*K)]
       ! KappaPerpConstHeat = 3.802064132280014164E-12
       !/
       KappaPerpConstHeat = 6*(sqrt(2*cProtonMass)/ &
            cElectronCharge**2)*(cTwoPi*cBoltzmann)**1.5 / &
            cElectronCharge**2*cEps**2*cBoltzmann                   / &
            cProtonMass*cBoltzmann
       !\
       ! Take into account the Coulomb logarithm::
       ! KappaPerpConstHeat = KappaPerpConstHeat/LogLambdaHeat
       !/
       KappaPerpConstHeat = KappaPerpConstHeat/LogLambdaHeat
       !\
       ! Compute Omega_i*Tau_ii2::
       ! Omega_i*Tau_ii2 = 2.552789521363004029E+30
       !/
       Omega_iTau_ii2Heat = (3*(sqrt(2/cProtonMass)/ &
            cElectronCharge)*(cTwoPi*cBoltzmann)**1.5* &
            cEps**2/cElectronCharge**2)**2
       !\
       ! Omega_i*Tau_ii2 = Omega_i*Tau_ii2*[B/(Ni*LogLambdaHeat)]^2*Ti^3
       !                        or
       ! Omega_i*Tau_ii2 = Omega_i*Tau_ii2*[mp*B/(rho*LogLambdaHeat)]^2*Ti^3
       !/
       Omega_iTau_ii2Heat = max(Omega_iTau_ii2Heat*cProtonMass**2* &
            (No2Si_V(UnitB_)**2*(BX2(i,j,k)+BY2(i,j,k)+BZ2(i,j,k)))     * &
            (No2Si_V(UnitTemperature_)*PressureUp/DensityUp)**3         / &
            (No2Si_V(UnitRho_)*DensityUp*LogLambdaHeat)**2,cTiny**2)
       !\
       ! Finally, compute KappaPerpHeat::
       ! KappaPerpHeat = KappaPerpConstHeat*Ti^2.5/(1+Omega_iTau_ii2Heat)
       !/
       !\
       ! Note:: for LogLambdaHeat = 20, (KappaPerpHeat/KappaParHeat) = 
       ! 3.273491465571024622E-30*(Ni/B)^2/Ti^3
       !/       
       KappaPerpHeat = KappaPerpConstHeat*(No2Si_V(UnitTemperature_)  * &
            (PressureUp+PressureDown)                          / &
            (DensityUp+DensityDown))**2.5             / &
            (cOne+Omega_iTau_ii2Heat)
       ! Dimensionalize KappaPerpHeat::
       KappaPerpHeat = KappaPerpHeat*CU_t*CU_temperature       / &
            (CU_energydens*CU_x**2)
    else
       !\
       ! Kappa0Heat is in units of [W/m/K]=[m^2*(J/m^3)/(s*K)]
       !/
      KappaParHeat  = Kappa0Heat*(No2Si_V(UnitTemperature_)   * &
            (PressureUp+PressureDown)                  / &
            (DensityUp+DensityDown))**ExponentHeat
       ! Dimensionalize KappaParHeat::
       KappaParHeat  = KappaParHeat*CU_t*CU_temperature/ &
            (CU_energydens*CU_x**2)
       ! Neglect KappaPerpHeat::
       KappaPerpHeat = cZero
    end if
  end subroutine compute_kappa_coefficient
  
end subroutine add_heat_flux
