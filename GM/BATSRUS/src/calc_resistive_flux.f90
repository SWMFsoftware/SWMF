!^CFG FILE DISSFLUX
!^CFG COPYRIGHT UM
Module ModResist
  use ModSize,     ONLY:nI,nJ,nK,gcn,nBLK
  implicit none
  
  logical :: UseResistivity=.false.
  logical :: DoInitEtaLocResist_B(nBLK)
  data DoInitEtaLocResist_B /nBLK*.true./
  logical, dimension (1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn,nBLK) :: &
       IsNotIgnitedYetResist_GB
  integer, parameter :: MaxGhostPE=(nI+2*gcn)*(nJ+2*gcn)*(nK+2*gcn)*nBLK
  data IsNotIgnitedYetResist_GB /MaxGhostPE*.true./
  real, dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn,nBLK):: &
       EtaLocResist_GB
  real :: EtaResist_G(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn)
  real :: Eta_GB(-1:nI+2,-1:nJ+2,-1:nK+2,nBlk)

  logical :: UseUserEta = .false., DoSetEta=.true.
end module ModResist

!==============================================================================

subroutine add_resistive_flux(DoResChangeOnly)
  use ModSize,     ONLY:nI,nJ,nK,gcn,nBLK
  use ModMain,     ONLY:nIFace,nJFace,nKFace,&
       iMinFaceY,iMaxFaceY,iMinFaceZ,iMaxFaceZ, &
       jMinFaceX,jMaxFaceX,jMinFaceZ,jMaxFaceZ, &
       kMinFaceX,kMaxFaceX,kMinFaceY,kMaxFaceY, &
       globalBLK,BLKtest,UseSpitzerForm, &
       X_, Y_, Z_
  use ModVarIndexes,ONLY:Bx_,By_,Bz_,&
       rho_,Energy_
  use ModGeometry, ONLY:x_BLK,y_BLK,z_BLK,dx_BLK, &
       dy_BLK,dz_BLK,fAx_BLK,fAy_BLK,fAz_BLK
  use ModParallel, ONLY:neiLeast,neiLwest,neiLsouth, &
       neiLnorth,neiLtop,neiLbot
  use ModAdvance,  ONLY:State_VGB, &
       B0xCell_BLK,B0yCell_BLK,B0zCell_BLK, & 
       VdtFace_x,VdtFace_y,VdtFace_z, &
       Flux_VX,Flux_VY,Flux_VZ
  use ModNumConst, ONLY:cOne,cTwo,cFour,cHalf, &
       cZero,cTiny,cHundred,cHundredth,cPi
  use ModPhysics,  ONLY: gm1
  use ModResist,   ONLY: EtaResist_G ,Eta_GB
  use ModMpi
  use ModCovariant, ONLY: UseCovariant, FaceArea2MinI_B,& !^CFG IF COVARIANT    
       FaceArea2MinJ_B,FaceArea2MinK_B, &                 !^CFG IF COVARIANT    
       FaceAreaI_DFB,FaceAreaJ_DFB,FaceAreaK_DFB          !^CFG IF COVARIANT  
  implicit none
  
  logical, intent(in):: DoResChangeOnly
  integer:: i,j,k
  real:: CristophCoefficient,EtaFluxCoefficient
  real:: DBXFace,DBYFace,DBZFACE,DB2Face
  real,dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn):: &
       BX,BY,BZ,B2
  
  real,dimension(1:nIFace,1:nJFace,1:nKFace,nBLK):: &
       CCfaceX_FB, CCfaceY_FB, CCfaceZ_FB 
  logical :: DoSetCoef_B(nBLK)=.true.
  real :: Area, InvDxyz
 !
  BX = State_VGB(Bx_,:,:,:,globalBLK)+B0xCell_BLK(:,:,:,globalBLK)
  BY = State_VGB(By_,:,:,:,globalBLK)+B0yCell_BLK(:,:,:,globalBLK)
  BZ = State_VGB(Bz_,:,:,:,globalBLK)+B0zCell_BLK(:,:,:,globalBLK)
  B2 = BX**2+BY**2+BZ**2   
  !\
  ! Compute the resistivity
  !/
  call compute_eta_coefficient(BX,BY,BZ,EtaResist_G)

  !!! Store it for plotting
   Eta_GB(:,:,:,GlobalBLK) = EtaResist_G

  !\
  ! Compute and add the x_resistive_flux to the x-face fluxes 
  !/
  CristophCoefficient = fAx_BLK(globalBLK)/dx_BLK(globalBLK)
  if(UseCovariant.and. DoSetCoef_B(globalBLK))then                !^CFG IF COVARIANT BEGIN
     DoSetCoef_B(globalBLK)=.false.

     do k=1,nK; do j=1,nJ; do i=1,nIFace 
      
        Area = sqrt(max(FaceAreaI_DFB(x_, i, j, k, globalBLK)**2  & 
             +          FaceAreaI_DFB(y_, i, j, k, globalBLK)**2  &
             +          FaceAreaI_DFB(z_, i, j, k, globalBLK)**2, & 
             FaceArea2MinI_B(globalBLK)) )                 
        InvDxyz  = 1.0/sqrt(&           
             ( x_BLK(i  , j, k, globalBLK)       &                         
             - x_BLK(i-1, j, k, globalBLK))**2 + &                         
             ( y_BLK(i  , j, k, globalBLK)       &                         
             - y_BLK(i-1, j, k, globalBLK))**2 + &                         
             ( z_BLK(i  , j, k, globalBLK)       &                         
             - z_BLK(i-1, j, k, globalBLK))**2 )                           
        CCfaceX_FB(i,j,k,globalBLK)= Area*InvDxyz 
     end do; end do; end do

     do k=1,nK; do j=1,nJFace; do i=1,nI
        Area = sqrt(max(FaceAreaJ_DFB(x_, i, j, k, globalBLK)**2  &
             +          FaceAreaJ_DFB(y_, i, j, k, globalBLK)**2  & 
             +          FaceAreaJ_DFB(z_, i, j, k, globalBLK)**2, & 
             FaceArea2MinJ_B(globalBLK)) )         
        InvDxyz  = 1.0/sqrt(&  
             ( x_BLK(i, j  , k, globalBLK)       &                         
             - x_BLK(i, j-1, k, globalBLK))**2 + &                         
             ( y_BLK(i, j  , k, globalBLK)       &                         
             - y_BLK(i, j-1, k, globalBLK))**2 + &                         
             ( z_BLK(i, j  , k, globalBLK)       &                         
             - z_BLK(i, j-1, k, globalBLK))**2 )                           
        CCfaceY_FB(i,j,k,globalBLK)= Area*InvDxyz 
     end do; end do; end do
     do k=1,nKFace; do j=1,nJ; do i=1,nI
        Area = sqrt(max(FaceAreaK_DFB(x_, i, j, k, globalBLK)**2  & 
             +          FaceAreaK_DFB(y_, i, j, k, globalBLK)**2  &  
             +          FaceAreaK_DFB(z_, i, j, k, globalBLK)**2, &    
             FaceArea2MinK_B(globalBLK)) )    
        InvDxyz  = 1.0/sqrt(&                
             ( x_BLK(i, j, k  , globalBLK)       &                         
             - x_BLK(i, j, k-1, globalBLK))**2 + &                         
             ( y_BLK(i, j, k  , globalBLK)       &                         
             - y_BLK(i, j, k-1, globalBLK))**2 + &                         
             ( z_BLK(i, j, k  , globalBLK)       &                         
             - z_BLK(i, j, k-1, globalBLK))**2 )                           
        CCfaceZ_FB(i,j,k,globalBLK)= Area*InvDxyz 
     end do; end do; end do
  end if                              !^CFG IF COVARIANT END
   
  if (.not.DoResChangeOnly) then
     do k=kMinFaceX,kMaxFaceX
        do j=jMinFaceX,jMaxFaceX
           do i=1,nIFace
              if(UseCovariant) CristophCoefficient = CCfaceX_FB(i,j,k,globalBLK) !^CFG IF COVARIANT
              call add_resistive_flux_x
           end do
        end do
     end do
  else if (neiLeast(globalBLK)==+1) then
     i=1
     do k=1,nK
        do j=1,nJ
           if(UseCovariant) CristophCoefficient = CCfaceX_FB(i,j,k,globalBLK) !^CFG IF COVARIANT
           call add_resistive_flux_x
        end do
     end do
  else if ( neiLwest(globalBLK)==+1) then
     i=nIFace
     do k=1,nK
        do j=1,nJ
           if(UseCovariant) CristophCoefficient = CCfaceX_FB(i,j,k,globalBLK) !^CFG IF COVARIANT
           call add_resistive_flux_x
        end do
     end do
  end if
  !\  
  ! Compute and add the y_resistive_flux to the y-face fluxes 
  !/
  CristophCoefficient = fAy_BLK(globalBLK)/dy_BLK(globalBLK)
  if (.not.DoResChangeOnly) then
     do k=kMinFaceY,kMaxFaceY
        do j=1,nJFace
           do i=iMinFaceY,iMaxFaceY
              if(UseCovariant) CristophCoefficient = CCfaceY_FB(i,j,k,globalBLK) !^CFG IF COVARIANT
              call add_resistive_flux_y
           end do
        end do
     end do
  else if(neiLsouth(globalBLK)==+1)then
     j=1
     do k=1,nK
        do i=1,nI
           if(UseCovariant) CristophCoefficient = CCfaceY_FB(i,j,k,globalBLK) !^CFG IF COVARIANT
           call add_resistive_flux_y
        end do
     end do
  else if (neiLnorth(globalBLK)==+1) then
     j=nJFace 
     do k=1,nK
        do i=1,nI
           if(UseCovariant) CristophCoefficient = CCfaceY_FB(i,j,k,globalBLK) !^CFG IF COVARIANT
           call add_resistive_flux_y
        end do
     end do
  end if
  !\
  ! Compute and add the z_resistive_flux to the z-face fluxes  
  !/
  CristophCoefficient = fAz_BLK(globalBLK)/dz_BLK(globalBLK)
  if (.not.DoResChangeOnly) then
     do k=1,nKFace
        do j=jMinFaceZ,jMaxFaceZ
           do i=iMinFaceZ,iMaxFaceZ
              if(UseCovariant) CristophCoefficient = CCfaceZ_FB(i,j,k,globalBLK) !^CFG IF COVARIANT
              call add_resistive_flux_z
           end do
        end do
     end do
  else if (neiLbot(globalBLK)==+1) then
     k=1
     do j=1,nJ
        do i=1,nI
           if(UseCovariant) CristophCoefficient = CCfaceZ_FB(i,j,k,globalBLK) !^CFG IF COVARIANT
           call add_resistive_flux_z
        end do
     end do
  else if ( neiLtop(globalBLK)==+1) then
     k=nKFace            
     do j=1,nJ
        do i=1,nI
           if(UseCovariant) CristophCoefficient = CCfaceZ_FB(i,j,k,globalBLK) !^CFG IF COVARIANT
           call add_resistive_flux_z
        end do
     end do
  end if

contains
  !============================================================================
  subroutine add_resistive_flux_x
    implicit none
    EtaFluxCoefficient = CristophCoefficient* &
         cHalf*(EtaResist_G(i-1,j,k)         + &
         EtaResist_G(i  ,j,k))
    DBXFace = BX(i,j,k)-BX(i-1,j,k)
    DBYFace = BY(i,j,k)-BY(i-1,j,k)
    DBZFace = BZ(i,j,k)-BZ(i-1,j,k)
    DB2Face = B2(i,j,k)-B2(i-1,j,k)
    !\
    ! Update the `resistive flux' in the induction equation
    ! Flux,ij|dis = -eta*DjBi
    ! Ignore the term `eta*DiBj' related with `eta*grad(divB)'
    !/
    ! Flux,xx|dis = -eta*DxBx
    Flux_VX(Bx_,i,j,k) = Flux_VX(Bx_,i,j,k)- &
         EtaFluxCoefficient*DBXFace
    ! Flux,yx|dis = -eta*DxBy
    Flux_VX(By_,i,j,k) = Flux_VX(By_,i,j,k)- &
         EtaFluxCoefficient*DBYFace
    ! Flux,zx|dis = -eta*DxBz
    Flux_VX(Bz_,i,j,k) = Flux_VX(Bz_,i,j,k)- &
         EtaFluxCoefficient*DBZFace
    !\
    ! Update the `resistive flux' in the energy equation
    ! Flux,j|dis = -(eta/mu)*Dj(B^2/2)
    ! Ignore the term `(eta/mu)*(B.grad)Bj' related with 
    ! `eta*B.grad(divB)'
    !/
    ! Flux,x|dis = -(eta/mu)*Dx(B^2/2)
    Flux_VX(Energy_,i,j,k) = Flux_VX(Energy_,i,j,k)-cHalf * &
         EtaFluxCoefficient*DB2Face
    ! Update the CFL condition
    VdtFace_x(i,j,k) = VdtFace_x(i,j,k)+cTwo* &
         EtaFluxCoefficient
  end subroutine add_resistive_flux_x

  !============================================================================

  subroutine add_resistive_flux_y
    implicit none
    EtaFluxCoefficient = CristophCoefficient* &
      cHalf*(EtaResist_G(i,j-1,k)         + &
             EtaResist_G(i,j  ,k))
    DBXFace = BX(i,j,k)-BX(i,j-1,k)
    DBYFace = BY(i,j,k)-BY(i,j-1,k)
    DBZFace = BZ(i,j,k)-BZ(i,j-1,k)
    DB2Face = B2(i,j,k)-B2(i,j-1,k)
    !\
    ! Update the `resistive flux' in the induction equation
    ! Flux,ij|dis = -eta*DjBi
    ! Ignore the term `eta*DiBj' related with `eta*grad(divB)'
    !/
    ! Flux,xy|dis = -eta*DyBx
    Flux_VY(Bx_,i,j,k) = Flux_VY(Bx_,i,j,k)- &
         EtaFluxCoefficient*DBXFace
    ! Flux,yy|dis = -eta*DyBy
    Flux_VY(By_,i,j,k) = Flux_VY(By_,i,j,k)- &
         EtaFluxCoefficient*DBYFace
    ! Flux,zy|dis = -eta*DyBz
    Flux_VY(Bz_,i,j,k) = Flux_VY(Bz_,i,j,k)- &
         EtaFluxCoefficient*DBZFace
    !\
    ! Update the `resistive flux' in the energy equation
    ! Flux,j|dis = -(eta/mu)*Dj(B^2/2)
    ! Ignore the term `(eta/mu)*(B.grad)Bj' related with 
    ! `eta*B.grad(divB)'
    !/
    ! Flux,y|dis = -(eta/mu)*Dy(B^2/2)
    Flux_VY(Energy_,i,j,k) = Flux_VY(Energy_,i,j,k)-cHalf * &
         EtaFluxCoefficient*DB2Face
    ! Update the CFL condition
    VdtFace_y(i,j,k) = VdtFace_y(i,j,k)+cTwo* &
         EtaFluxCoefficient
  end subroutine add_resistive_flux_y

  !============================================================================

  subroutine add_resistive_flux_z
    implicit none
    EtaFluxCoefficient = CristophCoefficient* &
      cHalf*(EtaResist_G(i,j,k-1)         + &
             EtaResist_G(i,j,k  ))
    DBXFace = BX(i,j,k)-BX(i,j,k-1)
    DBYFace = BY(i,j,k)-BY(i,j,k-1)
    DBZFace = BZ(i,j,k)-BZ(i,j,k-1)
    DB2Face = B2(i,j,k)-B2(i,j,k-1)
    !\
    ! Update the `resistive flux' in the induction equation
    ! Flux,ij|dis = -eta*DjBi
    ! Ignore the term `eta*DiBj' related with `eta*grad(divB)'
    !/
    ! Flux,xz|dis = -eta*DzBx
    Flux_VZ(Bx_,i,j,k) = Flux_VZ(Bx_,i,j,k)- &
         EtaFluxCoefficient*DBXFace
    ! Flux,yz|dis = -eta*DzBy
    Flux_VZ(By_,i,j,k) = Flux_VZ(By_,i,j,k)- &
         EtaFluxCoefficient*DBYFace
    ! Flux,zz|dis = -eta*DzBz
    Flux_VZ(Bz_,i,j,k) = Flux_VZ(Bz_,i,j,k)- &
         EtaFluxCoefficient*DBZFace
    !\
    ! Update the `resistive flux' in the energy equation
    ! Flux,j|dis = -(eta/mu)*Dj(B^2/2)
    ! Ignore the term `(eta/mu)*(B.grad)Bj' related with 
    ! `eta*B.grad(divB)'
    !/
    ! Flux,z|dis = -(eta/mu)*Dz(B^2/2)
    Flux_VZ(Energy_,i,j,k) = Flux_VZ(Energy_,i,j,k)-cHalf * &
         EtaFluxCoefficient*DB2Face
    ! Update the CFL condition
    VdtFace_z(i,j,k) = VdtFace_z(i,j,k)+cTwo* &
         EtaFluxCoefficient
  end subroutine add_resistive_flux_z

end subroutine add_resistive_flux

!==============================================================================

subroutine compute_eta_coefficient(BX,BY,BZ,Eta_G)
  use ModSize,     ONLY:nI,nJ,nK,gcn
  use ModProcMH,   ONLY:iProc
  use ModMain,     ONLY:globalBLK,BLKtest,UseSpitzerForm, &
       UseAnomResist,Time_Simulation
  use ModVarIndexes,ONLY:Bx_,By_,Bz_,rho_,P_
  use ModGeometry, ONLY:x_BLK,y_BLK,z_BLK,x1,x2,y1,y2, &
       dx_BLK,dy_BLK,dz_BLK
  use ModAdvance,  ONLY:State_VGB
  use ModConst,    ONLY:cBoltzmann,cProtonMass,cElectronMass, &
       cElectronCharge,cEps,cMu,cLightSpeed
  use ModNumConst, ONLY:cOne,cTwo,cFour,cHalf,cZero,cTiny, &
       cHundred,cHundredth,cPi
  use ModPhysics,  ONLY:TypeResist,Eta0Resist,Eta0AnomResist, &
       EtaAnomMaxResist,jCritResist,&
       No2Si_V, UnitT_, UnitX_, UnitRho_, UnitB_, UnitJ_, UnitTemperature_,&
       Alpha0Resist,yShiftResist,TimeInitRise, &
       TimeConstLev
  use ModResist
  use ModUser, ONLY : user_set_resistivity
  use ModMpi
  implicit none
  
  logical:: UseEtaLocDebug = .false.
  logical:: UseEtaAnomDebug = .false.
  integer:: i,j,k
  real:: CU_x,CU_t
  real:: xxx,yyy,zzz,scr,scr1,scr2
  real:: LogLambdaResist,EtaPerpConstResist
  real:: ElapsedTimeResist,TimeFactorResist
  real:: Eta0Resist_ND,Eta0AnomResist_ND,EtaAnomMaxResist_ND,jCritInv_ND
  real, intent(in), dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn):: &
       BX,BY,BZ
  real, intent(out), dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn):: &
       Eta_G
  !  real, dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn):: &
  !       EtaHallResist
  real, dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn):: &
       Omega_eTau_ei2Resist,EtaAnomLocResist,JmagResist
  real :: current_D(3), Jmag2

  !-------------------------------------------------------
  if(UseUserEta)then
     call user_set_resistivity(GlobalBlk, Eta_G)
     return
  end if

  !\
  ! Introduce some units for dimensionalization:: 
  ! This is not necessary, but makes things more clear!
  !/
  CU_x = No2Si_V(UnitX_) ! in [m]
  CU_t = No2Si_V(UnitT_) ! in [s]

  !\
  ! Dimensionalize:: Eta* is in units of [m^2/s]
  !/
  Eta0Resist_ND       = Eta0Resist*CU_t/CU_x**2

  if(UseAnomResist)then
     Eta0AnomResist_ND   = Eta0AnomResist*CU_t/CU_x**2
     EtaAnomMaxResist_ND = EtaAnomMaxResist*CU_t/CU_x**2
     jCritInv_ND         = cOne/(jCritResist/No2Si_V(UnitJ_))
  end if
  !\
  ! Select type of resistivity
  !/
  if (UseSpitzerForm) then
     !\
     ! Compute Spitzer-type, classical resistivity::
     ! Eta_G = EtaPerpConstResist*LogLambdaResist/Te^1.5, 
     ! where EtaPerpConstResist = 5.237943200000000000E+07
     !/ 
     ! Compute EtaPerpConstResist::
     EtaPerpConstResist = cHalf*sqrt(cElectronMass)    * &
          (cLightSpeed*cElectronCharge)**2/((cOne+cTwo)* &
          cEps*(cTwo*cPi*cBoltzmann)**(cOne+cHalf))
     !\
     ! Assume ln(Lambda) = 20 (solar corona) -- coded
     !                  or 30 (magnetosphere)
     ! _________________or 10 (lab. plasma)
     !/
     LogLambdaResist = cFour*(cOne+cFour) ! Coronal value
     !\
     ! Take into account the Coulomb logarithm::
     ! EtaPerpConstResist = EtaPerpConstResist*LogLambdaResist
     !/
     EtaPerpConstResist = EtaPerpConstResist*LogLambdaResist
     !\
     ! Finally, compute Eta_G::
     ! Eta_G = EtaPerpConstResist/Te^1.5
     !/
     Eta_G(:,:,:) = EtaPerpConstResist/(No2Si_V(UnitTemperature_)* &
          State_VGB(P_,:,:,:,globalBLK)/State_VGB(rho_,:,:,:,globalBLK))**(cOne+cHalf)
     !\
     ! Take into account the dependance of the B field::
     ! Eta_G = Eta_G*(1+Omega_eTau_ei2Resist), where
     ! Omega_eTau_ei2Resist = [B*mp/(rho*e*Eta_G)]^2
     !/
     Omega_eTau_ei2Resist(:,:,:) = (cProtonMass/cElectronCharge)**2* &
          (No2Si_V(UnitB_)**2*(BX**2+BY**2+BZ**2))   / &
          (No2Si_V(UnitRho_)*State_VGB(rho_,:,:,:,globalBLK)*Eta_G)**2
     ! Eta_G = Eta_G*(1+Omega_eTau_ei2Resist)
     Eta_G(:,:,:) = Eta_G(:,:,:)* &
          (cOne+Omega_eTau_ei2Resist)
     ! Dimensionalize Eta_G::
     Eta_G = Eta_G*CU_t/CU_x**2
     !\
     ! Compute the Hall resistivity, for whatever needed::
     ! EtaHallResist = Eta_G/sqrt( &
     ! max(Omega_eTau_ei2Resist,cTiny**2))
     !/
     !     EtaHallResist(:,:,:) = Eta_G(:,:,:)/ &
     !          sqrt(max(Omega_eTau_ei2Resist(:,:,:),cTiny**2))
     !\
     ! Add anomalous resistivity, if UseAnomResist == .true.
     !/
     if (UseAnomResist) &
          call add_anomalous_resistivity
  else
     select case(TypeResist)
     case('Constant','constant')
        !\
        ! Set uniform resistivity everywhere
        !/
        Eta_G(:,:,:) = Eta0Resist_ND
        !\
        ! Add anomalous resistivity, if UseAnomResist == .true.
        !/
        if (UseAnomResist) &
             call add_anomalous_resistivity
        
     case('Localized','localized')
        !\
        ! Compute localized in space, time-dependant resistivity
        !/
        if (DoInitEtaLocResist_B(globalBLK)) then
           DoInitEtaLocResist_B(globalBLK) = .false.
           do i=1-gcn,nI+gcn; do j=1-gcn,nJ+gcn; do k=1-gcn,nK+gcn
              xxx  = x_BLK(i,j,k,globalBLK)
              yyy  = y_BLK(i,j,k,globalBLK)+yShiftResist*y2
              scr1 = Alpha0Resist*xxx**2
              scr2 = Alpha0Resist*yyy**2
              if (scr1.le.(cOne+cHalf)*cHundred) then
                 scr1 = exp(-scr1)
              else
                 scr1 = cZero
              end if
              if (scr2.le.(cOne+cHalf)*cHundred) then
                 scr2 = exp(-scr2)
              else
                 scr2 = cZero
              end if
              EtaLocResist_GB(i,j,k,globalBLK) = Eta0Resist_ND*scr1*scr2
           end do; end do; end do
           if (iProc==0.and.globalBLK==BLKtest) then
              write(6,*) ''
              write(6,*) '>>>>>>>>>>>>>>>>> Localized Resistivity <<<<<<<<<<<<<<<<<'
              write(6,*) '>>>>>>>>>>>>>>>>>       Parameters      <<<<<<<<<<<<<<<<<'
              write(6,*) ''
              write(6,*) 'Eta0    = ',Eta0Resist
              write(6,*) 'Eta0_ND = ',Eta0Resist_ND
              write(6,*) ''
              write(6,*) 'TimeInitRise,TimeConstLev:: ', &
                          TimeInitRise,TimeConstLev
              write(6,*) 'Alpha0Resist,yShiftResist:: ', &
                          Alpha0Resist,yShiftResist
              write(6,*) 'minval(EtaLoc),maxval(EtaLoc):: ', &
                          minval(EtaLocResist_GB(:,:,:,:)),  &
                          maxval(EtaLocResist_GB(:,:,:,:))
              write(6,*) ''
              write(6,*) '>>>>>>>>>>>>>>>>>                       <<<<<<<<<<<<<<<<<'
              write(6,*) ''
           end if
        end if
        ElapsedTimeResist = Time_Simulation/No2Si_V(UnitT_)
        if (ElapsedTimeResist.le.TimeInitRise) then
           TimeFactorResist = ElapsedTimeResist/TimeInitRise
        else if (ElapsedTimeResist.le.TimeConstLev) then
           TimeFactorResist = cOne
        else
           TimeFactorResist = cZero
        endif
        if (UseEtaLocDebug.and.iProc==0.and.globalBLK==BLKtest) then
           write(6,*) 'Elapsed Time T = ',ElapsedTimeResist
           write(6,*) 'minval(TimeFactor*EtaLoc),maxval(TimeFactor*EtaLoc):: ', &
                       minval(TimeFactorResist*EtaLocResist_GB(:,:,:,:)),       &
                       maxval(TimeFactorResist*EtaLocResist_GB(:,:,:,:))
        end if
        Eta_G(:,:,:) = TimeFactorResist*EtaLocResist_GB(:,:,:,globalBLK)
        !\
        ! Add anomalous resistivity, if UseAnomResist == .true.
        !/
        if (UseAnomResist) &
             call add_anomalous_resistivity
        
     case('Hall_localized')
        !\
        ! Compute localized in space, time-dependant resistivity
        !/
        if (DoInitEtaLocResist_B(globalBLK)) then
           DoInitEtaLocResist_B(globalBLK) = .false.
           do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
              call get_current(i,j,k,GlobalBlk,current_D)
              Jmag2 = sum(current_D**2)
              if(Jmag2 > 4.0 )then
                 xxx  = x_BLK(i,j,k,globalBLK)
                 zzz  = z_BLK(i,j,k,globalBLK)
                 scr  = sqrt(xxx**2/2.0+zzz**2)
                 if (scr < 150.0) then
                    scr = 1.0/ cosh(scr)
                 else
                    scr = 0.0
                 end if
                 EtaLocResist_GB(i,j,k,globalBLK) = Eta0Resist_ND*(1+Alpha0Resist*scr)
              else
                 EtaLocResist_GB(i,j,k,globalBLK) = Eta0Resist_ND
              end if
           end do; end do; end do
           if (iProc==0.and.globalBLK==1) then
              write(*,*) ''
              write(*,*) '>>>>>>>>>>>>>>>>> Localized Resistivity <<<<<<<<<<<<<<<<<'
              write(*,*) '>>>>>>>>>>>>>>>>>       Parameters      <<<<<<<<<<<<<<<<<'
              write(*,*) ''
              write(*,*) 'Eta0    = ',Eta0Resist
              write(*,*) 'Eta0_ND = ',Eta0Resist_ND
              write(*,*) ''
              write(*,*) 'TimeInitRise,TimeConstLev:: ', &
                   TimeInitRise,TimeConstLev
              write(*,*) 'Alpha0Resist,yShiftResist:: ', &
                   Alpha0Resist,yShiftResist
              write(*,*) 'minval(EtaLoc),maxval(EtaLoc):: ', &
                   minval(EtaLocResist_GB(:,:,:,:)),  &
                   maxval(EtaLocResist_GB(:,:,:,:))
              write(*,*) ''
              write(*,*) '>>>>>>>>>>>>>>>>>                       <<<<<<<<<<<<<<<<<'
              write(*,*) ''
           end if
        end if
        TimeFactorResist = cOne
        Eta_G(:,:,:) = EtaLocResist_GB(:,:,:,globalBLK)
        !\
        ! Add anomalous resistivity, if UseAnomResist == .true.
        !/
        if (UseAnomResist) call add_anomalous_resistivity
        
     case default
        call stop_mpi('Unknown TypeResist::'//TypeResist)
     end select
  end if

contains
  !============================================================================
  subroutine add_anomalous_resistivity
    implicit none
    real :: current_D(3)

    !\
    ! Compute the magnitude of the current density |J| and eta_anomalous
    !/
    do i=0,nI+1; do j=0,nJ+1; do k=0,nK+1

       call get_current(i,j,k,GlobalBlk,current_D)
       JmagResist(i,j,k) = sqrt(sum(current_D**2))
       !\
       ! Compute the anomalous resistivity:: 
       ! EtaAnomLocResist(i,j,k) = Eta0AnomResist_ND*(|J|/Jcrit-1) 
       !/
       EtaAnomLocResist(i,j,k) = &
            min(EtaAnomMaxResist_ND,max(cZero, &
            Eta0AnomResist_ND*(JcritInv_ND*JmagResist(i,j,k)-cOne)))

       !\
       ! Add the anomalous to the background resistivity
       !/
       Eta_G(i,j,k) = Eta_G(i,j,k)+EtaAnomLocResist(i,j,k)

    end do; end do; end do
    if (UseEtaAnomDebug.and.iProc==0.and.globalBLK==BLKtest) then
       write(*,*) ''
       write(*,*) 'min(EtaAnomLocResist) ,max(EtaAnomLocResist) :: ', &
                   minval(EtaAnomLocResist), maxval(EtaAnomLocResist)
       write(*,*) 'min(jMagResist) ,max(jMagResist) :: ', &
                   minval(jMagResist), maxval(jMagResist)
       write(*,*) ''
    end if
  end subroutine add_anomalous_resistivity
  
end subroutine compute_eta_coefficient
