!^CFG FILE DISSFLUX
!^CFG COPYRIGHT UM
Module ModResist
  use ModSize,     ONLY:nI,nJ,nK,gcn,nBLK
  implicit none
  
  logical :: DoInitEtaLocResist_B(nBLK)
  data DoInitEtaLocResist_B /nBLK*.true./
  logical, dimension (1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn,nBLK) :: &
       IsNotIgnitedYetResist_GB
  integer, parameter :: MaxGhostPE=(nI+2*gcn)*(nJ+2*gcn)*(nK+2*gcn)*nBLK
  data IsNotIgnitedYetResist_GB /MaxGhostPE*.true./
  real, dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn,nBLK):: &
       EtaLocResist_GB,JcritInvResist_GB
  
end module ModResist

subroutine add_resistive_flux(DoResChangeOnly)
  use ModSize,     ONLY:nI,nJ,nK,gcn,nBLK
  use ModMain,     ONLY:nIFace,nJFace,nKFace,&
       iMinFaceY,iMaxFaceY,iMinFaceZ,iMaxFaceZ, &
       jMinFaceX,jMaxFaceX,jMinFaceZ,jMaxFaceZ, &
       kMinFaceX,kMaxFaceX,kMinFaceY,kMaxFaceY, &
       globalBLK,BLKtest,UseSpitzerForm
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
  use ModMpi
  implicit none
  
  logical, intent(in):: DoResChangeOnly
  integer:: i,j,k
  real:: CristophCoefficient,EtaFluxCoefficient
  real:: DBXFace,DBYFace,DBZFace,DB2Face
  real,dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn):: &
       BX,BY,BZ,B2,EtaPerpResist
  !
  BX(:,:,:) = State_VGB(Bx_,:,:,:,globalBLK)+B0xCell_BLK(:,:,:,globalBLK)
  BY(:,:,:) = State_VGB(By_,:,:,:,globalBLK)+B0yCell_BLK(:,:,:,globalBLK)
  BZ(:,:,:) = State_VGB(Bz_,:,:,:,globalBLK)+B0zCell_BLK(:,:,:,globalBLK)
  B2(:,:,:) = BX**2+BY**2+BZ**2   
  !\
  ! Compute the resistivity
  !/
  call compute_eta_coefficient(BX,BY,BZ,EtaPerpResist)
  !\
  ! Compute and add the x_resistive_flux to the x-face fluxes 
  !/
  CristophCoefficient = fAx_BLK(globalBLK)/dx_BLK(globalBLK)
  if (.not.DoResChangeOnly) then
     do k=kMinFaceX,kMaxFaceX
        do j=jMinFaceX,jMaxFaceX
           do i=1,nIFace
              call add_resistive_flux_x
           end do
        end do
     end do
  else if (neiLeast(globalBLK)==+1) then
     i=1
     do k=1,nK
        do j=1,nJ
           call add_resistive_flux_x
        end do
     end do
  else if ( neiLwest(globalBLK)==+1) then
     i=nIFace
     do k=1,nK
        do j=1,nJ
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
              call add_resistive_flux_y
           end do
        end do
     end do
  else if(neiLsouth(globalBLK)==+1)then
     j=1
     do k=1,nK
        do i=1,nI
           call add_resistive_flux_y
        end do
     end do
  else if (neiLnorth(globalBLK)==+1) then
     j=nJFace 
     do k=1,nK
        do i=1,nI
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
              call add_resistive_flux_z
           end do
        end do
     end do
  else if (neiLbot(globalBLK)==+1) then
     k=1
     do j=1,nJ
        do i=1,nI
           call add_resistive_flux_z
        end do
     end do
  else if ( neiLtop(globalBLK)==+1) then
     k=nKFace            
     do j=1,nJ
        do i=1,nI
           call add_resistive_flux_z
        end do
     end do
  end if
Contains
  
  subroutine add_resistive_flux_x
    implicit none
    EtaFluxCoefficient = CristophCoefficient* &
      cHalf*(EtaPerpResist(i-1,j,k)         + &
             EtaPerpResist(i  ,j,k))
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
  
  subroutine add_resistive_flux_y
    implicit none
    EtaFluxCoefficient = CristophCoefficient* &
      cHalf*(EtaPerpResist(i,j-1,k)         + &
             EtaPerpResist(i,j  ,k))
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

  subroutine add_resistive_flux_z
    implicit none
    EtaFluxCoefficient = CristophCoefficient* &
      cHalf*(EtaPerpResist(i,j,k-1)         + &
             EtaPerpResist(i,j,k  ))
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

subroutine compute_eta_coefficient(BX,BY,BZ,EtaPerpResist)
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
       EtaAnomMaxResist,unitSI_t,unitSI_x,unitSI_rho,unitSI_B, &
       unitSI_temperature,Alpha0Resist,yShiftResist,TimeInitRise, &
       TimeConstLev,ThresholdFactorResist
  use ModResist
  use ModMpi
  implicit none
  
  logical:: UseEtaLocDebug = .false.
  logical:: UseEtaAnomDebug = .false.
  integer:: i,j,k
  real:: CU_x,CU_t
  real:: xxx,yyy,scr1,scr2
  real:: LogLambdaResist,EtaPerpConstResist
  real:: ElapsedTimeResist,TimeFactorResist
  real:: UdriftResist,UthresholdResist
  real:: Eta0Resist_ND,Eta0AnomResist_ND,EtaAnomMaxResist_ND
  real, intent(in), dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn):: &
       BX,BY,BZ
  real, intent(out), dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn):: &
       EtaPerpResist
  !  real, dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn):: &
  !       EtaHallResist
  real, dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn):: &
       Omega_eTau_ei2Resist,EtaAnomLocResist,JmagResist
  !\
  ! Introduce some units for dimensionalization:: 
  ! This is not necessary, but makes things more clear!
  !/
  CU_x = unitSI_x ! in [m]
  CU_t = unitSI_t ! in [s]
  !\
  ! Dimensionalize:: Eta* is in units of [m^2/s]
  !/
  Eta0Resist_ND       = Eta0Resist*CU_t/CU_x**2
  Eta0AnomResist_ND   = Eta0AnomResist*CU_t/CU_x**2
  EtaAnomMaxResist_ND = EtaAnomMaxResist*CU_t/CU_x**2
  !\
  ! Select type of resistivity
  !/
  if (UseSpitzerForm) then
     !\
     ! Compute Spitzer-type, classical resistivity::
     ! EtaPerpResist = EtaPerpConstResist*LogLambdaResist/Te^1.5, 
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
     ! Finally, compute EtaPerpResist::
     ! EtaPerpResist = EtaPerpConstResist/Te^1.5
     !/
     EtaPerpResist(:,:,:) = EtaPerpConstResist/(unitSI_temperature* &
          State_VGB(P_,:,:,:,globalBLK)/State_VGB(rho_,:,:,:,globalBLK))**(cOne+cHalf)
     !\
     ! Take into account the dependance of the B field::
     ! EtaPerpResist = EtaPerpResist*(1+Omega_eTau_ei2Resist), where
     ! Omega_eTau_ei2Resist = [B*mp/(rho*e*EtaPerpResist)]^2
     !/
     Omega_eTau_ei2Resist(:,:,:) = (cProtonMass/cElectronCharge)**2* &
          (unitSI_B**2*(BX(:,:,:)**2+BY(:,:,:)**2+BZ(:,:,:)**2))   / &
          (unitSI_rho*State_VGB(rho_,:,:,:,globalBLK)*EtaPerpResist(:,:,:))**2
     ! EtaPerpResist = EtaPerpResist*(1+Omega_eTau_ei2Resist)
     EtaPerpResist(:,:,:) = EtaPerpResist(:,:,:)* &
          (cOne+Omega_eTau_ei2Resist(:,:,:))
     ! Dimensionalize EtaPerpResist::
     EtaPerpResist(:,:,:) = EtaPerpResist(:,:,:)*CU_t/CU_x**2
     !\
     ! Compute the Hall resistivity, for whatever needed::
     ! EtaHallResist = EtaPerpResist/sqrt( &
     ! max(Omega_eTau_ei2Resist,cTiny**2))
     !/
     !     EtaHallResist(:,:,:) = EtaPerpResist(:,:,:)/ &
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
        EtaPerpResist(:,:,:) = Eta0Resist_ND
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
        ElapsedTimeResist = Time_Simulation/unitSI_t
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
        EtaPerpResist(:,:,:) = TimeFactorResist*EtaLocResist_GB(:,:,:,globalBLK)
        !\
        ! Add anomalous resistivity, if UseAnomResist == .true.
        !/
        if (UseAnomResist) &
             call add_anomalous_resistivity
        
     case default
        call stop_mpi('Unknown TypeResist::'//TypeResist)
     end select
  end if
Contains
  
  subroutine add_anomalous_resistivity
    implicit none
    !\
    ! Compute the magnitude of the current density:: |J|
    !/
    JmagResist(0:nI+1,0:nJ+1,0:nK+1) = sqrt((cHalf*(                                &
         (BZ(0:nI+1,1:nJ+2,0:nK+1)-BZ(0:nI+1,-1:nJ,0:nK+1))/dy_BLK(globalBLK)     - &
         (BY(0:nI+1,0:nJ+1,1:nK+2)-BY(0:nI+1,0:nJ+1,-1:nK))/dz_BLK(globalBLK)))**2+ &
         (cHalf*(                                                                   &
         (BX(0:nI+1,0:nJ+1,1:nK+2)-BX(0:nI+1,0:nJ+1,-1:nK))/dz_BLK(globalBLK)     - &
         (BZ(1:nI+2,0:nJ+1,0:nK+1)-BZ(-1:nI,0:nJ+1,0:nK+1))/dx_BLK(globalBLK)))**2+ &
         (cHalf*(                                                                   &
         (BY(1:nI+2,0:nJ+1,0:nK+1)-BY(-1:nI,0:nJ+1,0:nK+1))/dx_BLK(globalBLK)     - &
         (BX(0:nI+1,1:nJ+2,0:nK+1)-BX(0:nI+1,-1:nJ,0:nK+1))/dy_BLK(globalBLK)))**2)
    !\
    ! Loop over to see where the anomalous resistivity ignites::
    !/
    do i=1-gcn,nI+gcn; do j=1-gcn,nJ+gcn; do k=1-gcn,nK+gcn
       !\
       ! Compute the threshold velocity:: 
       ! UthresholdResist = sqrt(cProtonMass/cElectronMass*ThresholdFactorResist*Te)
       !/
       UthresholdResist = sqrt(cProtonMass/cElectronMass*ThresholdFactorResist* &
            State_VGB(P_,i,j,k,globalBLK)/State_VGB(rho_,i,j,k,globalBLK))
       !\
       ! Compute the drift velocity:: UdriftResist = |J|/rho
       !/
       UdriftResist = JmagResist(i,j,k)/State_VGB(rho_,i,j,k,globalBLK) 
       if (UdriftResist.ge.UthresholdResist.and. &
          IsNotIgnitedYetResist_GB(i,j,k,globalBLK)) then
          IsNotIgnitedYetResist_GB(i,j,k,globalBLK) = .false.
          JcritInvResist_GB(i,j,k,globalBLK) = cOne/JmagResist(i,j,k)
       end if
       if (UdriftResist.lt.UthresholdResist) then
          IsNotIgnitedYetResist_GB(i,j,k,globalBLK) = .true.
          JcritInvResist_GB(i,j,k,globalBLK) = cZero
       end if
       !\
       ! Compute the anomalous resistivity:: 
       ! EtaAnomLocResist(i,j,k) = Eta0AnomResist_ND*(|J|/Jcrit-1) 
       !/
       EtaAnomLocResist(i,j,k) = Eta0AnomResist_ND* &
            (JcritInvResist_GB(i,j,k,globalBLK)*JmagResist(i,j,k)-cOne)
       !\
       ! If EtaAnomLocResist(i,j,k)<0, then EtaAnomLocResist(i,j,k)=0
       !/
       if (EtaAnomLocResist(i,j,k).lt.cZero) &
           EtaAnomLocResist(i,j,k) = cZero
       !\
       ! If EtaAnomLocResist(i,j,k)>EtaAnomMaxResist_ND, then 
       ! EtaAnomLocResist(i,j,k)=EtaAnomMaxResist_ND 
       !/
       if (EtaAnomLocResist(i,j,k).ge.EtaAnomMaxResist_ND) &
           EtaAnomLocResist(i,j,k) = EtaAnomMaxResist_ND
       !\
       ! Add the anomalous to the background resistivity
       !/
       EtaPerpResist(i,j,k) = EtaPerpResist(i,j,k)+EtaAnomLocResist(i,j,k)
    end do; end do; end do
    if (UseEtaAnomDebug.and.iProc==0.and.globalBLK==BLKtest) then
       write(6,*) ''
       write(6,*) 'min(EtaAnomLocResist) ,max(EtaAnomLocResist) :: ', &
                   minval(EtaAnomLocResist(:,:,:)),    &
                   maxval(EtaAnomLocResist(:,:,:))
       write(6,*) 'min(JcritInvResist_GB),max(JcritInvResist_GB):: ', &
                   minval(JcritInvResist_GB(:,:,:,:)), &
                   maxval(JcritInvResist_GB(:,:,:,:))
       write(6,*) ''
    end if
  end subroutine add_anomalous_resistivity
  
end subroutine compute_eta_coefficient
