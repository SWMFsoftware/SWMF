!^CFG COPYRIGHT UM
subroutine calc_sources

  use ModProcMH
  use ModMain
  use ModVarIndexes
  use ModGeometry, ONLY : dx_BLK, dy_BLK, dz_BLK, R_BLK,&
       body_BLK, Rmin_BLK, vInv_CB
  use ModGeometry, ONLY : R2_BLK                        !^CFG IF SECONDBODY
  use ModAdvance
  use ModParallel, ONLY : NOBLK, neiLEV, &
       neiLtop,neiLbot,neiLeast,neiLwest,neiLnorth,neiLsouth
  use ModPhysics
  use ModNumConst
  use ModResistivity,ONLY : UseResistivity, Eta_GB   !^CFG IF DISSFLUX
  use ModUser,     ONLY : user_calc_sources
  use ModCoordTransform
  use ModHallResist, ONLY: &
       UseHallResist, HallHyperFactor, calc_hyper_resistivity 
  use ModGrayDiffusion, ONLY: UseGrayDiffusion, calc_source_gray_diffusion
  use ModMultiFluid
  use ModPointImplicit, ONLY: UsePointImplicit, UsePointImplicit_B
  use ModMultiIon, ONLY: multi_ion_sources

  use ModCovariant, ONLY: UseCovariant 

  implicit none

  integer :: i, j, k, iDim

  real :: Coef

  ! Variable for div B diffusion
  real :: Dr

  real,dimension(3)::CurlB0CrossB_D

  ! Variables needed for Boris source terms also used for div(u)
  real :: FullB_DC(nDim,nI,nJ,nK),FullBx, FullBy, FullBz, Ux, Uy, Uz, RhoInv
  real :: E_D(3), DivE

  ! Variables needed for Joule heating
  real :: Current_D(3)

  integer:: iBlock

  logical :: DoTest, DoTestMe
  character(len=*), parameter :: NameSub = 'calc_sources'
  !---------------------------------------------------------------------------
  iBlock = GlobalBlk

  if(iProc==PROCtest .and. iBlock==BLKtest)then
     call set_oktest(NameSub, DoTest, DoTestMe)
  else
     DoTest=.false.; DoTestMe=.false.
  end if

  Source_VC   = 0.0

  ! Calculate source terms for ion pressure
  if(UseNonconservative)then
     do iFluid = 1, nFluid
        call select_fluid
        ! Adiabatic heating: -(g-1)*P*Div(U)
        do k=1,nK; do j=1,nJ; do i=1,nI
           Source_VC(iP,i,j,k) = -(g-1)*State_VGB(iP,i,j,k,iBlock)*&
                vInv_CB(i,j,k,iBlock)*&
                (uDotArea_XI(i+1,j,k,iFluid) - uDotArea_XI(i,j,k,iFluid) &
                +uDotArea_YI(i,j+1,k,iFluid) - uDotArea_YI(i,j,k,iFluid) &
                +uDotArea_ZI(i,j,k+1,iFluid) - uDotArea_ZI(i,j,k,iFluid))
        end do; end do; end do

        if(DoTestMe.and.VarTest==iP)call write_source('After p div U')
     end do

     ! Joule heating: dP/dt += (gamma-1)*eta*j**2
     if(UseResistFlux .or. UseResistivity)then  !^CFG IF DISSFLUX BEGIN
        do k=1,nK; do j=1,nJ; do i=1,nI           
           call get_current(i,j,k,iBlock,Current_D)
           Source_VC(P_,i,j,k) = Source_VC(P_,i,j,k) + &
                (g-1) * Eta_GB(i,j,k,iBlock) * sum(Current_D**2)
        end do; end do; end do

        if(DoTestMe.and.VarTest==P_)call write_source('After eta j')

     end if                 !^CFG END DISSFLUX
  end if

  if(UseElectronPressure)then
     ! Adiabatic heating for electron pressure: -(g-1)*Pe*Div(U)
     do k=1,nK; do j=1,nJ; do i=1,nI
        Source_VC(Pe_,i,j,k) = -(g-1)*State_VGB(Pe_,i,j,k,iBlock)*&
             vInv_CB(i,j,k,iBlock)*&
             (uDotArea_XI(i+1,j,k,eFluid_) - uDotArea_XI(i,j,k,eFluid_) &
             +uDotArea_YI(i,j+1,k,eFluid_) - uDotArea_YI(i,j,k,eFluid_) &
             +uDotArea_ZI(i,j,k+1,eFluid_) - uDotArea_ZI(i,j,k,eFluid_))
     end do; end do; end do

     if(DoTestMe.and.VarTest==Pe_)call write_source('After Pe div Ue')

  end if


  if(UseDivbSource)then
     if(UseCovariant)then   
        call calc_divb_source_covar
     else                   
        call calc_divb_source
     end if                 

     if(DoTestMe)write(*,*)'divb=',DivB1_GB(iTest,jTest,kTest,BlkTest)
     if(DoTestMe.and.VarTest>=RhoUx_.and.VarTest<=RhoUz_)&
          call write_source('After B0B1 source')

     ! Add contributions to other source terms
     do k=1,nK; do j=1,nJ; do i=1,nI
        RhoInv=cOne/State_VGB(rho_,i,j,k,iBlock)
        Source_VC(Bx_:Bz_,i,j,k)    = Source_VC(Bx_:Bz_,i,j,k) &
             -DivB1_GB(i,j,k,iBlock)* &
             State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock)*RhoInv

        if(.not. IsMhd) CYCLE

        Source_VC(rhoUx_:rhoUz_,i,j,k) = Source_VC(rhoUx_:rhoUz_,i,j,k) -&
             DivB1_GB(i,j,k,iBlock)* &
             State_VGB(Bx_:Bz_,i,j,k,iBlock)
        Source_VC(Energy_,i,j,k)     = Source_VC(Energy_,i,j,k) &
             -DivB1_GB(i,j,k,iBlock)* &
             sum(State_VGB(Bx_:Bz_,i,j,k,iBlock)*&
             State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock))*RhoInv
     end do;end do;end do

     if (UseB0Source) then
        do k=1,nK; do j=1,nJ; do i=1,nI
           Source_VC(rhoUx_:rhoUz_,i,j,k)=Source_VC(rhoUx_:rhoUz_,i,j,k) - &
                State_VGB(Bx_:Bz_,i,j,k,iBlock)*DivB0_CB(i,j,k,iBlock)+&
                cross_product(&
                State_VGB(Bx_:Bz_,i,j,k,iBlock),&
                CurlB0_DCB(:,i,j,k,iBlock))
        end do; end do; end do
     end if
  else
     if(UseB)call calc_divb(iBlock)
  end if

  if(UseCurlB0)then
     do k=1,nK; do j=1,nJ; do i=1,nI
        if(R_BLK(i,j,k,iBlock)<rCurrentFreeB0)CYCLE
        CurlB0CrossB_D = cross_product(&
             CurlB0_DCB(:,i,j,k,iBlock),&
             State_VGB(Bx_:Bz_,i,j,k,iBlock)+&
             B0_DGB(:,i,j,k,iBlock))
        Source_VC(rhoUx_:rhoUz_,i,j,k)= Source_VC(rhoUx_:rhoUz_,i,j,k) +&
             CurlB0CrossB_D
        Source_VC(Energy_,i,j,k)     = Source_VC(Energy_,i,j,k)        +&
             sum(CurlB0CrossB_D*State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock))&
             /State_VGB(rho_,i,j,k,iBlock)
     end do;end do;end do

     if(DoTestMe .and. &
          (VarTest==Energy_.or.VarTest>=RhoUx_.and.VarTest<=RhoUz_))&
          call write_source('After curl B0')
  end if

  if(boris_correction &                             !^CFG IF BORISCORR BEGIN
       .and. boris_cLIGHT_factor < 0.9999 & 
       .and. index(test_string,'nodivE')<1) then

     coef= (boris_cLIGHT_factor**2 - 1.0)*inv_c2LIGHT
     FullB_DC = State_VGB(Bx_:Bz_,1:nI,1:nJ,1:nK,iBlock)
     if(UseB0)FullB_DC = FullB_DC + B0_DGB(:,1:nI,1:nJ,1:nK,iBlock) 
     do k=1,nK; do j=1,nJ; do i=1,nI
        E_D = cross_product(FullB_DC(:,i,j,k),&
               State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock))/&
               State_VGB(rho_,i,j,k,iBlock)
        ! Calculate divergence of electric field 
        DivE = vInv_CB(i,j,k,iBlock)*&
             (EDotFA_X(i+1,j,k)-EDotFA_X(i,j,k)+&
             EDotFA_Y(i,j+1,k) -EDotFA_Y(i,j,k)+&
             EDotFA_Z(i,j,k+1) -EDotFA_Z(i,j,k))

        Source_VC(rhoUx_:rhoUz_,i,j,k) = Source_VC(rhoUx_:rhoUz_,i,j,k) &
             + Coef*DivE*E_D 

        if(DoTestMe.and.VarTest>=RhoUx_.and.VarTest<=RhoUz_) &
             call write_source('After E div E')

     end do; end do; end do
  end if                                                 !^CFG END BORISCORR

  ! These source terms apply to all the fluids
  do iFluid = 1, nFluid
     call select_fluid
     ! Add gravity and/or centrifugal force
     if(UseGravity.or.UseRotatingFrame) then
        Source_VC(iRhoUx,:,:,:) = Source_VC(iRhoUx,:,:,:) + &
             State_VGB(iRho,1:nI,1:nJ,1:nK,iBlock)* &
             fbody_x_BLK(:,:,:,iBlock)
        Source_VC(iRhoUy,:,:,:) = Source_VC(iRhoUy,:,:,:) + &
             State_VGB(iRho,1:nI,1:nJ,1:nK,iBlock)* &
             fbody_y_BLK(:,:,:,iBlock)
        Source_VC(iRhoUz,:,:,:) = Source_VC(iRhoUz,:,:,:) + &
             State_VGB(iRho,1:nI,1:nJ,1:nK,iBlock)* &
             fbody_z_BLK(:,:,:,iBlock)
        Source_VC(Energy_,:,:,:) = Source_VC(Energy_,:,:,:) + &
             (State_VGB(iRhoUx,1:nI,1:nJ,1:nK,iBlock)* &
             fbody_x_BLK(:,:,:,iBlock) + & 
             State_VGB(iRhoUy,1:nI,1:nJ,1:nK,iBlock)* &
             fbody_y_BLK(:,:,:,iBlock) + &
             State_VGB(iRhoUz,1:nI,1:nJ,1:nK,iBlock)* &
             fbody_z_BLK(:,:,:,iBlock)) 

        if(DoTestMe.and. &
             (VarTest==Energy_ .or. VarTest>=iRhoUx .and. VarTest<=iRhoUz)) &
             call write_source('After gravity')
     end if

     ! Add Coriolis forces
     if(UseRotatingFrame)then
        select case(TypeCoordSystem)
        case('HGC','HGR')
           ! This is a special case since Omega is parallel with the Z axis
           do k=1,nK; do j=1,nJ; do i=1,nI
              Source_VC(iRhoUx,i,j,k) = Source_VC(iRhoUx,i,j,k) + &
                   2*OmegaBody*State_VGB(iRhoUy,i,j,k,iBlock)
              Source_VC(iRhoUy,i,j,k) = Source_VC(iRhoUy,i,j,k) - &
                   2*OmegaBody*State_VGB(iRhoUx,i,j,k,iBlock)
           end do; end do; end do
        case default
           call stop_mpi(NameSub // &
                ' Coriolis force is not implemented for'// &
                ' TypeCoordSystem='//TypeCoordSystem)
        end select
        if(DoTestMe.and.VarTest>=iRhoUx .and. VarTest<=iRhoUy) &
             call write_source('After Coriolis')
     end if
  end do

  ! Explicit evaluation of multi-ion source for development purposes only
  if(UseMultiIon .and. &
       .not. (UsePointImplicit .and. UsePointImplicit_B(iBlock)) ) then
     call multi_ion_sources
     if(DoTestMe) call write_source('After MultiIon sources')
  end if

  if(UseHallResist .and. HallHyperFactor > 0.0) then
     call calc_hyper_resistivity(iBlock)
     if(DoTestMe) call write_source('After HyperResist')
  end if

  if(UseGrayDiffusion) call calc_source_gray_diffusion(iBlock)

  if(UseUserSource)then
     call user_calc_sources
     if(DoTestMe) call write_source('After user sources')
  end if

  if(DoTestMe) call write_source('final')

contains
  !===========================================================================
  subroutine calc_divb_source

    ! Variables needed for div B source terms
    real:: DxInvHalf, DyInvHalf, DzInvHalf, DivBInternal_C(1:nI,1:nJ,1:nK)
    real:: dB1nEast, dB1nWest, dB1nSouth, dB1nNorth, dB1nTop, dB1nBot

    DxInvHalf = 0.5/Dx_BLK(iBlock)
    DyInvHalf = 0.5/Dy_BLK(iBlock)
    DzInvHalf = 0.5/Dz_BLK(iBlock)

    do k=1,nK; do j=1,nJ; do i=1,nI
       dB1nEast = DxInvHalf*&
            (RightState_VX(Bx_,i,j,k)-LeftState_VX(Bx_,i,j,k))

       dB1nWest = DxInvHalf*&
            (RightState_VX(Bx_,i+1,j,k)-LeftState_VX(Bx_,i+1,j,k))

       dB1nSouth = DyInvHalf* &
            (RightState_VY(By_,i,j,k)-LeftState_VY(By_,i,j,k))

       dB1nNorth = DyInvHalf* &
            (RightState_VY(By_,i,j+1,k)-LeftState_VY(By_,i,j+1,k))

       dB1nBot = DzInvHalf * &
            (RightState_VZ(Bz_,i,j,k)-LeftState_VZ(Bz_,i,j,k))

       dB1nTop = DzInvHalf * &
            (RightState_VZ(Bz_,i,j,k+1)-LeftState_VZ(Bz_,i,j,k+1))

       DivBInternal_C(i,j,k) = 2*(&
            DxInvHalf*(LeftState_VX(Bx_,i+1,j,k) -RightState_VX(Bx_,i,j,k))+&
            DyInvHalf*(LeftState_VY(By_,i,j+1,k) -RightState_VY(By_,i,j,k))+&
            DzInvHalf*(LeftState_VZ(Bz_,i,j,k+1) -RightState_VZ(Bz_,i,j,k)))

       DivB1_GB(i,j,k,iBlock)  = DivBInternal_C(i,j,k) + &
            dB1nEast + dB1nWest + dB1nSouth + dB1nNorth + dB1nTop + dB1nBot

       if(.not.IsMhd) CYCLE

       Source_VC(rhoUx_:rhoUz_,i,j,k) = Source_VC(rhoUx_:rhoUz_,i,j,k) &
            -B0_DX(:,i,j,k)*dB1nEast    &
            -B0_DX(:,i+1,j,k)*dB1nWest  &
            -B0_DY(:,i,j,k)*dB1nSouth   &
            -B0_DY(:,i,j+1,k)*dB1nNorth &
            -B0_DZ(:,i,j,k)*dB1nBot     &
            -B0_DZ(:,i,j,k+1)*dB1nTop               
    end do; end do; end do
    if((.not.IsMhd).or.(.not.UseB0))return
    do k=1,nK; do j=1,nJ; do i=1,nI
       Source_VC(rhoUx_:rhoUz_,i,j,k) = Source_VC(rhoUx_:rhoUz_,i,j,k)  - &
            DivBInternal_C(i,j,k)*B0_DGB(:,i,j,k,iBlock)
    end do; end do; end do
  end subroutine calc_divb_source
  !===========================================================================
  subroutine calc_divb_source_covar

    use ModCovariant, ONLY: FaceAreaI_DFB, FaceAreaJ_DFB, FaceAreaK_DFB
    real :: FaceArea_D(3), vInvHalf
    real :: B1nJumpL, B1nJumpR, DivBInternal_C(1:nI,1:nJ,1:nK)
    integer :: i,j,k

    do k=1,nK; do j=1,nJ; do i=1,nI
       VInvHalf=vInv_CB(i,j,k,iBlock)*cHalf
       FaceArea_D=FaceAreaI_DFB(:,i,j,k,iBlock)
       B1nJumpL =VInvHalf*&
            (FaceArea_D(1)*(RightState_VX(Bx_,i,j,k)-LeftState_VX(Bx_,i,j,k))+&
            FaceArea_D(2)*(RightState_VX(By_,i,j,k)-LeftState_VX(By_,i,j,k))+&
            FaceArea_D(3)*(RightState_VX(Bz_,i,j,k)-LeftState_VX(Bz_,i,j,k)))
       DivBInternal_C(i,j,k)=&
            -(FaceArea_D(1)*RightState_VX(Bx_,i,j,k)+&
            FaceArea_D(2)*RightState_VX(By_,i,j,k)+&
            FaceArea_D(3)*RightState_VX(Bz_,i,j,k))

       FaceArea_D=FaceAreaI_DFB(:,i+1,j,k,iBlock)
       B1nJumpR =  VInvHalf*&
            (FaceArea_D(1)*(RightState_VX(Bx_,i+1,j,k)-LeftState_VX(Bx_,i+1,j,k))+&
            FaceArea_D(2)*(RightState_VX(By_,i+1,j,k)-LeftState_VX(By_,i+1,j,k))+&
            FaceArea_D(3)*(RightState_VX(Bz_,i+1,j,k)-LeftState_VX(Bz_,i+1,j,k)))
       DivBInternal_C(i,j,k)=DivBInternal_C(i,j,k)+&
            (FaceArea_D(1)*LeftState_VX(Bx_,i+1,j,k)+&
            FaceArea_D(2)*LeftState_VX(By_,i+1,j,k)+&
            FaceArea_D(3)*LeftState_VX(Bz_,i+1,j,k))

       DivB1_GB(i,j,k,iBlock)  = B1nJumpL + B1nJumpR

       if(.not.IsMhd) CYCLE

       Source_VC(rhoUx_:rhoUz_,i,j,k) = Source_VC(rhoUx_:rhoUz_,i,j,k) &
            -B0_DX(:,i,j,k)*B1nJumpL   &
            -B0_DX(:,i+1,j,k)*B1nJumpR

    end do; end do; end do

    do k=1,nK; do j=1,nJ; do i=1,nI 
       VInvHalf=vInv_CB(i,j,k,iBlock)*cHalf
       FaceArea_D=FaceAreaJ_DFB(:,i,j,k,iBlock)
       B1nJumpL = VInvHalf*&
            (FaceArea_D(1)*(RightState_VY(Bx_,i,j,k)-LeftState_VY(Bx_,i,j,k))+&
            FaceArea_D(2)*(RightState_VY(By_,i,j,k)-LeftState_VY(By_,i,j,k))+&
            FaceArea_D(3)*(RightState_VY(Bz_,i,j,k)-LeftState_VY(Bz_,i,j,k)))
       DivBInternal_C(i,j,k)=DivBInternal_C(i,j,k)-&
            (FaceArea_D(1)*RightState_VY(Bx_,i,j,k)+&
            FaceArea_D(2)*RightState_VY(By_,i,j,k)+&
            FaceArea_D(3)*RightState_VY(Bz_,i,j,k))

       FaceArea_D=FaceAreaJ_DFB(:,i,j+1,k,iBlock)
       B1nJumpR = VInvHalf*&
            (FaceArea_D(1)*(RightState_VY(Bx_,i,j+1,k)-LeftState_VY(Bx_,i,j+1,k))+&
            FaceArea_D(2)*(RightState_VY(By_,i,j+1,k)-LeftState_VY(By_,i,j+1,k))+&
            FaceArea_D(3)*(RightState_VY(Bz_,i,j+1,k)-LeftState_VY(Bz_,i,j+1,k)))
       DivBInternal_C(i,j,k)=DivBInternal_C(i,j,k)+&
            (FaceArea_D(1)*LeftState_VY(Bx_,i,j+1,k)+&
            FaceArea_D(2)*LeftState_VY(By_,i,j+1,k)+&
            FaceArea_D(3)*LeftState_VY(Bz_,i,j+1,k))

       DivB1_GB(i,j,k,iBlock)  = DivB1_GB(i,j,k,iBlock) &
            + B1nJumpL + B1nJumpR

       if(.not.IsMhd) CYCLE

       Source_VC(rhoUx_:rhoUz_,i,j,k) = Source_VC(rhoUx_:rhoUz_,i,j,k)&
            -B0_DY(:,i,j,k)*B1nJumpL &
            -B0_DY(:,i,j+1,k)*B1nJumpR
    end do; end do; end do

    do k=1,nK; do j=1,nJ; do i=1,nI 
       VInvHalf=vInv_CB(i,j,k,iBlock)*cHalf
       FaceArea_D=FaceAreaK_DFB(:,i,j,k,iBlock)
       B1nJumpL = VInvHalf*&
            (FaceArea_D(1)*(RightState_VZ(Bx_,i,j,k)-LeftState_VZ(Bx_,i,j,k))+&
            FaceArea_D(2)*(RightState_VZ(By_,i,j,k)-LeftState_VZ(By_,i,j,k))+&
            FaceArea_D(3)*(RightState_VZ(Bz_,i,j,k)-LeftState_VZ(Bz_,i,j,k)))
       DivBInternal_C(i,j,k)=DivBInternal_C(i,j,k)-&
            (FaceArea_D(1)*RightState_VZ(Bx_,i,j,k)+&
            FaceArea_D(2)*RightState_VZ(By_,i,j,k)+&
            FaceArea_D(3)*RightState_VZ(Bz_,i,j,k))

       FaceArea_D=FaceAreaK_DFB(:,i,j,k+1,iBlock)
       B1nJumpR = VInvHalf*&
            (FaceArea_D(1)*(RightState_VZ(Bx_,i,j,k+1)-LeftState_VZ(Bx_,i,j,k+1))+&
            FaceArea_D(2)*(RightState_VZ(By_,i,j,k+1)-LeftState_VZ(By_,i,j,k+1))+&
            FaceArea_D(3)*(RightState_VZ(Bz_,i,j,k+1)-LeftState_VZ(Bz_,i,j,k+1)))
       DivBInternal_C(i,j,k)=(DivBInternal_C(i,j,k)+&
            (FaceArea_D(1)*LeftState_VZ(Bx_,i,j,k+1)+&
            FaceArea_D(2)*LeftState_VZ(By_,i,j,k+1)+&
            FaceArea_D(3)*LeftState_VZ(Bz_,i,j,k+1)))*&
            vInv_CB(i,j,k,iBlock)

       DivB1_GB(i,j,k,iBlock)  = DivB1_GB(i,j,k,iBlock) &
            + B1nJumpL + B1nJumpR

       if(.not.IsMhd) CYCLE

       Source_VC(rhoUx_:rhoUz_,i,j,k) = Source_VC(rhoUx_:rhoUz_,i,j,k)&
            -B0_DZ(:,i,j,k)*B1nJumpL &
            -B0_DZ(:,i,j,k+1)*B1nJumpR
    end do; end do; end do

    do k=1,nK; do j=1,nJ; do i=1,nI 
       DivB1_GB(i,j,k,iBlock)  = DivB1_GB(i,j,k,iBlock)+&
            DivBInternal_C(i,j,k)
    end do; end do; end do
    if((.not.IsMhd).or.(.not.UseB0))return
    do k=1,nK; do j=1,nJ; do i=1,nI 
       Source_VC(rhoUx_:rhoUz_,i,j,k) = Source_VC(rhoUx_:rhoUz_,i,j,k) &
            -DivBInternal_C(i,j,k)*B0_DGB(:,i,j,k,iBlock)            
    end do; end do; end do
  end subroutine calc_divb_source_covar
  !===========================================================================
 
  subroutine write_source(String)
    character(len=*), intent(in) :: String
    write(*,'(a,es13.5)') NameSub//": "//String//" S(VarTest)=",&
         Source_VC(VarTest,iTest,jTest,kTest) 
  end subroutine write_source

end subroutine calc_sources

!=============================================================================

subroutine calc_divb(iBlock)

  use ModMain, ONLY : &              
       UseDivbDiffusion,&            !^CFG IF DIVBDIFFUSE
       nI,nJ,nK,test_string
  use ModVarIndexes, ONLY : Bx_,By_,Bz_
  use ModAdvance, ONLY : DivB1_GB,State_VGB, &
       LeftState_VX,RightState_VX,&
       LeftState_VY,RightState_VY,&
       LeftState_VZ,RightState_VZ
  use ModNumConst
  use ModGeometry, ONLY : dx_BLK,dy_BLK,dz_BLK 
  implicit none

  integer, intent(in) :: iBlock
  !---------------------------------------------------------------------------

  !\
  ! Calculate div B for a block and store result into DivB1_GB
  !/
  if(index(test_string,'DIVB_CD')>0)then
     ! Use central differencing if test string contains DIVB_CD
     DivB1_GB(1:nI,1:nJ,1:nK,iBlock) = cHalf*(&
          (State_VGB(Bx_, 2:nI+1, 1:nJ  , 1:nK  ,iBlock)- &
          State_VGB(Bx_, 0:nI-1, 1:nJ  , 1:nK  ,iBlock))/dx_BLK(iBlock)+&
          (State_VGB(By_, 1:nI  , 2:nJ+1, 1:nK  ,iBlock)- &
          State_VGB(By_, 1:nI  , 0:nJ-1, 1:nK  ,iBlock))/dy_BLK(iBlock)+&
          (State_VGB(Bz_, 1:nI  , 1:nJ  , 2:nK+1,iBlock)- &
          State_VGB(Bz_, 1:nI  , 1:nJ  , 0:nK-1,iBlock))/dz_BLK(iBlock))
  else
     ! Compute divB using averaged and conservatively corrected 
     ! left and right values
     
     DivB1_GB(1:nI,1:nJ,1:nK,iBlock) = &
          cHalf *( &
          ((LeftState_VX(Bx_,2:nI+1,1:nJ,1:nK)+    &
            RightState_VX(Bx_,2:nI+1,1:nJ,1:nK))-   &
           (LeftState_VX(Bx_,1:nI,1:nJ,1:nK)+    &
            RightState_VX(Bx_,1:nI,1:nJ,1:nK))  )/dx_BLK(iBlock)+  &
          ((LeftState_VY(By_,1:nI,2:nJ+1,1:nK)+    &
            RightState_VY(By_,1:nI,2:nJ+1,1:nK))-   &
           (LeftState_VY(By_,1:nI,1:nJ,1:nK)+    &
            RightState_VY(By_,1:nI,1:nJ,1:nK))  )/dy_BLK(iBlock)+  &
          ((LeftState_VZ(Bz_,1:nI,1:nJ,2:nK+1)+    &
            RightState_VZ(Bz_,1:nI,1:nJ,2:nK+1))-   &
           (LeftState_VZ(Bz_,1:nI,1:nJ,1:nK)+    &
            RightState_VZ(Bz_,1:nI,1:nJ,1:nK))  )/dz_BLK(iBlock) )
  endif                                         
end subroutine calc_divb

!==============================================================================

subroutine get_current(i,j,k,iBlock,Current_D)

  ! Calculate the current in a cell of a block

  use ModAdvance,  ONLY: State_VGB, Bx_, By_, Bz_
  use ModGeometry, ONLY: True_Cell, Dx_BLK, Dy_BLK, Dz_BLK
  use ModCovariant,ONLY: UseCovariant                

  implicit none
  integer, intent(in) :: i,j,k,iBlock
  real,    intent(out):: Current_D(3)

  real :: DxInvHalf, DyInvHalf, DzInvHalf
  !----------------------------------------------------------------------------

  ! Exclude cells next to the body because they produce incorrect currents
  if(.not.all(True_Cell(i-1:i+1,j-1:j+1,k-1:k+1,iBlock)))then
     Current_D = 0.0
     RETURN
  endif

  if(UseCovariant)then                               
     call covariant_curlb(i,j,k,iBlock,Current_D,.true.)
     RETURN
  end if                                             

  DxInvHalf = 0.5/Dx_BLK(iBlock)
  DyInvHalf = 0.5/Dy_BLK(iBlock)
  DzInvHalf = 0.5/Dz_BLK(iBlock)

  Current_D(1) = &
       (State_VGB(Bz_,i,j+1,k,iBlock) &
       -State_VGB(Bz_,i,j-1,k,iBlock))*DyInvHalf - &
       (State_VGB(By_,i,j,k+1,iBlock) &
       -State_VGB(By_,i,j,k-1,iBlock))*DzInvHalf

  Current_D(2) = &
       (State_VGB(Bx_,i,j,k+1,iBlock) &
       -State_VGB(Bx_,i,j,k-1,iBlock))*DzInvHalf- &
       (State_VGB(Bz_,i+1,j,k,iBlock) &
       -State_VGB(Bz_,i-1,j,k,iBlock))*DxInvHalf

  Current_D(3) = &
       (State_VGB(By_,i+1,j,k,iBlock) &
       -State_VGB(By_,i-1,j,k,iBlock))*DxInvHalf- &
       (State_VGB(Bx_,i,j+1,k,iBlock) &
       -State_VGB(Bx_,i,j-1,k,iBlock))*DyInvHalf

end subroutine get_current
