!^CFG COPYRIGHT UM
subroutine calc_sources

  use ModProcMH
  use ModMain
  use ModVarIndexes
  use ModGeometry,      ONLY: dx_BLK, dy_BLK, dz_BLK, R_BLK,&
       body_BLK, Rmin_BLK, vInv_CB, TypeGeometry, y_BLK, true_cell
  use ModGeometry,      ONLY: R2_BLK                   !^CFG IF SECONDBODY
  use ModAdvance
  use ModParallel,      ONLY: NOBLK, neiLEV, &
       neiLtop,neiLbot,neiLeast,neiLwest,neiLnorth,neiLsouth
  use ModPhysics
  use ModNumConst
  use ModResistivity,   ONLY: UseResistivity, Eta_GB   !^CFG IF DISSFLUX
  use ModUser,          ONLY: user_calc_sources
  use ModCoordTransform
  use ModImplicit,      ONLY: UseFullImplicit            !^CFG IF IMPLICIT
  use ModRadDiffusion,  ONLY: calc_source_rad_diffusion  !^CFG IF IMPLICIT
  use ModMultiFluid
  use ModPointImplicit, ONLY: UsePointImplicit, UsePointImplicit_B
  use ModMultiIon,      ONLY: multi_ion_source_expl, multi_ion_source_impl
  use ModCovariant,     ONLY: UseCovariant 
  use ModCurrent,       ONLY: get_current
  use ModWaves,         ONLY: UseWavePressure, GammaWave, DivU_C

  implicit none

  integer :: i, j, k, iDim, iVar
  logical :: UseRzGeometry
  real :: Pe, Pwave, DivU
  real :: Coef, GammaMinus1

  ! Variable for div B diffusion
  real :: Dr

  real,dimension(3)::CurlB0CrossB_D

  ! Variables needed for Boris source terms also used for div(u)
  real :: FullB_DC(nDim,nI,nJ,nK),FullBx, FullBy, FullBz, Ux, Uy, Uz, RhoInv
  real :: E_D(3), DivE

  ! Variables needed for Joule heating
  real :: Current_D(3), JouleHeating, HeatExchange

  ! Varibles needed for anisotropic pressure
  real :: b_D(3), GradU_DD(3,3), bDotGradparU, AnisoRelaxation
  logical :: IsFirehose, IsMirror 

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

  UseRzGeometry = TypeGeometry == 'rz'

  Source_VC   = 0.0

  ! Calculate source terms for ion pressure
  if(UseNonconservative)then
     do iFluid = 1, nFluid
        call select_fluid

        if(UseAnisopressure)then
           ! Source terms for anisotropic pressure equations
           do k=1,nK; do j=1,nJ; do i=1,nI
              if(.not.true_cell(i,j,k,iBlock)) CYCLE
              ! Calculate bDotGradparU = b dot (b matmul GradU)

              ! Calculate unit vector parallel with full B field
              b_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)
              if(UseB0) b_D = b_D + B0_DGB(:,i,j,k,iBlock)
              b_D = b_D/sqrt(max(1e-30, sum(B_D**2)))

              ! Calculate gradient tensor of velocity
              GradU_DD(x_,:) = (0.5*(LeftState_VX(Ux_:Uz_,i+1,j,k) &
                   + RightState_VX(Ux_:Uz_,i+1,j,k))               &
                   - 0.5*(LeftState_VX(Ux_:Uz_,i,j,k)              &
                   + RightState_VX(Ux_:Uz_,i,j,k)))/dx_BLK(iBlock)
              if(nJ > 1) then
                 GradU_DD(y_,:) = &
                      ( 0.5*(LeftState_VY(Ux_:Uz_,i,j+1,k)  &
                      +     RightState_VY(Ux_:Uz_,i,j+1,k)) &
                      - 0.5*(LeftState_VY(Ux_:Uz_,i,j,k)    &
                      +     RightState_VY(Ux_:Uz_,i,j,k)))/dy_BLK(iBlock)
              else
                 GradU_DD(y_,:) = 0.0
              end if
              if(nK > 1) then
                 GradU_DD(z_,:) = &
                      ( 0.5*(LeftState_VZ(Ux_:Uz_,i,j,k+1) &
                      +     RightState_VZ(Ux_:Uz_,i,j,k+1)) &
                      - 0.5*(LeftState_VZ(Ux_:Uz_,i,j,k) &
                      +     RightState_VZ(Ux_:Uz_,i,j,k)))/dz_BLK(iBlock)
              else
                 GradU_DD(z_,:) = 0.0
              end if

              ! Calculate b.grad u.b
              bDotGradparU = dot_product(b_D, matmul(b_D, GradU_DD))

              ! p parallel: -2*Ppar*b.(b.(Grad U))
              Source_VC(Ppar_,i,j,k) = Source_VC(Ppar_,i,j,k) &
                   - 2*State_VGB(Ppar_,i,j,k,iBlock)*bDotGradparU
              ! p perpendicular: Pperp*b.(b.(GradU))
              Source_VC(Pperp_,i,j,k) = Source_VC(Pperp_,i,j,k) &
                   + State_VGB(Pperp_,i,j,k,iBlock)*bDotGradparU
           end do; end do; end do

           if(DoTestMe .and. (VarTest == Ppar_ .or. VarTest == Pperp_)) &
                call write_source('After bDotGradparU')

           ! For Pperp equation, GammaPerp = 2, thus -Pperp*Div(U)
           GammaMinus1 = 1
        else
           GammaMinus1 = g - 1
        end if

        ! Adiabatic heating: -(g-1)*P*Div(U)
        do k=1,nK; do j=1,nJ; do i=1,nI
           if(.not.true_cell(i,j,k,iBlock)) CYCLE
           DivU = uDotArea_XI(i+1,j,k,iFluid) - uDotArea_XI(i,j,k,iFluid)
           if(nJ > 1) DivU = DivU + &
                uDotArea_YI(i,j+1,k,iFluid) - uDotArea_YI(i,j,k,iFluid)
           if(nK > 1) DivU = DivU + &
                uDotArea_ZI(i,j,k+1,iFluid) - uDotArea_ZI(i,j,k,iFluid)
           DivU = vInv_CB(i,j,k,iBlock)*DivU
           Source_VC(iP,i,j,k) = Source_VC(iP,i,j,k) &
                - GammaMinus1*State_VGB(iP,i,j,k,iBlock)*DivU
        end do; end do; end do

        if(DoTestMe .and. VarTest==iP)call write_source('After p div U')

     end do
  end if

  ! Joule heating: dP/dt += (gamma-1)*eta*j**2    !^CFG IF DISSFLUX BEGIN
  if(UseResistivity .and. .not.UseMultiIon .and. &
       (UseElectronPressure .or. UseNonConservative))then  

     do k=1,nK; do j=1,nJ; do i=1,nI
        if(.not.true_cell(i,j,k,iBlock)) CYCLE
        call get_current(i,j,k,iBlock,Current_D)
        JouleHeating = (g-1) * Eta_GB(i,j,k,iBlock) * sum(Current_D**2)
        if(UseElectronPressure) then
           ! For single ion fluid the ion-electron collision results in a 
           ! heat exchange term for the electron pressure 
           ! See eq. 4.124c in Schunk and Nagy.
           HeatExchange =0.0
           if(.not.UseMultiSpecies)then
              ! Explicit heat exchange
              HeatExchange = (g-1) * Eta_GB(i,j,k,iBlock) * &
                   3*State_VGB(Rho_,i,j,k,iBlock)*(1./IonMassPerCharge**2)

              ! Point-implicit correction for stability: H' = H/(1+dt*H)
              if(UseAnisoPressure)then
                 HeatExchange = &
                   HeatExchange/(1 + Cfl*HeatExchange*time_BLK(i,j,k,iBlock)) &
                   *(2./3.*State_VGB(Pperp_,i,j,k,iBlock) &
                   + 1./3.*State_VGB(Ppar_,i,j,k,iBlock) &
                   - State_VGB(Pe_,i,j,k,iBlock))
              else
                 HeatExchange = HeatExchange / &
                      (1 + Cfl*HeatExchange*time_BLK(i,j,k,iBlock)) &
                      *(State_VGB(P_,i,j,k,iBlock) &
                      - State_VGB(Pe_,i,j,k,iBlock))
              end if

           end if

           Source_VC(Pe_,i,j,k) = Source_VC(Pe_,i,j,k) &
                + JouleHeating + HeatExchange

           ! Heat exchange applies to ions too
           Source_VC(P_,i,j,k) = Source_VC(P_,i,j,k) - HeatExchange

           if(UseAnisoPressure) then
              ! Heat exchange for parallel ion pressure
              Source_VC(Ppar_,i,j,k) = Source_VC(Ppar_,i,j,k) - HeatExchange
              
              ! Relaxation term due to collisions
              AnisoRelaxation = Eta_GB(i,j,k,iBlock) &
                   *State_VGB(Rho_,i,j,k,iBlock)/IonMassPerCharge**2 &
                   *(State_VGB(Ppar_,i,j,k,iBlock) &
                   - State_VGB(Pperp_,i,j,k,iBlock))
              Source_VC(Pperp_,i,j,k) = Source_VC(Pperp_,i,j,k) &
                   + 1./3.*AnisoRelaxation
              Source_VC(Ppar_,i,j,k)  = Source_VC(Ppar_,i,j,k)  &
                   - 2./3.*AnisoRelaxation
           end if

           ! Remove Joule heating and apply heat exchange to ion energy
           Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) &
                - inv_gm1*(JouleHeating + HeatExchange)
        else
           Source_VC(P_,i,j,k) = Source_VC(P_,i,j,k) + JouleHeating
        end if
     end do; end do; end do

     if(DoTestMe.and.VarTest==P_)call write_source('After eta j')

  end if                                       !^CFG END DISSFLUX


  if(UseAnisoPressure .and. TauWaveParticle > -1.0)then
     ! "artificial" pressure relaxation for anisotropic pressure 
     ! due to wave-particle interaction 
     ! only done in unstable regions
     do k=1,nK; do j=1,nJ; do i=1,nI
        if(.not.true_cell(i,j,k,iBlock)) CYCLE
        B_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)
        if(UseB0) B_D = B_D + B0_DGB(:,i,j,k,iBlock)         
        IsFirehose = (State_VGB(Ppar_,i,j,k,iBlock) &
             - State_VGB(Pperp_,i,j,k,iBlock)) > sum(B_D**2)
        IsMirror = State_VGB(Pperp_,i,j,k,iBlock)/State_VGB(Ppar_,i,j,k,iBlock) &
             > 1 + sum(B_D**2)/2./State_VGB(Pperp_,i,j,k,iBlock)
        if(IsFirehose .or. IsMirror)then           
           ! Point-implicit for stability
           AnisoRelaxation = 1./(TauWaveParticle + Cfl*time_BLK(i,j,k,iBlock)) &
                *(State_VGB(Ppar_,i,j,k,iBlock) - State_VGB(Pperp_,i,j,k,iBlock))
           Source_VC(Pperp_,i,j,k) = Source_VC(Pperp_,i,j,k) + AnisoRelaxation/3.
           Source_VC(Ppar_,i,j,k)  = Source_VC(Ppar_,i,j,k) - 2./3*AnisoRelaxation
        end if
     end do; end do; end do  
  end if


  if(UseWavePressure)then
     do k = 1, nK; do j = 1, nJ; do i = 1, nI
        if(.not.true_cell(i,j,k,iBlock)) CYCLE
        DivU            = uDotArea_XI(i+1,j,k,1) - uDotArea_XI(i,j,k,1)
        if(nJ > 1) DivU = DivU + uDotArea_YI(i,j+1,k,1) - uDotArea_YI(i,j,k,1)
        if(nK > 1) DivU = DivU + uDotArea_ZI(i,j,k+1,1) - uDotArea_ZI(i,j,k,1)
        DivU = vInv_CB(i,j,k,iBlock)*DivU

        ! Store div U so it can be used in ModWaves
        DivU_C(i,j,k) = DivU

        Pwave = (GammaWave - 1) &
             *sum(State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock))

        do iVar = WaveFirst_, WaveLast_
           Source_VC(iVar,i,j,k) = Source_VC(iVar,i,j,k) &
                - DivU*(GammaWave - 1)*State_VGB(iVar,i,j,k,iBlock)
        end do
        Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) + DivU*Pwave

        ! Add "geometrical source term" p/r to the radial momentum equation
        ! The "radial" direction is along the Y axis
        ! NOTE: here we have to use signed radial distance!
        if(UseRzGeometry) Source_VC(RhoUy_,i,j,k) = Source_VC(RhoUy_,i,j,k) &
             + Pwave/y_BLK(i,j,k,iBlock)
     end do; end do; end do
  end if

  if(UseElectronPressure)then
     do k = 1, nK; do j = 1, nJ; do i = 1, nI
        if(.not.true_cell(i,j,k,iBlock)) CYCLE
        DivU = uDotArea_XI(i+1,j,k,eFluid_) - uDotArea_XI(i,j,k,eFluid_)
        if(nJ > 1) DivU = DivU &
             + uDotArea_YI(i,j+1,k,eFluid_) - uDotArea_YI(i,j,k,eFluid_)
        if(nK > 1) DivU = DivU &
             + uDotArea_ZI(i,j,k+1,eFluid_) - uDotArea_ZI(i,j,k,eFluid_)
        DivU = vInv_CB(i,j,k,iBlock)*DivU

        Pe = State_VGB(Pe_,i,j,k,iBlock)

        ! Adiabatic heating for electron pressure: -(g-1)*Pe*Div(U)
        Source_VC(Pe_,i,j,k) = Source_VC(Pe_,i,j,k) - (g-1)*Pe*DivU

        ! The energy equation contains the work of electron pressure
        ! -u.grad Pe = -div(u Pe) + Pe div u
        ! The -div(u Pe) is implemented as a flux in ModFaceFlux. 
        ! Here we add the Pe div u source term
        Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) + Pe*DivU

        ! Add "geometrical source term" p/r to the radial momentum equation
        ! The "radial" direction is along the Y axis
        ! NOTE: here we have to use signed radial distance!
        if(UseRzGeometry) Source_VC(RhoUy_,i,j,k) = Source_VC(RhoUy_,i,j,k) &
             + Pe/y_BLK(i,j,k,iBlock)
     end do; end do; end do

     if(DoTestMe.and.VarTest==Pe_)call write_source('After Pe div Ue')
  end if

  if(TypeGeometry == 'rz')then
     if(UseB)call stop_mpi('RZ geometry is not implemented for MHD')

     ! Add "geometrical source term" p/r to the radial momentum equation.
     ! The azimuthal component of the velocity is assumed to be zero.
     ! The axis is along X, the "radial" direction is along the Y axis
     ! NOTE: here we have to use signed radial distance!

     do k=1,nK; do j=1, nJ; do i=1, nI
        if(.not.true_cell(i,j,k,iBlock)) CYCLE
        Source_VC(iRhoUy_I,i,j,k) = Source_VC(iRhoUy_I,i,j,k) &
             + State_VGB(iP_I,i,j,k,iBlock) / y_BLK(i,j,k,iBlock)
     end do; end do; end do
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
        if(.not.true_cell(i,j,k,iBlock)) CYCLE
        RhoInv=cOne/State_VGB(rho_,i,j,k,iBlock)
        Source_VC(Bx_:Bz_,i,j,k)    = Source_VC(Bx_:Bz_,i,j,k) &
             -DivB1_GB(i,j,k,iBlock)* &
             State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock)*RhoInv

        if(.not. IsMhd) CYCLE

        Source_VC(rhoUx_:rhoUz_,i,j,k) = Source_VC(RhoUx_:RhoUz_,i,j,k) -&
             DivB1_GB(i,j,k,iBlock)* &
             State_VGB(Bx_:Bz_,i,j,k,iBlock)
        Source_VC(Energy_,i,j,k)     = Source_VC(Energy_,i,j,k) &
             -DivB1_GB(i,j,k,iBlock)* &
             sum(State_VGB(Bx_:Bz_,i,j,k,iBlock)*&
             State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock))*RhoInv
     end do;end do;end do

     if (UseB0Source) then
        do k=1,nK; do j=1,nJ; do i=1,nI
           if(.not.true_cell(i,j,k,iBlock)) CYCLE
           Source_VC(RhoUx_:RhoUz_,i,j,k) = Source_VC(rhoUx_:rhoUz_,i,j,k) &
                - State_VGB(Bx_:Bz_,i,j,k,iBlock)*DivB0_CB(i,j,k,iBlock) &
                + cross_product( &
                State_VGB(Bx_:Bz_,i,j,k,iBlock), CurlB0_DCB(:,i,j,k,iBlock))
        end do; end do; end do
     end if
  else
     if(UseB)call calc_divb(iBlock)
  end if

  if(UseCurlB0)then
     do k=1,nK; do j=1,nJ; do i=1,nI
        if(.not.true_cell(i,j,k,iBlock)) CYCLE
        if(R_BLK(i,j,k,iBlock)<rCurrentFreeB0)CYCLE
        CurlB0CrossB_D = cross_product(&
             CurlB0_DCB(:,i,j,k,iBlock),&
             State_VGB(Bx_:Bz_,i,j,k,iBlock) + B0_DGB(:,i,j,k,iBlock))
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

     Coef = (boris_cLIGHT_factor**2 - 1.0)*inv_c2LIGHT
     FullB_DC = State_VGB(Bx_:Bz_,1:nI,1:nJ,1:nK,iBlock)
     if(UseB0)FullB_DC = FullB_DC + B0_DGB(:,1:nI,1:nJ,1:nK,iBlock) 
     do k=1,nK; do j=1,nJ; do i=1,nI
        if(.not.true_cell(i,j,k,iBlock)) CYCLE
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
              if(.not.true_cell(i,j,k,iBlock)) CYCLE
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

  if(UseMultiIon)then
     ! Add momentum source terms containing the gradient of electron pressure
     call multi_ion_source_expl(iBlock)

     ! Add stiff momentum source terms (uPlus - Uion) and artificial friction
     ! Explicit evaluation of these source terms is for code development only
     if (.not. (UsePointImplicit .and. UsePointImplicit_B(iBlock)) ) &
          call multi_ion_source_impl

     if(DoTestMe) call write_source('After MultiIon sources')
  end if

  ! Add JxB term for nonconservativ MHD scheme (like LFM)
  if(UseB .and. .not.IsMhd .and. .not.UseMultiIon)then
     call multi_ion_source_expl(iBlock)

     if(DoTestMe) call write_source('After JxB term')
  end if


  !^CFG IF  IMPLICIT BEGIN
  if(UseRadDiffusion .and. UseFullImplicit) &
       call calc_source_rad_diffusion(iBlock)
  !^CFG END IMPLICIT

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
       if(.not.true_cell(i,j,k,iBlock)) CYCLE
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

    if((.not.IsMhd).or.(.not.UseB0)) RETURN

    do k=1,nK; do j=1,nJ; do i=1,nI
       if(.not.true_cell(i,j,k,iBlock)) CYCLE
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
    !------------------------------------------------------------------------
    
    do k=1,nK; do j=1,nJ; do i=1,nI
       if(.not.true_cell(i,j,k,iBlock)) CYCLE
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
       DivBInternal_C(i,j,k) = DivBInternal_C(i,j,k)+&
            (FaceArea_D(1)*LeftState_VX(Bx_,i+1,j,k)+&
            FaceArea_D(2)*LeftState_VX(By_,i+1,j,k)+&
            FaceArea_D(3)*LeftState_VX(Bz_,i+1,j,k))

       DivB1_GB(i,j,k,iBlock)  = B1nJumpL + B1nJumpR

       if(.not.IsMhd) CYCLE

       Source_VC(rhoUx_:rhoUz_,i,j,k) = Source_VC(rhoUx_:rhoUz_,i,j,k) &
            - B0_DX(:,i,j,k)*B1nJumpL   &
            - B0_DX(:,i+1,j,k)*B1nJumpR

    end do; end do; end do

    do k=1,nK; do j=1,nJ; do i=1,nI 
       if(.not.true_cell(i,j,k,iBlock)) CYCLE
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
       if(.not.true_cell(i,j,k,iBlock)) CYCLE
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
       if(.not.true_cell(i,j,k,iBlock)) CYCLE
       DivB1_GB(i,j,k,iBlock)  = DivB1_GB(i,j,k,iBlock)+&
            DivBInternal_C(i,j,k)
    end do; end do; end do

    if((.not.IsMhd).or.(.not.UseB0))RETURN

    do k=1,nK; do j=1,nJ; do i=1,nI 
       if(.not.true_cell(i,j,k,iBlock)) CYCLE
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

  ! Calculate div B for a block and store result into DivB1_GB
  ! Compute divB using averaged and conservatively corrected 
  ! left and right values

  use ModMain,       ONLY: nI,nJ,nK
  use ModVarIndexes, ONLY: Bx_,By_,Bz_
  use ModGeometry,   ONLY: dx_BLK, dy_BLK, dz_BLK 
  use ModAdvance,    ONLY: DivB1_GB,State_VGB, &
       LeftState_VX, RightState_VX, &
       LeftState_VY, RightState_VY, &
       LeftState_VZ, RightState_VZ

  implicit none

  integer, intent(in) :: iBlock

  integer:: i, j, k
  real   :: DivB, InvDx, InvDy, InvDz
  !---------------------------------------------------------------------------

  InvDx            = 1/dx_BLK(iBlock)
  if(nJ > 1) InvDy = 1/dy_BLK(iBlock)
  if(nK > 1) InvDz = 1/dz_BLK(iBlock)

  do k=1, nK; do j=1, nJ; do i=1, nI
     DivB = InvDx* &
          (  LeftState_VX(Bx_,i+1,j,k)  &
          + RightState_VX(Bx_,i+1,j,k)  &
          -  LeftState_VX(Bx_,i,j,k)    &
          - RightState_VX(Bx_,i,j,k) )
  
     if(nJ > 1) DivB = DivB + InvDy* &
          (  LeftState_VY(By_,i,j+1,k)   &
          + RightState_VY(By_,i,j+1,k)   &
          -  LeftState_VY(By_,i,j,k)     &
          - RightState_VY(By_,i,j,k) )
  
     if(nK > 1) DivB = DivB + InvDz* &
          (  LeftState_VZ(Bz_,i,j,k+1)    &
          + RightState_VZ(Bz_,i,j,k+1)    &
          -  LeftState_VZ(Bz_,i,j,k)      &
          - RightState_VZ(Bz_,i,j,k) )

     DivB1_GB(i,j,k,iBlock) = 0.5*DivB

  end do; end do; end do

end subroutine calc_divb

