!^CFG COPYRIGHT UM
!==============================================================================
module ModCalcSource

  implicit none
  save

  private !except

  ! Public methods
  public :: calc_sources
  public :: get_tesi_c

contains

  !============================================================================

  subroutine calc_sources

    use ModProcMH
    use ModMain
    use ModVarIndexes
    use ModGeometry,      ONLY: dx_BLK, dy_BLK, dz_BLK, R_BLK,&
         vInv_CB, y_BLK, true_cell
    use ModAdvance
    use ModPhysics
    use ModUser,          ONLY: user_calc_sources
    use ModCoordTransform
    use ModImplicit,      ONLY: UseFullImplicit            !^CFG IF IMPLICIT
    use ModRadDiffusion,  ONLY: calc_source_rad_diffusion  !^CFG IF IMPLICIT
    use ModMultiFluid
    use ModPointImplicit, ONLY: UsePointImplicit, UsePointImplicit_B
    use ModMultiIon,      ONLY: multi_ion_source_expl, multi_ion_source_impl
    use ModWaves,         ONLY: UseWavePressure, GammaWave, DivU_C
    use ModCoronalHeating,ONLY: UseCoronalHeating, get_block_heating, &
         CoronalHeating_C, UseAlfvenWaveDissipation, WaveDissipation_VC, &
         QeByQtotal
    use ModRadiativeCooling, ONLY: RadCooling_C,UseRadCooling, &
         get_radiative_cooling, add_chromosphere_heating
    use ModChromosphere,  ONLY: DoExtendTransitionRegion, extension_factor, &
         UseChromosphereHeating
    use ModFaceFlux,      ONLY: Pe_G
    use ModHallResist,    ONLY: UseBiermannBattery, IonMassPerCharge_G
    use BATL_lib, ONLY: IsCartesian, IsRzGeometry

    integer :: i, j, k, iVar
    real :: Pe, DivU
    real :: Coef

    ! Variable for div B diffusion

    real :: CurlB0CrossB_D(3)

    ! Variables needed for Boris source terms also used for div(u)
    real :: FullB_DC(nDim,nI,nJ,nK), RhoInv
    real :: E_D(3), DivE

    ! Varibles needed for anisotropic pressure
    real :: b_D(3), GradU_DD(3,3), bDotGradparU

    ! Electron temperature in K:
    real :: TeSi_C(nI,nJ,nK)

    integer:: iBlock

    logical :: DoTest, DoTestMe

    character(len=*), parameter :: NameSub = 'ModCalcSource::calc_sources'
    !--------------------------------------------------------------------------
    iBlock = GlobalBlk

    if(iProc==PROCtest .and. iBlock==BLKtest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    Source_VC = 0.0

    ! Calculate source terms for ion pressure
    if(UseNonconservative .or. UseAnisoPressure)then
       do iFluid = 1, nFluid
          call select_fluid

          if(UseAnisoPressure)then
             ! Source terms for anisotropic pressure equations
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
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

                ! p parallel: -2*ppar*b.(b.(Grad U))
                Source_VC(Ppar_,i,j,k) = Source_VC(Ppar_,i,j,k) &
                     - 2*State_VGB(Ppar_,i,j,k,iBlock)*bDotGradparU
                ! p : 2/3*(pperp - ppar)*b.(b.(GradU))
                !     = (p - ppar)*b.(b.(GradU)) 
                Source_VC(p_,i,j,k) = Source_VC(p_,i,j,k) &
                     + (State_VGB(p_,i,j,k,iBlock) -  &
                     State_VGB(Ppar_,i,j,k,iBlock))*bDotGradparU
             end do; end do; end do

             if(DoTestMe .and. (VarTest == Ppar_ .or. VarTest == p_)) &
                  call write_source('After bDotGradparU')

          end if

          ! Adiabatic heating: -(g-1)*P*Div(U)
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.true_cell(i,j,k,iBlock)) CYCLE

             DivU = uDotArea_XI(i+1,j,k,iFluid) - uDotArea_XI(i,j,k,iFluid)
             if(nJ > 1) DivU = DivU + &
                  uDotArea_YI(i,j+1,k,iFluid) - uDotArea_YI(i,j,k,iFluid)
             if(nK > 1) DivU = DivU + &
                  uDotArea_ZI(i,j,k+1,iFluid) - uDotArea_ZI(i,j,k,iFluid)
             DivU = vInv_CB(i,j,k,iBlock)*DivU
             if(UseAnisoPressure)then
                Source_VC(iP,i,j,k) = Source_VC(iP,i,j,k) &
                     - (State_VGB(iP,i,j,k,iBlock) &
                     - State_VGB(Ppar_,i,j,k,iBlock)/3.0)*DivU
             else
                Source_VC(iP,i,j,k) = Source_VC(iP,i,j,k) &
                     - gm1*State_VGB(iP,i,j,k,iBlock)*DivU
             end if
          end do; end do; end do

          if(DoTestMe .and. VarTest==iP)call write_source('After p div U')

       end do
    end if


    ! Joule heating: dP/dt += (gamma-1)*eta*j**2    !^CFG IF DISSFLUX BEGIN
    ! Heat exchange between electrons and ions (mult-ion is not coded).
    !if(.not.UseMultiIon .and. UseResistivity .and. &
    !     (UseElectronPressure .or. UseNonConservative))then  
    !   call calc_resistivity_source(iBlock)
    !   if(DoTestMe.and.VarTest==P_)call write_source('After resistive src')
    !end if                                       !^CFG END DISSFLUX

    if(UseWavePressure)then
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(.not.true_cell(i,j,k,iBlock)) CYCLE

          DivU            = uDotArea_XI(i+1,j,k,1) - uDotArea_XI(i,j,k,1)
          if(nJ > 1) DivU = DivU + uDotArea_YI(i,j+1,k,1) -uDotArea_YI(i,j,k,1)
          if(nK > 1) DivU = DivU + uDotArea_ZI(i,j,k+1,1) -uDotArea_ZI(i,j,k,1)
          DivU = vInv_CB(i,j,k,iBlock)*DivU

          ! Store div U so it can be used in ModWaves
          DivU_C(i,j,k) = DivU

          do iVar = WaveFirst_, WaveLast_
             Source_VC(iVar,i,j,k) = Source_VC(iVar,i,j,k) &
                  - DivU*(GammaWave - 1)*State_VGB(iVar,i,j,k,iBlock)
          end do

          ! Add "geometrical source term" p/r to the radial momentum equation
          ! The "radial" direction is along the Y axis
          ! NOTE: here we have to use signed radial distance!
          if(IsRzGeometry) &
               Source_VC(RhoUy_,i,j,k) = Source_VC(RhoUy_,i,j,k) &
               + sum(State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock)) &
               *(GammaWave - 1)/y_BLK(i,j,k,iBlock)
       end do; end do; end do
    end if

    if(UseCoronalHeating .and. DoExtendTransitionRegion .or. UseRadCooling) &
         call get_tesi_c(iBlock, TeSi_C)

    if(UseCoronalHeating)then
       call get_block_heating(iBlock)

       if(UseChromosphereHeating.and. DoExtendTransitionRegion)then
          call add_chromosphere_heating(TeSi_C, iBlock)
          do k=1,nK; do j=1,nJ; do i=1,nI
             CoronalHeating_C(i,j,k) = &
                  CoronalHeating_C(i,j,k)/extension_factor(TeSi_C(i,j,k))
          end do; end do; end do
       end if

       if(UseAlfvenWaveDissipation)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             Source_VC(WaveFirst_:WaveLast_,i,j,k) = &
                  Source_VC(WaveFirst_:WaveLast_,i,j,k) &
                  - WaveDissipation_VC(:,i,j,k)
             Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) &
                  - sum(WaveDissipation_VC(:,i,j,k))
          end do; end do; end do
       end if

       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(UseElectronPressure)then
             Source_VC(p_,i,j,k) = Source_VC(p_,i,j,k) &
                  + CoronalHeating_C(i,j,k)*gm1*(1.0-QeByQtotal)
             Source_VC(Pe_,i,j,k) = Source_VC(Pe_,i,j,k) &
                  + CoronalHeating_C(i,j,k)*gm1*QeByQtotal
          else
             Source_VC(p_,i,j,k) = Source_VC(p_,i,j,k) &
                  + CoronalHeating_C(i,j,k)*gm1
          end if
          Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) &
               + CoronalHeating_C(i,j,k)
       end do; end do; end do

    end if

    if(UseRadCooling)then
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          call get_radiative_cooling(i, j, k, iBlock, TeSi_C(i,j,k), &
               RadCooling_C(i,j,k))

          if(UseElectronPressure)then
             Source_VC(Pe_,i,j,k) = Source_VC(Pe_,i,j,k) &
                  + RadCooling_C(i,j,k)*gm1
          else
             Source_VC(p_,i,j,k)  = Source_VC(p_,i,j,k) &
                  + RadCooling_C(i,j,k)*gm1
          end if
          Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) &
               + RadCooling_C(i,j,k)
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
          Source_VC(Pe_,i,j,k) = Source_VC(Pe_,i,j,k) - gm1*Pe*DivU

          if(.not.UseMultiIon)then
             ! Add "geometrical source term" p/r to the radial momentum
             ! equation. The "radial" direction is along the Y axis
             ! NOTE: here we have to use signed radial distance!
             if(IsRzGeometry) Source_VC(RhoUy_,i,j,k) = &
                  Source_VC(RhoUy_,i,j,k) + Pe/y_BLK(i,j,k,iBlock)
          end if
       end do; end do; end do
       if(DoTestMe.and.VarTest==Pe_)call write_source('After Pe div Ue')
    end if

    if(IsRzGeometry)then
       ! The following geometrical source terms are added for the MHD equations
       ! Source[mr]  =(p+B^2/2-Bphi**2+mphi**2/rho)/radius
       ! Source[mphi]=(-mphi*mr/rho+Bphi*Br)/radius(if no angular momentum fix)
       ! Source[Bphi]=((Bphi*mr-Br*mphi)/rho)/radius

       ! The cylindrical coordinates are renamed to X,Y,Z
       ! The cylindrical axis is along X
       ! The "radial" direction is along the Y axis
       ! The azimuthal direction is along the Z axis

       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(.not.true_cell(i,j,k,iBlock)) CYCLE

          ! Source[mr] = (p+mphi**2/rho)/radius
          Source_VC(iRhoUy_I,i,j,k) = Source_VC(iRhoUy_I,i,j,k) &
               + (State_VGB(iP_I,i,j,k,iBlock) &
               +  State_VGB(iRhoUz_I,i,j,k,iBlock)**2 &
               /  State_VGB(iRho_I,i,j,k,iBlock)) &
               / y_BLK(i,j,k,iBlock)

          ! Source[mphi] = (-mphi*mr/rho)/radius
          Source_VC(iRhoUz_I,i,j,k) = Source_VC(iRhoUz_I,i,j,k) &
               - State_VGB(iRhoUz_I,i,j,k,iBlock) &
               * State_VGB(iRhoUy_I,i,j,k,iBlock) &
               /State_VGB(iRho_I,i,j,k,iBlock)/y_BLK(i,j,k,iBlock)

          if(UseB)then
             ! Source[mr] = (B^2/2-Bphi**2)/radius
             Source_VC(RhoUy_,i,j,k) = Source_VC(RhoUy_,i,j,k) &
                  + (0.5*sum(State_VGB(Bx_:Bz_,i,j,k,iBlock)**2) &
                  -  State_VGB(Bz_,i,j,k,iBlock)**2) / y_BLK(i,j,k,iBlock)

             ! Source[mphi]=Bphi*Br/radius
             Source_VC(RhoUz_,i,j,k) = Source_VC(RhoUz_,i,j,k) &
                  + State_VGB(Bz_,i,j,k,iBlock)*State_VGB(By_,i,j,k,iBlock) &
                  / y_BLK(i,j,k,iBlock)

             ! Source[Bphi]=((Bphi*mr-Br*mphi)/rho)/radius
             Source_VC(Bz_,i,j,k) = Source_VC(Bz_,i,j,k) &
                  + (State_VGB(Bz_,i,j,k,iBlock) &
                  *   State_VGB(RhoUy_,i,j,k,iBlock) &
                  -  State_VGB(By_,i,j,k,iBlock) &
                  *   State_VGB(RhoUz_,i,j,k,iBlock))&
                  /State_VGB(Rho_,i,j,k,iBlock)/y_BLK(i,j,k,iBlock)
          end if
          if(UseB .and. UseB0)then
             ! Source[mr] = (B0.B1 - 2 B0phi * Bphi)/radius
             Source_VC(RhoUy_,i,j,k) = Source_VC(RhoUy_,i,j,k) &
                  + (sum(State_VGB(Bx_:Bz_,i,j,k,iBlock) &
                  *      B0_DGB(:,i,j,k,iBlock)) &
                  - 2.0*State_VGB(Bz_,i,j,k,iBlock)*B0_DGB(z_,i,j,k,iBlock)) &
                  / y_BLK(i,j,k,iBlock)

             ! Source[mphi] = (B0phi * Br + Bphi * B0r)/radius
             Source_VC(RhoUz_,i,j,k) = Source_VC(RhoUz_,i,j,k) &
                  + (B0_DGB(z_,i,j,k,iBlock)*State_VGB(By_,i,j,k,iBlock) &
                  +  B0_DGB(y_,i,j,k,iBlock)*State_VGB(Bz_,i,j,k,iBlock)) &
                  / y_BLK(i,j,k,iBlock)

             ! Source[Bphi]=((B0phi * mr - B0r * mphi)/rho)/radius
             Source_VC(Bz_,i,j,k) = Source_VC(Bz_,i,j,k) &
                  + (B0_DGB(z_,i,j,k,iBlock)*State_VGB(RhoUy_,i,j,k,iBlock) &
                  -  B0_DGB(y_,i,j,k,iBlock)*State_VGB(RhoUz_,i,j,k,iBlock))&
                  /State_VGB(Rho_,i,j,k,iBlock)/y_BLK(i,j,k,iBlock)
          end if
       end do; end do; end do

       ! For now, no Hall MHD implementation for rz-geometry
       if(UseB .and. UseBiermannBattery .and. &
            (UseElectronPressure .or. ElectronPressureRatio > 0.0 .or. &
            .not.UseIdealEos))then

          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.true_cell(i,j,k,iBlock)) CYCLE

             ! Source[Bphi] = [ 1/(q_e*n_e) * (dP_e/dZ) ] / radius
             Source_VC(Bz_,i,j,k) = Source_VC(Bz_,i,j,k) &
                  + IonMassPerCharge_G(i,j,k)/State_VGB(Rho_,i,j,k,iBlock) &
                  /y_Blk(i,j,k,iBlock) &
                  *0.5*(Pe_G(i+1,j,k) - Pe_G(i-1,j,k))/Dx_Blk(iBlock)
          end do; end do; end do
       end if
    end if


    ! We consider two cases: curl(B0) is zero analytically or non-zero
    ! These are distinguished by UseCurlB0 being true or failes.
    !
    ! Momentum equation has the Lorentz force J x B
    !     = (curl B1) x (B1 + B0)     if curl B0 = 0
    !     = (curl B1+B0) x (B1 + B0)  if curl B0 is not 0
    !
    ! Conservative fluxes add the divergence of the Maxwell tensor
    !     div(B1^2 + B1.B0 - B1 B1 - B1 B0 - B0 B1)
    ! 
    ! Deviations between these two are
    !   -B1 div(B1)       - usual div B source
    !   -B1 div(B0)       - div(B0) source
    !   -curl(B0) x B1    - remove this if curl B0 = 0
    !   +curl(B0) x B0    - add this if curl B0 is not 0

    if(UseB .and. UseDivbSource)then
       if(IsCartesian)then
          call calc_divb_source
       else
          call calc_divb_source_gencoord
       end if

       if(DoTestMe)write(*,*)'divb=',DivB1_GB(iTest,jTest,kTest,BlkTest)
       if(DoTestMe.and.VarTest>=RhoUx_.and.VarTest<=RhoUz_)&
            call write_source('After B0B1 source')

       ! Add contributions to other source terms
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(.not.true_cell(i,j,k,iBlock)) CYCLE

          RhoInv = 1.0/State_VGB(rho_,i,j,k,iBlock)
          Source_VC(Bx_:Bz_,i,j,k) = Source_VC(Bx_:Bz_,i,j,k) &
               -DivB1_GB(i,j,k,iBlock)* &
               State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock)*RhoInv

          if(.not. IsMhd) CYCLE

          ! -B1 div(B1)       - usual div B source

          Source_VC(rhoUx_:rhoUz_,i,j,k) = Source_VC(RhoUx_:RhoUz_,i,j,k) &
               - DivB1_GB(i,j,k,iBlock)*State_VGB(Bx_:Bz_,i,j,k,iBlock)
          Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) &
               - DivB1_GB(i,j,k,iBlock)* &
               sum(State_VGB(Bx_:Bz_,i,j,k,iBlock)*&
               State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock))*RhoInv
       end do; end do; end do

       if(DoTestMe)call write_source('After divb source')

       if (UseB0Source) then

          !   -B1 div(B0)     - div(B0) source
          ! -curl(B0) x B1    - remove this term (in case curl B0 should be 0) 
          !                     have to undo this if curl B0 is actually not 0

          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.true_cell(i,j,k,iBlock)) CYCLE
             Source_VC(RhoUx_:RhoUz_,i,j,k) = Source_VC(rhoUx_:rhoUz_,i,j,k) &
                  - State_VGB(Bx_:Bz_,i,j,k,iBlock)*DivB0_CB(i,j,k,iBlock) &
                  - cross_product( &
                  CurlB0_DCB(:,i,j,k,iBlock), State_VGB(Bx_:Bz_,i,j,k,iBlock))
          end do; end do; end do

          if(DoTestMe.and.VarTest>=RhoUx_.and.VarTest<=RhoUz_)then 
             write(*,*)'DivB0_CB  =',DivB0_CB(iTest,jTest,kTest,BlkTest)
             write(*,*)'CurlB0_DCB=',CurlB0_DCB(:,iTest,jTest,kTest,BlkTest)
             call write_source('After B0 source')
          end if
       end if
    else
       if(UseB)call calc_divb(iBlock)
    end if

    if(UseB .and. UseCurlB0)then

       !   +curl(B0) x B1    - undo source term above
       !   +curl(B0) x B0    - add this since curl B0 is not 0
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(.not.true_cell(i,j,k,iBlock)) CYCLE
          if(R_BLK(i,j,k,iBlock) < rCurrentFreeB0)CYCLE

          CurlB0CrossB_D = cross_product(&
               CurlB0_DCB(:,i,j,k,iBlock),&
               State_VGB(Bx_:Bz_,i,j,k,iBlock) + B0_DGB(:,i,j,k,iBlock))
          Source_VC(rhoUx_:rhoUz_,i,j,k) = Source_VC(rhoUx_:rhoUz_,i,j,k) &
               + CurlB0CrossB_D
          Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) &
               + sum(CurlB0CrossB_D*State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock))&
               /State_VGB(rho_,i,j,k,iBlock)
       end do; end do; end do

       if(DoTestMe .and. &
            (VarTest==Energy_.or.VarTest>=RhoUx_.and.VarTest<=RhoUz_))&
            call write_source('After curl B0')
    end if

    if(UseB .and. boris_correction &                   !^CFG IF BORISCORR BEGIN
         .and. boris_cLIGHT_factor < 0.9999 & 
         .and. index(test_string,'nodivE')<1) then

       Coef = (boris_cLIGHT_factor**2 - 1.0)*inv_c2LIGHT
       FullB_DC = State_VGB(Bx_:Bz_,1:nI,1:nJ,1:nK,iBlock)
       if(UseB0)FullB_DC = FullB_DC + B0_DGB(:,1:nI,1:nJ,1:nK,iBlock) 
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
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
    end if                                               !^CFG END BORISCORR

    ! These source terms apply to all the fluids
    do iFluid = 1, nFluid
       call select_fluid
       ! Add gravity and/or centrifugal force
       if(UseGravity .or. UseRotatingFrame) then

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
               (VarTest==Energy_ .or. VarTest>=iRhoUx.and.VarTest<=iRhoUz))then
             call write_source('After gravity')
          end if
       end if

       ! Add Coriolis forces
       if(UseRotatingFrame)then
          select case(TypeCoordSystem)
          case('HGC','HGR')
             ! This is a special case since Omega is parallel with the Z axis
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
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

    if(SignB_>1 .and. DoThinCurrentSheet)then
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(.not.true_cell(i,j,k,iBlock)) CYCLE

          ! Note that the velocity of the first (and only) fluid is used
          DivU            =        uDotArea_XI(i+1,j,k,1) -uDotArea_XI(i,j,k,1)
          if(nJ > 1) DivU = DivU + uDotArea_YI(i,j+1,k,1) -uDotArea_YI(i,j,k,1)
          if(nK > 1) DivU = DivU + uDotArea_ZI(i,j,k+1,1) -uDotArea_ZI(i,j,k,1)
          DivU = vInv_CB(i,j,k,iBlock)*DivU

          Source_VC(SignB_,i,j,k) = Source_VC(SignB_,i,j,k) &
               + State_VGB(SignB_,i,j,k,iBlock)*DivU
       end do; end do; end do
    end if

    if(UseUserSource)then
       call user_calc_sources
       if(DoTestMe) call write_source('After user sources')
    end if

    if(DoTestMe) call write_source('final')

  contains
    !==========================================================================
    subroutine calc_divb_source

      ! Variables needed for div B source terms
      real:: DxInvHalf, DyInvHalf, DzInvHalf, DivBInternal_C(1:nI,1:nJ,1:nK)
      real:: dB1nEast, dB1nWest, dB1nSouth, dB1nNorth, dB1nTop, dB1nBot
      !------------------------------------------------------------------------

      DxInvHalf = 0.5/Dx_BLK(iBlock)
      DyInvHalf = 0.5/Dy_BLK(iBlock)
      DzInvHalf = 0.5/Dz_BLK(iBlock)

      do k = 1, nK; do j = 1, nJ; do i = 1, nI
         if(.not.true_cell(i,j,k,iBlock)) CYCLE

         dB1nEast = DxInvHalf*&
              (RightState_VX(Bx_,i,j,k)-LeftState_VX(Bx_,i,j,k))

         dB1nWest = DxInvHalf*&
              (RightState_VX(Bx_,i+1,j,k)-LeftState_VX(Bx_,i+1,j,k))

         dB1nSouth = DyInvHalf* &
              (RightState_VY(By_,i,j,k)-LeftState_VY(By_,i,j,k))

         dB1nNorth = DyInvHalf* &
              (RightState_VY(By_,i,j+1,k)-LeftState_VY(By_,i,j+1,k))

         if(nK > 1)then
            dB1nBot = DzInvHalf * &
                 (RightState_VZ(Bz_,i,j,k)-LeftState_VZ(Bz_,i,j,k))

            dB1nTop = DzInvHalf * &
                 (RightState_VZ(Bz_,i,j,k+1)-LeftState_VZ(Bz_,i,j,k+1))
         end if

         DivBInternal_C(i,j,k) = 2*(&
              DxInvHalf*(LeftState_VX(Bx_,i+1,j,k) -RightState_VX(Bx_,i,j,k))+&
              DyInvHalf*(LeftState_VY(By_,i,j+1,k) -RightState_VY(By_,i,j,k)))

         if(nK > 1) DivBInternal_C(i,j,k) = DivBInternal_C(i,j,k) + &
              2*DzInvHalf*(LeftState_VZ(Bz_,i,j,k+1) -RightState_VZ(Bz_,i,j,k))

         DivB1_GB(i,j,k,iBlock)  = DivBInternal_C(i,j,k) + &
              dB1nEast + dB1nWest + dB1nSouth + dB1nNorth

         if(nK > 1) DivB1_GB(i,j,k,iBlock) = DivB1_GB(i,j,k,iBlock) &
              + dB1nTop + dB1nBot

         ! Momentum source term from B0 only needed for div(B^2/2 - BB) 
         ! discretization
         if(.not.(IsMhd.and.UseB0)) CYCLE

         Source_VC(RhoUx_:RhoUz_,i,j,k) = Source_VC(RhoUx_:RhoUz_,i,j,k) &
              -B0_DX(:,i,j,k)*dB1nEast    &
              -B0_DX(:,i+1,j,k)*dB1nWest  &
              -B0_DY(:,i,j,k)*dB1nSouth   &
              -B0_DY(:,i,j+1,k)*dB1nNorth

         if(nK > 1) &
              Source_VC(RhoUx_:RhoUz_,i,j,k) = Source_VC(RhoUx_:RhoUz_,i,j,k) &
              -B0_DZ(:,i,j,k)*dB1nBot     &
              -B0_DZ(:,i,j,k+1)*dB1nTop

      end do; end do; end do

      ! Momentum source term from B0 only needed for true MHD equations
      if(.not.(IsMhd .and. UseB0)) RETURN

      do k = 1, nK; do j = 1, nJ; do i = 1, nI
         if(.not.true_cell(i,j,k,iBlock)) CYCLE
         Source_VC(rhoUx_:rhoUz_,i,j,k) = Source_VC(rhoUx_:rhoUz_,i,j,k) &
              - DivBInternal_C(i,j,k)*B0_DGB(:,i,j,k,iBlock)
      end do; end do; end do

    end subroutine calc_divb_source
    !==========================================================================
    subroutine calc_divb_source_gencoord

      use BATL_lib, ONLY: FaceNormal_DDFB

      real :: FaceArea_D(3), vInvHalf
      real :: B1nJumpL, B1nJumpR, DivBInternal_C(1:nI,1:nJ,1:nK)
      integer :: i, j, k

      character(len=*), parameter:: NameSub = 'calc_divb_source_gencoord'
      !------------------------------------------------------------------------

      do k = 1, nK; do j = 1, nJ; do i = 1, nI
         if(.not.true_cell(i,j,k,iBlock)) CYCLE

         VInvHalf=vInv_CB(i,j,k,iBlock)*0.5
         FaceArea_D = FaceNormal_DDFB(:,1,i,j,k,iBlock)
         B1nJumpL =VInvHalf*&
              sum(FaceArea_D*(RightState_VX(Bx_:Bz_,i,j,k) &
              -               LeftState_VX(Bx_:Bz_,i,j,k)))
         DivBInternal_C(i,j,k) = &
              -sum(FaceArea_D*RightState_VX(Bx_:Bz_,i,j,k))

         FaceArea_D = FaceNormal_DDFB(:,1,i+1,j,k,iBlock)
         B1nJumpR =  VInvHalf*&
              sum(FaceArea_D*(RightState_VX(Bx_:Bz_,i+1,j,k) &
              -               LeftState_VX(Bx_:Bz_,i+1,j,k)))

         DivBInternal_C(i,j,k) = DivBInternal_C(i,j,k) &
              + sum(FaceArea_D*LeftState_VX(Bx_:Bz_,i+1,j,k))

         DivB1_GB(i,j,k,iBlock)  = B1nJumpL + B1nJumpR

         if(.not.IsMhd) CYCLE

         Source_VC(RhoUx_:RhoUz_,i,j,k) = Source_VC(RhoUx_:RhoUz_,i,j,k) &
              - B0_DX(:,i,j,k)*B1nJumpL   &
              - B0_DX(:,i+1,j,k)*B1nJumpR

      end do; end do; end do

      if(DoTestMe)write(*,*)NameSub,' after i divbint, divb1=', &
           DivBInternal_C(iTest,jTest,kTest), &
           DivB1_GB(iTest,jTest,kTest,BlkTest)

      do k = 1, nK; do j = 1, nJ; do i = 1, nI 
         if(.not.true_cell(i,j,k,iBlock)) CYCLE
         VInvHalf=vInv_CB(i,j,k,iBlock)*0.5
         FaceArea_D = FaceNormal_DDFB(:,2,i,j,k,iBlock)
         B1nJumpL = VInvHalf*&
              sum(FaceArea_D*(RightState_VY(Bx_:Bz_,i,j,k) &
              -               LeftState_VY(Bx_:Bz_,i,j,k)))
         DivBInternal_C(i,j,k) = DivBInternal_C(i,j,k) &
              - sum(FaceArea_D*RightState_VY(Bx_:Bz_,i,j,k))

         FaceArea_D =  FaceNormal_DDFB(:,2,i,j+1,k,iBlock)
         B1nJumpR = VInvHalf*&
              sum(FaceArea_D*(RightState_VY(Bx_:Bz_,i,j+1,k) &
              -               LeftState_VY(Bx_:Bz_,i,j+1,k)))

         DivBInternal_C(i,j,k) = DivBInternal_C(i,j,k) &
              + sum(FaceArea_D*LeftState_VY(Bx_:Bz_,i,j+1,k))

         DivB1_GB(i,j,k,iBlock)  = DivB1_GB(i,j,k,iBlock) &
              + B1nJumpL + B1nJumpR

         if(.not.IsMhd) CYCLE

         Source_VC(RhoUx_:RhoUz_,i,j,k) = Source_VC(RhoUx_:RhoUz_,i,j,k)&
              -B0_DY(:,i,j,k)*B1nJumpL &
              -B0_DY(:,i,j+1,k)*B1nJumpR

      end do; end do; end do

      if(DoTestMe)write(*,*)NameSub,' after j divbint, divb1=', &
           DivBInternal_C(iTest,jTest,kTest), &
           DivB1_GB(iTest,jTest,kTest,BlkTest)

      if(nK > 1)then
         do k = 1, nK; do j = 1, nJ; do i = 1, nI 
            if(.not.true_cell(i,j,k,iBlock)) CYCLE

            VInvHalf = vInv_CB(i,j,k,iBlock)*0.5
            FaceArea_D = FaceNormal_DDFB(:,3,i,j,k,iBlock)
            B1nJumpL = VInvHalf*&
                 sum(FaceArea_D*(RightState_VZ(Bx_:Bz_,i,j,k) &
                 -                LeftState_VZ(Bx_:Bz_,i,j,k)))

            DivBInternal_C(i,j,k) = DivBInternal_C(i,j,k) &
                 - sum(FaceArea_D*RightState_VZ(Bx_:Bz_,i,j,k))

            FaceArea_D = FaceNormal_DDFB(:,3,i,j,k+1,iBlock)
            B1nJumpR = VInvHalf*&
                 sum(FaceArea_D*(RightState_VZ(Bx_:Bz_,i,j,k+1) &
                 -               LeftState_VZ(Bx_:Bz_,i,j,k+1)))

            DivBInternal_C(i,j,k) = (DivBInternal_C(i,j,k) + &
                 sum(FaceArea_D*LeftState_VZ(Bx_:Bz_,i,j,k+1))) &
                 *vInv_CB(i,j,k,iBlock)

            DivB1_GB(i,j,k,iBlock)  = DivB1_GB(i,j,k,iBlock) &
                 + B1nJumpL + B1nJumpR

            if(.not.IsMhd) CYCLE

            Source_VC(rhoUx_:rhoUz_,i,j,k) = Source_VC(rhoUx_:rhoUz_,i,j,k)&
                 -B0_DZ(:,i,j,k)*B1nJumpL &
                 -B0_DZ(:,i,j,k+1)*B1nJumpR
         end do; end do; end do
      end if

      if(DoTestMe)write(*,*)NameSub,' after k divbint, divb1=', &
           DivBInternal_C(iTest,jTest,kTest), &
           DivB1_GB(iTest,jTest,kTest,BlkTest)

      do k = 1, nK; do j = 1, nJ; do i = 1, nI 
         if(.not.true_cell(i,j,k,iBlock)) CYCLE
         DivB1_GB(i,j,k,iBlock) = DivB1_GB(i,j,k,iBlock) +DivBInternal_C(i,j,k)
      end do; end do; end do

      if(DoTestMe)write(*,*)NameSub,' final divb1=', &
           DivB1_GB(iTest,jTest,kTest,BlkTest)

      if((.not.IsMhd).or.(.not.UseB0))RETURN

      do k = 1, nK; do j = 1, nJ; do i = 1, nI 
         if(.not.true_cell(i,j,k,iBlock)) CYCLE
         Source_VC(RhoUx_:RhoUz_,i,j,k) = Source_VC(RhoUx_:RhoUz_,i,j,k) &
              - DivBInternal_C(i,j,k)*B0_DGB(:,i,j,k,iBlock)            
      end do; end do; end do

    end subroutine calc_divb_source_gencoord
    !==========================================================================

    subroutine write_source(String)
      character(len=*), intent(in) :: String
      write(*,'(a,es13.5)') NameSub//": "//String//" S(VarTest)=",&
           Source_VC(VarTest,iTest,jTest,kTest) 
    end subroutine write_source

  end subroutine calc_sources

  !============================================================================

  subroutine calc_divb(iBlock)

    ! Calculate div B for a block and store result into DivB1_GB
    ! Compute divB using averaged and conservatively corrected 
    ! left and right values

    use ModMain,       ONLY: nI, nJ, nK
    use ModVarIndexes, ONLY: Bx_, By_, Bz_
    use ModGeometry,   ONLY: dx_BLK, dy_BLK, dz_BLK 
    use ModAdvance,    ONLY: DivB1_GB, &
         LeftState_VX, RightState_VX, &
         LeftState_VY, RightState_VY, &
         LeftState_VZ, RightState_VZ

    integer, intent(in) :: iBlock

    integer:: i, j, k
    real   :: DivB, InvDx, InvDy, InvDz
    !--------------------------------------------------------------------------

    InvDx            = 1/dx_BLK(iBlock)
    if(nJ > 1) InvDy = 1/dy_BLK(iBlock)
    if(nK > 1) InvDz = 1/dz_BLK(iBlock)

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
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

  !============================================================================

  subroutine get_tesi_c(iBlock, TeSi_C)

    use ModAdvance,    ONLY: UseElectronPressure, UseIdealEos
    use ModAdvance,    ONLY: State_VGB, p_, Pe_, Rho_
    use ModSize,       ONLY: nI, nJ, nK
    use ModPhysics,    ONLY: No2Si_V, UnitTemperature_, &
         AverageIonCharge, PePerPtotal
    use ModMultifluid, ONLY: MassIon_I
    use ModUser,       ONLY: user_material_properties

    integer, intent(in)  :: iBlock
    real,    intent(out) :: TeSi_C(1:nI, 1:nJ, 1:nK)

    integer:: i, j, k
    !--------------------------------------------------------------------------
    if(UseIdealEos)then
       if(UseElectronPressure)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             TeSi_C(i,j,k) = State_VGB(Pe_,i,j,k,iBlock) &
                  /State_VGB(Rho_,i,j,k,iBlock)
          end do; end do; end do
          TeSi_C = TeSi_C * No2Si_V(UnitTemperature_ ) * &
               MassIon_I(1)/AverageIonCharge
       else
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             TeSi_C(i,j,k) = State_VGB(p_,i,j,k,iBlock) &
                  /State_VGB(Rho_,i,j,k,iBlock)
          end do; end do; end do
          TeSi_C = TeSi_C * No2Si_V(UnitTemperature_ ) * &
               MassIon_I(1)/AverageIonCharge * PePerPtotal
       end if
    else
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          call user_material_properties( &
               State_VGB(:,i,j,k,iBlock), TeOut=TeSi_C(i,j,k))
       end do; end do; end do
    end if
  end subroutine get_tesi_c

end module ModCalcSource
