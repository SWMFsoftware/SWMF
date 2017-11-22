!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModUpdateState

  use BATL_lib, ONLY: &
       test_start, test_stop, iTest, jTest, kTest, iBlockTest, &
       iProcTest, iVarTest

  implicit none

  private ! except

  public:: update_state         ! call user_update_state or update_state_normal
  public:: update_state_normal  ! normal update of state variables
  public:: update_b0            ! update time varying B0 field
  public:: update_te0           ! update Te0 variable
  public:: update_check         ! check and correct update if necessary
  public:: fix_anisotropy       ! fix pressure anisotropy after update
  public:: select_conservative  ! select cells to use conservative update

contains
  !============================================================================

  subroutine update_state(iBlock)

    use ModProcMH
    use ModMain
    use ModAdvance
    use ModMultiFluid, ONLY: iFluid, nFluid
    use BATL_lib, ONLY: CellVolume_GB
    use ModHeatFluxCollisionless, ONLY: UseHeatFluxCollisionless, &
         update_heatflux_collisionless
    use ModUserInterface ! user_update_states
    use ModMessagePass, ONLY: fix_buffer_grid

    integer, intent(in) :: iBlock
    integer :: iVar

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'update_state'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    if(DoTest)then
       write(*,*)NameSub,' dt=',time_BLK(iTest,jTest,kTest,iBlock)*Cfl
       if(allocated(IsConserv_CB)) write(*,*)NameSub,' IsConserv=', &
            IsConserv_CB(iTest,jTest,kTest,iBlock)
       write(*,*)
       do iVar=1,nVar
          write(*,'(2x,2a,es23.15)')NameVar_V(iVar), '(TestCell)  =',&
               State_VGB(iVar,iTest,jTest,kTest,iBlockTest)
       end do
       do iFluid = 1, nFluid
          write(*,'(2x,a,i2,a,es23.15)') &
               'E(',iFluid,')=',Energy_GBI(iTest,jTest,kTest,iBlockTest,iFluid)
       end do
       write(*,*)'Fluxes and sources for ',NameVar_V(iVarTest)
       write(*,'(2x,a,2es23.15)') &
            'X fluxes L,R =',Flux_VX(iVarTest,iTest,jTest,kTest),&
            Flux_VX(iVarTest,iTest+1,jTest,kTest)
       write(*,'(2x,a,2es23.15)') &
            'Y fluxes L,R =',Flux_VY(iVarTest,iTest,jTest,kTest),&
            Flux_VY(iVarTest,iTest,jTest+1,kTest)
       write(*,'(2x,a,2es23.15)') &
            'Z fluxes L,R =',Flux_VZ(iVarTest,iTest,jTest,kTest),&
            Flux_VZ(iVarTest,iTest,jTest,kTest+1)
       write(*,'(2x,a,es23.15)')'source=',Source_VC(iVarTest,iTest,jTest,kTest)
       write(*,'(2x,a,es23.15)')'fluxes=',(Flux_VX(iVarTest,iTest,jTest,kTest) &
            -Flux_VX(iVarTest,iTest+1,jTest,kTest)                        &
            +Flux_VY(iVarTest,iTest,jTest,kTest)                          &
            -Flux_VY(iVarTest,iTest,jTest+1,kTest)                        &
            +Flux_VZ(iVarTest,iTest,jTest,kTest)                          &
            -Flux_VZ(iVarTest,iTest,jTest,kTest+1) )                      &
            /CellVolume_GB(iTest,jTest,kTest,iBlockTest)
    end if

    ! Note must copy state to old state only if iStage is 1.
    if(iStage==1) then
       StateOld_VGB(:,:,:,:,iBlock) = State_VGB(:,:,:,:,iBlock)
       EnergyOld_CBI(:,:,:,iBlock,:) = Energy_GBI(1:nI,1:nJ,1:nK,iBlock,:)
    end if

    if(UseUserUpdateStates)then
       call user_update_states(iBlock)
    else
       call update_state_normal(iBlock)
    end if

    if(Ehot_ > 1 .and. UseHeatFluxCollisionless) then
       call update_heatflux_collisionless(iBlock)
       if(UseBufferGrid) call fix_buffer_grid(iBlock)
    end if

    if(DoTest)then
       write(*,*)NameSub,' final:'
       do iVar=1,nVar
          write(*,'(2x,2a,es23.15)')NameVar_V(iVar),'(TestCell)  =',&
               State_VGB(iVar,iTest,jTest,kTest,iBlockTest)
       end do
       do iFluid = 1, nFluid
          write(*,'(2x,a,i2,a,es23.15)') &
               'E(',iFluid,')=',Energy_GBI(iTest,jTest,kTest,iBlockTest,iFluid)
       end do
    end if

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine update_state
  !============================================================================

  subroutine update_state_normal(iBlock)

    use ModProcMH
    use ModMain
    use ModAdvance
    use ModB0, ONLY: B0_DGB
    use ModPhysics
    use ModGeometry, ONLY: true_cell
    use ModSemiImplVar, ONLY: UseStableImplicit
    use ModVarIndexes, ONLY: pe_, p_
    use ModPointImplicit, ONLY: UsePointImplicit, UsePointImplicit_B, &
         IsDynamicPointImplicit, update_point_implicit
    use ModMultiIon, ONLY: multi_ion_source_impl, multi_ion_init_point_impl, &
         multi_ion_set_restrict, multi_ion_update, DoRestrictMultiIon
    use ModEnergy
    use ModWaves, ONLY: nWave, WaveFirst_,WaveLast_, &
         UseWavePressure, UseWavePressureLtd, UseAlfvenWaves, DoAdvectWaves, &
         update_wave_group_advection
    use ModResistivity, ONLY: UseResistivity, UseResistiveFlux, &
         calc_resistivity_source
    use ModFaceValue, ONLY: UseFaceIntegral4
    use BATL_lib, ONLY: CellVolume_GB
    use ModUserInterface ! user_calc_sources, user_init_point_implicit
    use ModMessagePass, ONLY: fix_buffer_grid
    use ModIonElectron, ONLY: ion_electron_source_impl, &
         ion_electron_init_point_impl, HypEDecay

    integer, intent(in) :: iBlock

    integer :: i, j, k

    ! These variables have to be double precision for accurate Boris scheme
    real:: DtLocal, DtFactor, Coeff

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'update_state_normal'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    ! Nothing to do if time step is zero
    if(time_accurate .and. Dt == 0.0) RETURN

    ! Add Joule heating: dPe/dt or dP/dt += (gamma-1)*eta*j**2
    ! also dE/dt += eta*j**2 for semi-implicit scheme (UseResistiveFlux=F)
    ! and heat exchange between electrons and ions (mult-ion is not coded).

    if(.not.UseMultiIon .and. UseResistivity .and. &
         (UseElectronPressure .or. UseNonConservative .or. &
         .not.UseResistiveFlux)) then

       call calc_resistivity_source(iBlock)
       if(DoTest)write(*,'(2x,2a,15es20.12)') &
            NameSub, ' after add_resistive_source          =', &
            State_VGB(iVarTest,iTest,jTest,kTest,iBlock),       &
            Energy_GBI(iTest,jTest,kTest,iBlock,:)
    end if

    ! Calculate partial step size compared to largest stable time step
    if(nStage==4)then
       ! Classical 4th order Runge-Kutta scheme
       select case(iStage)
       case(1)
          DtFactor = Cfl/2
       case(2)
          DtFactor = Cfl/2
       case(3)
          DtFactor = Cfl
       case(4)
          DtFactor = Cfl/6
       end select
    elseif(UseHalfStep)then
       DtFactor = (Cfl*iStage)/nStage
    else
       DtFactor = Cfl
    end if

    ! Modify electron pressure source term to electron entropy if necessary
    ! d(Se)/d(Pe) = Pe^(1/gammaE-1)/gammaE
    if(UseElectronPressure .and. UseElectronEntropy)then
       do k = 1,nK; do j = 1,nJ; do i = 1,nI
          Source_VC(Pe_,i,j,k) = Source_VC(Pe_,i,j,k)/GammaElectron &
               * State_VGB(Pe_,i,j,k,iBlock)**(1/GammaElectron - 1)
       end do; end do; end do
    end if
    do k = 1,nK; do j = 1,nJ; do i = 1,nI
       DtLocal = DtFactor*time_BLK(i,j,k,iBlock)
       Source_VC(:,i,j,k) = &
            DtLocal* (Source_VC(:,i,j,k) + &
            ( Flux_VX(:,i,j,k) - Flux_VX(:,i+1,j,k) &
            + Flux_VY(:,i,j,k) - Flux_VY(:,i,j+1,k) &
            + Flux_VZ(:,i,j,k) - Flux_VZ(:,i,j,k+1) ) &
            /CellVolume_GB(i,j,k,iBlock) )
    end do; end do; end do

    if(nOrder == 4 .and. UseFaceIntegral4 .and. nDim > 1)then
       ! Integrate fluxes in the transverse direction (eq. 20)
       ! <F> = F + Laplace_transverse(F)/24
       do k = 1,nK; do j = 1,nJ; do i = 1,nI
          Coeff = DtFactor*time_BLK(i,j,k,iBlock) &
               /  (24.0*CellVolume_GB(i,j,k,iBlock))
          ! Add f
          Source_VC(:,i,j,k) = Source_VC(:,i,j,k) + Coeff* &
               ( Flux_VX(:,i,j+1,k) &
               + Flux_VX(:,i,j-1,k) &
               - 2*(nDim-1)*Flux_VX(:,i,j,k) &
               - Flux_VX(:,i+1,j+1,k) &
               - Flux_VX(:,i+1,j-1,k) &
               + 2*(nDim-1)*Flux_VX(:,i+1,j,k) &
               + Flux_VY(:,i+1,j,k) &
               + Flux_VY(:,i-1,j,k) &
               - 2*(nDim-1)*Flux_VY(:,i,j,k) &
               - Flux_VY(:,i+1,j+1,k) &
               - Flux_VY(:,i-1,j+1,k) &
               + 2*(nDim-1)*Flux_VY(:,i,j+1,k) )
          if(nK == 1) CYCLE
          ! Remaining terms for 3D
          Source_VC(:,i,j,k) = Source_VC(:,i,j,k) + Coeff* &
               ( Flux_VX(:,i,j,k+1) &
               + Flux_VX(:,i,j,k-1) &
               - Flux_VX(:,i+1,j,k+1) &
               - Flux_VX(:,i+1,j,k-1) &
               + Flux_VY(:,i,j,k+1) &
               + Flux_VY(:,i,j,k-1) &
               - Flux_VY(:,i,j+1,k+1) &
               - Flux_VY(:,i,j+1,k-1) &
               + Flux_VZ(:,i+1,j,k) &
               + Flux_VZ(:,i-1,j,k) &
               + Flux_VZ(:,i,j+1,k) &
               + Flux_VZ(:,i,j-1,k) &
               - 4*Flux_VZ(:,i,j,k) &
               - Flux_VZ(:,i+1,j,k+1) &
               - Flux_VZ(:,i-1,j,k+1) &
               - Flux_VZ(:,i,j+1,k+1) &
               - Flux_VZ(:,i,j-1,k+1) &
               + 4*Flux_VZ(:,i,j,k+1) )

       end do; end do; end do
    end if

    if(UseMultiIon .and. DoRestrictMultiIon)call multi_ion_set_restrict(iBlock)

    if(DoTest)write(*,'(2x,2a,15es20.12)') &
         NameSub, ' original testvar and energy         =', &
         State_VGB(iVarTest,iTest,jTest,kTest,iBlock),       &
         Energy_GBI(iTest,jTest,kTest,iBlock,:)

    call update_explicit

    if(UseMultiIon .and. &
         (UseSingleIonVelocity .or. UseSingleIonTemperature)) then

       call fix_multi_ion_update(iBlock)
       call calc_energy_cell(iBlock)

       if(DoTest)write(*,'(2x,2a,15es20.12)') &
            NameSub, ' after fix multiion update           =', &
            State_VGB(iVarTest,iTest,jTest,kTest,iBlock),       &
            Energy_GBI(iTest,jTest,kTest,iBlock,:)

    end if

    if(UseMultiIon .and. IsMhd)then
       call multi_ion_update(iBlock, IsFinal = .false.)

       if(DoTest)write(*,'(2x,2a,15es20.12)')  &
            NameSub, ' after multiion update1              =', &
            State_VGB(iVarTest,iTest,jTest,kTest,iBlock),       &
            Energy_GBI(iTest,jTest,kTest,iBlock,:)

    end if

    ! The point implicit update and other stuff below are only done in
    ! the last stage except for ion-electron equations where the point
    ! implicit has to be done in every stage.
    if(.not.UseEfield .and. iStage < nStage)then
       if(UseBufferGrid) call fix_buffer_grid(iBlock)
       RETURN
    end if

    ! Add point implicit user or multi-ion source terms
    if (UsePointImplicit)then
       if(IsDynamicPointImplicit .or. UsePointImplicit_B(iBlock)) then
          if(UseEfield)then
             call update_point_implicit(iBlock, ion_electron_source_impl, &
                  ion_electron_init_point_impl)
          elseif(UseMultiIon .and. .not.UseSingleIonVelocity)then
             call update_point_implicit(iBlock, multi_ion_source_impl, &
                  multi_ion_init_point_impl)
          elseif(UseUserSource) then
             call update_point_implicit(iBlock, user_calc_sources, &
                  user_init_point_implicit)
          end if

          ! Make ion temperatures equal if requested
          if(UseMultiIon .and. &
               (UseSingleIonVelocity .or. UseSingleIonTemperature)) &
               call fix_multi_ion_update(iBlock)

          ! Make sure that energy is consistent
          call calc_energy_cell(iBlock)

          if(DoTest)write(*,'(2x,2a,15es20.12)') &
               NameSub, ' after point impl state              =', &
               State_VGB(iVarTest, iTest,jTest,kTest,iBlock),      &
               Energy_GBI(iTest,jTest,kTest,iBlock,:)
       end if
    end if

    if(UseMultiIon .and. IsMhd)then
       call multi_ion_update(iBlock, IsFinal = .true.)
       if(DoTest)write(*,'(2x,2a,15es20.12)') &
            NameSub, ' after multiion update2              =', &
            State_VGB(iVarTest,iTest,jTest,kTest,iBlock),       &
            Energy_GBI(iTest,jTest,kTest,iBlock,:)
    end if

    ! The parabolic div B decay is only done in the last stage.
    if(UseHyperbolicDivb .and. HypDecay > 0 .and. iStage == nStage) &
         State_VGB(Hyp_,1:nI,1:nJ,1:nK,iBlock) = &
         State_VGB(Hyp_,1:nI,1:nJ,1:nK,iBlock)*(1 - HypDecay)

    if(UseEfield .and. HypEDecay > 0 .and. iStage == nStage) &
         State_VGB(HypE_,1:nI,1:nJ,1:nK,iBlock) = &
         State_VGB(HypE_,1:nI,1:nJ,1:nK,iBlock)*(1 - HypEDecay)

    if(UseStableImplicit) call deduct_expl_source

    if(UseBufferGrid) call fix_buffer_grid(iBlock)

    if(DoTest)write(*,'(2x,2a,15es20.12)') &
         NameSub, ' final state                         =', &
         State_VGB(iVarTest,iTest,jTest,kTest,iBlock),       &
         Energy_GBI(iTest,jTest,kTest,iBlock,:)

    call test_stop(NameSub, DoTest, iBlock)
  contains
    !==========================================================================
    subroutine update_explicit

      use ModBorisCorrection, ONLY: mhd_to_boris, boris_to_mhd, &
           mhd_to_boris_simple, boris_simple_to_mhd

      ! Allocatable storage for classical 4th order Runge-Kutta scheme
      real, allocatable, save:: Rk4_VCB(:,:,:,:,:), Rk4_CBI(:,:,:,:,:)

      real, parameter:: cThird = 1./3.
      real:: Coeff1, Coeff2
      integer:: iFluid, iRho
      !------------------------------------------------------------------------

      if(boris_correction) then
         ! Convert StateOld_VGB and State_VGB from classical MHD variables
         ! to semi-relativistic MHD variable
         do k=1,nK; do j=1,nJ; do i=1,nI
            if(.not.true_cell(i,j,k,iBlock)) CYCLE

            if(UseB0)then
               call mhd_to_boris(StateOld_VGB(:,i,j,k,iBlock), &
                    EnergyOld_CBI(i,j,k,iBlock,1), B0_DGB(:,i,j,k,iBlock))
               ! State_VGB is not used in 1-stage and HalfStep schemes
               if(.not.UseHalfStep .and. nStage > 1) &
                    call mhd_to_boris(State_VGB(:,i,j,k,iBlock), &
                    Energy_GBI(i,j,k,iBlock,1), B0_DGB(:,i,j,k,iBlock))
            else
               call mhd_to_boris(StateOld_VGB(:,i,j,k,iBlock), &
                    EnergyOld_CBI(i,j,k,iBlock,1))
               ! State_VGB is not used in 1-stage and HalfStep schemes
               if(.not.UseHalfStep .and. nStage > 1) &
                    call mhd_to_boris(State_VGB(:,i,j,k,iBlock), &
                    Energy_GBI(i,j,k,iBlock,1))
            end if
         end do; end do; end do

         if(DoTest)write(*,'(2x,2a,15es20.12)') &
              NameSub, ' after mhd_to_boris                  =', &
              State_VGB(iVarTest,iTest,jTest,kTest,iBlock),       &
              Energy_GBI(iTest,jTest,kTest,iBlock,:)
      endif

      if(UseBorisSimple .and. IsMhd) then
         ! Update using simplified Boris correction, i.e. update
         !
         !    RhoUBorisSimple = (1+B^2/(rho*c^2)) * RhoU
         !
         ! instead of RhoU. See Gombosi et al JCP 2002, 177, 176 (eq. 38-39)
         do k=1,nK; do j=1,nJ; do i=1,nI
            if(.not.true_cell(i,j,k,iBlock)) CYCLE

            if(UseB0)then
               call mhd_to_boris_simple(StateOld_VGB(:,i,j,k,iBlock), &
                    B0_DGB(:,i,j,k,iBlock))
               ! State_VGB is not used in 1-stage and HalfStep schemes
               if(.not.UseHalfStep .and. nStage > 1) &
                    call mhd_to_boris_simple(State_VGB(:,i,j,k,iBlock), &
                    B0_DGB(:,i,j,k,iBlock))
            else
               call mhd_to_boris_simple(StateOld_VGB(:,i,j,k,iBlock))
               ! State_VGB is not used in 1-stage and HalfStep schemes
               if(.not.UseHalfStep .and. nStage > 1) &
                    call mhd_to_boris_simple(State_VGB(:,i,j,k,iBlock))
            end if
         end do; end do; end do

         if(DoTest)write(*,'(2x,2a,15es20.12)') &
              NameSub, ' after mhd_to_boris_simple           =', &
              State_VGB(iVarTest,iTest,jTest,kTest,iBlock),       &
              Energy_GBI(iTest,jTest,kTest,iBlock,:)
      end if

      if(UseElectronPressure .and. UseElectronEntropy)then
         ! Convert electron pressure to entropy
         ! Se = Pe^(1/GammaE)
         do k=1,nK; do j=1,nJ; do i=1,nI
            if(.not.true_cell(i,j,k,iBlock)) CYCLE

            StateOld_VGB(Pe_,i,j,k,iBlock) = &
                 StateOld_VGB(Pe_,i,j,k,iBlock)**(1/GammaElectron)
            ! State_VGB is not used in 1-stage and HalfStep schemes
            if(.not.UseHalfStep .and. nStage > 1) &
                 State_VGB(Pe_,i,j,k,iBlock) = &
                 State_VGB(Pe_,i,j,k,iBlock)**(1/GammaElectron)
         end do; end do; end do
      end if

      ! Now update State_VGB

      if(UseHalfStep .or. nStage == 1 .or. nStage == 4)then
         ! Update state variables starting from level n (=old) state
         do k=1,nK; do j=1,nJ; do i=1,nI
            State_VGB(:,i,j,k,iBlock) = &
                 StateOld_VGB(:,i,j,k,iBlock) + Source_VC(1:nVar,i,j,k)
         end do; end do; end do

         ! Update energy variables
         do iFluid = 1, nFluid; do k=1,nK; do j=1,nJ; do i=1,nI
            Energy_GBI(i,j,k,iBlock,iFluid) = &
                 EnergyOld_CBI(i,j,k,iBlock,iFluid) &
                 + Source_VC(nVar+iFluid,i,j,k)
         end do; end do; end do; end do
      else
         ! Update state variables starting from previous stage
         do k=1,nK; do j=1,nJ; do i=1,nI
            State_VGB(:,i,j,k,iBlock) = &
                 State_VGB(:,i,j,k,iBlock) + Source_VC(1:nVar,i,j,k)
         end do; end do; end do

         ! Update energy variables
         do iFluid = 1, nFluid; do k=1,nK; do j=1,nJ; do i=1,nI
            Energy_GBI(i,j,k,iBlock,iFluid) = &
                 Energy_GBI(i,j,k,iBlock,iFluid) + Source_VC(nVar+iFluid,i,j,k)
         end do; end do; end do; end do
      end if

      if(nStage == 4)then
         ! Classical 4th order Runge-Kutta scheme. Requires extra storage.
         if(.not.allocated(Rk4_VCB)) allocate( &
              Rk4_VCB(nVar,nI,nJ,nK,MaxBlock), &
              Rk4_CBI(nI,nJ,nK,MaxBlock,nFluid))

         select case(iStage)
         case(1)
            ! Rk4 = U1 = Un + Dt/2*Rn
            do k=1,nK; do j=1,nJ; do i=1,nI
               Rk4_VCB(:,i,j,k,iBlock) = State_VGB(:,i,j,k,iBlock)
            end do; end do; end do
            do iFluid = 1, nFluid; do k=1,nK; do j=1,nJ; do i=1,nI
               Rk4_CBI(i,j,k,iBlock,iFluid) = Energy_GBI(i,j,k,iBlock,iFluid)
            end do; end do; end do; end do
         case(2)
            ! U2 = Un + Dt/2*R1
            ! Rk4 = Rk4 + 2*U2 = 3*Un + Dt/2*Rn + Dt*R1
            do k=1,nK; do j=1,nJ; do i=1,nI
               Rk4_VCB(:,i,j,k,iBlock) = Rk4_VCB(:,i,j,k,iBlock) &
                    + 2*State_VGB(:,i,j,k,iBlock)
            end do; end do; end do
            do iFluid = 1, nFluid; do k=1,nK; do j=1,nJ; do i=1,nI
               Rk4_CBI(i,j,k,iBlock,iFluid) = Rk4_CBI(i,j,k,iBlock,iFluid) &
                    + 2*Energy_GBI(i,j,k,iBlock,iFluid)
            end do; end do; end do; end do
         case(3)
            ! U3 = Un + Dt*R2
            ! Rk4 = Rk4 + U3 - 4Un = Dt/2*Rn + Dt*R1 + Dt*R2
            do k=1,nK; do j=1,nJ; do i=1,nI
               Rk4_VCB(:,i,j,k,iBlock) = Rk4_VCB(:,i,j,k,iBlock) &
                    + State_VGB(:,i,j,k,iBlock) &
                    - 4*StateOld_VGB(:,i,j,k,iBlock)
            end do; end do; end do
            do iFluid = 1, nFluid; do k=1,nK; do j=1,nJ; do i=1,nI
               Rk4_CBI(i,j,k,iBlock,iFluid) = Rk4_CBI(i,j,k,iBlock,iFluid) &
                    + Energy_GBI(i,j,k,iBlock,iFluid) &
                    - 4*EnergyOld_CBI(i,j,k,iBlock,iFluid)
            end do; end do; end do; end do
         case(4)
            ! U4 = Un + Dt/6*R3
            ! Un+1 = U4 + Rk4/3 = Un + Dt/6*(Rn + 2*R1 + 2*R2 + R3)
            do k=1,nK; do j=1,nJ; do i=1,nI
               State_VGB(:,i,j,k,iBlock) = &
                    + State_VGB(:,i,j,k,iBlock) &
                    + cThird*Rk4_VCB(:,i,j,k,iBlock)
            end do; end do; end do
            do iFluid = 1, nFluid; do k=1,nK; do j=1,nJ; do i=1,nI
               Energy_GBI(i,j,k,iBlock,iFluid) = &
                    + Energy_GBI(i,j,k,iBlock,iFluid) &
                    + cThird*Rk4_CBI(i,j,k,iBlock,iFluid)
            end do; end do; end do; end do
         end select

      elseif(.not.UseHalfStep .and. iStage > 1)then
         ! Interpolation step for 2nd and 3rd order Runge-Kutta schemes

         ! Runge-Kutta scheme coefficients
         if (nStage==2) then
            Coeff1 = 0.5
         elseif (nStage==3) then
            if (iStage==2) then
               Coeff1 = 0.75
            elseif (iStage==3) then
               Coeff1 = 1./3.
            end if
         end if
         Coeff2 = 1 - Coeff1

         ! Interpolate state variables
         do k=1,nK; do j=1,nJ; do i=1,nI
            State_VGB(:,i,j,k,iBlock) = &
                 Coeff1*StateOld_VGB(:,i,j,k,iBlock) + &
                 Coeff2*State_VGB(:,i,j,k,iBlock)
         end do; end do; end do

         ! Interpolate energies
         do iFluid = 1, nFluid; do k=1,nK; do j=1,nJ; do i=1,nI
            Energy_GBI(i,j,k,iBlock,iFluid) = &
                 Coeff1*EnergyOld_CBI(i,j,k,iBlock,iFluid) + &
                 Coeff2*Energy_GBI(i,j,k,iBlock,iFluid)
         end do; end do; end do; end do

      endif

      if(DoTest)write(*,'(2x,2a,15es20.12)') &
           NameSub, ' after flux/source                   =', &
           State_VGB(iVarTest,iTest,jTest,kTest,iBlock),       &
           Energy_GBI(iTest,jTest,kTest,iBlock,:)

      if(boris_correction) then
         ! Convert StateOld_VGB and State_VGB back from
         ! semi-relativistic to classical MHD variables
         do k=1,nK; do j=1,nJ; do i=1,nI
            if(.not.true_cell(i,j,k,iBlock)) CYCLE

            if(UseB0)then
               call boris_to_mhd(StateOld_VGB(:,i,j,k,iBlock), &
                    EnergyOld_CBI(i,j,k,iBlock,1), B0_DGB(:,i,j,k,iBlock))
               call boris_to_mhd(State_VGB(:,i,j,k,iBlock), &
                    Energy_GBI(i,j,k,iBlock,1), B0_DGB(:,i,j,k,iBlock))
            else
               call boris_to_mhd(StateOld_VGB(:,i,j,k,iBlock), &
                    EnergyOld_CBI(i,j,k,iBlock,1))
               call boris_to_mhd(State_VGB(:,i,j,k,iBlock), &
                    Energy_GBI(i,j,k,iBlock,1))
            end if
         end do; end do; end do

         if(DoTest)write(*,'(2x,2a,15es20.12)') &
              NameSub, ' after boris_to_mhd                  =', &
              State_VGB(iVarTest,iTest,jTest,kTest,iBlock),       &
              Energy_GBI(iTest,jTest,kTest,iBlock,:)
      endif

      if(UseBorisSimple .and. IsMhd) then
         ! Convert mometum in StateOld_VGB and State_VGB back from
         ! enhanced momentum
         do k=1,nK; do j=1,nJ; do i=1,nI
            if(.not.true_cell(i,j,k,iBlock)) CYCLE

            if(UseB0)then
               call boris_simple_to_mhd(StateOld_VGB(:,i,j,k,iBlock), &
                    B0_DGB(:,i,j,k,iBlock))
               call boris_simple_to_mhd(State_VGB(:,i,j,k,iBlock), &
                    B0_DGB(:,i,j,k,iBlock))
            else
               call boris_simple_to_mhd(StateOld_VGB(:,i,j,k,iBlock))
               call boris_simple_to_mhd(State_VGB(:,i,j,k,iBlock))
            end if
         end do; end do; end do

         if(DoTest)write(*,'(2x,2a,15es20.12)') &
              NameSub, ' after boris_to_mhd                  =', &
              State_VGB(iVarTest,iTest,jTest,kTest,iBlock),       &
              Energy_GBI(iTest,jTest,kTest,iBlock,:)
      endif

      if(UseElectronPressure .and. UseElectronEntropy)then
         ! Convert electron entropy back to pressure
         ! Pe = Se^GammaE
         do k=1,nK; do j=1,nJ; do i=1,nI
            if(.not.true_cell(i,j,k,iBlock)) CYCLE

            StateOld_VGB(Pe_,i,j,k,iBlock) = &
                 StateOld_VGB(Pe_,i,j,k,iBlock)**GammaElectron
            State_VGB(Pe_,i,j,k,iBlock) = &
                 State_VGB(Pe_,i,j,k,iBlock)**GammaElectron
         end do; end do; end do
      end if

      if(UseMultiSpecies)then
         ! Fix negative species densities
         State_VGB(SpeciesFirst_:SpeciesLast_,1:nI,1:nJ,1:nK,iBlock) = &
              max(0.0,&
              State_VGB(SpeciesFirst_:SpeciesLast_,1:nI,1:nJ,1:nK,iBlock))

         if(DoReplaceDensity)then
            ! Add up species densities to total density
            do k=1,nK; do j=1,nJ; do i=1,nI
               State_VGB(Rho_,i,j,k,iBlock) = &
                    sum(State_VGB(SpeciesFirst_:SpeciesLast_,i,j,k,iBlock))
            end do; end do; end do
         end if

         if(DoTest)write(*,'(2x,2a,15es20.12)') &
              NameSub, ' after multispecies correct          =', &
              State_VGB(iVarTest,iTest,jTest,kTest,iBlock)

      end if

      if(any(RhoMin_I > 0.0))then
         do iFluid = 1, nFluid
            if(RhoMin_I(iFluid) < 0) CYCLE
            iRho = iRho_I(iFluid)
            do k=1,nK; do j=1,nJ; do i=1,nI
               State_VGB(iRho,i,j,k,iBlock) = max(RhoMin_I(iFluid), &
                    State_VGB(iRho,i,j,k,iBlock))
            end do; end do; end do
         end do

         if(DoTest)write(*,'(2x,2a,15es20.12)') &
              NameSub, ' after min density correct densities =', &
              State_VGB(iRho_I,iTest,jTest,kTest,iBlock)
      end if

      if( IsMhd .and. &
           ((nStage==1.and..not.time_accurate) &
           .or.(nStage==2.and.iStage==1.and.UseHalfStep)))then

         ! A desparate attempt to maintain positivity by adding dB^2/2 to the
         ! energy. This is fine for steady state, and is 2nd order accurate
         ! for half+full step method. But it cannot be used for RK schemes!

         do k=1,nK; do j=1,nJ; do i=1,nI
            Energy_GBI(i,j,k,iBlock,1) = Energy_GBI(i,j,k,iBlock,1) + 0.5*( &
                 Source_VC(Bx_,i,j,k)**2 + &
                 Source_VC(By_,i,j,k)**2 + &
                 Source_VC(Bz_,i,j,k)**2)
         end do; end do; end do

         if(DoTest)write(*,'(2x,2a,15es20.12)') &
              NameSub, ' after energy dB correct             =', &
              State_VGB(iVarTest,iTest,jTest,kTest,iBlock),       &
              Energy_GBI(iTest,jTest,kTest,iBlock,:)

      end if

      if(UseWavePressure)then
         if(DoAdvectWaves .and. iStage==nStage .and. nWave>2)&
              call update_wave_group_advection(iBlock)
         if(UseWavePressureLtd)then
            do k = 1, nK; do j = 1, nJ; do i = 1, nI
               State_VGB(Ew_,i,j,k,iBlock) = &
                    sum(State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock))
            end do; end do; end do
         end if
         ! Avoid negative wave pressure
         if(UseAlfvenWaves)then
            do k = 1, nK; do j = 1, nJ; do i = 1, nI
               State_VGB(WaveFirst_:WaveLast_,i,j,k, iBlock) = &
                    max(State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock), 0.0)
            end do; end do; end do
         end if
      end if

      ! Update energy or pressure based on UseConservative and IsConserv_CB
      call calc_energy_or_pressure(iBlock)

      if(DoTest)write(*,'(2x,2a,15es20.12)') &
           NameSub, ' after pressure/energy update        =', &
           State_VGB(iVarTest,iTest,jTest,kTest,iBlock),       &
           Energy_GBI(iTest,jTest,kTest,iBlock,:)

    end subroutine update_explicit
    !==========================================================================

    subroutine deduct_expl_source()
      integer:: iVarSemi_

      character(len=*), parameter:: NameSub = 'deduct_expl_source'
      !------------------------------------------------------------------------
      if(UseElectronPressure) then
         iVarSemi_  = pe_
      else
         iVarSemi_ = p_
      endif

      do k=1,nK; do j=1,nJ; do i=1,nI
         ! DtLocal = Cfl*time_BLK(i,j,k,iBlock)

         ! For the first iteration, dt = 0;
         ! if(DtLocal < 1e-15) CYCLE
         Source_VCB(iVarSemi_,i,j,k,iBlock) = &
              State_VGB(iVarSemi_,i,j,k,iBlock) - &
              StateOld_VGB(iVarSemi_,i,j,k,iBlock)
         State_VGB(iVarSemi_,i,j,k,iBlock) = &
              StateOld_VGB(iVarSemi_,i,j,k,iBlock)
      end do; end do; end do

    end subroutine deduct_expl_source
    !==========================================================================

  end subroutine update_state_normal
  !============================================================================

  subroutine update_te0

    use ModPhysics, ONLY: UnitTemperature_,Si2No_V
    use ModAdvance, ONLY: State_VGB,  nI, nJ, nK
    use ModMain,    ONLY: nBlock, Unused_B, UseERadInput
    use ModVarIndexes, ONLY: Te0_
    use ModUserInterface ! user_material_properties

    real:: Te0Si
    integer:: i, j, k, iBlock
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'update_te0'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    !\
    ! At the end of time step just calculated values of ERad are used to
    ! calculate Te (and accordingly B(Te)).
    !/
    UseERadInput = .true.

    do iBlock = 1, nBlock
       if(Unused_B(iBlock))CYCLE
       do k=1,nK; do j=1,nJ; do i=1,nI
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               i,j,k,iBlock, TeOut=Te0SI)
          State_VGB(Te0_,i,j,k,iBlock) = Te0SI * Si2No_V(UnitTemperature_)
       end do; end do; end do
    end do

    ! Reset UseERadInput
    UseERadInput = .false.

    call test_stop(NameSub, DoTest)
  end subroutine update_te0
  !============================================================================

  subroutine update_check

    ! Check updated values for allowed change in density or pressure

    use ModProcMH
    use ModMain
    use ModImplicit, ONLY: UsePartImplicit
    use ModAdvance
    use ModB0, ONLY: B0_DGB
    use ModPhysics
    use ModGeometry, ONLY : true_cell
    use ModNumConst, ONLY: cTiny
    use ModMpi
    use ModEnergy
    use ModMultiFluid, ONLY: IsMhd
    use ModMultiIon,   ONLY: DoRestrictMultiIon, IsMultiIon_CB
    use BATL_lib, ONLY: Xyz_DGB

    integer, parameter :: max_checks=25
    integer :: i,j,k, iVar, iBlock,num_checks, iError
    real :: time_fraction_rho, min_time_fraction_rho
    real :: time_fraction_p,   min_time_fraction_p
    real :: time_fraction, cell_time_fraction, report_tf, report_tf_all
    real, dimension(2) :: percent_chg_rho
    real, dimension(2) :: percent_chg_p
    real, dimension(4) :: PercentChangePE, PercentChangeMax

    real :: Value,B0_DC(3,nI,nJ,nK)
    integer :: i_D(3)
    logical :: update_check_done, IsNegative, DoStop
    logical :: DoTest1
    logical :: DoTest2
    logical :: DoTest3

    real    :: RhoChangeMax_I(2), pChangeMax_I(2)
    character(len=*), parameter:: format="(a,i5,a,f6.1,a,3es11.3)"

    integer :: iError1=-1
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'update_check'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    call test_start('convergence_history', DoTest1)
    call test_start('update_check_detail', DoTest2)
    call test_start('locations',           DoTest3)

    ! Check for allowable percentage changed from update
    if(time_accurate) then
       !\\\
       !    TIME ACCURATE  ===================================================
       !///
       report_tf = 1.
       do num_checks = 1,max_checks
          percent_chg_rho = 0.1
          percent_chg_p   = 0.1
          do iBlock = 1, nBlockMax
             if (Unused_B(iBlock)) CYCLE
             if (num_checks == 1) then
                do k=1,nK; do j=1,nJ; do i=1,nI
                   do iVar = 1, nVar
                      if (DefaultState_V(iVar) <= cTiny) CYCLE

                      ! For sake of backward compatibility
                      if(iVar == P_) CYCLE

                      ! Do not check multi-ion variables in regions that are
                      ! not truely multi-ion
                      if(UseMultiIon .and. iVar > p_ &
                           .and. DoRestrictMultiIon)then
                         if(.not.IsMultiIon_CB(i,j,k,iBlock)) CYCLE
                      end if

                      ! Do not check minor species if not necessary
                      if(UseMultiSpecies .and. &
                           iVar >= SpeciesFirst_ .and. iVar <= SpeciesLast_ &
                           .and. StateOld_VGB(iVar,i,j,k,iBlock) &
                           < SpeciesPercentCheck*0.01*&
                           StateOld_VGB(Rho_,i,j,k,iBlock)) CYCLE

                      percent_chg_rho(1) = &
                           max(percent_chg_rho(1), 100*abs( min(0.,&
                           (State_VGB(iVar,i,j,k,iBlock)- &
                           StateOld_VGB(iVar,i,j,k,iBlock)) &
                           /StateOld_VGB(iVar,i,j,k,iBlock) ) ) )
                      percent_chg_rho(2) = &
                           max(percent_chg_rho(2), 100*abs( max(0.,&
                           (State_VGB(iVar,i,j,k,iBlock)- &
                           StateOld_VGB(iVar,i,j,k,iBlock)) &
                           /StateOld_VGB(iVar,i,j,k,iBlock) ) ) )
                   end do
                end do; end do; end do
             end if
             percent_chg_p(1) = &
                  max(percent_chg_p(1), 100. * abs( min( 0., minval( &
                  (State_VGB(P_,1:nI,1:nJ,1:nK,iBlock)- &
                  StateOld_VGB(P_,1:nI,1:nJ,1:nK,iBlock)) &
                  /StateOld_VGB(P_,1:nI,1:nJ,1:nK,iBlock) ) ) ) )
             percent_chg_p(2) = &
                  max(percent_chg_p(2), 100. * abs( max( 0., maxval( &
                  (State_VGB(P_,1:nI,1:nJ,1:nK,iBlock)- &
                  StateOld_VGB(P_,1:nI,1:nJ,1:nK,iBlock)) &
                  /StateOld_VGB(P_,1:nI,1:nJ,1:nK,iBlock) ) ) ) )
          end do
          if(DoTest)then
             ! Find location of maximum change
             call MPI_allreduce(percent_chg_rho, RhoChangeMax_I, 2, &
                  MPI_REAL, MPI_MAX, iComm, iError)
             call MPI_allreduce(percent_chg_p, pChangeMax_I, 2, &
                  MPI_REAL, MPI_MAX, iComm, iError)

             do iBlock = 1, nBlockMax
                if (Unused_B(iBlock)) CYCLE
                do k=1,nK; do j=1,nJ; do i=1,nI
                   do iVar = 1, nVar
                      if (DefaultState_V(iVar) <= cTiny) CYCLE

                      if(UseMultiSpecies .and. &
                           iVar >= SpeciesFirst_ .and. iVar <= SpeciesLast_ &
                           .and. StateOld_VGB(iVar,i,j,k,iBlock) &
                           < SpeciesPercentCheck*0.01*&
                           StateOld_VGB(Rho_,i,j,k,iBlock)) CYCLE

                      if(iVar == p_)then
                         if(pChangeMax_I(1) > percent_max_p(1) .and. &
                              1e-4 > abs(pChangeMax_I(1) - 100. * abs( &
                              (   State_VGB(P_,i,j,k,iBlock)- &
                              StateOld_VGB (P_,i,j,k,iBlock)) &
                              /StateOld_VGB(P_,i,j,k,iBlock) ))) &
                              write(*,format)NameSub//' nStep=',n_step,&
                              ' max p drop=',pChangeMax_I(1),&
                              '% at x,y,z=',&
                              Xyz_DGB(:,i,j,k,iBlock)

                         if(pChangeMax_I(2) > percent_max_p(2) .and. &
                              1e-4 > abs(pChangeMax_I(2) - 100. * abs( &
                              (   State_VGB(P_,i,j,k,iBlock)- &
                              StateOld_VGB (P_,i,j,k,iBlock)) &
                              /StateOld_VGB(P_,i,j,k,iBlock)  ))) &
                              write(*,format)NameSub//' nStep=',n_step,&
                              ' max p increase=',&
                              pChangeMax_I(2), &
                              '% at x,y,z=',&
                              Xyz_DGB(:,i,j,k,iBlock)

                         CYCLE
                      end if

                      if(RhoChangeMax_I(1) > percent_max_rho(1) .and. &
                           1e-4 > abs(RhoChangeMax_I(1) - 100*abs( &
                           (State_VGB(iVar,i,j,k,iBlock)- &
                           StateOld_VGB(iVar,i,j,k,iBlock)) &
                           /StateOld_VGB(iVar,i,j,k,iBlock) ))) &
                           write(*,format)NameSub//' nStep=',n_step,&
                           ' max '//trim(NameVar_V(iVar))//' drop=', &
                           RhoChangeMax_I(1), &
                           '% at x,y,z=',&
                           Xyz_DGB(:,i,j,k,iBlock)

                      if(RhoChangeMax_I(2) > percent_max_rho(2) .and. &
                           1e-4 > abs(RhoChangeMax_I(2) - 100*abs( &
                           (State_VGB(iVar,i,j,k,iBlock)- &
                           StateOld_VGB(iVar,i,j,k,iBlock)) &
                           /StateOld_VGB(iVar,i,j,k,iBlock) ))) &
                           write(*,format)NameSub//' nStep=',n_step,&
                           ' max '//trim(NameVar_V(iVar))//' increase=',&
                           RhoChangeMax_I(2), &
                           '% at x,y,z=',&
                           Xyz_DGB(:,i,j,k,iBlock)
                   end do
                end do; end do; end do
             end do
          end if ! DoTest
          time_fraction_rho = 1.0 / maxval(percent_chg_rho/percent_max_rho)
          call MPI_allreduce(time_fraction_rho, min_time_fraction_rho, 1, &
               MPI_REAL, MPI_MIN, iComm, iError)
          time_fraction_p   = 1.0 / maxval(percent_chg_p  /percent_max_p  )
          call MPI_allreduce(time_fraction_p, min_time_fraction_p, 1, &
               MPI_REAL, MPI_MIN, iComm, iError)
          if (min_time_fraction_rho >= 1. .and. min_time_fraction_p >= 1.) EXIT

          if (num_checks == 1) then
             time_fraction = 1.
             if (min_time_fraction_rho < 1.) &
                  time_fraction = 0.9*min_time_fraction_rho
             if (min_time_fraction_p   < 1.) &
                  time_fraction = min(time_fraction, 0.75)
          else
             time_fraction = 0.5
          end if
          dt = dt * time_fraction
          report_tf = report_tf*time_fraction
          do iBlock = 1, nBlockMax
             if (Unused_B(iBlock)) CYCLE
             if(UseB0)then
                B0_DC=B0_DGB(:,1:nI,1:nJ,1:nK,iBlock)
             else
                B0_DC=0.00
             end if
             ! Fix the update in the cells
             do k=1,nK; do j=1,nJ; do i=1,nI
                call fix_update
             end do; end do; end do
          end do
       end do
       PercentChangePE(1:2) =  percent_chg_rho(1:2) - 0.1
       PercentChangePE(3:4) =  percent_chg_p(1:2) - 0.1

       ! The part implicit scheme can get here if all blocks become explicit
       ! due to time step reductions. To be able to recover the time step,
       ! increase fixed time step if there was no time step reduction above.
       if(UsePartImplicit .and. dt == DtFixed) &
            DtFixed = min(DtFixedOrig, DtFixed*1.05)

       if(DoTest) then
          if (iProc == 0 .and. report_tf < 1.) &
               write(*,'(a,a,i6,a,f12.8,a,f12.8)') 'update_check TA:', &
               ' nStep=',n_step,'     dt reduction=',report_tf,' dt=',dt
       end if
    else
       !\\\
       !    LOCAL TIMESTEPPING
       !///
       report_tf = 1.
       PercentChangePE = 0.
       do iBlock = 1, nBlockMax
          if (Unused_B(iBlock)) CYCLE
          if(UseB0)then
             B0_DC=B0_DGB(:,1:nI,1:nJ,1:nK,iBlock)
          else
             B0_DC=0.00
          end if
          do k=1,nK; do j=1,nJ; do i=1,nI
             cell_time_fraction = 1.
             do num_checks = 1, max_checks
                update_check_done = .true.
                percent_chg_rho = 0.1
                percent_chg_p   = 0.1
                if (num_checks == 1) then
                   do iVar = 1, nVar
                      if (DefaultState_V(iVar) <= cTiny) CYCLE

                      ! This is for backwards compatibility
                      if(iVar == P_) CYCLE

                      if(UseMultiSpecies .and. &
                           iVar >= SpeciesFirst_ .and. iVar <= SpeciesLast_ &
                           .and. StateOld_VGB(iVar,i,j,k,iBlock) < &
                           SpeciesPercentCheck*0.01*&
                           StateOld_VGB(Rho_,i,j,k,iBlock)) CYCLE

                      percent_chg_rho(1) = max(percent_chg_rho(1), &
                           0.1 + 100. * abs( min( 0., &
                           (State_VGB(iVar,i,j,k,iBlock)-&
                           StateOld_VGB(iVar,i,j,k,iBlock)) &
                           /StateOld_VGB(iVar,i,j,k,iBlock) ) ) )

                      percent_chg_rho(2) = max(percent_chg_rho(2), &
                           0.1 + 100. * abs( max( 0., &
                           (State_VGB(iVar,i,j,k,iBlock)-&
                           StateOld_VGB(iVar,i,j,k,iBlock)) &
                           /StateOld_VGB(iVar,i,j,k,iBlock) ) ) )
                   end do
                end if
                percent_chg_p(1) = 0.1 + 100. * abs( min( 0., &
                     (State_VGB(P_,i,j,k,iBlock)-&
                     StateOld_VGB(P_,i,j,k,iBlock)) &
                     /StateOld_VGB(P_,i,j,k,iBlock) ) )
                percent_chg_p(2) = 0.1 + 100. * abs( max( 0., &
                     (State_VGB(P_,i,j,k,iBlock)-&
                     StateOld_VGB(P_,i,j,k,iBlock)) &
                     /StateOld_VGB(P_,i,j,k,iBlock) ) )
                time_fraction_rho = 1/maxval(percent_chg_rho/percent_max_rho)
                time_fraction_p   = 1/maxval(percent_chg_p  /percent_max_p  )
                if (time_fraction_rho < 1. .or. time_fraction_p < 1.) then
                   if (num_checks == 1) then
                      time_fraction = 1.
                      if (time_fraction_rho < 1.) &
                           time_fraction = 0.9*time_fraction_rho
                      if (time_fraction_p   < 1.) &
                           time_fraction = min(time_fraction, 0.75)
                   else
                      time_fraction = 0.5
                   end if
                   update_check_done = .false.
                   cell_time_fraction = cell_time_fraction * time_fraction
                   if(DoTest2) then
                      write(*,*) &
                           'update_check LT: changing cell value, PE=',iProc, &
                           ' BLK=',iBlock,' i,j,k=',i,' ',j,' ',k, &
                           '  time_fraction=',time_fraction
                      write(*,*) &
                           iProc,' ',iBlock,' ',i,' ',j,' ',k,' OLD:  ', &
                           NameVar_V(1),'=',StateOld_VGB(1,i,j,k,iBlock),&
                           '    ',&
                           NameVar_V(nVar),' ', StateOld_VGB(nVar,i,j,k,iBlock)
                      write(*,*) &
                           iProc,' ',iBlock,' ',i,' ',j,' ',k,' BAD: ', &
                           NameVar_V(1),'=',State_VGB(1,i,j,k,iBlock), &
                           '   ', NameVar_V(nVar),State_VGB(nVar,i,j,k,iBlock)
                   end if
                   call fix_update
                   if(DoTest2) then
                      write(*,*) &
                           iProc,' ',iBlock,' ',i,' ',j,' ',k,' NEW: ', &
                           NameVar_V(Rho_),'=',State_VGB(Rho_,i,j,k,iBlock),&
                           '   ',NameVar_V(p_),'=', State_VGB(p_,i,j,k,iBlock)
                   end if
                end if
                if (update_check_done) EXIT
             end do
             PercentChangePE(1:2) = &
                  max(percent_chg_rho(1:2)-0.1, PercentChangePE(1:2))
             PercentChangePE(3:4) = &
                  max(percent_chg_p(1:2)-0.1, PercentChangePE(3:4))
             report_tf = min(report_tf, cell_time_fraction)
          end do; end do; end do
       end do
       call MPI_allreduce(report_tf, report_tf_all, 1, &
            MPI_REAL, MPI_MIN, iComm, iError)
       report_tf = report_tf_all
       if(DoTest) then
          if (iProc == 0 .and. report_tf < 1.) &
               write(*,'(a,a,i6,a,f12.8)') 'update_check LT:', &
               ' nStep=',n_step,' max dt reduction=',report_tf
       end if
    end if

    if(DoTest1.and.iStage==nStage) then
       call MPI_allreduce(PercentChangePE,  PercentChangeMax, 4, &
            MPI_REAL, MPI_MAX, iComm, iError)
       if(iProc==0) then
          write(*,*)'Maximum change in pressure on proc 0:',&
               - PercentChangeMax(3),' %,   ',  PercentChangeMax(4),' %'
          write(*,*) 'Maximum change in other positive variables:', &
               - PercentChangeMax(1),' %,   ',  PercentChangeMax(2),' %'
       end if
       if(DoTest3) then
          do iBlock = 1,nBlockMax
             if(Unused_B(iBlock))CYCLE
             do k=1,nK;do j=1,nJ; do i=1,nI
                if(.not.true_cell(i,j,k,iBlock))CYCLE
                if (abs(100. * abs( min( 0., &
                     (State_VGB(Rho_,i,j,k,iBlock)-&
                     StateOld_VGB(Rho_,i,j,k,iBlock)) &
                     /StateOld_VGB(Rho_,i,j,k,iBlock) ) )-&
                     PercentChangeMax(1))<cTiny*PercentChangeMax(1))&
                     write(*,*)'Maximum decrease in density at X Y Z=',&
                     Xyz_DGB(:,i,j,k,iBlock),&
                     ': rho_old = ',StateOld_VGB(Rho_,i,j,k,iBlock),&
                     ' rho_new = ',State_VGB(Rho_,i,j,k,iBlock)

                if (abs(100. * abs( max( 0., &
                     (State_VGB(Rho_,i,j,k,iBlock)-&
                     StateOld_VGB(Rho_,i,j,k,iBlock)) &
                     /StateOld_VGB(Rho_,i,j,k,iBlock) ) )-&
                     PercentChangeMax(2))<cTiny*PercentChangeMax(2))&
                     write(*,*)'Maximum increase in density at the point',&
                     Xyz_DGB(:,i,j,k,iBlock),&
                     'is: rho_old = ',&
                     StateOld_VGB(Rho_,i,j,k,iBlock),&
                     'rho_new=',State_VGB(Rho_,i,j,k,iBlock)

                if (abs(100. * abs( min( 0., &
                     (State_VGB(p_,i,j,k,iBlock)-&
                     StateOld_VGB(p_,i,j,k,iBlock)) &
                     /StateOld_VGB(p_,i,j,k,iBlock) ) )-&
                     PercentChangeMax(3))<cTiny*PercentChangeMax(3))&
                     write(*,*)'Maximum decrease in',NameVar_V(p_), &
                     'at the point',&
                     Xyz_DGB(:,i,j,k,iBlock),&
                     'is: valeu_old = ',StateOld_VGB(P_,i,j,k,iBlock),&
                     'value_new=',State_VGB(p_,i,j,k,iBlock)
                if (abs(100. * abs( max( 0., &
                     (State_VGB(p_,i,j,k,iBlock)-&
                     StateOld_VGB(p_,i,j,k,iBlock)) &
                     /StateOld_VGB(p_,i,j,k,iBlock) ) )-&
                     PercentChangeMax(4))<cTiny*PercentChangeMax(4))&
                     write(*,*)'Maximum increase in',NameVar_V(p_), &
                     'at the point',&
                     Xyz_DGB(:,i,j,k,iBlock),&
                     'is: value_old = ',StateOld_VGB(p_,i,j,k,iBlock),&
                     'value_new=',State_VGB(p_,i,j,k,iBlock)
             end do; end do; end do
          end do
       end if
    end if

    if(iProc == 0 .and. report_tf<1.)&
         call error_report('Time step reduction, min(factor)',&
         report_tf,iError1,.true.)

    ! Check for positivity of variables
    IsNegative = .false.
    do iBlock = 1, nBlockMax
       if (Unused_B(iBlock)) CYCLE
       do iVar = 1, nVar
          ! Ignore variables that do not have to be positive
          if(DefaultState_V(iVar) < cTiny) CYCLE

          ! Do not check species densities if check threshold is positive
          ! (i.e. minor species densities are allowed to go negative,
          !       and they are fixed to be zero in update_state)
          if(UseMultiSpecies .and. SpeciesPercentCheck > 0.0 &
               .and. iVar >= SpeciesFirst_ .and. iVar <= SpeciesLast_) CYCLE

          Value = minval(State_VGB(iVar,1:nI,1:nJ,1:nK,iBlock))

          if(Value < 0.0)then
             i_D = minloc(State_VGB(iVar,1:nI,1:nJ,1:nK,iBlock))
             i = i_D(1); j = i_D(2); k = i_D(3)
             write (*,'(a,3i3,2i5,i3,a,3f12.4,/,5x,a,a,es12.4)') &
                  ' I J K iBlock iProc iVar=',i_D,iBlock,iProc,iVar, &
                  ' X Y Z=', &
                  Xyz_DGB(:,i,j,k,iBlock), &
                  ' Var='//trim(NameVar_V(iVar)), &
                  ' Value=', State_VGB(iVar,i,j,k,iBlock)
             IsNegative = .true.
          end if
       end do
    end do
    if (IsNegative) then
       if(time_accurate) then
          write(*,'(a,i4,a,a,i6,a,f12.8,a,f12.8)') &
               'Negative updated value: PE=',iProc, &
               'update_check TA:',' nStep=',n_step, &
               '     dt reduction=',report_tf,' dt=',dt
       else
          write(*,'(a,i4,a,a,i6,a,f12.8)') &
               'Negative updated value: PE=',iProc, &
               'update_check LT:',' nStep=',n_step, &
               ' max dt reduction=',report_tf
       end if
    end if

    call MPI_allreduce(IsNegative, DoStop,1,MPI_LOGICAL,MPI_LOR,iComm,iError)
    if (DoStop) call stop_mpi('Stopping, negative density or pressure')

    call test_stop(NameSub, DoTest)
  contains
    !==========================================================================

    subroutine fix_update

      logical :: IsConserv

      real :: fullBx, fullBy, fullBz, fullBB, UdotBc2, rhoc2, gA2_Boris
      real :: rhoUx_Boris, rhoUy_Boris, rhoUz_Boris, E_Boris, &
           rhoUx_o_Boris, rhoUy_o_Boris, rhoUz_o_Boris, E_o_Boris

      logical, parameter :: DoTestCell = .false.
      character(len=*), parameter:: NameSub = 'fix_update'
      !------------------------------------------------------------------------
      ! DoTestCell = DoTestMe .and. iBlock==iBlockTest .and. &
      !     i==iTest .and. j==jTest .and. k==kTest

      if(allocated(IsConserv_CB))then
         IsConserv = IsConserv_CB(i,j,k,iBlock)
      else
         IsConserv = .not. UseNonConservative
      end if

      if(DoTestCell)then
         write(*,*)NameSub,' IsConserv    =',IsConserv
         write(*,*)NameSub,' initial state=',State_VGB(:,i,j,k,iBlock)
         write(*,*)NameSub,' old     state=',StateOld_VGB(:,i,j,k,iBlock)
      end if

      if(boris_correction .and. nFluid==1) then

         ! Convert old state
         fullBx = B0_DC(x_,i,j,k) + &
              StateOld_VGB(Bx_,i,j,k,iBlock)
         fullBy = B0_DC(y_,i,j,k) + &
              StateOld_VGB(By_,i,j,k,iBlock)
         fullBz = B0_DC(z_,i,j,k) + &
              StateOld_VGB(Bz_,i,j,k,iBlock)
         fullBB = fullBx**2 + fullBy**2 + fullBz**2
         rhoc2  = &
              StateOld_VGB(rho_,i,j,k,iBlock)*c2LIGHT
         UdotBc2= (StateOld_VGB(rhoUx_,i,j,k,iBlock)*fullBx + &
              StateOld_VGB(rhoUy_,i,j,k,iBlock)*fullBy + &
              StateOld_VGB(rhoUz_,i,j,k,iBlock)*fullBz)/rhoc2
         gA2_Boris=1.+fullBB/rhoc2

         ! rhoU_Boris = rhoU - ((U x B) x B)/c^2
         !            = rhoU + (U B^2 - B U.B)/c^2
         !            = rhoU*(1+BB/(rho*c2)) - B UdotB/c^2
         rhoUx_o_Boris = StateOld_VGB(rhoUx_,i,j,k,iBlock)*ga2_Boris - &
              fullBx*UdotBc2
         rhoUy_o_Boris = StateOld_VGB(rhoUy_,i,j,k,iBlock)*ga2_Boris - &
              fullBy*UdotBc2
         rhoUz_o_Boris = StateOld_VGB(rhoUz_,i,j,k,iBlock)*ga2_Boris - &
              fullBz*UdotBc2

         if(IsConserv)then
            ! e_boris = e + 0.5/c^2 * (V x B)^2
            E_o_Boris = EnergyOld_CBI(i,j,k,iBlock,1) + (0.5/c2LIGHT)*( &
                 ((StateOld_VGB(rhoUy_,i,j,k,iBlock)*fullBz     &
                 -StateOld_VGB(rhoUz_,i,j,k,iBlock)*fullBy)**2 &
                 +(StateOld_VGB(rhoUx_,i,j,k,iBlock)*fullBz     &
                 -StateOld_VGB(rhoUz_,i,j,k,iBlock)*fullBx)**2 &
                 +(StateOld_VGB(rhoUx_,i,j,k,iBlock)*fullBy     &
                 -StateOld_VGB(rhoUy_,i,j,k,iBlock)*fullBx)**2 &
                 )/StateOld_VGB(rho_,i,j,k,iBlock)**2                 )
         end if

         ! Convert current state
         fullBx = B0_DC(x_,i,j,k) + State_VGB(Bx_,i,j,k,iBlock)
         fullBy = B0_DC(y_,i,j,k) + State_VGB(By_,i,j,k,iBlock)
         fullBz = B0_DC(z_,i,j,k) + State_VGB(Bz_,i,j,k,iBlock)
         fullBB = fullBx**2 + fullBy**2 + fullBz**2
         rhoc2  = State_VGB(rho_,i,j,k,iBlock)*c2LIGHT
         UdotBc2= (State_VGB(rhoUx_,i,j,k,iBlock)*fullBx + &
              State_VGB(rhoUy_,i,j,k,iBlock)*fullBy + &
              State_VGB(rhoUz_,i,j,k,iBlock)*fullBz)/rhoc2
         gA2_Boris = 1 + fullBB/rhoc2

         ! rhoU_Boris = rhoU - ((U x B) x B)/c^2
         !            = rhoU + (U B^2 - B U.B)/c^2
         !            = rhoU*(1+BB/(rho*c2)) - B UdotB/c^2
         rhoUx_Boris = State_VGB(rhoUx_,i,j,k,iBlock)*ga2_Boris - fullBx*UdotBc2
         rhoUy_Boris = State_VGB(rhoUy_,i,j,k,iBlock)*ga2_Boris - fullBy*UdotBc2
         rhoUz_Boris = State_VGB(rhoUz_,i,j,k,iBlock)*ga2_Boris - fullBz*UdotBc2

         if(IsConserv)then
            ! e_boris = e + 0.5/c^2 * (V x B)^2
            E_Boris = Energy_GBI(i,j,k,iBlock,1) + (0.5/c2LIGHT)*( &
                 ((State_VGB(rhoUy_,i,j,k,iBlock)*fullBz     &
                 -State_VGB(rhoUz_,i,j,k,iBlock)*fullBy)**2 &
                 +(State_VGB(rhoUx_,i,j,k,iBlock)*fullBz     &
                 -State_VGB(rhoUz_,i,j,k,iBlock)*fullBx)**2 &
                 +(State_VGB(rhoUx_,i,j,k,iBlock)*fullBy     &
                 -State_VGB(rhoUy_,i,j,k,iBlock)*fullBx)**2 &
                 )/State_VGB(rho_,i,j,k,iBlock)**2                 )
         end if

         ! Interpolate
         ! For possible extension to multifluid
         ! State_VGB(iRho_I,i,j,k,iBlock) = &
         !     (   time_fraction) *   State_VGB(iRho_I,i,j,k,iBlock) + &
         !     (1.0-time_fraction) * StateOld_VGB(iRho_I,i,j,k,iBlock)
         State_VGB(rho_,i,j,k,iBlock) = &
              (   time_fraction) *   State_VGB(rho_,i,j,k,iBlock) + &
              (1.0-time_fraction) * StateOld_VGB(rho_,i,j,k,iBlock)
         rhoUx_Boris = &
              (   time_fraction) * rhoUx_Boris + &
              (1.0-time_fraction) * rhoUx_o_Boris
         rhoUy_Boris = &
              (   time_fraction) * rhoUy_Boris + &
              (1.0-time_fraction) * rhoUy_o_Boris
         rhoUz_Boris = &
              (   time_fraction) * rhoUz_Boris + &
              (1.0-time_fraction) * rhoUz_o_Boris
         State_VGB(Bx_,i,j,k,iBlock) = &
              (   time_fraction) *   State_VGB(Bx_,i,j,k,iBlock) + &
              (1.0-time_fraction) * StateOld_VGB(Bx_,i,j,k,iBlock)
         State_VGB(By_,i,j,k,iBlock) = &
              (   time_fraction) *   State_VGB(By_,i,j,k,iBlock) + &
              (1.0-time_fraction) * StateOld_VGB(By_,i,j,k,iBlock)
         State_VGB(Bz_,i,j,k,iBlock) = &
              (   time_fraction) *   State_VGB(Bz_,i,j,k,iBlock) + &
              (1.0-time_fraction) * StateOld_VGB(Bz_,i,j,k,iBlock)

         ! Convert Back
         fullBx = B0_DC(x_,i,j,k) + State_VGB(Bx_,i,j,k,iBlock)
         fullBy = B0_DC(y_,i,j,k) + State_VGB(By_,i,j,k,iBlock)
         fullBz = B0_DC(z_,i,j,k) + State_VGB(Bz_,i,j,k,iBlock)
         fullBB = fullBx**2 + fullBy**2 + fullBz**2
         rhoc2  = State_VGB(rho_,i,j,k,iBlock)*c2LIGHT
         UdotBc2= (rhoUx_Boris*fullBx + rhoUy_Boris*fullBy + rhoUz_Boris*fullBz)&
              /rhoc2
         gA2_Boris= 1.0/(1.0+fullBB/rhoc2)

         ! rhoU = 1/(rho c^2 + B^2) * (I rho c^2 + B B) * rhoU_Boris
         !      = 1/[1+BB/(rho c^2)]* (rhoU_Boris + (rhoUBorisdotB/(rho c^2) * B)

         State_VGB(rhoUx_,i,j,k,iBlock) =  gA2_Boris * (rhoUx_Boris + UdotBc2*fullBx)
         State_VGB(rhoUy_,i,j,k,iBlock) =  gA2_Boris * (rhoUy_Boris + UdotBc2*fullBy)
         State_VGB(rhoUz_,i,j,k,iBlock) =  gA2_Boris * (rhoUz_Boris + UdotBc2*fullBz)

         if(IsConserv)then
            E_boris= &
                 (   time_fraction) * E_Boris + &
                 (1.0-time_fraction) * E_o_Boris

            ! E = E_boris - 0.5/c^2 * (V x B)^2
            Energy_GBI(i,j,k,iBlock,1) = E_Boris - (0.5/c2LIGHT)*( &
                 ((State_VGB(rhoUy_,i,j,k,iBlock)*fullBz     &
                 -State_VGB(rhoUz_,i,j,k,iBlock)*fullBy)**2 &
                 +(State_VGB(rhoUx_,i,j,k,iBlock)*fullBz     &
                 -State_VGB(rhoUz_,i,j,k,iBlock)*fullBx)**2 &
                 +(State_VGB(rhoUx_,i,j,k,iBlock)*fullBy     &
                 -State_VGB(rhoUy_,i,j,k,iBlock)*fullBx)**2 &
                 )/State_VGB(rho_,i,j,k,iBlock)**2               )

            if((nStage==1.and..not.time_accurate).or.&
                 (nStage>1.and.iStage==1)) &
                 Energy_GBI(i,j,k,iBlock,1) =  Energy_GBI(i,j,k,iBlock,1) - &
                 (0.5/time_fraction - 0.5)*&
                 ((State_VGB(Bx_,i,j,k,iBlock) - &
                 StateOld_VGB(Bx_,i,j,k,iBlock))**2 +&
                 (State_VGB(By_,i,j,k,iBlock) - &
                 StateOld_VGB(By_,i,j,k,iBlock))**2 +&
                 (State_VGB(Bz_,i,j,k,iBlock) - &
                 StateOld_VGB(Bz_,i,j,k,iBlock))**2 )

            call calc_pressure(i,i,j,j,k,k,iBlock,1,1)

            ! For multifluid update all other energies and
            ! call calc_pressure_point
         else
            ! For possible extension to multifluid
            ! State_VGB(iP_I,i,j,k,iBlock) = &
            !     (   time_fraction) *   State_VGB(iP_I,i,j,k,iBlock) + &
            !     (1.0-time_fraction) * StateOld_VGB(iP_I,i,j,k,iBlock)
            ! call calc_energy_point(i,j,k,iBlock)

            State_VGB(p_,i,j,k,iBlock) = &
                 (   time_fraction) *   State_VGB(p_,i,j,k,iBlock) + &
                 (1.0-time_fraction) * StateOld_VGB(p_,i,j,k,iBlock)

            call calc_energy(i,i,j,j,k,k,iBlock,1,1)
         end if
      else                      ! Non-Boris interpolation

         State_VGB(1:nVar,i,j,k,iBlock) = &
              (   time_fraction) *   State_VGB(1:nVar,i,j,k,iBlock) + &
              (1.0-time_fraction) * StateOld_VGB(1:nVar,i,j,k,iBlock)
         if(IsConserv)then
            Energy_GBI(i,j,k,iBlock,:) = &
                 (   time_fraction) *   Energy_GBI(i,j,k,iBlock,:) + &
                 (1.0-time_fraction) * EnergyOld_CBI(i,j,k,iBlock,:)

            if(IsMhd .and. (nStage==1.and..not.time_accurate).or.&
                 (nStage>1.and.iStage==1)) then
               Energy_GBI(i,j,k,iBlock,1) = &
                    Energy_GBI(i,j,k,iBlock,1) - &
                    (0.5/time_fraction - 0.5)*&
                    sum((State_VGB(Bx_:Bz_,i,j,k,iBlock) - &
                    StateOld_VGB(Bx_:Bz_,i,j,k,iBlock))**2)
            end if

            call calc_pressure_point(i,j,k,iBlock)
         else
            call calc_energy_point(i,j,k,iBlock)
         end if
      end if
      time_BLK(i,j,k,iBlock) = time_BLK(i,j,k,iBlock)*time_fraction

      if(DoTestCell)write(*,*)NameSub,' final state=',State_VGB(:,i,j,k,iBlock)

    end subroutine fix_update
    !==========================================================================

  end subroutine update_check
  !============================================================================

  subroutine select_conservative

    ! Set the global variable IsConserv_CB

    use ModNumConst
    use ModMain
    use ModAdvance
    use ModB0, ONLY: B0_DGB
    use ModGeometry
    use BATL_lib, ONLY: Xyz_DGB

    integer :: iBlock, iCrit, i, j, k

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'select_conservative'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    call timing_start('nonconservative')

    if(DoTest)write(*,*) NameSub,': starting with ',&
         'UseNonConservative, nConservCrit=',UseNonConservative, nConservCrit

    if(.not.allocated(IsConserv_CB))then
       allocate(IsConserv_CB(nI,nJ,nK,MaxBlock))
       if(DoTest)write(*,*) NameSub,': allocated IsConserv_CB'
    end if

    if(nConservCrit < 1)then
       ! There are no criteria so use fully non-conservative
       IsConserv_CB = .false.
       if(DoTest)write(*,*) NameSub,': set IsConserv_CB = F'
       RETURN
    endif

    if(any(TypeConservCrit_I == 'p' .or. TypeConservCrit_I == 'gradp' &
         .or. TypeConservCrit_I == 'jumpp') )then

       if(DoTest)write(*,*) NameSub, ': Apply physics based criteria'

       ! These all have to be true to use non-conservative,
       ! so any of the criteria can switch to conservative
       IsConserv_CB = .false.

       do iBlock = 1, nBlock
          if( Unused_B(iBlock) ) CYCLE

          if(UseElectronPressure)then
             do k = MinK,MaxK; do j = MinJ,MaxJ; do i = MinI,MaxI
                State_VGB(p_,i,j,k,iBlock) = State_VGB(p_,i,j,k,iBlock) &
                     + State_VGB(Pe_,i,j,k,iBlock)
             end do; end do; end do
          end if

          do iCrit = 1, nConservCrit
             select case(TypeConservCrit_I(iCrit))
             case('p')
                if(UseB0)then
                   do k=1,nK; do j=1,nJ; do i=1,nI
                      IsConserv_CB(i,j,k,iBlock) = &
                           IsConserv_CB(i,j,k,iBlock) .or. &
                           State_VGB(P_,i,j,k,iBlock) > pCoeffConserv * &
                           (Energy_GBI(i,j,k,iBlock,1) + 0.5 * &
                           ((State_VGB(Bx_,i,j,k,iBlock) &
                           + B0_DGB(x_,i,j,k,iBlock))**2 &
                           +(State_VGB(By_,i,j,k,iBlock) &
                           + B0_DGB(y_,i,j,k,iBlock))**2 &
                           +(State_VGB(Bz_,i,j,k,iBlock) &
                           + B0_DGB(z_,i,j,k,iBlock))**2 &
                           -State_VGB(Bx_,i,j,k,iBlock)**2 &
                           -State_VGB(By_,i,j,k,iBlock)**2 &
                           -State_VGB(Bz_,i,j,k,iBlock)**2 &
                           ))
                   end do;end do;end do
                else
                   do k=1,nK; do j=1,nJ; do i=1,nI
                      IsConserv_CB(i,j,k,iBlock) = &
                           IsConserv_CB(i,j,k,iBlock) .or. &
                           State_VGB(P_,i,j,k,iBlock) > pCoeffConserv * &
                           Energy_GBI(i,j,k,iBlock,1)
                   end do;end do;end do
                end if
             case('gradp')
                ! Switch to conservative if gradient of pressure is large
                do k=1,nK; do j=1,nJ; do i=1,nI
                   IsConserv_CB(i,j,k,iBlock) = IsConserv_CB(i,j,k,iBlock) &
                        .or. &
                        (abs(State_VGB(P_,i+1,j,k,iBlock) &
                        -    State_VGB(P_,i-1,j,k,iBlock))  &
                        +abs(State_VGB(P_,i,j+1,k,iBlock) &
                        -    State_VGB(P_,i,j-1,k,iBlock))  &
                        +abs(State_VGB(P_,i,j,k+1,iBlock) &
                        -    State_VGB(P_,i,j,k-1,iBlock))) &
                        > GradPCoeffConserv * min(    &
                        State_VGB(P_,i,j,k,iBlock),   &
                        State_VGB(P_,i+1,j,k,iBlock), &
                        State_VGB(P_,i-1,j,k,iBlock), &
                        State_VGB(P_,i,j+1,k,iBlock), &
                        State_VGB(P_,i,j-1,k,iBlock), &
                        State_VGB(P_,i,j,k+1,iBlock), &
                        State_VGB(P_,i,j,k-1,iBlock))
                end do; end do; end do
             case('jumpp')
                ! Switch to conservative if pressure jump is large
                do k=1,nK; do j=1,nJ; do i=1,nI
                   IsConserv_CB(i,j,k,iBlock) = IsConserv_CB(i,j,k,iBlock) &
                        .or. &
                        maxval(State_VGB(P_,i-2:i+2,j,k,iBlock)) &
                        > GradPCoeffConserv* &
                        minval(State_VGB(P_,i-2:i+2,j,k,iBlock)) .or. &
                        nJ > 1 .and. &
                        maxval(State_VGB(P_,i,j-2:j+2,k,iBlock)) &
                        > GradPCoeffConserv* &
                        minval(State_VGB(P_,i:i,j-2:j+2,k,iBlock)) .or. &
                        nK > 1 .and. &
                        maxval(State_VGB(P_,i,j,k-2:k+2,iBlock)) &
                        > GradPCoeffConserv* &
                        minval(State_VGB(P_,i:i,j,k-2:k+2,iBlock))
                end do; end do; end do
             case default
                CYCLE
             end select

             if(DoTest .and. iBlock==iBlockTest)&
                  write(*,*) NameSub, ': TypeCrit, IsConserv=',&
                  TypeConservCrit_I(iCrit), &
                  IsConserv_CB(iTest,jTest,kTest,iBlock)

          end do ! iCrit

          if(UseElectronPressure)then
             do k = MinK,MaxK; do j = MinJ,MaxJ; do i = MinI,MaxI
                State_VGB(p_,i,j,k,iBlock) = State_VGB(p_,i,j,k,iBlock) &
                     - State_VGB(Pe_,i,j,k,iBlock)
             end do; end do; end do
          end if

       end do ! iBlock

    else
       ! If there are no physics based criteria we start from
       ! the assumption of conservative everywhere
       IsConserv_CB = .true.

       if(DoTest .and. iProc==iProcTest)&
            write(*,*) NameSub, ': default IsConserv is true'
    endif

    do iBlock = 1, nBlock
       if( Unused_B(iBlock) ) CYCLE

       ! Apply geometry based criteria
       ! Any of these can switch from conservative to non-conservative
       do iCrit = 1, nConservCrit
          select case(TypeConservCrit_I(iCrit))
          case('r')
             ! Switch to non-conservative inside radius rConserv
             IsConserv_CB(:,:,:,iBlock) = IsConserv_CB(:,:,:,iBlock) .and. &
                  R_BLK(1:nI,1:nJ,1:nK,iBlock) > rConserv
          case('parabola')
             ! Switch to non-conservative behind parabola inside the bow shock
             IsConserv_CB(:,:,:,iBlock) = IsConserv_CB(:,:,:,iBlock) .and. &
                  Xyz_DGB(x_,1:nI,1:nJ,1:nK,iBlock) > xParabolaConserv - &
                  ( Xyz_DGB(y_,1:nI,1:nJ,1:nK,iBlock)**2 &
                  + Xyz_DGB(z_,1:nI,1:nJ,1:nK,iBlock)**2 ) / yParabolaConserv
          case default
             CYCLE
          end select
          if(DoTest.and.iBlock==iBlockTest)&
               write(*,*) NameSub, ': TypeCrit, IsConserv=',&
               TypeConservCrit_I(iCrit), IsConserv_CB(iTest,jTest,kTest,iBlock)
       end do
    end do

    call timing_stop('nonconservative')

    call test_stop(NameSub, DoTest)
  end subroutine select_conservative
  !============================================================================

  subroutine fix_anisotropy

    ! Calculate the pressure anisotropy relaxation term for anisotropic MHD.
    ! Correct the parallel pressure based on the firehose, mirror and proton
    ! cyclotron instability criteria in unstable regions and the global
    ! relaxation, if present, in the whole domain.
    ! The instability that changes Ppar most is applied.
    !
    ! If UseConstantTau = true, use TauInstability read from PARAM.in as the
    ! contant relaxation time, same for different instabilities.
    ! TauGlobal is a constant read from PARAM.in.

    use ModVarIndexes, ONLY: Bx_, Bz_, Ppar_, p_
    use ModMain,    ONLY: nI, nJ, nK, nBlock, Unused_B, UseB0, &
         time_accurate, Cfl, dt
    use ModB0,      ONLY: B0_DGB
    use ModAdvance, ONLY: State_VGB, time_BLK
    use ModPhysics, ONLY: UseConstantTau, TauInstability, &
         IonMassPerCharge, TauGlobal
    use ModGeometry, ONLY: true_cell

    ! Variables for anisotropic pressure
    real:: B_D(3), B2, p, Ppar, Pperp, Dp, DtCell
    real:: InvGyroFreq, PparOverLimit, Deltapf, Deltapm

    integer:: i, j, k, iBlock

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'fix_anisotropy'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       do k=1,nK; do j=1,nJ; do i=1,nI
          if(.not.true_cell(i,j,k,iBlock)) CYCLE

          ! Avoid Pperp < 0
          State_VGB(Ppar_,i,j,k,iBlock) = &
               min(3*State_VGB(p_,i,j,k,iBlock),State_VGB(Ppar_,i,j,k,iBlock))

          ! Do not apply the relaxation term in this case
          if(UseConstantTau .and. TauInstability < 0.0) CYCLE

          B_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)
          if(UseB0) B_D = B_D + B0_DGB(:,i,j,k,iBlock)
          B2     = sum(B_D**2)
          Ppar   = State_VGB(Ppar_,i,j,k,iBlock)
          p = State_VGB(p_,i,j,k,iBlock)
          Pperp  = (3*p - Ppar)/2.
          if(.not. time_accurate)then
             DtCell = Cfl*time_BLK(i,j,k,iBlock)
          else
             DtCell = dt
          end if

          InvGyroFreq = IonMassPerCharge/max(1e-8, sqrt(B2))

          ! Find the instability that changes ppar the most
          Dp = 0.0

          ! This is what global relaxation would do
          if(TauGlobal > 0) Dp = DtCell*(p - Ppar)/(DtCell + TauGlobal)

          ! Check for the firehose, mirror and ion cyclotron instabilities
          ! Limit anisotropy to instability criteria in unstable regions
          if(Ppar - Pperp > B2)then
             ! firehose
             ! by how much the instability limit is exceeded:
             ! ppar - ppar_marginalstable
             PparOverLimit = Ppar - p - 2/3.*B2

             ! Calc firehose relaxation time based on the maximum
             ! growth rate calculated from eqn (2) of Hall [1981]
             ! with theta = 0 and ppar < 4*pperp
             ! MaxGrowthRate =
             !    0.5*GyroFreq*Delta pf/sqrt(ppar*(pperp-ppar/4))
             ! where Delta pf = ppar-pperp-B^2 = 3/2*PparOverLimit
             if(.not. UseConstantTau)then
                Deltapf = 3/2.*PparOverLimit
                TauInstability = 2.0*InvGyroFreq* &
                     sqrt(max(3.0*Ppar*(Pperp-0.25*Ppar),1e-8))/Deltapf
             end if
             Dp = min(Dp, -DtCell*PparOverLimit/(DtCell + TauInstability))

          else
             if(Pperp**2 > Ppar*Pperp + 0.5*B2*Ppar)then
                ! mirror
                ! ppar_marginalstable - ppar
                PparOverLimit = (B2 + 6.0*p &
                     - sqrt(B2**2 + 12.0*B2*p + 9.0*p**2))/3. - Ppar

                ! Calc mirror relaxation time based on the maximum
                ! growth rate from eqn (7) of Southwood [1993],
                ! with the wavelength at maximum growth from eqn (21)
                ! of Hall [1980]
                ! MaxGrowthRate =
                !    4/3/sqrt(5)*GyroFreq*sqrt(2*Delta pm/ppar)
                ! where Delta pm = pperp-ppar-B^2*ppar/(2*pperp)
                if(.not. UseConstantTau)then
                   Deltapm = Pperp - Ppar - 0.5*B2*Ppar/Pperp
                   TauInstability = 0.75*InvGyroFreq*sqrt(2.5*Ppar/Deltapm)
                end if
                Dp = max(Dp, DtCell*PparOverLimit/(DtCell + TauInstability))
             end if
             if(Pperp > Ppar + 0.3*sqrt(0.5*B2*Ppar))then
                ! ion cyclotron
                ! ppar_marginalstable - ppar
                PparOverLimit = (sqrt(0.01*B2 + 2.0*p) &
                     - 0.1*sqrt(B2))**2/2. - Ppar

                ! Estimate ion cyclotron relaxation time from
                ! observations in the magnetosphere and theories
                if(.not. UseConstantTau) &
                     TauInstability = 100*InvGyroFreq

                Dp = max(Dp, DtCell*PparOverLimit/(DtCell + TauInstability))
             end if
          end if
          State_VGB(Ppar_,i,j,k,iBlock) = Ppar + Dp
       end do; end do; end do
    end do

    call test_stop(NameSub, DoTest)
  end subroutine fix_anisotropy
  !============================================================================

  subroutine update_b0

    use ModMain,          ONLY: nBlock, Unused_B, UseNewMagnetogram, &
         time_simulation, NameThisComp, t_max, tMagnetogram, DoThreads_B
    use ModPhysics,       ONLY: ThetaTilt
    use ModAdvance,       ONLY: Bx_, By_, Bz_, State_VGB
    use ModGeometry,      ONLY: true_cell, body_BLK
    use CON_axes,         ONLY: get_axes
    use ModNumConst,      ONLY: cRadToDeg
    use ModIO,            ONLY: iUnitOut, write_prefix
    use ModEnergy,        ONLY: calc_energy_ghost
    use ModB0,            ONLY: B0_DGB, set_b0_cell, set_b0_reschange
    use ModFieldLineThread, ONLY: UseFieldLineThreads, set_threads
    use ModMagnetogram,   ONLY: update_magnetogram
    use ModMessagePass,   ONLY: exchange_messages

    integer :: iBlock

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'update_b0'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    ! Update ThetaTilt
    if(NameThisComp=='GM') &
         call get_axes(Time_Simulation, MagAxisTiltGsmOut=ThetaTilt)

    if (DoTest) then
       if(NameThisComp=='GM')then
          call write_prefix; write(iUnitOut,*) &
               "update_b0 at tSimulation, TiltGsm=", &
               Time_Simulation, ThetaTilt*cRadToDeg
       else
          call write_prefix; write(iUnitOut,*) &
               "update_b0 at tSimulation=",Time_Simulation
       end if
    end if
    call timing_start(NameSub)
    if(UseNewMagnetogram)&
         call update_magnetogram(time_simulation, t_max, tMagnetogram)
    do iBlock=1,nBlock
       if(Unused_B(iBlock)) CYCLE

       ! Save total magnetic field into Bx, By ,Bz
       State_VGB(Bx_:Bz_,:,:,:,iBlock) = State_VGB(Bx_:Bz_,:,:,:,iBlock) &
            + B0_DGB(:,:,:,:,iBlock)

       call set_b0_cell(iBlock)

       ! Split total B again using new B0
       State_VGB(Bx_:Bz_,:,:,:,iBlock) = State_VGB(Bx_:Bz_,:,:,:,iBlock) &
            - B0_DGB(:,:,:,:,iBlock)

       ! Set B1 to 0 inside bodies
       if(Body_BLK(iBlock))then
          where(.not.true_cell(:,:,:,iBlock))
             State_VGB(Bx_,:,:,:,iBlock)=0.0
             State_VGB(By_,:,:,:,iBlock)=0.0
             State_VGB(Bz_,:,:,:,iBlock)=0.0
          end where
       end if

       ! Recalculate energy
       call calc_energy_ghost(iBlock)

    end do

    ! Recalculate B0 face values at resolution changes
    call set_b0_reschange
    if(UseFieldLineThreads)then
       DoThreads_B = .true.
       call set_threads
       call exchange_messages
    end if
    call timing_stop(NameSub)

    call test_stop(NameSub, DoTest)
  end subroutine update_b0
  !============================================================================

  subroutine fix_multi_ion_update(iBlock)

    ! This subroutine sets the ion velocities of the individual fluids equal to
    ! their mass density weighted average, and/or the ion temperature of the
    ! individual fluids equal to their number density weighted average.

    use ModSize,       ONLY: nI, nJ, nK
    use ModAdvance,    ONLY: State_VGB, &
         UseSingleIonVelocity, UseSingleIonTemperature
    use ModMultiFluid, ONLY: nIonFluid, IonFirst_, IonLast_, &
         iRho_I, iRhoUx_I, iRhoUy_I, iRhoUz_I, iP_I, MassIon_I, MassFluid_I
    use ModGeometry,   ONLY: true_cell

    integer, intent(in) :: iBlock

    integer :: i, j, k, iFluid
    real :: RhoInv, Ux, Uy, Uz, Temp, NumDens_I(nIonFluid)

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'fix_multi_ion_update'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    if(DoTest)then
       write(*,*) NameSub,': ion fluid velocities and temperatures before fix:'
       do iFluid = IonFirst_, IonLast_
          write(*,*)'iIonFluid,ux,uy,uz,T=', &
               State_VGB(iRhoUx_I(iFluid):iRhoUz_I(iFluid),&
               iTest,jTest,kTest,iBlock) &
               / State_VGB(iRho_I(iFluid),iTest,jTest,kTest,iBlock), &
               State_VGB(iP_I(iFluid),iTest,jTest,kTest,iBlock)    &
               *MassFluid_I(iFluid) &
               /State_VGB(iRho_I(iFluid),iTest,jTest,kTest,iBlock)
       end do
    end if

    if(UseSingleIonVelocity) then
       do k=1, nK; do j=1, nJ; do i=1, nI
          if(.not.true_cell(i,j,k,iBlock)) CYCLE

          ! Calcualate average velocity from total momentum and density
          RhoInv= 1/sum(State_VGB(iRho_I(IonFirst_:IonLast_),i,j,k,iBlock))
          Ux = RhoInv*sum(State_VGB(iRhoUx_I(IonFirst_:IonLast_),i,j,k,iBlock))
          Uy = RhoInv*sum(State_VGB(iRhoUy_I(IonFirst_:IonLast_),i,j,k,iBlock))
          Uz = RhoInv*sum(State_VGB(iRhoUz_I(IonFirst_:IonLast_),i,j,k,iBlock))

          ! Reset the momentum of all ion fluids
          State_VGB(iRhoUx_I(IonFirst_:IonLast_),i,j,k,iBlock) = &
               Ux*State_VGB(iRho_I(IonFirst_:IonLast_),i,j,k,iBlock)
          State_VGB(iRhoUy_I(IonFirst_:IonLast_),i,j,k,iBlock) = &
               Uy*State_VGB(iRho_I(IonFirst_:IonLast_),i,j,k,iBlock)
          State_VGB(iRhoUz_I(IonFirst_:IonLast_),i,j,k,iBlock) = &
               Uz*State_VGB(iRho_I(IonFirst_:IonLast_),i,j,k,iBlock)

       end do; end do; end do
    end if

    if(UseSingleIonTemperature) then
       do k=1, nK; do j=1, nJ; do i=1, nI
          if(.not.true_cell(i,j,k,iBlock)) CYCLE

          ! Number density
          NumDens_I = &
               State_VGB(iRho_I(IonFirst_:IonLast_),i,j,k,iBlock)/MassIon_I

          ! Average temperature = sum(p)/sum(n) = sum(p)/sum(rho/M)
          Temp = sum(State_VGB(iP_I(IonFirst_:IonLast_),i,j,k,iBlock)) &
               / sum(NumDens_I)

          ! Reset the pressure of all ion fluids
          State_VGB(iP_I(IonFirst_:IonLast_),i,j,k,iBlock) = Temp*NumDens_I

       end do; end do; end do
    end if

    if(DoTest)then
       write(*,*) NameSub,': ion fluid velocities and temperatures after fix:'
       do iFluid = IonFirst_, IonLast_
          write(*,*)'iIonFluid,ux,uy,uz,T=', &
               State_VGB(iRhoUx_I(iFluid):iRhoUz_I(iFluid),&
               iTest,jTest,kTest,iBlock) &
               / State_VGB(iRho_I(iFluid),iTest,jTest,kTest,iBlock), &
               State_VGB(iP_I(iFluid),iTest,jTest,kTest,iBlock)    &
               *MassFluid_I(iFluid) &
               /State_VGB(iRho_I(iFluid),iTest,jTest,kTest,iBlock)
       end do
    end if

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine fix_multi_ion_update
  !============================================================================

end module ModUpdateState
!==============================================================================
