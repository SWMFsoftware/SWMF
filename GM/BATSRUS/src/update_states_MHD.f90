!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine update_states_MHD(iBlock)

  use ModProcMH
  use ModMain
  use ModAdvance
  use ModB0, ONLY: B0_DGB
  use ModPhysics
  use ModGeometry, ONLY: true_cell
  use ModKind, ONLY: Real8_
  use ModSemiImplVar, ONLY: UseStableImplicit
  use ModVarIndexes, ONLY: pe_, p_
  use ModPointImplicit, ONLY: UsePointImplicit, UsePointImplicit_B, &
       update_point_implicit
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
  use ModIonElectron, ONLY: ion_electron_source_impl, ion_electron_init_point_impl,&
       HypEDecay

  implicit none

  integer, intent(in) :: iBlock

  integer :: i, j, k

  ! These variables have to be double precision for accurate Boris scheme
  real:: DtLocal, DtFactor, Coeff
  logical :: DoTest, DoTestMe
  character(len=*), parameter :: NameSub = 'update_states_mhd'
  !--------------------------------------------------------------------------
  if(iBlock==BLKtest .and. iProc==PROCtest)then
     call set_oktest(NameSub,DoTest,DoTestMe)
  else
     DoTest=.false.; DoTestMe=.false.
  endif

  !\
  ! Update the new solution state and calc residuals for the mth stage.
  ! Note must copy state to old state only if m is 1.
  !/

  if(iStage==1) then
     StateOld_VCB(1:nVar,1:nI,1:nJ,1:nK,iBlock) = & 
          State_VGB(1:nVar,1:nI,1:nJ,1:nK,iBlock)
     EnergyOld_CBI(:,:,:,iBlock,:) = Energy_GBI(1:nI,1:nJ,1:nK,iBlock,:)
  end if

  ! Nothing to do if time step is zero
  if(time_accurate .and. Dt == 0.0) RETURN

  ! Add Joule heating: dPe/dt or dP/dt += (gamma-1)*eta*j**2
  ! also dE/dt += eta*j**2 for semi-implicit scheme (UseResistiveFlux=F)
  ! and heat exchange between electrons and ions (mult-ion is not coded).

  if(.not.UseMultiIon .and. UseResistivity .and. &
       (UseElectronPressure .or. UseNonConservative .or. &
       .not.UseResistiveFlux)) then

     call calc_resistivity_source(iBlock)   
     if(DoTestMe)write(*,*) NameSub, ' after add_resistive_source=', &
          State_VGB(VarTest,iTest,jTest,kTest,iBlock), &
          Energy_GBI(iTest,jTest,kTest,iBlock,:)
  end if

  !Get Residual.
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

  if(DoTestMe)write(*,*) NameSub, ' original testvar and energy=', &
       State_VGB(VarTest, iTest, jTest, kTest, iBlock), &
       Energy_GBI(iTest,jTest,kTest,iBlock,:)

  call update_explicit

  if(UseMultiIon .and. (UseSingleIonVelocity .or. UseSingleIonTemperature) &
       .and. .not.UsePointImplicit) then

     call fix_multi_ion_update(iBlock)
     call calc_energy_cell(iBlock)

     if(DoTestMe)write(*,*) NameSub, ' after fix multiion update=', &
          State_VGB(VarTest,iTest,jTest,kTest,iBlock), &
          Energy_GBI(iTest,jTest,kTest,iBlock,:)

  end if

  if(UseMultiIon .and. IsMhd)then
     call multi_ion_update(iBlock, IsFinal = .false.)

     if(DoTestMe)write(*,*) NameSub, ' after multiion update1=', &
          State_VGB(VarTest,iTest,jTest,kTest,iBlock), &
          Energy_GBI(iTest,jTest,kTest,iBlock,:)

  end if

  ! The point implicit update and other stuff below are only done in last stage
  ! except for ion-electron equations where the point implicit has to be done
  ! in every stage.
  if(.not.UseEfield .and. iStage < nStage)then
     if(UseBufferGrid) call fix_buffer_grid(iBlock)
     RETURN
  end if

  ! Add point implicit user or multi-ion source terms
  if (UsePointImplicit .and. UsePointImplicit_B(iBlock))then
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

     if(DoTestMe)write(*,*) NameSub, ' after point impl state=', &
          State_VGB(VarTest, iTest,jTest,kTest,iBlock), &
          Energy_GBI(iTest,jTest,kTest,iBlock,:)
  end if

  if(UseMultiIon .and. IsMhd)then
     call multi_ion_update(iBlock, IsFinal = .true.)
     if(DoTestMe)write(*,*) NameSub,' after multiion update2=', &
          State_VGB(VarTest,iTest,jTest,kTest,iBlock), &
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

  if(DoTestMe)write(*,*) NameSub, ' final state=', &
       State_VGB(VarTest,iTest,jTest,kTest,iBlock), &
       Energy_GBI(iTest,jTest,kTest,iBlock,:)

contains

  subroutine update_explicit

    use ModBorisCorrection, ONLY: mhd_to_boris, boris_to_mhd, &
         mhd_to_boris_simple, boris_simple_to_mhd

    ! Allocatable storage for classical 4th order Runge-Kutta scheme
    real, allocatable, save:: Rk4_VCB(:,:,:,:,:), Rk4_CBI(:,:,:,:,:)

    real, parameter:: cThird = 1./3.
    real:: Coeff1, Coeff2
    integer:: iFluid, iRho
    !--------------------------------------------------------------------------

    if(boris_correction) then
       ! Convert StateOld_VCB and State_VGB from classical MHD variables
       ! to semi-relativistic MHD variable
       do k=1,nK; do j=1,nJ; do i=1,nI
          if(.not.true_cell(i,j,k,iBlock)) CYCLE

          if(UseB0)then
             call mhd_to_boris(StateOld_VCB(:,i,j,k,iBlock), &
                  EnergyOld_CBI(i,j,k,iBlock,1), B0_DGB(:,i,j,k,iBlock))
             ! State_VGB is not used in 1-stage and HalfStep schemes
             if(.not.UseHalfStep .and. nStage > 1) &
                  call mhd_to_boris(State_VGB(:,i,j,k,iBlock), &
                  Energy_GBI(i,j,k,iBlock,1), B0_DGB(:,i,j,k,iBlock))
          else
             call mhd_to_boris(StateOld_VCB(:,i,j,k,iBlock), &
                  EnergyOld_CBI(i,j,k,iBlock,1))
             ! State_VGB is not used in 1-stage and HalfStep schemes
             if(.not.UseHalfStep .and. nStage > 1) &
                  call mhd_to_boris(State_VGB(:,i,j,k,iBlock), &
                  Energy_GBI(i,j,k,iBlock,1))
          end if
       end do; end do; end do
          
       if(DoTestMe)write(*,*) NameSub, ' after mhd_to_boris=', &
            State_VGB(VarTest,iTest,jTest,kTest,iBlock), &
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
             call mhd_to_boris_simple(StateOld_VCB(:,i,j,k,iBlock), &
                  B0_DGB(:,i,j,k,iBlock))
             ! State_VGB is not used in 1-stage and HalfStep schemes
             if(.not.UseHalfStep .and. nStage > 1) &
                  call mhd_to_boris_simple(State_VGB(:,i,j,k,iBlock), &
                  B0_DGB(:,i,j,k,iBlock))
          else
             call mhd_to_boris_simple(StateOld_VCB(:,i,j,k,iBlock))
             ! State_VGB is not used in 1-stage and HalfStep schemes
             if(.not.UseHalfStep .and. nStage > 1) &
                  call mhd_to_boris_simple(State_VGB(:,i,j,k,iBlock))
          end if
       end do; end do; end do
          
       if(DoTestMe)write(*,*) NameSub, ' after mhd_to_boris_simple=', &
            State_VGB(VarTest,iTest,jTest,kTest,iBlock), &
            Energy_GBI(iTest,jTest,kTest,iBlock,:)
    end if

    if(UseHalfStep .or. nStage == 1 .or. nStage == 4)then
       ! Update state variables starting from level n (=old) state
       do k=1,nK; do j=1,nJ; do i=1,nI
          State_VGB(:,i,j,k,iBlock) = &
               StateOld_VCB(:,i,j,k,iBlock) + Source_VC(1:nVar,i,j,k)
       end do; end do; end do

       ! Update energy variables
       do iFluid = 1, nFluid; do k=1,nK; do j=1,nJ; do i=1,nI
          Energy_GBI(i,j,k,iBlock,iFluid) = &
               EnergyOld_CBI(i,j,k,iBlock,iFluid)+Source_VC(nVar+iFluid,i,j,k)
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
                  - 4*StateOld_VCB(:,i,j,k,iBlock)
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
               Coeff1*StateOld_VCB(:,i,j,k,iBlock) + &
               Coeff2*State_VGB(:,i,j,k,iBlock)
       end do; end do; end do

       ! Interpolate energies
       do iFluid = 1, nFluid; do k=1,nK; do j=1,nJ; do i=1,nI
          Energy_GBI(i,j,k,iBlock,iFluid) = &
               Coeff1*EnergyOld_CBI(i,j,k,iBlock,iFluid) + &
               Coeff2*Energy_GBI(i,j,k,iBlock,iFluid)
       end do; end do; end do; end do

    endif

    if(DoTestMe)write(*,*) NameSub, ' after flux/source=', &
         State_VGB(VarTest,iTest,jTest,kTest,iBlock), &
         Energy_GBI(iTest,jTest,kTest,iBlock,:)

    if(boris_correction) then
       ! Convert StateOld_VCB and State_VGB back from 
       ! semi-relativistic to classical MHD variables
       do k=1,nK; do j=1,nJ; do i=1,nI
          if(.not.true_cell(i,j,k,iBlock)) CYCLE

          if(UseB0)then
             call boris_to_mhd(StateOld_VCB(:,i,j,k,iBlock), &
                  EnergyOld_CBI(i,j,k,iBlock,1), B0_DGB(:,i,j,k,iBlock))
             call boris_to_mhd(State_VGB(:,i,j,k,iBlock), &
                  Energy_GBI(i,j,k,iBlock,1), B0_DGB(:,i,j,k,iBlock))
          else
             call boris_to_mhd(StateOld_VCB(:,i,j,k,iBlock), &
                  EnergyOld_CBI(i,j,k,iBlock,1))
             call boris_to_mhd(State_VGB(:,i,j,k,iBlock), &
                  Energy_GBI(i,j,k,iBlock,1))
          end if
       end do; end do; end do
          
       if(DoTestMe)write(*,*) NameSub, ' after boris_to_mhd=', &
            State_VGB(VarTest,iTest,jTest,kTest,iBlock), &
            Energy_GBI(iTest,jTest,kTest,iBlock,:)
    endif

    if(UseBorisSimple) then
       ! Convert mometum in StateOld_VCB and State_VGB back from 
       ! enhanced momentum
       do k=1,nK; do j=1,nJ; do i=1,nI
          if(.not.true_cell(i,j,k,iBlock)) CYCLE

          if(UseB0)then
             call boris_simple_to_mhd(StateOld_VCB(:,i,j,k,iBlock), &
                  B0_DGB(:,i,j,k,iBlock))
             call boris_simple_to_mhd(State_VGB(:,i,j,k,iBlock), &
                  B0_DGB(:,i,j,k,iBlock))
          else
             call boris_simple_to_mhd(StateOld_VCB(:,i,j,k,iBlock))
             call boris_simple_to_mhd(State_VGB(:,i,j,k,iBlock))
          end if
       end do; end do; end do
          
       if(DoTestMe)write(*,*) NameSub, ' after boris_to_mhd=', &
            State_VGB(VarTest,iTest,jTest,kTest,iBlock), &
            Energy_GBI(iTest,jTest,kTest,iBlock,:)
    endif

    if(UseMultiSpecies)then
       ! Fix negative species densities
       State_VGB(SpeciesFirst_:SpeciesLast_,1:nI,1:nJ,1:nK,iBlock) = max(0.0,&
            State_VGB(SpeciesFirst_:SpeciesLast_,1:nI,1:nJ,1:nK,iBlock))

       if(DoReplaceDensity)then
          ! Add up species densities to total density
          do k=1,nK; do j=1,nJ; do i=1,nI
             State_VGB(Rho_,i,j,k,iBlock) = &
                  sum(State_VGB(SpeciesFirst_:SpeciesLast_,i,j,k,iBlock))
          end do; end do; end do
       end if

       if(DoTestMe)write(*,*) NameSub, ' after multispecies correct=', &
            State_VGB(VarTest,iTest,jTest,kTest,iBlock)

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

       if(DoTestMe)write(*,*) NameSub, &
            ' after min density correct densities=', &
            State_VGB(iRho_I,iTest,jTest,kTest,iBlock)
    end if

    if( IsMhd .and. &
         ((nStage==1.and..not.time_accurate) &
         .or.(nStage==2.and.iStage==1.and.UseHalfStep)))then

       ! A desparate attempt to maintain positivity by adding dB^2/2 to the 
       ! energy. This is fine for steady state, and is 2nd order accurate 
       ! for half+full step method. But it cannot be used for RK schemes!

       do k=1,nK; do j=1,nJ; do i=1,nI
          Energy_GBI(i,j,k,iBlock,1) = Energy_GBI(i,j,k,iBlock,1) + cHalf*( &
               Source_VC(Bx_,i,j,k)**2 + &
               Source_VC(By_,i,j,k)**2 + &
               Source_VC(Bz_,i,j,k)**2)
       end do; end do; end do

       if(DoTestMe)write(*,*) NameSub, ' after energy dB correct=', &
            State_VGB(VarTest,iTest,jTest,kTest,iBlock), &
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

!!! The electron entropy and the Boris corrections do not work for RK update!!!

    if(UseElectronPressure .and. UseElectronEntropy)then
       ! Fix update to use the electron entropy equation
       ! Se = Pe^(1/GammaE), dSe/Pe = 1/GammaE*Se/Pe
       ! Pe_n+1 = (Se_n+1)^GammaE = (Se_n + R(Se))^GammaE 
       !        = (Pe_n^(1/GammaE + R(Se))^GammaE
       do k=1,nK; do j=1,nJ; do i=1,nI
          if(.not.true_cell(i,j,k,iBlock)) CYCLE

          State_VGB(Pe_,i,j,k,iBlock) = &
               (StateOld_VCB(Pe_,i,j,k,iBlock)**(1/GammaElectron) &
               + Source_VC(Pe_,i,j,k))**GammaElectron
       end do; end do; end do
    end if

    ! Update energy or pressure based on UseConservative and IsConserv_CB
    call calc_energy_or_pressure(iBlock)

    if(DoTestMe)write(*,*) NameSub, ' after pressure/energy update=', &
         State_VGB(VarTest,iTest,jTest,kTest,iBlock), &
         Energy_GBI(iTest,jTest,kTest,iBlock,:)

  end subroutine update_explicit
  !----------------------------------------------------------------------

  subroutine deduct_expl_source()
    integer:: iVarSemi_
    character(len=*), parameter :: NameSub = 'deduct_expl_source'

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
            StateOld_VCB(iVarSemi_,i,j,k,iBlock)
       State_VGB(iVarSemi_,i,j,k,iBlock) = &
            StateOld_VCB(iVarSemi_,i,j,k,iBlock)
    end do; end do; end do
  end subroutine deduct_expl_source
  !----------------------------------------------------------------------

end subroutine update_states_mhd

!==============================================================================

subroutine fix_anisotropy

  !\
  ! Calculate the pressure anisotropy relaxation term for anisotropic MHD.
  ! Correct the parallel pressure based on the firehose, mirror and proton 
  ! cyclotron instability criteria in unstable regions and the global relaxation, 
  ! if present, in the whole domain. The one that changes Ppar most is applied.
  !
  ! If UseConstantTau = true, use TauInstability read from PARAM.in as the 
  ! contant relaxation time, same for different instabilities. 
  ! TauGlobal is a constant read from PARAM.in.
  !/

  use ModVarIndexes, ONLY: Bx_, Bz_, Ppar_, p_
  use ModMain,    ONLY: nI, nJ, nK, nBlock, Unused_B, UseB0, &
       time_accurate, Cfl, dt
  use ModB0,      ONLY: B0_DGB
  use ModAdvance, ONLY: State_VGB, time_BLK
  use ModPhysics, ONLY: UseConstantTau, TauInstability, IonMassPerCharge, TauGlobal
  use ModGeometry,ONLY: true_cell

  implicit none

  ! Variables for anisotropic pressure
  real:: B_D(3), B2, p, Ppar, Pperp, Dp, DtCell
  real:: InvGyroFreq, PparOverLimit, Deltapf, Deltapm

  integer:: i, j, k, iBlock
  !---------------------------------------------------------------------------
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
  use ModMessagePass, ONLY: exchange_messages
  implicit none

  character(len=*), parameter :: NameSub = 'update_b0'
  logical :: DoTest, DoTestMe
  integer :: iBlock
  !==========================================================================

  call set_oktest(NameSub, DoTest, DoTestMe)

  ! Update ThetaTilt
  if(NameThisComp=='GM') &
       call get_axes(Time_Simulation, MagAxisTiltGsmOut=ThetaTilt)

  if (DoTestMe) then
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

end subroutine update_b0

!===========================================================================

subroutine fix_multi_ion_update(iBlock)

  ! This subroutine sets the ion velocities of the individual fluids equal to
  ! their mass density weighted average, and/or the ion temperature of the 
  ! individual fluids equal to their number density weighted average.

  use ModSize,       ONLY: nI, nJ, nK
  use ModAdvance,    ONLY: State_VGB, &
       UseSingleIonVelocity, UseSingleIonTemperature
  use ModMain,       ONLY: iTest, jTest, kTest, BlkTest, ProcTest
  use ModMultiFluid, ONLY: nIonFluid, IonFirst_, IonLast_, &
       iRho_I, iRhoUx_I, iRhoUy_I, iRhoUz_I, iP_I, MassIon_I, MassFluid_I
  use ModGeometry,   ONLY: true_cell
  use ModProcMH,     ONLY: iProc

  integer, intent(in) :: iBlock

  integer :: i, j, k, iFluid
  real :: RhoInv, Ux, Uy, Uz, Temp, NumDens_I(nIonFluid)

  logical :: DoTest, DoTestMe
  character(len=*), parameter :: NameSub = 'fix_multi_ion_update'
  !----------------------------------------------------------------------
  if (iBlock == BlkTest .and. iProc == ProcTest) then 
     call set_oktest(NameSub, DoTest, DoTestMe)
  else
     DoTest = .false. ; DoTestMe = .false.
  end if

  if(DoTestMe)then
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

  if(DoTestMe)then
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

end subroutine fix_multi_ion_update

