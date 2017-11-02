!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModTimeStepControl

  implicit none

  SAVE

  private ! except

  ! Public methods
  public:: read_time_step_control_param
  public:: calc_timestep
  public:: set_global_timestep
  public:: control_time_step

  logical, public:: UseTimeStepControl   = .false.
  real,    public:: TimeStepControlInit  = 1.0
  logical, public:: UseMaxTimeStep      = .false. ! calculate max time step
  real,    public:: DtMin = -1.0, DtMax = -1.0    ! value of min/max time steps

  ! Local variables
  integer             :: nVarControl = -1
  integer, allocatable:: iVarControl_I(:)
  real,    allocatable:: VarRatio_I(:)

  real :: &
       RejectStepLevel1   = 0.3, &
       RejectStepLevel2   = 3.0, &
       RejectStepFactor   = 0.50, &
       ReduceStepLevel1   = 0.6, &
       ReduceStepLevel2   = 1.5, &
       ReduceStepFactor   = 0.95, &
       IncreaseStepLevel1 = 0.8, &
       IncreaseStepLevel2 = 1.2, &
       IncreaseStepFactor = 1.05

  ! Minimum time step checking
  logical:: DoCheckTimeStep = .false.
  integer:: nCheckTimeStep = 0, iCheckTimeStep = 0
  real   :: TimeStepMin = 0.0

contains

  !===========================================================================
  subroutine read_time_step_control_param(NameCommand)

    use ModMain,       ONLY: NameVarLower_V
    use ModReadParam
    use ModVarIndexes, ONLY: nVar
    use ModUtilities,  ONLY: split_string

    character(len=*), intent(in):: NameCommand

    integer :: iControl, iVar
    character(len=100):: NameVarControl
    character(len=20) :: NameVarControl_I(nVar)

    character(len=*), parameter:: NameSub='read_time_step_control_param'
    !------------------------------------------------------------------------
    select case(NameCommand)
    case("#CONTROLTIMESTEP", "#TIMESTEPCONTROL")
       call read_var('UseTimeStepControl', UseTimeStepControl)
    case("#CONTROLINIT")
       call read_var('TimeStepControlInit', TimeStepControlInit)
    case("#CONTROLDECREASE")
       call read_var('RejectStepLevel1' ,  RejectStepLevel1)
       call read_var('ReduceStepLevel1' ,  ReduceStepLevel1)
       call read_var('IncreaseStepLevel1', IncreaseStepLevel1)
    case("#CONTROLINCREASE")
       call read_var('RejectStepLevel2' ,  RejectStepLevel2)
       call read_var('ReduceStepLevel2' ,  ReduceStepLevel2)
       call read_var('IncreaseStepLevel2', IncreaseStepLevel2)
    case("#CONTROLFACTOR")
       call read_var('RejectStepFactor',   RejectStepFactor)
       call read_var('ReduceStepFactor',   ReduceStepFactor)
       call read_var('IncreaseStepFactor', IncreaseStepFactor)
    case("#CONTROLVAR")
       call read_var('NameVarControl', NameVarControl, IsLowerCase=.true.)
       ! Split list of variables and match with NameVar_V
       call split_string(NameVarControl, nVar, NameVarControl_I, nVarControl)
       if(allocated(iVarControl_I)) deallocate(iVarControl_I, VarRatio_I)
       allocate( iVarControl_I(nVarControl), VarRatio_I(nVarControl) )
       do iControl=1, nVarControl
          do iVar = 1, nVar
             if(NameVarLower_V(iVar) == NameVarControl_I(iControl)) EXIT
          end do
          if(iVar > nVar)call stop_mpi(NameSub//' invalid NameVarControl='//&
               NameVarControl_I(iControl))
          iVarControl_I(iControl) = iVar
       end do
    case("#CHECKTIMESTEP")
       call read_var('DoCheckTimeStep', DoCheckTimeStep)
       call read_var('nCheckTimeStep',  nCheckTimeStep)
       call read_var('TimeStepMin'   ,  TimeStepMin)
    case default
       call stop_mpi(NameSub//' invalid NameCommand='//NameCommand)
    end select

  end subroutine read_time_step_control_param

  !===========================================================================

  subroutine calc_timestep(iBlock)

    use ModVarIndexes, ONLY: p_, WaveFirst_, WaveLast_
    use ModSize, ONLY: nI, nJ, nK
    use ModMain, ONLY: UseDtFixed, DtFixed, Dt_BLK, Cfl, &
         iTest, jTest, kTest, BlkTest, UseDtLimit, DtLimit
    use ModAdvance, ONLY : VdtFace_x, VdtFace_y, VdtFace_z, time_BLK, &
         DoFixAxis, rFixAxis, r2FixAxis, State_VGB, &
         UseElectronPressure
    use ModGeometry, ONLY: true_cell, true_BLK, rMin_BLK
    use ModCoronalHeating, ONLY: UseCoronalHeating, get_block_heating, &
         CoronalHeating_C, UseAlfvenWaveDissipation, WaveDissipation_VC
    use ModRadiativeCooling, ONLY: UseRadCooling, &
         get_radiative_cooling, add_chromosphere_heating
    use ModChromosphere, ONLY: DoExtendTransitionRegion, extension_factor, &
         UseChromosphereHeating, get_tesi_c, TeSi_C
    use ModPhysics, ONLY: InvGammaMinus1
    use BATL_lib, ONLY: CellVolume_GB, CoordMin_DB, CoordMax_DB, &
         IsCylindricalAxis, IsLatitudeAxis, r_, Lat_
    use ModNumConst, ONLY: cHalfPi
    use ModCoarseAxis,ONLY: UseCoarseAxis, calc_coarse_axis_timestep,&
         NorthHemiSph_, SouthHemiSph_

    integer, intent(in) :: iBlock

    integer :: i, j, k, Di, Dk
    real:: Vdt

    ! Variables for time step control due to loss terms
    real :: Einternal, Source, Dt_loss, Coef

    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'calc_timestep'
    !--------------------------------------------------------------------------

    if(iBlock==BLKtest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    endif

    if(DoTestMe)write(*,*) NameSub,' starting, true_BLK=', true_BLK(iBlock)

    ! Calculate time step limit based on maximum speeds across 6 faces
    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(.not. true_cell(i,j,k,iBlock)) then
          time_BLK(i,j,k,iBlock) = 0
       else
          Vdt = max(VdtFace_x(i,j,k),VdtFace_x(i+1,j,k))
          if(nJ > 1) Vdt = Vdt + max(VdtFace_y(i,j,k), VdtFace_y(i,j+1,k))
          if(nK > 1) Vdt = Vdt + max(VdtFace_z(i,j,k), VdtFace_z(i,j,k+1))
          time_BLK(i,j,k,iBlock) = CellVolume_GB(i,j,k,iBlock) / Vdt
       end if
    end do; end do; end do

    if(DoFixAxis)then
       ! Use the time step in the super cell of the cell just outside.
       ! In time accurate this removes the time step constraints from supercell
       ! In local time stepping mode it increases the time step
       if(IsCylindricalAxis)then
          if(CoordMin_DB(r_,iBlock) <= 0.0)then
             Di = 1; if(r2FixAxis > 0.0) Di = 2
             do j = 1, nJ; do i = 1, nI
                time_BLK(1:Di,j,k,iBlock) = time_BLK(Di+1,j,k,iBlock)
             end do; end do
          end if
       elseif(IsLatitudeAxis .and. rMin_Blk(iBlock) < rFixAxis)then
          Dk = 1; if(rMin_Blk(iBlock) < r2FixAxis) Dk = 2
          if(CoordMax_DB(Lat_,iBlock) > cHalfPi-1e-8)then
             do k = nK+1-Dk, nK
                time_BLK(1:nI,1:nJ,k,iBlock) = time_BLK(1:nI,1:nJ,nK-Dk,iBlock)
             end do
          end if
          if(CoordMin_DB(Lat_,iBlock) < -cHalfPi+1e-8)then
             do k = 1, Dk
                time_BLK(1:nI,1:nJ,k,iBlock) = time_BLK(1:nI,1:nJ,Dk+1,iBlock)
             end do
          end if
       end if
    elseif(UseCoarseAxis)then
       if(CoordMax_DB(Lat_,iBlock) > cHalfPi-1e-8)then
          call calc_coarse_axis_timestep(iBlock,NorthHemiSph_)
       elseif(CoordMin_DB(Lat_,iBlock) < -cHalfPi+1e-8)then
          call calc_coarse_axis_timestep(iBlock,SouthHemiSph_)
       end if
    end if

    ! Time step restriction due to point-wise loss terms
    ! (only explicit source terms)
    if(UseAlfvenWaveDissipation .or.UseRadCooling)then
       if(UseRadCooling .or. DoExtendTransitionRegion) &
            call get_tesi_c(iBlock, TeSi_C)

       if(UseCoronalHeating)then
          call get_block_heating(iBlock)

          if(UseChromosphereHeating .and. DoExtendTransitionRegion)then
             call add_chromosphere_heating(TeSi_C, iBlock)
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                CoronalHeating_C(i,j,k) = &
                     CoronalHeating_C(i,j,k)/extension_factor(TeSi_C(i,j,k))
             end do; end do; end do
          end if
       end if

       if(UseAlfvenWaveDissipation)then
          if(DoExtendTransitionRegion)then
             ! Does not work together with UseChromosphereHeating
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                Coef = extension_factor(TeSi_C(i,j,k))
                WaveDissipation_VC(:,i,j,k) = WaveDissipation_VC(:,i,j,k)/Coef
                CoronalHeating_C(i,j,k) = CoronalHeating_C(i,j,k)/Coef
             end do; end do; end do
          end if

          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not. true_cell(i,j,k,iBlock)) CYCLE

             Dt_loss = 0.5*minval(State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock)&
                  /WaveDissipation_VC(:,i,j,k))
             ! The following prevents the wave energies from becoming negative
             ! due to too large loss terms.
             time_BLK(i,j,k,iBlock) = min(time_BLK(i,j,k,iBlock), Dt_loss)
          end do; end do; end do
       end if

       ! No time step restriction yet if the electron pressure is used
       if(UseRadCooling .and. .not.UseElectronPressure)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             call get_radiative_cooling(i, j, k, iBlock, TeSi_C(i,j,k), Source)

             if(UseCoronalHeating) Source = Source + CoronalHeating_C(i,j,k)

             ! Only limit for losses
             if(Source >= -1e-30) CYCLE

             Einternal = InvGammaMinus1*State_VGB(p_,i,j,k,iBlock)

             Dt_loss = 0.5*Einternal/abs(Source)
             ! The following prevents the pressure from becoming negative
             ! due to too large loss terms.
             time_BLK(i,j,k,iBlock) = min(time_BLK(i,j,k,iBlock), Dt_loss)
          end do; end do; end do
       end if
    end if

    if(DoTestMe)then
       write(*,*)NameSub,' VdtFace_x(iTest:iTest+1)=', &
            VdtFace_x(Itest:Itest+1,jTest,kTest)
       if(nJ>1) write(*,*) NameSub,' VdtFace_y(jTest:jTest+1)=', &
            VdtFace_y(Itest,Jtest:jTest+1,kTest)
       if(nK>1) write(*,*) NameSub,' VdtFace_z(kTest:kTest+1)=', &
            VdtFace_z(Itest,Jtest,kTest:kTest+1)
       write(*,*) NameSub,' time_BLK=',time_BLK(Itest,Jtest,Ktest,iBlock)
    end if

    !\
    ! Compute maximum stable time step for this solution block 
    !/

    if(true_BLK(iBlock)) then
       Dt_BLK(iBlock) = minval(time_BLK(:,:,:,iBlock))
    else
       ! If the block has no true cells, set Dt_BLK=1.0E20
       Dt_BLK(iBlock) = min(1e20, &
            minval(time_BLK(:,:,:,iBlock), &
            MASK=true_cell(1:nI,1:nJ,1:nK,iBlock)))

       if(DoTestMe)write(*,*) NameSub,' minval(time_BLK,MASK), Dt_BLK=',&
            minval(time_BLK(:,:,:,iBlock), &
            MASK=true_cell(1:nI,1:nJ,1:nK,iBlock)), Dt_BLK(iBlock)
    end if

    if(DoTestMe)write(*,*)NameSub,' Dt_BLK, loc=',Dt_BLK(iBlock),&
         minloc(time_BLK(:,:,:,iBlock),&
         MASK=true_cell(1:nI,1:nJ,1:nK,iBlock))

    ! Reset time_BLK for fixed time step (but Dt_BLK is kept!)
    if(UseDtFixed) &
         time_BLK(:,:,:,iBlock) = DtFixed

    ! Limit local time step so that Cfl*time_BLK <= DtLimit,
    if(UseDtLimit) &
         time_BLK(:,:,:,iBlock) = min(DtLimit/Cfl, time_BLK(:,:,:,iBlock))

    if(DoTestMe .and. UseDtFixed) &
         write(*,*) NameSub,' after UseDtFixed, time_BLK =', &
         time_BLK(Itest,Jtest,Ktest,iBlock)

    if(DoTestMe .and. UseDtLimit) &
         write(*,*) NameSub,' after limiting, time_BLK =', &
         time_BLK(Itest,Jtest,Ktest,iBlock)

    ! Set time step to zero inside body.
    if(.not.true_BLK(iBlock)) then
       where (.not.true_cell(1:nI,1:nJ,1:nK,iBlock))&
            time_BLK(:,:,:,iBlock) = 0.0
    end if

  end subroutine calc_timestep

  !============================================================================

  subroutine set_global_timestep(TimeSimulationLimit)

    use ModProcMH
    use ModMain
    use ModAdvance,  ONLY: time_BLK, State_VGB, rho_, Bx_, Bz_, P_, &
         iTypeAdvance_B, ExplBlock_
    use ModB0,       ONLY: B0_DGB
    use ModGeometry, ONLY: true_cell, true_BLK
    use ModImplicit, ONLY: UsePartImplicit
    use ModPhysics,  ONLY: No2Si_V, Si2No_V, No2Io_V, &
         UnitX_, UnitU_, UnitT_, UnitB_, UnitRho_, UnitP_, Gamma
    use ModNumConst
    use ModMpi
    use BATL_lib,    ONLY: Xyz_DGB, CellSize_DB, &
         MaxNode, nNode, iNode_B, iTimeLevel_A, nTimeLevel

    real, intent(in) :: TimeSimulationLimit !Simulation time not to be exceeded

    integer :: iBlock
    integer :: iError, Ijk_D(3), i, j, k
    real    :: DtMinPe, Cmax, Cmax_C(nI,nJ,nK)

    real :: DtDim

    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'set_global_timestep'
    !--------------------------------------------------------------------------
    ! time_BLK is already set in calc_timestep, 
    ! and Dt=DtLimit is set in MH_set_parameters
    if(UseDtLimit) RETURN

    call set_oktest(NameSub, DoTest, DoTestMe)

    ! if calc_timestep is tested, test this routine too
    if(.not.DoTest) call set_oktest('calc_timestep', DoTest, DoTestMe)

    if(DoTestMe)write(*,*) NameSub, &
         ' starting with TimeSimulationLimit, dt_BLK, time_BLK=', &
         TimeSimulationLimit, dt_BLK(BlkTest), &
         time_BLK(iTest,jTest,kTest,BlkTest)

    if(UseMaxTimeStep)then
       if(.not.allocated(iTimeLevel_A)) allocate(iTimeLevel_A(MaxNode))
       iTimeLevel_A(1:nNode) = 0
       nTimeLevel = 0
    end if

    if(UseDtFixed)then
       Dt = DtFixed
    elseif(n_step < 1 .or. n_step == 1 .and. TimeSimulationLimit > 0.0)then
       Dt    = 0.0
       DtMin = 0.0
       DtMax = 0.0
    else
       !\
       ! Impose global time step for time-accurate calculations as required
       !/
       if(UsePartImplicit)then
          ! Implicit blocks are not taken into account for partially implicit
          ! run
          DtMinPE = minval(Dt_BLK(1:nBlock),&
               MASK=iTypeAdvance_B(1:nBlock) == ExplBlock_)
          if(UseMaxTimeStep) DtMax = maxval(Dt_BLK(1:nBlock),&
               MASK=iTypeAdvance_B(1:nBlock) == ExplBlock_)
       else
          DtMinPE = minval(Dt_BLK(1:nBlock), MASK=.not.Unused_B(1:nBlock))
          if(UseMaxTimeStep) DtMax = min(DtLimit/Cfl, &
               maxval(Dt_BLK(1:nBlock), MASK=.not.Unused_B(1:nBlock)))
       end if

       ! Set Dt to minimum time step over all the PE-s
       call MPI_allreduce(DtMinPE, DtMin, 1, MPI_REAL, MPI_MIN, iComm, iError)
       Dt = DtMin

       if(DoTestMe .and. DtMinPE == Dt)then
          do iBlock = 1, nBlock
             if(Unused_B(iBlock)) CYCLE
             if(Dt_BLK(iBlock) /= Dt) CYCLE
             write(*,*) NameSub, ' Dt=',Dt,'=', Dt*No2Si_V(UnitT_),&
                  ' s  is controlled by block with ',&
                  'iProc, iBlock, true_BLK, true_cell =', &
                  iProc, iBlock, true_BLK(iBlock), &
                  all(true_cell(1:nI,1:nJ,1:nK,iBlock)), &
		  any(true_cell(1:nI,1:nJ,1:nK,iBlock))
             write(*,*) NameSub, ' X,Y,Z coordinates of (1,1,1) cell center=',&
                  Xyz_DGB(:,1,1,1,iBlock)
             write(*,*) NameSub, ' cell size in normalized and SI units:', &
                  CellSize_DB(:,iBlock), ', ', &
                  CellSize_DB(:,iBlock)*No2Si_V(UnitX_),' m'
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                Cmax_C(i,j,k) = Gamma*State_VGB(P_,i,j,k,iBlock)
                if(UseB)then
                   if(UseB0)then
                      Cmax_C(i,j,k) = Cmax_C(i,j,k) +           &
                           sum((State_VGB(Bx_:Bz_,i,j,k,iBlock) &
                           +    B0_DGB(:,i,j,k,iBlock))**2    )
                   else
                      Cmax_C(i,j,k) = Cmax_C(i,j,k) + &
                           sum(State_VGB(Bx_:Bz_,i,j,k,iBlock)**2) 
                   end if
                end if
                CMax_C(i,j,k) = Cmax_C(i,j,k)/State_VGB(rho_,i,j,k,iBlock)
             end do; end do; end do
             Ijk_D = maxloc(Cmax_C, MASK=true_cell(1:nI,1:nJ,1:nK,iBlock))
             i=Ijk_D(1); j=Ijk_D(2); k=Ijk_D(3)
             Cmax = sqrt(Cmax_C(i,j,k))
             write(*,*) NameSub, ' Cmax=',Cmax*No2Si_V(UnitU_),&
                  ' m/s is reached at X,Y,Z=', Xyz_DGB(:,i,j,k,iBlock)
             if(UseB0)write(*,*) NameSub,' B0=',&
                  B0_DGB(:,i,j,k,iBlock)*No2Si_V(UnitB_), ' T'
             if(UseB)write(*,*) NameSub, ' B1=',&
                  State_VGB(Bx_:Bz_,i,j,k,iBlock)*No2Si_V(UnitB_), ' T'
             write(*,*) NameSub,' Rho=',&
                  State_VGB(Rho_,i,j,k,iBlock)*No2Si_V(UnitRho_), ' kg/m3'
             write(*,*) NameSub,' p=',&
                  State_VGB(p_,i,j,k,iBlock)*No2Si_V(UnitP_), ' Pa'
             EXIT
          end do
       end if

       if(UseMaxTimeStep)then
          ! Calculate largest time step
          if(nProc > 1) call MPI_allreduce(&
               MPI_IN_PLACE, DtMax, 1, MPI_REAL, MPI_MAX, iComm, iError)

          if(DoTestMe) write(*,*)NameSub,' DtMin, DtMax=', DtMin, DtMax

          ! Make DtMax a power of 2 multiple of DtMin
          nTimeLevel = floor(log(DtMax/DtMin)/log(2.0))
          DtMax = DtMin * 2**nTimeLevel

          if(DoTestMe)write(*,*)NameSub,' DtMax after rounding=',DtMax

          Dt = DtMax
       end if

    end if

    ! Limit Dt such that the simulation time cannot exceed TimeSimulationLimit.
    ! If statement avoids real overflow when TimeSimulationLimit = Huge(0.0)
    if(TimeSimulationLimit > 0.0 .and. &
         Time_Simulation + Cfl*Dt*No2Si_V(UnitT_) > TimeSimulationLimit)then
       Dt = (TimeSimulationLimit - Time_Simulation)*Si2No_V(UnitT_)/Cfl
       if(UseMaxTimeStep)then
          DtMax = Dt
          DtMin = min(DtMin, Dt)
          nTimeLevel = ceiling(log(DtMax/DtMin)/log(2.0))
          DtMin = DtMax/2**nTimeLevel
       end if
    end if

    do iBlock = 1, nBlock
       if (Unused_B(iBlock)) CYCLE

       if(UseMaxTimeStep .and. Dt_BLK(iBlock) > 0)then
          ! Make block time step power of 2 multiple of DtMin that is < Dt_BLK
          ! Limit by DtMax in case DtMax was limited by TimeSimulationLimit
          if(Dt_BLK(iBlock) >= DtMax)then
             Dt_BLK(iBlock) = DtMax
          else
             ! Time level of this block
             iTimeLevel_A(iNode_B(iBlock)) = &
                  ceiling(log(DtMax/Dt_BLK(iBlock))/log(2.0))
             ! Time step rounded to power of 2 fraction of DtMax
             Dt_BLK(iBlock) = DtMax / 2**iTimeLevel_A(iNode_B(iBlock))
          end if
          time_BLK(:,:,:,iBlock) = Dt_BLK(iBlock)
       else
          time_BLK(:,:,:,iBlock) = Dt
       end if

       ! Reset time step to zero inside body.
       if(.not.true_BLK(iBlock))then
          where(.not.true_cell(1:nI,1:nJ,1:nK,iBlock)) &
               time_BLK(:,:,:,iBlock) = 0.0
       end if

    end do

    ! Collect time level information from all processors
    if(UseMaxTimeStep .and. nProc > 1) call MPI_allreduce(MPI_IN_PLACE, &
         iTimeLevel_A, nNode, MPI_INTEGER, MPI_SUM, iComm, iError)

    ! Set global time step to the actual time step used
    Dt = Cfl*Dt

    ! Check if time step has been too small for nCheckTimeStep steps
    if(DoTestMe)write(*,*) NameSub,' DoCheckTimeStep=', DoCheckTimeStep

    if(DoCheckTimeStep)then
       DtDim = Dt*No2Io_V(UnitT_)
       if(DoTestMe)write(*,*) NameSub,' checking DtDim=', DtDim
       if(DtDim < TimeStepMin)then
          iCheckTimeStep = iCheckTimeStep + 1
          if(DoTestMe)write(*,*) NameSub,' iCheckTimeStep=', iCheckTimeStep
          if(iCheckTimeStep >= nCheckTimeStep)then
             if(iProc==0)then
                write(*,*) NameSub,' ERROR: time step=', &
                     DtDim,' is smaller than TimeStepMin=', TimeStepMin, &
                     ' for the last ',iCheckTimeStep, ' time steps'
                write(*,*) NameSub,' saving output files before abort'
             end if
             ! call BATS_save_files('FINAL')
             if(iProc==0) call stop_mpi(NameSub//' time steps are too small')
          end if
       else
          ! Reset counter for small time steps
          iCheckTimeStep = 0
       end if
    end if

    if(DoTestMe)write(*,*) NameSub,' finished with Dt, dt_BLK, time_BLK=', &
         Dt, dt_BLK(BlkTest), time_BLK(iTest,jTest,kTest,BlkTest)

  end subroutine set_global_timestep

  !===========================================================================

  subroutine control_time_step

    use ModMain,     ONLY: nBlock, nI, nJ, nK, Unused_B, Dt, Cfl, CflOrig, &
         DtFixed, DtFixedOrig, UseDtFixed, Time_Simulation, &
         DtLimit, DtLimitOrig, UseDtLimit, UseLocalTimeStep
    use ModAdvance,  ONLY: Rho_, p_, &
         State_VGB, StateOld_VGB, Energy_GBI, EnergyOld_CBI, time_BLK
    use ModPhysics,  ONLY: No2Si_V, UnitT_
    use ModGeometry, ONLY: true_cell
    use ModProcMH,   ONLY: nProc, iComm
    use ModMpi

    integer:: iBlock, i, j, k, iError
    real   :: RelativeChangeMin,  RelativeChangeMax, Tmp, Factor

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub='control_time_step'
    !-------------------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)

    ! Use density and pressure as control variables by default
    if(nVarControl < 0)then
       nVarControl = 2
       allocate( iVarControl_I(nVarControl), VarRatio_I(nVarControl) )
       iVarControl_I = (/Rho_, p_/)
    end if

    ! Calculate the largest relative drop in the control variables
    RelativeChangeMin = 1e+6
    RelativeChangeMax = 1e-6

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(true_cell(i,j,k,iBlock))then
             VarRatio_I = State_VGB(iVarControl_I,i,j,k,iBlock) &
                  /    StateOld_VGB(iVarControl_I,i,j,k,iBlock)
             RelativeChangeMin = min(RelativeChangeMin, minval(VarRatio_I))
             RelativeChangeMax = max(RelativeChangeMax, maxval(VarRatio_I))
          end if
       end do; end do; end do
    end do

    if(nProc > 1)then
       ! Take minimum and maximum over all the processors
       Tmp = RelativeChangeMin
       call MPI_allreduce(Tmp, RelativeChangeMin, 1, MPI_REAL, MPI_MIN, &
            iComm, iError)
       Tmp = RelativeChangeMax
       call MPI_allreduce(Tmp, RelativeChangeMax, 1, MPI_REAL, MPI_MAX, &
            iComm, iError)
    end if

    ! Figure out time step reduction factor
    if(       RelativeChangeMin < RejectStepLevel1  &
         .or. RelativeChangeMax > RejectStepLevel2 )then
       ! Redo step if change is outside [RejectStepLevel1, RejectStepLevel2]
       Time_Simulation = Time_Simulation - Dt*No2Si_V(UnitT_)
       Dt = 0.0
       ! Do not use previous step in BDF2 scheme
!!! n_prev = -1
       ! Reset the state variable, the energy and set time_BLK variable to 0
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          State_VGB(:,1:nI,1:nJ,1:nK,iBlock)  &
               = StateOld_VGB(:,1:nI,1:nJ,1:nK,iBlock)
          Energy_GBI(1:nI,1:nJ,1:nK,iBlock,:) = EnergyOld_CBI(:,:,:,iBlock,:)
          time_BLK(1:nI,1:nJ,1:nK,iBlock)     = 0.0
       end do
       ! Reduce next time step
       Factor = RejectStepFactor

       !if (RelativeChangeMax > 100.0) then
       !   Factor = 0.1
       !endif

    elseif(   RelativeChangeMin < ReduceStepLevel1 &
         .or. RelativeChangeMax > ReduceStepLevel2 )then
       ! Reduce next time step if change exceeded ReduceStepLevel
       Factor = ReduceStepFactor

    elseif(    RelativeChangeMin > IncreaseStepLevel1  &
         .and. RelativeChangeMax < IncreaseStepLevel2 )then
       ! Increase next time step if change is within IncreaseStepLevel
       Factor = IncreaseStepFactor
    else
       Factor = 1.0
    end if

    if(UseDtFixed)then
       ! Do not exceed DtFixedOrig
       DtFixed = min(DtFixedOrig, DtFixed*Factor)
    elseif(UseDtLimit .or. UseLocalTimeStep)then
       ! Do not exceed DtLimitOrig
       DtLimit = min(DtLimitOrig, DtLimit*Factor)
       Cfl     = min(CflOrig, Cfl*Factor)
    else
       ! Do not exceed CflOrig
       Cfl     = min(CflOrig, Cfl*Factor)
    end if

    if(DoTestMe)then
       write(*,*) NameSub,': RelativeChangeMin,Max,Factor=', &
            RelativeChangeMin, RelativeChangeMax, Factor
       if(UseDtFixed)then
          write(*,*) NameSub,': Dt, DtFixed, Cfl=',&
               Dt*No2Si_V(UnitT_), DtFixed*No2Si_V(UnitT_), Cfl
       elseif(UseDtLimit .or. UseLocalTimeStep)then
          write(*,*) NameSub,': Dt, DtLimit, Cfl=',&
               Dt*No2Si_V(UnitT_), DtLimit*No2Si_V(UnitT_), Cfl
       else
          write(*,*) NameSub,': Dt, Cfl=', Dt*No2Si_V(UnitT_), Cfl
       end if
    end if

  end subroutine control_time_step

end module ModTimeStepControl
