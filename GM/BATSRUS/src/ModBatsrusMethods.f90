!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

! This file contains the top level methods for BATSRUS

subroutine BATS_setup

  use ModMpi
  use ModProcMH
  use ModMain
  use ModCT, ONLY : DoInitConstrainB
  use ModIO
  use ModAMR,      ONLY: nRefineLevelIC
  use ModAdvance,  ONLY: iTypeAdvance_B, iTypeAdvance_BP, ExplBlock_
  use ModParallel, ONLY: init_mod_parallel

  implicit none

  ! Local variables

  character(len=*), parameter :: NameSub = 'BATS_setup'
  integer :: iError 

  !---------------------------------------------------------------------------

  ! Allocate and initialize variables dependent on number of PEs
  call init_mod_parallel

  if(.not.IsStandAlone)call write_progress(0)

  call grid_setup   ! Restart reads info integer only

  call set_initial_conditions ! Restart reads all real data

  call find_test_cell

  call initialize_files

  call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------

  call write_runtime_values

  IsRestartCoupler = restart
  restart = .false.

contains
  !===========================================================================

  subroutine grid_setup

    ! Set up problem geometry, blocks, and grid structure.

    use ModIO, ONLY: restart
    use ModRestartFile, ONLY: NameRestartInDir, &
         UseRestartInSeries, string_append_iter

    use ModMain, ONLY: iteration_number
    use BATL_lib, ONLY: init_grid_batl, read_tree_file,set_amr_criteria,&
         set_amr_geometry, nBlock, Unused_B, init_amr_criteria, &
         nInitialAmrLevel
    use ModBatlInterface, ONLY: set_batsrus_grid
    ! Dummy variables, to avoid array size issues with State_VGB in
    ! set_amr_criteria
    use ModAdvance, ONLY : nVar, State_VGB
    use ModLoadBalance, ONLY: load_balance

    !LOCAL VARIABLES:
    character(len=*), parameter :: NameSubSub = NameSub//'::grid_setup'
    character(len=100) :: NameFile

    integer:: iBlock
    !--------------------------------------------------------------------------

    call set_batsrus_grid

    if (.not.restart) then
       ! Create initial solution block geometry.
       ! Set all arrays for AMR
       call init_amr_criteria

       ! Perform initial refinement of mesh and solution blocks.
       do nRefineLevel = 1, nInitialAmrLevel

          if (iProc == 0 .and. lVerbose > 0) then
             call write_prefix; write (iUnitOut,*) NameSub, &
                  ' starting initial refinement level, nBlockAll =', &
                  nRefineLevel, nBlockAll
          end if
          call set_amr_criteria(nVar, State_VGB,TypeAmrIn='geo')
          call init_grid_batl
          call set_batsrus_grid
          ! need to update node information, maybe not all
          ! of load balancing
          do iBlock = 1, nBlock
             if(Unused_B(iBlock)) CYCLE
             call set_amr_geometry(iBlock)
          end do
       end do
    else
       ! Read initial solution block geometry from octree restart file.

       NameFile = trim(NameRestartInDir)//'octree.rst'
       if (UseRestartInSeries) &
            call string_append_iter(NameFile,iteration_number)
       call read_tree_file(NameFile)
       call init_grid_batl
       call set_batsrus_grid

       ! Set all arrays for AMR
       call init_amr_criteria

    end if

    ! Set initial block types
    where(.not.Unused_B) iTypeAdvance_B = ExplBlock_
    call MPI_ALLGATHER(iTypeAdvance_B, MaxBlock, MPI_INTEGER, &
         iTypeAdvance_BP, MaxBlock, MPI_INTEGER, iComm, iError)

    ! Move coordinates around except for restart because the
    ! coordinate info is in the .rst files and not in the octree (1st).
    ! Do not move the data, it is not yet set. There are new blocks.
    call load_balance(DoMoveCoord=.not.restart, DoMoveData=.false., &
         IsNewBlock=.true.)
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       call set_amr_geometry(iBlock)
    end do

    if (iProc == 0.and.lVerbose>0)then
       call write_prefix; write (iUnitOut,*) '    total blocks = ',nBlockALL
    end if

    nRefineLevel = nInitialAmrLevel

    if(DoSetLevels) call set_levels
  end subroutine grid_setup

  !===========================================================================

  subroutine set_initial_conditions

    use ModIO,          ONLY: restart
    use ModIO,          ONLY: restart_Bface
    use ModRestartFile, ONLY: read_restart_files
    use ModMessagePass, ONLY: exchange_messages
    use ModMain,        ONLY: UseB0
    use ModB0,          ONLY: set_b0_reschange
    use ModFieldLineThread, ONLY: UseFieldLineThreads, set_threads
    use ModAMR,         ONLY: do_amr
    use ModLoadBalance, ONLY: load_balance

    use ModUserInterface ! user_initial_perturbation, user_action

    ! Set intial conditions for solution in each block.

    !LOCAL VARIABLES:
    character(len=*), parameter :: NameSubSub = &
         NameSub//'::set_initial_conditions'
    integer :: iLevel, iBlock
    !-------------------------------------------------------------------------
    if(.not.restart .and. nRefineLevelIC>0)then
       call timing_start('amr_ics')
       do iLevel=1, nRefineLevelIC
          call timing_start('amr_ics_set')
          do iBlock = 1, nBlockMax
             call set_ICs(iBlock)
          end do
          call timing_stop('amr_ics_set')

          ! Allow the user to add a perturbation and use that 
          ! for physical refinement.
          if (UseUserPerturbation)then
             call timing_start('amr_ics_perturb')
             call user_initial_perturbation
             call timing_stop('amr_ics_perturb')
          end if

          call timing_start('amr_ics_amr')
          ! Do physics based AMR with the message passing
          call do_amr(.true.,'phy')

          call timing_stop('amr_ics_amr')
       end do

       ! Move coordinates, but not data (?), there are new blocks
       call timing_start('amr_ics_balance')
       call load_balance(DoMoveCoord=.true.,DoMoveData=.false., &
            IsNewBlock=.true.)

       call timing_stop('amr_ics_balance')

       call timing_stop('amr_ics')
    end if
    ! nRefineLevelIC has done its work, reset to zero
    nRefineLevelIC = 0

    ! Read initial data from restart files as necessary.
    if(restart)then
       call user_action('reading restart files')
       call read_restart_files
    end if

    do iBlock = 1, nBlockMax
       ! Initialize solution blocks
       call set_ICs(iBlock)
    end do

    call user_action('initial condition done')

    ! Allow the user to add a perturbation to the initial condition.
    if (UseUserPerturbation) then
       call user_initial_perturbation
       UseUserPerturbation=.false.
    end if

    if (restart) then
       if(iProc==0)then
          call write_prefix; write(iUnitOut,*)&
               NameSub,' restarts at n_step,Time_Simulation=',&
               n_step,Time_Simulation
       end if
       ! Load balance for the inner blocks:
       call load_balance(DoMoveCoord=.true., DoMoveData=.true., &
            IsNewBlock=.true.)


       ! Redo the AMR level constraints for fixed body level
       ! The coordinates of the blocks are only known now
       if(DoSetLevels) call set_levels
    end if

    ! Fix face centered B0 at resolution changes
    if(UseB0)call set_b0_reschange
    if(UseFieldLineThreads)call set_threads

    ! Ensure zero divergence for the CT scheme
    if(UseConstrainB)then
       if(restart_Bface)then
          DoInitConstrainB=.false.
       else
          call BATS_init_constrain_b
       end if
    end if

    call exchange_messages

  end subroutine set_initial_conditions

  !============================================================================

  subroutine initialize_files
    use ModSatelliteFile, ONLY: set_satellite_file_status, nSatellite

    ! Local variables
    character(len=*), parameter :: NameSubSub = NameSub//'::initialize_files'
    integer :: iSat

    if (iProc == 0) then
       do iSat = 1, nSatellite
          call set_satellite_file_status(iSat,'open')
          call set_satellite_file_status(iSat,'close')
       end do
    end if

    plot_type(restart_)='restart'
    plot_type(logfile_)='logfile'

  end subroutine initialize_files

end subroutine BATS_setup

!==============================================================================

subroutine BATS_init_session

  use ModMain, ONLY: iSignRotationIC, UseUserPerturbation, &
       UseRadDiffusion, UseHeatConduction, UseIonHeatConduction, &
       UseProjection, UseConstrainB, UseParticles
  use ModCT,   ONLY: DoInitConstrainB
  use ModHallResist, ONLY: UseHallResist, init_hall_resist, UseBiermannBattery
  use ModImplicit, ONLY: UseFullImplicit, UseSemiImplicit, TypeSemiImplicit
  use ModRadDiffusion, ONLY: init_rad_diffusion
  use ModHeatConduction, ONLY: init_heat_conduction
  use ModParticleFieldLine, ONLY: init_particle_line
  use ModRestartFile, ONLY: UseRestartOutSeries
  use ModMessagePass, ONLY: exchange_messages
  use ModUserInterface ! user_initial_perturbation
  use ModLocalTimestep, ONLY: UseLocalTimeStepNew
  use ModProcMH, ONLY: iProc
  use BATL_lib, ONLY: init_amr_criteria
  use ModLoadBalance, ONLY: load_balance, select_stepping

  implicit none

  ! Local variables
  character(len=*), parameter :: NameSub = 'BATS_init_session '
  !---------------------------------------------------------------------------
  
  if(UseLocalTimeStepNew)then
     if(iProc == 0)write(*,*) &
          NameSub,': redo load balance for time accurate local timestepping'
     ! move coordinates and data, there are no new blocks  
     call load_balance(DoMoveCoord=.true., DoMoveData=.true., &
          IsNewBlock=.false.)
     UseLocalTimeStepNew = .false.
  end if

  ! Find the test cell defined by #TESTIJK or #TESTXYZ commands
  call find_test_cell

  !\
  ! Allow the user to add a perturbation to the initial condition.
  !/
  if (UseUserPerturbation) then
     call user_initial_perturbation
     UseUserPerturbation=.false.
  end if

  ! Set number of explicit and implicit blocks
  ! Partially implicit/local selection will be done in each time step
  call select_stepping(DoPartSelect=.false.)

  ! Transform velocities from a rotating system to the HGI system if required
  if(iSignRotationIC /= 0)then
     ! add or subtract rho*(omega x r) to the momentum of all fluids for all blocks
     call add_rotational_velocity(iSignRotationIC, 0) 
     iSignRotationIC = 0
  end if

  ! Ensure zero divergence for the CT scheme
  if(UseConstrainB .and. DoInitConstrainB)&
       call BATS_init_constrain_b

  if(UseHallResist .or. UseBiermannBattery)call init_hall_resist

  if(UseHeatConduction .or. UseIonHeatConduction) &
       call init_heat_conduction

  if(UseParticles) &
       call init_particle_line

  if(UseSemiImplicit)then
     select case(TypeSemiImplicit)
     case('radiation', 'radcond', 'cond')
        call init_rad_diffusion
     end select
  elseif(UseFullImplicit.and.UseRadDiffusion)then
     call init_rad_diffusion
  end if

  ! Make sure that ghost cells are up to date
  call exchange_messages

  if(UseProjection)call project_B

  call BATS_save_files('INITIAL')

  ! save initial restart series
  if (UseRestartOutSeries) call BATS_save_files('RESTART')

  ! Set all arrays for AMR
  call init_amr_criteria

end subroutine BATS_init_session

!============================================================================

subroutine BATS_advance(TimeSimulationLimit)

  ! Advance solution with one time step

  use ModKind
  use ModProcMH
  use ModMain
  use ModIO, ONLY: iUnitOut, write_prefix, save_plots_amr
  use ModAmr, ONLY: DnAmr, DoAmr, automatic_refinement, do_amr, DtAmr
  use ModPhysics, ONLY : No2Si_V, UnitT_
  use ModAdvance, ONLY: UseNonConservative, nConservCrit, UseAnisoPressure, &
       UseElectronPressure
  use ModPartSteady, ONLY: UsePartSteady, IsSteadyState, &
       part_steady_select, part_steady_switch
  use ModImplicit, ONLY: UseImplicit, n_prev, UseSemiImplicit
  use ModSemiImplicit, ONLY: advance_semi_impl
  use ModIeCoupling, ONLY: apply_ie_velocity
  use ModTimeStepControl, ONLY: UseTimeStepControl, control_time_step, &
       UseMaxTimeStep
  use ModParticleFieldLine, ONLY: advect_particle_line
  use ModLaserHeating,    ONLY: add_laser_heating
  use ModVarIndexes, ONLY: Te0_
  use ModMessagePass, ONLY: exchange_messages
  use ModTimeStepControl, ONLY: set_global_timestep
  use ModB0, ONLY: DoUpdateB0, DtUpdateB0
  use ModResistivity, ONLY: UseResistivity, UseHeatExchange, calc_heat_exchange
  use ModMultiFluid, ONLY: UseMultiIon
  use ModLocalTimeStep, ONLY: advance_localstep, UseLocalTimeStep
  use ModPartImplicit, ONLY: advance_part_impl
  use ModHeatConduction, ONLY: calc_ei_heat_exchange
  use ModFieldLineThread, ONLY: UseFieldLineThreads, advance_threads, Enthalpy_
  use ModLoadBalance, ONLY: load_balance_blocks
  implicit none

  !INPUT ARGUMENTS:
  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded

  ! Local variables
  character(len=*), parameter :: NameSub = 'BATS_advance'

  logical :: DoTest, DoTestMe

  real :: AmrTime = 0
  !---------------------------------------------------------------------------
  ! Check if time limit is reached (to be safe)
  if(Time_Simulation >= TimeSimulationLimit) RETURN

  ! Check if steady state is achieved
  if(.not.time_accurate .and. UsePartSteady .and. IsSteadyState)then
     ! Create stop condition for stand alone mode
     nIter = iteration_number
     ! There is nothing to do, simply return
     RETURN
  end if

  call set_oktest(NameSub, DoTest, DoTestMe)

  ! We are advancing in time
  time_loop = .true.

  ! Some files should be saved at the beginning of the time step
  call BATS_save_files('BEGINSTEP')

  n_step = n_step + 1
  iteration_number = iteration_number+1

  if(time_accurate) call set_global_timestep(TimeSimulationLimit)

  ! Select block types and load balance blocks
  call load_balance_blocks

  ! Switch off steady blocks to reduce calculation
  if(UsePartSteady) call part_steady_switch(.true.)

  ! Calculate time step dt
  !if (time_accurate) call set_global_timestep(TimeSimulationLimit)

  call timing_start('advance')

  if(UseNonConservative .and. nConservCrit > 0)&
       call select_conservative

  if(UseImplicit.and.nBlockImplALL>0)then
     call advance_part_impl
  elseif(UseLocalTimeStep .and. n_step > 1 .and. time_accurate) then
     call advance_localstep(TimeSimulationLimit)
  else
     call advance_expl(.true., -1)
  endif

  ! Adjust Time_Simulation to match TimeSimulationLimit if it is very close
  if(  time_accurate .and. &
       Time_Simulation < TimeSimulationLimit .and. &
       TimeSimulationLimit - Time_Simulation <= 1e-6*Dt*No2Si_V(UnitT_))then

     if(iProc == 0 .and. lVerbose > 0)then
        call write_prefix; write(iUnitOut,*) NameSub, &
             ': adjusting Time_Simulation=', Time_Simulation,&
             ' to TimeSimulationLimit=', TimeSimulationLimit,&
             ' with Dt=', Dt
     end if

     Time_Simulation = TimeSimulationLimit
  end if

  if(UseIM)call apply_im_pressure
  
  if(.not.UseMultiIon)then
     if(UseHeatConduction .and. UseElectronPressure)then
        if(.not.UseSemiImplicit)call calc_ei_heat_exchange
     elseif(UseResistivity .and. UseHeatExchange .and. UseElectronPressure)then
        call calc_heat_exchange
     end if
  end if

  if(UseAnisoPressure)call fix_anisotropy

  if(UseIE)call apply_ie_velocity

  if(UseDivBDiffusion)call clean_divb

  if(UseLaserHeating) call add_laser_heating

  ! Calculate temperature at the end of time step
  if(Te0_>1)call update_te0

  if(UseFieldLineThreads)call advance_threads(Enthalpy_)
  call exchange_messages

  if(UseSemiImplicit .and. (Dt>0 .or. .not.time_accurate)) &
       call advance_semi_impl

  if(UseTimeStepControl .and. time_accurate .and. Dt>0) call control_time_step

  if(UsePartSteady) then
     ! Select steady and unsteady blocks
     if(.not. (Time_Accurate .and. Time_Simulation == 0.0))then
        call timing_start('part_steady')
        call part_steady_select
        call timing_stop('part_steady')
     end if

     ! Switch on steady blocks to be included in AMR, plotting, etc.
     call part_steady_switch(.false.)
  end if

  call advect_all_points

  if(UseParticles) &
       call advect_particle_line


  call timing_stop('advance')

!  if(time_accurate)&
!       call update_lagrangian_grid(&
!       Time_Simulation - Dt*No2Si_V(UnitT_),Time_Simulation)

  if(DoTest)write(*,*)NameSub,' iProc,new n_step,Time_Simulation=',&
       iProc,n_step,Time_Simulation

  if (DoUpdateB0) then
     ! dB0/dt term is added at the DtUpdateB0 frequency

     if ( int(Time_Simulation/DtUpdateB0) >  &
          int((Time_Simulation - Dt*No2Si_V(UnitT_))/DtUpdateB0)) &
          call update_b0
  end if

  if (UseProjection) call project_B


  ! AmrTime is the time to do AMR.
  if(DoAmr .and. AmrTime < DtAmr) AmrTime = DtAmr
  if ( DoAmr .and. ((DnAmr > 0 .and. mod(n_step, DnAmr) == 0) .or. &
       ( DtAmr >0 .and. time_simulation  >AmrTime)))then
     if(DtAmr >0) AmrTime = DtAmr + AmrTime
     call timing_start(NameThisComp//'_amr')
     if(iProc==0 .and. lVerbose > 0 .and. (DnAmr > 1 .or. DtAmr >0))then
        call write_prefix; write(iUnitOut,*) &
             '>>>>>>>>>>>>>>>>>>>> AMR <<<<<<<<<<<<<<<<<<<<'
        if (time_accurate) call write_timeaccurate
     end if

     if (.not. automatic_refinement) nRefineLevel = nRefineLevel + 1

     ! BDF2 scheme should not use previous step after AMR
     n_prev = -100

     ! Do AMR without full initial message passing
     call do_amr(.false.,'all')

     ! Output timing after AMR.
     call timing_stop(NameThisComp//'_amr')
     if(iProc == 0 .and. lVerbose > 0 .and. (DnAmr > 1 .or. DtAmr >0))then
        call timing_show(NameThisComp//'_amr',1)
        call timing_show('load_balance',1)
        call write_prefix; write(iUnitOut,*) &
             '>>>>>>>>>>>>>>>>>>>> AMR <<<<<<<<<<<<<<<<<<<<'
     end if

     if (UseProjection) call project_B

     ! Write plotfiles after AMR if required
     if(save_plots_amr)call BATS_save_files('AMRPLOTS')

  else
     ! If AMR is done, then the plotting of BATS_save_files('NORMAL')
     ! is called in ModAMR to save the AMR criteria.  
     call BATS_save_files('NORMAL')
  end if

end subroutine BATS_advance

!============================================================================

subroutine BATS_init_constrain_b
  use ModProcMH
  use ModMain, ONLY: lVerbose, x_, y_, z_, nBlock
  use ModCT, ONLY : DoInitConstrainB
  use ModNumConst, ONLY: cTiny
  use ModAdvance, ONLY : Bx_, Bz_, State_VGB, tmp1_BLK
  use ModIO, ONLY: write_prefix, iUnitOut
  use BATL_lib, ONLY: Xyz_DGB, message_pass_cell
  implicit none

  ! Local variables
  character(len=*), parameter :: NameSub ='BATS_init_constrain_b '
  real, external :: maxval_loc_abs_BLK
  integer :: iBlock
  integer :: iLoc_I(5)  ! full location index
  real    :: divbmax_now
  !--------------------------------------------------------------------------
  DoInitConstrainB=.false.

  call message_pass_cell(3,State_VGB(Bx_:Bz_,:,:,:,:), nWidthIn=1, &
       nProlongOrderIn=1, DoSendCornerIn=.false., DoRestrictFaceIn=.true.)

  do iBlock=1, nBlock
     ! Estimate Bface from the centered B values
     call Bcenter2Bface(iBlock)
     ! Calculate energy (it is not set in set_ICs)
     ! because the projection scheme will need it
!!! call calc_energy(iBlock)
  end do

  call proj_get_divb(tmp1_BLK)
  divbmax_now=maxval_loc_abs_BLK(nProc,tmp1_BLK,iLoc_I)
  if(iProc == 0.and.lVerbose>0)then
     call write_prefix; write(iUnitOut,*)
     call write_prefix; write(iUnitOut,*) NameSub, &
          ' maximum of |div B| before projection=',divbmax_now
     call write_prefix; write(iUnitOut,*)
  end if
  if(divbmax_now>cTiny)then
     if(iProc==iLoc_I(5))then
        call write_prefix; write(iUnitOut,*) NameSub, &
             ' divB,loc,x,y,z=',divbmax_now,iLoc_I,&
             Xyz_DGB(:,iLoc_I(x_),iLoc_I(y_),iLoc_I(z_),iLoc_I(4))
     end if

     if(iProc == 0.and.lVerbose>0)then
        call write_prefix; write(iUnitOut,*)
        call write_prefix; write(iUnitOut,*) &
             NameSub,' projecting B for CT scheme...'
     end if

     ! Do the projection with UseConstrainB true
     call project_B

     ! Check and report the accuracy of the projection
     call proj_get_divb(tmp1_BLK)
     divbmax_now=maxval_loc_abs_BLK(nProc,tmp1_BLK,iLoc_I)
     if(iProc == 0.and.lVerbose>0)then
        call write_prefix; write(iUnitOut,*)
        call write_prefix; write(iUnitOut,*) NameSub, &
             ' maximum of |div B| after projection=',divbmax_now
        call write_prefix; write(iUnitOut,*)
     end if
     if(iProc==iLoc_I(5).and.divbmax_now>cTiny)then
        call write_prefix; write(iUnitOut,*) NameSub, &
             ' divB,loc,x,y,z=',divbmax_now,iLoc_I,&
             Xyz_DGB(:,iLoc_I(x_),iLoc_I(y_),iLoc_I(z_),iLoc_I(4))
     end if
  end if

end subroutine BATS_init_constrain_b

!===========================================================================

subroutine BATS_save_files(TypeSaveIn)

  use ModProcMH
  use ModMain
  use ModIO
  use ModUtilities, ONLY : upper_case
  use ModMessagePass, ONLY: exchange_messages

  implicit none

  character(len=*), intent(in) :: TypeSaveIn

  character(len=len(TypeSaveIn)) :: TypeSave
  logical :: DoExchangeAgain, DoAssignNodeNumbers, IsFound, DoSaveRestartTmp
  integer :: iFile

  character(len=*), parameter :: NameSub='BATS_save_files'
  !--------------------------------------------------------------------------

  DoExchangeAgain     = .false.
  DoAssignNodeNumbers = .true.
  TypeSave = TypeSaveIn
  call upper_case(TypeSave)

  select case(TypeSave)
  case('INITIAL')
     ! Do not save current step or time
     n_output_last = n_step

     ! Initialize last save times
     where(dt_output>0.) &
          t_output_last=int(time_simulation/dt_output)

     ! DoSaveInitial may be set to true in the #SAVEINITIAL command
     if(DoSaveInitial .or. (time_accurate .and. time_simulation == 0.0))then
        if(DoSaveInitial)then
           ! Save all (except restart files)
           n_output_last = -1
           t_output_last = -1.0
        else
           ! Save only those with a positive time frequency
           where(dt_output>0.)
              n_output_last = -1
              t_output_last = -1.0
           end where
        end if
        ! Do not save restart file in any case
        n_output_last(restart_) = n_step
        call save_files
     end if
     ! Set back to default value (for next session)
     DoSaveInitial = .false.
  case('FINAL')
     save_restart_file = .false.
     call save_files_final
  case('FINALWITHRESTART')
     call save_files_final
  case('NORMAL', 'BEGINSTEP', 'BEFOREAMR')
     call save_files
  case('AMRPLOTS')
     do iFile=plot_+1, plot_+nPlotFile
        call save_file
     end do
     if(DoExchangeAgain) &
          call exchange_messages(DoResChangeOnlyIn=.true.)
  case('RESTART')
     DoSaveRestartTmp = save_restart_file
     save_restart_file = .true.
     iFile = restart_
     call save_file
     save_restart_file = DoSaveRestartTmp
  case default
     call stop_mpi(NameSub//' ERROR incorrect TypeSave='//TypeSave)
  end select

contains
  !==========================================================================
  subroutine save_files

    logical :: DoSave = .false.
    integer :: t_output_current

    do iFile = 1, nFile
       ! We want to use the IE magnetic perturbations that were passed 
       ! in the last coupling together with the current GM perturbations.
       if( (iFile==magfile_ .or. iFile==maggridfile_) &
            .neqv. TypeSave == 'BEGINSTEP') CYCLE

       if(dn_output(ifile)>=0)then
          if(dn_output(ifile)==0)then
             call save_file
          else if(mod(n_step,dn_output(ifile))==0)then
             call save_file
          end if
       else if(time_accurate .and. dt_output(ifile)>0.)then
          t_output_current = int(time_simulation/dt_output(ifile))
          DoSave = .false.
          if(t_output_current>t_output_last(ifile)) DoSave = .true.
          if(DoSave)then
             t_output_last(ifile)=t_output_current
             call save_file
          end if
       end if
    end do
    ! If message passing with corners was done in save_file for tecplot plots,
    ! then do exchange_messages over again to get expected values 
    ! in ghost cells.

    if(DoExchangeAgain)then
       if(iProc==0.and.lVerbose>0)then
          call write_prefix; write(iUnitOut,*)&
               'Calling exchange_messages to reset ghost cells ...'
       end if
       call exchange_messages(DoResChangeOnlyIn=.true.)
    end if

  end subroutine save_files

  !==========================================================================

  subroutine save_file

    use ModRestartFile, ONLY: write_restart_files
    use ModSatelliteFile, ONLY: &
         nSatellite, set_satellite_file_status, set_satellite_flags, &
         TimeSatStart_I, TimeSatEnd_I, iCurrent_satellite_position
    use ModGroundMagPerturb, ONLY: DoSaveMags, &
         DoSaveGridmag, write_magnetometers, DoWriteIndices, &
         write_geoindices
    use ModMessagePass, ONLY: exchange_messages

    integer :: iSat

    ! Backup location for the Time_Simulation variable.
    ! Time_Simulation is used in steady-state runs as a loop parameter
    ! in the save_files subroutine, where set_satellite_flags and 
    ! write_logfile are called with different Time_Simulation values
    ! spanning all the satellite trajectory cut. Old Time_Simulation value
    ! is saved here before and it is restored after the loop.
    !
    real :: tSimulationBackup = 0.0
    !---------------------------------------------------------------------

    if(n_step<=n_output_last(ifile) .and. dn_output(ifile)/=0 &
         .and. (n_step/=0 .or. ifile/=restart_) ) return

    if(ifile==restart_) then
       ! Case for restart file
       if(.not.save_restart_file)return
       call write_restart_files

    elseif(ifile==logfile_) then
       ! Case for logfile 

       if(.not.save_logfile)return
       call timing_start('save_logfile')
       call write_logfile(0,ifile)
       call timing_stop('save_logfile')

    elseif(ifile>plot_ .and. ifile<=plot_+nplotfile) then
       ! Case for plot files
       IsFound=.false.

       if(TypeSave /= 'BEFOREAMR' .and. .not.DoExchangeAgain .and. ( &
            index(plot_type(iFile),'lin')==1 .or. &
            index(plot_type(iFile),'pnt')==1 .or. &
            index(plot_type(iFile),'eqr')==1 .or. &
            index(plot_type(iFile),'eqb')==1 .or. &
            index(plot_type(iFile),'ieb')==1 .or. &
            index(plot_type(iFile),'lcb')==1 .or. &
            index(plot_type(iFile),'los')==1 .or. &
            index(plot_type(iFile),'sph')==1 .or. &
            plot_form(iFile) == 'tec')) then

          if(iProc==0.and.lVerbose>0)then
             call write_prefix; write(iUnitOut,*)&
                  ' Message passing for plot files ...'
          end if
          call exchange_messages(UseOrder2In=.true.)
          DoExchangeAgain = .true.
       end if

       if(index(plot_type(iFile),'los')>0) then
          IsFound = .true.
          call write_plot_los(iFile)
       end if

       if(index(plot_type(iFile),'pnt')>0) then
          IsFound = .true.
          call write_plot_particle(iFile)
       end if

       if(index(plot_type(iFile),'rfr')>0) then
          IsFound = .true.
          call write_plot_radiowave(iFile)
       end if

       if(index(plot_type(iFile),'lin')>0) then
          IsFound = .true.
          call write_plot_line(iFile)
       end if

       if(  index(plot_type(iFile),'eqr')>0 .or. &
            index(plot_type(iFile),'eqb')>0) then
          IsFound = .true.
          call plot_ray_equator(iFile)
       end if

       if(index(plot_type(iFile),'ieb')>0) then
          IsFound = .true.
          call ieb_plot(iFile)
       end if

       if(index(plot_type(iFile),'lcb')>0) then
          IsFound = .true.
          call lcb_plot(iFile)
       end if

       if(index(plot_type(iFile),'buf')>0)then
          IsFound = .true.
          if(TypeSaveIn/='INITIAL')call plot_buffer(iFile)
       end if

       if(plot_type(ifile)/='nul' .and. .not.IsFound ) then
          ! Assign node numbers for tec plots
          if( index(plot_form(ifile),'tec')>0 .and. DoAssignNodeNumbers)then
             call assign_node_numbers
             DoAssignNodeNumbers = .false.
          end if


          if(  index(plot_type(ifile),'ray')>0 .or. &
               index(plot_vars(ifile),'status')>0) call ray_trace

          call timing_start('save_plot')
          call write_plot_common(ifile)
          call timing_stop('save_plot')
       end if
    elseif(iFile > Satellite_ .and. iFile <= Satellite_ + nSatellite) then

       ! Case for satellite files
       iSat = iFile - Satellite_
       call timing_start('save_satellite')
       if(iProc==0)call set_satellite_file_status(iSat,'append')
       !
       ! Distinguish between time_accurate and .not. time_accurate:
       !

       if (time_accurate) then
          call set_satellite_flags(iSat)
          ! write one line for a single trajectory point
          call write_logfile(iSat,ifile)
       else
          tSimulationBackup = Time_Simulation    ! Save ...
          Time_Simulation = TimeSatStart_I(iSat)
          do while (Time_Simulation <= TimeSatEnd_I(iSat))
             call set_satellite_flags(iSat)
             ! write for ALL the points of trajectory cut
             call write_logfile(iSat,ifile)  
             Time_Simulation = Time_Simulation + dt_output(iSat+Satellite_)
          end do
          Time_Simulation = tSimulationBackup    ! ... Restore
          icurrent_satellite_position(iSat) = 1
          if(iProc==0)call set_satellite_file_status(iSat,'close')
       end if
       call timing_stop('save_satellite')

    elseif(ifile == magfile_) then
       !Cases for magnetometer files
       if(.not.DoSaveMags) RETURN
       if(time_accurate) then  
          call timing_start('save_magnetometer')
          call write_magnetometers('stat')   
          call timing_stop('save_magnetometer')  
       end if

    elseif(ifile == maggridfile_) then
       !Case for grid magnetometer files
       if(.not. DoSaveGridmag) RETURN
       if(time_accurate) then  
          call timing_start('grid_magnetometer')
          call write_magnetometers('grid')   
          call timing_stop('grid_magnetometer')  
       end if

    elseif(ifile == indexfile_) then
       ! Write geomagnetic index file.
       if(time_accurate .and. DoWriteIndices) call write_geoindices
    end if

    n_output_last(ifile) = n_step

    if(iProc==0 .and. lVerbose>0 .and. &
         ifile /= logfile_ .and. iFile /= magfile_ .and. iFile /= indexfile_ &
         .and. iFile /= maggridfile_ &
         .and. (iFile <= satellite_ .or. iFile > satellite_ + nSatellite))then
       if(time_accurate)then
          call write_prefix; 
          write(iUnitOut,'(a,i2,a,a,a,i7,a,i4,a,i2.2,a,i2.2,a)') &
               'saved ifile=',ifile,' type=',plot_type(ifile),&
               ' at n_step=',n_step,' time=', &
               int(                            Time_Simulation/3600.),':', &
               int((Time_Simulation-(3600.*int(Time_Simulation/3600.)))/60.),':', &
               int( Time_Simulation-(  60.*int(Time_Simulation/  60.))), &
               ' h:m:s'
       else
          call write_prefix; write(iUnitOut,'(a,i2,a,a,a,i7)') &
               'saved ifile=',ifile,' type=',plot_type(ifile), &
               ' at n_step=',n_step
       end if
    end if

  end subroutine save_file

  !===========================================================================

  subroutine save_files_final
    use ModSatelliteFile, ONLY: set_satellite_file_status, nSatellite
    use ModGroundMagPerturb, ONLY: finalize_magnetometer

    integer :: iSat
    !-----------------------------------------------------------------------
    do iFile = 1, plot_ + nPlotFile
       call save_file
    end do

    ! Close files
    if (iProc==0) then
       do iSat = 1, nSatellite
          call set_satellite_file_status(iSat,'close')
       end do
    end if

    call finalize_magnetometer

    if (save_logfile.and.iProc==0.and.unit_log>0) close(unit_log)

  end subroutine save_files_final

end subroutine BATS_save_files

!=============================================================================
subroutine BATSRUS_finalize

  use ModAdvance,      ONLY: clean_mod_advance
  use ModGeometry,     ONLY: clean_mod_geometry
  use ModNodes,        ONLY: clean_mod_nodes
  use ModCT,           ONLY: clean_mod_ct
  use ModRaytrace,     ONLY: clean_mod_raytrace
  use ModPartImplicit, ONLY: clean_mod_part_impl
  use ModSemiImplicit, ONLY: clean_mod_semi_impl
  use ModIeCoupling,   ONLY: clean_mod_ie_coupling
  use BATL_lib,        ONLY: clean_batl

  implicit none

  integer:: iError
  !---------------------------------------------------------------------------
  call clean_batl
  call clean_mod_advance
  call clean_mod_ct
  call clean_mod_part_impl
  call clean_mod_semi_impl
  call clean_mod_geometry
  call clean_mod_nodes
  call clean_mod_raytrace
  call clean_mod_ie_coupling

  ! call clean_mod_boundary_cells !!! to be implemented
  ! call clean_mod_resistivity !!! to be implemented

  call error_report('PRINT',0.,iError,.true.)

end subroutine BATSRUS_finalize
