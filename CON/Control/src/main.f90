!^CFG COPYRIGHT UM
!============================================================
!           SWMF: Space Weather Modeling Framework          |
!                   University of Michigan                  |
!============================================================
program SWMF

  use CON_world, world_init => init, world_setup => setup, world_clean => clean
  use CON_comp_param
  use CON_wrapper
  use CON_coupler, ONLY: set_grid_pointer_all, DnCouple_II, DtCouple_II, &
       check_couple_symm
  use CON_io, ONLY : read_inputs, DnShowProgressShort, DnShowProgressLong, &
       DoSaveRestart, DnSaveRestart, DtSaveRestart, &
       CON_save_restart => save_restart

  use CON_time
  use CON_axes, ONLY: DtUpdateB0, set_axes
  use CON_variables, ONLY: DoCheckStopFile
  use ModMpi
  use ModIoUnit, ONLY: UNITTMP_
  
  use ModMain
  use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK
  use ModCT, ONLY : DoInitConstrainB               !^CFG IF CONSTRAINB
  use ModImplicit, ONLY : UsePartImplicit,n_prev   !^CFG IF IMPLICIT
  use ModIO
  use ModAMR, ONLY : dn_refine,initial_refine_levels,nRefineLevelIC,automatic_refinement
  use ModPhysics, ONLY : unitSI_t
  use ModNumConst
  implicit none

  character (len=*), parameter :: NameSub = 'main'

  !\
  ! Local variable definitions.
  !/
  integer :: iComm, iProc, nProc, lComp, iComp
  integer :: idepth, isession=1, ifile, nIterExpect, nIterExpectTime
  integer :: iBlock, nBlockMoved
  integer :: iError

  logical :: last_session, stop_now
  logical :: IsFound
  logical :: Ray_uninitialized=.true.
  logical :: RCM_uninitialized=.true.

  real    :: divbmax_now
  real, external :: maxval_loc_abs_BLK
  real*8 :: time_start

  real*8, external :: timing_func_d

  logical :: DoExchangeAgain

  logical :: DoIM        !^CFG  IF RCM

  logical :: oktest, oktest_me
  integer :: iLoc_I(5)  ! full location index

  !---------------------------------------------------------------------------

  !================================================================
  ! INITIALIZATIONS INITIALIZATIONS INITIALIZATIONS INITIALIZATIONS
  !================================================================

  !\
  ! Initialize control component (MPI)
  !/
  call world_init

  !\
  ! Delete SWMF.SUCCESS and SWMF.STOP files if found
  !/
  if(is_proc0())then

     inquire(file='SWMF.SUCCESS',EXIST=IsFound)
     if(IsFound)then
        open(UNITTMP_, file = 'SWMF.SUCCESS')
        close(UNITTMP_,STATUS = 'DELETE')
     end if

     inquire(file='SWMF.STOP',EXIST=IsFound)
     if (IsFound) then 
        open(UNITTMP_, file = 'SWMF.STOP')
        close(UNITTMP_, STATUS = 'DELETE')
     endif

  end if

  !\
  ! Read component information from LAYOUT.in
  !/
  call world_setup

  !\
  ! Initialize inside time loop indicator and other parameters.
  !/
  time_loop  = .false.
  time_start = MPI_WTIME()

  ! local MPI parameters
  iProc = i_proc()
  nProc = n_proc()
  iComm = i_comm()


  ! read and store version name and number of registered components
  ! initialize the MPI parameters for the registered components
  do lComp = 1, n_comp()
     iComp=i_comp(lComp)
     call set_param_comp(iComp,"VERSION")
     if(.not.use_comp(iComp))then
        if(is_proc0())then
           write(*,'(a)')'CON_main SWMF_ERROR registered component '// &
                NameComp_I(iComp)//' is OFF!'
           write(*,'(a)')'Compile in a working component or remove it from '//&
                NameMapFile
        end if
        ! stop in a clean fashion (without abort)
        call world_clean
     end if
     call set_param_comp(iComp,'MPI')     ! Initialize MPI parameters
     call set_param_comp(iComp,'STDOUT')  ! Set prefix string for STDOUT
  end do
  ! Show framework and component information
  if(is_proc0())call show_all_comp
  
  ! Set individual grid pointers
  call set_grid_pointer_all

  ! Allocate and initialize variables dependent on number of PEs
  call allocate_vars  

  !\
  ! read input parameters
  !/
  call read_inputs(.true.,last_session)

!!! Temporary solutions
  ! This may be needed somewhere
  RealWorldTime       = TimeStart % Time
  time_simulation     = tSimulation
  write(*,*)'time_simulation=',time_simulation
  write(*,*)'TimeStart=',TimeStart

  save_restart_file   = DoSaveRestart
  dn_output(restart_) = DnSaveRestart
  dt_output(restart_) = DtSaveRestart
  UseIonosphere = use_comp(IE_)
  UseIM         = use_comp(IM_)
  call check_couple_symm(IE_,GM_,NameSub)
  dn_couple_ionosphere = DnCouple_II(IE_,GM_)
  dt_couple_ionosphere = DtCouple_II(IE_,GM_)
  call check_couple_symm(IM_,GM_,NameSub)
  dnCoupleIM = DnCouple_II(IM_,GM_)
  dtCoupleIM = DtCouple_II(IM_,GM_)
!!! 

  !!! This belongs to GM only !!!
  call write_progress(0)

  call timing_start('SWMF')
  call timing_start('setup')

  call set_oktest('main',oktest,oktest_me)

  call grid_setup              ! Restart reads integer only (no X_BLK or dx_BLK)

  call set_initial_conditions  ! Restart reads all real data

  call find_test_cell

  oktest_me = oktest .and. iProc==ProcTest ! Correct oktest_me based on ProcTest

  call initialize_files

  call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------

  call timing_stop('setup')

  call timing_show('setup',1)

  call write_runtime_values

  !\
  ! Begin sessions
  !/
  SESSIONLOOP: do 

     if(oktest_me)write(*,*)'main: starting SESSIONLOOP with last_session=',&
          last_session

     if (use_comp(UA_) .and. is_proc(UA_)) &
          call UA_init_session(iSession, tSimulation)

     if(UseIonosphere .and. is_proc(IE_)) call IE_session_init(iSession)

     !^CFG IF USERFILES BEGIN
     !\
     ! Allow the user to load a restart file and then add a perturbation to it
     ! before restarting a run. 
     !/
     if (time_accurate .and. UseUserPerturbation .and. restart) &
          call user_initial_perturbation
     !^CFG END USERFILES

     ! Set number of explicit and implicit blocks !^CFG IF  IMPLICIT BEGIN
     ! Partially implicit/local selection will be done in each time step
     call select_stepping(.false.)                !^CFG END IMPLICIT 

     ! Ensure zero divergence for the CT scheme   !^CFG IF CONSTRAINB
     if(UseConstrainB .and. DoInitConstrainB)&    !^CFG IF CONSTRAINB
          call init_constrain_b                   !^CFG IF CONSTRAINB

     ! Make sure that ghost cells are up to date
     call exchange_messages

     if (use_comp(IE_)) &                        !^CFG  IF IONOSPHERE
          call magneto_iono_coupling             !^CFG  IF IONOSPHERE

     !^CFG  IF RCM BEGIN
     if (use_comp(IM_).and.RCM_uninitialized) then

        !^CFG  IF RAYTRACE BEGIN
        if(.not.UseRaytrace) UseRaytrace=.true.
        if (Ray_uninitialized) then
           if(oktest_me) &
                write(*,*)'main: executing first ray_trace for first rcm.'
           call timing_start('Ray Tracing')
           call ray_trace(.true.)
           call timing_stop('Ray Tracing')
           Ray_uninitialized = .false.
        end if
        !^CFG END RAYTRACE

        if(oktest_me)write(*,*)'main: executing first rcm at start of session.'
        call timing_start('IM_0')
        call IM(0)
        call timing_stop('IM_0')

        call timing_start('MH_volume')
        call couple_gm_im
        call timing_stop('MH_volume')

        call couple_ie_im

        call timing_start('IM_1')
        call IM(1)
        call timing_stop('IM_1')

        call timing_start('IM_pressure_reset')
        call couple_im_gm
        call timing_stop('IM_pressure_reset')

        call timing_start('MH_volume')
        call couple_gm_im
        call timing_stop('MH_volume')
        call couple_ie_im
        RCM_uninitialized = .false.
     end if
     !^CFG END RCM

     !^CFG  IF RAYTRACE BEGIN
     if (UseRaytrace.and.Ray_uninitialized) then
        if(oktest_me)&
             write(*,*)'main: executing first ray_trace at start of session.'
        call timing_start('Ray Tracing')
        call ray_trace(.true.)
        call timing_stop('Ray Tracing')
        Ray_uninitialized = .false.
     end if
     !^CFG END RAYTRACE

     !=================================================================
     ! COMPUTATIONS COMPUTATIONS COMPUTATIONS COMPUTATIONS COMPUTATIONS
     !=================================================================

     if(nIteration==0)then
        if(dn_timing > -3)call timing_report_total
        if(is_proc0())write(*,*)'Resetting timing counters after setup.'
        call timing_reset('#all',3)
     end if

     time_loop = .true.

     !\
     ! Begin time step (iterations) loop.
     !/
     TIMELOOP: do

        if (use_comp(UA_) .and. is_proc(UA_)) call UA_run(tSimulation, 1.0)

        if(oktest_me)write(*,*)'main: Starting TIMELOOP at n_step=',n_step

        !^CFG IF IMPLICIT BEGIN

        ! Select blocks for partially local time stepping
        if(UsePartLocal)call select_stepping(.true.)

        ! Select and load balance blocks for partially implicit scheme
        if(UsePartImplicit)then                      
           ! Redo load balancing for partially implicit scheme
           call load_balance(.true.,.true.,nBlockMoved)
           if(nBlockMoved>0)then
              if(is_proc0().and.lVerbose>0)write(*,*)&
                   'load_balance finished: nBlockMoved=',nBlockMoved
              call find_neighbors
              call analyze_neighbors   !^CFG IF DEBUGGING
              call find_test_cell
           end if
        end if

        !^CFG END IMPLICIT

        if (UseProjection) call project_B    !^CFG IF PROJECTION

        call save_files

        if(MaxIteration >= 0 .and. nIteration >= MaxIteration) exit TIMELOOP
        if(DoTimeAccurate .and. tSimulationMax > cZero &
             .and. tSimulation >= tSimulationMax) &
             exit TIMELOOP

        stop_now=.false.
        if(is_proc0())then
           if(CpuTimeMax > cZero .and. MPI_WTIME()-time_start >= CpuTimeMax)&
                then
              write(*,*)'CPU time exceeded:',CpuTimeMax,MPI_WTIME()-time_start
              stop_now=.true.
           end if
           if(.not.stop_now .and. DoCheckStopfile) then
              inquire(file='SWMF.STOP',exist=stop_now)
              if (stop_now) &
                   write(*,*)'SWMF.STOP file exists: recieved stop signal'
           end if
        end if

        call MPI_BCAST(stop_now,1,MPI_LOGICAL,0,iComm,iError)
        if(stop_now)exit SESSIONLOOP

        nStep = nStep + 1
        nIteration = nIteration+1
        n_step = n_step + 1                     ! GM_run
        iteration_number = iteration_number + 1 ! GM_run
        call timing_step(nStep)

        ! Advance solution

        if (time_accurate) call set_global_timestep

        call timing_start('advance')
        if(UseImplicit.and.nBlockImplALL>0)then !^CFG IF IMPLICIT BEGIN
           call advance_impl
        else                                    !^CFG END IMPLICIT
           call advance_expl(.true., .true.)
        endif                                   !^CFG IF IMPLICIT
        call timing_stop('advance')

        if (time_accurate)then
           tSimulation = tSimulation + dt*unitSI_t
           TimeCurrent % Time =  TimeStart % Time + tSimulation
           call time_real_to_int(TimeCurrent)
           Time_Simulation = Time_Simulation + dt*unitSI_t  ! GM_run
           RealWorldTime = TimeCurrent % Time               ! GM_run
        end if

        if(mod(nStep,DnShowProgressShort)==0 .and. is_proc0())then
           write(*,*)'Short progress report at nStep,tSimulation=',&
                nStep,tSimulation
        end if

        if(mod(nStep,DnShowProgressLong)==0 .and. is_proc0())then
           write(*,*)'Long progress report at nStep=',&
                nStep,tSimulation
        end if


        !^CFG  IF RAYTRACE BEGIN
        !\
        ! Periodically perform raytracing calculations as required.
        !/
        if (UseRaytrace) then
           if (mod(n_step,dn_raytrace)==0)then
              call timing_start('Ray Tracing')
              call ray_trace(check_rayloop)
              call timing_stop('Ray Tracing')
           end if
        end if
        !^CFG END RAYTRACE

        !^CFG  IF RCM BEGIN
        !\
        ! Periodically perform IM coupling as required.
        !/
        if (UseIM) then
           DoIM=.false.
           if(time_accurate)then
              if(  int((Time_Simulation+0.5*dtCoupleIM)/dtCoupleIM) > &
                   int((Time_Simulation+0.5*dtCoupleIM-dt*unitSI_t)/dtCoupleIM))then
                 DoIM=.true.
                 if(is_proc0())then
                    write(*,*)' '
                    write(*,*)' Starting ray_trace for coupling at T=',Time_Simulation
                 end if
              end if
           else
              if(dnCoupleIM > 0)then
                 if(mod(nStep,dnCoupleIM)==0)then
                    DoIM=.true.
                 end if
              end if
           end if
           if (DoIM) then
              call timing_start('Ray Tracing')
              call ray_trace(check_rayloop)
              call timing_stop('Ray Tracing')

              call timing_start('MH_volume')
              call couple_gm_im
              call timing_stop('MH_volume')

              call couple_ie_im

              call timing_start('IM_2')
              call IM(2)
              call timing_stop('IM_2')

              call timing_start('IM_pressure_reset')
              call couple_im_gm
              call timing_stop('IM_pressure_reset')
           end if
        end if
        !^CFG END RCM

        !^CFG IF IONOSPHERE BEGIN
        !\
        ! Periodically perform ionosphere/magnetosphere coupling 
        ! calculations as required.
        !/
        if (UseIonosphere) then
           if (dn_couple_ionosphere > 0) then
              if (mod(nStep,dn_couple_ionosphere)==0) &
                   call magneto_iono_coupling
           else
              if (time_accurate) then
                 if(abs(dt_couple_ionosphere)<1.E-8) &
                      call stop_mpi('dt_couple_ionosphere<=1.E-8, fix it.')
                 if (int(Time_Simulation/dt_couple_ionosphere) > &
                      int((Time_Simulation-dt*unitSI_t)/dt_couple_ionosphere))&
                      call magneto_iono_coupling
              endif
           end if
        end if
        !^CFG END IONOSPHERE

        if (time_accurate .and. DtUpdateB0 > cZero) then
           if (int(Time_Simulation/DtUpdateB0) >  &
                int((Time_Simulation - dt*unitSI_t)/DtUpdateB0)) then
              if (is_proc0().and.lVerbose>0) &
                   write(*,*) "=> Updating B0 at Simulation time",&
                   Time_Simulation
              call timing_start('update_B0')
              do iBlock=1,nBlock
                 if(unusedBLK(iBlock)) CYCLE
                 call set_b0(iBlock)
              end do
              call timing_stop('update_B0')
           endif
        endif

        if ( dn_refine > 0 .and. mod(n_step,dn_refine)==0 )  &
             call amr_refinement

     end do TIMELOOP

     !\
     ! Check if there is anything else to do
     !/

     if(last_session)then
        EXIT SESSIONLOOP
     else
        if(is_proc0().and.lVerbose>=0) &
             write(*,*)'----- End of Session   ',isession,' ------'
        isession=isession+1
        if (dn_timing > -2) call timing_report
        call read_inputs(.false., last_session)
!!! Temporary solutions
        save_restart_file   = DoSaveRestart
        dn_output(restart_) = DnSaveRestart
        dt_output(restart_) = DtSaveRestart
        UseIonosphere = use_comp(IE_)
        UseIM         = use_comp(IM_)
        call check_couple_symm(IE_,GM_,NameSub)
        dn_couple_ionosphere = DnCouple_II(IE_,GM_)
        dt_couple_ionosphere = DtCouple_II(IE_,GM_)
        call check_couple_symm(IM_,GM_,NameSub)
        dnCoupleIM = DnCouple_II(IM_,GM_)
        dtCoupleIM = DtCouple_II(IM_,GM_)
!!!
        call timing_reset_all
        if(is_proc0().and.lVerbose>=0)&
             write(*,*)'----- Starting Session ',isession,' ------'
        call find_test_cell
        call set_oktest('main',oktest,oktest_me)
     end if

  end do SESSIONLOOP

  if(is_proc0())then
     if(lVerbose>=0)then
        write(*,*)
        write(*,'(a)')'    Finished Numerical Simulation'
        write(*,'(a)')'    -----------------------------'
        if (time_accurate) call write_timeaccurate
     end if
     if (dn_timing > -2) call timing_report
  end if


  !============================================================
  ! DATA OUTPUT DATA OUTPUT DATA OUTPUT DATA OUTPUT DATA OUTPUT
  !============================================================

  time_loop = .false.

  call save_files_final
  if(is_proc('IE'))call IE_finalize

  if(is_proc0().and.lVerbose>0)then
     write(*,*)
     write(*,'(a)')'    Finished Saving Output Files'
     write(*,'(a)')'    ----------------------------'
  end if

  !============================================================
  ! END END END END END END END END END END END END END END END
  !============================================================

  call timing_stop('SWMF')

  if(dn_timing > -3)call timing_report_total

  call error_report('PRINT',0.,iError,.true.)

  if(is_proc0())then
     open(UNITTMP_, file = 'SWMF.SUCCESS')
     close(UNITTMP_)
  end if

  call world_clean

contains

  !===========================================================================

  subroutine grid_setup

    !\
    ! Set up problem geometry, blocks, and grid structure.
    !/
    logical :: local_refine(nBLK)

    !--------------------------------------------------------------------------

    call set_root_block_geometry
    call build_octree_roots   ! Initialize octree data structure.
    call find_neighbors       ! Get initial neighbor information.

    if (.not.restart) then
       ! Create initial solution block geometry.

       ! Perform initial refinement of mesh and solution blocks.
       do idepth = 1, initial_refine_levels
          if (is_proc0().and.lVerbose>0)&
               write (*,*) 'Starting initial refinement level ',idepth
          call specify_initial_refinement(local_refine, idepth)
          call refine_grid(local_refine)
          call fixRefinementLevels
       end do
    else
       ! Read initial solution block geometry from octree restart file.

       ! Read restart header file only if old type.
       if(.not.restart_reals)call read_restart_header  
       call read_octree_file     ! Read octree restart file.

    end if
    ! number blocks and balance load
    call number_soln_blocks

    ! Move coordinates around except for restart from new restart files 
    ! which have coordinate info in the .rst files and not in the octree.
    call load_balance(.not.(restart .and. restart_reals),.false.,nBlockMoved)

    call find_neighbors

    if (is_proc0().and.lVerbose>0)&
         write (*,*) '    total blocks = ',nBlockALL

    idepth = initial_refine_levels

    if(DoSetLevels) call set_levels

    call analyze_neighbors    !^CFG IF DEBUGGING 

  end subroutine grid_setup

  !===========================================================================

  subroutine set_initial_conditions

    !\
    ! Set intial conditions for solution in each block.
    !/

    integer :: iLevel
    logical :: restart_read

    !-------------------------------------------------------------------------
    restart_read = .false.

    if(.not.restart .and. nRefineLevelIC>0)then
       do iLevel=1,nRefineLevelIC
          do globalBLK = 1, nBlockMax
             call set_ICs
          end do
          !^CFG IF USERFILES BEGIN
          !\
          ! Allow the user to load a restart file and then add a perturbation to it
          ! before restarting a run. 
          !/
          if (UseUserPerturbation) &
               call user_initial_perturbation
          !^CFG END USERFILES
          call amr_physics
          call fixRefinementLevels
          call number_soln_blocks
       end do
       call load_balance(.true.,.false.,nBlockMoved)
       call find_neighbors
       call analyze_neighbors      !^CFG IF DEBUGGING
       call find_test_cell
    end if

    do globalBLK = 1, nBlockMax
       !\
       ! Read initial data for solution block
       ! from restart file as necessary.
       !/
       if (restart .and. .not.unusedBLK(globalBLK)) then
          call timing_start('read_restart')
          call read_restart_file
          call timing_stop('read_restart')

          if(restart_reals)call fix_block_geometry(globalBLK)

          if (.not.restart_read) then
             restart_read = .true.
             call write_progress(1)

             ! Now we have time_simulation read from the first restart file
             if (Time_Simulation > cZero) then
                tSimulation = Time_Simulation
                TimeCurrent % Time =  TimeStart % Time + Time_Simulation
                RealWorldTime = TimeCurrent % Time
                call time_real_to_int(TimeCurrent)
                call time_int_to_string(TimeCurrent)
                if (time_accurate .and. DtUpdateB0 > cZero) then
                   if (is_proc0().and.lVerbose>0) &
                        write(*,*) "=> Reinitalizing axes at simulation time",&
                        Time_Simulation
                   call set_axes(Time_Simulation,.true.)
                endif
             end if
          end if
       end if

       !\
       ! Initialize solution block.
       !/
       call set_ICs

    end do ! Multi-block initialization loop.

    !^CFG IF USERFILES BEGIN
    !\
    ! Allow the user to load a restart file and then add a perturbation to it
    ! before restarting a run. 
    !/
    if (UseUserPerturbation) then
       call user_initial_perturbation
       UseUserPerturbation=.false.
    end if
    !^CFG END USERFILES   
 

    if (restart) then
       call MPI_BCAST(n_step,1,MPI_INTEGER,0,iComm,iError)
       call MPI_BCAST(Time_Simulation,1,MPI_REAL,0,iComm,iError)
!!!       call MPI_BCAST(tSimulation,1,MPI_DOUBLE_PRECISION,0,iComm,iError)
    end if

    !^CFG IF CONSTRAINB BEGIN
    ! Ensure zero divergence for the CT scheme
    if(UseConstrainB)then
       if(restart_Bface)then
          DoInitConstrainB=.false.
       else
          call init_constrain_b
       end if
    end if
    !^CFG END CONSTRAINB

  end subroutine set_initial_conditions

  subroutine initialize_files

    logical :: delete_file

    if (save_satellite_data .and. is_proc0()) &
         call open_satellite_output_files

    n_output_last=n_step
    where(dt_output>0.)
       t_output_last=int(time_simulation/dt_output)
    end where
    plot_type(restart_)='restart'
    plot_type(logfile_)='logfile'

  end subroutine initialize_files

  !===========================================================================

  subroutine save_files

    DoExchangeAgain = .false.
    do ifile=1,nfile
       if(dn_output(ifile)>=0)then
          if(dn_output(ifile)==0)then
             call save_file
          else if(mod(n_step,dn_output(ifile))==0)then
             call save_file
          end if
       else if(time_accurate .and. dt_output(ifile)>0.)then
          if(int(time_simulation/dt_output(ifile))>t_output_last(ifile))then
             t_output_last(ifile)=int(time_simulation/dt_output(ifile))
             call save_file
          end if
       end if
    end do
    ! If message passing with corners was done in save_file for tecplot plots,
    !   then do exchange_messages over again to get expected values in ghost cells.
    if(DoExchangeAgain)then
       if(is_proc0().and.lVerbose>0)write(*,*)&
            '  Calling exchange_messages to reset ghost cells ...'
       call exchange_messages
    end if

  end subroutine save_files

  !===========================================================================

  subroutine save_file
    use ModParallel, ONLY : UsePlotMessageOptions
    integer :: iFileLoop
    logical :: IsLogicalTmp
    if(n_step<=n_output_last(ifile) .and. dn_output(ifile)/=0) return

    if(ifile==restart_) then
       ! Case for restart file
       if(.not.save_restart_file)return

       !!! Added for CON
       call CON_save_restart

       call timing_start('save_restart')
       call write_octree_file
       if(is_proc0())call write_restart_header
       do globalBLK = 1,nBlockMax
          if (.not.unusedBLK(globalBLK)) call write_restart_file
       end do
       call timing_stop('save_restart')
       !^CFG  IF IONOSPHERE BEGIN
       if (UseIonosphere.and.is_proc0()) then
          call timing_start('save_iono_restart')
          call ionosphere(n_step,iono_save_restart)
          call timing_stop('save_iono_restart')
       end if
       !^CFG END IONOSPHERE

    elseif(ifile==logfile_) then
       ! Case for logfile 
       if(.not.save_logfile)return
       call timing_start('save_logfile')
       call write_logfile(0,ifile)
       call timing_stop('save_logfile')

    elseif(ifile>plot_ .and. ifile<=plot_+nplotfile) then
       ! Case for plot files
       IsFound=.false.
       !^CFG IF NOT SIMPLE BEGIN
       if(index(plot_type(ifile),'los')>0 .and. (.not. IsFound)) then
          IsFound = .true.
          call write_plot_lineofsight(ifile)
       end if
       !^CFG END SIMPLE 
       !^CFG IF IONOSPHERE BEGIN
       if(index(plot_type(ifile),'ion')>0 .and. (.not. IsFound)) then
          IsFound=.true. 
          if (is_proc0() .and. UseIonosphere) then
             call timing_start('save_iono')
             call ionosphere_write_output(ifile)
             call timing_stop('save_iono')
          endif  
       end if
       !^CFG END IONOSPHERE
       if(plot_type(ifile)/='nul' .and. (.not. IsFound) ) then
          ! Do message passing with corners once for plots
          if(.not.DoExchangeAgain)then
             if(  index(plot_form(ifile),'tec')>0 .or. &
                  index(plot_form(ifile),'tp1')>0 ) then
                if(is_proc0().and.lVerbose>0)write(*,*)&
                     '  Message passing for plot files ...'
                UsePlotMessageOptions = .true.
                call exchange_messages
                UsePlotMessageOptions = .false.
                DoExchangeAgain = .true.
                if(index(plot_form(ifile),'tp1')>0)call assign_node_numbers
             end if
          end if
          !^CFG  IF RAYTRACE BEGIN
          ! Do raytrace once if any plots need it
          do iFileLoop=1,ifile
             if(iFileLoop<ifile)then
                if(index(plot_type(iFileLoop),'ray')>0) EXIT
             else
                if(index(plot_type(ifile),'ray')>0) then
                   if(is_proc0().and.lVerbose>0)write(*,*)&
                        '  Ray tracing for plot files ...'
                   call timing_start('Ray Tracing')
                   call ray_trace(check_rayloop)
                   call timing_stop('Ray Tracing')
                end if
             end if
          end do
          !^CFG END RAYTRACE
          call timing_start('save_plot')
          call write_plot_common(ifile)
          call timing_stop('save_plot')
       end if
    elseif(ifile>satellite_ .and. ifile<=satellite_+nsatellite) then
       ! Case for satellite files
       if(.not.save_satellite_data)return
       call set_satellite_positions(IsLogicalTmp)
       call set_satellite_flags
       call write_logfile(ifile-satellite_,ifile)

    end if

    n_output_last(ifile)=n_step

    if(is_proc0() .and. lVerbose>0 .and. (ifile /= logfile_ .and. &
         (.not. (ifile > satellite_ .and. &
         ifile<=satellite_+maxsatellitefile))))then
       if(time_accurate)then
          write(*,'(a,i2,a,a,a,i7,a,i4,a,i2.2,a,i2.2,a)') &
               'saved ifile=',ifile,' type=',plot_type(ifile),&
               ' at n_step=',n_step,' time=', &
               int(                            Time_Simulation/3600.),':', &
               int((Time_Simulation-(3600.*int(Time_Simulation/3600.)))/60.),':', &
               int( Time_Simulation-(  60.*int(Time_Simulation/  60.))), &
               ' h:m:s'
       else
          write(*,'(a,i2,a,a,a,i7)') &
               'saved ifile=',ifile,' type=',plot_type(ifile),&
               ' at n_step=',n_step
       end if
    end if

  end subroutine save_file

  !===========================================================================

  subroutine save_files_final

    DoExchangeAgain = .false.
    do ifile=1,plot_+nplotfile
       call save_file
    end do

    call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------

    !^CFG  IF IONOSPHERE BEGIN
    !\
    ! Deallocate magnetosphere FAC memory as necessary.
    !/
    if (UseIonosphere) call magnetosphere_deallocate
    !^CFG END IONOSPHERE

    !\
    ! Close files
    !/
    if (save_satellite_data .and. is_proc0()) &
         call close_satellite_output_files

    if (save_logfile.and.is_proc0().and.unit_log>0) close(unit_log)

  end subroutine save_files_final

  !===========================================================================

  subroutine amr_refinement
    use ModAdvance, ONLY : tmp1_BLK

    !\
    ! Adaptive Mesh Refinement (AMR):
    ! Refine and coarsen the solution grid (on the fly) as needed.
    !/
    !\
    ! Output timing before AMR.
    !/
    if (dn_timing > -2) call timing_report
    call timing_reset_all
    call timing_start('amr')
    if(is_proc0().and.lVerbose>0)then
       write(*,*) '>>>>>>>>>>>>>>>>>>>> AMR <<<<<<<<<<<<<<<<<<<<'
       if (time_accurate) call write_timeaccurate
    end if
    !\
    ! Write plotfiles before AMR?
    !/
    if(save_plots_amr)then
       do ifile=plot_+1,plot_+nplotfile
          call save_file
       end do
    end if

    !\
    ! Perform the AMR.
    !/
    if (.not. automatic_refinement) idepth = idepth + 1

    ! BDF2 scheme should not use previous step after AMR  !^CFG IF IMPLICIT
    n_prev = -100                                         !^CFG IF IMPLICIT

    if(UseConstrainB)call b_face_fine_pass     !^CFG IF CONSTRAINB

    call amr(idepth)
    call analyze_neighbors                     !^CFG IF DEBUGGING    
    call find_test_cell

    !^CFG IF CONSTRAINB BEGIN
    if(UseConstrainB)then
       !Check for divb
       call proj_get_divb(tmp1_BLK)

       divbmax_now=maxval_loc_abs_BLK(nProc,tmp1_BLK,iLoc_I)
       if(is_proc0().and.lVerbose>0)then
          write(*,*)
          write(*,*)'Maximum of |div B| after AMR=',divbmax_now
          write(*,*)
       end if
       if(iProc==iLoc_I(5).and.divbmax_now>cTiny)write(*,*)&
            'divB,loc,x,y,z=',divbmax_now,iLoc_I,&
            x_BLK(iLoc_I(x_),iLoc_I(y_),iLoc_I(z_),iLoc_I(4)),&
            y_BLK(iLoc_I(x_),iLoc_I(y_),iLoc_I(z_),iLoc_I(4)),&
            z_BLK(iLoc_I(x_),iLoc_I(y_),iLoc_I(z_),iLoc_I(4))
    end if
    !^CFG END CONSTRAINB

    !^CFG  IF RAYTRACE BEGIN
    if(UseRaytrace)then
       call timing_start('Ray Tracing')
       call ray_trace(.true.)
       call timing_stop('Ray Tracing')
    end if
    !^CFG END RAYTRACE

    !\
    ! Output timing after AMR.
    !/
    call timing_stop('amr')
    call timing_show('amr',1)
    call timing_show('load_balance',1)
    if(is_proc0().and.lVerbose>0)&
         write(*,*) '>>>>>>>>>>>>>>>>>>>>     <<<<<<<<<<<<<<<<<<<<'

  end subroutine amr_refinement

  !^CFG IF CONSTRAINB BEGIN
  subroutine init_constrain_b
    use ModAdvance, ONLY : Bx_BLK,By_BLK,Bz_BLK,tmp1_BLK

    ! This should be a use ModInterface in GM/IH
    interface
       subroutine message_pass_dir(&
            idirmin,idirmax,width,sendghostcells,prolongorder,nvar,sol_BLK,&
            sol2_BLK,sol3_BLK,sol4_BLK,sol5_BLK,sol6_BLK,sol7_BLK,sol8_BLK,sol9_BLK,&
            restrictface)

         use ModSize
         implicit none

         integer, intent(in) :: idirmin,idirmax, width
         logical, intent(in) :: sendghostcells
         integer, intent(in) :: prolongorder, nvar
         real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), intent(inout) :: sol_BLK
         real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), optional, intent(inout) :: &
              sol2_BLK,sol3_BLK,sol4_BLK,sol5_BLK,sol6_BLK,sol7_BLK,sol8_BLK,sol9_BLK
         logical, optional, intent(in):: restrictface

       end subroutine message_pass_dir
    end interface


    DoInitConstrainB=.false.

    call message_pass_dir(1,3,1,.false.,1,3,Bx_BLK,By_BLK,Bz_BLK, &
         restrictface=.true.)

    do globalBLK=1,nBlockMax
       ! Estimate Bface from the centered B values
       call Bcenter2Bface
       ! Calculate energy (it is not set in set_ICs)
       ! because the projection scheme will need it
       call correctE
    end do

    call proj_get_divb(tmp1_BLK)
    divbmax_now=maxval_loc_abs_BLK(nProc,tmp1_BLK,iLoc_I)
    if(is_proc0().and.lVerbose>0)then
       write(*,*)
       write(*,*)'Maximum of |div B| before projection=',divbmax_now
       write(*,*)
    end if
    if(divbmax_now>cTiny)then
       if(iProc==iLoc_I(5))write(*,*)&
            'divB,loc,x,y,z=',divbmax_now,iLoc_I,&
            x_BLK(iLoc_I(x_),iLoc_I(y_),iLoc_I(z_),iLoc_I(4)),&
            y_BLK(iLoc_I(x_),iLoc_I(y_),iLoc_I(z_),iLoc_I(4)),&
            z_BLK(iLoc_I(x_),iLoc_I(y_),iLoc_I(z_),iLoc_I(4))

       if(is_proc0().and.lVerbose>0)then
          write(*,*)
          write(*,*)'Projecting B for CT scheme...'
       end if

       ! Do the projection with UseConstrainB true
       call project_B

       ! Check and report the accuracy of the projection
       call proj_get_divb(tmp1_BLK)
       divbmax_now=maxval_loc_abs_BLK(nProc,tmp1_BLK,iLoc_I)
       if(is_proc0().and.lVerbose>0)then
          write(*,*)
          call timing_show('projection',1)
          write(*,*)'Maximum of |div B| after projection=',divbmax_now
          write(*,*)
       end if
       if(iProc==iLoc_I(5).and.divbmax_now>cTiny)write(*,*)&
            'divB,loc,x,y,z=',divbmax_now,iLoc_I,&
            x_BLK(iLoc_I(x_),iLoc_I(y_),iLoc_I(z_),iLoc_I(4)),&
            y_BLK(iLoc_I(x_),iLoc_I(y_),iLoc_I(z_),iLoc_I(4)),&
            z_BLK(iLoc_I(x_),iLoc_I(y_),iLoc_I(z_),iLoc_I(4))
    end if

  end subroutine init_constrain_b
  !^CFG END CONSTRAINB

end program SWMF

!^CFG IF IMPLICIT BEGIN
!===========================================================================
subroutine select_stepping(DoPartSelect)

  ! Set logical arrays for implicit blocks, 
  ! set number of implicit and explicit blocks,
  ! and if DoPartSelect is true then select explicit and implicit blocks
  ! based on the stepping selection criteria.

  use CON_world
  use ModMpi

  use ModMain
  use ModGeometry, ONLY : Rmin_BLK
  use ModImplicit, ONLY : UseFullImplicit,UsePartImplicit, &
       implicitBLK,ImplCritType,explCFL,Rimplicit
  use ModParallel, ONLY : implicitBlock_BP
  implicit none

  logical, intent(in) :: DoPartSelect

  integer :: nBlockExpl, nBlockImpl
  integer :: iComm, iError
  logical :: oktest, oktest_me
  !---------------------------------------------------------------------------
  call set_oktest('select_stepping',oktest,oktest_me)

  iComm = i_comm()

  if(oktest_me)write(*,*) 'select_stepping starting with ',&
       'UseFullImplicit, UsePartImplicit, UsePartLocal, DoPartSelect=',&
       UseFullImplicit, UsePartImplicit, UsePartLocal, DoPartSelect

  if(((UsePartLocal .or. UsePartImplicit) .and. .not. DoPartSelect) &
       .or. .not. (UseImplicit .or. UsePartLocal))then
     nBlockExplALL    = nBlockALL
     nBlockImplALL    = 0
     implicitBLK(1:nBlock)           = .false.
     implicitBlock_BP(1:nBlockMax,:) = .false.

  elseif(UseFullImplicit)then
     nBlockExplALL = 0
     nBlockImplALL = nBlockALL
     implicitBLK(1:nBlock)           = .not.unusedBLK(1:nBlock)

     call MPI_ALLGATHER(implicitBLK,      nBLK, MPI_LOGICAL, &
          implicitBlock_BP, nBLK, MPI_LOGICAL, &
          iComm, iError)
  else

     if(is_proc0().and.lVerbose>0)&
          write(*,*)'select_stepping: ImplCritType=',ImplCritType

     ! Select implicitly treated blocks
     select case(ImplCritType)
     case('dt')
        ! Just checking
        if(.not.time_accurate)call stop_mpi(&
             'ImplCritType=dt is only valid in time_accurate mode')

        ! Set implicitBLK based on the time step.
        do globalBLK=1,nBlockMax
           if(unusedBLK(globalBLK))then
              implicitBLK(globalBLK)=.false.
           else
              ! Obtain the time step based on CFL condition

              ! For first iteration calculate dt_BLK when inside time loop,
              ! otherwise use the available dt_BLK from previous time step,
              ! or from the restart file, or simply 0 set in read_inputs.
              ! The latter two choices will be overruled later anyways.
              if(iteration_number==0 .and. time_loop)then
                 ! For first iteration in the time loop
                 ! calculate stable time step
                 call calc_facevalues(.false.)
                 call calc_facefluxes(.false.)
                 call calc_timestep
              end if
              ! If the smallest allowed timestep is below the fixed DtFixed
              ! then only implicit scheme will work
              implicitBLK(globalBLK) = dt_BLK(globalBLK)*explCFL <= DtFixed
           end if
        end do

        if(oktest_me)write(*,*)&
             'SELECT: unused,implicit,dt_BLK,explCFL,dt=',&
             unusedBLK(BLKtest),implicitBLK(BLKtest),dt_BLK(BLKtest),&
             explCFL,dt

     case('r','R')
        ! implicitly treated blocks are within Rimplicit and not unused
        implicitBLK(1:nBlockMax) = &
             Rmin_BLK(1:nBlockMax) <= Rimplicit .and. .not.unusedBLK(1:nBlockMax)
     case('test')
        implicitBLK(1:nBlockMax) = .false.
        if(i_proc()==PROCtest) implicitBLK(BLKtest) = .true.
     end select

     nBlockImpl = count(implicitBLK(1:nBlock))
     nBlockExpl = count(.not.(unusedBLK(1:nBlock).or.implicitBLK(1:nBlock)))

     call MPI_allreduce(nBlockImpl, nBlockImplALL, 1, MPI_INTEGER, MPI_SUM, &
          iComm, iError)

     call MPI_allreduce(nBlockExpl, nBlockExplALL, 1, MPI_INTEGER, MPI_SUM, &
          iComm, iError)

     call MPI_ALLGATHER(implicitBLK,      nBLK, MPI_LOGICAL, &
          implicitBlock_BP, nBLK, MPI_LOGICAL, &
          iComm, iError)
     if(is_proc0().and.lVerbose>0)oktest_me=.true. ! report for part implicit
  end if
  if(oktest_me)write(*,*)'select_stepping finished with ',&
       'nBlockExplALL, nBlockImplALL=',nBlockExplALL, nBlockImplALL

end subroutine select_stepping
!^CFG END IMPLICIT
!===========================================================================

subroutine SWMF_version(IsOn,NameVersion,Version)

  use CON_variables,  ONLY: VersionSwmf
  use CON_comp_param, ONLY: lNameVersion

  logical, intent(out)                      :: IsOn
  character (len=lNameVersion), intent(out) :: NameVersion
  real, intent(out)                         :: Version

  IsOn         =.true.
  NameVersion  ='SWMF by Univ. of Michigan'
  Version      = VersionSwmf

end subroutine SWMF_version

subroutine gettimestring
  use ModMain

  write(TimeH4,'(i4.4)') &
       int(                            Time_Simulation/3600.)
  write(TimeM2,'(i2.2)') &
       int((Time_Simulation-(3600.*int(Time_Simulation/3600.)))/60.)
  write(TimeS2,'(i2.2)') &
       int( Time_Simulation-(  60.*int(Time_Simulation/  60.)))

end subroutine gettimestring
!==============================================================================
subroutine show_all_comp

  use CON_world
  use CON_comp_param
  use CON_wrapper, ONLY : get_version_comp
  implicit none

  integer, parameter :: lWidth = 77

  logical                      :: IsOn, IsFound
  character (LEN=lNameVersion) :: NameVersion
  real                         :: Version

  integer :: lComp, iComp, iProc0Comp, nProcComp, iStrideComp
  !---------------------------------------------------------------------------
  ! Show framework version

  ! set code version 
  call SWMF_version(IsOn, NameVersion, Version)

  write(*,'(a)')'#'//repeat('=',lWidth)//'#'
  write(*,'(2a)') &
       '# ID  Version                                               ',&
       'nproc proc0 stride#'
  write(*,'(a)')'#'//repeat('-',lWidth)//'#'
  write(*,'(a,a,f5.2,3i6,a)')&
       '# CON ',NameVersion//' version',Version,n_proc(),0,1,' #'
  write(*,'(a)')'#'//repeat('-',lWidth)//'#'

  ! Show registered components
  do lComp = 1,n_comp()
     iComp = i_comp(lComp)

     call get_comp_info(iComp,&
          NameVersion=NameVersion,Version=Version,&
          iProcZero=iProc0Comp, nProc=nProcComp, iProcStride=iStrideComp)

     write(*,'(a,f5.2,3i6,a)')&
          '# '//NameComp_I(iComp)//'  '//NameVersion//' version',Version,&
          nProcComp,iProc0Comp,iStrideComp,' #'
  end do

  ! Show unregistered but compiled (ON) components
  IsFound = .false.
  do iComp=1,MaxComp
     if(use_comp(iComp)) CYCLE ! registered component
     call get_version_comp(iComp,IsOn,NameVersion,Version)
     if(.not.IsOn) CYCLE       ! empty component
     if(.not.IsFound)then
        write(*,'(a)')'#'//repeat('-',lWidth)//'#'
        IsFound=.true.
     end if
     write(*,'(a,f5.2,a)')&
          '# '//NameComp_I(iComp)//'  '//NameVersion//' version',Version,&
          '    not registered #'
  end do
  write(*,'(a)')'#'//repeat('=',lWidth)//'#'
end subroutine show_all_comp
