!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!               Space Weather Modeling Framework (SWMF)                !
!    Center for Space Environment Modeling, The University of Michigan !
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOI
! !TITLE:
! !AUTHORS:
! !AFFILIATION:
! !DATE:
! !INTRODUCTION:
!EOI
!-------------------------------------------------------------------------

subroutine BATS_setup
  use ModMpi
  use ModProcMH
  use ModIoUnit, ONLY: UNITTMP_
  use ModMain
  use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK
  use ModCT, ONLY : DoInitConstrainB               !^CFG IF CONSTRAINB
  use ModImplicit, ONLY : UsePartImplicit,n_prev   !^CFG IF IMPLICIT
  use ModIO
  use ModAMR, ONLY : &
       dn_refine,initial_refine_levels,nRefineLevelIC,nRefineLevel,&
       automatic_refinement
  use ModPhysics, ONLY : unitSI_t
  use ModNumConst

  implicit none

  ! Local variables

  character(len=*), parameter :: NameSub = 'BATS_setup'
  integer :: iError 
  !---------------------------------------------------------------------------

  ! Allocate and initialize variables dependent on number of PEs
  call allocate_vars  

  if(.not.IsStandAlone)call write_progress(0)

  call grid_setup   ! Restart reads integer only (no X_BLK or dx_BLK)

  call set_initial_conditions ! Restart reads all real data

  call find_test_cell

  call initialize_files

  call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------

  call write_runtime_values

contains
  !===========================================================================

  subroutine grid_setup
    !\
    ! Set up problem geometry, blocks, and grid structure.
    !/

    !LOCAL VARIABLES:
    integer :: nBlockMoved
    character(len=*), parameter :: NameSubSub = NameSub//'::grid_setup'
    logical :: local_refine(nBLK)

    !--------------------------------------------------------------------------

    call set_root_block_geometry
    call build_octree_roots   ! Initialize octree data structure.
    call find_neighbors       ! Get initial neighbor information.

    if (.not.restart) then
       ! Create initial solution block geometry.

       ! Perform initial refinement of mesh and solution blocks.
       do nRefineLevel = 1, initial_refine_levels
          if (iProc == 0.and.lVerbose>0) then
             call write_prefix; write (iUnitOut,*) NameSub, &
                  ' starting initial refinement level ',nRefineLevel
          end if
          call specify_initial_refinement(local_refine, nRefineLevel)
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

    if (iProc == 0.and.lVerbose>0)then
       call write_prefix; write (iUnitOut,*) '    total blocks = ',nBlockALL
    end if

    nRefineLevel = initial_refine_levels

    if(DoSetLevels) call set_levels

    call analyze_neighbors    !^CFG IF DEBUGGING 

  end subroutine grid_setup
  !===========================================================================

  subroutine set_initial_conditions

    !\
    ! Set intial conditions for solution in each block.
    !/

    use ModCompatibility, ONLY: calculate_dipole_tilt, SetDipoleTilt

    !LOCAL VARIABLES:
    character(len=*), parameter :: NameSubSub = &
         NameSub//'::set_initial_conditions'
    integer :: iLevel, iError
    integer :: nBlockMoved
    logical :: restart_read = .false. ! This is needed for UseNewAxes=F only

    !-------------------------------------------------------------------------
    if(.not.restart .and. nRefineLevelIC>0)then
       do iLevel=1,nRefineLevelIC
          do globalBLK = 1, nBlockMax
             call set_ICs
          end do
          !^CFG IF USERFILES BEGIN
          !\
          ! Allow the user to add a perturbation and use that for
          ! physical refinement.
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

          ! For sake of backwards compatibility
          if(.not.UseNewAxes .and. .not.restart_read .and. DoUpdateB0)then
             ! Now we have time_simulation read from the first restart file
             restart_read = .true.
             call calculate_dipole_tilt
          end if

       end if

       !\
       ! Initialize solution block.
       !/
       call set_ICs

    end do ! Multi-block initialization loop.

    !^CFG IF USERFILES BEGIN
    !\
    ! Allow the user to add a perturbation to the initial condition.
    !/
    if (UseUserPerturbation) then
       call user_initial_perturbation
       UseUserPerturbation=.false.
    end if
    !^CFG END USERFILES   

    if (restart) then
       if(iProc==0)then
          write(*,*)NameSub,' restarts at n_step,Time_Simulation=',&
               n_step,Time_Simulation
          !if(.not.IsStandAlone)then
          !   call get_physics(tSimulationOut=tSimulation)
          !   if(abs(tSimulation-Time_Simulation)>0.001) &
          !        write(*,*)NameSub,' WARNING Time_Simulation differs from ',&
          !        'tSimulation = ',tSimulation,' !!!'
          !end if
       end if
       ! ???n_step is already known, BCAST maybe for backwards compatibility???
       call MPI_BCAST(n_step,1,MPI_INTEGER,0,iComm,iError)
       ! Overwrite time_simulation with binary value read from restart file
       call MPI_BCAST(Time_Simulation,1,MPI_REAL,0,iComm,iError)
    end if

    !^CFG IF CONSTRAINB BEGIN
    ! Ensure zero divergence for the CT scheme
    if(UseConstrainB)then
       if(restart_Bface)then
          DoInitConstrainB=.false.
       else
          call BATS_init_constrain_b
       end if
    end if
    !^CFG END CONSTRAINB

  end subroutine set_initial_conditions

  !============================================================================

  subroutine initialize_files

    ! Local variables
    character(len=*), parameter :: NameSubSub = NameSub//'::initialize_files'
    logical :: delete_file

    if (save_satellite_data .and. iProc == 0) &
         call open_satellite_output_files

    plot_type(restart_)='restart'
    plot_type(logfile_)='logfile'

  end subroutine initialize_files

end subroutine BATS_setup

!=============================================================================

subroutine BATS_init_session

  use ModMain, ONLY: UseProjection                 !^CFG IF PROJECTION
  use ModMain, ONLY: UseConstrainB                 !^CFG IF CONSTRAINB
  use ModCT, ONLY : DoInitConstrainB               !^CFG IF CONSTRAINB

  implicit none

  ! Local variables
  character(len=*), parameter :: NameSub = 'BATS_init_session '
  !--------------------------------------------------------------------------
  ! Find the test cell defined by #TESTIJK or #TESTXYZ commands
  call find_test_cell

  ! Set number of explicit and implicit blocks !^CFG IF  IMPLICIT BEGIN
  ! Partially implicit/local selection will be done in each time step
  call select_stepping(.false.)                !^CFG END IMPLICIT 

  ! Ensure zero divergence for the CT scheme   !^CFG IF CONSTRAINB
  if(UseConstrainB .and. DoInitConstrainB)&    !^CFG IF CONSTRAINB
       call BATS_init_constrain_b              !^CFG IF CONSTRAINB

  ! Make sure that ghost cells are up to date
  call exchange_messages

  if(UseProjection)call project_B              !^CFG IF PROJECTION

  call BATS_save_files('INITIAL')

end subroutine BATS_init_session

!===========================================================================

subroutine BATS_advance(TimeSimulationLimit)
  !\
  ! Advance solution with one time step
  !/
  use ModKind
  use ModProcMH
  use ModMain
  use ModIO, ONLY: iUnitOut, write_prefix, save_plots_amr
  use ModAmr, ONLY: dn_refine
  use ModPhysics, ONLY: UnitSI_t
  use ModAdvance, ONLY: UseNonConservative, nConservCrit

  implicit none

  !INPUT ARGUMENTS:
  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded

  ! Local variables
  character(len=*), parameter :: NameSub = 'BATS_advance'

  integer      :: iBlock

  logical :: DoTest, DoTestMe
  !-------------------------------------------------------------------------
  call set_oktest(NameSub,DoTest,DoTestMe)

  ! We are advancing in time
  time_loop = .true.

  call BATS_select_blocks              !^CFG IF IMPLICIT

  n_step = n_step + 1
  iteration_number = iteration_number+1

  if (time_accurate) then
     call set_global_timestep
     dt = min(dt,(TimeSimulationLimit-Time_Simulation)/UnitSI_t)
  end if

  ! Calculate unsplit dB0Dt term for every time step
  if(DoUpdateB0 .and. .not.DoSplitDb0Dt)then
     call timing_start('update_B0')
     call calc_db0_dt(dt)
     call timing_stop('update_B0')
  end if

  call timing_start('advance')

  if(UseNonConservative .and. nConservCrit > 0)&
       call select_conservative

  if(UseImplicit.and.nBlockImplALL>0)then !^CFG IF IMPLICIT BEGIN
     call advance_impl
  else                                    !^CFG END IMPLICIT
     call advance_expl(.true.)
  endif                                   !^CFG IF IMPLICIT  
  
  if(UseIM)call apply_im_pressure         !^CFG IF RCM

  Time_Simulation = Time_Simulation + dt*UnitSI_t

  call exchange_messages
  
  call advect_all_points
  
  call timing_stop('advance')

  if(time_accurate)&
       call update_lagrangian_grid(&
       Time_Simulation - dt*UnitSI_t,Time_Simulation)

  if(DoTest)write(*,*)NameSub,' iProc,new n_step,Time_Simulation=',&
       iProc,n_step,Time_Simulation

  if (DoUpdateB0) then
     ! Unsplit dB0/dBt term is added every time step
     ! Split dB0/dt term is added at the dt_updateB0 frequency
     if (.not.DoSplitDb0Dt .or. &
          int(Time_Simulation/dt_UpdateB0) >  &
          int((Time_Simulation - dt*unitSI_t)/dt_UpdateB0)) &
          call update_b0
  end if

  if ( dn_refine > 0 .and. mod(n_step,dn_refine)==0 )then

     !\
     ! Output time before AMR.
     !/
     !if(IsStandAlone)then
     !   Due to the ordering this reset and report does not work well
     !   It may not be so important anyway to reset timing before the amr
     !   if (dn_timing > -2) call timing_report
     !   call timing_reset_all
     !end if
     call timing_start(NameThisComp//'_amr')
     if(iProc==0 .and.lVerbose>0)then
        call write_prefix; write(iUnitOut,*) &
             '>>>>>>>>>>>>>>>>>>>> AMR <<<<<<<<<<<<<<<<<<<<'
        if (time_accurate) call write_timeaccurate
     end if
     !\
     ! Write plotfiles before AMR?
     !/
     if(save_plots_amr)call BATS_save_files('PLOTS')
     call BATS_amr_refinement

     !\
     ! Output timing after AMR.
     !/
     call timing_stop(NameThisComp//'_amr')
     if(iProc == 0.and.lVerbose>0)then
        call timing_show(NameThisComp//'_amr',1)
        call timing_show('load_balance',1)
        call write_prefix; write(iUnitOut,*) &
             '>>>>>>>>>>>>>>>>>>>> AMR <<<<<<<<<<<<<<<<<<<<'
     end if

  end if

  if (UseProjection) call project_B    !^CFG IF PROJECTION

  call BATS_save_files('NORMAL')

end subroutine BATS_advance

!=============================================================================
subroutine BATS_amr_refinement
  !\
  ! Adaptive Mesh Refinement (AMR):
  ! Refine and coarsen the solution grid (on the fly) as needed.
  !/

  use ModProcMH
  use ModIO, ONLY: iUnitOut, write_prefix
  use ModMain, ONLY: lVerbose, x_, y_, z_
  use ModMain, ONLY: UseConstrainB                 !^CFG IF CONSTRAINB
  use ModCT, ONLY : DoInitConstrainB               !^CFG IF CONSTRAINB
  use ModImplicit, ONLY : UsePartImplicit,n_prev   !^CFG IF IMPLICIT
  use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK
  use ModAMR, ONLY : nRefineLevel, automatic_refinement
  use ModNumConst, ONLY: cTiny
  use ModAdvance, ONLY : tmp1_BLK

  implicit none

  !LOCAL VARIABLES:
  character(len=*), parameter :: NameSub = 'BATS_amr_refinement '
  real    :: divbmax_now
  real, external :: maxval_loc_abs_BLK
  integer :: ifile
  integer :: iLoc_I(5)  ! full location index
  !----------------------------------------------------------------------------

  !\
  ! Perform the AMR.
  !/
  if (.not. automatic_refinement) nRefineLevel = nRefineLevel + 1

  ! BDF2 scheme should not use previous step after AMR  !^CFG IF IMPLICIT
  n_prev = -100                                         !^CFG IF IMPLICIT

  if(UseConstrainB)call b_face_fine_pass     !^CFG IF CONSTRAINB

  call amr(nRefineLevel)
  call analyze_neighbors                     !^CFG IF DEBUGGING    
  call find_test_cell

  !^CFG IF CONSTRAINB BEGIN
  if(UseConstrainB)then
     !Check for divb
     call proj_get_divb(tmp1_BLK)

     divbmax_now=maxval_loc_abs_BLK(nProc,tmp1_BLK,iLoc_I)
     if(iProc == 0.and.lVerbose>0)then
        call write_prefix; write(iUnitOut,*)
        call write_prefix; write(iUnitOut,*) NameSub, &
             ' maximum of |div B| after AMR=',divbmax_now
        call write_prefix; write(iUnitOut,*)
     end if
     if(iProc==iLoc_I(5).and.divbmax_now>cTiny)write(*,*)&
          NameSub,' WARNING divB,loc,x,y,z=',divbmax_now,iLoc_I,&
          x_BLK(iLoc_I(x_),iLoc_I(y_),iLoc_I(z_),iLoc_I(4)),&
          y_BLK(iLoc_I(x_),iLoc_I(y_),iLoc_I(z_),iLoc_I(4)),&
          z_BLK(iLoc_I(x_),iLoc_I(y_),iLoc_I(z_),iLoc_I(4))
  end if
  !^CFG END CONSTRAINB

end subroutine BATS_amr_refinement

!^CFG IF CONSTRAINB BEGIN
!=============================================================================

subroutine BATS_init_constrain_b
  use ModProcMH
  use ModMain, ONLY: lVerbose, x_, y_, z_, globalBLK, nBlockMax
  use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK
  use ModCT, ONLY : DoInitConstrainB
  use ModNumConst, ONLY: cTiny
  use ModAdvance, ONLY : Bx_,By_,Bz_,State_VGB,tmp1_BLK
  use ModInterface, ONLY: message_pass_dir
  use ModIO, ONLY: write_prefix, iUnitOut
  implicit none

  ! Local variables
  character(len=*), parameter :: NameSub ='BATS_init_constrain_b '
  real, external :: maxval_loc_abs_BLK
  integer :: iLoc_I(5)  ! full location index
  real    :: divbmax_now
  !---------------------------------------------------------------------------
  DoInitConstrainB=.false.

  call message_pass_dir(1,3,1,.false.,1,3,Sol_VGB=State_VGB(Bx_:Bz_,:,:,:,:), &
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
          x_BLK(iLoc_I(x_),iLoc_I(y_),iLoc_I(z_),iLoc_I(4)),&
          y_BLK(iLoc_I(x_),iLoc_I(y_),iLoc_I(z_),iLoc_I(4)),&
          z_BLK(iLoc_I(x_),iLoc_I(y_),iLoc_I(z_),iLoc_I(4))
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
             x_BLK(iLoc_I(x_),iLoc_I(y_),iLoc_I(z_),iLoc_I(4)),&
             y_BLK(iLoc_I(x_),iLoc_I(y_),iLoc_I(z_),iLoc_I(4)),&
             z_BLK(iLoc_I(x_),iLoc_I(y_),iLoc_I(z_),iLoc_I(4))
     end if
  end if

end subroutine BATS_init_constrain_b

!^CFG END CONSTRAINB
!^CFG IF IMPLICIT BEGIN
!=============================================================================
subroutine BATS_select_blocks

  use ModProcMH
  use ModMain, ONLY: UsePartLocal, lVerbose
  use ModImplicit, ONLY : UsePartImplicit
  use ModIO, ONLY: write_prefix, iUnitOut
  implicit none

  !LOCAL VARIABLES:
  character(len=*), parameter :: NameSub = 'BATS_select_blocks'
  integer :: iError
  integer :: nBlockMoved
  !----------------------------------------------------------------------------

  ! Select blocks for partially local time stepping

  if(UsePartLocal)call select_stepping(.true.)

  ! Select and load balance blocks for partially implicit scheme
  if(UsePartImplicit)then                      
     ! Redo load balancing for partially implicit scheme
     call load_balance(.true.,.true.,nBlockMoved)
     if(nBlockMoved>0)then
        if(iProc == 0.and.lVerbose>0)then
           call write_prefix; write(iUnitOut,*)&
                'load_balance finished: nBlockMoved=',nBlockMoved
        end if
        call find_neighbors
        call analyze_neighbors   !^CFG IF DEBUGGING
        call find_test_cell
     end if
  end if

end subroutine BATS_select_blocks
!^CFG END IMPLICIT
!===========================================================================

subroutine BATS_save_files(TypeSaveIn)

  use ModProcMH
  use ModMain
  use ModIO
  use ModUtilities, ONLY : upper_case
  implicit none
  character(len=*), intent(in) :: TypeSaveIn

  character(len=len(TypeSaveIn)) :: TypeSave
  logical :: DoExchangeAgain, DoAssignNodeNumbers, IsFound
  integer :: iFile
  
  character(len=*), parameter :: NameSub='BATS_save_files'
  logical :: DoSaveRestartTmp
  logical :: IsTimeAccuratePrevious = .false.
  !--------------------------------------------------------------------------

  DoExchangeAgain     = .false.
  DoAssignNodeNumbers = .true.
  TypeSave = TypeSaveIn
  call upper_case(TypeSave)
  select case(TypeSave)
  case('INITIAL')
     ! Do not save for current time step (unless next if statement is true)
     n_output_last=n_step
     if(time_accurate .and. .not.IsTimeAccuratePrevious)then

        ! Save plot and log files at the beginning of a time accurate session

        where(dt_output>0.)
           ! The -1 allows a plot file to be written at the beginning 
           ! of the first time accurate session
           t_output_last=int(time_simulation/dt_output)-1
           n_output_last=n_step-1
        end where

        DoSaveRestartTmp = save_restart_file
        save_restart_file = .false.
        call save_files
        save_restart_file = DoSaveRestartTmp
     else
        ! Do not save for current time
        where(dt_output>0.) t_output_last=int(time_simulation/dt_output)
     end if
     IsTimeAccuratePrevious = time_accurate
  case('FINAL')
     save_restart_file = .false.
     call save_files_final
  case('FINALWITHRESTART')
     call save_files_final
  case('NORMAL')
     call save_files
  case('PLOT','PLOTS')
     do iFile=plot_+1, plot_+nPlotFile
        call save_file
     end do
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

  subroutine save_files

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
    ! then do exchange_messages over again to get expected values 
    ! in ghost cells.

    if(DoExchangeAgain)then
       if(iProc==0.and.lVerbose>0)then
          call write_prefix; write(iUnitOut,*)&
               'Calling exchange_messages to reset ghost cells ...'
       end if
       call exchange_messages
    end if

  end subroutine save_files

  !===========================================================================

  subroutine save_file
    use ModParallel, ONLY : UsePlotMessageOptions
    integer :: iFileLoop, iSat

    if(n_step<=n_output_last(ifile) .and. dn_output(ifile)/=0) return

    if(ifile==restart_) then
       ! Case for restart file
       if(.not.save_restart_file)return
       call timing_start('save_restart')
       call write_octree_file
       if(iProc==0)call write_restart_header
       do globalBLK = 1,nBlockMax
          if (.not.unusedBLK(globalBLK)) call write_restart_file
       end do
       if(iProc==0)call save_advected_points
       call timing_stop('save_restart')

    elseif(ifile==logfile_) then
       ! Case for logfile 

       if(.not.save_logfile)return
       call timing_start('save_logfile')
       call write_logfile(0,ifile)
       call timing_stop('save_logfile')

    elseif(ifile>plot_ .and. ifile<=plot_+nplotfile) then
       ! Case for plot files
       IsFound=.false.
       if(index(plot_type(ifile),'ion')>0)then
          if (iProc == 0) then
             call write_myname
             write(*,*)NameSub//' WARNING: only IE can write ion files!'
          end if
          RETURN
       end if

       if(.not.DoExchangeAgain .and. ( &
            index(plot_type(iFile),'lin')==1 .or. &    !^CFG IF RAYTRACE
            index(plot_type(iFile),'los')==1 .or. &    !^CFG IF NOT SIMPLE
            plot_form(iFile) == 'tec')) then

          if(iProc==0.and.lVerbose>0)then
             call write_prefix; write(iUnitOut,*)&
                  ' Message passing for plot files ...'
          end if
          UsePlotMessageOptions = .true.
          call exchange_messages
          UsePlotMessageOptions = .false.
          DoExchangeAgain = .true.
       end if

       !^CFG IF SIMPLE BEGIN
       if(index(plot_type(iFile),'los')>0) then
          IsFound = .true.
          call write_plot_los(iFile)
       end if
       !^CFG END SIMPLE

       !^CFG IF RAYTRACE BEGIN
       if(index(plot_type(iFile),'lin')>0) then
          IsFound = .true.
          call write_plot_line(iFile)
       end if
       !^CFG END RAYTRACE

       if(plot_type(ifile)/='nul' .and. .not.IsFound ) then
          ! Assign node numbers for tec plots
          if( index(plot_form(ifile),'tec')>0 .and. DoAssignNodeNumbers)then
             call assign_node_numbers
             DoAssignNodeNumbers = .false.
          end if

          !^CFG IF RAYTRACE BEGIN
          if(  index(plot_type(ifile),'ray')>0 .or. &
               index(plot_vars(ifile),'status')>0) call ray_trace
          !^CFG END RAYTRACE
          call timing_start('save_plot')
          call write_plot_common(ifile)
          call timing_stop('save_plot')
       end if
    elseif(ifile>satellite_ .and. ifile<=satellite_+nsatellite) then

       ! Case for satellite files
       if(.not.save_satellite_data)return
       iSat=ifile-satellite_
       call timing_start('save_satellite')
       if(iSat==1)call set_satellite_flags

       call write_logfile(iSat,ifile)
       call timing_stop('save_satellite')
    end if

    n_output_last(ifile)=n_step

    if(iProc==0 .and. lVerbose>0 .and. (ifile /= logfile_ .and. &
         (.not. (ifile > satellite_ .and. &
         ifile<=satellite_+maxsatellitefile))))then
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

    do ifile=1,plot_+nplotfile
       call save_file
    end do

    !\
    ! Close files
    !/
    if (save_satellite_data .and. iProc==0) &
         call close_satellite_output_files

    if (save_logfile.and.iProc==0.and.unit_log>0) close(unit_log)

  end subroutine save_files_final

end subroutine BATS_save_files


