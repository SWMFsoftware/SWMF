!  Copyright (C) 2002 Regents of the University of Michigan
! portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModMain
  use SP_ModProc, ONLY: iProc
  use SP_ModSize, ONLY: nDim, nLat, nLon, nNode, nParticleMax, &
       Particle_, OriginLat_, OriginLon_
  use SP_ModWrite, ONLY: set_write_param, write_output, finalize_write,&
       NamePlotDir
  use SP_ModReadMhData, ONLY: &
       set_read_mh_data_param, init_read_mh_data, read_mh_data, &
       finalize_read_mh_data, offset, DoReadMhData
  use SP_ModRestart, ONLY: save_restart=>write_restart, read_restart
  use SP_ModGrid, ONLY: nVar, copy_old_state,&
       LagrID_,X_,Y_,Z_,Rho_, Bx_,By_,Bz_,B_, Ux_,Uy_,Uz_, T_, BOld_, RhoOld_,&
       Wave1_, Wave2_, Length_, nBlock, &
       Proc_, Block_, nParticle_B, Shock_, ShockOld_, DLogRho_,&
       LatMin, LatMax, LonMin, LonMax, &
       RMin, RBufferMin, RBufferMax, RMax, ROrigin, &
       iShock_IB, iGridGlobal_IA, iNode_II, iNode_B, State_VIB, &
       FootPoint_VB, TypeCoordSystem,&
       set_grid_param, init_grid, get_node_indexes
  use SP_ModAdvance, ONLY: StartTime, iStartTime_I, &
       TimeGlobal, iIterGlobal, DoTraceShock, UseDiffusion, &
       set_momentum_param, advance
  use SP_ModDistribution, ONLY:  offset
  use ModKind, ONLY: Real8_

  implicit none

  SAVE

  private ! except

  real :: DataInputTime

  !\
  ! Stopping conditions. These variables are only used in stand alone mode.
  real   :: TimeMax = -1.0, CpuTimeMax = -1.0
  integer::nIterMax = -1
  logical:: UseStopFile = .true.
  logical:: IsLastRead=.false.
  ! Indicator of stand alone mode
  logical:: IsStandAlone=.false.
  !/

  !\
  !Timing variables
  logical:: UseTiming = .true.
  integer:: nTiming = -2
  integer:: nTimingDepth = -1
  character(len=10):: TimingStyle = 'cumm'
  !/

  ! Methods and variables from this module 
  public:: &
       read_param, initialize, finalize, run, check, save_restart, &
       TimeGlobal, iIterGlobal, DataInputTime, DoRestart,     & 
       UseTiming, nTiming, nTimingDepth, TimingStyle,         &
       IsLastRead, UseStopFile, CpuTimeMax, TimeMax, nIterMax,&
       IsStandAlone, copy_old_state

  ! Methods and variables from ModSize
  public:: &
       nDim, nLat, nLon, nNode, nParticleMax,&
       Particle_, OriginLat_, OriginLon_

  ! Methods and variables from ModGrid
  public:: &
       nVar, &
       LagrID_,X_,Y_,Z_,Rho_, Bx_,By_,Bz_,B_, Ux_,Uy_,Uz_, T_, RhoOld_, BOld_,&
       Wave1_, Wave2_, Length_, nBlock, &
       Proc_, Block_, nParticle_B, Shock_, ShockOld_,&
       LatMin, LatMax, LonMin, LonMax, &
       RMin, RBufferMin, RBufferMax, RMax, ROrigin,&
       iShock_IB, iGridGlobal_IA, iNode_II, iNode_B, State_VIB, &
       FootPoint_VB, TypeCoordSystem,& 
       get_node_indexes

  ! Methods and variables from ModReadMhData
  public:: &
       DoReadMhData, offset
  !\
  ! Logicals for actions
  !----------------------------------------------------------------------------
  ! run the component
  logical:: DoRun = .true.
  ! restart the run 
  logical:: DoRestart = .false.
  ! perform initialization
  logical:: DoInit = .true.
  !/
contains

  subroutine read_param(TypeAction)
    use ModTimeConvert, ONLY: time_int_to_real
    ! Read input parameters for SP component
    use ModReadParam, ONLY: &
         read_var, read_line, read_command, i_session_read, read_echo_set
    character (len=*), intent(in)     :: TypeAction ! What to do  

    ! aux variables 
    integer:: nParticleCheck, nLonCheck, nLatCheck
    logical:: DoEcho
    ! The name of the command
    character (len=100) :: NameCommand
    character (len=*), parameter :: NameSub='SP:read_param'
    !--------------------------------------------------------------------------
    ! Read the corresponding section of input file
    do
       if(.not.read_line() ) then
          IsLastRead = .true.
          EXIT
       end if
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case('#RESTART')
          call read_var('DoRestart',DoRestart)
       case('#CHECKGRIDSIZE')
          call read_var('nParticleMax',nParticleCheck)
          call read_var('nLon',     nLonCheck)
          call read_var('nLat',     nLatCheck)
          if(iProc==0.and.any(&
               (/nParticleMax,     nLon,     nLat/) /= &
               (/nParticleCheck,nLonCheck,nLatCheck/)))then
             write(*,*)'Code is compiled with nParticleMax,nLon,nLat=',&
                  (/nParticleMax, nLon, nLat/)
             call CON_stop(&
                  'Change nParticle,nLon,nLat with Config.pl -g & recompile!')
          end if
       case('#NSTEP')
          call read_var('nStep',iIterGlobal)
       case('#TIMESIMULATION')
          call read_var('tSimulation',TimeGlobal)
          DataInputTime = TimeGlobal
       case('#GRID', '#ORIGIN')
          call set_grid_param(NameCommand)
       case('#DORUN')
          call read_var('DoRun',DoRun)
       case('#SAVEPLOT','#USEDATETIME','#SAVEINITIAL')
          call set_write_param(NameCommand)
       case('#READMHDATA','#MHDATA')
          call set_read_mh_data_param(NameCommand)
       case('#COORDSYSTEM','#COORDINATESYSTEM')
          call read_var('TypeCoordSystem',TypeCoordSystem,IsUpperCase=.true.)
       case('#INJECTION')
          if(i_session_read() /= 1)CYCLE
          call set_momentum_param
       case('#TIMING')
          if(i_session_read() /= 1)CYCLE
          call read_var('UseTiming',UseTiming)
          if(.not.UseTiming)&
               CYCLE
          call read_var('DnTiming',nTiming)
          call read_var('nDepthTiming',nTimingDepth)
          call read_var('TypeTimingReport',TimingStyle)
       case('#TEST')
          ! various test modes allowing to disable certain features
          call read_var('DoTraceShock', DoTraceShock)
          call read_var('UseDiffusion', UseDiffusion)
       case('#END')
          call check_stand_alone
          IsLastRead=.true.
          EXIT
       case('#RUN')
          call check_stand_alone
          IsLastRead=.false.
          EXIT
       case('#STOP')
          call check_stand_alone
          call read_var('MaxIteration',nIterMax)
          call read_var('tSimulationMax',TimeMax)
       case('#CPUTIMEMAX')
          call check_stand_alone
          call read_var('CpuTimeMax',CpuTimeMax)
       case('#CHECKSTOPFILE')
          call check_stand_alone
          call read_var('UseStopFile',UseStopFile)
       case('#DESCRIPTION')
          call check_stand_alone
       case('#ECHO')
          call check_stand_alone
          call read_var('DoEcho', DoEcho)
          if(iProc==0)call read_echo_set(DoEcho)
       case("#STARTTIME", "#SETREALTIME")
          if(i_session_read() /= 1)CYCLE
          call check_stand_alone
          call read_var('iYear'  ,iStartTime_I(1))
          call read_var('iMonth' ,iStartTime_I(2))
          call read_var('iDay'   ,iStartTime_I(3))
          call read_var('iHour'  ,iStartTime_I(4))
          call read_var('iMinute',iStartTime_I(5))
          call read_var('iSecond',iStartTime_I(6))
          iStartTime_I(7) = 0
          call time_int_to_real(iStartTime_I, StartTime)
       case default
          call CON_stop(NameSub//': Unknown command '//NameCommand)
       end select
    end do
  contains
    subroutine check_stand_alone
      ! certain options are only available for stand alone mode;
      ! check whether the mode is used and stop the code if it's no the case
      !------------------------------------------------------------------------
      if(IsStandAlone)&
           RETURN
      call CON_stop(NameSub//': command '//trim(NameCommand)//&
           ' is only allowed in stand alone mode, correct PARAM.in')
    end subroutine check_stand_alone

  end subroutine read_param

  !============================================================================

  subroutine initialize
    use SP_ModDistribution, ONLY: init_distribution_function
    ! initialize the model
    character(LEN=*),parameter:: NameSub='SP:initialize'
    !--------------------------------------------------------------------------
    if(DoInit)then
       DoInit=.false.
    else
       RETURN
    end if
    call init_grid(DoRestart .or. DoReadMhData)
    call init_read_mh_data ! if input files are used, TimeGlobal is set here
    call init_distribution_function 
    if(DoRestart) call read_restart
    DataInputTime = TimeGlobal
  end subroutine initialize

  !============================================================================

  subroutine finalize
    ! finalize the model
    character(LEN=*),parameter:: NameSub='SP:finalize'
    !------------------------------------------------------------------------
    write(*,*)'Finalize'
    call finalize_write
    call finalize_read_mh_data
  end subroutine finalize

  !============================================================================

  subroutine run(TimeInOut, TimeLimit)
    use SP_ModGrid, ONLY: get_other_state_var
    ! advance the solution in time
    real, intent(inout):: TimeInOut
    real, intent(in)   :: TimeLimit
    logical, save:: IsFirstCall = .true.
    !------------------------------
    !\
    ! write the initial background state to the output file
    !/
    if(IsFirstCall)then
       ! print the initial state
       call write_output(IsInitialOutputIn = .true.)
       IsFirstCall = .false.
    end if

    !\
    ! May need to read background data from files
    !/
    if(DoReadMhData)then
       !\
       ! copy some variables from the previous time step
       !/
       call copy_old_state
       !\
       ! Read the background data from file
       !/
       call read_mh_data(DataInputTime)
       !Read from file: State_VIB(0:nRead,::) for the time moment
       !DataInputTime
       TimeInOut = DataInputTime
    else
       !Received from coupler: : State_VIB(0:nRead,::) for the 
       !time moment DataInputTime
       TimeInOut = TimeLimit
    end if
    !\
    ! recompute the derived components of state vector, e.g. 
    ! magnitude of magnetic field and velocity etc.
    !/
    call get_other_state_var
    !\
    ! if no new background data loaded, don't advance in time
    !/
    if(DataInputTime <= TimeGlobal) RETURN
    call lagr_time_derivative
    !\
    ! run the model
    !/
    if(DoRun) call advance(min(DataInputTime,TimeLimit))

    ! update time & iteration counters
    iIterGlobal = iIterGlobal + 1
    TimeGlobal = min(DataInputTime,TimeLimit)
    call write_output
  contains
    !=====================================================================
    subroutine lagr_time_derivative
      integer:: iBlock, iParticle
      !----------------------------------------------------------------------
      do iBlock = 1, nBlock
         do iParticle = 1, nParticle_B(  iBlock)            
            ! divergence of plasma velocity

            State_VIB(DLogRho_,iParticle,iBlock) = log(&
                 State_VIB(Rho_,iParticle,iBlock)/&
                 State_VIB(RhoOld_,iParticle,iBlock))
         end do
         ! location of shock
         call get_shock_location(iBlock)
      end do
    end subroutine lagr_time_derivative
  end subroutine run
  !=================
  subroutine get_shock_location(iBlock)
    use SP_ModAdvance, ONLY: nWidth
    use SP_ModGrid,    ONLY: R_
    integer, intent(in) :: iBlock
    ! find location of a shock wave on every field line
    !--------------------------------------------------------------------------
    ! loop variable
    integer:: iSearchMin, iSearchMax
    integer:: iShockCandidate
    !--------------------------------------------------------------------------
    if(.not.DoTraceShock)then
       iShock_IB(Shock_, iBlock) = 1
       RETURN
    end if

    ! shock front is assumed to be location of max gradient log(Rho1/Rho2);
    ! shock never moves back
    iSearchMin = max(iShock_IB(ShockOld_, iBlock), 1 + nWidth )
    iSearchMax = nParticle_B(iBlock) - nWidth - 1
    iShockCandidate = iSearchMin - 1 + maxloc(&
         State_VIB(DLogRho_,iSearchMin:iSearchMax,iBlock),&
         1, MASK = State_VIB(R_,iSearchMin:iSearchMax,iBlock) > 1.2)

    if(State_VIB(DLogRho_,iShockCandidate,iBlock) > 0.0)&
         iShock_IB(Shock_, iBlock) = iShockCandidate
  end subroutine get_shock_location
  !============================================================================
  subroutine check
    use ModUtilities, ONLY: make_dir
    character(LEN=*),parameter:: NameSub='SP:check'
    !--------------------------------------------------------------------------
    ! Make output and check input directories
    if(iProc==0) call make_dir(NamePlotDir)
    !\
    ! Initialize timing
    !/
    if(iProc==0)then
       call timing_active(UseTiming)
       call timing_step(0)
       call timing_depth(nTimingDepth)
       call timing_report_style(TimingStyle)
    end if
  end subroutine check
end module SP_ModMain
