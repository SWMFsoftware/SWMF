!  Copyright (C) 2002 Regents of the University of Michigan
! portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==================================================================
module SP_ModMain
  use SP_ModProc,    ONLY: iProc
  use SP_ModSize,    ONLY: nDim, nLat, nLon, nNode, nParticleMax
  use SP_ModPlot,    ONLY: save_plot_all, NamePlotDir
  use SP_ModReadMhData, ONLY: read_mh_data, DoReadMhData
  use SP_ModRestart, ONLY: save_restart, read_restart
  use SP_ModGrid,    ONLY: copy_old_state, LagrID_, X_,  Y_, Z_,  &
       Rho_, Bx_, By_, Bz_, Ux_, Uy_, Uz_, T_, Wave1_, Wave2_,    &
       Length_, nBlock, nParticle_B, Shock_, ShockOld_, DLogRho_, &
       RhoOld_, iShock_IB, iNode_B, State_VIB, FootPoint_VB
       
  use SP_ModAdvance, ONLY: DoTraceShock, UseDiffusion, advance
  use SP_ModTime,    ONLY: SPTime, DataInputTime, iIter
  implicit none
  SAVE
  private ! except
  !\
  ! Grid size, boundaries, coordinates
  ! Starting position of field lines in Rs
  real         :: ROrigin = 2.5
  ! Size of angular grid, latitude and longitude, at origin 
  ! surface R=ROrigin
  real         :: LonMin = 0.0, LonMax = 360.0 
  real         :: LatMin = -70.0, LatMax = 70.0
  ! Lower/Upper boundary of the domain in Rs
  real         :: RScMin=-1.0, RIhMax = -1.0
  ! Boundaries of the buffer layer between SC and IH Rs
  real         :: RIhMin=-1.0, RScMax=-1.0
  !\
  ! Stopping conditions. These variables are only used in stand alone mode.
  real   :: TimeMax  = -1.0, CpuTimeMax = -1.0
  integer:: nIterMax = -1
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
       SPTime, iIter, DataInputTime, DoRestart,     & 
       UseTiming, nTiming, nTimingDepth, TimingStyle,         &
       IsLastRead, UseStopFile, CpuTimeMax, TimeMax, nIterMax,&
       IsStandAlone, copy_old_state

  ! Methods and variables from ModSize
  public:: &
       nDim, nLat, nLon, nNode, nParticleMax

  ! Methods and variables from ModGrid
  public:: &
       LagrID_,X_, Y_, Z_, Rho_, Bx_, Bz_, Ux_, Uz_, T_, &
       Wave1_, Wave2_, Length_, nBlock, nParticle_B, Shock_,   &
       ShockOld_, RScMin, RIhMin, RScMax, RIhMax,          &
       iShock_IB, iNode_B, State_VIB, FootPoint_VB

  ! Methods and variables from ModReadMhData
  public:: DoReadMhData
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

  subroutine read_param
    use SP_ModGrid    , ONLY: read_param_grid=>read_param
    use SP_ModUnit    , ONLY: read_param_unit=>read_param
    use SP_ModTime    , ONLY: read_param_time=>read_param
    use SP_ModDistribution, ONLY:read_param_dist=>read_param
    use SP_ModAdvance , ONLY: read_param_adv =>read_param
    use SP_ModPlot    , ONLY: read_param_plot=>read_param
    use SP_ModReadMHData, ONLY:read_param_mhdata=>read_param
    ! Read input parameters for SP component
    use ModReadParam, ONLY: &
         read_var, read_line, read_command, i_session_read, read_echo_set
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
          !\
          ! read parameters for each module
          !/
       case('#ORIGIN')
          if(IsStandAlone)CYCLE
          call read_var('ROrigin', ROrigin)
          call read_var('LonMin', LonMin)
          call read_var('LatMin', LatMin)
          call read_var('LonMax', LonMax)
          call read_var('LatMax', LatMax)
       case('#COMPDOMAINS','#GRID')
          if(IsStandAlone)CYCLE
          call read_var('RScMin',RScMin)
          call read_var('RIhMin', RIhMin)
          call read_var('RScMax', RScMax)
          call read_var('RIhMax',RIhMax)
       case('#COORDSYSTEM', '#COORDINATESYSTEM',&
            '#CHECKGRIDSIZE')
          if(i_session_read() /= 1)CYCLE
          call read_param_grid(NameCommand)
       case('#PARTICLEENERGYUNIT')
          if(i_session_read() /= 1)CYCLE
          call read_param_unit(NameCommand)
       case('#MOMENTUMGRID')
          if(i_session_read() /= 1)CYCLE
          call read_param_dist(NameCommand)
       case('#INJECTION','#CFL')
          call read_param_adv(NameCommand)
       case('#SAVEPLOT','#USEDATETIME','#SAVEINITIAL','#NTAG')
          call read_param_plot(NameCommand)
       case('#READMHDATA','#MHDATA')
          call read_param_mhdata(NameCommand)
       case('#DORUN')
          call read_var('DoRun',DoRun)
       case('#TIMING')
          call read_var('UseTiming',UseTiming)
          if(.not.UseTiming) CYCLE
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
          call read_var('nIterMax',nIterMax)
          call read_var('TimeMax' ,TimeMax)
       case('#CPUTIMEMAX')
          call check_stand_alone
          call read_var('CpuTimeMax',CpuTimeMax)
       case('#CHECKSTOPFILE')
          call check_stand_alone
          call read_var('UseStopFile',UseStopFile)
       case('#ECHO')
          call check_stand_alone
          call read_var('DoEcho', DoEcho)
          if(iProc==0)call read_echo_set(DoEcho)
       case('#NSTEP','#TIMESIMULATION')
          if(i_session_read() /= 1)CYCLE
          call read_param_time(NameCommand)
       case("#STARTTIME", "#SETREALTIME")
          if(i_session_read() /= 1)CYCLE
          call check_stand_alone
          call read_param_time(NameCommand)
       case default
          call CON_stop(NameSub//': Unknown command '//NameCommand)
       end select
    end do
  contains
    subroutine check_stand_alone
      ! certain options are only available for stand alone mode;
      ! check whether the mode is used and stop the code if it's no the case
      !------------------------------------------------------------------------
      if(IsStandAlone)RETURN
      call CON_stop(NameSub//': command '//trim(NameCommand)//&
           ' is only allowed in stand alone mode, correct PARAM.in')
    end subroutine check_stand_alone

  end subroutine read_param
  !============================================================================
  subroutine initialize
    use SP_ModGrid        , ONLY: init_grid=>init
    use SP_ModUnit        , ONLY: init_unit=>init 
    use SP_ModDistribution, ONLY: init_dist=>init  
    use SP_ModAdvance     , ONLY: init_advance=>init
    use SP_ModPlot        , ONLY: init_plot=>init
    use SP_ModReadMhData  , ONLY: init_mhdata=>init
    ! initialize the model
    character(LEN=*),parameter:: NameSub='SP:initialize'
    !--------------------------------------------------------------------
    call init_grid
    call init_unit
    call init_dist
    call init_advance
    call init_plot
    call init_mhdata 
    if(DoRestart) call read_restart
    if((.not.IsStandAlone).and.(.not.DoRestart).and.(.not.DoReadMhData))&
         call get_origin_points
    DoInit=.false.
  contains
    subroutine get_origin_points
      use ModNumConst,       ONLY: cDegToRad
      use ModCoordTransform, ONLY: rlonlat_to_xyz
      use SP_ModGrid,        ONLY: Block_, Proc_, iGridGlobal_IA, iNode_II 
      integer:: iLat, iLon, iNode, iBlock      
      !Sell size on the origin surface, per line
      real         ::  DLon, DLat
      !--------------------------
      ! convert angels from degrees to radians
      LonMax = LonMax*cDegToRad
      LonMin = LonMin*cDegToRad
      ! angular grid's step
      DLon = (LonMax - LonMin)/nLon
      
      ! convert angels from degrees to radians
      LatMax = LatMax*cDegToRad
      LatMin = LatMin*cDegToRad
      ! angular grid's step
      DLat = (LatMax - LatMin)/nLat
      
      
      do iLat = 1, nLat
         do iLon = 1, nLon
            iNode = iNode_II(iLon, iLat)
            iBlock = iGridGlobal_IA(Block_, iNode)
            if(iProc == iGridGlobal_IA(Proc_, iNode))then
               nParticle_B(iBlock) = 1
               call rlonlat_to_xyz((/ROrigin, LonMin + (iLon - 0.5)*DLon, &
                    LatMin + (iLat - 0.5)*DLat/), State_VIB(X_:Z_,1,iBlock))
            end if
         end do
      end do
    end subroutine get_origin_points
  end subroutine initialize
  !============================================================================
  subroutine finalize
    use SP_ModPlot,       ONLY: finalize_plot=>finalize
    use SP_ModReadMhData, ONLY: finalize_read=>finalize
    ! finalize the model
    character(LEN=*),parameter:: NameSub='SP:finalize'
    !------------------------------------------------------------------------
    call finalize_plot
    call finalize_read
  end subroutine finalize
  !============================================================================

  subroutine run(TimeLimit)
    use SP_ModGrid, ONLY: get_other_state_var
    ! advance the solution in time
    real, intent(in)   :: TimeLimit
    logical, save:: IsFirstCall = .true.
    !------------------------------
    !\
    ! write the initial background state to the output file
    !/
    if(IsFirstCall)then
       ! print the initial state
       call save_plot_all(IsInitialOutputIn = .true.)
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
       call read_mh_data()
       !Read from file: State_VIB(0:nMHData,::) for the time moment
       !DataInputTime
    end if
    !\
    ! recompute the derived components of state vector, e.g. 
    ! magnitude of magnetic field and velocity etc.
    !/
    call get_other_state_var
    !\
    ! if no new background data loaded, don't advance in time
    !/
    if(DataInputTime <= SPTime) RETURN
    call lagr_time_derivative
    !\
    ! run the model
    !/
    if(DoRun) call advance(min(DataInputTime,TimeLimit))

    ! update time & iteration counters
    iIter = iIter + 1
    SPTime = min(DataInputTime,TimeLimit)
    call save_plot_all
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
    use SP_ModGrid,    ONLY: R_, NoShock_
    integer, intent(in) :: iBlock
    !\
    ! find location of a shock wave on a given line (block)
    !/
    !\
    ! Do not search too close to the Sun
    integer         :: iShockMin
    real, parameter :: RShockMin = 1.20  !*RSun
    ! Do not search too close to the heliosphere boundary
    integer:: iShockMax
    ! Misc
    integer:: iShockCandidate
    !---------------------------------------------------------------------
    if(.not.DoTraceShock)then
       iShock_IB(Shock_, iBlock) = NoShock_
       RETURN
    end if

    ! shock front is assumed to be location of max gradient log(Rho1/Rho2);
    ! shock never moves back
    iShockMin = max(iShock_IB(ShockOld_, iBlock), 1 + nWidth )
    iShockMax = nParticle_B(iBlock) - nWidth - 1
    iShockCandidate = iShockMin - 1 + maxloc(&
         State_VIB(DLogRho_,iShockMin:iShockMax,iBlock),&
         1, MASK = State_VIB(R_,iShockMin:iShockMax,iBlock) > RShockMin)

    if(State_VIB(DLogRho_,iShockCandidate,iBlock) > 0.0)&
         iShock_IB(Shock_, iBlock) = iShockCandidate
  end subroutine get_shock_location
  !=======================================================================
  subroutine check
    use ModUtilities, ONLY: make_dir
    character(LEN=*),parameter:: NameSub='SP:check'
    !---------------------------------------------------------------------
    if(.not.IsStandAlone)then
       ! check if the domains for SC and IH are correctly set
       if(RScMin < 0.0.or. RIhMin < 0.0 .or.RScMax < 0.0 .or.RIhMax < 0.0)&
            call CON_stop(NameSub//&
            ': RScMin, RScMax, RIhMin, RIhMax must be set to positive values')
       if(any((/RIhMax, RScMax/) < RIhMin) .or. RIhMax < RScMax)&
            call CON_stop(NameSub//&
           ': inconsistent values of RIhMin, RScMax, RIhMax')
       if(.not.(DoRestart.or.DoReadMhData))then
          ! check consistency of the origin point parameters
          if(LonMax <= LonMin .or. LatMax <= LatMin)&
               call CON_stop(NameSub//': Origin surface grid is inconsistent')
          if(ROrigin <= RScMin)call CON_stop(NameSub//&
               ': ROrigin, if set, must be greater than RScMin')
          if(any((/RIhMin, RScMax, RIhMax/) < ROrigin))&
               call CON_stop(NameSub//&
               ': inconsistent values of ROrigin, RIhMin, RScMax, RIhMax')
       end if
    else
       if(DoRestart)call CON_stop(NameSub//&
            'Restart save/read are not implemented for the stand-alone mode')
    end if
    ! Make output directory
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
