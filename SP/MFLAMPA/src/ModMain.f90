!  Copyright (C) 2002 Regents of the University of Michigan
! portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModMain
  use SP_ModProc, ONLY: iProc

  use SP_ModSize, ONLY: &
       nDim, nLat, nLon, nNode, nParticleMax, &
       Particle_, OriginLat_, OriginLon_
  
  use SP_ModWrite, ONLY: &
       set_write_param, write_output, NamePlotDir

  use SP_ModReadMhData, ONLY: &
       set_read_mh_data_param, init_read_mh_data, read_mh_data, DoReadMhData

  use SP_ModRestart, ONLY: &
       save_restart=>write_restart, read_restart

  use SP_ModGrid, ONLY: &
       nVar, copy_old_state,&
       LagrID_,X_,Y_,Z_,Rho_, Bx_,By_,Bz_,B_, Ux_,Uy_,Uz_, T_, BOld_, RhoOld_,&
       Wave1_, Wave2_, Length_, nBlock, &
       Proc_, Block_, nParticle_B, Shock_, ShockOld_,&
       LatMin, LatMax, LonMin, LonMax, &
       RMin, RBufferMin, RBufferMax, RMax, ROrigin, &
       iShock_IB, iGridGlobal_IA, iNode_II, iNode_B, State_VIB, &
       FootPoint_VB, TypeCoordSystem,&
       set_grid_param, init_grid, get_node_indexes
  
  use SP_ModAdvance, ONLY: &
       TimeGlobal, iIterGlobal, DoTraceShock, UseDiffusion, &
       Distribution_IIB, advance, set_momentum_param

  implicit none

  SAVE

  private ! except

  real :: DataInputTime

  !\
  ! Stopping conditions. These variables are only used in stand alone mode.
  real    :: TimeMax = -1.0, CpuTimeMax = -1.0
  integer ::nIterMax = -1
  logical :: UseStopFile = .true.
  logical :: IsLastRead=.false.
  !/

  !Timing variables
  logical:: UseTiming = .true.
  integer:: nTiming = -2
  integer:: nTimingDepth = -1
  character(len=10):: TimingStyle = 'cumm'

  ! Methods and variables from this module 
  public:: &
       read_param, initialize, run, check, save_restart,      &
       TimeGlobal, iIterGlobal, DataInputTime, DoRestart,     & 
       UseTiming, nTiming, nTimingDepth, TimingStyle,         &
       IsLastRead, UseStopFile, CpuTimeMax, TimeMax, nIterMax,&
       copy_old_state, offset

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
       Distribution_IIB, FootPoint_VB, TypeCoordSystem,& 
       get_node_indexes

  ! Methods and variables from ModWrite

  ! Methods and variables from ModAdvance

  ! Methods and variables from ModReadMhData
  public:: &
       DoReadMhData

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
    ! Read input parameters for SP component
    use ModReadParam, ONLY: read_var, read_line, read_command, i_session_read
    character (len=*), intent(in)     :: TypeAction ! What to do  

    ! aux variables 
    integer:: nParticleCheck, nLonCheck, nLatCheck
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
       case('#SAVEPLOT')
          call set_write_param
       case('#READMHDATA','#MHDATA')
          call set_read_mh_data_param(NameCommand)
       case('#COORDSYSTEM',"#COORDINATESYSTEM")
          call read_var('TypeCoordSystem',TypeCoordSystem,IsUpperCase=.true.)
       case('#INJECTION')
          call set_momentum_param
       case("#TIMING")
          if(i_session_read() /= 1)&
               CYCLE
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
       case("#END")
          IsLastRead=.true.
          EXIT
       case("#RUN")
          IsLastRead=.false.
          EXIT
       case("#STOP")
          call read_var('MaxIteration',nIterMax)
          call read_var('tSimulationMax',TimeMax)
       case("#CPUTIMEMAX")
          call read_var('CpuTimeMax',CpuTimeMax)
       case("#CHECKSTOPFILE")
          call read_var('UseStopFile',UseStopFile)
       case default
          call CON_stop(NameSub//': Unknown command '//NameCommand)
       end select
    end do
  end subroutine read_param

  !============================================================================

  subroutine initialize
    use SP_ModAdvance, ONLY: init_distribution_function
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
       !Reaceived from coupler: : State_VIB(0:nRead,::) for the 
       !time moment DataInputTime
       TimeInOut = TimeLimit
    end if
    !\
    ! recompute the derived components of state vector, e.g. 
    ! magnitude of magnetic field and velocity etc.
    !/
    call get_other_state_var
    call fix_grid_consistency   !????

    !\
    ! if no new background data loaded, don't advance in time
    !/
    if(DataInputTime <= TimeGlobal) RETURN
    !\
    ! run the model
    !/
    if(DoRun) call advance(min(DataInputTime,TimeLimit))

    ! update time & iteration counters
    iIterGlobal = iIterGlobal + 1
    TimeGlobal = min(DataInputTime,TimeLimit)
    call write_output(IsInitialOutputIn=.not.DoRun)
  contains
    !=====================================================================
    subroutine fix_grid_consistency
      use SP_ModGrid, ONLY: DLogRho_
      integer:: iBlock, iParticle
      !----------------------------------------------------------------------
      do iBlock = 1, nBlock
         do iParticle = 1, nParticle_B(  iBlock)            
            ! divergence of plasma velocity
            if(DataInputTime > TimeGlobal) then
               State_VIB(DLogRho_,iParticle,iBlock) = log(&
                    State_VIB(Rho_,iParticle,iBlock)/&
                    State_VIB(RhoOld_,iParticle,iBlock))
            end if
         end do
         ! location of shock
         iShock_IB(:, iBlock) = max( iShock_IB(:, iBlock), 1)
      end do
    end subroutine fix_grid_consistency
  end subroutine run
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
  !============================================================================
  subroutine offset(iBlock, iOffset, AlphaIn)
    ! shift in the data arrays is required if the grid point(s) is  
    ! appended or removed at the foot point of the magnetic field line
    !SHIFTED ARE:  State_VIB((/RhoOld_,BOld_/),:,:), Distribution_IIB,
    !ShockOld, nParticle_B
    integer, intent(in)        :: iBlock
    integer, intent(in)        :: iOffset
    real, optional, intent(in) :: AlphaIn
    real :: Alpha
    character(len=*), parameter :: NameSub = "SP: offset"
    !------------
    Alpha = 0
    if(present(AlphaIn))Alpha=AlphaIn
    if(iOffset==0)then
       RETURN
    elseif(iOffset==1)then
       State_VIB((/RhoOld_,BOld_/),2:nParticle_B(iBlock)+1,iBlock) &
            = State_VIB((/RhoOld_,BOld_/),1:nParticle_B(iBlock),iBlock)
       Distribution_IIB(:,2:nParticle_B( iBlock)+iOffset, iBlock)&
            = Distribution_IIB(:,1:nParticle_B( iBlock), iBlock)
       State_VIB((/RhoOld_, BOld_/), 1, iBlock) = &
            (Alpha + 1)*State_VIB((/RhoOld_, BOld_/), 2, iBlock) &
            -Alpha     * State_VIB((/RhoOld_, BOld_/), 3, iBlock)
       Distribution_IIB(:,1,iBlock) = Distribution_IIB(:,2,iBlock) + &
            Alpha*(Distribution_IIB(:,2,iBlock) - Distribution_IIB(:,3,iBlock))
    elseif(iOffset < 0)then
       State_VIB((/RhoOld_,BOld_/),1:nParticle_B(iBlock)+iOffset,iBlock) &
            =  State_VIB((/RhoOld_,BOld_/),1-iOffset:nParticle_B(iBlock),&
            iBlock)
       Distribution_IIB(:,1:nParticle_B( iBlock)+iOffset, iBlock)&
            = Distribution_IIB(:,1-iOffset:nParticle_B( iBlock), iBlock)
    else
       call CON_stop('No algorithm for iOffset >1 in '//NameSub)
    end if
    iShock_IB(ShockOld_, iBlock) = &
         iShock_IB(ShockOld_, iBlock) + iOffset
    nParticle_B(iBlock) = nParticle_B( iBlock) + iOffset
  end subroutine offset
end module SP_ModMain
