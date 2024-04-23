!  Copyright (C) 2002 Regents of the University of Michigan
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModMain

  use PT_ModProc,       ONLY: iProc
  use PT_ModReadMhData, ONLY: DoReadMhData
  use ModUtilities,     ONLY: CON_stop

  implicit none

  SAVE
  private ! except

  ! Indicator of stand alone mode
  logical:: IsStandAlone = .false.

  ! Stopping conditions. These variables are only used in stand alone mode.
  real    :: TimeMax     = -1.0, CpuTimeMax = -1.0
  integer :: nIterMax    = -1
  logical :: UseStopFile = .true.
  logical :: IsLastRead  = .false.

  ! Logicals for actions
  !----------------------
  ! run the component
  logical :: DoRun       = .true.
  ! restart the run
  logical :: DoRestart   = .false.
  ! Methods and variables from ModReadMhData
  public  :: DoReadMhData
  logical, public :: DoTraceShock = .true.
  ! Methods and variables from this module
  public:: read_param, initialize, finalize, run, check, DoRestart,        &
       IsLastRead, UseStopFile, CpuTimeMax, TimeMax, nIterMax, IsStandAlone
contains
  !============================================================================
  subroutine read_param
    use PT_ModGrid,          ONLY: read_param_grid       => read_param
    use PT_ModOriginPoints,  ONLY: read_param_origin     => read_param
    use PT_ModReadMHData,    ONLY: read_param_mhdata     => read_param
    use PT_ModTime,          ONLY: read_param_time       => read_param

    ! Read input parameters for SP component
    use ModReadParam, ONLY: &
         read_var, read_line, read_command, i_session_read, read_echo_set

    ! aux variables
    integer:: nParticleCheck, nLonCheck, nLatCheck
    logical:: DoEcho

    ! The name of the command
    character (len=100) :: NameCommand
    ! Read the corresponding section of input file
    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    do
       if(.not.read_line() ) then
          IsLastRead = .true.
          EXIT
       end if
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case('#RESTART')
          call read_var('DoRestart', DoRestart)
          ! read parameters for each module
       case('#ORIGIN')
          if(IsStandAlone) CYCLE
          call read_param_origin
       case('#COORDSYSTEM', '#COORDINATESYSTEM', '#TESTPOS', &
            '#CHECKGRIDSIZE', '#GRIDNODE')
          ! Currently we do not need '#DOSMOOTH'
          if(i_session_read() /= 1) CYCLE
          call read_param_grid(NameCommand)
       case('#READMHDATA','#MHDATA')
          call read_param_mhdata(NameCommand)
       case('#TRACESHOCK')
          call read_var('DoTraceShock', DoTraceShock)
       case('#DORUN')
          call read_var('DoRun', DoRun)
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
          call read_var('nIterMax', nIterMax)
          call read_var('TimeMax' , TimeMax)
       case('#CPUTIMEMAX')
          call check_stand_alone
          call read_var('CpuTimeMax', CpuTimeMax)
       case('#CHECKSTOPFILE')
          call check_stand_alone
          call read_var('UseStopFile',UseStopFile)
       case('#ECHO')
          call check_stand_alone
          call read_var('DoEcho', DoEcho)
          if(iProc==0)call read_echo_set(DoEcho)
       case("#STARTTIME",'#NSTEP','#TIMESIMULATION')
          if(i_session_read() /= 1)CYCLE
          call read_param_time(NameCommand)
       case("#SETREALTIME")
          if(i_session_read() /= 1)CYCLE
          call check_stand_alone
          call read_param_time(NameCommand)
       case("#TIMEACCURATE")
          call check_stand_alone
          call read_param_time(NameCommand)
       case default
          call CON_stop(NameSub//': Unknown command '//NameCommand)
       end select
    end do
  contains
    !==========================================================================
    subroutine check_stand_alone

      ! certain options are only available for stand alone mode;
      ! check whether the mode is used and stop the code if it's no the case
      !------------------------------------------------------------------------
      if(IsStandAlone)RETURN
      call CON_stop(NameSub//': command '//trim(NameCommand)//&
           ' is only allowed in stand alone mode, correct PARAM.in')
    end subroutine check_stand_alone
    !==========================================================================
  end subroutine read_param
  !============================================================================
  subroutine initialize
    use PT_ModReadMhData,    ONLY: init_mhdata     => init
    ! initialize the model
    character(len=*), parameter:: NameSub = 'initialize'
    !--------------------------------------------------------------------------
    if(iProc==0)then
       write(*,'(a)')'PT: '
       write(*,'(a)')'PT: initialize'
       write(*,'(a)')'PT: '
    end if
    call init_mhdata
  end subroutine initialize
  !============================================================================
  subroutine finalize
    use PT_ModReadMhData, ONLY: finalize_mhdata    => finalize
    ! use SP_ModTurbulence, ONLY: finalize_turbulence => finalize

    ! finalize the model
    character(len=*), parameter:: NameSub = 'finalize'
    !--------------------------------------------------------------------------
    ! if(IsStandAlone)call stand_alone_final_restart
    call finalize_mhdata
  end subroutine finalize
  !============================================================================
  subroutine run(TimeLimit)

    use PT_ModGrid,          ONLY: get_other_state_var, copy_old_state,  &
         get_shock_location
    use PT_ModReadMhData,    ONLY: read_mh_data
    use PT_ModTime,          ONLY: SPTime, DataInputTime, iIter, IsSteadyState
    ! advance the solution in time
    real, intent(in)   :: TimeLimit
    logical, save :: IsFirstCall = .true.
    real :: Dt ! time increment in the current call

    ! write the initial background state to the output file
    !--------------------------------------------------------------------------
    if(IsFirstCall)then
       ! recompute the derived components of state vector, e.g.
       ! magnitude of magnetic field and velocity etc. Smooth if needed.
       if(.not.DoReadMhData) call get_other_state_var
       IsFirstCall = .false.
    end if
    if(IsSteadyState)then
       ! call iterate_steady_state
    else
       ! May need to read background data from files
       if(DoReadMhData)then
          ! copy some variables from the previous time step
          call copy_old_state
          ! Read the background data from file
          call read_mh_data()
          ! Read from file: MHData_VIB(0:nMHData,::) for the time moment
          ! DataInputTime
       end if
       ! recompute the derived components of state vector, e.g.
       ! magnitude of magnetic field and velocity etc. Smooth if needed.
       call get_other_state_var
       ! if no new background data loaded, don't advance in time
       if(DataInputTime <= SPTime) RETURN
       if(DoTraceShock) call get_shock_location
       ! run the model
       ! if(DoRun) call advance(min(DataInputTime, TimeLimit))
    end if
    ! update time & iteration counters
    iIter = iIter + 1
    Dt = min(DataInputTime, TimeLimit) - SPTime
    SPTime = SPTime + Dt
    ! call save_plot_all

    ! save restart if needed
    ! call check_save_restart(Dt)
  end subroutine run
  !============================================================================
  subroutine check

    use ModUtilities, ONLY: make_dir
    ! Make output directory
    character(len=*), parameter:: NameSub = 'check'
    !--------------------------------------------------------------------------

  end subroutine check
  !============================================================================
end module PT_ModMain
!==============================================================================
