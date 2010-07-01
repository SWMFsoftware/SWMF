!^CMP COPYRIGHT UM
!
!BOP
!
!MODULE: CON_io - Input and output methods for CON
!
!DESCRIPTION:
! Input and output methods used by CON
!INTERFACE:
module CON_io

  !USES:
  use ModIoUnit
  use CON_world
  use CON_comp_param
  use CON_wrapper
  use CON_planet, ONLY: NamePlanet, read_planet_var, check_planet_var
  use CON_time
  use CON_axes
  use ModReadParam
  use CON_variables
  use ModUtilities, ONLY: DoFlush, flush_unit, fix_dir_name, check_dir, &
       split_string, lower_case

  implicit none

  private ! except

  !PUBLIC MEMBER FUNCTIONS:
  public :: read_inputs  !  read, broadcast and distribute parameters
  public :: save_restart ! save restart information

  !REVISION HISTORY:
  ! 07/22/03 Gabor Toth <gtoth@umich.edu> - initial prototype based on
  !                                         BATSRUS -r v7_5_0_FRAMEWORK
  ! 07/31/03 Aaron Ridley and G.Toth - added commands to read physical params
  ! 08/20/03 G.Toth - added save_restart, Protex prolog, and made it private
  ! 08/25/03 G.Toth - save_restart calls save_restart_comp for all components
  ! 05/20/04 G.Toth - added #CYCLE command
  ! 09/01/05 G.Toth - replaced call CON_stop with call world_clean and RETURN
  !                   check first read condition for each command separately.
  !EOP

  character (len=*), parameter :: NameMod='CON_io'

  ! Name of the input parameter file
  integer, parameter :: lNameFile=100
  character (LEN=*), parameter :: NameParamFile="PARAM.in"

  ! writing to STDOUT vs. log files
  logical                     :: UseStdout=.true.
  character(len=lNameFile)    :: NameStdoutDir='STDOUT/'
  character(len=*), parameter :: NameStdoutExt='.log'

  ! Saving restart information
  character(len=lNameFile)    :: NameRestartFile='RESTART.out'

  ! How often shall we save restart files?
  ! Default: every 100000 time steps or 1e30 seconds (ie at the end of the run)
  type(FreqType), public :: SaveRestart = FreqType(.true.,100000,1E30,-1,-1.0)

  ! Showing progress
  integer, public :: DnShowProgressShort = 10
  integer, public :: DnShowProgressLong  = 100

contains

  !BOP =======================================================================
  !IROUTINE: read_inputs - read, broadcast and distribute parameters
  !INTERFACE:
  subroutine read_inputs(IsLastRead)

    !USES:
    use CON_coupler, ONLY: &
         Couple_CC, MaxCouple, nCouple, iCompCoupleOrder_II, DoCoupleOnTime_C
    use CON_physics

    implicit none

    character (LEN=*), parameter:: NameSub = NameMod//'::read_inputs'
    real :: VersionRead

    !OUTPUT ARGUMENTS:
    logical, intent(out) :: IsLastRead ! True if last session was read

    !DESCRIPTION:
    ! Reads the parameters for control, and passes parameters 
    ! between 
    ! \begin{verbatim}
    ! #BEGIN_COMP PM
    ! ....
    ! #END_COMP PM
    ! \end{verbatim}
    ! to the physical component PM.
    !
    ! The parameter file consists of one or more sessions. 
    ! The \#RUN command means that the session should be executed.
    ! If the file ends or an \#END command is found, the last session
    ! is executed.
    ! The ordering of the commands is more or less arbitrary
    ! with the exception of \#PLANET (\#MOON, \#COMET) that must preceed 
    ! commands that idealize or overwrite the real parameters of the planet.
    !EOP

    ! True if first session will be read
    logical :: IsFirstRead=.true. 

    ! Logical for echo of commands
    logical :: DoEcho=.false.

    ! One line of input and command name
    character (len=lStringLine) :: StringLine, NameCommand

    ! Current line number and total number of lines
    ! Must be SAVE-d or initialized
    integer :: iLine=0

    ! Temporary variables
    integer :: nByteRealRead
    real    :: FracSecond ! Default precision is sufficient for reading

    integer            :: nDepthTiming  = -1
    character (len=10) :: TypeTimingReport  = 'cumu'
    logical            :: UseTimingAll = .false.

    ! Timing parameters
    character (lNameVersion) :: NameVersion
    real               :: VersionNumber
    logical            :: IsOn

    ! Beginning line number and number of lines for component parameters
    integer :: iLineModule, nLineModule

    ! Unit number for STDOUT from components
    integer :: iUnitOut

    ! Names, indexes and logicals for components
    integer :: lComp, iComp, iComp1, iComp2, nName, iCouple
    logical :: UseComp
    character (len=lNameComp) :: NameComp, NameComp1, NameComp2
    character (len=lNameComp) :: NameSourceTarget_I(2)
    character (len=lStringLine) :: NameSourceTarget

    !-------------------------------------------------------------------------
    if(is_proc0())write(*,'(a,i3)')NameSub//': iSession=',iSession

    !\
    ! Proc 0 reads the file and broadcasts to all PE-s
    !/
    if(IsFirstRead) call read_file(NameParamFile,i_comm())

    !\
    ! Read input data from string array stored in ModReadParam
    !/
    call read_init('  ',iSession,iLine)

    IsLastRead=.true.
    do
       if( .not. read_line(StringLine, iLine) )then
          IsLastRead=.true.
          EXIT
       end if
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case("#ECHO")

          call read_var('DoEcho',DoEcho)
          if(is_proc0())call read_echo_set(DoEcho)

       case("#FLUSH")
          call read_var('DoFlush', DoFlush)

       case("#STDOUT")

          call read_var('UseStdout',UseStdout)
          call set_stdout

       case("#STDOUTDIR")

          call read_var('NameStdoutDir',NameStdoutDir)
          call fix_dir_name(NameStdoutDir)
          if(is_proc0())call check_dir(NameStdoutDir)

       case("#DESCRIPTION")
          call read_var('StringDescription',StringDescription)

       case("#BEGIN_COMP")

          ! Obtain and check module name following #BEGIN_COMP ...
          NameComp=StringLine(13:14)
          if(.not.use_comp(NameComp))then
             if(is_proc0()) write(*,*) NameSub// &
                  ' SWMF_ERROR unregistered componet name '//trim(NameComp)
             iErrorSwmf = 11
             RETURN
          end if

          ! Find #END_COMP MODULENAME and read module parameters
          iLineModule = iLine
          MODULELINES: do
             if(.not.read_line(StringLine, iLine)) then
                if(is_proc0()) write(*,*) NameSub// &
                     ' SWMF_ERROR: could not find #END_COMP '//trim(NameComp)
                iErrorSwmf = 12
                RETURN
             end if
             if( StringLine(1:12) == '#END_COMP '//NameComp ) then
                nLineModule = iLine - iLineModule

                ! Decide where to write STDOUT
                if(UseStdout)then
                   iUnitOut=STDOUT_
                else
                   call get_comp_info(NameComp,iUnitOut=iUnitOut)
                end if
                ! Initialize ModReadParam
                call read_init(NameComp,iSession,iLineModule,iLine-1,iUnitOut)
                
                ! Echo component input on root PE of component
                if(is_proc0())        call read_echo_set(.false.)
                if(is_proc0(NameComp))call read_echo_set(DoEcho)

                ! Set the parameters for the component
                call set_param_comp(NameComp,'READ')

                if(is_proc0(NameComp))then
                   if(     UseStdout)call flush_unit(STDOUT_)
                   if(.not.UseStdout)call flush_unit(iUnitOut)
                end if

                ! reset ModReadParam to read the CON parameters
                call read_init('  ',iSession,iLine)

                ! Echo CON input on root PE of CON
                if(is_proc0(NameComp))call read_echo_set(.false.)
                if(is_proc0())then
                   call read_echo_set(DoEcho)
                   if(DoEcho)write(*,'(a)')trim(StringLine)
                end if

                EXIT MODULELINES
             endif
          end do MODULELINES

       case("#END")

          IsLastRead=.true.
          EXIT

       case("#RUN")

          IsLastRead=.false.
          EXIT

       case("#TIMEACCURATE")
          call read_var('DoTimeAccurate',DoTimeAccurate)
          if(.not. IsStandAlone .and. .not. DoTimeAccurate)then
             if(is_proc0())write(*,*)NameSub// &
                  ' SWMF_ERROR: steady-state mode is available only'// &
                  ' when the SWMF runs as stand-alone'
             if(UseStrict)then
                call world_clean
                iErrorSwmf = 20
                RETURN
             end if
             DoTimeAccurate = .true.
          endif

       case("#SAVERESTART")

          call read_var('DoSaveRestart',SaveRestart % DoThis)
          if(SaveRestart % DoThis)then
             call read_var('DnSaveRestart',SaveRestart % Dn)
             call read_var('DtSaveRestart',SaveRestart % Dt)
             SaveRestart % nNext = SaveRestart % Dn
             SaveRestart % tNext = SaveRestart % Dt
          endif

       case("#RESTARTFILE")

          call read_var('NameRestartFile',NameRestartFile)

       case("#STOP")

          if(IsStandAlone)then
             call read_var('MaxIteration'  ,MaxIteration)
             call read_var('tSimulationMax',tSimulationMax)
          else if(is_proc0()) then
             write(*,*)NameSub// &
                  ' SWMF_WARNING: stop condition has been set externally'
             write(*,*)NameSub// &
                  ' SWMF_WARNING: tSimulationMax = ',tSimulationMax
          end if

       case("#STRICT")

          call read_var('UseStrict',UseStrict)

       case("#VERBOSE")

          call read_var('lVerbose',lVerbose)

       case("#TEST")

          call read_var('StringTest',StringTest)

       case("#TIMING")

          call read_var('UseTiming',UseTiming)
          if(UseTiming)then
             call read_var('DnTiming'        ,DnTiming)
             call read_var('nDepthTiming'    ,nDepthTiming)
             call read_var('TypeTimingReport',TypeTimingReport)
             UseTimingAll = index(TypeTimingReport,'all') > 0
             TypeTimingReport = TypeTimingReport(1:4)
          end if

       case("#CPUTIMEMAX")

          call read_var('CpuTimeMax',CpuTimeMax)

       case("#CHECKSTOPFILE")

          call read_var('DoCheckStopFile',DoCheckStopFile)

       case("#CHECKSTOP")

          call read_var('DoCheckStop',CheckStop % DoThis)
          if(CheckStop % DoThis)then
             call read_var('DnCheckStop',CheckStop % Dn)
             call read_var('DtCheckStop',CheckStop % Dt)
          end if

       case("#PROGRESS")

          call read_var('DnShowProgressShort',DnShowProgressShort)
          call read_var('DnShowProgressLong', DnShowProgressLong)

       case("#COMPONENT")

          call read_var('NameComp',NameComp)
          call read_var('UseComp',UseComp)
          call put_comp_info(NameComp,Use=UseComp)

       case("#COUPLEORDER")

          call read_var('nCouple',nCouple)
          if(nCouple > MaxCouple)then
             if(is_proc0()) write(*,*)NameSub// &
                  ' SWMF_ERROR: nCouple > MaxCouple=',MaxCouple
             if(UseStrict)then
                iErrorSwmf = 13
                RETURN
             end if
             nCouple = MaxCouple
          end if
          do iCouple = 1,nCouple
             ! Read in a string with the names of source and target components
             call read_var('NameSource NameTarget',NameSourceTarget)
             ! Split the string
             call split_string(NameSourceTarget,2,&
                  NameSourceTarget_I,nName)
             ! Check the number of names
             if(nName /= 2)then
                if(is_proc0()) write(*,*) NameSub//' SWMF_ERROR: '// &
                     ' cannot find 2 component IDs in NameSourceTarget='// &
                     trim(NameSourceTarget)
                iErrorSwmf = 14
                RETURN
             end if
             if(.not.is_valid_comp_name(NameSourceTarget_I(1))) then
                if(is_proc0()) write(*,*) NameSub//' SWMF_ERROR: '// &
                     ' invalid source component name='//NameSourceTarget_I(1)
                iErrorSwmf = 15
                RETURN
             end if
             if(.not.is_valid_comp_name(NameSourceTarget_I(2))) then
                if(is_proc0()) write(*,*) NameSub//' SWMF_ERROR: '// &
                     ' invalid target component name='//NameSourceTarget_I(2)
                iErrorSwmf = 16
                RETURN
             end if

             ! Store source and target ID-s.
             iCompCoupleOrder_II(1,iCouple) = i_comp(NameSourceTarget_I(1))
             iCompCoupleOrder_II(2,iCouple) = i_comp(NameSourceTarget_I(2))
          end do

       case("#COUPLE1","#COUPLE2","#COUPLE1SHIFT","#COUPLE2SHIFT")

          call read_var('NameSource',NameComp1)
          iComp1=i_comp(NameComp1)
          if(.not.use_comp(iComp1)) then
             if(is_proc0()) write(*,*) NameSub//' SWMF_ERROR: '// &
                  NameComp1//' is OFF or not registered in '//NameMapFile
             iErrorSwmf = 17
             RETURN
          end if
          call read_var('NameTarget',NameComp2)
          iComp2=i_comp(NameComp2)
          if(.not.use_comp(iComp2))then
             if(is_proc0()) write(*,*) NameSub//' SWMF_ERROR: '// &
                  NameComp2//' is OFF or not registered in '//NameMapFile
             iErrorSwmf = 18
             RETURN
          end if

          call read_var('DnCouple',Couple_CC(iComp1,iComp2) % Dn)
          call read_var('DtCouple',Couple_CC(iComp1,iComp2) % Dt)
          Couple_CC(iComp1,iComp2) % DoThis = &
               Couple_CC(iComp1,iComp2) % Dn >= 0 .or. &
               Couple_CC(iComp1,iComp2) % Dt >= 0.0

          if(NameCommand(1:8)=="#COUPLE2")then
             ! Symmetric coupling
             Couple_CC(iComp2,iComp1)=Couple_CC(iComp1,iComp2)
          endif
          if(NameCommand(9:13)=="SHIFT")then
             ! Shift Comp1->Comp2 coupling
             call read_var('nNext12',Couple_CC(iComp1,iComp2) % nNext)
             call read_var('tNext12',Couple_CC(iComp1,iComp2) % tNext)
          end if
          if(NameCommand=="#COUPLE2SHIFT")then
             ! Shift Comp2->Comp1 coupling
             call read_var('nNext21',Couple_CC(iComp2,iComp1) % nNext)
             call read_var('tNext21',Couple_CC(iComp2,iComp1) % tNext)
          end if

       case("#COUPLETIME")
          call read_var('NameComp',NameComp)
          call read_var('DoCoupleOnTime',DoCoupleOnTime_C(i_comp(NameComp)))
       case("#CYCLE")
          call read_var('NameComp',NameComp)
          call read_var('DnRun',DnRun_C(i_comp(NameComp)))
       case("#VERSION")
          VersionRead = VersionSwmf ! set default value
          call read_var('Version',VersionRead)
          if(abs(VersionRead-VersionSwmf)>0.005.and.is_proc0())&
               write(*,'(a,f6.3,a,f6.3)') &
               NameSub//': SWMF_WARNING version in file is ',&
               VersionRead,' but SWMF version is ',VersionSwmf

       case("#PRECISION")
          call read_var('nByteReal',nByteRealRead)
          if(nByteReal/=nByteRealRead)then
             if(is_proc0()) then
                write(*,'(a,i1,a)') NameSub//' WARNING: '//&
                     'SWMF was compiled with ',nByteReal,' byte reals'
                write(*,'(a,i1,a)') NameSub//' WARNING: '// &
                     'requested precision is ',nByteRealRead,' byte reals'
             end if
             if(UseStrict)then
                if(is_proc0())write(*,*) &
                     'Change precision and recompile or set #STRICT .false.'
                iErrorSwmf = 19
                RETURN
             end if
          end if

       case("#NSTEP")
          if(.not.is_first_read())then
             if(UseStrict)RETURN
             CYCLE
          end if
          call read_var('nStep',nStep)

       case("#TIMESIMULATION")
          if(.not.is_first_read())then
             if(UseStrict)RETURN
             CYCLE
          end if
          if(IsStandAlone)then
             call read_var('tSimulation',tSimulation)
          else
             write(*,*)NameSub// &
                  ' SWMF_WARNING: simulation time has been set externally'
             write(*,*)NameSub// &
                  ' SWMF_WARNING: tSimulation = ',TimeStart % String
          end if

       case("#STARTTIME", "#SETREALTIME")
          if(.not.is_first_read())then
             if(UseStrict)RETURN
             CYCLE
          end if
          if(IsStandAlone)then
             call read_var('iYear'  ,TimeStart % iYear)
             call read_var('iMonth' ,TimeStart % iMonth)
             call read_var('iDay'   ,TimeStart % iDay)
             call read_var('iHour'  ,TimeStart % iHour)
             call read_var('iMinute',TimeStart % iMinute)
             call read_var('iSecond',TimeStart % iSecond)
             FracSecond = TimeStart % FracSecond ! set default value
             call read_var('FracSecond',FracSecond)
             TimeStart % FracSecond = FracSecond
             call time_int_to_real(TimeStart)
          else
             write(*,*)NameSub// &
                  ' SWMF_WARNING: start time has been set externally'
             write(*,*)NameSub// &
                  ' SWMF_WARNING: TimeStart = ',TimeStart % String
          end if

       case('#PLANET','#MOON','#COMET', &
            '#IDEALAXES','#ROTATIONAXIS','#MAGNETICAXIS','#MAGNETICCENTER',&
            '#ROTATION','#NONDIPOLE','#DIPOLE','#UPDATEB0')
          if(.not.is_first_read())then
             if(UseStrict)RETURN
             CYCLE
          end if

          call read_planet_var(NameCommand)

       case("#ROTATEHGR")
          if(.not.is_first_read())then
             if(UseStrict)RETURN
             CYCLE
          end if
          call read_var('dLongitudeHgr', dLongitudeHgrDeg)
          dLongitudeHgr = dLongitudeHgrDeg * cDegToRad

       case("#ROTATEHGI")
          if(.not.is_first_read())then
             if(UseStrict)RETURN
             CYCLE
          end if
          call read_var('dLongitudeHgi', dLongitudeHgiDeg)
          dLongitudeHgi = dLongitudeHgiDeg * cDegToRad

       case default
          if(is_proc0()) write(*,*) NameSub,' ERROR: Invalid command ',&
               trim(NameCommand),' at line',iLine,' in PARAM.in'
          if(UseStrict)then
             iErrorSwmf = 30
             RETURN
          end if
       end select
    end do
    ! end reading parameters

    !^CMP IF IE BEGIN
    !^CMP IF UA BEGIN
    ! Check if UA is used but IE is not
    if(use_comp(UA_) .and. .not. use_comp(IE_)) then
       write(*,*) NameSub//' SWMF_ERROR: '//&
            'UA is used without IE. This combination is not allowed!'
       call world_clean
       iErrorSwmf = 31
       RETURN
    end if
    ! Check if UA and IE are coupled
    if(use_comp(UA_) .and. .not. Couple_CC(IE_,UA_) % DoThis) then
       write(*,*) NameSub//' SWMF_ERROR: '//&
            'UA is used without IE-->UA coupling! Not allowed!'
       call world_clean
       iErrorSwmf = 32
       RETURN
    end if

    !^CMP END UA
    !^CMP END IE

    ! Switch off couplings for unused/switched off components
    do iComp1=1,MaxComp;
       if(use_comp(iComp1)) CYCLE
       Couple_CC(iComp1,:) % DoThis = .false.
       Couple_CC(:,iComp1) % DoThis = .false.
    end do

    ! Adjust nNext and tNext fields
    if(is_proc0()) then
       call check_freq('SWMF save restart', SaveRestart, DoTimeAccurate)
       call check_freq('SWMF check stop', CheckStop, DoTimeAccurate)
    end if
    call adjust_freq(SaveRestart, nStep  , tSimulation, DoTimeAccurate)
    call adjust_freq(CheckStop  , nStep+1, tSimulation+cTiny, DoTimeAccurate)

    do iComp1=1,MaxComp
       do iComp2=1,MaxComp

          if(is_proc0())call check_freq(&
               'SWMF couple '//NameComp_I(iComp1)//'->'//NameComp_I(iComp2), &
               Couple_CC(iComp1,iComp2), DoTimeAccurate)

          call adjust_freq( Couple_CC(iComp1,iComp2), &
               nStep+1,tSimulation+cTiny, DoTimeAccurate)
       end do
    end do

    ! Check and set some planet variables (e.g. DoUpdateB0)
    call check_planet_var(is_proc0(), DoTimeAccurate)

    ! Initialize axes
    call init_axes(TimeStart % Time)

    do lComp=1,n_comp()
       iComp = i_comp(lComp)
       if(use_comp(iComp))then
          ! Set component name so that i_line_command() works
          call read_init(NameComp_I(iComp),iSession,iLine)
          call set_param_comp(iComp,'CHECK')
       end if
    end do
    ! reset ModReadParam to read the CON parameters
    call read_init('  ',iSession,iLine)

    if(UseTiming)then
       call timing_version(IsOn,NameVersion,VersionNumber)
       if(.not.IsOn)then
          if(is_proc0()) &
               write(*,'(a)')NameSub//' WARNING: TIMING module is OFF'
          if(UseStrict)then
             iErrorSwmf = 33
             RETURN
          end if
          if(is_proc0()) write(*,*)NameSub//': setting UseTiming=.false.'
          UseTiming=.false.
       end if
    end if

    ! Initialize timing for CON
    if(is_proc0() .or. UseTimingAll)then
       call timing_active(UseTiming)
       call timing_comp_proc('  ',i_proc())
    end if
    ! Initialize timing for components
    do lComp=1,n_comp()
       iComp = i_comp(lComp)
       if(.not.use_comp(iComp))CYCLE
       if(is_proc0(iComp) .or. (UseTimingAll .and. is_proc(iComp)))then
          call timing_active(UseTiming)
          call timing_comp_proc(NameComp_I(iComp),i_proc())
       end if
    end do

    if(IsFirstRead)call timing_step(0)
    call timing_depth(nDepthTiming)
    call timing_report_style(TypeTimingReport)

    IsFirstRead = .false.

  contains
    !==========================================================================
    logical function is_first_read()

      if(.not.IsFirstRead)then
         if(is_proc0()) write(*,*) NameSub,' ERROR: command ',&
                  trim(NameCommand),&
                  ' can be used in first session only at line',&
                  iLine,' in PARAM.in'
         if(UseStrict) iErrorSwmf = 40
      end if
      is_first_read = IsFirstRead

    end function is_first_read
    !==========================================================================
  end subroutine read_inputs

  !BOP ========================================================================
  !IROUTINE: set_stdout - redirect standard output of components
  !INTERFACE:
  subroutine set_stdout
    !DESCRIPTION:
    ! If UseStdout is true tell the components to use STDOUT with
    ! a prefix string.
    ! If UseStdout is false, open an output file for each component and 
    ! each PE used by the component in the directory NameStdoutDir.
    ! Tell each component to write STDOUT into the file, which has 
    ! a unique unit number.
    ! Empty files will be deleted if the run ends in a normal fashion.
    !EOP
    character (len=*), parameter :: NameSub=NameMod//'::set_stdout'

    integer :: lComp, iComp, iUnitOut
    character (len=lNameFile) :: NameFile
    !-------------------------------------------------------------------------
    if(.not.UseStdout .and. is_proc0())call check_dir(NameStdoutDir)
    do lComp=1,n_comp()
       iComp=i_comp(lComp)

       if(.not.is_proc(iComp)) CYCLE

       if(.not.UseStdout)then
          ! Get the unit number from the registry
          call get_comp_info(iComp,iUnitOut=iUnitOut)
          ! Check if unit has been set
          if(iUnitOut==STDOUT_)then
             ! Construct the unique file name for each component and PE
             write(NameFile,'(a,i4.4,a)') trim(NameStdoutDir)// &
                  NameComp_I(iComp),i_proc(),NameStdoutExt
             ! Get a new unit number
             iUnitOut=io_unit_new()
             ! Open the file
             open(iUnitOut,FILE=NameFile,STATUS='unknown')
             ! Store the unit number into the registry
             call put_comp_info(iComp,iUnitOut=iUnitOut)
          end if
          ! Tell the component to write into a file
          call set_param_comp(iComp,"FILEOUT")
          ! Tell the timing utility to use the file
          call timing_iounit(iUnitOut)
       else
          ! Tell the component to write to STDOUT
          call set_param_comp(iComp,"STDOUT")
          ! Tell the timing utilty to use STDOUT
          call timing_iounit(STDOUT_)
       end if
    end do
  end subroutine set_stdout

  !BOP =======================================================================
  !IROUTINE: save_restart - save restart information for SWMF
  !INTERFACE:
  subroutine save_restart

    !DESCRIPTION:
    ! Save restart information for all components.
    ! Then save restart information for CON, such as code version,
    ! precision, name of planet, problem description,
    ! start date and time, simulation time, number of time steps,
    ! etc. The file is in the \#COMMAND parameters format, and
    ! it should simply be included into the input parameter file
    ! with the \#INCLUDE file command.
    !EOP

    character(len=*), parameter :: NameSub=NameMod//'::save_restart'
    
    integer :: lComp, iComp, iError
    !------------------------------------------------------------------------

    if(lVerbose>0 .and. is_proc0()) &
         write(*,*)NameSub,' is called at nStep,tSimulation=',&
         nStep,tSimulation

    do lComp = 1,n_comp()
       iComp = i_comp(lComp)
       call save_restart_comp(iComp,tSimulation)
    end do

    ! Ensure that all components have written restart state before 
    ! writing the CON restart file (RESTART.out)
    call MPI_barrier(i_comm(), iError)

    if(.not.is_proc0()) RETURN

    open(UNITTMP_,FILE=NameRestartFile)

    write(UNITTMP_,'(a)')'#DESCRIPTION'
    write(UNITTMP_,'(a)')trim(StringDescription)
    write(UNITTMP_,*)
    write(UNITTMP_,'(a)')'#PLANET'
    write(UNITTMP_,'(a25,a15)') NamePlanet,         '     NamePlanet'
    write(UNITTMP_,*)
    if(i_line_command("#IDEALAXES", iSessionIn=1) > 0)then
       write(UNITTMP_,'(a)')'#IDEALAXES'
       write(UNITTMP_,*)
    end if
    write(UNITTMP_,'(a)')'#STARTTIME'
    write(UNITTMP_,'(i8,a32)')TimeStart % iYear,         'iYear      '
    write(UNITTMP_,'(i8,a32)')TimeStart % iMonth,        'iMonth     '
    write(UNITTMP_,'(i8,a32)')TimeStart % iDay,          'iDay       '
    write(UNITTMP_,'(i8,a32)')TimeStart % iHour,         'iHour      '
    write(UNITTMP_,'(i8,a32)')TimeStart % iMinute,       'iMinute    '
    write(UNITTMP_,'(i8,a32)')TimeStart % iSecond,       'iSecond    '
    write(UNITTMP_,'(f15.12,a25)')TimeStart % FracSecond,'FracSecond '
    write(UNITTMP_,*)
    write(UNITTMP_,'(a)')'#NSTEP'
    write(UNITTMP_,'(i8,a32)')nStep,                     'nStep      '
    write(UNITTMP_,*)
    write(UNITTMP_,'(a)')'#TIMESIMULATION'
    write(UNITTMP_,'(es15.8,a25)')tSimulation,           'tSimulation'
    write(UNITTMP_,*)
    write(UNITTMP_,'(a)')'#VERSION'
    write(UNITTMP_,'(f5.2,a35)')VersionSwmf,             'VersionSwmf'
    write(UNITTMP_,*)
    write(UNITTMP_,'(a)')'#PRECISION'
    write(UNITTMP_,'(i1,a39)')nByteReal,                 'nByteReal  '
    write(UNITTMP_,*)
    if(dLongitudeHgr /= 0.0)then
       write(UNITTMP_,'(a)')'#ROTATEHGR'
       write(UNITTMP_,'(es15.8,a27)')dLongitudeHgrDeg,   'dLongitudeHgr'
       write(UNITTMP_,*)
    endif
    if(dLongitudeHgi /= 0.0)then
       write(UNITTMP_,'(a)')'#ROTATEHGI'
       write(UNITTMP_,'(es15.8,a27)')dLongitudeHgiDeg,   'dLongitudeHgi'
       write(UNITTMP_,*)
    endif

    close(UNITTMP_)

  end subroutine save_restart

end module CON_io
