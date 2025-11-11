!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!
!
!
!
!
!
! Contains general subroutines which mimic the behaviour of the
! physics components. Based on the real wrapper and the real physical
! components. Selectt INTVERION = Stubs in Makefile.def and you can select
! the Empty version for all the components.
!
! This is useful for testing and profiling the control component.
! For example one can set a fraction of the realistic CPU times per
! time step and measure performance with various layouts, coupling
! strategies and session models. Only one command is recognized by
! the stub components:
! \begin{verbatim}
! #TIMESTEP
! 1.0        DtRun (the typical time step of the component)
! 0.12       DtCpu (the CPU time needed for 1 time step)
! \end{verbatim}
! Both run time and CPU time are in seconds.

module CON_wrapper
  !
  !
  use CON_world
  use CON_comp_param
  use CON_comp_info
  use ModMpi, ONLY: MPI_WTIME
  use ModIoUnit, ONLY: UNITTMP_
  use CON_time, ONLY: nStep, tSimulation
  use ModReadParam
  use ModKind
  use ModUtilities, ONLY: sleep

  implicit none

  private   ! except

  public :: set_param_comp       ! set parameters for registered component
  interface set_param_comp
     module procedure set_param_name, set_param_id
  end interface

  public :: get_version_comp     ! get version from (unregistered) component
  interface get_version_comp
     module procedure get_version_name, get_version_id
  end interface

  public :: init_session_comp    ! initialize component for session
  interface init_session_comp
     module procedure init_session_comp_id
  end interface

  public :: run_comp             ! run component (advance/solve)
  interface run_comp
     module procedure run_comp_id
  end interface

  public :: save_restart_comp    ! save restart info for component
  interface save_restart_comp
     module procedure save_restart_comp_id
  end interface

  public :: finalize_comp       ! finalize component
  interface finalize_comp
     module procedure finalize_comp_id
  end interface

  logical, public  :: IsNewInputIe = .false. ! Store if IE got new input

  real,    public  :: TimeNewInputSp=-2.0    ! Last time when SP got input

  real,    public  :: DtCpu_C(MaxComp)=1.0   ! CPU time needed for a step

  ! revision history:
  ! 30Jul03 - Gabor Toth <gtoth@umich.edu> - initial prototype/prolog/code
  !                                          based on the actual wrapper
  ! 20Jul04 - Gabor Toth added exceptions for SP coupling (experimental)

  character(len=*),parameter :: NameMod='CON_wrapper_stub'

  real :: TimeSp = -1.0                      ! The actual time SP has reached

  real,    dimension(MaxComp) :: DtRun_C=1.0
  real(Real8_)                :: Cpu0
  integer, dimension(MaxComp) :: nRestart_C=0

  character (len=100) :: NameFile, NameCommand

contains
  !============================================================================

  !
  subroutine set_param_name(NameComp,TypeAction)

    !
    character(len=*), intent(in)           :: NameComp, TypeAction

    ! This routine calls set\_param\_id with the component ID defined by the
    ! name of the component.

    ! revision history:
    ! 29Aug03 - G. Toth <gtoth@umich.edu> - initial prototype

    character(len=*), parameter:: NameSub = 'set_param_name'
    !--------------------------------------------------------------------------
    call set_param_id(i_comp(NameComp),TypeAction)

  end subroutine set_param_name
  !============================================================================

  !
  subroutine set_param_id(iComp,TypeAction)

    use CON_comp_info

    !
    integer, intent(in)          :: iComp
    character(len=*), intent(in) :: TypeAction

    ! This routine calls $**$\_set\_param for component $**$
    ! defined with the component ID.

    ! revision history:
    ! 29Aug03 - G. Toth <gtoth@umich.edu> - initial prototype/prolog/code

    type(CompInfoType) :: CompInfo
    integer :: iUnitOut
    character(len=*), parameter:: NameSub = 'set_param_id'
    !--------------------------------------------------------------------------
    if(.not.use_comp(iComp)) call CON_stop(NameSub//' '//TypeAction// &
         ' SWMF_ERROR component '//NameComp_I(iComp)//' is not used')

    ! Check action and return if possible
    select case(TypeAction)
    case('VERSION','GRID')
       ! These actions should be done on all PE-s
    case('MPI','READ','CHECK','STDOUT','FILEOUT')
       if(.not.is_proc(iComp)) RETURN
    case default
       call CON_stop(NameSub//' SWMF_ERROR: invalid TypeAction='//TypeAction)
    end select

    ! Check out component information
    call get_comp_info(iComp,CompInfo=CompInfo,iUnitOut=iUnitOut)

    select case(TypeAction)
    case('VERSION')
       call put(CompInfo,&
            NameVersion=NameComp_I(iComp)//' stub',&
            Version=1.0,&
            Use=.true.)
    case('GRID','MPI','CHECK','STDOUT','FILEOUT')
       write(iUnitOut,*)NameComp_I(iComp),': set_param_comp iProc, Action=',&
            i_proc(),' ',TypeAction
    case('READ')
       do
          if(.not.read_line()) EXIT
          if(.not.read_command(NameCommand)) CYCLE
          select case(NameCommand)
          case('#TIMESTEP')
             call read_var('DtRun',DtRun_C(iComp))
             call read_var('DtCpu',DtCpu_C(iComp))
          case default
             write(iUnitOut,*)NameComp_I(iComp),&
                  ': stub cannot interpret command NameCommand=',NameCommand
          end select
       end do
    case default
       call CON_stop(NameSub//' SWMF_ERROR: invalid TypeAction='//TypeAction)
    end select

    ! Store returned component information
    select case(TypeAction)
    case('VERSION')
       call put_comp_info(iComp,Use=CompInfo%Use,Version=CompInfo%Version,&
            NameVersion=CompInfo%NameVersion)
    end select

  end subroutine set_param_id
  !============================================================================

  !
  ! Obtain ON/OFF status, version name and version number
  ! for component specified by its name!

  subroutine get_version_name(NameComp,IsOn,NameVersion,Version)

    character (len=*),            intent(in)  :: NameComp
    logical,                      intent(out) :: IsOn
    character (len=lNameVersion), intent(out) :: NameVersion
    real,                         intent(out) :: Version

    ! revision history:
    ! 30Aug03 - G. Toth <gtoth@umich.edu> - initial prototype/prolog/code

    character(len=*), parameter:: NameSub = 'get_version_name'
    !--------------------------------------------------------------------------
    if(.not.is_valid_comp_name(NameComp)) call CON_stop(NameSub// &
         ' SWMF_ERROR '//NameComp//' is not a valid component name')

    call get_version_id(i_comp(NameComp),IsOn,NameVersion,Version)

  end subroutine get_version_name
  !============================================================================

  !
  ! Obtain ON/OFF status, version name and version number
  ! for component specified by its ID!

  subroutine get_version_id(iComp,IsOn,NameVersion,Version)

    integer,                      intent(in)  :: iComp
    logical,                      intent(out) :: IsOn
    character (len=lNameVersion), intent(out) :: NameVersion
    real,                         intent(out) :: Version

    ! revision history:
    ! 21Jul03 - G. Toth <gtoth@umich.edu> - use temporary CompInfo
    ! 15Jul03 - G. Toth <gtoth@umich.edu> - initial prototype/prolog/code

    type(CompInfoType) :: CompInfo

    character(len=*), parameter:: NameSub = 'get_version_id'
    !--------------------------------------------------------------------------
    IsOn=.true.
    NameVersion=NameComp_I(iComp)//' stub'
    Version=1.0

  end subroutine get_version_id
  !============================================================================

  subroutine init_session_comp_id(iComp, iSession, TimeSimulation)

    integer,  intent(in) :: iComp            ! component ID
    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    ! Initialize component for the session. Session numbers start from 1.

    ! revision history:
    ! 18Aug03 - G. Toth <gtoth@umich.edu> O. Volberg <volov@umich.edu>
    !         - initial prototype/prolog/code
    integer :: iUnitOut

    character(len=*), parameter:: NameSub = 'init_session_comp_id'
    !--------------------------------------------------------------------------
    call check_i_comp(iComp,NameSub)

    if(.not.use_comp(iComp) .or. .not.is_proc(iComp)) RETURN

    call get_comp_info(iComp,iUnitOut=iUnitOut)

    ! Set initial SP time
    if(iComp == SP_)then
       write(iUnitOut,*)'SP: ',NameSub,' old TimeSp=',TimeSp
       if(TimeSp < 0.0) TimeSp = TimeSimulation
       write(iUnitOut,*)'SP: ',NameSub,' new TimeSp=',TimeSp
    end if

    write(iUnitOut,*)NameComp_I(iComp),': ',NameSub,&
         ' iProc,iSession,tSim=',&
         i_proc(),iSession,TimeSimulation

  end subroutine init_session_comp_id
  !============================================================================

  !
  subroutine finalize_comp_id(iComp, TimeSimulation)

    integer,  intent(in) :: iComp            ! component ID
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    ! Finalize the component at the end of the run, save plots, restart info,
    ! report errors, etc.

    ! revision history:
    ! 18Aug03 - G. Toth <gtoth@umich.edu> O. Volberg <volov@umich.edu>
    !         - initial prototype/prolog/code
    integer :: iUnitOut
    character(len=*), parameter:: NameSub = 'finalize_comp_id'
    !--------------------------------------------------------------------------

    call check_i_comp(iComp,NameSub)

    if(.not.use_comp(iComp) .or. .not.is_proc(iComp)) RETURN

    call get_comp_info(iComp,iUnitOut=iUnitOut)

    write(iUnitOut,*)NameComp_I(iComp),': ',NameSub,' iProc,TimeSimulation=',&
         i_proc(),TimeSimulation

  end subroutine finalize_comp_id
  !============================================================================

  subroutine save_restart_comp_id(iComp, TimeSimulation)

    integer,  intent(in) :: iComp            ! component ID
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    ! Save restart information for the component.

    ! revision history:
    ! 18Aug03 - G. Toth <gtoth@umich.edu> O. Volberg <volov@umich.edu>
    !         - initial prototype/prolog/code
    integer :: iUnitOut
    character(len=*), parameter:: NameSub = 'save_restart_comp_id'
    !--------------------------------------------------------------------------

    call check_i_comp(iComp,NameSub)

    if(.not.use_comp(iComp) .or. .not.is_proc(iComp)) RETURN

    call get_comp_info(iComp,iUnitOut=iUnitOut)

    write(iUnitOut,*)NameComp_I(iComp),': ',NameSub,' iProc,TimeSimulation=',&
         i_proc(),TimeSimulation

    if(is_proc0(iComp))then
       nRestart_C(iComp)=nRestart_C(iComp)+1
       write(NameFile,'(a,i2.2)') &
            NameComp_I(iComp)//'restart_',nRestart_C(iComp)
       open(UNITTMP_,file=NameFile)
       write(UNITTMP_,*)nStep,'   nStep'
       write(UNITTMP_,*)tSimulation,'   tSimulation'
       close(UNITTMP_)
    end if

  end subroutine save_restart_comp_id
  !============================================================================

  subroutine run_comp_id(iComp, TimeSimulation, TimeSimulationLimit)

    integer,  intent(in) :: iComp               ! component ID
    real,     intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded

    real,     intent(inout) :: TimeSimulation   ! current time of component

    ! Do one time step or solve for the component.
    ! The component should update TimeSimulation.
    ! The updated TimeSimulation should not exceed TimeSimulationLimit.

    ! revision history:
    ! 30Aug03 - G. Toth <gtoth@umich.edu> initial prototype/prolog/code

    integer :: iUnitOut, iError

    character(len=*), parameter:: NameSub = 'run_comp_id'
    !--------------------------------------------------------------------------
    call check_i_comp(iComp,NameSub)

    if(.not.use_comp(iComp) .or. .not.is_proc(iComp)) RETURN

    call get_comp_info(iComp,iUnitOut=iUnitOut)

    write(iUnitOut,*)NameComp_I(iComp),': ',NameSub,&
         ' iProc,nStep,tSim,Limit =',&
         i_proc(),nStep,TimeSimulation, TimeSimulationLimit

    TimeSimulation = min(TimeSimulationLimit, TimeSimulation + DtRun_C(iComp))

    ! SP can only solve up to the last coupling time TimeNewInputSp
    if(iComp==SP_)then
       ! write(iUnitOut,*)'SP: ',NameSub,' TimeNewInputSp=',TimeNewInputSp
       ! write(iUnitOut,*)'SP: ',NameSub,' TimeSp        =',TimeSp
       if(TimeNewInputSp > TimeSp)then
          ! Advance SP time by one step
          TimeSp = min(TimeNewInputSp, TimeSp + DtRun_C(iComp))
          write(iUnitOut,*)'SP: ',NameSub,' advancing to TimeSp=',TimeSp
          ! Tell the framework if we are done by setting a fake TimeSimulation
          TimeSimulation = TimeSimulationLimit + (TimeSp - TimeNewInputSp)
       else
          write(iUnitOut,*)'SP: ',NameSub,' cannot solve yet iProc=',i_proc()
          ! Pretend that we are done
          TimeSimulation = TimeSimulationLimit
          RETURN
       end if
    end if

    ! IE only solves if there is new info
    if(iComp==IE_)then
       if(.not.IsNewInputIe)then
          write(iUnitOut,*)'IE: ',NameSub,' no need to solve iProc=',i_proc()
          RETURN
       else
          write(iUnitOut,*)'IE: ',NameSub,' solve iProc=',i_proc()
          IsNewInputIe = .false.
       end if
    end if

    ! Synchronize PE-s used by the component
    call timing_start(NameComp_I(iComp)//'_barrier')
    call MPI_barrier(i_comm(iComp),iError)
    call timing_stop(NameComp_I(iComp)//'_barrier')

    ! Work on this time step
    call timing_start(NameComp_I(iComp)//'_run')
    call sleep(DtCpu_C(iComp))
    call timing_stop(NameComp_I(iComp)//'_run')

  end subroutine run_comp_id
  !============================================================================

end module CON_wrapper
!==============================================================================
