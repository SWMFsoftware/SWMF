!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module CON_wrapper

  ! General subroutines which call the specific subroutines
  ! in the component wrappers. The components are identified by their
  ! name or the component ID.

  use CON_world
  use CON_comp_param
  use CON_comp_info
  use ModUtilities, ONLY: CON_set_do_test, CON_stop
  use omp_lib

  use CZ_wrapper, ONLY: CZ_set_param, CZ_init_session, &       !^CMP IF CZ
       CZ_run, CZ_finalize, CZ_save_restart                    !^CMP IF CZ
  use EE_wrapper, ONLY: EE_set_param, EE_init_session, &       !^CMP IF EE
       EE_run, EE_finalize, EE_save_restart                    !^CMP IF EE
  use GM_wrapper, ONLY: GM_set_param, GM_init_session, &       !^CMP IF GM
       GM_run, GM_finalize, GM_save_restart                    !^CMP IF GM
  use IE_wrapper, ONLY: IE_set_param, IE_init_session, &       !^CMP IF IE
       IE_run, IE_finalize, IE_save_restart                    !^CMP IF IE
  use IH_wrapper, ONLY: IH_set_param, IH_init_session, &       !^CMP IF IH
       IH_run, IH_finalize, IH_save_restart                    !^CMP IF IH
  use IM_wrapper, ONLY: IM_set_param, IM_init_session, &       !^CMP IF IM
       IM_run, IM_finalize, IM_save_restart                    !^CMP IF IM
  use OH_wrapper, ONLY: OH_set_param, OH_init_session, &       !^CMP IF OH
       OH_run, OH_finalize, OH_save_restart                    !^CMP IF OH
  use PS_wrapper, ONLY: PS_set_param, PS_init_session, &       !^CMP IF PS
       PS_run, PS_finalize, PS_save_restart                    !^CMP IF PS
  use PC_wrapper, ONLY: PC_set_param, PC_init_session, &       !^CMP IF PC
       PC_run, PC_finalize, PC_save_restart                    !^CMP IF PC
  use PT_wrapper, ONLY: PT_set_param, PT_init_session, &       !^CMP IF PT
       PT_run, PT_finalize, PT_save_restart                    !^CMP IF PT
  use PW_wrapper, ONLY: PW_set_param, PW_init_session, &       !^CMP IF PW
       PW_run, PW_finalize, PW_save_restart                    !^CMP IF PW
  use RB_wrapper, ONLY: RB_set_param, RB_init_session, &       !^CMP IF RB
       RB_run, RB_finalize, RB_save_restart                    !^CMP IF RB
  use SC_wrapper, ONLY: SC_set_param, SC_init_session, &       !^CMP IF SC
       SC_run, SC_finalize, SC_save_restart                    !^CMP IF SC
  use SP_wrapper, ONLY: SP_set_param, SP_init_session, &       !^CMP IF SP
       SP_run, SP_finalize, SP_save_restart                    !^CMP IF SP
  use UA_wrapper, ONLY: UA_set_param, UA_init_session, &       !^CMP IF UA
       UA_run, UA_finalize, UA_save_restart                    !^CMP IF UA

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

  ! revision history:
  ! 09Jul03 - Gabor Toth <gtoth@umich.edu> - initial prototype/prolog/code
  ! 20Aug03 - O. Volberg and G. Toth - added new methods to initialize,
  !                                    run session, finalize and save restart.

  integer :: nThread, MaxThread
  character(len=*),parameter :: NameMod='CON_wrapper'

contains
  !============================================================================
  subroutine set_param_name(NameComp,TypeAction)

    !
    character(len=*), intent(in)           :: NameComp, TypeAction

    ! This routine calls set\_param\_id with the component ID defined by the
    ! name of the component.

    ! revision history:
    ! 21Jul03 - G. Toth <gtoth@umich.edu> - Pass CompInfo for safety
    ! 19Jul03 - G. Toth <gtoth@umich.edu> - simplified and generalized version
    ! 09Jul03 - G. Toth <gtoth@umich.edu> - initial prototype/prolog/code

    character(len=*), parameter:: NameSub = 'set_param_name'
    !--------------------------------------------------------------------------
    call set_param_id(i_comp(NameComp),TypeAction)

  end subroutine set_param_name
  !============================================================================
  subroutine set_param_id(iComp,TypeAction)

    integer, intent(in)          :: iComp
    character(len=*), intent(in) :: TypeAction

    ! This routine calls $**$\_set\_param for component $**$
    ! defined with the component ID.

    ! revision history:
    ! 19Jul03 - G. Toth <gtoth@umich.edu> - simplified and generalized version
    ! 09Jul03 - G. Toth <gtoth@umich.edu> - initial prototype/prolog/code
    ! 03Jun08 - R Oran  <oran@umich.edu>  - add component OH

    logical :: DoTest, DoTestMe
    type(CompInfoType) :: CompInfo

    character(len=*), parameter:: NameSub = 'set_param_id'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTestMe)write(*,*) NameSub,' starting for Comp, TypeAction=',&
         NameComp_I(iComp),' ',TypeAction

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
    call get_comp_info(iComp,CompInfo=CompInfo)

    !$ call omp_set_num_threads(CompInfo%nThread)

    select case(iComp)
    case(CZ_)                                     !^CMP IF CZ
       call CZ_set_param(CompInfo,TypeAction)     !^CMP IF CZ
    case(EE_)                                     !^CMP IF EE
       call EE_set_param(CompInfo,TypeAction)     !^CMP IF EE
    case(GM_)                                     !^CMP IF GM
       call GM_set_param(CompInfo,TypeAction)     !^CMP IF GM
    case(IE_)                                     !^CMP IF IE
       call IE_set_param(CompInfo,TypeAction)     !^CMP IF IE
    case(IH_)                                     !^CMP IF IH
       call IH_set_param(CompInfo,TypeAction)     !^CMP IF IH
    case(IM_)                                     !^CMP IF IM
       call IM_set_param(CompInfo,TypeAction)     !^CMP IF IM
    case(OH_)                                     !^CMP IF OH
       call OH_set_param(CompInfo,TypeAction)     !^CMP IF OH
    case(PC_)                                     !^CMP IF PC
       call PC_set_param(CompInfo,TypeAction)     !^CMP IF PC
    case(PS_)                                     !^CMP IF PS
       call PS_set_param(CompInfo,TypeAction)     !^CMP IF PS
    case(PT_)                                     !^CMP IF PT
       call PT_set_param(CompInfo,TypeAction)     !^CMP IF PT
    case(PW_)                                     !^CMP IF PW
       call PW_set_param(CompInfo,TypeAction)     !^CMP IF PW
    case(RB_)                                     !^CMP IF RB
       call RB_set_param(CompInfo,TypeAction)     !^CMP IF RB
    case(SC_)                                     !^CMP IF SC
       call SC_set_param(CompInfo,TypeAction)     !^CMP IF SC
    case(SP_)                                     !^CMP IF SP
       call SP_set_param(CompInfo,TypeAction)     !^CMP IF SP
    case(UA_)                                     !^CMP IF UA
       call UA_set_param(CompInfo,TypeAction)     !^CMP IF UA
    case default
       call CON_stop(NameSub//' '//TypeAction//&
            ' SWMF_ERROR: not implemented for component'//NameComp_I(iComp))
    end select

    ! Reset the number of threads
    !$ MaxThread = omp_get_max_threads()
    !$ call omp_set_num_threads(MaxThread)

    ! Store returned component information
    select case(TypeAction)
    case('VERSION')
       call put_comp_info(iComp,Use=CompInfo%Use,Version=CompInfo%Version,&
            NameVersion=CompInfo%NameVersion)
    end select

    if(DoTestMe)write(*,*) NameSub,' finished TypeAction=',TypeAction, &
         ' for ',NameComp_I(iComp)

  end subroutine set_param_id
  !============================================================================
  subroutine get_version_name(NameComp, IsOn, NameVersion)

    ! Obtain ON/OFF status, version name and version number
    ! for component specified by its name!

    character (len=*),            intent(in)  :: NameComp
    logical,                      intent(out) :: IsOn
    character (len=lNameVersion), intent(out) :: NameVersion

    ! revision history:
    ! 15Jul03 - G. Toth <gtoth@umich.edu> - initial prototype/prolog/code

    character(len=*), parameter:: NameSub = 'get_version_name'
    !--------------------------------------------------------------------------
    if(.not.is_valid_comp_name(NameComp)) call CON_stop(NameSub// &
         ' SWMF_ERROR '//NameComp//' is not a valid component name')

    call get_version_id(i_comp(NameComp), IsOn, NameVersion)

  end subroutine get_version_name
  !============================================================================
  subroutine get_version_id(iComp, IsOn, NameVersion)

    ! Obtain ON/OFF status, version name and version number
    ! for component specified by its ID!

    use CON_comp_info

    integer,                      intent(in)  :: iComp
    logical,                      intent(out) :: IsOn
    character (len=lNameVersion), intent(out) :: NameVersion

    ! revision history:
    ! 21Jul03 - G. Toth <gtoth@umich.edu> - use temporary CompInfo
    ! 15Jul03 - G. Toth <gtoth@umich.edu> - initial prototype/prolog/code
    ! 03Jun08 - R. Oran <oran@umich.edu>  - added component OH

    logical :: DoTest, DoTestMe
    type(CompInfoType) :: CompInfo

    character(len=*), parameter:: NameSub = 'get_version_id'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTestMe)write(*,*) NameSub,' starting for ',NameComp_I(iComp)

    call check_i_comp(iComp,NameSub)

    select case(iComp)
    case(EE_)                                    !^CMP IF EE
       call EE_set_param(CompInfo,'VERSION')     !^CMP IF EE
    case(GM_)                                    !^CMP IF GM
       call GM_set_param(CompInfo,'VERSION')     !^CMP IF GM
    case(IE_)                                    !^CMP IF IE
       call IE_set_param(CompInfo,'VERSION')     !^CMP IF IE
    case(IH_)                                    !^CMP IF IH
       call IH_set_param(CompInfo,'VERSION')     !^CMP IF IH
    case(OH_)                                    !^CMP IF OH
       call OH_set_param(CompInfo,'VERSION')     !^CMP IF OH
    case(IM_)                                    !^CMP IF IM
       call IM_set_param(CompInfo,'VERSION')     !^CMP IF IM
    case(PC_)                                    !^CMP IF PC
       call PC_set_param(CompInfo,'VERSION')     !^CMP IF PC
    case(PS_)                                    !^CMP IF PS
       call PS_set_param(CompInfo,'VERSION')     !^CMP IF PS
    case(PT_)                                    !^CMP IF PT
       call PT_set_param(CompInfo,'VERSION')     !^CMP IF PT
    case(PW_)                                    !^CMP IF PW
       call PW_set_param(CompInfo,'VERSION')     !^CMP IF PW
    case(RB_)                                    !^CMP IF RB
       call RB_set_param(CompInfo,'VERSION')     !^CMP IF RB
    case(SC_)                                    !^CMP IF SC
       call SC_set_param(CompInfo,'VERSION')     !^CMP IF SC
    case(SP_)                                    !^CMP IF SP
       call SP_set_param(CompInfo,'VERSION')     !^CMP IF SP
    case(UA_)                                    !^CMP IF UA
       call UA_set_param(CompInfo,'VERSION')     !^CMP IF UA
    case(CZ_)                                    !^CMP IF CZ
       call CZ_set_param(CompInfo,'VERSION')     !^CMP IF CZ
    case default
       call put(CompInfo, Use=.false., NameVersion='not implemented')
    end select

    call get(CompInfo, Use=IsOn, NameVersion=NameVersion)

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
    ! 03Jun   - R. Oran <oran@umich.edu> - added component OH
    logical :: DoTest, DoTestMe

    character(len=*), parameter:: NameSub = 'init_session_comp_id'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    call check_i_comp(iComp,NameSub)

    if(DoTestMe)write(*,*) NameSub,' starting for ',NameComp_I(iComp),&
         ': iSession, Time=',iSession, TimeSimulation

    if(.not.use_comp(iComp) .or. .not.is_proc(iComp)) RETURN

    ! Set the number of threads
    !$ call get_comp_info(iComp, nThread=nThread)
    !$ call omp_set_num_threads(nThread)

    select case(iComp)
    case(EE_)                                            !^CMP IF EE
       call EE_init_session(iSession,TimeSimulation)     !^CMP IF EE
    case(GM_)                                            !^CMP IF GM
       call GM_init_session(iSession,TimeSimulation)     !^CMP IF GM
    case(IE_)                                            !^CMP IF IE
       call IE_init_session(iSession,TimeSimulation)     !^CMP IF IE
    case(IH_)                                            !^CMP IF IH
       call IH_init_session(iSession,TimeSimulation)     !^CMP IF IH
    case(OH_)                                            !^CMP IF OH
       call OH_init_session(iSession,TimeSimulation)     !^CMP IF OH
    case(IM_)                                            !^CMP IF IM
       call IM_init_session(iSession,TimeSimulation)     !^CMP IF IM
    case(PC_)                                            !^CMP IF PC
       call PC_init_session(iSession,TimeSimulation)     !^CMP IF PC
    case(PS_)                                            !^CMP IF PS
       call PS_init_session(iSession,TimeSimulation)     !^CMP IF PS
    case(PT_)                                            !^CMP IF PT
       call PT_init_session(iSession,TimeSimulation)     !^CMP IF PT
    case(PW_)                                            !^CMP IF PW
       call PW_init_session(iSession,TimeSimulation)     !^CMP IF PW
    case(RB_)                                            !^CMP IF RB
       call RB_init_session(iSession,TimeSimulation)     !^CMP IF RB
    case(SC_)                                            !^CMP IF SC
       call SC_init_session(iSession,TimeSimulation)     !^CMP IF SC
    case(SP_)                                            !^CMP IF SP
       call SP_init_session(iSession,TimeSimulation)     !^CMP IF SP
    case(UA_)                                            !^CMP IF UA
       call UA_init_session(iSession,TimeSimulation)     !^CMP IF UA
    case(CZ_)                                            !^CMP IF CZ
       call CZ_init_session(iSession,TimeSimulation)     !^CMP IF CZ
    case default
       call CON_stop(NameSub//' SWMF_ERROR incorrect iComp value')
    end select

    ! Reset the number of threads
    !$ MaxThread = omp_get_max_threads()
    !$ call omp_set_num_threads(MaxThread)

    if(DoTestMe)write(*,*) NameSub,' finished for ',NameComp_I(iComp)

  end subroutine init_session_comp_id
  !============================================================================
  subroutine finalize_comp_id(iComp, TimeSimulation)

    integer,  intent(in) :: iComp            ! component ID
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    ! Finalize the component at the end of the run, save plots, restart info,
    ! report errors, etc.

    ! revision history:
    ! 18Aug03 - G. Toth <gtoth@umich.edu> O. Volberg <volov@umich.edu>
    !         - initial prototype/prolog/code
    ! 03Jun08 - R. Oran <oran@umich.edu> - added component OH
    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'finalize_comp_id'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    call check_i_comp(iComp,NameSub)

    if(DoTestMe)write(*,*) NameSub,' starting for ',NameComp_I(iComp),&
         ': Time=', TimeSimulation

    if(.not.use_comp(iComp) .or. .not.is_proc(iComp)) RETURN

    ! Set the number of threads
    !$ call get_comp_info(iComp, nThread=nThread)
    !$ call omp_set_num_threads(nThread)

    select case(iComp)
    case(EE_)                               !^CMP IF EE
       call EE_finalize(TimeSimulation)     !^CMP IF EE
    case(GM_)                               !^CMP IF GM
       call GM_finalize(TimeSimulation)     !^CMP IF GM
    case(IE_)                               !^CMP IF IE
       call IE_finalize(TimeSimulation)     !^CMP IF IE
    case(IH_)                               !^CMP IF IH
       call IH_finalize(TimeSimulation)     !^CMP IF IH
    case(OH_)                               !^CMP IF OH
       call OH_finalize(TimeSimulation)     !^CMP IF OH
    case(IM_)                               !^CMP IF IM
       call IM_finalize(TimeSimulation)     !^CMP IF IM
    case(PC_)                               !^CMP IF PC
       call PC_finalize(TimeSimulation)     !^CMP IF PC
    case(PS_)                               !^CMP IF PS
       call PS_finalize(TimeSimulation)     !^CMP IF PS
    case(PT_)                               !^CMP IF PT
       call PT_finalize(TimeSimulation)     !^CMP IF PT
    case(PW_)                               !^CMP IF PW
       call PW_finalize(TimeSimulation)     !^CMP IF PW
    case(RB_)                               !^CMP IF RB
       call RB_finalize(TimeSimulation)     !^CMP IF RB
    case(SC_)                               !^CMP IF SC
       call SC_finalize(TimeSimulation)     !^CMP IF SC
    case(SP_)                               !^CMP IF SP
       call SP_finalize(TimeSimulation)     !^CMP IF SP
    case(UA_)                               !^CMP IF UA
       call UA_finalize(TimeSimulation)     !^CMP IF UA
    case(CZ_)                               !^CMP IF CZ
       call CZ_finalize(TimeSimulation)     !^CMP IF CZ
    case default
       call CON_stop(NameSub//' SWMF_ERROR incorrect iComp value')
    end select

    ! Reset the number of threads
    !$ MaxThread = omp_get_max_threads()
    !$ call omp_set_num_threads(MaxThread)

    if(DoTestMe)write(*,*) NameSub,' finished for ',NameComp_I(iComp)

  end subroutine finalize_comp_id
  !============================================================================
  subroutine save_restart_comp_id(iComp, TimeSimulation)

    integer,  intent(in) :: iComp            ! component ID
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    ! Save restart information for the component.

    ! revision history:
    ! 18Aug03 - G. Toth <gtoth@umich.edu> O. Volberg <volov@umich.edu>
    !         - initial prototype/prolog/code
    ! 03Jun08 - R. Oran <oran@umich.edu> - added component OH
    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'save_restart_comp_id'
    !--------------------------------------------------------------------------

    call check_i_comp(iComp,NameSub)

    if(.not.use_comp(iComp) .or. .not.is_proc(iComp)) RETURN

    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTestMe)write(*,*) NameSub,' starting for ',NameComp_I(iComp),&
         ': Time=', TimeSimulation

    ! Set the number of threads
    !$ call get_comp_info(iComp, nThread=nThread)
    !$ call omp_set_num_threads(nThread)

    select case(iComp)
    case(EE_)                                   !^CMP IF EE
       call EE_save_restart(TimeSimulation)     !^CMP IF EE
    case(GM_)                                   !^CMP IF GM
       call GM_save_restart(TimeSimulation)     !^CMP IF GM
    case(IE_)                                   !^CMP IF IE
       call IE_save_restart(TimeSimulation)     !^CMP IF IE
    case(IH_)                                   !^CMP IF IH
       call IH_save_restart(TimeSimulation)     !^CMP IF IH
    case(OH_)                                   !^CMP IF OH
       call OH_save_restart(TimeSimulation)     !^CMP IF OH
    case(IM_)                                   !^CMP IF IM
       call IM_save_restart(TimeSimulation)     !^CMP IF IM
    case(PC_)                                   !^CMP IF PC
       call PC_save_restart(TImeSimulation)     !^CMP IF PC
    case(PS_)                                   !^CMP IF PS
       call PS_save_restart(TimeSimulation)     !^CMP IF PS
    case(PT_)                                   !^CMP IF PT
       call PT_save_restart(TimeSimulation)     !^CMP IF PT
    case(PW_)                                   !^CMP IF PW
       call PW_save_restart(TimeSimulation)     !^CMP IF PW
    case(RB_)                                   !^CMP IF RB
       call RB_save_restart(TimeSimulation)     !^CMP IF RB
    case(SC_)                                   !^CMP IF SC
       call SC_save_restart(TimeSimulation)     !^CMP IF SC
    case(SP_)                                   !^CMP IF SP
       call SP_save_restart(TimeSimulation)     !^CMP IF SP
    case(UA_)                                   !^CMP IF UA
       call UA_save_restart(TimeSimulation)     !^CMP IF UA
    case(CZ_)                                   !^CMP IF CZ
       call CZ_save_restart(TimeSimulation)     !^CMP IF CZ
    case default
       call CON_stop(NameSub//' SWMF_ERROR incorrect iComp value')
    end select

    ! Reset the number of threads
    !$ MaxThread = omp_get_max_threads()
    !$ call omp_set_num_threads(MaxThread)

    if(DoTestMe)write(*,*) NameSub,' finished for ',NameComp_I(iComp)

  end subroutine save_restart_comp_id
  !============================================================================
  subroutine run_comp_id(iComp, TimeSimulation, TimeSimulationLimit)

    integer,  intent(in) :: iComp               ! component ID
    real,     intent(in) :: TimeSimulationLimit ! simulation time limit
    real,     intent(inout) :: TimeSimulation   ! current time of component

    ! Do one time step or solve for the component.
    ! The component should update TimeSimulation.
    ! The updated TimeSimulation should not exceed TimeSimulationLimit.

    ! revision history:
    ! 18Aug03 - G. Toth <gtoth@umich.edu> O. Volberg <volov@umich.edu>
    !         - initial prototype/prolog/code
    ! 03Jun08 - R. Oran <oran@umich.edu - added component OH
    logical :: DoTest, DoTestMe

    character(len=*), parameter:: NameSub = 'run_comp_id'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    call check_i_comp(iComp,NameSub)

    if(.not.use_comp(iComp) .or. .not.is_proc(iComp)) RETURN

    if(DoTestMe)write(*,*) NameSub,' starting for ',NameComp_I(iComp),&
         ': Time, Limit=', TimeSimulation, TimeSimulationLimit

    ! Set the number of threads
    !$ call get_comp_info(iComp, nThread=nThread)
    !$ call omp_set_num_threads(nThread)

    call timing_start(NameComp_I(iComp)//'_run')
    select case(iComp)
    case(EE_)                                               !^CMP IF EE
       call EE_run(TimeSimulation, TimeSimulationLimit)     !^CMP IF EE
    case(GM_)                                               !^CMP IF GM
       call GM_run(TimeSimulation, TimeSimulationLimit)     !^CMP IF GM
    case(IE_)                                               !^CMP IF IE
       call IE_run(TimeSimulation, TimeSimulationLimit)     !^CMP IF IE
    case(IH_)                                               !^CMP IF IH
       call IH_run(TimeSimulation, TimeSimulationLimit)     !^CMP IF IH
    case(OH_)                                               !^CMP IF OH
       call OH_run(TimeSimulation, TimeSimulationLimit)     !^CMP IF OH
    case(IM_)                                               !^CMP IF IM
       call IM_run(TimeSimulation, TimeSimulationLimit)     !^CMP IF IM
    case(PC_)                                               !^CMP IF PC
       call PC_run(TimeSimulation, TimeSimulationLimit)     !^CMP IF PC
    case(PS_)                                               !^CMP IF PS
       call PS_run(TimeSimulation, TimeSimulationLimit)     !^CMP IF PS
    case(PT_)                                               !^CMP IF PT
       call PT_run(TimeSimulation, TimeSimulationLimit)     !^CMP IF PT
    case(PW_)                                               !^CMP IF PW
       call PW_run(TimeSimulation, TimeSimulationLimit)     !^CMP IF PW
    case(RB_)                                               !^CMP IF RB
       call RB_run(TimeSimulation, TimeSimulationLimit)     !^CMP IF RB
    case(SC_)                                               !^CMP IF SC
       call SC_run(TimeSimulation, TimeSimulationLimit)     !^CMP IF SC
    case(SP_)                                               !^CMP IF SP
       call SP_run(TimeSimulation, TimeSimulationLimit)     !^CMP IF SP
    case(UA_)                                               !^CMP IF UA
       call UA_run(TimeSimulation, TimeSimulationLimit)     !^CMP IF UA
    case(CZ_)                                               !^CMP IF CZ
       call CZ_run(TimeSimulation, TimeSimulationLimit)     !^CMP IF CZ
    case default
       call CON_stop(NameSub//' SWMF_ERROR incorrect iComp value')
    end select
    call timing_stop(NameComp_I(iComp)//'_run')

    ! Reset the number of threads
    !$ MaxThread = omp_get_max_threads()
    !$ call omp_set_num_threads(MaxThread)

    if(DoTestMe)write(*,*) NameSub,' finished for ',NameComp_I(iComp)

  end subroutine run_comp_id
  !============================================================================
end module CON_wrapper
!==============================================================================
