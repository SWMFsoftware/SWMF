!^CMP COPYRIGHT UM
!
!QUOTE: \clearpage
!
!BOP -------------------------------------------------------------------
!
!QUOTE: \section{CON/Interface: between CON and Components and between Components}
!
!MODULE: CON_wrapper - Wrapper for component specific subroutines
!
!DESCRIPTION:
!
! Contains general subroutines which call the specific subroutines
! in the component wrappers. The components are identified by their 
! name or the component ID. 
!
!INTERFACE:

module CON_wrapper
  !
  !USES:
  !
  use CON_world
  use CON_comp_param
  use CON_comp_info

  implicit none

  private   ! except

  !PUBLIC MEMBER FUNCTIONS:

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

  !REVISION HISTORY:
  ! 09Jul03 - Gabor Toth <gtoth@umich.edu> - initial prototype/prolog/code
  ! 20Aug03 - O. Volberg and G. Toth - added new methods to initialize,
  !                                    run session, finalize and save restart.
  !EOP ___________________________________________________________________

  character(len=*),parameter :: NameMod='CON_wrapper'

contains

  !BOP -------------------------------------------------------------------
  !IROUTINE: set_param_name - Set parameters for component given by name
  !
  !INTERFACE:
  subroutine set_param_name(NameComp,TypeAction)

    implicit none

    !INPUT PARAMETERS:
    !
    character(len=*), intent(in)           :: NameComp, TypeAction

    !DESCRIPTION:
    ! This routine calls set\_param\_id with the component ID defined by the 
    ! name of the component.

    !REVISION HISTORY:
    ! 21Jul03 - G. Toth <gtoth@umich.edu> - Pass CompInfo for safety
    ! 19Jul03 - G. Toth <gtoth@umich.edu> - simplified and generalized version
    ! 09Jul03 - G. Toth <gtoth@umich.edu> - initial prototype/prolog/code
    !EOP ___________________________________________________________________

    character(len=*),parameter :: NameSub=NameMod//'.set_param_name'

    !-------------------------------------------------------------------
    call set_param_id(i_comp(NameComp),TypeAction)

  end subroutine set_param_name

  !BOP -------------------------------------------------------------------
  !IROUTINE: set_param_id - Set parameters for component given by ID
  !
  !INTERFACE:
  subroutine set_param_id(iComp,TypeAction)

    implicit none

    !INPUT PARAMETERS:
    !
    integer, intent(in)          :: iComp
    character(len=*), intent(in) :: TypeAction

    !DESCRIPTION:
    ! This routine calls $**$\_set\_param for component $**$ 
    ! defined with the component ID.

    !REVISION HISTORY:
    ! 19Jul03 - G. Toth <gtoth@umich.edu> - simplified and generalized version
    ! 09Jul03 - G. Toth <gtoth@umich.edu> - initial prototype/prolog/code
    ! 03Jun08 - R Oran  <oran@umich.edu>  - add component OH
    ! 16Sep09 - R. ORan <oran@umich.edu> - add component LC
    !EOP ___________________________________________________________________

    character(len=*),parameter :: NameSub=NameMod//'.set_param_id'
    logical :: DoTest, DoTestMe
    type(CompInfoType) :: CompInfo
    !-------------------------------------------------------------------
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

    select case(iComp)
    case(GM_)                                     !^CMP IF GM
       call GM_set_param(CompInfo,TypeAction)     !^CMP IF GM
    case(IE_)                                     !^CMP IF IE
       call IE_set_param(CompInfo,TypeAction)     !^CMP IF IE
    case(IH_)                                     !^CMP IF IH
       call IH_set_param(CompInfo,TypeAction)     !^CMP IF IH
    case(OH_)                                     !^CMP IF OH                  
       call OH_set_param(CompInfo,TypeAction)     !^CMP IF OH  
    case(IM_)                                     !^CMP IF IM
       call IM_set_param(CompInfo,TypeAction)     !^CMP IF IM
    case(LA_)                                     !^CMP IF LA
       call LA_set_param(CompInfo,TypeAction)     !^CMP IF LA
    case(PS_)                                     !^CMP IF PS
       call PS_set_param(CompInfo,TypeAction)     !^CMP IF PS
    case(PW_)                                     !^CMP IF PW
       call PW_set_param(CompInfo,TypeAction)     !^CMP IF PW
    case(RB_)                                     !^CMP IF RB
       call RB_set_param(CompInfo,TypeAction)     !^CMP IF RB
    case(SC_)                                     !^CMP IF SC
       call SC_set_param(CompInfo,TypeAction)     !^CMP IF SC
    case(LC_)                                     !^CMP IF LC
       call LC_set_param(CompInfo,TypeAction)     !^CMP IF LC
    case(SP_)                                     !^CMP IF SP
       call SP_set_param(CompInfo,TypeAction)     !^CMP IF SP
    case(UA_)                                     !^CMP IF UA
       call UA_set_param(CompInfo,TypeAction)     !^CMP IF UA
    case default
       call CON_stop(NameSub//' '//TypeAction//&
            ' SWMF_ERROR: not implemented for component'//NameComp_I(iComp))
    end select

    ! Store returned component information
    select case(TypeAction)
    case('VERSION')
       call put_comp_info(iComp,Use=CompInfo%Use,Version=CompInfo%Version,&
            NameVersion=CompInfo%NameVersion)
    end select

    if(DoTestMe)write(*,*) NameSub,' finished TypeAction=',TypeAction, &
         ' for ',NameComp_I(iComp)

  end subroutine set_param_id

  !BOP -------------------------------------------------------------------
  !IROUTINE: get_version_name - call **_get_version for named component
  !
  !DESCRIPTION:
  ! Obtain ON/OFF status, version name and version number 
  ! for component specified by its name! 

  !INTERFACE:
  subroutine get_version_name(NameComp,IsOn,NameVersion,Version)

    implicit none

    !INPUT PARAMETERS:
    character (len=*),            intent(in)  :: NameComp
    !OUTPUT PARAMETERS:
    logical,                      intent(out) :: IsOn
    character (len=lNameVersion), intent(out) :: NameVersion
    real,                         intent(out) :: Version

    !REVISION HISTORY:
    ! 15Jul03 - G. Toth <gtoth@umich.edu> - initial prototype/prolog/code
    !EOP ___________________________________________________________________

    character(len=*),parameter :: NameSub=NameMod//'::get_version_name'
    !-------------------------------------------------------------------
    if(.not.is_valid_comp_name(NameComp)) call CON_stop(NameSub// &
         ' SWMF_ERROR '//NameComp//' is not a valid component name')

    call get_version_id(i_comp(NameComp),IsOn,NameVersion,Version)

  end subroutine get_version_name

  !BOP -------------------------------------------------------------------
  !IROUTINE: get_version_id - call **_version for component given by ID
  !
  !DESCRIPTION:
  ! Obtain ON/OFF status, version name and version number 
  ! for component specified by its ID! 

  !INTERFACE:
  subroutine get_version_id(iComp,IsOn,NameVersion,Version)

    use CON_comp_info

    !INPUT PARAMETERS:
    integer,                      intent(in)  :: iComp
    !OUTPUT PARAMETERS:
    logical,                      intent(out) :: IsOn
    character (len=lNameVersion), intent(out) :: NameVersion
    real,                         intent(out) :: Version

    !REVISION HISTORY:
    ! 21Jul03 - G. Toth <gtoth@umich.edu> - use temporary CompInfo
    ! 15Jul03 - G. Toth <gtoth@umich.edu> - initial prototype/prolog/code
    ! 03Jun08 - R. Oran <oran@umich.edu>  - added component OH
    !EOP ___________________________________________________________________

    character(len=*),parameter :: NameSub=NameMod//'::get_version_id'
    logical :: DoTest, DoTestMe
    type(CompInfoType) :: CompInfo

    !-------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTestMe)write(*,*) NameSub,' starting for ',NameComp_I(iComp)

    call check_i_comp(iComp,NameSub)

    select case(iComp)
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
    case(LA_)                                    !^CMP IF LA
       call LA_set_param(CompInfo,'VERSION')     !^CMP IF LA
    case(PS_)                                    !^CMP IF PS
       call PS_set_param(CompInfo,'VERSION')     !^CMP IF PS
    case(PW_)                                    !^CMP IF PW
       call PW_set_param(CompInfo,'VERSION')     !^CMP IF PW
    case(RB_)                                    !^CMP IF RB
       call RB_set_param(CompInfo,'VERSION')     !^CMP IF RB
    case(SC_)                                    !^CMP IF SC
       call SC_set_param(CompInfo,'VERSION')     !^CMP IF SC
    case(LC_)                                    !^CMP IF LC
       call LC_set_param(CompInfo,'VERSION')     !^CMP IF LC
    case(SP_)                                    !^CMP IF SP
       call SP_set_param(CompInfo,'VERSION')     !^CMP IF SP
    case(UA_)                                    !^CMP IF UA
       call UA_set_param(CompInfo,'VERSION')     !^CMP IF UA 
    case default
       call put(CompInfo,Use=.false.,NameVersion='not implemented',Version=0.0)
    end select

    call get(CompInfo,Use=IsOn,NameVersion=NameVersion,Version=Version)
    
  end subroutine get_version_id

  !BOP -------------------------------------------------------------------
  !IROUTINE: init_session_comp_id - call **_init_session for component given by ID
  !INTERFACE:
  subroutine init_session_comp_id(iComp, iSession, TimeSimulation)

    !INPUT PARAMETERS:
    integer,  intent(in) :: iComp            ! component ID
    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    !DESCRIPTION:
    ! Initialize component for the session. Session numbers start from 1.

    !REVISION HISTORY:
    ! 18Aug03 - G. Toth <gtoth@umich.edu> O. Volberg <volov@umich.edu> 
    !         - initial prototype/prolog/code
    ! 03Jun   - R. Oran <oran@umich.edu> - added component OH
    !EOP ___________________________________________________________________
    character(len=*), parameter :: NameSub = NameMod//'::init_session_comp_id'
    logical :: DoTest, DoTestMe
    !-------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    call check_i_comp(iComp,NameSub)

    if(DoTestMe)write(*,*) NameSub,' starting for ',NameComp_I(iComp),&
         ': iSession, Time=',iSession, TimeSimulation

    if(.not.use_comp(iComp) .or. .not.is_proc(iComp)) RETURN

    select case(iComp)
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
    case(LA_)                                            !^CMP IF LA
       call LA_init_session(iSession,TimeSimulation)     !^CMP IF LA
    case(PS_)                                            !^CMP IF PS
       call PS_init_session(iSession,TimeSimulation)     !^CMP IF PS
    case(PW_)                                            !^CMP IF PW
       call PW_init_session(iSession,TimeSimulation)     !^CMP IF PW
    case(RB_)                                            !^CMP IF RB
       call RB_init_session(iSession,TimeSimulation)     !^CMP IF RB
    case(SC_)                                            !^CMP IF SC
       call SC_init_session(iSession,TimeSimulation)     !^CMP IF SC
    case(LC_)                                            !^CMP IF LC
       call LC_init_session(iSession,TimeSimulation)     !^CMP IF LC
    case(SP_)                                            !^CMP IF SP
       call SP_init_session(iSession,TimeSimulation)     !^CMP IF SP
    case(UA_)                                            !^CMP IF UA
       call UA_init_session(iSession,TimeSimulation)     !^CMP IF UA
    case default
       call CON_stop(NameSub//' SWMF_ERROR incorrect iComp value')
    end select

    if(DoTestMe)write(*,*) NameSub,' finished for ',NameComp_I(iComp)

  end subroutine init_session_comp_id

  !BOP -------------------------------------------------------------------
  !IROUTINE: finalize_comp_id - call **_finalize for component given by ID
  !
  !INTERFACE:
  subroutine finalize_comp_id(iComp, TimeSimulation)

    !INPUT PARAMETERS:
    integer,  intent(in) :: iComp            ! component ID
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    !DESCRIPTION:
    ! Finalize the component at the end of the run, save plots, restart info,
    ! report errors, etc.

    !REVISION HISTORY:
    ! 18Aug03 - G. Toth <gtoth@umich.edu> O. Volberg <volov@umich.edu> 
    !         - initial prototype/prolog/code
    ! 03Jun08 - R. Oran <oran@umich.edu> - added component OH
    !EOP ___________________________________________________________________
    character(len=*), parameter :: NameSub = NameMod//'::finalize_comp_id'
    logical :: DoTest, DoTestMe
    !-------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)
    
    call check_i_comp(iComp,NameSub)

    if(DoTestMe)write(*,*) NameSub,' starting for ',NameComp_I(iComp),&
         ': Time=', TimeSimulation

    if(.not.use_comp(iComp) .or. .not.is_proc(iComp)) RETURN

    select case(iComp)
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
    case(LA_)                               !^CMP IF LA
       call LA_finalize(TimeSimulation)     !^CMP IF LA
    case(PS_)                               !^CMP IF PS
       call PS_finalize(TimeSimulation)     !^CMP IF PS
    case(PW_)                               !^CMP IF PW
       call PW_finalize(TimeSimulation)     !^CMP IF PW
    case(RB_)                               !^CMP IF RB
       call RB_finalize(TimeSimulation)     !^CMP IF RB
    case(SC_)                               !^CMP IF SC
       call SC_finalize(TimeSimulation)     !^CMP IF SC
    case(LC_)                               !^CMP IF LC
       call LC_finalize(TimeSimulation)     !^CMP IF LC
    case(SP_)                               !^CMP IF SP
       call SP_finalize(TimeSimulation)     !^CMP IF SP
    case(UA_)                               !^CMP IF UA
       call UA_finalize(TimeSimulation)     !^CMP IF UA
    case default
       call CON_stop(NameSub//' SWMF_ERROR incorrect iComp value')
    end select

    if(DoTestMe)write(*,*) NameSub,' finished for ',NameComp_I(iComp)

  end subroutine finalize_comp_id

  !BOP -------------------------------------------------------------------
  !IROUTINE: save_restart_comp_id - call **_save_restart for component given by ID

  !INTERFACE:
  subroutine save_restart_comp_id(iComp, TimeSimulation)

    !INPUT PARAMETERS:
    integer,  intent(in) :: iComp            ! component ID
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    !DESCRIPTION:
    ! Save restart information for the component.

    !REVISION HISTORY:
    ! 18Aug03 - G. Toth <gtoth@umich.edu> O. Volberg <volov@umich.edu> 
    !         - initial prototype/prolog/code
    ! 03Jun08 - R. Oran <oran@umich.edu> - added component OH
    !EOP ___________________________________________________________________
    character(len=*), parameter :: NameSub = NameMod//'::save_restart_comp_id'
    logical :: DoTest, DoTestMe
    !-------------------------------------------------------------------
    
    call check_i_comp(iComp,NameSub)

    if(.not.use_comp(iComp) .or. .not.is_proc(iComp)) RETURN

    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTestMe)write(*,*) NameSub,' starting for ',NameComp_I(iComp),&
         ': Time=', TimeSimulation

    select case(iComp)
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
    case(LA_)                                   !^CMP IF LA
       call LA_save_restart(TimeSimulation)     !^CMP IF LA
    case(PS_)                                   !^CMP IF PS
       call PS_save_restart(TimeSimulation)     !^CMP IF PS
    case(PW_)                                   !^CMP IF PW
       call PW_save_restart(TimeSimulation)     !^CMP IF PW
    case(RB_)                                   !^CMP IF RB
       call RB_save_restart(TimeSimulation)     !^CMP IF RB
    case(SC_)                                   !^CMP IF SC
       call SC_save_restart(TimeSimulation)     !^CMP IF SC
    case(LC_)                                   !^CMP IF LC
       call LC_save_restart(TImeSimulation)     !^CMP IF LC
    case(SP_)                                   !^CMP IF SP
       call SP_save_restart(TimeSimulation)     !^CMP IF SP
    case(UA_)                                   !^CMP IF UA
       call UA_save_restart(TimeSimulation)     !^CMP IF UA
    case default
       call CON_stop(NameSub//' SWMF_ERROR incorrect iComp value')
    end select

    if(DoTestMe)write(*,*) NameSub,' finished for ',NameComp_I(iComp)

  end subroutine save_restart_comp_id

  !BOP -------------------------------------------------------------------
  !IROUTINE: run_comp_id - call **_run for component given by ID
  !INTERFACE:
  subroutine run_comp_id(iComp, TimeSimulation, TimeSimulationLimit)

    !INPUT PARAMETERS:
    integer,  intent(in) :: iComp               ! component ID
    real,     intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded

    !INPUT/OUTPUT ARGUMENTS:
    real,     intent(inout) :: TimeSimulation   ! current time of component

    !DESCRIPTION:
    ! Do one time step or solve for the component.
    ! The component should update TimeSimulation.
    ! The updated TimeSimulation should not exceed TimeSimulationLimit.

    !REVISION HISTORY:
    ! 18Aug03 - G. Toth <gtoth@umich.edu> O. Volberg <volov@umich.edu> 
    !         - initial prototype/prolog/code
    ! 03Jun08 - R. Oran <oran@umich.edu - added component OH
    !EOP ___________________________________________________________________
    character(len=*), parameter :: NameSub = NameMod//'::run_comp_id'
    logical :: DoTest, DoTestMe
    !-------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    call check_i_comp(iComp,NameSub)

    if(.not.use_comp(iComp) .or. .not.is_proc(iComp)) RETURN

    if(DoTestMe)write(*,*) NameSub,' starting for ',NameComp_I(iComp),&
         ': Time, Limit=', TimeSimulation, TimeSimulationLimit

    call timing_start(NameComp_I(iComp)//'_run')
    select case(iComp)
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
    case(LA_)                                               !^CMP IF LA
       call LA_run(TimeSimulation, TimeSimulationLimit)     !^CMP IF LA
    case(PS_)                                               !^CMP IF PS
       call PS_run(TimeSimulation, TimeSimulationLimit)     !^CMP IF PS
    case(PW_)                                               !^CMP IF PW
       call PW_run(TimeSimulation, TimeSimulationLimit)     !^CMP IF PW
    case(RB_)                                               !^CMP IF RB
       call RB_run(TimeSimulation, TimeSimulationLimit)     !^CMP IF RB
    case(SC_)                                               !^CMP IF SC
       call SC_run(TimeSimulation, TimeSimulationLimit)     !^CMP IF SC
    case(LC_)                                               !^CMP IF LC
       call LC_run(TimeSimulation, TimeSimulationLimit)     !^CMP IF LC
    case(SP_)                                               !^CMP IF SP
       call SP_run(TimeSimulation, TimeSimulationLimit)     !^CMP IF SP
    case(UA_)                                               !^CMP IF UA
       call UA_run(TimeSimulation, TimeSimulationLimit)     !^CMP IF UA
    case default
       call CON_stop(NameSub//' SWMF_ERROR incorrect iComp value')
    end select
    call timing_stop(NameComp_I(iComp)//'_run')

    if(DoTestMe)write(*,*) NameSub,' finished for ',NameComp_I(iComp)

  end subroutine run_comp_id

end module CON_wrapper
