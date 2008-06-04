!^CMP COPYRIGHT UM
!
!BOP
!
!MODULE: CON_couple_all - provide access to all the couplers
!
!DESCRIPTION:
! This module uses all the pairwise couplers. It also provides public
! methods for initializing the couplers, or accessing the couplers
! based on the two component IDs.
!
!INTERFACE:
module CON_couple_all

  !USES:
  use CON_comp_param
  use CON_world, ONLY: use_comp, is_proc, i_proc
  !^CMP IF GM BEGIN
  use CON_couple_ih_gm        !^CMP IF IH
  use CON_couple_gm_ie        !^CMP IF IE
  use CON_couple_gm_im        !^CMP IF IM
  use CON_couple_gm_pw        !^CMP IF PW
  use CON_couple_gm_rb        !^CMP IF RB
  !^CMP END GM
  !^CMP IF IE BEGIN
  use CON_couple_ie_im        !^CMP IF IM
  use CON_couple_ie_pw        !^CMP IF PW
  use CON_couple_ie_rb        !^CMP IF RB
  use CON_couple_ie_ua        !^CMP IF UA
  use CON_couple_ie_ps
  !^CMP END IE
  !^CMP IF IH BEGIN
  use CON_couple_ih_sc        !^CMP IF SC
  use CON_couple_ih_oh        !^CMP IF OH
  !^CMP END IH
  use CON_couple_mh_sp        !^CMP IF SP

  implicit none

  private   ! except

  !PUBLIC MEMBER FUNCTIONS:

  public :: couple_all_init ! initialize all couplers
  public :: couple_two_comp     ! couple 2 components based on their IDs

  !REVISION HISTORY:
  ! 27Aug03 - G. Toth <gtoth@umich.edu> initial prototype/prolog/code
  ! 14Jan05 - G. Toth commented out _swmf version of GM-IE couplers.
  ! 03Jun08 - R. Oran <oran@umich.edu> added ih_oh coupling
  !EOP
  character(len=*), parameter :: NameMod='CON_couple_all'

contains

  !BOP -------------------------------------------------------------------
  !IROUTINE: couple_all_init - initialize all the couplers
  !INTERFACE:
  subroutine couple_all_init
    !EOP
    !BOC
    !                                                     ^CMP IF GM BEGIN
    if(use_comp(GM_).and.use_comp(IE_))call couple_gm_ie_init  !^CMP IF IE
    if(use_comp(GM_).and.use_comp(IM_))call couple_gm_im_init  !^CMP IF IM
    if(use_comp(GM_).and.use_comp(PW_))call couple_gm_pw_init  !^CMP IF PW
    if(use_comp(GM_).and.use_comp(RB_))call couple_gm_rb_init  !^CMP IF RB
    if(use_comp(IH_).and.use_comp(GM_))call couple_ih_gm_init  !^CMP IF IH
    !                                                     ^CMP END GM
    !                                                     ^CMP IF IE BEGIN
    if(use_comp(IE_).and.use_comp(IM_))call couple_ie_im_init  !^CMP IF IM
    if(use_comp(IE_).and.use_comp(PW_))call couple_ie_pw_init  !^CMP IF PW
    if(use_comp(IE_).and.use_comp(RB_))call couple_ie_rb_init  !^CMP IF RB
    if(use_comp(IE_).and.use_comp(UA_))call couple_ie_ua_init  !^CMP IF UA
    if(use_comp(IE_).and.use_comp(PS_))call couple_ie_ps_init
    !                                                     ^CMP END IE
    !                                                     ^CMP IF IH BEGIN
    if(use_comp(IH_).and.use_comp(SC_))call couple_ih_sc_init  !^CMP IF SC
    if(use_comp(IH_).and.use_comp(OH_))call couple_oh_ih_init  !^CMP IF OH
    !                                                     ^CMP END IH
    if((&                                                 !^CMP IF SP BEGIN
         use_comp(IH_).or.&                               !^CMP IF IH
         use_comp(SC_).or.&                               !^CMP IF SC
         .false.).and.use_comp(SP_))call couple_mh_sp_init !^CMP END SP
    !EOC
  end subroutine couple_all_init

  !BOP =======================================================================
  !IROUTINE: couple_two_comp - call couple_**_** for components given by IDs
  !INTERFACE:
  subroutine couple_two_comp(iCompSource, iCompTarget, TimeSimulation)

    !INPUT PARAMETERS:
    integer,  intent(in) :: iCompSource, iCompTarget ! component IDs
    real,     intent(in) :: TimeSimulation           ! coupling simulation time

    !DESCRIPTION:
    ! Couple two components given with their IDs. The simulation time
    ! is shared at the time of coupling. Call the appropriate coupling.
    ! Stop with an error for invalid component pairs.

    !REVISION HISTORY:
    ! 27Aug03 - G. Toth <gtoth@umich.edu> initial prototype/prolog/code
    ! 03Jun08 - R. Oran <oran@umich.edu> add IH OH coupling
    !EOP

    character(len=*), parameter :: NameSub = NameMod//'::couple_two_comp'

    logical :: DoTest,DoTestMe
    !-------------------------------------------------------------------
    
    call check_i_comp(iCompSource,NameSub//': source')
    call check_i_comp(iCompTarget,NameSub//': target')

    !\
    ! Return if any component is not used or if the PE is
    ! used by neither components.
    !/
    if(.not.(use_comp(iCompSource))) RETURN
    if(.not.(use_comp(iCompTarget))) RETURN
    if(.not.(is_proc(iCompSource).or.is_proc(iCompTarget))) RETURN

    call CON_set_do_test(NameSub,DoTest,DoTestMe)
    if(DoTest)write(*,*)NameSub,': coupling iProc=',i_proc(),' ',&
         NameComp_I(iCompSource),' --> ',NameComp_I(iCompTarget)


    call timing_start(NameComp_I(iCompSource)//'_'//NameComp_I(iCompTarget)//'_couple')

    select case(iCompSource)
    case(SC_)                                 !^CMP IF SC BEGIN
       select case(iCompTarget)
       case(IH_)                              !^CMP IF IH
          call couple_sc_ih(TimeSimulation)   !^CMP IF IH
       case(SP_)                              !^CMP IF SP
          call couple_sc_sp(TimeSimulation)   !^CMP IF SP
       case default                           
          call error
       end select                             !^CMP END SC
    case(IH_)                                 !^CMP IF IH BEGIN
       select case(iCompTarget)
       case(GM_)                                   !^CMP IF GM
          call couple_ih_gm(TimeSimulation)        !^CMP IF GM
       case(SC_)                                   !^CMP IF SC
          call couple_ih_sc(TimeSimulation)        !^CMP IF SC
       case(OH_)                                   !^CMP IF OH
          call couple_ih_oh(TimeSimulation)        !^CMP IF OH
       case(SP_)                                   !^CMP IF SP
          call couple_ih_sp(TimeSimulation)        !^CMP IF SP
       case default
          call error
       end select                             !^CMP END IH
    case(OH_)
       select case(iCompTarget)               !^CMP IF OH BEGIN
       case(IH_)                                   !^CMP IF IH
          call couple_oh_ih(TimeSimulation)        !^CMP IF IH
       case default
          call error
       end select                                  !^CMP END OH
    case(GM_)                                 !^CMP IF GM BEGIN
       select case(iCompTarget)
       case(IE_)                                   !^CMP IF IE
          call couple_gm_ie(TimeSimulation)        !^CMP IF IE
       case(IM_)                                   !^CMP IF IM
          call couple_gm_im(TimeSimulation)        !^CMP IF IM
       case(PW_)                                   !^CMP IF PW
          call couple_gm_pw(TimeSimulation)        !^CMP IF PW
       case(RB_)                                   !^CMP IF RB
          call couple_gm_rb(TimeSimulation)        !^CMP IF RB
       case default
          call error
       end select                             !^CMP END GM
    case(UA_)                                 !^CMP IF UA BEGIN
       select case(iCompTarget)
       case(IE_)                                   !^CMP IF IE
          call couple_ua_ie(TimeSimulation)        !^CMP IF IE
       case default
          call error
       end select                             !^CMP END UA 
    case(IE_)                                 !^CMP IF IE BEGIN
       select case(iCompTarget)
       case(GM_)                                   !^CMP IF GM
          call couple_ie_gm(TimeSimulation)        !^CMP IF GM
       case(IM_)                                   !^CMP IF IM
          call couple_ie_im(TimeSimulation)        !^CMP IF IM
       case(PS_)
          call couple_ie_ps(TimeSimulation)
       case(PW_)                                   !^CMP IF PW
          call couple_ie_pw(TimeSimulation)        !^CMP IF PW
       case(RB_)                                   !^CMP IF RB
          call couple_ie_rb(TimeSimulation)        !^CMP IF RB
       case(UA_)                                   !^CMP IF UA
          call couple_ie_ua(TimeSimulation)        !^CMP IF UA
       case default
          call error
       end select                             !^CMP END IE
    case(IM_)                                 !^CMP IF IM BEGIN
       select case(iCompTarget)
       case(GM_)                                   !^CMP IF GM
          call couple_im_gm(TimeSimulation)        !^CMP IF GM
       case default
          call error
       end select                             !^CMP END IM
    case(PW_)                                 !^CMP IF PW BEGIN
       select case(iCompTarget)
       case(GM_)                                   !^CMP IF GM
          call couple_pw_gm(TimeSimulation)        !^CMP IF GM
       case default
          call error
       end select                             !^CMP END PW
    case default
       call CON_stop(NameSub//&
            ' SWMF_ERROR: no coupling implemented from source '// &
            NameComp_I(iCompSource))
    end select
    call timing_stop(NameComp_I(iCompSource)//'_'//NameComp_I(iCompTarget)//'_couple')
    if(DoTest)write(*,*)NameSub,': finished coupling iProc=',i_proc(),' ',&
         NameComp_I(iCompSource),' --> ',NameComp_I(iCompTarget)

  contains

    subroutine error
      call CON_stop(NameSub//' SWMF_ERROR: no coupling implemented for '//&
           NameComp_I(iCompSource)//' --> '//&
           NameComp_I(iCompTarget))
    end subroutine error

  end subroutine couple_two_comp

end module CON_couple_all
