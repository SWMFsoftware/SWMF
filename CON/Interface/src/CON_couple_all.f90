!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!
!
!
! This module uses all the pairwise couplers. It also provides public
! methods for initializing the couplers, or accessing the couplers
! based on the two component IDs.
!
module CON_couple_all

  use CON_comp_param
  use CON_world,   ONLY: use_comp, is_proc, i_proc, lComp_I
  use CON_coupler, ONLY: &
       iCompSourceCouple, iCompTargetCouple, &
       iVar_V, iVar_VCC, nVarCouple, nVarCouple_CC, &
       DoCoupleVar_V, DoCoupleVar_VCC, &
       nVarBuffer, nVarBuffer_CC, &
       iVarSource_V, iVarSource_VCC, &
       iVarTarget_V, iVarTarget_VCC
  use CON_bline, ONLY: UseBLine_C, IsSource4BL_C
  use ModUtilities, ONLY: CON_set_do_test, CON_stop

  !^CMP IF GM BEGIN
  use CON_couple_ih_gm        !^CMP IF IH
  use CON_couple_sc_gm        !^CMP IF SC
  use CON_couple_gm_ie        !^CMP IF IE
  use CON_couple_gm_im        !^CMP IF IM
  use CON_couple_gm_ps        !^CMP IF PS
  use CON_couple_gm_pt        !^CMP IF PT
  use CON_couple_gm_pc        !^CMP IF PC
  use CON_couple_gm_pw        !^CMP IF PW
  use CON_couple_gm_rb        !^CMP IF RB
  use CON_couple_ua_gm        !^CMP IF UA
  !^CMP END GM
  !^CMP IF IE BEGIN
  use CON_couple_ie_im        !^CMP IF IM
  use CON_couple_ie_pw        !^CMP IF PW
  use CON_couple_ie_rb        !^CMP IF RB
  use CON_couple_ie_ua        !^CMP IF UA
  use CON_couple_ie_ps        !^CMP IF PS
  !^CMP END IE
  !^CMP IF IH BEGIN
  use CON_couple_ih_sc        !^CMP IF SC
  use CON_couple_ih_oh        !^CMP IF OH
  use CON_couple_ih_pt        !^CMP IF PT
  !^CMP END IH
  !^CMP IF SC BEGIN
  use CON_couple_ee_sc        !^CMP IF EE
  use CON_couple_sc_pt        !^CMP IF PT
  !^CMP END SC
  use CON_couple_mh_sp        !^CMP IF SP
  !^CMP IF OH BEGIN
  use CON_couple_oh_pt	      !^CMP IF PT
  !^CMP END OH

  implicit none

  private   ! except

  public :: couple_all_init ! initialize all couplers
  public :: couple_two_comp ! couple 2 components based on their IDs

  ! revision history:
  ! 27Aug03 - G. Toth <gtoth@umich.edu> initial prototype/prolog/code
  ! 14Jan05 - G. Toth commented out _swmf version of GM-IE couplers.
  ! 03Jun08 - R. Oran <oran@umich.edu> added IH-OH coupling
  ! 23Jul08 - A. Ridley added IM-IE coupling
  ! 24Sep13 - G. Toth added GM-PT coupling
  ! 18Feb14 - L. Daldorff added GM-PC coupling
  ! 12Mar15 - A. Michael added OH-PT coupling
  ! 15Sep20 - V. Tenishev added IH-PT coupling
  character(len=*), parameter :: NameMod='CON_couple_all'

contains
  !============================================================================

  subroutine couple_all_init
    !                                                     ^CMP IF GM BEGIN
    !--------------------------------------------------------------------------
    if(use_comp(GM_).and.use_comp(IE_))call couple_gm_ie_init  !^CMP IF IE
    if(use_comp(GM_).and.use_comp(IM_))call couple_gm_im_init  !^CMP IF IM
    if(use_comp(GM_).and.use_comp(PS_))call couple_gm_ps_init  !^CMP IF PS
    if(use_comp(GM_).and.use_comp(PW_))call couple_gm_pw_init  !^CMP IF PW
    if(use_comp(GM_).and.use_comp(RB_))call couple_gm_rb_init  !^CMP IF RB
    if(use_comp(IH_).and.use_comp(GM_))call couple_ih_gm_init  !^CMP IF IH
    if(use_comp(SC_).and.use_comp(GM_))call couple_sc_gm_init  !^CMP IF IH
    if(use_comp(GM_).and.use_comp(PT_))call couple_gm_pt_init  !^CMP IF PT
    if(use_comp(GM_).and.use_comp(PC_))call couple_gm_pc_init  !^CMP IF PC
    if(use_comp(UA_).and.use_comp(GM_))call couple_ua_gm_init  !^CMP IF UA
    !                                                     ^CMP END GM
    !                                                     ^CMP IF IE BEGIN
    if(use_comp(IE_).and.use_comp(IM_))call couple_ie_im_init  !^CMP IF IM
    if(use_comp(IE_).and.use_comp(PW_))call couple_ie_pw_init  !^CMP IF PW
    if(use_comp(IE_).and.use_comp(RB_))call couple_ie_rb_init  !^CMP IF RB
    if(use_comp(IE_).and.use_comp(UA_))call couple_ie_ua_init  !^CMP IF UA
    if(use_comp(IE_).and.use_comp(PS_))call couple_ie_ps_init  !^CMP IF PS
    !                                                     ^CMP END IE
    !                                                     ^CMP IF IH BEGIN
    if(use_comp(IH_).and.use_comp(SC_))call couple_ih_sc_init  !^CMP IF SC
    if(use_comp(IH_).and.use_comp(OH_))call couple_ih_oh_init  !^CMP IF OH
    !                                                     ^CMP END IH
    !                                                     ^CMP IF SC BEGIN
    if(use_comp(SC_).and.use_comp(EE_))call couple_ee_sc_init  !^CMP IF EE
    !                                                     ^CMP END SC
    if(UseBLine_C(SP_).and.(&                             !^CMP IF SP BEGIN
         IsSource4Bl_C(OH_).or.&                          !^CMP IF OH
         IsSource4Bl_C(IH_).or.&                          !^CMP IF IH
         IsSource4Bl_C(SC_).or.&                          !^CMP IF SC
         .false.))call couple_mh_sp_init                  !^CMP END SP
    if(UseBLine_C(PT_).and.(&                             !^CMP IF PT BEGIN
         IsSource4Bl_C(OH_).or.&                          !^CMP IF OH
         IsSource4Bl_C(IH_).or.&                          !^CMP IF IH
         IsSource4Bl_C(SC_).or.&                          !^CMP IF SC
         .false.))then
       call couple_mh_sp_init
    else
       if(use_comp(SC_).and.use_comp(PT_))call couple_sc_pt_init !^CMP IF SC
       if(use_comp(IH_).and.use_comp(PT_))call couple_ih_pt_init !^CMP IF IH
    end if                                                !^CMP END PT
    !	 				 		       ^CMP IF OH BEGIN
    if(use_comp(OH_).and.use_comp(PT_))call couple_oh_pt_init  !^CMP IF PT
    !					    		       ^CMP END OH
  end subroutine couple_all_init
  !============================================================================

  subroutine couple_two_comp(iCompSource, iCompTarget, TimeSimulation)

    integer,  intent(in) :: iCompSource, iCompTarget ! component IDs
    real,     intent(in) :: TimeSimulation           ! coupling simulation time

    ! Couple two components given with their IDs. The simulation time
    ! is shared at the time of coupling. Call the appropriate coupling.
    ! Stop with an error for invalid component pairs.

    ! revision history:
    ! 27Aug03 - G. Toth <gtoth@umich.edu> initial prototype/prolog/code
    ! 03Jun08 - R. Oran <oran@umich.edu> add IH OH coupling
    ! 22Dec11 - R. Oran   added capability to use global mpi coupler.

    ! indexes by component list
    integer:: lCompSource, lCompTarget

    logical :: DoTest, DoTestMe

    character(len=*), parameter:: NameSub = 'couple_two_comp'
    !--------------------------------------------------------------------------
    call check_i_comp(iCompSource,NameSub//': source')
    call check_i_comp(iCompTarget,NameSub//': target')

    ! Return if any component is not used
    if(.not.(use_comp(iCompSource))) RETURN
    if(.not.(use_comp(iCompTarget))) RETURN

    ! Return if the PE is not used by either components
    if(.not.(is_proc(iCompSource).or.is_proc(iCompTarget))) RETURN

    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTest)write(*,*)NameSub,': coupling iProc=',i_proc(),' ',&
         NameComp_I(iCompSource),' --> ',NameComp_I(iCompTarget)

    call timing_start(NameComp_I(iCompSource)// &
         '_'//NameComp_I(iCompTarget)//'_couple')

    iVar_V        = iVar_VCC(:,iCompSource,iCompTarget)
    nVarCouple    = nVarCouple_CC(iCompSource,iCompTarget)
    DoCoupleVar_V = DoCoupleVar_VCC(:,iCompSource,iCompTarget)

    ! Make the component indexes public through CON_coupler
    iCompSourceCouple = iCompSource
    iCompTargetCouple = iCompTarget

    lCompSource = lComp_I(iCompSource)
    lCompTarget = lComp_I(iCompTarget)

    if(allocated(iVarSource_VCC))then
       nVarBuffer    = nVarBuffer_CC(lCompSource,lCompTarget)
       iVarSource_V  = iVarSource_VCC(:,lCompSource,lCompTarget)
       iVarTarget_V  = iVarTarget_VCC(:,lCompSource,lCompTarget)
    end if

    select case(iCompSource)
    case(EE_)                                 !^CMP IF EE BEGIN
       select case(iCompTarget)
       case(SC_)                              !^CMP IF SC
          call couple_ee_sc(TimeSimulation)   !^CMP IF SC
       case default
          call error
       end select                             !^CMP END EE
    case(SC_)                                 !^CMP IF SC BEGIN
       select case(iCompTarget)
       case(GM_)                              !^CMP IF GM
          call couple_sc_gm(TimeSimulation)   !^CMP IF GM
       case(IH_)                              !^CMP IF IH
          call couple_sc_ih(TimeSimulation)   !^CMP IF IH
       case(SP_)                              !^CMP IF SP
          call couple_sc_sp(TimeSimulation)   !^CMP IF SP
       case(EE_)                              !^CMP IF EE
          call couple_sc_ee(TimeSimulation)   !^CMP IF EE
       case(PT_)                              !^CMP IF PT BEGIN
          if(UseBLine_C(PT_))then
             call couple_sc_sp(TimeSimulation)
          else
             call couple_sc_pt(TimeSimulation)
          end if                              !^CMP END PT
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
       case(PT_)                              !^CMP IF PT BEGIN
          if(UseBLine_C(PT_))then
             call couple_ih_sp(TimeSimulation)
          else
             call couple_ih_pt(TimeSimulation)
          end if                              !^CMP END PT
       case(SP_)                                   !^CMP IF SP
          call couple_ih_sp(TimeSimulation)        !^CMP IF SP
       case default
          call error
       end select                             !^CMP END IH
    case(OH_)
       select case(iCompTarget)               !^CMP IF OH BEGIN
       case(IH_)                                   !^CMP IF IH
          call couple_oh_ih(TimeSimulation)        !^CMP IF IH
       case(SP_)                                   !^CMP IF SP
          call couple_oh_sp(TimeSimulation)        !^CMP IF SP
       case(PT_)				   !^CMP IF PT
          if(UseBLine_C(PT_))then
             call couple_oh_sp(TimeSimulation)
          else
             call couple_oh_pt(TimeSimulation)
          end if                                   !^CMP IF PT
       case default
          call error
       end select                                  !^CMP END OH
    case(GM_)                                 !^CMP IF GM BEGIN
       select case(iCompTarget)
       case(IE_)                                   !^CMP IF IE
          call couple_gm_ie(TimeSimulation)        !^CMP IF IE
       case(IM_)                                   !^CMP IF IM
          call couple_gm_im(TimeSimulation)        !^CMP IF IM
       case(PC_)                                   !^CMP IF PC
          call couple_gm_pc(TimeSimulation)        !^CMP IF PC
       case(PT_)                                   !^CMP IF PT
          call couple_gm_pt(TimeSimulation)        !^CMP IF PT
       case(PW_)                                   !^CMP IF PW
          call couple_gm_pw(TimeSimulation)        !^CMP IF PW
       case(RB_)                                   !^CMP IF RB
          call couple_gm_rb(TimeSimulation)        !^CMP IF RB
       case(SC_)                                   !^CMP IF SC
          call couple_gm_sc(TimeSimulation)        !^CMP IF SC
       case default
          call error
       end select                             !^CMP END GM
    case(UA_)                                 !^CMP IF UA BEGIN
       select case(iCompTarget)
       case(GM_)                                   !^CMP IF GM
          call couple_ua_gm(TimeSimulation)        !^CMP IF GM
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
       case(PS_)                                   !^CMP IF PS
          call couple_ie_ps(TimeSimulation)        !^CMP IF PS
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
       case(IE_)                                   !^CMP IF IE
          call couple_im_ie(TimeSimulation)        !^CMP IF IE
       case(GM_)                                   !^CMP IF GM
          call couple_im_gm(TimeSimulation)        !^CMP IF GM
       case default
          call error
       end select                             !^CMP END IM
    case(PS_)                                 !^CMP IF PS BEGIN
       select case(iCompTarget)
       case(GM_)                                   !^CMP IF GM
          call couple_ps_gm(TimeSimulation)        !^CMP IF GM
       case default
          call error
       end select                             !^CMP END PS
    case(PC_)                                 !^CMP IF PC BEGIN
       select case(iCompTarget)
       case(GM_)                                   !^CMP IF GM
          call couple_pc_gm(TimeSimulation)        !^CMP IF GM
       case default
          call error
       end select                             !^CMP END PC
    case(PW_)                                 !^CMP IF PW BEGIN
       select case(iCompTarget)
       case(GM_)                                   !^CMP IF GM
          call couple_pw_gm(TimeSimulation)        !^CMP IF GM
       case default
          call error
       end select                             !^CMP END PW
    case(PT_)
       select case(iCompTarget)               !^CMP IF PT BEGIN
       case(OH_)                                   !^CMP IF OH
          call couple_pt_oh(TimeSimulation)        !^CMP IF OH
       case default
          call error
       end select                             !^CMP END PT
    case default
       call CON_stop(NameSub//&
            ' SWMF_ERROR: no coupling implemented from source '// &
            NameComp_I(iCompSource))
    end select
    call timing_stop(NameComp_I(iCompSource)// &
         '_'//NameComp_I(iCompTarget)//'_couple')

    if(DoTest)write(*,*)NameSub,': finished coupling iProc=',i_proc(),' ',&
         NameComp_I(iCompSource),' --> ',NameComp_I(iCompTarget)

  contains
    !==========================================================================

    subroutine error
      !------------------------------------------------------------------------
      call CON_stop(NameSub//' SWMF_ERROR: no coupling implemented for '//&
           NameComp_I(iCompSource)//' --> '//&
           NameComp_I(iCompTarget))
    end subroutine error
    !==========================================================================

  end subroutine couple_two_comp
  !============================================================================

end module CON_couple_all
!==============================================================================
