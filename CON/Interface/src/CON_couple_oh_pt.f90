!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE OH
!^CMP FILE PT

!
! Couple OH and PT components both ways.
!
module CON_couple_oh_pt

  use CON_coupler

  use CON_couple_points, ONLY: &
       couple_points_init, couple_points, CouplePointsType
  use OH_wrapper, ONLY: &
       OH_get_grid_info, OH_find_points, OH_get_for_pt, OH_put_from_pt, &
       OH_get_for_pt_dt
  use PT_wrapper, ONLY: &
       PT_get_grid_info, PT_find_points, PT_get_for_oh, PT_put_from_oh, &
       PT_put_from_oh_dt,PT_divu_coupling_state

  implicit none
  save

  private ! except

  public :: couple_oh_pt_init ! initialize both couplings
  public :: couple_oh_pt      ! couple OH to PT
  public :: couple_pt_oh      ! couple PT to OH

  ! revision history:
  ! 03/16/2015 A.Michael and G.Toth - initial version

  ! Router communicator info
  type(CouplePointsType) :: CouplerOhToPt, CouplerPtToOh

contains
  !============================================================================

  subroutine couple_oh_pt_init

    ! Initialize OH->PT coupler.
    ! This subroutine should be called from all PE-s

    logical :: DoTest, DoTestMe,flag

    character(len=*), parameter:: NameSub = 'couple_oh_pt_init'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTestMe)write(*,*) NameSub,' starting'

    CouplerOhToPt%iCompTarget = PT_
    CouplerOhToPt%iCompSource = OH_

    ! OH sends all its variables to PT. Take information from Grid_C
    CouplerOhToPt%NameVar = Grid_C(OH_)%NameVar
    CouplerOhToPt%nVar    = Grid_C(OH_)%nVar

    ! Add divu to the list of communicated variables if needed
    call PT_divu_coupling_state(flag)

    if (flag) then
      CouplerOhToPt%nVar=CouplerOhToPt%nVar+1
      CouplerOhToPt%NameVar=trim(CouplerOhToPt%NameVar)//" divu"
    endif

    call couple_points_init(CouplerOhToPt)

    CouplerPtToOh%iCompSource = PT_
    CouplerPtToOh%iCompTarget = OH_

    ! 5 source terms per fluid. Number of fluids is from number of MHD vars.
    CouplerPtToOh%nVar    = 5 * (CouplerOhToPt%nVar / 5)

    ! charge exchange sources
    if(CouplerPtToOh%nVar == 5)then
       CouplerPtToOh%NameVar = 'Srho Smx Smy Smz Se'
    else
       CouplerPtToOh%NameVar = 'Srho Smx Smy Smz Se Srho2 Smx2 Smy2 Smz2 Se2'
    end if

    call couple_points_init(CouplerPtToOh)

    if(DoTestMe)write(*,*) NameSub,' finished'

  end subroutine couple_oh_pt_init
  !============================================================================

  subroutine couple_oh_pt(tSimulation)
    use CON_transfer_data, ONLY: transfer_real
    real, intent(in) :: tSimulation

    ! Couple between two components:
    !    Ourer Heliosphere          (OH) source
    !    Particle Tracker           (PT) target
    !
    ! Send information from OH to PT.

    logical :: DoTest, DoTestMe
    real:: DtSi ! time step in SI units

    character(len=*), parameter:: NameSub = 'couple_oh_pt'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTest)write(*,*) NameSub,' starting iProc=', CouplerOhToPt%iProcWorld

    if(IsTightCouple_CC(OH_,PT_)) then
       if(is_proc(OH_)) call OH_get_for_pt_dt(DtSi)
       call transfer_real(OH_,PT_,DtSi)
       if(is_proc(PT_)) call PT_put_from_oh_dt(DtSi)
    endif

    call couple_points(CouplerOhToPt, OH_get_grid_info, OH_find_points, &
         OH_get_for_pt, PT_get_grid_info, PT_put_from_oh)

    if(DoTest)write(*,*) NameSub,' finished iProc=', CouplerOhToPt%iProcWorld

  end subroutine couple_oh_pt
  !============================================================================

  subroutine couple_pt_oh(tSimulation)

    ! List of variables to pass
    real, intent(in) :: tSimulation

    ! The first time there is no back coupling
    ! This is actually wrong for restarts!!!
    logical:: IsFirstTime = .true.

    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'couple_pt_oh'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTest)write(*,*)NameSub,' starting iProc=',CouplerPTtoOH%iProcWorld

    call couple_points(CouplerPtToOh, PT_get_grid_info,  PT_find_points , &
         PT_get_for_oh, OH_get_grid_info, OH_put_from_pt)

    if(DoTest) write(*,*) NameSub,' finished, iProc=', i_proc()

  end subroutine couple_pt_oh
  !============================================================================

end module CON_couple_oh_pt
!==============================================================================
