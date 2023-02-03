!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE SC
!^CMP FILE PT

!
! Couple SC and PT components: SC -> PT.
!
module CON_couple_sc_pt

  use CON_coupler

  use CON_couple_points, ONLY: &
       couple_points_init, couple_points, CouplePointsType
  use SC_wrapper, ONLY: &
       SC_get_grid_info, SC_find_points, SC_get_for_pt,  &
       SC_get_for_pt_dt
  use PT_wrapper, ONLY: &
       PT_get_grid_info, PT_find_points, PT_put_from_sc, &
       PT_put_from_sc_dt,PT_divu_coupling_state

  implicit none
  save

  private ! except

  public :: couple_sc_pt_init ! initialize both couplings
  public :: couple_sc_pt      ! couple SC to PT

  ! revision history:
  ! 03/16/2015 A.Michael and G.Toth - initial version

  ! Router communicator info
  type(CouplePointsType) :: CouplerScToPt

contains
  !============================================================================

  subroutine couple_sc_pt_init

    ! Initialize SC->PT coupler.
    ! This subroutine should be called from all PE-s

    logical :: DoTest, DoTestMe,flag

    character(len=*), parameter:: NameSub = 'couple_sc_pt_init'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTestMe)write(*,*) NameSub,' starting'

    CouplerScToPt%iCompTarget = PT_
    CouplerScToPt%iCompSource = SC_

    ! SC sends all its variables to PT. Take information from Grid_C
    CouplerScToPt%NameVar = Grid_C(SC_)%NameVar
    CouplerScToPt%nVar    = Grid_C(SC_)%nVar

    ! Add divu to the list of communicated variables if needed
    call PT_divu_coupling_state(flag)

    if (flag) then
      CouplerScToPt%nVar=CouplerScToPt%nVar+1
      CouplerScToPt%NameVar=trim(CouplerScToPt%NameVar)//" divu"
    endif

    call couple_points_init(CouplerScToPt)

    if(DoTestMe)write(*,*) NameSub,' finished'

  end subroutine couple_sc_pt_init
  !============================================================================

  subroutine couple_sc_pt(tSimulation)
    use CON_transfer_data, ONLY: transfer_real
    real, intent(in) :: tSimulation

    ! Couple between two components:
    !    Solar Corona               (SC) source
    !    Particle Tracker           (PT) target
    !
    ! Send information from SC to PT.

    logical :: DoTest, DoTestMe
    real:: DtSi ! time step in SI units

    character(len=*), parameter:: NameSub = 'couple_sc_pt'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTest)write(*,*) NameSub,' starting iProc=', CouplerScToPt%iProcWorld

    if(IsTightCouple_CC(SC_,PT_)) then
       if(is_proc(SC_)) call SC_get_for_pt_dt(DtSi)
       call transfer_real(SC_,PT_,DtSi)
       if(is_proc(PT_)) call PT_put_from_sc_dt(DtSi)
    endif

    call couple_points(CouplerScToPt, SC_get_grid_info, SC_find_points, &
         SC_get_for_pt, PT_get_grid_info, PT_put_from_sc)

    if(DoTest)write(*,*) NameSub,' finished iProc=', CouplerScToPt%iProcWorld

  end subroutine couple_sc_pt
  !============================================================================

end module CON_couple_sc_pt
!==============================================================================
