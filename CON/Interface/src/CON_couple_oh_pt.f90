!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE OH
!^CMP FILE PT

!BOP
!MODULE: CON_couple_oh_pt - couple OH and PT components
!
!DESCRIPTION:
! Couple OH and PT components both ways. 
!
!INTERFACE:
module CON_couple_oh_pt

  !USES:
  use CON_coupler

  use CON_couple_points, ONLY: &
       couple_points_init, couple_points, CouplePointsType
  use OH_wrapper, ONLY: &
       OH_get_grid_info, OH_find_points, OH_get_for_pt, OH_put_from_pt, &
       OH_get_for_pt_dt
  use PT_wrapper, ONLY: &
       PT_get_grid_info, PT_find_points, PT_get_for_oh, PT_put_from_oh, &
       PT_put_from_oh_dt

  implicit none
  save

  private ! except

  !PUBLIC MEMBER FUNCTIONS:

  public :: couple_oh_pt_init ! initialize both couplings
  public :: couple_oh_pt      ! couple OH to PT
  public :: couple_pt_oh      ! couple PT to OH

  !REVISION HISTORY:
  ! 03/16/2015 A.Michael and G.Toth - initial version
  !EOP

  ! Router communicator info
  type(CouplePointsType) :: CouplerOhToPt, CouplerPtToOh

contains

  !BOP =======================================================================
  !IROUTINE: couple_oh_pt_init - initialize OH-PT couplings
  !INTERFACE:
  subroutine couple_oh_pt_init

    !DESCRIPTION:
    ! Initialize OH->PT coupler.
    ! This subroutine should be called from all PE-s
    !EOP

    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'couple_oh_pt_init'
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTestMe)write(*,*) NameSub,' starting'

    CouplerOhToPt%iCompTarget = PT_
    CouplerOhToPt%iCompSource = OH_
    
    ! OH sends all its variables to PT. Take information from Grid_C
    CouplerOhToPt%NameVar = Grid_C(OH_)%NameVar
    CouplerOhToPt%nVar    = Grid_C(OH_)%nVar

    call couple_points_init(CouplerOhToPt)

    CouplerPtToOh%iCompSource = PT_
    CouplerPtToOh%iCompTarget = OH_
    CouplerPtToOh%NameVar     = 'Srho Smx Smy Smz Se' ! charge exchange sources
    CouplerPtToOh%nVar        = 5

    call couple_points_init(CouplerPtToOh)

    if(DoTestMe)write(*,*) NameSub,' finished'

  end subroutine couple_oh_pt_init

  !BOP =======================================================================
  !IROUTINE: couple_oh_pt - couple OH to PT
  !INTERFACE:
  subroutine couple_oh_pt(tSimulation)
    use CON_transfer_data, ONLY: transfer_real
    !INPUT ARGUMENT:
    real, intent(in) :: tSimulation

    !DESCRIPTION:
    ! Couple between two components:\\
    !    Ourer Heliosphere          (OH) source\\
    !    Particle Tracker           (PT) target
    !
    ! Send information from OH to PT. 
    !EOP

    logical :: DoTest, DoTestMe
    real:: DtSi ! time step in SI units
    character (len=*), parameter:: NameSub = 'couple_oh_pt'
    !-------------------------------------------------------------------------
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

  !=======================================================================
  subroutine couple_pt_oh(tSimulation)

    ! List of variables to pass
    real, intent(in) :: tSimulation

    ! The first time there is no back coupling
    ! This is actually wrong for restarts!!!
    logical:: IsFirstTime = .true.

    logical :: DoTest, DoTestMe
    character (len=*), parameter :: NameSub='couple_pt_oh'
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTest)write(*,*)NameSub,' starting iProc=',CouplerPTtoOH%iProcWorld

    call couple_points(CouplerPtToOh, PT_get_grid_info,  PT_find_points , &
         PT_get_for_oh, OH_get_grid_info, OH_put_from_pt)

    if(DoTest) write(*,*) NameSub,' finished, iProc=', i_proc()

  end subroutine couple_pt_oh

end module CON_couple_oh_pt
