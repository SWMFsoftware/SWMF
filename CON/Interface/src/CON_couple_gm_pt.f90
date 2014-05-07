!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE GM
!^CMP FILE PT

!BOP
!MODULE: CON_couple_gm_pt - couple GM and PT components
!
!DESCRIPTION:
! Couple GM and PT components both ways. 
!
!INTERFACE:
module CON_couple_gm_pt

  !USES:
  use CON_coupler, ONLY: GM_, PT_, Grid_C, lNameVar
  use CON_couple_points, ONLY: &
       couple_points_init, couple_points, CouplePointsType

  use GM_wrapper, ONLY: GM_get_grid_info, GM_find_points, GM_get_for_pt
  use PT_wrapper, ONLY: PT_get_grid_info, PT_put_from_gm

  implicit none
  save

  private ! except

  !PUBLIC MEMBER FUNCTIONS:

  public :: couple_gm_pt_init ! initialize both couplings
  public :: couple_gm_pt      ! couple GM to PT

  !REVISION HISTORY:
  ! 09/24/2013 G.Toth <gtoth@umich.edu> - initial version
  !EOP

  type(CouplePointsType) :: Coupler
  
contains

  !BOP =======================================================================
  !IROUTINE: couple_gm_pt_init - initialize GM-PT couplings
  !INTERFACE:
  subroutine couple_gm_pt_init

    !DESCRIPTION:
    ! Initialize GM->PT coupler.
    ! This subroutine should be called from all PE-s
    !EOP

    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'couple_gm_pt_init'
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTestMe)write(*,*) NameSub,' starting'

    Coupler%iCompTarget = PT_
    Coupler%iCompSource = GM_
    
    ! GM sends all its variables to PT. Take information from Grid_C
    Coupler%NameVar = Grid_C(Coupler%iCompSource)%NameVar
    Coupler%nVar    = Grid_C(Coupler%iCompSource)%nVar

    call couple_points_init(Coupler)

    if(DoTestMe)write(*,*) NameSub,' finished'

  end subroutine couple_gm_pt_init

  !BOP =======================================================================
  !IROUTINE: couple_gm_pt - couple GM to PT
  !INTERFACE:
  subroutine couple_gm_pt(tSimulation)
    
    !INPUT ARGUMENT:
    real, intent(in) :: tSimulation

    !DESCRIPTION:
    ! Couple between two components:\\
    !    Global Magnetosphere       (GM) source\\
    !    Particle Tracker           (PT) target
    !
    ! Send information from GM to PT. 
    !EOP

    logical :: DoTest, DoTestMe
    character (len=*), parameter:: NameSub = 'couple_gm_pt'
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTest)write(*,*) NameSub,' starting iProc=', Coupler%iProcWorld

    call couple_points(Coupler, GM_get_grid_info, GM_find_points, &
         GM_get_for_pt, PT_get_grid_info, PT_put_from_gm)

    if(DoTest)write(*,*) NameSub,' finished iProc=', Coupler%iProcWorld

  end subroutine couple_gm_pt

end module CON_couple_gm_pt
