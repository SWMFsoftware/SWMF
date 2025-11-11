!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE UA
!^CMP FILE GM

module CON_couple_gm_ua

  ! UA provides source terms for GM where possible

  use CON_coupler

  use CON_couple_points

  use GM_wrapper, ONLY: GM_find_points, GM_put_from_ua, GM_get_grid_info
  use UA_wrapper, ONLY: UA_find_points, UA_get_for_gm, UA_get_grid_info

  implicit none
  save

  private ! except

  public :: couple_ua_gm_init
  public :: couple_ua_gm      ! couple UA to GM

  ! Router communicator info
  type(CouplePointsType) :: CouplerUAtoGM

contains
  !============================================================================
  subroutine couple_ua_gm_init

    logical :: DoTest, DoTestMe

    ! This subroutine should be called from all PE-s

    character(len=*), parameter:: NameSub = 'couple_ua_gm_init'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    ! Initialize point couplers

    CouplerUAtoGM%iCompSource = UA_
    CouplerUAtoGM%iCompTarget = GM_
    CouplerUAtoGM%nVar        = 5

    call couple_points_init(CouplerUAtoGM)

  end subroutine couple_ua_gm_init
  !============================================================================
  subroutine couple_ua_gm(tSimulation)

    real, intent(in) :: tSimulation

    logical :: DoTest, DoTestMe

    character(len=*), parameter:: NameSub = 'couple_ua_gm'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTest) write(*,*) NameSub,' starting iProc=', CouplerUAtoGM%iProcWorld

    call couple_points(CouplerUAtoGM, UA_get_grid_info, UA_find_points, &
         UA_get_for_GM, GM_get_grid_info, GM_put_from_UA)

    if(DoTest) write(*,*) NameSub,' finished, iProc=', i_proc()

  end subroutine couple_ua_gm
  !============================================================================
end module CON_couple_gm_ua
!==============================================================================
