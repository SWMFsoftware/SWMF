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
       OH_get_grid_info, OH_find_points, OH_get_for_pt, OH_put_from_pt
  use PT_wrapper, ONLY: &
       PT_get_grid_info, PT_find_points, PT_get_for_oh, PT_put_from_oh

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

    integer :: nFluid

    integer :: l, i, i0

    logical :: DoTest, DoTestMe

    character(len=*), parameter:: NameSub = 'couple_oh_pt_init'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTestMe)write(*,*) NameSub,' starting'

    CouplerOhToPt%iCompTarget = PT_
    CouplerOhToPt%iCompSource = OH_

    ! Take information from both Grid_C(OH_) and Grid_C(PT_)
    CouplerOhToPt%NameVar = trim(Grid_C(OH_)%NameVar)//' '//Grid_C(PT_)%NameVar
    CouplerOhToPt%nVar    = Grid_C(OH_)%nVar + Grid_C(PT_)%nVar

    call couple_points_init(CouplerOhToPt)

    CouplerPtToOh%iCompSource = PT_
    CouplerPtToOh%iCompTarget = OH_

    ! Find the number of fluids from MHD variable names
    l = len(Grid_C(OH_)%NameVar)
    nFluid = 0
    i = 1
    i0 = 1
    do while(i0 < l .and. i > 0)
       ! Assume the fluid density variable name is '*rho*'
       i = index(Grid_C(OH_)%NameVar(i0:), 'rho')
       if(i > 0) then
          i0 = i0 + i + 3
          nFluid = nFluid + 1
       end if
    end do

    ! 5 source terms per fluid
    CouplerPtToOh%nVar    = 5 * nFluid

    ! charge exchange sources
    if(nFluid == 1)then
       CouplerPtToOh%NameVar = 'Srho Smx Smy Smz Se'
    else if(nFluid == 2) then
       CouplerPtToOh%NameVar = 'Srho Smx Smy Smz Se Srho2 Smx2 Smy2 Smz2 Se2'
    else if(nFluid == 3) then
       CouplerPtToOh%NameVar = 'Srho Smx Smy Smz Se Srho2 Smx2 Smy2 Smz2 Se2 '&
            //'Srho3 Smx3 Smy3 Smz3 Se3'
    else
       write(*,*)NameSub, ' Error!'
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
