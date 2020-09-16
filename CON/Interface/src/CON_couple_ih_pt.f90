!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE IH
!^CMP FILE PT

!BOP
!MODULE: CON_couple_ih_pt - couple IH and PT components
!
!DESCRIPTION:
! Couple IH and PT components: IH -> PT. 
!
!INTERFACE:
module CON_couple_ih_pt

  !USES:
  use CON_coupler

  use CON_couple_points, ONLY: &
       couple_points_init, couple_points, CouplePointsType
  use IH_wrapper, ONLY: &
       IH_get_grid_info, IH_find_points, IH_get_for_pt,  &
       IH_get_for_pt_dt
  use PT_wrapper, ONLY: &
       PT_get_grid_info, PT_find_points, PT_put_from_ih, &
       PT_put_from_ih_dt

  implicit none
  save

  private ! except

  !PUBLIC MEMBER FUNCTIONS:

  public :: couple_ih_pt_init ! initialize both couplings
  public :: couple_ih_pt      ! couple IH to PT

  !REVISION HISTORY:
  ! 03/16/2015 A.Michael and G.Toth - initial version
  !EOP

  ! Router communicator info
  type(CouplePointsType) :: CouplerIhToPt

contains

  !BOP =======================================================================
  !IROUTINE: couple_ih_pt_init - initialize IH-PT couplings
  !INTERFACE:
  subroutine couple_ih_pt_init

    !DESCRIPTION:
    ! Initialize IH->PT coupler.
    ! This subroutine should be called from all PE-s
    !EOP

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub = 'couple_ih_pt_init'
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTestMe)write(*,*) NameSub,' starting'

    CouplerIhToPt%iCompTarget = PT_
    CouplerIhToPt%iCompSource = IH_
    
    ! IH sends all its variables to PT. Take information from Grid_C
    CouplerIhToPt%NameVar = Grid_C(IH_)%NameVar
    CouplerIhToPt%nVar    = Grid_C(IH_)%nVar

    call couple_points_init(CouplerIhToPt)

    if(DoTestMe)write(*,*) NameSub,' finished'

  end subroutine couple_ih_pt_init

  !BOP =======================================================================
  !IROUTINE: couple_ih_pt - couple IH to PT
  !INTERFACE:
  subroutine couple_ih_pt(tSimulation)
    use CON_transfer_data, ONLY: transfer_real
    !INPUT ARGUMENT:
    real, intent(in) :: tSimulation

    !DESCRIPTION:
    ! Couple between two components:\\
    !    Inner Heliosphere          (IH) source\\
    !    Particle Tracker           (PT) target
    !
    ! Send information from IH to PT. 
    !EOP

    logical :: DoTest, DoTestMe
    real:: DtSi ! time step in SI units
    character(len=*), parameter :: NameSub = 'couple_ih_pt'
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTest)write(*,*) NameSub,' starting iProc=', CouplerIhToPt%iProcWorld

    if(IsTightCouple_CC(IH_,PT_)) then 
       if(is_proc(IH_)) call IH_get_for_pt_dt(DtSi)
       call transfer_real(IH_,PT_,DtSi)
       if(is_proc(PT_)) call PT_put_from_ih_dt(DtSi)
    endif

    call couple_points(CouplerIhToPt, IH_get_grid_info, IH_find_points, &
         IH_get_for_pt, PT_get_grid_info, PT_put_from_ih)

    if(DoTest)write(*,*) NameSub,' finished iProc=', CouplerIhToPt%iProcWorld

  end subroutine couple_ih_pt

  !=======================================================================
end module CON_couple_ih_pt
