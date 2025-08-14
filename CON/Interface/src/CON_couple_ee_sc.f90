!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE EE
!^CMP FILE SC

! Couple EE and SC components both ways.
! EE overwrites SC in the whole EE domain. Extra SC variables are solved
! for in SC.
! SC provides boundary conditions for EE where possible.

module CON_couple_ee_sc

  use CON_coupler

  use CON_couple_points

  use EE_wrapper
  use SC_wrapper

  implicit none
  save

  private ! except

  public :: couple_ee_sc_init ! initialize both couplings
  public :: couple_ee_sc      ! couple EE to SC
  public :: couple_sc_ee      ! couple SC to EE

  ! revision history:
  ! 09/24/2013 G.Toth <gtoth@umich.edu> - initial version

  ! Router communicator info
  type(CouplePointsType) :: CouplerEEtoSC,  CouplerSCtoEE

contains
  !============================================================================
  subroutine couple_ee_sc_init

    logical :: DoTest, DoTestMe

    ! This subroutine should be called from all PE-s

    character(len=*), parameter:: NameSub = 'couple_ee_sc_init'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    ! Set up variable coupling information
    call set_couple_var_info(EE_, SC_)

    ! Initialize point couplers

    CouplerEEtoSC%iCompSource = EE_
    CouplerEEtoSC%iCompTarget = SC_
    CouplerEEtoSC%NameVar     = NameVarBuffer
    CouplerEEtoSC%nVar        = nVarBuffer

    call couple_points_init(CouplerEEtoSC)

    ! Set up variable coupling information
    call set_couple_var_info(SC_, EE_)

    CouplerSCtoEE%iCompSource = SC_
    CouplerSCtoEE%iCompTarget = EE_
    CouplerSCtoEE%NameVar     = NameVarBuffer
    CouplerSCtoEE%nVar        = nVarBuffer

    call couple_points_init(CouplerSCtoEE)

  end subroutine couple_ee_sc_init
  !============================================================================
  subroutine couple_ee_sc(tSimulation)

    real, intent(in) :: tSimulation

    logical :: DoTest, DoTestMe

    ! Name of this interface
    character(len=*), parameter:: NameSub = 'couple_ee_sc'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTest)write(*,*)NameSub,' starting iProc=', i_proc()

    call couple_points(CouplerEEtoSC, EE_get_grid_info, EE_find_points, &
         EE_get_for_sc, SC_get_grid_info, SC_put_from_ee)

    if(DoTest) write(*,*) NameSub,' finished, iProc=', i_proc()

  end subroutine couple_ee_sc
  !============================================================================
  subroutine couple_sc_ee(tSimulation)

    ! List of variables to pass
    real, intent(in) :: tSimulation

    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'couple_sc_ee'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTest)write(*,*)NameSub,' starting iProc=',CouplerSCtoEE%iProcWorld

    call couple_points(CouplerSCtoEE, SC_get_grid_info,  SC_find_points , &
         SC_get_for_ee, EE_get_grid_info, EE_put_from_sc)

    if(DoTest) write(*,*) NameSub,' finished, iProc=', i_proc()

  end subroutine couple_sc_ee
  !============================================================================
end module CON_couple_ee_sc
!==============================================================================
