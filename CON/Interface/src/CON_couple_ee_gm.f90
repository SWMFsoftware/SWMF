!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE EE
!^CMP FILE GM

!
! Couple EE and GM components both ways. Here EE is just an arbitrary BATS-R-US
! model that runs on the same grid as the GM model.
! The goal is to solve different sets of equations with differnt time steps
! and pass the State_VGB pointers back-and-force so that various source terms
! can be calculated.
!
module CON_couple_ee_gm

  use CON_coupler

  use EE_wrapper, ONLY: EE_use_pointer
  use GM_wrapper, ONLY: GM_use_pointer

  implicit none
  save

  private ! except

  public :: couple_ee_gm_init ! initialize both couplings
  public :: couple_ee_gm      ! couple EE to GM
  public :: couple_gm_ee      ! couple GM to EE

  ! revision history:
  ! 05/09/2016 G.Toth snf Y. Shou - initial version

contains
  !============================================================================
  subroutine couple_ee_gm_init

    ! This subroutine should be called from all PE-s
    !--------------------------------------------------------------------------
    if(is_proc(GM_) .neqv. is_proc(EE_)) &
         call CON_stop('EE and GM should run on same processors')

  end subroutine couple_ee_gm_init
  !============================================================================
  subroutine couple_ee_gm(tSimulation)

    real, intent(in) :: tSimulation
    !--------------------------------------------------------------------------
    call GM_use_pointer(EE_, tSimulation)

  end subroutine couple_ee_gm
  !============================================================================
  subroutine couple_gm_ee(tSimulation)

    ! List of variables to pass
    real, intent(in) :: tSimulation
    !--------------------------------------------------------------------------
    call EE_use_pointer(GM_, tSimulation)

  end subroutine couple_gm_ee
  !============================================================================
end module CON_couple_ee_gm
!==============================================================================
