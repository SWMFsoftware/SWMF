! !  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
! !  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE EE
!^CMP FILE SC
!BOP
!MODULE: CON_couple_ee_sc - couple EE and SC components
!
!DESCRIPTION:
! Couple EE and SC components both ways
!
!INTERFACE:
module CON_couple_ee_sc

  !USES:
  use CON_coupler

  implicit none

  private !except

  !PUBLIC MEMBER FUNCTIONS:

  public:: couple_ee_sc_init
  public:: couple_ee_sc, couple_sc_ee

  !REVISION HISTORY:
  !EOP

contains
  !============================================================================
  subroutine couple_ee_sc_init
    !--------------------------------------------------------------------------
  end subroutine couple_ee_sc_init

  !============================================================================
  subroutine couple_ee_sc(tSimulation)

    real, intent(in) :: tSimulation ! simulation time at coupling
    !--------------------------------------------------------------------------
  end subroutine couple_ee_sc

  !============================================================================
  subroutine couple_sc_ee(tSimulation)

    real, intent(in) :: tSimulation ! simulation time at coupling
    !--------------------------------------------------------------------------
  end subroutine couple_sc_ee

end module CON_couple_ee_sc

