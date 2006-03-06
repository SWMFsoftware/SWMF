!^CFG COPYRIGHT UM
subroutine stand_alone

  ! Set defaults for the stand alone code.
  ! This version behaves like the IH component of SWMF.

  use ModMain,   ONLY: NameThisComp, UseRotatingBc

  implicit none

  NameThisComp  = 'IH'
  UseRotatingBc = .true.

end subroutine stand_alone
