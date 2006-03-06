!^CFG COPYRIGHT UM
subroutine stand_alone

  ! Set defaults for the stand alone code.
  ! This version behaves like the SC component of SWMF.

  use ModMain,   ONLY: NameThisComp, UseRotatingBc

  implicit none

  NameThisComp  = 'SC'
  UseRotatingBc = .false.

end subroutine stand_alone
