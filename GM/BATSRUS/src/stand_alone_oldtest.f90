!^CFG COPYRIGHT UM
subroutine stand_alone

  ! Set defaults for the standalone code.
  ! This version behaves like BATSRUS before the merge with old directories.

  use ModProcMH, ONLY: iProc
  use ModMain,   ONLY: NameThisComp, UseNewParam, UseNewAxes, UseRotatingBc,&
       time_accurate

  implicit none

  NameThisComp  = 'GM'
  UseNewParam   = .false.
  UseNewAxes    = .false.
  time_accurate = .false.
  UseRotatingBc = .false.

end subroutine stand_alone
