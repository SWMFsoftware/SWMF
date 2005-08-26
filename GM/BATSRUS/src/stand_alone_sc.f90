!^CFG COPYRIGHT UM
subroutine stand_alone

  ! Set defaults for the stand alone code.
  ! This version behaves like the SC component of SWMF.

  use ModProcMH, ONLY: iProc
  use ModMain,   ONLY: NameThisComp, UseRotatingBc
  use ModIO,     ONLY: NamePlotDir, NameRestartInDir, NameRestartOutDir

  implicit none

  NameThisComp  = 'SC'
  UseRotatingBc = .false.

  NamePlotDir       = 'SC/'//NamePlotDir
  NameRestartInDir  = 'SC/'//NameRestartInDir
  NameRestartOutDir = 'SC/'//NameRestartOutDir

end subroutine stand_alone
