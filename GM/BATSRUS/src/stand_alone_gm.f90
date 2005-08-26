!^CFG COPYRIGHT UM
subroutine stand_alone

  ! Set defaults for the stand alone code.
  ! This version behaves like the GM component of SWMF.

  use ModProcMH, ONLY: iProc
  use ModMain,   ONLY: NameThisComp, UseRotatingBc
  use ModIO,     ONLY: NamePlotDir, NameRestartInDir, NameRestartOutDir

  implicit none

  NameThisComp  = 'GM'
  UseRotatingBc = .true.

  NamePlotDir       = 'GM/'//NamePlotDir
  NameRestartInDir  = 'GM/'//NameRestartInDir
  NameRestartOutDir = 'GM/'//NameRestartOutDir

end subroutine stand_alone
