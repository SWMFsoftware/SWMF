!^CFG COPYRIGHT UM
subroutine stand_alone

  ! Set defaults for the standalone code.
  ! This version behaves like BATSRUS before the merge but uses new directories

  use ModProcMH, ONLY: iProc
  use ModMain,   ONLY: NameThisComp, UseNewParam, UseNewAxes, UseRotatingBc,&
       time_accurate
  use ModIO,     ONLY: NamePlotDir, NameRestartInDir, NameRestartOutDir

  implicit none

  NameThisComp  = 'GM'
  UseNewParam   = .false.
  UseNewAxes    = .false.
  time_accurate = .false.
  UseRotatingBc = .false.

  NamePlotDir       = 'GM/'//NamePlotDir
  NameRestartInDir  = 'GM/'//NameRestartInDir
  NameRestartOutDir = 'GM/'//NameRestartOutDir

end subroutine stand_alone
