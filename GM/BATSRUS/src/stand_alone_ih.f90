!^CFG COPYRIGHT UM
subroutine stand_alone

  ! Set defaults for the stand alone code.
  ! This version behaves like the IH component of SWMF.

  use ModProcMH, ONLY: iProc
  use ModMain,   ONLY: NameThisComp, UseNewParam, UseNewAxes, UseCorotation,&
       time_accurate
  use ModIO,     ONLY: NamePlotDir, NameRestartInDir, NameRestartOutDir

  implicit none

  NameThisComp  = 'IH'
  UseNewParam   = .true.
  UseNewAxes    = .true.
  time_accurate = .true.
  UseCorotation = .true.

  NamePlotDir       = 'IH/'//NamePlotDir
  NameRestartInDir  = 'IH/'//NameRestartInDir
  NameRestartOutDir = 'IH/'//NameRestartOutDir

end subroutine stand_alone
