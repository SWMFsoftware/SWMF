!^CFG COPYRIGHT UM
subroutine stand_alone

  ! Set defaults for the stand alone code.
  ! This version behaves like the IH component of SWMF.

  use ModProcMH, ONLY: iProc
  use ModMain,   ONLY: NameThisComp, UseNewParam, UseNewAxes, UseCorotation,&
       time_accurate
  use ModIO,     ONLY: NamePlotDir, NameRestartInDir, NameRestartOutDir

  implicit none

  NameThisComp  = 'SC'
  UseNewParam   = .true.
  UseNewAxes    = .false.
  time_accurate = .false.
  UseCorotation = .false.

  NamePlotDir       = 'SC/'//NamePlotDir
  NameRestartInDir  = 'SC/'//NameRestartInDir
  NameRestartOutDir = 'SC/'//NameRestartOutDir

end subroutine stand_alone
