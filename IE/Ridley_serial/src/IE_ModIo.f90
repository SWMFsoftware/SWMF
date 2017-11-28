!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module IE_ModIo

  ! Unit number for temporary files
  use ModIoUnit, ONLY: iUnit => UNITTMP_, STDOUT_
  use ModUtilities, ONLY: open_file, close_file, check_dir, make_dir

  implicit none
  save

  integer            :: iUnitOut=STDOUT_, unitlog=-1
  integer, parameter :: lStringPrefix=6
  character (len=lStringPrefix) :: StringPrefix = ''

  ! Maximum number of plot files
  integer, parameter :: MaxFile = 10

  ! Actual number of plot files
  integer :: nFile=0

  character(LEN=3)   :: plot_form(MaxFile)
  character(LEN=100) :: plot_vars(MaxFile)
  character(len=100) :: NameIonoDir='ionosphere/'

  logical            :: DoRestart = .false.
  character(len=100) :: NameRestartInDir  = 'IE/restartIN/'
  character(len=100) :: NameRestartOutDir = 'IE/restartOUT/'

  character(len=100) :: NameFile

  integer            :: dn_output(MaxFile)=-1
  real               :: dt_output(MaxFile)=-1.0
  integer            :: t_output_last(MaxFile)=-1

  ! Plot file name string logicals.
  logical :: IsPlotName_e = .false.

  ! Log file name string logicals.
  logical :: IsLogName_e = .false.

contains

  !===========================================================================
  subroutine write_prefix

    if(iUnitOut==STDOUT_)write(*,'(a)',ADVANCE='NO')trim(StringPrefix)

  end subroutine write_prefix

end module IE_ModIo
