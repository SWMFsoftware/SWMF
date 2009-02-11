!^CFG COPYRIGHT UM
Module ModIO

  use ModSize
  use ModIoUnit, ONLY: unit_tmp => UNITTMP_, STDOUT_
  use ModVarIndexes, ONLY: nVar

  implicit none
  SAVE

  ! All constants and variables related to Input/Output for GM

  ! Maximum number of output files and output variables
  ! note that:
  !     maxfile > MaxPlotFile + MaxSatelliteFile + extras
  ! is required. MaxSatelliteFile is defined in ModSatelliteFile.f90
  integer, parameter :: MaxPlotFile=25
  integer, parameter :: MaxFile = 350
  integer, parameter :: nPlotvarLosMax=10
  integer, parameter :: nPlotRfrFreqMax=20
  integer, parameter :: nPlotvarMax = max(30,nVar+10) ! Max number of plot vars
  integer, parameter :: MaxLine=20          ! Max number of lines/plot file

  ! Named indexes for output files
  integer, parameter :: &
       restart_=1, logfile_=2, plot_=2, satellite_ = plot_+MaxPlotFile

  ! I/O
  integer             :: iUnitOut = STDOUT_
  character (len=7)   :: StringPrefix=''
  character (len=100) :: NamePlotDir="GM/IO2/"

  ! Generic file name variable
  character (len=80) :: filename

  ! The largest time unit used in the plot file names in time-accurate runs
  character (len=20) :: NameMaxTimeUnit = 'hour'
  ! Simulation time or physical date and time used in file names
  character (len=14) :: StringDateOrTime


  logical :: restart=.false.        ! read restart file
  logical :: restart_Bface =.false. ! Bface restarted ^CFG IF CONSTRAINB

  logical :: DoSaveInitial = .false.

  logical :: save_restart_file=.true., &
       save_plots_amr=.false.,save_logfile=.false.,save_binary=.true.

  ! Unit numbers for the log file and the second SPH plot file
  integer :: unit_log = -1, unit_tmp2 = -1

  ! variables for the line of sight integration plots
  character (LEN=10) :: TypeLosImage
  integer :: n_pix_r(maxfile)
  real :: r_size_image(maxfile), xoffset(maxfile), yoffset(maxfile)
  real :: radius_occult(maxfile), mu_los
  real :: offset_angle(maxfile)
  real, dimension(3,maxfile) :: ObsPos_DI

  ! Variables for radiowave image
  ! ObsPos_DI is borrowed from the LOS plot
  integer, dimension(MaxFile) :: n_Pix_X, n_Pix_Y
  real,    dimension(MaxFile) :: X_Size_Image, Y_Size_Image
  ! String read from PARAM.in, like '1500kHz, 11MHz, 42.7MHz, 1.08GHz':
  character(len=100), dimension(MaxFile) :: StringRadioFrequency_I  
  real, dimension(MaxFile,nPlotRfrFreqMax) :: RadioFrequency_II


  ! Variables for field/stream/current line files
  logical :: IsSingleLine_I(MaxPlotFile)      ! One subfile for the plot file?
  integer :: nLine_I(MaxPlotFile)             ! Number of lines for a plot file
  character :: NameLine_I(MaxPlotFile)                 ! Name of vector field
  real      :: XyzStartLine_DII(3,MaxLine,MaxPlotFile) ! Starting positions
  logical   :: IsParallelLine_II(MaxLine,MaxPlotFile)  ! Parallel/anti-parallel

  ! Actual number of output files and plot files
  ! note that nfile is not the number of output files but rather the 
  ! index of the maximum file number.  The array FileUsed contains a 
  ! value that tells whether or not each file is used or not.
  integer :: nFile=0, nPlotFile=0

  ! Saving frequencies and the last saved time step and snapshot number
  real,    dimension(maxfile) :: dt_output=-1.
  integer, dimension(maxfile) :: dn_output=-1, &
       n_output_last=-1, t_output_last=-1

  ! Frequency of writing progress reports in terms of time steps
  integer :: dn_progress1=10, dn_progress2=100

  character (LEN=10) :: plot_type(maxfile), plot_type1
  character (LEN=3)  :: plot_form(maxfile)
  character (LEN=3) :: log_form
  character (LEN=10) :: TypeIdlFile_I(maxfile)

  ! x1, x2, y1, y2, z1, z2 limits for plotting
  real, dimension(6,MaxFile) :: plot_range 

  ! x, y, z point for arbitrary slice plotting
  real, dimension(3,MaxFile) :: plot_point=-99999.

  ! x, y, z normal vector for arbitrary slice plotting
  real, dimension(3,MaxFile) :: plot_normal

  ! dx resolution for equidistant plotting
  real, dimension(3,MaxFile) :: plot_dx

  ! variables to plot
  character (len=500) :: plot_vars(MaxFile), plot_vars1
  character (len=50)  :: plot_pars(MaxFile), plot_pars1

  ! variables to put in log file
  character (len=500) :: log_vars, log_R_str

  ! variables to control time output format 
  character (len=100) :: log_time

  ! dimensionalize the output
  logical :: plot_dimensional(MaxFile)    

  ! Plot variable names and units defined in the user module
  character(len=10), dimension(nPlotVarMax) :: &
       NameVarUserTec_I, NameUnitUserTec_I, NameUnitUserIdl_I

contains

  !===========================================================================
  subroutine write_prefix

    use ModMain, ONLY: IsStandAlone

    if(IsStandAlone) RETURN
    if(iUnitOut==STDOUT_)write(*,'(a)',ADVANCE='NO')trim(StringPrefix)

  end subroutine write_prefix

  !===========================================================================
  subroutine write_myname

    use ModMain, ONLY: NameThisComp, IsStandAlone

    if(IsStandAlone) RETURN
    if(len_trim(NameThisComp)>0) &
         write(*,'(a)',ADVANCE='NO')NameThisComp//':'

  end subroutine write_myname

end module ModIO
