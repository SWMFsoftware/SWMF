!  Copyright (C) 2002 Regents of the University of Michigan
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!This code is a copyright protected software (c) 2002- University of Michigan
Module ModIO

  use ModSize
  use ModIoUnit, ONLY: UnitTmp_, STDOUT_
  use ModVarIndexes, ONLY: nVar

  implicit none
  SAVE

  ! All constants and variables related to Input/Output for GM

  ! Maximum number of output files and output variables
  ! note that:
  !     MaxFile > MaxPlotFile + MaxSatelliteFile + extras
  ! is required. MaxSatelliteFile is defined in ModSatelliteFile.f90
  integer, parameter :: MaxPlotFile=25
  integer, parameter :: MaxFile = 350
  integer, parameter :: nPlotvarLosMax=20
  integer, parameter :: nPlotRfrFreqMax=20
  integer, parameter :: nPlotvarMax = max(30,nVar+10) ! Max number of plot vars
  integer, parameter :: MaxLine=20          ! Max number of lines/plot file
  integer, parameter :: lNameLogVar = 20    ! Max length of NameLogVar

  ! Named indexes for output files
  integer, parameter :: &
       restart_=1, logfile_=2, magfile_=3, indexfile_=4, maggridfile_=5, &
       plot_=5, satellite_ = plot_+MaxPlotFile
  
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
  logical :: IsRestartCoupler=.false. ! Informs coupler if restarted

  logical :: DoSaveInitial = .false.

  logical :: save_restart_file=.true., &
       save_plots_amr=.false.,save_logfile=.false.,save_binary=.true.

  ! Unit numbers for the log file
  integer :: unit_log = -1

  ! variables for the line of sight integration plots
  character (LEN=10) :: TypeLosImage
  integer :: n_pix_r(MaxFile)
  real :: r_size_image(MaxFile), xoffset(MaxFile), yoffset(MaxFile)
  real :: radius_occult(MaxFile), mu_los
  real :: offset_angle(MaxFile)
  real, dimension(3,MaxFile) :: ObsPos_DI
  character (LEN=20) :: NameLosTable(MaxFile)
  logical:: UseLosSimple = .false. !!! experiment with simple LOS algorithm

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
  real,    dimension(MaxFile) :: dt_output=-1.
  integer, dimension(MaxFile) :: dn_output=-1, &
       n_output_last=-1, t_output_last=-1

  ! Frequency of writing progress reports in terms of time steps
  integer :: dn_progress1=10, dn_progress2=100

  character(LEN=10) :: plot_type(MaxFile), plot_type1
  character(LEN=3)  :: plot_form(MaxFile)
  character(LEN=3)  :: log_form
  character(LEN=10) :: TypeFile_I(MaxFile)

  character(LEN=3)  :: TypeCoordPlot_I(MaxFile) = '???'

  ! x1, x2, y1, y2, z1, z2 limits for plotting
  real, dimension(6,MaxFile) :: plot_range=0. 
  ! plot range for the current file
  real :: PlotRange_I(6)

  ! x, y, z point for arbitrary slice plotting
  real, dimension(3,MaxFile) :: plot_point=-99999.

  ! x, y, z normal vector for arbitrary slice plotting
  real, dimension(3,MaxFile) :: plot_normal

  ! dx resolution for equidistant plotting
  real, dimension(3,MaxFile) :: plot_dx

  ! variables to plot
  character (len=500) :: plot_vars(MaxFile), plot_vars1
  character (len=500) :: plot_pars(MaxFile), plot_pars1

  ! variables to put in log file
  character (len=500) :: log_vars, log_R_str

  ! variables to control time output format 
  character (len=100) :: log_time

  ! dimensionalize the output
  logical :: plot_dimensional(MaxFile)    

  ! Plot variable names and units defined in the user module
  character(len=10), dimension(nPlotVarMax) :: &
       NameVarUserTec_I, NameUnitUserTec_I, NameUnitUserIdl_I

  ! Plot file name string logicals.
  !  One of the three must be true
  !  Time accurate flag can override
  logical :: IsPlotName_n = .true.   !  true if time accurate false
  logical :: IsPlotName_t = .true.   ! false if time accurate false
  logical :: IsPlotName_e = .false.  ! false if time accurate false

  ! Log file name string logicals.
  !  One of the two must be true
  !  Time accurate flag false will automatically set IsLogName_n true
  !     when file is opened.
  logical :: IsLogName_n = .true.
  logical :: IsLogName_e = .false.

  logical :: DoSaveOneTecFileOrig = .false., DoSaveOneTecFile

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
