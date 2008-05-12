Module ModIO

  use ModSize

  ! All constants and variables related to Input/Output

  character (len=*), parameter :: IO_dir="./IO2/"
  character (len=*), parameter :: restartIN_dir ="./restartIN/"
  character (len=*), parameter :: restartOUT_dir="./restartOUT/"
  character (len=80) :: filename

  logical :: restart                      ! Use a restart file?
  logical :: restart_Bface=.false.        ! Bface from restart file?
  logical :: restart_ghost=.true.         ! ghost cells from restart file?
  logical :: restart_reals=.false.        ! Real numbers only in restart file?

  logical :: save_restart_file=.true., save_satellite_data=.false., &
       save_plots_amr=.false.,save_logfile=.false.,save_binary=.false.

  ! Unit number names
  ! unit_in,unit_in+1,...,unit_in+10 can be used for nested input files
  integer, parameter :: unit_in=10, unit_out=50
  integer, parameter :: &
       unit_ion=30, &
       unit_log=31, &
       unit_solarwind=33, &
       unit_satellite_=33, &
       unit_tmp=99

  ! Maximum number of output files and output variables
  ! note that:
  !     maxfile > maxplotfile + maxsatellitefile + extras
  ! is required
  integer, parameter :: maxplotfile=15
  integer, parameter :: maxsatellitefile=10
  integer, parameter :: maxfile = 30

  ! Ifile index names
  integer, parameter :: restart_=1, &
       logfile_=2, plot_=2, satellite_ = plot_+maxplotfile

  ! Actual number of output files and plot files
  ! note that nfile is not the number of output files but rather the 
  ! index of the maximum file number.  The array FileUsed contains a 
  ! value that tells whether or not each file is used or not.
  integer :: nfile=0, nplotfile=0, nsatellite=0

  ! Saving frequencies and the last saved time step and snapshot number
  real,    dimension(maxfile) :: dt_output=-1.
  integer, dimension(maxfile) :: dn_output=-1, n_output_last=0, t_output_last=0

  integer :: dn_progress1=10, dn_progress2=100

  character (LEN=10) :: plot_type(maxfile), plot_type1
  character (LEN=3)  :: plot_form(maxfile)
  character (LEN=3) :: log_form
  logical :: write_speed_files

  ! x1, x2, y1, y2, z1, z2 limits for plotting
  real, dimension(6,maxfile) :: plot_range 

  ! dx resolution for equidistant plotting
  real, dimension(3,maxfile) :: plot_dx

  ! variables to plot
  character (len=100) :: plot_vars(maxfile), plot_vars1
  character (len=50)  :: plot_pars(maxfile), plot_pars1

  ! variables to put in log file
  character (len=100) :: log_vars, log_R_str
  real :: R_log

  ! variables to control time output format 
  character (len=100) :: log_time, sat_time(maxfile)

  ! variables to write to the satellite files
  character (len=100) :: satellite_vars(maxfile)
  
  ! dimensionalize the output
  logical :: plot_dimensional(maxfile)    

  !\
  ! Variables for the satellite locations
  !/
  logical, dimension(maxsatellitefile,nBLK)                   :: SatelliteInBLK=.false.
  logical, dimension(maxsatellitefile)                        :: TrackSatellite = .false.
  logical, dimension(maxsatellitefile)                        :: TrackSatelliteOLD = .false.
  logical, dimension(maxsatellitefile)                        :: UseSatelliteFile = .true.
  logical, dimension(maxsatellitefile)                        :: Satellite_first_write = .true.
  integer, parameter                                          :: Max_Satellite_Npts = 5000
  integer, dimension(maxsatellitefile)                        :: Satellite_Npts, icurrent_satellite_position=1
  integer, dimension(maxsatellitefile)                        :: iPEsatellite, iBLKsatellite
  real,    dimension(maxsatellitefile, Max_Satellite_Npts, 3) :: XSatellite_traj
  real,    dimension(maxsatellitefile, 3)                     :: XSatellite
  real*8,  dimension(maxsatellitefile, Max_Satellite_Npts)    :: Satellite_Time
  character (len=50)                                          :: Satellite_name(maxsatellitefile)
  integer, dimension(maxsatellitefile)                        :: Satellite_Coordinate_System

end module ModIO
