!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!=====================================
! IONO: Stand-Alone Ionosphere Model |
!=====================================
program IONO
  use ModProcIE
  use ModIonosphere
  use ModMain
  use ModIO
  use ModPhysics
  implicit none
  
  integer :: n_start, n_end, dummy, ierr
  real    :: dt_dummy
  
  ! One line of input
  character (len=100) :: line

  character (LEN=*), parameter :: inputfile="PARAM.in"
  integer, parameter :: unitin = 24

  logical :: done

  ! Initialize Ionosphere
  call ionosphere(0,1)
  
  !\\\
  !
  ! Solar EUV depends on the day of year and the solar EUV
  !   intensity depends on the F10.7 flux.
  !
  ! Here are some reasonable values for these constants.
  !
  conductance_model = 5
  f107_flux = 150.0
  Start_Time_Array(1) = 1982
  Start_Time_Array(2) = 06
  Start_Time_Array(3) = 21
  Start_Time_Array(4) = 00
  Start_Time_Array(5) = 00
  Start_Time_Array(6) = 00
  Start_Time_Array(7) = 00
  n_start = 2000
  n_end = 2000
  StarLightPedConductance = 0.1
  PolarCapPedConductance = 0.50

  done = .false.

  open(20,file='IONO.in',status="old",iostat = ierr)
  if (ierr /=0) done = .true.

  write(*,*) done

  do while (.not.done)
     read(20,'(a)',iostat = ierr) line

     if (ierr /=0) then
        done = .true.
     else

        if(index(line,'#END')>0)then
           done = .true.
        else

           if(index(line,'#SETREALTIME')>0)then

              read(20,*) Start_Time_Array(1)
              read(20,*) Start_Time_Array(2)
              read(20,*) Start_Time_Array(3)
              read(20,*) Start_Time_Array(4)
              read(20,*) Start_Time_Array(5)
              read(20,*) Start_Time_Array(6)

           endif

           if(index(line,'#IONOSPHERE')>0)then

              read(20,*) UseIonosphere
              read(20,*) conductance_model
              read(20,*) UseFullCurrent
              read(20,*) f107_flux
              read(20,*) StarLightPedConductance
              read(20,*) PolarCapPedConductance
              read(20,*) dn_couple_ionosphere

           endif

        endif

     endif

  enddo

  Time_Array = Start_Time_Array

  UseUAM = .false.                              !^CFG  IF TIEGCM
  UseIM  = .false.
  AlignDipoleCorotation = .true.
  SetDipoleTilt = .false.
  iProc = 0
  problem_type = problem_earth
  Max_DoY_Tilt = 23.5 * 2.0 * iono_pi / 360.0
  Magnetic_Pole_Colat = 11.0 * 2.0 * iono_pi / 360.0
  Magnetic_Pole_East_Lon = (-70.9 + 360.0) * iono_pi / 180.0
  rot_period_dim = 24.0
  rBody = 6378.00E03

  UseFakeRegion2 = .true.
             
  call Calculate_Dipole_Tilt

  !
  !///
  
  !\
  ! Peform ionosphere calculations.
  !/
  
  IONO_NORTH_nMagBndPts = 100
  IONO_SOUTH_nMagBndPts = 100
  
  plot_vars(1) = 'minimum'
  plot_form(1) = 'idl'

  do n_step=n_start,n_end
     
     ! Read in restart file
     call ionosphere(n_step,3)
     
     ! Solve for potential
     call ionosphere(n_step,6)
     
     !\
     ! Save ionosphere solution in output data file.
     !/

     call ionosphere_write_output(1)

     call increment_real_world_time(dt_dummy)
     
  enddo
  
  write(*,'(/,1X,''IONO: Execution complete.'',/)')
  
end program IONO
