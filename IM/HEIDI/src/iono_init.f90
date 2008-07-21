
!*************************************************************************
!
! INITIALIZATION Routines
!
!*************************************************************************

!-------------------------------------------------------------------------
! ionosphere_fine_grid
!
!
!
!-------------------------------------------------------------------------

subroutine ionosphere_fine_grid
  !\
  ! This routine sets the fine grid meshes for the
  ! northern and southern hemispheres.
  !/
  use ModIonosphere
  use ModPhysics
  implicit none

  integer :: i,j
  real :: dTheta_l, dPsi_l
        
  IONO_PI = 2.00*asin(1.00)

! These are factors taken out of a run......

  SSPearth = 1.0E+3*50.0
  IONO_Ref_Density = 1.0e06*CON_mp*5.0
  IONO_Ref_SoundSpeed = SSPearth 

  Bdp = -0.58528E+04

  IONO_Bdp = Bdp* &
       sqrt(IONO_MU*IONO_Ref_Density* &
       IONO_Ref_SoundSpeed*IONO_Ref_SoundSpeed)
  
  IONO_Radius = 6372.0 * 1.0e3 + 110.0 * 1.0e3

  dTheta_l = 0.50*IONO_PI/real(IONO_nTheta-1)
  dPsi_l = dTheta_l

  do j = 1, IONO_nPsi
     IONO_NORTH_Theta(1,j) = IONO_Theta_0
     IONO_NORTH_Psi(1,j) = real(j-1)*dPsi_l
     IONO_NORTH_X(1,j) = sin(IONO_NORTH_Theta(1,j))* &
                         cos(IONO_NORTH_Psi(1,j))
     IONO_NORTH_Y(1,j) = sin(IONO_NORTH_Theta(1,j))* &
                         sin(IONO_NORTH_Psi(1,j))
     IONO_NORTH_Z(1,j) = cos(IONO_NORTH_Theta(1,j))

     do i = 1, IONO_nTheta
        IONO_NORTH_Theta(i,j) = IONO_PI*(0.5*real(i-1)/real(IONO_nTheta))
        IONO_NORTH_Psi(i,j) = IONO_NORTH_Psi(1,j)
        IONO_NORTH_X(i,j) = sin(IONO_NORTH_Theta(i,j))* &
                            cos(IONO_NORTH_Psi(i,j))
        IONO_NORTH_Y(i,j) = sin(IONO_NORTH_Theta(i,j))* &
                            sin(IONO_NORTH_Psi(i,j))
        IONO_NORTH_Z(i,j) = cos(IONO_NORTH_Theta(i,j))
     end do  

     IONO_NORTH_Theta(IONO_nTheta,j) = 0.50*IONO_PI
     IONO_NORTH_Psi(IONO_nTheta,j) = IONO_NORTH_Psi(1,j)
     IONO_NORTH_X(IONO_nTheta,j) = sin(IONO_NORTH_Theta(IONO_nTheta,j))* &
                                   cos(IONO_NORTH_Psi(IONO_nTheta,j))
     IONO_NORTH_Y(IONO_nTheta,j) = sin(IONO_NORTH_Theta(IONO_nTheta,j))* &
                                   sin(IONO_NORTH_Psi(IONO_nTheta,j))
     IONO_NORTH_Z(IONO_nTheta,j) = cos(IONO_NORTH_Theta(IONO_nTheta,j))
  end do

  do j = 1, IONO_nPsi
     IONO_SOUTH_Theta(1,j) = IONO_PI-0.50*IONO_PI
     IONO_SOUTH_Psi(1,j) = real(j-1)*dPsi_l
     IONO_SOUTH_X(1,j) = sin(IONO_SOUTH_Theta(1,j))* &
                         cos(IONO_SOUTH_Psi(1,j))
     IONO_SOUTH_Y(1,j) = sin(IONO_SOUTH_Theta(1,j))* &
                         sin(IONO_SOUTH_Psi(1,j))
     IONO_SOUTH_Z(1,j) = cos(IONO_SOUTH_Theta(1,j))

     do i = 1, IONO_nTheta
        IONO_SOUTH_Theta(i,j) = IONO_PI*(0.5+0.5*real(i-1)/real(IONO_nTheta))
        IONO_SOUTH_Psi(i,j) = IONO_SOUTH_Psi(1,j)
        IONO_SOUTH_X(i,j) = sin(IONO_SOUTH_Theta(i,j))* &
                            cos(IONO_SOUTH_Psi(i,j))
        IONO_SOUTH_Y(i,j) = sin(IONO_SOUTH_Theta(i,j))* &
                            sin(IONO_SOUTH_Psi(i,j))
        IONO_SOUTH_Z(i,j) = cos(IONO_SOUTH_Theta(i,j))
     end do  

     IONO_SOUTH_Theta(IONO_nTheta,j) = IONO_PI-IONO_Theta_0
     IONO_SOUTH_Psi(IONO_nTheta,j) = IONO_SOUTH_Psi(1,j)
     IONO_SOUTH_X(IONO_nTheta,j) = sin(IONO_SOUTH_Theta(IONO_nTheta,j))* &
                                   cos(IONO_SOUTH_Psi(IONO_nTheta,j))
     IONO_SOUTH_Y(IONO_nTheta,j) = sin(IONO_SOUTH_Theta(IONO_nTheta,j))* &
                                   sin(IONO_SOUTH_Psi(IONO_nTheta,j))
     IONO_SOUTH_Z(IONO_nTheta,j) = cos(IONO_SOUTH_Theta(IONO_nTheta,j))
  end do

!
! dPsi is going to be defined as 2*dPsi for all of the points:
!

  do j = 2, IONO_nPsi-1
     dPsi_North(j) = IONO_NORTH_Psi(1,j+1)-IONO_NORTH_Psi(1,j-1)
  enddo
  dPsi_North(1)    = IONO_NORTH_Psi(1,2) -                                    &
       (IONO_NORTH_Psi(1,IONO_nPsi-1) - 2.0*IONO_PI)
  dPsi_North(IONO_nPsi) = dPsi_North(1)

!
! dTheta is going to be defined as 2*dTheta for all but the top and bottom
!   points.  As these locations, we switch to 1st order
!

  do i = 2, IONO_nTheta-1
     dTheta_North(i)        = IONO_NORTH_Theta(i+1,1) - IONO_NORTH_Theta(i-1,1)
  enddo
  dTheta_North(1)           = IONO_NORTH_Theta(2,1) - IONO_NORTH_Theta(1,1)
  dTheta_North(IONO_nTheta) =                                                &
       IONO_NORTH_Theta(IONO_nTheta,1)-IONO_NORTH_Theta(IONO_nTheta-1,1)

!
! dPsi is going to be defined as 2*dPsi for all of the points:
!

  do j = 2, IONO_nPsi-1
     dPsi_South(j)      = IONO_SOUTH_Psi(1,j+1)-IONO_SOUTH_Psi(1,j-1)
  enddo
  dPsi_South(1)         = IONO_SOUTH_Psi(1,2)-                               &
       (IONO_SOUTH_Psi(1,IONO_nPsi-1)-2.0*IONO_PI)
  dPsi_South(IONO_nPsi) = dPsi_South(1)

!
! dTheta is going to be defined as 2*dTheta for all but the top and bottom
!   points.  As these locations, we switch to 1st order
!

  do i = 2, IONO_nTheta-1
     dTheta_South(i)        = IONO_SOUTH_Theta(i+1,1) - IONO_SOUTH_Theta(i-1,1)
  enddo
  dTheta_South(1)           = IONO_SOUTH_Theta(2,1) - IONO_SOUTH_Theta(1,1)
  dTheta_South(IONO_nTheta) = IONO_SOUTH_Theta(IONO_nTheta,1) -              &
       IONO_SOUTH_Theta(IONO_nTheta-1,1)

end subroutine ionosphere_fine_grid

!-------------------------------------------------------------------------
! ionosphere_init
!
!
!
!-------------------------------------------------------------------------

subroutine ionosphere_init(year, day, ut)
  !\
  ! This routine initializes the fine grid 
  ! ionospheric solutions for the
  ! northern and southern hemispheres.
  !/
  use ModIonosphere
  use ModPhysics
  use ModMain
  implicit none

  integer, intent(in) :: year, day
  real, intent(in) :: ut

  IONO_NORTH_PHI = 0.00
  IONO_NORTH_JR = 0.00 
  IONO_NORTH_Jx = 0.00
  IONO_NORTH_Jy = 0.00
  IONO_NORTH_Jz = 0.00
  IONO_NORTH_Ex = 0.00
  IONO_NORTH_Ey = 0.00
  IONO_NORTH_Ez = 0.00
  IONO_NORTH_ETh = 0.00
  IONO_NORTH_EPs = 0.00
  IONO_NORTH_Ux = 0.00
  IONO_NORTH_Uy = 0.00
  IONO_NORTH_Uz = 0.00

  IONO_SOUTH_PHI = 0.00
  IONO_SOUTH_JR = 0.00
  IONO_SOUTH_Jx = 0.00
  IONO_SOUTH_Jy = 0.00
  IONO_SOUTH_Jz = 0.00
  IONO_SOUTH_Ex = 0.00
  IONO_SOUTH_Ey = 0.00
  IONO_SOUTH_Ez = 0.00
  IONO_SOUTH_ETh = 0.00
  IONO_SOUTH_EPs = 0.00
  IONO_SOUTH_Ux = 0.00
  IONO_SOUTH_Uy = 0.00
  IONO_SOUTH_Uz = 0.00  

  IONO_NORTH_PHI_BC = 0.00
  IONO_NORTH_ETh_BC = 0.00
  IONO_NORTH_EPs_BC = 0.00
  IONO_NORTH_UR_BC = 0.00 
  IONO_NORTH_UTh_BC = 0.00 
  IONO_NORTH_UPs_BC = 0.00

  IONO_SOUTH_PHI_BC = 0.00
  IONO_SOUTH_ETh_BC = 0.00
  IONO_SOUTH_EPs_BC = 0.00
  IONO_SOUTH_UR_BC = 0.00
  IONO_SOUTH_UTh_BC = 0.00 
  IONO_SOUTH_UPs_BC = 0.00

  IONO_NORTH_TGCM_JR = 0.00
  IONO_SOUTH_TGCM_JR = 0.00

  IONO_Old_Seconds = 0.00

  new_conductance = .true.

  start_time_array(1) = year
  start_time_array(2) = 1
  start_time_array(3) = day
  start_time_array(4) = 0.0
  start_time_array(5) = 0.0
  start_time_array(6) = ut
  start_time_array(7) = 0.0

end subroutine ionosphere_init

!*************************************************************************
!
! INPUT/OUTPUT Routines
!
!*************************************************************************

!-------------------------------------------------------------------------
! write_timegcm_file
!
!
!
!-------------------------------------------------------------------------

subroutine write_timegcm_file(iter, phi_north, phi_south,   &
     eflux_north, eflux_south, avee_north, avee_south)

  use ModIonosphere
  use ModMain
  implicit none

  integer, intent(in) :: iter
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::             &
       phi_north, phi_south, eflux_north, eflux_south,      &
       avee_north, avee_south
  real psi_offset

  integer :: Iunit, i, j
  character (len=13), Parameter :: iono_dir="./ionosphere/"
  character (len=4), Parameter :: IO_ext=".dat"

  write(*,*) '=> Writing output datafiles for TIMEGCM.'

  Iunit = 23
  open(unit=Iunit,file=iono_dir//"MHD_to_TIMEGCM"//IO_ext, &
       status="unknown")

!  Write Year, month, day, hour, minute, second

  write(Iunit,fmt="(6I5)") (Start_Time_Array(i),i=1,6)

! Write number of thetas and psis

  write(Iunit,fmt="(2I10)") IONO_nTheta*2, IONO_nPsi

! Silliness to occur here:
!   This ionospheric code starts with noon at 0.0 longitude.  This is 
!   different than almost all other ionospheres, so we need to shift
!   the output by 180.0 degrees.
!   Also, we want 0.0 and 360.0 to both be printed - but the new ones,
!   so we have to ignore the old 0.0 and repeat the old 180.0.

  do i = IONO_nTheta, 1, -1
     do j = IONO_nPsi/2+1, IONO_nPsi
        psi_offset = (360.00/(2.00*IONO_PI))*IONO_SOUTH_Psi(i,j) - 180.0
        if (psi_offset < 0) psi_offset = psi_offset + 360.0
        write(Iunit,fmt="(1P,5(E13.5))")  &
             (360.00/(2.00*IONO_PI))*IONO_SOUTH_Theta(i,j), &
             psi_offset, &
             1.0e-03*phi_south(i,j),eflux_south(i,j), &
             avee_south(i,j)
     end do
     do j = 2, IONO_nPsi/2+1
        psi_offset = (360.00/(2.00*IONO_PI))*IONO_SOUTH_Psi(i,j) - 180.0
        if (psi_offset < 0) psi_offset = psi_offset + 360.0
        write(Iunit,fmt="(1P,5(E13.5))")  &
             (360.00/(2.00*IONO_PI))*IONO_SOUTH_Theta(i,j), &
             psi_offset, &
             1.0e-03*phi_south(i,j),eflux_south(i,j), &
             avee_south(i,j)
     end do
  end do
  do i = IONO_nTheta, 1, -1
     do j = IONO_nPsi/2+1, IONO_nPsi
        psi_offset = (360.00/(2.00*IONO_PI))*IONO_SOUTH_Psi(i,j) - 180.0
        if (psi_offset < 0) psi_offset = psi_offset + 360.0
        write(Iunit,fmt="(1P,5(E13.5))")  &
             (360.00/(2.00*IONO_PI))*IONO_NORTH_Theta(i,j), &
             psi_offset, &
             1.0e-03*phi_north(i,j),eflux_north(i,j), &
             avee_north(i,j)
     end do
     do j = 2, IONO_nPsi/2+1
        psi_offset = (360.00/(2.00*IONO_PI))*IONO_SOUTH_Psi(i,j) - 180.0
        if (psi_offset < 0) psi_offset = psi_offset + 360.0
        write(Iunit,fmt="(1P,5(E13.5))")  &
             (360.00/(2.00*IONO_PI))*IONO_NORTH_Theta(i,j), &
             psi_offset, &
             1.0e-03*phi_north(i,j),eflux_north(i,j), &
             avee_north(i,j)
     end do
  end do

  close(UNIT=Iunit)
 
end subroutine write_timegcm_file

!-------------------------------------------------------------------------
! read_timegcm_file
!
!
!
!-------------------------------------------------------------------------

subroutine read_timegcm_file

  use ModIonosphere
  implicit none

  integer :: Iunit, i, j, iter, nTheta, nPsi
  character (len=13), Parameter :: iono_dir="./"
  character (len=4), Parameter :: IO_ext=".dat"
  character (len=100) :: TIMEGCM_file
  real :: dum_lat, dum_lon, min_lat, max_lat, max_fac
  integer, dimension(1:6) :: time_array
  real :: tgcm_time

  min_lat = 45.0
  max_lat = 85.0
  max_fac = 0.0

  write(*,*) '=> Reading output file from the TIMEGCM.'

  TIMEGCM_file = "./ionosphere/TIMEGCM_to_MHD.dat"

  Iunit = 24
  open(unit=Iunit,file=TIMEGCM_file, status="old")

  read(Iunit,fmt="(6I5)") (time_array(i),i=1,6)

  read(Iunit,fmt="(2I10)") nTheta, nPsi

  do i = IONO_nTheta,1,-1
     do j = 1, IONO_nPsi
        read(Iunit,fmt="(3(E13.5))")  &
             dum_lat, dum_lon, IONO_SOUTH_TGCM_JR(i,j)

        if (dum_lat.lt.(90+min_lat)) then
           if (abs(iono_south_tgcm_jr(i,j)) > max_fac) &
                max_fac = abs(iono_south_tgcm_jr(i,j))
           IONO_SOUTH_TGCM_JR(i,j) = 0.0
        endif

        if (dum_lat.gt.(90+max_lat)) then
           IONO_SOUTH_TGCM_JR(i,j) = 0.0
        endif

     end do
  end do

  if (min_lat > 0.0) &
       write(6,*) "=> Maximum abs(FAC) below ",min_lat," (SH) is ",max_fac

  max_fac = 0.0
  do i = IONO_nTheta,1,-1
     do j = 1, IONO_nPsi
        read(Iunit,fmt="(3(E13.5))")  &
             dum_lat, dum_lon, IONO_NORTH_TGCM_JR(i,j)

        if (dum_lat.gt.(90-min_lat)) then
           if (abs(iono_north_tgcm_jr(i,j)) > max_fac) &
                max_fac = abs(iono_south_tgcm_jr(i,j))
           IONO_NORTH_TGCM_JR(i,j) = 0.0
        endif

        if (dum_lat.lt.(90-max_lat)) then
           IONO_NORTH_TGCM_JR(i,j) = 0.0
        endif
     end do
  end do

  if (min_lat > 0.0) &
       write(6,*) "=> Maximum abs(FAC) below ",min_lat," (NH) is ",max_fac

  close(UNIT=Iunit)
 
end subroutine read_timegcm_file

!-------------------------------------------------------------------------
! read_rcm_file
!
!
!
!-------------------------------------------------------------------------

subroutine read_rcm_file

  use ModIonosphere
  implicit none

  integer :: Iunit, i, j, iter, nTheta, nPsi
  character (len=13), Parameter :: iono_dir="../d.timegcm/"
  character (len=4), Parameter :: IO_ext=".dat"
  character (len=100) :: TIMEGCM_file
  real :: dum_lat, dum_lon, min_lat, max_lat
  integer, dimension(1:6) :: time_array
  real :: tgcm_time

  min_lat = 5.0
  max_lat = 180.0 - min_lat

  write(*,*) '=> Reading output file from the TIMEGCM.'

  Iunit = 23
  open(unit=Iunit,file=iono_dir//"TIMEGCM_to_MHD_filemanager"//IO_ext, &
       status="old")
  read(Iunit,fmt="(A100)") TIMEGCM_file
  close(UNIT=Iunit)

  open(unit=Iunit,file=iono_dir//TIMEGCM_file, status="old")

  read(Iunit,fmt="(6I5)") (time_array(i),i=1,6)

  read(Iunit,fmt="(2I10)") nTheta, nPsi

  do i = 1, IONO_nTheta
     do j = 1, IONO_nPsi
        read(Iunit,fmt="(3(E13.5))")  &
             dum_lat, dum_lon, IONO_NORTH_TGCM_JR(i,j)
        if (dum_lat.lt.min_lat) IONO_NORTH_TGCM_JR(i,j) = 0.0
     end do
  end do

  do i = 1, IONO_nTheta
     do j = 1, IONO_nPsi
        read(Iunit,fmt="(3(E13.5))")  &
             dum_lat, dum_lon, IONO_SOUTH_TGCM_JR(i,j)
        if (dum_lat.gt.max_lat) IONO_SOUTH_TGCM_JR(i,j) = 0.0
     end do
  end do

  close(UNIT=Iunit)
 
end subroutine read_rcm_file

!-------------------------------------------------------------------------
! ionosphere_write_output
!
!
!
!-------------------------------------------------------------------------

subroutine ionosphere_write_output(ifile,time,IO_prefix,IO_suffix)
  !\
  ! This routine writes out the fine grid 
  ! ionospheric solutions for the
  ! northern and southern hemispheres to
  ! a data file.
  !/
  use ModIonosphere
  use ModMain
  use ModIO
  implicit none

  integer, intent(in) :: ifile
  real, intent(in) :: time
  character (len=5), intent(in) :: IO_prefix
  character (len=3), intent(in) :: IO_suffix

  integer, parameter :: tec_type = 1
  integer, parameter :: idl_type = 2
  integer, parameter :: all_vars = 1
  integer, parameter :: min_vars = 2
  integer, parameter :: uam_vars = 3

  integer :: output_type, variables, year, iError

  integer :: Iunit, i, j
!  character (len=13), Parameter :: iono_dir="./ionosphere/"
  character (len=12), Parameter :: iono_dir="_ionosphere/"
  character (len=2), Parameter :: iono_dir_path="./"
  character (len=4) :: IO_ext=".dat"
  character (len=7) :: NUMiter
  character (len=2) :: syear,smonth, sday, shour, sminute, ssecond
  character (len=4) :: sms

  real :: day,hour,minute

  variables = min_vars
  output_type = idl_type

  if (plot_vars(ifile) == 'minimum') then
     variables = min_vars
  elseif (plot_vars(ifile) == 'maximum') then
     variables = all_vars
  elseif (plot_vars(ifile) == 'uam') then
     variables = uam_vars
  endif

  if (plot_form(ifile) == 'idl') then
     output_type = idl_type
  elseif (plot_form(ifile) == 'tec') then
     output_type = tec_type
  endif

  Iunit = 23

  if (output_type == idl_type) then
     IO_ext = ".idl"
  elseif (output_type == tec_type) then
     IO_ext = ".tec"
  endif

!  day  = int(time)
!  hour = int((time - day)*24.0)
!  minute = int(((time - day)*24.0 - hour)*60.0)
!!! time is in seconds from the beginning of the run, not decimal days
  day = int(time/86400.)
  hour = int((time - day)/3600.)
  minute = int((time - day - hour*3600.)/60.)

  if (day < 10) then
     write (sday,'(a1,i1)') '0',int(day)
  else
     write (sday,'(i2)') int(day)
  endif
  
  if (hour < 10) then
     write (shour,'(a1,i1)') '0',int(hour)
  else
     write (shour,'(i2)') int(hour)
  endif
  
  if (minute < 10) then
     write (sminute,'(a1,i1)') '0',int(minute)
  else
     write (sminute,'(i2)') int(minute)
  endif
  
!     open(unit=Iunit,file=iono_dir//"it"//&
!          syear//smonth//sday//"_"&
!          //shour//sminute//ssecond//"_"&
!          //sms//IO_ext, &
!          status="unknown")

!!! For use with day, hour, minute output format:
!  open(unit=Iunit,file=iono_dir//"it"//&
!       sday//shour//sminute//IO_ext, &
!       status="unknown", iostat=iError)
!  if (iError /= 0) then
!     write(*,*) "Error in opening file : ",iono_dir//"it"//&
!       sday//shour//sminute//IO_ext
!     stop
!  endif

!!! For use with mram-given output suffix:

  write(*,*) "file : ",iono_dir_path//IO_prefix//iono_dir//"it_"//&
       IO_suffix//IO_ext

  open(unit=Iunit,file=iono_dir_path//IO_prefix//iono_dir//"it_"//&
       IO_suffix//IO_ext, &
       status="unknown")

  if (output_type == idl_type) then

     write(Iunit, '(a)')  'TITLE='
     write(Iunit, '(a)')  '  "BATSRUS: Ionospheric Potential Solution"'

     write(Iunit, *) ' '

     write(Iunit, '(a)')  'NUMERICAL VALUES'
     if (variables == min_vars) then
        write(Iunit, '(I5,a)')           6, ' nvars'
     elseif (variables == all_vars) then
        write(Iunit, '(I5,a)')           18, ' nvars'
     elseif (variables == uam_vars) then
        write(Iunit, '(I5,a)')           7, ' nvars'
     endif
     write(Iunit, '(I5,a)') IONO_nTheta, ' nTheta'
     write(Iunit, '(I5,a)')   IONO_nPsi, ' nPhi'

     write(Iunit, *) ' '

     write(Iunit, '(a)') 'VARIABLE LIST'

     if (variables == min_vars) then

        write(Iunit, '(I5,a)')  1, ' Theta [deg]'
        write(Iunit, '(I5,a)')  2, ' Psi [deg]'
        write(Iunit, '(I5,a)')  3, ' SigmaH [mhos]'
        write(Iunit, '(I5,a)')  4, ' SigmaP [mhos]'
        write(Iunit, '(I5,a)')  5, ' Jr [mA/m^2]'
        write(Iunit, '(I5,a)')  6, ' Phi [kV]'

     elseif (variables == all_vars) then

        write(Iunit, '(I5,a)')  1, ' X [Re]'
        write(Iunit, '(I5,a)')  2, ' Y [Re]'
        write(Iunit, '(I5,a)')  3, ' Z [Re]'
        write(Iunit, '(I5,a)')  4, ' Theta [deg]'
        write(Iunit, '(I5,a)')  5, ' Psi [deg]'
        write(Iunit, '(I5,a)')  6, ' SigmaH [mhos]'
        write(Iunit, '(I5,a)')  7, ' SigmaP [mhos]'
        write(Iunit, '(I5,a)')  8, ' Jr [mA/m^2]'
        write(Iunit, '(I5,a)')  9, ' Phi [kV]'
        write(Iunit, '(I5,a)') 10, ' Ex [mV/m]'
        write(Iunit, '(I5,a)') 11, ' Ey [mV/m]'
        write(Iunit, '(I5,a)') 12, ' Ez [mV/m]'
        write(Iunit, '(I5,a)') 13, ' Jx [microA/m2]'
        write(Iunit, '(I5,a)') 14, ' Jy [microA/m2]'
        write(Iunit, '(I5,a)') 15, ' Jz [microA/m2]'
        write(Iunit, '(I5,a)') 16, ' Ux [km/s]'
        write(Iunit, '(I5,a)') 17, ' Uy [km/s]'
        write(Iunit, '(I5,a)') 18, ' Uz [km/s]'

     elseif (variables == uam_vars) then

        write(Iunit, '(I5,a)')  1, ' Theta [deg]'
        write(Iunit, '(I5,a)')  2, ' Psi [deg]'
        write(Iunit, '(I5,a)')  3, ' SigmaH [mhos]'
        write(Iunit, '(I5,a)')  4, ' SigmaP [mhos]'
        write(Iunit, '(I5,a)')  5, ' Jr [mA/m^2]'
        write(Iunit, '(I5,a)')  6, ' Jr(NW) [mA/m^2]'
        write(Iunit, '(I5,a)')  7, ' Phi [kV]'

     endif

     write(Iunit, *) ' '

     write(Iunit, '(a)') 'TIME'
     write(Iunit, '(I5,a)')  0, ' Year'
     write(Iunit, '(I5,a)')  0, ' Month'
     write(Iunit, '(I5,a)')  int(day), ' Day'
     write(Iunit, '(I5,a)')  int(hour), ' Hour'
     write(Iunit, '(I5,a)')  int(minute), ' Minute'
     write(Iunit, '(I5,a)')  0, ' Second'
     write(Iunit, '(I5,a)')  0, ' Millisecond'

     write(Iunit, *) ' '
     write(Iunit, '(a)') 'BEGIN NORTHERN HEMISPHERE'

  elseif (output_type == tec_type) then

     write(Iunit, *)  'TITLE="BATSRUS: Ionospheric Potential Solution"'

     if (variables == min_vars) then

        write(Iunit, *)  'VARIABLES= "Theta [deg]","Psi [deg]"'
        write(Iunit, *)  ' "SigmaH [S]","SigmaP [S]"'
        write(Iunit, *)  ' "JR [`mA/m^2]","PHI [kV]"'

     elseif (variables == all_vars) then

        write(Iunit, *)  'VARIABLES= "X [R]","Y [R]","Z [R]"'
        write(Iunit, *)  ' "Theta [deg]","Psi [deg]"'
        write(Iunit, *)  ' "SigmaH [S]","SigmaP [S]"'
        write(Iunit, *)  ' "JR [`mA/m^2]","PHI [kV]"'
        write(Iunit, *)  ' "Ex [mV/m]","Ey [mV/m]","Ez [mV/m]"'
        write(Iunit, *)  ' "Jx [`mA/m^2]","Jy [`mA/m^2]","Jz [`mA/m^2]"'
        write(Iunit, *)  ' "Ux [km/s]","Uy [km/s]","Uz [km/s]"'

     elseif (variables == uam_vars) then


        write(Iunit, *)  'VARIABLES= "Theta [deg]","Psi [deg]"'
        write(Iunit, *)  ' "SigmaH [S]","SigmaP [S]"'
        write(Iunit, *)  ' "JR [`mA/m^2]","JR (NW) [`mA/m^2]","PHI [kV]"'

     endif

     write(Iunit, *)  'ZONE T="Northern Hemisphere" ','I= ',IONO_nTheta, &
          ' J= ',IONO_nPsi,' F=POINT'

  endif

  if (variables == min_vars) then
     do j = 1, IONO_nPsi
        do i = 1, IONO_nTheta
           write(Iunit,fmt="(6(E13.5))")  &
                (360.00/(2.00*IONO_PI))*IONO_NORTH_Theta(i,j), &
                (360.00/(2.00*IONO_PI))*IONO_NORTH_Psi(i,j), &
                IONO_NORTH_SigmaH(i,j),IONO_NORTH_SigmaP(i,j), &
                1.0e06*IONO_NORTH_JR(i,j),   &
                1.0e-03*IONO_NORTH_PHI(i,j)
        end do
     end do
  elseif (variables == all_vars) then
     do j = 1, IONO_nPsi
        do i = 1, IONO_nTheta
           write(Iunit,fmt="(18(E13.5))")  &
                IONO_NORTH_X(i,j),IONO_NORTH_Y(i,j),IONO_NORTH_Z(i,j), &
                (360.00/(2.00*IONO_PI))*IONO_NORTH_Theta(i,j), &
                (360.00/(2.00*IONO_PI))*IONO_NORTH_Psi(i,j), &
                IONO_NORTH_SigmaH(i,j),IONO_NORTH_SigmaP(i,j), &
                1.0e06*IONO_NORTH_JR(i,j),1.0e-03*IONO_NORTH_PHI(i,j), &
                1.0e03*IONO_NORTH_Ex(i,j),1.0e03*IONO_NORTH_Ey(i,j), &
                1.0e03*IONO_NORTH_Ez(i,j), &
                1.0e06*IONO_NORTH_Jx(i,j),1.0e06*IONO_NORTH_Jy(i,j), &
                1.0e06*IONO_NORTH_Jz(i,j), &
                1.0e-03*IONO_NORTH_Ux(i,j),1.0e-03*IONO_NORTH_Uy(i,j), &
                1.0e-03*IONO_NORTH_Uz(i,j)
        end do
     end do
  elseif (variables == uam_vars) then
     do j = 1, IONO_nPsi
        do i = 1, IONO_nTheta
           write(Iunit,fmt="(6(E13.5))")  &
                (360.00/(2.00*IONO_PI))*IONO_NORTH_Theta(i,j), &
                (360.00/(2.00*IONO_PI))*IONO_NORTH_Psi(i,j), &
                IONO_NORTH_SigmaH(i,j),IONO_NORTH_SigmaP(i,j), &
                1.0e06*IONO_NORTH_JR(i,j), 1.0e06*IONO_NORTH_TGCM_JR(i,j), &
                1.0e-03*IONO_NORTH_PHI(i,j)
        end do
     end do
  endif

  if (output_type == idl_type) then
     write(Iunit, *) ' '
     write(Iunit, '(a)') 'BEGIN SOUTHERN HEMISPHERE'
  elseif (output_type == tec_type) then

     write(Iunit, *)  'ZONE T="Southern Hemisphere" ','I= ',IONO_nTheta, &
          ' J= ',IONO_nPsi,' F=POINT'
  endif

  if (variables == min_vars) then
     do j = 1, IONO_nPsi
        do i = 1, IONO_nTheta
           write(Iunit,fmt="(6(E13.5))")  &
                (360.00/(2.00*IONO_PI))*IONO_SOUTH_Theta(i,j), &
                (360.00/(2.00*IONO_PI))*IONO_SOUTH_Psi(i,j), &
                IONO_SOUTH_SigmaH(i,j),IONO_SOUTH_SigmaP(i,j), &
                1.0e06*IONO_SOUTH_JR(i,j),   &
                1.0e-03*IONO_SOUTH_PHI(i,j)
        end do
     end do
  elseif (variables == all_vars) then
     do j = 1, IONO_nPsi
        do i = 1, IONO_nTheta
           write(Iunit,fmt="(18(E13.5))")  &
                IONO_SOUTH_X(i,j),IONO_SOUTH_Y(i,j),IONO_SOUTH_Z(i,j), &
                (360.00/(2.00*IONO_PI))*IONO_SOUTH_Theta(i,j), &
                (360.00/(2.00*IONO_PI))*IONO_SOUTH_Psi(i,j), &
                IONO_SOUTH_SigmaH(i,j),IONO_SOUTH_SigmaP(i,j), &
                1.0e06*IONO_SOUTH_JR(i,j),1.0e-03*IONO_SOUTH_PHI(i,j), &
                1.0e03*IONO_SOUTH_Ex(i,j),1.0e03*IONO_SOUTH_Ey(i,j), &
                1.0e03*IONO_SOUTH_Ez(i,j), &
                1.0e06*IONO_SOUTH_Jx(i,j),1.0e06*IONO_SOUTH_Jy(i,j), &
                1.0e06*IONO_SOUTH_Jz(i,j), &
                1.0e-03*IONO_SOUTH_Ux(i,j),1.0e-03*IONO_SOUTH_Uy(i,j), &
                1.0e-03*IONO_SOUTH_Uz(i,j)
        end do
     end do
  elseif (variables == uam_vars) then
     do j = 1, IONO_nPsi
        do i = 1, IONO_nTheta
           write(Iunit,fmt="(6(E13.5))")  &
                (360.00/(2.00*IONO_PI))*IONO_SOUTH_Theta(i,j), &
                (360.00/(2.00*IONO_PI))*IONO_SOUTH_Psi(i,j), &
                IONO_SOUTH_SigmaH(i,j),IONO_SOUTH_SigmaP(i,j), &
                1.0e06*IONO_SOUTH_JR(i,j), 1.0e06*IONO_SOUTH_TGCM_JR(i,j),  &
                1.0e-03*IONO_SOUTH_PHI(i,j)
        end do
     end do
  endif

  close(UNIT=Iunit)

end subroutine ionosphere_write_output

!-------------------------------------------------------------------------
! ionosphere_read_restart_file
!
!
!
!-------------------------------------------------------------------------

subroutine ionosphere_read_restart_file(iter)
  !\
  ! This routine creates an ionospheric 
  ! restart solution file.
  !/
  use ModIonosphere
  implicit none

  integer, intent(in) :: iter

  integer :: Iunit, i, j, i2, j2, nhemi, ierror
  real :: xx,yy
  character (len=13), Parameter :: iono_dir="./ionosphere/"
  character (len=4), Parameter :: IO_ext=".rst"
  character (len=7) :: NUMiter

  write(*,*) '=> Reading restart file for ionosphere.'

  if (iter < 10) then
     write (NUMiter,'(a6,i1)') 'n00000',iter
  elseif (iter < 100) then
     write (NUMiter,'(a5,i2)') 'n0000',iter
  elseif (iter < 1000) then
     write (NUMiter,'(a4,i3)') 'n000',iter
  elseif (iter < 10000) then
     write (NUMiter,'(a3,i4)') 'n00',iter
  elseif (iter < 100000) then
     write (NUMiter,'(a2,i5)') 'n0',iter
  elseif (iter < 1000000) then
     write (NUMiter,'(a1,i6)') 'n',iter
  else
     write(*,*) "ionosphere_restart: iter = ",iter, &
                " Error, too many iterations for filename"
     stop
  end if

  Iunit = 23
  open(unit=Iunit,file=iono_dir//"ionosphere"//"_"//NUMiter//IO_ext, &
       status="old",iostat=ierror)

  if (ierror.eq.0) then

     read(Iunit,*) IONO_Radius, IONO_Height, &
          IONO_Ref_Density, IONO_Ref_SoundSpeed, &
          IONO_Bdp, IONO_Radius_Mag_Boundary, &
          IONO_NORTH_Theta_Max, IONO_SOUTH_Theta_Min

     do j = 1, IONO_nPsi
        do i = 1, IONO_nTheta
           read(Iunit,*) i2,j2,xx,yy,IONO_NORTH_JR(i,j),                  &
                SAVE_NORTH_SigmaH(i,j),SAVE_NORTH_SigmaP(i,j)
        end do
     end do

     do j = 1, IONO_nPsi
        do i = 1, IONO_nTheta
           read(Iunit,*) i2,j2,xx,yy,IONO_SOUTH_JR(i,j),                  &
                SAVE_SOUTH_SigmaH(i,j),SAVE_SOUTH_SigmaP(i,j) 
        end do
     end do

     close(UNIT=Iunit)

     skip_ionosphere_restart = .false.

  else

     write(6,*) "==================================================="
     write(6,*) "= Error in finding ionosphere restart file.       ="
     write(6,*) "==================================================="

     skip_ionosphere_restart = .true.

  endif

end subroutine ionosphere_read_restart_file

!-------------------------------------------------------------------------
! ionosphere_write_restart_file
!
!
!
!-------------------------------------------------------------------------

subroutine ionosphere_write_restart_file(iter)
  !\
  ! This routine reads in an ionospheric 
  ! restart solution file.
  !/
  use ModIonosphere
  implicit none

  integer, intent(in) :: iter

  integer :: Iunit, i, j
  character (len=13), Parameter :: iono_dir="./ionosphere/"
  character (len=4), Parameter :: IO_ext=".rst"
  character (len=7) :: NUMiter

  write(*,*) '=> Writing restart file for ionosphere.'

  if (iter < 10) then
     write (NUMiter,'(a6,i1)') 'n00000',iter
  elseif (iter < 100) then
     write (NUMiter,'(a5,i2)') 'n0000',iter
  elseif (iter < 1000) then
     write (NUMiter,'(a4,i3)') 'n000',iter
  elseif (iter < 10000) then
     write (NUMiter,'(a3,i4)') 'n00',iter
  elseif (iter < 100000) then
     write (NUMiter,'(a2,i5)') 'n0',iter
  elseif (iter < 1000000) then
     write (NUMiter,'(a1,i6)') 'n',iter
  else
     write(*,*) "ionosphere_restart: iter = ",iter, &
                " Error, too many iterations for filename"
     stop
  end if

  Iunit = 23
  open(unit=Iunit,file=iono_dir//"ionosphere_"//NUMiter//IO_ext,              &
       status="unknown")

  write(Iunit,*) IONO_Radius, IONO_Height,                                    &
                 IONO_Ref_Density, IONO_Ref_SoundSpeed,                       &
                 IONO_Bdp, IONO_Radius_Mag_Boundary,                          &
                 IONO_NORTH_Theta_Max, IONO_SOUTH_Theta_Min

  do j = 1, IONO_nPsi
     do i = 1, IONO_nTheta
        write(Iunit,*) i,j,                                                   &
                       IONO_NORTH_Theta(i,j),IONO_NORTH_Psi(i,j),             &
                       IONO_NORTH_JR(i,j),SAVE_NORTH_SigmaH(i,j),             &
                       SAVE_NORTH_SigmaP(i,j)
     end do
  end do

  do j = 1, IONO_nPsi
     do i = 1, IONO_nTheta
        write(Iunit,*) i,j,                                                   &
                       IONO_SOUTH_Theta(i,j),IONO_SOUTH_Psi(i,j),             &
                       IONO_SOUTH_JR(i,j),SAVE_SOUTH_SigmaH(i,j),             &
                       SAVE_SOUTH_SigmaP(i,j)
     end do
  end do

  close(UNIT=Iunit)

end subroutine ionosphere_write_restart_file

 
