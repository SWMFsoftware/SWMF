!-------------------------------------------------------------------------
! ionosphere_write_output
!
!
!
!-------------------------------------------------------------------------

subroutine ionosphere_write_output(time)
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

  real, intent(in) :: time

  integer, parameter :: tec_type = 1
  integer, parameter :: idl_type = 2
  integer, parameter :: all_vars = 1
  integer, parameter :: min_vars = 2
  integer, parameter :: uam_vars = 3

  integer :: output_type, variables, year

  integer :: Iunit, i, j
  character (len=13), Parameter :: iono_dir="./ionosphere/"
  character (len=4) :: IO_ext=".dat"
  character (len=7) :: NUMiter
  character (len=2) :: syear,smonth, sday, shour, sminute, ssecond
  character (len=4) :: sms

  integer :: day, hour, minute

  variables = min_vars
  output_type = idl_type

  Iunit = 23

  if (output_type == idl_type) then
     IO_ext = ".idl"
  elseif (output_type == tec_type) then
     IO_ext = ".tec"
  endif

!  year = mod(time_array(1),100)
!  if (year < 10) then
!     write (syear,'(a1,i1)') '0',year
!  else
!     write (syear,'(i2)') year
!  endif

!  if (time_array(2) < 10) then
!     write (smonth,'(a1,i1)') '0',time_array(2)
!  else
!     write (smonth,'(i2)') time_array(2)
!  endif
  
  day  = int(time)
  hour = int((time - day)*24.0)
  minute = int(((time - day)*24.0 - hour)*60.0)

  if (day < 10) then
     write (sday,'(a1,i1)') '0',day
  else
     write (sday,'(i2)') day
  endif
  
  if (hour < 10) then
     write (shour,'(a1,i1)') '0',hour
  else
     write (shour,'(i2)') hour
  endif
  
  if (minute < 10) then
     write (sminute,'(a1,i1)') '0',minute
  else
     write (sminute,'(i2)') minute
  endif
  
!     open(unit=Iunit,file=iono_dir//"it"//&
!          syear//smonth//sday//"_"&
!          //shour//sminute//ssecond//"_"&
!          //sms//IO_ext, &
!          status="unknown")

  open(unit=Iunit,file=iono_dir//"it"//&
       sday//shour//sminute//IO_ext, &
       status="unknown",iostat=iError)

  if (iError /= 0) then
     write(*,*) "Error in opening file : ",iono_dir//"it"//&
       sday//shour//sminute//IO_ext
     stop
  endif

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
     write(Iunit, '(I5,a)')  day, ' Day'
     write(Iunit, '(I5,a)')  hour, ' Hour'
     write(Iunit, '(I5,a)')  minute, ' Minute'
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
