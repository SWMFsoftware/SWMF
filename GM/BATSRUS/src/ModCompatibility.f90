!^CFG COPYRIGHT UM
module ModCompatibility

  use ModNumConst

  implicit none
  save

  integer, dimension(7) :: Time_Array

  logical :: SetDipoleTilt=.false.
  logical :: AlignDipoleCorotation=.true.
  logical :: SetCorotationTilt=.false.
  real    :: CorotationTilt=0.0, CorotationTiltDeg=0.0
  real    :: CorotationLon=0.0, CorotationLonDeg=0.0

  real    :: Magnetic_Pole_Colat    = 11.0 * 2.0 * cPi / 360.0
  real    :: Magnetic_Pole_East_Lon = (-70.9 + 360.0) * cPi / 180.0
  real    :: Max_DoY_Tilt           = 23.5 * 2.0 * cPi / 360.0

  integer :: CoordinateSystem

  integer, parameter :: GSM = 1
  integer, parameter :: GSE = 2
  integer, parameter :: GEO = 3

contains

  subroutine read_compatible_command(NameCommand)

    use ModProcMH, ONLY : iProc
    use ModReadParam
    use ModPhysics, ONLY: Bdp_dim, rot_period_dim
    use ModMain, ONLY: ThetaTiltDeg, Dt_UpdateB0, UseRotatingBc, UseStrict

    character (len=*), intent(in) :: NameCommand

    ! local variables

    character (len=10) :: AxesType_Array(2)
    integer            :: NAxesTypes
    ! character (len=50) :: CorotationType
    character (len=50) :: CoordinateSystemString, AxesType

    select case(NameCommand)

    case("#AXES")
       call read_var('CoordinateSystemString',CoordinateSystemString) 
       if (index(CoordinateSystemString,'GSM') <= 0 .and. &
            index(CoordinateSystemString,'gsm') <= 0) then
          write(*,*) "----------------------------------------"
          write(*,*) "              Warning:                  "
          write(*,*) "BATSRUS only works in GSM at this point."
          write(*,*) "----------------------------------------"
          if(UseStrict)call stop_mpi('Correct PARAM.in')              
          CoordinateSystem = GSM
       else                                                          
          CoordinateSystem = GSM
       endif

       call read_var('AxesType',AxesType)

       call split_str(AxesType,2,AxesType_Array,NAxesTypes)
       if (NAxesTypes == 1) AxesType_Array(2) = AxesType_Array(1)

       ! Determine Dipole Information First

       if (index(AxesType_Array(1),'none') > 0) then
          Bdp_dim = 0.0
          SetDipoleTilt = .true.
          ThetaTiltDeg = 0.0
          dt_UpdateB0 = -1.0
       elseif (index(AxesType_Array(1),'planet') > 0) then
          ! This means that the code will automatically set
          ! the dipole tilt to whatever it needs to be:
          SetDipoleTilt = .false.
          ThetaTiltDeg = 0.0  ! Loaded for corotation "ideal" case
          ! only. Overwritten later.
          dt_UpdateB0 = 0.0
       elseif (index(AxesType_Array(1),'ideal') > 0) then
          SetDipoleTilt = .true.
          ThetaTiltDeg = 0.0
          dt_UpdateB0 = -1.0
       elseif (index(AxesType_Array(1),'user') > 0) then
          call read_var('Bdp_dim'   ,Bdp_dim)
          call read_var('ThetaTilt' ,ThetaTiltDeg)
          SetDipoleTilt = .true.
          dt_UpdateB0   = -1.0
       endif

       if (index(AxesType_Array(2),'none') > 0 ) then
          UseRotatingBc = .false.
       elseif (index(AxesType_Array(2),'planet') > 0 ) then
          UseRotatingBc = .true.
          AlignDipoleCorotation = .false.
          SetCorotationTilt = .false.
       elseif (index(AxesType_Array(2),'ideal') > 0 ) then
          UseRotatingBc = .true.
          AlignDipoleCorotation = .true.
          CorotationTiltDeg = ThetaTiltDeg
          CorotationLonDeg = 0.0
          SetCorotationTilt = .true.
       elseif (index(AxesType_Array(2),'user') > 0 ) then
          UseRotatingBc = .true.
          call read_var('CorotationTilt',CorotationTiltDeg)
          call read_var('CorotationLon',CorotationLonDeg)
          call read_var('rot_period',rot_period_dim)
          SetCorotationTilt = .true.
       endif
    case("#DIPOLE")
       call read_var('Bdp_dim'      ,Bdp_dim)
       call read_var('SetDipoleTilt',SetDipoleTilt)
       call read_var('ThetaTiltDeg' ,ThetaTiltDeg)
       call read_var('dt_UpdateB0'  ,dt_UpdateB0)
    case("#COROTATION")
       call read_var('UseRotatingBc',UseRotatingBc)
       if (UseRotatingBc) then
          call read_var('AlignDipoleCorotation',AlignDipoleCorotation)
          if (.not.AlignDipoleCorotation) then
             call read_var('SetCorotationTilt',SetCorotationTilt)
             if (SetCorotationTilt) then
                call read_var('CorotationTilt',CorotationTiltDeg)
             endif
          else
             SetCorotationTilt = .true.
             CorotationTiltDeg = 0.0
          endif
       endif
    case default
       if(iProc==0) then
          write(*,*) 'Invalid #COMMAND, unknown. '//NameCommand
          if(UseStrict)call stop_mpi('Correct PARAM.in!')
       end if
    end select

  end subroutine read_compatible_command

  !-------------------------------------------------------------------------
  subroutine Rotate_into_Magneto(xyz, angle)
    implicit none

    real, dimension(1:3) :: xyz, temp_xyz
    real                 :: angle, sinangle, cosangle

    sinangle = sin(angle)
    cosangle = cos(angle)

    temp_xyz(1) = xyz(1)*cosangle - xyz(3)*sinangle
    temp_xyz(2) = xyz(2)
    temp_xyz(3) = xyz(1)*sinangle + xyz(3)*cosangle

    xyz = temp_xyz

  end subroutine Rotate_into_Magneto
  !-------------------------------------------------------------------------
  !
  !-------------------------------------------------------------------------
  subroutine Rotate_into_Iono(xyz, angle)
    implicit none

    real, dimension(1:3) :: xyz
    real                 :: angle, angle_tmp

    angle_tmp = -angle
    call Rotate_into_Magneto(xyz, angle_tmp)

  end subroutine Rotate_into_Iono
  !-------------------------------------------------------------------------
  !
  !-------------------------------------------------------------------------
  subroutine GSM_to_GEO(xyz)
    use ModNumConst
    use ModPhysics, ONLY : rot_period_dim
    implicit none

    real, dimension(1:3) :: xyz, xyz_out
    real                 :: Theta, Phi, Psi, Xtemp0, Ztemp0, Xtemp1

    if (.not.SetCorotationTilt) then

       Theta = -CorotationTilt
       Phi   = -Magnetic_Pole_Colat
       Psi   = -(Time_Array(4) * 3600.0 + &
            Time_Array(5) * 60.0 + &
            Time_Array(6) + &
            Time_Array(7)/1000.0) * cTwoPi/ (rot_period_dim * 3600.0) + &
            Magnetic_Pole_East_Lon

    else

       Theta = 0.0
       Phi   = -CorotationTilt
       Psi   = -CorotationLon

    endif

    Xtemp0     = xyz(1) * cos(Theta) - xyz(3) * sin(Theta)
    Ztemp0     = xyz(1) * sin(Theta) + xyz(3) * cos(Theta)

    Xtemp1     = Xtemp0 * cos(Psi) - xyz(2) * sin(Psi)
    xyz_out(2) = Xtemp0 * sin(Psi) + xyz(2) * cos(Psi)

    xyz_out(1) = Xtemp1 * cos(Phi) - Ztemp0 * sin(Phi)
    xyz_out(3) = Xtemp1 * sin(Phi) + Ztemp0 * cos(Phi)

    xyz = xyz_out

  end subroutine GSM_to_GEO
  !-------------------------------------------------------------------------
  !
  !-------------------------------------------------------------------------
  subroutine GEO_to_GSM(xyz)
    use ModNumConst
    use ModPhysics, ONLY : rot_period_dim
    implicit none

    real, dimension(1:3) :: xyz, xyz_out
    real                 :: Theta, Phi, Psi, lon, Xtemp0, Xtemp1, Ztemp0

    if (.not.SetCorotationTilt) then

       Theta = CorotationTilt
       Phi   = Magnetic_Pole_Colat

       Psi   = ((Time_Array(4) * 3600.0 + &
            Time_Array(5) * 60.0 + &
            Time_Array(6) + &
            Time_Array(7)/1000.0) * cTwoPi / (rot_period_dim * 3600.0) + &
            Magnetic_Pole_East_Lon)

    else

       Phi   = CorotationTilt
       Psi   = CorotationLon
       Theta = 0.0

    endif

    Xtemp0     = xyz(1) * cos(Phi) - xyz(3) * sin(Phi)
    Ztemp0     = xyz(1) * sin(Phi) + xyz(3) * cos(Phi)

    Xtemp1     = Xtemp0 * cos(Psi) - xyz(2) * sin(Psi)
    xyz_out(2) = Xtemp0 * sin(Psi) + xyz(2) * cos(Psi)

    xyz_out(1) = Xtemp1 * cos(Theta) - Ztemp0 * sin(Theta)
    xyz_out(3) = Xtemp1 * sin(Theta) + Ztemp0 * cos(Theta)

    xyz = xyz_out

  end subroutine GEO_to_GSM
  !-------------------------------------------------------------------------
  !\
  ! This program calculates cartesian corotation velocity v_phi as a
  ! function of the cartesian coordinates "location"
  !/
  !-------------------------------------------------------------------------
  subroutine calc_corotation_velocities(iter,time_now,location,v_phi)
    use ModPhysics, ONLY : OMEGAbody
    use ModNumConst
    implicit none

    integer,intent(in) :: iter
    integer, intent(in) :: time_now
    real,dimension(1:3),intent(in) :: location  
    real,dimension(1:3),intent(out) :: v_phi

    ! ------------------------------------------------------------------
    ! Need to find corotation velocity.  This should be the velocity
    ! around the geographic pole.  Here it is the rotation around the
    ! geomagnetic pole, which is WRONG, but close.
    ! ------------------------------------------------------------------

    call GSM_to_GEO(location)
    ! Find the rotation velocity in the XYZ rotated frame

    v_phi(1)  = -OMEGAbody*location(2)
    v_phi(2)  = OMEGAbody*location(1) 
    v_phi(3)  = cZero

    ! Rotate into real coordinate system

    call GEO_to_GSM(v_phi)

  end subroutine calc_corotation_velocities
  !----------------------------------------------------------------------------
  !  
  !  This routine calculates the dipole tilt angles as well as the
  !  corotation tilt angles needed in the code.  Here is how the angles are
  !  defined:
  !
  !  For and ideal case, the angles (is GSM) can be described as the following:
  !
  !  For the Dipole:
  !    THETAtilt   - The angle of the dipole in the X-Z plane
  !
  !  For the Corotation:
  !    CorotationTilt - The angle of the corotation tilt in the X-Z plane
  !    CorotationLon  - The angle of the corotation tilt around the Z-axis.
  !
  !  When we want to do planets, though, the angles are set as:
  !
  !  For the Dipole:
  !    THETAtilt - Because the dipole is most likely not aligned with the 
  !                rotation axis, and the rotation axis changes its angle
  !                with respect to the sun as a function of the day of the
  !                year, THETAtilt has a time of day (ToD) dependence and
  !                a day of year (DoY) dependence.  For these to work 
  !                correctly, the programmer has to know a number of different
  !                quantities, such as the offset of the magnetic and geographic
  !                poles (Magnetic_Pole_Colat), the longitude of the magnetic
  !                pole in geographic coordinates (Magnetic_Pole_East_Lon),
  !                the maximum tilt of the planet towards the sun (Max_DoY_Tilt),
  !                and one of the two days of the year in which the planet
  !                is pointing vertical with respect to the sun in the X-Z plane.
  !                The proper tilt in the X-Z plane can then be determined.
  !
  !  For the Corotation:
  !    CorotationTilt - This is set as the Day of the Year Tilt.
  !    In 'rotate.f90', the time of day angle is also set.  This is done there
  !    so that the corotation axis can move from time step to time step, unlike
  !    the dipole.  This angle uses the same angles as described above.
  !
  !----------------------------------------------------------------------------

  subroutine Calculate_Dipole_Tilt

    use ModProcMH, ONLY: iProc
    use ModMain
    use ModPhysics
    use ModNumConst
    implicit none

    real :: DoY, ToD, THETAtilt_DoY, THETAtilt_ToD

    ! Make sure that time_array is set
    call increment_real_world_time(time_simulation)

    DoY = 0.0
    ToD = 0.0

    if (.not.SetDipoleTilt) then

       select case(problem_type)                      !^CFG IF NOT SIMPLE

       case(problem_earth)                            !^CFG IF NOT SIMPLE

          ! This is an angle:

          ToD = (Time_Array(4) * 3600.0 + &
               Time_Array(5) * 60.0 + &
               Time_Array(6) + &
               Time_Array(7)/1000.0) * cTwoPi / (rot_period_dim * 3600.0) + &
               Magnetic_Pole_East_Lon 

          ! This is an angle:

          DoY = (jday(Time_Array(1),  &
               Time_Array(2),  &
               Time_Array(3)) - &
               jday(Time_Array(1), 3, 21)) * cTwoPi / 365.0

       end select                                    !^CFG IF NOT SIMPLE

       THETAtilt_DoY = -Max_DoY_Tilt * sin(DoY)
       THETAtilt_ToD = -Magnetic_Pole_Colat * sin(ToD - cHalfPi)

       THETAtilt = (THETAtilt_DoY + THETAtilt_ToD)

    else

       ! Convert to radians :
       THETAtilt = cTwoPi*THETAtiltDeg/360.00
       THETAtilt_DoY = THETAtilt
       THETAtilt_ToD = 0.0

    endif

    cosTHETAtilt = cos(THETAtilt)
    sinTHETAtilt = sin(THETAtilt)

    if (AlignDipoleCorotation) then

       ! This mean that if the user selects "planet ideal" then the corotation
       ! will update with the tilt angle.

       CorotationTilt = Thetatilt
       CorotationLon  = 0.0

    else

       if (.not.SetCorotationTilt) then
          CorotationTilt = THETAtilt_DoY
          CorotationLon  = ToD
       else
          CorotationTilt = CorotationTiltDeg * cPi / 180.0
          CorotationLon  = CorotationLonDeg * cPi / 180.0
       endif

    endif

    if (iProc == 0 .and. lVerbose>0) then
       write(*,*)"=> Setting Dipole Tilt to",ThetaTilt*180/cPi," deg."
       write(*,*)"=> Setting Corotation Tilt to",CorotationTilt*180/cPi," deg."
       write(*,*)"=> Setting Corotation Longitude to",&
            CorotationLon*180/cPi," deg."
    endif

  end subroutine Calculate_Dipole_Tilt

  !----------------------------------------------------------------
  ! These routines are for keeping track of the real world time in
  ! the code.  This is mainly used in MHD for knowing when to couple
  ! to other codes, setting the dipole tilt, etc.  In the ionosphere
  ! these routines are needed to calculate the solar zenith angles
  ! for the conductances.
  !----------------------------------------------------------------

  subroutine time_int_to_real(iTime, Time)

    implicit none

    integer, dimension(1:7), intent(in) :: iTime
    real*8, intent(out) :: Time

    integer, dimension(1:12) :: nDaysInMonth
    integer :: iYear, nLeapYears, iMonth, iDay, iHour, iMinute, iSecond, i

    nDaysInMonth(1) = 31
    nDaysInMonth(2) = 28
    nDaysInMonth(3) = 31
    nDaysInMonth(4) = 30
    nDaysInMonth(5) = 31
    nDaysInMonth(6) = 30
    nDaysInMonth(7) = 31
    nDaysInMonth(8) = 31
    nDaysInMonth(9) = 30
    nDaysInMonth(10) = 31
    nDaysInMonth(11) = 30
    nDaysInMonth(12) = 31

    if (mod(iTime(1),4).eq.0) nDaysInMonth(2) = 29

    Time = 0.0;

    if (iTime(1) > 1900) then
       iYear = iTime(1) - 1965
    else 
       if (iTime(1) > 65) then
          iYear = iTime(1) - 65 
       else 
          iYear = iTime(1) + 100 - 65
       endif
    endif
    nLeapYears = iYear/4

    iMonth = iTime(2) - 1

    iDay = 0

    do i=1, iMonth
       iDay = iDay + nDaysInMonth(i)
    enddo

    iDay = iDay + iTime(3) - 1
    iHour = iTime(4)
    iMinute = iTime(5)
    iSecond = iTime(6)

    Time = (dble(iSecond) * dble(1.0)) +                  &
         (dble(iMinute) * dble(60.0)) +                       &
         (dble(iHour) * dble(60.0*60.0)) +                 &
         (dble(iDay) * dble(24.0*60.0*60.0)) +             &
         (dble(nLeapYears) * dble(24.0*60.0*60.0)) +            &
         (dble(iYear) * dble(365.0*24.0*60.0*60.0)) +      &
         iTime(7)/1000.0

  end subroutine time_int_to_real


  subroutine time_real_to_int(Time, iTime)

    implicit none

    real*8, intent(in) :: Time
    integer, dimension(1:7), intent(out) :: iTime

    integer, dimension(1:12) :: nDaysInMonth
    integer :: iYear, nLeapYears, iMonth, iDay, iHour, iMinute, iSecond
    real*8 :: SecondsPerYear = 31536000.0
    real*8 :: SecondsPerDay = 86400.0
    real*8 :: SecondsPerHour = 3600.0
    real*8 :: SecondsPerMinute = 60.0
    real*8 :: TimeRemaining

    nDaysInMonth(1) = 31
    nDaysInMonth(2) = 28
    nDaysInMonth(3) = 31
    nDaysInMonth(4) = 30
    nDaysInMonth(5) = 31
    nDaysInMonth(6) = 30
    nDaysInMonth(7) = 31
    nDaysInMonth(8) = 31
    nDaysInMonth(9) = 30
    nDaysInMonth(10) = 31
    nDaysInMonth(11) = 30
    nDaysInMonth(12) = 31

    iYear = int(Time/SecondsPerYear)
    nLeapYears = iYear/4
    iDay = int((Time - (dble(iYear)*SecondsPerYear))/SecondsPerDay)

    if (iDay.le.nLeapYears) then
       iYear = int((Time - (dble(nLeapYears)*SecondsPerDay))/SecondsPerYear)
       nLeapYears = iYear/4
       iDay = int((Time - (dble(iYear)*SecondsPerYear))/SecondsPerDay)
       if (iDay.le.nLeapYears) then
          iYear = int((Time - (dble(nLeapYears)*SecondsPerDay))/SecondsPerYear)
          nLeapYears = iYear/4
          iDay = int((Time - (dble(iYear)*SecondsPerYear))/SecondsPerDay)
       endif
    endif

    if (mod((iYear+65),4).eq.0) nDaysInMonth(2) = nDaysInMonth(2) + 1

    iDay = iDay - nLeapYears

    TimeRemaining = Time - dble(iYear)*SecondsPerYear
    TimeRemaining = TimeRemaining - dble(iDay+nLeapYears)*SecondsPerDay

    iHour = int(TimeRemaining/SecondsPerHour)
    TimeRemaining = TimeRemaining - dble(iHour)*SecondsPerHour

    iMinute = int(TimeRemaining/SecondsPerMinute)
    TimeRemaining = TimeRemaining - dble(iMinute)*SecondsPerMinute

    iSecond = int(TimeRemaining)

    iMonth = 1;

    do while (iDay.ge.nDaysInMonth(iMonth))
       iDay = iDay - nDaysInMonth(iMonth)
       iMonth = iMonth + 1
    end do

    iTime(1) = iYear + 1965
    iTime(2) = iMonth
    iTime(3) = iDay + 1
    iTime(4) = iHour
    iTime(5) = iMinute
    iTime(6) = iSecond
    iTime(7) = (TimeRemaining - iSecond)*1000

  end subroutine time_real_to_int

  integer function jday(year, mon, day) result(Julian_Day)

    implicit none

    integer :: i
    integer, dimension(1:12) :: nDaysInMonth
    integer :: year, mon, day

    nDaysInMonth(1) = 31
    nDaysInMonth(2) = 28
    nDaysInMonth(3) = 31
    nDaysInMonth(4) = 30
    nDaysInMonth(5) = 31
    nDaysInMonth(6) = 30
    nDaysInMonth(7) = 31
    nDaysInMonth(8) = 31
    nDaysInMonth(9) = 30
    nDaysInMonth(10) = 31
    nDaysInMonth(11) = 30
    nDaysInMonth(12) = 31

    if (mod(year,4).eq.0) nDaysInMonth(2) = nDaysInMonth(2) + 1
    Julian_Day = 0
    do i = 1, mon-1
       Julian_Day = Julian_Day + nDaysInMonth(i)
    enddo
    Julian_Day = Julian_Day + day

  end function jday


  subroutine increment_real_world_time(dt_in)

    use ModMain, ONLY : StartTime

    implicit none

    real, intent(in) :: dt_in

    call time_real_to_int(StartTime + dble(dt_in), Time_Array)

  end subroutine increment_real_world_time

end module ModCompatibility
