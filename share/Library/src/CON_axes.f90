!^CFG COPYRIGHT UM
!BOP
!MODULE: CON_axes - coordinate system initialization and setting
!INTERFACE:
module CON_axes

  !DESCRIPTION:
  ! CON uses GSE coordinates for planetary data, because it is convenient 
  !     as well as it is inertial except for orbital motion. It connects
  !     the planet and the Sun with the X axis, which makes it the 
  !     ideal choice for describing the whole space weather simulations.
  !
  ! \bigskip
  !
  ! {\bf Coordinate system definitions for SWMF}
  !
  ! \bigskip
  !
  ! Coordinate systems with their origin in the center of the planet:
  ! \begin{verbatim}
  ! GEI (Geocentric Equatorial Inertial)
  ! PEI (Planetocentric Equatorial Inertial)
  !
  !   Z parallel with the rotation axis.
  !   X points towards the vernal equinox.
  !   Y completes the right handed coordinate system.
  !   GEI orbits around the Sun.
  !   GEI does not rotate except for the precession and nutation of the 
  !       rotation axis of the Earth.
  !   Inertial forces are negligible.
  !
  ! GSE (Geocentric Solar Ecliptic) 
  ! PSO (Planetocentric Solar Orbital)
  !
  !   X towards the Sun (S)
  !   Z orthogonal to Orbital/Ecliptic plane pointing "North"
  !   Y opposite of orbital velocity
  !
  !   GSE is rotating around the Z axis with the orbital angular speed.
  !   GSE orbits around the Sun. 
  !   The inertial forces can be neglected.

  ! GSM (Geocentric Solar Magnetic) 
  ! PSM (Planetocentric Solar Magnetic)
  !
  !   X axis points towards the Sun (as GSE)
  !   Z axis is rotated such that the magnetic axis lies in the X-Z plane
  !   Y completes the right handed coordinate system
  !
  !   GSM is rotating around the Z axis with the orbital angular speed.
  !   GSM is rotating around the X axis back and forth with the projection
  !       of the magnetic axis motion, which depends on the rotational angular
  !       speed and on the angle between the magnetic and rotational axes.
  !   GSM orbits around the Sun.
  !   Inertial forces can be neglected if the magnetic and rotational axes are
  !   (almost) aligned and/or the rotation speed is slow relative to dynamical
  !   time scales.

  ! SMG (Solar MaGnetic Coordinates)
  !
  !   Z is the magnetic axis pointing "North".
  !   Y is orthogonal to the direction to the Sun.
  !   X completes the right handed cooridinate system with X 
  !     pointing "towards" the Sun.
  !
  !   SMG wobbles around due to the rotation of the magnetic axis.
  !   SMG orbits around the Sun.
  !   SMG differs from GSM in a rotation around the Y axis.
  !
  !   Inertial forces can be neglected if the magnetic and rotational axes are
  !   (almost) aligned and/or the rotation speed is slow relative to dynamical
  !   time scales.

  ! GEO (GEOgraphic) or PGR (PlanetoGRaphic)
  !
  !   Z is the rotation axis pointing "North".
  !   X goes through the 0 meridian which is defined for the planet.
  !   Y completes the right handed coordinate system.
  !
  !   GEO is a corotating coordinate system.
  !   GEO rotates around the Z axis with the inertial angular speed 
  !       of the planet (which is NOT 2*Pi/(24*3600.0) for the Earth).
  !       For the Earth the 0 meridian goes through Greenich. For other planets
  !       the 0 meridian is defined as the half plane which is AngleEquinox 
  !       away from the direction of equinox at the time of equinox.
  !   GEO orbits around the Sun.
  !   Inertial forces may or may not be negligible.

  ! MAG (Magnetic coordinates)
  !
  !   Z is the magnetic axis pointing "North".
  !   Y axis is orthogonal to rotational axis pointing towards Omega x Z
  !   X completes the right handed coordinate system.
  !
  !   MAG rotates around the rotational axis which does not coincide with 
  !       any of the principal axes.
  !   MAG orbits around the Sun.
  !
  !   Inertial forces may or may not be negligible.

  ! HGI (HelioGraphic Inertial coordinates)
  !   Z is the rotation axis of the Sun pointing "North".
  !   X axis is the intersection of the ecliptic and solar equatorial planes. 
  !     which is at around 75 degrees ecliptic longitude.
  !   Y axis completes the right handed coordinate system.
  !
  !   HGI is a truly inertial system.

  ! HGR (HelioGRaphic coordinates)
  !   
  !   Z is the rotation axis of the Sun pointing "North".
  !   X axis rotates with the Carrington rotation with a 27.2753 day period.
  !     The X axis coincided with the X axis of the HGI system on 
  !     09 Nov. 1853 00:00:00
  !   Y axis completes the right handed coordinate system.
  !
  ! HGR is a rotating system, inertial forces should be taken into account. 
  !
  !\end{verbatim}
  !TODO:
  ! Generalize transformations to and from heliocentric coordinate systems 
  ! for non-Earth planets.

  !USES:

  use ModKind
  use ModConst
  use ModCoordTransform
  use ModTimeConvert, ONLY : time_int_to_real,time_real_to_int
  use CON_planet
  use CON_geopack, ONLY: &
       HgiGse_DD, CON_recalc, CON_sun, SunEMBDistance, JulianDay

  !REVISION HISTORY:
  ! 01Aug03 - Gabor Toth and Aaron Ridley  - initial version
  ! 14Aug03 - Gabor Toth <gtoth@umich.edu> - major revision and extension
  ! 23Mar04 - Gabor Toth eliminated the use of CON_time to make 
  !                      CON_axes more indepenedent
  ! 17Jan05 - Ofer Cohen and G. Toth merged in GEOPACK and added functions
  !                      angular_velocity and transform_velocity
  !EOP

  implicit none

  save

  character(len=*), parameter, private :: NameMod='CON_axes'

  integer, parameter, private :: x_=1, y_=2, z_=3

  ! Difference between 01/01/1965 00:00:00 and 09/11/1853 00:00:00 in seconds
  real(Real8_), parameter :: tStartCarrington   = -3.5074080D+09
  real(Real8_), parameter :: CarringtonRotation = cTwoPi/(27.2753D0*24*3600)

  ! Position and Velocity of Planet in HGI
  real :: XyzPlanet_D(3), vPlanetHgi_D(3)

  ! Initial time in 8 byte real
  real(Real8_) :: tStart = -1.0

  ! Rotational axis in GSE and GSM
  real    :: RotAxis_D(3)      ! Permanent Cartesian components in GSE
  real    :: RotAxisGsm_D(3)   ! Changing  Cartesian components in GSM

  ! Magnetic axis in GEO, GEI and GSE
  real    :: MagAxisGeo_D(3)                         ! Permanent vector in GEO
  real    :: MagAxis0Gei_D(3)  ! Starting position of the magnetix axis in GEI
  real    :: MagAxis_D(3)      ! Current  position of the magnetix axis in GSE
  real    :: MagAxisGsm_D(3)   ! Current  position of the magnetix axis in GSM
  real    :: MagAxisTiltGsm    ! Current  tilt  in GSM

  ! Logical tells if the time independent axis parameters have been set
  logical :: DoInitializeAxes=.true.

  ! Coordinate transformation matrices connecting all the systems
  real, dimension(3,3) :: &
       SmgGsm_DD, &            ! vSmg_D = matmul(SmgGsm_DD,vGsm_D)
       GsmGse_DD, &            ! vGsm_D = matmul(GsmGse_DD,vGse_D)
       GseGei_DD, &            ! vGse_D = matmul(GseGei_DD,vGei_D)
       GeiGeo_DD, &            ! vGei_D = matmul(GeiGeo_DD,vGeo_D)
       MagGeo_DD, &            ! vMag_D = matmul(MagGeo_DD,vGeo_D)
       HgrHgi_DD, &            ! vHgr_D = matmul(HgrHgi_DD,vHgi_D)
       HgrGse_DD               ! vHgr_D = matmul(HgrGse_DD,vGse_D)

  ! Remaining coordinate transformation matrices to convert to/from GSE
  real, dimension(3,3) :: &
       SmgGse_DD, GeoGse_DD, MagGse_DD

contains

  !BOP ========================================================================
  !IROUTINE: init_axes - initialize the axes
  !INTERFACE:
  subroutine init_axes(tStartIn)

    !INPUT ARGUMENTS:
    real(Real8_) :: tStartIn

    !DESCRIPTION:
    ! Set the direction and position of the rotation axis in GSE
    ! Set the initial direction and position of the magnetic axis in 
    ! GSE, GSM, GEI and GEO systems.
    !
    ! Calculate conversion matrices between MAG-GEO-GEI-GSE systems.
    !EOP

    character(len=*), parameter :: NameSub=NameMod//'::init_axes'

    logical :: DoTest, DoTestMe
    !-------------------------------------------------------------------------

    if (.not.DoInitializeAxes) return

    call CON_set_do_test(NameSub, DoTest,DoTestMe)

    tStart = tStartIn

    call time_int_to_real(TimeEquinox)

    ! Get GSE position for the rotational axis
    if(.not.UseSetRotAxis)then
       if(UseRealRotAxis .or. UseRealMagAxis)then
          RotAxisTheta = TiltRotation
          RotAxisPhi   = mod( &
               cHalfPi - OmegaOrbit*(tStart - TimeEquinox % Time), &
               real(cTwoPi,kind=Real8_))
       else
          ! Rotational axis must be aligned with magnetic axis
          if(UseSetMagAxis)then
             RotAxisTheta = MagAxisTheta
             RotAxisPhi   = MagAxisPhi
          else
             call CON_stop(NameSub// &
                  ' SWMF_ERROR both rotation and magnetic axes'//&
                  ' are aligned with the other one?!')
          end if
       end if
    end if

    if(DoTestMe)then
       write(*,*)'tStart,TimeEquinox=',tStart,TimeEquinox
       write(*,*)'RotAxisTheta,RotAxisPhi=',&
            RotAxisTheta*cRadToDeg, RotAxisPhi*cRadToDeg
    end if

    ! Using the RotAxisTheta and RotAxisPhi 
    ! set the GseGei matrix to convert between GSE and  GEI systems
    call set_gse_gei_matrix

    ! Calculate initial position for the magnetic axis in GSE and GEI systems
    if(UseRealMagAxis)then
       ! Cartesian coordinates of the magnetic axis unit vector in GEO
       call dir_to_xyz(MagAxisThetaGeo,MagAxisPhiGeo,MagAxis_D)

       if(DoTestMe)then
          write(*,*)'MagAxisThetaGeo,MagAxisPhiGeo=',&
               MagAxisThetaGeo*cRadToDeg,MagAxisPhiGeo*cRadToDeg
          write(*,*)'MagAxisGeo_D=',MagAxis_D
       end if

       ! GEO --> GEI
       call set_gei_geo_matrix(0.0)
       MagAxis0Gei_D = matmul(GeiGeo_DD,MagAxis_D)

       ! GEI --> GSE
       MagAxis_D = matmul(GseGei_DD,MagAxis0Gei_D)

       ! Cartesian vector to spherical direction
       call xyz_to_dir(MagAxis_D,MagAxisTheta,MagAxisPhi)

       if(DoTestMe)then
          write(*,*)'UseRealMagAxis:'
          write(*,*)'MagAxisGei_D=',MagAxis0Gei_D
          write(*,*)'GseGei_DD='
          call show_rot_matrix(GseGei_DD)
          write(*,*)'MagAxisGse_D=',MagAxis_D
          write(*,*)'MagAxisTheta,MagAxisPhi=',&
               MagAxisTheta*cRadToDeg,MagAxisPhi*cRadToDeg
       end if

    else
       if(.not.UseSetMagAxis)then
          ! Must be aligned with rotational axis
          MagAxisTheta = RotAxisTheta
          MagAxisPhi   = RotAxisPhi
       end if
       ! Convert direction to Cartesian coordinates in GSE
       call dir_to_xyz(MagAxisTheta,MagAxisPhi,MagAxis_D)

       ! Calculate the GEI position too 
       ! (in case mag axis is not aligned and rotates)
       call set_gei_geo_matrix(0.0)
       MagAxis0Gei_D = matmul(MagAxis_D,GseGei_DD)

       if(DoTestMe)then
          write(*,*)'Aligned=',.not.UseSetMagAxis,' Set=',UseSetMagAxis
          write(*,*)'MagAxisGei_D=',MagAxis0Gei_D
          write(*,*)'MagAxisGse_D=',MagAxis_D
          write(*,*)'MagAxisTheta,MagAxisPhi=',&
               MagAxisTheta*cRadToDeg,MagAxisPhi*cRadToDeg
       end if

    end if

    if(.not.(UseRealRotAxis .or. UseSetRotAxis) .and. UseRealMagAxis)then
       ! Rotation axis is aligned.
       ! The angles are reset now because we needed the real rotational axis 
       ! to find the real magnetic axis. We do not need that anymore.
       RotAxisTheta = MagAxisTheta
       RotAxisPhi   = MagAxisPhi
       ! The "permanent" matrices are also recalculated
       call set_gse_gei_matrix
       call set_gei_geo_matrix(0.0)
    endif

    ! Obtain the cartesian components of the rotational axis (in GSE)
    call dir_to_xyz(RotAxisTheta,RotAxisPhi,RotAxis_D)

    if(.not.(UseRealRotAxis.and.UseRealMagAxis))then
       ! Recalculate the magnetic axis direction in GEO

       MagAxisGeo_D = matmul(MagAxis0Gei_D,GeiGeo_DD)

       call xyz_to_dir(MagAxisGeo_D, MagAxisThetaGeo, MagAxisPhiGeo)

       if(DoTestMe)write(*,*)'Final MagAxisThetaGeo, MagAxisPhiGeo=',&
            MagAxisThetaGeo*cRadToDeg, MagAxisPhiGeo*cRadToDeg
    else
       ! Set the magnetic direction in Cartesian GEO coordinates
       call dir_to_xyz(MagAxisThetaGeo, MagAxisPhiGeo, MagAxisGeo_D)
    end if

    ! From MagAxisThetaGeo and MagAxisPhiGeo obtain the MAG-GEO matrix
    ! This matrix does not change with simulation time.
    call set_mag_geo_matrix

    ! Set the time dependent axes for the initial time
    call set_axes(0.0,.true.)

    if(DoTestMe)then
       write(*,*)'Final rotation axis:'
       write(*,*)'RotAxisTheta,RotAxisPhi=',&
            RotAxisTheta*cRadToDeg, RotAxisPhi*cRadToDeg
       write(*,*)'RotAxisGse_D=',RotAxis_D
       write(*,*)'RotAxisGsm_D=',RotAxisGsm_D
    end if

    DoInitializeAxes=.false.

    ! Set Hgi-Gse matrix and planet velocity 
    call set_hgi_gse_v_planet

  contains

    !=========================================================================

    subroutine set_gse_gei_matrix

      ! The GseGei_DD matrix converts between GSE and GEI with two rotations:
      !
      !   rotate around X_GEI with RotAxisTheta      so that Z_GEI->Z_GSE
      !   rotate around Z_GSE with RotAxisPhi + Pi/2 so that X_GEI->X_GSE
      !
      ! The GseGei_DD matrix changes at the order of TimeSimulation/TimeOrbit.
      ! For usual simulations that change can be safely neglected.

      !-----------------------------------------------------------------------

      GseGei_DD = matmul(&
           rot_matrix_z(RotAxisPhi + cHalfPi), &
           rot_matrix_x(RotAxisTheta) &
           )

    end subroutine set_gse_gei_matrix

    !=========================================================================

    subroutine set_mag_geo_matrix

      ! The first rotation is around the Z_GEO axis with MagAxisPhiGeo,
      ! which rotates Y_GEO into Y_MAG.
      ! The second rotation is around the Y_MAG axis with MagAxisThetaGeo,
      ! which rotates Z_GEO into Z_MAG.
      !
      ! This matrix only changes with the slow motion of the magnetix axis
      ! relative to the Earth.

      MagGeo_DD = matmul( &
           rot_matrix_y(-MagAxisThetaGeo), &
           rot_matrix_z(-MagAxisPhiGeo))

    end subroutine set_mag_geo_matrix

    !=========================================================================

    subroutine set_hgi_gse_d_planet(tSimulation)
      ! Calculate HgiGse matrix from CON_recalc in CON_geopack
      real, intent(in) :: tSimulation

      integer :: iTime_I(1:7)
      integer::iYear,iMonth,iDay,iHour,iMin,iSec,jDay
      real :: GSTime, SunLongitude, Obliq
      
      !-----------------------------------------------------------------------
      call time_real_to_int(tStart + tSimulation, iTime_I)
      iYear=iTime_I(1);iMonth=iTime_I(2);iDay=iTime_I(3)
      iHour=iTime_I(4);iMin=iTime_I(5);iSec=iTime_I(6)
      call CON_recalc(iYear,iMonth,iDay,iHour,iMin,iSec)
      jDay = JulianDay(iYear,iMonth,iDay)
      call CON_sun(iYear,jDay,iHour,iMin,iSec,GSTime,SunLongitude,Obliq)

    end subroutine set_hgi_gse_d_planet

    !=========================================================================

    subroutine set_hgi_gse_v_planet
      ! Caculate vPlanet in HGI system
      real, dimension(3) :: XyzPlus_D
      !----------------------------------------------------------------------
      ! Calculate planet position for TimeSim and TimeSim+dt
      call set_hgi_gse_d_planet(1.0)
      XyzPlus_D = matmul(HgiGse_DD, (/-cAU*SunEMBDistance, 0.0, 0.0/))
      
      call set_hgi_gse_d_planet(0.0)
      XyzPlanet_D = matmul(HgiGse_DD, (/-cAU*SunEMBDistance, 0.0, 0.0/))
      
      ! Finite difference velocity with the 1 second time perturbation
      vPlanetHgi_D = XyzPlus_D - XyzPlanet_D
      
    end subroutine set_hgi_gse_v_planet

    !=========================================================================

  end subroutine init_axes

  !============================================================================

  subroutine set_gei_geo_matrix(TimeSim)

    ! The rotation is around the Z axis, which is the rotational axis
    !
    ! This matrix only changes due to the precession of Earth.

    real, intent(in) :: TimeSim

    real :: AlphaEquinox

    !-------------------------------------------------------------------------
    if(.not.UseRotation)then
       ! If the planet does not rotate we may take GEI=GEO
       GeiGeo_DD = cUnit_DD
       RETURN
    end if

    AlphaEquinox = (TimeSim + tStart - TimeEquinox % Time) &
         * OmegaPlanet + AngleEquinox

    GeiGeo_DD = rot_matrix_z(AlphaEquinox)

  end subroutine set_gei_geo_matrix

  !BOP ========================================================================
  !IROUTINE: set_axes - set time dependent axes and transformation matrices
  !INTERFACE:
  subroutine set_axes(TimeSim,DoSetAxes)

    !INPUT ARGUMENTS:
    real,              intent(in) :: TimeSim
    logical, optional, intent(in) :: DoSetAxes

    !DESCRIPTION:
    ! The magnetic axis as well as the corotating GEO and MAG frames
    ! are rotating around the rotational axis with OmegaPlanet 
    ! angular speed. Calculate the position of the axes and the
    ! transformation matrices for the given simulation time TimeSim.
    !
    ! When the optional DoSetAxes argument is present (its value is ignored), 
    ! the magnetic axis and the related variables are always set. 
    ! This is needed for the initial setting.
    !
    ! Otherwise the update is based on a number of parameters.
    !
    ! If the planet does not rotate, or the magnetic axis is aligned with 
    ! the rotation axis, no calculation is performed.
    !
    ! The last simulation time with which an update was done
    ! is stored into TimeSimLast. 
    !
    ! If DoUpdateB0 == .false. the magnetic axis is taken to be fixed.
    ! If DtUpdateB0 <= 0.001 (sec) the magnetic axis is updated
    !      if TimeSim differs from TimeSimLast.
    ! If DtUpdateB0 >  0.001 (sec) then the magnetic axis is updated
    !      if int(TimeSim/DtUpdateB0) differs from int(TimeSimLast/DtUpdateB0)
    !
    !EOP

    real :: MagAxisGei_D(3)

    character(len=*), parameter :: NameSub=NameMod//'::set_axes'

    real :: TimeSimLast = -1000.0 ! Last simulation time

    logical :: DoTest, DoTestMe
    !-------------------------------------------------------------------------
    if(.not.present(DoSetAxes))then
       ! If magnetic axis does not move, no need to update
       if(.not.DoUpdateB0) RETURN

       ! If DtUpdateB0 is more than 0.001 update if int(time/DtUpdateB0) differ
       if(DtUpdateB0 > 0.001)then
          if(int(TimeSim/DtUpdateB0) == int(TimeSimLast/DtUpdateB0)) RETURN
       end if

       ! If DtUpdateB0 is less than 1 msec update unless time is the same
       if(abs(TimeSim - TimeSimLast) < cTiny) RETURN

    end if

    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTestMe)then
       write(*,*) NameSub,'UseAlignedAxes,UseRotation,DoUpdateB0=',&
            UseAlignedAxes,UseRotation,DoUpdateB0
       write(*,*) NameSub,'DtUpdateB0,TimeSim,TimeSimLast=',&
            DtUpdateB0,TimeSim,TimeSimLast
    end if

    ! Remember the simulation time
    TimeSimLast = TimeSim

    ! Rotate MagAxis0Gei around Z axis to get current position in GEI 
    MagAxisGei_D = matmul(rot_matrix_z(OmegaPlanet*TimeSim),MagAxis0Gei_D)

    ! Transform from GEI to GSE
    MagAxis_D = matmul(GseGei_DD,MagAxisGei_D)

    ! Set the angles in GSE
    call xyz_to_dir(MagAxis_D, MagAxisTheta, MagAxisPhi)

    ! Set the transformation matrices

    ! Calculate the rotation matrix to convert between GSE and GSM systems.
    ! This is a rotation around the shared X axis with the angle between the
    ! Z_GSE axis and the magnetic axis projected onto the Y-Z plane.
    !
    ! This matrix changes with simulation time unless 
    !    UseRotation=.false. or UseAlignedAxes=.true.

    GsmGse_DD = rot_matrix_x( atan2(MagAxis_D(y_),MagAxis_D(z_)) )

    ! Calculate the rotation matrix to convert between SMG and GSM systems.
    ! This is a rotation around the Y axis with the magnetic tilt_GSM,
    ! which is -asin(MagAxis_D(x_)
    !
    ! This matrix changes with simulation time unless 
    !    UseRotation=.false. or UseAlignedAxes=.true.

    MagAxisTiltGsm   = -asin(MagAxis_D(x_))
    SmgGsm_DD = rot_matrix_y( MagAxisTiltGsm )

    ! SMG-GSE transformation matrix

    SmgGse_DD = matmul(SmgGsm_DD, GsmGse_DD)

    ! Calculate GSM coordinates and tilt of the magnetic axis.
    ! and calculate the rotation axis in GSM coordinates.
    ! These are useful to obtain the dipole field and the corotation velocity
    ! in the GSM system.
    MagAxisGsm_D     = matmul(GsmGse_DD,MagAxis_D)
    RotAxisGsm_D     = matmul(GsmGse_DD,RotAxis_D)

    ! Now calculate the transformation matrices for the rotating systems
    call set_gei_geo_matrix(TimeSim)

    GeoGse_DD = transpose(matmul(GseGei_DD,GeiGeo_DD))
    MagGse_DD = matmul(MagGeo_DD,GeoGse_DD)

    if(DoTestMe)then
       write(*,*)NameSub,' new MagAxis_D     =',MagAxis_D
       write(*,*)NameSub,' new MagAxisTiltGsm=',MagAxisTiltGsm*cRadToDeg
       write(*,*)NameSub,' new RotAxisGsm_D  =',RotAxisGsm_D
    end if

    ! Calculate HgiHgr_DD  and HgrGse_DD Matrices
    HgrHgi_DD = rot_matrix_z( &
         (TimeSim + tStart - tStartCarrington)*CarringtonRotation )
    HgrGse_DD = matmul(HgrHgi_DD, HgiGse_DD)

  end subroutine set_axes

  !BOP ========================================================================
  !IROUTINE: get_axes - get parameters of axes at a given time
  !INTERFACE:
  subroutine get_axes(TimeSim, &
       MagAxisTiltGsmOut, RotAxisGsmOut_D, RotAxisGseOut_D)

    !INPUT ARGUMENTS:
    real, intent(in) :: TimeSim
    !OUTPUT ARGUMENTS:
    real, intent(out), optional :: MagAxisTiltGsmOut
    real, intent(out), optional :: RotAxisGsmOut_D(3)
    real, intent(out), optional :: RotAxisGseOut_D(3)
    !DESCRIPTION:
    ! Provides various information about the rotation and magnetic axes
    ! through the optional output arguments.
    !EOP
    character(len=*), parameter :: NameSub=NameMod//'::get_axes'
    !--------------------------------------------------------------------------
    ! Set time independent information
    if(DoInitializeAxes)&
         call CON_stop(NameSub//' ERROR: init_axes has not been called')

    ! Set time dependent information (TimeSim is cashed)
    call set_axes(TimeSim)

    if (present(MagAxisTiltGsmOut)) MagAxisTiltGsmOut = MagAxisTiltGsm
    if (present(RotAxisGsmOut_D))   RotAxisGsmOut_D   = RotAxisGsm_D
    if (present(RotAxisGseOut_D))   RotAxisGseOut_D   = RotAxis_D

  end subroutine get_axes

  !BOP ========================================================================
  !IROUTINE: transform_matrix - return transform matrix between 2 coord systems
  !INTERFACE:
  function transform_matrix(TimeSim,TypeCoordIn,TypeCoordOut) result(Rot_DD)

    !INPUT ARGUMENTS:
    real,             intent(in) :: TimeSim      ! Simulation time
    character(len=*), intent(in) :: TypeCoordIn  ! Type of input coord. system
    character(len=*), intent(in) :: TypeCoordOut ! Type of output coord. system

    !RETURN VALUE:
    real :: Rot_DD(3,3)

    !DESCRIPTION:
    ! Calculate the transformation matrix between two coordinate systems.
    ! One should store the transformation matrix and reuse it, because
    ! this general routine is not very efficient. Typical usage:
    ! \begin{verbatim}
    ! real :: IeUa_DD(3,3)
    ! ! Obtain the transformation matrix for the current time
    ! IeUa_DD = transform_matrix(TimeSimulation,'GEO','SMG')
    ! ! transform vectors in UA (GEO system) to IE (SMG system):
    ! VecIe_D = matmul(IeUa_DD,VecUa_D)
    ! ...
    ! \end{verbatim}
    !EOP
    character(len=*), parameter :: NameSub=NameMod//'::transform_matrix'

    real :: InGse_DD(3,3), OutGse_DD(3,3)
    !------------------------------------------------------------------------
    if(TypeCoordIn == TypeCoordOut)then
       Rot_DD = cUnit_DD
       RETURN
    end if

    ! Set time dependent information
    call set_axes(TimeSim)

    select case(TypeCoordIn)
    case('GSE')
       InGse_DD = cUnit_DD
    case('GSM')
       InGse_DD = GsmGse_DD
    case('SMG')
       InGse_DD = SmgGse_DD
    case('MAG')
       InGse_DD = MagGse_DD
    case('GEO')
       InGse_DD = GeoGse_DD
    case('GEI')
       InGse_DD = transpose(GseGei_DD)
    case('HGI')
       InGse_DD = HgiGse_DD
    case('HGR')
       InGse_DD = HgrGse_DD
    case default
       call CON_stop(NameSub//' unknown TypeCoordIn='//TypeCoordIn)
    end select

    select case(TypeCoordOut)
    case('GSE')
       OutGse_DD = cUnit_DD
    case('GSM')
       OutGse_DD = GsmGse_DD
    case('SMG')
       OutGse_DD = SmgGse_DD
    case('MAG')
       OutGse_DD = MagGse_DD
    case('GEO')
       OutGse_DD = GeoGse_DD
    case('GEI')
       OutGse_DD = transpose(GseGei_DD)  
    case('HGI')
       OutGse_DD = HgiGse_DD
    case('HGR')
       OutGse_DD = HgrGse_DD
    case default
       call CON_stop(NameSub//' unknown TypeCoordOut='//TypeCoordOut)
    end select

    Rot_DD = matmul(OutGse_DD,transpose(InGse_DD))
    
  end function transform_matrix

  !BOP ========================================================================
  !IROUTINE: angular_velocity - get angular velocity between two coord systems
  !INTERFACE:

  function angular_velocity(TimeSim, NameCoord1, NameCoord2In, iFrameIn) &
       result(Omega_D)

    !INPUT ARGUMENTS:
    real,                       intent(in) :: TimeSim      ! Simulation time
    character(len=*),           intent(in) :: NameCoord1   ! Name of 1st coord. system
    character(len=*), optional, intent(in) :: NameCoord2In ! Name of 2nd coord. system
    integer, optional,intent(in) :: iFrameIn               ! Frame for result
    
    !RETURN VALUE:
    real :: Omega_D(3) ! Angular velocity components
    
    !DESCRIPTION:
    ! This subroutine calculates the 3 components of the angular velocity between
    ! two coordinate systems from the transformation matrix between them.
    ! If the second frame is not present in the argument list, the result is
    ! the angular velocity of the first frame relative to an inertial frame.
    ! The angular velocity is given in the moving frame.
    ! When both frames are given, the relative angular rotation is returned.
    ! If iFrame is presemt. it defines whether the output angular velocity is with 
    ! respect to the first (iFrame=1) or second (iFrame=2) system.
    ! If the iFrame argument is not present, the result is in the first frame.
    !EOP
    
    ! Local variables
    character (len=3) :: NameCoord2
    integer ::  iFrame
    real    ::  dTimeSim
    real, dimension(3,3) ::  Rot_DD, RotMid_DD, RotPlus_DD, RotMinus_DD, dRot_DD
    
    character (len=*), parameter :: NameSub = NameMod // '::get_omega'
    !--------------------------------------------------------------------------
    ! Check optional arguments and set defaults
    if(present(NameCoord2In))then
       NameCoord2 = NameCoord2In
       if(present(iFrameIn))then
          if(iFrameIn /= 1 .and. iFrameIn /=2)then
             write(*,*) NameSub, ' ERROR iFrame = ',iFrame
             call CON_stop(NameSub // ': invalid value for iFrame = 1 or 2')
          end if
          iFrame = iFrameIn
       else
          iFrame = 2 ! Default is to provide Omega_D in the output coord. system
       end if
    else
       if(NameCoord1(1:1) == 'H')then
          ! For heliocentric coordinate systems set the inertial frame to HGI
          NameCoord2 = 'HGI'
       else
          ! For geocentric systems GSE is assumed to be inertial ! Otherwise use GEI!
          NameCoord2 = 'GSE'
       end if
       iFrame = 1
    end if

    ! Determine the perturbation of time
    if(precision(TimeSim) >= 12) then 
       dTimeSim = max(1.0, 1e-6*TimeSim)
    else
       dTimeSim = max(1.0, 1e-4*TimeSim)
    end if
    
    if(NameCoord1 == NameCoord2)then
       ! Nothing to do
       Omega_D = 0.0
       RETURN
    end if
       
    RotMinus_DD = transform_matrix(TimeSim-dTimeSim,NameCoord1,NameCoord2)
    RotPlus_DD  = transform_matrix(TimeSim+dTimeSim,NameCoord1,NameCoord2)
    dRot_DD = (RotPlus_DD-RotMinus_DD)/(cTwo*dTimeSim)

    RotMid_DD = transform_matrix(TimeSim,NameCoord1,NameCoord2)
    Rot_DD  = matmul(transpose(RotMid_DD), dRot_DD)

    Omega_D = (/ Rot_DD(2,3), Rot_DD(3,1), Rot_DD(1,2) /)
    
    if(iFrame == 2) Omega_D = matmul(RotMid_DD, Omega_D)
    where(abs(Omega_D) < 1e-12) Omega_D = 0.00
    
  end function angular_velocity

  !BOP ========================================================================
  !IROUTINE: transform_velocity - transforms velocity between two coord systems
  !INTERFACE:
  
  function transform_velocity(TimeSim, v1_D, Position_D, &
       NameCoord1, NameCoord2) result(v2_D)

    !INPUT ARGUMENTS:
    real,             intent(in) :: TimeSim       ! Simulation time
    real,             intent(in) :: v1_D(3)       ! Velocity in 1st system
    real,             intent(in) :: Position_D(3) ! Position in 1st system
    character(len=3), intent(in) :: NameCoord1    ! Name of 1st coord. system
    character(len=3), intent(in) :: NameCoord2    ! Name of 2nd coord. system

    !RETURN VALUE:
    real :: v2_D(3)                                        ! v2 components

    !DESCRIPTION:
    ! This function transforms the velocity vector from one coordinate system 
    ! to another. The input position and velocity should be in SI units and 
    ! the output velocity vector is also in SI units.
    ! If the two systems have the same name, then the input and output 
    ! velocity vectors are the same.
    !EOP

    ! Local variables
    character (len=3) :: NameCoord1Last = 'XXX', NameCoord2Last = 'XXX' 
    real :: TimeSimLast = -1.0

    real, dimension(3)   :: v1Total_D, Omega12_D, vPlanet1_D
    real, dimension(3,3) :: Transform12_DD

    character (len=*), parameter :: NameSub = NameMod // '::get_omega'
    !--------------------------------------------------------------------------

    !If NameCoord1 is the same as NameCoord2 there is no transformation.
    if(NameCoord1 == NameCoord2)then
       v2_D = v1_D
       RETURN
    end if

    if(.not.(TimeSim == TimeSimLast .and. NameCoord1 == NameCoord1Last &
         .and. NameCoord1 == NameCoord2Last) )then

       ! Store current time and coordinate system names
       TimeSimLast = TimeSim
       NameCoord1Last = NameCoord1
       NameCoord2Last = NameCoord2

       ! Get transformation matrix and angular velocity between frames
       Transform12_DD = transform_matrix(TimeSim,NameCoord1,NameCoord2)
       Omega12_D      = angular_velocity(TimeSim, NameCoord1, NameCoord2)

       if(NameCoord1(1:1) == 'H' .eqv. NameCoord2(1:1) == 'H')then 
          ! Both helio-centric or both planet-centric, no planet speed added
          vPlanet1_D = 0.0

       else if(NameCoord1(1:1) /= 'H' .and. NameCoord2(1:1) == 'H')then
          ! Planet-centric --> helio-centric: add planet speed
          vPlanet1_D = +matmul(&
               transform_matrix(TimeSim,'HGI',NameCoord1), vPlanetHgi_D)
       else
          ! Helio-centric --> planet-centric: subtract planet speed
          vPlanet1_D = -matmul( &
               transform_matrix(TimeSim,'HGI',NameCoord1), vPlanetHgi_D)
       end if
    end if

    v1Total_D = v1_D + cross_product(Omega12_D, Position_D) + vPlanet1_D
    v2_D = matmul(Transform12_DD, v1Total_D)

  end function transform_velocity

  !BOP ========================================================================
  !IROUTINE: test_axes - test the CON_axes class
  !INTERFACE:
  subroutine test_axes

    !DESCRIPTION:
    ! Do some self consistency checks. Stop with an error message if
    ! test fails. Otherwise write out success.
    !EOP
    real :: MagAxisTilt
    real :: RotAxisGsm_D(3), RotAxisGeo_D(3), Rot_DD(3,3)
    real:: Omega_D(3),v1_D(3),v2_D(3)
    !-------------------------------------------------------------------------

    if(.not.DoInitializeAxes) write(*,*)'test failed: DoInitializeAxes=',&
         DoInitializeAxes,' should be true'

    call time_int_to_real(TimeEquinox)
    if(TimeEquinox % Time <= 0.0) write(*,*)'test failed: TimeEquinox =',&
         TimeEquinox,' should have a large positive double in the %Time field'

    write(*,'(a)')'Testing init_axes'
    call init_axes(TimeEquinox % Time)

    if(tStart /= TimeEquinox % Time)write(*,*)'test init_axes failed: ',&
         'tStart=',tStart,' should be equal to TimeEquinox % Time=',&
         TimeEquinox % Time

    if(DoInitializeAxes) write(*,*)'test init_axes failed: DoInitializeAxes=',&
         DoInitializeAxes,' should be fales'

    write(*,'(a)')'Testing get_axes'

    call get_axes(0.0, MagAxisTilt, RotAxisGsm_D)

    if(abs(MagAxisTilt*cRadToDeg - 8.0414272433221718) > 0.00001)write(*,*) &
         'test get_axes failed: MagAxisTilt =',MagAxisTilt*cRadToDeg,&
         ' should be 8.0414272433221718 degrees within round off error'

    if(maxval(abs(RotAxisGsm_D - (/0.0, 0.131054271126, 0.991375195382/))) &
         > 0.00001)write(*,*) 'test get_axes failed: RotAxisGsm_D =',&
         RotAxisGsm_D,' should be equal to ', &
         '(/0.0, 0.131054271126, 0.991375195382/) within round off errors'

    write(*,'(a)')'Testing transform_matrix'

    Rot_DD = transform_matrix(0.0,'GSM','GEO')
    RotAxisGeo_D = matmul(Rot_DD,RotAxisGsm_D)

    if(maxval(abs(RotAxisGeo_D - (/0.,0.,1./)))>0.0001) &
         write(*,*)'test transform_matrix failed: RotAxisGeo_D=',&
         RotAxisGeo_D,' should be (/0.,0.,1./) within round off errors'

    write(*,'(a)')'Testing show_rot_matrix'
    write(*,'(a)')'GeoGsm_DD='; call show_rot_matrix(Rot_DD)
    
    write(*,'(a)')'Testing angular_velocity'
    
    Omega_D = angular_velocity(0.0,'GSE','SMG')
    
    if(maxval(abs(Omega_D - (/1.01128925E-05,9.55662927E-06,-1.42873158E-06/)))>1e-10) &
         write(*,*)'test angular_velocity failed: Omega_D = ',Omega_D,&
         ' should be (/1.0112892E-05,9.5566292E-06,-1.4287315E-06/)within round off errors'   
    
    write(*,'(a)')'Testing transform_velocity'
    
    v2_D = transform_velocity(0.0, (/0.0, 0.0, 0.0/), (/cAU, 0.0, 0.0/), 'HGI', 'HGR')
    
    if(maxval(abs(v2_D -(/-1.99357971E+05,-3.45466075E+05,4.02175140E-12/)))>10) &
         write(*,*)'test angular_velocity failed: v2_D = ',v2_D,&
         ' should be (/-1.9935797E+05,-3.4546607E+05,4.0217514E-12/)within round off errors'

  end subroutine test_axes

end module CON_axes
