!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module BATL_geometry

  use BATL_size, ONLY: MaxDim, nDim

  implicit none

  SAVE

  private ! except

  public:: init_geometry  ! initialize the module
  public:: clean_geometry ! clean up storage
  public:: xyz_to_coord   ! convert XYZ coordinates to generalized coordinates
  public:: coord_to_xyz   ! convert generalized coordinates to XYZ coordinates
  public:: radius_to_gen  ! convert radial coordinate to generalized coordinate
  public:: gen_to_radius  ! convert generalized coordinate to radial coordinate
  public:: test_geometry  ! unit test

  character(len=20), public:: TypeGeometry = 'cartesian'

  ! Cartesian, cylindrical or spherical coordinates
  logical, public:: IsCartesianGrid   = .true.  ! Cartesian grid (possibly RZ)
  logical, public:: IsCartesian       = .true.  ! Normal Cartesian geometry
  logical, public:: IsRotatedCartesian= .false. ! Rotated Cartesian grid
  logical, public:: IsRzGeometry      = .false. ! RZ geometry (x is symm. axis)
  logical, public:: IsRoundCube       = .false. ! square/cube stretched
  logical, public:: IsCylindrical     = .false. ! cylindrical: r, phi, z
  logical, public:: IsSpherical       = .false. ! spherical: r, theta, phi
  logical, public:: IsRLonLat         = .false. ! spherical: r, lon, lat
  logical, public:: IsCylindricalAxis = .false. ! r=0 boundary for cylindrical
  logical, public:: IsSphericalAxis   = .false. ! theta=0 and pi boundaries
  logical, public:: IsLatitudeAxis    = .false. ! |lat|=pi/2 boundaries
  logical, public:: IsAnyAxis         = .false. ! true if any of the above 3 is
  logical, public:: IsLogRadius       = .false. ! logarithmic radial coordinate
  logical, public:: IsGenRadius       = .false. ! stretched radial coordinate

  ! Periodicity of the domain per dimension
  logical, public:: IsPeriodic_D(MaxDim) = .false.

  ! Index names for coordinates limited by nDim
  integer, parameter, public:: Dim1_=1, Dim2_=min(nDim,2), Dim3_=min(nDim,3)

  ! Index names for Cartesian components (not limited by nDim)
  integer, parameter, public:: x_=1, y_=2, z_=3

  ! The following index names will be set in init_geometry
  integer, public:: r_=-1, Phi_=-1, Theta_=-1, Lon_=-1, Lat_=-1

  ! General radial coordinates
  integer, public:: nRgen = -1    ! number of elements in LogRgen_I
  real,    public, allocatable:: LogRgen_I(:)  ! array of log(r) values

  ! Rotation matrix from generalized to X,Y,Z coordinates
  real, public:: GridRot_DD(MaxDim,MaxDim)

contains

  !=========================================================================

  subroutine init_geometry(TypeGeometryIn, IsPeriodicIn_D, RgenIn_I)

    use ModNumConst,       ONLY: i_DD
    use ModCoordTransform, ONLY: rot_matrix_z

    character(len=*), optional, intent(in):: TypeGeometryIn
    logical,          optional, intent(in):: IsPeriodicIn_D(nDim)
    real,             optional, intent(in):: RgenIn_I(:)

    ! Initialize geometry for BATL
    !
    ! TypeGeometry can be
    !    'cartesian'
    !    'rotatedcartesian'
    !    'rz'
    !    'roundcube'
    !    'cylindrical'
    !    'cylindrical_lnr'
    !    'cylindrical_genr'
    !    'spherical'
    !    'spherical_lnr'
    !    'spherical_genr'
    !    'rlonlat'
    !    'rlonlat_lnr'
    !    'rlonlat_genr'
    !
    ! IsPeriodic_D defines periodicity for each dimension
    !
    ! RgenIn_I defines a mapping from a general index space to the radial 
    ! coordinate. The index space is mapped to the 0-1 interval so the 
    ! first element RgenIn_I corresponds to 0.0, and the last element to 1.0. 
    ! The interpolation is done in log(R), so we get a logarithmic radial 
    ! grid within each interval.

    character(len=*), parameter:: NameSub = 'init_geometry'
    !-----------------------------------------------------------------------

    TypeGeometry = 'cartesian'
    if(present(TypeGeometryIn)) TypeGeometry = TypeGeometryIn

    IsPeriodic_D = .false.
    if(present(IsPeriodicIn_D)) IsPeriodic_D(1:nDim) = IsPeriodicIn_D

    ! Logicals are useful for efficient code
    IsCartesian        = TypeGeometry(1:9)  == 'cartesian'
    IsRotatedCartesian = TypeGeometry(1:16) == 'rotatedcartesian'
    IsRzGeometry       = TypeGeometry(1:2)  == 'rz'
    IsSpherical        = TypeGeometry(1:3)  == 'sph'
    IsRLonLat          = TypeGeometry(1:3)  == 'rlo'
    IsCylindrical      = TypeGeometry(1:3)  == 'cyl'
    IsRoundCube        = TypeGeometry(1:5)  == 'round'

    IsLogRadius   = index(TypeGeometry,'lnr')  > 0
    IsGenRadius   = index(TypeGeometry,'genr') > 0

    ! Grid is Cartesian (even in RZ geometry)
    IsCartesianGrid = IsCartesian .or. IsRzGeometry

    ! Set up a rotation matrix
    if(IsRotatedCartesian)then
       ! Rotate around the Z axis with atan(3/4)
       GridRot_DD = rot_matrix_z(0.6, 0.8)
    else
       ! Just in case it would be used
       GridRot_DD = i_DD
    end if

    r_ = -1; Phi_ = -1; Theta_ = -1; Lon_ = -1; Lat_ = -1
    if(IsRzGeometry)then
       r_ = 2
    elseif(IsCylindrical)then
       r_ = 1; Phi_ = 2
    elseif(IsSpherical)then
       r_ = 1; Theta_ = 2; Phi_ = 3
    elseif(IsRLonLat)then
       r_ = 1; Phi_ = 2; Theta_ = 3; Lon_ =2; Lat_ = 3
    end if

    nRgen = -1
    if(allocated(LogRgen_I)) deallocate(LogRgen_I)
    if(IsGenRadius)then
       if(.not.present(RgenIn_I)) call CON_stop(NameSub// &
            ': RgenIn_I argument is missing for TypeGeometry=' //TypeGeometry)
 
       ! Store general radial coordinate table
       nRgen = size(RgenIn_I)
       allocate(LogRgen_I(nRgen))
       LogRgen_I = alog(RgenIn_I)
    end if

  end subroutine init_geometry

  !=========================================================================

  subroutine clean_geometry

    ! Release storage

    nRgen = -1
    if(allocated(LogRgen_I)) deallocate(LogRgen_I)

  end subroutine clean_geometry

  !=========================================================================
  
  subroutine xyz_to_coord(XyzIn_D, CoordOut_D)

    use ModCoordTransform, ONLY: atan2_check, xyz_to_sph
    use ModNumConst,       ONLY: cHalfPi

    real, intent(in) :: XyzIn_D(MaxDim)
    real, intent(out):: CoordOut_D(MaxDim)

    real:: x, y, r2

    character(len=*), parameter:: NameSub = 'BATL_geometry::xyz_to_coord'
    !----------------------------------------------------------------------

    if(IsCartesianGrid)then
       CoordOut_D = XyzIn_D
       RETURN
    elseif(IsRotatedCartesian)then
       CoordOut_D = matmul(XyzIn_D, GridRot_DD) 
    elseif(IsCylindrical)then
       x = XyzIn_D(1); y = XyzIn_D(2)
       CoordOut_D(1) = sqrt(x**2 + y**2)
       CoordOut_D(2) = atan2_check(y, x)
       CoordOut_D(3) = XyzIn_D(3)
    elseif(IsSpherical)then
       call xyz_to_sph(XyzIn_D, CoordOut_D)
    elseif(IsRLonLat)then
       ! Use xyz to r,theta,phi conversion
       call xyz_to_sph(XyzIn_D, &
            CoordOut_D(r_), CoordOut_D(Theta_), CoordOut_D(Phi_))
       ! Convert colatitude to latitude
       CoordOut_D(Lat_) = cHalfPi - CoordOut_D(Theta_)
    elseif(IsRoundCube)then
       r2 = sum(XyzIn_D**2)
       if(r2 > 0.0)then
          CoordOut_D = sqrt(r2/maxval(XyzIn_D**2)) * XyzIn_D
       else
          CoordOut_D = 0.0
       end if
    else
       call CON_stop(NameSub// &
            ' not yet implemented for TypeGeometry='//TypeGeometry)
    end if

    if(IsLogRadius)then
       CoordOut_D(1) = log(CoordOut_D(1))
    elseif(IsGenRadius)then
       call radius_to_gen(CoordOut_D(1))
    end if

  end subroutine xyz_to_coord

  !=========================================================================
  
  subroutine coord_to_xyz(CoordIn_D, XyzOut_D)

    use ModCoordTransform, ONLY: sph_to_xyz
    use ModNumConst,       ONLY: cHalfPi

    real, intent(in) :: CoordIn_D(MaxDim)
    real, intent(out):: XyzOut_D(MaxDim)

    real:: r, r2, Phi, Coord_D(MaxDim)

    character(len=*), parameter:: NameSub = 'BATL_geometry::coord_to_xyz'
    !----------------------------------------------------------------------

    if(IsCartesianGrid)then
       XyzOut_D = CoordIn_D
       RETURN
    elseif(IsRotatedCartesian)then
       XyzOut_D = matmul(GridRot_DD, CoordIn_D)
       RETURN
    endif

    Coord_D = CoordIn_D
    if(IsLogRadius)then
       Coord_D(1) = exp(Coord_D(1))
    elseif(IsGenRadius)then
       call gen_to_radius(Coord_D(1))
    end if

    if(IsCylindrical)then
       r = Coord_D(r_); 
       Phi = Coord_D(Phi_)
       XyzOut_D(1) = r*cos(Phi)
       XyzOut_D(2) = r*sin(Phi)
       XyzOut_D(3) = Coord_D(3)
    elseif(IsSpherical)then
       call sph_to_xyz(Coord_D, XyzOut_D)
    elseif(IsRLonLat)then
       call sph_to_xyz(Coord_D(r_), cHalfPi - Coord_D(Lat_), Coord_D(Phi_), &
            XyzOut_D)
    elseif(IsRoundCube)then
       r2 = sum(CoordIn_D**2)
       if(r2 > 0.0)then
          XyzOut_D = maxval(abs(CoordIn_D))/sqrt(r2) * CoordIn_D
       else
          XyzOut_D = 0.0
       end if
    else
       call CON_stop(NameSub// &
            ' not yet implemented for TypeGeometry='//TypeGeometry)
    end if

  end subroutine coord_to_xyz

  !============================================================================

  subroutine radius_to_gen(r)

    use ModInterpolate, ONLY: find_cell

    ! Convert true radial coordinate to general coordinate
    real, intent(inout):: r

    integer:: i
    real:: dCoord

    character(len=*), parameter:: NameSub='BATL_geometry::radius_to_gen'
    !-------------------------------------------------------------------------
    ! interpolate the logarithm of r
    call find_cell(0, nRgen-1, alog(r), &
         i, dCoord, LogRgen_I, DoExtrapolate=.true.)

    ! r is the real coordinate scaled to the 0-1 interval
    r = (i + dCoord)/(nRgen - 1)

  end subroutine radius_to_gen

  !============================================================================

  subroutine gen_to_radius(r)

    use ModInterpolate, ONLY: linear

    ! Convert generalized radial coordinate to true radial coordinate
    real, intent(inout):: r

    character(len=*), parameter:: NameSub='BATL_geometry::gen_to_radius'
    !-------------------------------------------------------------------------

    ! interpolate the LogRgen_I array for the general coordinate
    r = exp(linear(LogRgen_I, 0, nRgen-1, r*(nRgen-1), DoExtrapolate=.true.))

  end subroutine gen_to_radius

  !============================================================================

  subroutine test_geometry

    use BATL_mpi, ONLY: iProc

    logical:: IsPeriodicTest_D(MaxDim)
    real:: Xyz_D(MaxDim), Coord_D(MaxDim), Good_D(MaxDim)

    real:: r, GenR
    real:: Rgen_I(5) = (/ 1.0, 1.2, 5.0, 25.0, 100.0 /)

    logical:: DoTestMe
    character(len=*), parameter :: NameSub = 'test_geometry'
    !--------------------------------------------------------------------------
    DoTestMe = iProc == 0

    if(DoTestMe) write(*,*)'Starting ',NameSub
    if(DoTestMe) write(*,*)'Testing init_geometry for Cartesian'
    call init_geometry

    if(TypeGeometry /= 'cartesian') &
         write(*,*)'ERROR: init_geometry failed, ', &
         'TypeGeometry=', TypeGeometry, ' should be Cartesian by default'

    if(.not.IsCartesian .or. IsRotatedCartesian .or. &
         IsRzGeometry .or. IsSpherical .or. IsCylindrical)&
         write(*,*)'ERROR: init_geometry failed for Cartesian grid, ', &
         'IsCartesian, IsRzGeometry, IsSpherical, IsCylindrical=', &
         IsCartesian, IsRzGeometry, IsSpherical, IsCylindrical

    if(IsLogRadius .or. IsGenRadius) &
         write(*,*)'ERROR: init_geometry failed for Cartesian grid, ', &
         'IsLogRadius, IsGenRadius =', IsLogRadius, IsGenRadius

    if(any(IsPeriodic_D)) &
         write(*,*)'ERROR: init_geometry failed, ', &
         'IsPeriodic_D =', IsPeriodic_D, ' should be all false'

    if(DoTestMe) write(*,*)'Testing xyz_to_coord for Cartesian'
    Xyz_D = (/1., 2., 3./)
    call xyz_to_coord(Xyz_D, Coord_D)
    if(any(Coord_D /= Xyz_D)) &
         write(*,*)'ERROR: xyz_to_coord failed for Cartesian grid, ', &
         'Xyz_D =', Xyz_D, ' Coord_D =', Coord_D

    if(DoTestMe) write(*,*)'Testing coord_to_xyz for Cartesian'
    call coord_to_xyz(Coord_D, Xyz_D)
    if(any(Coord_D /= Xyz_D)) &
         write(*,*)'ERROR: coord_to_xyz failed for Cartesian grid, ', &
         'Xyz_D =', Xyz_D, ' Coord_D =', Coord_D

    if(nDim == 2)then
       if(DoTestMe) write(*,*)'Testing init_geometry for RZ geometry'
       IsPeriodicTest_D = (/.true., .false., .false./)

       call init_geometry('rz', IsPeriodicIn_D = IsPeriodicTest_D(1:nDim))

       if(TypeGeometry /= 'rz') &
            write(*,*)'ERROR: init_geometry failed, ', &
            'TypeGeometry=', TypeGeometry, ' should be rz'

       if(.not.IsRzGeometry .or. IsCartesian.or.IsSpherical.or.IsCylindrical)&
            write(*,*)'ERROR: init_geometry failed for RZ grid, ', &
            'IsCartesian, IsRzGeometry, IsSpherical, IsCylindrical=', &
            IsCartesian, IsRzGeometry, IsSpherical, IsCylindrical

       if(IsLogRadius .or. IsGenRadius) &
            write(*,*)'ERROR: init_geometry failed for RZ grid, ', &
            'IsLogRadius, IsGenRadius =', IsLogRadius, IsGenRadius

       if(any(IsPeriodic_D(1:nDim) .neqv. IsPeriodicTest_D(1:nDim))) &
            write(*,*)'ERROR: init_geometry failed, ', &
            'for TypeGeometry=', TypeGeometry,         &
            'IsPeriodic_D =', IsPeriodic_D(1:nDim),    &
            ' should be ', IsPeriodicTest_D(1:nDim)

       if(DoTestMe) write(*,*)'Testing xyz_to_coord for RZ geometry'
       Xyz_D = (/1., 2., 3./)
       call xyz_to_coord(Xyz_D, Coord_D)
       if(any(Coord_D /= Xyz_D)) &
            write(*,*)'ERROR: xyz_to_coord failed for RZ grid, ', &
            'Xyz_D =', Xyz_D, ' Coord_D =', Coord_D

       if(DoTestMe) write(*,*)'Testing coord_to_xyz for RZ geometry'
       call coord_to_xyz(Coord_D, Xyz_D)
       if(any(Coord_D /= Xyz_D)) &
            write(*,*)'ERROR: coord_to_xyz failed for RZ grid, ', &
            'Xyz_D =', Xyz_D, ' Coord_D =', Coord_D

    end if

    if(nDim == 1) RETURN

    if(DoTestMe) write(*,*)'Testing init_geometry for cylindrical_lnr'
    IsPeriodicTest_D = (/.false., .true., .true./)

    call init_geometry('cylindrical_lnr', &
         IsPeriodicIn_D = IsPeriodicTest_D(1:nDim))

    if(TypeGeometry /= 'cylindrical_lnr') &
         write(*,*)'ERROR: init_geometry failed, ', &
         'TypeGeometry=', TypeGeometry, ' should be cylindrical_lnr'

    if(.not.IsCylindrical .or. IsCartesian.or.IsRzGeometry.or.IsSpherical)&
         write(*,*)'ERROR: init_geometry failed for cylindrical_lnr, ', &
         'IsCartesian, IsRzGeometry, IsSpherical, IsCylindrical=', &
         IsCartesian, IsRzGeometry, IsSpherical, IsCylindrical

    if(.not.IsLogRadius .or. IsGenRadius) &
         write(*,*)'ERROR: init_geometry failed for cylindrical_lnr, ',&
         'IsLogRadius, IsGenRadius =', IsLogRadius, IsGenRadius

    if(any(IsPeriodic_D(1:nDim) .neqv. IsPeriodicTest_D(1:nDim))) &
         write(*,*)'ERROR: init_geometry failed, ', &
         'for TypeGeometry=', TypeGeometry,         &
         'IsPeriodic_D =', IsPeriodic_D(1:nDim),    &
         ' should be ', IsPeriodicTest_D(1:nDim)

    if(DoTestMe) write(*,*)'Testing xyz_to_coord for cylindrical_lnr'
    Xyz_D = (/1., 2., 3./)
    Good_D = (/log(sqrt(5.)), atan2(2.,1.), 3./)
    call xyz_to_coord(Xyz_D, Coord_D)
    if(any(abs(Coord_D - Good_D) > 1e-6)) &
         write(*,*)'ERROR: xyz_to_coord failed for cylindrical_lnr, ', &
         'Xyz_D =', Xyz_D, ' Coord_D =', Coord_D,' should be ', Good_D

    if(DoTestMe) write(*,*)'Testing init_geometry for roundcube'
    IsPeriodicTest_D = (/.false., .false., .false./)

    call init_geometry('roundcube', &
         IsPeriodicIn_D = IsPeriodicTest_D(1:nDim))

    if(TypeGeometry /= 'roundcube') &
         write(*,*)'ERROR: init_geometry failed, ', &
         'TypeGeometry=', TypeGeometry, ' should be roundcube'

    if(.not.IsRoundCube .or. IsCartesian.or.IsRzGeometry &
         .or.IsCylindrical.or.IsSpherical)&
         write(*,*)'ERROR: init_geometry failed for roundcube, ', &
         'IsRoundCube,IsCartesian,IsRzGeometry,IsSpherical,IsCylindrical=', &
         IsRoundCube, IsCartesian, IsRzGeometry, IsSpherical, IsCylindrical

    if(IsLogRadius .or. IsGenRadius) &
         write(*,*)'ERROR: init_geometry failed for roundcube, ',&
         'IsLogRadius, IsGenRadius =', IsLogRadius, IsGenRadius

    if(any(IsPeriodic_D(1:nDim) .neqv. IsPeriodicTest_D(1:nDim))) &
         write(*,*)'ERROR: init_geometry failed, ', &
         'for TypeGeometry=', TypeGeometry,         &
         'IsPeriodic_D =', IsPeriodic_D(1:nDim),    &
         ' should be ', IsPeriodicTest_D(1:nDim)

    if(DoTestMe) write(*,*)'Testing xyz_to_coord for cubsedsphere'
    Xyz_D  = (/3., 4., nDim-2.0/)
    Good_D = sqrt(sum(Xyz_D**2))/4.0 * Xyz_D
    call xyz_to_coord(Xyz_D, Coord_D)
    if(any(abs(Coord_D - Good_D) > 1e-6)) &
         write(*,*)'ERROR: xyz_to_coord failed for cubsedsphere, ', &
         'Xyz_D =', Xyz_D, ' Coord_D =', Coord_D,' should be ', Good_D

    if(DoTestMe) write(*,*)'Testing coord_to_xyz for cubsedsphere'
    Good_D = Xyz_D
    call coord_to_xyz(Coord_D, Xyz_D)
    if(any(abs(Xyz_D - Good_D) > 1e-6)) &
         write(*,*)'ERROR: coord_to_xyz failed for cubsedsphere, ', &
         'Coord_D =', Coord_D, ' Xyz_D =', Xyz_D,' should be ', Good_D

    if(nDim < 3) RETURN

    if(DoTestMe) write(*,*)'Testing init_geometry for spherical_genr'
    IsPeriodicTest_D = (/.false., .true., .false./)

    call init_geometry('spherical_genr', IsPeriodicIn_D = IsPeriodicTest_D, &
         RgenIn_I = Rgen_I )

    if(TypeGeometry /= 'spherical_genr') &
         write(*,*)'ERROR: init_geometry failed, ', &
         'TypeGeometry=', TypeGeometry, ' should be spherical_genr'

    if(.not.IsSpherical .or. IsCartesian.or.IsRzGeometry.or.IsCylindrical)&
         write(*,*)'ERROR: init_geometry failed for spherical_genr, ', &
         'IsCartesian, IsRzGeometry, IsCylindrical, IsSpherical=', &
         IsCartesian, IsRzGeometry, IsCylindrical, IsSpherical

    if(.not.IsGenRadius .or. IsLogRadius) &
         write(*,*)'ERROR: init_geometry failed for spherical_genr, ',&
         'IsLogRadius, IsGenRadius =', IsLogRadius, IsGenRadius

    if(any(IsPeriodic_D(1:nDim) .neqv. IsPeriodicTest_D(1:nDim))) &
         write(*,*)'ERROR: init_geometry failed, ', &
         'for TypeGeometry=', TypeGeometry,         &
         'IsPeriodic_D =', IsPeriodic_D(1:nDim),    &
         ' should be ', IsPeriodicTest_D(1:nDim)

    if(nRgen /= size(Rgen_I)) &
         write(*,*)'ERROR: init_geometry failed, ', &
         'for TypeGeometry=', TypeGeometry,         &
         'nRgen=', nRgen,' should be ',size(Rgen_I)
    
    if(.not.allocated(LogRgen_I)) &
         write(*,*)'ERROR: init_geometry failed, ', &
         'for TypeGeometry=', TypeGeometry,         &
         'LogRgen_I is not allocated'

    if(any(abs(exp(LogRgen_I) - Rgen_I) > 1e-6)) &
         write(*,*)'ERROR: init_geometry failed, ', &
         'for TypeGeometry=', TypeGeometry,         &
         'exp(LogRgen_I) =', exp(LogRgen_I),' should be ',Rgen_I

    if(DoTestMe) write(*,*)'Testing radius_to_gen and gen_to_radius'
    r = sqrt(Rgen_I(2)*Rgen_I(3))
    GenR = r
    call radius_to_gen(GenR)
    if(abs(GenR - 1.5/4) > 1e-6) &
         write(*,*)'ERROR: radius_to_gen failed for spherical_genr, ', &
         'r=', r,' GenR =', GenR, ' should be ', 1.5/4

    ! Test conversion back
    call gen_to_radius(GenR)
    if(abs(GenR - r) > 1e-6) &
         write(*,*)'ERROR: gen_to_radius failed for spherical_genr, ', &
         'Orig r=', r,' new r =', GenR

    r = 1.0/Rgen_I(2)**2
    GenR = r
    call radius_to_gen(GenR)
    if(abs(GenR + 2.0/4) > 1e-6) &
         write(*,*)'ERROR: radius_to_gen failed for spherical_genr, ', &
         'r=', r,' GenR =', GenR, ' should be ', -2.0/4

    ! Test conversion back
    call gen_to_radius(GenR)
    if(abs(GenR - r) > 1e-6) &
         write(*,*)'ERROR: gen_to_radius failed for spherical_genr, ', &
         'Orig r=', r,' new r =', GenR

    r = 1600.0
    GenR = r
    call radius_to_gen(GenR)
    if(abs(GenR - (1+2./4)) > 1e-6) &
         write(*,*)'ERROR: radius_to_gen failed for spherical_genr, ', &
         'r=', r,' GenR =', GenR, ' should be ', 1 + 2./4.

    ! Test conversion back
    call gen_to_radius(GenR)
    if(abs(GenR - r) > 1e-6) &
         write(*,*)'ERROR: gen_to_radius failed for spherical_genr, ', &
         'Orig r=', r,' new r =', GenR

    if(DoTestMe) write(*,*)'Testing xyz_to_coord for spherical_genr'
    Xyz_D = (/9., 12., 20./)
    Good_D = (/0.75, atan2(15.,20.),  atan2(12., 9.)/)
    call xyz_to_coord(Xyz_D, Coord_D)
    if(any(abs(Coord_D - Good_D) > 1e-6)) &
         write(*,*)'ERROR: xyz_to_coord failed for spherical_genr, ', &
         'Xyz_D =', Xyz_D, ' Coord_D =', Coord_D,' should be ', Good_D

    if(DoTestMe) write(*,*)'Testing coord_to_xyz for spherical_genr'
    Good_D = Xyz_D
    call coord_to_xyz(Coord_D, Xyz_D)
    if(any(abs(Xyz_D - Good_D) > 1e-6)) &
         write(*,*)'ERROR: coord_to_xyz failed for spherical, ', &
         'Coord_D =', Coord_D, ' Xyz_D =', Xyz_D,' should be ', Good_D

  end subroutine test_geometry

end module BATL_geometry
