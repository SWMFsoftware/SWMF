module BATL_geometry

  use BATL_size, ONLY: MaxDim, nDim

  implicit none

  SAVE

  private ! except

  public:: init_geometry  ! initialize the module
  public:: xyz_to_coord   ! convert XYZ coordinates to generalized coordinates
  public:: coord_to_xyz   ! convert generalized coordinates to XYZ coordinates
  public:: radius_to_gen  ! convert radial coordinate to generalized coordinate
  public:: gen_to_radius  ! convert generalized coordinate to radial coordinate
  public:: test_geometry  ! unit test

  character(len=20), public:: TypeGeometry = 'cartesian'

  ! Cartesian, cylindrical or spherical coordinates
  logical, public:: IsCartesian       = .true.
  logical, public:: IsRzGeometry      = .false.
  logical, public:: IsSpherical       = .false.
  logical, public:: IsSphericalAxis   = .false. ! theta=0 and theta=pi boundary
  logical, public:: IsCylindrical     = .false.
  logical, public:: IsCylindricalAxis = .false. ! r=0 boundary
  logical, public:: IsCubedSphere     = .false.
  logical, public:: IsLogRadius       = .false. ! logarithmic radial coordinate
  logical, public:: IsGenRadius       = .false. ! stretched radial coordinate

  ! Periodicity of the domain per dimension
  logical, public:: IsPeriodic_D(MaxDim) = .false.

  ! Index names for coordinates
  integer, parameter, public:: x_=1, y_=min(nDim,2), z_=min(nDim,3)

  ! The following index names will be set in init_geometry
  integer, public:: r_=-1, Phi_=-1, Theta_=-1, Lon_=-1, Lat_=-1

contains

  !=========================================================================

  subroutine init_geometry(TypeGeometryIn, IsPeriodicIn_D)

    character(len=*), optional, intent(in):: TypeGeometryIn
    logical,          optional, intent(in):: IsPeriodicIn_D(nDim)
    !-----------------------------------------------------------------------

    TypeGeometry = 'cartesian'
    if(present(TypeGeometryIn)) TypeGeometry = TypeGeometryIn

    IsPeriodic_D = .false.
    if(present(IsPeriodicIn_D)) IsPeriodic_D(1:nDim) = IsPeriodicIn_D

    ! Logicals are useful for efficient code
    IsCartesian   = TypeGeometry(1:9)  == 'cartesian'
    IsRzGeometry  = TypeGeometry(1:2)  == 'rz'
    IsSpherical   = TypeGeometry(1:3)  == 'sph'
    IsCylindrical = TypeGeometry(1:3)  == 'cyl'
    IsCubedSphere = TypeGeometry(1:3)  == 'cub'

    IsLogRadius   = index(TypeGeometry,'lnr')  > 0
    IsGenRadius   = index(TypeGeometry,'genr') > 0

    r_ = -1; Phi_ = -1; Theta_ = -1
    if(IsRzGeometry)then
       r_ = 2
    elseif(IsCylindrical)then
       r_ = 1; Phi_ = 2
    elseif(IsSpherical)then
       r_ = 1; Theta_ = 2; Phi_ = 3
    end if

  end subroutine init_geometry

  !=========================================================================
  
  subroutine xyz_to_coord(XyzIn_D, CoordOut_D)

    use ModCoordTransform, ONLY: atan2_check, xyz_to_sph

    real, intent(in) :: XyzIn_D(MaxDim)
    real, intent(out):: CoordOut_D(MaxDim)

    real:: x, y, r2

    character(len=*), parameter:: NameSub = 'BATL_geometry::xyz_to_coord'
    !----------------------------------------------------------------------

    if(IsCartesian .or. IsRzGeometry)then
       CoordOut_D = XyzIn_D
       RETURN
    elseif(IsCylindrical)then
       x = XyzIn_D(1); y = XyzIn_D(2)
       CoordOut_D(1) = sqrt(x**2 + y**2)
       CoordOut_D(2) = atan2_check(y, x)
       CoordOut_D(3) = XyzIn_D(3)
    elseif(IsSpherical)then
       call xyz_to_sph(XyzIn_D, CoordOut_D)
    elseif(IsCubedSphere)then
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

    real, intent(in) :: CoordIn_D(MaxDim)
    real, intent(out):: XyzOut_D(MaxDim)

    real:: r, r2, Phi, Coord_D(MaxDim)

    character(len=*), parameter:: NameSub = 'BATL_geometry::coord_to_xyz'
    !----------------------------------------------------------------------

    if(IsCartesian .or. IsRzGeometry)then
       XyzOut_D = CoordIn_D
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
    elseif(IsCubedSphere)then
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

    ! Convert true radial coordinate to general coordinate
    real, intent(inout):: r

    character(len=*), parameter:: NameSub='BATL_geometry::radius_to_gen'

    call CON_stop(NameSub//' not yet implemented')
  end subroutine radius_to_gen

  !============================================================================

  subroutine gen_to_radius(r)

    ! Convert generalized radial coordinate to true radial coordinate
    real, intent(inout):: r

    character(len=*), parameter:: NameSub='BATL_geometry::gen_to_radius'

    call CON_stop(NameSub//' not yet implemented')

  end subroutine gen_to_radius

  !============================================================================

  subroutine test_geometry

    use BATL_mpi, ONLY: iProc

    logical:: IsPeriodicTest_D(MaxDim)
    real:: Xyz_D(MaxDim), Coord_D(MaxDim), Good_D(MaxDim)

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

    if(.not.IsCartesian .or. IsRzGeometry .or. IsSpherical .or. IsCylindrical)&
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

    if(DoTestMe) write(*,*)'Testing init_geometry for cubedsphere'
    IsPeriodicTest_D = (/.false., .false., .false./)

    call init_geometry('cubedsphere', &
         IsPeriodicIn_D = IsPeriodicTest_D(1:nDim))

    if(TypeGeometry /= 'cubedsphere') &
         write(*,*)'ERROR: init_geometry failed, ', &
         'TypeGeometry=', TypeGeometry, ' should be cubedsphere'

    if(.not.IsCubedSphere .or. IsCartesian.or.IsRzGeometry &
         .or.IsCylindrical.or.IsSpherical)&
         write(*,*)'ERROR: init_geometry failed for cubedsphere, ', &
         'IsCubedSphere,IsCartesian,IsRzGeometry,IsSpherical,IsCylindrical=', &
         IsCubedSphere, IsCartesian, IsRzGeometry, IsSpherical, IsCylindrical

    if(IsLogRadius .or. IsGenRadius) &
         write(*,*)'ERROR: init_geometry failed for cubedsphere, ',&
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

    call init_geometry('spherical_genr', IsPeriodicIn_D = IsPeriodicTest_D)

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

    ! genr is not yet implemented, so switch to plain spherical
    call init_geometry('spherical')

    if(DoTestMe) write(*,*)'Testing xyz_to_coord for spherical'
    Xyz_D = (/9., 12., 20./)
    Good_D = (/25., atan2(15.,20.),  atan2(12., 9.)/)
    call xyz_to_coord(Xyz_D, Coord_D)
    if(any(abs(Coord_D - Good_D) > 1e-6)) &
         write(*,*)'ERROR: xyz_to_coord failed for spherical, ', &
         'Xyz_D =', Xyz_D, ' Coord_D =', Coord_D,' should be ', Good_D

    if(DoTestMe) write(*,*)'Testing coord_to_xyz for spherical'
    Good_D = Xyz_D
    call coord_to_xyz(Coord_D, Xyz_D)
    if(any(abs(Xyz_D - Good_D) > 1e-6)) &
         write(*,*)'ERROR: coord_to_xyz failed for spherical, ', &
         'Coord_D =', Coord_D, ' Xyz_D =', Xyz_D,' should be ', Good_D


end subroutine test_geometry

end module BATL_geometry
