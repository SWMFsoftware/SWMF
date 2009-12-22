module BATL_geometry

  use BATL_size, ONLY: MaxDim, nDim

  implicit none

  SAVE

  private ! except

  public:: init_geometry
  public:: xyz_to_coord

  character(len=20), public:: TypeGeometry = 'cartesian'

  ! Cartesian, cylindrical or spherical coordinates
  logical, public:: IsCartesian   = .true.
  logical, public:: IsRzGeometry  = .false.
  logical, public:: IsSpherical   = .false.
  logical, public:: IsCylindrical = .false.

  ! Periodicity of the domain per dimension
  logical, public:: IsPeriodic_D(MaxDim) = .false.

  ! Index names for coordinates
  integer, parameter, public:: x_=1, y_=min(nDim,2), z_=min(nDim,3), r_=y_ 

contains

  subroutine init_geometry(TypeGeometryIn, IsPeriodicIn_D)

    character(len=*), optional, intent(in):: TypeGeometryIn
    logical,          optional, intent(in):: IsPeriodicIn_D(nDim)
    !-----------------------------------------------------------------------

    TypeGeometry = 'cartesian'
    if(present(TypeGeometryIn)) TypeGeometry = TypeGeometryIn

    IsCartesian   = TypeGeometry(1:9)  == 'cartesian'
    IsRzGeometry  = TypeGeometry(1:2)  == 'rz'
    IsSpherical   = TypeGeometry(1:9)  == 'spherical'
    IsCylindrical = TypeGeometry(1:11) == 'cylindrical'

    ! Set default for periodicity
    IsPeriodic_D = .false.
    if(IsCylindrical) IsPeriodic_D = (/ .false., .true., .false./)
    if(IsSpherical)   IsPeriodic_D = (/ .false., .false., .true./)

    ! Set periodicity if argument is provided
    if(present(IsPeriodicIn_D)) IsPeriodic_D(1:nDim) = IsPeriodicIn_D

  end subroutine init_geometry

  !=========================================================================
  
  subroutine xyz_to_coord(XyzIn_D, CoordOut_D)

    real, intent(in) :: XyzIn_D(MaxDim)
    real, intent(out):: CoordOut_D(MaxDim)

    character(len=*), parameter:: NameSub = 'BATL_geometry::xyz_to_coord'
    !----------------------------------------------------------------------

    if(IsCartesian .or. IsRzGeometry)then
       CoordOut_D = XyzIn_D
    else
       call CON_stop(NameSub// &
            ' not yet implemented for TypeGeometry='//TypeGeometry)
    end if

  end subroutine xyz_to_coord

end module BATL_geometry
