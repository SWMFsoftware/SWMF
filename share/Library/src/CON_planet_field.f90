!^CFG COPYRIGHT UM
!BOP
!MODULE: CON_planet_field - provide value and mapping of magnetic field
!INTERFACE:
module CON_planet_field

  !DESCRIPTION:
  ! This class provides the magnetic field of the planet
  ! for an arbitrary spatial position at an arbitrary time. 
  ! It also provides the mapping from an arbitrary point to a given 
  ! radial distance.
  ! The position as well as the magnetic field can be represented as 3 
  ! scalars or 1 three-element array.
  ! The coordinate system and the normalization of the coordinates
  ! and the magnetic field can be given with string input arguments.

  !USES:
  use CON_planet
  use CON_axes

  implicit none

  private ! except

  !PUBLIC MEMBER FUNCTIONS
  public :: get_planet_field  ! Get planet field at some time and place
  public :: map_planet_field  ! Map planet field from a point to a radius 
  public :: test_planet_field ! Test the methods in this module

  !REVISION HISTORY:
  ! 11Aug03 - Gabor Toth <gtoth@umich.edu> - initial prototype/prolog/code
  !EOP -------------------------------------------------------------------

  interface get_planet_field
     module procedure &
          get_planet_field11, &
          get_planet_field13, &
          get_planet_field31, &
          get_planet_field33
  end interface

  interface map_planet_field
     module procedure &
          map_planet_field11, &
          map_planet_field33
  end interface

  integer, parameter :: x_=1, y_=2, z_=3

  character(len=*), parameter :: NameMod = 'CON_planet_field'


contains

  !IROUTINE: get_planet_field - get planet field at some time and position
  !INTERFACE:
  subroutine get_planet_field11(TimeSim, XyzIn_D, TypeCoord, b_D)

    !INPUT ARGUMENTS:
    real,              intent(in) :: TimeSim      ! simulation time
    real,              intent(in) :: XyzIn_D(3)   ! spatial position
    character(len=*),  intent(in) :: TypeCoord    ! type of coordinates

    !OUTPUT ARGUMENTS:
    real,              intent(out):: b_D(3)       ! magnetic field

    !DESCRIPTION:
    ! This is the fundamental subroutine that provides the magnetic
    ! field at a given position at a given simulation time. 
    ! If called repeatedly, the subroutine remembers the last simulation time
    ! argument, so it does not recalculate the position of the magnetic axis.
    ! The position may be normalized with the radius of the planet.
    ! The coordinate system and normalization information 
    ! for the position is given by the string TypeCoord.
    ! The first 3 characters should contain the coordinate system.
    ! This may be followed (after some spaces) by the characters "NORM" 
    ! in all capitals. For example "MAG", "GSM NORM", "GSE NORMALIZED" etc.

    !EOP
    character(len=*), parameter :: NameSub=NameMod//'::get_planet_field'

    real :: Xyz_D(3)     ! Normalized (and rotated) position
    real :: Dipole_D(3)  ! Dipole moment
    real :: r, r2, rInv, r2Inv, r3Inv, Term1
    character (len=3) :: NameCoordSystem

!!! logical :: DoTest, DoTestMe
    !------------------------------------------------------------------------

!!! call CON_set_do_test(NameSub,DoTest,DoTestMe)

    if(TypeBField == 'NONE')then
       b_D = 0
       RETURN
    end if

    if(len(TypeCoord)<3)call CON_stop(NameSub//&
         ' SWMF_ERROR: coordinate type should be at least 3 characters,'// &
         ' TypeCoord='//TypeCoord)

    ! Normalize position if necessary
    if(index(TypeCoord,"NORM")>0)then
       Xyz_D = XyzIn_D
    else
       Xyz_D = XyzIn_D / RadiusPlanet
    end if

    ! radial distance squared
    r2 = sum(Xyz_D**2)
    if(r2 < 1E-12) then
       ! return zero field if very small
       b_D = 0
       RETURN
    end if

    ! Update axes
    call set_axes(TimeSim)

    ! The coord system name is stored in the first 3 characters
    NameCoordSystem = TypeCoord(1:3)

    ! Calculate magnetic field
    select case(TypeBField)
    case('DIPOLE')
       ! Various powers of radial distance
       r2Inv = 1/r2
       r     = sqrt(r2)
       rInv  = 1/r
       r3Inv = rInv*r2Inv

       select case(NameCoordSystem)
       case('GSM')
          ! Dipole_D has X and Z components only
          Dipole_D = DipoleStrength*MagAxisGsm_D
          Term1 = sum(Dipole_D(x_:z_:2)*Xyz_D(x_:z_:2))*3*r2Inv
          b_D = (Term1*Xyz_D - Dipole_D)*r3Inv
       case('GSE')
          ! Dipole_D has X,Y and Z components in general
          Dipole_D = DipoleStrength*MagAxis_D
          Term1 = sum(Dipole_D*Xyz_D)*3*r2Inv
          b_D = (Term1*Xyz_D - Dipole_D)*r3Inv
       case('GEO')
          ! Dipole_D has X,Y and Z components in general
          Dipole_D = DipoleStrength*MagAxisGeo_D
          Term1 = sum(Dipole_D*Xyz_D)*3*r2Inv
          b_D = (Term1*Xyz_D - Dipole_D)*r3Inv
       case('MAG','SMG')
          ! Dipole is aligned with the Z axis
          Term1      = DipoleStrength*Xyz_D(3)*3*r2Inv
          b_D(x_:y_) = Term1*Xyz_D(x_:y_)*r3Inv
          b_D(z_)    = (Term1*Xyz_D(z_)-DipoleStrength)*r3Inv
       case default
          call CON_stop(NameSub// &
               ' SWMF_ERROR: unimplemented NameCoordSystem='//NameCoordSystem)
       end select
       !case('QUADRUPOLE','OCTUPOLE')
       !   ! Transform to MAG system
       !   Xyz_D = coord_transform(TimeSim,Xyz_D,NameCoordSystem,'MAG')
       !
       !   ! Dipole is aligned with Z
       !   Term1      = DipoleStrength*Xyz_D(3)*3*r2Inv
       !   b_D(x_:y_) = Term1*Xyz_D(1:2)*r3Inv
       !   b_D(z_)    = (Term1*Xyz_D(3)-DipoleStrength)*r3Inv
       !
       !   ! Add quadrupole terms
       !
       !   if(TypeBField == 'OCTUPOLE')then
       !      ! Add octupole terms
       !   end if
       !
       !   ! Transform the magnetic field back to the input coordinate system
       !   b_D = coord_transform(TimeSim,b_D,'MAG',NameCoordSystem)
    case default
       call CON_stop(NameSub//' SWMF_ERROR: unimplemented TypeBField='//&
            TypeBField)
    end select

  end subroutine get_planet_field11

  !============================================================================

  subroutine get_planet_field13(TimeSim, XyzIn_D, TypeCoord, Bx, By, Bz)

    real,              intent(in) :: TimeSim      ! simulation time
    real,              intent(in) :: XyzIn_D(3)   ! spatial position
    character(len=*),  intent(in) :: TypeCoord    ! type of coordinates
    real,              intent(out):: Bx, By, Bz   ! magnetic field

    real :: b_D(3)

    call get_planet_field(TimeSim, XyzIn_D, TypeCoord, b_D)

    Bx = b_D(x_)
    By = b_D(y_)
    Bz = b_D(z_)

  end subroutine get_planet_field13

  !============================================================================

  subroutine get_planet_field31(TimeSim, x, y, z, TypeCoord, b_D)

    real,              intent(in) :: TimeSim      ! simulation time
    real,              intent(in) :: x, y, z      ! spatial position
    character(len=*),  intent(in) :: TypeCoord    ! type of coordinates
    real,              intent(out):: b_D(3)       ! magnetic field

    call get_planet_field(TimeSim, (/x, y, z/), TypeCoord, b_D)

  end subroutine get_planet_field31

  !============================================================================

  subroutine get_planet_field33(TimeSim, x, y, z, TypeCoord, Bx, By, Bz)

    real,              intent(in) :: TimeSim      ! simulation time
    real,              intent(in) :: x, y, z      ! spatial position
    character(len=*),  intent(in) :: TypeCoord    ! type of coordinates
    real,              intent(out):: Bx, By, Bz   ! magnetic field

    real :: b_D(3)

    call get_planet_field(TimeSim, (/x, y, z/), TypeCoord, b_D)

    Bx = b_D(x_)
    By = b_D(y_)
    Bz = b_D(z_)

  end subroutine get_planet_field33

  !BOP ========================================================================
  !IROUTINE: map_planet_field - map planet field from a position to some radius
  !INTERFACE:
  subroutine map_planet_field11(TimeSim, XyzIn_D, TypeCoord, &
       rMapIn, XyzMap_D, iHemisphere)

    !INPUT ARGUMENTS:
    real,              intent(in) :: TimeSim      ! simulation time
    real,              intent(in) :: XyzIn_D(3)   ! spatial position
    character(len=*),  intent(in) :: TypeCoord    ! type of coordinates
    real,              intent(in) :: rMapIn       ! radial distance to map to

    !OUTPUT ARGUMENTS:
    real,              intent(out):: XyzMap_D(3)  ! mapped position
    integer,           intent(out):: iHemisphere  ! which hemisphere

    !DESCRIPTION:
    ! Map the planet field from the input position to the mapping radius.
    ! The coordinate system and normalization used for the input position
    ! will also be used for the output position. The routine also returns
    ! which hemisphere the point is mapping to: +1 for north, -1 for south.
    ! If the point does not map to the defined radius at all, 0 is returned,
    ! and the output position is set to a radial projection of the input
    ! position to the magnetic equator.

    !PARAMETERS:

    ! If the  normalized input or mapping radius is less than this value
    ! an error occurs
    real, parameter :: rNormLimit = 0.9

    ! If difference between the normalized input and mapping radii 
    ! is less than DrNormLimit a trivial mapping is done
    real, parameter :: DrNormLimit = 0.0001

    !EOP

    character(len=*), parameter :: NameSub=NameMod//'::map_planet_field'
    real             :: Xyz_D(3)        ! Normalized and rotated position
    character(len=3) :: NameCoordSystem ! Input/Output coordinate system

    ! Temporary variables for the analytic mapping
    real :: rMap, rMap2, rMap3, r, r3, XyRatio, XyMap2

!!! logical :: DoTest, DoTestMe
    !------------------------------------------------------------------------

!!! call CON_set_do_test(NameSub,DoTest,DoTestMe)

    if(TypeBField == 'NONE') call CON_stop(NameSub// &
         ' SWMF_ERROR: the planet has no magnetic field')

    if(len(TypeCoord)<3)call CON_stop(NameSub//&
         ' SWMF_ERROR: coordinate type should be at least 3 characters,'// &
         ' TypeCoord='//TypeCoord)

    ! Normalize position and mapping radius if necessary
    if(index(TypeCoord,"NORM")>0)then
       Xyz_D = XyzIn_D
       rMap  = rMapIn
    else
       Xyz_D = XyzIn_D / RadiusPlanet
       rMap  = rMapIn  / RadiusPlanet
    end if

    ! Check if the mapping radius is outside of the planet
    if ( rMap < rNormLimit ) then
       write(*,*)NameSub,' mapping radius and coord type =',rMapIn,TypeCoord
       write(*,*)NameSub,' normalized mapping radius rMap=',rMap
       call CON_stop(NameSub// &
            ' SWMF_ERROR mapping radius is less than planet radius')
    end if

    ! Convert input position into MAG or SMG system
    NameCoordSystem = TypeCoord(1:3)
    select case(NameCoordSystem)
    case('MAG','SMG')
       ! There is nothing to do
    case ('GEO')
       ! Convert into MAG
       Xyz_D = matmul(MagGeo_DD,Xyz_D)
    case('GSM')
       ! Convert into SMG
       call set_axes(TimeSim)
       Xyz_D = matmul(SmgGsm_DD,Xyz_D)
    case('GSE')
       ! Convert into SMG
       call set_axes(TimeSim)
       Xyz_D = matmul(SmgGse_DD,Xyz_D)
    case default
       call CON_stop(NameSub//' SWMF_ERROR unimplemented NameCoordSystem='//&
            NameCoordSystem)
    end select

    ! In MAG/SMG coordinates the hemisphere depends on the sign of Z
    iHemisphere = sign(1.0,Xyz_D(3))

    ! Normalized radial distance
    r = sqrt(sum(Xyz_D**2))

    ! Check if the point is outside the planet
    if ( r < rNormLimit ) then
       write(*,*)NameSub,' input position and coord type=',XyzIn_D,TypeCoord
       write(*,*)NameSub,' normalized radius r=',r
       call CON_stop(NameSub// &
            ' SWMF_ERROR radial distance is less than planet radius')
    end if

    ! Check if the mapping radius differs from the radius of input position
    if( abs(r-rMap) < DrNormLimit ) then
       ! Trivial mapping
       XyzMap_D = XyzIn_D
       ! The hemisphere has been established already
       RETURN
    end if

    ! Find the mapped position
    select case(TypeBField)
    case('DIPOLE')
       ! Solution of the vector potential equation
       ! The vector potential is proportional to (x^2+y^2)/r^3
       ! so sqrt(xMap^2+yMap^2)/sqrt(x^2+y^2) = sqrt(rMap^3/r^3)

       ! Calculate powers of the radii
       rMap2 = rMap**2
       rMap3 = rMap2*rMap
       r3    = r**3

       ! This is the ratio of input and mapped X and Y components
       XyRatio = sqrt(rMap3/r3)

       ! Calculate the X and Y components of the mapped position
       XyzMap_D(1:2) = XyRatio*Xyz_D(1:2)

       ! The squared distance of the mapped position from the magnetic axis
       XyMap2 = XyzMap_D(1)**2 + XyzMap_D(2)**2

       ! Check if there is a mapping at all
       if(rMap2 < XyMap2)then
          ! The point does not map to the given radius
          iHemisphere = 0

          ! Put mapped point to the magnetic equator
          XyzMap_D(1:2) = (rMap/sqrt(sum(Xyz_D(1:2)**2)))*Xyz_D(1:2)
          XyzMap_D(3) = 0
       else
          ! Calculate the Z component of the mapped position
          ! Select the same hemisphere as for the input position
          XyzMap_D(3) = iHemisphere*sqrt(rMap2 - XyMap2)
       end if
    case default
       call CON_stop(NameSub//' unimplemented TypeBField='//TypeBField)
    end select

    ! Convert position back to the input coordinate system
    select case(NameCoordSystem)
    case('GEO')
       ! Convert back from MAG
       XyzMap_D = matmul(XyzMap_D,MagGeo_DD)
    case('GSM')
       ! Convert back from SMG
       XyzMap_D = matmul(XyzMap_D,SmgGsm_DD)
    case('GSE')
       ! Convert back from SMG
       XyzMap_D = matmul(XyzMap_D,SmgGse_DD)
    end select

    ! Undo the normalization
    if( index(TypeCoord,"NORM")<=0 ) &
         Xyzmap_D = Xyzmap_D*RadiusPlanet

  end subroutine map_planet_field11

  !BOP ========================================================================
  !IROUTINE: map_planet_field33 - map planet field from a position to some radius
  !INTERFACE:
  subroutine map_planet_field33(TimeSim, xIn, yIn, zIn, TypeCoord, &
       rMap, xMap, yMap, zMap, iHemisphere)

    !INPUT ARGUMENTS:
    real,             intent(in) :: TimeSim       ! simulation time
    real,             intent(in) :: xIn, yIn, zIn ! spatial position
    character(len=*), intent(in) :: TypeCoord     ! type of coordinates
    real,             intent(in) :: rMap          ! radial distance to map to

    !OUTPUT ARGUMENTS:
    real,              intent(out):: xMap, yMap, zMap  ! mapped position
    integer,           intent(out):: iHemisphere       ! mapped hemisphere

    !DESCRIPTION:
    ! Interface to the map\_planet\_field11 routine with 3 scalars for both
    ! input and output positions

    !LOCAL VARIABLES:
    real :: XyzIn_D(3), XyzMap_D(3)
    !EOP
    !-------------------------------------------------------------------------
    !BOC
    XyzIn_D(1)=xIn; XyzIn_D(2)=yIn; XyzIn_D(3)=zIn

    call map_planet_field(TimeSim, XyzIn_D, TypeCoord, rMap, &
         XyzMap_D, iHemisphere)

    xMap=XyzMap_D(1); yMap=XyzMap_D(2); zMap=XyzMap_D(3)
    !EOC
  end subroutine map_planet_field33

  !BOP =======================================================================
  !IROUTINE: test_planet_field - test methods in CON_planet_field
  !INTERFACE:
  subroutine test_planet_field
    !DESCRIPTION:
    ! Test the methods in this class.
    !EOP

    real :: TimeSim
    real :: xSmg_D(3), xGsm_D(3), xGse_D(3), bSmg_D(3), bGsm_D(3), bGse_D(3)
    real :: x_D(3), rMap, xMap_D(3)
    integer :: iHemisphere
    !------------------------------------------------------------------------

    call init_axes(TimeEquinox % Time)

    write(*,*)
    write(*,*)'TEST GET_PLANET_FIELD'
    write(*,*)

    xSmg_D = (/1.0, 1.0, 0.1/)
    write(*,*)'Location xSmg_D = ',xSmg_D
    call get_planet_field(0.0,xSmg_D,'SMG NORM',bSmg_D)
    write(*,*)'Field    bSmg_D = ',bSmg_D
    write(*,*)
    call get_planet_field(0.0,xSmg_D*RadiusPlanet,'SMG',bSmg_D)
    write(*,*)'Field    bSmg_D = ',bSmg_D
    write(*,*)
    xGsm_D = matmul(xSmg_D,SmgGsm_DD)
    write(*,*)'Location xGsm_D =',xGsm_D
    call get_planet_field(0.0,xGsm_D,'GSM NORM',bGsm_D)
    write(*,*)'Field    bGsm_D = ',bGsm_D
    write(*,*)'Rotated  bGsm_D = ',matmul(SmgGsm_DD,bGsm_D)
    write(*,*)
    xGse_D = matmul(xGsm_D,GsmGse_DD)
    write(*,*)'Location xGse_D =',xGse_D
    call get_planet_field(0.0,xGse_D,'GSE NORM',bGse_D)
    write(*,*)'Field    bGse_D = ',bGse_D
    write(*,*)'Rotated  bGse_D = ',matmul(SmgGsm_DD,matmul(GsmGse_DD,bGse_D))

    TimeSim = 6*3600
    write(*,*)'Test at time=', TimeSim
    call get_planet_field(TimeSim,xSmg_D,'SMG NORM',bSmg_D)
    write(*,*)'Field    bSmg_D = ',bSmg_D
    call get_planet_field(TimeSim,xGsm_D,'GSM NORM',bGsm_D)
    write(*,*)'Field    bGsm_D = ',bGsm_D

    write(*,*)
    write(*,*)'TEST MAP_PLANET_FIELD'
    write(*,*)
    X_D=(/6.0, -8.0, -0.0001/)
    rMap  = 1.0
    write(*,*)'x_D,rMap=',x_D,rMap
    call map_planet_field(0.0,x_D,'SMG NORM',rMap,xMap_D,iHemisphere)
    write(*,*)'xMap_D,iHemisphere=',xMap_D,iHemisphere

  end subroutine test_planet_field

end module CON_planet_field
