!^CFG COPYRIGHT UM
!
!QUOTE: \clearpage
!
!BOP
!
!QUOTE: \section{Planet}
!
!MODULE: CON_planet - planet parameters shared by components
!INTERFACE:
module CON_planet

  !DESCRIPTION:
  ! Physical information about the planet. The planet is described 
  ! with its name. Default values can be set with {\bf planet\_init}.
  ! Simplifying assumptions, such as no rotation, aligned magnetic
  ! and rotational axes etc. can be made.
  !
  ! This is a public class. The variables should be modified by CON only.
  ! Components can only access the data through the inquiry methods 
  ! via the {\bf CON\_physics} class.

  !USES:
  use ModConst
  use ModTimeConvert, ONLY: TimeType

  !REVISION HISTORY:
  ! 01Aug03 - Aaron Ridly <ridley@umich.edu> and 
  !           Gabor Toth <gtoth@umich.edu>   - initial prototype/prolog/code
  ! 23Mar   - added get_planet subroutine for OO type access
  !EOP

  implicit none

  save

  character(len=*), parameter, private :: NameMod='CON_planet'

  integer, parameter    :: lNamePlanet = 40
  character (len=lNamePlanet) :: NamePlanet = 'EARTH'

  real               :: RadiusPlanet     = rEarth
  real               ::   MassPlanet     = mEarth
  real               ::  OmegaPlanet     = OmegaEarth
  real               :: TiltRotation     = TiltEarth
  real               :: IonosphereHeight = IonoHeightEarth
  ! Orbital data
  real           :: OmegaOrbit   = cTwoPi/OrbitalPeriodEarth

  ! Default values are for Earth
  type(TimeType) :: TimeEquinox  = TimeType(&
       iYearEquinoxEarth, iMonthEquinoxEarth, iDayEquinoxEarth, &
       iHourEquinoxEarth, iMinuteEquinoxEarth, iSecondEquinoxEarth, &
       FracSecondEquinoxEarth, &
       0.0_Real8_, '20000320073500')

  real           :: AngleEquinox = AngleEquinoxEarth

  ! Magnetic field type and strength in teslas
  integer, parameter          :: lTypeBField = 40
  character (len=lTypeBField) :: TypeBField = 'DIPOLE'
  real                        :: DipoleStrength = DipoleStrengthEarth
  real    :: MagAxisThetaGeo   = bAxisThetaEarth     ! Permanent theta  in GEO
  real    :: MagAxisPhiGeo     = bAxisPhiEarth       ! Permanent phi    in GEO

  ! Orientation of the axes
  real    :: RotAxisTheta      ! Permanent theta angle in GSE
  real    :: RotAxisPhi        ! Permanent phi   angle in GSE
  real    :: MagAxisTheta      ! Current   theta angla in GSE
  real    :: MagAxisPhi        ! Current   phi   angla in GSE

  ! Optional changes relative to the "real" planet
  logical :: UseRotation     = .true.
  logical :: UseAlignedAxes  = .false.
  logical :: UseRealRotAxis  = .true.
  logical :: UseSetRotAxis   = .false.
  logical :: UseRealMagAxis  = .true.
  logical :: UseSetMagAxis   = .false.

  ! Frequency of updating the magnetic field information
  logical :: DoUpdateB0      = .true.
  real    :: DtUpdateB0      = 0.0001

  ! A primary axis is set to the true value
  ! a secondary axis is aligned with the primary axis
  logical :: IsRotAxisPrimary = .true., IsMagAxisPrimary = .true.

contains

  !BOP ========================================================================
  !IROUTINE: is_planet_init - initialize parameters if planet is known
  !INTERFACE:
  function is_planet_init(NamePlanetIn) result(IsKnown)

    !INPUT ARGUMENTS:
    character(len=*), intent(in) :: NamePlanetIn
    
    !RETURN VALUE:
    logical :: IsKnown

    !DESCRIPTION:
    ! Initialize parameters for the planet identified by its name and
    ! return true if the planet is known. If it is not known return false. 
    ! Store the name in either case. The planet data can be initialized at most
    ! once.
    !EOP
    character (len=*), parameter :: NameSub = NameMod//'::is_planet_init'

    logical :: IsInitialized = .false.
    !-------------------------------------------------------------------------

    if(IsInitialized)then
       if(NamePlanet == NamePlanetIn) RETURN
       call CON_stop(NameSub//&
            ' SWMF_error attempt to change planet name from '// &
            trim(NamePlanet)//' to '//NamePlanetIn)
    end if

    NamePlanet    = NamePlanetIn
    IsInitialized = .true.

    select case(NamePlanet)
    case('EARTH')
       ! This is the default so everything is already set
    case default
       IsKnown = .false.
       RETURN
    end select

    IsKnown = .true.

  end function is_planet_init

  !===========================================================================

  subroutine read_planet_var(NameCommand)

    use ModUtilities, ONLY: upper_case
    use ModReadParam, ONLY: read_var, lStringLine

    character (len=*), intent(in) :: NameCommand

    character (len=*), parameter :: NameSub = NameMod//'::read_planet_var'

    ! Planet related temporary variables
    character (len=lNamePlanet) :: NamePlanetIn
    character (len=lStringLine) :: NamePlanetCommands=''
    logical :: UseNonDipole
    !-------------------------------------------------------------------------

    select case(NameCommand)
    case("#PLANET")

       if (NamePlanetCommands /= '') &
            call CON_stop(NameSub// &
            ' SWMF_ERROR #PLANET should precede '// &
            NamePlanetCommands)

       call read_var('NamePlanet',NamePlanetIn)

       call upper_case(NamePlanetIn)
       if ( .not. is_planet_init(NamePlanetIn) ) then

          call read_var('RadiusPlanet', RadiusPlanet)
          call read_var('MassPlanet',   MassPlanet)
          call read_var('OmegaPlanet',  OmegaPlanet)
          call read_var('TiltRotation', TiltRotation)
          TiltRotation = TiltRotation * cDegToRad
          call read_var('TypeBField',   TypeBField)

          call upper_case(TypeBField)

          select case(TypeBField)

          case('NONE')

             MagAxisTheta   = cZero
             MagAxisPhi     = cZero
             UseSetMagAxis  = .true.
             UseRealMagAxis = .false.

          case('DIPOLE','QUADRUPOLE','OCTUPOLE')

             call read_var('MagAxisThetaGeo', MagAxisThetaGeo)
             MagAxisThetaGeo = MagAxisThetaGeo * cDegToRad
             call read_var('MagAxisPhiGeo',   MagAxisPhiGeo)
             MagAxisPhiGeo = MagAxisPhiGeo * cDegToRad
             call read_var('DipoleStrength',DipoleStrength)

             if (TypeBField == 'QUADRUPOLE') then
                call CON_stop(NameSub// &
                     ' SWMF_ERROR quadrupole field unimplemented')
             endif

             if (TypeBField == 'OCTUPOLE') then
                call CON_stop(NameSub// &
                     ' SWMF_ERROR octupole field unimplemented')
             endif

          case default
             call CON_stop(NameSub// &
                  ' SWMF_ERROR unimplemented TypeBField='// &
                  TypeBField)

          end select

       endif

    case('#IDEALAXES')
       ! This is a short version of setting one axis parallel with Z
       ! and the other one aligned with it

       NamePlanetCommands = '#IDEALAXES ' // NamePlanetCommands

       UseRealRotAxis   = .false.
       IsRotAxisPrimary = .true.
       UseSetRotAxis    = .true.
       RotAxisTheta     = cZero
       RotAxisPhi       = cZero
       UseRealMagAxis   = .false.
       IsMagAxisPrimary = .false.
       UseSetMagAxis    = .false.

    case('#ROTATIONAXIS')

       NamePlanetCommands = '#ROTATIONAXIS ' // NamePlanetCommands
       UseRealRotAxis = .false.

       call read_var('IsRotAxisPrimary', IsRotAxisPrimary)
       if (IsRotAxisPrimary) then

          UseSetRotAxis = .true.

          call read_var('RotAxisTheta', RotAxisTheta)
          if(RotAxisTheta < cZero)call CON_stop(NameSub// &
               'Negative tilt should be entered as Phi=180.0')
          RotAxisTheta = cRadToDeg * RotAxisTheta

          call read_var('RotAxisPhi', RotAxisPhi)
          RotAxisPhi = cRadToDeg * RotAxisPhi
       else
          if(.not.IsMagAxisPrimary)call CON_stop(NameSub//&
               ': SWMF_ERROR ',&
               'Either rotation or magnetic axis must be primary')
       end if

    case('#MAGNETICAXIS')

       NamePlanetCommands = '#MAGNETICAXIS ' // NamePlanetCommands
       UseRealMagAxis = .false.

       call read_var('IsMagAxisPrimary', IsMagAxisPrimary)
       if (IsMagAxisPrimary) then

          UseSetMagAxis = .true.

          call read_var('MagAxisTheta', MagAxisTheta)
          if(MagAxisTheta < 0.0)call CON_stop(NameSub// &
               'Negative tilt should be entered as Phi=180.0')
          MagAxisTheta = cRadToDeg * MagAxisTheta

          call read_var('MagAxisPhi', MagAxisPhi)
          MagAxisPhi = cRadToDeg * MagAxisPhi
       else
          if(.not.IsRotAxisPrimary)call CON_stop(NameSub//&
               ': SWMF_ERROR ',&
               'Either rotation or magnetic axis must be primary')
       end if

    case('#ROTATION')

       NamePlanetCommands = '#ROTATION ' // NamePlanetCommands

       call read_var('UseRotation', UseRotation)
       if (.not.UseRotation) then
          OmegaPlanet  = cZero
       else
          call read_var('Rotation period [hours]',  OmegaPlanet)
          OmegaPlanet = cTwoPi / (3600 * OmegaPlanet)
       endif

    case('#NONDIPOLE')

       NamePlanetCommands = '#NONDIPOLE ' // NamePlanetCommands

       call read_var('UseNonDipole',UseNonDipole)
       if (.not.UseNonDipole) then
          TypeBField = 'DIPOLE'
       else
          call CON_stop(NameSub// &
               ' SWMF_ERROR Nondipole mags unimplemented')
       endif

    case('#DIPOLE')

       NamePlanetCommands = '#DIPOLE ' // NamePlanetCommands

       call read_var('DipoleStrength',DipoleStrength)

    end select

  end subroutine read_planet_var

  !==========================================================================

  subroutine check_planet_var(IsProc0)

    logical, intent(in) :: IsProc0

    character (len=*), parameter :: NameSub=NameMod//'::check_planet_var'

    ! The rotation and magnetic axes are aligned if any of them is not a 
    ! primary axis.
    UseAlignedAxes = (.not. IsRotAxisPrimary) .or. (.not. IsMagAxisPrimary)

    ! Warn if setting is unphysical
    if(UseSetMagAxis .and. UseRealRotAxis .and. IsProc0) &
         write(*,*)NameSub,' WARNING: magnetic axis is explicitly set ',&
         'while rotation axis is calculated from real time ?!'

  end subroutine check_planet_var

  !==========================================================================

  subroutine get_planet( &
       NamePlanetOut, RadiusPlanetOut, IonosphereHeightOut, &
       UseRotationOut, DipoleStrengthOut, DoUpdateB0Out, DtUpdateB0Out)

    character(len=*), optional, intent(out) :: NamePlanetOut
    real,             optional, intent(out) :: RadiusPlanetOut
    real,             optional, intent(out) :: IonosphereHeightOut
    logical,          optional, intent(out) :: UseRotationOut
    real,             optional, intent(out) :: DipoleStrengthOut
    logical,          optional, intent(out) :: DoUpdateB0Out
    real,             optional, intent(out) :: DtUpdateB0Out
    !-----------------------------------------------------------------------
    if(present(NamePlanetOut))      NamePlanetOut     = NamePlanet
    if(present(RadiusPlanetOut))    RadiusPlanetOut   = RadiusPlanet
    if(present(IonosphereHeightOut))IonosphereHeightOut = IonosphereHeight
    if(present(UseRotationOut))     UseRotationOut    = UseRotation
    if(present(DipoleStrengthOut))  DipoleStrengthOut = DipoleStrength
    if(present(DoUpdateB0Out))      DoUpdateB0Out     = DoUpdateB0
    if(present(DtUpdateB0Out))      DtUpdateB0Out     = DtUpdateB0

  end subroutine get_planet

end module CON_planet
