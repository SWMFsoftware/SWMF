!^CFG COPYRIGHT UM
!
!QUOTE: \clearpage
!
!BOP
!
!QUOTE: \section{Shared Methods for Components and CON}
!
!MODULE: CON_physics - physics parameters accessible by components
!INTERFACE:
module CON_physics

  !USES:

  use CON_time
  use CON_axes
  use CON_planet
  use CON_planet_field

  implicit none

  private ! except
  
  !PUBLIC DATA MEMBERS:
  public :: lNamePlanet      ! length of NamePlanet string
  public :: lTypeBField      ! length of TypeBField string

  !PUBLIC MEMBER FUNCTIONS:
  public :: get_physics      ! get physics parameters
  public :: get_planet_field ! get planet field at a point and time
  public :: map_planet_field ! map planet field to some radius
  public :: get_axes         ! get information about (coordinate) axes
  public :: transform_matrix ! return transformation matrix between two systems
  public :: time_real_to_int ! transform double precision seconds into date
  public :: time_int_to_real ! transform date into double precision seconds
  public :: n_day_of_year    ! retrun day of year calculated from date

  !REVISION HISTORY:
  ! 14Aug03 - Gabor Toth <gtoth@umich.edu> - initial prototype/prolog/code
  !EOP

  save

  character(len=*), parameter, private :: NameMod='CON_physics'

contains
  !BOP
  !IROUTINE: get_physics - get physics parameters
  !INTERFACE:
  subroutine get_physics(&
       NamePlanetOut,RadiusPlanetOut,IonosphereHeightOut,&
       DipoleStrengthOut,TypeBFieldOut,&
       DtUpdateB0Out,&
       UseRotationOut,&
       DoTimeAccurateOut,tSimulationOut,TimeStartOut,TimeCurrentOut,&
       tStartOut, tCurrentOut, nStepOut)

    !OUTPUT ARGUMENTS:
    character(len=*), optional, intent(out) :: NamePlanetOut
    real,             optional, intent(out) :: RadiusPlanetOut
    real,             optional, intent(out) :: IonosphereHeightOut
    character(len=*), optional, intent(out) :: TypeBFieldOut
    real,             optional, intent(out) :: DipoleStrengthOut
    real,             optional, intent(out) :: DtUpdateB0Out
    logical,          optional, intent(out) :: UseRotationOut
    logical,          optional, intent(out) :: DoTimeAccurateOut
    real,             optional, intent(out) :: tSimulationOut
    type(TimeType),   optional, intent(out) :: TimeStartOut
    real(Real8_),     optional, intent(out) :: tStartOut
    type(TimeType),   optional, intent(out) :: TimeCurrentOut
    real(Real8_),     optional, intent(out) :: tCurrentOut
    integer,          optional, intent(out) :: nStepOut
    !EOP
    !-------------------------------------------------------------------------

    if(present(NamePlanetOut))     NamePlanetOut     = NamePlanet
    if(present(RadiusPlanetOut))   RadiusPlanetOut   = RadiusPlanet
    if(present(IonosphereHeightOut))IonosphereHeightOut = IonosphereHeight
    if(present(TypeBFieldOut))     TypeBFieldOut     = TypeBField
    if(present(DipoleStrengthOut)) DipoleStrengthOut = DipoleStrength
    if(present(DtUpdateB0Out))     DtUpdateB0Out     = DtUpdateB0
    if(present(UseRotationOut))    UseRotationOut    = UseRotation
    if(present(DoTimeAccurateOut)) DoTimeAccurateOut = DoTimeAccurate
    if(present(tSimulationOut))    tSimulationOut    = tSimulation
    if(present(tSimulationOut))    tSimulationOut    = tSimulation
    if(present(TimeCurrentOut))    TimeCurrentOut    = TimeCurrent
    if(present(tCurrentOut))       tCurrentOut       = TimeCurrent % Time
    if(present(TimeStartOut))      TimeStartOut      = TimeStart
    if(present(tStartOut))         tStartOut         = TimeStart % Time
    if(present(nStepOut))          nStepOut          = nStep

  end subroutine get_physics

end module CON_physics
