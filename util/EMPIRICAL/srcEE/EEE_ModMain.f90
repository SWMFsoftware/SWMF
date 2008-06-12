!^CFG COPYRIGHT UM
!==============================================================================
Module EEE_ModMain
  use EEE_ModGetStateInit, ONLY: EEE_get_state_init
  use EEE_ModGetStateBC,   ONLY: EEE_get_state_BC
  use EEE_ModGetB0,        ONLY: EEE_get_B0
  implicit none
  save

  private

  public :: EEE_initialize
  public :: EEE_set_parameters
  public :: EEE_get_state_init
  public :: EEE_get_state_BC
  public :: EEE_get_B0

contains

  subroutine EEE_initialize(BodyNDim,BodyTDim,gamma)
    use EEE_ModCommonVariables
    implicit none

    real, intent(in) :: BodyNDim,BodyTDim,gamma

    integer :: iComm, iError
    !--------------------------------------------------------------------------

    iComm = MPI_COMM_WORLD
    call MPI_COMM_RANK(iComm,iProc,iError)

    g = gamma
    inv_g = 1.0/g
    gm1 = g - 1.0
    inv_gm1 = 1.0/(g - 1.0)

    ! assume MassIon_I(1) = 1.0
    No2Si_V(UnitX_)   = rSun
    No2Si_V(UnitU_)   = sqrt(g*cBoltzmann*BodyTDim/cProtonMass)
    No2Si_V(UnitRho_) = 1000000*cProtonMass*BodyNDim

    !\
    ! Set other normalizing SI variables from the independent ones.
    !
    ! For sake of convenience
    !  units of B are chosen to satisfy v_A = B/sqrt(rho)       (mu = 1)
    !  units of n are chosen to satisfy  n  = rho/(ionmass/amu) (mp = 1)
    !  units of T are chosen to satisfy  T  = p/n               (kBoltzmann = 1)
    !
    ! Note that No2Si_V(UnitN_) is NOT EQUAL TO 1/No2Si_V(UnitX_)^3 !!!
    !/
    No2Si_V(UnitT_)          = No2Si_V(UnitX_)/No2Si_V(UnitU_)         ! s
    No2Si_V(UnitN_)          = No2Si_V(UnitRho_)/cProtonMass           ! #/m^3
    No2Si_V(UnitP_)          = No2Si_V(UnitRho_)*No2Si_V(UnitU_)**2    ! Pa
    No2Si_V(UnitB_)          = No2Si_V(UnitU_) &
         *sqrt(cMu*No2Si_V(UnitRho_))                                  ! T
    No2Si_V(UnitRhoU_)       = No2Si_V(UnitRho_)*No2Si_V(UnitU_)       ! kg/m^2/s
    No2Si_V(UnitEnergyDens_) = No2Si_V(UnitP_)                         ! J/m^3
    No2Si_V(UnitPoynting_)   = No2Si_V(UnitEnergyDens_)*No2Si_V(UnitU_)! J/m^2/s
    No2Si_V(UnitJ_)          = No2Si_V(UnitB_)/( No2Si_V(UnitX_)*cMu ) ! A/m^2
    No2Si_V(UnitElectric_)   = No2Si_V(UnitU_)*No2Si_V(UnitB_)         ! V/m
    No2Si_V(UnitTemperature_)= No2Si_V(UnitP_) &
         /( No2Si_V(UnitN_)*cBoltzmann )                               ! K 
    No2Si_V(UnitDivB_)       = No2Si_V(UnitB_)/No2Si_V(UnitX_)         ! T/m
    No2Si_V(UnitAngle_)      = 1.0                                     ! radian

    !\
    ! Set inverse conversion SI -> normalized
    !/
    Si2No_V = 1.0/No2Si_V

    ! As a default use SI units, so below only the differences need to be set
    Io2Si_V = 1.0
    No2Io_V = No2Si_V

    Io2Si_V(UnitX_)           = rSun                      ! R
    Io2Si_V(UnitRho_)         = 1.0E+3                    ! g/cm^3
    Io2Si_V(UnitN_)           = 1.0E+6                    ! #/cm^3
    Io2Si_V(UnitU_)           = 1.0E+3                    ! km/s
    Io2Si_V(UnitP_)           = 1.0E-1                    ! dyne/cm^2
    Io2Si_V(UnitB_)           = 1.0E-4                    ! Gauss
    Io2Si_V(UnitRhoU_)        = 1.0E+1                    ! g/cm^2/s
    Io2Si_V(UnitEnergydens_)  = 1.0E-1                    ! erg/cm^3
    Io2Si_V(UnitJ_)           = 1.0E-6                    ! uA/m^2
    Io2Si_V(UnitDivB_)        = 1.0E-2                    ! Gauss/cm
    Io2Si_V(UnitAngle_)       = cRadToDeg                 ! degrees

    ! Calculate the remaining unit conversions
    Si2Io_V = 1/Io2Si_V
    No2Io_V = No2Si_V*Si2Io_V
    Io2No_V = 1/No2Io_V

    Gbody  = -cGravitation*mSun*(Si2No_V(UnitU_)**2 * Si2No_V(UnitX_))

  end subroutine EEE_initialize

  !============================================================================

  subroutine EEE_set_parameters(NameCommand)
    use ModReadParam,     ONLY: read_var
    use EEE_ModGL98,      ONLY: set_parameters_GL98
    use EEE_ModTD99,      ONLY: set_parameters_TD99
    use EEE_ModArch,      ONLY: set_parameters_arch
    use EEE_ModShearFlow, ONLY: set_parameters_shearflow
    implicit none

    character(len=*), intent(in) :: NameCommand
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case("#TD99FLUXROPE")
       call set_parameters_TD99
    case("#GL98FLUXROPE")
       call set_parameters_GL98
    case("#ARCH")
       call set_parameters_arch
    case("#SHEARFLOW")
       call set_parameters_shearflow
    end select

  end subroutine EEE_set_parameters

end Module EEE_ModMain
