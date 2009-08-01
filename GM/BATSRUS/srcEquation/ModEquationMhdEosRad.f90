module ModVarIndexes
  use ModSingleFluid, Redefine1 => Erad_, Redefine2 => ExtraEint_
  implicit none

  save

  ! This equation module contains the standard MHD equations with
  ! an additional energy deficit for the equation of state and radiation
  character (len=*), parameter :: NameEquation='MHD+eos+radiation'

  integer, parameter :: nVar = 10

  ! Named indexes for State_VGB and other variables
  ! These indexes should go subsequently, from 1 to nVar+1.
  ! The energy is handled as an extra variable, so that we can use
  ! both conservative and non-conservative scheme and switch between them.
  integer, parameter :: &
       Rho_       = 1,          &
       RhoUx_     = 2, Ux_ = 2, &
       RhoUy_     = 3, Uy_ = 3, &
       RhoUz_     = 4, Uz_ = 4, &
       Bx_        = 5,          &
       By_        = 6,          &
       Bz_        = 7,          &
       Erad_      = 8,          &
       ExtraEint_ = 9,          &
       p_         = nVar,       &
       Energy_    = nVar+1

  ! This allows to calculate RhoUx_ as RhoU_+x_ and so on.
  integer, parameter :: U_ = Ux_ - 1, RhoU_ = RhoUx_ - 1, B_ = Bx_ - 1

  ! These arrays are useful for multifluid
  integer, parameter :: iRho_I(nFluid)   = (/Rho_/)
  integer, parameter :: iRhoUx_I(nFluid) = (/RhoUx_/)
  integer, parameter :: iRhoUy_I(nFluid) = (/RhoUy_/)
  integer, parameter :: iRhoUz_I(nFluid) = (/RhoUz_/)
  integer, parameter :: iP_I(nFluid)     = (/p_/)

  ! The default values for the state variables:
  ! Variables which are physically positive should be set to 1,
  ! variables that can be positive or negative should be set to 0:
  real, parameter :: DefaultState_V(nVar+nFluid) = (/ & 
       1.0, & ! Rho_
       0.0, & ! RhoUx_
       0.0, & ! RhoUy_
       0.0, & ! RhoUz_
       0.0, & ! Bx_
       0.0, & ! By_
       0.0, & ! Bz_
       0.0, & ! Erad_
       0.0, & ! ExtraEint_
       1.0, & ! p_
       1.0 /) ! Energy_

  ! The names of the variables used in i/o
  character(len=*), parameter :: NameVar_V(nVar+nFluid) = (/ &
       'Rho ', & ! Rho_
       'Mx  ', & ! RhoUx_
       'My  ', & ! RhoUy_
       'Mz  ', & ! RhoUz_
       'Bx  ', & ! Bx_
       'By  ', & ! By_
       'Bz  ', & ! Bz_
       'Erad', & ! Erad_
       'EInt', & ! ExtraEint_
       'P   ', & ! p_
       'E   '/)  ! Energy_

  ! The space separated list of nVar conservative variables for plotting
  character(len=*), parameter :: NameConservativeVar = &
       'rho mx my mz bx by bz Erad EInt e'

  ! The space separated list of nVar primitive variables for plotting
  character(len=*), parameter :: NamePrimitiveVar = &
       'Rho Ux Uy Uz bx by bz Erad EInt p'

  ! The space separated list of nVar primitive variables for TECplot output
  character(len=*), parameter :: NamePrimitiveVarTec = &
       '"`r", "U_x", "U_y", "U_z", "B_x", "B_y", "B_z", "Erad", "EInt", "p"'

  ! Names of the user units for IDL and TECPlot output
  character(len=20) :: &
       NameUnitUserIdl_V(nVar+nFluid) = '', NameUnitUserTec_V(nVar+nFluid) = ''

  ! The user defined units for the variables
  real :: UnitUser_V(nVar+nFluid) = 1.0

  ! The only scalar to be advected is the radiation energy density
  integer, parameter :: ScalarFirst_ = Erad_, ScalarLast_ = ExtraEint_

  ! There are no multi-species
  logical, parameter :: UseMultiSpecies = .false.

  ! Declare the following variables to satisfy the compiler
  integer, parameter :: SpeciesFirst_ = 1, SpeciesLast_ = 1
  real               :: MassSpecies_V(SpeciesFirst_:SpeciesLast_)

contains

  subroutine init_mod_equation

    call init_mhd_variables

    ! Set the unit and unit name for the wave energy variable
    UnitUser_V(Erad_)        = UnitUser_V(Energy_)
    NameUnitUserTec_V(Erad_) = NameUnitUserTec_V(Energy_)
    NameUnitUserIdl_V(Erad_) = NameUnitUserIdl_V(Energy_)
    UnitUser_V(ExtraEint_)        = UnitUser_V(Energy_)
    NameUnitUserTec_V(ExtraEint_) = NameUnitUserTec_V(Energy_)
    NameUnitUserIdl_V(ExtraEint_) = NameUnitUserIdl_V(Energy_)

  end subroutine init_mod_equation

end module ModVarIndexes
