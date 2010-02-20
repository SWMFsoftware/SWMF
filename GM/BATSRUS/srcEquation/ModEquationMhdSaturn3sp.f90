module ModVarIndexes

  use ModSingleFluid
  use ModExtraVariables

  implicit none

  save

  ! This equation module contains the standard MHD equations with
  ! three species for Saturn.  1 - solar wind protons+ionosphere, 
  ! 2 - water group plasma from the rings and Enceladus, 
  ! 3 - nitrogen group plasma from Titan
  character (len=*), parameter :: NameEquation= &
       'Saturn MHD 3 Species (Saturn3sp), Hansen, May, 2007'

  ! The variables numbered from 1 to nVar are:
  !
  ! 1. defined in set_ICs.
  ! 2. prolonged and restricted in AMR
  ! 3. saved into the restart file
  ! 4. sent and recieved in the exchange message
  ! 5. filled in the outer ghostcells by the program set_outer_BCs
  ! 5. integrated by subroutine integrate_all for saving to logfile
  ! 6. should be updated by advance_*

  integer, parameter :: nVar = 11

  ! Named indexes for State_VGB and other variables
  ! These indexes should go subsequently, from 1 to nVar+1.
  ! The energy is handled as an extra variable, so that we can use
  ! both conservative and non-conservative scheme and switch between them.
  integer, parameter :: &
       Rho_    = 1,    &
       RhoH_   = 2,    &
       RhoH2O_ = 3,    &
       RhoN_   = 4,    &
       RhoUx_  = 5,    &
       RhoUy_  = 6,    &
       RhoUz_  = 7,    &
       Bx_     = 8,    &
       By_     = 9,    &
       Bz_     = 10,    &
       p_      = nVar, &
       Energy_ = nVar+1

  ! This allows to calculate RhoUx_ as rhoU_+x_ and so on.
  integer, parameter :: RhoU_ = RhoUx_-1, B_ = Bx_-1

  ! These arrays are useful for multifluid
  integer, parameter :: iRho_I(nFluid)   = (/Rho_/)
  integer, parameter :: iRhoUx_I(nFluid) = (/RhoUx_/)
  integer, parameter :: iRhoUy_I(nFluid) = (/RhoUy_/)
  integer, parameter :: iRhoUz_I(nFluid) = (/RhoUz_/)
  integer, parameter :: iP_I(nFluid)     = (/p_/)

  ! The default values for the state variables:
  ! Variables which are physically positive should be set to 1,
  ! variables that can be positive or negative should be set to 0:
  real, parameter :: DefaultState_V(nVar+1) = (/ & 
       1.0, & ! Rho_
       1.0, & ! RhoH_
       1.0, & ! RhoH2O_
       1.0, & ! RhoN_
       0.0, & ! RhoUx_
       0.0, & ! RhoUy_
       0.0, & ! RhoUz_
       0.0, & ! Bx_
       0.0, & ! By_
       0.0, & ! Bz_
       1.0, & ! p_
       1.0 /) ! Energy_

  ! The names of the variables used in i/o
  character(len=6) :: NameVar_V(nVar+1) = (/ &
       'Rho   ', & ! Rho_
       'RhoH  ', & ! RhoH_
       'RhoH2O', & ! RhoH2O_
       'RhoN  ', & ! RhoN_
       'Mx    ', & ! RhoUx_
       'My    ', & ! RhoUy_
       'Mz    ', & ! RhoUz_
       'Bx    ', & ! Bx_
       'By    ', & ! By_
       'Bz    ', & ! Bz_
       'p     ', & ! p_
       'e     ' /) ! Energy_

  ! The space separated list of nVar conservative variables for plotting
  character(len=*), parameter :: NameConservativeVar = &
       'rho RhoH RhoH2O RhoN mx my mz bx by bz e'

  ! The space separated list of nVar primitive variables for plotting
  character(len=*), parameter :: NamePrimitiveVar = &
       'rho RhoH RhoH2O RhoN ux uy uz bx by bz p'

  ! The space separated list of nVar primitive variables for TECplot output
  character(len=*), parameter :: NamePrimitiveVarTec = &
       '"`r", "RhoH", "RhoH2O", "RhoN", "U_x", "U_y", "U_z", "B_x", "B_y", "B_z", "p"'

  ! Names of the user units for IDL and TECPlot output
  character(len=20) :: &
       NameUnitUserIdl_V(nVar+1) = '', NameUnitUserTec_V(nVar+1) = ''

  ! The user defined units for the variables
  real :: UnitUser_V(nVar+1) = 1.0

  ! Primitive variable names
  integer, parameter :: U_ = RhoU_, Ux_ = RhoUx_, Uy_ = RhoUy_, Uz_ = RhoUz_

  ! There are three extra scalars
  integer, parameter :: ScalarFirst_ = RhoH_, ScalarLast_ = RhoN_

  ! Species
  logical, parameter :: UseMultiSpecies = .true.
  integer, parameter :: SpeciesFirst_   = ScalarFirst_
  integer, parameter :: SpeciesLast_    = ScalarLast_

  ! Molecular mass of solarwind, H2O and N species in AMU:
  real :: MassSpecies_V(SpeciesFirst_:SpeciesLast_) = (/ 1.0, 16.6, 14.0 /)

contains

  subroutine init_mod_equation

    call init_mhd_variables

  end subroutine init_mod_equation

end module ModVarIndexes
