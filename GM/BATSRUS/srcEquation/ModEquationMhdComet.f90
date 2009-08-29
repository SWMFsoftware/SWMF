module ModVarIndexes

  use ModSingleFluid
  use ModExtraVariables

  implicit none

  save

  ! This equation module contains the MHD equations with species for Comets
  character (len=*), parameter :: NameEquation='Cometary MHD'

  ! The variables numbered from 1 to nVar are:
  !
  ! 1. defined in set_ICs.
  ! 2. prolonged and restricted in AMR
  ! 3. saved into the restart file
  ! 4. sent and recieved in the exchange message
  ! 5. filled in the outer ghostcells by the program set_outer_BCs
  ! 5. integrated by subroutine integrate_all for saving to logfile
  ! 6. should be updated by advance_*

  integer, parameter :: nVar = 14  !8 + 6 ion species

  ! Named indexes for State_VGB and other variables
  ! These indexes should go subsequently, from 1 to nVar+1.
  ! The energy is handled as an extra variable, so that we can use
  ! both conservative and non-conservative scheme and switch between them.
  integer, parameter :: &
       Rho_     = 1,    &
       RhoUx_   = 2,    &
       RhoUy_   = 3,    &
       RhoUz_   = 4,    &
       Bx_      = 5,    &
       By_      = 6,    &
       Bz_      = 7,    &
       RhoH2Op_ = 8,    &
       RhoHp_   = 9,    &
       RhoH3Op_ = 10,   &
       RhoOHp_  = 11,   &
       RhoOp_   = 12,   &
       RhoCOp_  = 13,   &
       p_       = nVar, &
       Energy_  = nVar+1

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
       0.0, & ! RhoUx_
       0.0, & ! RhoUy_
       0.0, & ! RhoUz_
       0.0, & ! Bx_
       0.0, & ! By_
       0.0, & ! Bz_
       1.0, & ! H2Op_
       1.0, & ! Hp_
       1.0, & ! H3Op_
       1.0, & ! OHp_
       1.0, & ! Op_
       1.0, & ! COp_
       1.0, & ! p_
       1.0 /) ! Energy_

  ! The names of the variables used in i/o
  character(len=*), parameter :: NameVar_V(nVar+1) = (/ &
       'Rho  ', & ! Rho_
       'RhoUx', & ! RhoUx_
       'RhoUy', & ! RhoUy_
       'RhoUz', & ! RhoUz_
       'Bx   ', & ! Bx_
       'By   ', & ! By_
       'Bz   ', & ! Bz_
       'H2Op ', & ! H2Op_
       'Hp   ', & ! Hp_
       'H3Op ', & ! H3Op_
       'OHp  ', & ! OHp_
       'Op   ', & ! Op_
       'COp  ', & ! COp_
       'p    ', & ! p_
       'e    ' /) ! Energy_

  ! The space separated list of nVar conservative variables for plotting
  character(len=*), parameter :: NameConservativeVar = &
       'rho mx my mz bx by bz H2Op Hp H3Op OHp Op COp e'

  ! The space separated list of nVar primitive variables for plotting
  character(len=*), parameter :: NamePrimitiveVar = &
       'rho ux uy uz bx by bz H2Op Hp H3Op OHp Op COp p'

  ! The space separated list of nVar primitive variables for TECplot output
  character(len=*), parameter :: NamePrimitiveVarTec = &
       '"`r", "U_x", "U_y", "U_z", "B_x", "B_y", "B_z", "H2Op", "Hp", "H3Op", "OHp", "Op", "COp", "p"'

  ! Names of the user units for IDL and TECPlot output
  character(len=20) :: &
       NameUnitUserIdl_V(nVar+1) = '', NameUnitUserTec_V(nVar+1) = ''

  ! The user defined units for the variables
  real :: UnitUser_V(nVar+1) = 1.0

  ! Primitive variable names
  integer, parameter :: U_ = RhoU_, Ux_ = RhoUx_, Uy_ = RhoUy_, Uz_ = RhoUz_

  ! There are no extra scalars
  integer, parameter :: ScalarFirst_ = RhoH2Op_, ScalarLast_ = RhoCOp_

  ! Species
  logical, parameter :: UseMultiSpecies = .true.
  integer, parameter :: SpeciesFirst_   = ScalarFirst_
  integer, parameter :: SpeciesLast_    = ScalarLast_

  ! Molecular mass of species H, O2, O, CO2 in AMU:
  real, parameter :: MassSpecies_V(SpeciesFirst_:SpeciesLast_) = &
       (/18.0,1.0,19.0,17.0,16.0,28.0/)
contains

  subroutine init_mod_equation

    call init_mhd_variables

  end subroutine init_mod_equation

end module ModVarIndexes
