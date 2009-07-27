module ModVarIndexes

  use ModExtraVariables

  implicit none

  save

  ! This equation module contains the standard MHD equations.
  character (len=*), parameter :: NameEquation='MHD with 2 ions'

  integer, parameter :: nVar = 18

  integer, parameter :: nFluid    = 3
  integer, parameter :: IonFirst_ = 2
  integer, parameter :: IonLast_  = 3
  logical, parameter :: IsMhd     = .true.
  real               :: MassFluid_I(2:3) = (/ 1.0, 16.0 /)

  character (len=3), parameter :: NameFluid_I(nFluid)= (/ 'All', 'Hp ', 'Op '/)

  ! Named indexes for State_VGB and other variables
  ! These indexes should go subsequently, from 1 to nVar+nFluid.
  ! The energies are handled as an extra variable, so that we can use
  ! both conservative and non-conservative scheme and switch between them.
  integer, parameter :: &
       Rho_       =  1,          &
       RhoUx_     =  2, Ux_ = 2, &
       RhoUy_     =  3, Uy_ = 3, &
       RhoUz_     =  4, Uz_ = 4, &
       Bx_        =  5, &
       By_        =  6, &
       Bz_        =  7, &
       p_         =  8, &
       HpRho_     =  9, &
       HpRhoUx_   = 10, &
       HpRhoUy_   = 11, &
       HpRhoUz_   = 12, &
       HpP_       = 13, &
       OpRho_     = 14, &
       OpRhoUx_   = 15, &
       OpRhoUy_   = 16, &
       OpRhoUz_   = 17, &
       OpP_       = 18, &
       Energy_    = nVar+1, &
       HpEnergy_  = nVar+2, &
       OpEnergy_  = nVar+3

  ! This allows to calculate RhoUx_ as RhoU_+x_ and so on.
  integer, parameter :: U_ = Ux_ - 1, RhoU_ = RhoUx_-1, B_ = Bx_-1

  ! These arrays are useful for multifluid
  integer, parameter :: iRho_I(nFluid)   = (/Rho_,   HpRho_,   OpRho_/)
  integer, parameter :: iRhoUx_I(nFluid) = (/RhoUx_, HpRhoUx_, OpRhoUx_/)
  integer, parameter :: iRhoUy_I(nFluid) = (/RhoUy_, HpRhoUy_, OpRhoUy_/)
  integer, parameter :: iRhoUz_I(nFluid) = (/RhoUz_, HpRhoUz_, OpRhoUz_/)
  integer, parameter :: iP_I(nFluid)     = (/p_,     HpP_,     OpP_/)

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
       1.0, & ! p_
       1.0, & ! HpRho_
       0.0, & ! HpRhoUx_
       0.0, & ! HpRhoUy_
       0.0, & ! HpRhoUz_
       1.0, & ! HpP_
       1.0, & ! OpRho_
       0.0, & ! OpRhoUx_
       0.0, & ! OpRhoUy_
       0.0, & ! OpRhoUz_
       1.0, & ! OpP_
       1.0, & ! Energy_
       1.0, & ! HpEnergy_
       1.0  /)! OpEnergy_

  ! The names of the variables used in i/o
  character(len=*), parameter :: NameVar_V(nVar+nFluid) = (/ &
       'Rho  ', & ! Rho_
       'Mx   ', & ! RhoUx_
       'My   ', & ! RhoUy_
       'Mz   ', & ! RhoUz_
       'Bx   ', & ! Bx_
       'By   ', & ! By_
       'Bz   ', & ! Bz_
       'P    ', & ! p_
       'HpRho', & ! HpRho_
       'HpMx ', & ! HpRhoUx_
       'HpMy ', & ! HpRhoUy_
       'HpMz ', & ! HpRhoUz_
       'HpP  ', & ! HpP_
       'OpRho', & ! OpRho_
       'OpMx ', & ! OpRhoUx_
       'OpMy ', & ! OpRhoUy_
       'OpMz ', & ! OpRhoUz_
       'OpP  ', & ! OpP_
       'E    ', & ! Energy_
       'HpE  ', & ! HpEnergy_
       'OpE  ' /) ! OpEnergy_

  ! The space separated list of nVar conservative variables for plotting
  character(len=*), parameter :: NameConservativeVar = &
       'Rho Mx My Mz Bx By Bz E '// &
       'HpRho HpMx HpMy HpMz HpE '// &
       'OpRho OpMx OpMy OpMz OpE'

  ! The space separated list of nVar primitive variables for plotting
  character(len=*), parameter :: NamePrimitiveVar = &
       'Rho Ux Uy Uz Bx By Bz P '// &
       'HpRho HpUx HpUy HpUz HpP '// &
       'OpRho OpUx OpUy OpUz OpP'

  ! The space separated list of nVar primitive variables for TECplot output
  character(len=*), parameter :: NamePrimitiveVarTec = &
       '"`r", "U_x", "U_y", "U_z", "B_x", "B_y", "B_z", "p", ' // &
       '"`r^H^+", "U_x^H^+", "U_y^H^+", "U_z^H^+", "P^H^+", ' // &
       '"`r^O^+", "U_x^O^+", "U_y^O^+", "U_z^O^+", "P^O^+"'

  ! Names of the user units for IDL and TECPlot output
  character(len=20) :: &
       NameUnitUserIdl_V(nVar+nFluid) = '', NameUnitUserTec_V(nVar+nFluid) = ''

  ! The user defined units for the variables
  real :: UnitUser_V(nVar+nFluid) = 1.0

  ! There are no extra scalars
  integer, parameter :: ScalarFirst_ = 2, ScalarLast_ = 1

  ! There are no multi-species
  logical, parameter :: UseMultiSpecies = .false.

  ! Declare the following variables to satisfy the compiler
  integer, parameter :: SpeciesFirst_ = 1, SpeciesLast_ = 1
  real               :: MassSpecies_V(SpeciesFirst_:SpeciesLast_)

contains

  subroutine init_mod_equation

    call init_mhd_variables

  end subroutine init_mod_equation

end module ModVarIndexes
