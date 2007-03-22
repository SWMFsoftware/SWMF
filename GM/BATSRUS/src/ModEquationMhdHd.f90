module ModVarIndexes

  use ModFluid
  implicit none

  save

  ! This equation module contains the standard MHD equations.
  character (len=*), parameter :: NameEquation='MHD and HD'

  ! The variables numbered from 1 to nVar are:
  !
  ! 1. defined in set_ICs.
  ! 2. prolonged and restricted in AMR
  ! 3. saved into the restart file
  ! 4. sent and recieved in the exchange message
  ! 5. filled in the outer ghostcells by set_outer_bcs
  ! 5. integrated by subroutine integrate_all for saving to logfile
  ! 6. should be updated by advance_*

  integer, parameter :: nVar = 13

  integer, parameter :: nFluid = 2

  character (len=3), parameter :: NameFluid_I(nFluid) = (/ 'Ion'    , 'Neu' /)
  character (len=7), parameter :: TypeFluid_I(nFluid) = (/ 'ion    ', 'neutral' /)

  ! Named indexes for State_VGB and other variables
  ! These indexes should go subsequently, from 1 to nVar+1.
  ! The energy is handled as an extra variable, so that we can use
  ! both conservative and non-conservative scheme and switch between them.
  integer, parameter :: &
       NeuRho_    = Bz_+1,  &
       Energy_    = nVar+1, &
       NeuEnergy_ = nVar+2

  integer, parameter :: iVarFluid_I(nFluid)  = (/ 0, NeuRho_ - 1 /)

  ! The default values for the state variables:
  ! Variables which are physically positive should be set to 1,
  ! variables that can be positive or negative should be set to 0:
  real, parameter :: DefaultState_V(nVar+nFluid) = (/ & 
       1.0, & ! Rho_
       0.0, & ! RhoUx_
       0.0, & ! RhoUy_
       0.0, & ! RhoUz_
       1.0, & ! p_
       0.0, & ! Bx_
       0.0, & ! By_
       0.0, & ! Bz_
       1.0, & ! NeuRho_
       0.0, & ! NeuRhoUx_
       0.0, & ! NeuRhoUy_
       0.0, & ! NeuRhoUz_
       1.0, & ! NeuP_
       1.0, & ! Energy_
       1.0 /) ! NeuEnergy_

  ! The names of the variables used in i/o
  character(len=*), parameter :: NameVar_V(nVar+nFluid) = (/ &
       'Rho     ', & ! Rho_
       'RhoUx   ', & ! RhoUx_
       'RhoUy   ', & ! RhoUy_
       'RhoUz   ', & ! RhoUz_
       'p       ', & ! p_
       'Bx      ', & ! Bx_
       'By      ', & ! By_
       'Bz      ', & ! Bz_
       'NeuRho  ', & ! NeuRho_
       'NeuRhoUx', & ! NeuRhoUx_
       'NeuRhoUy', & ! NeuRhoUy_
       'NeuRhoUz', & ! NeuRhoUz_
       'NeuP    ', & ! NeuP_
       'e       ', & ! Energy_
       'NeuE    ' /) ! NeuEnergy_

  ! The space separated list of nVar conservative variables for plotting
  character(len=*), parameter :: NameConservativeVar = &
       'Rho RhoUx RhoUy RhoUz Bx By Bz e NeuRho NeuRhoUx NeuRhoUy NeuRhoUz NeuE'

  ! The space separated list of nVar primitive variables for plotting
  character(len=*), parameter :: NamePrimitiveVar = &
       'Rho Ux Uy Uz Bx By Bz p NeuRho NeuUx NeuUy NeuUz NeuP'

  ! The space separated list of nVar primitive variables for TECplot output
  character(len=*), parameter :: NamePrimitiveVarTec = &
       '"`r", "U_x", "U_y", "U_z", "B_x", "B_y", "B_z", "p" "`r^n", "U^n_x", "U^n_y", "U^n_z", "P^n" '

  ! Names of the user units for IDL and TECPlot output
  character(len=20) :: &
       NameUnitUserIdl_V(nVar+1) = '', NameUnitUserTec_V(nVar+1) = ''

  ! The user defined units for the variables
  real :: UnitUser_V(nVar+nFluid) = 1.0

  ! Named indexes for corrected fluxes
  integer, parameter :: Vdt_ = nVar+1
  integer, parameter :: BnL_ = nVar+2
  integer, parameter :: BnR_ = nVar+3
  integer, parameter :: nCorrectedFaceValues = BnR_

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
