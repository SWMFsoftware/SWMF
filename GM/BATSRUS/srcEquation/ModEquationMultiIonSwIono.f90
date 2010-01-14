module ModVarIndexes

  use ModExtraVariables

  implicit none

  save

  ! This equation module declares two Ions (solar wind H and ionosphere H)
  ! for the purposes of tracking plasma entry into the magnetosphere in
  ! multifluid MHD.  It requires changes to the #MAGNETOSPHERE command
  ! in order to set the Iono species as the dominant outflow species at
  ! the inner boundary.  
  ! Based on ModEquationMultiIon and ModEquationMhdSwIono(Hansen, 2006)
  ! Welling, 2010.

  ! This equation module contains the standard MHD equations.
  character (len=*), parameter :: NameEquation='Multi-ion (Sw and Iono H+)'

  integer, parameter :: nVar = 13

  ! There are two ion fluids but no total ion fluid
  integer, parameter :: nFluid    = 2
  integer, parameter :: IonFirst_ = 1
  integer, parameter :: IonLast_  = 2
  logical, parameter :: IsMhd     = .false.
  real               :: MassFluid_I(nFluid) = (/ 1.0, 1.0 /) !both are H+.

  character (len=3), parameter :: NameFluid_I(nFluid) = (/ 'Sw', 'Iono' /)

  ! Named indexes for State_VGB and other variables
  ! These indexes should go subsequently, from 1 to nVar+nFluid.
  ! The energies are handled as an extra variable, so that we can use
  ! both conservative and non-conservative scheme and switch between them.
  integer, parameter ::            &
       Rho_         =  1,          &
       RhoUx_       =  2, Ux_ = 2, &
       RhoUy_       =  3, Uy_ = 3, &
       RhoUz_       =  4, Uz_ = 4, &
       Bx_          =  5,          &
       By_          =  6,          &
       Bz_          =  7,          &
       p_           =  8,          &
       IonoRho_     =  9,          &
       IonoRhoUx_   = 10,          &
       IonoRhoUy_   = 11,          &
       IonoRhoUz_   = 12,          &
       IonoP_       = 13,          &
       Energy_      = nVar+1,      &
       IonoEnergy_  = nVar+2

  ! This allows to calculate RhoUx_ as RhoU_+x_ and so on.
  integer, parameter :: U_ = Ux_ - 1, RhoU_ = RhoUx_-1, B_ = Bx_-1

  ! These arrays are useful for multifluid
  integer, parameter :: iRho_I(nFluid)   = (/Rho_,   IonoRho_/)
  integer, parameter :: iRhoUx_I(nFluid) = (/RhoUx_, IonoRhoUx_/)
  integer, parameter :: iRhoUy_I(nFluid) = (/RhoUy_, IonoRhoUy_/)
  integer, parameter :: iRhoUz_I(nFluid) = (/RhoUz_, IonoRhoUz_/)
  integer, parameter :: iP_I(nFluid)     = (/p_,     IonoP_/)

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
       1.0, & ! IonoRho_
       0.0, & ! IonoRhoUx_
       0.0, & ! IonoRhoUy_
       0.0, & ! IonoRhoUz_
       1.0, & ! IonoP_
       1.0, & ! Energy_
       1.0 /) ! IonoEnergy_

  ! The names of the variables used in i/o
  character(len=*), parameter :: NameVar_V(nVar+nFluid) = (/ &
       'Rho  ',   & ! Rho_
       'Mx   ',   & ! RhoUx_
       'My   ',   & ! RhoUy_
       'Mz   ',   & ! RhoUz_
       'Bx   ',   & ! Bx_
       'By   ',   & ! By_
       'Bz   ',   & ! Bz_
       'P    ',   & ! p_
       'IonoRho', & ! IonoRho_
       'IonoMx ', & ! IonoRhoUx_
       'IonoMy ', & ! IonoRhoUy_
       'IonoMz ', & ! IonoRhoUz_
       'IonoP  ', & ! IonoP_
       'SwE  ',   & ! Energy_
       'IonoE  ' /) ! IonoEnergy_

  ! The space separated list of nVar conservative variables for plotting
  character(len=*), parameter :: NameConservativeVar = &
       'Rho Mx My Mz Bx By Bz E IonoRho IonoMx IonoMy IonoMz IonoE'

  ! The space separated list of nVar primitive variables for plotting
  character(len=*), parameter :: NamePrimitiveVar = &
       'Rho Ux Uy Uz Bx By Bz P IonoRho IonoUx IonoUy IonoUz IonoP'

  ! The space separated list of nVar primitive variables for TECplot output
  character(len=*), parameter :: NamePrimitiveVarTec = &
       '"`r^Sw^+", "U_x^Sw^+", "U_y^Sw^+", "U_z^Sw^+", "B_x", "B_y", ' // &
       '"B_z", "p", "`r^Iono^+", "U_x^Iono^+", "U_y^Iono^+", ' // &
       '"U_z^Iono^+", "P^Iono^+"'

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
