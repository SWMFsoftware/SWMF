module ModVarIndexes

  use ModExtraVariables

  implicit none

  save
  ! This equation module declares two Ions (solar wind H and ionosphere H)
  ! for the purposes of tracking plasma entry into the magnetosphere in
  ! multifluid MHD.  It requires changes to the #MAGNETOSPHERE command
  ! in order to set the Iono species as the dominant outflow species at
  ! the inner boundary.  
  ! Based on ModEquationMhdIons and ModEquationMhdSwIono(Hansen, 2006)
  ! Welling, 2010.

  ! This equation module contains the standard MHD equations.
  character (len=*), parameter :: NameEquation='MHD with SW and Iono Hydrogen'

  integer, parameter :: nVar = 18

  integer, parameter :: nFluid    = 3
  integer, parameter :: IonFirst_ = 2
  integer, parameter :: IonLast_  = 3
  logical, parameter :: IsMhd     = .true.
  real               :: MassFluid_I(2:3) = (/ 1.0, 1.0 /)

  character (len=4), parameter :: NameFluid_I(nFluid)= (/ 'All', 'Sw ', 'Iono '/)

  ! Named indexes for State_VGB and other variables
  ! These indexes should go subsequently, from 1 to nVar+nFluid.
  ! The energies are handled as an extra variable, so that we can use
  ! both conservative and non-conservative scheme and switch between them.
  integer, parameter :: &
       Rho_       =  1,          &
       RhoUx_     =  2, Ux_ = 2, &
       RhoUy_     =  3, Uy_ = 3, &
       RhoUz_     =  4, Uz_ = 4, &
       Bx_        =  5,     &
       By_        =  6,     &
       Bz_        =  7,     &
       p_         =  8,     &
       SwRho_     =  9,     &
       SwRhoUx_   = 10,     &
       SwRhoUy_   = 11,     &
       SwRhoUz_   = 12,     &
       SwP_       = 13,     &
       IonoRho_     = 14,   &
       IonoRhoUx_   = 15,   &
       IonoRhoUy_   = 16,   &
       IonoRhoUz_   = 17,   &
       IonoP_       = 18,   &
       Energy_    = nVar+1, &
       SwEnergy_  = nVar+2, &
       IonoEnergy_  = nVar+3

  ! This allows to calculate RhoUx_ as RhoU_+x_ and so on.
  integer, parameter :: U_ = Ux_ - 1, RhoU_ = RhoUx_-1, B_ = Bx_-1

  ! These arrays are useful for multifluid
  integer, parameter :: iRho_I(nFluid)   = (/Rho_,   SwRho_,   IonoRho_/)
  integer, parameter :: iRhoUx_I(nFluid) = (/RhoUx_, SwRhoUx_, IonoRhoUx_/)
  integer, parameter :: iRhoUy_I(nFluid) = (/RhoUy_, SwRhoUy_, IonoRhoUy_/)
  integer, parameter :: iRhoUz_I(nFluid) = (/RhoUz_, SwRhoUz_, IonoRhoUz_/)
  integer, parameter :: iP_I(nFluid)     = (/p_,     SwP_,     IonoP_/)

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
       1.0, & ! SwRho_
       0.0, & ! SwRhoUx_
       0.0, & ! SwRhoUy_
       0.0, & ! SwRhoUz_
       1.0, & ! SwP_
       1.0, & ! IonoRho_
       0.0, & ! IonoRhoUx_
       0.0, & ! IonoRhoUy_
       0.0, & ! IonoRhoUz_
       1.0, & ! IonoP_
       1.0, & ! Energy_
       1.0, & ! SwEnergy_
       1.0  /)! IonoEnergy_

  ! The names of the variables used in i/o
  character(len=7) :: NameVar_V(nVar+nFluid) = (/ &
       'Rho    ', & ! Rho_
       'Mx     ', & ! RhoUx_
       'My     ', & ! RhoUy_
       'Mz     ', & ! RhoUz_
       'Bx     ', & ! Bx_
       'By     ', & ! By_
       'Bz     ', & ! Bz_
       'P      ', & ! p_
       'SwRho  ', & ! SwRho_
       'SwMx   ', & ! SwRhoUx_
       'SwMy   ', & ! SwRhoUy_
       'SwMz   ', & ! SwRhoUz_
       'SwP    ', & ! SwP_
       'IonoRho', & ! IonoRho_
       'IonoMx ', & ! IonoRhoUx_
       'IonoMy ', & ! IonoRhoUy_
       'IonoMz ', & ! IonoRhoUz_
       'IonoP  ', & ! IonoP_
       'E      ', & ! Energy_
       'SwE    ', & ! SwEnergy_
       'IonoE  ' /) ! IonoEnergy_

  ! The space separated list of nVar conservative variables for plotting
  character(len=*), parameter :: NameConservativeVar = &
       'Rho Mx My Mz Bx By Bz E '// &
       'SwRho SwMx SwMy SwMz SwE '// &
       'IonoRho IonoMx IonoMy IonoMz IonoE'

  ! The space separated list of nVar primitive variables for plotting
  character(len=*), parameter :: NamePrimitiveVar = &
       'Rho Ux Uy Uz Bx By Bz P '// &
       'SwRho SwUx SwUy SwUz SwP '// &
       'IonoRho IonoUx IonoUy IonoUz IonoP'

  ! The space separated list of nVar primitive variables for TECplot output
  character(len=*), parameter :: NamePrimitiveVarTec = &
       '"`r", "U_x", "U_y", "U_z", "B_x", "B_y", "B_z", "p", ' // &
       '"`r^SwH^+", "U_x^SwH^+", "U_y^SwH^+", "U_z^SwH^+", "P^SwH^+", ' // &
       '"`r^IonoH^+", "U_x^IonoH^+", "U_y^IonoH^+", "U_z^IonoH^+", "P^IonoH^+"'

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
