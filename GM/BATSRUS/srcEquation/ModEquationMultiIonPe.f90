!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModVarIndexes

  use ModExtraVariables, &
       Redefine => iPparIon_I, Redefine2 => Pe_

  implicit none

  save

  ! This equation module contains the standard MHD equations.
  character (len=*), parameter :: NameEquation='MultiIon: H+, O+ and Pe'

  integer, parameter :: nVar = 14

  ! There are two ion fluids but no total ion fluid
  integer, parameter :: nFluid    = 2
  integer, parameter :: IonFirst_ = 1
  integer, parameter :: IonLast_  = 2
  logical, parameter :: IsMhd     = .false.
  real               :: MassFluid_I(nFluid) = (/ 1.0, 16.0 /)

  character (len=2), parameter :: NameFluid_I(nFluid) = (/ 'Hp', 'Op' /)

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
       Pe_        =  8, &
       p_         =  9, &
       OpRho_     = 10, &
       OpRhoUx_   = 11, &
       OpRhoUy_   = 12, &
       OpRhoUz_   = 13, &
       OpP_       = 14, &
       Energy_    = nVar+1, &
       OpEnergy_  = nVar+2

  ! This allows to calculate RhoUx_ as RhoU_+x_ and so on.
  integer, parameter :: U_ = Ux_ - 1, RhoU_ = RhoUx_-1, B_ = Bx_-1

  ! These arrays are useful for multifluid
  integer, parameter :: iRho_I(nFluid)   = (/Rho_,   OpRho_/)
  integer, parameter :: iRhoUx_I(nFluid) = (/RhoUx_, OpRhoUx_/)
  integer, parameter :: iRhoUy_I(nFluid) = (/RhoUy_, OpRhoUy_/)
  integer, parameter :: iRhoUz_I(nFluid) = (/RhoUz_, OpRhoUz_/)
  integer, parameter :: iP_I(nFluid)     = (/p_,     OpP_/)

  integer, parameter :: iPparIon_I(IonFirst_:IonLast_) = 1

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
       1.0, & ! Pe_
       1.0, & ! p_
       1.0, & ! OpRho_
       0.0, & ! OpRhoUx_
       0.0, & ! OpRhoUy_
       0.0, & ! OpRhoUz_
       1.0, & ! OpP_
       1.0, & ! Energy_
       1.0 /) ! OpEnergy_

  ! The names of the variables used in i/o
  character(len=5) :: NameVar_V(nVar+nFluid) = (/ &
       'HpRho', & ! Rho_
       'HpMx ', & ! RhoUx_
       'HpMy ', & ! RhoUy_
       'HpMz ', & ! RhoUz_
       'Bx   ', & ! Bx_
       'By   ', & ! By_
       'Bz   ', & ! Bz_
       'Pe   ', & ! Pe_
       'HpP  ', & ! p_
       'OpRho', & ! OpRho_
       'OpMx ', & ! OpRhoUx_
       'OpMy ', & ! OpRhoUy_
       'OpMz ', & ! OpRhoUz_
       'OpP  ', & ! OpP_
       'HpE  ', & ! Energy_
       'OpE  ' /) ! OpEnergy_

  ! The space separated list of nVar conservative variables for plotting
  character(len=*), parameter :: NameConservativeVar = &
       'HpRho HpMx HpMy HpMz Bx By Bz Pe HpE OpRho OpMx OpMy OpMz OpE'

  ! The space separated list of nVar primitive variables for plotting
  character(len=*), parameter :: NamePrimitiveVar = &
       'HpRho HpUx HpUy HpUz Bx By Bz Pe HpP OpRho OpUx OpUy OpUz OpP'

  ! The space separated list of nVar primitive variables for TECplot output
  character(len=*), parameter :: NamePrimitiveVarTec = &
       '"`r^H^+", "U_x^H^+", "U_y^H^+", "U_z^H^+", "B_x", "B_y", "B_z", ' // &
       '"p_e", "p", "`r^O^+", "U_x^O^+", "U_y^O^+", "U_z^O^+", "P^O^+"'

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

end module ModVarIndexes
