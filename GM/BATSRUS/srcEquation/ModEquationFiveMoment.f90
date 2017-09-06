!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModVarIndexes

  use ModExtraVariables, &
       Redefine1 => iPparIon_I, &
       Redefine2 => Ex_,  Redefine3 => Ey_, Redefine4 => Ez_, &
       Redefine5 => HypE_,Redefine6 => nElectronFluid

  implicit none

  save

  ! This equation module contains the standard MHD equations.
  character (len=*), parameter :: NameEquation='Five Moment Closure'

  integer, parameter :: nVar = 17

  ! There are two ion fluids but no total ion fluid
  integer, parameter :: nFluid         = 2
  integer, parameter :: nElectronFluid = 1
  integer, parameter :: IonFirst_      = 1
  integer, parameter :: IonLast_       = 2
  logical, parameter :: IsMhd          = .false.
  real               :: MassFluid_I(nFluid) = (/ 1.0, 1.0/1836.0 /)

  character (len=3), parameter :: NameFluid_I(nFluid) = (/ 'Ion', 'El ' /)

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
       Ex_        =  8, &
       Ey_        =  9, &
       Ez_        = 10, &
       HypE_      = 11, &
       p_         = 12, &
       eRho_      = 13, &
       eRhoUx_    = 14, &
       eRhoUy_    = 15, &
       eRhoUz_    = 16, &
       eP_        = 17, &
       Energy_    = nVar+1, &
       eEnergy_   = nVar+2

  ! This allows to calculate RhoUx_ as RhoU_+x_ and so on.
  integer, parameter :: U_ = Ux_ - 1, RhoU_ = RhoUx_-1, B_ = Bx_-1

  ! These arrays are useful for multifluid
  integer, parameter :: iRho_I(nFluid)   = (/Rho_,   eRho_/)
  integer, parameter :: iRhoUx_I(nFluid) = (/RhoUx_, eRhoUx_/)
  integer, parameter :: iRhoUy_I(nFluid) = (/RhoUy_, eRhoUy_/)
  integer, parameter :: iRhoUz_I(nFluid) = (/RhoUz_, eRhoUz_/)
  integer, parameter :: iP_I(nFluid)     = (/p_,     eP_/)

  integer, parameter :: iPparIon_I(IonFirst_:IonLast_) = (/1,2/)

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
       0.0, & ! Ex_
       0.0, & ! Ey_
       0.0, & ! Ez_
       0.0, & ! HypE_
       1.0, & ! p_
       1.0, & ! eRho_
       0.0, & ! eRhoUx_
       0.0, & ! eRhoUy_
       0.0, & ! eRhoUz_
       1.0, & ! eP_
       1.0, & ! Energy_
       1.0 /) ! eEnergy_

  ! The names of the variables used in i/o
  character(len=5) :: NameVar_V(nVar+nFluid) = (/ &
       'Rho  ', & ! Rho_
       'Mx   ', & ! RhoUx_
       'My   ', & ! RhoUy_
       'Mz   ', & ! RhoUz_
       'Bx   ', & ! Bx_
       'By   ', & ! By_
       'Bz   ', & ! Bz_
       'Ex   ', & ! Ex_
       'Ey   ', & ! Ey_
       'Ez   ', & ! Ez_
       'HypE ', & ! HypE_
       'P    ', & ! p_
       'ElRho', & ! eRho_
       'ElMx ', & ! eRhoUx_
       'ElMy ', & ! eRhoUy_
       'ElMz ', & ! eRhoUz_
       'ElP  ', & ! eP_
       'E    ', & ! Energy_
       'ElE  ' /) ! eEnergy_

  ! The space separated list of nVar conservative variables for plotting
  character(len=*), parameter :: NameConservativeVar = &
       'Rho Mx My Mz Bx By Bz Ex Ey Ez HypE E ElRho ElMx ElMy ElMz ElE'

  ! The space separated list of nVar primitive variables for plotting
  character(len=*), parameter :: NamePrimitiveVar = &
       'Rho Ux Uy Uz Bx By Bz Ex Ey Ez HypE P ElRho ElUx ElUy ElUz ElP'

  ! The space separated list of nVar primitive variables for TECplot output
  character(len=*), parameter :: NamePrimitiveVarTec = &
       '"`r^H^+", "U_x^H^+", "U_y^H^+", "U_z^H^+", "B_x", "B_y", "B_z", ' // &
       '"E_x", "E_y", "E_z", "HypE", "p", ' // &
       '"`r^e^-", "U_x^e^-", "U_y^e^-", "U_z^e^-", "P^e^-"'

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
