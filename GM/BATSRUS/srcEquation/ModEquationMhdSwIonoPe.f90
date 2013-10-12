!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModVarIndexes

  use ModExtraVariables, Redefine => Pe_

  implicit none

  save

  ! This equation file declares three ion fluids: solar wind H+, ionospheric
  ! H+, and ionospheric O+ along with ion electron pressure.  This allows for
  ! thorough investigations of each populations entry and heating mechanisms
  ! within the terrestrial magnetosphere.  Solar wind values default to 
  ! first fluid; user must specify inner boundary densities using 
  ! #MAGNETOSPHERE command.
  character (len=*), parameter :: NameEquation = &
       'MHD with SW and Iono H+, Iono O+, and electron pressure'

  integer, parameter :: nVar = 24

  integer, parameter :: nFluid    = 4
  integer, parameter :: IonFirst_ = 2
  integer, parameter :: IonLast_  = 4
  logical, parameter :: IsMhd     = .true.
  real               :: MassFluid_I(2:4) = (/ 1.0, 1.0, 16.0 /)

  character (len=3), parameter :: NameFluid_I(nFluid)= &
       (/ 'All', 'Sw ', 'Hp ', 'Op '/)

  ! Named indexes for State_VGB and other variables
  ! These indexes should go subsequently, from 1 to nVar+nFluid.
  ! The energies are handled as an extra variable, so that we can use
  ! both conservative and non-conservative scheme and switch between them.
  integer, parameter ::   &
       Rho_     =  1,          &
       RhoUx_   =  2, Ux_ = 2, &
       RhoUy_   =  3, Uy_ = 3, &
       RhoUz_   =  4, Uz_ = 4, &
       Bx_      =  5, &
       By_      =  6, &
       Bz_      =  7, &
       Pe_      =  8, &
       p_       =  9, & 
       SwRho_   = 10, &
       SwRhoUx_ = 11, &
       SwRhoUy_ = 12, &
       SwRhoUz_ = 13, &
       SwP_     = 14, &
       HpRho_   = 15, &
       HpRhoUx_ = 16, &
       HpRhoUy_ = 17, &
       HpRhoUz_ = 18, &
       HpP_     = 19, &
       OpRho_   = 20, &
       OpRhoUx_ = 21, &
       OpRhoUy_ = 22, &
       OpRhoUz_ = 23, &
       OpP_     = 24, &
       Energy_  = nVar+1, &
       SwEnergy_= nVar+2, &
       HpEnergy_= nVar+3, &
       OpEnergy_= nVar+4

  ! This allows to calculate RhoUx_ as RhoU_+x_ and so on.
  integer, parameter :: U_ = Ux_ - 1, RhoU_ = RhoUx_-1, B_ = Bx_-1

  ! These arrays are useful for multifluid
  integer, parameter :: &
       iRho_I(nFluid)   = (/Rho_,   SwRho_,   HpRho_,   OpRho_/), &
       iRhoUx_I(nFluid) = (/RhoUx_, SwRhoUx_, HpRhoUx_, OpRhoUx_/), &
       iRhoUy_I(nFluid) = (/RhoUy_, SwRhoUy_, HpRhoUy_, OpRhoUy_/), &
       iRhoUz_I(nFluid) = (/RhoUz_, SwRhoUz_, HpRhoUz_, OpRhoUz_/), &
       iP_I(nFluid)     = (/p_,     SwP_,     HpP_,     OpP_/)

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
       1.0, & ! SwRho_
       0.0, & ! SwRhoUx_
       0.0, & ! SwRhoUy_
       0.0, & ! SwRhoUz_
       1.0, & ! SwP_
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
       1.0, & ! SwEnergy_
       1.0, & ! HpEnergy_
       1.0  /)! OpEnergy_

  ! The names of the variables used in i/o
  character(len=5) :: NameVar_V(nVar+nFluid) = (/ &
       'Rho  ', & ! Rho_
       'Mx   ', & ! RhoUx_
       'My   ', & ! RhoUy_
       'Mz   ', & ! RhoUz_
       'Bx   ', & ! Bx_
       'By   ', & ! By_
       'Bz   ', & ! Bz_
       'Pe   ', & ! Pe_
       'P    ', & ! p_
       'SwRho', & ! HpRho_
       'SwMx ', & ! HpRhoUx_
       'SwMy ', & ! HpRhoUy_
       'SwMz ', & ! HpRhoUz_
       'SwP  ', & ! HpP_
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
       'SwE  ', & ! HpEnergy_
       'HpE  ', & ! HpEnergy_
       'OpE  ' /) ! OpEnergy_

  ! The space separated list of nVar conservative variables for plotting
  character(len=*), parameter :: NameConservativeVar = &
       'Rho Mx My Mz Bx By Bz Pe E '// &
       'SwRho SwMx SwMy SwMz SwE '// &
       'HpRho HpMx HpMy HpMz HpE '// &
       'OpRho OpMx OpMy OpMz OpE'

  ! The space separated list of nVar primitive variables for plotting
  character(len=*), parameter :: NamePrimitiveVar = &
       'Rho Ux Uy Uz Bx By Bz Pe P '// &
       'SwRho SwUx SwUy SwUz SwP '// &
       'HpRho HpUx HpUy HpUz HpP '// &
       'OpRho OpUx OpUy OpUz OpP'

  ! The space separated list of nVar primitive variables for TECplot output
  character(len=*), parameter :: NamePrimitiveVarTec = &
       '"`r", "U_x", "U_y", "U_z", "B_x", "B_y", "B_z", "p_e", "p", ' // &
       '"`r^Hsw^+", "U_x^Hsw^+", "U_y^Hsw^+", "U_z^Hsw^+", "P^Hsw^+", ' // &
       '"`r^Hiono^+", "U_x^Hiono^+", "U_y^Hiono^+", "U_z^Hiono^+", "P^Hiono^+", ' // &
       '"`r^Oiono^+", "U_x^Oiono^+", "U_y^Oiono^+", "U_z^Oiono^+", "P^Oiono^+"'

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
