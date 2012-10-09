module ModVarIndexes

  use ModExtraVariables, Redefine => Pe_

  implicit none

  save

  ! This equation file declares three ion fluids: solar wind H+, ionospheric
  ! H+, and ionospheric O+ along with ion electron pressure.  This allows for
  ! thorough investigations of each populations entry and heating mechanisms
  ! within the terrestrial magnetosphere.
  character (len=*), parameter :: NameEquation = &
       'MHD with SW and Iono H+, Iono O+, and electron pressure'

  integer, parameter :: nVar = 24

  integer, parameter :: nFluid    = 4
  integer, parameter :: IonFirst_ = 2
  integer, parameter :: IonLast_  = 4
  logical, parameter :: IsMhd     = .true.
  real               :: MassFluid_I(2:4) = (/ 1.0, 1.0, 16.0 /)

  character (len=6), parameter :: NameFluid_I(nFluid)= &
       (/ 'All   ', 'HpSw  ', 'HpIono', 'OpIono'/)

  ! Named indexes for State_VGB and other variables
  ! These indexes should go subsequently, from 1 to nVar+nFluid.
  ! The energies are handled as an extra variable, so that we can use
  ! both conservative and non-conservative scheme and switch between them.
  integer, parameter ::   &
       Rho_         =  1,          &
       RhoUx_       =  2, Ux_ = 2, &
       RhoUy_       =  3, Uy_ = 3, &
       RhoUz_       =  4, Uz_ = 4, &
       Bx_          =  5, &
       By_          =  6, &
       Bz_          =  7, &
       Pe_          =  8, &
       p_           =  9, & 
       HpSwRho_     = 10, &
       HpSwRhoUx_   = 11, &
       HpSwRhoUy_   = 12, &
       HpSwRhoUz_   = 13, &
       HpSwP_       = 14, &
       HpIonoRho_   = 15, &
       HpIonoRhoUx_ = 16, &
       HpIonoRhoUy_ = 17, &
       HpIonoRhoUz_ = 18, &
       HpIonoP_     = 19, &
       OpIonoRho_   = 20, &
       OpIonoRhoUx_ = 21, &
       OpIonoRhoUy_ = 22, &
       OpIonoRhoUz_ = 23, &
       OpIonoP_     = 24, &
       Energy_        = nVar+1, &
       HpSwEnergy_    = nVar+2, &
       HpIonoEnergy_  = nVar+3, &
       OpIonoEnergy_  = nVar+4

  ! This allows to calculate RhoUx_ as RhoU_+x_ and so on.
  integer, parameter :: U_ = Ux_ - 1, RhoU_ = RhoUx_-1, B_ = Bx_-1

  ! These arrays are useful for multifluid
  integer, parameter :: &
       iRho_I(nFluid)   = (/Rho_,   HpSwRho_,   HpIonoRho_,   OpIonoRho_/), &
       iRhoUx_I(nFluid) = (/RhoUx_, HpSwRhoUx_, HpIonoRhoUx_, OpIonoRhoUx_/), &
       iRhoUy_I(nFluid) = (/RhoUy_, HpSwRhoUy_, HpIonoRhoUy_, OpIonoRhoUy_/), &
       iRhoUz_I(nFluid) = (/RhoUz_, HpSwRhoUz_, HpIonoRhoUz_, OpIonoRhoUz_/), &
       iP_I(nFluid)     = (/p_,     HpSwP_,     HpIonoP_,     OpIonoP_/)

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
       1.0, & ! HpSwRho_
       0.0, & ! HpSwRhoUx_
       0.0, & ! HpSwRhoUy_
       0.0, & ! HpSwRhoUz_
       1.0, & ! HpSwP_
       1.0, & ! HpIonoRho_
       0.0, & ! HpIonoRhoUx_
       0.0, & ! HpIonoRhoUy_
       0.0, & ! HpIonoRhoUz_
       1.0, & ! HpIonoP_
       1.0, & ! OpIonoRho_
       0.0, & ! OpIonoRhoUx_
       0.0, & ! OpIonoRhoUy_
       0.0, & ! OpIonoRhoUz_
       1.0, & ! OpIonoP_
       1.0, & ! Energy_
       1.0, & ! HpSwEnergy_
       1.0, & ! HpIonoEnergy_
       1.0  /)! OpIonoEnergy_

  ! The names of the variables used in i/o
  character(len=9) :: NameVar_V(nVar+nFluid) = (/ &
       'Rho      ', & ! Rho_
       'Mx       ', & ! RhoUx_
       'My       ', & ! RhoUy_
       'Mz       ', & ! RhoUz_
       'Bx       ', & ! Bx_
       'By       ', & ! By_
       'Bz       ', & ! Bz_
       'Pe       ', & ! Pe_
       'P        ', & ! p_
       'HpSwRho  ', & ! HpRho_
       'HpSwMx   ', & ! HpRhoUx_
       'HpSwMy   ', & ! HpRhoUy_
       'HpSwMz   ', & ! HpRhoUz_
       'HpSwP    ', & ! HpP_
       'HpIonoRho', & ! HpRho_
       'HpIonoMx ', & ! HpRhoUx_
       'HpIonoMy ', & ! HpRhoUy_
       'HpIonoMz ', & ! HpRhoUz_
       'HpIonoP  ', & ! HpP_
       'OpIonoRho', & ! OpIonoRho_
       'OpIonoMx ', & ! OpIonoRhoUx_
       'OpIonoMy ', & ! OpIonoRhoUy_
       'OpIonoMz ', & ! OpIonoRhoUz_
       'OpIonoP  ', & ! OpIonoP_
       'E        ', & ! Energy_
       'HpSwE    ', & ! HpEnergy_
       'HpIonoE  ', & ! HpIonoEnergy_
       'OpIonoE  ' /) ! OpEnergy_

  ! The space separated list of nVar conservative variables for plotting
  character(len=*), parameter :: NameConservativeVar = &
       'Rho Mx My Mz Bx By Bz Pe E '// &
       'HpSwRho HpSwMx HpSwMy HpSwMz HpSwE '// &
       'HpIonoRho HpIonoMx HpIonoMy HpIonoMz HpIonoE '// &
       'OpIonoRho OpIonoMx OpIonoMy OpIonoMz OpIonoE'

  ! The space separated list of nVar primitive variables for plotting
  character(len=*), parameter :: NamePrimitiveVar = &
       'Rho Ux Uy Uz Bx By Bz Pe P '// &
       'HpSwRho HpSwUx HpSwUy HpSwUz HpSwP '// &
       'HpIonoRho HpIonoUx HpIonoUy HpIonoUz HpIonoP '// &
       'OpIonoRho OpIonoUx OpIonoUy OpIonoUz OpIonoP'

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
