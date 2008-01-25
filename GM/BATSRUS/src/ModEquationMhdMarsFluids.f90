module ModVarIndexes

  implicit none

  save

  ! This equation module contains the standard MHD equations.
  character (len=*), parameter :: NameEquation='Multi-fluid MHD for Mars'

  integer, parameter :: nVar = 23

  integer, parameter :: nFluid = 4
  integer, parameter :: nIonFluid = 4
  logical, parameter :: UseMultiIon = .true.
  real               :: MassFluid_I(nFluid) = (/ 1.0, 32.0, 16.0, 44.0 /)

  character (len=6), parameter :: NameFluid_I(nFluid) = &
       (/ 'All  ', 'O2p ', 'Op  ', 'CO2p' /)
  character (len=7), parameter :: TypeFluid_I(nFluid) = &
       (/ 'ion ', 'ions', 'ions', 'ions' /)

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
       O2pRho_     =  9, &
       O2pRhoUx_   = 10, &
       O2pRhoUy_   = 11, &
       O2pRhoUz_   = 12, &
       O2pP_       = 13, &
       OpRho_     =  14, &
       OpRhoUx_   = 15, &
       OpRhoUy_   = 16, &
       OpRhoUz_   = 17, &
       OpP_       = 18, &
       CO2pRho_     =  19, &
       CO2pRhoUx_   = 20, &
       CO2pRhoUy_   = 21, &
       CO2pRhoUz_   = 22, &
       CO2pP_       = 23, &
       Energy_    = nVar+1, &
       O2pEnergy_ = nVar+2,&
       OpEnergy_ = nVar+3,&
       CO2pEnergy_ = nVar+4


  ! This allows to calculate RhoUx_ as RhoU_+x_ and so on.
  integer, parameter :: U_ = Ux_ - 1, RhoU_ = RhoUx_-1, B_ = Bx_-1

  ! These arrays are useful for multifluid
  integer, parameter :: &
       iRho_I(nFluid)   = (/Rho_,   O2pRho_, OpRho_,  CO2pRho_ /), &
       iRhoUx_I(nFluid) = (/RhoUx_, O2pRhoUx_,  OpRhoUx_, CO2pRhoUx_ /), &
       iRhoUy_I(nFluid) = (/RhoUy_, O2pRhoUy_,  OpRhoUy_,  CO2pRhoUy_ /), &
       iRhoUz_I(nFluid) = (/RhoUz_, O2pRhoUz_, OpRhoUz_, CO2pRhoUz_ /), &
       iP_I(nFluid)     = (/p_,     O2pP_, OpP_ , CO2pP_ /)

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
       1.0, & ! O2pRho_
       0.0, & ! O2pRhoUx_
       0.0, & ! O2pRhoUy_
       0.0, & ! O2pRhoUz_
       1.0, & ! O2pP_
       1.0, & ! OpRho_
       0.0, & ! OpRhoUx_
       0.0, & ! OpRhoUy_
       0.0, & ! OpRhoUz_
       1.0, & ! OpP_
       1.0, & ! CO2pRho_
       0.0, & ! CO2pRhoUx_
       0.0, & ! CO2pRhoUy_
       0.0, & ! CO2pRhoUz_
       1.0, & ! CO2pP_
       1.0, & ! Energy_
       1.0, & ! O2pEnergy_
       1.0, & ! OpEnergy_
       1.0 /) ! CO2pEnergy_

  ! The names of the variables used in i/o
  character(len=*), parameter :: NameVar_V(nVar+nFluid) = (/ &
       'Rho    ', & ! Rho_
       'Mx     ', & ! RhoUx_
       'My     ', & ! RhoUy_
       'Mz     ', & ! RhoUz_
       'Bx     ', & ! Bx_
       'By     ', & ! By_
       'Bz     ', & ! Bz_
       'P      ', & ! p_
       'O2pRho ', & ! O2pRho_
       'O2pMx  ', & ! O2pRhoUx_
       'O2pMy  ', & ! O2pRhoUy_
       'O2pMz  ', & ! O2pRhoUz_
       'O2pP   ', & ! O2pP_
       'OpRho  ', & ! OpRho_
       'OpMx   ', & ! OpRhoUx_
       'OpMy   ', & ! OpRhoUy_
       'OpMz   ', & ! OpRhoUz_
       'OpP    ', & ! OpP_
       'CO2pRho', & ! CO2pRho_
       'CO2pMx ', & ! CO2pRhoUx_
       'CO2pMy ', & ! CO2pRhoUy_
       'CO2pMz ', & ! CO2pRhoUz_
       'CO2pP  ', & ! CO2pP_
       'HpE    ', & ! Energy_
       'O2pE   ', & ! O2pEnergy_
       'OpE    ', & ! OpEnergy_
       'CO2pE  ' /) ! CO2pEnergy_

  ! The space separated list of nVar conservative variables for plotting
  character(len=*), parameter :: NameConservativeVar = &
       'Rho Mx My Mz Bx By Bz E O2pRho O2pMx O2pMy O2pMz O2pE OpRho OpMx OpMy OpMz CO2pE CO2pRho CO2pMx CO2pMy CO2pMz CO2pE '

  ! The space separated list of nVar primitive variables for plotting
  character(len=*), parameter :: NamePrimitiveVar = &
       'Rho Ux Uy Uz Bx By Bz P O2pRho O2pUx O2pUy O2pUz O2pP OpRho OpUx OpUy OpUz OpP CO2pRho CO2pUx CO2pUy CO2pUz CO2pP '

  ! The space separated list of nVar primitive variables for TECplot output
  character(len=*), parameter :: NamePrimitiveVarTec = &
       '"`r^H^+", "U_x^H^+", "U_y^H^+", "U_z^H^+", "B_x", "B_y", "B_z", ' // &
       '"p", "`r^O2^+", "U_x^O2^+", "U_y^O2^+", "U_z^O2^+", "P^O2^+"' // &
       '"`r^O^+", "U_x^O^+", "U_y^O^+", "U_z^O^+", "P^O^+"'// &
       '"`r^CO2^+", "U_x^CO2^+", "U_y^CO2^+", "U_z^CO2^+", "P^CO2^+"'

  ! Names of the user units for IDL and TECPlot output
  character(len=20) :: &
       NameUnitUserIdl_V(nVar+nFluid) = '', NameUnitUserTec_V(nVar+nFluid) = ''

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
