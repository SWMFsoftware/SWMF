!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModVarIndexes

  use ModExtraVariables, &
       Redefine1 => nWave, &
       Redefine2 => WaveFirst_, &
       Redefine3 => WaveLast_, &
       Redefine4 => Erad_, &
       Redefine5 => Pe_, &
       Redefine6 => Ehot_

  implicit none

  save

  ! This equation module contains the standard MHD equations with wave energy
  character (len=*), parameter :: &
       NameEquation='Multi-fluid MHD + Alfven waves + electron pressure'

  ! loop variable for implied do-loop over spectrum
  integer, private :: iWave

  integer, parameter :: nFluid    = 3
  integer, parameter :: IonFirst_ = 2
  integer, parameter :: IonLast_  = 3
  logical, parameter :: IsMhd     = .true.
  real :: MassFluid_I(2:nFluid) = (/ 1.0, 4.0 /)

  character(len=4), parameter :: NameFluid_I(nFluid) = &
       (/ 'All ', 'Hp  ', 'He2p' /)

  ! Number of wave bins in spectrum
  integer, parameter :: nWave = 2

  integer, parameter :: nVar = 10 + (nFluid - 1)*5 + nWave

  ! Named indexes for State_VGB and other variables
  ! These indexes should go subsequently, from 1 to nVar+1.
  ! The energy is handled as an extra variable, so that we can use
  ! both conservative and non-conservative scheme and switch between them.
  integer, parameter :: &
       Rho_        = 1,                  &
       RhoUx_      = 2,                  &
       RhoUy_      = 3,                  &
       RhoUz_      = 4,                  &
       Bx_         = 5,                  &
       By_         = 6,                  &
       Bz_         = 7,                  &
       p_          = 8,                  &
       HpRho_      = 9,                  &
       HpRhoUx_    = 10,                 &
       HpRhoUy_    = 11,                 &
       HpRhoUz_    = 12,                 &
       HpP_        = 13,                 &
       He2pRho_    = 14,                 &
       He2pRhoUx_  = 15,                 &
       He2pRhoUy_  = 16,                 &
       He2pRhoUz_  = 17,                 &
       He2pP_      = 18,                 &
       Pe_         = 19,                 &
       Ehot_       = 20,                 &
       WaveFirst_  = 21,                 &
       WaveLast_   = WaveFirst_+nWave-1, &
       Energy_     = nVar+1,             &
       HpEnergy_   = nVar+2,             &
       He2pEnergy_ = nVar+3

  ! This is for backward compatibility with single group radiation
  integer, parameter :: Erad_ = WaveFirst_

  ! This allows to calculate RhoUx_ as rhoU_+x_ and so on.
  integer, parameter :: RhoU_ = RhoUx_-1, B_ = Bx_-1

  ! These arrays are useful for multifluid
  integer, parameter :: iRho_I(nFluid)   = (/Rho_,   HpRho_,   He2pRho_/)
  integer, parameter :: iRhoUx_I(nFluid) = (/RhoUx_, HpRhoUx_, He2pRhoUx_/)
  integer, parameter :: iRhoUy_I(nFluid) = (/RhoUy_, HpRhoUy_, He2pRhoUy_/)
  integer, parameter :: iRhoUz_I(nFluid) = (/RhoUz_, HpRhoUz_, He2pRhoUz_/)
  integer, parameter :: iP_I(nFluid)     = (/p_,     HpP_,     He2pP_/)

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
       1.0, & ! He2pRho_
       0.0, & ! He2pRhoUx_
       0.0, & ! He2pRhoUy_
       0.0, & ! He2pRhoUz_
       1.0, & ! He2pP_        
       1.0, & ! Pe_
       0.0, & ! Ehot_
       (1.0, iWave=WaveFirst_,WaveLast_), &
       1.0, & ! Energy_
       1.0, & ! HpEnergy_
       1.0 /) ! He2pEnergy_

  ! The names of the variables used in i/o
  character(len=7) :: NameVar_V(nVar+nFluid) = (/ &
       'Rho    ', & ! Rho_
       'Mx     ', & ! RhoUx_
       'My     ', & ! RhoUy_
       'Mz     ', & ! RhoUz_
       'Bx     ', & ! Bx_
       'By     ', & ! By_
       'Bz     ', & ! Bz_
       'p      ', & ! p_
       'HpRho  ', & ! HpRho_
       'HpMx   ', & ! HpRhoUx_
       'HpMy   ', & ! HpRhoUy_
       'HpMz   ', & ! HpRhoUz_
       'HpP    ', & ! HpP_             
       'He2pRho', & ! He2pRho_
       'He2pMx ', & ! He2pRhoUx_
       'He2pMy ', & ! He2pRhoUy_
       'He2pMz ', & ! He2pRhoUz_
       'He2pP  ', & ! He2pP_             
       'Pe     ', & ! Pe_
       'Ehot   ', & ! Ehot_
       ('I??    ', iWave=WaveFirst_,WaveLast_), &
       'E      ', & ! Energy_
       'HpE    ', & ! HpEnergy_
       'He2pE  ' /) ! He2pEnergy_

  ! The space separated list of nVar conservative variables for plotting
  character(len=*), parameter :: NameConservativeVar = &
       'rho mx my mz bx by bz E '// &
       'HpRho HpMx HpMy HpMz HpE '// &
       'He2pRho He2pMx He2pMy He2pMz He2pE '// &
       'pe ehot ew'

  ! The space separated list of nVar primitive variables for plotting
  character(len=*), parameter :: NamePrimitiveVar = &
       'rho ux uy uz bx by bz p '// &
       'HpRho HpUx HpUy HpUz HpP '// &
       'He2pRho He2pUx He2pUy He2pUz He2pP '// &
       'pe ehot I(02)'

  ! The space separated list of nVar primitive variables for TECplot output
  character(len=*), parameter :: NamePrimitiveVarTec = &
       '"`r", "U_x", "U_y", "U_z", "B_x", "B_y", "B_z", "p", '// &
       '"`r^H^+", "U_x^H^+", "U_y^H^+", "U_z^H^+", "P^H^+", '// &
       '"`r^He2^+", "U_x^He2^+", "U_y^He2^+", "U_z^He2^+", "P^He2^+", '// &
       '"p_e", "Ehot", "Ew"'

  ! Names of the user units for IDL and TECPlot output
  character(len=20) :: &
       NameUnitUserIdl_V(nVar+nFluid) = '', NameUnitUserTec_V(nVar+nFluid) = ''

  ! The user defined units for the variables
  real :: UnitUser_V(nVar+nFluid) = 1.0

  ! Primitive variable names
  integer, parameter :: U_ = RhoU_, Ux_ = RhoUx_, Uy_ = RhoUy_, Uz_ = RhoUz_

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


