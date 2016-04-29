!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!#NOTPUBLIC  email:rubinmar@umich.edu  expires:12/31/2099
module ModVarIndexes

  use ModExtraVariables, &
       Redefine1 => Pe_, &
       Redefine2 => iPparIon_I

  implicit none

  save

  ! This equation module contains the standard MHD equations 
  ! with electron pressure
  character (len=*), parameter :: NameEquation= &
       '3-fluid + Pe MHD for Comets'

  integer, parameter :: nVar = 24

  integer, parameter :: nFluid    = 4
  integer, parameter :: IonFirst_ = 2        ! First individual ion fluid
  integer, parameter :: IonLast_  = 4        ! Last individual ion fluid
  logical, parameter :: IsMhd     = .true.   ! First total ion fluid obeys MHD
  real               :: MassFluid_I(2:nFluid) = (/ 1.0, 1.0, 17.0 /)

  ! Fluids: total fluid, solar wind protons, cometary protons, water ions
  character (len=4), parameter :: NameFluid_I(nFluid) = &
       (/ 'All ', 'Sw  ', 'Hp  ', 'H2Op' /)

  ! Named indexes for State_VGB and other variables
  ! These indexes should go subsequently, from 1 to nVar+nFluid.
  ! The energy is handled as an extra variable, so that we can use
  ! both conservative and non-conservative scheme and switch between them.
  integer, parameter :: &
       Rho_       =  1, &
       RhoUx_     =  2, Ux_ = 2, &
       RhoUy_     =  3, Uy_ = 3, &
       RhoUz_     =  4, Uz_ = 4, &
       Bx_        =  5, &
       By_        =  6, &
       Bz_        =  7, &
       Pe_        =  8, &
       p_         =  9, &
       SwRho_     = 10, &
       SwRhoUx_   = 11, &
       SwRhoUy_   = 12, &
       SwRhoUz_   = 13, &
       SwP_       = 14, &
       HpRho_     = 15, &
       HpRhoUx_   = 16, &
       HpRhoUy_   = 17, &
       HpRhoUz_   = 18, &
       HpP_       = 19, &
       H2OpRho_   = 20, &
       H2OpRhoUx_ = 21, &
       H2OpRhoUy_ = 22, &
       H2OpRhoUz_ = 23, &
       H2OpP_     = 24, &
       Energy_    = nVar+1, &
       SwEnergy_ = nVar+2, &
       HpEnergy_  = nVar+3, &
       H2OpEnergy_= nVar+4

  ! This allows to calculate RhoUx_ as RhoU_+x_ and so on.
  integer, parameter :: U_ = Ux_ - 1, RhoU_ = RhoUx_-1, B_ = Bx_-1

  ! These arrays are useful for multifluid
  integer, parameter :: &
       iRho_I(nFluid)  =(/Rho_,   SwRho_,   HpRho_,   H2OpRho_ /) ,&
       iRhoUx_I(nFluid)=(/RhoUx_, SwRhoUx_, HpRhoUx_, H2OpRhoUx_ /),&
       iRhoUy_I(nFluid)=(/RhoUy_, SwRhoUy_, HpRhoUy_, H2OpRhoUy_ /),&
       iRhoUz_I(nFluid)=(/RhoUz_, SwRhoUz_, HpRhoUz_, H2OpRhoUz_ /),&
       iP_I(nFluid)    =(/p_,     SwP_,     HpP_,     H2OpP_ /)

  integer, parameter :: iPparIon_I(IonFirst_:IonLast_) = (/1,2,3/)

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
       1.0, & ! H2OpRho_
       0.0, & ! H2OpRhoUx_
       0.0, & ! H2OpRhoUy_
       0.0, & ! H2OpRhoUz_
       1.0, & ! H2OpP_
       1.0, & ! Energy_
       1.0, & ! SwEnergy_
       1.0, & ! HpEnergy_
       1.0 /) ! H2OpEnergy_

  ! The names of the variables used in i/o
  character(len=7) :: NameVar_V(nVar+nFluid) = (/ &
       'Rho    ', & ! Rho_
       'Mx     ', & ! RhoUx_
       'My     ', & ! RhoUy_
       'Mz     ', & ! RhoUz_
       'Bx     ', & ! Bx_
       'By     ', & ! By_
       'Bz     ', & ! Bz_
       'Pe     ', & ! Pe_
       'p      ', & ! p_
       'SwRho  ', & ! SwRho_
       'SwMx   ', & ! SwRhoUx_
       'SwMy   ', & ! SwRhoUy_
       'SwMz   ', & ! SwRhoUz_
       'SwP    ', & ! SwP_
       'HpRho  ', & ! HpRho_
       'HpMx   ', & ! HpRhoUx_
       'HpMy   ', & ! HpRhoUy_
       'HpMz   ', & ! HpRhoUz_
       'HpP    ', & ! HpP_
       'H2OpRho', & ! H2OpRho_
       'H2OpMx ', & ! H2OpRhoUx_
       'H2OpMy ', & ! H2OpRhoUy_
       'H2OpMz ', & ! H2OpRhoUz_
       'H2OpP  ', & ! H2OpP_
       'E      ', & ! Energy_
       'SwE    ', & ! SwEnergy_
       'HpE    ', & ! HpEnergy_
       'H2OpE  ' /) ! H2OpEnergy_

  ! The space separated list of nVar conservative variables for plotting
  character(len=*), parameter :: NameConservativeVar = &
       'Rho Mx My Mz Bx By Bz E '// &
       'SwRho SwMx SwMy SwMz SwE '// &
       'HpRho HpMx HpMy HpMz HpE '// &
       'H2OpRho H2OpMx H2OpMy H2OpMz H2OpE '

  ! The space separated list of nVar primitive variables for plotting
  character(len=*), parameter :: NamePrimitiveVar = &
       'Rho Ux Uy Uz Bx By Bz Pe P ' // &
       'SwRho SwUx SwUy SwUz SwP '// &
       'HpRho HpUx HpUy HpUz HpP '// &
       'H2OpRho H2OpUx H2OpUy H2OpUz H2OpP '

  ! The space separated list of nVar primitive variables for TECplot output
  character(len=*), parameter :: NamePrimitiveVarTec = &
       '"`r", "U_x", "U_y", "U_z", "B_x", "B_y", "B_z", "p_e", "p",' // &
       '"`r^SW^+", "U_x^SW^+", "U_y^SW^+", "U_z^SW^+", "P^SW^+"'// &
       '"`r^H^+", "U_x^H^+", "U_y^H^+", "U_z^H^+", "P^H^+"'// &
       '"`r^H2O^+", "U_x^H2O^+", "U_y^H2O^+", "U_z^H2O^+", "P^H2O^+"'

  ! Names of the user units for IDL and TECPlot output
  character(len=20) :: &
       NameUnitUserIdl_V(nVar+nFluid) = '', NameUnitUserTec_V(nVar+nFluid) = ''

  ! The user defined units for the variables
  real :: UnitUser_V(nVar+nFluid) = 1.0

  ! There are no extra scalars (Pe has its own flux)
  integer, parameter :: ScalarFirst_ = 2, ScalarLast_ = 1

  ! There are no multi-species
  logical, parameter :: UseMultiSpecies = .false.

  ! Declare the following variables to satisfy the compiler
  integer, parameter :: SpeciesFirst_ = 1, SpeciesLast_ = 1
  real               :: MassSpecies_V(SpeciesFirst_:SpeciesLast_)

end module ModVarIndexes

