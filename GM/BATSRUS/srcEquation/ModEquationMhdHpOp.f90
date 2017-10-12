!  Copyright (C) 2002 Regents of the University of Michigan
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModVarIndexes

  use ModSingleFluid
  use ModExtraVariables, &
       Redefine1 => SpeciesFirst_, &
       Redefine2 => SpeciesLast_,  &
       Redefine3 => MassSpecies_V

  implicit none

  save

  ! This equation module contains the standard MHD equations with
  ! two species for Earth.  1 - solar wind protons, 2 - ionospheric plasma
  character (len=*), parameter :: NameEquation= 'MHD 2 Species (Hp+Op), Toth, 2017'

  ! Number of variables without energy:
  integer, parameter :: nVar = 10

  ! Named indexes for State_VGB and other variables
  ! These indexes should go subsequently, from 1 to nVar+1.
  ! The energy is handled as an extra variable, so that we can use
  ! both conservative and non-conservative scheme and switch between them.
  integer, parameter :: &
       Rho_   = 1,    &
       RhoHp_ = 2,    &
       RhoOp_ = 3,    &
       RhoUx_ = 4,    &
       RhoUy_ = 5,    &
       RhoUz_ = 6,    &
       Bx_    = 7,    &
       By_    = 8,    &
       Bz_    = 9,    &
       p_     = nVar, &
       Energy_= nVar+1

  ! This allows to calculate RhoUx_ as rhoU_+x_ and so on.
  integer, parameter :: RhoU_ = RhoUx_-1, B_ = Bx_-1

  ! These arrays are useful for multifluid
  integer, parameter :: iRho_I(nFluid)   = (/Rho_/)
  integer, parameter :: iRhoUx_I(nFluid) = (/RhoUx_/)
  integer, parameter :: iRhoUy_I(nFluid) = (/RhoUy_/)
  integer, parameter :: iRhoUz_I(nFluid) = (/RhoUz_/)
  integer, parameter :: iP_I(nFluid)     = (/p_/)

  ! The default values for the state variables:
  ! Variables which are physically positive should be set to 1,
  ! variables that can be positive or negative should be set to 0:
  real, parameter :: DefaultState_V(nVar+1) = (/ & 
       1.0, & ! Rho_
       1.0, & ! RhoHp_
       1.0, & ! RhoOp_
       0.0, & ! RhoUx_
       0.0, & ! RhoUy_
       0.0, & ! RhoUz_
       0.0, & ! Bx_
       0.0, & ! By_
       0.0, & ! Bz_
       1.0, & ! p_
       1.0 /) ! Energy_

  ! The names of the variables used in i/o
  character(len=5) :: NameVar_V(nVar+1) = (/ &
       'Rho  ', & ! Rho_
       'RhoHp', & ! RhoHp_
       'RhoOp', & ! RhoOp_
       'Mx   ', & ! RhoUx_
       'My   ', & ! RhoUy_
       'Mz   ', & ! RhoUz_
       'Bx   ', & ! Bx_
       'By   ', & ! By_
       'Bz   ', & ! Bz_
       'p    ', & ! p_
       'e    ' /) ! Energy_


  ! Primitive variable names
  integer, parameter :: U_ = RhoU_, Ux_ = RhoUx_, Uy_ = RhoUy_, Uz_ = RhoUz_

  ! There are two extra scalars
  integer, parameter :: ScalarFirst_ = RhoHp_, ScalarLast_ = RhoOp_

  ! Species
  integer, parameter :: SpeciesFirst_   = ScalarFirst_
  integer, parameter :: SpeciesLast_    = ScalarLast_

  ! Molecular mass of solarwind and ionosphere species in AMU:
  real :: MassSpecies_V(SpeciesFirst_:SpeciesLast_) = (/ 1.0, 16.0 /)

end module ModVarIndexes
