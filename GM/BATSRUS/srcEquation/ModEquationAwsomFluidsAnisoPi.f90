!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModVarIndexes

  use ModExtraVariables, &
       Redefine1 => nWave, &
       Redefine2 => WaveFirst_, &
       Redefine3 => WaveLast_, &
       Redefine4 => Erad_, &
       Redefine5 => Pe_, &
       Redefine6 => Ehot_, &
       Redefine7 => iPparIon_I, &
       Redefine8 => Ppar_

  implicit none

  save

  ! This equation module contains the standard MHD equations with wave energy
  character (len=*), parameter :: NameEquation = &
       'Multi-fluid MHD + Alfven waves + isotropic Pe + anisotropic Pi'

  ! loop variable for implied do-loop over spectrum
  integer, private :: iWave

  integer, parameter :: nFluid    = 2
  integer, parameter :: IonFirst_ = 1
  integer, parameter :: IonLast_  = 2
  logical, parameter :: IsMhd     = .false.
  real :: MassFluid_I(nFluid) = (/ 1.0, 4.0 /)

  character(len=4), parameter :: NameFluid_I(nFluid) = &
       (/ 'Hp  ', 'He2p' /)

  ! Number of wave bins in spectrum
  integer, parameter :: nWave = 2

  integer, parameter :: nVar = nFluid*6 + 5 + nWave

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
       Pe_         = 8,                  &
       Ehot_       = 9,                  &
       WaveFirst_  = 10,                 &
       WaveLast_   = WaveFirst_+nWave-1, &
       Ppar_       = WaveLast_+1,        &
       p_          = WaveLast_+2,        &
       He2pRho_    = p_+1,               &
       He2pRhoUx_  = p_+2,               &
       He2pRhoUy_  = p_+3,               &
       He2pRhoUz_  = p_+4,               &
       He2pPpar_   = p_+5,               &
       He2pP_      = p_+6,               &
       Energy_     = nVar+1,             &
       He2pEnergy_ = nVar+2

  ! This is for backward compatibility with single group radiation
  integer, parameter :: Erad_ = WaveFirst_

  ! This allows to calculate RhoUx_ as rhoU_+x_ and so on.
  integer, parameter :: RhoU_ = RhoUx_-1, B_ = Bx_-1

  ! These arrays are useful for multifluid
  integer, parameter :: iRho_I(nFluid)   = (/Rho_,   He2pRho_/)
  integer, parameter :: iRhoUx_I(nFluid) = (/RhoUx_, He2pRhoUx_/)
  integer, parameter :: iRhoUy_I(nFluid) = (/RhoUy_, He2pRhoUy_/)
  integer, parameter :: iRhoUz_I(nFluid) = (/RhoUz_, He2pRhoUz_/)
  integer, parameter :: iP_I(nFluid)     = (/p_,     He2pP_/)

  integer, parameter :: iPparIon_I(IonFirst_:IonLast_) = (/Ppar_, He2pPpar_/)

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
       0.0, & ! Ehot_
       (1.0, iWave=WaveFirst_,WaveLast_), &
       1.0, & ! Ppar_
       1.0, & ! p_
       1.0, & ! He2pRho_
       0.0, & ! He2pRhoUx_
       0.0, & ! He2pRhoUy_
       0.0, & ! He2pRhoUz_
       1.0, & ! He2pPpar_
       1.0, & ! He2pP_
       1.0, & ! Energy_
       1.0 /) ! He2pEnergy_

  ! The names of the variables used in i/o
  character(len=8) :: NameVar_V(nVar+nFluid) = (/ &
       'Rho     ', & ! Rho_
       'Mx      ', & ! RhoUx_
       'My      ', & ! RhoUy_
       'Mz      ', & ! RhoUz_
       'Bx      ', & ! Bx_
       'By      ', & ! By_
       'Bz      ', & ! Bz_
       'Pe      ', & ! Pe_
       'Ehot    ', & ! Ehot_
       ('I??     ', iWave=WaveFirst_,WaveLast_), &
       'Ppar    ', & ! Ppar_
       'p       ', & ! p_
       'He2pRho ', & ! He2pRho_
       'He2pMx  ', & ! He2pRhoUx_
       'He2pMy  ', & ! He2pRhoUy_
       'He2pMz  ', & ! He2pRhoUz_
       'He2pPpar', & ! He2pPpar_
       'He2pP   ', & ! He2pP_             
       'E       ', & ! Energy_
       'He2pE   ' /) ! He2pEnergy_


  ! Primitive variable names
  integer, parameter :: U_ = RhoU_, Ux_ = RhoUx_, Uy_ = RhoUy_, Uz_ = RhoUz_

  ! There are no extra scalars
  integer, parameter :: ScalarFirst_ = 2, ScalarLast_ = 1


end module ModVarIndexes


