module ModVarIndexes

  use ModExtraVariables

  implicit none

  save

  ! This equation module contains the standard MHD equations.
  character (len=*), parameter :: NameEquation='Multi-fluid MHD for Mars with Hpsw and Hpion separate'

  integer, parameter :: nVar = 33

  integer, parameter :: nFluid    = 6
  integer, parameter :: IonFirst_ = 2        ! First individual ion fluid
  integer, parameter :: IonLast_  = 6        ! Last individual ion fluid
  logical, parameter :: IsMhd     = .true.   ! First total ion fluid obeys MHD
  real               :: MassFluid_I(2:nFluid) = (/ 1.0, 32.0, 16.0, 44.0, 1.0 /)

  character (len=5), parameter:: NameFluid_I(nFluid) = &
       (/ 'All  ', 'Hpsw ', 'O2p  ', 'Op   ', 'CO2p ', 'Hp   ' /)

  ! Named indexes for State_VGB and other variables
  ! These indexes should go subsequently, from 1 to nVar+nFluid.
  ! The energies are handled as an extra variable, so that we can use
  ! both conservative and non-conservative scheme and switch between them.
  integer, parameter :: &
       Rho_        =  1,          &
       RhoUx_      =  2, Ux_ = 2, &
       RhoUy_      =  3, Uy_ = 3, &
       RhoUz_      =  4, Uz_ = 4, &
       Bx_         =  5, &
       By_         =  6, &
       Bz_         =  7, &
       p_          =  8, &
       HpswRho_    =  9, &
       HpswRhoUx_  = 10, &
       HpswRhoUy_  = 11, &
       HpswRhoUz_  = 12, &
       HpswP_      = 13, &
       O2pRho_     = 14, &
       O2pRhoUx_   = 15, &
       O2pRhoUy_   = 16, &
       O2pRhoUz_   = 17, &
       O2pP_       = 18, &
       OpRho_      = 19, &
       OpRhoUx_    = 20, &
       OpRhoUy_    = 21, &
       OpRhoUz_    = 22, &
       OpP_        = 23, &
       CO2pRho_    = 24, &
       CO2pRhoUx_  = 25, &
       CO2pRhoUy_  = 26, &
       CO2pRhoUz_  = 27, &
       CO2pP_      = 28, &
       HpRho_      = 29, &
       HpRhoUx_    = 30, &
       HpRhoUy_    = 31, &
       HpRhoUz_    = 32, &
       HpP_        = 33, &
       Energy_     = nVar+1, &
       HpswEnergy_ = nVar+2, &
       O2pEnergy_  = nVar+3,&
       OpEnergy_   = nVar+4,&
       CO2pEnergy_ = nVar+5,&
       HpEnergy_   = nVar+6

  
  ! This allows to calculate RhoUx_ as RhoU_+x_ and so on.
  integer, parameter :: U_ = Ux_ - 1, RhoU_ = RhoUx_-1, B_ = Bx_-1

  ! These arrays are useful for multifluid
  integer, parameter :: &
       iRho_I(nFluid)  =(/Rho_,   HpswRho_,   O2pRho_,   OpRho_,   CO2pRho_,   HpRho_ /) ,&
       iRhoUx_I(nFluid)=(/RhoUx_, HpswRhoUx_, O2pRhoUx_, OpRhoUx_, CO2pRhoUx_, HpRhoUx_/),&
       iRhoUy_I(nFluid)=(/RhoUy_, HpswRhoUy_, O2pRhoUy_, OpRhoUy_, CO2pRhoUy_, HpRhoUy_/),&
       iRhoUz_I(nFluid)=(/RhoUz_, HpswRhoUz_, O2pRhoUz_, OpRhoUz_, CO2pRhoUz_, HpRhoUz_/),&
       iP_I(nFluid)    =(/p_,     HpswP_,     O2pP_,     OpP_ ,    CO2pP_,     HpP_ /)

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
       1.0, & ! HpswRho_
       0.0, & ! HpswRhoUx_
       0.0, & ! HpswRhoUy_
       0.0, & ! HpswRhoUz_
       1.0, & ! HpswP_
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
       1.0, & ! HpRho_ 
       0.0, & ! HpUx_   
       0.0, & ! HpUy_   
       0.0, & ! HpUz_   
       1.0, & ! HpP_      
       1.0, & ! Energy_
       1.0, & ! HpswEnergy_
       1.0, & ! O2pEnergy_
       1.0, & ! OpEnergy_
       1.0, & ! CO2pEnergy_
       1.0/)  ! HpEnergy


  ! The names of the variables used in i/o
  character(len=8) :: NameVar_V(nVar+nFluid) = (/ &
       'Rho     ', & ! Rho_
       'Mx      ', & ! RhoUx_
       'My      ', & ! RhoUy_
       'Mz      ', & ! RhoUz_
       'Bx      ', & ! Bx_
       'By      ', & ! By_
       'Bz      ', & ! Bz_
       'P       ', & ! p_
       'HpswRho ', & ! HpswRho_
       'HpswMx  ', & ! HpswRhoUx_
       'HpswMy  ', & ! HpswRhoUy_
       'HpswMz  ', & ! HpswRhoUz_
       'HpswP   ', & ! HpswP_
       'O2pRho  ', & ! O2pRho_
       'O2pMx   ', & ! O2pRhoUx_
       'O2pMy   ', & ! O2pRhoUy_
       'O2pMz   ', & ! O2pRhoUz_
       'O2pP    ', & ! O2pP_
       'OpRho   ', & ! OpRho_
       'OpMx    ', & ! OpRhoUx_
       'OpMy    ', & ! OpRhoUy_
       'OpMz    ', & ! OpRhoUz_
       'OpP     ', & ! OpP_
       'CO2pRho ', & ! CO2pRho_
       'CO2pMx  ', & ! CO2pRhoUx_
       'CO2pMy  ', & ! CO2pRhoUy_
       'CO2pMz  ', & ! CO2pRhoUz_
       'CO2pP   ', & ! CO2pP_
       'HpRho   ', & ! HpRho
       'HpUx    ', & ! HpUx_
       'HpUy    ', & ! HpUy_
       'HpUz    ', & ! HpUz_
       'HpP     ', &  ! HpP_
       'E       ', & ! Energy_
       'HpswE   ', & ! HpswEnergy_
       'O2pE    ', & ! O2pEnergy_
       'OpE     ', & ! OpEnergy_
       'CO2pE   ', & ! CO2pEnergy_
       'HpE     '/)  ! HpEnergy_ 

  ! The space separated list of nVar conservative variables for plotting
  character(len=*), parameter :: NameConservativeVar = &
       'Rho Mx My Mz Bx By Bz E '// &
       'HpswRho HpswMx HpswMy HpswMz HpswE '// &
       'O2pRho O2pMx O2pMy O2pMz O2pE '// &
       'OpRho OpMx OpMy OpMz OpE '// &
       'CO2pRho CO2pMx CO2pMy CO2pMz CO2pE '//&
       'HpRho HpMx HpMy HpMz HpE'

  ! The space separated list of nVar primitive variables for plotting
  character(len=*), parameter :: NamePrimitiveVar = &
       'Rho Ux Uy Uz Bx By Bz P ' // &
       'HpswRho HpswUx HpswUy HpswUz HpswP '// &
       'O2pRho O2pUx O2pUy O2pUz O2pP '// &
       'OpRho OpUx OpUy OpUz OpP '// &
       'CO2pRho CO2pUx CO2pUy CO2pUz CO2pP '//&
       'HpRho HpUx HpUy HpUz HpP'

  ! The space separated list of nVar primitive variables for TECplot output
  character(len=*), parameter :: NamePrimitiveVarTec = &
       '"`r", "U_x", "U_y", "U_z", "B_x", "B_y", "B_z", "p", ' // &
       '"`r^H^+sw", "U_x^H^+sw", "U_y^H^+sw", "U_z^H^+sw", "P^H^+sw", ' // &
       '"`r^O2^+", "U_x^O2^+", "U_y^O2^+", "U_z^O2^+", "P^O2^+" ' // &
       '"`r^O^+", "U_x^O^+", "U_y^O^+", "U_z^O^+", "P^O^+" '// &
       '"`r^CO2^+", "U_x^CO2^+", "U_y^CO2^+", "U_z^CO2^+", "P^CO2^+" '//&
       '"`r^H^+", "U_x^H^+", "U_y^H^+", "U_z^H^+", "P^H^+"'

  ! Names of the user units for IDL and TECPlot output
  character(len=80) :: &
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
