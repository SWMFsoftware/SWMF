!^FILE NOT ALWAVES
Module ModVarIndexes
  use ModNumConst

  implicit none

  character (len=*), parameter :: NameEquation='MHD_modif'

! Define the numbur of cell centered variables (8 for a single fluid MHD,
! 5 for hydrodynamics.
! The variables numbered from 1 to nCentered are:
!
! 1. defined in set_ICs.
! 2. prolonged and restricted in AMR
! 3. saved and downloaded in the restart file
! 4. sent and recieved in the exchange message
! 5. filled in the outer ghostcells by the program fill_outer_ghostsells
! 5. integrated by subroutine integrate_all for saving to 
! 6. should be updated by advance_*

  integer, parameter :: nVar=9

! The use of these indexes in the arrays, which are not named as
! CellCenteredVars... can result in an error
  ! Named indexes for CellCentered variables
  integer, parameter:: &
       rho_=1,  &
       rhoUx_=2,&
       rhoUy_=3,&
       rhoUz_=4,&
       Bx_=5,   &
       By_=6,   &
       Bz_=7,   &
       EnergyRL_=8,&
       P_=nVar
  !These indexes should go subsequently, from 1 to nVar.

  integer, parameter::rhoU_=1,B_=4
! This allows to calculate rhoUx_ as rhoU_+x_ and so on.



  !This is an important step: define the default set of the cell centered
  !variables.
  real,dimension(nVar)::DefaultState_V
  data DefaultState_V/&
       cOne,& !rho_
       cZero,&!rhoUx_
       cZero,&!rhoUy_
       cZero,&!rhoUz_
       cZero,&!Bx_
       cZero,&!By_
       cZero,&!Bz_
       cZero,&!EnergyRL_
       cOne/  !P_
  !For the cell centered vars from the first up to nVar 'th
  !the components of the Default array MUST be defined following STRICTLY the
  !following rules: the components which are POSITIVE according to their 
  !physical sense, such as density, pressure and so on,
  !\
  ! MUST be initialized with cOne 
  !/
  !(or something greater than cTiny). If this rule is ignored, then some 
  !idle opertaions in the ghost cells can result in dividing per zero or cal-
  !culating the sqrt(negative). On the other hand, the things which can be both
  !positive and negative, such as all the components of the vector variables
  !(velocity, magnetic field, some definitions of entropy),
  !\
  !MUST be initialized with cZero !
  !/
  !(or smth less then cTiny).Otherwise the code (see correct_monotone_restrict
  !in the exchange_messages) will interprete the physically reasonable negative
  !values for such the variable as the discretization error and will "correct"
  !it by applying some limiters, what will result in the loss in accuracy
  !at least. The array DefaultState_V is only used for initializing the 
  !ghostcells and the unused blocks, in the ideal world and with all the 
  !restrictions listed above the particular values of its components 
  !(say, the change cOne=>cTwo) should make absolutely no difference on
  !any results.

  character(len=*),parameter,dimension(nVar)::NameVar_V = (/ &
       'rho  ', & !rho_
       'rhoUx', & !rhoUx_
       'rhoUy', & !rhoUy_
       'rhoUz', & !rhoUz_
       'Bx   ', & !Bx_
       'By   ', & !By_
       'Bz   ', & !Bz_
       'RLEn ', & !EnergyRL_
       'P    ' /) !P_
  !These are the character names of the variables used in i/o

! Additional convention (related to the program correct_monotone_restrict:
! The variables from rhoUx_ to B_+ndim are not assumed to be definite positive, all
! the others are assumed to be definite positive.

  !Conservative or non-conservative equation for energy or both!
  integer,parameter::     Energy_=nVar+1
  ! The possibility to use both conservative and non-conservative approach for energy
  ! and to switch between them results in a following important restriction. There 
  ! are only two choices for the parameter Energy_, which, according to the !
  ! definition is always the LAST variable component in the arrays for fluxes (see
  ! ModAdvance). These two options are as follows:
  ! OPTION 1. Energy_=nVar - fully conservative approach
  ! (is not used in BATSRUS for a very long time and the use of this option
  ! results in a lot of conflicts between the available subroutines)
  ! OPTION 2. Energy_=nCellsCenteredVars+1 - default in BATSRUS. In this case
  ! the nVar 'th variable should satisfy the following definition:
  ! the equation for this "energy-like" variable is the linear combination
  ! of the conservation laws system, the coefficient for the energy equation in
  ! this combination being essentially non-zero. Examples: pressure of ideal gas 
  ! (the corresponding coefficient being (\gamma-1) in this case),density of entropy
  ! (coefficient is 1/Temperature).



  !Corrected Fluxes:
 
  integer,parameter :: Vdt_=nVar+1
  integer,parameter :: BnL_=nVar+2
  integer,parameter :: BnR_=nVar+3
  integer,parameter :: nCorrectedFaceValues = BnR_

  !The fluxes calculated in the finer blocks are send to the coarser neighboring !
  !block to mantain the scheme conservation property. nCorrectedFaceValues of fluxes
  !are sent recieved, among them are:
  !1. nVar-1 componenents of fluxes and Energy_. So the flux of 
  !the "energy-like" variable is NOT SENT. Also C_max*FaceArea is sent and right 
  !and left normal components of B_ vector at face are sent. If magnetic field is
  ! not involved in the model one can define B_=rhoU_ and so on, hence the normal
  ! components of the velocity will be sent.

  ! Primitive variable names
  integer, parameter :: U_=rhoU_, Ux_=rhoUx_, Uy_=rhoUy_, Uz_=rhoUz_

end Module ModVarIndexes
