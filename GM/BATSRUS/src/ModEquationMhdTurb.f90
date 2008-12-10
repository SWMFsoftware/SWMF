module ModVarIndexes

  use ModSingleFluid

  implicit none

  save

  ! This equation module contains the standard MHD equations plus one
  ! extra wave energy of a single frequency w, Iw that carries the extra energy.
  character (len=*), parameter :: NameEquation='Corona MHD Turb'

  ! The variables numbered from 1 to nVar are:
  !
  ! 1. defined in set_ICs.
  ! 2. prolonged and restricted in AMR
  ! 3. saved into the restart file
  ! 4. sent and recieved in the exchange message
  ! 5. filled in the outer ghostcells by the program set_outer_BCs
  ! 5. integrated by subroutine integrate_all for saving to logfile
  ! 6. should be updated by advance_*

  integer, parameter :: nVar = 59

  ! Named indexes for State_VGB and other variables
  ! These indexes should go subsequently, from 1 to nVar+1.
  ! The energy is handled as an extra variable, so that we can use
  ! both conservative and non-conservative scheme and switch between them.
  integer, parameter :: &
       Rho_   = 1,    &
       RhoUx_ = 2,    &
       RhoUy_ = 3,    &
       RhoUz_ = 4,    &
       Bx_    = 5,    &
       By_    = 6,    &
       Bz_    = 7,    &
       Ew_    = 8,    &
       I01_   = 9,  I02_= 10, I03_= 11, I04_= 12, I05_= 13, &  
       I06_   = 14, I07_= 15, I08_= 16, I09_= 17, I10_= 18, & 
       I11_   = 19,  I12_= 20, I13_= 21, I14_= 22, I15_= 23, &
       I16_   = 24, I17_= 25, I18_= 26, I19_= 27, I20_= 28, &
       I21_   = 29,  I22_= 30, I23_= 31, I24_= 32, I25_= 33, &
       I26_   = 34, I27_= 35, I28_= 36, I29_= 37, I30_= 38, &
       I31_   = 39,  I32_= 40, I33_= 41, I34_= 42, I35_= 43, &
       I36_   = 44, I37_= 45, I38_= 46, I39_= 47, I40_= 48, &
       I41_   = 49, I42_= 50, I43_= 51, I44_= 52, I45_= 53, &
       I46_   = 54, I47_= 55, I48_= 56, I49_= 57, I50_= 58, & 
       p_     = nVar, &
       Energy_= nVar+1  

 ! Modified for wave turbulance (Rona Oran, Nov 2008). Here we add 50 new variables Ixx where xx=01.....50.
 ! The 50 new variables represent the wave energy I of 50 frequency groups of Alfven waves. These are scalar 
 ! quanteties which are advected with the velocity U of the MHD solution and they do not affect the MHD variables.
 ! In the current verion each wave energy is advected according to the equation:
 ! dI/dt+div(U I)=0, where d stands for partial derivative. In LaTex:
 ! \frac{\partial I}{\partial t}+\nabla(\vec{U} I)=0
 ! In A later developement the full wave energy equation will become:
 ! dI/dt+div[(U+Vg*b)I-0.5(divU)dI/dlogw=gamma I
 ! where U and b are vectors, I and gamma are functions of w
 ! In this equation the wave energy advects at the velocity of the medium plus the Alfven waves group velocity
 ! in the direction of the magnetic field. The third term describes evolution of the energy spectrum with 
 ! expansion/contraction of the plasma and the right hand side describes dissipation of wave energy.
 ! NOTE such an equation exists for each wave mode under considiration.
 ! The last equation in LaTex:
 ! \frac{\partial I}{\partial t}+\nabla((\vec{U}+V_g\vec{b})I)+\frac{1}{2}(\nabla\vec{U})\frac{\partial I}{\partial log w}=\gamma I

 ! The freqquency groups range from 1 s^{-1} to 1/100000 s^{-1} in a logarithmic scale. This is in order to cover the range of 
 ! Alfven wave frequancy observed at the chromosphere and photosphere by Hinode.
 ! The energy spectrum can be varied is defined in the ModUser file (ModUserTurbCorona.f90), see additional detailed therein.

  ! This allows to calculate rhoUx_ as rhoU_+x_ and so on.
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
       0.0, & ! RhoUx_
       0.0, & ! RhoUy_
       0.0, & ! RhoUz_
       0.0, & ! Bx_
       0.0, & ! By_
       0.0, & ! Bz_
       0.0, & ! Ew_ 
       1.0,  1.0, 1.0, 1.0, 1.0 ,& ! Iw(1 -5)
       1.0,  1.0, 1.0, 1.0, 1.0 ,& ! Iw(6 -10)    
       1.0,  1.0, 1.0, 1.0, 1.0 ,& ! Iw(11 -15)    
       1.0,  1.0, 1.0, 1.0, 1.0 ,& ! Iw(16 -20)    
       1.0,  1.0, 1.0, 1.0, 1.0 ,& ! Iw(21 -25)    
       1.0,  1.0, 1.0, 1.0, 1.0 ,& ! Iw(26 -30)    
       1.0,  1.0, 1.0, 1.0, 1.0 ,& ! Iw(31 -35)    
       1.0,  1.0, 1.0, 1.0, 1.0 ,& ! Iw(36 -40)    
       1.0,  1.0, 1.0, 1.0, 1.0 ,& ! Iw(41 -45)    
       1.0,  1.0, 1.0, 1.0, 1.0 ,& ! Iw(46 -50)    
       1.0, & ! p_
       1.0 /) ! Energy_ 
 
  ! The names of the variables used in i/o
  character(len=*), parameter :: NameVar_V(nVar+1) = (/ &
       'Rho', & ! Rho_
       'Mx ', & ! RhoUx_
       'My ', & ! RhoUy_
       'Mz ', & ! RhoUz_
       'Bx ', & ! Bx_
       'By ', & ! By_
       'Bz ', & ! Bz_
       'Ew ', & ! Ew_  
       'I01', 'I02', 'I03', 'I04', 'I05', &  
       'I06', 'I07', 'I08', 'I09', 'I10', & ! I01_ t0 I10_ 
       'I11', 'I12', 'I13', 'I14', 'I15', & 
       'I16', 'I17', 'I18', 'I19', 'I20', & ! I11_ t0 I20_  
       'I21', 'I22', 'I23', 'I24', 'I25', & 
       'I26', 'I27', 'I28', 'I29', 'I30', & ! I21_ t0 I30_  
       'I31', 'I32', 'I33', 'I34', 'I35', & 
       'I36', 'I37', 'I38', 'I39', 'I40', & ! I31_ t0 I40_  
       'I41', 'I42', 'I43', 'I44', 'I45', & 
       'I46', 'I47', 'I48', 'I49', 'I50', & ! I41_ t0 I50_
       'p  ', & ! p_                                                                                          
       'e  '/) ! Energy_        
  
  ! The space separated list of nVar conservative variables for plotting
  character(len=*), parameter :: NameConservativeVar = &
       'rho mx my mz bx by bz Ew' // &
       ' I01 I02 I03 I04 I05 I06 I07 I08 I09 I10' // &
       ' I11 I12 I13 I14 I15 I16 I17 I18 I19 I20' //& 
       ' I21 I22 I23 I24 I25 I26 I27 I28 I29 I30' //& 
       ' I31 I32 I33 I34 I35 I36 I37 I38 I39 I40' //& 
       ' I41 I42 I43 I44 I45 I46 I47 I48 I49 I50' //&
       ' e'

  ! The space separated list of nVar primitive variables for plotting
  character(len=*), parameter :: NamePrimitiveVar = &
       'rho ux uy uz bx by bz Ew' // &
        ' I01 I02 I03 I04 I05 I06 I07 I08 I09 I10' //&
        ' I11 I12 I13 I14 I15 I16 I17 I18 I19 I20' //&
        ' I21 I22 I23 I24 I25 I26 I27 I28 I29 I30' //&
        ' I31 I32 I33 I34 I35 I36 I37 I38 I39 I40' //& 
        ' I41 I42 I43 I44 I45 I46 I47 I48 I49 I50' //&
        ' p'
  
! The space separated list of nVar primitive variables for TECplot output
  character(len=*), parameter :: NamePrimitiveVarTec = &
       '"`r", "U_x", "U_y", "U_z", "B_x", "B_y", "B_z", "E_w", ' //&
       '"I01", "I02", "I03", "I04", "I05", "I06", "I07", "I08", "I09", "I10", '// & 
       '"I11", "I12", "I13", "I14", "I15", "I16", "I17", "I18", "I19", "I20", '//&
       '"I21", "I22", "I23", "I24", "I25", "I26", "I27", "I28", "I29", "I30", '//&
       '"I31", "I32", "I33", "I34", "I35", "I36", "I37", "I38", "I39", "I40", '//&
       '"I41", "I42", "I43", "I44", "I45", "I46", "I47", "I48", "I49", "I50"  '//&
       '"p" '

  ! Names of the user units for IDL and TECPlot output
  character(len=20) :: &
       NameUnitUserIdl_V(nVar+1) = '', NameUnitUserTec_V(nVar+1) = ''

  ! The user defined units for the variables
  real :: UnitUser_V(nVar+1) = 1.0

  ! Primitive variable names
  integer, parameter :: U_ = RhoU_, Ux_ = RhoUx_, Uy_ = RhoUy_, Uz_ = RhoUz_

  ! Specify scalar to be advected
  integer, parameter :: ScalarFirst_ = Ew_, ScalarLast_ = I50_

  ! There are no multi-species
  logical, parameter :: UseMultiSpecies = .false.

  ! Declare the following variables to satisfy the compiler
  integer, parameter :: SpeciesFirst_ = 1, SpeciesLast_ = 1
  real               :: MassSpecies_V(SpeciesFirst_:SpeciesLast_)

contains
  subroutine init_mod_equation
    ! Initialize usre units and names for the MHD variables

    call init_mhd_variables

    ! Set the unit and unit name for the wave energy variable
    UnitUser_V(Ew_)        = UnitUser_V(Energy_)
    NameUnitUserTec_V(Ew_) = NameUnitUserTec_V(Energy_)
    NameUnitUserIdl_V(Ew_) = NameUnitUserIdl_V(Energy_)
   
    UnitUser_V(I01_:I50_)        = UnitUser_V(Energy_)
    NameUnitUserTec_V(I01_:I50_) = NameUnitUserTec_V(Energy_)
    NameUnitUserIdl_V(I01_:I50_) = NameUnitUserIdl_V(Energy_)
    
   end subroutine init_mod_equation

end module ModVarIndexes


