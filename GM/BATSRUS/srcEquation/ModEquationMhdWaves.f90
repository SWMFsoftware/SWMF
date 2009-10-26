module ModVarIndexes

  use ModSingleFluid
  use ModExtraVariables,        &
       Redefine1 => ExtraEInt_, &
       Redefine2 => nWave,      &
       Redefine3 => WaveFirst_, &
       Redefine4 => WaveLast_

  implicit none

  save

  ! This equation module contains the standard MHD equations plus one
  ! extra wave energy of a single frequency w, Iw that carries the extra energy.
  character (len=*), parameter :: NameEquation='MHD Waves'



  ! loop variable over spectrum
  integer, private :: iWave

  ! Number of frequency bins in spectrum
  integer, parameter :: nWave =1
  integer, parameter :: nVar = 9 + nWave
  
  ! Array of strings for plotting
  character(len=3)   :: NameWaveVar_I(nWave)

  ! Named indexes for State_VGB and other variables
  ! These indexes should go subsequently, from 1 to nVar+1.
  ! The energy is handled as an extra variable, so that we can use
  ! both conservative and non-conservative scheme and switch between them.

  ! The variables numbered from 1 to nVar are:
  !
  ! 1. defined in set_ICs.
  ! 2. prolonged and restricted in AMR
  ! 3. saved into the restart file
  ! 4. sent and recieved in the exchange message
  ! 5. filled in the outer ghostcells by the program set_outer_BCs
  ! 5. integrated by subroutine integrate_all for saving to logfile
  ! 6. should be updated by advance_*

  integer, parameter :: &
       Rho_   = 1,    &
       RhoUx_ = 2,    &
       RhoUy_ = 3,    &
       RhoUz_ = 4,    &
       Bx_    = 5,    &
       By_    = 6,    &
       Bz_    = 7,    &
       ExtraEInt_    = 8,    &
       WaveFirst_ = 9, &
       WaveLast_  = WaveFirst_+nWave-1, &
       p_     = nVar, &
       Energy_= nVar+1  

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
       0.0, & ! ExtraEInt_ 
       (0.0, iWave=WaveFirst_,WaveLast_), & 
       1.0, & ! p_
       1.0 /) ! Energy_ 
 
  ! The names of the variables used in i/o
  character(len=3) :: NameVar_V(nVar+1) = (/ &
       'Rho', & ! Rho_
       'Mx ', & ! RhoUx_
       'My ', & ! RhoUy_
       'Mz ', & ! RhoUz_
       'Bx ', & ! Bx_
       'By ', & ! By_
       'Bz ', & ! Bz_
       'Ew ', & ! ExtraEInt_  
       ('I??', iWave=1,nWave), & ! Waves to be reset
       'p  ', & ! p_
       'e  '/)  ! Energy_        
  
  ! The space separated list of nVar conservative variables for plotting
  character(len=*), parameter :: NameConservativeVarPref = &
       'rho mx my mz bx by bz Ew' 
     
  character(len=*), parameter :: NameConservativeVarSuff = &
       ' e'
 
  character(len=4*nWave+len(NameConservativeVarPref) + &
       len(NameConservativeVarSuff))      :: NameConservativeVar = &
       trim(NameConservativeVarPref)         ! modified by init_mod_equation

  ! The space separated list of nVar primitive variables for plotting
  character(len=*), parameter :: NamePrimitiveVarPref = &
       'rho ux uy uz bx by bz Ew' 
  character(len=*), parameter :: NamePrimitiveVarSuff = &
       ' p'
  character(len=4*nWave+len(NamePrimitiveVarPref)+len(NamePrimitiveVarSuff)) :: &
       NamePrimitiveVar = trim(NamePrimitiveVarPref)

  ! The space separated list of nVar primitive variables for TECplot output
  character(len=*), parameter :: NamePrimitiveVarTecPref = &
       '"`r", "U_x", "U_y", "U_z", "B_x", "B_y", "B_z", "E_w",'

  character(len=*), parameter :: NamePrimitiveVarTecSuff = &
       ' "p"'

  character(len=7*nWave+len(NamePrimitiveVarTecPref) + &
       len(NamePrimitiveVarTecSuff))     :: NamePrimitiveVarTec = &
       trim(NamePrimitiveVarTecPref)

       
  ! Names of the user units for IDL and TECPlot output
  character(len=20) :: &
       NameUnitUserIdl_V(nVar+1) = '', NameUnitUserTec_V(nVar+1) = ''

  ! The user defined units for the variables
  real :: UnitUser_V(nVar+1) = 1.0

  ! Primitive variable names
  integer, parameter :: U_ = RhoU_, Ux_ = RhoUx_, Uy_ = RhoUy_, Uz_ = RhoUz_

  ! Specify scalar to be advected
  integer, parameter :: ScalarFirst_ = ExtraEInt_, ScalarLast_ = WaveLast_

  ! There are no multi-species
  logical, parameter :: UseMultiSpecies = .false.

  ! Declare the following variables to satisfy the compiler
  integer, parameter :: SpeciesFirst_ = 1, SpeciesLast_ = 1
  real               :: MassSpecies_V(SpeciesFirst_:SpeciesLast_)

contains

  subroutine init_mod_equation
    ! Initialize user units and names for the MHD variables

    call init_mhd_variables

    ! Create variable name strings for plotting
    
    do iWave = 1, nWave
       write(NameWaveVar_I(iWave),'(a,i2.2)') 'I',iWave
       NameVar_V(WaveFirst_+iWave-1) = NameWaveVar_I(iWave)
       NameConservativeVar=trim(NameConservativeVar)//&
            ' '//trim(NameWaveVar_I(iWave))
       NamePrimitiveVar=trim(NamePrimitiveVar)//&
            ' '//trim(NameWaveVar_I(iWave))
       NamePrimitiveVarTec=trim(NamePrimitiveVarTec)//' "'//trim(NameWaveVar_I(iWave))//'",'
    end do
    NameConservativeVar=trim(NameConservativeVar)//trim(NameConservativeVarSuff)
    NamePrimitiveVar=trim(NamePrimitiveVar)//trim(NamePrimitiveVarSuff)
    NamePrimitiveVarTec=trim(NamePrimitiveVarTec)//trim(NamePrimitiveVarSuff)

  end subroutine init_mod_equation

end module ModVarIndexes


