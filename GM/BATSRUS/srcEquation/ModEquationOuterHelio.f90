module ModVarIndexes

   implicit none

   save

   ! This equation module contains the standard MHD equations + 4 neutrals
   character (len=*), parameter :: NameEquation='MHD and four neutrals'

   ! The variables numbered from 1 to nVar are:
   !
   ! 1. defined in set_ICs.
   ! 2. prolonged and restricted in AMR
   ! 3. saved into the restart file
   ! 4. sent and recieved in the exchange message
   ! 5. filled in the outer ghostcells by set_outer_bcs
   ! 5. integrated by subroutine integrate_all for saving to logfile
   ! 6. should be updated by advance_*

   integer, parameter :: nVar = 28

   ! 1 ion fluid and 4 neutral fluids
   integer, parameter :: nFluid    = 5
   integer, parameter :: IonFirst_ = 1
   integer, parameter :: IonLast_  = 1
   logical, parameter :: IsMhd     = .true.
   real               :: MassFluid_I(nFluid) = (/ 1.0, 1.0, 1.0, 1.0, 1.0 /)

   character (len=3), parameter :: &
        NameFluid_I(nFluid) = (/'Ion', 'Neu', 'Ne2', 'Ne3', 'Ne4'/)

   ! Named indexes for State_VGB and other variables
   ! These indexes should go subsequently, from 1 to nVar+nFluid.
   ! The energies are handled as an extra variable, so that we can use
   ! both conservative and non-conservative scheme and switch between  them.
   integer, parameter :: &
        Rho_       =  1,          &
        RhoUx_     =  2, Ux_ = 2, &
        RhoUy_     =  3, Uy_ = 3, &
        RhoUz_     =  4, Uz_ = 4, &
        Bx_        =  5, &
        By_        =  6, &
        Bz_        =  7, &
        p_         =  8, &
        NeuRho_    =  9, &
        NeuRhoUx_  = 10, NeuUx_ = 10, &
        NeuRhoUy_  = 11, NeuUy_ = 11, &
        NeuRhoUz_  = 12, NeuUz_ = 12, &
        NeuP_      = 13, &
        Ne2Rho_    = 14, &
        Ne2RhoUx_  = 15, Ne2Ux_ = 15, &
        Ne2RhoUy_  = 16, Ne2Uy_ = 16, &
        Ne2RhoUz_  = 17, Ne2Uz_ = 17, &
        Ne2P_      = 18, &
        Ne3Rho_    = 19, &
        Ne3RhoUx_  = 20, Ne3Ux_ = 20, &
        Ne3RhoUy_  = 21, Ne3Uy_ = 21, &
        Ne3RhoUz_  = 22, Ne3Uz_ = 22, &
        Ne3P_      = 23, &
        Ne4Rho_    = 24, &
        Ne4RhoUx_  = 25, Ne4Ux_ = 25, &
        Ne4RhoUy_  = 26, Ne4Uy_ = 26, &
        Ne4RhoUz_  = 27, Ne4Uz_ = 27, &
        Ne4P_      = 28, &
        Energy_    = nVar+1, &
        NeuEnergy_ = nVar+2, &
        Ne2Energy_ = nVar+3, &
        Ne3Energy_ = nVar+4, &
        Ne4Energy_ = nVar+5 

   ! This allows to calculate RhoUx_ as RhoU_+x_ and so on.
   integer, parameter :: U_ = Ux_ - 1, RhoU_ = RhoUx_-1, B_ = Bx_-1

   ! These arrays are useful for multifluid
   integer, parameter :: &
        iRho_I(nFluid)   = (/Rho_,   NeuRho_,  Ne2Rho_,   Ne3Rho_, Ne4Rho_ /),     &
        iRhoUx_I(nFluid) = (/RhoUx_, NeuRhoUx_,  Ne2RhoUx_, Ne3RhoUx_, Ne4RhoUx_ /), &
        iRhoUy_I(nFluid) = (/RhoUy_, NeuRhoUy_,  Ne2RhoUy_, Ne3RhoUy_, Ne4RhoUy_ /), &
        iRhoUz_I(nFluid) = (/RhoUz_, NeuRhoUz_,  Ne2RhoUz_, Ne3RhoUz_, Ne4RhoUz_ /), &
        iP_I(nFluid)     = (/p_,     NeuP_,      Ne2P_,     Ne3P_ , Ne4P_/)

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
        1.0, & ! NeuRho_
        0.0, & ! NeuRhoUx_
        0.0, & ! NeuRhoUy_
        0.0, & ! NeuRhoUz_
        1.0, & ! NeuP_
        1.0, & ! Ne2Rho_
        0.0, & ! Ne2RhoUx_
        0.0, & ! Ne2RhoUy_
        0.0, & ! Ne2RhoUz_
        1.0, & ! Ne2P_
        1.0, & ! Ne3Rho_
        0.0, & ! Ne3RhoUx_
        0.0, & ! Ne3RhoUy_
        0.0, & ! Ne3RhoUz_
        1.0, & ! Ne3P_
        1.0, & ! Ne4Rho_
        0.0, & ! Ne4RhoUx_
        0.0, & ! Ne4RhoUy_
        0.0, & ! Ne4RhoUz_
        1.0, & ! Ne4P_
        1.0, & ! Energy_
        1.0, & ! NeuEnergy_
        1.0, & ! Neu2Energy_
        1.0, & ! Neu3Energy_
        1.0 /) ! Neu4Energy__

   ! The names of the variables used in i/o
   character(len=*), parameter :: NameVar_V(nVar+nFluid) = (/ &
        'Rho   ', & ! Rho_
        'Mx    ', & ! RhoUx_
        'My    ', & ! RhoUy_
        'Mz    ', & ! RhoUz_
        'Bx    ', & ! Bx_
        'By    ', & ! By_
        'Bz    ', & ! Bz_
        'p     ', & ! p_
        'NeuRho', & ! NeuRho_
        'NeuMx ', & ! NeuRhoUx_
        'NeuMy ', & ! NeuRhoUy_
        'NeuMz ', & ! NeuRhoUz_
        'NeuP  ', & ! NeuP_
        'Ne2Rho', & ! Ne2Rho_
        'Ne2Mx ', & ! Ne2RhoUx_
        'Ne2My ', & ! Ne2RhoUy_
        'Ne2Mz ', & ! Ne2RhoUz_
        'Ne2P  ', & ! Ne2P_
        'Ne3Rho', & ! Ne3Rho_
        'Ne3Mx ', & ! Ne3RhoUx_
        'Ne3My ', & ! Ne3RhoUy_
        'Ne3Mz ', & ! Ne3RhoUz_
        'Ne3P  ', & ! Ne3P_
        'Ne4Rho', & ! Ne3Rho_
        'Ne4Mx ', & ! Ne3RhoUx_
        'Ne4My ', & ! Ne3RhoUy_
        'Ne4Mz ', & ! Ne3RhoUz_
        'Ne4P  ', & ! Ne3P_
        'e     ', & ! Energy_
        'NeuE  ', & ! NeuEnergy_
        'Ne2E  ', & ! Ne2Energy_
        'Ne3E  ', & ! Ne3Energy_
        'Ne4E  ' /) ! Ne4Energy_

   ! The space separated list of nVar conservative variables for  plotting
   character(len=*), parameter :: NameConservativeVar = &
        'Rho Mx My Mz Bx By Bz e NeuRho NeuMx NeuMy ' // &
        'NeuMz NeuE Ne2Rho Ne2Mx Ne2My Ne2Mz Ne2E Ne3Rho ' // &
        'Ne3Mx Ne3My Ne3Mz Ne3E Ne4Rho Ne4Mx Ne4My Ne4Mz Ne4E'

   ! The space separated list of nVar primitive variables for plotting
   character(len=*), parameter :: NamePrimitiveVar = &
        'Rho Ux Uy Uz Bx By Bz p NeuRho NeuUx NeuUy NeuUz ' // &
        'NeuP Ne2Rho Ne2Ux Ne2Uy Ne2Uz Ne2P Ne3Rho Ne3Ux ' //  &
        'Ne3Uy Ne3Uz Ne3P Ne4Rho Ne4Ux Ne4Uy Ne4Uz Ne4P'

   ! The space separated list of nVar primitive variables for TECplot output
   character(len=*), parameter :: NamePrimitiveVarTec = &
        '"`r", "U_x", "U_y", "U_z", "B_x", "B_y", "B_z", "p", ' // &
        '"NeuRho", "NeuUx", "NeuUy", "NeuUz", "NeuP", "Ne2Rho", ' // &
        '"Ne2Ux", "Ne2Uy", "Ne2Uz", "Ne2P", "Ne3Rho", "Ne3Ux", ' // &
        '"Ne3Uy", "Ne3Uz", "Ne3P", "Ne4Rho", "Ne4Ux", ' // &
        '"Ne4Uy", "Ne4Uz", "Ne4P"'

   ! Names of the user units for IDL and TECPlot output
   character(len=20) :: &
        NameUnitUserIdl_V(nVar+nFluid) = '', &
        NameUnitUserTec_V(nVar +nFluid) = ''

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
