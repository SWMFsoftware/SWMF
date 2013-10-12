!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!#NOTPUBLIC  email:mopher@bu.edu  expires:12/31/2099
module ModVarIndexes

   use ModExtraVariables

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

   integer, parameter :: nVar = 38 !C.P. edited

 
   ! 1 tot ion fluid, 2 single ion fluids and 4 neutral fluids
   integer, parameter :: nFluid    = 7
   integer, parameter :: IonFirst_ = 2
   integer, parameter :: IonLast_  = 3
   logical, parameter :: IsMhd     = .true.
   real               :: MassFluid_I(2:nFluid) = (/ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 /)


   ! All is total ion fluid, SWH is the Solar wind hydrogen fluid, Pu3
   ! are pick up ions produced in region 3 (see mod user), 
   ! and Neu, Ne# are neutrals produced in the correspondign region

   character (len=3), parameter :: &
        NameFluid_I(nFluid) = (/'All', 'SWH', 'Pu3', 'Neu', 'Ne2', 'Ne3', 'Ne4'/) 

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
        P_      =  8, &
        SWHRho_       =  9,          &
        SWHRhoUx_     =  10, SWHUx_ = 10, &
        SWHRhoUy_     =  11, SWHUy_ = 11, &
        SWHRhoUz_     =  12, SWHUz_ = 12, &
        SWHP_         =  13, &
        Pu3Rho_    = 14, &
        Pu3RhoUx_  = 15, Pu3Ux_ = 15, &
        Pu3RhoUy_  = 16, Pu3Uy_ = 16, &
        Pu3RhoUz_  = 17, Pu3Uz_ = 17, &
        Pu3P_      = 18

   integer, parameter:: &
        NeuRho_    = 19, &
        NeuRhoUx_  = 20, NeuUx_ = 20, &
        NeuRhoUy_  = 21, NeuUy_ = 21, &
        NeuRhoUz_  = 22, NeuUz_ = 22, &
        NeuP_      = 23, &
        Ne2Rho_    = 24, &
        Ne2RhoUx_  = 25, Ne2Ux_ = 25, &
        Ne2RhoUy_  = 26, Ne2Uy_ = 26, &
        Ne2RhoUz_  = 27, Ne2Uz_ = 27, &
        Ne2P_      = 28, &
        Ne3Rho_    = 29, &
        Ne3RhoUx_  = 30, Ne3Ux_ = 30, &
        Ne3RhoUy_  = 31, Ne3Uy_ = 31, &
        Ne3RhoUz_  = 32, Ne3Uz_ = 32, &
        Ne3P_      = 33, &
        Ne4Rho_    = 34, &
        Ne4RhoUx_  = 35, Ne4Ux_ = 35, &
        Ne4RhoUy_  = 36, Ne4Uy_ = 36, &
        Ne4RhoUz_  = 37, Ne4Uz_ = 37, &
        Ne4P_      = 38, &   
        Energy_    = nVar+1, &
        SWHEnergy_ = nVar+2, &
        Pu3Energy_ = nVar+3, &
        NeuEnergy_ = nVar+4, &
        Ne2Energy_ = nVar+5, &
        Ne3Energy_ = nVar+6, &
        Ne4Energy_ = nVar+7
	

    ! is this name dependent on the first 
   ! This allows to calculate RhoUx_ as RhoU_+x_ and so on.
   integer, parameter :: U_ = Ux_ - 1, RhoU_ = RhoUx_-1, B_ = Bx_-1

   ! These arrays are useful for multifluid
   integer, parameter :: &
        iRho_I(nFluid)   = (/ Rho_,  SWHRho_,   Pu3Rho_,   NeuRho_,    Ne2Rho_,   Ne3Rho_,   Ne4Rho_ /),   &
        iRhoUx_I(nFluid) = (/ RhoUx_, SWHRhoUx_, Pu3RhoUx_, NeuRhoUx_,  Ne2RhoUx_, Ne3RhoUx_, Ne4RhoUx_ /), &
        iRhoUy_I(nFluid) = (/ RhoUy_, SWHRhoUy_, Pu3RhoUy_, NeuRhoUy_,  Ne2RhoUy_, Ne3RhoUy_, Ne4RhoUy_ /), &
        iRhoUz_I(nFluid) = (/ RhoUz_, SWHRhoUz_, Pu3RhoUz_, NeuRhoUz_,  Ne2RhoUz_, Ne3RhoUz_, Ne4RhoUz_ /), &
        iP_I(nFluid)     = (/ P_,  SWHP_,     Pu3P_,     NeuP_,      Ne2P_,     Ne3P_,     Ne4P_ /)

   ! The default values for the state variables:
   ! Variables which are physically positive should be set to 1,
   ! variables that can be positive or negative should be set to 0:
   real, parameter :: DefaultState_V(nVar+nFluid) = (/ &
        1.0,           & ! Rho_
        0.0, 0.0, 0.0, & ! RhoUx_ .. RhoUz_
        0.0, 0.0, 0.0, & ! Bx_ .. Bz_
        1.0,           & ! P_
        1.0,           & ! SWHRho_
        0.0, 0.0, 0.0, & ! SWHRhoUx_ .. SWHRhoUz_
        1.0,           & ! SWHP_
        1.0,           & ! Pu3Rho_
        0.0, 0.0, 0.0, & ! Pu3RhoUx_ .. Pu3RhoUz_
        1.0,           & ! Pu3P_
        1.0,           & ! NeuRho_
        0.0, 0.0, 0.0, & ! NeuRhoUx_ .. NeuRhoUz_
        1.0,           & ! NeuP_
        1.0,           & ! Ne2Rho_
        0.0, 0.0, 0.0, & ! Ne2RhoUx_ .. Ne2RhoUz_
        1.0,           & ! Ne2P_
        1.0,           & ! Ne3Rho_
        0.0, 0.0, 0.0, & ! Ne3RhoUx_ .. Ne3RhoUz_
        1.0,           & ! Ne3P_
        1.0,           & ! Ne4Rho_
        0.0, 0.0, 0.0, & ! Ne4RhoUx_ .. Ne4RhoUz_
        1.0,           & ! Ne4P_
        1.0,           & ! Energy_
        1.0,           & ! SWHEnergy_
        1.0,           & ! Pu3Energy_
        1.0,           & ! NeuEnergy_
        1.0,           & ! Neu2Energy_
        1.0,           & ! Neu3Energy_
        1.0 /)           ! Neu4Energy_
            
   ! The names of the variables used in i/o
   character(len=6):: NameVar_V(nVar+nFluid) = (/ &
        'Rho   ', & ! Rho_
        'Mx    ', 'My    ', 'Mz    ', & ! RhoUx_  RhoUz_
        'Bx    ', 'By    ', 'Bz    ', & ! Bx_  Bz_
        'p     ', & ! p_
        'SWHRho', & ! SWHRho_
        'SWHMx ', 'SWHMy ', 'SWHMz ', & ! SWHRhoUx_ SWHRhoUz_
        'SWHp  ', & ! SWHp_
        'Pu3Rho', & ! Pu3Rho_
        'Pu3Mx ', 'Pu3My ', 'Pu3Mz ', & ! Pu3RhoUx_ Pu3RhoUz_
        'Pu3P  ', & ! Pu3P_
        'NeuRho', & ! NeuRho_
        'NeuMx ', 'NeuMy ', 'NeuMz ', & ! NeuRhoUx_ NeuRhoUz_
        'NeuP  ', & ! NeuP_
        'Ne2Rho', & ! Ne2Rho_
        'Ne2Mx ', 'Ne2My ', 'Ne2Mz ', & ! Ne2RhoUx_ Ne2RhoUz_
        'Ne2P  ', & ! Ne2P_
        'Ne3Rho', & ! Ne3Rho_
        'Ne3Mx ', 'Ne3My ', 'Ne3Mz ', & ! Ne3RhoUx_ Ne3RhoUz_
        'Ne3P  ', & ! Ne3P_
        'Ne4Rho', & ! Ne4Rho_
        'Ne4Mx ', 'Ne4My ', 'Ne4Mz ', & ! Ne4RhoUx_ Ne4RhoUz_
        'Ne4P  ', & ! Ne4P_
        'E     ', & ! Energy_
        'SWHE  ', & ! SWHEnergy_
        'Pu3E  ', & ! Pu3Energy_
        'NeuE  ', & ! NeuEnergy_
        'Ne2E  ', & ! Ne2Energy_
        'Ne3E  ', & ! Ne3Energy_
        'Ne4E  ' /) ! Ne4Energy_

   ! The space separated list of nVar conservative variables for  plotting
   character(len=*), parameter :: NameConservativeVar = &
        'Rho Mx My Mz Bx By Bz E ' // &
        'SWHRho SWHMx SWHMy SWHMz SWHE Pu3Rho Pu3Mx Pu3My Pu3Mz Pu3E ' // &
        'NeuRho NeuMx NeuMy NeuMz NeuE Ne2Rho Ne2Mx Ne2My Ne2Mz Ne2E ' // &
        'Ne3Rho Ne3Mx Ne3My Ne3Mz Ne3E Ne4Rho Ne4Mx Ne4My Ne4Mz Ne4E' 
	

   ! The space separated list of nVar primitive variables for plotting
   character(len=*), parameter :: NamePrimitiveVar = &
        'Rho Ux Uy Uz Bx By Bz P ' // &
        'SWHRho SWHUx SWHUy SWHUz SWHP Pu3Rho Pu3Ux Pu3Uy Pu3Uz Pu3P ' // & 
        'NeuRho NeuUx NeuUy NeuUz NeuP Ne2Rho Ne2Ux Ne2Uy Ne2Uz Ne2P ' // &
        'Ne3Rho Ne3Ux Ne3Uy Ne3Uz Ne3P Ne4Rho Ne4Ux Ne4Uy Ne4Uz Ne4P'
	

   ! The space separated list of nVar primitive variables for TECplot output
   character(len=*), parameter :: NamePrimitiveVarTec = &
        '"`r", "U_x", "U_y", "U_z", "B_x", "B_y", "B_z", "P", ' // &
        '"SWHRho", "SWHUx", "SWHUy", "SWHUz", "SWHP", ' // &
        '"Pu3Rho", "Pu3Ux", "Pu3Uy", "Pu3Uz", "Pu3P", ' // &
        '"NeuRho", "NeuUx", "NeuUy", "NeuUz", "NeuP", "Ne2Rho", ' // &
        '"Ne2Ux", "Ne2Uy", "Ne2Uz", "Ne2P", "Ne3Rho", "Ne3Ux", ' // &
        '"Ne3Uy", "Ne3Uz", "Ne3P", "Ne4Rho", "Ne4Ux", ' // &
        '"Ne4Uy", "Ne4Uz", "Ne4P" '

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
