!^CFG COPYRIGHT UM
module CRASH_ModEos

  ! Equation Of State (EOS) for ionized plasma
  !
  ! Thermodynamic variables and other notations
  !
  !        Rho - the mass density
  !        E - internal energy of the unit of mass
  !        e, i - electron, ion
  !        V, vol - volume or volumetric
  !        \left(\frac{\partial...}{\partial...}\right)_V - thermodynamical
  !             derivative at constant volume
  !        Te, Ti - electron and ion temperature
  !        iMaterial - integer index of the material:
  !        iMaterial=0 - xenon
  !        iMaterial=1 - beryllium   
  !        iMaterial=2 - plastic (or polyimide, if more than one plastic is used
  !        iMaterial=3 - gold
  !        iMaterial=4 - Acrylic (acronim is Ay_, as long as Ac and Ar are both in use.
  !        iMaterial=90 - plasma with eos E=aT^4/4; p=aT^4/12; C_V=aT^3 
  !
  ! In the initial CRASH treatment of materials,
  !
  ! 1. We can treat any given material as having a single average Z.
  !
  ! 2. Mixtures can be treated either as an average material or as mixtures.
  !
  ! 3. If mixtures are treated as mixtures, then collisional rates should be 
  !    calculated using the "effective Z", which is the average of Z squared 
  !    divided by the average Z.
  !
  ! 4. The average Z can be determined from the Saha equation, but must not 
  !    exceed the nuclear charge of the material in question.
  !
  ! 5. We do not need to account for electron degeneracy in the initial model.
  !
  ! 6. In our regime of interest, the electrons behave as an ideal gas in an 
  !    ion-sphere environment within which Coulomb interactions do affect the 
  !    electron pressure and internal energy. The electron pressure and 
  !    internal energy are best calculated using equations 3.47 through 3.50 
  !    in R. P. Drake, High Energy Density Physics
  !
  ! 7. The ion pressure is the ideal gas pressure. The ion internal energy 
  !    includes the particle energy of random motion and the energy of 
  !    ionization. The model of eqs. 3.74 through 3.76 in the mentioned book 
  !    is acceptable. Alternatively, a more complex model using actual 
  !    ionization energies would be acceptable.
  !!WARNING!!!
  !Correction in this item. Since the ionization partition function is controlled
  !by the electron temperature, the ionization energy as well as not mentioned
  !excitation energy are both included to the ELECTRON ENERGY DENSITY. See detail
  !in HEDP.pdf. To make this document go to util/CRASH/doc/Tex directory and
  !make PDF
  !
  ! 8. The materials that matter are
  !    - Beryllium
  !    - Xenon
  !    - Polyimide (C_22 H_10 N_2 O_5)
  !
  ! Error flag (iError) values:
  ! 0 - OK 
  ! 2 - relativistic temperature
  !\
  !!   WARNING !!!
  !You cannot use total pressure and total energy density as input or output
  !parameters, if the electron temperature is not equal to ion temperature.
  !In this case ONLY electron energy density and electron pressure may be 
  !used.
  !/

  use CRASH_ModPolyimide
  use CRASH_ModAcrylic
  use CRASH_ModStatSum
  use CRASH_ModAtomicMass
  use CRASH_ModPowerLawEos
  use CRASH_ModFermiGas, ONLY: UseFermiGas, LogGeMinBoltzmann, LogGeMinFermi
  use CRASH_ModMultiGroup, ONLY: meshhv, abscon, nGroup, &
       OpacityPlanck_I, OpacityRosseland_I, opacys

  implicit none

  private !Except

  integer, public, parameter:: Xe_=0      ! Xenon
  integer, public, parameter:: Be_=1      ! Beryllium
  integer, public, parameter:: Plastic_=2 ! Polyimide (C_22 H_10 N_2 O_5)
  integer, public, parameter:: Au_=3      ! Gold
  integer, public, parameter:: Ay_=4      ! Acrylic

  public:: cAtomicMass_I, cAPolyimide ! inherited from ModPolyimide
  public:: cAAcrylic

  public:: UsePreviousTe ! inherited from CRASH_ModStatSum
  public:: eos, read_eos_parameters, fix_hyades_state
  interface eos
     module procedure eos_material
     module procedure eos_mixture
  end interface

  ! Local variables

  ! test material with the EOS e \propto T^4, p \propto T^4
  integer, parameter:: Test_ = 90 
  integer, parameter :: nZ_I(Xe_:Ay_)=(/54, 4, 6, 79, 6/)

  real, parameter, dimension(Xe_:Ay_),public :: cAtomicMassCRASH_I=&
       (/cAtomicMass_I(54),   &!  Xe
         cAtomicMass_I(4),    &!  Be
         cAPolyimide,         &!  Pl
         cAtomicMass_I(79),   &!  Au
         cAAcrylic/)

  character(LEN=2), public ::&
       NameMaterial_I(Xe_:Ay_) = (/'Xe','Be','Pl','Au','Ay'/)

  logical, public, dimension(Xe_:Ay_) :: &
       UseEosTable_I = .false., UseOpacityTable_I = .false.

  integer :: P_      =  1, &
             E_      =  2, &
             Pe_     =  3, &
             Ee_     =  4, &
             Cv_     =  5, &
             Cve_    =  6, &
             Gamma_  =  7, &
             GammaE_ =  8, &
             TeTi_   =  9, &
             Cond_   = 10, &
             Z_      = 11, &
             Z2_     = 12

  integer :: nVarEos =12   
  
  character(LEN=100):: NameVarEos = &
       'P E Pe Ee Cv Cve Gamma GammaE TeTi Cond Z Z2'


contains

  !============================================================================
  subroutine read_name_var_eos
    use ModReadParam, ONLY: read_var
    use ModUtilities, ONLY: split_string, lower_case
    
    integer, parameter:: MaxString = 200
    character(LEN=20):: NameVar_I(MaxString)
    integer:: iVar
    !---------------------

    call read_var('NameVarEos', NameVarEos)
    call split_string(NameVarEos, MaxString,  NameVar_I, nVarEos)
    
    !Reset named indices
    P_ = -1; E_ = -1; Pe_ = -1; Ee_ = -1
    Cv_ = -1; Cve_ = -1; Gamma_ = -1; GammaE_ = -1
    Cond_ = -1; TeTi_ = -1; Z_ = -1; Z2_ = -1

    do iVar = 1, nVarEos
       call lower_case(NameVar_I(iVar))
       select case(trim(NameVar_I(iVar)))
       case('p')
          P_ = iVar
       case('e')
          E_ = iVar
       case('pe')
          Pe_ = iVar
       case('ee')
          Ee_ = iVar
       case('cv')
          Cv_ = iVar
       case('cve')
          Cve_= iVar
       case('gamma')
          Gamma_= iVar
       case('gammae')
          GammaE_ = iVar
       case('cond')
          Cond_ = iVar
       case('teti', 'tite')
          TeTi_ = iVar
       case('z')
          Z_ = iVar
       case('z2')
          Z2_= iVar
       end select
    end do
  end subroutine read_name_var_eos
  !===============================
  subroutine read_eos_parameters

    ! Usage (with default values shown):
    !
    ! #EOS
    ! T                     UseFermiGas
    ! 4.0                   LogGeMinBoltzmann
    ! 0.0                   LogGeMinFermi
    ! 
    ! Recommended value for the last parameter: -4.0

    use ModReadParam, ONLY: read_var
    !-----------------------------------------------------------------------
    ! For now. But it should/could read other things
    call read_var('UseFermiGas',         UseFermiGas      )
    call read_var('LogGeMinBoltzmann',   LogGeMinBoltzmann)
    call read_var('LogGeMinFermi',       LogGeMinFermi    )

  end subroutine read_eos_parameters

  !============================================================================

  subroutine eos_material(iMaterial,Rho,&
       TeIn, eTotalIn, pTotalIn, eElectronIn, pElectronIn,   &
       TeOut, eTotalOut, pTotalOut, GammaOut, CvTotalOut,    &
       eElectronOut, pElectronOut, GammaEOut, CvElectronOut, &
       OpacityPlanckOut_I, OpacityRosselandOut_I,            &
       HeatCond, TeTiRelax, Ne, iError)
    use ModLookupTable, ONLY: i_lookup_table, interpolate_lookup_table
    ! Eos function for single material

    integer, intent(in):: iMaterial     ! index of material
    real,    intent(in):: Rho           ! mass density [kg/m^3]
    !\
    !!   WARNING !!!
    !You cannot use total pressure and total energy density as input or output
    !parameters, if the electron temperature is not equal to ion temperature.
    !In this case ONLY electron energy density and electron pressure may be 
    !used.
    !/

    ! One of the following five energetic input parameters must be present
    real,    optional, intent(in)  :: TeIn         ! temperature SI[K]
    real,    optional, intent(in)  :: eTotalIn     ! internal energy density
    real,    optional, intent(in)  :: pTotalIn     ! pressure
    real,    optional, intent(in)  :: eElectronIn  ! internal energu density of electrons
    real,    optional, intent(in)  :: pElectronIn  ! pressure of electrons
    
    ! One or more of the output parameters can be present
    real,    optional, intent(out) :: TeOut        ! temperature
    real,    optional, intent(out) :: pTotalOut    ! pressure
    real,    optional, intent(out) :: eTotalOut    ! internal energy density
    real,    optional, intent(out) :: GammaOut     ! polytropic index
    real,    optional, intent(out) :: CvTotalOut   ! specific heat / unit volume
                                                   ! Electrons !!!!!!   
    real,    optional, intent(out) :: pElectronOut ! pressure
    real,    optional, intent(out) :: eElectronOut ! internal energy density
    real,    optional, intent(out) :: GammaEOut    ! polytropic index
    real,    optional, intent(out) :: CvElectronOut! specific heat / unit volume
    real,    optional, intent(out) :: Ne           ! electron concentration [m-3]

    real,    optional, intent(out), &              ! Opacities
                   dimension(nGroup) :: OpacityPlanckOut_I, OpacityRosselandOut_I

    real,    optional, intent(out) :: HeatCond     ! electron heat conductivity (SI)
    real,    optional, intent(out) :: TeTiRelax    ! electron-ion interaction rate (SI)
    integer, optional, intent(out) :: iError       ! error flag

    real   :: Natomic
    real   :: Te, Value_V(1:nVarEos)
    integer:: iTable
    character(LEN=*), parameter:: NameSub = 'eos_material'
    !-------------------------------------------------------------------------
    if(iMaterial == Test_)then
       call eos_esimt4(TeIn, eTotalIn, pTotalIn, &
            TeOut, eTotalOut, pTotalOut, GammaOut, CvTotalOut)
       if(present(iError))iError = 0
       RETURN
    end if

    !Get the atomic concentration
    Natomic = Rho /  ( cAtomicMass * cAtomicMassCRASH_I(iMaterial) )

    if(UseEosTable_I(iMaterial))then
       iTable = i_lookup_table(NameMaterial_I(iMaterial)//'_eos')
       
       if(present(TeIn))then

          Te = TeIn * cKToEV
          call interpolate_lookup_table(iTable, Te, Natomic, Value_V, DoExtrapolate=.false.)

       elseif(present(eTotalIn))then

          ! Get an energy per the atomic cell, express in eV
          ! Find temperature from dentity and internal energy
          call interpolate_lookup_table(iTable, E_, eTotalIn/ (cEV * Natomic), Natomic, &
               Value_V, Arg1Out = Te, DoExtrapolate=.false.)
       

       elseif(present(pTotalIn))then
          ! Divide pressure by Na , express in eV
          !Find temperature from dentity and pressure
          call interpolate_lookup_table(iTable, P_,  pTotalIn / (cEV * Natomic), Natomic, &
               Value_V, Arg1Out = Te, DoExtrapolate=.false.)
          
     
       elseif(present(eElectronIn))then
          ! Get an energy per the atomic cell, express in eV
          ! Find temperature from dentity and internal energy
          call interpolate_lookup_table(iTable, Ee_, eElectronIn/ (cEV * Natomic), Natomic, &
               Value_V, Arg1Out = Te, DoExtrapolate=.false.)


       elseif(present(pElectronIn))then

          ! Divide pressure by Na , express in eV
          !Find temperature from dentity and pressure
          call interpolate_lookup_table(iTable, Pe_,  pElectronIn / (cEV * Natomic), Natomic, &
               Value_V, Arg1Out = Te, DoExtrapolate=.false.)
          
      
       else
          call CON_stop(NameSub// &
               ': none of Te, eTotal, or pTotal is among the input parameters')
       end if

       if(present(TeOut))      TeOut     = Te*cEVToK
       if(present(eTotalOut))  eTotalOut = Natomic*cEV*Value_V(E_)
       if(present(pTotalOut))  pTotalOut = Natomic*cEV*Value_V(P_)
       if(present(eElectronOut)) eElectronOut = Natomic*cEV*Value_V(Ee_)
       if(present(pElectronOut)) pElectronOut = Natomic*cEV*Value_V(Pe_)
       if(present(GammaEOut))  GammaEOut = Value_V(GammaE_)
       if(present(GammaOut))   GammaOut  = Value_V(Gamma_)
       if(present(CvTotalOut)) CvTotalOut = (Natomic*cBoltzmann)*Value_V(Cv_)
       if(present(CvElectronOut)) CvElectronOut = (Natomic*cBoltzmann)*Value_V(Cve_)
   
       if(present(HeatCond))   HeatCond  = Value_V(Cond_)
       if(present(TeTiRelax))  TeTiRelax = Value_V(TeTi_)
       if(present(Ne))         Ne = Value_V(Z_) * NAtomic
       
       return
    end if

    if(iMaterial == Plastic_)then
       call set_mixture(nPolyimide, nZPolyimide_I, CPolyimide_I)
       
    elseif(iMaterial == Ay_)then
       call set_mixture(nAcrylic, nZAcrylic_I, cAcrylic_I)

    else
       call set_element(nZ_I(iMaterial))

    end if
    call eos_generic(Natomic, &
         TeIn, eTotalIn, pTotalIn, eElectronIn, pElectronIn,   &
         TeOut, eTotalOut, pTotalOut, GammaOut, CvTotalOut,    &
         eElectronOut, pElectronOut, GammaEOut, CvElectronOut, & 
         OpacityPlanckOut_I, OpacityRosselandOut_I,            &
         HeatCond, TeTiRelax, Ne, iError)

  end subroutine eos_material

  !============================================================================
  !Cannot be used for mixed-cell simulations if gold and/or Acrylic is used
  subroutine eos_mixture(RhoToARatio_I,&
       TeIn, eTotalIn, pTotalIn, eElectronIn, pElectronIn,   &
       TeOut, eTotalOut, pTotalOut, GammaOut, CvTotalOut,    &
       eElectronOut, pElectronOut, GammaEOut, CvElectronOut, &
       OpacityPlanckOut_I, OpacityRosselandOut_I,            & 
       HeatCond, TeTiRelax, Ne, iError)
    !\
    !!   WARNING !!!
    !You cannot use total pressure and total energy density as input or output
    !parameters, if the electron temperature is not equal to ion temperature.
    !In this case ONLY electron energy density and electron pressure may be 
    !used.
    !/

    ! Eos function for mixed material
    real, intent(in) :: RhoToARatio_I(Xe_:Plastic_) ! Mass densities/A

    ! One of the following five energetic input parameters must be present
    real,    optional, intent(in)  :: TeIn         ! temperature
    real,    optional, intent(in)  :: eTotalIn     ! internal energy density
    real,    optional, intent(in)  :: pTotalIn     ! pressure
    real,    optional, intent(in)  :: eElectronIn  ! internal energu density of electrons
    real,    optional, intent(in)  :: pElectronIn  ! pressure of electrons

    ! One or more of the output parameters can be present
    real,    optional, intent(out) :: TeOut        ! temperature
    real,    optional, intent(out) :: pTotalOut    ! pressure
    real,    optional, intent(out) :: eTotalOut    ! internal energy density 
    real,    optional, intent(out) :: GammaOut     ! polytropic index
    real,    optional, intent(out) :: CvTotalOut   ! specific heat per volume
    ! Electrons !!!!!!   
    real,    optional, intent(out) :: pElectronOut ! pressure
    real,    optional, intent(out) :: eElectronOut ! internal energy density
    real,    optional, intent(out) :: GammaEOut    ! polytropic index
    real,    optional, intent(out) :: CvElectronOut! specific heat / unit volume
    real,    optional, intent(out) :: Ne           ! electron concentration, [m-3]

    real,    optional, intent(out), &              !Opacities m^-1
                   dimension(nGroup) :: OpacityPlanckOut_I, OpacityRosselandOut_I


    real,    optional, intent(out) :: HeatCond     ! electron heat conductivity (SI)
    real,    optional, intent(out) :: TeTiRelax    ! electron-ion interaction rate (SI)
    integer, optional, intent(out) :: iError       ! error flag

    real :: RhoToATotal, Natomic

    integer, parameter :: nAll = 1 + 1 + nPolyimide   

    integer, parameter :: nZAll_I(nAll) = (/54 , &  !Xe
         4 , &  !Be
         6 , &  !C
         1 , &  !H
         7 , &  !N
         8 /)   !O
    real :: ConcentrationAll_I(nAll)
    !-------------------------------------------------------------------------
    RhoToATotal = sum( RhoToARatio_I ) 

    !Relative atomic concentrations of Xe, Be and polyimide:
    ConcentrationAll_I( 1:3 ) = RhoToARatio_I / RhoToATotal

    !Specify concentrations for C, H, N, O

    ConcentrationAll_I( 3:6 ) = ConcentrationAll_I( 3 ) * CPolyimide_I

    call set_mixture(nAll, nZAll_I, ConcentrationAll_I)

    !Get the atomic concentration
    Natomic = RhoToATotal / cAtomicMass 

    call eos_generic(Natomic, &
         TeIn, eTotalIn, pTotalIn, eElectronIn, pElectronIn,   &
         TeOut, eTotalOut, pTotalOut, GammaOut, CvTotalOut,    &
         eElectronOut, pElectronOut, GammaEOut, CvElectronOut, &
         OpacityPlanckOut_I, OpacityRosselandOut_I,            & 
         HeatCond, TeTiRelax, Ne,  iError)

  end subroutine eos_mixture

  !============================================================================

  subroutine eos_generic(Natomic, &
       TeIn, eTotalIn, pTotalIn, eElectronIn, pElectronIn,   &
       TeOut, eTotalOut, pTotalOut, GammaOut, CvTotalOut,    &
       eElectronOut, pElectronOut, GammaEOut, CvElectronOut, & 
       OpacityPlanckOut_I, OpacityRosselandOut_I,            &
       HeatCond, TeTiRelax, Ne, iError)
    use CRASH_ModTransport, ONLY: electron_heat_conductivity, te_ti_relaxation
    use CRASH_ModPartition, ONLY: zAv
    !\
    !!   WARNING !!!
    !You cannot use total pressure and total energy density as input or output
    !parameters, if the electron temperature is not equal to ion temperature.
    !In this case ONLY electron energy density and electron pressure may be 
    !used.
    !/

    real,              intent(in)  :: Natomic      ! Atomic concentration
    real,    optional, intent(in)  :: TeIn         ! temperature
    real,    optional, intent(in)  :: eTotalIn     ! internal energy density
    real,    optional, intent(in)  :: PtotalIn     ! pressure
    real,    optional, intent(in)  :: eElectronIn  ! internal energu density of electrons
    real,    optional, intent(in)  :: pElectronIn  ! pressure of electrons


    real,    optional, intent(out) :: TeOut        ! temperature
    real,    optional, intent(out) :: pTotalOut    ! pressure
    real,    optional, intent(out) :: eTotalOut    ! internal energy density 
    real,    optional, intent(out) :: GammaOut     ! polytropic index
    real,    optional, intent(out) :: CvTotalOut   ! Specific heat per volume
    ! Electrons !!!!!!   
    real,    optional, intent(out) :: pElectronOut ! pressure
    real,    optional, intent(out) :: eElectronOut ! internal energy density
    real,    optional, intent(out) :: GammaEOut    ! polytropic index
    real,    optional, intent(out) :: CvElectronOut! specific heat / unit volume
    real,    optional, intent(out) :: Ne           ! electron concentration, [m-3]

    real,    optional, intent(out), &              ! Opacities m^-1
                   dimension(nGroup) :: OpacityPlanckOut_I, OpacityRosselandOut_I


    real,    optional, intent(out) :: HeatCond     ! electron heat conductivity (SI)
    real,    optional, intent(out) :: TeTiRelax    ! electron-ion interaction rate (SI)
    integer, optional, intent(out) :: iError       ! error flag

    real :: ePerAtom, pPerAtom,TeInEV  !All in eV

    character (len=*), parameter:: NameSub='CRASH_ModEos::eos'
    !----------------------------------------------------------------------!
    if(present(TeIn))then

       TeInEV = TeIn * cKToEV
       call set_ionization_equilibrium(TeInEV, Natomic, iError)

    elseif(present(eTotalIn))then

       ! Get an energy per the atomic cell, express in eV
       ePerAtom = eTotalIn/ (cEV * Natomic)

       ! Find temperature from dentity and internal energy

       call set_temperature(ePerAtom, Natomic, iError)

    elseif(present(pTotalIn))then
       ! Divide pressure by Na , express in eV
       pPerAtom = pTotalIn / (cEV * Natomic)

       !Find temperature from dentity and pressure
       call pressure_to_temperature(pPerAtom, Natomic, iError)
    elseif(present(eElectronIn))then

       ! Get an energy per the atomic cell, express in eV
       ePerAtom = eElectronIn/ (cEV * Natomic)

       ! Find temperature from dentity and internal energy

       call u_e_to_temperature(ePerAtom, Natomic, iError)

    elseif(present(pElectronIn))then
       ! Divide pressure by Na , express in eV
       pPerAtom = pElectronIn / (cEV * Natomic)

       !Find temperature from dentity and pressure
       call pressure_e_to_temperature(pPerAtom, Natomic, iError)
    else
       call CON_stop(NameSub// &
            ': none of Te, eTotal, or pTotal is among the input parameters')
    end if

    if(present(TeOut))      TeOut     = Te*cEVToK
    if(present(eTotalOut))  eTotalOut = Natomic*cEV*internal_energy()
    if(present(pTotalOut))  pTotalOut = pressure()
    if(present(GammaOut))   call get_gamma(GammaSOut=GammaOut)
    if(present(eElectronOut)) eElectronOut = Natomic*cEV*internal_energy_e()
    if(present(pElectronOut))  pElectronOut = pressure_e()
    if(present(GammaEOut))   call get_gamma(GammaSeOut=GammaEOut)
    if(present(OpacityPlanckOut_I).or.present(OpacityRosselandOut_I))then
       call meshhv
       call abscon
       call opacys(TRadIn = Te)
       if(present(OpacityPlanckOut_I)) &
            OpacityPlanckOut_I = OpacityPlanck_I(1:nGroup)*100.0       ![m^-1]
       if(present(OpacityRosselandOut_I)) &
            OpacityRosselandOut_I = OpacityRosseland_I(1:nGroup)*100.0 ![m^-1]
    end if

    if(present(HeatCond))   HeatCond = electron_heat_conductivity()
    if(present(TeTiRelax))  TeTiRelax = te_ti_relaxation()
    if(present(CvTotalOut)) CvTotalOut = (Natomic*cBoltzmann)*heat_capacity()
    if(present(CvElectronOut)) CvElectronOut = (Natomic*cBoltzmann)*heat_capacity_e()
    if(present(Ne)) Ne = NAtomic * zAv
  end subroutine eos_generic
  !=========================
  subroutine fix_hyades_state(iMaterial, StateCgs_V, PMinSi)
    use ModConst
    integer,intent(in)         :: iMaterial
    real   ,intent(inout)      :: StateCgs_V(4) !Rho[Cgs], P[Cgs], Te[KeV], Ti[Kev]
    real, OPTIONAL, intent(in) :: PMinSi

    real:: DensitySi, NAtomicSi, PressureSi, TeSi, TiSi
    !---------------------------------
    DensitySi  = 1.0e3 * StateCgs_V(1)
    NAtomicSi  = DensitySi/(cAtomicMassCRASH_I(iMaterial)*cAtomicMass)
    TeSi       = 1.0e3 * StateCgs_V(3) * ceVToK
    TiSi       = 1.0e3 * StateCgs_V(4) * ceVToK
    call eos(iMaterial, DensitySi, TeIn = TeSi, pTotalOut = PressureSi)
    PressureSi = PressureSi + (TiSi - TeSi) * nAtomicSi * cBoltzmann
    StateCgs_V(2) = 10.0 * PressureSi
    if(present(PMinSi))then
       if(PressureSi < PMinSi)then
          PressureSi = PMinSi
          StateCgs_V(2) = 10.0 * PressureSi
          call eos(iMaterial, DensitySi, pTotalIn = PressureSi, TeOut = TeSi)
          StateCgs_V(3:4) = TeSi * cKToeV * 1.0e-3
       end if
    end if
  end subroutine fix_hyades_state
end module CRASH_ModEos
