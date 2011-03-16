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
  use CRASH_ModPolyimide,ONLY:cAPolyimide
  use CRASH_ModAcrylic, ONLY:cAAcrylic
  use CRASH_ModAtomicDataMix, ONLY: nMixMax
  use CRASH_ModStatSum
  use CRASH_ModAtomicMass
  use CRASH_ModPowerLawEos
  use CRASH_ModFermiGas, ONLY: UseFermiGas, LogGeMinBoltzmann, LogGeMinFermi
  use CRASH_ModMultiGroup, ONLY: meshhv, abscon, nGroup, &
       OpacityPlanck_I, OpacityRosseland_I, opacys
  use ModLookupTable, ONLY: MaxTable
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

  !This the main eos function, which may be implemented both via 
  !internal or external  eos tables and via the built-in EOS model 
  public:: eos

  interface eos
     module procedure eos_material
     module procedure eos_mixture
  end interface

  ! Local variables

  ! test material with the EOS e \propto T^4, p \propto T^4
  integer, parameter:: Test_ = 90 
  integer, parameter:: nMaterialMax = 5
  integer           :: nZ_I(0:nMaterialMax-1)=&
                                   (/54,&  !Xenon 
                                      4,&  !Beryllium
                                     -4,&  !Minus number of elements in polyimide
                                     79,&  !Gold
                                     -3/)  !Minus number of elements in acrylic

  real, dimension(0:nMaterialMax-1),public :: cAtomicMassCRASH_I=&
       (/cAtomicMass_I(54),          &!  Xe
         cAtomicMass_I( 4),          &!  Be
         (cAtomicMass_I(6) * 22.0 + cAtomicMass_I(1) * 10.0 + cAtomicMass_I(7) *  2.0 &
         + cAtomicMass_I(8) *  5.0) / (22.0 + 10.0 + 2.0 +5.0),                          &!  Pl
         cAtomicMass_I(79),          &!  Au
         (cAtomicMass_I(6) * 5.0 + cAtomicMass_I(8) * 2.0 + cAtomicMass_I(1) * 8.0 ) / 15/) !Ay
 
  character(LEN=2), public ::&
       NameMaterial_I(Xe_:Ay_) = (/'Xe','Be','Pl','Au','Ay'/)

  !Arrays for mixtures
  integer, dimension(1:nMixMax, 0:nMaterialMax-1):: nZMix_II=reshape((/&
                                   1, 0, 0, 0, 0, 0, &
                                   1, 0, 0, 0, 0, 0, &
                                   6, 1, 7, 8, 0, 0, &
                                   1, 0, 0, 0, 0, 0, &
                                   6, 8, 1, 0, 0, 0/),(/6,5/))

  real, dimension(1:nMixMax, 0:nMaterialMax-1) :: cMix_II=reshape((/&
                                   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                                   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                                   22.0/39, 10.0/39, 2.0/39, 5.0/39,0.0,0.0,&
                                   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                                   5.0/15, 2.0/15, 8.0/15,0.0,0.0,0.0/),(/6,5/))


  ! The logicals determine if we use or not tabulated EOS
  ! and opacities  
  
  logical, public, dimension(0:nMaterialMax-1) :: &
       UseEosTable_I = .false., UseOpacityTable_I = .false.
  
  !Subroutine which may be used to set/reset UseEosTable_I

  public:: read_if_use_eos_table, read_if_use_opac_table
 
  !The columns in the EOS table
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

  !The number of columns in the EOS table
  integer :: nVarEos =12   
  
  !The variable names in the EOS table
  character(LEN=100):: NameVarEos = &
       'P E Pe Ee Cv Cve Gamma GammaE TeTi Cond Z Z2'

  !The following subroutine may be used to reset the
  !list of variables tabulated in the EOS tables.
  !Accordingly the named indexes are resen as well as
  !nVarEos. The list may be shortened but not extended.
  !The named indices for the excluded variables are set to
  !-1.

  public:: read_name_var_eos 

  !The number of columns in the Opac table
  integer :: nVarOpac =2   
  
  !The variable names in the opacity table
  character(LEN=100):: NameVarOpac = &
       'Opac(1) Ross(1)'
  

  !The following subroutine:
  !1. Checks if the tables are available for the
  !   materials for which UseEosTable_I is True
  !2. If the table is described in the PARAM.in file and should be
  !   made, it is made  
  !3. If the table is not described in the PARAM.in file form it with
  !   the default ranges and fill in using the built-in EOS model

  public:: check_eos_table

  !The following subroutine:
  !1. Checks if the opacity tables are available for the
  !   materials for which UseOpacityTable_I is True
  !2. If the table is described in the PARAM.in file and should be
  !   made, it is made  
  !3. If the table is not described in the PARAM.in file form it with
  !   the default ranges and fill in using the built-in opacity model

  public:: check_opac_table

  !Arrays which relate the iTable for the EOS table with 
  !the material number:
  integer:: iMaterial4EosTable_I(MaxTable) = -1
  integer:: iTableEos4Material_I(0:nMaterialMax-1) = -1

  !Arrays which relate the iTable for the opacity table with 
  !the material number:
  integer:: iMaterial4OpacTable_I(MaxTable) = -1
  integer:: iTableOpac4Material_I(0:nMaterialMax-1) = -1

  !\
  ! Miscellaneous subroutnies (probably, redundant)
  !/
  public:: read_eos_parameters, fix_hyades_state

  !\
  ! Defaultparameters for EOS tables:
  !/
  integer, parameter :: Min_=1, Max_=2, &
                        IndexDefaultEos_I(2) = (/501, 501/)
  integer, public    :: IndexDefaultOpac_I(2)= (/201, 201/)
  real,dimension(Min_:Max_,0:nMaterialMax-1), parameter::&
       TeDefaultEos_II = reshape(&     ! original minimum
       (/1.0e-2, 1.0e+3, & !Xe_     1e-2
         1.0e-3, 2.0e+3, & !Be_     1e-3
         1.0e-3, 1.0e+2, & !Plastic 1e-3
         1.0e-3, 1.0e+2, & !Au_     1e-3
         1.0e-3, 1.0e+2  & !Ay_     1e-3
         /), (/2,5/)),    &
       TeDefaultOpac_II = reshape(&     ! original minimum
       (/3.0e-2, 1.0e+3, & !Xe_     1e-2
         3.0e-2, 2.0e+3, & !Be_     1e-3
         5.0e-2, 1.0e+2, & !Plastic 1e-3
         2.0e-1, 1.0e+2, & !Au_     1e-3
         5.0e-2, 1.0e+2  & !Ay_     1e-3
         /), (/2,5/)),    &
       NaDefault_II = reshape(&
       (/1.0e+24, 1.0e+29, & !Xe_
         1.0e+23, 2.0e+29, & !Be_
         1.0e+24, 1.5e+29, & !Plastic
         1.0e+24, 1.2e+29, & !Au_
         1.0e+24, 1.5e+29  & !Ay_
         /), (/2,5/))
  !Note that at 1 Atm and at the room temperature
  !In gas: N~3.10^{25} m-3
  !In solids: N<10^{29} m-3
  ! electron temperature of 1e-3 eV is approximately 11.6 K
contains

  !============================================================================
  ! This subroutine may be used to exclude undesired columns 
  ! from the eos lookup tables
  ! ================================
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
       case default
          call CON_stop(NameVar_I(iVar)//' is not allowed in the eos lookup tables')
       end select
    end do
  end subroutine read_name_var_eos
  !===============================
  subroutine read_if_use_eos_table(nMaterial)
    !Usage
    !#USEEOSTABLE
    !T                      Use Eos Table for Xe
    !T                      Use Eos Table for Be
    !T                      Use Eos Table for Pl
    !T                      Use Eos Table for Au
    !T                      Use Eos Table for Ay
    use ModReadParam, ONLY: read_var
    !--------------------
    integer,intent(in) :: nMaterial !The number of materials

    integer:: iMaterial
    !------------------!
    do iMaterial = 0, nMaterial-1
       call read_var(&
            'Use EOS table for '//NameMaterial_I(iMaterial),&
            UseEosTable_I(iMaterial))
    end do
       
  end subroutine read_if_use_eos_table
  !===============================
  subroutine read_if_use_opac_table(nMaterial)
    !Usage
    !#USEOPACTABLE
    !T                      Use Opac Table for Xe
    !T                      Use Opac Table for Be
    !T                      Use Opac Table for Pl
    !T                      Use Opac Table for Au
    !T                      Use Opac Table for Ay
    use ModReadParam, ONLY: read_var
    !-------------------
    integer, intent(in) :: nMaterial !The number of materials

    integer:: iMaterial
    !------------------!
    do iMaterial = 0, nMaterial-1
       call read_var(&
            'Use Opac table for '//NameMaterial_I(iMaterial),&
            UseOpacityTable_I(iMaterial))
    end do
       
  end subroutine read_if_use_opac_table
  !============================
  subroutine check_eos_table(iComm,Save,TypeFile)
    use ModLookupTable, ONLY: i_lookup_table, make_lookup_table
    use ModLookupTable, ONLY: init_lookup_table, get_name_description

    integer, optional, intent(in) :: iComm
    logical, optional, intent(in) :: Save 
    character(LEN=*), optional, intent(in)::TypeFile
    integer:: iMaterial, iTable, iPosition
    character(LEN=5)::TypeFileHere
    character(LEN=100)::NameDescription
    !-------------------
    
    do iMaterial = 0,nMaterialMax-1
       if(.not.UseEosTable_I(iMaterial))CYCLE
       iTable =  i_lookup_table(NameMaterial_I(iMaterial)//'_eos')

       if(iTable<0)then
          !Declare the table with the default parameters
          if(.not.present(Save))then
             call init_lookup_table(&
               NameTable = NameMaterial_I(iMaterial)//'_eos', &
               NameCommand = 'make'                         , &
               NameVar = 'logTe logNa '//NameVarEos         , &
               nIndex_I = IndexDefaultEos_I                 , &
               IndexMin_I = (/TeDefaultEos_II(Min_, iMaterial), &
                              NaDefault_II(Min_, iMaterial)/)  , &
               IndexMax_I = (/TeDefaultEos_II(Max_, iMaterial), &
                              NaDefault_II(Max_, iMaterial)/)    )
          else
             TypeFileHere = 'real8'
             if(present(TypeFile))TypeFileHere = TypeFile
             write(NameDescription,'(a,e13.7)')&
                  'CRASH EOS for '//NameMaterial_I(iMaterial)//&
                  'Atomic Mass = ',&
                  cAtomicMassCrash_I(iMaterial)
             call init_lookup_table(&
               NameTable = NameMaterial_I(iMaterial)//'_eos', &
               NameCommand = 'save'                         , &
               NameVar = 'logTe logNa '//NameVarEos         , &
               nIndex_I = IndexDefaultEos_I                 , &
               IndexMin_I = (/TeDefaultEos_II(Min_, iMaterial), &
                              NaDefault_II(Min_, iMaterial)/)  , &
               IndexMax_I = (/TeDefaultEos_II(Max_, iMaterial), &
                              NaDefault_II(max_, iMaterial)/)  , &
               NameFile   = &
                 NameMaterial_I(iMaterial)//'_eos_CRASH.dat', &
               TypeFile = TypeFileHere                      , &
               StringDescription = trim(NameDescription))
          end if
          iTable =  i_lookup_table(NameMaterial_I(iMaterial)//'_eos')
       else
          call get_name_description(iTable, NameDescription)
          if(index(NameDescription,'Mass =')>0)then
             !Read the atomic weight from the table
             iPosition = index(NameDescription,'=')
             read(&
                  NameDescription(iPosition+1:len_trim(NameDescription)),&
                  *)cAtomicMassCRASH_I(iMaterial)
          end if
       end if

       !The table is at least declared. Check if it is filled in
       
       iTableEos4Material_I(iMaterial) = iTable
       iMaterial4EosTable_I(iTable)    = iMaterial
       
       !Temporary disable the use of table in EOS
       UseEosTable_I(iMaterial) = .false.
       
       !Fill in the table
       call make_lookup_table(iTable, calc_eos_table, iComm)
       
       !Recover the true value for UseEos:
       UseEosTable_I(iMaterial) = .true.
    end do
  end subroutine check_eos_table
  !============================
  subroutine calc_eos_table(iTable, Arg1, Arg2, Value_V)
    integer, intent(in):: iTable
    real, intent(in)   :: Arg1, Arg2
    real, intent(out)  :: Value_V(:)
    
    real:: ValueTmp_V(-1:nVarEos)
    real:: Rho, Te
    integer::iMaterial
    !-----------------
    iMaterial = iMaterial4EosTable_I(iTable)
    
    Rho = Arg2 * cAtomicMass * cAtomicMassCRASH_I(iMaterial)
    Te  = Arg1 * cEvToK
    
    ValueTmp_V = 0.0
    call eos_material(iMaterial,Rho,&
       TeIn = Te, &
       eTotalOut = ValueTmp_V(E_)    ,&
       pTotalOut = ValueTmp_V(P_)    ,&
       GammaOut  = ValueTmp_V(Gamma_),&
       CvTotalOut= ValueTmp_V(Cv_)   ,&
       eElectronOut = ValueTmp_V(Ee_),&
       pElectronOut = ValueTmp_V(Pe_),&
       GammaEOut    = ValueTmp_V(GammaE_),&
       CvElectronOut = ValueTmp_V(Cve_)  ,& 
       HeatCond      = ValueTmp_V(Cond_) ,&
       TeTiRelax    = ValueTmp_V(TeTi_)  ,&
       zAverageOut   = ValueTmp_V(Z_)    ,&
       z2AverageOut  = ValueTmp_V(Z2_)   )

    ValueTmp_V(-1:0) = 0.0
    ValueTmp_V(E_)   =  ValueTmp_V(E_)/(cEV * Arg2)
    ValueTmp_V(P_)   =  ValueTmp_V(P_)/(cEV * Arg2)
    ValueTmp_V(Ee_)  =  ValueTmp_V(Ee_)/(cEV * Arg2)
    ValueTmp_V(Pe_)  =  ValueTmp_V(Pe_)/(cEV * Arg2)
    ValueTmp_V(Cv_)  =  ValueTmp_V(Cv_)/(cBoltzmann * Arg2)
    ValueTmp_V(Cve_) =  ValueTmp_V(Cve_)/(cBoltzmann * Arg2)

    Value_V(1:nVarEos) = ValueTmp_V(1:nVarEos)
  end subroutine calc_eos_table
   !============================
  subroutine check_opac_table(iComm,Save,TypeFile)
    use ModLookupTable, ONLY: i_lookup_table, make_lookup_table
    use ModLookupTable, ONLY: init_lookup_table
    
    integer, optional, intent(in) :: iComm
    logical, optional, intent(in) :: Save 
    character(LEN=*), optional, intent(in)::TypeFile
    integer:: iMaterial, iTable
    character(LEN=5)::TypeFileHere
    !-------------------
    nVarOpac = 2 * nGroup

    !Construct NameVarOpac

    NameVarOpac = ''

    if(nGroup==1)then
       NameVarOpac = 'Planck Ross'
    elseif(2 <= nGroup.and.nGroup <= 9)then
       write(NameVarOpac,'(a,i1,a,i1,a)')&
            'Planck(', nGroup, ')  Ross(',nGroup,')'
    elseif(10 <= nGroup.and.nGroup <= 99)then
       write(NameVarOpac,'(a,i2,a,i2,a)')&
            'Planck(', nGroup, ')  Ross(',nGroup,')'
    elseif(100 <= nGroup.and.nGroup <= 999)then
       write(NameVarOpac,'(a,i3,a,i3,a)')&
            'Planck(', nGroup, ')  Ross(',nGroup,')'
    elseif(1000 <= nGroup.and.nGroup <= 9999)then
       write(NameVarOpac,'(a,i4,a,i4,a)')&
            'Planck(', nGroup, ')  Ross(',nGroup,')'
    else
       call CON_stop(&
            'The opacity table cannot be set with the declared nGroup=')
    end if
    do iMaterial = 0,nMaterialMax-1
       if(.not.UseOpacityTable_I(iMaterial))CYCLE
       iTable =  i_lookup_table(NameMaterial_I(iMaterial)//'_opac')

       if(iTable<0)then
          !Declare the table with the default parameters
          if(.not.present(Save))then
             call init_lookup_table(&
               NameTable = NameMaterial_I(iMaterial)//'_opac', &
               NameCommand = 'make'                          , &
               NameVar = 'logRho logTe '//NameVarOpac        , &
               nIndex_I = IndexDefaultOpac_I                 , &
               IndexMin_I = (/NaDefault_II(Min_, iMaterial)*  &
                   cAtomicMass * cAtomicMassCRASH_I(iMaterial),  &
                              TeDefaultOpac_II(Min_, iMaterial)/)  , &
               IndexMax_I = (/NaDefault_II(Max_, iMaterial)*  &
                   cAtomicMass * cAtomicMassCRASH_I(iMaterial),  &
                             TeDefaultOpac_II(Max_, iMaterial)/)    )
          else
             TypeFileHere = 'real8'
             if(present(TypeFile))TypeFileHere = TypeFile
             call init_lookup_table(&
               NameTable = NameMaterial_I(iMaterial)//'_opac', &
               NameCommand = 'save'                          , &
               NameVar = 'logRho logTe '//NameVarOpac        , &
               nIndex_I = IndexDefaultOpac_I                 , &
               IndexMin_I = (/NaDefault_II(Min_, iMaterial)*  &
                   cAtomicMass * cAtomicMassCRASH_I(iMaterial),  &
                              TeDefaultOpac_II(Min_, iMaterial)/)  , &
               IndexMax_I = (/NaDefault_II(Max_, iMaterial)*  &
                   cAtomicMass * cAtomicMassCRASH_I(iMaterial),  &
                             TeDefaultOpac_II(Max_, iMaterial)/), &
               NameFile   = &
                 NameMaterial_I(iMaterial)//'_opac_CRASH.dat', &
               TypeFile = TypeFileHere                      , &
               StringDescription = &
                 'CRASH Opacity for '//NameMaterial_I(iMaterial)  )
          end if
          iTable =  i_lookup_table(NameMaterial_I(iMaterial)//'_opac')
       end if

       !The table is at least declared. Check if it is filled in
       
       iTableOpac4Material_I(iMaterial) = iTable
       iMaterial4OpacTable_I(iTable)    = iMaterial
       
       !Temporary disable the use of table in EOS
       UseOpacityTable_I(iMaterial) = .false.
       
       !Fill in the table
       call make_lookup_table(iTable, calc_opac_table, iComm)
       
       !Recover the true value for UseEos:
       UseOpacityTable_I(iMaterial) = .true.
    end do
  end subroutine check_opac_table
  !============================
  subroutine calc_opac_table(iTable, Arg1, Arg2, Value_V)
    integer, intent(in):: iTable
    real, intent(in)   :: Arg1, Arg2
    real, intent(out)  :: Value_V(:)
    
    real:: Rho, Te, PlanckTmp_I(nGroup), RosselandTmp_I(nGroup)
    integer::iMaterial
    !-----------------
    iMaterial = iMaterial4OpacTable_I(iTable)
    
    Rho = Arg1 
    Te  = Arg2 * cEvToK
    
    Value_V = 0.0
    call eos_material(iMaterial,Rho,&
       TeIn = Te, &
       OpacityPlanckOut_I = PlanckTmp_I, &
       OpacityRosselandOut_I = RosselandTmp_I )
!       OpacityPlanckOut_I = Value_V(1:nGroup),&
!       OpacityRosselandOut_I = Value_V(1+nGroup:2*nGroup) )

    Value_V(1:nGroup) = PlanckTmp_I/Arg1
    Value_V(1+nGroup:2*nGroup) = RosselandTmp_I/Arg1
!    Value_V(:) = Value_V(:)/Arg1
  end subroutine calc_opac_table
  !============================
  !\
  ! Beginning of EOS functions
  !/
  !For two different kinds of input parameters for
  !to characterize the material (either the material number,
  !or the vector of their contents in a mixture) we apply
  !eos_material or eos_mixture. Both subroutine call eos_generic
  !and have very similar set of optional input and output parameters.
  !=================================================================

  subroutine eos_material(iMaterial,Rho,&
       TeIn, eTotalIn, pTotalIn, eElectronIn, pElectronIn,   &
       TeOut, eTotalOut, pTotalOut, GammaOut, CvTotalOut,    &
       eElectronOut, pElectronOut, GammaEOut, CvElectronOut, &
       OpacityPlanckOut_I, OpacityRosselandOut_I,            &
       HeatCond, TeTiRelax, Ne, zAverageOut, z2AverageOut, iError)
    use ModLookupTable, ONLY: interpolate_lookup_table
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
    real,    optional, intent(out) :: zAverageOut  ! <z>
    real,    optional, intent(out) :: z2AverageOut ! <z^2>

    real,    optional, intent(out), &              ! Opacities
                   dimension(nGroup) :: OpacityPlanckOut_I, OpacityRosselandOut_I

    real,    optional, intent(out) :: HeatCond     ! electron heat conductivity (SI)
    real,    optional, intent(out) :: TeTiRelax    ! electron-ion interaction rate (SI)
    integer, optional, intent(out) :: iError       ! error flag

    real   :: Natomic
    real   :: TeEV, Value_V(1:nVarEos), Opacity_V(1:2*nGroup)
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
       iTable = iTableEos4Material_I(iMaterial)
       
       if(present(TeIn))then

          TeEV = TeIn * cKToEV
          call interpolate_lookup_table(iTable, TeEV, Natomic, Value_V, DoExtrapolate=.false.)

       elseif(present(eTotalIn))then

          ! Get an energy per the atomic cell, express in eV
          ! Find temperature from dentity and internal energy
          call interpolate_lookup_table(iTable, E_, eTotalIn/ (cEV * Natomic), Natomic, &
               Value_V, Arg1Out = TeEV, DoExtrapolate=.false.)
       

       elseif(present(pTotalIn))then
          ! Divide pressure by Na , express in eV
          !Find temperature from dentity and pressure
          call interpolate_lookup_table(iTable, P_,  pTotalIn / (cEV * Natomic), Natomic, &
               Value_V, Arg1Out = TeEV, DoExtrapolate=.false.)
          
     
       elseif(present(eElectronIn))then
          ! Get an energy per the atomic cell, express in eV
          ! Find temperature from dentity and internal energy
          call interpolate_lookup_table(iTable, Ee_, eElectronIn/ (cEV * Natomic), Natomic, &
               Value_V, Arg1Out = TeEV, DoExtrapolate=.false.)


       elseif(present(pElectronIn))then

          ! Divide pressure by Na , express in eV
          !Find temperature from dentity and pressure
          call interpolate_lookup_table(iTable, Pe_,  pElectronIn / (cEV * Natomic), Natomic, &
               Value_V, Arg1Out = TeEV, DoExtrapolate=.false.)
          
      
       else
          call CON_stop(NameSub// &
               ': none of Te, eTotal, or pTotal is among the input parameters')
       end if

       if(present(TeOut))      TeOut     = TeEV*cEVToK
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
       if(present(zAverageOut))zAverageOut = Value_V(Z_)
       if(present(z2AverageOut))z2AverageOut = Value_V(Z2_)
       if(present(OpacityPlanckOut_I).or.&
            present(OpacityRosselandOut_I))then
          if(UseOpacityTable_I(iMaterial))then
             iTable = iTableOpac4Material_I(iMaterial)
             call interpolate_lookup_table(iTable, Rho, TeEV, Opacity_V, DoExtrapolate=.false.)
             if(present(OpacityPlanckOut_I))OpacityPlanckOut_I=Opacity_V(1:nGroup) * Rho
             if(present(OpacityRosselandOut_I))&
                  OpacityRosselandOut_I=Opacity_V(nGroup+1:2*nGroup) * Rho
             return
          end if
          !Else: we need to calculate opacities only
       else
          !Opacities are not needed, all the other parameters are calculated
          return
       end if
    end if

    if(nZ_I(iMaterial)<0)then
       call set_mixture(-nZ_I(iMaterial),&
            nZMix_II(1:-nZ_I(iMaterial),iMaterial),&
            cMix_II(1:-nZ_I(iMaterial), iMaterial) )
    else
       call set_element(nZ_I(iMaterial))
       
    end if
    if(UseEosTable_I(iMaterial))then
       !We need only opacities
       call eos_generic(Natomic, &
                  TeIn = TeEV * cEvToK, &
                  OpacityPlanckOut_I = OpacityPlanckOut_I, &
                  OpacityRosselandOut_I = OpacityRosselandOut_I , &
                  iError= iError)
    else
       call eos_generic(Natomic, &
            TeIn, eTotalIn, pTotalIn, eElectronIn, pElectronIn,   &
            TeOut, eTotalOut, pTotalOut, GammaOut, CvTotalOut,    &
            eElectronOut, pElectronOut, GammaEOut, CvElectronOut, & 
            OpacityPlanckOut_I, OpacityRosselandOut_I,            &
            HeatCond, TeTiRelax, Ne, zAverageOut, z2AverageOut, iError)
    end if
  end subroutine eos_material

  !============================================================================
  !Cannot be used for mixed-cell simulations if gold and/or Acrylic is used
  subroutine eos_mixture(RhoToARatio_I,&
       TeIn, eTotalIn, pTotalIn, eElectronIn, pElectronIn,   &
       TeOut, eTotalOut, pTotalOut, GammaOut, CvTotalOut,    &
       eElectronOut, pElectronOut, GammaEOut, CvElectronOut, &
       OpacityPlanckOut_I, OpacityRosselandOut_I,            & 
       HeatCond, TeTiRelax, Ne, zAverageOut, z2AverageOut, iError)
    !\
    !!   WARNING !!!
    !You cannot use total pressure and total energy density as input or output
    !parameters, if the electron temperature is not equal to ion temperature.
    !In this case ONLY electron energy density and electron pressure may be 
    !used.
    !/
    use CRASH_ModPolyimide
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
    real,    optional, intent(out) :: zAverageOut  ! <z>
    real,    optional, intent(out) :: z2AverageOut ! <z^2>

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
         HeatCond, TeTiRelax, Ne, zAverageOut, z2AverageOut,  iError)

  end subroutine eos_mixture

  !============================================================================

  subroutine eos_generic(Natomic, &
       TeIn, eTotalIn, pTotalIn, eElectronIn, pElectronIn,   &
       TeOut, eTotalOut, pTotalOut, GammaOut, CvTotalOut,    &
       eElectronOut, pElectronOut, GammaEOut, CvElectronOut, & 
       OpacityPlanckOut_I, OpacityRosselandOut_I,            &
       HeatCond, TeTiRelax, Ne, zAverageOut, z2AverageOut, iError)
    use CRASH_ModTransport, ONLY: electron_heat_conductivity, te_ti_relaxation
    use CRASH_ModPartition, ONLY: zAv, Z2
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
    real,    optional, intent(out) :: zAverageOut  ! <z>
    real,    optional, intent(out) :: z2AverageOut ! <z^2>

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
    if(present(zAverageOut)) zAverageOut = zAv
    if(present(z2AverageOut))z2AverageOut = Z2
  end subroutine eos_generic
  !====================
  !End of EOS functions
  !=============================
  !\
  ! Miscellaneous subroutnies (probably, redundant)
  !/
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
