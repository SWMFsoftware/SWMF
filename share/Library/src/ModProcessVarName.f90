! ^CFG COPYRIGHT UM
! ======================================================
module ModProcessVarName

  implicit none

  private
  public :: process_var_name, nVarMax

  integer,parameter :: nVarMax = 100   ! maximum number of state variables
  integer,parameter :: nSubstance = 27 ! number of distinct fluids/species

  ! Number of state variables associated with each substance to be standarized
  integer,parameter :: nVarPerSubstance = 7

  ! Number of allowed alternative names for each variable
  integer  :: nSynonym = 4

  ! State variables not associated with a specific fluid/ specie
  integer,parameter  :: nVarExtra = 7

  ! Named indices for species/fluids
  integer, parameter :: &
       H_    = 1,  &
       Hp_   = 2,  &
       HpSw_ = 3,  &
       H2p_  = 4,  &
       O_    = 5,  &
       Op_   = 6,  &
       O2p_  = 7,  & 
       He_   = 8,  &
       OHp_  = 9,  &
       N_    = 10, &
       COp_  = 11, &
       CO2p_ = 12, &
       H2O_  = 13, &
       H2Op_ = 14, &
       H3Op_ = 15, &
       Mp_   = 16, &
       Lp_   = 17, &
       MHCp_ = 18, &
       HHCp_ = 19, &
       HNIp_ = 20, &
       Sw_   = 21, &
       Iono_ = 22, &
       Neu1_ = 23, &
       Neu2_ = 24, &
       Neu3_ = 25, &
       Neu4_ = 26, &
       Main_ = 27 ! main component, MHD/HD

  ! named indices for basic state variables associated with a substance
  integer,parameter :: &
       Rho_   = 1, &
       RhoUx_ = 2, &
       RhoUy_ = 3, &
       RhoUz_ = 4, &
       p_     = 5, &
       Ppar_  = 6, &
       Energy_= 7

  ! string array containing basic state variables associated with a substance
  character(len = 6) :: NameSubstanceVar_I(nVarPerSubstance) = (/ &
       'Rho   ', &
       'Mx    ', &
       'My    ', &
       'Mz    ', &
       'P     ', &
       'Ppar  ', &
       'Energy'  /)

  ! string arrays containing variables not associated with any substance
  character(len=4) :: NameVarExtra_I(nVarExtra) = (/ &
       'bx  ', &
       'by  ', &
       'bz  ', &
       'pe  ', &
       'ew  ', &
       'eint', &
       'hyp ' /)

 character(len=4) :: NameVarExtraStandardized_I(nVarExtra) = (/ &
       'Bx  ', &
       'By  ', &
       'Bz  ', &
       'Pe  ', &
       'Ew  ', &
       'EInt', &
       'Hyp ' /)

  ! Array storing standarized variable names for all species / fluids
  character(len = 20),allocatable :: SubstanceStandardName_II(:,:)

  ! Array storing all possible names 
  character(len = 20),allocatable :: Dictionary_III(:, :, :)
  ! -------------------------------------------------------------------------
contains

  subroutine process_var_name(nVarName, NameVar_V,  &
       nDensity, nSpeed, nP, nPpar, nWave, nMaterial)

    use ModUtilities,  ONLY: lower_case

    integer,intent(in)                :: nVarName
    character(len=*), intent(inout)   :: NameVar_V(nVarName)
    integer,intent(out)               :: nDensity, nSpeed, nP, nPpar
    integer,intent(out)               :: nWave, nMaterial
   ! DESCRIPTION:
    ! ------------
    ! 1. Creates standard names and a dictionary for each standard name.
    ! The dictionary only contains the basic hydro quantities for 
    ! different substances. Other quantities, e.g. magnetic field, 
    ! that are not associated with a specific substance are
    ! stored separately. This allows us to avoid searching the complete 
    ! dictionary when it is not needed.
    ! The dictionary is a string array:
    !    Dictionary_III(nSubstance, nVarPerSubstance, nSynonym)
    ! where nSubstance is the number of possible species/ fluids,
    !    nVarPerSubstance 
    !              enumarates the variables associated with each substance.
    !    nSynonym  is the number of alternative names representing the same
    !              physical quantity, used by different ModEquation files.
    !
    ! 2. Look up the elements of NameVar_V and replace them with standard names
    ! The look up procedure in the dictionary is done by 
    !    call find_substance_replace_name
    ! Once a specific NameVarIn is found to be identical to a dictionary item:
    ! it is replaced with SubstanceStandardName_II(iSubstance,iVarPerSubstance)
    !
    ! 3. The number of fluids and species found are returned by 
    !      nDensity and nSpeed.
 
    integer                   :: nDistinctSubstanceVar_I(nVarPerSubstance)
    character(len=15)                 :: NameVarIn
    integer                           :: iName, iVar, iSubstanceFound = 0
    logical                           :: IsFoundVar 

    character(len=*), parameter:: NameSub = 'process_var_name'
    ! ------------------------------------------------------------------------
    nDistinctSubstanceVar_I(:) = 0
    nWave = 0
    nMaterial = 0

    ! create standard names and dictionary arrays
    allocate(SubstanceStandardName_II(nSubstance, nVarPerSubstance))
    allocate(Dictionary_III(nSubstance, nVarPerSubstance, nSynonym))
    call create_dictionary

    ! Look up each var name
    NAMELOOP: do iName = 1, nVarName
       ! init search
       IsFoundVar = .false.
       NameVarIn = NameVar_V(iName)
       call lower_case(NameVarIn)

       ! Don't look up in dictionary for: bx, by, bz, EInt, ew, pe, hyp
       do iVar = 1, nVarExtra
          if(NameVarIn == NameVarExtra_I(iVar)) then
             NameVar_V(iName) = NameVarExtraStandardized_I(iVar)
             IsFoundVar = .true.
             CYCLE NAMELOOP
          end if
       end do
    
       ! check dictionary ( loop over density, momentum. pressure, energy)
       do iVar = 1, nVarPerSubstance 
          call find_substance_replace_name
          if(IsFoundVar) then
             ! Count how many distinct substance variables are present
             nDistinctSubstanceVar_I(iVar) = &
                  nDistinctSubstanceVar_I(iVar) +1 
             CYCLE NAMELOOP
          end if
       end do
       
       ! variable name may correspond to numbered wave/material
       ! These names are created  in BATSRUS:MH_set_parameters 
       ! and need not be changed
       if (lge(NameVarIn, 'i01') .and. lle(NameVarIn, 'i99')) then
          nWave = nWave + 1
          IsFoundVar = .true.
          CYCLE NAMELOOP
       end if
       
       if (lge(NameVarIn, 'm01') .and. lle(NameVarIn, 'm99')) then          
          nMaterial = nMaterial + 1
          IsFoundVar = .true.
          CYCLE NAMELOOP
       end if
 
       if(.not. IsFoundVar) then 
          write(*,*) 'ERROR: Var name not in dictionary: ',NameVarIn
          write(*,*) 'Please use standard variable names in ModEuqation '// &
               'and recompile:'
          !write(*,*) SubstanceStandardName_II
          write(*,*) ''
          call CON_stop(NameSub//': unknown variable '//NameVarIn)
       end if

    end do NAMELOOP
   
    nDensity = nDistinctSubstanceVar_I(Rho_)
    nSpeed   = nDistinctSubstanceVar_I(RhoUx_)
    nP       = nDistinctSubstanceVar_I(P_)
    nPpar    = nDistinctSubstanceVar_I(Ppar_)

    deallocate(Dictionary_III)
    deallocate(SubstanceStandardName_II)

  contains
    !==========================================================================
    subroutine find_substance_replace_name

      ! lookup var name in dictionary, replace with standard name

      implicit none

      integer             :: iSubstance, iSynonym
      character(len=15)   :: DictionaryItem
      !----------------------------------------------------------------
      do iSubstance = 1, nSubstance 
         do iSynonym = 2, nSynonym
            DictionaryItem = Dictionary_III(iSubstance, iVar, iSynonym)
            if(len_trim(DictionaryItem) > 0) then
               if(NameVarIn ==  DictionaryItem) then
                  iSubstanceFound = iSubstance
                  IsFoundVar = .true.
                  NameVar_V(iName) = &
                       SubstanceStandardName_II(iSubstanceFound, iVar)
                  RETURN
               end if
            end if
         end do
      end do
    end subroutine find_substance_replace_name

  end subroutine process_var_name
  ! =========================================================================  
  subroutine create_standard_name

    implicit none

    integer   :: iVar, iSubstance
    
    character(len = 6) :: NameSubstance_I(nSubstance)
    ! ---------------------------------------------------------------------

    NameSubstance_I = (/ &
          'H   ',  &
          'Hp  ',  &
          'HpSw',  &
          'H2p ',  &
          'O   ',  &
          'Op  ',  &
          'O2p ',  & 
          'He  ',  &
          'OHp ',  &
          'N   ',  &
          'COp ',  &
          'CO2p',  &
          'H2O ',  &
          'H2Op',  &
          'H3Op',  &
          'Mp  ',  &
          'Lp_ ',  &
          'MHCp',  &
          'HHCp',  &
          'HNIp',  &
          'Sw  ',  &
          'Iono',  &
          'Neu1',  &
          'Neu2',  &
          'Neu3',  &
          'Neu4',  &
          '    '  /) ! main component, MHD / HD 
          
 
    ! loop over all possible species/fluids to fill in Name arrays
    do iSubstance = 1, nSubstance
       do iVar = 1, nVarPerSubstance
          SubstanceStandardName_II(iSubstance,iVar) = &
              ''//trim(NameSubstance_I(iSubstance))//NameSubstanceVar_I(iVar)
       
       end do
    end do
    
  end subroutine create_standard_name
  ! =========================================================================
  subroutine create_dictionary

    Dictionary_III(:,:,:) = ''

    ! first page in dictionary is a 2 by 2 array of standard names
    call create_standard_name
    Dictionary_III(:,:,1) = SubstanceStandardName_II

    !\
    ! fill in alternative names
    !/
    Dictionary_III(Main_, Rho_,      2) = 'rho'
    Dictionary_III(Main_, RhoUx_,    2) = 'mx'
    Dictionary_III(Main_, RhoUx_,    3) = 'rhoux'
    Dictionary_III(Main_, RhoUy_,    2) = 'my'
    Dictionary_III(Main_, RhoUy_,    3) = 'rhouy'
    Dictionary_III(Main_, RhoUz_,    2) = 'mz'
    Dictionary_III(Main_, RhoUz_,    3) = 'rhouz'
    Dictionary_III(Main_, Energy_,   2) = 'e'
    Dictionary_III(Main_, p_,        2) = 'p'
    Dictionary_III(Main_, Ppar_,     2) = 'ppar'

    ! H atoms
    Dictionary_III(H_, Rho_,   2) = 'rhoh'

    ! H+ ions
    Dictionary_III(Hp_, Rho_,   2) = 'hprho'
    Dictionary_III(Hp_, Rho_,   3) = 'h1p'
    Dictionary_III(Hp_, Rho_,   4) = 'hp'
    Dictionary_III(Hp_, RhoUx_, 2) = 'hpmx'
    Dictionary_III(Hp_, RhoUx_, 3) = 'hpux'
    Dictionary_III(Hp_, RhoUy_, 2) = 'hpmy'
    Dictionary_III(Hp_, RhoUy_, 3) = 'hpuy'
    Dictionary_III(Hp_, RhoUz_, 2) = 'hpmz'
    Dictionary_III(Hp_, RhoUz_, 3) = 'hpuz'
    Dictionary_III(Hp_, p_,     2) = 'hpp'
    Dictionary_III(Hp_, Energy_,2) = 'hpe'

    ! H+ ions in solar wind
    Dictionary_III(HpSw_, Rho_,   2) = 'hpswrho'
    Dictionary_III(HpSw_, RhoUx_, 2) = 'hpswmx'
    Dictionary_III(HpSw_, RhoUy_, 2) = 'hpswmy'
    Dictionary_III(HpSw_, RhoUz_, 2) = 'hpswmz'
    Dictionary_III(HpSw_, p_,     2) = 'hpswp'
    Dictionary_III(HpSw_, Energy_,2) = 'hpswe'

    ! H2+ ions
    Dictionary_III(H2p_, Rho_,    2) = 'h2p'

    ! He atoms
    Dictionary_III(He_, Rho_,     2) = 'rhohe'

    ! O atoms
    Dictionary_III(O_, Rho_,      2) = 'rhoo'

    ! O+ ions
    Dictionary_III(Op_, Rho_,   2) = 'oprho'
    Dictionary_III(Op_, Rho_,   3) = 'op'
    Dictionary_III(Op_, RhoUx_, 2) = 'opmx'
    Dictionary_III(Op_, RhoUy_, 2) = 'opmy'
    Dictionary_III(Op_, RhoUz_, 2) = 'opmz'
    Dictionary_III(Op_, p_,     2) = 'opp'
    Dictionary_III(Op_, Energy_,2) = 'ope'

    ! O2+ ions
    Dictionary_III(O2p_, Rho_,   2) = 'o2p'
    Dictionary_III(O2p_, Rho_,   3) = 'o2prho'   
    Dictionary_III(O2p_, RhoUx_, 2) = 'o2pmx'
    Dictionary_III(O2p_, RhoUy_, 2) = 'o2pmy'
    Dictionary_III(O2p_, RhoUz_, 2) = 'o2pmz'
    Dictionary_III(O2p_, p_,     2) = 'o2pp'
    Dictionary_III(O2p_, Energy_,2) = 'o2pe'
   
    ! CO+ ions
    Dictionary_III(COp_, Rho_,   2) = 'cop'
   
    ! CO2+ ions
    Dictionary_III(CO2p_, Rho_,   2) = 'co2prho'
    Dictionary_III(CO2p_, Rho_,   3) = 'co2p'
    Dictionary_III(CO2p_, RhoUx_, 2) = 'co2pmx'
    Dictionary_III(CO2p_, RhoUy_, 2) = 'co2pmy'
    Dictionary_III(CO2p_, RhoUz_, 2) = 'co2pmz'
    Dictionary_III(CO2p_, p_,     2) = 'co2pp'
    Dictionary_III(CO2p_, Energy_,2) = 'co2pe'

    ! H2O molecules
    Dictionary_III(H2O_, Rho_,    2) = 'rhoh2o'
    
    ! H2O+ ions
    Dictionary_III(H2Op_, Rho_,   2) = 'h2op'
    Dictionary_III(H2Op_, Rho_,   3) = 'rhoh2op'

    ! H3O+ ions
    Dictionary_III(H3Op_, Rho_,   2) = 'h3op'

    ! OH+ ions
    Dictionary_III(OHp_, Rho_,   2) = 'ohp'
   
    ! Titan ions
    Dictionary_III(N_,    Rho_,   2) = 'rhon'
    Dictionary_III(Mp_,   Rho_,   2) = 'mp'
    Dictionary_III(Lp_,   Rho_,   2) = 'lp'
    Dictionary_III(MHCp_, Rho_,   2) = 'mhcp'
    Dictionary_III(HHCp_, Rho_,   2) = 'hhcp'
    Dictionary_III(HNIp_, Rho_,   2) = 'hnip'

    ! solar wind
    Dictionary_III(Sw_, Rho_,   2) = 'swrho'
    Dictionary_III(Sw_, Rho_,   3) = 'rhosw'
    Dictionary_III(Sw_, RhoUx_, 2) = 'swmx'
    Dictionary_III(Sw_, RhoUy_, 2) = 'swmy'
    Dictionary_III(Sw_, RhoUz_, 2) = 'swmz'
    Dictionary_III(Sw_, p_,     2) = 'swp'
    Dictionary_III(Sw_, Energy_,2) = 'swe'

    ! ionosphere
    Dictionary_III(Iono_, Rho_,   2) = 'ionorho'
    Dictionary_III(Iono_, Rho_,   3) = 'rhoion'
    Dictionary_III(Iono_, RhoUx_, 2) = 'ionomx'
    Dictionary_III(Iono_, RhoUy_, 2) = 'ionomy'
    Dictionary_III(Iono_, RhoUz_, 2) = 'ionomz'
    Dictionary_III(Iono_, p_,     2) = 'ionop'
    Dictionary_IIi(Iono_, Energy_,2) = 'ionoe'

    ! Outer Heliosphere Pop1 / arbitrary neutral
    Dictionary_III(Neu1_, Rho_,   2) = 'neurho'
    Dictionary_III(Neu1_, RhoUx_, 2) = 'neumx'
    Dictionary_III(Neu1_, RhoUy_, 2) = 'neumy'
    Dictionary_III(Neu1_, RhoUz_, 2) = 'neumz'
    Dictionary_III(Neu1_, p_,     2) = 'neup'
    Dictionary_III(Neu1_, Energy_,2) = 'neue'
   
    ! Outer Heliosphere Pop2 / arbitrary neutral
    Dictionary_III(Neu2_, Rho_,   2) = 'ne2rho'
    Dictionary_III(Neu2_, RhoUx_, 2) = 'ne2mx'
    Dictionary_III(Neu2_, RhoUy_, 2) = 'ne2my'
    Dictionary_III(Neu2_, RhoUz_, 2) = 'ne2mz'
    Dictionary_III(Neu2_, p_,     2) = 'ne2p'
    Dictionary_III(Neu2_, Energy_,2) = 'ne2e'

    ! Outer Heliosphere Pop3 / arbitrary neutral
    Dictionary_III(Neu3_, Rho_,   2) = 'ne3rho'
    Dictionary_III(Neu3_, RhoUx_, 2) = 'ne3mx'
    Dictionary_III(Neu3_, RhoUy_, 2) = 'ne3my'
    Dictionary_III(Neu3_, RhoUz_, 2) = 'ne3mz'
    Dictionary_III(Neu3_, p_,     2) = 'ne3p'
    Dictionary_III(Neu3_, Energy_,2) = 'ne3e'

    ! Outer Heliosphere Pop4 / arbitrary neutral
    Dictionary_III(Neu4_, Rho_,   2) = 'ne4rho'
    Dictionary_III(Neu4_, RhoUx_, 2) = 'ne4mx'
    Dictionary_III(Neu4_, RhoUy_, 2) = 'ne4my'
    Dictionary_III(Neu4_, RhoUz_, 2) = 'ne4mz'
    Dictionary_III(Neu4_, p_,     2) = 'ne4p'
    Dictionary_III(Neu4_, Energy_,2) = 'ne4e'
    
  end subroutine create_dictionary
  
  ! =========================================================================
 
end module ModProcessVarName
