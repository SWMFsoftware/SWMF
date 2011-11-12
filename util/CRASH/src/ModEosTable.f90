!^CFG COPYRIGHT UM
module CRASH_ModEosTable
  use CRASH_ModEos
  use ModLookupTable, ONLY: MaxTable
  use CRASH_ModMultiGroup, ONLY: nGroup, &
      OpacityPlanck_I, OpacityRosseland_I
  use ModConst, ONLY: cAtomicMass, cEvToK, cEv, cBoltzmann

  implicit none

  !The following subroutine:
  !1. Checks if the tables are available for these 
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

  ! Arrays which relate the iTable for the EOS table with 
  ! the material number:
  integer:: iMaterial4EosTable_I(MaxTable) = -1

  ! Arrays which relate the iTable for the opacity table with 
  ! the material number:
  integer:: iMaterial4OpacTable_I(MaxTable) = -1

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

  !\
  ! Defaultparameters for EOS tables:
  !/
  integer, parameter :: Min_=1, Max_=2, &
                        IndexDefaultEos_I(2) = (/201, 201/)
  integer, public    :: IndexDefaultOpac_I(2)= (/201, 201/)
  real,dimension(Min_:Max_,0:nMaterialMax-1), parameter::&
       TeDefaultEos_II = reshape(&     ! original minimum
       (/1.0e-2, 1.0e+3, & !Xe_     1e-2
         1.0e-3, 4.0e+3, & !Be_     1e-3
         1.0e-3, 1.0e+2, & !Plastic 1e-3
         1.0e-3, 1.0e+2, & !Au_     1e-3
         1.0e-3, 1.0e+2  & !Ay_     1e-3
         /), (/2,5/)),    &
       TeDefaultOpac_II = reshape(&     ! original minimum
       (/3.0e-2, 1.0e+3, & !Xe_     1e-2
         3.0e-2, 4.0e+3, & !Be_     1e-3
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
                  ' Atomic Mass = ',&
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
    call eos(iMaterial,Rho,&
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
    call eos(iMaterial,Rho,&
       TeIn = Te, &
       OpacityPlanckOut_I = PlanckTmp_I, &
       OpacityRosselandOut_I = RosselandTmp_I )


    Value_V(1:nGroup) = PlanckTmp_I/Arg1
    Value_V(1+nGroup:2*nGroup) = RosselandTmp_I/Arg1

  end subroutine calc_opac_table
  !============================

end module CRASH_ModEosTable
