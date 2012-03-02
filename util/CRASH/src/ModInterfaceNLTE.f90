!^CFG COPYRIGHT UM
module CRASH_ModInterfaceNLTE
  use CRASH_ModMultiGroup, ONLY:nGroup, EnergyGroup_I
  use CRASH_M_EOS,   ONLY: UseNLTE=>UseCrashEos
  implicit none
  PRIVATE !Except
  public:: UseNLTE  !Logical describing if the nonLTE is used 
  public:: read_nlte,check_nlte !reads UseNlte, makes a check if UseNlte==.true.
  public:: NLTE_EOS !Full list of the eos function parameters (no pIn)
contains
  subroutine read_nlte
    !Use in PARAM.in:
    !#NLTE
    !T/F               UseNLTE
    !
    use ModReadParam,  ONLY: read_var
    !----------------------
    call read_var('UseNLTE',UseNLTE)
  end subroutine read_nlte
  !====================
  subroutine check_nlte
    use CRASH_M_EOS,   ONLY: SetOptions
    use CRASH_M_expTab,ONLY: exp_tab8
    use ModConst,            ONLY: cHPlanckEV
    use CRASH_M_NLTE,only : ng_rad
    use M_RADIOM, only : prep_projE, prepCorrUbar
    logical,save:: DoInit = .true.
    !---------------------
    if(.not.DoInit)return
    DoInit = .false.
    !Initialize NLTE calculations
    call exp_tab8()
    call setoptions(.false., .false., .true.)
    
    !What else?
   
   
    !\
    ! Coefficients for transforming from the user defined grid to
    ! the refined logrithmic-uniform internal fixed grid
    !/ 
    call prep_projE(EnergyGroup_I(0:nGroup),nGroup)

    !\
    ! Initialize and calculate some internal arrays
    !/
    call prepCorrUbar()
   
    ng_rad=nGroup

  end subroutine check_nlte
  !==========================
  subroutine NLTE_EOS(& !Full list of the eos function parameters (no pIn)
       iMaterialIn,Rho,&
       TeIn, eTotalIn, eElectronIn,   &
       EoBIn_I,                                              &
       TeOut, eTotalOut, pTotalOut, GammaOut, CvTotalOut,    &
       eElectronOut, pElectronOut, GammaEOut, CvElectronOut, &
       OpacityPlanckOut_I, OpacityRosselandOut_I,            &
       HeatCond, TeTiRelax, Ne, zAverageOut, z2AverageOut)

    use CRASH_M_EOS,   ONLY: iMaterial, set_kbr
    use CRASH_M_NLTE,only : ng_rad,EoB, NLTE=>NLTE_EOS, setErad 
    use CRASH_ModEos,ONLY: eos, cAtomicMassCRASH_I, &
                           nZMix_II, cMix_II
    use CRASH_M_localProperties,only : atoNum,atoMass
    use ModConst
    ! Eos function for single material

    integer, intent(in):: iMaterialIn     ! index of material
    real,    intent(in):: Rho             ! mass density [kg/m^3]
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
    real,    optional, intent(in)  :: eElectronIn  ! internal energu density of electrons

    ! E-over-B ratio for all known groups
    real,    optional, intent(in)  :: EoBIn_I(1:ng_rad)

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
    
    real:: Tz, NAtomic, Te, EIn, TzSi    !in eV, cm-3, eV, erg/cm3 
    !---------------
    !Set iMaterial and dependent variables

    iMaterial = iMaterialIn

    !Calculate atomic density
    NAtomic = Rho/( cAtomicMassCRASH_I(iMaterial)*cAtomicMass ) & !In 1/m3
         * 1.0e-6                                   ![cm-3] !Convert units

    atomass = cAtomicMassCRASH_I(iMaterial)
   
    atonum  = sum(nZMix_II(:,iMaterial)*&
                      cMix_II(:,iMaterial))
    
    if(present(EoBIn_I))then
       EoB(1:nGroup) = EoBIn_I
    else
       EoB(1:nGroup)=0.0  !Zero radiation energy
    end if

    call set_kbr(NAtom=NAtomic)
    call setErad(eg_o_bg= EoB(1:nGroup),&
                 hnug=EnergyGroup_I(0:nGroup),&          
                 ng=nGroup)

    if(present(TeIn))then
       !Convert temperature to eV
       Te=TeIn * cKToeV
       if(present(EElectronOut).or.present(PElectronOut))then
          !get Tz
          call NLTE(Natom=NAtomic,&
               Te_in=Te,             &
               Zbar_out=zAverageOut, &
               Tz_out=Tz,            &
               Ee_out=EElectronOut,  &
               Pe_out=PElectronOut)
       else
          call NLTE(Natom=NAtomic,   &
               Te_in=Te,             &
               Zbar_out=zAverageOut, &
               Tz_out=Tz,            &
               Et_out=Ein,           &
               Pt_out=PTotalOut)
          !We use a fake parameter for ETotalOut,
          !which is needed in case no energetic output
          !parameters are present
          if(present(ETotalOut))ETotalOut=EIn
       end if
    elseif(present(EElectronIn))then
       !Convert J/m3 = 10^7erg/10^6cm3=10 erg/cm3
       EIn = EElectronIn * 10.0

       !Get Tz
       call NLTE(Natom=NAtomic,&
         Ee_in=EIn,           &
         Zbar_out=zAverageOut,&
         Tz_out=Tz,           &
         Te_out=Te,           &
         Pe_out=PElectronOut)
       
    elseif(present(ETotalIn))then
       !Convert J/m3 = 10^7erg/10^6cm3=10 erg/cm3
        EIn = ETotalIn * 10.0

       !Get Tz

       call NLTE(Natom=NAtomic,&
         Et_in=EIn,           &
         Zbar_out=zAverageOut,&
         Tz_out=Tz,           &
         Te_out=Te,           &
         Pt_out=PTotalOut)
    else
       call CON_stop(&
            'Stop NLTE_eos: TeIn or EElectronIn or ETotalIn must be present')
    end if
    !\
    ! CONVERT
    !/

    !eV to K
    TzSi = Tz*ceVToK
    
    if(present(TeOut))then
       TeOut = Te*ceVToK
    end if
    !erg/cm3=0.1 J/m3
    if(present(EElectronOut))then
       EElectronOut = EElectronOut*0.10
    end if
    if(present(ETotalOut   ))then
       ETotalOut    = ETotalOut   *0.10
    end if
    if(present(PElectronOut))then
       PElectronOut = PElectronOut*0.10
    end if
    if(present(PTotalOut   ))then
       PTotalOut    = PTotalOut   *0.10
    end if

    if(&
         present(GammaOut).or.      &
         present(GammaEOut).or.     &
         present(CvTotalOut).or.    &
         present(CvElectronOut).or. &
         present(OpacityPlanckOut_I).or. &
         present(OpacityRosselandOut_I).or. &
         present(HeatCond).or.      &
         present(TeTiRelax).or.     &
         present(Ne).or.            &
         present(zAverageOut).or.   & 
         present(z2AverageOut) )    &
         call eos(&
         iMaterial=iMaterialIn,     &
         Rho=Rho,                   &
         TeIn=TzSi,                   &
         GammaOut=GammaOut,         &
         CvTotalOut=CvTotalOut,     &
         GammaEOut=GammaEOut,       &
         CvElectronOut=CvElectronOut, &
         OpacityPlanckOut_I=OpacityPlanckOut_I,       &
         OpacityRosselandOut_I=OpacityRosselandOut_I, &
         HeatCond=HeatCond,           &
         TeTiRelax=TeTiRelax,         &
         Ne=Ne,                       &
         zAverageOut=zAverageOut,   &
         z2AverageOut=z2AverageOut)
    !TBD  Correct heat conduction and TeTiRelax with (Te/Tz)factors
  end subroutine NLTE_EOS
end module CRASH_ModInterfaceNLTE
