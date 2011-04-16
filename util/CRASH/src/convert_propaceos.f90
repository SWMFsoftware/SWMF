program PROPACEOS
  use  CRASH_ModAtomicNotation
  use CRASH_ModAtomicMass,ONLY: cAtomicMass_I
  use ModConst
  use ModPlotFile
  use CRASH_ModPartition,ONLY: &
       Z2PerA, Z2, ZAv, CoulombLog, Na, Te
  use CRASH_ModAtomicDataMix,ONLY: nMix,Concentration_I, nZ_I
  use CRASH_ModPartition,ONLY: Population_II,iZMin_I
  use CRASH_ModTransport,ONLY: electron_heat_conductivity, te_ti_relaxation
  implicit none

  integer,parameter::iUnit = 11, nDensity=201, nTemperature = 201
  integer,parameter::nFrequency=30
  real,  &
       dimension(nFrequency+1) :: hNu_I

  real,  &
       dimension(nTemperature,nDensity) :: &
       zAvr_II, ETotal_II, dETotalOverDT_II, dETotalOverDRho_II, EIon_II, EElectron_II, &
       dEIonOverDT_II, dEElectronOverDT_II, pTotal_II,&
       pIon_II, pElectron_II, dPIonOverDT_II, dPElectronOverDT_II,&
       dPTotalOverRho_II, dPElectronOverRho_II
  
  real,  &
       dimension(nTemperature,nDensity) :: &
       RossOpacity_II, &
       PlanckOpacity_II, PlanckOpacEms_II
  
  real,  &
       dimension(nTemperature) :: &
       Temperature_I
  real,  &
       dimension(nDensity) :: &
       Density_I, Rho_I     
  
  real,  &
       dimension(nFrequency,nTemperature,nDensity) :: &
       PlanckOpacity_III, RossOpacity_III, PlanckOpacEms_III

  real, allocatable:: Value_VII(:,:,:) 

  character(len=80)  :: StringHeader  
  integer::iString, iRho, iTe
  character(LEN=6)::NameFile
  character(LEN=2)::NameMaterial
  integer:: iMaterial
  real::AtomicMass


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
  character(LEN=100):: NameVarEos = &
       'P E Pe Ee Cv Cve Gamma GammaE TeTi Cond Z Z2'
  
  character(LEN=100)::NameDescription
  !Commented out:
  ! logical,parameter:: DoHeaderOnly=.false., DoCut=.false.
  ! integer,parameter:: iStringStart = -1, iStringLast = -1
 
  !--------------

  Population_II= 0.0; Concentration_I = 0.0
  nZ_I = 0

  write(*,*)'Enter the name of element (e.g. Ar ) - note than Ar.prp should be present'

  read(*,'(a)')NameMaterial
  if(len_trim(NameMaterial)==1)NameMaterial=trim(NameMaterial)//'_'
  iMaterial = i_material(NameMaterial)
  if(iMaterial>0)then
     write(*,*)'iMaterial =', iMaterial
     nZ_I(1) = iMaterial
     AtomicMass = cAtomicMass_I(iMaterial)
     write(*,*)'Atomic weight =',AtomicMass
  else
     call CON_stop('The electron heat conductivity is not supported for mixtures')
  end if
     
  NameFile = NameMaterial//'.prp'
  open(11,file=NameFile,status='old')
  
  !Only prints out the header lines with the numbers. 
  !May be used to find unneeded pieces
  !Commented out!!!!! 
  !if(DoHeaderOnly)then
  !   iString = 0
  !   do 
  !      read(11,'(a)',err=1111,end=1111)StringHeader
  !      iString = iString + 1
  !      if(index(StringHeader,'*')>0&
  !           .and..not.index(StringHeader,'Mean')>0)&
  !           write(*,'(i6,a)')iString,':'//StringHeader
  !
  !   end do
  ! 1111 continue
  !   close (11)
  !   stop
  !end if
  !May be used to delete unneeded piece

  !if(DoHeaderOnly)then
  !   open(12, file = NameMaterial//'.new',status='replace') 
  !   iString = 0
  !   do 
  !      read(11,'(a)',err=2222,end=2222)StringHeader
  !      iString = iString + 1
  !      if(iString<iStringStart.or.iString>iStringLast)&
  !           write(12,'(a)')trim(StringHeader)
  !
  !   end do
!!!!!2222 continue
  !   close (11)
  !   close (12)
  !   stop
  !end if


  do iString =1,24
     read(11,'(a)')StringHeader
     write(*,'(a)')StringHeader
  end do

     call read_eosopa_file_main ( iUnit, &  
          nDensity, nTemperature, &
          nFrequency, &                                        
          hNu_I, zAvr_II, &
          ETotal_II, dETotalOverDT_II, dETotalOverDRho_II, EIon_II, &
          EElectron_II, dEIonOverDT_II, dEElectronOverDT_II, &
          pIon_II, pElectron_II, dPIonOverDT_II, dPElectronOverDT_II, &
          RossOpacity_II, &
          PlanckOpacity_II, &
          PlanckOpacEms_II, &
          Temperature_I, Density_I, &
          PlanckOpacity_III, &
          RossOpacity_III, PlanckOpacEms_III)

  close(11)
  !Rescale Density array

  Density_I = 1.0e6 * Density_I  ! 1/cm3 = 10+6 1/m3
  Rho_I     = cAtomicMass * AtomicMass * Density_I  !Comes in kg/m3 

  !============================Opacities===============================
  !Save opacities
  !Swap indexes (Rho<=>Te), merge to a sigle table
  !Convert units: cm2/g=0.10 m2/kg
  
  allocate(Value_VII(2*nFrequency,nDensity,nTemperature))
  do iTe = 1, nTemperature
     do iRho = 1, nDensity
        Value_VII(1:nFrequency, iRho, iTe) = PlanckOpacity_III(:,iTe,iRho) * 0.10
        Value_VII(1+nFrequency:2*nFrequency, iRho, iTe) = RossOpacity_III(:,iTe,iRho) * 0.10
     end do
  end do

  call save_plot_file( &
         NameMaterial//'_opac_PRISM.dat',                                &
         TypeFileIn     = 'real8',                     &
         StringHeaderIn = 'PROPACEOS Opacity for'//NameMaterial, &
         NameVarIn      = 'logRho logTe Planck01 Planck02 Planck03 Planck04 Planck05'//&
                                      ' Planck06 Planck07 Planck08 Planck09 Planck10'//&
                                      ' Planck11 Planck12 Planck13 Planck14 Planck15'//&
                                      ' Planck16 Planck17 Planck18 Planck19 Planck20'//&
                                      ' Planck21 Planck22 Planck23 Planck24 Planck25'//&
                                      ' Planck26 Planck27 Planck28 Planck29 Planck30'//&
                                      ' Ross01 Ross02 Ross03 Ross04 Ross05'//&
                                      ' Ross06 Ross07 Ross08 Ross09 Ross10'//&
                                      ' Ross11 Ross12 Ross13 Ross14 Ross15'//&
                                      ' Ross16 Ross17 Ross18 Ross19 Ross20'//&
                                      ' Ross21 Ross22 Ross23 Ross24 Ross25'//&
                                      ' Ross26 Ross27 Ross28 Ross29 Ross30', &
         CoordMinIn_D   = (/Rho_I(1), Temperature_I(1)/),             &                             
         CoordMaxIn_D   = (/Rho_I(nDensity), Temperature_I(nTemperature)/),             &
         VarIn_VII      = Value_VII)
  deallocate(Value_VII)
  
  !=========================EOS====================
  !Convert energy in J/g to J/kg, 1/g=10+3 1/kg
  
  ETotal_II = 1.0e3* ETotal_II
  EElectron_II = 1.0e3* EElectron_II

  !Convert from J/kg to J/m3
  do iRho = 1, nDensity
     ETotal_II(:,iRho) = ETotal_II(:,iRho) * Rho_I(iRHo)
     EElectron_II(:,iRho) = EElectron_II(:,iRho) * Rho_I(iRHo)
  end do

  !Convert pressure from dyne/cm2 to Pa: dyne/cm2=0.10 * Pa:
  pElectron_II = pElectron_II *0.10
  pIon_II      = pIon_II      *0.10


  pTotal_II =  pIon_II +  pElectron_II


  !Convert specific heat from J/g/eV to J/m^3/K
  !Convert energy in J/g to J/kg, 1/g=10+3 1/kg
  !1/eV = 1/ceVToK 
  dETotalOverDT_II    = 1.0e3 * dETotalOverDT_II   /cEVToK
  dEElectronOverDT_II = 1.0e3 * dEElectronOverDT_II/cEVToK

  !Convert from J/kg to J/m3
  do iRho = 1, nDensity
     dETotalOverDT_II(:,iRho) = dETotalOverDT_II(:,iRho) * Rho_I(iRHo)
     dEElectronOverDT_II(:,iRho) = dEElectronOverDT_II(:,iRho) * Rho_I(iRHo)
  end do

  !Pressure derivatives
  dPTotalOverRho_II(:,1) = Rho_I(1)/ pTotal_II(:,1) * &
       (pTotal_II(:,2) - pTotal_II(:,1))/(Rho_I(2)-Rho_I(1))

  dPElectronOverRho_II(:,1) = Rho_I(1)/ pElectron_II(:,1) * &
       (pElectron_II(:,2) - pElectron_II(:,1))/(Rho_I(2)-Rho_I(1))

  do iRho = 2, nDensity -1
     dPTotalOverRho_II(:,iRho) = Rho_I(iRho)/ pTotal_II(:,iRho) * &
          (pTotal_II(:,iRho+1) - pTotal_II(:,iRho-1))/&
          (Rho_I(iRho+1)-Rho_I(iRho -1))

     dPElectronOverRho_II(:,iRho) = Rho_I(iRho)/ pElectron_II(:,iRho) * &
          (pElectron_II(:,iRho+1) - pElectron_II(:,iRho-1))/&
          (Rho_I(iRho+1)-Rho_I(iRho -1))
  end do

  dPTotalOverRho_II(:,nDensity) = Rho_I(nDensity)/ pTotal_II(:,nDensity) * &
       (pTotal_II(:,nDensity) - pTotal_II(:,nDensity-1))/&
       (Rho_I(nDensity)-Rho_I(nDensity-1))

  dPElectronOverRho_II(:,nDensity) = Rho_I(nDensity)/ pElectron_II(:,nDensity) * &
       (pElectron_II(:,nDensity) - pElectron_II(:,nDensity-1))/&
       (Rho_I(nDensity)-Rho_I(nDensity-1))

  !dyne/cm2 /eV = 0.1/ceVToK 1/K
  dPIonOverDT_II = dPIonOverDT_II * 0.10/cEVToK
  dPElectronOverDT_II = dPElectronOverDT_II * 0.10/cEVToK


  allocate( Value_VII(nVarEos, nTemperature, nDensity) ) 
  Value_VII = 0.0
  do iRho = 1, nDensity
     do iTe = 1, nTemperature
        Value_VII(E_, iTe, iRho) = ETotal_II( iTe, iRho)/(cEV * Density_I(iRho))

        Value_VII(P_, iTe, iRho) = pTotal_II( iTe, iRho)/(cEV * Density_I(iRho))

        Value_VII(Ee_, iTe, iRho) = ETotal_II( iTe, iRho)/(cEV * Density_I(iRho))

        Value_VII(Pe_, iTe, iRho) = pElectron_II( iTe, iRho)/(cEV * Density_I(iRho))

        Value_VII(Cv_, iTe, iRho) = dETotalOverDT_II(iTe, iRho)/(cBoltzmann * Density_I(iRho))
        Value_VII(Cve_, iTe, iRho) = dEElectronOverDT_II(iTe, iRho)/(cBoltzmann * Density_I(iRho))

        Value_VII(Z_,  iTe, iRho) = zAvr_II(iTe, iRho)
        Value_VII(Z2_,  iTe, iRho) = zAvr_II(iTe, iRho)**2
        ZAv =  Value_VII(Z2_,  iTe, iRho)
        Z2  =  ZAv **2
        Z2PerA = Z2/AtomicMass
        Na     = Density_I(iRho)
        Te     = Temperature_I(iTe)
        nMix = 1
        iZMin_I = 0
        Concentration_I = 0.0; Concentration_I(1) = 1.0

        !Concentration of neutrals is switched off at 2 eV 
        Population_II(0,1) = max(0.0, min(1.0,2.0 - Te))
        Value_VII(Cond_,  iTe, iRho) = electron_heat_conductivity()
        Value_VII(TeTi_,  iTe, iRho) = te_ti_relaxation()

        !Use formula (24) from HEDP.pdf:
        
         Value_VII(GammaE_,  iTe, iRho) = &
              dPElectronOverRho_II(iTe, iRho) + &
              dPElectronOverDT_II( iTe, iRho)**2&
              *Temperature_I(iTe) * cEVToK /(pElectron_II(iTe, iRho)*&
                                    dEElectronOverDT_II(iTe, iRho))
         
         Value_VII(Gamma_,  iTe, iRho) = &
              dPTotalOverRho_II(iTe, iRho) + &
              (dPElectronOverDT_II( iTe, iRho) +&
               dPIonOverDT_II(iTe, iRho))**2&
              *Temperature_I(iTe) * cEVToK /(pTotal_II(iTe, iRho)*&
                                    dETotalOverDT_II(iTe, iRho))

        ! GammaSe = compressibility_at_const_te_e() + &
        !    d_pressure_e_over_d_te()**2 * (Na * Te * cEV) /(heat_capacity_e() * pressure_e())
        ! GammaS =  compressibility_at_const_te  () + &
        !    d_pressure_over_d_te  ()**2 * (Na * Te * cEV) /(heat_capacity  () * pressure  ())
        ! Where:
        !  compressibility_at_const_t = 
        !(Na /P) \left(\partial P / \partial Na \right)_{T=const}
        !   the temperuture derivative of the specific pressure (P/Na)
        ! (1/Na)*(\partial P/\partial T)_{\rho=const} 
        ! Heat capacity - dimensionless, per atom
     end do
  end do

  write(NameDescription,'(a,e13.7)')&
       'PROPACEOS EOS for '//NameMaterial//&
       'Atomic Mass = ',AtomicMass
  call save_plot_file( &
       NameMaterial//'_eos_PRISM.dat',                                &
       TypeFileIn     = 'real8',                     &
       StringHeaderIn = trim(NameDEscription), &
       NameVarIn      = 'logTe logNa '//NameVarEos, &
       CoordMinIn_D   = (/Temperature_I(1),Density_I(1)/),             &                             
       CoordMaxIn_D   = (/Temperature_I(nTemperature),Density_I(nDensity)/),             &
       VarIn_VII      = Value_VII)
  open(11,file='PARAM.'//NameMaterial//'.PRISM',status='replace')
  write(11,'(a)')'----------------EOS TABLE-------------'
  write(11,'(a)')'    '
  write(11,'(a)')'#LOOKUPTABLE'
  write(11,'(a)')NameMaterial//'_eos                     NameTable'
  write(11,'(a)')'load                  NameCommand'
  write(11,'(a)')'Tables/'//NameMaterial//'_eos_PRISM.dat'
  write(11,'(a)')'real8                 TypeFile'
  write(11,'(a)')'EOS(Te, Na) for '//NameMaterial
  write(11,'(a)')'logTe logNa P E Pe Ee Cv Cve Gamma GammaE TeTi Cond Z Z2'
  write(11,'(a)')'201                   nIndex1'
  write(11,'(e13.7,a)')Temperature_I(1),'                       Index1Min (eV)'
  write(11,'(e13.7,a)')Temperature_I(nTemperature),&
                                        '                       Index1Max (eV)'
  write(11,'(a)')'201                   nIndex2'
  write(11,'(e13.7,a)')Density_I(1),    '                       Index2Min (m-3)'
  write(11,'(e13.7,a)')Density_I(nDensity),&
                                        '                       Index2Max (m-3)'
  write(11,'(a)')'  '
  write(11,'(a)')'----------------OPACITY TABLE----------'
  write(11,'(a)')'  '
  write(11,'(a)')'#LOOKUPTABLE'
  write(11,'(a)')NameMaterial//'_opac                   NameTable'
  write(11,'(a)')'load                  NameCommand'
  write(11,'(a)')'Tables/'//NameMaterial//'_opac_PRISM.dat'
  write(11,'(a)')'real8                 TypeFile'
  write(11,'(a)')'Opacity(rho,Te) for '//NameMaterial
  write(11,'(a)')'logrho logT Planck(30) Ross(30)'
  write(11,'(a)')'201                   nIndex1'
  write(11,'(e13.7,a)')Rho_I(1),        '                       Index1Min (kg/m3)'
  write(11,'(e13.7,a)')Rho_I(nDensity), '                       Index1Max (kg/m3)'
  write(11,'(a)')'201                   nIndex2'
  write(11,'(e13.7,a)')Temperature_I(1),'                       Index2Min (eV)'
  write(11,'(e13.7,a)')Temperature_I(nTemperature),&
                                        '                       Index2Max (eV)'
  write(11,'(a)')'  '
  write(11,'(a)')'#END'
  write(11,'(a)')'  '
  close(11)



end program PROPACEOS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                                  +
! + Copyright (c) Prism Computational Sciences, Inc., 1998-2000.     +
! +                                                                  +
! + File:         read.f90                          +
! +                                                                  +
! + Written by:   J. J. MacFarlane                                   +
! +                                                                  +
! + Created:      7/1/99                                             +
! + Modified:                                                        +
! + Modified:                                                        +
! +                                                                  +
! + Description:  Reads in equation of state and opacity data.       +
! +               The opacity data file contains:                    +
! +               electron and ion internal energies and derivatives,+
! +               electron and ion pressure and derivatives,         +
! +               mean charge states, and                            +
! +               multigroup opacities.                              +
! +                                                                  +
! + Note: All real variables and arrays in output other than hnu_in  +
! +       are small_real_kind                                        +
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine read_eosopa_file_main ( iUnit, &  
     nDensity, nTemperature, &
     nFrequency, &                                        
     hNu_I, zAvr_II, &
     ETotal_II, dETotalOverDT_II, dETotalOverDRho_II, EIon_II, &
     EElectron_II, dEIonOverDT_II, dEElectronOverDT_II, &
     pIon_II, pElectron_II, dPIonOverDT_II, dPElectronOverDT_II, &
     RossOpacity_II, &
     PlanckOpacity_II, &
     PlanckOpacEms_II, &
     Temperature_I, Density_I, &
     PlanckOpacity_III, &
     RossOpacity_III, PlanckOpacEms_III)


  ! ... use statements


  ! ... implicit statements
  implicit none

  !
  ! ... input variables:
  integer, &
       intent(in) :: iUnit
  integer, &
       intent(in) :: nDensity
  integer, &
       intent(in) :: nTemperature
  integer, &
       intent(in) :: nFrequency     
 

  !
  ! ... output variables

  real, intent(out), &
       dimension(nFrequency+1) :: hNu_I

  real, intent(out), &
       dimension(nTemperature,nDensity) :: &
       zAvr_II, ETotal_II, dETotalOverDT_II, dETotalOverDRho_II, EIon_II, EElectron_II, &
       dEIonOverDT_II, dEElectronOverDT_II, &
       pIon_II, pElectron_II, dPIonOverDT_II, dPElectronOverDT_II

  real, intent(out), &
       dimension(nTemperature,nDensity) :: &
       RossOpacity_II, &
       PlanckOpacity_II, PlanckOpacEms_II

  real, intent(out), &
       dimension(nTemperature) :: &
       Temperature_I
  real, intent(out), &
       dimension(nDensity) :: &
       Density_I     

  real, intent(out), &
       dimension(nFrequency,nTemperature,nDensity) :: &
       PlanckOpacity_III, RossOpacity_III, PlanckOpacEms_III


  ! ... type declaration statements:

  integer :: num_freq_grps_local
  integer,parameter::TableInID_=3

  character(len=80)  :: StringHeader


  integer :: i, nfgp1, l, it, id, ig, m, &
       n_temp_opc_tbl, n_dens_opc_tbl, &
       n_temp_eos_tbl, n_dens_eos_tbl

  real :: rho0es
  
  logical,parameter:: DoReadOpacity = .true.

  ! *******************************************************************
  !
  !                          begin execution
  !
  ! *******************************************************************

  ! ... initializations

  hNu_I(:) = 0.0

  zAvr_II(:,:) = 0.0
  ETotal_II(:,:) = 0.0
  dETotalOverDT_II(:,:) = 0.0
  dETotalOverDRho_II(:,:) = 0.0
  EIon_II(:,:) = 0.0
  EElectron_II(:,:) = 0.0
  dEIonOverDT_II(:,:) = 0.0
  dEElectronOverDT_II(:,:) = 0.0
  pIon_II(:,:) = 0.0
  pElectron_II(:,:) = 0.0
  dPIonOverDT_II(:,:) = 0.0
  dPElectronOverDT_II(:,:) = 0.0
  RossOpacity_II(:,:) = 0.0
  PlanckOpacity_II(:,:) = 0.0
  PlanckOpacEms_II(:,:) = 0.0
  Temperature_I(:) = 0.0
  Density_I(:) = 0.0
  

  if ( ( TableInID_ .eq. 1 ) .or. & 
       ( TableInID_ .eq. 2 ) .or. & 
       ( TableInID_ .eq. 3 ) .or. &
       ( TableInID_ .eq. 4 ) .or. &
       ( TableInID_ .eq. 5 ) .or. &
       ( TableInID_ .eq. 6 ) .or. &
       ( TableInID_ .eq. 7 ) ) then

     ! ...    format is already checked in get_opacity_table_grid


     ! ...    read temperature and density grid parameters,
     !        and number of frequency groups
     !        ---------------------------------------------

     ! ...    first, the EOS grid parameters



     read (iUnit,*) n_temp_eos_tbl
     write(*,*)'n_temp_eos_tbl=',n_temp_eos_tbl
     
     read (iUnit,*) (Temperature_I(it),it=1,n_temp_eos_tbl)
     write(*,*)'Temperature range:',Temperature_I(1),Temperature_I(n_temp_eos_tbl)
     read (iUnit,*) n_dens_eos_tbl
     write(*,*)'n_dens_eos_tbl=',n_dens_eos_tbl
     read (iUnit,*) (Density_I(it),it=1,n_dens_eos_tbl)
      write(*,*)'Density range:',Density_I(1),Density_I(n_dens_eos_tbl)
     read (iUnit,*) rho0es
     !write(*,*)rho0es
     



     if ( DoReadOpacity ) then
        do i=1,4
           read (iUnit,802) StringHeader
           write(*,*)StringHeader
        enddo
     end if



     ! ...    next, the opacity grid parameters

     if ( DoReadOpacity ) then


        read (iUnit,*) n_temp_opc_tbl
        write(*,*) 'n_temp_opc_tbl=',n_temp_opc_tbl
        read (iUnit,*) (Temperature_I(it),it=1,n_temp_opc_tbl)
        write(*,*)'Temperature range:',Temperature_I(1),Temperature_I(n_temp_opc_tbl)
        read (iUnit,*) n_dens_opc_tbl
        write(*,*)'n_dens_opc_tbl=',n_dens_opc_tbl
        read (iUnit,*) (Density_I(it),it=1,n_dens_opc_tbl)
        write(*,*)'Density range:',Density_I(1),Density_I(n_dens_opc_tbl)
        read (iUnit,*) num_freq_grps_local
        write(*,*)'nFrequency=',num_freq_grps_local
        
     end if



     ! ...    radiation transport is assumed to be a multi-group model;
     !        read in the frequency group boundaries (eV)
     !        ---------------------------------------------------------

     if ( DoReadOpacity ) then

        nfgp1 = num_freq_grps_local + 1


        read (iUnit,802) StringHeader
        write(*,*)StringHeader
        read (iUnit,*) (hNu_I(l),l=1,nfgp1)
        write(*,*)'Frequency range:', hNu_I(1),hNu_I(nfgp1)
       
     end if

  else


     return

  endif


  ! ... read in the charge states, ion and electron specific energies,
  !     ion and electron pressures, and their derivatives wrt temperature;
  !     store in log form
  !     ------------------------------------------------------------------

 

  ! Zbar (esu)
  !write(*,*)'Start to read zAvr'
  read (iUnit,802) StringHeader
  read (iUnit,*) ((zAvr_II(l,m),  l=1,n_temp_eos_tbl), &
       m=1,n_dens_eos_tbl)
  write(*,*)StringHeader
 

  if ( DoReadOpacity ) then
     ! ... added frequency integrated mean opacities. prp format ID = 3 (November 2003 PRW)
     if ( TableInID_ .gt. 2 ) then
        ! Frequency integrated Rosseland group opacity (cm2/g)
   !     write(*,*)'Start to read integral Rosseland opac'
        read (iUnit,802) StringHeader
        read (iUnit,*) ((RossOpacity_II(l,m), l=1,n_temp_opc_tbl), &
             m=1,n_dens_opc_tbl)
        write(*,*)StringHeader


        ! Frequency integrated Planck group absorption opacity (cm2/g)
        write(*,*)'Start to read integral Planck opac'
        read (iUnit,802) StringHeader
        read (iUnit,*) ((PlanckOpacity_II(l,m), l=1,n_temp_opc_tbl), &
             m=1,n_dens_opc_tbl)
        write(*,*)StringHeader

        ! Frequency integrated Planck group emission opacity (cm2/g)
        write(*,*)'Start to read integral Planck emissivity'
        read (iUnit,802) StringHeader
        read (iUnit,*) ((PlanckOpacEms_II(l,m), l=1,n_temp_opc_tbl), &
             m=1,n_dens_opc_tbl)      
        write(*,*)StringHeader
         
     end if
     ! ... end add

  end if


  ! Total internal energy (J/g)
  !write(*,*)'Start to read internal energy density'
  read (iUnit,802) StringHeader
  read (iUnit,*) ((ETotal_II(l,m), l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
  write(*,*)StringHeader

  if ( TableInID_ .le. 5 ) then
     ! write(*,*)'Start to read internal energy density T-derivative'
     ! Total energy T-derivative (J/g/eV)         
     read (iUnit,802) StringHeader
     read (iUnit,*) ((dETotalOverDT_II(l,m),l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)StringHeader

    ! write(*,*)'Start to read internal energy density scaled rho-derivative'

     ! Scaled total energy rho-derivative (1/eV)
     read (iUnit,802) StringHeader
     read (iUnit,*) ((dETotalOverDRho_II(l,m),l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)StringHeader

     !write(*,*)'Start to read ion energy density'
     ! Ion internal energy (J/g)
     read (iUnit,802) StringHeader
     read (iUnit,*) ((EIon_II(l,m), l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)StringHeader

    ! write(*,*)'Start to read electron energy density'
     ! Electron internal energy (J/g)
     read (iUnit,802) StringHeader
     read (iUnit,*) ((EElectron_II(l,m), l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)StringHeader

     ! write(*,*)'Start to read ion energy density T-derivative'
     ! Ion energy T-derivative (J/g/eV)
     read (iUnit,802) StringHeader
     read (iUnit,*) ((dEIonOverDT_II(l,m),l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)StringHeader

     ! write(*,*)'Start to read electron energy density T-derivative'
     ! Electron energy T-derivative (J/g/eV)
     read (iUnit,802) StringHeader
     read (iUnit,*) ((dEElectronOverDT_II(l,m),l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)StringHeader

     !write(*,*)'Start to read ion pressure'
     ! Ion pressure (erg/cm**3)
     read (iUnit,802) StringHeader
     read (iUnit,*) ((pIon_II(l,m), l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)StringHeader

    ! write(*,*)'Start to read electron pressure'
     ! Electron pressure (erg/cm**3)
     read (iUnit,802) StringHeader
     read (iUnit,*) ((pElectron_II(l,m), l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)StringHeader

     !write(*,*)'Start to read ion pressure T-derivative'
     ! Ion pressure T-derivative (erg/cm**3/eV)
     read (iUnit,802) StringHeader
     read (iUnit,*) ((dPIonOverDT_II(l,m),l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)StringHeader

     !write(*,*)'Start to read electron pressure T-derivative'
     ! Electron pressure T-derivative (erg/cm**3/eV)
     read (iUnit,802) StringHeader
     read (iUnit,*) ((dPElectronOverDT_II(l,m),l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)StringHeader

  endif



  !     ----------------------------
  ! ... now the multigroup opacities
  !     ----------------------------


  RossOpacity_III       = 0.
  PlanckOpacity_III = 0.
  PlanckOpacEms_III= 0.

  ! ... read in the rosseland, planck emission, and planck absorption tables;
  !     values are in (cm**2/g)

  if ( DoReadOpacity ) then



     do id=1, n_dens_opc_tbl
        do it=1, n_temp_opc_tbl
           read (iUnit,802) StringHeader
!           write(*,*)'Read Ross Opac:'//StringHeader
           read (iUnit,*) (RossOpacity_III(ig,it,id),ig=1,num_freq_grps_local)
           read (iUnit,802) StringHeader
!           write(*,*)'Read Planck Opac Emis:'//StringHeader
           read (iUnit,*) (PlanckOpacEms_III(ig,it,id),ig=1,num_freq_grps_local)
           read (iUnit,802) StringHeader
!           write(*,*)'Read Planck Opac Abs'//StringHeader
           read (iUnit,*) (PlanckOpacity_III(ig,it,id),ig=1,num_freq_grps_local)

        enddo
     enddo
!     read (iUnit,802) StringHeader
!     write(*,*)StringHeader
     
  end if
  return

  ! ... format statement to read the header
802 format (a80)

end subroutine read_eosopa_file_main
!============================================================================
! The following subroutines are here so that we can use SWMF library routines
! Also some features available in SWMF mode only require empty subroutines
! for compilation of the stand alone code.
!============================================================================
subroutine CON_stop(StringError)
  implicit none
  character (len=*), intent(in) :: StringError
  write(*,*)StringError
  stop
end subroutine CON_stop
!============================================================================
subroutine CON_set_do_test(String,DoTest,DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe
end subroutine CON_set_do_test
