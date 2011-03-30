program conversion
  use  CRASH_ModAtomicNotation
  use CRASH_ModAtomicMass
  implicit none
  integer,parameter::iUnit = 11, nDensity=201, nTemperature = 201
  integer,parameter::nFrequency=30
  real,  &
       dimension(nFrequency+1) :: hNu_I

  real,  &
       dimension(nTemperature,nDensity) :: &
       zAvr_II, ETotal_II, dETotalOverDT_II, dETotalOverDRho_II, EIon_II, EElectron_II, &
       dEIonOverDT_II, dEElectronOverDT_II, &
       pIon_II, pElectron_II, dPIonOverDT_, dPElectronOverDT_II
  
  real,  &
       dimension(nTemperature,nDensity) :: &
       RossOpacity_II, &
       PlanckOpacity_II, PlanckOpacEms_II
  
  real,  &
       dimension(nTemperature) :: &
       Temperature_I
  real,  &
       dimension(nDensity) :: &
       Density_I     
  
  real,  &
       dimension(nFrequency,nTemperature,nDensity) :: &
       PlanckOpacity_III, RossOpacity_III, PlanckOpacEms_III
  character(len=80)  :: header  
  integer::iString
  character(LEN=6)::NameFile
  character(LEN=2)::NameElement
  integer:: iMaterial
  real::AtomicWeight
 
  !--------------
  write(*,*)'Enter the name of element (e.g. Ar ) - note than Ar.prp should be present'

  read(*,'(a)')NameElement
  if(len_trim(NameElement)==1)NameElement=trim(NameElement)//'_'
  iMaterial = i_material(NameElement)
  if(iMaterial>0)then
     AtomicWeight = cAtomicMass_I(iMaterial)
     write(*,*)'Atomic weight =',AtomicWeight
  end if
     
  NameFile = NameElement//'.prp'
  open(11,file=NameFile,status='old')
  do iString =1,24
     read(11,'(a)')header
  !   write(*,*)header
  end do

     call read_eosopa_file_main ( iUnit, &  
          nDensity, nTemperature, &
          nFrequency, &                                        
          hNu_I, zAvr_II, &
          ETotal_II, dETotalOverDT_II, dETotalOverDRho_II, EIon_II, &
          EElectron_II, dEIonOverDT_II, dEElectronOverDT_II, &
          pIon_II, pElectron_II, dPIonOverDT_, dPElectronOverDT_II, &
          RossOpacity_II, &
          PlanckOpacity_II, &
          PlanckOpacEms_II, &
          Temperature_I, Density_I, &
          PlanckOpacity_III, &
          RossOpacity_III, PlanckOpacEms_III)



  close(11)
end program

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
     pIon_II, pElectron_II, dPIonOverDT_, dPElectronOverDT_II, &
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
       pIon_II, pElectron_II, dPIonOverDT_, dPElectronOverDT_II

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

  character(len=80)  :: header


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
  dPIonOverDT_(:,:) = 0.0
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
           read (iUnit,802) header
      !     write(*,*)header
        enddo
     end if



     ! ...    next, the opacity grid parameters

     if ( DoReadOpacity ) then


        read (iUnit,*) n_temp_opc_tbl
       ! write(*,*) 'n_temp_opc_tbl=',n_temp_opc_tbl
        read (iUnit,*) (Temperature_I(it),it=1,n_temp_opc_tbl)
       ! write(*,*)'Temperature range:',Temperature_I(1),Temperature_I(n_temp_opc_tbl)
        read (iUnit,*) n_dens_opc_tbl
       !  write(*,*)'n_dens_opc_tbl=',n_dens_opc_tbl
        read (iUnit,*) (Density_I(it),it=1,n_dens_opc_tbl)
       ! write(*,*)'Density range:',Density_I(1),Density_I(n_dens_opc_tbl)
        read (iUnit,*) num_freq_grps_local
        write(*,*)'nFrequency=',num_freq_grps_local
        
     end if



     ! ...    radiation transport is assumed to be a multi-group model;
     !        read in the frequency group boundaries (eV)
     !        ---------------------------------------------------------

     if ( DoReadOpacity ) then

        nfgp1 = num_freq_grps_local + 1


        read (iUnit,802) header
       ! write(*,*)header
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
  read (iUnit,802) header
  read (iUnit,*) ((zAvr_II(l,m),  l=1,n_temp_eos_tbl), &
       m=1,n_dens_eos_tbl)
  write(*,*)header
 

  if ( DoReadOpacity ) then
     ! ... added frequency integrated mean opacities. prp format ID = 3 (November 2003 PRW)
     if ( TableInID_ .gt. 2 ) then
        ! Frequency integrated Rosseland group opacity (cm2/g)
  !      write(*,*)'Start to read integral Rosseland opac'
        read (iUnit,802) header
        read (iUnit,*) ((RossOpacity_II(l,m), l=1,n_temp_opc_tbl), &
             m=1,n_dens_opc_tbl)
        write(*,*)header


        ! Frequency integrated Planck group absorption opacity (cm2/g)
 !       write(*,*)'Start to read integral Planck opac'
        read (iUnit,802) header
        read (iUnit,*) ((PlanckOpacity_II(l,m), l=1,n_temp_opc_tbl), &
             m=1,n_dens_opc_tbl)
        write(*,*)header

        ! Frequency integrated Planck group emission opacity (cm2/g)
 !       write(*,*)'Start to read integral Planck emissivity'
        read (iUnit,802) header
        read (iUnit,*) ((PlanckOpacEms_II(l,m), l=1,n_temp_opc_tbl), &
             m=1,n_dens_opc_tbl)      
         write(*,*)header
         
     end if
     ! ... end add

  end if


  ! Total internal energy (J/g)
 ! write(*,*)'Start to read internal energy density'
  read (iUnit,802) header
  read (iUnit,*) ((ETotal_II(l,m), l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
  write(*,*)header

  if ( TableInID_ .le. 5 ) then
   !   write(*,*)'Start to read internal energy density T-derivative'
     ! Total energy T-derivative (J/g/eV)         
     read (iUnit,802) header
     read (iUnit,*) ((dETotalOverDT_II(l,m),l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)header

    ! write(*,*)'Start to read internal energy density scaled rho-derivative'

     ! Scaled total energy rho-derivative (1/eV)
     read (iUnit,802) header
     read (iUnit,*) ((dETotalOverDRho_II(l,m),l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)header

     !write(*,*)'Start to read ion energy density'
     ! Ion internal energy (J/g)
     read (iUnit,802) header
     read (iUnit,*) ((EIon_II(l,m), l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)header

    ! write(*,*)'Start to read electron energy density'
     ! Electron internal energy (J/g)
     read (iUnit,802) header
     read (iUnit,*) ((EElectron_II(l,m), l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)header

     ! write(*,*)'Start to read ion energy density T-derivative'
     ! Ion energy T-derivative (J/g/eV)
     read (iUnit,802) header
     read (iUnit,*) ((dEIonOverDT_II(l,m),l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)header

     ! write(*,*)'Start to read electron energy density T-derivative'
     ! Electron energy T-derivative (J/g/eV)
     read (iUnit,802) header
     read (iUnit,*) ((dEElectronOverDT_II(l,m),l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)header

     !write(*,*)'Start to read ion pressure'
     ! Ion pressure (erg/cm**3)
     read (iUnit,802) header
     read (iUnit,*) ((pIon_II(l,m), l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)header

    ! write(*,*)'Start to read electron pressure'
     ! Electron pressure (erg/cm**3)
     read (iUnit,802) header
     read (iUnit,*) ((pElectron_II(l,m), l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)header

     !write(*,*)'Start to read ion pressure T-derivative'
     ! Ion pressure T-derivative (erg/cm**3/eV)
     read (iUnit,802) header
     read (iUnit,*) ((dPIonOverDT_(l,m),l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)header

     !write(*,*)'Start to read electron pressure T-derivative'
     ! Electron pressure T-derivative (erg/cm**3/eV)
     read (iUnit,802) header
     read (iUnit,*) ((dPElectronOverDT_II(l,m),l=1,n_temp_eos_tbl), m=1,n_dens_eos_tbl)
     write(*,*)header

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
           read (iUnit,802) header
!           write(*,*)'Read Ross Opac:'//header
           read (iUnit,*) (RossOpacity_III(ig,it,id),ig=1,num_freq_grps_local)
           read (iUnit,802) header
!           write(*,*)'Read Planck Opac Emis:'//header
           read (iUnit,*) (PlanckOpacEms_III(ig,it,id),ig=1,num_freq_grps_local)
           read (iUnit,802) header
!           write(*,*)'Read Planck Opac Abs'//header
           read (iUnit,*) (PlanckOpacity_III(ig,it,id),ig=1,num_freq_grps_local)

        enddo
     enddo
!     read (iUnit,802) header
!     write(*,*)header
     
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
