Module ModReadArtep
  use ModIoUnit, ONLY: io_unit_new
  implicit none
  SAVE
  real, allocatable:: Value_VII(:,:,:) 
  integer::line
contains
  subroutine read_artep
    integer:: iFile, iRho, iTe
    real:: Aux_V(7), Aux
    !---------------

    allocate(Value_VII(60,201,201))

    iFile = io_unit_new()
    open(iFile, file='xz_opac.artep', form='formatted',status='old')
    line=0
    do iTe = 1, 201
       do iRho = 1,201
          line = line+1;write(*,*)line
          read(iFile,*)Aux_V,Value_VII(1:30, iRho, iTe), Aux, Value_VII(31:60, iRho, iTe)
       end do
    end do
    close(iFile)
  end subroutine read_artep
end Module ModReadArtep
program save_eos_table
  use CRASH_ModStatSum
  use CRASH_ModPartition, ONLY: UseCoulombCorrection
  use CRASH_ModPolyimide
  use CRASH_ModFermiGas
  use CRASH_ModMultiGroup
  use CRASH_ModExcitation
  use CRASH_ModAtomicDataMix
  use CRASH_ModExcitationData,ONLY : n_ground, cExcitationN_III
  use CRASH_ModIonMix
  use ModReadArtep
  use ModPlotFile
  use ModMpi
  use CRASH_ModAtomicMass,ONLY: cAtomicMass_I
  use ModConst
  implicit none
  integer:: iError, iMaterial
  real:: AtomicMass
  !----------------
  call read_artep
  AtomicMass = cAtomicMass_I(54) * cAtomicMass
  Value_VII = 0.10 * Value_VII  !cm^2/g = 1e-4 m^2/(1e-3 kg) = 0.1 m2/kg
  call save_plot_file( &
         'Xe_opac.dat',                                &
         TypeFileIn     = 'real8',                     &
         StringHeaderIn = 'ARTEP Opacity for Xenon', &
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
                                      ' Ross26 Ross27 Ross28 Ross29 Ross30 EvMin EvMax', &
         CoordMinIn_D   = (/log10(AtomicMass*1.0e+24), log10(0.030)/),       &                             
         CoordMaxIn_D   = (/log10(AtomicMass*1.0e+29),          3.0/),       &
         ParamIn_I      = (/0.1,2.0e+4/),                                    &
         VarIn_VII      = Value_VII)
  stop
end program save_eos_table
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
