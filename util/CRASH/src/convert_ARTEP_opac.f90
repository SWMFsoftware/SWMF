module ModReadArtep
  use ModIoUnit, ONLY: io_unit_new
  implicit none

  integer, parameter:: nRho = 201, nTe = 201, nGroup = 30

  SAVE
  real, allocatable:: Value_VII(:,:,:) 
  integer::line
contains
  !==================================================================================
  subroutine read_artep
    integer:: iFile, iRho, iTe
    real:: Aux_V(7), Aux
    !--------------------------------------------------------------------------------

    allocate(Value_VII(2*nGroup,nRho,nTe))

    iFile = io_unit_new()
    open(iFile, file='xz_opac.artep', form='formatted',status='old')
    line=0
    do iTe = 1, nTe
       do iRho = 1, nRho
          line = line+1;write(*,*)line
          read(iFile,*) Aux_V, Value_VII(1:nGroup,iRho,iTe), &
               Aux,            Value_VII(nGroup+1:2*nGroup,iRho,iTe)
       end do
    end do
    close(iFile)
  end subroutine read_artep

end Module ModReadArtep
!====================================================================================
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
  !----------------------------------------------------------------------------------
  call read_artep
  AtomicMass = cAtomicMass_I(54) * cAtomicMass
  Value_VII = 0.10 * Value_VII  !cm^2/g = 1e-4 m^2/(1e-3 kg) = 0.1 m2/kg
  call save_plot_file( &
         'Xe_opac.dat',                                                      &
         TypeFileIn     = 'real8',                                           &
         StringHeaderIn = 'ARTEP Opacity for Xenon',                         &
         NameVarIn      = 'logRho logTe Planck(30) Ross(30) EvMin EvMax',    &
         CoordMinIn_D   = (/log10(AtomicMass*1.0e+24), log10(0.03)/),        &
         CoordMaxIn_D   = (/log10(AtomicMass*1.0e+29),       3.0  /),        &
         ParamIn_I      = (/0.1, 2e4/),                                      &
         VarIn_VII      = Value_VII)

end program save_eos_table
!====================================================================================
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
