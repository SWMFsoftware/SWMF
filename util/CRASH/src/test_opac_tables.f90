
program save_opac_table
  use CRASH_ModStatSum
  use CRASH_ModPartition, ONLY: UseCoulombCorrection
  use CRASH_ModPolyimide
  use CRASH_ModFermiGas
  use CRASH_ModMultiGroup
  use CRASH_ModExcitation
  use CRASH_ModAtomicDataMix
  use CRASH_ModExcitationData,ONLY : n_ground, cExcitationN_III
  use CRASH_ModIonMix
  use CRASH_ModEos, ONLY:UseEosTable_I, UseOpacityTable_I,&
                         check_opac_table,IndexDefaultOpac_I
  use ModMpi
  implicit none
  integer:: iError, iMaterial
  !----------------
  UseCoulombCorrection = .false.
  UseExcitation = .true.
  UseEosTable_I = .false.
  UseOpacityTable_I = .true.
  !Set reduced table izes
  IndexDefaultOpac_I = (/9,9/)
  call set_multigroup(10,0.1/cHPlanckEV,20000.0/cHPlanckEV)
  
  call MPI_Init(iError)
  call check_opac_table(MPI_COMM_WORLD,.true.,'ascii')
  call MPI_Finalize(iError)
end program save_opac_table
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
