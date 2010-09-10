
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
  use CRASH_ModEos, ONLY:UseEosTable_I, check_eos_table
  use ModMpi
  implicit none
  integer:: iError
  !----------------
  UseCoulombCorrection = .true.
  UseExcitation = .true.
  UseEosTable_I = .true.
  
  call MPI_Init(iError)
  call check_eos_table(MPI_COMM_WORLD)!.true.)
  call MPI_Finalize(iError)
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
