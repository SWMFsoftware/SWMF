program abs
  use CRASH_ModStatSum
  use CRASH_ModStatSumMix
  use CRASH_ModPolyimide
  use CRASH_ModFermiGas
  use CRASH_ModMultiGroup
  implicit NONE
  real:: vTe = 10.0 !eV
  real::NaTrial = 1.0e22
  integer:: iPlot,iError
  !---------------
  call set_mixture(nPolyimide, nZPolyimide_I, CPolyimide_I)
  UsePreviousTe = .false.
  call set_ionization_equilibrium(vTe,NaTrial*1000000.0,iError)
  call set_default_multigroup
  call meshhv
  call abscon
  do iPlot = 1, nPhoton
     write(*,*)PhotonEnergy_I(iPlot),abscfs(iPlot)
  end do

end program abs
!============================================================================
! The following subroutines are here so that we can use SWMF library routines
! Also some features available in SWMF mode only require empty subroutines
! for compilation of the stand alone code.
!============================================================================
subroutine CON_stop(StringError)
  implicit none
  character (len=*), intent(in) :: StringError
end subroutine CON_stop
!============================================================================
subroutine CON_set_do_test(String,DoTest,DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe
end subroutine CON_set_do_test

