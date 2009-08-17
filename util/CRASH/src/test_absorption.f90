program abs
  use CRASH_ModStatSum
  use CRASH_ModStatSumMix
  use CRASH_ModPolyimide
  use CRASH_ModFermiGas
  use CRASH_ModMultiGroup
  use CRASH_ModExcitation
  use CRASH_ModAtomicDataMix
  implicit NONE
  real:: vTe = 10.0 !eV
  real::NaTrial = 1.0e22
  integer:: iPlot,iError
  integer::iMix,iZ
  !---------------
  UseExcitation = .true.
  call set_mixture(nPolyimide, nZPolyimide_I, CPolyimide_I)
  UsePreviousTe = .false.

  DoNotAddLineCore = .false.
  UseBremsstrahlung = .false.
  UsePhotoionization = .false.
 
  UseCoulombCorrection = .true.
  call set_ionization_equilibrium(vTe,NaTrial*1000000.0,iError)
  open(24,file='report.txt')
  do iMix = 1,nMix
     write(24,*)'iMix,nZ_I(iMix),iZMin_I(iMix), iZMax_I(iMix)',iMix,nZ_I(iMix),iZMin_I(iMix), iZMax_I(iMix)
     do iZ = 0, nZ_I(iMix)
        write(24,*)' iZ = ', iZ
        write(24,*)Partition_III(:,iZ,iMix)
     end do
  end do
  close(24)
  call set_default_multigroup
  call meshhv
  call abscon
  
  open(24,file='../doc/Polymide_Abs.dat')
  write(24,'(a,i6,a)') &
       'Photon energy [eV]  Absorbtion Coeff cm-1, in ',nPhoton,' points' 
  do iPlot = 1, nPhoton
     if(PhotonEnergy_I(iPlot)< 0.01*Te.or.PhotonEnergy_I(iPlot)>100.0*Te)&
          CYCLE
     write(24,*)log10(PhotonEnergy_I(iPlot)),log10(abscfs(iPlot))
  end do
  close(24)
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

