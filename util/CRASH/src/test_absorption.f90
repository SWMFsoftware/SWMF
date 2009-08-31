program abs
  use CRASH_ModStatSum
  use CRASH_ModPartition
  use CRASH_ModPolyimide
  use CRASH_ModFermiGas
  use CRASH_ModMultiGroup
  use CRASH_ModExcitation
  use CRASH_ModAtomicDataMix
  use CRASH_ModExcitationData,ONLY : n_ground, cExcitationN_III
  use ModConst
  implicit NONE

  integer, parameter :: unit = 24

  real:: vTe = 10.0 !eV
  real::NaTrial = 1.0e22
  integer:: iPlot,iError
  integer :: iL, iN, iZ, iMix
  integer :: nGround

  real,dimension(1:nZMax) :: IonizPotential_I
  !---------------

  open(unit,file='../doc/excited_levels.tex',status='replace')
  write(unit,'(a)')'\newcolumntype{x}[1]{>{\centering\hspace{0pt}}p{#1}}'
  write(unit,'(a)')'\begin{tabular}'//&
       '{|x{1cm}|x{1cm}||x{2cm}||x{2cm}|x{2cm}|x{2cm}|x{2cm}|x{2cm}|}'
  write(unit,'(a)')'\hline'
  write(unit,'(a)')'\multirow{2}{*}{i} & \multirow{2}{*}{n} & '//&
       '\multirow{2}{*}{formula} & '//&
       '\multicolumn{5}{c|}{database, for different values of $l$} \tabularnewline'
  write(unit,'(a)')'\cline{4-8}'
  write(unit,'(a)')' & & & s & p & d & f & g \tabularnewline'
  write(unit,'(a)')'\hline'
  write(unit,'(a)')'\hline'

  call get_ioniz_potential(7, IonizPotential_I(1:7))

  do iZ = 0, 6
     write(unit,'(a,i3,a)') '\multirow{4}{*}{', iZ, '} '

     do iN = 2, 5

        nGround = n_ground(iZ, 7)

        write(unit,'(a,i3,a,f6.1,5(a,f6.1),a)') ' & ', iN, ' & ', &
             max(IonizPotential_I(iZ+1) - &
             cRyToEV * (real(iZ+1)/iN)**2, 0.0), &
             (' & ', cExcitationN_III(iL, iN, iZ), iL = 0, 4), &
             ' \tabularnewline'

     end do

     write(unit,'(a)')'\hline'

  end do

  write(unit,'(a)')'\end{tabular}'
  close(unit)



  UseExcitation = .true.
  call set_mixture(nPolyimide, nZPolyimide_I, CPolyimide_I)
  UsePreviousTe = .false.

  DoNotAddLineCore = .false.
  UseBremsstrahlung = .false.
  UsePhotoionization = .false.
 
  UseCoulombCorrection = .true.
  call set_ionization_equilibrium(vTe,NaTrial*1000000.0,iError)
  open(unit,file='report.txt')
  do iMix = 1,nMix
     write(unit,*)'iMix,nZ_I(iMix),iZMin_I(iMix), iZMax_I(iMix)',iMix,nZ_I(iMix),iZMin_I(iMix), iZMax_I(iMix)
     do iZ = 0, nZ_I(iMix)
        write(unit,*)' iZ = ', iZ
        write(unit,*)Partition_III(:,iZ,iMix)
     end do
  end do
  close(unit)
  call set_multigroup(100, 0.1/cHPlanckEV, 1000.0/cHPlanckEV)
  call meshhv
  call abscon
  
  open(unit,file='../doc/polyimide_absorption.dat')
  write(unit,'(a,i6,a)') &
       'Photon energy [eV]  Absorbtion Coeff cm-1, in ',nPhoton,' points' 
  do iPlot = 1, nPhoton
     if(PhotonEnergy_I(iPlot)< 0.01*Te.or.PhotonEnergy_I(iPlot)>100.0*Te)&
          CYCLE
     write(unit,*)log10(PhotonEnergy_I(iPlot)),&
          log10(AbsorptionCoefficient_I(iPlot))
  end do
  close(unit)
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

