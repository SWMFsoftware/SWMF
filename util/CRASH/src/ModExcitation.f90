!^CFG COPYRIGHT UM

module CRASH_ModExcitation
  use CRASH_ModAtomicMass,ONLY : nZMax
  implicit none
  PRIVATE !Except

  integer, parameter, public :: nMixMax = 6

  !\
  ! The logical to handle whether excitation levels should
  ! be accounted for
  !/
  logical,public :: UseExcitation = .false.

  !Excitation energy of ion of ionization state i, averaged by possible
  !excitation levels [eV]
  real,dimension(0:nZMax,nMixMax),public :: ExcitationEnergyAv_II = 0.0

  !\logarithm of the statistical sub-sum over the excitated states, for
  ! a given sort of ions.
  real,dimension(0:nZMax,nMixMax),public :: LogGi_II = 0.0
  
  !\
  !  Maximum principal quantum number for the excitation states.
  !/
  integer,parameter:: nExcitation = 10
  
  !\
  ! Partition function for the excited states
  !/
  real,dimension(nExcitation,0:nZMax,nMixMax)::Partition_III = 0.0

  !The principal quantum number of the outermost electron in bounded with an ion
  !with I electrons
  integer,parameter :: nGround_I(100) = (/ 1, 1, &                                                       !  2
       2, 2, 2, 2, 2, 2, 2, 2, &                                                                         !  8
       3, 3, 3, 3, 3, 3, 3, 3, &                                                                         !  8
       4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &                                           ! 18
       5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &                                           ! 18
       6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, & ! 32
       7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7 /)

  public :: get_excitation_levels
contains

  !Fill in arrays ExcitationEnergyAv_II and LogGi_II
  subroutine get_excitation_levels(iMix, iZMin, iZMax, nZ, IonizPotential_I, TeInv)
    integer :: iMix
    integer :: iZMin, iZMax
    integer :: nZ
    integer :: iZ, iN  !Loop variables

    real, dimension(nZMax),intent(in) :: IonizPotential_I
    real :: TeInv

    real    :: Gi, GiInv
    integer :: nGround
    real    :: ExcitationEnergy
    !------------

    if (UseExcitation) then
       Partition_III(:,:,iMix) = 0.0
       ExcitationEnergyAv_II(:,iMix) = 0.0

       do iZ = iZMin, iZMax
      
          nGround = nGround_I(nZ - iZ)

          Partition_III(nGround,iZ,iMix) = 2.0 * nGround*nGround
          Gi = Partition_III(nGround,iZ,iMix) 
          do iN = nGround+1, nExcitation

             ExcitationEnergy = IonizPotential_I(iZ+1) * (1.0 - (real(nGround) / iN)**2)
             Partition_III(iN,iZ,iMix) = 2.0 * iN*iN * exp(-ExcitationEnergy * TeInv)

             Gi = Gi + Partition_III(iN,iZ,iMix)
             ExcitationEnergyAv_II(iZ,iMix) = ExcitationEnergyAv_II(iZ,iMix) + &
                  Partition_III(iN,iZ,iMix) * ExcitationEnergy

          end do

          LogGi_II(iZ,iMix) = log(Gi); GiInv = 1.0/Gi

          !Normalize to obtain the average excitation energy
          ExcitationEnergyAv_II(iZ,iMix) = ExcitationEnergyAv_II(iZ,iMix) * GiInv
          Partition_III(nGround:nExcitation,iZ,iMix) = Partition_III(nGround:nExcitation,iZ,iMix) * GiInv
       end do

    end if

  end subroutine get_excitation_levels
    
end module CRASH_ModExcitation
