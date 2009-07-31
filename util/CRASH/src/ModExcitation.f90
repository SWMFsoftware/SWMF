!^CFG COPYRIGHT UM

module CRASH_ModExcitation
  use CRASH_ModAtomicDataMix
  use CRASH_ModExcitationData,ONLY : n_ground
  implicit none
  PRIVATE !Except


  !Excitation energy of ion of ionization state i, averaged by possible
  !excitation levels [eV]
  real,dimension(0:nZMax,nMixMax),public :: ExtraEnergyAv_II = 0.0

  !The average value of -V\frac{\partial \Delta E}{\partial V}
  !for given i and iMix [eV]
  real,dimension(0:nZMax,nMixMax),public :: VirialCoeffAv_II = 0.0

  !The average value of V^2\frac{\partial^2 \Delta E}{\partial V^2}
  !for given i and iMix [eV]
  real,dimension(0:nZMax,nMixMax),public :: VirialCoeff2Av_II = 0.0

  !The variance of the extra energy, E_x, for given i and iMix [eV^2]
  real,dimension(0:nZMax,nMixMax),public :: Cov2ExtraEnergy_II = 0.0

  !The covariance between the extra energy and -V\frac{\partial \Delta E}{\partial V}
  !for given i and iMix [eV^2]
  real,dimension(0:nZMax,nMixMax),public :: CovExtraEnergyVirial_II = 0.0

  !The variance of -V\frac{\partial \Delta E}{\partial V}
  !for given i and iMix [eV^2]
  real,dimension(0:nZMax,nMixMax),public :: Cov2VirialCoeff_II = 0.0


  !\logarithm of the statistical sub-sum over the excitated states, for
  ! a given sort of ions.
  real,dimension(0:nZMax,nMixMax),public :: LogGi_II = 0.0


  !\
  ! Partition function for the excited states
  !/
  real,dimension(0:nExcitation-1,nExcitation,0:nZMax,nMixMax) :: Partition_IIII = 0.0


  public :: get_excitation_levels
  public :: nMixMax
  real :: PowerOfRIono

  real,public :: IonizationPotentialLowering_I(0:nZMax) = 0.0
contains

  !Fill in arrays ExtraEnergyAv_II and LogGi_II
  subroutine get_excitation_levels(iMix, iZMin, iZMax, nZ, TeInv, rIonoSphereInv)
    integer,intent(in) :: iMix
    integer,intent(in) :: iZMin, iZMax
    integer,intent(in) :: nZ
    real,intent(in) :: TeInv
    real,intent(in) :: rIonoSphereInv

    integer :: iZ, iN, iL  !Loop variables

    real    :: DeltaEnergy
    real    :: Gi, GiInv
    integer :: nGround
  
    real,parameter :: IndexPerThree    = iPressureEffectIndex/3.0
    real,parameter :: IndexSecondDeriv = IndexPerThree * (1.0 + IndexPerThree)

    !Energy correction due to the pressure ionization effect, averaged by
    !possible excitation levels [eV]
    real,dimension(0:nZMax,nMixMax) :: DeltaEnergyAv_II = 0.0

    !Energies of excited levels relative to the ground level, E_x = E_{exc} + \Delta E,
    !where E_{exc} stands for the excitation energy, and \Delta E is the correction
    !accounting for the pressure ionization
    real,dimension(0:nExcitation-1,nExcitation,0:nZMax,nMixMax) :: ExtraEnergy_IIII = 0.0

    ! -V \partial E_x/\partial V
    real,dimension(0:nExcitation-1,nExcitation,0:nZMax,nMixMax) :: VirialCoeff_IIII = 0.0
    !------------

    if (.not.UseExcitation) return


    Partition_IIII     (:,:,:,iMix) = 0.0
    DeltaEnergyAv_II       (:,iMix) = 0.0
    ExtraEnergyAv_II       (:,iMix) = 0.0
    VirialCoeffAv_II       (:,iMix) = 0.0
    VirialCoeff2Av_II      (:,iMix) = 0.0
    Cov2ExtraEnergy_II     (:,iMix) = 0.0
    CovExtraEnergyVirial_II(:,iMix) = 0.0
    Cov2VirialCoeff_II     (:,iMix) = 0.0

    PowerOfRIono = rIonoSphereInv**iPressureEffectIndex

    do iZ = iZMin, iZMax

       nGround = n_ground(iZ, nZ)
       Gi = 0.0

       do iN = nGround, nExcitation
          do iL = 0, iN-1

             if(UsePressureIonization.and.(ExcitationEnergy_IIII(iL,iN,iZ,iMix) + &
                  VirialCoeff4Energy_IIII(iL,iN,iZ,iMix) * PowerOfRIono >= &
                  IonizPotential_II(iZ+1,iMix) - IonizationPotentialLowering_I(iZ)))CYCLE

             DeltaEnergy = VirialCoeff4Energy_IIII(iL,iN,iZ,iMix) * PowerOfRIono
             ExtraEnergy_IIII(iL,iN,iZ,iMix) = ExcitationEnergy_IIII(iL,iN,iZ,iMix) + &
                  DeltaEnergy

             Partition_IIII(iL,iN,iZ,iMix) = Degeneracy_IIII(iL,iN,iZ,iMix)*&
                  exp(-ExtraEnergy_IIII(iL,iN,iZ,iMix) * TeInv)


             Gi = Gi + Partition_IIII(iL,iN,iZ,iMix)

             DeltaEnergyAv_II(iZ,iMix) = DeltaEnergyAv_II(iZ,iMix) + &
                  Partition_IIII(iL,iN,iZ,iMix) * DeltaEnergy
             ExtraEnergyAv_II(iZ,iMix) = ExtraEnergyAv_II(iZ,iMix) + &
                  Partition_IIII(iL,iN,iZ,iMix) * ExtraEnergy_IIII(iL,iN,iZ,iMix)

             VirialCoeff_IIII(iL,iN,iZ,iMix) = IndexPerThree * DeltaEnergy

          end do
       end do


       if(Gi==0.0)then
          LogGi_II(iZ,iMix) = -30.0
          GiInv = 0.0
          CYCLE
       else
          LogGi_II(iZ,iMix) = log(Gi)
          GiInv = 1.0/Gi
       end if
       
       !Normalize to obtain the partial partition function and average energies
       Partition_IIII(:,nGround:nExcitation,iZ,iMix) = Partition_IIII(:,nGround:nExcitation,iZ,iMix) * GiInv
       ExtraEnergyAv_II(iZ,iMix) = ExtraEnergyAv_II(iZ,iMix) * GiInv
       DeltaEnergyAv_II(iZ,iMix) = DeltaEnergyAv_II(iZ,iMix) * GiInv

       VirialCoeffAv_II (iZ,iMix) = IndexPerThree    * DeltaEnergyAv_II(iZ,iMix)
       VirialCoeff2Av_II(iZ,iMix) = IndexSecondDeriv * DeltaEnergyAv_II(iZ,iMix)



       do iN = nGround, nExcitation
          do iL = 0, iN-1
             if (Partition_IIII(iL,iN,iZ,iMix) == 0.0)CYCLE

             Cov2ExtraEnergy_II(iZ,iMix) = Cov2ExtraEnergy_II(iZ,iMix) + &
                  Partition_IIII(iL,iN,iZ,iMix) * &
                  (ExtraEnergy_IIII(iL,iN,iZ,iMix) - ExtraEnergyAv_II(iZ,iMix))**2

             CovExtraEnergyVirial_II(iZ,iMix) = CovExtraEnergyVirial_II(iZ,iMix) + &
                  Partition_IIII(iL,iN,iZ,iMix) * &
                  (ExtraEnergy_IIII(iL,iN,iZ,iMix) - ExtraEnergyAv_II(iZ,iMix)) * &
                  VirialCoeff_IIII(iL,iN,iZ,iMix)

             Cov2VirialCoeff_II(iZ,iMix) = Cov2VirialCoeff_II(iZ,iMix) + &
                  Partition_IIII(iL,iN,iZ,iMix) * &
                  (VirialCoeff_IIII(iL,iN,iZ,iMix) - VirialCoeffAv_II(iZ,iMix))**2

          end do
       end do

    end do
  end subroutine get_excitation_levels

end module CRASH_ModExcitation
