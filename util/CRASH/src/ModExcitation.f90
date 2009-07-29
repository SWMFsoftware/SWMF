!^CFG COPYRIGHT UM

module CRASH_ModExcitation
  use CRASH_ModAtomicDataMix
  use CRASH_ModExcitationData,ONLY : n_ground
  implicit none
  PRIVATE !Except


  !Excitation energy of ion of ionization state i, averaged by possible
  !excitation levels [eV]
  real,dimension(0:nZMax,nMixMax),public :: ExtraEnergyAv_II = 0.0

  !The average value of -V\frac{\partial \Delta E}{\partial V} for given i and iMix
  real,dimension(0:nZMax,nMixMax),public :: VirialCoeffAv_II = 0.0


  !\logarithm of the statistical sub-sum over the excitated states, for
  ! a given sort of ions.
  real,dimension(0:nZMax,nMixMax),public :: LogGi_II = 0.0
  real,dimension(0:nZMax,nMixMax) :: GiInv_II = 0.0


  !\
  ! Partition function for the excited states
  !/
  real,dimension(0:nExcitation-1,nExcitation,0:nZMax,nMixMax)::Partition_IIII = 0.0

  real,dimension(0:nExcitation-1,nExcitation,0:nZMax,nMixMax) :: ExtraEnergy_IIII = 0.0
 
  public :: get_excitation_levels
  public :: nMixMax
  real ::PowerOfRIono
contains



  logical function is_dead_orbit(iL,iN,iZ,iMix,nZ,rIonoSphereInv)
    use ModConst
    integer,intent(in) :: iL, iN, iZ, iMix
    integer,intent(in) :: nZ
    real,intent(in)    :: rIonoSphereInv
    !--------------------------------
    integer :: nGround
    
  

    nGround = n_ground(iZ, nZ)
    if (iN == nGround) then
       is_dead_orbit = .false.
       return
    end if
    
    is_dead_orbit = (ExcitationEnergy_IIII(iL,iN,iZ,iMix) + &
         VirialCoeff4Energy_IIII(iL,iN,iZ,iMix) * PowerOfRIono >= &
         IonizPotential_II(iZ+1,iMix) - 1.8 * (2.0*iZ + 1) * &
         cRyToEv * rIonoSphereInv)

  end function is_dead_orbit

  !Fill in arrays ExtraEnergyAv_II and LogGi_II
  subroutine get_excitation_levels(iMix, iZMin, iZMax, nZ, TeInv, rIonoSphereInv)
    integer :: iMix
    integer :: iZMin, iZMax
    integer :: nZ
    integer :: iZ, iN, iL  !Loop variables

    real :: TeInv
    real :: rIonoSphereInv
    real :: eMadelungPerTe

    real    :: Gi, GiInv
    integer :: nGround
  
    real,parameter:: IndexPerThree = iPressureEffectIndex/3.0
    !------------

    if (.not.UseExcitation) return

    Partition_IIII(:,:,:,iMix) = 0.0
    ExtraEnergyAv_II(:,iMix) = 0.0
    VirialCoeffAv_II(:,iMix) = 0.0

    PowerOfRIono =    rIonoSphereInv**iPressureEffectIndex

    do iZ = iZMin, iZMax
       
       nGround = n_ground(iZ, nZ)
       Gi = 0.0

       do iN = nGround, nExcitation
          do iL = 0, iN-1
             if (is_dead_orbit(iL, iN, iZ, iMix, nZ, rIonoSphereInv))CYCLE
                             
                ExtraEnergy_IIII(iL,iN,iZ,iMix) = ExcitationEnergy_IIII(iL,iN,iZ,iMix) + &
                     VirialCoeff4Energy_IIII(iL,iN,iZ,iMix) * PowerOfRIono
                
                Partition_IIII(iL,iN,iZ,iMix) = Degeneracy_IIII(iL,iN,iZ,iMix)*&
                     exp(-ExtraEnergy_IIII(iL,iN,iZ,iMix) * TeInv)
                
                
                Gi = Gi + Partition_IIII(iL,iN,iZ,iMix)
                
                ExtraEnergyAv_II(iZ,iMix) = ExtraEnergyAv_II(iZ,iMix) + &
                     Partition_IIII(iL,iN,iZ,iMix) * ExtraEnergy_IIII(iL,iN,iZ,iMix)

                VirialCoeffAv_II(iZ,iMix) = VirialCoeffAv_II(iZ,iMix) + &
                      Partition_IIII(iL,iN,iZ,iMix) *  PowerOfRIono *  IndexPerThree * &
                      VirialCoeff4Energy_IIII(iL,iN,iZ,iMix)
          end do
       end do
       
       if(Gi==0.0)then
          LogGi_II(iZ,iMix) = -30.0
          GiInv_II(iZ,iMix) = 0.0
          GiInv = 0.0
       else
          LogGi_II(iZ,iMix) = log(Gi)
          GiInv = 1.0/Gi
          GiInv_II(iZ,iMix) = GiInv
       end if
       
       !Normalize to obtain the average excitation energy
       ExtraEnergyAv_II(iZ,iMix) = ExtraEnergyAv_II(iZ,iMix) * GiInv
       Partition_IIII(:,nGround:nExcitation,iZ,iMix) = Partition_IIII(:,nGround:nExcitation,iZ,iMix) * GiInv
       VirialCoeffAv_II(iZ,iMix) = VirialCoeffAv_II(iZ,iMix) * GiInv_II(iZ,iMix)
    end do
  end subroutine get_excitation_levels

end module CRASH_ModExcitation
