module CRASH_ModExcitationData
  use CRASH_ModIonization,ONLY : get_ioniz_potential
  implicit none
  PRIVATE !except

  !public methods
  public:: n_ground  !Principal quantum number as a function of the charge state and the element number
  public:: get_excitation_energy
  public:: get_virial_coeff4_energy
  public:: get_degeneracy
contains
  !=====================================================================================
  integer function n_ground(iZ,nZ)
    integer,intent(in)::iZ,nZ
    !The principal quantum number of the outermost electron in bounded with an ion
    !with I electrons
    integer,parameter :: nGround0_I(100) = (/ &
         1, 1, &                                                                    !  2
         2, 2, 2, 2, 2, 2, 2, 2, &                                                  !  8
         3, 3, 3, 3, 3, 3, 3, 3, &                                                  !  8
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &                    ! 18
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &                    ! 18
         6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &                          !
         6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &                          ! 32
         7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7 /)                                !
    integer,parameter :: nGround_I(100) = (/ &
         1, 1, &                                                                    !  2
         2, 2, 2, 2, 2, 2, 2, 2, &                                                  !  8
         3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &                    ! 18
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,       &                    ! 
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,       &                    ! 32
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,       &                    !
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,       &                    !
         5, 5, 5, 5, 5, 5, 5, 5    /)                                               ! 40
    !-------------------------------------------------------------------------------!
    if(iZ==0)then 
       n_ground = nGround0_I(nZ)
    else
       n_ground = nGround_I(nZ - iZ)
    end if
  end function n_ground
  !======================================================================================
  !Fill in the array ExcitationEnergy_III with tabulated or calculted excitation energies
  !for quantum states defined by n, l of atoms or ions of a particular element defined
  !by their charge state
  subroutine get_excitation_energy(nExcitation, nZ, ExcitationEnergy_III)
    integer :: nExcitation
    integer :: nZ
    real,dimension(0:nExcitation-1,nExcitation,0:nZ-1),intent(out) :: ExcitationEnergy_III
    real,dimension(1:nZ) :: IonizPotential_I

    integer :: iZ, iN, nGround
    !--------------
    ExcitationEnergy_III = 0.0

    call get_ioniz_potential(nZ, IonizPotential_I(1:nZ))

    do iZ = 0, nZ-1
       nGround = n_ground(iZ, nZ)

       do iN = nGround+1, nExcitation
          ExcitationEnergy_III(0:iN-1,iN,iZ) = IonizPotential_I(iZ+1) * &
               (1.0 - (real(nGround) / iN)**2.0)
       end do
    end do
  end subroutine get_excitation_energy
  !======================================================================================
  subroutine get_virial_coeff4_energy(nExcitation, nZ, VirialCoeff4Energy_III)
    integer :: nExcitation
    integer :: nZ
    real,dimension(0:nExcitation-1,nExcitation,0:nZ-1),intent(out) :: VirialCoeff4Energy_III
    real,dimension(1:nZ) :: IonizPotential_I

    integer :: iZ, iN, iL, nGround
    !--------------
    VirialCoeff4Energy_III = 0.0

    call get_ioniz_potential(nZ, IonizPotential_I(1:nZ))

    do iZ = 0, nZ-1
       nGround = n_ground(iZ, nZ)

       do iN = nGround, nExcitation
          do iL = 0, iN-1
             VirialCoeff4Energy_III(iL,iN,iZ) = IonizPotential_I(iZ+1) * &
                  real(nGround)**2 * real(iN - iL)**2 / real(iZ+1)**2
          end do
       end do
    end do
  end subroutine get_virial_coeff4_energy
  !======================================================================================
  subroutine get_degeneracy(nExcitation, nZ, Degeneracy_III)
    integer :: nExcitation
    integer :: nZ
    integer,dimension(0:nExcitation-1,nExcitation,0:nZ-1),intent(out) :: Degeneracy_III
    integer :: iMix

    integer :: iZ, iN, iL, nGround
    !--------------
    Degeneracy_III = 0

    do iZ = 0, nZ-1
       nGround = n_ground(iZ, nZ)

       do iN = nGround, nExcitation
          do iL = 0, iN-1
             Degeneracy_III(iZ,iN,iL) = 2 * (2*iL + 1)
          end do
       end do
    end do
  end subroutine get_degeneracy
  !======================================================================================
end module CRASH_ModExcitationData
