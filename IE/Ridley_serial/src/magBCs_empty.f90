!^CFG COPYRIGHT UM
subroutine ionosphere_magBCs(PHI_BC, ETh_BC, EPs_BC,                         &
                             UR_BC, UTh_BC, UPs_BC, Radius_BC,               &
                             PHI, X, Y, Z,                                   &
                             Theta_BC, Psi_BC, Radius, nTheta, nPsi,         &
                             dTheta, dPsi)

  use ModIonosphere
  implicit none

  integer :: nTheta, nPsi
  real :: Radius, Radius_BC
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::                              &
                  PHI_BC, ETh_BC, EPs_BC,                                    &
                  UR_BC, UTh_BC, UPs_BC,                                     &
                  PHI, X, Y, Z, Theta_BC, Psi_BC

  real, dimension(1:IONO_nTheta,1:IONO_nPsi,2) ::                            &
                  Theta_temp, Psi_temp, phi_temp

  real, dimension(1:IONO_nTheta) :: dTheta
  real, dimension(1:IONO_nPsi)   :: dPsi

  ! This routine does nothing.

end subroutine ionosphere_magBCs
