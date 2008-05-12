!======================================
!                                     |
!    Module for Ionosphere Model      |
!                                     |
!======================================

module ModIonosphere

  !\
  ! Ionosphere array parameters
  !/
  integer, parameter :: IONO_nTheta = 65
  integer, parameter :: IONO_nPsi = 4*(IONO_nTheta-1)+1
  integer, parameter :: IONO_nRadial = 21
  integer, parameter :: IONO_NGMAX = 6
  integer :: IONO_NORTH_nMagBndPts, IONO_SOUTH_nMagBndPts

  !\
  ! Ionospheric Model Parameters
  !/

  integer :: ionosphere_type
  logical :: UseIonosphere
  logical :: TIMEGCM_Coupled
  logical :: RCM_Coupled
  logical :: new_conductance
  real    :: f107_flux, PolarCapPedConductance, StarLightPedConductance
  real    :: Hall_to_Ped_Ratio

  !\
  ! Some time dependent stuff
  !/

  real :: IONO_Old_Seconds

  !\
  ! IONOSPHERE RESTART ERROR stuff
  !/
  logical :: skip_ionosphere_restart

  !\
  ! Ionosphere solution parameters
  !/
  integer, parameter :: IONO_NPRE = 10, &
                        IONO_NPOST = 10, &
                        IONO_NEXACT = 25

  real, parameter :: IONO_OMEGA = 1.00,                       &
                     IONO_TOLER = 5.0e-05,                    &
                     IONO_MU = 1.256637e-06,                  &
                     IONO_Theta_0 = 0.0001,                   &
                     IONO_Min_EFlux = 0.1e-16,                &     ! W/m2
                     IONO_Min_Ave_E = 0.5,                    &     ! keV
                     Polar_Rain = 0.1e-2                            ! W/m2


  integer, parameter :: IONO_Model_No_Hall = 1,               &
                        IONO_Model_With_Hall = 2,             &
                        IONO_Model_With_Simple_Aurora = 3,    &
                        IONO_Model_With_Complex_Aurora = 4

  real :: IONO_PI, IONO_Bdp, IONO_Radius_Mag_Boundary, &
          IONO_NORTH_Theta_Max, IONO_SOUTH_Theta_Min, &
          IONO_Radius, IONO_Height, IONO_Ref_Density, &
          IONO_Ref_SoundSpeed

  real, dimension(1:360,1:90) :: &
       hal_a0_up,ped_a0_up,      &
       hal_a0_do,ped_a0_do,      &
       hal_a1_up,ped_a1_up,      &
       hal_a1_do,ped_a1_do,      &
       hal_a2_up,ped_a2_up,      &
       hal_a2_do,ped_a2_do

  real, dimension(1:360) :: cond_mlts
  real, dimension(1:90) :: cond_lats

  integer :: i_cond_nmlts, i_cond_nlats

  !\
  ! Ionosphere solution array definitions
  !/
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::  &
    IONO_NORTH_PHI,IONO_SOUTH_PHI,               & !Ionospheric potential
    IONO_NORTH_X,IONO_NORTH_Y,IONO_NORTH_Z,      & !Ionospheric coordinates
    IONO_NORTH_Theta,IONO_NORTH_Psi,             & !
    IONO_SOUTH_X,IONO_SOUTH_Y,IONO_SOUTH_Z,      & !
    IONO_SOUTH_Theta,IONO_SOUTH_Psi,             & !
    IONO_NORTH_JR,                               & !Ionospheric current
    IONO_NORTH_JTh,IONO_NORTH_JPs,               &
    IONO_NORTH_Jx,IONO_NORTH_Jy,IONO_NORTH_Jz,   &
    IONO_SOUTH_JR,                               &
    IONO_SOUTH_JTh,IONO_SOUTH_JPs,               &
    IONO_SOUTH_Jx,IONO_SOUTH_Jy,IONO_SOUTH_Jz,   & !
    IONO_NORTH_TGCM_JR,                          & !TIMEGCM current
    IONO_SOUTH_TGCM_JR,                          & !
    IONO_NORTH_RCM_JR,                           & !RCM current
    IONO_SOUTH_RCM_JR,                           & !
    IONO_NORTH_AMIE_JR,                          & !AMIE current
    IONO_SOUTH_AMIE_JR,                          & !
    IONO_NORTH_Ex,IONO_NORTH_Ey,IONO_NORTH_Ez,   & !Ionospheric electric field
    IONO_NORTH_ETh,IONO_NORTH_EPs,               & !
    IONO_SOUTH_Ex,IONO_SOUTH_Ey,IONO_SOUTH_Ez,   & !
    IONO_SOUTH_ETh,IONO_SOUTH_EPs,               & !
    IONO_NORTH_Ux,IONO_NORTH_Uy,IONO_NORTH_Uz,   & !Ionospheric flow velocities
    IONO_NORTH_UTh,IONO_NORTH_UPs,               & !
    IONO_SOUTH_Ux,IONO_SOUTH_Uy,IONO_SOUTH_Uz,   & !
    IONO_SOUTH_UTh,IONO_SOUTH_UPs,               & !
    IONO_NORTH_EFlux,IONO_NORTH_Ave_E,           & !Ionospheric particle inputs
    IONO_SOUTH_EFlux,IONO_SOUTH_Ave_E,           & !
    IONO_NORTH_Sigma0,                           & !Ionospheric conductances
    IONO_NORTH_SigmaH,IONO_NORTH_SigmaP,         &
    IONO_NORTH_SigmaThTh,IONO_NORTH_SigmaThPs,   & !
    IONO_NORTH_SigmaPsPs
    
    
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::  &
    IONO_SOUTH_Sigma0,                           &
    IONO_SOUTH_SigmaH,IONO_SOUTH_SigmaP,         & !
    IONO_SOUTH_SigmaThTh,IONO_SOUTH_SigmaThPs,   & !
    IONO_SOUTH_SigmaPsPs,                        & !
    IONO_NORTH_dSigmaThTh_dTheta,                &
    IONO_NORTH_dSigmaThPs_dTheta,                & !Ionospheric cond. gradients
    IONO_NORTH_dSigmaPsPs_dTheta,                &
    IONO_NORTH_dSigmaThTh_dPsi,                  & !
    IONO_NORTH_dSigmaThPs_dPsi,                  &
    IONO_NORTH_dSigmaPsPs_dPsi,                  & !
    IONO_SOUTH_dSigmaThTh_dTheta,                &
    IONO_SOUTH_dSigmaThPs_dTheta,                & !
    IONO_SOUTH_dSigmaPsPs_dTheta,                &
    IONO_SOUTH_dSigmaThTh_dPsi,                  & !
    IONO_SOUTH_dSigmaThPs_dPsi,                  &
    IONO_SOUTH_dSigmaPsPs_dPsi,                  &
    SAVE_NORTH_SigmaH,SAVE_NORTH_SigmaP,         &
    SAVE_SOUTH_SigmaH,SAVE_SOUTH_SigmaP

  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::  &
    IONO_NORTH_PHI_BC,IONO_SOUTH_PHI_BC,         & !Magnetosphere bound pot
    IONO_NORTH_ETh_BC,IONO_NORTH_EPs_BC,         & !Magnetosphere bound E field
    IONO_SOUTH_ETh_BC,IONO_SOUTH_EPs_BC,         & !
    IONO_NORTH_UTh_BC,IONO_NORTH_UPs_BC,         & !Magnetosphere bound flow
    IONO_NORTH_UR_BC,                            & !Velocities
    IONO_SOUTH_UTh_BC,IONO_SOUTH_UPs_BC,         &
    IONO_SOUTH_UR_BC

  real, dimension(1:IONO_nTheta) :: dTheta_North, dTheta_South
  real, dimension(1:IONO_nPsi)   :: dPsi_North, dPsi_South

  real, dimension(1:IONO_nTheta,1:IONO_nPsi,1:IONO_nRadial) ::  &
    IONO_NORTH_UTh_velB0,IONO_NORTH_UPs_velB0,                & !Magnetosphere
    IONO_NORTH_UR_velB0,                                      & !B0 convection
    IONO_SOUTH_UTh_velB0,IONO_SOUTH_UPs_velB0,                & !velocities
    IONO_SOUTH_UR_velB0

  real, dimension(1:IONO_nTheta,1:IONO_nPsi,1:IONO_NGMAX) ::  &
    IONO_U,IONO_RES,IONO_RHO,IONO_RHS,                        & !Multigrid solution residuals and RHS
    IONO_X,IONO_Y,                                            & !Multigrid solution coordinates
    IONO_S0,IONO_SH,IONO_SP,                                  & !Multigrid conductances
    IONO_Sxx,IONO_Sxy,IONO_Syy,                               & !
    IONO_dSxxdx,IONO_dSxydx,IONO_dSyydx,                      & !
    IONO_dSxxdy,IONO_dSxydy,IONO_dSyydy

  real, dimension(:), allocatable ::   &
    MAG_NORTH_X,MAG_NORTH_Y,MAG_NORTH_Z,MAG_NORTH_R,          & !Magnetospheric coordinates
    MAG_NORTH_Theta,MAG_NORTH_Psi,                            & !
    MAG_SOUTH_X,MAG_SOUTH_Y,MAG_SOUTH_Z,MAG_SOUTH_R,          & !
    MAG_SOUTH_Theta,MAG_SOUTH_Psi,                            & !
    MAG_NORTH_JR,MAG_NORTH_Jx,MAG_NORTH_Jy,MAG_NORTH_Jz,      & !Magnetospheric current
    MAG_SOUTH_JR,MAG_SOUTH_Jx,MAG_SOUTH_Jy,MAG_SOUTH_Jz

  !\
  ! Magnetosphere inner boundary current solution variable definitions.
  !/
  integer :: nMagBndPts_North, &
             nMagBndPts_South
  integer, dimension(:), allocatable :: nMagBndPts_North_PE, &
                                        nMagBndPts_South_PE
  real, dimension(:), allocatable ::   &
    Xmag_North,Ymag_North,Zmag_North,         & !Magnetospheric coordinates
    Xmag_South,Ymag_South,Zmag_South,         & !
    JXmag_North,JYmag_North,JZmag_North,      & !Magnetospheric current
    JXmag_South,JYmag_South,JZmag_South

end module ModIonosphere
