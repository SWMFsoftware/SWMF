!^CFG COPYRIGHT UM
Module ModAdvance
  use ModSize
  use ModVarIndexes
  use ModIO,         ONLY: iUnitOut, write_prefix
  use ModProcMH,     ONLY: iProc
  
  implicit none
  save

  ! Logical parameter indicating static vs. dynamic allocation
  logical, parameter :: IsDynamicAdvance = .false.

  ! Update check parameters
  logical :: UseUpdateCheck
  real :: percent_max_rho(2), percent_max_p(2)


  ! Numerical flux type and optimization
  character (len=10) :: FluxType

  !\
  ! Conservative/Non-conservative parameters
  !/
  logical :: UseNonConservative

  ! Number and type of criteria
  integer :: nConservCrit
  character (len=10), dimension(:), allocatable :: TypeConservCrit_I

  ! Geometrical parameters
  real    :: rConserv, xParabolaConserv, yParabolaConserv 

  ! Physics based parameters (to locate shocks)
  real    :: pCoeffConserv, GradPCoeffConserv

  ! Cells selected to be updated with conservative equations
  logical, dimension(:,:,:,:), allocatable :: IsConserv_CB


  !\
  ! Block cell-centered MHD solution
  !/
  real, dimension(nVar,1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn,nBLK) :: &
       State_VGB
  real, dimension(1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn,nBLK) :: &
       E_BLK
  common /MHDSolnBlock1/ &
       State_VGB,E_BLK 
  !\
  ! Block cell-centered MHD solution old state
  !/
  real,  dimension(nVar,1:nI, 1:nJ, 1:nK,nBLK) :: &
       StateOld_VCB
  real,  dimension(1:nI, 1:nJ, 1:nK,nBLK) :: &
       E_o_BLK,time_BLK
  common /MHDSolnBlock2/ &
       StateOld_VCB,E_o_BLK,time_BLK  

  !\
  ! Block cell-centered intrinsic magnetic field, time, and temporary storage
  !/
  real,  dimension(1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn,nBLK) :: &
       B0xCell_BLK,B0yCell_BLK,B0zCell_BLK, &
       tmp1_BLK, tmp2_BLK

  ! Array for storing dB0/dt derivatives
  real, allocatable :: Db0Dt_CDB(:,:,:,:,:)


  !\
  ! Block cell-centered body forces & heat sources
  !/
  real, dimension(1:nI, 1:nJ, 1:nK,nBLK) :: &
       fbody_x_BLK,fbody_y_BLK,fbody_z_BLK,qheat_BLK

  !\
  ! Local cell-centered source terms and divB.
  !/
  real, dimension(Energy_,1:nI,1:nJ,1:nK)::Source_VC
  real, dimension(1:nI,1:nJ,1:nK) :: Theat0
  real,dimension(0:nI+1,0:nJ+1,0:nK+1,nBLK):: DivB1_GB
  real, dimension( 0:nI+1, 0:nJ+1, 0:nK+1) :: &
       gradX_rho,gradX_Ux,gradX_Uy,gradX_Uz,gradX_Bx,gradX_By,gradX_Bz,gradX_p,&
       gradX_VAR,&
       gradY_rho,gradY_Ux,gradY_Uy,gradY_Uz,gradY_Bx,gradY_By,gradY_Bz,gradY_p,&
       gradY_VAR,&
       gradZ_rho,gradZ_Ux,gradZ_Uy,gradZ_Uz,gradZ_Bx,gradZ_By,gradZ_Bz,gradZ_p,&
       gradZ_VAR

  !\
  ! Block face-centered intrinsic magnetic field array definitions.
  !/
  real,  dimension(2-gcn:nI+gcn,0:nJ+1,0:nK+1,nBLK) :: &
       B0xFace_x_BLK,B0yFace_x_BLK,B0zFace_x_BLK 
  real,  dimension(0:nI+1,2-gcn:nJ+gcn,0:nK+1,nBLK) :: &
       B0xFace_y_BLK,B0yFace_y_BLK,B0zFace_y_BLK
  real,  dimension(0:nI+1,0:nJ+1,2-gcn:nK+gcn,nBLK) :: &
       B0xFace_z_BLK,B0yFace_z_BLK,B0zFace_z_BLK
!^CFG IF NOT SIMPLE BEGIN
  real,dimension(3,3,1:nI,1:nJ,1:nK,nBLK) :: B0SourceMatrix_DDCB
!^CFG END SIMPLE

  !\
  ! X Face local MHD solution array definitions.
  !/
  integer, parameter :: nFaceValueVars = nVar
  real, dimension(nFaceValueVars,2-gcn:nI+gcn,0:nJ+1,0:nK+1) ::     &
       LeftState_VX,      &  ! Face Left  X
       RightState_VX          ! Face Right X 
  real, dimension(2-gcn:nI+gcn,0:nJ+1,0:nK+1) ::     &
       EDotFA_X,                              & !^CFG IF BORISCORR      
       VdtFace_x,UDotFA_X    ! V/dt Face  X
  real, dimension(Energy_,0:nI+1,2-gcn:nJ+gcn,0:nK+1) ::     &
       Flux_VX         ! Face Flux  X 
  common /MHDSolnFaceX/   &
       LeftState_VX,      &
       RightState_VX,      &
       Flux_VX,      & 
       EDotFA_X,                              & !^CFG IF BORISCORR      
       VdtFace_x,UDotFA_X 

  !\
  ! Y Face local MHD solution array definitions.
  !/
  real, dimension(nFaceValueVars,0:nI+1,2-gcn:nJ+gcn,0:nK+1) ::     &
       LeftState_VY,      &  ! Face Left  Y
       RightState_VY          ! Face Right Y

  real, dimension(0:nI+1,2-gcn:nJ+gcn,0:nK+1) ::     &
       EDotFA_Y,                              & !^CFG IF BORISCORR      
       VdtFace_y,UDotFA_Y    ! V/dt Face  Y
  real, dimension(Energy_,0:nI+1,2-gcn:nJ+gcn,0:nK+1) ::     &
       Flux_VY         ! Face Flux  Y 

  common /MHDSolnFaceY/   &
       LeftState_VY,      &
       RightState_VY,      &
       Flux_VY,      &
       EDotFA_Y,                              & !^CFG IF BORISCORR   
       VdtFace_y,UDotFA_Y                                          

  !\
  ! Z Face local MHD solution array definitions.
  !/
  real, dimension(nFaceValueVars,0:nI+1,0:nJ+1,2-gcn:nK+gcn) ::     &
       LeftState_VZ,      &  ! Face Left  Z
       RightState_VZ          ! Face Right Z

  real, dimension(0:nI+1,0:nJ+1,2-gcn:nK+gcn) ::     &
       EDotFA_Z,                              & !^CFG IF BORISCORR
       VdtFace_z,UDotFA_Z    !V/dt Face Z
  real, dimension(Energy_,0:nI+1,2-gcn:nJ+gcn,0:nK+1) ::     &
       Flux_VZ         ! Face Flux  Z 
  common /MHDSolnFaceZ/   &
       LeftState_VZ,      &
       RightState_VZ,      & 
       Flux_VZ,      &
       EDotFA_Z,                              & !^CFG IF BORISCORR
       VdtFace_z,UDotFA_Z

  !\
  ! The number of the face variables, which are corrected at the
  ! resolution changes
  !/

  !\
  !  Face conservative or corrected flux.
  !/
  real, dimension(nCorrectedFaceValues,1:nJ,1:nK,1:2,nBLK) :: &
          CorrectedFlux_VXB
  real, dimension(nCorrectedFaceValues,1:nI,1:nK,1:2,nBLK) :: &
          CorrectedFlux_VYB
  real, dimension(nCorrectedFaceValues,1:nI,1:nJ,1:2,nBLK) :: &
          CorrectedFlux_VZB

contains

  !============================================================================

  subroutine init_mod_advance

    if(IsDynamicAdvance .and. iProc==0)then
       call write_prefix
       write(iUnitOut,'(a)') 'init_mod_advance allocated arrays'
    end if

  end subroutine init_mod_advance

  !============================================================================

  subroutine clean_mod_advance

    if(IsDynamicAdvance .and. iProc==0)then
       call write_prefix
       write(iUnitOut,'(a)') 'clean_mod_advance deallocated arrays'
    end if

  end subroutine clean_mod_advance

  !============================================================================

end Module ModAdvance
