!^CFG COPYRIGHT UM
!
!-------------------------------------------------------------------------
Module ModMappingParam
  use ModNumConst

  real,parameter::IONO_PI=cPi
  real, parameter :: IONO_TOLER = 5.0e-05
  real :: IONO_Radius_Mag_Boundary, &
       IONO_Radius, IONO_Height, rMap, DipoleSign

end Module ModMappingParam

!============================================================================

subroutine set_mapping_param

  use ModPhysics,  ONLY: unitSI_x,Rbody
  use CON_physics, ONLY: get_planet
  use ModMappingParam
  implicit none
  real::DipoleStrength

  call get_planet(&
       RadiusPlanetOut      = IONO_Radius,&
       IonosphereHeightOut  = IONO_Height,&
       DipoleStrengthOut    = DipoleStrength)

  ! final radius is fixed, store it
  rMap=(IONO_Radius+ IONO_Height)/UnitSI_x
  IONO_Radius_Mag_Boundary = Rbody*UnitSI_x
  DipoleSign=sign(cOne,DipoleStrength)

end subroutine set_mapping_param

!============================================================================

module ModInnerBC

  use ModMappingParam
  implicit none
  save

  !\
  ! Ionosphere array parameters
  !/
  integer :: IONO_nTheta = -1
  integer :: IONO_nPsi   = -1

  !\
  ! Ionosphere solution array definitions
  !/
  real, dimension(:,:), allocatable ::     &
       IONO_NORTH_Theta,IONO_NORTH_Psi,             & !
       IONO_SOUTH_Theta,IONO_SOUTH_Psi


  real, dimension(:,:), allocatable ::     &
       IONO_NORTH_PHI_BC,IONO_SOUTH_PHI_BC,         & !Magnetospheric bound pot
       IONO_NORTH_ETh_BC,IONO_NORTH_EPs_BC,         & !Magnetospheric bound
       IONO_SOUTH_ETh_BC,IONO_SOUTH_EPs_BC,         & !Electric fields
       IONO_NORTH_UTh_BC,IONO_NORTH_UPs_BC,         & !Magnetosphere bound flow
       IONO_NORTH_UR_BC,                            & !Velocities
       IONO_SOUTH_UTh_BC,IONO_SOUTH_UPs_BC,         &
       IONO_SOUTH_UR_BC

end module ModInnerBC

!============================================================================

subroutine init_inner_bc_arrays(iSize,jSize)

  use CON_coupler
  use ModInnerBc,ONLY:&
       IONO_nTheta, IONO_nPsi,&            
       IONO_NORTH_Theta, IONO_NORTH_Psi, & 
       IONO_SOUTH_Theta, IONO_SOUTH_Psi, &   
       IONO_NORTH_PHI_BC,IONO_SOUTH_PHI_BC,         & 
       IONO_NORTH_ETh_BC,IONO_NORTH_EPs_BC,& 
       IONO_SOUTH_ETh_BC,IONO_SOUTH_EPs_BC,&                  
       IONO_NORTH_UR_BC, &
       IONO_NORTH_UTh_BC, &
       IONO_NORTH_UPs_BC, &
       IONO_SOUTH_UR_BC, &
       IONO_SOUTH_UTh_BC, &
       IONO_SOUTH_UPs_BC

  implicit none

  integer, intent(in) :: iSize,jSize
  character(len=*), parameter :: NameSub='GM_allocate_iono_arrays'
  integer :: nCells_D(2),iError, i
  ! Allocate and calculate coordinates
  
  nCells_D=ncells_decomposition_d(IE_)
  if(  iSize /= nCells_D(1)+1 .or. &
       jSize /= nCells_D(2)+1 ) then
     
     write(*,*)NameSub//' grid sizes do not agree iSize,jSize,nCells=',&
          iSize,jSize,nCells_D(1:2)
     call CON_stop(NameSub//' SWMF_ERROR')
  end if
  
  IONO_nTheta = iSize
  IONO_nPsi   = jSize
  
  allocate(&
       IONO_NORTH_Theta(iSize,jSize), IONO_NORTH_Psi(iSize,jSize), &
       IONO_NORTH_Phi_BC(iSize,jSize),                             &
       IONO_SOUTH_Theta(iSize,jSize), IONO_SOUTH_Psi(iSize,jSize), &
       IONO_SOUTH_Phi_BC(iSize,jSize),                             &
       IONO_NORTH_ETh_BC(IONO_nTheta,IONO_nPsi),                   &
       IONO_NORTH_EPs_BC(IONO_nTheta,IONO_nPsi),                   &
       IONO_NORTH_UR_BC(IONO_nTheta,IONO_nPsi),                    &
       IONO_NORTH_UTh_BC(IONO_nTheta,IONO_nPsi),                   &
       IONO_NORTH_UPs_BC(IONO_nTheta,IONO_nPsi),                   &
       IONO_SOUTH_ETh_BC(IONO_nTheta,IONO_nPsi),                   &
       IONO_SOUTH_EPs_BC(IONO_nTheta,IONO_nPsi),                   &
       IONO_SOUTH_UR_BC(IONO_nTheta,IONO_nPsi),                    &
       IONO_SOUTH_UTh_BC(IONO_nTheta,IONO_nPsi),                   &
       IONO_SOUTH_UPs_BC(IONO_nTheta,IONO_nPsi),                   &
       STAT = iError)
  call check_allocate(iError,NameSub//' IONO_NORTH/SOUTH_ThetaPsiPhi')
  IONO_NORTH_Phi_BC=cZero
  IONO_SOUTH_Phi_BC=cZero
  IONO_NORTH_ETh_BC=cZero
  IONO_NORTH_EPs_BC=cZero
  IONO_NORTH_UR_BC=cZero
  IONO_NORTH_UTh_BC=cZero
  IONO_NORTH_UPs_BC=cZero
  IONO_SOUTH_ETh_BC=cZero
  IONO_SOUTH_EPs_BC=cZero
  IONO_SOUTH_UR_BC=cZero
  IONO_SOUTH_UTh_BC=cZero
  IONO_SOUTH_UPs_BC=cZero
  
  do i=1,iSize
     IONO_NORTH_Theta(i,:) = Grid_C(IE_) % Coord1_I(i)
     IONO_SOUTH_Theta(i,:) = Grid_C(IE_) % Coord1_I(i+iSize-1)
     IONO_NORTH_Psi(i,:)   = Grid_C(IE_) % Coord2_I
     IONO_SOUTH_Psi(i,:)   = Grid_C(IE_) % Coord2_I
  end do

end subroutine init_inner_bc_arrays
