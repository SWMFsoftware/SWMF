!^CFG COPYRIGHT UM
subroutine calc_corotation_velocities(Xyz_D, uRot_D)
  !-------------------------------------------------------------------------
  !\
  ! This routine calculates cartesian corotation velocity uRot_D as a
  ! function of the cartesian coordinates Xyz_D
  !/
  !-------------------------------------------------------------------------

  use CON_axes,          ONLY: get_axes
  use ModCoordTransform, ONLY: cross_product
  use ModMain,           ONLY: Time_Simulation, TypeCoordSystem
  use ModPhysics,        ONLY: OmegaBody
  use ModNumConst
  implicit none

  real, intent(in) :: Xyz_D(3)
  real, intent(out):: uRot_D(3)

  real    :: Omega_D(3)
  logical :: IsUninitialized = .true.

  !------------------------------------------------------------------------
  select case(TypeCoordSystem)
  case('HGI')
     ! In the HGI system the Solar angular velocity vector points towards +Z
     Omega_D = (/ 0., 0., OmegaBody /)
  case('GSE')
     if(IsUninitialized)then
        call get_axes(Time_Simulation,RotAxisGseOut_D=Omega_D)
        Omega_D = OmegaBody * Omega_D
        IsUninitialized = .false.
     end if
  case('GSM')
     ! GSM system, Omega_D may be changing
     call get_axes(Time_Simulation,RotAxisGsmOut_D=Omega_D)
     Omega_D = OmegaBody*Omega_D
  end select

  ! The corotation velocity is u = Omega x R

  uRot_D = cross_product(Omega_D, Xyz_D)

end subroutine calc_corotation_velocities
!=============================================================================
subroutine transform_to_hgi

  ! Transform velocities from rotating frame to the HGI frame
  ! u' = u + Omega x R, 
  ! where Omega is the angular velocity of the rotating frame
  ! Since Omega = (0,0,OmegaBody)
  ! ux = ux - OmegaBody*y
  ! uy = uy + OmegaBody*x

  use ModMain,     ONLY: nI, nJ, nK, nBlock, UnusedBlk
  use ModAdvance,  ONLY: State_VGB, Rho_, RhoUx_, RhoUy_
  use ModGeometry, ONLY: x_BLK, y_BLK, true_cell
  use ModPhysics,  ONLY: OmegaBody
  implicit none
  integer :: i,j,k,iBlock
  !---------------------------------------------------------------------------
  
  do iBlock=1, nBlock
     if(UnusedBlk(iBlock))CYCLE
     do k=-1,nK+2; do j=-1,nJ+2; do i=-1,nI+2
        if(.not.True_Cell(i,j,k,iBlock)) CYCLE
        State_VGB(RhoUx_,i,j,k,iBlock) = State_VGB(RhoUx_,i,j,k,iBlock) - &
             State_VGB(Rho_,i,j,k,iBlock)*OmegaBody*y_BLK(i,j,k,iBlock)

        State_VGB(RhoUy_,i,j,k,iBlock) = State_VGB(RhoUy_,i,j,k,iBlock) + &
             State_VGB(Rho_,i,j,k,iBlock)*OmegaBody*x_BLK(i,j,k,iBlock)

     end do; end do; end do
  end do

end subroutine transform_to_hgi
