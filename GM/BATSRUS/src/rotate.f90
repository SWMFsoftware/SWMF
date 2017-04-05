!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!This code is a copyright protected software (c) 2002- University of Michigan
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

  real, save:: Omega_D(3)
  logical   :: IsUninitialized = .true.

  !------------------------------------------------------------------------
  select case(TypeCoordSystem)
  case('HGI')
     ! In the HGI system the Solar angular velocity vector points towards +Z
     Omega_D = (/ 0., 0., OmegaBody /)
  case('GSE')
     if(IsUninitialized)then
        call get_axes(Time_Simulation, RotAxisGseOut_D=Omega_D)
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
!==========================================================================
subroutine add_rotational_velocity(iSign, iBlock)

  use ModSize,     ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nBlock, x_, y_
  use ModMain,     ONLY: Unused_B
  use ModAdvance,  ONLY: State_VGB
  use ModGeometry, ONLY: true_cell
  use ModPhysics,  ONLY: OmegaBody
  use ModMultiFluid, ONLY: iRho_I, iRhoUx_I, iRhoUy_I
  use BATL_lib,    ONLY: Xyz_DGB

  implicit none

  integer, intent(in):: iSign
  integer, intent(in):: iBlock

  ! Transform velocities between inertial and rotating frames
  ! where Omega is the angular velocity of the rotating frame
  ! Since Omega = (0,0,OmegaBody)
  ! ux = ux - iSign*OmegaBody*y
  ! uy = uy + iSign*OmegaBody*x
  ! iSign=+1: from rotating to inertial frame
  ! iSign=-1: from inertial to rotating frame
  !
  ! If iBlock is 0 or negative, do all blocks in State_VGB,
  ! otherwise do only block iBlock

  integer :: i, j, k, iBlk, iBlockFirst, iBlockLast
  !---------------------------------------------------------------------------      
  if(iBlock > 0)then
     iBlockFirst = 1; iBlockLast = nBlock
  else
     iBlockFirst = iBlock; iBlockLast = iBlock
  end if

  do iBlk = iBlockFirst, iBlockLast
     if(Unused_B(iBlk))CYCLE
     do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
        if(.not.true_Cell(i,j,k,iBlk)) CYCLE
        State_VGB(iRhoUx_I,i,j,k,iBlk) = State_VGB(iRhoUx_I,i,j,k,iBlk) - &
             iSign*State_VGB(iRho_I,i,j,k,iBlk)*OmegaBody*Xyz_DGB(y_,i,j,k,iBlk)
        
        State_VGB(iRhoUy_I,i,j,k,iBlk) = State_VGB(iRhoUy_I,i,j,k,iBlk) + &
             iSign*State_VGB(iRho_I,i,j,k,iBlk)*OmegaBody*Xyz_DGB(x_,i,j,k,iBlk)

     end do; end do; end do
  end do
  
end subroutine add_rotational_velocity

