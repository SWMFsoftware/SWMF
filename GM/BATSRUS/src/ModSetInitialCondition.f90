!^CFG COPYRIGHT UM
subroutine set_ics

  use ModMain
  use ModAdvance
  use ModGeometry, ONLY : x2,y2,z2,x_BLK,y_BLK,z_BLK,R_BLK,true_cell
  use ModIO, ONLY : restart
  use ModImplicit,ONLY: UsePartImplicit             !^CFG IF IMPLICIT
  use ModPhysics
  use ModNumConst
  use ModUser, ONLY: user_set_ics

  implicit none

  real   :: SinSlope, CosSlope, Rot_II(2,2)
  integer:: i, j, k, iBlock, iVar
  !----------------------------------------------------------------------------
  iBlock = GlobalBlk

  B0xCell_BLK(:,:,:,iBlock) = 0.00
  B0yCell_BLK(:,:,:,iBlock) = 0.00
  B0zCell_BLK(:,:,:,iBlock) = 0.00
  B0xFace_x_BLK(:,:,:,iBlock) = 0.00
  B0yFace_x_BLK(:,:,:,iBlock) = 0.00
  B0zFace_x_BLK(:,:,:,iBlock) = 0.00
  B0xFace_y_BLK(:,:,:,iBlock) = 0.00
  B0yFace_y_BLK(:,:,:,iBlock) = 0.00
  B0zFace_y_BLK(:,:,:,iBlock) = 0.00
  B0xFace_z_BLK(:,:,:,iBlock) = 0.00
  B0yFace_z_BLK(:,:,:,iBlock) = 0.00
  B0zFace_z_BLK(:,:,:,iBlock) = 0.00

  time_BLK(:,:,:,iBlock) = 0.00

  fbody_x_BLK(:,:,:,iBlock) = 0.00
  fbody_y_BLK(:,:,:,iBlock) = 0.00
  fbody_z_BLK(:,:,:,iBlock) = 0.00

  Flux_VX = cZero
  Flux_VY = cZero
  Flux_VZ = cZero

  call init_conservative_facefluxes(iBlock)

  if(unusedBLK(iBlock))then  
     do iVar=1,nVar
        State_VGB(iVar,:,:,:,iBlock) = DefaultState_V(iVar)
     end do
  else
     !\
     ! If used, initialize solution variables and parameters.
     !/
     call set_b0(iBlock)
     call body_force_averages

     if(.not.restart)then

        if(UseShockTube)then
           ! Calculate sin and cos from the tangent = ShockSlope
           SinSlope=ShockSlope/sqrt(cOne+ShockSlope**2)
           CosSlope=      cOne/sqrt(cOne+ShockSlope**2)
           ! Set rotational matrix
           Rot_II = reshape( (/CosSlope, SinSlope, CosSlope, -SinSlope/), &
                (/2,2/) )
        end if

        ! Loop through all the cells
        do k=-1,nK+2; do j=-1,nJ+2; do i=-1,nI+2
           if(.not.true_cell(i,j,k,iBlock))then
              State_VGB(:,i,j,k,iBlock)   = CellState_VI(:,body1_)
           elseif(.not.UseShockTube)then
              State_VGB(:,i,j,k,iBlock)   = CellState_VI(:,1)
           else
              if(x_BLK(i,j,k,iBlock) < -ShockSlope*y_BLK(i,j,k,iBlock))then
                 ! Set all variables first
                 State_VGB(:,i,j,k,iBlock)   = ShockLeftState_V

                 ! Rotate vector variables
                 State_VGB(Ux_:Uy_,i,j,k,iBlock) = &
                      matmul(Rot_II,ShockLeftState_V(Ux_:Uy_))
                 State_VGB(Bx_:By_,i,j,k,iBlock) = &
                      matmul(Rot_II,ShockLeftState_V(Bx_:By_))
              else
                 ! Set all variables first
                 State_VGB(:,i,j,k,iBlock)   = ShockRightState_V
                 ! Set vector variables
                 State_VGB(Ux_:Uy_,i,j,k,iBlock) = &
                      matmul(Rot_II,ShockRightState_V(Ux_:Uy_))
                 State_VGB(Bx_:By_,i,j,k,iBlock) = &
                      matmul(Rot_II,ShockRightState_V(Bx_:By_))
              end if
              ! Convert velocity to momentum
              State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = &
                   State_VGB(Rho_,i,j,k,iBlock)*State_VGB(Ux_:Uz_,i,j,k,iBlock)
           end if

        end do; end do; end do
        !\
        ! Initialize solution quantities
        !/
        if(UseConstrainB)call constrain_ics(iBlock) !^CFG IF CONSTRAINB

        if(UseUserICs) call user_set_ics

     end if ! not restart

  end if ! unusedBLK
  !\
  ! Compute energy from set values above.
  !/
  call calc_energy(iBlock)

end subroutine set_ics

