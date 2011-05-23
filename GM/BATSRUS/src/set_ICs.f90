!^CFG COPYRIGHT UM
subroutine set_ics

  use ModMain
  use ModAdvance
  use ModGeometry, ONLY : x2,y2,z2,x_BLK,y_BLK,z_BLK,R_BLK,true_cell
  use ModIO, ONLY : restart
  use ModImplicit,ONLY: UsePartImplicit             !^CFG IF IMPLICIT
  use ModPhysics
  use ModUser, ONLY: user_set_ics
  use ModMultiFluid
  use ModEnergy, ONLY: calc_energy_ghost
  use ModConserveFlux, ONLY: init_cons_flux

  implicit none

  real   :: SinSlope, CosSlope, Rot_II(2,2)
  real   :: ShockLeft_V(nVar), ShockRight_V(nVar)
  integer:: i, j, k, iBlock, iVar

  character(len=*), parameter:: NameSub = 'set_ics'
  logical :: DoTest, DoTestMe
  !----------------------------------------------------------------------------
  iBlock = GlobalBlk

  if(iProc == ProcTest .and. iBlock == BlkTest)then
     call set_oktest(NameSub, DoTest, DoTestMe)
  else
     DoTest = .false.; DoTestMe = .false.
  end if

  if(UseB0)then
     B0_DGB(:,:,:,:,iBlock) = 0.00
     B0ResChange_DXSB(:,:,:,:,iBlock) = 0.00
     B0ResChange_DYSB(:,:,:,:,iBlock) = 0.00
     B0ResChange_DZSB(:,:,:,:,iBlock) = 0.00
  end if

  time_BLK(:,:,:,iBlock) = 0.00

  if(UseGravity .or. UseRotatingFrame)then
     fbody_x_BLK(:,:,:,iBlock) = 0.00
     fbody_y_BLK(:,:,:,iBlock) = 0.00
     fbody_z_BLK(:,:,:,iBlock) = 0.00
  end if

  Flux_VX = 0.0
  Flux_VY = 0.0
  Flux_VZ = 0.0

  call init_cons_flux(iBlock)

  if(unusedBLK(iBlock))then  
     do iVar=1,nVar
        State_VGB(iVar,:,:,:,iBlock) = DefaultState_V(iVar)
     end do
  else
     !\
     ! If used, initialize solution variables and parameters.
     !/
     if(UseB0)call set_b0(iBlock)
     call body_force_averages

     if(.not.restart)then

        if(UseShockTube)then
           ! Calculate sin and cos from the tangent = ShockSlope
           SinSlope=ShockSlope/sqrt(cOne+ShockSlope**2)
           CosSlope=      cOne/sqrt(cOne+ShockSlope**2)
           ! Set rotational matrix
           Rot_II = reshape( (/CosSlope, SinSlope, -SinSlope, CosSlope/), &
                (/2,2/) )
           ! calculate normalized left and right states
           ShockLeft_V  = ShockLeftState_V /UnitUser_V(1:nVar)
           ShockRight_V = ShockRightState_V/UnitUser_V(1:nVar)

           ! fix the units for the velocities
           do iFluid = 1, nFluid
              call select_fluid
              ShockLeft_V(iUx:iUz)  = ShockLeftState_V(iUx:iUz) *Io2No_V(Ux_)
              ShockRight_V(iUx:iUz) = ShockRightState_V(iUx:iUz)*Io2No_V(Ux_)
           end do

        end if

        ! Loop through all the cells
        do k=-1,nK+2; do j=-1,nJ+2; do i=-1,nI+2
           if(.not.true_cell(i,j,k,iBlock))then
              State_VGB(:,i,j,k,iBlock)   = CellState_VI(:,body1_)
           elseif(.not.UseShockTube)then
              State_VGB(:,i,j,k,iBlock)   = CellState_VI(:,1)
           else
              if( (x_BLK(i,j,k,iBlock)-ShockPosition) &
                   < -ShockSlope*y_BLK(i,j,k,iBlock))then
                 ! Set all variables first
                 State_VGB(:,i,j,k,iBlock)   = ShockLeft_V

                 ! Rotate vector variables
                 do iFluid = 1, nFluid
                    call select_fluid
                    State_VGB(iUx:iUy,i,j,k,iBlock) = &
                         matmul(Rot_II,ShockLeft_V(iUx:iUy))
                 end do
                 if(UseB) State_VGB(Bx_:By_,i,j,k,iBlock) = &
                      matmul(Rot_II,ShockLeft_V(Bx_:By_))
              else
                 ! Set all variables first
                 State_VGB(:,i,j,k,iBlock)   = ShockRight_V
                 ! Set vector variables
                 do iFluid = 1, nFluid
                    call select_fluid
                    State_VGB(iUx:iUy,i,j,k,iBlock) = &
                         matmul(Rot_II,ShockRight_V(iUx:iUy))
                 end do
                 if(UseB) State_VGB(Bx_:By_,i,j,k,iBlock) = &
                      matmul(Rot_II,ShockRight_V(Bx_:By_))
              end if
              ! Convert velocity to momentum
              do iFluid = 1, nFluid
                 call select_fluid
                 State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = &
                      State_VGB(iRho,i,j,k,iBlock) * &
                      State_VGB(iUx:iUz,i,j,k,iBlock)
              end do
              if(.not.UseB0)CYCLE
              ! Remove B0 from B (if any)
              State_VGB(Bx_:Bz_,i,j,k,iBlock)=State_VGB(Bx_:Bz_,i,j,k,iBlock) - &
                   B0_DGB(:,i,j,k,iBlock)
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
  call calc_energy_ghost(iBlock)

  if(DoTestMe)write(*,*) &
       NameSub, 'State(test)=',State_VGB(:,iTest,jTest,kTest,BlkTest)

end subroutine set_ics
