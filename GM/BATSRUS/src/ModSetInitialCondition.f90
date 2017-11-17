!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModSetInitialCondition

  use BATL_lib, ONLY: &
       test_start, test_stop, iBlockTest, iTest, jTest, kTest

  implicit none

  private ! except

  public:: set_initial_condition   ! set initial condition for one block
  public:: add_rotational_velocity ! transform between rotating/inertial frames

contains
  !============================================================================

  subroutine set_initial_condition(iBlock)

    use ModMain
    use ModAdvance
    use ModB0, ONLY: B0_DGB, set_b0_cell, subtract_b0
    use ModGeometry, ONLY: true_cell
    use ModIO, ONLY : restart
    use ModPhysics
    use ModUserInterface ! user_set_ics
    use ModConstrainDivB, ONLY: constrain_ics
    use ModMultiFluid
    use ModEnergy, ONLY: calc_energy_ghost
    use ModConserveFlux, ONLY: init_cons_flux
    use ModRestartFile, ONLY: UseRestartWithFullB
    use ModBoundaryGeometry, ONLY: iBoundary_GB
    use BATL_lib, ONLY: Xyz_DGB

    integer, intent(in) :: iBlock

    real   :: SinSlope, CosSlope, Rot_II(2,2)
    real   :: ShockLeft_V(nVar), ShockRight_V(nVar)
    integer:: i, j, k, iVar, iBoundary

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'set_initial_condition'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    time_BLK(:,:,:,iBlock) = 0.00

    Flux_VX = 0.0
    Flux_VY = 0.0
    Flux_VZ = 0.0

    call init_cons_flux(iBlock)

    if(Unused_B(iBlock))then
       do iVar = 1, nVar
          State_VGB(iVar,:,:,:,iBlock) = DefaultState_V(iVar)
       end do
    else

       !\
       ! If used, initialize solution variables and parameters.
       !/
       if(UseB0) call set_b0_cell(iBlock)

       ! Subtract B0 from Full B0+B1 from restart to obtain B1
       if(UseB0 .and. restart .and. UseRestartWithFullB) call subtract_b0(iBlock)

       if(.not.restart)then

          if(UseShockTube)then
             ! Calculate sin and cos from the tangent = ShockSlope
             SinSlope=ShockSlope/sqrt(1.0 + ShockSlope**2)
             CosSlope=       1.0/sqrt(1.0 + ShockSlope**2)
             ! Set rotational matrix
             Rot_II = reshape( (/CosSlope, SinSlope, -SinSlope, CosSlope/), &
                  (/2,2/) )
             ! calculate normalized left and right states
             ShockLeft_V  = ShockLeftState_V /UnitUser_V(1:nVar)
             ShockRight_V = ShockRightState_V/UnitUser_V(1:nVar)

             ! fix the units for the velocities
             do iFluid = 1, nFluid
                call select_fluid
                ShockLeft_V(iUx:iUz) = ShockLeftState_V(iUx:iUz) *Io2No_V(UnitU_)
                ShockRight_V(iUx:iUz)= ShockRightState_V(iUx:iUz)*Io2No_V(UnitU_)
             end do

          end if

          ! Loop through all the cells
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             if(.not.true_cell(i,j,k,iBlock))then
                iBoundary = iBoundary_GB(i,j,k,iBlock)

                State_VGB(1:nVar,i,j,k,iBlock) = FaceState_VI(1:nVar,iBoundary)

                ! Convert velocity to momentum
                do iFluid = 1, nFluid
                   call select_fluid
                   State_VGB(iRhoUx,i,j,k,iBlock) = &
                        FaceState_VI(iUx,iBoundary)*FaceState_VI(iRho,iBoundary)
                   State_VGB(iRhoUy,i,j,k,iBlock) = &
                        FaceState_VI(iUy,iBoundary)*FaceState_VI(iRho,iBoundary)
                   State_VGB(iRhoUz,i,j,k,iBlock) = &
                        FaceState_VI(iUz,iBoundary)*FaceState_VI(iRho,iBoundary)
                end do

             elseif(.not.UseShockTube)then
                State_VGB(:,i,j,k,iBlock)   = CellState_VI(:,Coord1MaxBc_)
             else
                if( (Xyz_DGB(x_,i,j,k,iBlock)-ShockPosition) &
                     < -ShockSlope*Xyz_DGB(y_,i,j,k,iBlock))then
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
                State_VGB(Bx_:Bz_,i,j,k,iBlock) = &
                     State_VGB(Bx_:Bz_,i,j,k,iBlock) - B0_DGB(:,i,j,k,iBlock)
             end if

          end do; end do; end do

          if(UseConstrainB)call constrain_ics(iBlock)

          if(UseUserICs) call user_set_ics(iBlock)

          if(iSignRotationIC /= 0) &
               call add_rotational_velocity(iSignRotationIC, iBlock)

       end if ! not restart

    end if ! Unused_B

    !\
    ! Compute energy from set values above.
    !/
    call calc_energy_ghost(iBlock)

    if(DoTest)write(*,*) &
         NameSub, 'State(test)=',State_VGB(:,iTest,jTest,kTest,iBlockTest)

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine set_initial_condition
  !============================================================================

  subroutine add_rotational_velocity(iSign, iBlockIn)

    use ModSize,     ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nBlock, x_, y_
    use ModMain,     ONLY: Unused_B
    use ModAdvance,  ONLY: State_VGB
    use ModGeometry, ONLY: true_cell
    use ModPhysics,  ONLY: OmegaBody
    use ModMultiFluid, ONLY: iRho_I, iRhoUx_I, iRhoUy_I
    use BATL_lib,    ONLY: Xyz_DGB

    integer, intent(in):: iSign
    integer, optional, intent(in):: iBlockIn

    ! Transform velocities between inertial and rotating frames
    ! where Omega is the angular velocity of the rotating frame
    ! Since Omega = (0,0,OmegaBody)
    ! ux = ux - iSign*OmegaBody*y
    ! uy = uy + iSign*OmegaBody*x
    ! iSign=+1: from rotating to inertial frame
    ! iSign=-1: from inertial to rotating frame
    !
    ! If iBlockIn is present, do that block, otherwise do all blocks.

    integer :: i, j, k, iBlock, iBlockFirst, iBlockLast
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'add_rotational_velocity'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    if(present(iBlockIn))then
       iBlockFirst = iBlockIn; iBlockLast = iBlockIn
    else
       iBlockFirst = 1; iBlockLast = nBlock
    end if

    do iBlock = iBlockFirst, iBlockLast
       if(Unused_B(iBlock))CYCLE
       do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
          if(.not.true_Cell(i,j,k,iBlock)) CYCLE
          State_VGB(iRhoUx_I,i,j,k,iBlock) = State_VGB(iRhoUx_I,i,j,k,iBlock) &
               - iSign*State_VGB(iRho_I,i,j,k,iBlock) &
	       *OmegaBody*Xyz_DGB(y_,i,j,k,iBlock)

          State_VGB(iRhoUy_I,i,j,k,iBlock) = State_VGB(iRhoUy_I,i,j,k,iBlock) &
	       + iSign*State_VGB(iRho_I,i,j,k,iBlock) &
	       *OmegaBody*Xyz_DGB(x_,i,j,k,iBlock)

       end do; end do; end do
    end do

    call test_stop(NameSub, DoTest)
  end subroutine add_rotational_velocity
  !============================================================================

end module ModSetInitialCondition
!==============================================================================
