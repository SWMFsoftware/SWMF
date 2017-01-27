!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine set_ics(iBlock)

  use ModMain
  use ModAdvance
  use ModB0, ONLY: B0_DGB, set_b0_cell, subtract_b0
  use ModGeometry, ONLY: true_cell, R2_BLK
  use ModIO, ONLY : restart
  use ModPhysics
  use ModUserInterface ! user_set_ics
  use ModMultiFluid
  use ModEnergy, ONLY: calc_energy_ghost
  use ModConserveFlux, ONLY: init_cons_flux
  use ModRestartFile, ONLY: UseRestartWithFullB
  use BATL_lib, ONLY: Xyz_DGB
  use ModBoundaryGeometry, ONLY: iBoundary_GB

  implicit none

  integer, intent(in) :: iBlock

  real   :: SinSlope, CosSlope, Rot_II(2,2)
  real   :: ShockLeft_V(nVar), ShockRight_V(nVar)
  integer:: i, j, k, iVar, iBoundary

  character(len=*), parameter:: NameSub = 'set_ics'
  logical :: DoTest, DoTestMe
  !----------------------------------------------------------------------------

  if(iProc == ProcTest .and. iBlock == BlkTest)then
     call set_oktest(NameSub, DoTest, DoTestMe)
  else
     DoTest = .false.; DoTestMe = .false.
  end if

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
              ShockLeft_V(iUx:iUz)  = ShockLeftState_V(iUx:iUz) *Io2No_V(UnitU_)
              ShockRight_V(iUx:iUz) = ShockRightState_V(iUx:iUz)*Io2No_V(UnitU_)
           end do

        end if

        ! Loop through all the cells
        do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
           if(.not.true_cell(i,j,k,iBlock))then
              iBoundary = iBoundary_GB(i,j,k,iBlock)
              State_VGB(1:nVar,i,j,k,iBlock) = CellState_VI(1:nVar,iBoundary)
           elseif(.not.UseShockTube)then
              State_VGB(:,i,j,k,iBlock)   = CellState_VI(:,1)  
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

        if(iSignRotationIC /= 0) call add_rotational_velocity(iSignRotationIC, iBlock)

     end if ! not restart

  end if ! Unused_B

  !\
  ! Compute energy from set values above.
  !/
  call calc_energy_ghost(iBlock)

  if(DoTestMe)write(*,*) &
       NameSub, 'State(test)=',State_VGB(:,iTest,jTest,kTest,BlkTest)

end subroutine set_ics
