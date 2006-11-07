!^CFG COPYRIGHT UM
subroutine set_ICs
  use ModMain
  use ModAdvance
  use ModGeometry, ONLY : x2,y2,z2,x_BLK,y_BLK,z_BLK,R_BLK,true_cell
  use ModIO, ONLY : restart
  use ModImplicit,ONLY: UsePartImplicit             !^CFG IF IMPLICIT
  use ModPhysics
  use ModNumConst
  use ModUser, ONLY: user_set_ICs
  implicit none

  real :: Rmax, SinSlope, CosSlope
  real :: B4, dB4dx, zeta4, q4, epsi4, plobe, &
       XFace, YFace, ZFace
  integer::i,j,k,iVar

  B0xCell_BLK(:,:,:,globalBLK) = 0.00
  B0yCell_BLK(:,:,:,globalBLK) = 0.00
  B0zCell_BLK(:,:,:,globalBLK) = 0.00
  B0xFace_x_BLK(:,:,:,globalBLK) = 0.00
  B0yFace_x_BLK(:,:,:,globalBLK) = 0.00
  B0zFace_x_BLK(:,:,:,globalBLK) = 0.00
  B0xFace_y_BLK(:,:,:,globalBLK) = 0.00
  B0yFace_y_BLK(:,:,:,globalBLK) = 0.00
  B0zFace_y_BLK(:,:,:,globalBLK) = 0.00
  B0xFace_z_BLK(:,:,:,globalBLK) = 0.00
  B0yFace_z_BLK(:,:,:,globalBLK) = 0.00
  B0zFace_z_BLK(:,:,:,globalBLK) = 0.00

  time_BLK(:,:,:,globalBLK) = 0.00

  fbody_x_BLK(:,:,:,globalBLK) = 0.00
  fbody_y_BLK(:,:,:,globalBLK) = 0.00
  fbody_z_BLK(:,:,:,globalBLK) = 0.00

  Flux_VX = cZero
  Flux_VY = cZero
  Flux_VZ = cZero

  call init_conservative_facefluxes(globalBLK)

  if(unusedBLK(globalBLK))then  
     do iVar=1,nVar
        State_VGB(iVar,:,:,:,globalBLK)   = DefaultState_V(iVar)
     end do
  else


     !\
     ! If used, initialize solution variables and parameters.
     !/
     call set_b0(globalBLK)
     call body_force_averages

     if(.not.restart)then
        !\
        ! Initialize solution quantities
        !/
        do k=1-gcn,nK+gcn;do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
           if(true_cell(i,j,k,globalBLK))then
              State_VGB(:,i,j,k,globalBLK)   = CellState_VI(:,1)
           else
              State_VGB(:,i,j,k,globalBLK)   = CellState_VI(:,body1_)
           end if
        end do;end do;end do
     end if

     if(UseUserICs .and. .not.restart) call user_set_ICs

  end if ! unusedBLK
  ! Eliminate B1 around body if necessary
  if(.not.restart)then
     if(UseConstrainB)then  !^CFG IF CONSTRAINB BEGIN
        call constrain_ICs(globalBLK)
     else                   !^CFG END CONSTRAINB
        if(index(test_string,'DriftIn')>0)then
           where(x_BLK(:,:,:,globalBLK)<10.)
              State_VGB(Bx_,:,:,:,globalBLK)=0.0
              State_VGB(By_,:,:,:,globalBLK)=0.0
              State_VGB(Bz_,:,:,:,globalBLK)=0.0
              State_VGB(P_,:,:,:,globalBLK)=&
                   State_VGB(P_,:,:,:,globalBLK)+ &
                   0.5*(SW_Bx**2+SW_By**2+SW_Bz**2)
           elsewhere
              State_VGB(Bx_,:,:,:,globalBLK)=SW_Bx
              State_VGB(By_,:,:,:,globalBLK)=SW_By
              State_VGB(Bz_,:,:,:,globalBLK)=SW_Bz
           end where
        end if
     end if                 !^CFG IF CONSTRAINB
  end if
  !\
  ! Compute energy from set values above.
  !/
  call calc_energy(globalBLK)

end subroutine set_ICs

