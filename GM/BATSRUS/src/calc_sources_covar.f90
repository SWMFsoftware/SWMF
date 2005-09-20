!^CFG COPYRIGHT UM
!^CFG FILE NOT CARTESIAN

subroutine calc_sources
  use ModMain
  use ModProcMH
  use ModCovariant
  use ModGeometry, ONLY : VolumeInverse_I, x_BLK, y_BLK, z_BLK,&
       body_BLK, RMin_BLK, TypeGeometry
  use ModAdvance
  use ModNumConst
  use ModPhysics, ONLY:gm1
  implicit none

  integer :: i,j,k

  logical :: oktest,oktest_me

  real:: B1nJump,DivBInternal_C(1:nI,1:nJ,1:nK)!,DivU_C(1:nI,1:nJ,1:nK)

  ! Variables needed for Boris source terms also used for div(u)

  real ::  B1Cell_D(3)

  real,dimension(3):: FaceArea_D
  real:: VInvHalf,RhoInv
  integer::iVolumeCounter


  !---------------------------------------------------------------------------
  if(iProc==PROCtest .and. globalBLK==BLKtest)then
     call set_oktest('calc_sources',oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  end if

  Source_VC   = cZero

  iVolumeCounter= nIJK*(globalBLK-1)
  if(UseNonconservative)then
     do k=1,nK; do j=1,nJ; do i=1,nI
        iVolumeCounter = iVolumeCounter + 1
        Source_VC(P_,i,j,k) = -gm1*State_VGB(P_,i,j,k,globalBLK)*&
             VolumeInverse_I(iVolumeCounter)*&
             (UDotFA_X(i+1,j,k)-UDotFA_X(i,j,k)+&
             UDotFA_Y(i,j+1,k) -UDotFA_Y(i,j,k)+&
             UDotFA_Z(i,j,k+1) -UDotFA_Z(i,j,k))
     end do; end do; end do
  end if

  if(UseDivbSource)then

     iVolumeCounter= nIJK*(globalBLK-1)
     do k=1,nK; do j=1,nJ; do i=1,nI
        iVolumeCounter=iVolumeCounter+1
        VInvHalf=VolumeInverse_I(iVolumeCounter)*cHalf
        call calc_faceareaI(i,j,k,globalBLK,FaceArea_D)
        B1nJump =VInvHalf*&
             (FaceArea_D(1)*(RightState_VX(Bx_,i,j,k)-LeftState_VX(Bx_,i,j,k))+&
             FaceArea_D(2)*(RightState_VX(By_,i,j,k)-LeftState_VX(By_,i,j,k))+&
             FaceArea_D(3)*(RightState_VX(Bz_,i,j,k)-LeftState_VX(Bz_,i,j,k)))
        DivBInternal_C(i,j,k)=&
             -(FaceArea_D(1)*RightState_VX(Bx_,i,j,k)+&
             FaceArea_D(2)*RightState_VX(By_,i,j,k)+&
             FaceArea_D(3)*RightState_VX(Bz_,i,j,k))
        Source_VC(rhoUx_,i,j,k) = -B0xFace_x_BLK(i,j,k,globalBLK)*B1nJump
        Source_VC(rhoUy_,i,j,k) = -B0yFace_x_BLK(i,j,k,globalBLK)*B1nJump
        Source_VC(rhoUz_,i,j,k) = -B0zFace_x_BLK(i,j,k,globalBLK)*B1nJump
        DivB1_GB(i,j,k,globalBLK)  = B1nJump
        call calc_faceareaI(i+1,j,k,globalBLK,FaceArea_D)
        B1nJump =  VInvHalf*&
             (FaceArea_D(1)*(RightState_VX(Bx_,i+1,j,k)-LeftState_VX(Bx_,i+1,j,k))+&
             FaceArea_D(2)*(RightState_VX(By_,i+1,j,k)-LeftState_VX(By_,i+1,j,k))+&
             FaceArea_D(3)*(RightState_VX(Bz_,i+1,j,k)-LeftState_VX(Bz_,i+1,j,k)))
        DivBInternal_C(i,j,k)=DivBInternal_C(i,j,k)+&
             (FaceArea_D(1)*LeftState_VX(Bx_,i+1,j,k)+&
             FaceArea_D(2)*LeftState_VX(By_,i+1,j,k)+&
             FaceArea_D(3)*LeftState_VX(Bz_,i+1,j,k))
        Source_VC(rhoUx_,i,j,k) = Source_VC(rhoUx_,i,j,k)&
             -B0xFace_x_BLK(i+1,j,k,globalBLK)*B1nJump
        Source_VC(rhoUy_,i,j,k) = Source_VC(rhoUy_,i,j,k)&
             -B0yFace_x_BLK(i+1,j,k,globalBLK)*B1nJump
        Source_VC(rhoUz_,i,j,k) = Source_VC(rhoUz_,i,j,k)&
             -B0zFace_x_BLK(i+1,j,k,globalBLK)*B1nJump
        DivB1_GB(i,j,k,globalBLK)  = DivB1_GB(i,j,k,globalBLK)+B1nJump
     end do; end do; end do

     iVolumeCounter=nIJK*(globalBLK-1)
     do k=1,nK; do j=1,nJ; do i=1,nI 
        iVolumeCounter=iVolumeCounter+1
        VInvHalf=VolumeInverse_I(iVolumeCounter)*cHalf
        call calc_faceareaJ(i,j,k,globalBLK,FaceArea_D)
        B1nJump = VInvHalf*&
             (FaceArea_D(1)*(RightState_VY(Bx_,i,j,k)-LeftState_VY(Bx_,i,j,k))+&
             FaceArea_D(2)*(RightState_VY(By_,i,j,k)-LeftState_VY(By_,i,j,k))+&
             FaceArea_D(3)*(RightState_VY(Bz_,i,j,k)-LeftState_VY(Bz_,i,j,k)))
        DivBInternal_C(i,j,k)=DivBInternal_C(i,j,k)-&
             (FaceArea_D(1)*RightState_VY(Bx_,i,j,k)+&
             FaceArea_D(2)*RightState_VY(By_,i,j,k)+&
             FaceArea_D(3)*RightState_VY(Bz_,i,j,k))
        Source_VC(rhoUx_,i,j,k) = Source_VC(rhoUx_,i,j,k)&
             -B0xFace_y_BLK(i,j,k,globalBLK)*B1nJump
        Source_VC(rhoUy_,i,j,k) = Source_VC(rhoUy_,i,j,k)&
             -B0yFace_y_BLK(i,j,k,globalBLK)*B1nJump
        Source_VC(rhoUz_,i,j,k) = Source_VC(rhoUz_,i,j,k)&
             -B0zFace_y_BLK(i,j,k,globalBLK)*B1nJump
        DivB1_GB(i,j,k,globalBLK)  = DivB1_GB(i,j,k,globalBLK)+B1nJump

        call calc_faceareaJ(i,j+1,k,globalBLK,FaceArea_D)
        B1nJump = VInvHalf*&
             (FaceArea_D(1)*(RightState_VY(Bx_,i,j+1,k)-LeftState_VY(Bx_,i,j+1,k))+&
             FaceArea_D(2)*(RightState_VY(By_,i,j+1,k)-LeftState_VY(By_,i,j+1,k))+&
             FaceArea_D(3)*(RightState_VY(Bz_,i,j+1,k)-LeftState_VY(Bz_,i,j+1,k)))
        DivBInternal_C(i,j,k)=DivBInternal_C(i,j,k)+&
             (FaceArea_D(1)*LeftState_VY(Bx_,i,j+1,k)+&
             FaceArea_D(2)*LeftState_VY(By_,i,j+1,k)+&
             FaceArea_D(3)*LeftState_VY(Bz_,i,j+1,k))
        Source_VC(rhoUx_,i,j,k) = Source_VC(rhoUx_,i,j,k)&
             -B0xFace_y_BLK(i,j+1,k,globalBLK)*B1nJump
        Source_VC(rhoUy_,i,j,k) = Source_VC(rhoUy_,i,j,k)&
             -B0yFace_y_BLK(i,j+1,k,globalBLK)*B1nJump
        Source_VC(rhoUz_,i,j,k) = Source_VC(rhoUz_,i,j,k)&
             -B0zFace_y_BLK(i,j+1,k,globalBLK)*B1nJump
        DivB1_GB(i,j,k,globalBLK)  = DivB1_GB(i,j,k,globalBLK)+B1nJump
     end do; end do; end do

     iVolumeCounter=nIJK*(globalBLK-1)
     do k=1,nK; do j=1,nJ; do i=1,nI 
        iVolumeCounter=iVolumeCounter+1
        VInvHalf=VolumeInverse_I(iVolumeCounter)*cHalf
        call calc_faceareaK(i,j,k,globalBLK,FaceArea_D)
        B1nJump = VInvHalf*&
             (FaceArea_D(1)*(RightState_VZ(Bx_,i,j,k)-LeftState_VZ(Bx_,i,j,k))+&
             FaceArea_D(2)*(RightState_VZ(By_,i,j,k)-LeftState_VZ(By_,i,j,k))+&
             FaceArea_D(3)*(RightState_VZ(Bz_,i,j,k)-LeftState_VZ(Bz_,i,j,k)))
        DivBInternal_C(i,j,k)=DivBInternal_C(i,j,k)-&
             (FaceArea_D(1)*RightState_VZ(Bx_,i,j,k)+&
             FaceArea_D(2)*RightState_VZ(By_,i,j,k)+&
             FaceArea_D(3)*RightState_VZ(Bz_,i,j,k))

        Source_VC(rhoUx_,i,j,k) = Source_VC(rhoUx_,i,j,k)&
             -B0xFace_z_BLK(i,j,k,globalBLK)*B1nJump
        Source_VC(rhoUy_,i,j,k) = Source_VC(rhoUy_,i,j,k)&
             -B0yFace_z_BLK(i,j,k,globalBLK)*B1nJump
        Source_VC(rhoUz_,i,j,k) = Source_VC(rhoUz_,i,j,k)&
             -B0zFace_z_BLK(i,j,k,globalBLK)*B1nJump
        DivB1_GB(i,j,k,globalBLK)  = DivB1_GB(i,j,k,globalBLK)+B1nJump


        call calc_faceareaK(i,j,k+1,globalBLK,FaceArea_D)
        B1nJump = VInvHalf*&
             (FaceArea_D(1)*(RightState_VZ(Bx_,i,j,k+1)-LeftState_VZ(Bx_,i,j,k+1))+&
             FaceArea_D(2)*(RightState_VZ(By_,i,j,k+1)-LeftState_VZ(By_,i,j,k+1))+&
             FaceArea_D(3)*(RightState_VZ(Bz_,i,j,k+1)-LeftState_VZ(Bz_,i,j,k+1)))
        DivBInternal_C(i,j,k)=(DivBInternal_C(i,j,k)+&
             (FaceArea_D(1)*LeftState_VZ(Bx_,i,j,k+1)+&
             FaceArea_D(2)*LeftState_VZ(By_,i,j,k+1)+&
             FaceArea_D(3)*LeftState_VZ(Bz_,i,j,k+1)))*&
             VolumeInverse_I(iVolumeCounter)

        Source_VC(rhoUx_,i,j,k) = Source_VC(rhoUx_,i,j,k)&
             -B0xFace_z_BLK(i,j,k+1,globalBLK)*B1nJump
        Source_VC(rhoUy_,i,j,k) = Source_VC(rhoUy_,i,j,k)&
             -B0yFace_z_BLK(i,j,k+1,globalBLK)*B1nJump
        Source_VC(rhoUz_,i,j,k) = Source_VC(rhoUz_,i,j,k)&
             -B0zFace_z_BLK(i,j,k+1,globalBLK)*B1nJump
        DivB1_GB(i,j,k,globalBLK)  = DivB1_GB(i,j,k,globalBLK)+B1nJump
     end do; end do; end do
     do k=1,nK; do j=1,nJ; do i=1,nI 
        Source_VC(rhoUx_,i,j,k) = Source_VC(rhoUx_,i,j,k)-DivBInternal_C(i,j,k)*&
             B0xCell_BLK(i,j,k,globalBLK)
        Source_VC(rhoUy_,i,j,k) = Source_VC(rhoUy_,i,j,k)-DivBInternal_C(i,j,k)*&
             B0yCell_BLK(i,j,k,globalBLK)
        Source_VC(rhoUz_,i,j,k) = Source_VC(rhoUz_,i,j,k)-DivBInternal_C(i,j,k)*&
             B0zCell_BLK(i,j,k,globalBLK)               
        DivB1_GB(i,j,k,globalBLK)  = DivB1_GB(i,j,k,globalBLK)+&
             DivBInternal_C(i,j,k)
     end do; end do; end do

     ! Add contributions to other source terms
     do k=1,nK; do j=1,nJ; do i=1,nI
        Source_VC(rhoUx_:rhoUz_,i,j,k) =  Source_VC(rhoUx_:rhoUz_,i,j,k) -&
             DivB1_GB(i,j,k,globalBLK)* &
             State_VGB(Bx_:Bz_,i,j,k,globalBLK)
        RhoInv=cOne/State_VGB(rho_,i,j,k,globalBLK)
        Source_VC(Bx_:Bz_,i,j,k)    = Source_VC(Bx_:Bz_,i,j,k)-DivB1_GB(i,j,k,globalBLK)* &
             State_VGB(rhoUx_:rhoUz_,i,j,k,globalBLK)*RhoInv
        Source_VC(Energy_,i,j,k)     = Source_VC(Energy_,i,j,k)-DivB1_GB(i,j,k,globalBLK)* &
             sum(State_VGB(Bx_:Bz_,i,j,k,globalBLK)*&
             State_VGB(rhoUx_:rhoUz_,i,j,k,globalBLK))*&
             RhoInv
     end do;end do;end do

     if (UseB0Source) then
        do k=1,nK; do j=1,nJ; do i=1,nI
           Source_VC(rhoUx_:rhoUz_,i,j,k)=Source_VC(rhoUx_:rhoUz_,i,j,k) + &
                matmul( State_VGB(Bx_:Bz_,i,j,k,globalBLK),&
                B0SourceMatrix_DDCB(:,:,i,j,k,globalBLK)) 
        end do; end do; end do
     end if

  end if


  if(UseGravity.or.UseRotatingFrame) then
     Source_VC(rhoUx_,:,:,:) = Source_VC(rhoUx_,:,:,:) + &
          State_VGB(rho_,1:nI,1:nJ,1:nK,globalBLK)* &
          fbody_x_BLK(:,:,:,globalBLK)
     Source_VC(rhoUy_,:,:,:) = Source_VC(rhoUy_,:,:,:) + &
          State_VGB(rho_,1:nI,1:nJ,1:nK,globalBLK)* &
          fbody_y_BLK(:,:,:,globalBLK)
     Source_VC(rhoUz_,:,:,:) = Source_VC(rhoUz_,:,:,:) + &
          State_VGB(rho_,1:nI,1:nJ,1:nK,globalBLK)* &
          fbody_z_BLK(:,:,:,globalBLK)
     Source_VC(Energy_,:,:,:) = Source_VC(Energy_,:,:,:) + &
          (State_VGB(rhoUx_,1:nI,1:nJ,1:nK,globalBLK)* &
          fbody_x_BLK(:,:,:,globalBLK) + & 
          State_VGB(rhoUy_,1:nI,1:nJ,1:nK,globalBLK)* &
          fbody_y_BLK(:,:,:,globalBLK) + &
          State_VGB(rhoUz_,1:nI,1:nJ,1:nK,globalBLK)* &
          fbody_z_BLK(:,:,:,globalBLK)) 
  end if

  if(UseUserSource) call user_calc_sources

end subroutine calc_sources
!=============================End covariant_calc_sources.f90
