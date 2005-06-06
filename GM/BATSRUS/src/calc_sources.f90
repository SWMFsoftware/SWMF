!^CFG COPYRIGHT UM
!^CFG FILE CARTESIAN
subroutine calc_sources
  use ModProcMH
  use ModMain
  use ModVarIndexes
  use ModGeometry, ONLY : dx_BLK, dy_BLK, dz_BLK, R_BLK,&
       body_BLK, Rmin_BLK, VolumeInverse_I
  use ModGeometry, ONLY : R2_BLK         !^CFG IF SECONDBODY
  use ModAdvance
  use ModParallel, ONLY : NOBLK, neiLEV, &
       neiLtop,neiLbot,neiLeast,neiLwest,neiLnorth,neiLsouth
  use ModPhysics
  use ModNumConst
  use ModGeometry,ONLY :  fAx_BLK, fAy_BLK, fAz_BLK
  implicit none

  real:: B1nJump,DivBInternal

  integer :: i,j,k,iDim

  logical :: oktest,oktest_me

  real :: dr

  ! Variables needed for Boris source terms also used for div(u)

  real :: fullBx, fullBy, fullBz, Ux, Uy, Uz, coef, RhoInv
  real :: E_D(3), DivE

  ! Variables used by B0 source terms
  ! inverse of Dx, Dy, Dz, divergence B0
  real :: DxInv, DyInv, DzInv, DivB0
  ! non-linear difference terms, central B
  real :: dB0dBx, dB0dBy, dB0dBz, BxCell, ByCell, BzCell

  !---------------------------------------------------------------------------
  if(iProc==PROCtest .and. globalBLK==BLKtest)then
     call set_oktest('calc_sources',oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  end if

  Source_VC   = cZero

  ! Calculate source term for pressure: -(g-1)*P*Div(U)
  if(UseNonconservative)then
     do k=1,nK; do j=1,nJ; do i=1,nI
        Source_VC(P_,i,j,k) = -(g-1)*State_VGB(P_,i,j,k,globalBLK)*&
             VolumeInverse_I(globalBLK)*&
             (UDotFA_X(i+1,j,k)-UDotFA_X(i,j,k)+&
             UDotFA_Y(i,j+1,k) -UDotFA_Y(i,j,k)+&
             UDotFA_Z(i,j,k+1) -UDotFA_Z(i,j,k))
     end do; end do; end do
  end if
  ! No source terms should be calculated before DivbDiffusion (except SP)!

                        

  if(UseDivbDiffusion)then               !^CFG IF DIVBDIFFUSE BEGIN
     ! Apply Timur Linde's diffusive source terms
     call calc_divB 
     ! S[B]=coef_i*grad(div B)), where coef_i=divb_diffcoeff*(dx_i)^2/dt       
     Source_VC(Bx_,:,:,:) =(dx_BLK(globalBLK)*divb_diffcoeff*0.5)*&
          (DivB1_GB(2:nI+1,1:nJ  ,1:nK,globalBLK  )-&
          DivB1_GB(0:nI-1,1:nJ  ,1:nK,globalBLK  ))

     Source_VC(By_,:,:,:) =(dy_BLK(globalBLK)*divb_diffcoeff*0.5)*&
          (DivB1_GB(1:nI  ,2:nJ+1,1:nK,globalBLK  )-&
          DivB1_GB(1:nI  ,0:nJ-1,1:nK,globalBLK  ))

     Source_VC(Bz_,:,:,:) =(dz_BLK(globalBLK)*divb_diffcoeff*0.5)*&
          (DivB1_GB(1:nI  ,1:nJ  ,2:nK+1,globalBLK)-&
          DivB1_GB(1:nI  ,1:nJ  ,0:nK-1,globalBLK))


     ! Set diffusive source to zero at far field BC-s
     if(any(neiLEV(:,globalBLK)==NOBLK))then
        if(neiLeast(globalBLK)==NOBLK)then
           Source_VC(Bx_:Bz_,1,:,:)=cZero
        endif
        if(neiLwest(globalBLK)==NOBLK)then
           Source_VC(Bx_:Bz_,nI,:,:)=cZero
        endif
        if(neiLsouth(globalBLK)==NOBLK)then
           Source_VC(Bx_:Bz_,:,1,:)=cZero
        endif
        if(neiLnorth(globalBLK)==NOBLK)then
           Source_VC(Bx_:Bz_,:,nJ,:)=0.
        endif
        if(neiLbot(globalBLK)==NOBLK)then
           Source_VC(Bx_:Bz_,:,:,1)=cZero
        endif
        if(neiLtop(globalBLK)==NOBLK)then
           Source_VC(Bx_:Bz_,:,:,nK)=cZero
        endif
     end if

     if(body_BLK(globalBLK))then
        dr=sqrt(dx_BLK(globalBLK)**2+dy_BLK(globalBLK)**2+dz_BLK(globalBLK)**2)
        where(R_BLK(1:nI,1:nJ,1:nK,globalBLK)<Rbody+dr)
           Source_VC(Bx_,:,:,:)=0.0
           Source_VC(By_,:,:,:)=0.0
           Source_VC(Bz_,:,:,:)=0.0
        end where
        !^CFG IF SECONDBODY BEGIN
        if(UseBody2)then
           where(R2_BLK(1:nI,1:nJ,1:nK,globalBLK)<RBody2+dr)
              Source_VC(Bx_,:,:,:)=0.0
              Source_VC(By_,:,:,:)=0.0
              Source_VC(Bz_,:,:,:)=0.0
           end where
        end if
        !^CFG END SECONDBODY
     endif

     ! Divide by dt where possible

     where(time_BLK(1:nI,1:nJ,1:nK,globalBLK)>0.)
        Source_VC(Bx_,:,:,:)=Source_VC(Bx_,:,:,:)/cfl/time_BLK(1:nI,1:nJ,1:nK,globalBLK)
        Source_VC(By_,:,:,:)=Source_VC(By_,:,:,:)/cfl/time_BLK(1:nI,1:nJ,1:nK,globalBLK)
        Source_VC(Bz_,:,:,:)=Source_VC(Bz_,:,:,:)/cfl/time_BLK(1:nI,1:nJ,1:nK,globalBLK)
     endwhere

     ! S[E] = B.(coef_i.grad(div B)))

     Source_VC(Energy_,:,:,:) = &
          Source_VC(Bx_,:,:,:)*State_VGB(Bx_,1:nI,1:nJ,1:nK,globalBLK)+ &
          Source_VC(By_,:,:,:)*State_VGB(By_,1:nI,1:nJ,1:nK,globalBLK)+ &
          Source_VC(Bz_,:,:,:)*State_VGB(Bz_,1:nI,1:nJ,1:nK,globalBLK)
  end if                                     !^CFG END DIVBDIFFUSE
  if(UseDivbSource)then


     do k=1,nK; do j=1,nJ; do i=1,nI
        B1nJump = cHalf* VolumeInverse_I(globalBLK)*&
             fAx_BLK(globalBLK)*(RightState_VX(Bx_,i,j,k)-LeftState_VX(Bx_,i,j,k))
        Source_VC(rhoUx_,i,j,k) = -B0xFace_x_BLK(i,j,k,globalBLK)*B1nJump
        Source_VC(rhoUy_,i,j,k) = -B0yFace_x_BLK(i,j,k,globalBLK)*B1nJump
        Source_VC(rhoUz_,i,j,k) = -B0zFace_x_BLK(i,j,k,globalBLK)*B1nJump
        DivB1_GB(i,j,k,globalBLK)  = B1nJump
        B1nJump = cHalf* VolumeInverse_I(globalBLK)*&
             fAx_BLK(globalBLK)*(RightState_VX(Bx_,i+1,j,k)-LeftState_VX(Bx_,i+1,j,k))

        Source_VC(rhoUx_,i,j,k) = Source_VC(rhoUx_,i,j,k)&
             -B0xFace_x_BLK(i+1,j,k,globalBLK)*B1nJump
        Source_VC(rhoUy_,i,j,k) = Source_VC(rhoUy_,i,j,k)&
             -B0yFace_x_BLK(i+1,j,k,globalBLK)*B1nJump
        Source_VC(rhoUz_,i,j,k) = Source_VC(rhoUz_,i,j,k)&
             -B0zFace_x_BLK(i+1,j,k,globalBLK)*B1nJump
        DivB1_GB(i,j,k,globalBLK)  = DivB1_GB(i,j,k,globalBLK)+B1nJump
!     end do; end do; end do
!     do k=1,nK; do j=1,nJ; do i=1,nI 
        B1nJump = cHalf*VolumeInverse_I(globalBLK)*&
             fAy_BLK(globalBLK)*(RightState_VY(By_,i,j,k)-LeftState_VY(By_,i,j,k))

        Source_VC(rhoUx_,i,j,k) = Source_VC(rhoUx_,i,j,k)&
             -B0xFace_y_BLK(i,j,k,globalBLK)*B1nJump
        Source_VC(rhoUy_,i,j,k) = Source_VC(rhoUy_,i,j,k)&
             -B0yFace_y_BLK(i,j,k,globalBLK)*B1nJump
        Source_VC(rhoUz_,i,j,k) = Source_VC(rhoUz_,i,j,k)&
             -B0zFace_y_BLK(i,j,k,globalBLK)*B1nJump
        DivB1_GB(i,j,k,globalBLK)  = DivB1_GB(i,j,k,globalBLK)+B1nJump

        B1nJump = cHalf* VolumeInverse_I(globalBLK)*&
             fAy_BLK(globalBLK)*(RightState_VY(By_,i,j+1,k)-LeftState_VY(By_,i,j+1,k))

        Source_VC(rhoUx_,i,j,k) = Source_VC(rhoUx_,i,j,k)&
             -B0xFace_y_BLK(i,j+1,k,globalBLK)*B1nJump
        Source_VC(rhoUy_,i,j,k) = Source_VC(rhoUy_,i,j,k)&
             -B0yFace_y_BLK(i,j+1,k,globalBLK)*B1nJump
        Source_VC(rhoUz_,i,j,k) = Source_VC(rhoUz_,i,j,k)&
             -B0zFace_y_BLK(i,j+1,k,globalBLK)*B1nJump
        DivB1_GB(i,j,k,globalBLK)  = DivB1_GB(i,j,k,globalBLK)+B1nJump
!     end do; end do; end do
!     do k=1,nK; do j=1,nJ; do i=1,nI
        B1nJump = cHalf* VolumeInverse_I(globalBLK)*&
             fAz_BLK(globalBLK)*(RightState_VZ(Bz_,i,j,k)-LeftState_VZ(Bz_,i,j,k))

        Source_VC(rhoUx_,i,j,k) = Source_VC(rhoUx_,i,j,k)&
             -B0xFace_z_BLK(i,j,k,globalBLK)*B1nJump
        Source_VC(rhoUy_,i,j,k) = Source_VC(rhoUy_,i,j,k)&
             -B0yFace_z_BLK(i,j,k,globalBLK)*B1nJump
        Source_VC(rhoUz_,i,j,k) = Source_VC(rhoUz_,i,j,k)&
             -B0zFace_z_BLK(i,j,k,globalBLK)*B1nJump
        DivB1_GB(i,j,k,globalBLK)  = DivB1_GB(i,j,k,globalBLK)+B1nJump

        B1nJump = cHalf*VolumeInverse_I(globalBLK)*&
             fAz_BLK(globalBLK)*(RightState_VZ(Bz_,i,j,k+1)-LeftState_VZ(Bz_,i,j,k+1))

        Source_VC(rhoUx_,i,j,k) = Source_VC(rhoUx_,i,j,k)&
             -B0xFace_z_BLK(i,j,k+1,globalBLK)*B1nJump
        Source_VC(rhoUy_,i,j,k) = Source_VC(rhoUy_,i,j,k)&
             -B0yFace_z_BLK(i,j,k+1,globalBLK)*B1nJump
        Source_VC(rhoUz_,i,j,k) = Source_VC(rhoUz_,i,j,k)&
             -B0zFace_z_BLK(i,j,k+1,globalBLK)*B1nJump
        DivB1_GB(i,j,k,globalBLK)  = DivB1_GB(i,j,k,globalBLK)+B1nJump
!     end do; end do; end do
!     do k=1,nK; do j=1,nJ; do i=1,nI
        DivBInternal = VolumeInverse_I(globalBLK) *(&
             fAx_BLK(globalBLK)*(LeftState_VX(Bx_,i+1,j,k) -RightState_VX(Bx_,i,j,k))+&
             fAy_BLK(globalBLK)*(LeftState_VY(By_,i,j+1,k) -RightState_VY(By_,i,j,k))+&
             fAz_BLK(globalBLK)*(LeftState_VZ(Bz_,i,j,k+1) -RightState_VZ(Bz_,i,j,k)))
        Source_VC(rhoUx_,i,j,k) = Source_VC(rhoUx_,i,j,k)-DivBInternal*&
             B0xCell_BLK(i,j,k,globalBLK)
        Source_VC(rhoUy_,i,j,k) = Source_VC(rhoUy_,i,j,k)-DivBInternal*&
             B0yCell_BLK(i,j,k,globalBLK)
        Source_VC(rhoUz_,i,j,k) = Source_VC(rhoUz_,i,j,k)-DivBInternal*&
             B0zCell_BLK(i,j,k,globalBLK)               
        DivB1_GB(i,j,k,globalBLK)=DivB1_GB(i,j,k,globalBLK)+&
             DivBInternal
     end do; end do; end do

     if(oktest_me)write(*,*)'divb=',DivB1_GB(iTest,jTest,kTest,BlkTest)
     if(oktest_me.and.VarTest>=RhoUx_.and.VarTest<=RhoUz_)&
          call write_source('After B0B1 source')

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
        do k=1,nK; do j=1,nJ; do i=1,nI;do iDim=1,nDim
           Source_VC(rhoU_+iDim,i,j,k)=Source_VC(rhoU_+iDim,i,j,k) + &
                sum(State_VGB(Bx_:Bz_,i,j,k,globalBLK)*&
                B0SourceMatrix_DDCB(:,iDim,i,j,k,globalBLK)) 
        end do; end do; end do;end do
     end if
  else
     call calc_divB
  end if
  
  if(boris_correction .and. boris_cLIGHT_factor < 0.9999 & !^CFG IF BORISCORR BEGIN
       .and. index(test_string,'nodivE')<1) then
     
     coef= (boris_cLIGHT_factor**2 - 1.0)*inv_c2LIGHT
     do k=1,nK; do j=1,nJ; do i=1,nI
        fullBx = B0xCell_BLK(i,j,k,globalBLK)+State_VGB(Bx_,i,j,k,globalBLK)
        fullBy = B0yCell_BLK(i,j,k,globalBLK)+State_VGB(By_,i,j,k,globalBLK)
        fullBz = B0zCell_BLK(i,j,k,globalBLK)+State_VGB(Bz_,i,j,k,globalBLK)
        Ux     = State_VGB(rhoUx_,i,j,k,globalBLK)/State_VGB(rho_,i,j,k,globalBLK)
        Uy     = State_VGB(rhoUy_,i,j,k,globalBLK)/State_VGB(rho_,i,j,k,globalBLK)
        Uz     = State_VGB(rhoUz_,i,j,k,globalBLK)/State_VGB(rho_,i,j,k,globalBLK)
        E_D(x_) = fullBy*Uz - fullBz*Uy
        E_D(y_) = fullBz*Ux - fullBx*Uz
        E_D(z_) = fullBx*Uy - fullBy*Ux

     ! Calculate divergence of electric field 
        DivE = VolumeInverse_I(globalBLK)*&
             (EDotFA_X(i+1,j,k)-EDotFA_X(i,j,k)+&
             EDotFA_Y(i,j+1,k) -EDotFA_Y(i,j,k)+&
             EDotFA_Z(i,j,k+1) -EDotFA_Z(i,j,k))

        Source_VC(rhoUx_:rhoUz_,i,j,k) = Source_VC(rhoUx_:rhoUz_,i,j,k) + coef*DivE*E_D 
     end do; end do; end do
  end if                                                 !^CFG END BORISCORR
  
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
  
  if(UseUserSource) call user_calc_sources          !^CFG IF USERFILES
Contains
 
  subroutine write_source(String)
    character(len=*) :: String
    write(*,'(a,a)',advance='no')String," S=",Source_VC(VarTest,iTest,jTest,kTest) 
  end subroutine write_source
  !========================================================================= 
end subroutine calc_sources
subroutine calc_divb
  use ModMain, ONLY : &
       UseDivbDiffusion,&            !^CFG IF DIVBDIFFUSE
       nI,nJ,nK,globalBLK,test_string
  use ModVarIndexes, ONLY : Bx_,By_,Bz_
  use ModAdvance, ONLY : DivB1_GB,State_VGB, &
       LeftState_VX,RightState_VX,&
       LeftState_VY,RightState_VY,&
       LeftState_VZ,RightState_VZ
  use ModNumConst
  use ModGeometry, ONLY : dx_BLK,dy_BLK,dz_BLK,&
       fAx_BLK,fAy_BLK,fAz_BLK,VolumeInverse_I 
  implicit none

  !\
  ! Calculate div B for a block and store result into SdivB
  !/
  if(UseDivBDiffusion)then           !^CFG IF DIVBDIFFUSE BEGIN 
     ! Compute divB using central differencing because
     ! SdivB is needed in the first layer of ghost cells 
     !for diffusive div B control.

     DivB1_GB(0:nI+1,0:nJ+1,0:nK+1,globalBLK) = cHalf*(&
          (State_VGB(Bx_, 1:nI+2, 0:nJ+1, 0:nK+1,globalBLK)- &
          State_VGB(Bx_,-1:nI  , 0:nJ+1, 0:nK+1,globalBLK))/dx_BLK(globalBLK)+&
          (State_VGB(By_, 0:nI+1, 1:nJ+2, 0:nK+1,globalBLK)- &
          State_VGB(By_, 0:nI+1,-1:nJ  , 0:nK+1,globalBLK))/dy_BLK(globalBLK)+&
          (State_VGB(Bz_, 0:nI+1, 0:nJ+1, 1:nK+2,globalBLK)- &
          State_VGB(Bz_, 0:nI+1, 0:nJ+1,-1:nK  ,globalBLK))/dz_BLK(globalBLK))
  else                               !^CFG END DIVBDIFFUSE
     if(index(test_string,'DIVB_CD')>0)then
        ! Use central differencing if test string contains DIVB_CD
        DivB1_GB(1:nI,1:nJ,1:nK,globalBLK) = cHalf*(&
             (State_VGB(Bx_, 2:nI+1, 1:nJ  , 1:nK  ,globalBLK)- &
             State_VGB(Bx_, 0:nI-1, 1:nJ  , 1:nK  ,globalBLK))/dx_BLK(globalBLK)+&
             (State_VGB(By_, 1:nI  , 2:nJ+1, 1:nK  ,globalBLK)- &
             State_VGB(By_, 1:nI  , 0:nJ-1, 1:nK  ,globalBLK))/dy_BLK(globalBLK)+&
             (State_VGB(Bz_, 1:nI  , 1:nJ  , 2:nK+1,globalBLK)- &
             State_VGB(Bz_, 1:nI  , 1:nJ  , 0:nK-1,globalBLK))/dz_BLK(globalBLK))
     else
        ! Compute divB using averaged and conservatively corrected 
        ! left and right values

        DivB1_GB(1:nI,1:nJ,1:nK,globalBLK) = &
             cHalf * VolumeInverse_I(globalBLK) *( &
             fAx_BLK(globalBLK)*((LeftState_VX(Bx_,2:nI+1,1:nJ,1:nK)+    &
             RightState_VX(Bx_,2:nI+1,1:nJ,1:nK))-   &
             (LeftState_VX(Bx_,1:nI,1:nJ,1:nK)+    &
             RightState_VX(Bx_,1:nI,1:nJ,1:nK)))+  &
             fAy_BLK(globalBLK)*((LeftState_VY(By_,1:nI,2:nJ+1,1:nK)+    &
             RightState_VY(By_,1:nI,2:nJ+1,1:nK))-   &
             (LeftState_VY(By_,1:nI,1:nJ,1:nK)+    &
             RightState_VY(By_,1:nI,1:nJ,1:nK)))+  &
             fAz_BLK(globalBLK)*((LeftState_VZ(Bz_,1:nI,1:nJ,2:nK+1)+    &
             RightState_VZ(Bz_,1:nI,1:nJ,2:nK+1))-   &
             (LeftState_VZ(Bz_,1:nI,1:nJ,1:nK)+    &
             RightState_VZ(Bz_,1:nI,1:nJ,1:nK))))
     endif
  end if   !^CFG IF DIVBDIFFUSE
end subroutine calc_divb
