!^CFG COPYRIGHT UM
!^CFG FILE NOT COVARIANT
subroutine calc_sources
  use ModProcMH
  use ModMain
  use ModVarIndexes
  use ModGeometry, ONLY : dx_BLK, dy_BLK, dz_BLK, R_BLK,&
       body_BLK, Rmin_BLK, vInv_CB
  use ModGeometry, ONLY : R2_BLK                        !^CFG IF SECONDBODY
  use ModAdvance
  use ModParallel, ONLY : NOBLK, neiLEV, &
       neiLtop,neiLbot,neiLeast,neiLwest,neiLnorth,neiLsouth
  use ModPhysics
  use ModNumConst
  use ModResist,   ONLY : EtaResist_G                    !^CFG IF DISSFLUX
  use ModUser,     ONLY : user_calc_sources
  use ModHallResist,ONLY: UseHallResist, ResistDiag
  implicit none

  integer :: i, j, k, iDim

  logical :: DoTest, DoTestMe

  real :: Coef

  ! Variables needed for div B source terms
  real:: DxInvHalf, DyInvHalf, DzInvHalf, B1nJump, DivBInternal

  ! Variable for div B diffusion
  real :: Dr

  ! Variables needed for Boris source terms also used for div(u)
  real :: FullBx, FullBy, FullBz, Ux, Uy, Uz, RhoInv
  real :: E_D(3), DivE

  ! Variables needed for Joule heating
  real :: Current_D(3)

  !---------------------------------------------------------------------------
  if(iProc==PROCtest .and. globalBLK==BLKtest)then
     call set_oktest('calc_sources', DoTest, DoTestMe)
  else
     DoTest=.false.; DoTestMe=.false.
  end if

  Source_VC   = cZero

  ! Calculate source terms for pressure
  if(UseNonconservative)then
     ! Adiabatic heating: -(g-1)*P*Div(U)
     do k=1,nK; do j=1,nJ; do i=1,nI
        Source_VC(P_,i,j,k) = -(g-1)*State_VGB(P_,i,j,k,globalBLK)*&
             vInv_CB(i,j,k,globalBLK)*&
             (UDotFA_X(i+1,j,k)-UDotFA_X(i,j,k)+&
             UDotFA_Y(i,j+1,k) -UDotFA_Y(i,j,k)+&
             UDotFA_Z(i,j,k+1) -UDotFA_Z(i,j,k))
     end do; end do; end do

     ! Joule heating: dP/dt += (gamma-1)*eta*j**2
     if(UseResistFlux)then  !^CFG IF DISSFLUX BEGIN
        do k=1,nK; do j=1,nJ; do i=1,nI           
           call get_current(i,j,k,GlobalBlk,Current_D)
           Source_VC(P_,i,j,k) = Source_VC(P_,i,j,k) + &
                (g-1) * EtaResist_G(i,j,k) * sum(Current_D**2)
        end do; end do; end do
     end if                 !^CFG END DISSFLUX
     if(UseHallResist)then
        do k=1,nK; do j=1,nJ; do i=1,nI
           call get_current(i,j,k,GlobalBlk,Current_D)
           Source_VC(P_,i,j,k) = Source_VC(P_,i,j,k) + &
                (g-1) * ResistDiag * sum(Current_D**2)
        end do; end do; end do
     end if
  end if


  if(UseDivbSource)then

     DxInvHalf = 0.5/Dx_BLK(GlobalBlk)
     DyInvHalf = 0.5/Dy_BLK(GlobalBlk)
     DzInvHalf = 0.5/Dz_BLK(GlobalBlk)

     do k=1,nK; do j=1,nJ; do i=1,nI
        B1nJump = DxInvHalf*&
             (RightState_VX(Bx_,i,j,k)-LeftState_VX(Bx_,i,j,k))

        Source_VC(rhoUx_,i,j,k) = -B0xFace_x_BLK(i,j,k,globalBLK)*B1nJump
        Source_VC(rhoUy_,i,j,k) = -B0yFace_x_BLK(i,j,k,globalBLK)*B1nJump
        Source_VC(rhoUz_,i,j,k) = -B0zFace_x_BLK(i,j,k,globalBLK)*B1nJump
        DivB1_GB(i,j,k,globalBLK)  = B1nJump

        B1nJump = DxInvHalf*&
             (RightState_VX(Bx_,i+1,j,k)-LeftState_VX(Bx_,i+1,j,k))

        Source_VC(rhoUx_,i,j,k) = Source_VC(rhoUx_,i,j,k)&
             -B0xFace_x_BLK(i+1,j,k,globalBLK)*B1nJump
        Source_VC(rhoUy_,i,j,k) = Source_VC(rhoUy_,i,j,k)&
             -B0yFace_x_BLK(i+1,j,k,globalBLK)*B1nJump
        Source_VC(rhoUz_,i,j,k) = Source_VC(rhoUz_,i,j,k)&
             -B0zFace_x_BLK(i+1,j,k,globalBLK)*B1nJump
        DivB1_GB(i,j,k,globalBLK)  = DivB1_GB(i,j,k,globalBLK)+B1nJump

        B1nJump = DyInvHalf* &
             (RightState_VY(By_,i,j,k)-LeftState_VY(By_,i,j,k))

        Source_VC(rhoUx_,i,j,k) = Source_VC(rhoUx_,i,j,k)&
             -B0xFace_y_BLK(i,j,k,globalBLK)*B1nJump
        Source_VC(rhoUy_,i,j,k) = Source_VC(rhoUy_,i,j,k)&
             -B0yFace_y_BLK(i,j,k,globalBLK)*B1nJump
        Source_VC(rhoUz_,i,j,k) = Source_VC(rhoUz_,i,j,k)&
             -B0zFace_y_BLK(i,j,k,globalBLK)*B1nJump
        DivB1_GB(i,j,k,globalBLK)  = DivB1_GB(i,j,k,globalBLK)+B1nJump

        B1nJump = DyInvHalf* &
             (RightState_VY(By_,i,j+1,k)-LeftState_VY(By_,i,j+1,k))

        Source_VC(rhoUx_,i,j,k) = Source_VC(rhoUx_,i,j,k)&
             -B0xFace_y_BLK(i,j+1,k,globalBLK)*B1nJump
        Source_VC(rhoUy_,i,j,k) = Source_VC(rhoUy_,i,j,k)&
             -B0yFace_y_BLK(i,j+1,k,globalBLK)*B1nJump
        Source_VC(rhoUz_,i,j,k) = Source_VC(rhoUz_,i,j,k)&
             -B0zFace_y_BLK(i,j+1,k,globalBLK)*B1nJump
        DivB1_GB(i,j,k,globalBLK)  = DivB1_GB(i,j,k,globalBLK)+B1nJump

        B1nJump = DzInvHalf * &
             (RightState_VZ(Bz_,i,j,k)-LeftState_VZ(Bz_,i,j,k))

        Source_VC(rhoUx_,i,j,k) = Source_VC(rhoUx_,i,j,k)&
             -B0xFace_z_BLK(i,j,k,globalBLK)*B1nJump
        Source_VC(rhoUy_,i,j,k) = Source_VC(rhoUy_,i,j,k)&
             -B0yFace_z_BLK(i,j,k,globalBLK)*B1nJump
        Source_VC(rhoUz_,i,j,k) = Source_VC(rhoUz_,i,j,k)&
             -B0zFace_z_BLK(i,j,k,globalBLK)*B1nJump
        DivB1_GB(i,j,k,globalBLK)  = DivB1_GB(i,j,k,globalBLK)+B1nJump

        B1nJump = DzInvHalf * &
             (RightState_VZ(Bz_,i,j,k+1)-LeftState_VZ(Bz_,i,j,k+1))

        Source_VC(rhoUx_,i,j,k) = Source_VC(rhoUx_,i,j,k)&
             -B0xFace_z_BLK(i,j,k+1,globalBLK)*B1nJump
        Source_VC(rhoUy_,i,j,k) = Source_VC(rhoUy_,i,j,k)&
             -B0yFace_z_BLK(i,j,k+1,globalBLK)*B1nJump
        Source_VC(rhoUz_,i,j,k) = Source_VC(rhoUz_,i,j,k)&
             -B0zFace_z_BLK(i,j,k+1,globalBLK)*B1nJump
        DivB1_GB(i,j,k,globalBLK)  = DivB1_GB(i,j,k,globalBLK)+B1nJump

        DivBInternal = 2*(&
             DxInvHalf*(LeftState_VX(Bx_,i+1,j,k) -RightState_VX(Bx_,i,j,k))+&
             DyInvHalf*(LeftState_VY(By_,i,j+1,k) -RightState_VY(By_,i,j,k))+&
             DzInvHalf*(LeftState_VZ(Bz_,i,j,k+1) -RightState_VZ(Bz_,i,j,k)))
        Source_VC(rhoUx_,i,j,k) = Source_VC(rhoUx_,i,j,k)-DivBInternal*&
             B0xCell_BLK(i,j,k,globalBLK)
        Source_VC(rhoUy_,i,j,k) = Source_VC(rhoUy_,i,j,k)-DivBInternal*&
             B0yCell_BLK(i,j,k,globalBLK)
        Source_VC(rhoUz_,i,j,k) = Source_VC(rhoUz_,i,j,k)-DivBInternal*&
             B0zCell_BLK(i,j,k,globalBLK)               
        DivB1_GB(i,j,k,globalBLK)=DivB1_GB(i,j,k,globalBLK)+&
             DivBInternal
     end do; end do; end do

     if(DoTestMe)write(*,*)'divb=',DivB1_GB(iTest,jTest,kTest,BlkTest)
     if(DoTestMe.and.VarTest>=RhoUx_.and.VarTest<=RhoUz_)&
          call write_source('After B0B1 source')

     ! Add contributions to other source terms
     do k=1,nK; do j=1,nJ; do i=1,nI
        Source_VC(rhoUx_:rhoUz_,i,j,k) =  Source_VC(rhoUx_:rhoUz_,i,j,k) -&
             DivB1_GB(i,j,k,globalBLK)* &
             State_VGB(Bx_:Bz_,i,j,k,globalBLK)
        RhoInv=cOne/State_VGB(rho_,i,j,k,globalBLK)
        Source_VC(Bx_:Bz_,i,j,k)    = Source_VC(Bx_:Bz_,i,j,k) &
             - DivB1_GB(i,j,k,globalBLK)* &
             State_VGB(rhoUx_:rhoUz_,i,j,k,globalBLK)*RhoInv
        Source_VC(Energy_,i,j,k)     = Source_VC(Energy_,i,j,k) &
             -DivB1_GB(i,j,k,globalBLK)* &
             sum(State_VGB(Bx_:Bz_,i,j,k,globalBLK)*&
             State_VGB(rhoUx_:rhoUz_,i,j,k,globalBLK))*RhoInv
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
        FullBx = B0xCell_BLK(i,j,k,globalBLK)+State_VGB(Bx_,i,j,k,globalBLK)
        FullBy = B0yCell_BLK(i,j,k,globalBLK)+State_VGB(By_,i,j,k,globalBLK)
        FullBz = B0zCell_BLK(i,j,k,globalBLK)+State_VGB(Bz_,i,j,k,globalBLK)
        Ux = State_VGB(rhoUx_,i,j,k,globalBLK)/State_VGB(rho_,i,j,k,globalBLK)
        Uy = State_VGB(rhoUy_,i,j,k,globalBLK)/State_VGB(rho_,i,j,k,globalBLK)
        Uz = State_VGB(rhoUz_,i,j,k,globalBLK)/State_VGB(rho_,i,j,k,globalBLK)
        E_D(x_) = FullBy*Uz - FullBz*Uy
        E_D(y_) = FullBz*Ux - FullBx*Uz
        E_D(z_) = FullBx*Uy - FullBy*Ux

        ! Calculate divergence of electric field 
        DivE = vInv_CB(i,j,k,globalBLK)*&
             (EDotFA_X(i+1,j,k)-EDotFA_X(i,j,k)+&
             EDotFA_Y(i,j+1,k) -EDotFA_Y(i,j,k)+&
             EDotFA_Z(i,j,k+1) -EDotFA_Z(i,j,k))

        Source_VC(rhoUx_:rhoUz_,i,j,k) = Source_VC(rhoUx_:rhoUz_,i,j,k) &
             + Coef*DivE*E_D 
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

  ! Add Coriolis forces here
  if(UseRotatingFrame)then
     select case(TypeCoordSystem)
     case('HGC','HGR')
        ! This is a special case since Omega is parallel with the Z axis
        do k=1,nK; do j=1,nJ; do i=1,nI
           Source_VC(rhoUx_,i,j,k) = Source_VC(rhoUx_,i,j,k) + &
                cTwo*OmegaBody*State_VGB(rhoUy_,i,j,k,globalBLK)
           Source_VC(rhoUy_,i,j,k) = Source_VC(rhoUy_,i,j,k) - &
                cTwo*OmegaBody*State_VGB(rhoUx_,i,j,k,globalBLK)
        end do; end do; end do
     case default
        call stop_mpi('ERROR in calc_sources: '// &
             'Coriolis force is not implemented for '// &
             'TypeCoordSystem=',TypeCoordSystem)
     end select
  end if

  if(UseUserSource) call user_calc_sources

contains
  !===========================================================================
  subroutine write_source(String)
    character(len=*) :: String
    write(*,'(a,a)',advance='no')String," S=",Source_VC(VarTest,iTest,jTest,kTest) 
  end subroutine write_source

end subroutine calc_sources

!=============================================================================

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
       fAx_BLK,fAy_BLK,fAz_BLK,vInv_CB 
  implicit none

  !\
  ! Calculate div B for a block and store result into SdivB
  !/
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
          cHalf * vInv_CB(:,:,:,globalBLK) *( &
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
end subroutine calc_divb

!==============================================================================

subroutine get_current(i,j,k,iBlock,Current_D)

  ! Calculate the current in a cell of a block

  use ModAdvance,  ONLY: State_VGB, Bx_, By_, Bz_
  use ModGeometry, ONLY: True_Cell, Dx_BLK, Dy_BLK, Dz_BLK

  implicit none
  integer, intent(in) :: i,j,k,iBlock
  real,    intent(out):: Current_D(3)

  real :: DxInvHalf, DyInvHalf, DzInvHalf
  !----------------------------------------------------------------------------

  ! Exclude cells next to the body because they produce incorrect currents
  if(.not.all(True_Cell(i-1:i+1,j-1:j+1,k-1:k+1,iBlock)))then
     Current_D = 0.0
     RETURN
  endif

  DxInvHalf = 0.5/Dx_BLK(iBlock)
  DyInvHalf = 0.5/Dy_BLK(iBlock)
  DzInvHalf = 0.5/Dz_BLK(iBlock)

  Current_D(1) = &
       (State_VGB(Bz_,i,j+1,k,iBlock) &
       -State_VGB(Bz_,i,j-1,k,iBlock))*DyInvHalf - &
       (State_VGB(By_,i,j,k+1,iBlock) &
       -State_VGB(By_,i,j,k-1,iBlock))*DzInvHalf

  Current_D(2) = &
       (State_VGB(Bx_,i,j,k+1,iBlock) &
       -State_VGB(Bx_,i,j,k-1,iBlock))*DzInvHalf- &
       (State_VGB(Bz_,i+1,j,k,iBlock) &
       -State_VGB(Bz_,i-1,j,k,iBlock))*DxInvHalf

  Current_D(3) = &
       (State_VGB(By_,i+1,j,k,iBlock) &
       -State_VGB(By_,i-1,j,k,iBlock))*DxInvHalf- &
       (State_VGB(Bx_,i,j+1,k,iBlock) &
       -State_VGB(Bx_,i,j-1,k,iBlock))*DyInvHalf

end subroutine get_current
