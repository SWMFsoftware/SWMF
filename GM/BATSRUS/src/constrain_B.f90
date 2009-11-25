!^CFG COPYRIGHT UM
!^CFG FILE CONSTRAINB

! A flux averaged constrained transport scheme for block AMR grid. See 
!
! G. Toth, 2000, Journal of Computational Physics, 161, 605-652
!
! G. Toth and P. L. Roe, 2002, Journal of Computational Phys, 180, 736-750

subroutine get_VxB(iBlock)

  ! Calculate VxB from fluxes following Balsara and Spicer

  use ModMain, ONLY : nI,nJ,nK,iTest,jTest,kTest,BLKtest
  use ModVarIndexes, ONLY : Bx_,By_,Bz_
  use ModAdvance, ONLY : Flux_VX,Flux_VY,Flux_VZ
  use ModGeometry, ONLY : fAx_BLK,fAy_BLK,fAz_BLK
  use ModCT, ONLY : VxB_x,VxB_y,VxB_z
  implicit none

  integer, intent(in) :: iBlock

  logical :: oktest, oktest_me
  !-------------------------------------------------------------------------
  if(iBlock==BLKtest)then
     call set_oktest('get_vxb',oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  end if

  ! VxB_x=(fy+fy-fz-fz)/4
  VxB_x(1:nI,1:nJ+1,1:nK+1,iBlock)= 0.25*(                &
          (Flux_VY(Bz_,1:nI,1:nJ+1,0:nK  )                     &
          +Flux_VY(Bz_,1:nI,1:nJ+1,1:nK+1))/fAy_BLK(iBlock) &
         -(Flux_VZ(By_,1:nI,0:nJ  ,1:nK+1)                     &
          +Flux_VZ(By_,1:nI,1:nJ+1,1:nK+1))/fAz_BLK(iBlock))

  ! VxB_y=(fz+fz-fx-fx)/4
  VxB_y(1:nI+1,1:nJ,1:nK+1,iBlock)= 0.25*(                &
          (Flux_VZ(Bx_,0:nI  ,1:nJ,1:nK+1)                     &
          +Flux_VZ(Bx_,1:nI+1,1:nJ,1:nK+1))/fAz_BLK(iBlock) &
         -(Flux_VX(Bz_,1:nI+1,1:nJ,0:nK  )                     &
          +Flux_VX(Bz_,1:nI+1,1:nJ,1:nK+1))/fAx_BLK(iBlock)) 
     
  ! VxB_z=(fx+fx-fy-fy)/4
  VxB_z(1:nI+1,1:nJ+1,1:nK,iBlock)= 0.25*(                &
          (Flux_VX(By_,1:nI+1,0:nJ  ,1:nK)                     &
          +Flux_VX(By_,1:nI+1,1:nJ+1,1:nK))/fAx_BLK(iBlock) &
         -(Flux_VY(Bx_,0:nI  ,1:nJ+1,1:nK)                     &
          +Flux_VY(Bx_,1:nI+1,1:nJ+1,1:nK))/fAy_BLK(iBlock))

  if(oktest_me)then
     write(*,*)'get_vxb: final VxB (edge centered)'
     write(*,*)'VxB_xLL,LR,RL,RR=',&
          VxB_x(iTest,jTest:jTest+1,kTest:kTest+1,BlkTest)
     write(*,*)'VxB_yLL,LR,RL,RR=',&
          VxB_y(iTest:iTest+1,jTest,kTest:kTest+1,BlkTest)
     write(*,*)'VxB_zLL,LR,RL,RR=',&
          VxB_z(iTest:iTest+1,jTest:jTest+1,kTest,BlkTest)
  end if

end subroutine get_VxB

!=============================================================================

subroutine bound_VxB(iBlock)

  ! Apply boundary conditions on VxB 

  use ModSize
  use ModMain, ONLY : TypeBc_I
  use ModVarIndexes, ONLY : Bx_,By_,Bz_
  use ModAdvance, ONLY : Flux_VX,Flux_VY,Flux_VZ
  use ModParallel, ONLY : NOBLK,&
       neiLtop,neiLbot,neiLeast,neiLwest,neiLnorth,neiLsouth
  use ModGeometry, ONLY : fAx_BLK,fAy_BLK,fAz_BLK,true_cell,body_BLK
  use ModCT, ONLY : VxB_x,VxB_y,VxB_z
  use ModPhysics, ONLY: SW_UX,SW_UY,SW_UZ,SW_BX,SW_BY,SW_BZ
  implicit none

  integer, intent(in) :: iBlock

  integer, parameter :: VxB_BC_order=1

  integer:: i,j,k
  !-------------------------------------------------------------------------

  ! Apply continuous or fixed boundary conditions at outer boundaries
  if(neiLeast(iBlock)==NOBLK)then
     do k=1,nK+1; do j=1,nJ
        VxB_y(1,j,k,iBlock) = +Flux_VZ(Bx_,1,j,k)/fAz_BLK(iBlock)
     end do; end do
     do k=1,nK; do j=1,nJ+1
        VxB_z(1,j,k,iBlock) = -Flux_VY(Bx_,1,j,k)/fAy_BLK(iBlock)
     end do; end do
  end if
  if(neiLwest(iBlock)==NOBLK)then
     ! fixed inflow!
     !VxB_x(nI  ,:,:,iBlock)=SW_Uy*SW_Bz-SW_Uz*SW_Uy
     select case(TypeBc_I(west_))
     case('inflow','vary','fixed')
        VxB_y(nI+1,:,:,iBlock)=SW_Uz*SW_Bx-SW_Ux*SW_Bz
        VxB_z(nI+1,:,:,iBlock)=SW_Ux*SW_By-SW_Uy*SW_Bx
     case default
        ! continuous
        do k=1,nK+1; do j=1,nJ
           VxB_y(nI+1,j,k,iBlock) = +Flux_VZ(Bx_,nI,j,k)/fAz_BLK(iBlock)
        end do; end do
        do k=1,nK; do j=1,nJ+1
           VxB_z(nI+1,j,k,iBlock) = -Flux_VY(Bx_,nI,j,k)/fAy_BLK(iBlock)
        end do; end do
     end select
  end if
  if(neiLsouth(iBlock)==NOBLK)then
     do k=1,nK+1; do i=1,nI
        VxB_x(i,1,k,iBlock) = -Flux_VZ(By_,i,1,k)/fAz_BLK(iBlock)
     end do; end do
     do k=1,nK; do i=1,nI+1
        VxB_z(i,1,k,iBlock) = +Flux_VX(By_,i,1,k)/fAx_BLK(iBlock)
     end do; end do
  end if
  if(neiLnorth(iBlock)==NOBLK)then
     do k=1,nK+1; do i=1,nI
        VxB_x(i,nJ+1,k,iBlock) = -Flux_VZ(By_,i,nJ,k)/fAz_BLK(iBlock)
     end do; end do
     do k=1,nK; do i=1,nI+1
        VxB_z(i,nJ+1,k,iBlock) = +Flux_VX(By_,i,nJ,k)/fAx_BLK(iBlock)
     end do; end do
  end if
  if(neiLbot(iBlock)==NOBLK)then
     do j=1,nJ+1; do i=1,nI
        VxB_x(i,j,1,iBlock) = +Flux_VY(Bz_,i,j,1)/fAy_BLK(iBlock)
     end do; end do
     do j=1,nJ; do i=1,nI+1
        VxB_y(i,j,1,iBlock) = -Flux_VX(Bz_,i,j,1)/fAx_BLK(iBlock)
     end do; end do
  end if
  if(neiLtop(iBlock)==NOBLK)then
     do j=1,nJ+1; do i=1,nI
        VxB_x(i,j,nK+1,iBlock) = +Flux_VY(Bz_,i,j,nK)/fAy_BLK(iBlock)
     end do; end do
     do j=1,nJ; do i=1,nI+1
        VxB_y(i,j,nK+1,iBlock) = -Flux_VX(Bz_,i,j,nK)/fAx_BLK(iBlock)
     end do; end do
  end if

  !!! Set VxB to zero on the cell edges of the body cells !!!
  if(body_BLK(iBlock))then
     ! Apply inner boundary condition on the electric field
     ! Make sure that edges belonging to body ghost cells are also corrected
     do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
        if(.not.true_cell(i,j,k,iBlock))then
           VxB_x(i,j:j+1,k:k+1,iBlock)=0.0
           VxB_y(i:i+1,j,k:k+1,iBlock)=0.0
           VxB_z(i:i+1,j:j+1,k,iBlock)=0.0
        end if
     end do; end do; end do
  end if

end subroutine bound_VxB

!=============================================================================

subroutine constrain_B(iBlock)

  ! Use CT scheme for updating the B field so that div B is conserved

  use ModSize
  use ModMain, ONLY : Dt,BLKtest,iTest,jTest,kTest
  use ModGeometry, ONLY : dx_BLK,dy_BLK,dz_BLK
  use ModCT, ONLY : VxB_x,VxB_y,VxB_z,Bxface_BLK,Byface_BLK,Bzface_BLK
  implicit none

  integer, intent(in) :: iBlock

  real :: qdt
  logical :: oktest,oktest_me
  !-------------------------------------------------------------------------
  if(iBlock==BLKtest)then
     call set_oktest('constrain_b',oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  end if

  ! Calculate physical time step
  qdt=dt

  if(oktest_me)then
     write(*,*)'constrain_b: initial face centered B'
     write(*,*)'BxfaceL,R=',&
          BxFace_BLK(iTest:iTest+1,jTest,kTest,BlkTest)
     write(*,*)'ByfaceL,R=',&
          ByFace_BLK(iTest,jTest:jTest+1,kTest,BlkTest)
     write(*,*)'BzfaceL,BzfaceR=',&
          BzFace_BLK(iTest,jTest,kTest:kTest+1,BlkTest)
  end if

  ! dBx/dt=d(VxB_z)/dy-d(VxB_y)/dz
  Bxface_BLK(1:nI+1,1:nJ,1:nK,iBlock)=                &
       Bxface_BLK(1:nI+1,1:nJ,1:nK,iBlock) + qdt*(    &
       +(VxB_z(1:nI+1,2:nJ+1,1:nK  ,iBlock)           &
        -VxB_z(1:nI+1,1:nJ  ,1:nK  ,iBlock))          &
                                                    /dy_BLK(iBlock)  &
       -(VxB_y(1:nI+1,1:nJ  ,2:nK+1,iBlock)           &
        -VxB_y(1:nI+1,1:nJ  ,1:nK  ,iBlock))          &
                                                    /dz_BLK(iBlock))
  ! dBy/dt=d(VxB_x)/dz-d(VxB_z)/dx
  Byface_BLK(1:nI,1:nJ+1,1:nK,iBlock)=                &
       Byface_BLK(1:nI,1:nJ+1,1:nK,iBlock) + qdt*(    &
       +(VxB_x(1:nI  ,1:nJ+1,2:nK+1,iBlock)           &
        -VxB_x(1:nI  ,1:nJ+1,1:nK  ,iBlock))          &
                                                    /dz_BLK(iBlock)  &
       -(VxB_z(2:nI+1,1:nJ+1,1:nK  ,iBlock)           &
        -VxB_z(1:nI  ,1:nJ+1,1:nK  ,iBlock))          &
                                                    /dx_BLK(iBlock))

  ! dBz/dt=d(VxB_y)/dx-d(VxB_x)/dy
  Bzface_BLK(1:nI,1:nJ,1:nK+1,iBlock)=                &
       Bzface_BLK(1:nI,1:nJ,1:nK+1,iBlock) + qdt*(    &
       +(VxB_y(2:nI+1,1:nJ  ,1:nK+1,iBlock)           &
        -VxB_y(1:nI  ,1:nJ  ,1:nK+1,iBlock))          &
                                                    /dx_BLK(iBlock)  &
       -(VxB_x(1:nI  ,2:nJ+1,1:nK+1,iBlock)           &
        -VxB_x(1:nI  ,1:nJ  ,1:nK+1,iBlock))          &
                                                    /dy_BLK(iBlock))
  if(oktest_me)then
     write(*,*)'constrain_b: final face centered B'
     write(*,*)'BxfaceL,R=',&
          BxFace_BLK(iTest:iTest+1,jTest,kTest,BlkTest)
     write(*,*)'ByfaceL,R=',&
          ByFace_BLK(iTest,jTest:jTest+1,kTest,BlkTest)
     write(*,*)'BzfaceL,BzfaceR=',&
          BzFace_BLK(iTest,jTest,kTest:kTest+1,BlkTest)
  end if

end subroutine constrain_B

!==============================================================================

subroutine Bface2Bcenter(iBlock)

  use ModSize
  use ModVarIndexes, ONLY : Bx_,By_,Bz_
  use ModAdvance, ONLY : State_VGB
  use ModGeometry, ONLY : true_cell,body_BLK
  use ModCT, ONLY : Bxface_BLK,Byface_BLK,Bzface_BLK
  implicit none

  integer, intent(in) :: iBlock

  !---------------------------------------------------------------------------

  ! average in direction x (b->B)
  State_VGB(Bx_,1:nI,1:nJ,1:nK,iBlock)= 0.5*(      &
       Bxface_BLK(1:nI  ,1:nJ,1:nK,iBlock)+ &
       Bxface_BLK(2:nI+1,1:nJ,1:nK,iBlock))

  ! average in direction y (b->B)
  State_VGB(By_,1:nI,1:nJ,1:nK,iBlock)= 0.5*(      &
       Byface_BLK(1:nI,1:nJ  ,1:nK,iBlock)+ &
       Byface_BLK(1:nI,2:nJ+1,1:nK,iBlock))

  ! average in direction z (b->B)
  State_VGB(Bz_,1:nI,1:nJ,1:nK,iBlock)= 0.5*(      &
       Bzface_BLK(1:nI,1:nJ,1:nK  ,iBlock)+ &
       Bzface_BLK(1:nI,1:nJ,2:nK+1,iBlock))

  if(body_BLK(iBlock))then
     where(.not.true_cell(:,:,:,iBlock))
        State_VGB(Bx_,:,:,:,iBlock)=0.0
        State_VGB(By_,:,:,:,iBlock)=0.0
        State_VGB(Bz_,:,:,:,iBlock)=0.0
     end where
  end if

end subroutine Bface2Bcenter

!==============================================================================

subroutine Bcenter2Bface(iBlock)

  use ModSize
  use ModVarIndexes, ONLY : Bx_,By_,Bz_
  use ModAdvance, ONLY : State_VGB
  use ModCT, ONLY : Bxface_BLK,Byface_BLK,Bzface_BLK
  use ModMain, ONLY: UseConstrainB
  implicit none

  integer, intent(in) :: iBlock

  integer:: i,j,k
  !---------------------------------------------------------------------------

  ! Estimate BFace from Bcenter

  ! The conditional write statements avoid a compiler optimization bug
  ! in the NAGWare Fortran 95 compiler Release 4.0a(388).
  ! The condition is never true, because this routine is only called
  ! when UseConstrainB is true

  do k=1,nK; do j=1,nJ; do i=1,nI+1
     BxFace_BLK(i,j,k,iBlock)= 0.5*( &
          State_VGB(Bx_,i-1,j,k,iBlock)+ &
          State_VGB(Bx_,i  ,j,k,iBlock))
     if(.not.UseConstrainB)write(*,*)'!!!'
  end do; end do; end do
  do k=1,nK; do j=1,nJ+1; do i=1,nI
     ByFace_BLK(i,j,k,iBlock)= 0.5*( &
          State_VGB(By_,i,j-1,k,iBlock)+ &
          State_VGB(By_,i,j  ,k,iBlock))
     if(.not.UseConstrainB)write(*,*)'!!!'
  end do; end do; end do
  do k=1,nK+1; do j=1,nJ; do i=1,nI
     BzFace_BLK(i,j,k,iBlock)= 0.5*( &
          State_VGB(Bz_,i,j,k-1,iBlock)+ &
          State_VGB(Bz_,i,j,k  ,iBlock))
     if(.not.UseConstrainB)write(*,*)'!!!'
  end do; end do; end do

  call bound_Bface(iBlock)

end subroutine Bcenter2Bface

!==============================================================================
subroutine bound_Bface(iBlock)

  !!! Set Bface to zero on the cell faces of the body cells !!!
  ! Make sure that ghost cells inside the body are taken into account
  ! This may have to be generalized later

  use ModSize
  use ModMain,     ONLY: BLKtest
  use ModGeometry, ONLY: true_cell,body_BLK
  use ModCT,       ONLY: Bxface_BLK,Byface_BLK,Bzface_BLK
  implicit none

  integer, intent(in) :: iBlock

  integer :: i,j,k

  logical :: oktest, oktest_me
  !---------------------------------------------------------------------------

  if(iBlock==BLKtest)then
     call set_oktest('bound_Bface',oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  end if

  if(oktest_me)write(*,*)'bound_Bface, body_BLK=',body_BLK(iBlock)

  if(body_BLK(iBlock))then
     do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
        if(.not.true_cell(i,j,k,iBlock))then
           BxFace_BLK(i:i+1,j,k,iBlock)=0.0
           ByFace_BLK(i,j:j+1,k,iBlock)=0.0
           BzFace_BLK(i,j,k:k+1,iBlock)=0.0
        end if
     end do; end do; end do
  end if

end subroutine bound_Bface

!=============================================================================
subroutine prolong1_Bface(coarse_sol, iVar, iBLK, fine_sol)

  ! First order div B conserving prolongation for Bface

  use ModSize
  use ModVarIndexes, ONLY : Bx_,By_,Bz_
  implicit none

  integer, intent(in) :: iVar, iBLK
  real, dimension (1-gcn:nI+gcn,1-gcn:nJ+gcn, 1-gcn:nK+gcn) :: &
       coarse_sol,fine_sol

  integer :: i,j,k,i1,j1,k1,ishift,jshift,kshift
  !--------------------------------------------------------------------------

  call get_shifts(iBLK,ishift,jshift,kshift)

  ! Assign default solution state to fine block to get corners
  fine_sol=0.0

  !\
  ! Prolong coarse grid solution to finer block.
  !/

  select case(iVar)
  case(Bx_)
     !BxFace
     do i = 1+ishift, nI/2+ishift+1
        do j = 1+jshift, nJ/2+jshift
           do k = 1+kshift, nK/2+kshift
              i1 = 2*(i-ishift)-1
              j1 = 2*(j-jshift)-1
              k1 = 2*(k-kshift)-1

              fine_sol(i1,j1:j1+1,k1:k1+1) = coarse_sol(i,j,k)

              if(i1<nI+1)fine_sol(i1+1,j1:j1+1,k1:k1+1) = &
                   0.5*(coarse_sol(i,j,k)+coarse_sol(i+1,j,k))
           end do
        end do
     end do
  case(By_)
     !ByFace
     do i = 1+ishift, nI/2+ishift
        do j = 1+jshift, nJ/2+jshift+1
           do k = 1+kshift, nK/2+kshift
              i1 = 2*(i-ishift)-1
              j1 = 2*(j-jshift)-1
              k1 = 2*(k-kshift)-1

              fine_sol(i1:i1+1,j1,k1:k1+1) = coarse_sol(i,j,k)

              if(j1<nJ+1)fine_sol(i1:i1+1,j1+1,k1:k1+1) = &
                   0.5*(coarse_sol(i,j,k)+coarse_sol(i,j+1,k))
           end do
        end do
     end do
  case(Bz_)
     !BzFace
     do i = 1+ishift, nI/2+ishift
        do j = 1+jshift, nJ/2+jshift
           do k = 1+kshift, nK/2+kshift+1
              i1 = 2*(i-ishift)-1
              j1 = 2*(j-jshift)-1
              k1 = 2*(k-kshift)-1

              fine_sol(i1:i1+1,j1:j1+1,k1) = coarse_sol(i,j,k)

              if(k1<nK+1)fine_sol(i1:i1+1,j1:j1+1,k1+1) = &
                   0.5*(coarse_sol(i,j,k)+coarse_sol(i,j,k+1))
           end do
        end do
     end do
  end select

end subroutine prolong1_Bface

!=============================================================================

subroutine prolong_b_face(Bxf_c,Byf_c,Bzf_c,&
     BxFaceFine_XQS,ByFaceFine_YQS,BzFaceFine_ZQS,&
     IsFinerNei_E,iChild,iBlock,Bxf_f,Byf_f,Bzf_f)

  ! Second order div B conserving prolongation for Bface
  ! _c is coarse, _f is fine

  use ModSize
  use ModMain, ONLY : BLKtest, iTest, jTest, kTest
  use ModGeometry, ONLY: dx_BLK, dy_BLK, dz_BLK
  use ModAMR, ONLY: child2subface
  implicit none

  ! Coarse face centered B field components
  real, intent(in),  dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn) :: &
       Bxf_c, Byf_c, Bzf_c

  ! Did we have finer neighbors before prolongation
  logical, intent(in) :: IsFinerNei_E(east_:top_)

  ! Normal B components from finer neighbors 
  ! on the shared subfaces (index Q) on two sides (index S)
  real, intent(in) :: BxFaceFine_XQS(1:nJ,1:nK,4,2)
  real, intent(in) :: ByFaceFine_YQS(1:nI,1:nK,4,2)
  real, intent(in) :: BzFaceFine_ZQS(1:nI,1:nJ,4,2)

  ! The child index relative to the coarse parent
  integer, intent(in) :: iChild

  ! The block number of the fine block into which the prolongation is done
  integer, intent(in) :: iBlock

  ! Fine face centered B field components produced by prolongation
  real, intent(out), dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn) :: &
       Bxf_f, Byf_f, Bzf_f

  integer :: i,j,k,i1,j1,k1,ishift,jshift,kshift
  integer :: iP,iM,jP,jM,kP,kM
  real :: gradXl, gradXr, gradYl, gradYr, gradZl, gradZr
  real :: dBxdy, dBxdz, dBydx, dBydz, dBzdx, dBzdy
  real :: dBxdxx, dBydyy, dBzdzz, dBxdxyz, dBydxyz, dBzdxyz

  ! aspect ratio related constants
  real :: Dx,Dy,Dz,DxDy8Inv,DxDz8Inv,DyDx8Inv,DyDz8Inv,DzDx8Inv,DzDy8Inv
  real :: Dx2,Dy2,Dz2,Dx2Dxy2,Dy2Dxy2,Dx2Dxz2,Dz2Dxz2,Dy2Dyz2,Dz2Dyz2
  real :: Dy2DBxDxyz,Dz2DBxDxyz, Dx2DByDxyz,Dz2DByDxyz, Dx2DBzDxyz,Dy2DBzDxyz

  logical :: oktest, oktest_me
  !--------------------------------------------------------------------------

  call get_shifts(iChild,ishift,jshift,kshift)

  ! Calculate aspect ratios for a non-cubic cell
  Dx = dx_BLK(iBlock); Dy = dy_BLK(iBlock); Dz = dz_BLK(iBlock)
  DxDy8Inv = Dx/(8*Dy); DxDz8Inv = Dx/(8*Dz)
  DyDx8Inv = Dy/(8*Dx); DyDz8Inv = Dy/(8*Dz)
  DzDx8Inv = Dz/(8*Dx); DzDy8Inv = Dz/(8*Dy)

  Dx2 = Dx**2; Dy2 = Dy**2; Dz2 = Dz**2
  Dx2Dxy2 = Dx2/(Dx2+Dy2); Dy2Dxy2 = Dy2/(Dx2+Dy2)
  Dx2Dxz2 = Dx2/(Dx2+Dz2); Dz2Dxz2 = Dz2/(Dx2+Dz2)
  Dy2Dyz2 = Dy2/(Dy2+Dz2); Dz2Dyz2 = Dz2/(Dy2+Dz2)

  if(iBlock==BLKtest.and.&
       iShift<iTest.and.iTest<=iShift+nI/2.and. &
       jShift<jTest.and.jTest<=jShift+nJ/2.and. &
       kShift<kTest.and.kTest<=kShift+nK/2)then
     call set_oktest('prolong2_bface',oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  end if

  if(oktest_me)write(*,*)'prolong2_bface: iChild,iShift,jShift,kShift=',&
       iChild,iShift,jShift,kShift

  ! Assign default solution state to fine block to get corners
  Bxf_f=0.0; Byf_f=0.0; Bzf_f=0.0

  ! X faces
  do i = 1+ishift, nI/2+ishift+1
     do j = 1+jshift, nJ/2+jshift; jP=min(j+1,nJ); jM=max(j-1,1)
        do k = 1+kshift, nK/2+kshift; kP=min(k+1,nK); kM=max(k-1,1)
           i1 = 2*(i-ishift)-1; j1 = 2*(j-jshift)-1; k1 = 2*(k-kshift)-1
           ! Second order minmod limited interpolation on coarse cell faces

           gradYr = Bxf_c(i,jP ,k)-Bxf_c(i,j  ,k)
           gradYl = Bxf_c(i,j  ,k)-Bxf_c(i,jM ,k)
           dBxdy  = sign(0.25,gradyl)*&
                max(0.,min(abs(gradyl),sign(1.,gradyl)*gradyr))

           gradZr = Bxf_c(i,j,kP) - Bxf_c(i,j,k)
           gradZl = Bxf_c(i,j,k)  - Bxf_c(i,j,kM)
           dBxdz  = sign(0.25,gradzl)*&
                max(0.,min(abs(gradzl),sign(1.,gradzl)*gradzr))

           Bxf_f(i1,j1  ,k1  )   = Bxf_c(i,j,k) - dBxdy - dBxdz
           Bxf_f(i1,j1+1,k1  )   = Bxf_c(i,j,k) + dBxdy - dBxdz
           Bxf_f(i1,j1  ,k1+1)   = Bxf_c(i,j,k) - dBxdy + dBxdz
           Bxf_f(i1,j1+1,k1+1)   = Bxf_c(i,j,k) + dBxdy + dBxdz

        end do
     end do
  end do

  ! Y faces
  do i = 1+ishift, nI/2+ishift; iP=min(i+1,nI); iM=max(i-1,1)
     do j = 1+jshift, nJ/2+jshift+1
        do k = 1+kshift, nK/2+kshift; kP=min(k+1,nK); kM=max(k-1,1)
           i1 = 2*(i-ishift)-1; j1 = 2*(j-jshift)-1; k1 = 2*(k-kshift)-1

           gradXr = Byf_c(iP,j,k) - Byf_c(i ,j,k)
           gradXl = Byf_c(i ,j,k) - Byf_c(iM,j,k)
           dBydx  = sign(0.25,gradxl)*&
                max(0.,min(abs(gradxl),sign(1.,gradxl)*gradxr))

           gradZr = Byf_c(i,j,kP) - Byf_c(i,j,k)
           gradZl = Byf_c(i,j,k)  - Byf_c(i,j,kM)
           dBydz  = sign(0.25,gradzl)*&
                max(0.,min(abs(gradzl),sign(1.,gradzl)*gradzr))

           Byf_f(i1  ,j1,k1  )   = Byf_c(i,j,k) - dBydx - dBydz
           Byf_f(i1+1,j1,k1  )   = Byf_c(i,j,k) + dBydx - dBydz
           Byf_f(i1  ,j1,k1+1)   = Byf_c(i,j,k) - dBydx + dBydz
           Byf_f(i1+1,j1,k1+1)   = Byf_c(i,j,k) + dBydx + dBydz

        end do
     end do
  end do

  ! Z faces
  do i = 1+ishift, nI/2+ishift; iP=min(i+1,nI); iM=max(i-1,1)
     do j = 1+jshift, nJ/2+jshift; jP=min(j+1,nJ); jM=max(j-1,1)
        do k = 1+kshift, nK/2+kshift+1
           i1 = 2*(i-ishift)-1; j1 = 2*(j-jshift)-1; k1 = 2*(k-kshift)-1

           gradXr = Bzf_c(iP,j,k) - Bzf_c(i ,j,k)
           gradXl = Bzf_c(i ,j,k) - Bzf_c(iM,j,k)
           dBzdx  = sign(0.25,gradxl)*&
                max(0.,min(abs(gradxl),sign(1.,gradxl)*gradxr))

           gradYr = Bzf_c(i,jP,k) - Bzf_c(i,j ,k)
           gradYl = Bzf_c(i,j ,k) - Bzf_c(i,jM,k)
           dBzdy  = sign(0.25,gradyl)*&
                max(0.,min(abs(gradyl),sign(1.,gradyl)*gradyr))

           Bzf_f(i1  ,j1  ,k1)   = Bzf_c(i,j,k) - dBzdx - dBzdy
           Bzf_f(i1+1,j1  ,k1)   = Bzf_c(i,j,k) + dBzdx - dBzdy
           Bzf_f(i1  ,j1+1,k1)   = Bzf_c(i,j,k) - dBzdx + dBzdy
           Bzf_f(i1+1,j1+1,k1)   = Bzf_c(i,j,k) + dBzdx + dBzdy

        end do
     end do
  end do

  if(oktest_me)then
     ! Check if the interpolated fine B fluxes add up to the coarse B flux
     i1 = 2*(iTest-ishift)-1; j1 = 2*(jTest-jshift)-1; k1 = 2*(kTest-kshift)-1
     write(*,*)'Before correction'
     write(*,*)'Bx_c, avg Bx_f(-)=',Bxf_c(iTest,jTest,kTest),&
          sum(Bxf_f(i1,j1:j1+1,k1:k1+1))/4
     write(*,*)'Bx_c, avg Bx_f(+)=',Bxf_c(iTest+1,jTest,kTest),&
          sum(Bxf_f(i1+2,j1:j1+1,k1:k1+1))/4
     write(*,*)'By_c, avg By_f(-)=',Byf_c(iTest,jTest,kTest),&
          sum(Byf_f(i1:i1+1,j1,k1:k1+1))/4
     write(*,*)'By_c, avg By_f(+)=',Byf_c(iTest,jTest+1,kTest),&
          sum(Byf_f(i1:i1+1,j1+2,k1:k1+1))/4
     write(*,*)'Bz_c, avg Bz_f(-)=',Bzf_c(iTest,jTest,kTest),&
          sum(Bzf_f(i1:i1+1,j1:j1+1,k1))/4
     write(*,*)'Bz_c, avg Bz_f(+)=',Bzf_c(iTest,jTest,kTest+1),&
          sum(Bzf_f(i1:i1+1,j1:j1+1,k1+2))/4
     write(*,*)'divB_c=',&
          (Bxf_c(iTest+1,jTest,kTest)-Bxf_c(iTest,jTest,kTest))/Dx+&
          (Byf_c(iTest,jTest+1,kTest)-Byf_c(iTest,jTest,kTest))/Dy+&
          (Bzf_c(iTest,jTest,kTest+1)-Bzf_c(iTest,jTest,kTest))/Dz
  end if

  ! Correct normal components on faces which were shared with a finer block
  ! before the AMR so that we get a consistent normal flux
  if(IsFinerNei_E(east_).and.iShift==0)  Bxf_f(   1,1:nJ,1:nK)=&
       BxFaceFine_XQS(:,:,child2subface(iChild,east_),1)

  if(IsFinerNei_E(west_).and.iShift>0)   Bxf_f(nI+1,1:nJ,1:nK)=&
       BxFaceFine_XQS(:,:,child2subface(iChild,west_),2)

  if(IsFinerNei_E(south_).and.jShift==0) Byf_f(1:nI,   1,1:nK)=&
       ByFaceFine_YQS(:,:,child2subface(iChild,south_),1)

  if(IsFinerNei_E(north_).and.jShift>0)  Byf_f(1:nI,nJ+1,1:nK)=&
       ByFaceFine_YQS(:,:,child2subface(iChild,north_),2)

  if(IsFinerNei_E(bot_).and.kShift==0)   Bzf_f(1:nI,1:nJ,1)=&
       BzFaceFine_ZQS(:,:,child2subface(iChild,bot_),1)

  if(IsFinerNei_E(top_).and.kShift>0)    Bzf_f(1:nI,1:nJ,nK+1)=&
       BzFaceFine_ZQS(:,:,child2subface(iChild,top_),2)

  if(oktest_me)then
     ! Check if the corrected fine B fluxes add up to the coarse B flux
     i1 = 2*(iTest-ishift)-1; j1 = 2*(jTest-jshift)-1; k1 = 2*(kTest-kshift)-1
     write(*,*)'After correction'
     write(*,*)'Bx_c, avg Bx_f(-)=',Bxf_c(iTest,jTest,kTest),&
          sum(Bxf_f(i1,j1:j1+1,k1:k1+1))/4
     write(*,*)'Bx_c, avg Bx_f(+)=',Bxf_c(iTest+1,jTest,kTest),&
          sum(Bxf_f(i1+2,j1:j1+1,k1:k1+1))/4
     write(*,*)'By_c, avg By_f(-)=',Byf_c(iTest,jTest,kTest),&
          sum(Byf_f(i1:i1+1,j1,k1:k1+1))/4
     write(*,*)'By_c, avg By_f(+)=',Byf_c(iTest,jTest+1,kTest),&
          sum(Byf_f(i1:i1+1,j1+2,k1:k1+1))/4
     write(*,*)'Bz_c, avg Bz_f(-)=',Bzf_c(iTest,jTest,kTest),&
          sum(Bzf_f(i1:i1+1,j1:j1+1,k1))/4
     write(*,*)'Bz_c, avg Bz_f(+)=',Bzf_c(iTest,jTest,kTest+1),&
          sum(Bzf_f(i1:i1+1,j1:j1+1,k1+2))/4
  end if

  ! Do central faces of coarse cells according to Toth and Roe paper
  do i = 1+ishift, nI/2+ishift
     do j = 1+jshift, nJ/2+jshift
        do k = 1+kshift, nK/2+kshift
           i1 = 2*(i-ishift)-1; j1 = 2*(j-jshift)-1; k1 = 2*(k-kshift)-1

           ! Second order derivatives, correct for face areas

           ! dBxdxx = 1/8 * sum i j By + i k Bz
           dBxdxx = &
                DxDy8Inv*( &
                + Byf_f(i1  ,j1  ,k1  ) &
                - Byf_f(i1+1,j1  ,k1  ) &
                + Byf_f(i1  ,j1  ,k1+1) &
                - Byf_f(i1+1,j1  ,k1+1) &
                - Byf_f(i1  ,j1+2,k1  ) &
                + Byf_f(i1+1,j1+2,k1  ) &
                - Byf_f(i1  ,j1+2,k1+1) &
                + Byf_f(i1+1,j1+2,k1+1)) + &
                DxDz8Inv*(&
                + Bzf_f(i1  ,j1  ,k1  ) &
                - Bzf_f(i1+1,j1  ,k1  ) &
                + Bzf_f(i1  ,j1+1,k1  ) &
                - Bzf_f(i1+1,j1+1,k1  ) &
                - Bzf_f(i1  ,j1  ,k1+2) &
                + Bzf_f(i1+1,j1  ,k1+2) &
                - Bzf_f(i1  ,j1+1,k1+2) &
                + Bzf_f(i1+1,j1+1,k1+2))

           ! dBydyy = 1/8 * sum j i Bx + k j Bz
           dBydyy = &
                DyDx8Inv*( &
                + Bxf_f(i1  ,j1  ,k1  ) &
                - Bxf_f(i1  ,j1+1,k1  ) &
                + Bxf_f(i1  ,j1  ,k1+1) &
                - Bxf_f(i1  ,j1+1,k1+1) &
                - Bxf_f(i1+2,j1  ,k1  ) &
                + Bxf_f(i1+2,j1+1,k1  ) &
                - Bxf_f(i1+2,j1  ,k1+1) &
                + Bxf_f(i1+2,j1+1,k1+1)) + &
                DyDz8Inv*( &
                + Bzf_f(i1  ,j1  ,k1  ) &
                - Bzf_f(i1  ,j1+1,k1  ) &
                + Bzf_f(i1+1,j1  ,k1  ) &
                - Bzf_f(i1+1,j1+1,k1  ) &
                - Bzf_f(i1  ,j1  ,k1+2) &
                + Bzf_f(i1  ,j1+1,k1+2) &
                - Bzf_f(i1+1,j1  ,k1+2) &
                + Bzf_f(i1+1,j1+1,k1+2))

           ! dBzdzz = 1/8 * sum k i Bx + k j By
           dBzdzz = &
                DzDx8Inv*( &
                + Bxf_f(i1  ,j1  ,k1  ) &
                - Bxf_f(i1  ,j1  ,k1+1) &
                + Bxf_f(i1  ,j1+1,k1  ) &
                - Bxf_f(i1  ,j1+1,k1+1) &
                - Bxf_f(i1+2,j1  ,k1  ) &
                + Bxf_f(i1+2,j1  ,k1+1) &
                - Bxf_f(i1+2,j1+1,k1  ) &
                + Bxf_f(i1+2,j1+1,k1+1)) + &
                DzDy8Inv*( &
                + Byf_f(i1  ,j1  ,k1  ) &
                - Byf_f(i1  ,j1  ,k1+1) &
                + Byf_f(i1+1,j1  ,k1  ) &
                - Byf_f(i1+1,j1  ,k1+1) &
                - Byf_f(i1  ,j1+2,k1  ) &
                + Byf_f(i1  ,j1+2,k1+1) &
                - Byf_f(i1+1,j1+2,k1  ) &
                + Byf_f(i1+1,j1+2,k1+1))

           ! Third order derivatives
           ! These are zero unless a fine B face correction was done
           ! aspect ratios are taken into account

           ! dBxdxyz = 1/(8(dy^2+dz^2)) sum i j k Bx
           dBxdxyz= &
                -Bxf_f(i1  ,j1  ,k1  ) &
                +Bxf_f(i1  ,j1+1,k1  ) &
                +Bxf_f(i1  ,j1  ,k1+1) &
                -Bxf_f(i1  ,j1+1,k1+1) &
                +Bxf_f(i1+2,j1  ,k1  ) &
                -Bxf_f(i1+2,j1+1,k1  ) &
                -Bxf_f(i1+2,j1  ,k1+1) &
                +Bxf_f(i1+2,j1+1,k1+1)

           Dy2DBxDxyz = DzDx8Inv*Dy2Dyz2*DBxDxyz 
           Dz2DBxDxyz = DyDx8Inv*Dz2Dyz2*DBxDxyz

           ! dBydxyz = 1/(8(dx^2+dz^2)) sum i j k By
           dBydxyz= &
                -Byf_f(i1  ,j1  ,k1  ) &
                +Byf_f(i1  ,j1+2,k1  ) &
                +Byf_f(i1  ,j1  ,k1+1) &
                -Byf_f(i1  ,j1+2,k1+1) &
                +Byf_f(i1+1,j1  ,k1  ) &
                -Byf_f(i1+1,j1+2,k1  ) &
                -Byf_f(i1+1,j1  ,k1+1) &
                +Byf_f(i1+1,j1+2,k1+1)

           Dx2DByDxyz = DzDy8Inv*Dx2Dxz2*DByDxyz
           Dz2DByDxyz = DxDy8Inv*Dz2Dxz2*DByDxyz

           ! dBzdxyz = 1/(8(dy^2+dz^2)) sum i j k Bz
           dBzdxyz= &
                -Bzf_f(i1  ,j1  ,k1  ) &
                +Bzf_f(i1  ,j1+1,k1  ) &
                +Bzf_f(i1  ,j1  ,k1+2) &
                -Bzf_f(i1  ,j1+1,k1+2) &
                +Bzf_f(i1+1,j1  ,k1  ) &
                -Bzf_f(i1+1,j1+1,k1  ) &
                -Bzf_f(i1+1,j1  ,k1+2) &
                +Bzf_f(i1+1,j1+1,k1+2)

           Dx2DBzDxyz = DyDz8Inv*Dx2Dxy2*DBzDxyz
           Dy2DBzDxyz = DxDz8Inv*Dy2Dxy2*DBzDxyz

           ! Calculate internal fine solution

           ! Bx = (Bxp+Bxm)/2 + dBxdxx + k dz^2 dBydxyz + j dy^2 dBzdxyz
           Bxf_f(i1+1,j1  ,k1  )=0.5*(&
                Bxf_f(i1  ,j1  ,k1  )+ &
                Bxf_f(i1+2,j1  ,k1  ))  + dBxdxx - Dz2DByDxyz - Dy2DBzDxyz

           Bxf_f(i1+1,j1+1,k1  )=0.5*(&
                Bxf_f(i1  ,j1+1,k1  )+ &
                Bxf_f(i1+2,j1+1,k1  ))  + dBxdxx - Dz2DByDxyz + Dy2DBzDxyz

           Bxf_f(i1+1,j1  ,k1+1)=0.5*(&
                Bxf_f(i1  ,j1  ,k1+1)+ &
                Bxf_f(i1+2,j1  ,k1+1))  + dBxdxx + Dz2DByDxyz - Dy2DBzDxyz

           Bxf_f(i1+1,j1+1,k1+1)=0.5*(&
                Bxf_f(i1  ,j1+1,k1+1)+ &
                Bxf_f(i1+2,j1+1,k1+1))  + dBxdxx + Dz2DByDxyz + Dy2DBzDxyz

           ! By = (Byp+Bym)/2 + dBydyy + i dx^2 dBzdxyz + k dz^2 dBxdxyz
           Byf_f(i1  ,j1+1,k1  )=0.5*(&
                Byf_f(i1  ,j1  ,k1  )+ &
                Byf_f(i1  ,j1+2,k1  ))  + dBydyy - Dx2DBzDxyz - Dz2DBxDxyz

           Byf_f(i1+1,j1+1,k1  )=0.5*(&
                Byf_f(i1+1,j1  ,k1  )+ &
                Byf_f(i1+1,j1+2,k1  ))  + dBydyy + Dx2DBzDxyz - Dz2DBxDxyz

           Byf_f(i1  ,j1+1,k1+1)=0.5*(&
                Byf_f(i1  ,j1  ,k1+1)+ &
                Byf_f(i1  ,j1+2,k1+1))  + dBydyy - Dx2DBzDxyz + Dz2DBxDxyz

           Byf_f(i1+1,j1+1,k1+1)=0.5*(&
                Byf_f(i1+1,j1  ,k1+1)+ &
                Byf_f(i1+1,j1+2,k1+1))  + dBydyy + Dx2DBzDxyz + Dz2DBxDxyz

           ! Bz = (Bzp+Bzm)/2 + dBzdzz + j dy^2 dBxdxyz + i dx^2 dBydxyz
           Bzf_f(i1  ,j1  ,k1+1)=0.5*(&
                Bzf_f(i1  ,j1  ,k1  )+ &
                Bzf_f(i1  ,j1  ,k1+2))  + dBzdzz - Dy2DBxDxyz - Dx2DByDxyz

           Bzf_f(i1+1,j1  ,k1+1)=0.5*(&
                Bzf_f(i1+1,j1  ,k1  )+ &
                Bzf_f(i1+1,j1  ,k1+2))  + dBzdzz - Dy2DBxDxyz + Dx2DByDxyz

           Bzf_f(i1  ,j1+1,k1+1)=0.5*(&
                Bzf_f(i1  ,j1+1,k1  )+ &
                Bzf_f(i1  ,j1+1,k1+2))  + dBzdzz + Dy2DBxDxyz - Dx2DByDxyz

           Bzf_f(i1+1,j1+1,k1+1)=0.5*(&
                Bzf_f(i1+1,j1+1,k1  )+ &
                Bzf_f(i1+1,j1+1,k1+2))  + dBzdzz + Dy2DBxDxyz + Dx2DByDxyz

        end do
     end do
  end do

  if(oktest_me)then
     ! Check the divergence B condition for the prolonged cells
     i1 = 2*(iTest-iShift)-1; j1 = 2*(jTest-jshift)-1; k1 = 2*(kTest-kshift)-1
     do i=i1,i1+1; do j=j1,j1+1; do k=k1,k1+1
        write(*,*)'divB_f(',i,',',j,',',k,')=',&
             (Bxf_f(i+1,j,k)-Bxf_f(i,j,k))/Dx+&
             (Byf_f(i,j+1,k)-Byf_f(i,j,k))/Dy+&
             (Bzf_f(i,j,k+1)-Bzf_f(i,j,k))/Dz
     end do; end do; end do
  end if

end subroutine prolong_b_face
!============================================================================
subroutine assign_coarse_face_soln(sol_BLK,iVar) 
  use ModProcMH
  use ModSize
  use ModAMR, ONLY:local_cube,local_cubeBLK
  use ModMpi
  implicit none

  integer,intent(in)::iVar
  real, intent(inout), dimension (1-gcn:nI+gcn, &
       1-gcn:nJ+gcn, &
       1-gcn:nK+gcn,nBLK) :: sol_BLK

  real, dimension(1:nI/2+1, 1:nJ/2+1, 1:nK/2+1, 8) ::&
       restricted_soln_blks

  integer::remaining_PE,remaining_BLK,iCube

  integer,parameter :: isize=(nI/2+1)*(nJ/2+1)*(nK/2+1)
  integer :: iTag,iError, number_send_requests, send_requests(7)
  integer::receive_requests(7), number_receive_requests, &
       status(MPI_STATUS_SIZE, 7)

  number_send_requests = 0
  remaining_PE=local_cube(1)
  remaining_BLK=local_cubeBLK(1)

  do icube = 1, 8
     if (iProc == local_cube(icube)) then
        call restrict_Bface(sol_BLK(:,:,:,local_cubeBLK(icube)),&
             iVar, restricted_soln_blks(:,:,:,icube))

        if (icube > 1 .and. iProc .ne. remaining_PE) then
           itag = local_cubeBLK(icube)*100 + iVar
           number_send_requests = number_send_requests + 1
           call MPI_isend(restricted_soln_blks(1,1,1,icube), &
                isize, MPI_REAL, remaining_PE,itag, iComm,&
                send_requests(number_send_requests), iError)
        end if
     end if
  end do

  if (number_send_requests > 0) then
     call MPI_waitall(number_send_requests, &
          send_requests(1), &
          status(1,1), iError)
  end if

  number_receive_requests = 0

  if (iProc == remaining_PE) then ! remaining coarse block
     do icube = 2, 8
        if (local_cube(icube) .ne. remaining_PE) then
           itag = local_cubeBLK(icube)*100 + iVar
           number_receive_requests = number_receive_requests + 1
           call MPI_irecv(restricted_soln_blks(1,1,1,icube), &
                isize, MPI_REAL, local_cube(icube), itag, iComm,&
                receive_requests(number_receive_requests), iError)
        end if
     end do

     if (number_receive_requests > 0) then
        call MPI_waitall(number_receive_requests, &
             receive_requests(1), &
             status(1,1), iError)

     end if
     call assign_restricted_Bface(restricted_soln_blks, iVar, &
          sol_BLK(:,:,:,remaining_BLK))
  end if ! remaining coarse block

end subroutine assign_coarse_face_soln

!=============================================================================

subroutine restrict_Bface(fine_sol,iVar,coarse_sol)

  use ModSize
  use ModVarIndexes, ONLY : Bx_,By_,Bz_
  implicit none

  integer, intent(in) :: iVar
  real, intent(in) :: fine_sol(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn)
  real, intent(out):: coarse_sol(1:nI/2+1,1:nJ/2+1,1:nK/2+1)

  !---------------------------------------------------------------------------

  select case(iVar)
  case(Bx_)
     coarse_sol(1:nI/2+1,1:(nJ+1)/2,1:(nK+1)/2)=0.25*(&
          fine_sol(1:nI+1:2, 1:nJ  :2, 1:nK  :2)+&
          fine_sol(1:nI+1:2, 2:nJ+1:2, 1:nK  :2)+&
          fine_sol(1:nI+1:2, 1:nJ  :2, 2:nK+1:2)+&
          fine_sol(1:nI+1:2, 2:nJ+1:2, 2:nK+1:2))
  case(By_)
     coarse_sol(1:nI/2,1:nJ/2+1,1:(nK+1)/2)=0.25*(&
          fine_sol(1:nI:2, 1:nJ+1:2, 1:nK  :2)+&
          fine_sol(2:nI:2, 1:nJ+1:2, 1:nK  :2)+&
          fine_sol(1:nI:2, 1:nJ+1:2, 2:nK+1:2)+&
          fine_sol(2:nI:2, 1:nJ+1:2, 2:nK+1:2))
  case(Bz_)
     coarse_sol(1:nI/2,1:(nJ+1)/2,1:nK/2+1)=0.25*(&
          fine_sol(1:nI:2, 1:nJ  :2, 1:nK+1:2)+&
          fine_sol(2:nI:2, 1:nJ  :2, 1:nK+1:2)+&
          fine_sol(1:nI:2, 2:nJ+1:2, 1:nK+1:2)+&
          fine_sol(2:nI:2, 2:nJ+1:2, 1:nK+1:2))
  case default
     call stop_mpi('Invalid iVar in restrict_Bface')
  end select

end subroutine restrict_Bface

!==============================================================================

subroutine assign_restricted_Bface(r_sol,iVar,coarse_sol)

  use ModSize
  use ModVarIndexes, ONLY : Bx_,By_,Bz_
  implicit none

  real, intent(in)   :: r_sol(1:nI/2+1,1:nJ/2+1,1:nK/2+1,8)
  integer, intent(in):: iVar
  real, intent(out)  :: coarse_sol(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn)

  integer :: dI, dJ, dK

  !---------------------------------------------------------------------------

  ! Assign default for corners
  coarse_sol=0.0

  select case(iVar)
  case(Bx_)
     dI=1; dJ=0; dK=0;
  case(By_)
     dI=0; dJ=1; dK=0;
  case(Bz_)
     dI=0; dJ=0; dK=1;
  end select

  coarse_sol(1:nI/2      ,1:nJ/2      ,nK/2+1:nK+dK)= &  ! 001
       r_sol(1:nI/2      ,1:nJ/2      ,1:nK/2+dK   ,1)

  coarse_sol(nI/2+1:nI+dI,1:nJ/2      ,nK/2+1:nK+dK)= &  ! 101
       r_sol(1:nI/2+dI   ,1:nJ/2      ,1:nK/2+dK   ,2)

  coarse_sol(nI/2+1:nI+dI,1:nJ/2      ,1:nK/2      )= &  ! 100
       r_sol(1:nI/2+dI   ,1:nJ/2      ,1:nK/2      ,3)

  coarse_sol(1:nI/2      ,1:nJ/2      ,1:nK/2      )= &  ! 000
       r_sol(1:nI/2      ,1:nJ/2      ,1:nK/2      ,4)

  coarse_sol(1:nI/2      ,nJ/2+1:nJ+dJ,1:nK/2      )= &  ! 010
       r_sol(1:nI/2      ,1:nJ/2+dJ   ,1:nK/2      ,5)

  coarse_sol(nI/2+1:nI+dI,nJ/2+1:nJ+dJ,1:nK/2      )= &  ! 110
       r_sol(1:nI/2+dI   ,1:nJ/2+dJ   ,1:nK/2      ,6)

  coarse_sol(nI/2+1:nI+dI,nJ/2+1:nJ+dJ,nK/2+1:nK+dK)= &  ! 111
       r_sol(1:nI/2+dI   ,1:nJ/2+dJ   ,1:nK/2+dK   ,7)

  coarse_sol(1:nI/2      ,nJ/2+1:nJ+dJ,nK/2+1:nK+dK)= &  ! 011
       r_sol(1:nI/2      ,1:nJ/2+dJ   ,1:nK/2+dK   ,8)

end subroutine assign_restricted_Bface

!==============================================================================

subroutine constrain_ICs(iBlock)

  ! Initialize B field, Bface, and Bcenter for Constrained Transport

  use ModMain
  use ModVarIndexes
  use ModAdvance, ONLY : State_VGB
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK,body_BLK,true_cell
  use ModIO, ONLY : restart
  use ModPhysics, ONLY : SW_Bx,SW_By,SW_Bz
  use ModCT, ONLY : Bxface_BLK,Byface_BLK,Bzface_BLK
  implicit none

  integer, intent(in) :: iBlock
  !---------------------------------------------------------------------------

  if(unusedBLK(iBlock))then
     BxFace_BLK(:,:,:,iBlock)=0.0
     ByFace_BLK(:,:,:,iBlock)=0.0
     BzFace_BLK(:,:,:,iBlock)=0.0
  else

     if(.not.restart)then
        where(x_BLK(:,:,:,iBlock)<16.)
           ! Cancel B field at x<16Re to avoid non-zero initial divB
           ! x=16 is a good choice because it is a power of 2 so it is 
           ! a block boundary for all block sizes. 
           ! x=16 is larger than typical rBody.
           State_VGB(Bx_,:,:,:,iBlock)=0.0
           State_VGB(By_,:,:,:,iBlock)=0.0
           State_VGB(Bz_,:,:,:,iBlock)=0.0
           ! Balance total pressure
           State_VGB(P_,:,:,:,iBlock)=State_VGB(P_,:,:,:,iBlock)+ &
                0.5*(SW_Bx**2+SW_By**2+SW_Bz**2)
        elsewhere
           ! Use solar wind values ahead of the Earth
           State_VGB(Bx_,:,:,:,iBlock)=SW_Bx
           State_VGB(By_,:,:,:,iBlock)=SW_By
           State_VGB(Bz_,:,:,:,iBlock)=SW_Bz
        end where
     end if


     if(index(test_string,'testCTcoarse')>0)then
        State_VGB(Bx_,:,:,:,iBlock)=   x_BLK(:,:,:,iBlock)
        State_VGB(By_,:,:,:,iBlock)=   y_BLK(:,:,:,iBlock)
        State_VGB(Bz_,:,:,:,iBlock)=-2*z_BLK(:,:,:,iBlock)
        if(body_BLK(iBlock))then
           where(.not.true_cell(:,:,:,iBlock))
              State_VGB(Bx_,:,:,:,iBlock)=0.
              State_VGB(By_,:,:,:,iBlock)=0.
              State_VGB(Bz_,:,:,:,iBlock)=0.
           end where
        end if
     end if

  endif

end subroutine constrain_ICs

!=============================================================================
subroutine correct_VxB

  ! At refinement level changes correct VxB on coarse cell edge to the
  ! average of the 2 fine edge values. This requires the message passing
  ! of restricted faces and (in certain cases) edges

  ! For example the 4x4 block face orthogonal to Z is restricted like this:
  !
  !  +-Ex-+-Ex-+-Ex-+-Ex-+              +----Ex---+----Ex---+
  !  |    |    |    |    |              |         |         |
  !  Ey   Ey   Ey   Ey   Ey             |         |         |
  !  |    |    |    |    |              |         |         |
  !  +-Ex-+-Ex-+-Ex-+-Ex-+             Ey         Ey        Ey
  !  |    |    |    |    |              |         |         |
  !  Ey   Ey   Ey   Ey   Ey             |         |         |
  !  |    |    |    |    |              |         |         |
  !  +-Ex-+-Ex-+-Ex-+-Ex-+      --->    +----Ex---+----Ex---+
  !  |    |    |    |    |              |         |         |
  ! *Ey   Ey   Ey   Ey   Ey             |         |         |
  !  |    |    |    |    |              |         |         |
  !  +-Ex-+-Ex-+-Ex-+-Ex-+            *Ey         Ey        Ey
  !  |    |    |    |    |              |         |         |
  ! *Ey   Ey   Ey   Ey   Ey             |         |         |
  !  |    |    |    |    |              |         |         |
  !  +-Ex-+-Ex-+-Ex-+-Ex-+              +----Ex---+----Ex---+
  !
  ! Note that the Ex values are restricted from a 4x5 to 2x3,
  !     while the Ey values are restricted from a 5x4 to 3x2.
  !
  ! We denote Ex as orientation 1 and Ey as orientation 2.
  ! In general the array corresponding to orientation 1 is size 
  ! nxface*(nyface+1) while orientation 2 has size (nxface+1)*nyface,
  ! where nxface and nyface are half of the number of cells
  ! in the first and second directions (ordered as X,Y,Z).
  ! The two orientations are message passed in separate messages since
  ! their shapes and even sizes can be different.
  !
  ! In 3D it is possible that a coarse and a fine block share an
  ! edge without sharing a face. In that case the edge has to be
  ! message passed, but this only involves the orientation parallel
  ! to the edge. These messages and the corresponding variables are
  ! denoted by 3.
  !
  ! The restriction operator is a simple average from the two fine edge
  ! values to the one coarse edge value. See the *Ey-s in the figure.
  ! This conserves magnetic flux as required.
  !
  ! Each fine block can send 2 orientations per face, while a coarse block
  ! may receive 2 orientations per subface, i.e. 8 messages per face.
  ! Messages are tagged by the receiving block number and subface number,
  ! and indexed in the receive buffers accordingly.
  !
  ! A fine block has four edges per face. However, each edge belongs to two
  ! faces, so we can relate 2 edges to a face. In particular, we assign to
  ! the face orthogonal to X  the edges parallel to Y, to Y faces the Z edges
  ! and to Z faces the X edges. The two edges are indexed by 1 and 2 in the 
  ! order of increasing coordinates. Furthermore only one of the two edges 
  ! assigned to a fine face can be shared with a coarse block in a diagonal 
  ! direction. Thus only one edge message per face can be sent.
  !
  ! The coarse block can receive at most 2 subedges for the 2 edges belonging
  ! to a face, so at most 4 messages can arrive. Messages are tagged and 
  ! indexed in the receive buffer rbuf3 by the receiving block number, the
  ! edge and subedge indices. The subedges are also indexed by 1 and 2 
  ! with increasing coordinates. E.g. for a Z face the 2 edges parallel to X 
  ! are indexed and divided like this:
  !
  !  Y
  !  ^
  !  |
  !  |
  !  +---sub1---edge2-sub2---+
  !  |     |     |     |     |
  !  |     |     |     |     |
  !  |     |     |     |     |
  !  +-----+-----+-----+-----+
  !  |     |     |     |     |
  !  |     |     |     |     |
  !  |     |     |     |     |
  !  +-----+-----+-----+-----+
  !  |     |     |     |     |
  !  |     |     |     |     |
  !  |     |     |     |     |
  !  +-----+-----+-----+-----+
  !  |     |     |     |     |
  !  |     |     |     |     |
  !  |     |     |     |     |
  !  +---sub1---edge1-sub2---+ --> X

  use ModProcMH
  use ModMain
  use ModParallel, ONLY : neiLEV,neiBLK,neiPE, &
       BLKneighborPE,BLKneighborBLK,BLKneighborLEV,BLKneighborCHILD
  use ModCT, ONLY : VxB_x,VxB_y,VxB_z
  use ModAMR, ONLY : child2subface,child2subedge
  use ModMpi
  implicit none

  ! Local variables

  ! facedir=1,2,3 correspond to east-west, south-north, bot-top.
  integer :: isweep, facedir, edgedir, sidedir

  ! Index sent and received face (east=1,..,top=6) and recv edge (min=1,max=2)
  integer :: iface, rface, iedge

  ! Fixed coord index (1 or n?+1) for sent and received faces and edges
  integer :: isface, irface, isedge, iredge

  ! Restricted face sizes (n?/2)
  integer :: nxface, nyface

  ! subedge (1..2), subface (1..4) and child (1..8) index
  integer :: isubedge, isubface, ichild

  ! Block index (1..nBLK)
  integer :: iBLK

  ! Descriptors for neighbor face and neighbor subedges
  integer :: neiP,neiB,neiL,neiedgeP(2),neiedgeB(2)

  ! MPI variables
  integer :: itag, request, number_receive_requests, receive_requests(nBLK*72)
  integer :: status(MPI_STATUS_SIZE, nBLK*72)

  ! Maximum sizes of RESTRICTED VxB layers to be received for two orientations
  integer, parameter :: maxsize1= max(&
       nI*(nJ+2)/4,nI*(nK+2)/4,nJ*(nK+2)/4)

  integer, parameter :: maxsize2= max(&
       (nI+2)*nJ/4,(nI+2)*nK/4,(nJ+2)*nK/4)


  ! Maximum size of a restricted block edge, only aligned VxB orientation
  integer, parameter :: maxsize3=max(nI/2,nJ/2,nK/2)

  ! Receive buffers to hold 4 incoming RESTRICTED subface/subedge values
  ! for all blocks and 6 faces, and for the 2 orientations and the edge
  real, dimension(maxsize1,4  ,nBLK,6) :: rbuf1
  real, dimension(maxsize2,4  ,nBLK,6) :: rbuf2
  real, dimension(maxsize3,2,2,nBLK,6) :: rbuf3

  ! Actual size of messages for the two VxB orientations and the edge
  integer :: isize1, isize2, isize3

  ! Restricted values to be sent are stored in these buffers
  real, dimension(:,:), allocatable :: sbuf1, sbuf2
  real, dimension(:),   allocatable :: sbuf3

  logical :: oktest, oktest_me
  integer :: inow, jnow, iError

  !---------------------------------------------------------------------------

  if(index(test_string,'NOCORRECT_VXB')>0)return

  call set_oktest('correct_vxb',oktest, oktest_me)
  if(oktest)write(*,*)'correct_VxB me=',iProc

  select case(optimize_message_pass)
  case('dir')
     ! Send messages for two faces together
     do isweep=1,3
        call VxB_pass_faces(2*isweep-1,2*isweep)
     end do
  case('face','min')
     ! Send messages face by face
     do isweep=1,6
        call VxB_pass_faces(isweep,isweep)
     end do
  case default
     ! Send messages for all faces
     call VxB_pass_faces(1,6)
  end select

  if(oktest_me)write(*,*)'VxB_pass finished'

contains

  subroutine VxB_pass_faces(ifacemin,ifacemax)

    integer, intent(in):: ifacemin,ifacemax
    !------------------------------------------------------------------------

    if(oktest)write(*,*)&
         'VxB_pass_faces:me,ifacemin,ifacemax=',iProc,ifacemin,ifacemax

    ! Debug
    if(okdebug)then
       rbuf1  =0.; rbuf2=0.; rbuf3=0.
    end if

    number_receive_requests = 0
    receive_requests = MPI_REQUEST_NULL

    do iface=ifacemin,ifacemax

       ! Set index ranges for the face
       call setranges

       if(okdebug.and.oktest)then
          write(*,*)'setranges for receive done'
          write(*,*)'me,iface,nxface,nyface,isize1,isize2,isize3',&
               iProc,iface,nxface,nyface,isize1,isize2,isize3
          write(*,*)'facedir,edgedir,sidedir=',facedir,edgedir,sidedir
       end if

       do iBLK = 1,nBlockMax
          if(unusedBLK(iBLK))CYCLE

          ! Post non-blocking receive for opposite face of neighbor block
          neiL=neiLEV(rface,iBLK)

          if(neiL==0)then
             ! Check for shared edges
             do iedge=1,2
                if(.not.recv_edge()) CYCLE
                do isubedge=1,2
                   neiP=neiedgeP(isubedge)
                   if(neiP==iProc) CYCLE

                   ! Remote receive egde
                   neiB=neiedgeB(isubedge)
                   itag = 100*neiB+10*iface+3

                   if(oktest)write(*,*)&
                        'Remote recv edge:me,iBLK,itag,neiP,neiB,iedge,isubedge',&
                        iProc,iBLK,itag,neiP,neiB,iedge,isubedge

                   call MPI_irecv(rbuf3(1,isubedge,iedge,iBLK,iface),isize3,&
                        MPI_REAL,neiP,itag,iComm,request,iError)

                   number_receive_requests = number_receive_requests + 1
                   receive_requests(number_receive_requests) = request
                end do
             end do
          end if

          ! Check if neighboring block is finer
          if(neiL/=-1)CYCLE

          do isubface=1,4
             neiP=neiPE(isubface,rface,iBLK)
             if(neiP==iProc)CYCLE

             ! Remote receive
             neiB=neiBLK(isubface,rface,iBLK)
             itag = 100*neiB+10*iface
             if(oktest.and.okdebug)write(*,*)&
                  'Remote recv,me,iBLK,itag,neiP,neiB,isubface=',&
                  iProc,iBLK,itag,neiP,neiB,isubface

             call MPI_irecv(rbuf1(1,isubface,iBLK,iface),isize1,&
                  MPI_REAL,neiP,itag+1,iComm,request,iError)

             number_receive_requests = number_receive_requests + 1
             receive_requests(number_receive_requests) = request

             call MPI_irecv(rbuf2(1,isubface,iBLK,iface),isize2,&
                  MPI_REAL,neiP,itag+2,iComm,request,iError)

             number_receive_requests = number_receive_requests + 1
             receive_requests(number_receive_requests) = request
          end do ! isubface
       end do ! iBLK
    end do ! iface

    !\
    ! Wait for all receive commands to be posted for all processors
    !/
    call barrier_mpi

    if(oktest)write(*,*)'receives posted: me=',iProc

    !\
    ! Send blocking messages with Rsend (ready to receive)
    !/
    do iface=ifacemin,ifacemax

       ! Set index ranges for the face
       call setranges

       if(okdebug.and.oktest)write(*,*)&
            'setranges for send done: me, iface=',iProc, iface

       allocate(sbuf1(nxface,nyface+1),sbuf2(nxface+1,nyface),sbuf3(isize3))

       if(okdebug.and.oktest)write(*,*)'allocation done, me,iface=',&
            iProc,iface

       do iBLK=1,nBlockMax
          if(unusedBLK(iBLK))CYCLE

          ! Check if neighbouring block is coarser
          neiL=neiLEV(iface,iBLK)
          if(neiL/=1)CYCLE

          if(okdebug.and.oktest)write(*,*)&
               'sending: me, iface,iBLK,neiL=',iProc,iface,iBLK,neiL

          ! Restrict VxB
          select case(facedir)
          case(1)
             sbuf1=0.5*(VxB_y(isface,1:nJ-1:2,1:nK+1:2,iBLK)+&
                        VxB_y(isface,2:nJ  :2,1:nK+1:2,iBLK))
             sbuf2=0.5*(VxB_z(isface,1:nJ+1:2,1:nK-1:2,iBLK)+&
                        VxB_z(isface,1:nJ+1:2,2:nK  :2,iBLK))
          case(2)
             sbuf1=0.5*(VxB_x(1:nI-1:2,isface,1:nK+1:2,iBLK)+&
                        VxB_x(2:nI  :2,isface,1:nK+1:2,iBLK))

             sbuf2=0.5*(VxB_z(1:nI+1:2,isface,1:nK-1:2,iBLK)+&
                        VxB_z(1:nI+1:2,isface,2:nK  :2,iBLK))
          case(3)
             sbuf1=0.5*(VxB_x(1:nI-1:2,1:nJ+1:2,isface,iBLK)+&
                        VxB_x(2:nI  :2,1:nJ+1:2,isface,iBLK))
             sbuf2=0.5*(VxB_y(1:nI+1:2,1:nJ-1:2,isface,iBLK)+&
                        VxB_y(1:nI+1:2,2:nJ  :2,isface,iBLK))

          end select

          ! Check if any of the edges require message passing
          if(send_edge())then

             select case(facedir)
             case(1)
                ! X face passes Y edge
                sbuf3=sbuf1(:,isedge)
             case(2)
                ! Y face passes Z edge
                sbuf3=sbuf2(isedge,:)
             case(3)
                ! Z face passes X edge
                sbuf3=sbuf1(:,isedge)
             end select

             if(oktest)write(*,*)&
                  'Send edge: me,iBLK,iface,isedge,sbuf3',&
                  iProc,iBLK,iface,isedge,sbuf3

             if(neiP==iProc)then
                ! Local copy into appropriate subedge
                ichild=BLKneighborCHILD(0,0,0,1,iBLK)
                isubedge=child2subedge(ichild,iface)

                if(oktest)write(*,*)&
                     'Local copy edge,me,iBLK,neiB,iedge,isubedge=',&
                     iProc,iBLK,neiB,iedge,isubedge
                call buf2subedge(sbuf3,isize3,neiB)
             else
                ! Remote send of edge
                itag = 100*iBLK+10*iface+3

                if(oktest)write(*,*)&
                     'Remote send edge,me,iBLK,itag,neiP,neiB,isedge=',&
                     iProc,iBLK,itag,neiP,neiB,isedge

                call MPI_Rsend(sbuf3,isize3,&
                     MPI_REAL,neiP,itag,iComm,iError)
             end if
          end if

          neiP=neiPE(1,iface,iBLK)
          neiB=neiBLK(1,iface,iBLK)

          if(neiP==iProc)then
             ! Local copy into appropriate subface
             ! Subface index =1,2,3, or 4 with respect to the coarse neighbor
             ichild=BLKneighborCHILD(0,0,0,1,iBLK)
             isubface=child2subface(ichild,iface)
             call buf2subface(sbuf1,sbuf2,nxface,nyface,neiB)
          else
             ! Remote send of face
             itag = 100*iBLK+10*iface

             if(oktest.and.okdebug)write(*,*)&
                  'Remote send,me,iBLK,itag,neiP,neiB=',&
                  iProc,iBLK,itag,neiP,neiB

             call MPI_Rsend(sbuf1,isize1,&
                  MPI_REAL,neiP,itag+1,iComm,iError)
             call MPI_Rsend(sbuf2,isize2,&
                  MPI_REAL,neiP,itag+2,iComm,iError)
          end if
       end do ! iBLK

       deallocate(sbuf1,sbuf2,sbuf3)

       if(oktest_me)write(*,*)'messages sent, me, iface=',iProc,iface

    end do ! iface

    !\
    ! WAIT FOR ALL MESSAGES TO BE RECEIVED
    !/
    if (number_receive_requests > 0) &
         call MPI_waitall(number_receive_requests,receive_requests,status,iError)

    if(oktest_me)write(*,*)'messages received, me, facedir=',iProc, facedir

    ! Copy averaged VxB received from non-local finer neigbors
    ! and stored in the buffers into the coarse VxB

    do iface=ifacemin,ifacemax

       ! Set index ranges for the face
       call setranges

       do iBLK = 1,nBlockMax
          if(unusedBLK(iBLK))CYCLE

          neiL=neiLEV(rface,iBLK)
          if(neiL==0)then
             ! Check if remote edges were received
             do iedge=1,2
                if(.not.recv_edge())CYCLE

                if(oktest)write(*,*)'receive edge, me,iBLK,iface,iedge=',&
                     iProc,iBLK,iface,iedge

                do isubedge=1,2
                   if(neiedgeP(isubedge)==iProc)CYCLE

                   neiB=neiedgeB(isubedge)

                   if(oktest)write(*,*)'read buffer, me,isize3,neiB,rbuf3=',&
                        iProc,isize3,neiB,&
                        rbuf3(1:isize3,isubedge,iedge,neiB,iface)

                   call buf2subedge(rbuf3(1,isubedge,iedge,iBLK,iface),&
                        isize3,iBLK)
                end do
             enddo
          endif

          ! Check if neighboring block is finer
          if(neiL/=-1)CYCLE

          do isubface=1,4
             if(neiPE(isubface,rface,iBLK)==iProc) CYCLE

             neiB=neiBLK(isubface,rface,iBLK)
             if(okdebug.and.oktest)&
                  write(*,*)'buf2subface: me, isubface, iBLK, neiB=',&
                  iProc,isubface,iBLK,neiB

             call buf2subface(rbuf1(1,isubface,iBLK,iface),&
                  rbuf2(1,isubface,iBLK,iface),&
                  nxface,nyface,iBLK)
          end do
       end do ! iBLK
    end do ! iface

    if(oktest)write(*,*)'VxB_pass_faces finished: me, ifacemin, ifacemax=',&
         iProc,ifacemin,ifacemax

  end subroutine VxB_pass_faces

  !===========================================================================

  subroutine setranges

    ! Calculate directions for face, edge, and side direction for edge
    facedir=(iface+1)/2; edgedir=mod(facedir,3)+1; sidedir=mod(edgedir,3)+1

    ! Calculate the size of the subfaces for the 2 orientations of VxB
    ! and for the subedge associated with the face
    select case(facedir)
    case(1)
       nxface=nJ/2; nyface=nK/2
    case(2)
       nxface=nI/2; nyface=nK/2
    case(3)
       nxface=nI/2; nyface=nJ/2
    end select

    isize1=nxface*(nyface+1); isize2=(nxface+1)*nyface
    isize3=nIJK_D(edgedir)/2

    select case(iface)
    case(1,3,5)
       rface=iface+1; isface=1; irface=nIJK_D(facedir)+1
    case(2,4,6)
       rface=iface-1; isface=nIJK_D(facedir)+1; irface=1
    end select

  end subroutine setranges

  !===========================================================================
  subroutine buf2subface(qbuf1,qbuf2,qnxface,qnyface,qBLK)

    integer, intent(in) :: qnxface,qnyface,qBLK
    real, intent(inout) :: qbuf1(qnxface,qnyface+1),qbuf2(qnxface+1,qnyface)
    !-------------------------------------------------------------------------

    ! Assign VxB on a subface of receiving face

    select case(facedir)
    case(1)
       select case(isubface)
          ! Beware, case(2) and case(3) are swapped
       case(1)
          VxB_y(irface, 1:nJ/2  , 1:nK/2+1, qBLK)=qbuf1
          VxB_z(irface, 1:nJ/2+1, 1:nK/2  , qBLK)=qbuf2
       case(3)
          VxB_y(irface, nJ/2+1:nJ  , 1:nK/2+1, qBLK)=qbuf1
          VxB_z(irface, nJ/2+1:nJ+1, 1:nK/2  , qBLK)=qbuf2
       case(2)
          VxB_y(irface, 1:nJ/2   ,nK/2+1:nK+1, qBLK)=qbuf1
          VxB_z(irface, 1:nJ/2+1 ,nK/2+1:nK  , qBLK)=qbuf2
       case(4)
          VxB_y(irface,nJ/2+1:nJ  ,nK/2+1:nK+1,qBLK)=qbuf1
          VxB_z(irface,nJ/2+1:nJ+1,nK/2+1:nK  ,qBLK)=qbuf2
       end select
    case(2)
       select case(isubface)
          ! Beware, case(2) and case(3) are swapped
       case(1)
          VxB_x(1:nI/2  , irface, 1:nK/2+1, qBLK)=qbuf1
          VxB_z(1:nI/2+1, irface, 1:nK/2  , qBLK)=qbuf2
       case(3)
          VxB_x(nI/2+1:nI   ,irface, 1:nK/2+1, qBLK)=qbuf1
          VxB_z(nI/2+1:nI+1 ,irface, 1:nK/2  , qBLK)=qbuf2
       case(2)
          VxB_x(1:nI/2  , irface, nK/2+1:nK+1, qBLK)=qbuf1
          VxB_z(1:nI/2+1, irface, nK/2+1:nK  , qBLK)=qbuf2
       case(4)
          VxB_x(nI/2+1:nI  ,irface,nK/2+1:nK+1,qBLK)=qbuf1
          VxB_z(nI/2+1:nI+1,irface,nK/2+1:nK  ,qBLK)=qbuf2
       end select
    case(3)

       select case(isubface)
          ! Beware, case(2) and case(3) are not swapped
       case(1)
          VxB_x(1:nI/2  , 1:nJ/2+1, irface,qBLK)=qbuf1
          VxB_y(1:nI/2+1, 1:nJ/2  , irface,qBLK)=qbuf2
       case(2)
          VxB_x(nI/2+1:nI  , 1:nJ/2+1 ,irface, qBLK)=qbuf1
          VxB_y(nI/2+1:nI+1, 1:nJ/2   ,irface, qBLK)=qbuf2
       case(3)
          VxB_x(1:nI/2  , nJ/2+1:nJ+1, irface, qBLK)=qbuf1
          VxB_y(1:nI/2+1, nJ/2+1:nJ  , irface, qBLK)=qbuf2
       case(4)
          VxB_x(nI/2+1:nI  ,nJ/2+1:nJ+1,irface,qBLK)=qbuf1
          VxB_y(nI/2+1:nI+1,nJ/2+1:nJ  ,irface,qBLK)=qbuf2
       end select

    end select

  end subroutine buf2subface

  !===========================================================================
  subroutine buf2subedge(qbuf,qsize,qBLK)

    integer, intent(in) :: qsize,qBLK
    real, intent(inout) :: qbuf(qsize)
    !-------------------------------------------------------------------------

    ! Assign VxB on a subface of receiving subedge indexed by isubedge

    select case(facedir)
    case(1)
       if(isubedge==1)then
          VxB_y(irface,1:nJ/2        ,iredge,qBLK)=qbuf
       else
          VxB_y(irface,nJ/2+1:nJ,iredge,qBLK)=qbuf
       endif
    case(2)
       if(isubedge==1)then
          VxB_z(iredge,irface,1:nK/2        ,qBLK)=qbuf
       else
          VxB_z(iredge,irface,nK/2+1:nK,qBLK)=qbuf
       end if
    case(3)
       if(isubedge==1)then
          VxB_x(1:nI/2        ,iredge,irface,qBLK)=qbuf
       else
          VxB_x(nI/2+1:nI,iredge,irface,qBLK)=qbuf
       end if
    end select

  end subroutine buf2subedge

  !===========================================================================

  logical function send_edge()

    ! Check if any of the edges associated with the sending face are shared
    ! with a coarser block to which an edge message should be sent.
    ! We already know that the neighbor in the face direction is coarser,
    ! so the diagonal and the extra direction are checked. On the figures
    ! below the left case returns send_edge=.true., while the other cases
    ! return false. Here iface 3 (south, on the figure ==) of block iBLK (i).
    ! is checked for the east Z edge associated with the face (Z on the figure)
    ! The edge is orthogonal to the screen. Below each figure the caption
    ! indicates how the Z edge will get corrected, if it's needed at all.
    !
    ! +-----+                                 +-----+      
    ! |     |                                 |     |
    ! |     |                                 |     |
    ! |     +--+             +--+--+          |     +--+          +--+--+
    ! |     | i|             |  | i|          |     | i|          |  | i|
    ! |     |  |             |  |  |          |     |  |          |  |  |
    ! +-----Z==+--+       +---##Z==+--+       +--+--Z==+--+       +--Z==+--+
    ! |     |     |       |     |     |          |  |     |       |  |     |
    ! |     |     |       |     |     |          |  |     |       |  |     |
    ! |     |     |       |     |     |          +--+     |       +--+     |
    ! |     |     |       |     |     |             |     |          |     |
    ! |     |     |       |     |     |             |     |          |     |
    ! +-----+-----+       +-----+-----+             +-----+          +-----+
    !
    !   EDGE PASS           FACE PASS        SAME LEVEL          SAME LEVEL
    !
    ! When true is returned, all the necessary variables are set

    integer :: d_edge(3), d_side(3), q_edge
    !-------------------------------------------------------------------------

    if(index(test_string,'NOEDGEPASS')>0)then
       send_edge=.false. ; return
    end if

    ! Calculate edge and side face shifts
    d_edge=0; d_edge(facedir)=iface-rface
    d_side=0; 

    ! Check both edges
    do q_edge=-1,1,2
       d_edge(sidedir)=q_edge; d_side(sidedir)=q_edge

       if(BLKneighborLEV(d_edge(1),d_edge(2),d_edge(3),iBLK)==1 .and. &
          BLKneighborLEV(d_side(1),d_side(2),d_side(3),iBLK)==1)then

          if(q_edge==-1)then
             isedge=1;                   iredge=nIJK_D(sidedir)+1; iedge=2
          else
             isedge=nIJK_D(sidedir)/2+1; iredge=1;                 iedge=1
          end if

          neiB=BLKneighborBLK(d_edge(1),d_edge(2),d_edge(3),1,iBLK)
          neiP= BLKneighborPE(d_edge(1),d_edge(2),d_edge(3),1,iBLK)

          ! If one edge is to be sent, the other is definitely not
          send_edge=.true.
          return
       end if
    end do

    send_edge=.false.

  end function send_edge

  !===========================================================================

  logical function recv_edge()

    ! Check if the edge indexed with iedge associated with the recv face 
    ! is shared with a finer block from which an edge should be received.
    ! We already know that the neighbor in the face direction is equal,
    ! so the diagonal and the extra direction are checked. On the figures
    ! below the left case returns recv_edge=.true., while the other cases
    ! return false. Here iface=4 (north), and Z edge associated with Y-faces
    ! is orthogonal to the screen (see the Z character). Below each figure
    ! the caption shows how the Z edge is corrected, if it's needed at all.
    !
    ! +-----+             +-----+         +-----+-----+    +-----+-----+      
    ! |     |             |     |         |     |     |    |     |     |
    ! |     |             |     |         |     |     |    |     |     |
    ! |     +--+          |     +--+      |     |     |    |     |     |
    ! |     |  |          |     |  |      |     |     |    |     |     |
    ! |     |  |          |     |  |      |     |     |    |     |     |
    ! +=====Z--+--+       +=====Z--+      +=====Z-----+    +=====Z--+--+
    ! |     |     |       |     #  |      |     |     |    |     |  |
    ! |     |     |       |     #  |      |     |     |    |     |  |
    ! |  i  |     |       |  i  |--+      |  i  |     |    |  i  |--+
    ! |     |     |       |     |         |     |     |    |     |
    ! |     |     |       |     |         |     |     |    |     |
    ! +-----+-----+       +-----+         +-----+-----+    +-----+
    !
    !   EDGE PASS          FACE PASS        SAME LEVEL        SAME LEVEL
    !
    ! When true is returned, all the necessary variables are set

    integer :: d_edge(3), d_side(3)
    !-------------------------------------------------------------------------

    if(index(test_string,'NOEDGEPASS')>0)then
       recv_edge=.false. ; return
    end if

    ! Calculate shifts for edge and side face for edge iedge on face rface
    d_edge=0; d_edge(sidedir)=2*iedge-3; d_edge(facedir)=rface-iface
    d_side=0; d_side(sidedir)=2*iedge-3;

    if(BLKneighborLEV(d_edge(1),d_edge(2),d_edge(3),iBLK)==-1 .and. &
       BLKneighborLEV(d_side(1),d_side(2),d_side(3),iBLK)== 0)then

       if(iedge==1)then
          iredge=1
       else
          iredge=nIJK_D(sidedir)+1
       end if
       neiedgeP= BLKneighborPE(d_edge(1),d_edge(2),d_edge(3),1:2,iBLK)
       neiedgeB=BLKneighborBLK(d_edge(1),d_edge(2),d_edge(3),1:2,iBLK)

       recv_edge= .true.
    else
       recv_edge=.false.
    end if

  end function recv_edge

end subroutine correct_VxB

!==============================================================================
subroutine b_face_fine_pass

  ! Set B*FaceFine_*SB from finer face
  use ModProcMH
  use ModMain, ONLY : nBLock,unusedBLK,BLKtest
  use ModCT
  use ModAMR, ONLY : child2subface
  use ModParallel, ONLY : neiLEV,neiBLK,neiPE,BLKneighborCHILD
  use ModMpi
  implicit none

  integer :: iError
  integer :: iBlock, iTag, iProcNei, iBlockNei
  integer :: iSide, iFace, iFaceOther, iSubface, iChild, iSize

  integer :: number_receive_requests, request
  integer :: receive_requests(nBLK*24)
  integer :: status(MPI_STATUS_SIZE, nBLK*24)

  real, allocatable :: Buffer(:,:)

  logical :: oktest, oktest_me
  !------------------------------------------------------------------
  call set_oktest('b_face_fine_pass',oktest,oktest_me)

  ! Initialize counters for non-blocking receive requests
  number_receive_requests = 0
  receive_requests = MPI_REQUEST_NULL

  !\
  ! Non-blocking recieve messages from fine blocks
  ! or copy for local blocks
  !/
  do iBlock=1,nBlock
     if(unusedBLK(iBlock)) CYCLE
!!!     if(.not.refine_list(iBlock,iProc)) CYCLE

     do iFace=east_,top_
        if(neiLEV(iFace,iBlock)==-1)then
           do iSubFace=1,4
              iProcNei =neiPE(iSubFace,iFace,iBlock)
              if(iProcNei /= iProc)then
                 call recv_b_face_fine
              else
                 call copy_b_face_fine
              end if
           end do
        end if
     end do
  end do

  !\
  ! Wait for all receive commands to be posted for all processors
  !/
  call barrier_mpi

  !\
  ! Send blocking messages with Rsend (ready to receive)
  !/
  do iBlock=1,nBlock
     if(unusedBLK(iBlock)) CYCLE

     do iFace=east_,top_
        if(neiLEV(iFace,iBlock)/=1) CYCLE
        iBlockNei=neiBLK(1,iFace,iBlock)
        iProcNei =neiPE(1,iFace,iBlock)
!!!        if(.not.refine_list(iBlockNei,iProcNei)) CYCLE
        if(iProcNei==iProc) CYCLE ! local copy
        call send_b_face_fine
     end do

  end do

  !\
  ! WAIT FOR ALL MESSAGES TO BE RECEIVED
  !/
  if (number_receive_requests > 0) &
       call MPI_waitall(number_receive_requests,receive_requests,status,iError)

contains

  subroutine recv_b_face_fine

    ! write(*,*)'recv_b_face_fine: me,iBlock,iFace=',iProc,iBlock,iFace

    iTag=100*iBlock+10*iFace+iSubFace
    select case(iFace)
    case(east_,west_)
       iSize=nJ*nK
       iSide=iFace-east_+1
       call MPI_irecv(BxFaceFine_XQSB(1,1,iSubFace,iSide,iBlock), iSize, &
            MPI_REAL, iProcNei, iTag, iComm, request, iError)
    case(south_,north_)
       iSize=nI*nK
       iSide=iFace-south_+1
       call MPI_irecv(ByFaceFine_YQSB(1,1,iSubFace,iSide,iBlock), iSize, &
            MPI_REAL, iProcNei, iTag, iComm, request, iError)
    case(bot_,top_)
       iSize=nI*nJ
       iSide=iFace-bot_+1
       call MPI_irecv(BzFaceFine_ZQSB(1,1,iSubFace,iSide,iBlock), iSize, &
            MPI_REAL, iProcNei, iTag, iComm, request, iError)
    end select
    number_receive_requests = number_receive_requests + 1
    receive_requests(number_receive_requests) = request

  end subroutine recv_b_face_fine

  subroutine send_b_face_fine

    ! write(*,*)'send_b_face_fine: me,iBlock,iFace=',iProc,iBlock,iFace

    select case(iFace)
    case(east_)
       iSize=nJ*nK
       allocate(Buffer(nJ,nK))
       iFaceOther=west_
       Buffer=BxFace_BLK(1,1:nJ,1:nK,iBlock)
    case(west_)
       iSize=nJ*nK
       allocate(Buffer(nJ,nK))
       iFaceOther=east_
       Buffer=BxFace_BLK(nI+1,1:nJ,1:nK,iBlock)
    case(south_)
       iSize=nI*nK
       allocate(Buffer(nI,nK))
       iFaceOther=north_
       Buffer=ByFace_BLK(1:nI,1,1:nK,iBlock)
    case(north_)
       iSize=nI*nK
       allocate(Buffer(nI,nK))
       iFaceOther=south_
       Buffer=ByFace_BLK(1:nI,nJ+1,1:nK,iBlock)
    case(bot_)
       iSize=nI*nJ
       allocate(Buffer(nI,nJ))
       iFaceOther=top_
       Buffer=BzFace_BLK(1:nI,1:nJ,1,iBlock)
    case(top_)
       iSize=nI*nJ
       allocate(Buffer(nI,nJ))
       iFaceOther=bot_
       Buffer=BzFace_BLK(1:nI,1:nJ,nK+1,iBlock)
    end select

    iChild=BLKneighborCHILD(0,0,0,1,iBlock)
    iSubface=child2subface(iChild,iFace)
    iTag=100*iBlockNei+10*iFaceOther+iSubFace
    call MPI_Rsend(buffer, iSize, &
         MPI_REAL, iProcNei, iTag, iComm, iError)

    deallocate(Buffer)

  end subroutine send_b_face_fine

  subroutine copy_b_face_fine

    ! Copy fine normal B face component from the neighboring block

    iBlockNei=neiBLK(iSubFace,iFace,iBlock)

    if(oktest_me.and.(iBlock==BLKtest.or.iBlockNei==BLKtest))&
         write(*,*)'copy from iBlockNei=',iBlockNei,' to iBlock=',iBlock

    select case(iFace)
    case(east_)
       BxFaceFine_XQSB(:,:,iSubFace,1,iBlock)=&
            BxFace_BLK(nI+1,1:nJ,1:nK,iBlockNei)
    case(west_)
       BxFaceFine_XQSB(:,:,iSubFace,2,iBlock)=&
            BxFace_BLK(   1,1:nJ,1:nK,iBlockNei)
    case(south_)
       ByFaceFine_YQSB(:,:,iSubFace,1,iBlock)=&
            ByFace_BLK(1:nI,nJ+1,1:nK,iBlockNei)
    case(north_)
       ByFaceFine_YQSB(:,:,iSubFace,2,iBlock)=&
            ByFace_BLK(1:nI,   1,1:nK,iBlockNei)
    case(bot_)
       BzFaceFine_ZQSB(:,:,iSubFace,1,iBlock)=&
            BzFace_BLK(1:nI,1:nJ,nK+1,iBlockNei)
    case(top_)
       BzFaceFine_ZQSB(:,:,iSubFace,2,iBlock)=&
            BzFace_BLK(1:nI,1:nJ,   1,iBlockNei)
    end select

  end subroutine copy_b_face_fine

end subroutine b_face_fine_pass
