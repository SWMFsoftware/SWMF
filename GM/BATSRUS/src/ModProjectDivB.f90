!^CFG COPYRIGHT UM
!^CFG FILE PROJECTION

!=============================================================================
subroutine project_B

  ! Project B according to B'=B-grad phi, where Laplace phi=div B
  ! Solve the Poisson equation with an iterative scheme
  !
  ! See G. Toth, 2000, Journal of Computational Physics, 161, 605-652

  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,nBLK,Itest,Jtest,Ktest,BLKtest
  use ModVarIndexes,ONLY: Bx_,Bz_,P_
  use ModAdvance, ONLY : State_VGB
  use ModGeometry, ONLY : true_cell, x_BLK, y_BLK, z_BLK
  use ModProject
  use ModMain, ONLY : UseConstrainB                   !^CFG IF CONSTRAINB
  use ModCT, ONLY : Bxface_BLK,Byface_BLK,Bzface_BLK  !^CFG IF CONSTRAINB
  implicit none

  ! Local variables
  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK) :: &
       proj_divb, &                  ! original error to be projected out
       phi                           ! Scalar field in the Poisson problem


  integer :: info      ! error status of the iterative scheme
  integer :: nmatvec   ! number of matvex operations performed
  real ::    resid     ! residual after iterative solver stopped

  real :: divbmax_now  ! current value of max(|div B|)
  real :: pmin_old, pmin_new     ! original and projected value of min(p)


  integer :: loc(5)
  logical :: oktest, oktest_me

  ! Functions:

  real, external :: maxval_abs_BLK, minval_BLK, maxval_loc_abs_BLK, &
       minval_loc_BLK

  !---------------------------------------------------------------------------

  call set_oktest('project_B',oktest,oktest_me)

  call timing_start('project_B')

  if(oktest_me)then
     write(*,*)'Project_B, old B:', &
          State_VGB(Bx_:Bz_,Itest,Jtest,Ktest,BLKtest),  &
          true_cell(Itest,Jtest,Ktest,BLKtest)
     if(UseConstrainB)then                      !^CFG IF CONSTRAINB BEGIN
        write(*,*)'Project_B, Bxface,true_cell:',    &
             Bxface_BLK(Itest,Jtest,Ktest,BLKtest),  &
             Bxface_BLK(Itest+1,Jtest,Ktest,BLKtest),&
             true_cell(Itest-1,Jtest,Ktest,BLKtest), &
             true_cell(Itest+1,Jtest,Ktest,BLKtest)

        write(*,*)'Project_B, Byface,true_cell:',    &
             Byface_BLK(Itest,Jtest,Ktest,BLKtest),  &
             Byface_BLK(Itest,Jtest+1,Ktest,BLKtest),&
             true_cell(Itest,Jtest-1,Ktest,BLKtest), &
             true_cell(Itest,Jtest+1,Ktest,BLKtest)

        write(*,*)'Project_B, Bzface,true_cell:',    &
             Bzface_BLK(Itest,Jtest,Ktest,BLKtest),  &
             Bzface_BLK(Itest,Jtest,Ktest+1,BLKtest),&
             true_cell(Itest,Jtest,Ktest-1,BLKtest), &
             true_cell(Itest,Jtest,Ktest+1,BLKtest)

     end if                                     !^CFG END CONSTRAINB
  end if

  ! Determine the div B error
  call proj_get_divb(proj_divB)

  if(oktest_me)write(*,*)'proj_divb=',proj_divb(Itest,Jtest,Ktest,BLKtest)

  !DEBUG
  !call proj_bound(proj_divb)
  !call show_BLK('divb',proj_divb)
  !DEBUG
  !call stop_mpi('debug')

  if(oktest)then
     pmin_old=minval_BLK(nProc,State_VGB(P_,:,:,:,:))
     divbmax_now=maxval_abs_BLK(nProc,proj_divB)
     if(oktest_me)write(*,*)&
          'Project_B: old max(abs(divB)), min(p)=',divbmax_now,pmin_old
  endif

  ! Calculate and save divbmax, the maximum div B allowed
  if(divbmax<divbmin)then
     ! Default values if projection parameters were not set
     if(proj_divbcoeff<0.0)  proj_divbcoeff=0.1
     if(proj_divbconst<0.0)  proj_divbconst=0.0
     if(proj_matvecmax<0)proj_matvecmax=1000
     select case(proj_typestop)
     case('rel')
        divbmax=proj_divbcoeff
     case('max')
        if(proj_divbcoeff>0.0)then
           divbmax=proj_divbcoeff*maxval_abs_BLK(nProc,proj_divb)
        else
           divbmax=proj_divbconst
        endif
     end select
     if(oktest_me)write(*,*)'Allowed maximum for div B:',divbmax
     if(divbmax<divbmin)&
          call stop_mpi('Error in project_B: Too small value for divbmax')
  endif

  ! Solve the Poisson equation Laplace phi = div B
  call proj_poisson(proj_divb,divbmax,proj_typestop,proj_matvecmax, &
       info,nmatvec,resid,phi)

  if(oktest_me)write(*,*)'Project_B: Poisson solver info, nmatvec, resid:',&
       info, nmatvec, resid

  ! Do not do anything if the initial guess satisfied the stopping criterion
  if(info==3)goto 100

  ! Do not subtract grad(phi) from B if the iterations did not reduce error
  if(info<0)goto 100

  ! Get the ghost cell values for the solution phi
  call proj_boundphi(phi,1,3);

  !if(oktest_me)write(*,*)'proj_divb,lapl(phi)=',&
  !     proj_divb(Itest,Jtest,Ktest,BLKtest),&
  !     (phi(Itest+2,Jtest,Ktest,BLKtest)+phi(Itest-2,Jtest,Ktest,BLKtest) &
  !     +phi(Itest,Jtest+2,Ktest,BLKtest)+phi(Itest,Jtest-2,Ktest,BLKtest) &
  !     +phi(Itest,Jtest,Ktest-2,BLKtest)+phi(Itest,Jtest,Ktest-2,BLKtest) &
  !     -6*phi(Itest,Jtest,Ktest,BLKtest))/dx_BLK(BLKtest)**2/4

  ! Correct B field: B'=B-dB where dB=grad(phi) 
  call proj_correction(phi)




  ! Recalculate boundary conditions
  call exchange_messages

!!! test_string=test_string(1:50)

  ! Testing div B'=0
  if(oktest)then
     call proj_get_divB(proj_divB)

     if(oktest_me)write(*,*)'after proj_divb=',&
          proj_divb(Itest,Jtest,Ktest,BLKtest)

     divbmax_now=maxval_loc_abs_BLK(nProc,proj_divB,loc)
     if(abs(resid-divbmax_now)/(divbmax_now+1.e-12) > 0.1)then
        if(iProc==loc(5))write(*,*)&
             'Project_B: resid,max(abs(divB)),loc,X,Y,Z=',&
             resid, divbmax_now, loc, &
             x_BLK(loc(1),loc(2),loc(3),loc(4)),&
             y_BLK(loc(1),loc(2),loc(3),loc(4)),&
             z_BLK(loc(1),loc(2),loc(3),loc(4))
     else
        if(oktest_me)write(*,*)'Project_B: resid,max(abs(divB))',&
             resid, divbmax_now
     end if
     pmin_new   =minval_loc_BLK(nProc,State_VGB(P_,:,:,:,:),loc)
     if(pmin_new<0.5*pmin_old)then
        if(iProc==loc(5))write(*,*)'Project_B: min(p),loc,X,Y,Z=',&
             pmin_new,loc,x_BLK(loc(1),loc(2),loc(3),loc(4)),&
             y_BLK(loc(1),loc(2),loc(3),loc(4)),&
             z_BLK(loc(1),loc(2),loc(3),loc(4))
     else
        if(oktest_me)write(*,*)'Project_B: new min(p)',pmin_new
     endif

     if(oktest_me)then
        write(*,*)'Project_B, new B, true_cell:', &
             State_VGB(Bx_:Bz_,Itest,Jtest,Ktest,BLKtest),  &
             true_cell(Itest,Jtest,Ktest,BLKtest)
        if(UseConstrainB)then                      !^CFG IF CONSTRAINB BEGIN
           write(*,*)'Project_B, Bxface,true_cell:',    &
                Bxface_BLK(Itest,Jtest,Ktest,BLKtest),  &
                Bxface_BLK(Itest+1,Jtest,Ktest,BLKtest),&
                true_cell(Itest-1,Jtest,Ktest,BLKtest), &
                true_cell(Itest+1,Jtest,Ktest,BLKtest)

           write(*,*)'Project_B, Byface,true_cell:',    &
                Byface_BLK(Itest,Jtest,Ktest,BLKtest),  &
                Byface_BLK(Itest,Jtest+1,Ktest,BLKtest),&
                true_cell(Itest,Jtest-1,Ktest,BLKtest), &
                true_cell(Itest,Jtest+1,Ktest,BLKtest)

           write(*,*)'Project_B, Bzface,true_cell:',    &
                Bzface_BLK(Itest,Jtest,Ktest,BLKtest),  &
                Bzface_BLK(Itest,Jtest,Ktest+1,BLKtest),&
                true_cell(Itest,Jtest,Ktest-1,BLKtest), &
                true_cell(Itest,Jtest,Ktest+1,BLKtest)

        end if                                     !^CFG END CONSTRAINB
     end if
  endif

  !if(oktest_me)write(*,*)'final proj_divb=',&
  !        proj_divb(Itest,Jtest,Ktest,BLKtest)

100 call timing_stop('project_B')

end subroutine project_B

!=============================================================================
subroutine proj_get_divB(proj_divB)
  ! Calculate div B using simple finite differences
  ! Do corrections for mesh refinement
  use ModMain, ONLY : nI,nJ,nK,nBLK,nBlock,unusedBLK
  use ModVarIndexes, ONLY : Bx_,By_,Bz_
  use ModAdvance, ONLY : State_VGB
  use ModGeometry, ONLY : dx_BLK, dy_BLK, dz_BLK
  use ModProject
  use ModMain, ONLY : UseConstrainB                       !^CFG IF CONSTRAINB
  use ModCT, ONLY : Bxface_BLK,Byface_BLK,Bzface_BLK      !^CFG IF CONSTRAINB
  implicit none

  ! Argument

  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), intent(out) :: proj_divb

  ! Local variables
  integer :: iBLK, i, j, k
  real    :: DxInvHalf, DyInvHalf, DzInvHalf
  !---------------------------------------------------------------------------

  proj_divB=0.0

  if(UseConstrainB)then                      !^CFG IF CONSTRAINB BEGIN
     do iBLK=1,nBlock
        if(unusedBLK(iBLK)) CYCLE

        proj_divB(1:nI,1:nJ,1:nK,iBLK)= &
           (Bxface_BLK(2:nI+1,1:nJ  ,1:nK  ,iBLK)              &
           -Bxface_BLK(1:nI  ,1:nJ  ,1:nK  ,iBLK))/dx_BLK(iBLK)&
          +(Byface_BLK(1:nI  ,2:nJ+1,1:nK  ,iBLK)              &
           -Byface_BLK(1:nI  ,1:nJ  ,1:nK  ,iBLK))/dy_BLK(iBLK)&
          +(Bzface_BLK(1:nI  ,1:nJ  ,2:nK+1,iBLK)              &
           -Bzface_BLK(1:nI  ,1:nJ  ,1:nK  ,iBLK))/dz_BLK(iBLK)
     end do
  else                                       !^CFG END CONSTRAINB
     do iBLK=1,nBlock
        if(unusedBLK(iBLK))CYCLE
        DxInvHalf = 0.5/dx_BLK(iBLK);
        DyInvHalf = 0.5/dy_BLK(iBLK);
        DzInvHalf = 0.5/dz_BLK(iBLK);
        do k=1,nK; do j=1,nJ; do i=1,nI
           proj_divb(i,j,k,iBLK) = &
                DxInvHalf* &
                (State_VGB(Bx_,i+1,j,k,iBLK)-State_VGB(Bx_,i-1,j,k,iBLK)) + &
                DyInvHalf* &
                (State_VGB(By_,i,j+1,k,iBLK)-State_VGB(By_,i,j-1,k,iBLK)) + &
                DzInvHalf* &
                (State_VGB(Bz_,i,j,k+1,iBLK)-State_VGB(Bz_,i,j,k-1,iBLK))
       end do; end do; end do
    end do
  end if                                     !^CFG IF CONSTRAINB

end subroutine proj_get_divB

!=============================================================================
! This is just an interface to the various iterative schemes to
! Solve grad div phi = rhs 
subroutine proj_poisson(rhs,tolerance,typestop,matvecmax,&
     info,nmatvec,resid,phi)
  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,nBLK
  use ModProject
  implicit none

  ! Arguments

  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), intent(inout):: rhs
  !    on input: the right hand side of the Poisson equation
  !    on output: the residual

  character (len=3), intent(in) :: typestop
  real, intent(in) :: tolerance
  ! The required accuracy of the solution is given by tolerance and typestop:
  !    typestop='max': resid=max(abs(Laplace(phi)-rhs)) < tolerance
  !    typestop='abs': resid=sum((Laplace(phi)-rhs)**2) < tolerance
  !    typestop='rel': resid=sum((Laplace(phi)-rhs)**2) < tolerance*resid_0
  !       where resid_0 is the initial residuum, which is simply sum(rhs**2)

  integer, intent(in) :: matvecmax
  ! At most matvecmax matrix vector multiplications can be performed

  integer, intent(out) :: info
  ! The success of the iteration is given by returned value of info:
  !     abs(info)=  0 - solution found satisfying given tolerance.
  !                 1 - iteration aborted due to division by very small value.
  !                 2 - no convergence within maximum number of iterations.
  !                 3 - initial guess satisfies the stopping criterion.
  !    sign(info)=  + - residual decreased
  !                 - - residual did not reduce

  integer, intent(out) :: nmatvec
  ! The total number of matvec operations done during the iterations

  real, intent(out)    :: resid
  ! The residual after the iterations

  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), intent(out):: phi
  !    The solution

  ! Local variables

  integer :: ierror=-1

  logical :: oktest, oktest_me
  !---------------------------------------------------------------------------

  call set_oktest('proj_poisson',oktest,oktest_me)
  if(oktest_me)write(*,*)'Proj_Poisson solver called with ',proj_method

  ! Initialize parameters and phi for the iterative solvers
  nmatvec =matvecmax
  resid=tolerance

  ! Solve the Poisson equation
  select case(proj_method)
  case('cg')
     call       proj_cg(rhs,phi,nmatvec,resid,typestop,info)

  case('bicgstab')
     call proj_bicgstab(rhs,phi,nmatvec,resid,typestop,info)
  case default
     call stop_mpi('Error in Proj_Poisson: Unknown type of iterative method')
  end select

  if(oktest_me)write(*,*)'Poisson info,nmatvec,resid',info,nmatvec,resid

  if(info/=0.and.info/=3)then
     if(iProc==0)then
        if(ierror<1)then
           ! Print additional information
           write(*,"(a,i2,a,i5)")'info=',info,' nmatvec=',nmatvec
           select case(abs(info))
           case(1)
              write(*,*)'Breakdown due to division by a small value'
           case(2)
              write(*,*)'No convergence within maximum number of iterations'
           end select
           if(info>0)write(*,*)'The residual decreased'
           if(info<0)write(*,*)'The residual did not decrease'
        end if
        call error_report('Poisson solver failed, resid',resid,ierror,.true.)
     end if
  end if

end subroutine proj_poisson

!=============================================================================
! Calculate Laplace phi
subroutine proj_matvec(phi,laplace_phi)
  use ModMain, ONLY : nBLK,nBlock,unusedBLK,nI,nJ,nK
  use ModGeometry, ONLY : dx_BLK,dy_BLK,dz_BLK,true_cell,body_BLK
  use ModProject
  use ModMain, ONLY : UseConstrainB          !^CFG IF CONSTRAINB
  use ModCT                                  !^CFG IF CONSTRAINB
  implicit none

  ! Arguments

  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), intent(inout) :: phi

  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), intent(out) :: laplace_phi

  ! Local variables

  integer :: idim, iBLK
  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK) :: dphi,ddphi
  integer :: i,j,k
  real, dimension(-1:1,-1:1,-1:1) :: phiC
  !---------------------------------------------------------------------------

  call proj_boundphi(phi,1,3)

  if(UseConstrainB)then                      !^CFG IF CONSTRAINB BEGIN
     do iBLK=1,nBlock
        if(unusedBLK(iBLK))CYCLE

        if(body_BLK(iBLK))then
           ! If some cells are inside the body, phi should have zero gradient
           ! accross the face between the body and the physical cell so that
           ! there is correction done for Bface on that face
           do k=1,nK; do j=1,nJ; do i=1,nI

              if(.not.true_cell(i,j,k,iBLK)) CYCLE

              where(true_cell(i-1:i+1,j-1:j+1,k-1:k+1,iBLK))
                 phiC=phi(i-1:i+1,j-1:j+1,k-1:k+1,iBLK)
              elsewhere
                 phiC=phi(i,j,k,iBLK)
              end where

              laplace_phi(i,j,k,iBLK)= &
                   (phiC(1,0,0)+phiC(-1,0,0)-2*phiC(0,0,0))/dx_BLK(iBLK)**2 +&
                   (phiC(0,1,0)+phiC(0,-1,0)-2*phiC(0,0,0))/dy_BLK(iBLK)**2 +&
                   (phiC(0,0,1)+phiC(0,0,-1)-2*phiC(0,0,0))/dz_BLK(iBLK)**2

           end do; end do; end do
        else
           laplace_phi(1:nI,1:nJ,1:nK,iBLK)= &
                (phi(2:nI+1,1:nJ,1:nK,iBLK)+phi(0:nI-1,1:nJ,1:nK,iBLK) &
                -2*phi(1:nI,1:nJ,1:nK,iBLK))/dx_BLK(iBLK)**2 +         &
                (phi(1:nI,2:nJ+1,1:nK,iBLK)+phi(1:nI,0:nJ-1,1:nK,iBLK) &
                -2*phi(1:nI,1:nJ,1:nK,iBLK))/dy_BLK(iBLK)**2 +         &
                (phi(1:nI,1:nJ,2:nK+1,iBLK)+phi(1:nI,1:nJ,0:nK-1,iBLK) &
                -2*phi(1:nI,1:nJ,1:nK,iBLK))/dz_BLK(iBLK)**2
        endif

     end do

     RETURN
  end if                                     !^CFG END CONSTRAINB

  call proj_gradient(1,phi,dphi)

  call proj_boundphi(dphi,1,1)

  call proj_gradient(1,dphi,laplace_phi)

  do idim=2,3
     call proj_gradient(idim,phi,dphi)

     call proj_boundphi(dphi,idim,idim)

     call proj_gradient(idim,dphi,ddphi)

     call add_BLK(laplace_phi,ddphi)

  end do

end subroutine proj_matvec

!=============================================================================
! Calculate gradient of phi in direction idim for real cells only
subroutine proj_gradient(idim,phi,dphi)
  use ModMain, ONLY : nI,nJ,nK,nBLK,nBlock,unusedBLK
  use ModGeometry, ONLY : fAx_BLK,fAy_BLK,fAz_BLK,vInv_CB
  implicit none

  ! Arguments
  integer, intent(in) :: idim
  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), intent(in) :: phi

  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), intent(out) :: dphi

  ! Local variables
  integer :: iBLK

  logical :: oktest, oktest_me

  !---------------------------------------------------------------------------

  call set_oktest('proj_gradient',oktest,oktest_me)

  if(oktest_me)write(*,*)'Gradient, idim=',idim 
  !call show_BLK('phi',phi)

  !dphi=0.0

  do iBLK=1,nBlock
     if(unusedBLK(iBLK))CYCLE
     select case(idim)
     case(1)
        dphi(1:nI,1:nJ,1:nK,iBLK) = 0.5*vInv_CB(:,:,:,iBLK)*fAx_BLK(iBLK)* &
             (phi(2:nI+1,1:nJ,1:nK,iBLK)-phi(0:nI-1,1:nJ,1:nK,iBLK))
     case(2)
        dphi(1:nI,1:nJ,1:nK,iBLK) = 0.5*vInv_CB(:,:,:,iBLK)*fAy_BLK(iBLK)* &
             (phi(1:nI,2:nJ+1,1:nK,iBLK)-phi(1:nI,0:nJ-1,1:nK,iBLK))

     case(3)
        dphi(1:nI,1:nJ,1:nK,iBLK) = 0.5*vInv_CB(:,:,:,iBLK)*fAz_BLK(iBLK)* &
             (phi(1:nI,1:nJ,2:nK+1,iBLK)-phi(1:nI,1:nJ,0:nK-1,iBLK))
     end select
  end do ! All blocks are done

  !call show_BLK('dphi',dphi)

!!! CORRECT dphi FOR MESH REFINEMENT EFFECTS

  return
end subroutine proj_gradient

!=============================================================================
! Calculate boundary values for phi for dimensions idimmin..idimmax
subroutine proj_boundphi(phi,idimmin,idimmax)
  use ModMain, ONLY : nI,nJ,nK,nBLK,nBlock,unusedBLK,prolong_order, UseBatl
  use ModMain, ONLY : UseConstrainB                   !^CFG IF CONSTRAINB
  use ModGeometry, ONLY : body_BLK, true_cell
  use ModParallel, ONLY : NOBLK,neiLtop,neiLbot,neiLeast,neiLwest,neiLnorth,neiLsouth
  use ModMessagePass, ONLY: message_pass_dir
  use ModProject
  use BATL_lib, ONLY: message_pass_cell_scalar

  implicit none

  ! Arguments
  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), &
       intent(inout) :: phi
  integer, intent(in):: idimmin,idimmax

  ! Local variables
  integer :: iBLK

  !---------------------------------------------------------------------------
  if(UseBatl)then
     call message_pass_cell_scalar(Phi,nWidthIn=1, nProlongOrderIn=1, &
          DoSendCornerIn=.false., DoRestrictFaceIn = .true.)
  else
     call message_pass_dir(idimmin,idimmax,1,.false.,prolong_order,1,phi,&
          restrictface=.true.)
  end if

  do iBLK=1,nBlock
     if(unusedBLK(iBLK))CYCLE

     if(neiLeast(iBLK) ==NOBLK)phi(-1:0     ,:,:,iBLK)=0.0
     if(neiLwest(iBLK) ==NOBLK)phi(nI+1:nI+2,:,:,iBLK)=0.0
     if(neiLsouth(iBLK)==NOBLK)phi(:,-1:0     ,:,iBLK)=0.0
     if(neiLnorth(iBLK)==NOBLK)phi(:,nJ+1:nJ+2,:,iBLK)=0.0
     if(neiLbot(iBLK)  ==NOBLK)phi(:,:,-1:0     ,iBLK)=0.0
     if(neiLtop(iBLK)  ==NOBLK)phi(:,:,nK+1:nK+2,iBLK)=0.0

     if(body_BLK(iBLK))then
        where(.not.true_cell(:,:,:,iBLK))phi(:,:,:,iBLK)=0.0
     end if

     if(UseConstrainB)then                      !^CFG IF CONSTRAINB BEGIN
        ! Correct ghost cells at resolution changes to get consistent gradphi

        if(neiLeast(iBLK) ==-1)phi(0   ,1:nJ,1:nK,iBLK)=2*phi(0   ,1:nJ,1:nK,iBLK) &
                                                         -phi(1   ,1:nJ,1:nK,iBLK)
        if(neiLwest(iBLK) ==-1)phi(nI+1,1:nJ,1:nK,iBLK)=2*phi(nI+1,1:nJ,1:nK,iBLK) &
                                                         -phi(nI  ,1:nJ,1:nK,iBLK)
        if(neiLsouth(iBLK)==-1)phi(1:nI,0   ,1:nK,iBLK)=2*phi(1:nI,0   ,1:nK,iBLK) &
                                                         -phi(1:nI,1   ,1:nK,iBLK)
        if(neiLnorth(iBLK)==-1)phi(1:nI,nJ+1,1:nK,iBLK)=2*phi(1:nI,nJ+1,1:nK,iBLK) &
                                                         -phi(1:nI,nJ  ,1:nK,iBLK)
        if(neiLbot(iBLK)  ==-1)phi(1:nI,1:nJ,0   ,iBLK)=2*phi(1:nI,1:nJ,0   ,iBLK) &
                                                         -phi(1:nI,1:nJ,1   ,iBLK)
        if(neiLtop(iBLK)  ==-1)phi(1:nI,1:nJ,nK+1,iBLK)=2*phi(1:nI,1:nJ,nK+1,iBLK) &
                                                         -phi(1:nI,1:nJ,nK  ,iBLK)
     end if                                     !^CFG END CONSTRAINB

  end do

end subroutine proj_boundphi

!=============================================================================
! Correct B field by gradient of phi
subroutine proj_correction(phi)
  use ModMain, ONLY : nI,nJ,nK,nBLK,Itest,Jtest,Ktest,BLKtest, &
       nBlock,unusedBLK
  use ModVarIndexes, ONLY : Bx_,By_,Bz_
  use ModAdvance,    ONLY : State_VGB
  use ModGeometry,   ONLY : dx_BLK,dy_BLK,dz_BLK,true_cell
  use ModProject
  use ModMain, ONLY : UseConstrainB             !^CFG IF CONSTRAINB
  use ModCT                                     !^CFG IF CONSTRAINB
  use ModEnergy, ONLY: calc_energy_cell         !^CFG IF CONSTRAINB
  implicit none

  ! Arguments
  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), intent(inout) :: phi

  ! Local variables
  integer :: iBLK, i, j, k
  real    :: DxInvHalf, DyInvHalf, DzInvHalf

  logical :: oktest, oktest_me

  !---------------------------------------------------------------------------

  call set_oktest('proj_correction',oktest,oktest_me)

  if(UseConstrainB)then                      !^CFG IF CONSTRAINB BEGIN

     if(oktest_me)write(*,*)'proj_correction old Bzface=',&
          BzFace_BLK(Itest,Jtest,Ktest,BLKtest), &
          BzFace_BLK(Itest,Jtest,Ktest+1,BLKtest), &
          true_cell(Itest,Jtest,Ktest-1,BLKtest), &
          true_cell(Itest,Jtest,Ktest+1,BLKtest)

     do iBLK=1,nBlock
        if(unusedBLK(iBLK)) CYCLE

        BxFace_BLK(1:nI+1,1:nJ,1:nK,iBLK)=BxFace_BLK(1:nI+1,1:nJ,1:nK,iBLK) &
             -(phi(1:nI+1,1:nJ,1:nK,iBLK)-phi(0:nI,1:nJ,1:nK,iBLK)) &
             /dx_BLK(iBLK)

        ByFace_BLK(1:nI,1:nJ+1,1:nK,iBLK)=ByFace_BLK(1:nI,1:nJ+1,1:nK,iBLK) &
             -(phi(1:nI,1:nJ+1,1:nK,iBLK)-phi(1:nI,0:nJ,1:nK,iBLK)) &
             /dy_BLK(iBLK)

        BzFace_BLK(1:nI,1:nJ,1:nK+1,iBLK)=BzFace_BLK(1:nI,1:nJ,1:nK+1,iBLK) &
             -(phi(1:nI,1:nJ,1:nK+1,iBLK)-phi(1:nI,1:nJ,0:nK,iBLK)) &
             /dz_BLK(iBLK)

        if(oktest_me.and.BLKtest==iBLK)write(*,*)'before bound_Bface Bzface=',&
             BzFace_BLK(Itest,Jtest,Ktest,BLKtest), &
             BzFace_BLK(Itest,Jtest,Ktest+1,BLKtest)

        ! Make sure that the correction is NOT applied to the body boundaries
        call bound_Bface(iBLK)
        ! Recalculate the cell centered B
        call Bface2Bcenter(iBLK)
        ! Keep pressure and modify energy for sake of positivity
        call calc_energy_cell(iBLK)
     end do

     if(oktest_me)write(*,*)'proj_correction new Bzface=',&
          BzFace_BLK(Itest,Jtest,Ktest,BLKtest), &
          BzFace_BLK(Itest,Jtest,Ktest+1,BLKtest)

  else                                       !^CFG END CONSTRAINB
     do iBLK = 1, nBlock
        if(unusedBLK(iBLK)) CYCLE
        DxInvHalf = 0.5/dx_BLK(iBLK);
        DyInvHalf = 0.5/dy_BLK(iBLK);
        DzInvHalf = 0.5/dz_BLK(iBLK);
        do k=1,nK; do j=1,nJ; do i=1,nI
           if(.not.true_cell(i,j,k,iBLK)) CYCLE
           State_VGB(Bx_,i,j,k,iBLK) = State_VGB(Bx_,i,j,k,iBLK) - &
                DxInvHalf*(phi(i+1,j,k,iBLK)-phi(i-1,j,k,iBLK))

           State_VGB(By_,i,j,k,iBLK) = State_VGB(By_,i,j,k,iBLK) - &
                DyInvHalf*(phi(i,j+1,k,iBLK)-phi(i,j-1,k,iBLK))

           State_VGB(Bz_,i,j,k,iBLK) = State_VGB(Bz_,i,j,k,iBLK) - &
                DzInvHalf*(phi(i,j,k+1,iBLK)-phi(i,j,k-1,iBLK))
        end do; end do; end do
     end do
  end if                                     !^CFG IF CONSTRAINB

end subroutine proj_correction
