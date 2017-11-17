!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModProjectDivB

  use BATL_lib, ONLY: &
       test_start, test_stop, iTest, jTest, kTest, iBlockTest

  ! Parameters for projection scheme:
  !
  ! proj_method determines the iterative scheme to be used
  ! if proj_method='cg'       then Conjugate Gradient method is used
  !                           (this is only good for symmetric matrices! )
  !    proj_method='bicgstab' then BiConjugate Gradient method is used
  !
  ! proj_typestop determines the stopping condition for the iterations.
  ! if proj_typestop='rel' then sum((div B)**2) is reduced by a factor
  !    proj_divbcoeff
  ! if proj_typestop='max' then |div B| is kept below a limit determeined by
  !    proj_divbcoeff and
  !    proj_divbconst
  !    if proj_divbcoeff>0 then
  !       |div B|< proj_divbcoeff * divbmax_0
  !          where divbmax_0 is the max(|div B|) after the first time step
  !    otherwise
  !       |div B|< proj_divbconst
  !
  ! proj_matvecmax is an upper limit on matrix-vector products for iterations
  !
  ! proj_boundtype determines the boundary conditions for the Poisson equation
  ! if proj_boundtype='zero' : phi=0 in ghost cells outside of comput. domain

  use ModSize, ONLY: nI, nJ, nK, nBlock, &
       MinI, MaxI, MinJ, MaxJ, MinK, MaxK, MaxBlock

  use ModProcMH, ONLY: iProc, nProc, iComm
  use ModMain, ONLY: Unused_B
  use ModGeometry, ONLY: true_BLK, true_cell
  use ModMpi

  implicit none
  SAVE

  private ! except

  public:: project_divb            ! project to div B=0 with Poisson solve
  public:: read_project_divb_param ! read parameters for projection scheme
  public:: proj_get_divb           ! Calculate div B for projection or CT

  ! Maximum value for divB (absolute or relative)
  real, public :: DivbMax = -1.0

  ! Local variables

  character (len=10):: proj_method   ='cg        '
  character (len=3)::  proj_typestop ='rel'
  real ::              proj_divbcoeff=0.1
  real ::              proj_divbconst=0.0
  integer ::           proj_matvecmax=50

  ! Counter for matrix vector multiplications, and for errors of solver


  ! Minimum value for divB (a small number)

  real, parameter :: divbmin=1E-10

contains
  !============================================================================
  subroutine read_project_divb_param(NameCommand)

    use ModReadParam, ONLY: read_var

    character(len=*), intent(in):: NameCommand

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'read_project_divb_param'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    select case(NameCommand)
    case("#PROJECTION")
       call read_var('TypeProjectIter' ,proj_method)
       call read_var('TypeProjectStop' ,proj_typestop)
       call read_var('RelativeLimit'   ,proj_divbcoeff)
       call read_var('AbsoluteLimit'   ,proj_divbconst)
       call read_var('MaxMatvec'       ,proj_matvecmax)
       ! Make sure that DivbMax is recalculated
       DivbMax = -1.0
    case default
       call stop_mpi(NameSub//': unknown command='//NameCommand)
    end select

    call test_stop(NameSub, DoTest)
  end subroutine read_project_divb_param
  !============================================================================
  subroutine project_divb

    ! Project B according to B'=B-grad phi, where Laplace phi=div B
    ! Solve the Poisson equation with an iterative scheme
    !
    ! See G. Toth, 2000, Journal of Computational Physics, 161, 605-652

    use ModProcMH
    use ModVarIndexes, ONLY: Bx_,Bz_,P_
    use ModAdvance, ONLY: State_VGB
    use ModGeometry, ONLY: true_cell
    use ModMain, ONLY: UseConstrainB
    use ModConstrainDivB, ONLY: Bxface_BLK, Byface_BLK, Bzface_BLK
    use ModMessagePass, ONLY: exchange_messages
    use BATL_lib, ONLY: Xyz_DGB

    ! Local variables
    ! proj_divb: original error to be projected out
    ! phi:       Scalar field in the Poisson problem
    real, allocatable, save:: proj_divb(:,:,:,:), phi(:,:,:,:)

    integer :: info      ! error status of the iterative scheme
    integer :: nmatvec   ! number of matvex operations performed
    real ::    resid     ! residual after iterative solver stopped

    real :: divbmax_now  ! current value of max(|div B|)
    real :: pmin_old, pmin_new     ! original and projected value of min(p)

    integer :: loc(5)

    ! Functions:

    real, external :: maxval_abs_BLK, minval_BLK, maxval_loc_abs_BLK, &
         minval_loc_BLK

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'project_divb'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    call timing_start(NameSub)

    if(.not.allocated(proj_divb)) allocate( &
         proj_divb(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock), &
         phi(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))

    if(DoTest)then
       write(*,*) NameSub, ': old B:', &
            State_VGB(Bx_:Bz_,iTest,jTest,kTest,iBlockTest),  &
            true_cell(iTest,jTest,kTest,iBlockTest)
       if(UseConstrainB)then
          write(*,*) NameSub, ': Bxface,true_cell:',    &
               Bxface_BLK(iTest,jTest,kTest,iBlockTest),  &
               Bxface_BLK(iTest+1,jTest,kTest,iBlockTest),&
               true_cell(iTest-1,jTest,kTest,iBlockTest), &
               true_cell(iTest+1,jTest,kTest,iBlockTest)

          write(*,*) NameSub, ': Byface,true_cell:',    &
               Byface_BLK(iTest,jTest,kTest,iBlockTest),  &
               Byface_BLK(iTest,jTest+1,kTest,iBlockTest),&
               true_cell(iTest,jTest-1,kTest,iBlockTest), &
               true_cell(iTest,jTest+1,kTest,iBlockTest)

          write(*,*) NameSub, ': Bzface,true_cell:',    &
               Bzface_BLK(iTest,jTest,kTest,iBlockTest),  &
               Bzface_BLK(iTest,jTest,kTest+1,iBlockTest),&
               true_cell(iTest,jTest,kTest-1,iBlockTest), &
               true_cell(iTest,jTest,kTest+1,iBlockTest)

       end if
    end if

    ! Determine the div B error
    call proj_get_divb(proj_divB)

    if(DoTest)write(*,*)'proj_divb=',proj_divb(iTest,jTest,kTest,iBlockTest)

    ! DEBUG
    ! call proj_bound(proj_divb)
    ! call show_BLK('divb',proj_divb)
    ! DEBUG
    ! call stop_mpi('debug')

    if(DoTest)then
       pmin_old=minval_BLK(nProc,State_VGB(P_,:,:,:,:))
       divbmax_now=maxval_abs_BLK(nProc,proj_divB)
       if(DoTest) write(*,*) NameSub, &
            ': old max(abs(divB)), min(p)=',divbmax_now,pmin_old
    endif

    ! Calculate and save divbmax, the maximum div B allowed
    if(divbmax < divbmin)then
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
       if(DoTest)write(*,*)'Allowed maximum for div B:',divbmax
       if(divbmax<divbmin)&
            call stop_mpi(NameSub//': Too small value for divbmax')
    endif

    ! Solve the Poisson equation Laplace phi = div B
    call proj_poisson(proj_divb,divbmax,proj_typestop,proj_matvecmax, &
         info,nmatvec,resid,phi)

    if(DoTest)write(*,*) NameSub, ': Poisson solver info, nmatvec, resid:',&
         info, nmatvec, resid

    ! Do not do anything if the initial guess satisfied the stopping criterion
    if(info==3)GOTO 100

    ! Do not subtract grad(phi) from B if the iterations did not reduce error
    if(info<0)GOTO 100

    ! Get the ghost cell values for the solution phi
    call proj_boundphi(phi);

    ! if(DoTest)write(*,*)'proj_divb,lapl(phi)=',&
    !     proj_divb(iTest,jTest,kTest,iBlockTest),&
    !     (phi(iTest+2,jTest,kTest,iBlockTest)+phi(iTest-2,jTest,kTest,iBlockTest) &
    !     +phi(iTest,jTest+2,kTest,iBlockTest)+phi(iTest,jTest-2,kTest,iBlockTest) &
    !     +phi(iTest,jTest,kTest-2,iBlockTest)+phi(iTest,jTest,kTest-2,iBlockTest) &
    !     -6*phi(iTest,jTest,kTest,iBlockTest))/CellSize_DB(x_,iBlockTest)**2/4

    ! Correct B field: B'=B-dB where dB=grad(phi)
    call proj_correction(phi)

    ! Recalculate boundary conditions
    call exchange_messages

!!! StringTest=StringTest(1:50)

    ! Testing div B'=0
    if(DoTest)then
       call proj_get_divb(proj_divB)

       if(DoTest)write(*,*)'after proj_divb=',&
            proj_divb(iTest,jTest,kTest,iBlockTest)

       divbmax_now=maxval_loc_abs_BLK(nProc,proj_divB,loc)
       if(abs(resid-divbmax_now)/(divbmax_now+1.e-12) > 0.1)then
          if(iProc==loc(5))write(*,*)&
                NameSub, ': resid,max(abs(divB)),loc,X,Y,Z=',&
               resid, divbmax_now, loc, &
               Xyz_DGB(:,loc(1),loc(2),loc(3),loc(4))
       else
          if(DoTest)write(*,*) NameSub, ': resid,max(abs(divB))',&
               resid, divbmax_now
       end if
       pmin_new   =minval_loc_BLK(nProc,State_VGB(P_,:,:,:,:),loc)
       if(pmin_new<0.5*pmin_old)then
          if(iProc==loc(5))write(*,*) NameSub, ': min(p),loc,X,Y,Z=',&
               pmin_new,loc,Xyz_DGB(:,loc(1),loc(2),loc(3),loc(4))
       else
          if(DoTest)write(*,*) NameSub, ': new min(p)',pmin_new
       endif

       if(DoTest)then
          write(*,*) NameSub, ': new B, true_cell:', &
               State_VGB(Bx_:Bz_,iTest,jTest,kTest,iBlockTest),  &
               true_cell(iTest,jTest,kTest,iBlockTest)
          if(UseConstrainB)then
             write(*,*) NameSub, ': Bxface,true_cell:',    &
                  Bxface_BLK(iTest,jTest,kTest,iBlockTest),  &
                  Bxface_BLK(iTest+1,jTest,kTest,iBlockTest),&
                  true_cell(iTest-1,jTest,kTest,iBlockTest), &
                  true_cell(iTest+1,jTest,kTest,iBlockTest)

             write(*,*) NameSub, ': Byface,true_cell:',    &
                  Byface_BLK(iTest,jTest,kTest,iBlockTest),  &
                  Byface_BLK(iTest,jTest+1,kTest,iBlockTest),&
                  true_cell(iTest,jTest-1,kTest,iBlockTest), &
                  true_cell(iTest,jTest+1,kTest,iBlockTest)

             write(*,*) NameSub, ': Bzface,true_cell:',    &
                  Bzface_BLK(iTest,jTest,kTest,iBlockTest),  &
                  Bzface_BLK(iTest,jTest,kTest+1,iBlockTest),&
                  true_cell(iTest,jTest,kTest-1,iBlockTest), &
                  true_cell(iTest,jTest,kTest+1,iBlockTest)

          end if
       end if
    endif

    ! if(DoTest)write(*,*)'final proj_divb=',&
    !        proj_divb(iTest,jTest,kTest,iBlockTest)

100 call timing_stop(NameSub)

    call test_stop(NameSub, DoTest)
  end subroutine project_divb
  !============================================================================

  subroutine proj_get_divb(proj_divB)

    ! Calculate div B using simple finite differences
    ! Do corrections for mesh refinement

    use ModVarIndexes, ONLY: Bx_, By_, Bz_
    use ModAdvance, ONLY: State_VGB
    use ModMain, ONLY: UseConstrainB
    use ModConstrainDivB, ONLY: Bxface_BLK, Byface_BLK, Bzface_BLK
    use BATL_lib, ONLY:  nI,nJ,nK, MaxBlock, nBlock, Unused_B, &
         x_, y_, z_, CellSize_DB

    ! Argument

    real, intent(out) :: proj_divb(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

    ! Local variables
    integer :: iBlock, i, j, k
    real    :: DxInvHalf, DyInvHalf, DzInvHalf
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'proj_get_divb'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    proj_divB=0.0

    if(UseConstrainB)then
       do iBlock=1,nBlock
          if(Unused_B(iBlock)) CYCLE

          proj_divB(1:nI,1:nJ,1:nK,iBlock)= &
               (Bxface_BLK(2:nI+1,1:nJ  ,1:nK  ,iBlock)                      &
               -Bxface_BLK(1:nI  ,1:nJ  ,1:nK  ,iBlock))/CellSize_DB(x_,iBlock)&
               +(Byface_BLK(1:nI  ,2:nJ+1,1:nK  ,iBlock)                      &
               -Byface_BLK(1:nI  ,1:nJ  ,1:nK  ,iBlock))/CellSize_DB(y_,iBlock)&
               +(Bzface_BLK(1:nI  ,1:nJ  ,2:nK+1,iBlock)                      &
               -Bzface_BLK(1:nI  ,1:nJ  ,1:nK  ,iBlock))/CellSize_DB(z_,iBlock)
       end do
    else
       do iBlock=1,nBlock
          if(Unused_B(iBlock))CYCLE
          DxInvHalf = 0.5/CellSize_DB(x_,iBlock);
          DyInvHalf = 0.5/CellSize_DB(y_,iBlock);
          DzInvHalf = 0.5/CellSize_DB(z_,iBlock);
          do k=1,nK; do j=1,nJ; do i=1,nI
             proj_divb(i,j,k,iBlock) = &
                  DxInvHalf* &
                  (State_VGB(Bx_,i+1,j,k,iBlock)-State_VGB(Bx_,i-1,j,k,iBlock)) + &
                  DyInvHalf* &
                  (State_VGB(By_,i,j+1,k,iBlock)-State_VGB(By_,i,j-1,k,iBlock)) + &
                  DzInvHalf* &
                  (State_VGB(Bz_,i,j,k+1,iBlock)-State_VGB(Bz_,i,j,k-1,iBlock))
          end do; end do; end do
       end do
    end if

    call test_stop(NameSub, DoTest)
  end subroutine proj_get_divb
  !============================================================================

  ! This is just an interface to the various iterative schemes to
  ! Solve grad div phi = rhs
  subroutine proj_poisson(rhs,tolerance,typestop,matvecmax,&
       info,nmatvec,resid,phi)
    use ModProcMH
    use ModMain, ONLY:MaxBlock

    ! Arguments

    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock), intent(inout):: rhs
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

    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock), intent(out):: phi
    !    The solution

    ! Local variables

    integer :: ierror=-1

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'proj_poisson'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    if(DoTest)write(*,*)'Proj_Poisson solver called with ',proj_method

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

    if(DoTest)write(*,*)'Poisson info,nmatvec,resid',info,nmatvec,resid

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

    call test_stop(NameSub, DoTest)
  end subroutine proj_poisson
  !============================================================================

  ! Calculate Laplace phi
  subroutine proj_matvec(phi,laplace_phi)
    use ModMain, ONLY : MaxBlock,nBlock,Unused_B,nI,nJ,nK, x_, y_, z_
    use ModGeometry, ONLY : true_cell,body_BLK
    use ModMain, ONLY : UseConstrainB
!    use ModCT
    use BATL_lib, ONLY: CellSize_DB

    ! Arguments

    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock), intent(inout) :: phi

    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock), intent(out) :: laplace_phi

    ! Local variables

    integer :: idim, iBlock
    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock) :: dphi,ddphi
    integer :: i,j,k
    real :: phiC(-1:1,-1:1,-1:1), InvDx2, InvDy2, InvDz2
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'proj_matvec'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    call proj_boundphi(phi)

    if(UseConstrainB)then
       do iBlock=1,nBlock
          if(Unused_B(iBlock))CYCLE

          if(body_BLK(iBlock))then
             ! If some cells are inside the body, phi should have zero gradient
             ! accross the face between the body and the physical cell so that
             ! there is correction done for Bface on that face
             do k = 1, nK; do j = 1, nJ; do i = 1, nI

                if(.not.true_cell(i,j,k,iBlock)) CYCLE

                where(true_cell(i-1:i+1,j-1:j+1,k-1:k+1,iBlock))
                   phiC=phi(i-1:i+1,j-1:j+1,k-1:k+1,iBlock)
                elsewhere
                   phiC=phi(i,j,k,iBlock)
                end where

                laplace_phi(i,j,k,iBlock)= &
                     (phiC(1,0,0)+phiC(-1,0,0)-2*phiC(0,0,0)) &
                     /CellSize_DB(x_,iBlock)**2 + &
                     (phiC(0,1,0)+phiC(0,-1,0)-2*phiC(0,0,0)) &
                     /CellSize_DB(y_,iBlock)**2 + &
                     (phiC(0,0,1)+phiC(0,0,-1)-2*phiC(0,0,0)) &
                     /CellSize_DB(z_,iBlock)**2
             end do; end do; end do
          else
             InvDx2 = 1/CellSize_DB(x_,iBlock)**2
             InvDy2 = 1/CellSize_DB(y_,iBlock)**2
             InvDz2 = 1/CellSize_DB(z_,iBlock)**2
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                laplace_phi(i,j,k,iBlock)= &
                     (phi(i+1,j,k,iBlock) + phi(i-1,j,k,iBlock) &
                     -2*phi(i,j,k,iBlock))*InvDx2 + &
                     (phi(i,j+1,k,iBlock) + phi(i,j-1,k,iBlock) &
                     -2*phi(i,j,k,iBlock))*InvDy2 + &
                     (phi(i,j,k+1,iBlock) + phi(i,j,k-1,iBlock) &
                     -2*phi(i,j,k,iBlock))*InvDz2
             end do; end do; end do
          endif

       end do

       RETURN
    end if

    call proj_gradient(1,phi,dphi)

    call proj_boundphi(dphi)

    call proj_gradient(1,dphi,laplace_phi)

    do idim=2,3
       call proj_gradient(idim,phi,dphi)

       call proj_boundphi(dphi)

       call proj_gradient(idim,dphi,ddphi)

       call add_BLK(laplace_phi,ddphi)

    end do

    call test_stop(NameSub, DoTest)
  end subroutine proj_matvec
  !============================================================================

  subroutine proj_gradient(idim,phi,dphi)

    ! Calculate gradient of phi in direction idim

    use ModSize, ONLY: x_, y_, z_, &
         MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nI, nJ, nK, nBlock, MaxBlock
    use ModMain, ONLY : Unused_B
    use BATL_lib, ONLY: CellSize_DB

    ! Arguments
    integer, intent(in) :: idim
    real, intent(in) :: phi(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

    real, intent(out) :: dphi(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

    ! Local variables
    integer :: i, j, k, iBlock
    real:: Coef

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'proj_gradient'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    if(DoTest)write(*,*)'Gradient, idim=',idim
    ! call show_BLK('phi',phi)

    ! dphi=0.0

    do iBlock = 1, nBlock
       if(Unused_B(iBlock))CYCLE
       select case(idim)
       case(1)
          Coef = 0.5/CellSize_DB(x_,iBlock)
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             dphi(i,j,k,iBlock) = Coef*(phi(i+1,j,k,iBlock) - phi(i-1,j,k,iBlock))
          end do; end do; end do
       case(2)
          Coef = 0.5/CellSize_DB(y_,iBlock)
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             dphi(i,j,k,iBlock) = Coef*(phi(i,j+1,k,iBlock) - phi(i,j-1,k,iBlock))
          end do; end do; end do
       case(3)
          Coef = 0.5/CellSize_DB(z_,iBlock)
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             dphi(i,j,k,iBlock) = Coef*(phi(i,j,k+1,iBlock) - phi(i,j,k-1,iBlock))
          end do; end do; end do
       end select
    end do ! All blocks are done

    ! call show_BLK('dphi',dphi)

!!! CORRECT dphi FOR MESH REFINEMENT EFFECTS

    call test_stop(NameSub, DoTest)
  end subroutine proj_gradient
  !============================================================================

  subroutine proj_boundphi(phi)

    ! Calculate boundary values for phi for dimensions

    use ModMain, ONLY : MaxBlock, nBlock, Unused_B
    use ModGeometry, ONLY : body_BLK, true_cell
    use ModParallel, ONLY : &
         NOBLK,neiLtop,neiLbot,neiLeast,neiLwest,neiLnorth,neiLsouth
    use BATL_lib, ONLY: message_pass_cell, &
         MinI, MaxI, j0_, nJp1_, MinJ, MaxJ, k0_, nKp1_, MinK, MaxK

    ! Arguments
    real,    intent(inout) :: phi(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

    ! Local variables
    integer :: iBlock

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'proj_boundphi'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    call message_pass_cell(Phi, nWidthIn=1, nProlongOrderIn=1, &
         DoSendCornerIn=.false., DoRestrictFaceIn = .true.)

    do iBlock = 1, nBlock
       if(Unused_B(iBlock))CYCLE

       if(neiLeast(iBlock) ==NOBLK) phi(MinI:0   ,:,:,iBlock) = 0.0
       if(neiLwest(iBlock) ==NOBLK) phi(nI+1:MaxI,:,:,iBlock) = 0.0
       if(neiLsouth(iBlock)==NOBLK) phi(:,MinJ :j0_ ,:,iBlock) = 0.0
       if(neiLnorth(iBlock)==NOBLK) phi(:,nJp1_:MaxJ,:,iBlock) = 0.0

       if(nK>1)then
          if(neiLbot(iBlock)==NOBLK) phi(:,:,MinK :k0_ ,iBlock) = 0.0
          if(neiLtop(iBlock)==NOBLK) phi(:,:,nKp1_:MaxK,iBlock) = 0.0
       end if

       if(body_BLK(iBlock))then
          where(.not.true_cell(:,:,:,iBlock)) phi(:,:,:,iBlock)=0.0
       end if

       ! if(UseConstrainB)then
       !   ! Correct ghost cells at resolution changes to get consistent gradphi
       !
       !   if(neiLeast(iBlock) ==-1) &
       !        phi(0   ,1:nJ,1:nK,iBlock)=2*phi(0   ,1:nJ,1:nK,iBlock) &
       !        -phi(1   ,1:nJ,1:nK,iBlock)
       !   if(neiLwest(iBlock) ==-1) &
       !        phi(nI+1,1:nJ,1:nK,iBlock)=2*phi(nI+1,1:nJ,1:nK,iBlock) &
       !        -phi(nI  ,1:nJ,1:nK,iBlock)
       !   if(neiLsouth(iBlock)==-1) &
       !        phi(1:nI,0   ,1:nK,iBlock)=2*phi(1:nI,0   ,1:nK,iBlock) &
       !        -phi(1:nI,1   ,1:nK,iBlock)
       !   if(neiLnorth(iBlock)==-1) &
       !        phi(1:nI,nJ+1,1:nK,iBlock)=2*phi(1:nI,nJ+1,1:nK,iBlock) &
       !        -phi(1:nI,nJ  ,1:nK,iBlock)
       !   if(neiLbot(iBlock)  ==-1) &
       !        phi(1:nI,1:nJ,0   ,iBlock)=2*phi(1:nI,1:nJ,0   ,iBlock) &
       !        -phi(1:nI,1:nJ,1   ,iBlock)
       !   if(neiLtop(iBlock)  ==-1) &
       !        phi(1:nI,1:nJ,nK+1,iBlock)=2*phi(1:nI,1:nJ,nK+1,iBlock) &
       !        -phi(1:nI,1:nJ,nK  ,iBlock)
       ! end if

    end do

    call test_stop(NameSub, DoTest)
  end subroutine proj_boundphi
  !============================================================================

  subroutine proj_correction(phi)

    ! Correct B field by gradient of phi

    use ModMain, ONLY : nI,nJ,nK
    use ModVarIndexes, ONLY : Bx_,By_,Bz_
    use ModAdvance,    ONLY : State_VGB
    use ModGeometry,   ONLY : true_cell
    use ModMain, ONLY : UseConstrainB
    use ModConstrainDivB, ONLY: Bxface_BLK, Byface_BLK, Bzface_BLK, &
         Bface2Bcenter, bound_bface
    use ModEnergy, ONLY: calc_energy_cell
    use BATL_lib, ONLY: CellSize_DB, x_, y_, z_, nI, nJ, nK, nBlock,Unused_B

    ! Arguments
    real, intent(inout) :: phi(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

    ! Local variables
    integer :: iBlock, i, j, k
    real    :: DxInvHalf, DyInvHalf, DzInvHalf, DxInv, DyInv, DzInv

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'proj_correction'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    if(UseConstrainB)then

       if(DoTest)write(*,*)'proj_correction old Bzface=',&
            BzFace_BLK(iTest,jTest,kTest,iBlockTest), &
            BzFace_BLK(iTest,jTest,kTest+1,iBlockTest), &
            true_cell(iTest,jTest,kTest-1,iBlockTest), &
            true_cell(iTest,jTest,kTest+1,iBlockTest)

       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE

          DxInv = 1/CellSize_DB(x_,iBlock)
          do k = 1, nK; do j = 1, nJ; do i = 1, nI+1
             BxFace_BLK(i,j,k,iBlock) = BxFace_BLK(i,j,k,iBlock) &
                  - DxInv*(phi(i,j,k,iBlock) - phi(i-1,j,k,iBlock))
          end do; end do; end do

          DyInv = 1/CellSize_DB(y_,iBlock)
          do k = 1, nK; do j = 1, nJ+1; do i = 1, nI
             ByFace_BLK(i,j,k,iBlock) = ByFace_BLK(i,j,k,iBlock) &
                  - DyInv*(phi(i,j+1,k,iBlock) - phi(i,j-1,k,iBlock))
          end do; end do; end do

          if(nK > 1)then
             DzInv = 1/CellSize_DB(z_,iBlock)
             do k = 1, nK+1; do j = 1, nJ; do i = 1, nI
                BzFace_BLK(i,j,k,iBlock) = BzFace_BLK(i,j,k,iBlock) &
                     - DzInv*(phi(i,j,k+1,iBlock) - phi(i,j,k-1,iBlock))
             end do; end do; end do
          end if

          if(DoTest.and.iBlockTest==iBlock) &
               write(*,*)'before bound_Bface Bzface=',&
               BzFace_BLK(iTest,jTest,kTest,iBlockTest), &
               BzFace_BLK(iTest,jTest,kTest+1,iBlockTest)

          ! Make sure that the correction is NOT applied to the body boundaries
          call bound_Bface(iBlock)
          ! Recalculate the cell centered B
          call Bface2Bcenter(iBlock)
          ! Keep pressure and modify energy for sake of positivity
          call calc_energy_cell(iBlock)
       end do

       if(DoTest)write(*,*)'proj_correction new Bzface=',&
            BzFace_BLK(iTest,jTest,kTest,iBlockTest), &
            BzFace_BLK(iTest,jTest,kTest+1,iBlockTest)

    else
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          DxInvHalf = 0.5/CellSize_DB(x_,iBlock);
          DyInvHalf = 0.5/CellSize_DB(y_,iBlock);
          DzInvHalf = 0.5/CellSize_DB(z_,iBlock);
          do k=1,nK; do j=1,nJ; do i=1,nI
             if(.not.true_cell(i,j,k,iBlock)) CYCLE
             State_VGB(Bx_,i,j,k,iBlock) = State_VGB(Bx_,i,j,k,iBlock) - &
                  DxInvHalf*(phi(i+1,j,k,iBlock)-phi(i-1,j,k,iBlock))

             State_VGB(By_,i,j,k,iBlock) = State_VGB(By_,i,j,k,iBlock) - &
                  DyInvHalf*(phi(i,j+1,k,iBlock)-phi(i,j-1,k,iBlock))

             State_VGB(Bz_,i,j,k,iBlock) = State_VGB(Bz_,i,j,k,iBlock) - &
                  DzInvHalf*(phi(i,j,k+1,iBlock)-phi(i,j,k-1,iBlock))
          end do; end do; end do
       end do
    end if

    call test_stop(NameSub, DoTest)
  end subroutine proj_correction
  !============================================================================

  ! The CG-algorithm is implemented as shown on page 12 of the thesis
  ! "Preconditioning for sparse matrices with applications."
  ! Auke van der Ploeg, University of Groningen, 1994.
  ! Rewritten to F90 by G. Toth based on the F77 subroutine in src/conjgrad.f
  ! Rewritten by G. Toth for BATS-R-US code to be used with MPI, 2000.
  !
  ! This subroutine determines the solution of A.QX=RHS, where
  ! the matrix-vector multiplication with A is performed by
  ! the subroutine 'proj_matvec'.
  !
  ! !! If the matrix is not symmetric positive definite, CG is likely to fail.
  !
  subroutine proj_cg(rhs,qx,iter,tol,typestop,info)
    use ModProcMH
    use ModMain, ONLY:MaxBlock
    use ModAdvance, ONLY : tmp1_BLK,tmp2_BLK

    ! Arguments

    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock), intent(inout) :: rhs
    !        on input:  right-hand side vector.
    !        on output: residual vector.

    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock), intent(out):: qx
    !        on output: solution vector.

    integer, intent(inout) :: iter
    !       on input:  maximum number of iterations to be performed.
    !       on output: actual  number of iterations done.

    real, intent(inout) :: tol
    !       on input:  required (relative) 2-norm or maximum norm of residual
    !       on output: achieved (relative) 2-norm or maximum norm of residual

    character (len=3), intent(in) :: typestop
    !       Determine stopping criterion (||.|| denotes the 2-norm):
    !       typestop='rel'    -- relative stopping crit.: ||res|| <= tol*||res0||
    !       typestop='abs'    -- absolute stopping crit.: ||res|| <= tol
    !       typestop='max'    -- maximum  stopping crit.: max(abs(res)) <= tol

    integer, intent(out)   :: info
    !       Gives reason for returning:
    !     abs(info)=  0 - solution found satisfying given tolerance.
    !                 1 - iteration aborted due to division by very small value.
    !                 2 - no convergence within maximum number of iterations.
    !                 3 - initial guess satisfies the stopping criterion.
    !    sign(info)=  + - residual decreased
    !                 - - residual did not reduce

    ! Local variables

    integer ::itr,matv
    real :: rho,rhonew,res,res0,bet,alf,assumedzero

    ! external Functions:

    real :: maxval_abs_BLK

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'proj_cg'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    if(DoTest)write(*,*)'Proj_CG: tol, iter=',tol,iter

    ! Debug
    ! call show_BLK('rhs',rhs)

    assumedzero=1.E-16; itr=0; matv=0
    call set_BLK(qx,0.0)

    if (typestop/='rel'.and.typestop/='abs'.and.typestop/='max') &
         call stop_mpi('Error in CG: typestop='//typestop// &
         ' should be one of rel/abs/max.')

    if(DoTest) write(*,*)'n gives the number of CG-iterations.'

    ! rho=||rhs||
    rho=sqrt(dot_product_BLK(rhs,rhs))

    res0=rho
    if (typestop=='max') res0=maxval_abs_BLK(nProc,rhs)
    ! DEBUG
    ! write(*,*)'res0:',res0,' rhs:',rhs

    res=res0

    assumedzero = assumedzero*res0

    if (DoTest) then
       if (typestop=='max') then
          write(*,*)'n:',itr,' Maximum norm initial residual:',res0
       else
          write(*,*)'n:',itr,' 2-norm initial residual:',res0
       end IF
    end if

    if (res0<divbmin .or. (typestop/='rel'.and.res0<=tol)) then
       info = 3
    else
       ! Initialize rho and tmp1_BLK=Z
       rho=rho*rho
       call eq_BLK(tmp1_BLK,rhs)

       ! Do iteration
       do
          ! AZ=A.Z
          call proj_matvec(tmp1_BLK,tmp2_BLK)
          matv=matv+1

          ! Debug
          ! write(*,*)'Z =tmp1_BLK=',tmp1_BLK(:,:,:,1:2)
          ! write(*,*)'AZ=tmp2_BLK=',tmp2_BLK(:,:,:,1:2)
          ! call stop_mpi('Debug')

          ! alf=A.AZ
          alf=dot_product_BLK(tmp1_BLK,tmp2_BLK)

          if(DoTest)write(*,*)'alf=',alf
          ! Debug
          ! call show_BLK('Z',tmp1_BLK)
          ! call show_BLK('AZ',tmp2_BLK)
          ! call stop_mpi('Debug')
          if (abs(alf)<=assumedzero**2) then
             info = 1
             EXIT
          end if
          alf=rho/alf
          if(DoTest)write(*,*)'alf=',alf

          call add_times_BLK(qx,alf,tmp1_BLK)
          call add_times_BLK(rhs,-alf,tmp2_BLK)

          ! rhonew=||rhs||
          rhonew=sqrt(dot_product_BLK(rhs,rhs))
          if(DoTest)write(*,*)'rhonew=',rhonew

          select case(typestop)
          case('max')
             res=maxval_abs_BLK(nProc,rhs)

          case('rel')
             res=rhonew/res0

          case('abs')
             res=rhonew
          end select
          rhonew=rhonew*rhonew

          itr=itr+1
          if (DoTest) &
               write(*,*)'n:',itr,' ',typestop,'. norm of residual:',res

          if (res<=tol) then
             info = 0
             EXIT
          end if
          if (itr>=iter) then
             info = 2
             EXIT
          end if
          if (rho<=assumedzero**2) then
             info = 1
             EXIT
          end if

          bet=rhonew/rho
          if(DoTest)write(*,*)'bet=',bet

          call eq_plus_times_BLK(tmp1_BLK,rhs,bet,tmp1_BLK)

          if (DoTest) write(*,*)'alf,bet,rho,rhonew:',alf,bet,rho,rhonew
          rho=rhonew
       end do
    end if

    ! return number of iterations and achieved residual
    iter=itr
    tol =res
    if((typestop=='rel'.and.res>1.0).or.(typestop/='rel'.and.res>res0))info=-info

    ! report results
    if(DoTest)then
       write(*,*)'Total Number of CG-iterations:',itr
       write(*,*)'Number of matrix-vector mult.:',matv
       select case(abs(info))
       case(0)
          write(*,*)'Successful iteration, norm of res. is:',tol
       case(1)
          write(*,*)'Iteration aborted due to division by a'
          write(*,*)'very small value.'
       case(2)
          write(*,*)'Stopping crit. not fulfilled within given'
          write(*,*)'maximum number of iterations.'
       case(3)
          write(*,*)'Initial guess for the solution satisfies'
          write(*,*)'given stopping criterion.'
       case default
          write(*,*)'Impossible value for info:',info
       end select
       if(info<0)write(*,*)'The residual did not reduce'
    endif

    call test_stop(NameSub, DoTest)
  end subroutine proj_cg
  !============================================================================

  ! Simple BiCGstab(\ell=1) iterative method
  ! Modified by G.Toth from the \ell<=2 version written
  ! by M.A.Botchev, Jan.'98.
  ! Further modified for the BATS-R-US code (2000) by G. Toth
  !
  ! This is the "vanilla" version of BiCGstab(\ell) as described
  ! in PhD thesis of D.R.Fokkema, Chapter 3.  It includes two enhancements
  ! to BiCGstab(\ell) proposed by G.Sleijpen and H.van der Vorst in
  ! 1) G.Sleijpen and H.van der Vorst "Maintaining convergence
  !    properties of BiCGstab methods in finite precision arithmetic",
  !    Numerical Algorithms, 10, 1995, pp.203-223
  ! 2) G.Sleijpen and H.van der Vorst "Reliable updated residuals in
  !    hybrid BiCG methods", Computing, 56, 1996, 141-163
  !
  ! {{ This code is based on:
  ! subroutine bistbl v1.0 1995
  !
  ! Copyright (c) 1995 by D.R. Fokkema.
  ! Permission to copy all or part of this work is granted,
  ! provided that the copies are not made or distributed
  ! for resale, and that the copyright notice and this
  ! notice are retained.  }}
  !
  ! This subroutine determines the solution of A.QX=RHS, where
  ! the matrix-vector multiplication with A is performed by
  ! the subroutine 'proj_matvec'. For symmetric matrix use the more efficient
  ! proj_cg algorithm!
  !
  subroutine proj_bicgstab(rhs,qx,iter,tol,typestop,info)
    use ModProcMH
    use ModMain, ONLY:MaxBlock
    use ModAdvance, ONLY : tmp1_BLK,tmp2_BLK

    ! Arguments

    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock), intent(inout) :: rhs
    !        on input:  right-hand side vector.
    !        on output: residual vector.

    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock), intent(out):: qx
    !       on output: solution vector.

    integer, intent(inout) :: iter
    !       on input:  maximum number of iterations to be performed.
    !       on output: actual  number of iterations done.

    real, intent(inout) :: tol
    !       on input:  required (relative) 2-norm or maximum norm of residual
    !       on output: achieved (relative) 2-norm or maximum norm of residual

    character (len=3), intent(in) :: typestop
    !      Determine stopping criterion (||.|| denotes the 2-norm):
    !      typestop='rel'    -- relative stopping crit.: ||res|| <= tol*||res0||
    !      typestop='abs'    -- absolute stopping crit.: ||res|| <= tol
    !      typestop='max'    -- maximum  stopping crit.: max(abs(res)) <= tol

    ! NOTE for typestop='rel' and 'abs':
    !            To save computational work, the value of
    !            residual norm used to check the convergence inside the main
    !            iterative loop is computed from
    !            projections, i.e. it can be smaller than the true residual norm
    !            (it may happen when e.g. the 'matrix-free' approach is used).
    !            Thus, it is possible that the true residual does NOT satisfy
    !            the stopping criterion ('rel' or 'abs').
    !            The true residual norm (or residual reduction) is reported on
    !            output in parameter TOL -- this can be changed to save 1 MATVEC
    !            (see comments at the end of the subroutine)

    integer, intent(out)   :: info
    !       Gives reason for returning:
    !     abs(info)=  0 - solution found satisfying given tolerance.
    !                 1 - iteration aborted due to division by very small value.
    !                 2 - no convergence within maximum number of iterations.
    !                 3 - initial guess satisfies the stopping criterion.
    !    sign(info)=  + - residual decreased
    !                 - - residual did not reduce

    ! Local parameters

    integer, parameter :: qz_=1,zz_=3,y0_=5,yl_=6,qy_=7

    ! Local variables (2 vectors are needed in addition to tmp1_BLK=r and tmp2_BLK=u:

    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock):: &
         bicg_r1, bicg_u1

    real :: rwork(2,7)

    logical:: GoOn, rcmp, xpdt
    integer:: nmv
    real::    alpha, beta, omega, rho0, rho1, sigma
    real::    varrho, hatgamma
    real::    assumedzero, rnrm0, rnrm, rnrmMax0, rnrmMax
    real::    mxnrmx, mxnrmr, kappa0, kappal

    ! External function:

    real :: maxval_abs_BLK

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'proj_bicgstab'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    if(DoTest)write(*,*)'Proj_BiCGSTAB tol,iter:',tol,iter

    info = 0

    if (tol<=0.0) call stop_mpi('Error in Proj_BiCGSTAB: tolerance < 0')
    if (iter<=1)  call stop_mpi('Error in Proj_BiCGSTAB: maxmatvec < 2')

    !
    !     --- Initialize first residual
    !
    assumedzero = 1.e-16
    call eq_BLK(tmp1_BLK,rhs)
    call set_BLK(qx,0.0)

    nmv = 0
    !
    !     --- Initialize iteration loop
    !

    rnrm0 = sqrt( dot_product_BLK(tmp1_BLK,tmp1_BLK))

    rnrm = rnrm0
    if(DoTest) print *,'initial rnrm:',rnrm

    mxnrmx = rnrm0
    mxnrmr = rnrm0
    rcmp = .false.
    xpdt = .false.

    alpha = 0.0
    omega = 1.0
    sigma = 1.0
    rho0 =  1.0
    !
    !     --- Iterate
    !
    select case(typestop)
    case('rel')
       GoOn = rnrm>tol*rnrm0 .and. nmv<iter
       assumedzero = assumedzero*rnrm0
       rnrmMax = 0
       rnrmMax0 = 0

    case('abs')
       GoOn = rnrm>tol       .and. nmv<iter
       assumedzero = assumedzero*rnrm0
       rnrmMax = 0
       rnrmMax0 = 0

    case('max')
       rnrmMax0 = maxval_abs_BLK(nProc,tmp1_BLK)
       rnrmMax  = rnrmMax0
       if(DoTest) print *,'initial rnrmMax:',rnrmMax
       GoOn = rnrmMax>tol    .and. nmv<iter
       assumedzero = assumedzero*rnrmMax
    case default
       call stop_mpi('Error in Proj_BiCGSTAB: unknown typestop value')
    end select

    if (.not.GoOn) then
       if(DoTest) print *,'Proj_BiCGSTAB: nothing to do. info = ',info
       iter = nmv
       info = 3
       RETURN
    end if

    do while (GoOn)
       !
       !     =====================
       !     --- The BiCG part ---
       !     =====================
       !
       rho0 = -omega*rho0

       rho1 = dot_product_BLK(rhs,tmp1_BLK)

       if (abs(rho0)<assumedzero**2) then
          info = 1
          RETURN
       endif
       beta = alpha*(rho1/rho0)
       rho0 = rho1
       call eq_plus_times_BLK(tmp2_BLK,tmp1_BLK,-beta,tmp2_BLK)

       call proj_matvec(tmp2_BLK,bicg_u1)
       nmv = nmv+1

       ! DEBUG
       ! write(*,*)'u =',tmp2_BLK(1:nI,1:nJ,1:nK,iBlockTest)
       ! write(*,*)'u1=',bicg_u1(1:nI,1:nJ,1:nK,iBlockTest)

       sigma=dot_product_BLK(rhs,bicg_u1)

       if (abs(sigma)<assumedzero**2) then
          info = 1
          RETURN
       endif

       alpha = rho1/sigma
       call add_times_BLK(qx,alpha,tmp2_BLK)

       call add_times_BLK(tmp1_BLK,-alpha,bicg_u1)

       call proj_matvec(tmp1_BLK,bicg_r1)
       nmv = nmv+1

       rnrm = sqrt( dot_product_BLK(tmp1_BLK,tmp1_BLK) )

       mxnrmx = max (mxnrmx, rnrm)
       mxnrmr = max (mxnrmr, rnrm)

       ! DEBUG
       ! write(*,*)'rho0, rho1, beta, sigma, alpha, rnrm:',&
       !     rho0, rho1, beta, sigma, alpha, rnrm

       !
       !  ==================================
       !  --- The convex polynomial part ---
       !  ==================================
       !
       !    --- Z = R'R a 2 by 2 matrix
       ! i=1,j=0
       rwork(1,1) = dot_product_BLK(tmp1_BLK,tmp1_BLK)

       ! i=1,j=1
       rwork(2,1) = dot_product_BLK(bicg_r1,tmp1_BLK)
       rwork(1,2) = rwork(2,1)

       ! i=2,j=1
       rwork(2,2) = dot_product_BLK(bicg_r1,bicg_r1)

       !
       !   --- tilde r0 and tilde rl (small vectors)
       !
       rwork(1:2,zz_:zz_+1)   = rwork(1:2,qz_:qz_+1)
       rwork(1,y0_) = -1.0
       rwork(2,y0_) = 0.0

       rwork(1,yl_) = 0.0
       rwork(2,yl_) = -1.0
       !
       !   --- Convex combination
       !
       rwork(1:2,qy_) = rwork(1,yl_)*rwork(1:2,qz_) + &
            rwork(2,yl_)*rwork(1:2,qz_+1)

       kappal = sqrt( sum( rwork(1:2,yl_)*rwork(1:2,qy_) ) )

       rwork(1:2,qy_) = rwork(1,y0_)*rwork(1:2,qz_) + &
            rwork(2,y0_)*rwork(1:2,qz_+1)

       kappa0 = sqrt( sum( rwork(1:2,y0_)*rwork(1:2,qy_) ) )

       varrho = sum( rwork(1:2,yl_)*rwork(1:2,qy_) )
       varrho = varrho / (kappa0*kappal)

       hatgamma = sign(1.0,varrho)*max(abs(varrho),0.7) * (kappa0/kappal)

       rwork(1:2,y0_) = -hatgamma*rwork(1:2,yl_) + rwork(1:2,y0_)

       !
       !    --- Update
       !
       omega = rwork(2,y0_)

       call add_times_BLK(tmp2_BLK,-omega,bicg_u1)

       call add_times_BLK(qx,omega,tmp1_BLK)

       call add_times_BLK(tmp1_BLK,-omega,bicg_r1)

       rwork(1:2,qy_) = rwork(1,y0_)*rwork(1:2,qz_) + &
            rwork(2,y0_)*rwork(1:2,qz_+1)

       rnrm = sqrt( sum( rwork(1:2,y0_)*rwork(1:2,qy_) ) )

       select case(typestop)
       case('rel')
          GoOn = rnrm>tol*rnrm0 .and. nmv<iter
          if(DoTest) print *, nmv,' matvecs, ', ' ||rn||/||r0|| =',rnrm/rnrm0

       case('abs')
          GoOn = rnrm>tol       .and. nmv<iter
          if(DoTest) print *, nmv,' matvecs, ||rn|| =',rnrm

       case('max')
          rnrmMax = maxval_abs_BLK(nProc,tmp1_BLK)
          GoOn = rnrmMax>tol    .and. nmv<iter
          if(DoTest) print *, nmv,' matvecs, max(rn) =',rnrmMax
       end select

    end do
    !
    !     =========================
    !     --- End of iterations ---
    !     =========================

    select case(typestop)
    case('rel')
       if (rnrm>tol*rnrm0) info = 2
       tol = rnrm/rnrm0

    case('abs')
       if (rnrm>tol) info = 2
       tol = rnrm

    case('max')
       if (rnrmMax>tol) info = 2
       tol = rnrmMax
    end select

    if((typestop/='max'.and.rnrm>rnrm0).or.(typestop=='max'.and.rnrmMax&
         >rnrmMax0)) info=-info

    iter = nmv

    call test_stop(NameSub, DoTest)
  end subroutine proj_bicgstab
  !============================================================================
  subroutine set_BLK(qa,qb)

    ! Set qa=qb for all used blocks, where qb is a scalar

    ! Arguments

    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock), intent(out) :: qa
    real, intent(in) :: qb

    ! Local variables:
    integer:: iBlock

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'set_BLK'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    do iBlock=1,nBlock
       if(.not.Unused_B(iBlock)) then
          if(true_BLK(iBlock))then
             qa(1:nI,1:nJ,1:nK,iBlock)=qb
          else
             where(true_cell(1:nI,1:nJ,1:nK,iBlock)) &
                  qa(1:nI,1:nJ,1:nK,iBlock)=qb
          end if
       end if
    end do

    call test_stop(NameSub, DoTest)
  end subroutine set_BLK
  !============================================================================

  subroutine eq_BLK(qa,qb)

    ! Do qa=qb for all used blocks

    ! Arguments

    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock):: qa,qb
    intent(out) :: qa
    intent(in)  :: qb

    ! Local variables:
    integer:: iBlock

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'eq_BLK'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    do iBlock=1,nBlock
       if(.not.Unused_B(iBlock)) then
          if(true_BLK(iBlock))then
             qa(1:nI,1:nJ,1:nK,iBlock)= qb(1:nI,1:nJ,1:nK,iBlock)
          else
             where(true_cell(1:nI,1:nJ,1:nK,iBlock)) &
                  qa(1:nI,1:nJ,1:nK,iBlock)= qb(1:nI,1:nJ,1:nK,iBlock)
          end if

       end if
    end do

    call test_stop(NameSub, DoTest)
  end subroutine eq_BLK
  !============================================================================

  subroutine add_BLK(qa,qb)

    ! Do qa=qa+qb for all used blocks

    ! Arguments

    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock) :: qa, qb
    intent(inout) :: qa
    intent(in)    :: qb

    ! Local variables:
    integer:: iBlock

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'add_BLK'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    do iBlock=1,nBlock
       if(.not.Unused_B(iBlock)) then
          if(true_BLK(iBlock))then
             qa(1:nI,1:nJ,1:nK,iBlock)= &
                  qa(1:nI,1:nJ,1:nK,iBlock)+qb(1:nI,1:nJ,1:nK,iBlock)
          else
             where(true_cell(1:nI,1:nJ,1:nK,iBlock)) &
                  qa(1:nI,1:nJ,1:nK,iBlock)= &
                  qa(1:nI,1:nJ,1:nK,iBlock)+qb(1:nI,1:nJ,1:nK,iBlock)
          end if
       end if
    end do

    call test_stop(NameSub, DoTest)
  end subroutine add_BLK
  !============================================================================

  subroutine sub_BLK(qa,qb)

    ! Do qa=qa-qb for all used blocks

    ! Arguments

    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock) :: qa, qb
    intent(inout) :: qa
    intent(in)    :: qb

    ! Local variables:
    integer:: iBlock

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'sub_BLK'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    do iBlock=1,nBlock
       if(.not.Unused_B(iBlock)) then
          if(true_BLK(iBlock))then
             qa(1:nI,1:nJ,1:nK,iBlock)= &
                  qa(1:nI,1:nJ,1:nK,iBlock)-qb(1:nI,1:nJ,1:nK,iBlock)
          else
             where(true_cell(1:nI,1:nJ,1:nK,iBlock)) &
                  qa(1:nI,1:nJ,1:nK,iBlock)= &
                  qa(1:nI,1:nJ,1:nK,iBlock)-qb(1:nI,1:nJ,1:nK,iBlock)
          end if
       end if
    end do

    call test_stop(NameSub, DoTest)
  end subroutine sub_BLK
  !============================================================================

  subroutine eq_plus_BLK(qa,qb,qc)

    ! Do qa=qb+qc for all used blocks

    ! Arguments

    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock) :: qa,qb,qc
    intent(out) :: qa
    intent(in)  :: qb,qc

    ! Local variables:
    integer:: iBlock

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'eq_plus_BLK'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    do iBlock=1,nBlock
       if(.not.Unused_B(iBlock)) then
          if(true_BLK(iBlock))then
             qa(1:nI,1:nJ,1:nK,iBlock)= &
                  qb(1:nI,1:nJ,1:nK,iBlock)+qc(1:nI,1:nJ,1:nK,iBlock)
          else
             where(true_cell(1:nI,1:nJ,1:nK,iBlock)) &
                  qa(1:nI,1:nJ,1:nK,iBlock)= &
                  qb(1:nI,1:nJ,1:nK,iBlock)+qc(1:nI,1:nJ,1:nK,iBlock)
          end if
       end if
    end do

    call test_stop(NameSub, DoTest)
  end subroutine eq_plus_BLK
  !============================================================================

  subroutine add_times_BLK(qa,qb,qc)

    ! Do qa=qa+qb*qc for all used blocks, where qb is a scalar

    ! Arguments

    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock) :: qa,qc
    intent(inout) :: qa
    intent(in)    :: qc

    real, intent(in) :: qb

    ! Local variables:
    integer:: iBlock

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'add_times_BLK'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    do iBlock=1,nBlock
       if(.not.Unused_B(iBlock)) then
          if(true_BLK(iBlock))then
             qa(1:nI,1:nJ,1:nK,iBlock)= &
                  qa(1:nI,1:nJ,1:nK,iBlock)+qb*qc(1:nI,1:nJ,1:nK,iBlock)
          else
             where(true_cell(1:nI,1:nJ,1:nK,iBlock)) &
                  qa(1:nI,1:nJ,1:nK,iBlock)= &
                  qa(1:nI,1:nJ,1:nK,iBlock)+qb*qc(1:nI,1:nJ,1:nK,iBlock)
          end if
       end if
    end do

    call test_stop(NameSub, DoTest)
  end subroutine add_times_BLK
  !============================================================================

  subroutine eq_plus_times_BLK(qa,qb,qc,qd)

    ! Do qa=qb+qc*qd for all used blocks, where qc is a scalar

    ! Arguments

    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock) :: qa,qb,qd
    intent(inout) :: qa
    intent(in)    :: qb
    intent(inout) :: qd

    real, intent(in) :: qc

    ! Local variables:
    integer:: iBlock

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'eq_plus_times_BLK'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    do iBlock=1,nBlock
       if(.not.Unused_B(iBlock)) then
          if(true_BLK(iBlock))then
             qa(1:nI,1:nJ,1:nK,iBlock)= &
                  qb(1:nI,1:nJ,1:nK,iBlock)+qc*qd(1:nI,1:nJ,1:nK,iBlock)
          else
             where(true_cell(1:nI,1:nJ,1:nK,iBlock)) &
                  qa(1:nI,1:nJ,1:nK,iBlock)= &
                  qb(1:nI,1:nJ,1:nK,iBlock)+qc*qd(1:nI,1:nJ,1:nK,iBlock)
          end if
       end if
    end do

    call test_stop(NameSub, DoTest)
  end subroutine eq_plus_times_BLK
  !============================================================================

  real function dot_product_BLK(qa,qb)

    ! Return qa.qb=sum(qa*qb) for all used blocks

    ! Arguments

    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock), intent(in) :: qa,qb

    ! Local variables:
    real    :: qproduct, qproduct_all
    integer :: iBlock, iError

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'dot_product_BLK'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    qproduct=0.0

    do iBlock=1,nBlock
       if(.not.Unused_B(iBlock)) then
          if(true_BLK(iBlock)) then
             qproduct=qproduct + &
                  sum(qa(1:nI,1:nJ,1:nK,iBlock)*qb(1:nI,1:nJ,1:nK,iBlock))
          else
             qproduct=qproduct + &
                  sum(qa(1:nI,1:nJ,1:nK,iBlock)*qb(1:nI,1:nJ,1:nK,iBlock),&
                  MASK=true_cell(1:nI,1:nJ,1:nK,iBlock))
          end if
       end if
    end do

    if(nProc>1)then
       call MPI_allreduce(qproduct, qproduct_all, 1,  MPI_REAL, MPI_SUM, &
            iComm, iError)
       dot_product_BLK=qproduct_all
       if(DoTest)write(*,*)'me,product,product_all:',&
            iProc,qproduct,qproduct_all
    else
       dot_product_BLK=qproduct
       if(DoTest)write(*,*)'me,qproduct:',iProc,qproduct
    end if

    call test_stop(NameSub, DoTest)
  end function dot_product_BLK
  !============================================================================

  real function sum_BLK(qnum,qa)

    ! Return sum(qa) for all used blocks and true cells
    ! Do for each processor separately if qnum=1, otherwise add them all

    ! Arguments

    integer, intent(in) :: qnum
    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock), intent(in) :: qa

    ! Local variables:
    real    :: qsum, qsum_all
    integer :: iBlock, iError

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'sum_BLK'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    qsum=0.0

    do iBlock=1,nBlock
       if(.not.Unused_B(iBlock)) then
          if(true_BLK(iBlock)) then
             qsum=qsum + sum(qa(1:nI,1:nJ,1:nK,iBlock))
          else
             qsum=qsum + sum(qa(1:nI,1:nJ,1:nK,iBlock), &
                  MASK=true_cell(1:nI,1:nJ,1:nK,iBlock))
          end if
       end if
    end do

    if(qnum>1)then
       call MPI_allreduce(qsum, qsum_all, 1,  MPI_REAL, MPI_SUM, &
            iComm, iError)
       sum_BLK=qsum_all
       if(DoTest)write(*,*)'me,sum,sum_all:',iProc,qsum,qsum_all
    else
       sum_BLK=qsum
       if(DoTest)write(*,*)'me,qsum:',iProc,qsum
    end if

    call test_stop(NameSub, DoTest)
  end function sum_BLK
  !============================================================================

end module ModProjectDivB
!==============================================================================
