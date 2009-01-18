!^CFG COPYRIGHT UM
!^CFG FILE IMPLICIT
subroutine impl_matvec(qx,qy,n)

  ! Calculate qy= L.qx = (I - beta*dt*dR/dw).qx
  ! Do it matrix-free or with the matrix or preconditioned depending 
  ! on "JacobianType"
  !
  ! N contains the number of elements in qx and qy.

  use ModImplicit
  use ModMpi
  implicit none

  integer, intent(in):: n
  real, intent(in)   :: qx(n)
  real, intent(out)  :: qy(n)

  logical :: oktest, oktest_me
  !----------------------------------------------------------------------------

  call set_oktest('impl_matvec',oktest,oktest_me)

  call timing_start('impl_matvec')

  if(oktest_me.and.nImplBLK>0)&
       write(*,*)'matvec_impl initial n,sum(x**2),x(test)=',&
       n,sum(qx(1:n)**2),qx(implVARtest)

  nmatvec=nmatvec+1

  select case(JacobianType)
  case('free')
     call impl_matvec_free(qx,qy,n)
  case('with')
     call impl_matvec_with(qx,qy,n)
  case('prec')
     call impl_matvec_prec(qx,qy,n)
  case default
     call stop_mpi('Unknown value for JacobianType')
  end select

  if(oktest_me.and.nImplBLK>0)&
       write(*,*)'impl_matvec final sum(x**2),sum(y**2),y(test)=',&
       sum(qx(1:n)**2),sum(qy(1:n)**2),qy(implVARtest)

  !DEBUG
  !write(*,*)'implVARtest,maxval,loc(abs(qy))=',&
  !     implVARtest,maxval(abs(qy)),maxloc(abs(qy))

  call timing_stop('impl_matvec')

end subroutine impl_matvec

!=============================================================================
subroutine impl_matvec_free(qx,qy,nn)

  ! Calculate qy=L.qx for the iterative solver, matrix-free 
  ! where L= I - beta*dt*dR/dw   (dt=dt_implicit)
  !
  ! One sided derivative:
  !----------------------
  ! weps = w + eps*qx            ! perturb w
  !
  ! weps'=weps + R(weps,dtexpl)  ! advance weps
  !
  ! dR/dw.qx = (R(w+eps*qx)-R(w))/eps = [(weps'-weps) - (w'-w)]/eps/dtexpl
  !
  !                                   = (weps'-w')/eps/dtexpl - qx/dtexpl
  !
  ! L.qx = dx - beta*dt*dR/dw.qx 
  !      = (1 + beta*dtcoeff)*qx - beta*dtcoeff*(weps' - w')/eps
  !
  ! where w=w_k, w'=w+R_low, beta=ImplCoeff, eps=sqrt(JacobianEps)/||qx||
  ! instead of eps=(JacobianEps)^(1/2)*(w_k.qx)/(qx.qx) suggested by Keyes

  use ModProcMH
  use ModMain, ONLY : Itest,Jtest,Ktest,VARtest
  use ModImplicit
  use ModMpi
  implicit none

  integer, intent(in):: nn
  real, intent(in)   :: qx(nn)
  ! Sometimes this subroutine called with the same array in both arguments
  ! that's why the intent of qy cannot be set to out.
  real, intent(inout):: qy(nn)

  real               :: weps(nI,nJ,nK,nw,MaxImplBLK)
  integer:: n,i,j,k,iw,implBLK, iError
  real:: qeps, qxnrm, qxnrm_total, q1, q2

  logical :: oktest, oktest_me
  !----------------------------------------------------------------------------
  call set_oktest('impl_matvec_free', oktest, oktest_me)

  call timing_start('matvec_free')

  !if(UseSemiImplicit)then
  !   n=0
  !   do implBLK=1,nImplBLK; do k=1,nK; do j=1,nJ; do i=1,nI; do iw=1,nw
  !      n=n+1
  !      weps(i,j,k,iw,implBLK) = qx(n)*wnrm(iw)        
  !   enddo; enddo; enddo; enddo; enddo
  !   ! Advance weps
  !   call get_semi_impl_residual(weps)
  !   n=0
  !   q1=1.+ImplCoeff*dtcoeff
  !   q2=ImplCoeff*dtcoeff
  !   do implBLK=1,nImplBLK; do k=1,nK; do j=1,nJ; do i=1,nI; do iw=1,nw
  !      n=n+1
  !      qy(n) = q1*qx(n) - q2*weps(i,j,k,iw,implBLK)/wnrm(iw)
  !   enddo; enddo; enddo; enddo; enddo
  !
  !   call timing_stop('matvec_free')
  !   RETURN
  !end if

  qxnrm=sum(qx(1:nimpl)**2)
  call MPI_allreduce(qxnrm, qxnrm_total, 1, MPI_REAL, MPI_SUM,iComm,iError)

  if(oktest_me)write(*,*)'impl_matvec_free initial n,sum(x**2)=',nn,qxnrm_total

  qxnrm=sqrt(qxnrm_total/nimpl_total)

  if(qxnrm<smalldouble)qxnrm=1.0

  qeps=sqrt(JacobianEps)/qxnrm

  if(oktest_me)write(*,*)'qeps, qxnrm =',qeps,qxnrm

  n=0
  do implBLK=1,nImplBLK; do k=1,nK; do j=1,nJ; do i=1,nI; do iw=1,nw; 
     n=n+1
     weps(i,j,k,iw,implBLK)=w_k(i,j,k,iw,implBLK)+qeps*qx(n)*wnrm(iw)
     !DEBUG no perturbation
     !weps(i,j,k,iw,implBLK)=w_k(i,j,k,iw,implBLK)
  enddo; enddo; enddo; enddo; enddo

  if(oktest)then
     call MPI_allreduce(sum(weps(:,:,:,:,1:nImplBLK)**2),q1,1,&
          MPI_REAL,MPI_SUM,iComm,iError)
     if(oktest_me)write(*,*)'sum(weps**2)=',q1
     if(oktest_me.and.nImplBLK>0)write(*,*)'weps(test)=',&
          weps(Itest,Jtest,Ktest,VARtest,implBLKtest)
  endif

  ! Advance weps:low order,  no dt, don't subtract
  if(UseSemiImplicit)then
     call get_semi_impl_residual(weps)
  else
     call get_residual(.true.,.false.,.false.,weps,weps) 
  end if

  if(oktest)then
     call MPI_allreduce(sum(weps(:,:,:,:,1:nImplBLK)**2),q1,1,&
          MPI_REAL,MPI_SUM,iComm,iError)
     if(oktest_me)write(*,*)'after advance,sum(weps**2)=',q1
     if(oktest_me.and.nImplBLK>0)write(*,*)'after advance,weps(test)=',&
          weps(Itest,Jtest,Ktest,VARtest,implBLKtest)
  end if

  ! Calculate qy = L.qx = (1 + beta*dtcoeff)*qx - beta*dtcoeff*(weps' - w')/eps
  ! where weps = w + eps*qx, weps' = weps + dt*R(weps) and w' = w + dt*R(w)
  ! qy = qx + beta*dtcoeff*qx - beta*dtcoeff*(w + eps*qx + R(weps) - w - R(w))/eps
  !    = qx - beta*dtcoeff*(R(weps)-R(w))/eps = qx - beta*dt*dR/dU*qx

  q1=1.+ImplCoeff*dtcoeff
  q2=ImplCoeff*dtcoeff/qeps

  if(oktest_me)write(*,*)'dtcoeff,ImplCoeff,q1,q2=',dtcoeff,ImplCoeff,q1,q2

  n=0
  do implBLK=1,nImplBLK; do k=1,nK; do j=1,nJ; do i=1,nI; do iw=1,nw
     n=n+1
     qy(n)=q1*qx(n) - q2*(weps(i,j,k,iw,implBLK) &
          -w_k(i,j,k,iw,implBLK)-RES_impl(i,j,k,iw,implBLK))/wnrm(iw)
  enddo; enddo; enddo; enddo; enddo

  call timing_stop('matvec_free')

  if(oktest)then
     call MPI_allreduce(sum(qy(1:n)**2),q1,1,&
          MPI_REAL,MPI_SUM,iComm,iError)
     if(oktest_me.and.nImplBLK>0)&
          write(*,*)'y,x,weps,w_k,RES_impl(test)=',         &
          qy(implVARtest),qx(implVARtest),             &
          weps(Itest,Jtest,Ktest,VARtest,implBLKtest), &
          w_k(Itest,Jtest,Ktest,VARtest,implBLKtest),  &
          RES_impl(Itest,Jtest,Ktest,VARtest,implBLKtest)
     if(oktest_me)write(*,*)'impl_matvec_free final n,sum(qy**2)=',n,q1
  end if

end subroutine impl_matvec_free

!=============================================================================
subroutine impl_matvec_prec(qx,qy,n)

  ! Calculate qy=P_L.A.P_R.qx for the iterative solver, where 
  ! P_L and P_R are the left and right preconditioner matrices,
  ! L = I - beta*dt*dR/dw, and R is the residual from dw/dt = R(w).
  !
  ! The multiplication by L is done in a matrix free fashion.

  use ModImplicit
  use ModLinearSolver, ONLY: Lhepta, Uhepta
  implicit none

  integer, intent(in):: n
  real, intent(in)   :: qx(n)
  real, intent(out)  :: qy(n)

  integer :: implBLK

  logical :: oktest, oktest_me
  !-----------------------------------------------------------------------------

  call set_oktest('impl_matvec_prec',oktest,oktest_me)

  if(oktest_me)write(*,*)'impl_matvec_prec initial n,sum(x**2),sum(y**2)=',&
       nimpl,sum(qx(1:nimpl)**2),sum(qy(1:nimpl)**2)

  qy=qx

  ! qy = P_R.qx, where P_R = I, U^{-1}, or U^{-1}L^{-1}
  ! for left, symmetric and right preconditioning, respectively

  if(PrecondSide/='left')then
     do implBLK=1,nImplBLK
        if(PrecondSide=='right') &
             call Lhepta(nIJK,nw,nI,nI*nJ,&
             qy(nwIJK*(implBLK-1)+1) ,&
             MAT(1,1,1,1,1,1,implBLK),&   ! Main diagonal
             MAT(1,1,1,1,1,2,implBLK),&   ! -i
             MAT(1,1,1,1,1,4,implBLK),&   ! -j
             MAT(1,1,1,1,1,6,implBLK))    ! -k
        call Uhepta(.true.,nIJK,nw,nI,nI*nJ,&
             qy(nwIJK*(implBLK-1)+1) ,  &
             MAT(1,1,1,1,1,3,implBLK),  &   ! +i diagonal
             MAT(1,1,1,1,1,5,implBLK),  &   ! +j 
             MAT(1,1,1,1,1,7,implBLK))      ! +k
     end do
  end if

  call impl_matvec_free(qy,qy,n) ! qy = A.qy

  ! qy = P_L.qy, where P_L==U^{-1}.L^{-1}, L^{-1}, or I
  ! for left, symmetric, and right preconditioning, respectively

  if(PrecondSide/='right')then
     do implBLK=1,nImplBLK
        call Lhepta(nIJK,nw,nI,nI*nJ,&
             qy(nwIJK*(implBLK-1)+1) ,&
             MAT(1,1,1,1,1,1,implBLK),&   ! Main diagonal
             MAT(1,1,1,1,1,2,implBLK),&   ! -i
             MAT(1,1,1,1,1,4,implBLK),&   ! -j
             MAT(1,1,1,1,1,6,implBLK))    ! -k
        if(PrecondSide=='left') &
             call Uhepta(.true.,nIJK,nw,nI,nI*nJ,&
             qy(nwIJK*(implBLK-1)+1),   &
             MAT(1,1,1,1,1,3,implBLK),  &   ! +i diagonal
             MAT(1,1,1,1,1,5,implBLK),  &   ! +j
             MAT(1,1,1,1,1,7,implBLK))      ! +k
     end do
  end if

  if(oktest_me)write(*,*)'impl_matvec_prec final n,sum(x**2),sum(y**2)=',&
       n,sum(qx**2),sum(qy**2)

end subroutine impl_matvec_prec

!=============================================================================

subroutine impl_matvec_with(qx,qy,n)

  ! Calculate qy = L.qx, where L is stored in the block hepta diagonal 
  ! matrix MAT

  use ModImplicit
  implicit none

  integer, intent(in) :: n
  real, intent(in)    :: qx(n)
  real, intent(out)   :: qy(n)

  integer :: implBLK
  logical :: oktest, oktest_me
  !---------------------------------------------------------------------------

  call set_oktest('impl_matvec_with',oktest,oktest_me)

  if(oktest_me)write(*,*)'impl_matvec_with initial n,sum(x**2),sum(y**2)=',&
       n,sum(qx**2),sum(qy**2)

  qy = 0.

  do implBLK=1,nImplBLK
     call impl_matvec_with_BLK(MAT(1,1,1,1,1,1,implBLK),&
          qx(1+(implBLK-1)*nwIJK),&
          qy(1+(implBLK-1)*nwIJK))

  end do

  if(oktest_me)write(*,*)'impl_matvec_with final n,sum(x**2),sum(y**2)=',&
       n,sum(qx**2),sum(qy**2)

contains

  subroutine impl_matvec_with_BLK(JAC,qx,qy)

    ! JAC.qx = qy  for a single block

    implicit none

    real, intent(in) :: qx(nw,nIJK)
    real, intent(inout):: qy(nw,nIJK)
    real, intent(in) :: JAC(nw,nw,nIJK,nstencil)
    integer :: i,j,k

    do j=1,nIJK
       do i=1,nw
          do k=1,nw
             qy(i,j)=qy(i,j)+JAC(i,k,j,1)*qx(k,j)
          end do
       end do
       if(j>1)then
          do i=1,nw
             do k=1,nw
                qy(i,j)=qy(i,j)+JAC(i,k,j,2)*qx(k,j-1)
             end do
          end do
       end if
       if(j<nIJK)then
          do i=1,nw
             do k=1,nw
                qy(i,j)=qy(i,j)+JAC(i,k,j,3)*qx(k,j+1)
             end do
          end do
       end if
       if(j>nI)then
          do i=1,nw
             do k=1,nw
                qy(i,j)=qy(i,j)+JAC(i,k,j,4)*qx(k,j-nI)
             end do
          end do
       end if
       if(j<=nIJK-nI)then
          do i=1,nw
             do k=1,nw
                qy(i,j)=qy(i,j)+JAC(i,k,j,5)*qx(k,j+nI)
             end do
          end do
       end if
       if(j>nI*nJ)then
          do i=1,nw
             do k=1,nw
                qy(i,j)=qy(i,j)+JAC(i,k,j,6)*qx(k,j-nI*nJ)
             end do
          end do
       end if
       if(j<=nIJK-nI*nJ)then
          do i=1,nw
             do k=1,nw
                qy(i,j)=qy(i,j)+JAC(i,k,j,7)*qx(k,j+nI*nJ)
             end do
          end do
       end if
    end do
  end subroutine impl_matvec_with_BLK

end subroutine impl_matvec_with
